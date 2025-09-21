#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Deterministic screening script for DCA monomers (Windows-safe, pass-only output).

Run with NO arguments: `python screening_monomers.py`
Behavior equals:
  --input DCA_data.csv \
  --smiles_col isosmiles \
  --complexity_col complexity \
  --complexity_max 800 \
  --no-isomeric \
  --angle_min 90.0 \
  --max_perp 6.0
Output: only the passed rows to `DCA_screening_data.csv` with the SAME columns/order as the input.
Console: progress percentage ON by default (single updating line). No RDKit noise. Use `--no-progress` to disable.
"""

import argparse
import hashlib
import math
import struct
import sys
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
from rdkit import RDLogger

# Silence RDKit logs across versions
for ch in ("rdApp.info", "rdApp.warning", "rdApp.error", "rdApp.debug"):
    try:
        RDLogger.DisableLog(ch)
    except Exception:
        pass

# -----------------------------
# Config
# -----------------------------
EPS = 1e-9
ROUND_N = 6

@dataclass
class GeoResult:
    ok: bool
    theta_deg: Optional[float]
    max_perp: Optional[float]
    reason: str = ""

# -----------------------------
# Utilities
# -----------------------------

def canonicalize(smiles: str, isomeric: bool = False) -> Optional[str]:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None
        return Chem.MolToSmiles(mol, isomericSmiles=isomeric, canonical=True)
    except Exception:
        return None


def smiles_seed(cano: str) -> int:
    u32 = struct.unpack("<I", hashlib.sha256(cano.encode()).digest()[:4])[0]
    seed = int(u32 & 0x7FFFFFFF)
    if seed == 0:
        seed = 1
    return seed


def get_formal_charge_zero(mol: Chem.Mol) -> bool:
    try:
        chg = int(rdmolops.GetFormalCharge(mol))
    except Exception:
        chg = sum(int(a.GetFormalCharge()) for a in mol.GetAtoms())
    return chg == 0


def collect_carboxyl_c_atoms(mol: Chem.Mol) -> List[int]:
    patt = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    out: List[int] = []
    for m in mol.GetSubstructMatches(patt):
        c_idx = m[0]
        o_idx = m[2]
        o = mol.GetAtomWithIdx(o_idx)
        if o.GetTotalNumHs() >= 1:
            out.append(c_idx)
    return sorted(set(out))


def embed_deterministic_conf(mol: Chem.Mol, seed: int) -> Optional[Chem.Conformer]:
    molH = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3()
    try:
        ps.SetRandomSeed(int(seed))
    except Exception:
        ps.randomSeed = int(seed)
    ps.useSmallRingTorsions = True
    ps.useRandomCoords = False
    cid = AllChem.EmbedMolecule(molH, ps)
    if cid < 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(molH)
    except Exception:
        pass
    return molH.GetConformer(int(cid))


def pick_endpoints_topological(mol: Chem.Mol, carbox_cidx: List[int]) -> Optional[Tuple[int, int]]:
    if len(carbox_cidx) < 2:
        return None
    D = rdmolops.GetDistanceMatrix(mol)
    best = None
    best_val = -1.0
    for ii in range(len(carbox_cidx)):
        for jj in range(ii + 1, len(carbox_cidx)):
            a = carbox_cidx[ii]
            b = carbox_cidx[jj]
            d = float(D[a, b])
            if d > best_val + EPS or (abs(d - best_val) <= EPS and (a, b) < (best or (999999, 999999))):
                best_val = d
                best = (a, b)
    return best


def center_of_mass(conf: Chem.Conformer, mol: Chem.Mol) -> np.ndarray:
    coords = []
    masses = []
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        coords.append([p.x, p.y, p.z])
        masses.append(a.GetMass())
    coords_np = np.asarray(coords, dtype=float)
    masses_np = np.asarray(masses, dtype=float)
    w = masses_np / (masses_np.sum() + EPS)
    return (coords_np * w[:, None]).sum(axis=0)


def angle_at_com(conf: Chem.Conformer, mol: Chem.Mol, i: int, j: int) -> float:
    ci = conf.GetAtomPosition(i)
    cj = conf.GetAtomPosition(j)
    com = center_of_mass(conf, mol)
    vi = np.array([ci.x, ci.y, ci.z], dtype=float) - com
    vj = np.array([cj.x, cj.y, cj.z], dtype=float) - com
    ni = np.linalg.norm(vi) + EPS
    nj = np.linalg.norm(vj) + EPS
    cosang = float(np.clip(np.dot(vi, vj) / (ni * nj), -1.0, 1.0))
    deg = math.degrees(math.acos(cosang))
    return float(round(deg, ROUND_N))


def max_perp_thickness(conf: Chem.Conformer, mol: Chem.Mol, i: int, j: int) -> float:
    ci = conf.GetAtomPosition(i)
    cj = conf.GetAtomPosition(j)
    p0 = np.array([ci.x, ci.y, ci.z], dtype=float)
    p1 = np.array([cj.x, cj.y, cj.z], dtype=float)
    v = p1 - p0
    nv2 = max(float(np.dot(v, v)), EPS)
    maxd = 0.0
    for a in mol.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        r = np.array([p.x, p.y, p.z], dtype=float)
        t = float(np.dot(r - p0, v) / nv2)
        proj = p0 + t * v
        d = float(np.linalg.norm(r - proj))
        if d > maxd:
            maxd = d
    return float(round(maxd, ROUND_N))


def geometry_screen(smiles: str, angle_min: float, max_perp_limit: float, isomeric: bool) -> GeoResult:
    cano = canonicalize(smiles, isomeric=isomeric)
    if not cano:
        return GeoResult(False, None, None, "invalid_smiles")
    mol = Chem.MolFromSmiles(cano)
    if not mol:
        return GeoResult(False, None, None, "mol_parse_fail")
    if not get_formal_charge_zero(mol):
        return GeoResult(False, None, None, "non_neutral")
    carbox_c = collect_carboxyl_c_atoms(mol)
    if len(carbox_c) < 2:
        return GeoResult(False, None, None, "need_two_carboxyl")
    pair = pick_endpoints_topological(mol, carbox_c)
    if not pair:
        return GeoResult(False, None, None, "endpoint_pick_fail")
    i, j = pair
    seed = smiles_seed(cano)
    conf = embed_deterministic_conf(mol, seed)
    if conf is None:
        return GeoResult(False, None, None, "embed_fail")
    theta = angle_at_com(conf, mol, i, j)
    thickness = max_perp_thickness(conf, mol, i, j)
    angle_ok = (theta + EPS) >= angle_min
    thickness_ok = (thickness + EPS) <= max_perp_limit
    ok = angle_ok and thickness_ok
    reason = ""
    if not ok:
        reason = "angle_lt_min" if not angle_ok else "thickness_gt_max"
    return GeoResult(ok, theta, thickness, reason)


# -----------------------------
# Main
# -----------------------------

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="DCA_data.csv")
    ap.add_argument("--smiles_col", default="isosmiles")
    ap.add_argument("--complexity_col", default="complexity")
    ap.add_argument("--complexity_max", type=float, default=800.0)
    ap.add_argument("--isomeric", dest="isomeric", action="store_true")
    ap.add_argument("--no-isomeric", dest="isomeric", action="store_false")
    ap.set_defaults(isomeric=False)
    ap.add_argument("--angle_min", type=float, default=90.0)
    ap.add_argument("--max_perp", type=float, default=6.0)
    ap.add_argument("--out", default="DCA_screening_data.csv")
    ap.add_argument("--progress", dest="progress", action="store_true")
    ap.add_argument("--no-progress", dest="progress", action="store_false")
    ap.set_defaults(progress=True)
    args = ap.parse_args()

    try:
        df = pd.read_csv(args.input, low_memory=False)
    except FileNotFoundError:
        print("ERROR: input CSV '{}' not found.".format(args.input), file=sys.stderr)
        sys.exit(2)

    if args.smiles_col not in df.columns:
        print("ERROR: SMILES column '{}' not found in input.".format(args.smiles_col), file=sys.stderr)
        sys.exit(2)

    if args.complexity_col in df.columns:
        df[args.complexity_col] = pd.to_numeric(df[args.complexity_col], errors='coerce')
    else:
        df[args.complexity_col] = np.nan

    passed_idx: List[int] = []

    total = len(df)
    last_pct = -1

    for i, (idx, row) in enumerate(df.iterrows(), start=1):
        smi = str(row[args.smiles_col]) if pd.notna(row[args.smiles_col]) else None
        if not smi or smi.strip() == "":
            if args.progress and total > 0:
                pct = int(i * 100 / total)
                if pct != last_pct:
                    sys.stdout.write("\rProgress: {}% ({}/{})".format(pct, i, total))
                    sys.stdout.flush()
                    last_pct = pct
            continue

        cval = row[args.complexity_col]
        if pd.isna(cval) or float(cval) > float(args.complexity_max) + EPS:
            if args.progress and total > 0:
                pct = int(i * 100 / total)
                if pct != last_pct:
                    sys.stdout.write("\rProgress: {}% ({}/{})".format(pct, i, total))
                    sys.stdout.flush()
                    last_pct = pct
            continue

        g = geometry_screen(smi, args.angle_min, args.max_perp, args.isomeric)
        if g.ok:
            passed_idx.append(idx)

        if args.progress and total > 0:
            pct = int(i * 100 / total)
            if pct != last_pct:
                sys.stdout.write("\rProgress: {}% ({}/{})".format(pct, i, total))
                sys.stdout.flush()
                last_pct = pct

    if args.progress:
        sys.stdout.write("\r" + (" " * 60) + "\r")
        sys.stdout.flush()

    df_pass = df.loc[passed_idx]
    df_pass.to_csv(args.out, index=False)

    print("Saved {} passed molecules to {}".format(len(df_pass), args.out))


if __name__ == "__main__":
    main()
