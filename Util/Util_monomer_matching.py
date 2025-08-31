#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

def parse_mol2_atoms(filename):
    atoms = []
    bonds = []
    with open(filename, "r") as f:
        lines = f.readlines()

    atom_start = lines.index("@<TRIPOS>ATOM\n") + 1
    bond_start = lines.index("@<TRIPOS>BOND\n")

    for line in lines[atom_start:bond_start]:
        parts = line.split()
        if len(parts) >= 9:
            atom_id = int(parts[0])
            atom_name = parts[1]
            x, y, z = map(float, parts[2:5])
            atom_type = parts[5]
            charge = float(parts[-1])
            atoms.append({
                "id": atom_id,
                "name": atom_name,
                "x": x, "y": y, "z": z,
                "type": atom_type,
                "charge": charge
            })

    for line in lines[bond_start+1:]:
        parts = line.split()
        if len(parts) >= 4:
            bond_id = int(parts[0])
            a1, a2 = int(parts[1]), int(parts[2])
            btype = parts[3]
            bonds.append({"id": bond_id, "a1": a1, "a2": a2, "type": btype})

    return atoms, bonds, lines[:atom_start], lines[bond_start:bond_start+1]

def find_cooh_carbon(atoms, bonds):
    # Map atom IDs to element guess (by type/name)
    elem = {}
    for atom in atoms:
        name = atom["name"].upper()
        if name.startswith("C"):
            elem[atom["id"]] = "C"
        elif name.startswith("O"):
            elem[atom["id"]] = "O"
        elif name.startswith("H"):
            elem[atom["id"]] = "H"
        else:
            elem[atom["id"]] = "X"

    for atom in atoms:
        if elem[atom["id"]] != "C":
            continue
        has_O = False
        has_OH = False
        for b in bonds:
            if b["a1"] == atom["id"]:
                other = b["a2"]
            elif b["a2"] == atom["id"]:
                other = b["a1"]
            else:
                continue
            if elem.get(other) == "O":
                # check if O bonded to H
                bonded_H = any(
                    (bb["a1"] == other and elem.get(bb["a2"]) == "H") or
                    (bb["a2"] == other and elem.get(bb["a1"]) == "H")
                    for bb in bonds
                )
                if bonded_H:
                    has_OH = True
                else:
                    has_O = True
        if has_O and has_OH:
            return atom["id"]
    return None

def reorder_atoms(atoms, bonds, cooh_id):
    # Place COOH carbon first
    atoms_sorted = [a for a in atoms if a["id"] == cooh_id] + [a for a in atoms if a["id"] != cooh_id]

    # Remap IDs
    id_map = {old["id"]: new_id+1 for new_id, old in enumerate(atoms_sorted)}
    for i, atom in enumerate(atoms_sorted):
        atom["new_id"] = i+1

    # Update bond indices
    for b in bonds:
        b["a1"] = id_map[b["a1"]]
        b["a2"] = id_map[b["a2"]]

    return atoms_sorted, bonds

def write_mol2(filename, header_lines, atoms, bonds, bond_header):
    with open(filename, "w") as f:
        f.writelines(header_lines)
        for atom in atoms:
            f.write(f"{atom['new_id']} {atom['name']} {atom['x']:.4f} {atom['y']:.4f} {atom['z']:.4f} {atom['type']} 1 UNL1 {atom['charge']:.4f}\n")
        f.writelines(bond_header)
        for b in bonds:
            f.write(f"{b['id']} {b['a1']} {b['a2']} {b['type']}\n")

def main():
    atoms, bonds, header_lines, bond_header = parse_mol2_atoms("monomer.mol2")
    cooh_id = find_cooh_carbon(atoms, bonds)
    if cooh_id is None:
        sys.stderr.write("No -COOH groups found in the molecule!\n")
        sys.exit(1)
    atoms_sorted, bonds_updated = reorder_atoms(atoms, bonds, cooh_id)
    write_mol2("monomer_reorder.mol2", header_lines, atoms_sorted, bonds_updated, bond_header)

if __name__ == "__main__":
    main()

