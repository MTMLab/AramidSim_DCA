#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import subprocess
from pathlib import Path

def run(cmd: str):
    try:
        return subprocess.run(cmd, shell=True)
    except Exception:
        return None

def main():
    namelist_path = Path("name.txt")
    if not namelist_path.exists():
        print("Unable to open file")
        return

    names = [line.rstrip("\n") for line in namelist_path.open("r") if line.strip()]

    ikik = 0
    with Path("name_list.txt").open("w") as name_list_out:
        for name in names:
            Path("test_back.smi").write_text(name + "\n")
            run("cp test_back.smi test.smi")
            run("taskset --cpu-list 0 obabel -ismi test.smi -omol2 --gen3d --partialcharge eem -O monomer.mol2")
            run("python ../Util/Util_monomer_matching.py")
            run("taskset --cpu-list 0 obabel -imol2 monomer_reorder.mol2 -osmi -O test2.smi")
            system1, system2 = "", ""
            if Path("test2.smi").exists():
                tokens = Path("test2.smi").read_text().split()
                if len(tokens) >= 1:
                    system1 = tokens[0]
                if len(tokens) >= 2:
                    system2 = tokens[1]

            # Replace patterns (-COOH modifications)
            result_smiles = system1
            diamine_substitute = "(NC8=CC=C(N)C=C8)"

            # --- Step 1: Replace -COOH oxygen with diamine substitute ---
            cooh_replaced = False

            # Pattern 1: C(=O)O -> C(=O)[diamine]
            pos = result_smiles.find("C(=O)O")
            if pos != -1:
                result_smiles = (
                    result_smiles[:pos + 5]
                    + diamine_substitute
                    + result_smiles[pos + 5 + 1 :]
                )
                cooh_replaced = True

            # Pattern 2: C(O)=O -> C([diamine])=O
            if not cooh_replaced:
                pos = result_smiles.find("C(O)=O")
                if pos != -1:
                    result_smiles = (
                        result_smiles[:pos + 2]
                        + diamine_substitute
                        + result_smiles[pos + 2 + 1 :]
                    )
                    cooh_replaced = True

            # Pattern 3: OC(=O) -> [diamine]C(=O)
            if not cooh_replaced:
                pos = result_smiles.find("OC(=O)")
                if pos != -1:
                    result_smiles = (
                        result_smiles[:pos]
                        + diamine_substitute
                        + result_smiles[pos + 1 :]
                    )
                    cooh_replaced = True

            # Pattern 4: O=C(O) -> O=C([diamine])
            if not cooh_replaced:
                pos = result_smiles.find("O=C(O)")
                if pos != -1:
                    result_smiles = (
                        result_smiles[:pos + 4]
                        + diamine_substitute
                        + result_smiles[pos + 4 + 1 :]
                    )
                    cooh_replaced = True

            # Pattern 5: C(=O)(O) -> C(=O)([diamine])
            if not cooh_replaced:
                pos = result_smiles.find("C(=O)(O)")
                if pos != -1:
                    result_smiles = (
                        result_smiles[:pos + 6]
                        + diamine_substitute
                        + result_smiles[pos + 6 + 1 :]
                    )
                    cooh_replaced = True

            # --- Step 2: Replace remaining -COOH oxygen with Cl ---
            cooh_cl_replaced = False

            # Pattern 1: C(=O)O -> C(=O)Cl
            pos = result_smiles.find("C(=O)O")
            if pos != -1:
                result_smiles = (
                    result_smiles[:pos + 5]
                    + "Cl"
                    + result_smiles[pos + 5 + 1 :]
                )
                cooh_cl_replaced = True

            # Pattern 2: C(O)=O -> C(Cl)=O
            if not cooh_cl_replaced:
                pos = result_smiles.find("C(O)=O")
                if pos != -1:
                    result_smiles = (
                        result_smiles[:pos + 2]
                        + "Cl"
                        + result_smiles[pos + 2 + 1 :]
                    )
                    cooh_cl_replaced = True

            # Pattern 3: OC(=O) -> ClC(=O)
            if not cooh_cl_replaced:
                pos = result_smiles.find("OC(=O)")
                if pos != -1:
                    result_smiles = (
                        result_smiles[:pos]
                        + "Cl"
                        + result_smiles[pos + 1 :]
                    )
                    cooh_cl_replaced = True

            # Pattern 4: O=C(O) -> O=C(Cl)
            if not cooh_cl_replaced:
                pos = result_smiles.find("O=C(O)")
                if pos != -1:
                    result_smiles = (
                        result_smiles[:pos + 4]
                        + "Cl"
                        + result_smiles[pos + 4 + 1 :]
                    )
                    cooh_cl_replaced = True

            # Pattern 5: C(=O)(O) -> C(=O)(Cl)
            if not cooh_cl_replaced:
                pos = result_smiles.find("C(=O)(O)")
                if pos != -1:
                    result_smiles = (
                        result_smiles[:pos + 6]
                        + "Cl"
                        + result_smiles[pos + 6 + 1 :]
                    )
                    cooh_cl_replaced = True

            # Write final SMILES
            Path("monomer.smi").write_text(result_smiles + "\n")

            # Copy to structures/monomer_<ikik>.smi
            Path("structures").mkdir(exist_ok=True, parents=True)
            run(f"cp monomer.smi structures/monomer_{ikik}.smi")

            # Append to name_list.txt
            name_list_out.write(f"monomer_{ikik}\n")

            ikik += 1

if __name__ == "__main__":
    main()

