#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import shutil

def run_cmd(cmd: str):
    ret = os.system(cmd)
    if ret != 0:
        sys.stderr.write(f"Command failed: {cmd}\n")
        sys.exit(1)

def main():
    # Read monomer names
    if not os.path.exists("name.txt"):
        sys.stderr.write("Unable to open file: name.txt\n")
        sys.exit(1)

    with open("name.txt") as f:
        names = [line.strip() for line in f if line.strip()]

    # Backup E_h_bond.txt
    if os.path.exists("E_h_bond.txt"):
        shutil.move("E_h_bond.txt", "E_h_bond_back.txt")

    # Open summary output
    with open("E_h_bond.txt", "a") as out_summary:
        for name in names:
            # Generate run script
            with open("run", "w") as runf:
                runf.write("#!/bin/bash\n")
                runf.write("rm -rf lammps/*\n")
                runf.write(f"mkdir lammps/{name}\n")
                runf.write(f"cp back_S/* lammps/{name}\n")

                runf.write("if ! command -v conda &> /dev/null; then\n")
                runf.write("    echo \"❌ Conda not found. Please install Anaconda or Miniconda.\"\n")
                runf.write("    exit 1\nfi\n")

                runf.write("source \"$(conda info --base)/etc/profile.d/conda.sh\"\n")
                runf.write("conda activate AmberTools23 || {\n")
                runf.write("    echo \"❌ AmberTools23 environment not found.\"\n")
                runf.write("    exit 1\n}\n")

                # Monomer setup
                runf.write(f"cp structures/{name}.smi test.smi\n")
                runf.write("taskset --cpu-list 0 obabel -ismi test.smi -omol2 --gen3d --partialcharge eem -O monomer_com1.mol2\n")
                runf.write("antechamber -i monomer_com1.mol2 -fi mol2 -fo mol2 -o monomer_com2.mol2 -pf y -at gaff\n")
                runf.write("python ../Util/Util_monomer_reorder.py\n")
                runf.write("python ../Util/Util_make_lt_fiber.py\n")
                runf.write("cp monomer_reorder2.mol2 mol2tolt/test/monomer.mol2\n")

                runf.write("cp polymer_new.lt moltemplates/polymer.lt\n")
                runf.write("cp system_new.lt moltemplates/system.lt\n")

                runf.write("cd mol2tolt/ && ./run.sh && cp test/monomer.lt ../moltemplates/monomer_add.lt && cd ..\n")

                runf.write("cd moltemplates && ../../Util/moltemplate-master/moltemplate/scripts/moltemplate.sh system.lt > moltemplate.log 2>&1\n")
                runf.write("python ../../Util/Util_data_mol_modify.py\n")
                runf.write(f"cp system2.data ../lammps/{name}/system.data\ncd ..\n")

                # LAMMPS execution
                runf.write(f"cd lammps/{name}\n")
                runf.write("../../../Util/lammps-2Aug2023/src/lmp_serial -i run.in.npt2\n")
                runf.write("../../../Util/lammps-2Aug2023/src/lmp_serial -i run.in.npt2_pppm\n")
                runf.write("python ../../../Util/Util_Polymer_Output_Stretched.py\n")

            # Prepare and run script
            shutil.copy("run", "run_exe")
            run_cmd("chmod 770 run_exe")
            run_cmd("./run_exe")

            # Copy outputs
            run_cmd("cp output_all_back.txt output_all.txt")
            run_cmd(f"cp lammps/{name}/output_all.txt ./")

            # Read results
            try:
                with open("output_all.txt") as f:
                    vals = f.read().split()
                    EE, HH, Pi, HE = map(float, vals[:4])
            except Exception as e:
                sys.stderr.write(f"Error reading output_all.txt: {e}\n")
                sys.exit(1)

            # Append to summary
            out_summary.write(f"{name}    {EE}    {HH}    {Pi}    {HE}\n")

            # Cleanup
            os.remove("output_all.txt")

if __name__ == "__main__":
    main()

