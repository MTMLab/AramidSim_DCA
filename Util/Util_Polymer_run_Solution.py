#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import subprocess

def run_cmd(cmd: str):
    ret = os.system(cmd)
    if ret != 0:
        sys.stderr.write(f"Command failed: {cmd}\n")
        sys.exit(1)

def main():
    # Read names
    if not os.path.exists("name.txt"):
        sys.stderr.write("Unable to open file: name.txt\n")
        sys.exit(1)

    with open("name.txt") as f:
        names = [line.strip() for line in f if line.strip()]

    # Backup E_h_bond.txt
    if os.path.exists("E_h_bond.txt"):
        shutil.move("E_h_bond.txt", "E_h_bond_back.txt")

    # Open result file
    with open("E_h_bond.txt", "a") as out_summary:
        for name in names:
            # Generate run script
            with open("run", "w") as runf:
                runf.write("#!/bin/bash\n")
                runf.write("rm -rf lammps/*\n")
                runf.write(f"mkdir lammps/{name}\n")
                runf.write(f"cp back/* lammps/{name}\n")

                runf.write("if ! command -v conda &> /dev/null; then\n")
                runf.write("    echo \"❌ Conda not found. Please install Anaconda or Miniconda.\"\n")
                runf.write("    exit 1\n")
                runf.write("fi\n")

                runf.write("source \"$(conda info --base)/etc/profile.d/conda.sh\"\n")
                runf.write("conda activate AmberTools23 || {\n")
                runf.write("    echo \"❌ AmberTools23 environment not found.\"\n")
                runf.write("    exit 1\n}\n")

                # Polymer setup
                runf.write(f"cp structures/{name}.smi test.smi\n")
                runf.write("taskset --cpu-list 0 obabel -ismi test.smi -omol2 --gen3d --partialcharge eem -O monomer_com1.mol2\n")
                runf.write("antechamber -i monomer_com1.mol2 -fi mol2 -fo mol2 -o monomer_com2.mol2 -pf y -at gaff\n")
                runf.write("python ../Util/Util_monomer_reorder.py\n")
                runf.write("python ../Util/Util_make_lt_linear.py\n")
                runf.write("cp monomer_reorder2.mol2 mol2tolt/test/monomer.mol2\n")
                runf.write("cp polymer_new.lt moltemplates/polymer.lt\n")
                runf.write("cp system_new.lt moltemplates/system.lt\n")
                runf.write("cd mol2tolt/ && ./run.sh && cp test/monomer.lt ../moltemplates/monomer_add.lt && cd ..\n")
                runf.write("cd moltemplates && ../../Util/moltemplate-master/moltemplate/scripts/moltemplate.sh system.lt > moltemplate.log 2>&1\n")
                runf.write(f"cp system.data ../lammps/{name}/system.data\ncd ..\n")

                # Polymer mass
                runf.write("cd moltemplates && cp system.data ../ratio_calculation/polymer/system.data && cd ../ratio_calculation/polymer\n")
                runf.write("../../../Util/lammps-2Aug2023/src/lmp_serial -i run_polymer_mass.npt2\n")
                runf.write("python polymer_mass.py\n")
                runf.write("cp extracted_mass.txt ../extracted_mass.txt\ncd ../..\n")

                # Solvent setup
                runf.write("cp solvent/solvent.smi solvent.smi\n")
                runf.write("taskset --cpu-list 0 obabel -ismi solvent.smi -omol2 --gen3d --partialcharge eem -O solvent_com1.mol2\n")
                runf.write("antechamber -i solvent_com1.mol2 -fi mol2 -fo mol2 -o solvent_com2.mol2 -pf y -at gaff\n")
                runf.write("cp solvent_com2.mol2 mol2tolt_solvent/test/solvent.mol2\n")
                runf.write("cd mol2tolt_solvent/ && ./run.sh && cp test/solvent.lt ../moltemplates_solvent_single/solvent.lt && cd ..\n")
                runf.write("cd moltemplates_solvent_single && ../../Util/moltemplate-master/moltemplate/scripts/moltemplate.sh system.lt > moltemplate.log 2>&1\n")
                runf.write("cp system.data ../ratio_calculation/solvent/system.data\n")
                runf.write("cd ../ratio_calculation/solvent && ../../../Util/lammps-2Aug2023/src/lmp_serial -in run_solvent_mass.npt2\n")
                runf.write("python solvent_mass.py\n")
                runf.write("cp extracted_solvent_mass.txt ../extracted_solvent_mass.txt\ncd ..\n")

                # Solvent system assembly
                runf.write("python number_solvnet.py\npython factoring.py\npython solvent_lt.py\ncp system_solvent.lt ../system_solvent.lt\ncd ..\n")
                runf.write("cp system_solvent.lt moltemplates_solvent/system.lt\n")
                runf.write("cd mol2tolt_solvent/ && ./run.sh && cp test/solvent.lt ../moltemplates_solvent/solvent.lt && cd ..\n")
                runf.write("cd moltemplates_solvent && ../../Util/moltemplate-master/moltemplate/scripts/moltemplate.sh system.lt > moltemplate.log 2>&1\n")
                runf.write(f"cp system.data ../lammps/{name}/system_solvent.data\ncd ..\n")

                # Run LAMMPS
                runf.write(f"cd lammps/{name}\n")
                runf.write("python ../../group_polymer.py\npython ../../group_solvent.py\n")
                runf.write("../../../Util/lammps-2Aug2023/src/lmp_serial -in run_iso.in.npt2_wo_strain\n")
                runf.write("../../../Util/lammps-2Aug2023/src/lmp_serial -i run_iso.in.npt2_wo_strain_pppm\n")
                runf.write("python ../../../Util/Util_Polymer_Output_Solution.py\n")

            # Copy run → run_exe, chmod, execute
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
                    EE, HH, Pi, HE = map(float, vals[:5])
            except Exception as e:
                sys.stderr.write(f"Error reading output_all.txt: {e}\n")
                sys.exit(1)

            # Append results
            out_summary.write(f"{name}    {EE}    {HH}    {Pi}    {HE}\n")

            # Cleanup
            os.remove("output_all.txt")

if __name__ == "__main__":
    main()

