#!/usr/bin/env python3
import sys
import os

def main():
    # Constants
    namax = 300000
    
    # Initialize variables
    elem = [-1] * namax
    n_improper0 = 0
    n_improper_type0 = 0
    
    # Initialize lists and arrays
    nn = [[-1 for _ in range(4)] for _ in range(namax)]
    b_id0 = [[-1 for _ in range(4)] for _ in range(namax)]
    a_id0 = [[-1 for _ in range(5)] for _ in range(namax)]
    b_id1 = [[-1 for _ in range(4)] for _ in range(namax)]
    a_id1 = [[-1 for _ in range(5)] for _ in range(namax)]
    a_dih0 = [[-1 for _ in range(6)] for _ in range(namax)]
    a_dih1 = [[-1 for _ in range(6)] for _ in range(namax)]
    a_im0 = [[-1 for _ in range(6)] for _ in range(namax)]
    a_im1 = [[-1 for _ in range(6)] for _ in range(namax)]
    
    x = [[0.0 for _ in range(3)] for _ in range(namax)]
    len_matrix = [[0.0 for _ in range(3)] for _ in range(3)]
    x0 = [[0.0 for _ in range(3)] for _ in range(namax)]
    x1 = [[0.0 for _ in range(3)] for _ in range(namax)]
    
    b_coeff0 = [[0.0 for _ in range(5)] for _ in range(namax)]
    a_coeff0 = [[0.0 for _ in range(5)] for _ in range(namax)]
    d_coeff0 = [[0.0 for _ in range(5)] for _ in range(namax)]
    i_coeff0 = [[0.0 for _ in range(5)] for _ in range(namax)]
    b_coeff1 = [[0.0 for _ in range(5)] for _ in range(namax)]
    a_coeff1 = [[0.0 for _ in range(5)] for _ in range(namax)]
    d_coeff1 = [[0.0 for _ in range(5)] for _ in range(namax)]
    i_coeff1 = [[0.0 for _ in range(5)] for _ in range(namax)]
    
    # 1D arrays
    cnt = [0] * namax
    idb1 = [0] * namax
    idr1 = [0] * namax
    ide1 = [0] * namax
    id0 = [0] * namax
    id1 = [0] * namax
    atom_type0 = [0] * namax
    atom_type1 = [0] * namax
    mol_id0 = [0] * namax
    mol_id1 = [0] * namax
    
    hlen = [0.0] * 3
    cen = [0.0] * 3
    q0 = [0.0] * namax
    q1 = [0.0] * namax
    mass1 = [0.0] * namax
    mass0 = [0.0] * namax
    l_lo = [0.0] * 3
    l_hi = [0.0] * 3
    
    lj0_e = [0.0] * namax
    lj0_s = [0.0] * namax
    lj1_e = [0.0] * namax
    lj1_s = [0.0] * namax
    lj2_e = [0.0] * namax
    lj2_s = [0.0] * namax
    
    try:
        # Read input file
        with open("system.data", "r") as inputFile3333:
            lines = [line.strip() for line in inputFile3333.readlines()]
            
            # Parse file line by line
            line_idx = 0
            
            # Skip first line (header)
            line_idx += 1
            
            # Skip empty lines
            while line_idx < len(lines) and lines[line_idx] == "":
                line_idx += 1
            
            # Read counts - atoms
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                kkk = int(parts[0])
                line_idx += 1
            
            # bonds
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                n_bond0 = int(parts[0])
                line_idx += 1
            
            # angles
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                n_angle0 = int(parts[0])
                line_idx += 1
            
            # dihedrals
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                n_dihedral0 = int(parts[0])
                line_idx += 1
            
            # impropers
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                n_improper0 = int(parts[0])
                line_idx += 1
            
            # Skip empty line
            while line_idx < len(lines) and lines[line_idx] == "":
                line_idx += 1
            
            # Read type counts - atom types
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                n_type0 = int(parts[0])
                line_idx += 1
            
            # bond types
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                n_bond_type0 = int(parts[0])
                line_idx += 1
            
            # angle types
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                n_angle_type0 = int(parts[0])
                line_idx += 1
            
            # dihedral types
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                n_dihedral_type0 = int(parts[0])
                line_idx += 1
            
            # improper types
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                n_improper_type0 = int(parts[0])
                line_idx += 1
            
            # Skip empty lines
            while line_idx < len(lines) and lines[line_idx] == "":
                line_idx += 1
            
            # Read box dimensions - x
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                l_lo[0] = float(parts[0])
                l_hi[0] = float(parts[1])
                line_idx += 1
            
            # y
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                l_lo[1] = float(parts[0])
                l_hi[1] = float(parts[1])
                line_idx += 1
            
            # z
            if line_idx < len(lines):
                parts = lines[line_idx].split()
                l_lo[2] = float(parts[0])
                l_hi[2] = float(parts[1])
                line_idx += 1
            
            # Skip empty lines
            while line_idx < len(lines) and lines[line_idx] == "":
                line_idx += 1
            
            # Find "Masses" section
            while line_idx < len(lines) and lines[line_idx] != "Masses":
                line_idx += 1
            
            if line_idx < len(lines):
                line_idx += 1  # Skip "Masses" line
                
                # Skip empty lines
                while line_idx < len(lines) and lines[line_idx] == "":
                    line_idx += 1
                
                # Read masses
                for i in range(n_type0):
                    if line_idx < len(lines):
                        parts = lines[line_idx].split()
                        if len(parts) >= 2:
                            mass0[i] = float(parts[1])
                        line_idx += 1
            
            # Find "Atoms" section
            while line_idx < len(lines) and not lines[line_idx].startswith("Atoms"):
                line_idx += 1
            
            if line_idx < len(lines):
                line_idx += 1  # Skip "Atoms" line
                
                # Skip empty lines
                while line_idx < len(lines) and lines[line_idx] == "":
                    line_idx += 1
                
                # Read atoms
                for i in range(kkk):
                    if line_idx < len(lines):
                        parts = lines[line_idx].split()
                        if len(parts) >= 7:
                            id0[i] = int(parts[0])
                            # Skip mol_id (parts[1])
                            atom_type0[i] = int(parts[2])
                            q0[i] = float(parts[3])
                            x0[i][0] = float(parts[4])
                            x0[i][1] = float(parts[5])
                            x0[i][2] = float(parts[6])
                        line_idx += 1
            
            # Find "Bonds" section
            while line_idx < len(lines) and lines[line_idx] != "Bonds":
                line_idx += 1
            
            if line_idx < len(lines):
                line_idx += 1  # Skip "Bonds" line
                
                # Skip empty lines
                while line_idx < len(lines) and lines[line_idx] == "":
                    line_idx += 1
                
                # Read bonds
                for i in range(n_bond0):
                    if line_idx < len(lines):
                        parts = lines[line_idx].split()
                        if len(parts) >= 4:
                            b_id0[i][0] = int(parts[0])
                            b_id0[i][1] = int(parts[1])
                            b_id0[i][2] = int(parts[2])
                            b_id0[i][3] = int(parts[3])
                        line_idx += 1
            
            # Find "Angles" section
            while line_idx < len(lines) and lines[line_idx] != "Angles":
                line_idx += 1
            
            if line_idx < len(lines):
                line_idx += 1  # Skip "Angles" line
                
                # Skip empty lines
                while line_idx < len(lines) and lines[line_idx] == "":
                    line_idx += 1
                
                # Read angles
                for i in range(n_angle0):
                    if line_idx < len(lines):
                        parts = lines[line_idx].split()
                        if len(parts) >= 5:
                            a_id0[i][0] = int(parts[0])
                            a_id0[i][1] = int(parts[1])
                            a_id0[i][2] = int(parts[2])
                            a_id0[i][3] = int(parts[3])
                            a_id0[i][4] = int(parts[4])
                        line_idx += 1
            
            # Find "Dihedrals" section
            while line_idx < len(lines) and lines[line_idx] != "Dihedrals":
                line_idx += 1
            
            if line_idx < len(lines):
                line_idx += 1  # Skip "Dihedrals" line
                
                # Skip empty lines
                while line_idx < len(lines) and lines[line_idx] == "":
                    line_idx += 1
                
                # Read dihedrals
                for i in range(n_dihedral0):
                    if line_idx < len(lines):
                        parts = lines[line_idx].split()
                        if len(parts) >= 6:
                            a_dih0[i][0] = int(parts[0])
                            a_dih0[i][1] = int(parts[1])
                            a_dih0[i][2] = int(parts[2])
                            a_dih0[i][3] = int(parts[3])
                            a_dih0[i][4] = int(parts[4])
                            a_dih0[i][5] = int(parts[5])
                        line_idx += 1
            
            # Find "Impropers" section
            while line_idx < len(lines) and lines[line_idx] != "Impropers":
                line_idx += 1
            
            if line_idx < len(lines):
                line_idx += 1  # Skip "Impropers" line
                
                # Skip empty lines
                while line_idx < len(lines) and lines[line_idx] == "":
                    line_idx += 1
                
                # Read impropers
                for i in range(n_improper0):
                    if line_idx < len(lines):
                        parts = lines[line_idx].split()
                        if len(parts) >= 6:
                            a_im0[i][0] = int(parts[0])
                            a_im0[i][1] = int(parts[1])
                            a_im0[i][2] = int(parts[2])
                            a_im0[i][3] = int(parts[3])
                            a_im0[i][4] = int(parts[4])
                            a_im0[i][5] = int(parts[5])
                        line_idx += 1
        
    except FileNotFoundError:
        print("Error: system.data file not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)
    
    # Assign molecule IDs (divide atoms into 5 groups)
    for i in range(kkk // 5):
        mol_id0[i] = 1
    
    for i in range(kkk // 5, 2 * kkk // 5):
        mol_id0[i] = 2
    
    for i in range(2 * kkk // 5, 3 * kkk // 5):
        mol_id0[i] = 3
    
    for i in range(3 * kkk // 5, 4 * kkk // 5):
        mol_id0[i] = 4
    
    for i in range(4 * kkk // 5, kkk):
        mol_id0[i] = 5
    
    # Write output file
    try:
        with open("system2.data", "w") as outputFile5:
            outputFile5.write("LAMMPS Description\n")
            outputFile5.write("\n")
            outputFile5.write(f"{kkk} atoms\n")
            outputFile5.write(f"{n_bond0} bonds\n")
            outputFile5.write(f"{n_angle0} angles\n")
            outputFile5.write(f"{n_dihedral0} dihedrals\n")
            outputFile5.write(f"{n_improper0} impropers\n")
            outputFile5.write("\n")
            outputFile5.write(f"{n_type0} atom types\n")
            outputFile5.write(f"{n_bond_type0} bond types\n")
            outputFile5.write(f"{n_angle_type0} angle types\n")
            outputFile5.write(f"{n_dihedral_type0} dihedral types\n")
            outputFile5.write(f"{n_improper_type0} improper types\n")
            outputFile5.write("\n")
            outputFile5.write(f"{l_lo[0]} {l_hi[0]} xlo xhi\n")
            outputFile5.write(f"{l_lo[1]} {l_hi[1]} ylo yhi\n")
            outputFile5.write(f"{l_lo[2]} {l_hi[2]} zlo zhi\n")
            outputFile5.write("\n")
            outputFile5.write("Masses\n")
            outputFile5.write("\n")
            
            for i in range(n_type0):
                outputFile5.write(f"{i+1} {mass0[i]}\n")
            
            outputFile5.write("\n")
            outputFile5.write("Atoms # full\n")
            outputFile5.write("\n")
            
            for i in range(kkk):
                outputFile5.write(f"{id0[i]}   {mol_id0[i]}   {atom_type0[i]}   {q0[i]}   {x0[i][0]}   {x0[i][1]}   {x0[i][2]}\n")
            
            outputFile5.write("\n")
            outputFile5.write("Bonds\n")
            outputFile5.write("\n")
            
            for i in range(n_bond0):
                outputFile5.write(f"{b_id0[i][0]}  {b_id0[i][1]}  {b_id0[i][2]}  {b_id0[i][3]}\n")
            
            outputFile5.write("\n")
            outputFile5.write("Angles\n")
            outputFile5.write("\n")
            
            for i in range(n_angle0):
                outputFile5.write(f"{a_id0[i][0]}  {a_id0[i][1]}  {a_id0[i][2]}  {a_id0[i][3]}  {a_id0[i][4]}\n")
            
            outputFile5.write("\n")
            outputFile5.write("Dihedrals\n")
            outputFile5.write("\n")
            
            for i in range(n_dihedral0):
                outputFile5.write(f"{a_dih0[i][0]}  {a_dih0[i][1]}  {a_dih0[i][2]}  {a_dih0[i][3]}  {a_dih0[i][4]}  {a_dih0[i][5]}\n")
            
            outputFile5.write("\n")
            outputFile5.write("Impropers\n")
            outputFile5.write("\n")
            
            for i in range(n_improper0):
                outputFile5.write(f"{a_im0[i][0]}  {a_im0[i][1]}  {a_im0[i][2]}  {a_im0[i][3]}  {a_im0[i][4]}  {a_im0[i][5]}\n")
            
            outputFile5.write("\n")
    
    except Exception as e:
        print(f"Error: Could not write to system2.data - {e}")
        sys.exit(1)
    
    print("Successfully processed system.data and created system2.data")
    return 0

if __name__ == "__main__":
    main()
