#!/usr/bin/env python3
import sys
import os

def calculate_interval(filename="monomer_reorder2.mol2"):
    """
    Read MOL2 file and calculate interval based on x-coordinate range
    Returns range * 1.2 (20% margin)
    """
    try:
        with open(filename, 'r') as mol2_file:
            lines = mol2_file.readlines()
    except FileNotFoundError:
        return 10.0  # default value
    
    atom_section = False
    x_min = float('inf')
    x_max = float('-inf')
    atom_count = 0
    
    for line in lines:
        line = line.strip()
        
        # Check for @<TRIPOS>ATOM section
        if "@<TRIPOS>ATOM" in line:
            atom_section = True
            continue
        
        # Check for end of atom section
        if atom_section and line.startswith("@<TRIPOS>"):
            break
        
        # Parse atom coordinates
        if atom_section and line:
            try:
                parts = line.split()
                if len(parts) >= 6:
                    # MOL2 format: atom_id atom_name x y z atom_type ...
                    x = float(parts[2])
                    
                    if x < x_min:
                        x_min = x
                    if x > x_max:
                        x_max = x
                    atom_count += 1
            except (ValueError, IndexError):
                continue
    
    if atom_count == 0:
        return 10.0  # default value
    
    range_x = x_max - x_min
    interval = range_x * 1.2  # 20% margin
    
    return interval

def generate_polymer_structure():
    """
    Main function to generate polymer structure files
    """
    # Calculate interval from MOL2 file
    interval = calculate_interval("monomer_reorder2.mol2")
    
    # Generate polymer_new.lt file
    try:
        with open("polymer_new.lt", 'w') as output_file:
            # Header
            output_file.write("import monomer_add.lt\n")
            output_file.write("import PPTA.lt\n")
            output_file.write("\n")
            output_file.write("aramid inherits GAFF {\n")
            output_file.write("\n")
            output_file.write("   create_var {$mol}\n")
            output_file.write("\n")
            
            # Polymer chain generation
            d0 = 0
            
            # First chain (monomers 0-7)
            output_file.write("monomers[0] = new PPTA.move(0,0,0)\n")  # PPTA
            d0 = 12
            output_file.write(f"monomers[1] = new cation.move({d0},0,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[2] = new PPTA.move({d0},0,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[3] = new cation.move({d0},0,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[4] = new cation.move({d0},0,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[5] = new cation.move({d0},0,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[6] = new PPTA.move({d0},0,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[7] = new PPTA.move({d0},0,0)\n")  # PPTA
            d0 += 12
            
            # Second chain (monomers 8-15)
            output_file.write("monomers[8] = new cation.move(0,40,0)\n")  # cation
            d0 = interval
            output_file.write(f"monomers[9] = new PPTA.move({d0},40,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[10] = new cation.move({d0},40,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[11] = new PPTA.move({d0},40,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[12] = new cation.move({d0},40,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[13] = new cation.move({d0},40,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[14] = new PPTA.move({d0},40,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[15] = new PPTA.move({d0},40,0)\n")  # PPTA
            d0 += 12
            
            # Third chain (monomers 16-23)
            output_file.write("monomers[16] = new PPTA.move(0,80,0)\n")  # PPTA
            d0 = 12
            output_file.write(f"monomers[17] = new PPTA.move({d0},80,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[18] = new cation.move({d0},80,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[19] = new PPTA.move({d0},80,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[20] = new cation.move({d0},80,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[21] = new cation.move({d0},80,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[22] = new cation.move({d0},80,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[23] = new PPTA.move({d0},80,0)\n")  # PPTA
            d0 += 12
            
            # Fourth chain (monomers 24-31)
            output_file.write("monomers[24] = new PPTA.move(0,120,0)\n")  # PPTA
            d0 = 12
            output_file.write(f"monomers[25] = new PPTA.move({d0},120,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[26] = new cation.move({d0},120,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[27] = new cation.move({d0},120,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[28] = new PPTA.move({d0},120,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[29] = new cation.move({d0},120,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[30] = new cation.move({d0},120,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[31] = new PPTA.move({d0},120,0)\n")  # PPTA
            d0 += 12
            
            # Fifth chain (monomers 32-39)
            output_file.write("monomers[32] = new PPTA.move(0,160,0)\n")  # PPTA
            d0 = 12
            output_file.write(f"monomers[33] = new cation.move({d0},160,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[34] = new cation.move({d0},160,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[35] = new cation.move({d0},160,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[36] = new cation.move({d0},160,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[37] = new PPTA.move({d0},160,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[38] = new PPTA.move({d0},160,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[39] = new PPTA.move({d0},160,0)\n")  # PPTA
            d0 += 12
            
            # Bond structure definition
            output_file.write("   write('Data Bond List') {\n")
            
            # First chain bonds (0-7)
            bonds = [
                ("b1", "monomers[0]/atom2", "monomers[1]/atom1"),
                ("b2", "monomers[1]/atom2", "monomers[2]/atom1"),
                ("b3", "monomers[2]/atom2", "monomers[3]/atom1"),
                ("b4", "monomers[3]/atom2", "monomers[4]/atom1"),
                ("b5", "monomers[4]/atom2", "monomers[5]/atom1"),
                ("b6", "monomers[5]/atom2", "monomers[6]/atom1"),
                ("b7", "monomers[6]/atom2", "monomers[7]/atom1"),
                ("b8", "monomers[7]/atom2", "monomers[0]/atom1"),
                
                # Second chain bonds (8-15)
                ("b9", "monomers[8]/atom2", "monomers[9]/atom1"),
                ("b10", "monomers[9]/atom2", "monomers[10]/atom1"),
                ("b11", "monomers[10]/atom2", "monomers[11]/atom1"),
                ("b12", "monomers[11]/atom2", "monomers[12]/atom1"),
                ("b13", "monomers[12]/atom2", "monomers[13]/atom1"),
                ("b14", "monomers[13]/atom2", "monomers[14]/atom1"),
                ("b15", "monomers[14]/atom2", "monomers[15]/atom1"),
                ("b16", "monomers[15]/atom2", "monomers[8]/atom1"),
                
                # Third chain bonds (16-23)
                ("b17", "monomers[16]/atom2", "monomers[17]/atom1"),
                ("b18", "monomers[17]/atom2", "monomers[18]/atom1"),
                ("b19", "monomers[18]/atom2", "monomers[19]/atom1"),
                ("b20", "monomers[19]/atom2", "monomers[20]/atom1"),
                ("b21", "monomers[20]/atom2", "monomers[21]/atom1"),
                ("b22", "monomers[21]/atom2", "monomers[22]/atom1"),
                ("b23", "monomers[22]/atom2", "monomers[23]/atom1"),
                ("b24", "monomers[23]/atom2", "monomers[16]/atom1"),
                
                # Fourth chain bonds (24-31)
                ("b25", "monomers[24]/atom2", "monomers[25]/atom1"),
                ("b26", "monomers[25]/atom2", "monomers[26]/atom1"),
                ("b27", "monomers[26]/atom2", "monomers[27]/atom1"),
                ("b28", "monomers[27]/atom2", "monomers[28]/atom1"),
                ("b29", "monomers[28]/atom2", "monomers[29]/atom1"),
                ("b30", "monomers[29]/atom2", "monomers[30]/atom1"),
                ("b31", "monomers[30]/atom2", "monomers[31]/atom1"),
                ("b32", "monomers[31]/atom2", "monomers[24]/atom1"),
                
                # Fifth chain bonds (32-39)
                ("b33", "monomers[32]/atom2", "monomers[33]/atom1"),
                ("b34", "monomers[33]/atom2", "monomers[34]/atom1"),
                ("b35", "monomers[34]/atom2", "monomers[35]/atom1"),
                ("b36", "monomers[35]/atom2", "monomers[36]/atom1"),
                ("b37", "monomers[36]/atom2", "monomers[37]/atom1"),
                ("b38", "monomers[37]/atom2", "monomers[38]/atom1"),
                ("b39", "monomers[38]/atom2", "monomers[39]/atom1"),
                ("b40", "monomers[39]/atom2", "monomers[32]/atom1"),
            ]
            
            for bond_id, atom1, atom2 in bonds:
                output_file.write(f"    $bond:{bond_id}  $atom:{atom1} $atom:{atom2}\n")
            
            output_file.write("}\n")
            output_file.write("}\n")
        
        print("polymer_new.lt file generated successfully!")
        
    except IOError as e:
        return False
    
    # Generate system_new.lt file
    try:
        with open("system_new.lt", 'w') as output_file:
            output_file.write("import polymer.lt\n")
            output_file.write("polymers = new aramid\n")
            output_file.write("write_once(\"Data Boundary\") {\n")
            output_file.write(f"0.0    {d0}  xlo xhi\n")
            output_file.write("0.0    200  ylo yhi\n")
            output_file.write("0.0    40  zlo zhi\n")
            output_file.write("}\n")
        
        print("system_new.lt file generated successfully!")
        
    except IOError as e:
        return False
    
    return True

def main():
    """
    Main function
    """
    generate_polymer_structure()

if __name__ == "__main__":
    main()
