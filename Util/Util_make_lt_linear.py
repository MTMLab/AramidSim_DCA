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
            output_file.write("import H_head.lt\n")
            output_file.write("import H_tail.lt\n")
            output_file.write("import PPTA.lt\n")
            output_file.write("\n")
            output_file.write("aramid inherits GAFF {\n")
            output_file.write("\n")
            output_file.write("   create_var {$mol}\n")
            output_file.write("\n")
            
            # Polymer chain generation
            d0 = 0
            
            # Fiber chain
            output_file.write("monomers[0] = new H_head\n")
            d0 = 12
            output_file.write(f"monomers[1] = new PPTA.move({d0},0,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[2] = new cation.move({d0},0,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[3] = new PPTA.move({d0},0,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[4] = new cation.move({d0},0,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[5] = new PPTA.move({d0},0,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[6] = new PPTA.move({d0},0,0)\n")  # PPTA
            d0 += 12
            output_file.write(f"monomers[7] = new cation.move({d0},0,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[8] = new cation.move({d0},0,0)\n")  # cation
            d0 += interval
            output_file.write(f"monomers[9] = new H_tail.move({d0},0,0)\n")
            d0 += 12
            
            # Bond structure definition
            output_file.write("   write('Data Bond List') {\n")
            
            # Linear bonds connecting sequential monomers
            bonds = [
                ("b1", "monomers[0]/H_head", "monomers[1]/atom1"),
                ("b2", "monomers[1]/atom2", "monomers[2]/atom1"),
                ("b3", "monomers[2]/atom2", "monomers[3]/atom1"),
                ("b4", "monomers[3]/atom2", "monomers[4]/atom1"),
                ("b5", "monomers[4]/atom2", "monomers[5]/atom1"),
                ("b6", "monomers[5]/atom2", "monomers[6]/atom1"),
                ("b7", "monomers[6]/atom2", "monomers[7]/atom1"),
                ("b8", "monomers[7]/atom2", "monomers[8]/atom1"),
                ("b9", "monomers[8]/atom2", "monomers[9]/H_tail")
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
