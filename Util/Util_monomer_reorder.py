import re
import numpy as np
from typing import List, Tuple, Dict

class MOL2Processor:
    def __init__(self):
        # Element type constants
        self.aC = 0   # Carbon
        self.aCL = 1  # Chlorine
        self.aO = 2   # Oxygen
        self.aH = 3   # Hydrogen
        self.aN = 4   # Nitrogen
        
        # Data storage variables
        self.na = 0  # Number of atoms
        self.nb = 0  # Number of bonds
        self.atoms = []  # Atom information
        self.bonds = []  # Bond information
        self.coordinates = []  # Coordinate information
        self.charges = []  # Charge information
        self.elements = []  # Element types
        
        # Indices of special atoms
        self.ad0 = []  # Nitrogen bonded with 2 hydrogens
        self.ad0_cl = []  # Carbon bonded with both oxygen and chlorine
        self.ad0_del = []  # Hydrogens to be removed
        self.ad0_cl_del = []  # Chlorines to be removed
        self.b_del = []  # Bonds to be removed (hydrogen)
        self.b_cl_del = []  # Bonds to be removed (chlorine)

    def read_mol2_file(self, filename: str, read_charges: bool = True) -> None:
        """Read MOL2 file and store molecular information"""
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Find sections
        molecule_idx = -1
        atom_idx = -1
        bond_idx = -1
        
        for i, line in enumerate(lines):
            if '@<TRIPOS>MOLECULE' in line:
                molecule_idx = i
            elif '@<TRIPOS>ATOM' in line:
                atom_idx = i
            elif '@<TRIPOS>BOND' in line:
                bond_idx = i
        
        # Read molecule information
        if molecule_idx != -1:
            molecule_info = lines[molecule_idx + 2].split()
            self.na = int(molecule_info[0])
            self.nb = int(molecule_info[1])
        
        # Read atom information
        if atom_idx != -1:
            self.atoms = []
            self.coordinates = []
            self.charges = []
            self.elements = []
            
            for i in range(atom_idx + 1, atom_idx + 1 + self.na):
                if i < len(lines):
                    parts = lines[i].split()
                    if len(parts) >= 6:
                        atom_info = {
                            'id': int(parts[0]),
                            'name': parts[1],  # Atom name (C1, N2, etc.)
                            'x': float(parts[2]),
                            'y': float(parts[3]),
                            'z': float(parts[4]),
                            'type': parts[5],  # Element type (actual element symbol)
                            'subst_id': parts[6] if len(parts) > 6 else '1',
                            'subst_name': parts[7] if len(parts) > 7 else 'UNL1',
                            'charge': float(parts[8]) if len(parts) > 8 and read_charges else 0.0
                        }
                        self.atoms.append(atom_info)
                        self.coordinates.append([atom_info['x'], atom_info['y'], atom_info['z']])
                        self.charges.append(atom_info['charge'])
                        
                        # Use the second string (element type, not atom name) as in C++ code
                        # In MOL2, parts[5] is the element type
                        element_type = self.identify_element(parts[5])
                        self.elements.append(element_type)
        
        # Read bond information
        if bond_idx != -1:
            self.bonds = []
            for i in range(bond_idx + 1, bond_idx + 1 + self.nb):
                if i < len(lines):
                    parts = lines[i].split()
                    if len(parts) >= 4:
                        bond_info = {
                            'id': int(parts[0]),
                            'atom1': int(parts[1]),
                            'atom2': int(parts[2]),
                            'type': parts[3]
                        }
                        self.bonds.append(bond_info)

    def identify_element(self, atom_type: str) -> int:
        """Identify element from atom type"""
        atom_type = atom_type.upper().strip()
        
        # In MOL2 files, atom types are in the form C.2, N.ar, O.2, Cl, H, etc.
        # Extract the element symbol before the dot (.)
        base_element = atom_type.split('.')[0]
        
        # C++ order: N -> CL -> O -> H -> C -> others(10)
        if base_element.startswith('N'):
            return self.aN
        elif base_element.startswith('CL') or base_element.startswith('Cl'):
            return self.aCL
        elif base_element.startswith('O'):
            return self.aO
        elif base_element.startswith('H'):
            return self.aH
        elif base_element.startswith('C'):
            return self.aC
        else:
            return 10  # Unknown element

    def find_special_atoms(self) -> None:
        """Find special atoms (nitrogen bonded with 2 hydrogens, carbon bonded with oxygen-chlorine)"""
        # Find nitrogen bonded with 2 hydrogens
        for i, elem in enumerate(self.elements):
            if elem == self.aN:
                h_bonds = []
                h_count = 0
                
                for j, bond in enumerate(self.bonds):
                    if bond['atom1'] == i + 1 and self.elements[bond['atom2'] - 1] == self.aH:
                        h_count += 1
                        h_bonds.append((j, bond['atom2'] - 1))
                    elif bond['atom2'] == i + 1 and self.elements[bond['atom1'] - 1] == self.aH:
                        h_count += 1
                        h_bonds.append((j, bond['atom1'] - 1))
                
                if h_count == 2:
                    self.ad0.append(i)
                    # Set the first hydrogen as removal target
                    if h_bonds:
                        bond_idx, h_atom_idx = h_bonds[0]
                        self.ad0_del.append(h_atom_idx)
                        self.b_del.append(bond_idx)
                    break
        
        # Find carbon bonded with both oxygen and chlorine
        for i, elem in enumerate(self.elements):
            if elem == self.aC:
                has_oxygen = False
                has_chlorine = False
                cl_bond_idx = -1
                cl_atom_idx = -1
                
                for j, bond in enumerate(self.bonds):
                    if bond['atom1'] == i + 1:
                        if self.elements[bond['atom2'] - 1] == self.aO:
                            has_oxygen = True
                        elif self.elements[bond['atom2'] - 1] == self.aCL:
                            has_chlorine = True
                            cl_bond_idx = j
                            cl_atom_idx = bond['atom2'] - 1
                    elif bond['atom2'] == i + 1:
                        if self.elements[bond['atom1'] - 1] == self.aO:
                            has_oxygen = True
                        elif self.elements[bond['atom1'] - 1] == self.aCL:
                            has_chlorine = True
                            cl_bond_idx = j
                            cl_atom_idx = bond['atom1'] - 1
                
                if has_oxygen and has_chlorine:
                    self.ad0_cl.append(i)
                    self.ad0_cl_del.append(cl_atom_idx)
                    self.b_cl_del.append(cl_bond_idx)
                    break

    def redistribute_charges(self) -> None:
        """Redistribute charges of atoms to be removed to remaining hydrogens (same logic as C++ code)"""
        if not self.ad0_del or not self.ad0_cl_del:
            return
        
        # C++ code: Calculate total charge of atoms excluding those to be removed
        sum_q = 0.0
        for j in range(self.na):
            if j == self.ad0_del[0] or j == self.ad0_cl_del[0]:
                continue
            sum_q += self.charges[j]
        
        # Count remaining hydrogen atoms (excluding hydrogen to be removed)
        hydrogen_count = 0
        for j in range(self.na):
            if self.elements[j] == self.aH and j != self.ad0_del[0]:
                hydrogen_count += 1
        
        if hydrogen_count > 0:
            charge_adjustment = sum_q / hydrogen_count
            
            # Adjust charge for each hydrogen atom (subtract, not add)
            for j in range(self.na):
                if self.elements[j] == self.aH and j != self.ad0_del[0]:
                    self.charges[j] = self.charges[j] - charge_adjustment

    def adjust_molecule_orientation(self) -> None:
        """Adjust molecular orientation (axis with largest nitrogen-carbon distance becomes x-axis)"""
        if not self.ad0 or not self.ad0_cl:
            return
        
        n_coord = self.coordinates[self.ad0[0]]
        c_coord = self.coordinates[self.ad0_cl[0]]
        
        # Calculate distance for each axis
        dx = abs(n_coord[0] - c_coord[0])
        dy = abs(n_coord[1] - c_coord[1])
        dz = abs(n_coord[2] - c_coord[2])
        
        # Set the axis with the largest distance as x-axis
        coordinates_array = np.array(self.coordinates)
        
        if dy > dx and dy > dz:
            # Rotate y-axis to x-axis
            coordinates_array[:, [0, 1]] = coordinates_array[:, [1, 0]]
        elif dz > dx and dz > dy:
            # Rotate z-axis to x-axis
            coordinates_array[:, [0, 2]] = coordinates_array[:, [2, 0]]
        
        # Adjust so that nitrogen has smaller x value than carbon
        n_x = coordinates_array[self.ad0[0], 0]
        c_x = coordinates_array[self.ad0_cl[0], 0]
        
        if n_x > c_x:
            avg_x = np.mean(coordinates_array[:, 0])
            coordinates_array[:, 0] = 2 * avg_x - coordinates_array[:, 0]
        
        # Make all x coordinates positive
        min_x = np.min(coordinates_array[:, 0])
        coordinates_array[:, 0] = coordinates_array[:, 0] + abs(min_x)
        
        self.coordinates = coordinates_array.tolist()

    def write_output_file(self, output_filename: str, template_filename: str) -> None:
        """Save rearranged molecular structure as new MOL2 file"""
        # Read original format information from template file
        template_atoms = []
        template_bonds = []
        
        with open(template_filename, 'r') as f:
            lines = f.readlines()
        
        atom_idx = -1
        bond_idx = -1
        
        for i, line in enumerate(lines):
            if '@<TRIPOS>ATOM' in line:
                atom_idx = i
            elif '@<TRIPOS>BOND' in line:
                bond_idx = i
        
        # Read atom information from template
        if atom_idx != -1:
            for i in range(atom_idx + 1, atom_idx + 1 + self.na):
                if i < len(lines):
                    parts = lines[i].split()
                    template_atoms.append({
                        'name': parts[1],
                        'type': parts[5] if len(parts) > 5 else 'UNL',
                        'subst_id': parts[6] if len(parts) > 6 else '1',
                        'subst_name': parts[7] if len(parts) > 7 else 'UNL1'
                    })
        
        # Read bond information from template
        if bond_idx != -1:
            for i in range(bond_idx + 1, bond_idx + 1 + self.nb):
                if i < len(lines):
                    parts = lines[i].split()
                    template_bonds.append({
                        'type': parts[3] if len(parts) > 3 else '1'
                    })
        
        # Create new atom list excluding atoms to be removed
        new_atoms = []
        atom_mapping = {}  # Old index -> New index mapping
        new_idx = 1
        
        # First add nitrogen (N2000)
        if self.ad0:
            new_atoms.append({
                'id': new_idx,
                'name': 'N2000',
                'coord': self.coordinates[self.ad0[0]],
                'type': 'n',
                'charge': self.charges[self.ad0[0]]
            })
            atom_mapping[self.ad0[0]] = new_idx
            new_idx += 1
        
        # Next add carbon (C2000)
        if self.ad0_cl:
            template_type = template_atoms[self.ad0_cl[0]]['type'] if self.ad0_cl[0] < len(template_atoms) else 'c3'
            new_atoms.append({
                'id': new_idx,
                'name': 'C2000',
                'coord': self.coordinates[self.ad0_cl[0]],
                'type': template_type,
                'charge': self.charges[self.ad0_cl[0]]
            })
            atom_mapping[self.ad0_cl[0]] = new_idx
            new_idx += 1
        
        # Add remaining atoms (excluding removal targets)
        for i in range(self.na):
            if (i in self.ad0 or i in self.ad0_cl or 
                i in self.ad0_del or i in self.ad0_cl_del):
                continue
            
            template_info = template_atoms[i] if i < len(template_atoms) else {'name': f'X{i}', 'type': 'du'}
            new_atoms.append({
                'id': new_idx,
                'name': template_info['name'],
                'coord': self.coordinates[i],
                'type': template_info['type'],
                'charge': self.charges[i]
            })
            atom_mapping[i] = new_idx
            new_idx += 1
        
        # Create new bond list excluding bonds to be removed
        new_bonds = []
        bond_idx_new = 1
        
        for i, bond in enumerate(self.bonds):
            if i in self.b_del or i in self.b_cl_del:
                continue
            
            old_atom1 = bond['atom1'] - 1
            old_atom2 = bond['atom2'] - 1
            
            if old_atom1 in atom_mapping and old_atom2 in atom_mapping:
                bond_type = template_bonds[i]['type'] if i < len(template_bonds) else '1'
                new_bonds.append({
                    'id': bond_idx_new,
                    'atom1': atom_mapping[old_atom1],
                    'atom2': atom_mapping[old_atom2],
                    'type': bond_type
                })
                bond_idx_new += 1
        
        # Write output file
        with open(output_filename, 'w') as f:
            f.write("@<TRIPOS>MOLECULE\n")
            f.write("*****\n")
            f.write(f" {len(new_atoms)} {len(new_bonds)} 0 0 0\n")
            f.write("SMALL\n")
            f.write("GASTEIGER\n")
            f.write("\n")
            f.write("@<TRIPOS>ATOM\n")
            
            # Write atom information
            for atom in new_atoms:
                coord = atom['coord']
                f.write(f"      {atom['id']} {atom['name']}           "
                       f"{coord[0]:.4f}   {coord[1]:.4f}    {coord[2]:.4f} "
                       f"{atom['type']}  1  UNL1    {atom['charge']:.6f}\n")
            
            f.write("@<TRIPOS>BOND\n")
            
            # Write bond information
            for bond in new_bonds:
                f.write(f"     {bond['id']}    {bond['atom1']}    {bond['atom2']}    {bond['type']}\n")
            
            f.write("@<TRIPOS>SUBSTRUCTURE\n")
            f.write("      1 ***         1 TEMP              0 ****  ****    0 ROOT\n")

    def process_mol2_files(self, input_file1: str, input_file2: str, output_file: str) -> None:
        """Execute the entire processing workflow"""
        # Read coordinates and charges from first file
        self.read_mol2_file(input_file1, read_charges=True)
        
        # Find special atoms
        self.find_special_atoms()
        
        if not self.ad0 or not self.ad0_cl:
            print("Error: Required special atoms not found.")
            return
        
        # Redistribute charges
        self.redistribute_charges()
        
        # Adjust molecular orientation
        self.adjust_molecule_orientation()
        
        # Write output file
        self.write_output_file(output_file, input_file2)
        
        print("1 molecule converted to reordered format.")


def main():
    """Main function - usage example"""
    processor = MOL2Processor()
    
    try:
        # Execute file processing
        processor.process_mol2_files(
            input_file1="monomer_com1.mol2",  # Coordinate/charge data
            input_file2="monomer_com2.mol2",  # Format information
            output_file="monomer_reorder2.mol2"  # Output file
        )
        
    except FileNotFoundError as e:
        print(f"File not found: {e}")
    except Exception as e:
        print(f"Error during processing: {e}")


if __name__ == "__main__":
    main()
