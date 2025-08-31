#!/usr/bin/env python3

from rdkit import Chem
from rdkit import RDLogger
import os
import shutil
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json

# Hide RDKit warnings
RDLogger.DisableLog('rdApp.*')

# Fixed solvent SMILES
FIXED_SOLVENT_SMILES = "O=C1CCCN1C"

def clean_set_directory():
    """
    Clean the set directory by removing all files and subdirectories except structures folder,
    and ensure structures folder exists but is empty.
    """
    set_dir = os.path.join(os.getcwd(), "set")
    structures_dir = os.path.join(set_dir, "structures")
    
    try:
        # If set directory doesn't exist, create it
        if not os.path.exists(set_dir):
            os.makedirs(set_dir)
        
        # Get list of all items in set directory
        if os.path.exists(set_dir):
            items = os.listdir(set_dir)
            for item in items:
                item_path = os.path.join(set_dir, item)
                
                # If it's the structures directory, clean it but don't remove it
                if item == "structures":
                    if os.path.isdir(item_path):
                        # Remove all contents of structures directory
                        for sub_item in os.listdir(item_path):
                            sub_item_path = os.path.join(item_path, sub_item)
                            try:
                                if os.path.isfile(sub_item_path) or os.path.islink(sub_item_path):
                                    os.unlink(sub_item_path)
                                elif os.path.isdir(sub_item_path):
                                    shutil.rmtree(sub_item_path)
                            except Exception:
                                pass
                else:
                    # Remove all other files and directories
                    try:
                        if os.path.isfile(item_path) or os.path.islink(item_path):
                            os.unlink(item_path)
                        elif os.path.isdir(item_path):
                            shutil.rmtree(item_path)
                    except Exception:
                        pass
        
        # Ensure structures directory exists
        if not os.path.exists(structures_dir):
            os.makedirs(structures_dir)
        
    except Exception:
        return False
    
    return True

def clean_stretched_directory():
    """
    Clean the Stretched directory according to specific rules.
    """
    stretched_dir = os.path.join(os.getcwd(), "Stretched")
    
    try:
        if not os.path.exists(stretched_dir):
            return True
        
        # 1) Clean root Stretched directory - remove all files except specific ones
        files_to_keep = {"E_h_bond_back.txt", "output_all_back.txt"}
        items = os.listdir(stretched_dir)
        for item in items:
            item_path = os.path.join(stretched_dir, item)
            if os.path.isfile(item_path) and item not in files_to_keep:
                try:
                    os.unlink(item_path)
                except Exception:
                    pass
        
        # 2) Clean ./Stretched/structures - remove all files
        structures_dir = os.path.join(stretched_dir, "structures")
        if os.path.exists(structures_dir):
            for item in os.listdir(structures_dir):
                item_path = os.path.join(structures_dir, item)
                try:
                    if os.path.isfile(item_path) or os.path.islink(item_path):
                        os.unlink(item_path)
                    elif os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                except Exception:
                    pass
        else:
            os.makedirs(structures_dir)
        
        # 3) Clean ./Stretched/lammps - remove all files
        lammps_dir = os.path.join(stretched_dir, "lammps")
        if os.path.exists(lammps_dir):
            for item in os.listdir(lammps_dir):
                item_path = os.path.join(lammps_dir, item)
                try:
                    if os.path.isfile(item_path) or os.path.islink(item_path):
                        os.unlink(item_path)
                    elif os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                except Exception:
                    pass
        else:
            os.makedirs(lammps_dir)
        
        # 4) Clean ./Stretched/mol2tolt/test/ - remove all files
        test_dir = os.path.join(stretched_dir, "mol2tolt", "test")
        if os.path.exists(test_dir):
            for item in os.listdir(test_dir):
                item_path = os.path.join(test_dir, item)
                try:
                    if os.path.isfile(item_path) or os.path.islink(item_path):
                        os.unlink(item_path)
                    elif os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                except Exception:
                    pass
        else:
            os.makedirs(test_dir, exist_ok=True)
        
        # 5) Clean ./Stretched/moltemplates/ - remove all files except specific ones
        moltemplates_dir = os.path.join(stretched_dir, "moltemplates")
        files_to_keep_moltemplates = {"PPTA.lt", "PPTA_monomer.mol2"}
        if os.path.exists(moltemplates_dir):
            for item in os.listdir(moltemplates_dir):
                item_path = os.path.join(moltemplates_dir, item)
                if os.path.isfile(item_path) and item not in files_to_keep_moltemplates:
                    try:
                        os.unlink(item_path)
                    except Exception:
                        pass
                elif os.path.isdir(item_path):
                    try:
                        shutil.rmtree(item_path)
                    except Exception:
                        pass
        else:
            os.makedirs(moltemplates_dir)
        
    except Exception:
        return False
    
    return True

def clean_solution_directory():
    """
    Clean the Solution directory according to specific rules.
    """
    solution_dir = os.path.join(os.getcwd(), "Solution")
    
    try:
        if not os.path.exists(solution_dir):
            return True
        
        # 1) Clean root Solution directory - remove all files except specific ones
        files_to_keep = {"E_h_bond_back.txt", "output_all_back.txt", "group_polymer.py", "group_solvent.py"}
        items = os.listdir(solution_dir)
        for item in items:
            item_path = os.path.join(solution_dir, item)
            if os.path.isfile(item_path) and item not in files_to_keep:
                try:
                    os.unlink(item_path)
                except Exception:
                    pass
        
        # 2) Clean ./Solution/structures - remove all files
        structures_dir = os.path.join(solution_dir, "structures")
        if os.path.exists(structures_dir):
            for item in os.listdir(structures_dir):
                item_path = os.path.join(structures_dir, item)
                try:
                    if os.path.isfile(item_path) or os.path.islink(item_path):
                        os.unlink(item_path)
                    elif os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                except Exception:
                    pass
        else:
            os.makedirs(structures_dir)
        
        # 3) Clean ./Solution/solvent - remove all files
        solvent_dir = os.path.join(solution_dir, "solvent")
        if os.path.exists(solvent_dir):
            for item in os.listdir(solvent_dir):
                item_path = os.path.join(solvent_dir, item)
                try:
                    if os.path.isfile(item_path) or os.path.islink(item_path):
                        os.unlink(item_path)
                    elif os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                except Exception:
                    pass
        else:
            os.makedirs(solvent_dir)
        
        # 4) Clean ./Solution/lammps - remove all files
        lammps_dir = os.path.join(solution_dir, "lammps")
        if os.path.exists(lammps_dir):
            for item in os.listdir(lammps_dir):
                item_path = os.path.join(lammps_dir, item)
                try:
                    if os.path.isfile(item_path) or os.path.islink(item_path):
                        os.unlink(item_path)
                    elif os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                except Exception:
                    pass
        else:
            os.makedirs(lammps_dir)
        
        # 5) Clean ./Solution/mol2tolt/test/ - remove all files
        test_dir = os.path.join(solution_dir, "mol2tolt", "test")
        if os.path.exists(test_dir):
            for item in os.listdir(test_dir):
                item_path = os.path.join(test_dir, item)
                try:
                    if os.path.isfile(item_path) or os.path.islink(item_path):
                        os.unlink(item_path)
                    elif os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                except Exception:
                    pass
        else:
            os.makedirs(test_dir, exist_ok=True)
        
        # 6) Clean ./Solution/moltemplates/ - remove all files except specific ones
        moltemplates_dir = os.path.join(solution_dir, "moltemplates")
        files_to_keep_moltemplates = {"PPTA.lt", "PPTA_monomer.mol2", "H_head.lt", "H_tail.lt"}
        if os.path.exists(moltemplates_dir):
            for item in os.listdir(moltemplates_dir):
                item_path = os.path.join(moltemplates_dir, item)
                if os.path.isfile(item_path) and item not in files_to_keep_moltemplates:
                    try:
                        os.unlink(item_path)
                    except Exception:
                        pass
                elif os.path.isdir(item_path):
                    try:
                        shutil.rmtree(item_path)
                    except Exception:
                        pass
        else:
            os.makedirs(moltemplates_dir)
        
        # 7) Clean ./Solution/moltemplates_solvent/ - remove all files
        moltemplates_solvent_dir = os.path.join(solution_dir, "moltemplates_solvent")
        if os.path.exists(moltemplates_solvent_dir):
            for item in os.listdir(moltemplates_solvent_dir):
                item_path = os.path.join(moltemplates_solvent_dir, item)
                try:
                    if os.path.isfile(item_path) or os.path.islink(item_path):
                        os.unlink(item_path)
                    elif os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                except Exception:
                    pass
        else:
            os.makedirs(moltemplates_solvent_dir)
        
        # 8) Clean ./Solution/moltemplates_solvent_single/ - remove all files except system.lt
        moltemplates_solvent_single_dir = os.path.join(solution_dir, "moltemplates_solvent_single")
        files_to_keep_single = {"system.lt"}
        if os.path.exists(moltemplates_solvent_single_dir):
            for item in os.listdir(moltemplates_solvent_single_dir):
                item_path = os.path.join(moltemplates_solvent_single_dir, item)
                if os.path.isfile(item_path) and item not in files_to_keep_single:
                    try:
                        os.unlink(item_path)
                    except Exception:
                        pass
                elif os.path.isdir(item_path):
                    try:
                        shutil.rmtree(item_path)
                    except Exception:
                        pass
        else:
            os.makedirs(moltemplates_solvent_single_dir)
        
        # 9) Clean ./Solution/ratio_calculation/ - keep only specific files and folders
        ratio_calc_dir = os.path.join(solution_dir, "ratio_calculation")
        if os.path.exists(ratio_calc_dir):
            files_to_keep_ratio = {"factoring.py", "number_solvnet.py", "solvent_lt.py"}
            folders_to_keep_ratio = {"polymer", "solvent"}
            
            for item in os.listdir(ratio_calc_dir):
                item_path = os.path.join(ratio_calc_dir, item)
                
                if os.path.isfile(item_path) and item not in files_to_keep_ratio:
                    try:
                        os.unlink(item_path)
                    except Exception:
                        pass
                
                elif os.path.isdir(item_path) and item not in folders_to_keep_ratio:
                    try:
                        shutil.rmtree(item_path)
                    except Exception:
                        pass
            
            # Clean polymer folder - keep only specific files
            polymer_dir = os.path.join(ratio_calc_dir, "polymer")
            polymer_files_to_keep = {"polymer_mass.py", "run_polymer_mass.npt2", "system.in", "system.in.init", "system.in.settings"}
            if os.path.exists(polymer_dir):
                for item in os.listdir(polymer_dir):
                    item_path = os.path.join(polymer_dir, item)
                    if item not in polymer_files_to_keep:
                        try:
                            if os.path.isfile(item_path) or os.path.islink(item_path):
                                os.unlink(item_path)
                            elif os.path.isdir(item_path):
                                shutil.rmtree(item_path)
                        except Exception:
                            pass
            else:
                os.makedirs(polymer_dir)
            
            # Clean solvent folder - keep only specific files
            solvent_calc_dir = os.path.join(ratio_calc_dir, "solvent")
            solvent_files_to_keep = {"run_solvent_mass.npt2", "solvent_mass.py", "system.in", "system.in.init", "system.in.settings"}
            if os.path.exists(solvent_calc_dir):
                for item in os.listdir(solvent_calc_dir):
                    item_path = os.path.join(solvent_calc_dir, item)
                    if item not in solvent_files_to_keep:
                        try:
                            if os.path.isfile(item_path) or os.path.islink(item_path):
                                os.unlink(item_path)
                            elif os.path.isdir(item_path):
                                shutil.rmtree(item_path)
                        except Exception:
                            pass
            else:
                os.makedirs(solvent_calc_dir)
        else:
            # Create ratio_calculation directory structure if it doesn't exist
            os.makedirs(ratio_calc_dir)
            os.makedirs(os.path.join(ratio_calc_dir, "polymer"))
            os.makedirs(os.path.join(ratio_calc_dir, "solvent"))
        
        # 10) Clean ./Solution/mol2tolt_solvent/test/ - remove all files
        mol2tolt_solvent_test_dir = os.path.join(solution_dir, "mol2tolt_solvent", "test")
        if os.path.exists(mol2tolt_solvent_test_dir):
            for item in os.listdir(mol2tolt_solvent_test_dir):
                item_path = os.path.join(mol2tolt_solvent_test_dir, item)
                try:
                    if os.path.isfile(item_path) or os.path.islink(item_path):
                        os.unlink(item_path)
                    elif os.path.isdir(item_path):
                        shutil.rmtree(item_path)
                except Exception:
                    pass
        else:
            os.makedirs(mol2tolt_solvent_test_dir, exist_ok=True)
        
    except Exception:
        return False
    
    return True


def get_polymer_configuration():
    """
    Get fixed polymer configuration (binary type with fixed ratios)
    
    Returns:
        dict: Fixed polymer configuration dictionary
    """
    # Fixed configuration - no print output here
    polymer_type = "binary"
    ppta_count = 4
    oda_count = 0  # Binary type doesn't use ODA
    cation_count = 4
    
    return {
        'polymer_type': polymer_type,
        'ppta_count': ppta_count,
        'oda_count': oda_count,
        'cation_count': cation_count
    }

def save_polymer_config(config):
    """
    Save polymer configuration to JSON file
    
    Args:
        config (dict): Polymer configuration dictionary
        
    Returns:
        bool: True if successful, False otherwise
    """
    config_path = os.path.join(os.getcwd(), "set", "polymer_config.json")
    
    try:
        os.makedirs(os.path.dirname(config_path), exist_ok=True)
        with open(config_path, 'w', encoding='utf-8') as f:
            json.dump(config, f, indent=2, ensure_ascii=False)
        print(f"Polymer configuration saved to {config_path}")
        return True
    except Exception as e:
        print(f"Error saving polymer configuration: {e}")
        return False

def load_polymer_config():
    """
    Load polymer configuration from JSON file
    
    Returns:
        dict or None: Polymer configuration dictionary, or None if file not found
    """
    config_path = os.path.join(os.getcwd(), "set", "polymer_config.json")
    
    try:
        if os.path.exists(config_path):
            with open(config_path, 'r', encoding='utf-8') as f:
                config = json.load(f)
            return config
        else:
            print(f"Warning: Polymer configuration file not found at {config_path}")
            return None
    except Exception as e:
        print(f"Error reading polymer configuration: {e}")
        return None

def format_polymer_info(config):
    """
    Format polymer configuration into a readable string for result.txt
    
    Args:
        config (dict): Polymer configuration
        
    Returns:
        str: Formatted polymer information string
    """
    if config is None:
        return "N/A N/A N/A N/A"
    
    polymer_type = config.get('polymer_type', 'N/A')
    ppta_count = config.get('ppta_count', 'N/A')
    oda_count = config.get('oda_count', 0) if config.get('oda_count') is not None else 0
    cation_count = config.get('cation_count', 'N/A')
    
    return f"{polymer_type} {ppta_count} {oda_count} {cation_count}"

def display_configuration_summary(canonical, solvent_canonical, polymer_config):
    """
    Display simplified configuration summary
    
    Args:
        canonical (str): Canonical monomer SMILES
        solvent_canonical (str): Canonical solvent SMILES
        polymer_config (dict): Polymer configuration dictionary
    """
    print("\n" + "="*72)
    print("    Simulation System Configuration Summary")
    print("="*72)
    print(f"Monomer SMILES: {canonical}")
    print(f"Solvent SMILES: {solvent_canonical}")
    print(f"PPTA[TPC+PPD] Count: {polymer_config['ppta_count']}")
    print(f"[DCA+PPD] Count: {polymer_config['cation_count']}")
    print("="*72)

def check_existing_result(monomer_smiles, solvent_smiles):
    """
    Check for existing results - modified to only check monomer since solvent is fixed
    New format: DCA_SMILES stretched_interE #ofHbond pi_stacking_energy H_bond_interE solution_interE
    """
    result_file_path = os.path.join(os.getcwd(), "result.txt")
    
    if not os.path.exists(result_file_path):
        return None
    
    try:
        with open(result_file_path, 'r') as f:
            lines = f.readlines()
        
        # Parse each line to find matching system
        for line_num, line in enumerate(lines, 1):
            parts = line.strip().split()
            if len(parts) >= 6:  # New format has 6 columns (removed solvent)
                saved_monomer = parts[0]
                saved_stretched_interE = parts[1]
                saved_hbond_count = parts[2]
                saved_pi_stacking = parts[3]
                saved_hbond_interE = parts[4]
                saved_solution_interE = parts[5]
                
                # Canonicalize the saved SMILES for comparison
                try:
                    saved_monomer_canonical = canonicalize_smiles(saved_monomer)
                    
                    # Skip if canonicalization fails
                    if saved_monomer_canonical is None:
                        print(f"Warning: Invalid SMILES found in result.txt line {line_num}, skipping...")
                        continue
                    
                    # Compare canonical forms (only monomer since solvent is fixed)
                    if saved_monomer_canonical == monomer_smiles:
                        return {
                            'monomer_smiles': saved_monomer_canonical,
                            'stretched_interE': saved_stretched_interE,
                            'hbond_count': saved_hbond_count,
                            'pi_stacking_energy': saved_pi_stacking,
                            'hbond_interE': saved_hbond_interE,
                            'solution_interE': saved_solution_interE,
                            'line': line.strip(),
                            'original_monomer': saved_monomer
                        }
                except Exception as e:
                    print(f"Warning: Error processing SMILES in line {line_num}: {e}")
                    continue
    
    except Exception as e:
        print(f"Error reading result.txt: {e}")
        return None
    
    return None

def display_existing_result(result_data):
    """
    Display existing result - modified to only show monomer since solvent is fixed
    """
    print("\n" + "="*80)
    print("FOUND EXISTING SIMULATION RESULT!")
    print("="*80)
    
    # Show both original and canonical SMILES if different
    if 'original_monomer' in result_data and result_data['original_monomer'] != result_data['monomer_smiles']:
        print(f"Monomer SMILES (saved): {result_data['original_monomer']}")
        print(f"Monomer SMILES (canonical): {result_data['monomer_smiles']}")
    else:
        print(f"Monomer SMILES: {result_data['monomer_smiles']}")
    
    print(f"Solvent SMILES: {canonicalize_smiles(FIXED_SOLVENT_SMILES)} (fixed)")
    
    print("-"*80)
    print(f"Stretched interE [kcal/(mol·A^3)] : {result_data['stretched_interE']}")
    print(f"Number of H-bonds per Vol. : {result_data['hbond_count']}")
    print(f"pi-pi stacking energy [kcal/(mol·A^3)] : {result_data['pi_stacking_energy']}")
    print(f"H-bond interE [kcal/(mol·A^3)] : {result_data['hbond_interE']}")
    print(f"Solution interE [kcal/(mol·A^3)] : {result_data['solution_interE']}")
    print("="*80)
    
    # Check if plot exists and display path - only use monomer SMILES for filename
    plot_filename = f"{result_data['monomer_smiles']}.png"
    plot_path = os.path.join("./Result_plot/", plot_filename)
    
    # If not found, try with original SMILES
    if not os.path.exists(plot_path) and 'original_monomer' in result_data:
        plot_filename = f"{result_data['original_monomer']}.png"
        plot_path = os.path.join("./Result_plot/", plot_filename)
    
    # If still not found, try with sanitized filename
    if not os.path.exists(plot_path):
        sanitized_filename = sanitize_filename(result_data['monomer_smiles'])
        plot_filename = f"{sanitized_filename}.png"
        plot_path = os.path.join("./Result_plot/", plot_filename)
    
    if os.path.exists(plot_path):
        print(f"Plot available at: {plot_path}")
    else:
        print("No plot file found for this system")


def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=True)

def save_smiles(smiles, folder="set", filename="name.txt"):
    path = os.path.join(os.getcwd(), folder)
    os.makedirs(path, exist_ok=True)
    filepath = os.path.join(path, filename)
    try:
        with open(filepath, 'w') as f:
            f.write(smiles)
        print(f"Saved canonical SMILES to {filepath}")
        return filepath
    except Exception as e:
        print(f"Error saving SMILES: {e}")
        return None

def save_solvent_smiles(smiles, folder="Solution/solvent", filename="solvent.smi"):
    path = os.path.join(os.getcwd(), folder)
    os.makedirs(path, exist_ok=True)
    filepath = os.path.join(path, filename)
    try:
        with open(filepath, 'w') as f:
            f.write(smiles)
        print(f"Saved solvent SMILES to {filepath}")
    except Exception as e:
        print(f"Error saving solvent SMILES: {e}")

def run_simulation(name_file_path, monomer_smiles, solvent_smiles):
    set_dir = os.path.dirname(name_file_path)
    py_path = os.path.abspath(os.path.join(set_dir, "../Util/Util_Polymer_combine.py"))

    if not os.path.isfile(py_path):
        print(f"Error: Python script not found at {py_path}")
        return

    try:
        subprocess.run(["python3", py_path], cwd=set_dir, check=True)
        print("Initial simulation (combine) executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during initial simulation execution: {e}")

    copy_structures_to_targets(set_dir, monomer_smiles, solvent_smiles)

def copy_structures_to_targets(set_dir, monomer_smiles, solvent_smiles):
    source = os.path.join(set_dir, "structures")
    if not os.path.exists(source):
        print(f"Source directory {source} does not exist. Skipping copy.")
        return

    base_dir = os.getcwd()
    targets = ["Stretched", "Solution"]
    for target in targets:
        target_struct_path = os.path.join(base_dir, target, "structures")
        name_file_path = os.path.join(base_dir, target, "name.txt")

        if os.path.exists(target_struct_path):
            try:
                shutil.rmtree(target_struct_path)
                print(f"Removed existing {target_struct_path}")
            except Exception as e:
                print(f"Failed to remove {target_struct_path}: {e}")
                continue

        try:
            shutil.copytree(source, target_struct_path)
            print(f"Copied {source} → {target_struct_path}")
        except Exception as e:
            print(f"Failed to copy to {target_struct_path}: {e}")
            continue

        try:
            structure_files = sorted(os.listdir(source))
            with open(name_file_path, 'w') as f:
                for fname in structure_files:
                    name_without_ext = os.path.splitext(fname)[0]
                    f.write(name_without_ext + '\n')
            print(f"Wrote structure base names to {name_file_path}")
        except Exception as e:
            print(f"Failed to write structure filenames: {e}")

    run_final_stretch(base_dir)
    run_final_solution(base_dir, monomer_smiles, solvent_smiles)

def run_final_stretch(base_dir):
    stretch_dir = os.path.join(base_dir, "Stretched")
    py_path = os.path.abspath(os.path.join(stretch_dir, "../Util/Util_Polymer_run_Stretched.py"))

    if not os.path.isfile(py_path):
        print(f"Error: Final stretching Python script not found at {py_path}")
        return

    try:
        subprocess.run(["python3", py_path], cwd=stretch_dir, check=True)
        print("Final stretching simulation executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during final stretching simulation: {e}")

def run_final_solution(base_dir, monomer_smiles, solvent_smiles):
    solution_dir = os.path.join(base_dir, "Solution")
    py_path = os.path.abspath(os.path.join(solution_dir, "../Util/Util_Polymer_run_Solution.py"))

    if not os.path.isfile(py_path):
        print(f"Error: Final solution Python script not found at {py_path}")
        return

    try:
        subprocess.run(["python3", py_path], cwd=solution_dir, check=True)
        print("Final solution simulation executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error during final solution simulation: {e}")

    print_interaction_energies(base_dir, monomer_smiles, solvent_smiles)

# === Analysis Functions ===

def read_lammps_data(filename):
    """
    Read LAMMPS output file and return data
    
    Args:
        filename (str): Path to the file to read
    
    Returns:
        pandas.DataFrame: DataFrame containing TimeStep and Value columns
    """
    try:
        data = pd.read_csv(filename, sep=r'\s+', comment='#', 
                          names=['TimeStep', 'Value'])
        return data
    except FileNotFoundError:
        print(f"File not found: {filename}")
        return None
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

def sanitize_filename(filename):
    """
    Sanitize filename by replacing invalid characters with safe alternatives
    """
    replacements = {
        '/': '_','\\': '_',':': '_','*': '_','?': '_','"': '_','<': '_','>': '_','|': '_',
        '(': '[',')': ']','=': '-','+': 'plus','#': 'hash'
    }
    sanitized = filename
    for old_char, new_char in replacements.items():
        sanitized = sanitized.replace(old_char, new_char)
    sanitized = ''.join(c for c in sanitized if c.isalnum() or c in '-_[].')
    if len(sanitized) > 200:
        sanitized = sanitized[:200]
    return sanitized

def plot_combined_figure(stretched_data, solution_data, monomer_smiles=None, solvent_smiles=None, save_plot=True):
    """
    Display Stretched and Solution data in 2x1 subplots
    """
    polymer_config = load_polymer_config()
    polymer_info = format_polymer_info(polymer_config)
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['axes.unicode_minus'] = False
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 14))
    plt.subplots_adjust(hspace=0.4, top=0.92, bottom=0.08)
    
    color1 = 'tab:blue'
    color2 = 'tab:red'
    
    ax1_twin = ax1.twinx()
    ax1.plot(stretched_data['timesteps_ps'], stretched_data['inter_energy'], '-', color=color1, linewidth=2.5, alpha=0.8, label='Chain-Chain Interaction Energy')
    ax1.set_xlabel('Time [ps]', fontsize=13)
    ax1.set_ylabel('Chain-Chain interE [Kcal/mol·A^3]', color=color1, fontsize=12)
    ax1.tick_params(axis='y', labelcolor=color1, labelsize=10)
    ax1.tick_params(axis='x', labelsize=10)
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1_twin.plot(stretched_data['timesteps_ps'], stretched_data['density'], '-', color=color2, linewidth=2.5, alpha=0.8, label='Density')
    ax1_twin.set_ylabel('Density [g/cm³]', color=color2, fontsize=12)
    ax1_twin.tick_params(axis='y', labelcolor=color2, labelsize=10)
    ax1.set_title('Stretched State', fontsize=15, fontweight='bold', pad=15)
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1_twin.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=10)
    
    ax2_twin = ax2.twinx()
    ax2.plot(solution_data['timesteps_ps'], solution_data['inter_energy'], '-', color=color1, linewidth=2.5, alpha=0.8, label='Chain-Chain Interaction Energy')
    ax2.set_xlabel('Time [ps]', fontsize=13)
    ax2.set_ylabel('Chain-Chain interE [Kcal/mol·A^3]', color=color1, fontsize=12)
    ax2.tick_params(axis='y', labelcolor=color1, labelsize=10)
    ax2.tick_params(axis='x', labelsize=10)
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2_twin.plot(solution_data['timesteps_ps'], solution_data['density'], '-', color=color2, linewidth=2.5, alpha=0.8, label='Density')
    ax2_twin.set_ylabel('Density [g/cm³]', color=color2, fontsize=12)
    ax2_twin.tick_params(axis='y', labelcolor=color2, labelsize=10)
    ax2.set_title('Solution State', fontsize=15, fontweight='bold', pad=15)
    lines3, labels3 = ax2.get_legend_handles_labels()
    lines4, labels4 = ax2_twin.get_legend_handles_labels()
    ax2.legend(lines3 + lines4, labels3 + labels4, loc='upper right', fontsize=10)
    
    fig.suptitle('LAMMPS Time Series: Stretched vs Solution States\n', fontsize=17, fontweight='bold')
    
    if save_plot:
        output_dir = './Result_plot/'
        os.makedirs(output_dir, exist_ok=True)
        if monomer_smiles:
            # Use only monomer SMILES for filename since solvent is fixed
            filename = f"{monomer_smiles}.png"
        else:
            filename = 'stretched_vs_solution_combined.png'
        output_path = os.path.join(output_dir, filename)
        try:
            plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
            plt.show()
            return output_path
        except Exception as e:
            print(f"Warning: Could not save with original SMILES filename. Using sanitized version.")
            print(f"Error: {e}")
            sanitized_filename = sanitize_filename(monomer_smiles)
            filename = f"{sanitized_filename}.png"
            output_path = os.path.join(output_dir, filename)
            plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
            plt.show()
            return output_path
    plt.show()
    return None

def load_and_process_data(data_dir_name):
    """
    Load and process LAMMPS data from specific directory
    """
    data_dir = f"./{data_dir_name}/lammps/monomer_0/"
    output1_file = os.path.join(data_dir, "output1.txt")
    output5_file = os.path.join(data_dir, "output5.txt")
    
    inter_energy_data = read_lammps_data(output1_file)
    density_data = read_lammps_data(output5_file)
    
    if inter_energy_data is None or density_data is None:
        return None
    
    timesteps_raw = inter_energy_data['TimeStep'].values
    inter_energy = inter_energy_data['Value'].values
    density = density_data['Value'].values
    
    timesteps_ps = timesteps_raw * 0.5 / 1000  # ps units
    
    return {
        'timesteps_raw': timesteps_raw,
        'timesteps_ps': timesteps_ps,
        'inter_energy': inter_energy,
        'density': density,
        'name': data_dir_name
    }

def run_analysis(monomer_smiles=None, solvent_smiles=None):
    """
    Analyze simulation results and generate graphs
    """
    stretched_data = load_and_process_data("Stretched")
    if stretched_data is None:
        return None
    
    solution_data = load_and_process_data("Solution")
    if solution_data is None:
        return None
    
    try:
        plot_path = plot_combined_figure(stretched_data, solution_data, 
                                       monomer_smiles=monomer_smiles, 
                                       solvent_smiles=None,  # Not needed since solvent is fixed
                                       save_plot=True)
        return plot_path
    except Exception as e:
        return None

def read_last_energy(filepath):
    """
    Read the second column value from the last line of a file
    
    Args:
        filepath (str): Path to the file to read
        
    Returns:
        float or None: Second column value from the last line, None if failed
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
            if not lines:
                return None
            
            # Get the last line
            last_line = lines[-1].strip()
            if not last_line:
                # If last line is empty, check previous lines
                for line in reversed(lines):
                    if line.strip():
                        last_line = line.strip()
                        break
            
            parts = last_line.split()
            if len(parts) >= 2:
                return float(parts[1])
                
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return None

def read_stretched_data(base_dir):
    """
    Read all required data from Stretched directory
    
    Args:
        base_dir (str): Base directory path
        
    Returns:
        dict: Dictionary containing data from each file
    """
    stretched_dir = os.path.join(base_dir, "Stretched", "lammps", "monomer_0")
    
    data_files = {
        'stretched_interE': 'output1.txt',
        'hbond_count': 'output2.txt', 
        'pi_stacking_energy': 'output3.txt',
        'hbond_interE': 'output4.txt'
    }
    
    results = {}
    
    for key, filename in data_files.items():
        filepath = os.path.join(stretched_dir, filename)
        value = read_last_energy(filepath)
        results[key] = value
    
    return results

def read_first_energy(filepath):
    """
    Read the second column value from the first line of a file (for Solution data)
    
    Args:
        filepath (str): Path to the file to read
        
    Returns:
        float or None: Second column value from the first line, None if failed
    """
    try:
        with open(filepath, 'r') as f:
            line = f.readline()
            parts = line.strip().split()
            if len(parts) >= 2:
                return float(parts[1])
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    return None

def print_interaction_energies(base_dir, monomer_smiles, solvent_smiles):
    """
    Read all interaction energies and record to result.txt
    New format: DCA_SMILES stretched_interE #ofHbond pi_stacking_energy H_bond_interE solution_interE
    """
    # Read Stretched data
    stretched_data = read_stretched_data(base_dir)
    
    # Read Solution data (maintain existing approach)
    solution_ehbond = os.path.join(base_dir, "Solution", "E_h_bond.txt")
    solution_value = read_first_energy(solution_ehbond)  # Use existing function
    
    # Display 5 values in cmd window (maintain existing format)
    stretched_display = f"{stretched_data['stretched_interE']}" if stretched_data['stretched_interE'] is not None else "N/A"
    hbond_count_display = f"{stretched_data['hbond_count']}" if stretched_data['hbond_count'] is not None else "N/A"
    pi_stacking_display = f"{stretched_data['pi_stacking_energy']}" if stretched_data['pi_stacking_energy'] is not None else "N/A"
    hbond_interE_display = f"{stretched_data['hbond_interE']}" if stretched_data['hbond_interE'] is not None else "N/A"
    solution_display = f"{solution_value}" if solution_value is not None else "N/A"
    
    print(f"Stretched interE [kcal/(mol·A^3)] : {stretched_display}")
    print(f"Number of Hbond per Vol. : {hbond_count_display}")
    print(f"pi-pi stacking energy [kcal/(mol·A^3)] : {pi_stacking_display}")
    print(f"H-bond interE [kcal/(mol·A^3)] : {hbond_interE_display}")
    print(f"Solution interE [kcal/(mol·A^3)] : {solution_display}")
    
    # Convert each value to string for result.txt (handle None cases as "N/A")
    stretched_interE = f"{stretched_data['stretched_interE']:.6f}" if stretched_data['stretched_interE'] is not None else "N/A"
    hbond_count = f"{stretched_data['hbond_count']:.6f}" if stretched_data['hbond_count'] is not None else "N/A"
    pi_stacking = f"{stretched_data['pi_stacking_energy']:.6f}" if stretched_data['pi_stacking_energy'] is not None else "N/A"
    hbond_interE = f"{stretched_data['hbond_interE']:.6f}" if stretched_data['hbond_interE'] is not None else "N/A"
    solution_interE = f"{solution_value:.6f}" if solution_value is not None else "N/A"
    
    # Create result line for result.txt - only include monomer SMILES since solvent is fixed
    # New format: DCA_SMILES stretched_interE #ofHbond pi_stacking_energy H_bond_interE solution_interE
    result_line = f"{monomer_smiles} {stretched_interE} {hbond_count} {pi_stacking} {hbond_interE} {solution_interE}"
    
    # Append to result.txt
    result_file_path = os.path.join(os.getcwd(), "result.txt")
    result_relative_path = "./result.txt"
    try:
        with open(result_file_path, "a") as f:
            f.write(result_line + "\n")
        print(f"Results saved to {result_relative_path}")
    except Exception as e:
        print(f"Error writing to result.txt: {e}")
    
    # Generate analysis plots - only pass monomer SMILES
    plot_path = run_analysis(monomer_smiles=monomer_smiles, solvent_smiles=None)
    if plot_path:
        print(f"Results plots saved to {plot_path}")
    else:
        print("Error: Could not generate analysis plots.")

def main():
    # Get canonical form of the fixed solvent
    solvent_canonical = canonicalize_smiles(FIXED_SOLVENT_SMILES)
    if solvent_canonical is None:
        print(f"Error: Fixed solvent SMILES '{FIXED_SOLVENT_SMILES}' is invalid!")
        return
    
    while True:
        smiles_input = input("Enter a dicarboxylic acid(DCA) monomer SMILES string (or 'q' to quit): ")
        if smiles_input.lower() == 'q':
            break
        
        # Clean directories before processing new SMILES
        clean_set_directory()
        clean_stretched_directory()
        clean_solution_directory()

        canonical = canonicalize_smiles(smiles_input)

        if canonical is None:
            print("Invalid monomer SMILES. Please check the input.")
            continue

        # Get fixed polymer configuration
        polymer_config = get_polymer_configuration()
        
        # Display simplified configuration summary with fixed solvent
        display_configuration_summary(canonical, solvent_canonical, polymer_config)
        
        existing_result = check_existing_result(canonical, solvent_canonical)
        if existing_result:
            display_existing_result(existing_result)
            user_choice = input("\nDo you want to re-run the simulation anyway? (y/n): ").lower()
            if user_choice != 'y':
                continue
            else:
                print("\nRe-running simulation as requested...")
        
        name_file_path = save_smiles(canonical)
        if not name_file_path:
            continue

        # Save the fixed solvent SMILES
        save_solvent_smiles(solvent_canonical)
        
        if not save_polymer_config(polymer_config):
            print("Failed to save polymer configuration. Continuing anyway...")
        
        run_simulation(name_file_path, canonical, solvent_canonical)

if __name__ == "__main__":
    main()
