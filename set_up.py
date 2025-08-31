#!/usr/bin/env python3
import os
import subprocess

def run_command(command, cwd=None, description=""):
    """Execute a shell command with basic logging and error handling."""
    try:
        print(f"Running: {description if description else command}")
        subprocess.run(command, shell=True, check=True, cwd=cwd)
        print(f"✅ Success: {description if description else command}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error in {description if description else command}: {e}")
        return False

def get_user_choice(prompt):
    """Get user's y/n choice with input validation."""
    while True:
        choice = input(f"{prompt} (y/n): ").lower().strip()
        if choice in ['y', 'yes']:
            return True
        elif choice in ['n', 'no']:
            return False
        else:
            print("Please enter 'y' for yes or 'n' for no.")

# ✅ LAMMPS packages to enable (order matters)
LMP_PACKAGES_TO_ENABLE_IN_ORDER = [
    "body",
    "class2",
    "kspace",
    "dielectric",
    "dipole",
    "molecule",
    "extra-compute",
    "extra-fix",
    "extra-molecule",
    "extra-pair",
    "fep",
    "granular",
    "manifold",
    "manybody",
    "mc",
    "misc",
    "qeq",
    "reaction",
    "reaxff",
    "rigid",
]

def enable_only_selected_packages_in_order(lammps_src_path: str):
    print("\n=== Step 2A: Configure LAMMPS packages (enable only selected, in order) ===")

    # 1) Disable all packages
    if not run_command("make no-all", cwd=lammps_src_path, description="Disable all LAMMPS packages (make no-all)"):
        return False

    # 2) Enable only the desired packages in the specified order
    for pkg in LMP_PACKAGES_TO_ENABLE_IN_ORDER:
        if not run_command(f"make yes-{pkg}", cwd=lammps_src_path, description=f"Enable package in order: {pkg}"):
            return False

    # 3) Clean before build (optional)
    run_command("make clean-all", cwd=lammps_src_path, description="Clean previous LAMMPS builds (make clean-all)")
    return True

def install_moltemplate():
    """Install moltemplate if user chooses to."""
    print("\n=== Step 1: Installing moltemplate ===")
    
    if not get_user_choice("Do you want to install moltemplate?"):
        print("Skipping moltemplate installation.")
        return True
    
    base_util_path = os.path.join(os.getcwd(), "Util")
    moltemplate_path = os.path.join(base_util_path, "moltemplate-master")
    
    if os.path.exists(moltemplate_path):
        print(f"Found moltemplate directory: {moltemplate_path}")
        return run_command("pip3 install .", cwd=moltemplate_path, description="Installing moltemplate")
    else:
        print(f"Warning: moltemplate directory not found at {moltemplate_path}")
        return False

def install_lammps():
    """Configure and compile LAMMPS if user chooses to."""
    print("\n=== Step 2: Configuring & Compiling LAMMPS ===")
    
    if not get_user_choice("Do you want to configure and compile LAMMPS?"):
        print("Skipping LAMMPS installation.")
        return True
    
    base_util_path = os.path.join(os.getcwd(), "Util")
    lammps_src_path = os.path.join(base_util_path, "lammps-2Aug2023", "src")
    
    if not os.path.exists(lammps_src_path):
        print(f"Warning: LAMMPS src directory not found at {lammps_src_path}")
        return False
    
    print(f"Found LAMMPS src directory: {lammps_src_path}")

    # 2A) Package configuration (fixed order)
    if not enable_only_selected_packages_in_order(lammps_src_path):
        print("Aborting due to package configuration failure.")
        return False

    # Ask user which version(s) to compile
    compile_serial = get_user_choice("Do you want to compile LAMMPS serial version?")
    compile_mpi = get_user_choice("Do you want to compile LAMMPS MPI version?")
    
    success = True
    
    # 2B) Serial build
    if compile_serial:
        print("\nCompiling LAMMPS serial version...")
        if not run_command("make serial", cwd=lammps_src_path, description="Making LAMMPS (serial)"):
            success = False

    # 2C) MPI build
    if compile_mpi:
        print("\nCompiling LAMMPS MPI version...")
        if not run_command("make mpi", cwd=lammps_src_path, description="Making LAMMPS (mpi)"):
            success = False
    
    if not compile_serial and not compile_mpi:
        print("No LAMMPS version selected for compilation.")
    
    return success

def main():
    print("=== Interactive Setup Process ===")
    print("This script will help you install moltemplate and/or LAMMPS.")
    print("You can choose which components to install.\n")

    # Step 1: Install moltemplate (optional)
    moltemplate_success = install_moltemplate()
    
    # Step 2: Configure & Compile LAMMPS (optional)
    lammps_success = install_lammps()

    # Summary
    print("\n=== Setup Summary ===")
    if moltemplate_success:
        print("✅ Moltemplate step completed successfully")
    else:
        print("❌ Moltemplate step failed or was skipped")
    
    if lammps_success:
        print("✅ LAMMPS step completed successfully")
    else:
        print("❌ LAMMPS step failed or was skipped")
    
    print("\n=== Setup completed ===")

if __name__ == "__main__":
    main()
