# AramidSim_DCA : Polymer Simulation Interface Tool (DCA-based)

This repository provides a **Python-based simulation interface** for studying polymer systems derived from **dicarboxylic acids (DCAs)**.  
Given a **DCA monomer SMILES**, the tool automatically performs **LAMMPS simulations** in both **stretched (fiber)** and **solution** states using **NMP (N-Methyl-2-pyrrolidone)** as the fixed solvent, extracts interaction energies, and generates comparative plots.

The polymer chain consists of **8 monomer units**: 4 PPTA[TPC+PPD] units and 4 [DCA+PPD] units.

- **For academic use**: reproducible workflows for polymer–polymer & polymer–solvent interactions and structure–property relationship studies.  
- **For open-source users**: simple input/output interface, automatic result logging, and visualization.  

## Polymer System Configuration

The polymer chain structure consists of **8 monomer units**:
- **4 PPTA[TPC+PPD] units**: Para-aramid backbone segments
- **4 [DCA+PPD] units**: Dicarboxylic acid-based segments

### Simulation States
- **Stretched (Fiber) State**: 25 infinitely connected polymer chains arranged in a bundle to simulate fiber structure
- **Solution State**: H-terminated polymer chains dissolved in **NMP solvent** at 9 wt% concentration

---

## Features
- **End-to-End Automation**: DCA SMILES input → structure generation → stretched & solution simulations → results → figures  
- **Dual Simulation States**: 
  - **Stretched state**: 25 infinitely connected polymer chains forming a fiber bundle
  - **Solution state**: H-terminated polymer chains at 9 wt% concentration in NMP solvent
- **Fixed Solvent System**: Uses NMP (O=C1CCCN1C) as the standardized solvent for consistent comparisons
- **Pure Python Workflow**: Python orchestration + LAMMPS/Moltemplate backend  
- **Automated Analysis & Visualization**: Extracts chain–chain interaction energies (interE) and densities, with comparative plots  

---

## Repository Structure
```
.
└── AramidSim_DCA/
    ├── README.md                           # This documentation file
    ├── Simulation.py                       # Main simulation interface
    ├── set_up.py                           # Set up Moltemplate & LAMMPS and helper programs
    ├── result.txt                          # Accumulated results
    ├── Util/                               # Python helper modules & utilities
    ├── Stretched/                          # Stretched-state simulations
    ├── Solution/                           # Solution-state simulations
    ├── Result_plot/                        # Auto-generated figures
    └── set/                                # Stores polymer configuration & installation files
```

---

## Environment & Dependencies

This project was developed and tested in the following environment:

- **Python**: 3.12.7  
- **LAMMPS**: lammps-2Aug2023 (https://www.lammps.org)
- **Moltemplate**: 2.19.12 (https://www.moltemplate.org)
- **AmberTools**: 24.8 (https://ambermd.org)
- **Conda environment**: `Aramid` (with RDKit, numpy, pandas, matplotlib, scipy, etc.)  

### Installation (Recommended)
```bash
# Create environment
conda create -n Aramid_Interface python=3.12.7 -y
conda activate Aramid_Interface

# Install core dependencies
conda install -c conda-forge numpy scipy pandas matplotlib ambertools -y
pip install rdkit-pypi moltemplate
```

---

## Installation & Setup

### 1. Clone or download the repository
```bash
git clone <repository URL>
cd AramidSim_DCA
```

### 2. Install Moltemplate & LAMMPS
```bash
python set_up.py
```
- This script handles the installation and setup of Moltemplate and LAMMPS.  
- Files will be placed into the `set/` folder for internal use.  

If you **already have Moltemplate and LAMMPS** installed system-wide:  
- You may still need to run the setup script to ensure proper configuration.  

---

## Running a Simulation

### Run the main interface
```bash
python Simulation.py
```

### Example Input
```
Enter a dicarboxylic acid(DCA) monomer SMILES string (or 'q' to quit): O=C(O)c1cncc(C(=O)O)c1
```

**Note**: The solvent is fixed as NMP (N-Methyl-2-pyrrolidone, SMILES: O=C1CCCN1C) for all simulations.

### Example Output
```
Stretched interE [kcal/(mol·A^3)] : -0.131262
Number of H-bonds per Vol. : 0.006473
pi-pi stacking energy [kcal/(mol·A^3)] : -0.047482
H-bond interE [kcal/(mol·A^3)] : -0.118610
Solution interE [kcal/(mol·A^3)] : -0.212430

Plot available at: ./Result_plot/O=C(O)c1cncc(C(=O)O)c1.png
```

---

## Results

- **result.txt format**:
  ```
  Monomer_SMILES  stretched_interE  H-bonds_per_Vol  pi-pi_stacking_energy  H-bond_interE  solution_interE
  ```
  Example:
  ```
  O=C(O)c1cncc(C(=O)O)c1 -0.131262 0.006473 -0.047482 -0.118610 -0.212430
  ```

- **Example figure** (from `Result_plot/`):  
  - Filename format: `[DCA_SMILES].png`
  - Contains comparative plots of interaction energies between stretched and solution states

---

## Academic Notes

This tool provides a reproducible pipeline for computational polymer science:  
- Automated **polymer generation** from DCA SMILES  
- Dual-state **LAMMPS simulations** (fiber vs solution) with standardized NMP solvent
- Extraction of **interaction energies & densities**  
- Archiving in `result.txt` and comparative plots in `Result_plot/`  
- **Consistent solvent system** enables direct comparison across different DCA monomers

---

## Citation
If you use this tool in academic work, please cite as:

> *Polymer Simulation Interface Tool (DCA-based).* GitHub repository, 2025.  
> URL: `<repository URL>`  

### BibTeX
```bibtex
@misc{PolymerSimulationDCA2025,
  author       = {Your Name and Collaborators},
  title        = {Polymer Simulation Interface Tool (DCA-based)},
  year         = {2025},
  publisher    = {GitHub},
  journal      = {GitHub repository},
  howpublished = {\url{<repository URL>}},
}
```

---

## License
MIT License (or specify otherwise).
