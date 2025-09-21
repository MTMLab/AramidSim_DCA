# DCA_reaxff_conversion_and_simulation
This repository provides a minimal pipeline for performing tensile simulations with ReaxFF in LAMMPS, starting from NPT-equilibrated structures obtained using conventional force fields.
The core workflow involves automatically converting bond-based (Amber) data files into ReaxFF-compatible data files, followed by running the tensile simulation with run.ts (a LAMMPS input script).
ReaxFF naturally handles bond order changes (bond breaking and formation), it is well-suited for analyzing fracture and failure mechanisms.


## Quick Start

1) **Data Conversion**
```bash
python Amber_to_Reaxff.py
```

2) **Run Tensile Simulation (example)**
```bash
lmp_serial -in run.ts 
```
