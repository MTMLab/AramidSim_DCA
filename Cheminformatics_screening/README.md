# AramidSim_DCA_cheminformatics_screening

Deterministic screening of dicarboxylic acid (DCA) monomers for aramid polymer design.  
This repository provides `screening_monomers.py`, which implements the same filtering logic described in our related publication to ensure consistency and reproducibility.

---

## Screening Criteria

- **Aromatic core**: We required the presence of at least one benzene ring in the monomer structure. Aromatic rings impart rigidity and strength to aramid backbones (as seen in PPTA). This filter reduced the list to those compounds containing an aromatic ring and the two carboxylic acid groups.
- **Moderate structural complexity**: We eliminated molecules with a chemical complexity > **800** (PubChem complexity). Extremely complex molecules often indicate large sizes or many functional groups, which can complicate polymerization and introduce side reactions or processing issues. By capping complexity, we favored monomers that are synthetically accessible and more likely to polymerize cleanly.
- **Neutral charge**: Only neutral molecules were retained. Any candidate containing ionic groups (identified by “+” or “–” in the SMILES string) was removed. Charged functional groups could ionize during polymerization or fiber spinning, disrupting the polycondensation reaction or leading to undesirable salt formation.
- **Linear geometry**: To promote the formation of highly linear, crystalline polymer chains, we computed the angle formed by lines connecting the molecular center of mass to each of the two carboxylate carbon atoms. Monomers closer to **180°** are more linear; we required this angle to be **≥ 90°**.
- **Limited radius of gyration**: To avoid bulky substituents that hinder close packing, we excluded monomers whose maximum perpendicular thickness (atom radius of gyration about the principal axis connecting the two carboxyl carbons) exceeded **6 Å**.

---

## Quick Start

### 1) Install (recommended via conda)
```bash
# Create a fresh conda env (recommended)
conda create -n dca-screen python=3.10 -y
conda activate dca-screen

# Install RDKit via conda-forge
conda install -c conda-forge rdkit -y

```


### 2) Prepare input
Place an input CSV (default name: `DCA_data.csv`) with at least the following columns:
- `isosmiles` — SMILES string for each monomer
- `complexity` — PubChem complexity (numeric)


### 3) Run
With default parameters (no arguments):
```bash
python screening_monomers.py
```
This is equivalent to:
```bash
python screening_monomers.py \
  --input DCA_data.csv \
  --smiles_col isosmiles \
  --complexity_col complexity \
  --complexity_max 800 \
  --no-isomeric \
  --angle_min 90.0 \
  --max_perp 6.0 \
  --out DCA_screening_data.csv
```


### 4) Output
- **`DCA_screening_data.csv`** — only the rows that *passed* all screening criteria, preserving the same column order as the input.
- Console will print a summary like:
  ```
  Saved 3385 passed molecules to DCA_screening_data.csv
  ```

---

## Command-Line Options

```
--input            Path to input CSV (default: DCA_data.csv)
--smiles_col       Column name for SMILES (default: isosmiles)
--complexity_col   Column name for complexity (default: complexity)
--complexity_max   Max allowed complexity (default: 800.0)
--isomeric         Use isomeric SMILES (default: false)
--no-isomeric      Do not use isomeric SMILES (default)
--angle_min        Minimum CoM–C···CoM–C angle in degrees (default: 90.0)
--max_perp         Max perpendicular thickness in Å (default: 6.0)
--out              Output CSV path (default: DCA_screening_data.csv)
```

---

## Citing

If you use this code, please cite the associated paper (fill in your details here):
```
Author(s), Title, Journal (Year), DOI
```

---

## License

Specify your project license here (e.g., MIT).
