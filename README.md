# PDB Disulfide Bond Finder 

----

## Overview 

**Disulfide Bond Finder** is a Python script designed to identify potential disulfide bonds within protein structures from PDB files. It extracts cysteine residues and calculates distances and dihedral angles between them in order to predict potential disulfide bonds.

## Features 

- **Residue Identification**: Identifies cysteine residues in the protein structure.

- **Region Filtering**: Filters out potential mobile regions based on:
    - **B-factor**: for experimental structures (if B-factor >30 A²)
    - **pLDDT**: for AlphaFold models (if pLDDT <50)

- **Disulfide Bond Prediction**: Predicts potential disulfide bonds based on geometric criteria:
    - **Distance**: 1.5 - 2.5 Å
    - **Dihedral Angle**: 84° - 96°

## Requirements

- Python 3.X

- Biopython (Bio.PDB module)

Install Biopython via pip:

```bash
pip install biopython
```

Or via conda:

```bash
conda install -c conda-forge biopython
```

## Usage

Run the script from the command line with a PDB file as an argument:

```bash
python pdb_diS.py <example_protein.pdb>
```

If no argument is provided, the program will prompt to enter the file name:

```bash
python pdb_diS.py
```
```
Enter file name (PDB): example_protein.pdb
```

### Example Output

**When disulfide bonds are found:**

```
Potential disulfide bonds found in 'example_protein.pdb': 
* CYS 148 - CYS 191 | Distance: 2.04 Å | Angle: 91.12º
* CYS 182 - CYS 224 | Distance: 2.04 Å | Angle: 95.61º
* CYS 282 - CYS 332 | Distance: 2.03 Å | Angle: 85.46º
```
**When no disulfide bonds are found:**

```
Potential disulfide bonds NOT found in 'example_protein.pdb'.
```

## References 

- [Biopython Structural Bioinformatics Documentation](https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ)
