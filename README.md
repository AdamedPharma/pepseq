# pepseq
Chemoinformatic library to work flawlessly with **modified peptide** structures.

## Aim:
A **text format** to conveniently handle complicated modified peptide structures as standard formats like SMILES or Fasta sequence are either to simple or to complicated.

**New format** is a combination of:

- Easily readable 1aa code like : 
KYLDERAAQDFVQW

- 3aa code that allows to handle non-standard amino acids: 
Ala-Aib-Gly-Lys

- Common chemical names to define typical chemical substitution of amino acids, like: 
Thr(tBu), Ser(tBu), Phe(4-Cl)

- SMILES codes to describe unspecified modifications: 
KYLDCRAAQDFVQW<[5*]CCCC(N)=O>

## Installation
To use Pepseq, first install it using pip:
pip install pepseq

## Usage

python3.10 commands.py  pepseq-to-smiles 'CDEFGAHC' --out outpath.smi

python3.10 commands.py read-smiles outpath.smi --out smiles_reading_output

python3.10 commands.py  calculate-json-from  'CH3~SC{Cys(R1)}AFC~NH2' --mod-smiles '[*:1]CCC' --ketcher

python3.10 commands.py augment-db-json-command --sdf-path monomers.sdf

ls *smiles_reading_output*
smiles_reading_output.pepseq  smiles_reading_output.smi
## Testing

python3.10 -m pytest

### Arguments
out.smi <- modification SMILES


## Roadmap
- Add handling of more external modifications, e.g. staples (like ornithine), linkers (like PEG2), anchors (like C‐18 fatty di‐acid chain).
- Add decomposing of peptide modification where it is compounded, e.g. peptide stapled with ornithine which is linked through PEG2 to
C-18 fatty di-acid chain.

- Naming more atoms in residues SMARTS.

- Add more fixtures for testing.
## Documentation

to generate documentation

cd docs

make html

documentation in HTML format will be available in docs/build/html/index.html
