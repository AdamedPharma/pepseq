import pkgutil
import json
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import mol_json_to_nx


mol_N_C_smiles_val = (
    "CC(=O)N[C@H]1CSC(Br)CNP([Na])SC[C@@H](C(=O)N[C@"
    + "@H](C)C(=O)N[C@H]2CSSC[C@@H](C(=O)N[C@@H](CCCCN)C(N)=O)NC(=O)CNC"
    + "2=O)NC(=O)[C@H](CO)NC1=O"
)


residues_path = pkgutil.extend_path("residues.json", __name__)
with open(residues_path) as fp:
    residues_json = json.load(fp)


residues_graphs = [mol_json_to_nx(residue_json) for residue_json in residues_json]

tests = [(mol_N_C_smiles_val, residues_json)]
