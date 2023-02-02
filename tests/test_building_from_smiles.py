import json
import pkgutil

from BuildPeptideJSONFromSMILES import decompose_peptide_smiles

db_path = pkgutil.extend_path("Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(db_path)


mol_N_C_smiles_val = "CC(=O)N[C@H]1CSC(Br)CNP([Na])SC[C@@H](C(=O)N[C@@H](C)C(=O)N[C@H]2CSSC[C@@H](C(=O)N[C@@H](CCCCN)C(N)=O)NC(=O)CNC2=O)NC(=O)[C@H](CO)NC1=O"

correct_peptide_json = {
    "sequence": "CSCACGCK",
    "internal_modifications": [
        {
            1: [
                {"ResID": "5", "AtomName": "SG", "ResidueName": ""},
                {"ResID": "7", "AtomName": "SG", "ResidueName": ""},
            ]
        }
    ],
    "external_modifications": [
        {
            "smiles": "[1*]C(Br)CNP([2*])[Na]",
            "max_attachment_point_id": 2,
            "attachment_points_on_sequence": {
                1: {
                    "attachment_point_id": 1,
                    "ResID": "1",
                    "AtomName": "SG",
                    "ResidueName": "",
                },
                2: {
                    "attachment_point_id": 2,
                    "ResID": "3",
                    "AtomName": "SG",
                    "ResidueName": "",
                },
            },
        },
        {
            "smiles": "[1*]O",
            "max_attachment_point_id": 1,
            "attachment_points_on_sequence": {
                1: {
                    "attachment_point_id": 1,
                    "ResID": "8",
                    "AtomName": "CO",
                    "ResidueName": "",
                }
            },
        },
    ],
}


correct_peptide_json_NC = {
    "sequence": "CSCACGCK",
    "internal_modifications": [
        {
            1: [
                {"ResID": "5", "AtomName": "SG", "ResidueName": ""},
                {"ResID": "7", "AtomName": "SG", "ResidueName": ""},
            ]
        }
    ],
    "external_modifications": [
        {
            "smiles": "[2*]C(Br)CNP([3*])[Na]",
            "max_attachment_point_id": 3,
            "attachment_points_on_sequence": {
                2: {
                    "attachment_point_id": 2,
                    "ResID": "1",
                    "AtomName": "SG",
                    "ResidueName": "",
                },
                3: {
                    "attachment_point_id": 3,
                    "ResID": "3",
                    "AtomName": "SG",
                    "ResidueName": "",
                },
            },
        }
    ],
    "C_terminus": "NH2",
    "N_terminus": "Ac",
}


def test_decompose_peptide_smiles_db():
    peptide_json_NC = decompose_peptide_smiles(mol_N_C_smiles_val, db_json)
    assert peptide_json_NC == correct_peptide_json_NC