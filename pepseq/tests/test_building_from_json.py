import json
import pkgutil

import rdkit
from pepseq.BuildingModifiedPeptideFromPeptideJSON import (
    BuildingModifiedPeptideFromPeptideJSON,
)

db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


peptide_json = {
    "sequence": "CSCACGCK",
    "external_modifications": [
        {
            "smiles": "[1*]C([Br])CNP([Na])[2*]",
            "attachment_points_on_sequence": {
                1: {
                    "attachment_point_id": 1,
                    "ResID": "1",
                    "AtomName": "SG",
                    "ResidueName": "CYS",
                },
                2: {
                    "attachment_point_id": 2,
                    "ResID": "3",
                    "AtomName": "SG",
                    "ResidueName": "CYS",
                },
            },
        }
    ],
    "internal_modifications": {
        1: [
            {"ResID": 5, "AtomName": "SG", "ResidueName": "CYS"},
            {"ResID": 7, "AtomName": "SG", "ResidueName": "CYS"},
        ]
    },
}


def test_building():
    mol = BuildingModifiedPeptideFromPeptideJSON().execute(peptide_json, db_json)
    assert (
        rdkit.Chem.MolToSmiles(mol)
        == "[H]N[C@H]1CSC(Br)CNP([Na])SC[C@@H](C(=O)N[C@@H](C)C(=O)N["
        + "C@H]2CSSC[C@@H](C(=O)N[C@@H](CCCCN)C(=O)O)NC(=O)CNC2=O)N"
        + "C(=O)[C@H](CO)NC1=O"
    )
