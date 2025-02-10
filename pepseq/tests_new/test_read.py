import json
import pkgutil

from pepseq.read import from_json


db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


def test_from_json():
    """
    Test the `from_json` function to ensure it correctly parses a peptide JSON object
    and creates a peptide object with the expected attributes.
    The test uses a predefined peptide JSON object `peptide_json_fx_2` and compares
    the resulting peptide object's dictionary representation with the expected dictionary
    `peptide__dict__fx`.
    The expected dictionary includes:
    - `smiles`: The SMILES representation of the peptide.
    - `complete_smiles`: The complete SMILES representation of the peptide.
    - `peptide_json`: The original JSON object used to create the peptide.
    - `sequence`: The peptide sequence.
    - `length`: The length of the peptide.
    - `mw`: The molecular weight of the peptide.
    The test asserts that the dictionary representation of the peptide object matches
    the expected dictionary.
    """
    peptide_json_fx_2 = {
        "length": 17,
        "sequence": "CACDAPEPsEQCGCDEF",
        "internal_modifications": [],
        "C_terminus": "OH",
        "N_terminus": "H",
        "pepseq_format": "H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH",
        "symbols": [
            "H",
            "Cys(R1)",
            "A",
            "C",
            "D",
            "A",
            "P",
            "E",
            "P",
            "s",
            "E",
            "Q",
            "Cys(R2)",
            "G",
            "Cys(R3)",
            "D",
            "E",
            "F",
            "OH",
        ],
        "external_modifications": [
            {
                "smiles": "[*:1]CNCC[*:2]",
                "max_attachment_point_id": 2,
                "attachment_points_on_sequence": {
                    "2": {
                        "attachment_point_id": "2",
                        "ResID": "12",
                        "AtomName": "SG",
                        "ResidueName": "Cys",
                    },
                    "1": {
                        "attachment_point_id": "1",
                        "ResID": "1",
                        "AtomName": "SG",
                        "ResidueName": "Cys",
                    },
                },
            },
            {
                "smiles": "[*:3]CNCC",
                "max_attachment_point_id": 3,
                "attachment_points_on_sequence": {
                    "3": {
                        "attachment_point_id": "3",
                        "ResID": "14",
                        "AtomName": "SG",
                        "ResidueName": "Cys",
                    }
                },
            },
        ],
    }

    peptide = from_json(peptide_json_fx_2, db_json)
    smiles_fx = (
            "[H]N[C@H]1CSCNCCSC[C@@H](C(=O)NCC(=O)N[C@@H](CSCNCC)C"
            "(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H]"
            "(Cc2ccccc2)C(=O)O)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H]"
            "(CCC(=O)O)NC(=O)[C@@H](CO)NC(=O)[C@@H]2CCCN2C(=O)"
            "[C@H](CCC(=O)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](C)NC(=O)"
            "[C@H](CC(=O)O)NC(=O)[C@H](CS)NC(=O)[C@H](C)NC1=O"
    )
    complete_smiles_fx = (
            "[H]N[C@H]1CSCNCCSC[C@@H](C(=O)NCC(=O)N[C@@H](CSCNCC)"
            "C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H]"
            "(Cc2ccccc2)C(=O)O)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H]"
            "(CCC(=O)O)NC(=O)[C@@H](CO)NC(=O)[C@@H]2CCCN2C(=O)[C@H]"
            "(CCC(=O)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](C)NC(=O)[C@H]"
            "(CC(=O)O)NC(=O)[C@H](CS)NC(=O)[C@H](C)NC1=O"
    )

    peptide__dict__fx = {
        "smiles": smiles_fx,
        "complete_smiles": complete_smiles_fx,
        "peptide_json": {
            "length": 17,
            "sequence": "CACDAPEPsEQCGCDEF",
            "internal_modifications": [],
            "C_terminus": "OH",
            "N_terminus": "H",
            "pepseq_format": "H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH",
            "symbols": [
                "H",
                "Cys(R1)",
                "A",
                "C",
                "D",
                "A",
                "P",
                "E",
                "P",
                "s",
                "E",
                "Q",
                "Cys(R2)",
                "G",
                "Cys(R3)",
                "D",
                "E",
                "F",
                "OH",
            ],
            "external_modifications": [
                {
                    "smiles": "[*:1]CNCC[*:2]",
                    "max_attachment_point_id": 2,
                    "attachment_points_on_sequence": {
                        "2": {
                            "attachment_point_id": "2",
                            "ResID": "12",
                            "AtomName": "SG",
                            "ResidueName": "Cys",
                        },
                        "1": {
                            "attachment_point_id": "1",
                            "ResID": "1",
                            "AtomName": "SG",
                            "ResidueName": "Cys",
                        },
                    },
                },
                {
                    "smiles": "[*:3]CNCC",
                    "max_attachment_point_id": 3,
                    "attachment_points_on_sequence": {
                        "3": {
                            "attachment_point_id": "3",
                            "ResID": "14",
                            "AtomName": "SG",
                            "ResidueName": "Cys",
                        }
                    },
                },
            ],
        },
        "sequence": "CACDAPEPsEQCGCDEF",
        "length": 17,
        "mw": 1916.13,
    }

    assert peptide.__dict__ == peptide__dict__fx
