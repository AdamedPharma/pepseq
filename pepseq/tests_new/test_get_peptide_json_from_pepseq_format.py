import json
import pkgutil
from pepseq.get_peptide_json_from_pepseq_format import (
    get_single_modification_json,
    get_ext_mod_json,
    get_smiles_json,
    get_pep_json,
)

smiles2 = ["[*:1]CNCC[*:2]", "[*:3]CNCC"]
symbols_fx = [
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
    "A",
    "K",
    "Cys(R3)",
]


db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


def test_get_single_modification_json():
    """
    Test the `get_single_modification_json` function to ensure it correctly generates
    the JSON representation of a single modification based on the provided attachment
    points and SMILES strings.
    The test verifies the following:
    - The function correctly handles a modification with a single attachment point.
    - The function correctly handles a modification with multiple attachment points.
    Test data:
    - `attachment_points_on_sequence_fx`: A dictionary representing attachment points on a sequence.
    - `single_modification_json_fx2`: Expected JSON output for a modification with a single attachment point.
    - `single_modification_json_fx1`: Expected JSON output for a modification with multiple attachment points.
    Assertions:
    - The function's output for a single attachment point matches `single_modification_json_fx2`.
    - The function's output for multiple attachment points matches `single_modification_json_fx1`.
    """
    attachment_points_on_sequence_fx = {
        1: {
            "attachment_point_id": "1",
            "ResID": "1",
            "AtomName": "SG",
            "ResidueName": "Cys",
        },
        2: {
            "attachment_point_id": "2",
            "ResID": "12",
            "AtomName": "SG",
            "ResidueName": "Cys",
        },
        3: {
            "attachment_point_id": "3",
            "ResID": "15",
            "AtomName": "SG",
            "ResidueName": "Cys",
        },
    }
    single_modification_json_fx2 = {
        "smiles": "[*:3]CNCC",
        "max_attachment_point_id": 3,
        "attachment_points_on_sequence": {
            "3": {
                "attachment_point_id": "3",
                "ResID": "15",
                "AtomName": "SG",
                "ResidueName": "Cys",
            }
        },
    }

    single_modification_json_fx1 = {
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
    }
    res_fx2 = get_single_modification_json(attachment_points_on_sequence_fx, smiles2[1])

    res_fx1 = get_single_modification_json(attachment_points_on_sequence_fx, smiles2[0])

    assert res_fx2 == single_modification_json_fx2
    assert res_fx1 == single_modification_json_fx1
    return


def test_get_ext_mod_json():
    """
    Test the `get_ext_mod_json` function to ensure it returns the correct JSON structure
    for extended modifications.
    The function is expected to return a list of dictionaries, each representing a
     modification with its SMILES string, maximum attachment point ID, and attachment
     points on the sequence.
    The expected result `r_fx` contains two modifications:
    1. A modification with SMILES "[*:1]CNCC[*:2]", max attachment point ID 2, and
        attachment points on sequence with IDs 1 and 2.
    2. A modification with SMILES "[*:3]CNCC", max attachment point ID 3, and
        an attachment point on sequence with ID 3.
    The test asserts that the result from `get_ext_mod_json(symbols_fx, smiles2)`
     matches the expected result `r_fx`.
    """
    r = get_ext_mod_json(symbols_fx, smiles2)
    r_fx = [
        {
            "smiles": "[*:1]CNCC[*:2]",
            "max_attachment_point_id": 2,
            "attachment_points_on_sequence": {
                "1": {
                    "attachment_point_id": "1",
                    "ResID": "1",
                    "AtomName": "SG",
                    "ResidueName": "Cys",
                },
                "2": {
                    "attachment_point_id": "2",
                    "ResID": "12",
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
                    "ResID": "15",
                    "AtomName": "SG",
                    "ResidueName": "Cys",
                }
            },
        },
    ]
    assert r == r_fx
    return


def test_get_smiles_json():
    """
    Test the `get_smiles_json` function to ensure it correctly processes and returns
    the expected JSON structure for given SMILES strings and their attachment points.
    The test verifies that the function:
    - Correctly identifies and processes the attachment points on the sequence.
    - Returns the expected JSON structure with the correct SMILES strings and attachment points.
    Expected output:
    A list of dictionaries, each containing:
    - `smiles`: The SMILES string with attachment points.
    - `max_attachment_point_id`: The maximum attachment point ID.
    - `attachment_points_on_sequence`: A dictionary of attachment points with details such as
      attachment point ID, residue ID, atom name, and residue name.
    Asserts:
    - The output of `get_smiles_json` matches the expected JSON structure (`ext_mod_fx`).
    """
    ext_mod_fx = [
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
                    "ResID": "15",
                    "AtomName": "SG",
                    "ResidueName": "Cys",
                }
            },
        },
    ]

    ext_mod = get_smiles_json(symbols_fx, smiles2)
    assert ext_mod == ext_mod_fx
    return


def test_get_pep_json():
    """
    Test the `get_pep_json` function to ensure it correctly converts a peptide sequence
    in pepseq format to a JSON representation.
    The test uses a sample pepseq format string and compares the output of the
     `get_pep_json` function to an expected JSON structure.
    The expected JSON structure includes:
    - Length of the peptide sequence
    - The peptide sequence itself
    - Internal modifications (if any)
    - C-terminus and N-terminus information
    - The original pepseq format string
    - Symbols representing each part of the sequence
    - External modifications with details about attachment points
    The test asserts that the output of `get_pep_json` matches the expected JSON structure.
    """
    pepseq_format = "H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH"

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
    assert peptide_json_fx_2 == get_pep_json(
        pepseq_format, db_json, ["[*:1]CNCC[*:2]", "[*:3]CNCC"]
    )
