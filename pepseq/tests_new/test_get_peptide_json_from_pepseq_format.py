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
