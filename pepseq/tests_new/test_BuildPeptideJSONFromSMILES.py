import json
import pkgutil

from pepseq.BuildPeptideJSONFromSMILES import (
    from_smiles_to_pepseq_and_mod_smiles_strings,
    decompose_peptide_smiles_with_termini,
    decompose_peptide_smiles,
    from_smiles_to_pepseq_and_one_mod_smiles_strings,
)


db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


def test_from_smiles_to_pepseq_and_mod_smiles_strings():
    """
    Test the conversion of a SMILES string to a peptide sequence (pepseq) and
     modified SMILES strings (mod_smiles) using the functions
     from_smiles_to_pepseq_and_mod_smiles_strings and
     from_smiles_to_pepseq_and_one_mod_smiles_strings.
    The test verifies that the conversion functions correctly parse the input
     SMILES string and produce the expected peptide sequence and modified SMILES
     strings.
    Assertions:
        - The peptide sequence (pepseq) should match the expected sequence.
        - The list of modified SMILES strings (mod_smiles) should match the
          expected list of strings.
    """
    smiles = "".join(
        [
            "[H]N[C@H]1CSCNCCSC[C@@H](C(=O)NCC(=O)N[C@@H](CSCNCCSP)C(=O)N[C@@H]",
            "(CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc2ccccc2)C(=O)O)NC(=O)",
            "[C@H](CCC(N)=O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](CO)NC(=O)[C@@H]2C",
            "CCN2C(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](C)NC(=O)[C@H]",
            "(CC(=O)O)NC(=O)[C@H](CS)NC(=O)[C@H](C)NC1=O",
        ]
    )

    n_subst_limit = None

    pepseq, mod_smiles = from_smiles_to_pepseq_and_mod_smiles_strings(
        smiles, db_json, n_subst_limit=n_subst_limit
    )

    assert pepseq == "H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH"
    assert mod_smiles == ["C(C[*:2])NC[*:1]", "PSCCNC[*:3]"]

    pepseq, mod_smiles = from_smiles_to_pepseq_and_one_mod_smiles_strings(
        smiles, db_json, n_subst_limit=None
    )

    assert pepseq == "H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH"
    assert mod_smiles == ["C(C[*:2])NC[*:1]", "PSCCNC[*:3]"]


def test_decompose_peptide_smiles():
    """
    Test the decomposition of a peptide SMILES string into a JSON representation.
    This test verifies that the `decompose_peptide_smiles` and
     `decompose_peptide_smiles_with_termini` functions correctly convert a given
     peptide SMILES string into the expected JSON format, including internal and
     external modifications, as well as termini information.
    Fixtures:
        - smiles_fixture: A SMILES string representing the peptide.
        - peptide_json_fixture: The expected JSON representation of the peptide
          without termini information.
        - peptide_json_w_termini_fixture: The expected JSON representation of the
          peptide with termini information.
    Assertions:
        - The JSON output from `decompose_peptide_smiles` matches `peptide_json_fixture`.
        - The JSON output from `decompose_peptide_smiles_with_termini` matches
          `peptide_json_w_termini_fixture`.
    """
    smiles_fixture = "".join(
        [
            "[H]N[C@H]1CSCNCCSC[C@@H](C(=O)NCC(=O)N[C@@H](CSCNCCSP)C(=O)N[C@@H]",
            "(CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc2ccccc2)C(=O)O)NC(=O)",
            "[C@H](CCC(N)=O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](CO)NC(=O)[C@@H]2CC",
            "CN2C(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](C)NC(=O)[C@H]",
            "(CC(=O)O)NC(=O)[C@H](CS)NC(=O)[C@H](C)NC1=O",
        ]
    )

    peptide_json_fixture = {
        "sequence": "CACDAPEPsEQCGCDEF",
        "internal_modifications": [],
        "external_modifications": [
            {
                "smiles": "C(C[*:2])NC[*:1]",
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
                        "ResID": "12",
                        "AtomName": "SG",
                        "ResidueName": "",
                    },
                },
            },
            {
                "smiles": "PSCCNC[*:1]",
                "max_attachment_point_id": 1,
                "attachment_points_on_sequence": {
                    1: {
                        "attachment_point_id": 1,
                        "ResID": "14",
                        "AtomName": "SG",
                        "ResidueName": "",
                    }
                },
            },
            {
                "smiles": "O[*:1]",
                "max_attachment_point_id": 1,
                "attachment_points_on_sequence": {
                    1: {
                        "attachment_point_id": 1,
                        "ResID": "17",
                        "AtomName": "CO",
                        "ResidueName": "",
                    }
                },
            },
        ],
    }

    peptide_json_w_termini_fixture = {
        "sequence": "CACDAPEPsEQCGCDEF",
        "internal_modifications": [],
        "external_modifications": [
            {
                "smiles": "C(C[*:2])NC[*:1]",
                "max_attachment_point_id": 2,
                "attachment_points_on_sequence": {
                    2: {
                        "attachment_point_id": 2,
                        "ResID": "12",
                        "AtomName": "SG",
                        "ResidueName": "",
                    },
                    1: {
                        "attachment_point_id": 1,
                        "ResID": "1",
                        "AtomName": "SG",
                        "ResidueName": "",
                    },
                },
            },
            {
                "smiles": "PSCCNC[*:3]",
                "max_attachment_point_id": 3,
                "attachment_points_on_sequence": {
                    3: {
                        "attachment_point_id": 1,
                        "ResID": "14",
                        "AtomName": "SG",
                        "ResidueName": "",
                    }
                },
            },
        ],
        "C_terminus": "OH",
        "N_terminus": "H",
        "pepseq_format": "H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH",
    }
    n_subst_limit = None
    peptide_json = decompose_peptide_smiles(
        smiles_fixture, db_json, n_subst_limit=n_subst_limit
    )

    assert peptide_json == peptide_json_fixture
    peptide_json_w_termini = decompose_peptide_smiles_with_termini(
        smiles_fixture, db_json, n_subst_limit=n_subst_limit
    )
    assert peptide_json_w_termini == peptide_json_w_termini_fixture

    return
