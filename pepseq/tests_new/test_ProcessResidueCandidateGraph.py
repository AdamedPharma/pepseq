import json
import pkgutil
import rdkit

import importlib

from pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph import (
    add_connection_point_to_molecular_graph,
    process_external_modification,
    process_external_connections,
    full_decomposition,
    decompose_residues_internal,
)
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import (
    nx_to_mol,
    mol_json_to_nx,
    mol_json_to_mol,
)

from pepseq.tests_new.helpers import smiles_are_identical

from pepseq.BuildPeptideJSONFromSMILES import get_cx_smarts_db


db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


def load_tests(name):
    """
    Dynamically loads and yields test cases from a specified module.
    Args:
        name (str): The name of the module to import which contains the test cases.
    Yields:
        object: Each test case found in the `tests` variable of the imported module.
    """
    # Load module which contains test data
    tests_module = importlib.import_module(name)
    # Tests are to be found in the variable `tests` of the module
    for test in tests_module.tests:
        yield test


# content of test_example.py
def pytest_generate_tests(metafunc):
    """
    Automatically generates test cases for the given test function.
    Args:
        metafunc (Metafunc): The test function to generate test cases for.
    """
    for fixture in metafunc.fixturenames:
        print(fixture)
        if fixture.startswith("data_"):
            # Load associated test data
            tests = load_tests(fixture)
            metafunc.parametrize(fixture, tests)


def test_add_connection_point_to_molecular_graph():
    """
    Test the function `add_connection_point_to_molecular_graph` by adding
     a connection point to a molecular graph.
    This test uses a fixture for the point ID and the modification atom,
     and a predefined molecular graph in JSON format.
    It converts the modified graph to a SMILES string and compares it to an
     expected SMILES string to verify correctness.
    The test performs the following steps:
    1. Define the point ID and modification atom.
    2. Define the molecular graph in JSON format.
    3. Convert the JSON graph to a NetworkX graph.
    4. Add a connection point to the molecular graph.
    5. Convert the modified graph to a SMILES string.
    6. Compare the resulting SMILES string to the expected SMILES string.
    Asserts:
        - The SMILES string generated from the modified molecular graph matches the expected SMILES string.
    """
    point_id_fixture = 2
    mod_atom = 9

    mod_graph_copy_json = {
        "nodes_tuple": [
            [6, 0, 0, 4, False, 0, 8, 0, 8],
            [6, 0, 0, 4, False, 0, 9, 0, 9],
            [6, 0, 0, 4, False, 0, 6, 0, 6],
            [7, 0, 0, 4, False, 0, 7, 0, 7],
        ],
        "nodes_columns": [
            "atomic_num",
            "chiral_tag",
            "formal_charge",
            "hybridization",
            "is_aromatic",
            "isotope",
            "node_id",
            "num_explicit_hs",
        ],
        "edges_tuple": ((1, None, 8, 7), (1, None, 8, 9), (1, None, 6, 7)),
        "edges_columns": ["bond_type", "is_peptide_bond", "bond_start", "bond_end"],
    }

    smi = rdkit.Chem.MolToSmiles(
        nx_to_mol(
            add_connection_point_to_molecular_graph(
                mol_json_to_nx(mod_graph_copy_json), point_id_fixture, mod_atom
            )
        )
    )

    smi_fixture = "CNCC[*:2]"
    assert smiles_are_identical(smi, smi_fixture)


def test_process_external_modification():
    """
    Test the process_external_modification function.
    This test verifies that the process_external_modification function correctly processes
    external modifications on a molecular graph and returns the expected external modification
    details and attachment point ID.
    Fixtures:
    - G_mod_json_fixture: A dictionary representing the molecular graph in JSON format.
    - mod_bonds_fixture: A list of tuples representing the modification bonds.
    - res_matches_fixture: A dictionary mapping residue IDs to their corresponding details.
    - atom_names_dict_fixture: A dictionary mapping atom IDs to their names.
    - point_id: An integer representing the point ID.
    Expected Results:
    - external_modification: A dictionary containing the SMILES representation of the external
      modification, the maximum attachment point ID, and the attachment points on the sequence.
    - attachment_point_id: An integer representing the attachment point ID.
    Assertions:
    - The external_modification returned by the function should match the expected fixture.
    - The attachment_point_id returned by the function should match the expected fixture.
    """
    G_mod_json_fixture = {
        "nodes_tuple": [
            [6, 0, 0, 4, False, 0, 8, 0, 8],
            [6, 0, 0, 4, False, 0, 9, 0, 9],
            [6, 0, 0, 4, False, 0, 6, 0, 6],
            [7, 0, 0, 4, False, 0, 7, 0, 7],
        ],
        "nodes_columns": [
            "atomic_num",
            "chiral_tag",
            "formal_charge",
            "hybridization",
            "is_aromatic",
            "isotope",
            "node_id",
            "num_explicit_hs",
        ],
        "edges_tuple": ((1, None, 8, 7), (1, None, 8, 9), (1, None, 6, 7)),
        "edges_columns": ["bond_type", "is_peptide_bond", "bond_start", "bond_end"],
    }
    mod_bonds_fixture = [("Res_1_C", [(3, 6)]), ("Res_12_C", [(9, 10)])]
    res_matches_fixture = {
        "1": ("C", {3: "SG"}, (1, 2, 3, 4, 5, 0)),
        "12": ("C", {10: "SG"}, (12, 11, 10, 13, 14, 15)),
    }
    G_mod_fixture = mol_json_to_nx(G_mod_json_fixture)
    atom_names_dict_fixture = {
        0: "N",
        1: "CA",
        3: "SG",
        4: "CO",
        5: "O",
        10: "SG",
        12: "CA",
        13: "CO",
        14: "O",
        15: "N",
    }
    point_id = 2
    external_modification, attachment_point_id = process_external_modification(
        G_mod_fixture,
        mod_bonds_fixture,
        res_matches_fixture,
        atom_names_dict_fixture,
        point_id,
    )
    attachment_point_id_fixture = 4

    external_modification_fixture = {
        "smiles": "C(C[*:4])NC[*:3]",
        "max_attachment_point_id": 4,
        "attachment_points_on_sequence": {
            3: {
                "attachment_point_id": 3,
                "ResID": "1",
                "AtomName": "SG",
                "ResidueName": "",
            },
            4: {
                "attachment_point_id": 4,
                "ResID": "12",
                "AtomName": "SG",
                "ResidueName": "",
            },
        },
    }
    assert external_modification == external_modification_fixture
    assert attachment_point_id == attachment_point_id_fixture


def test_process_external_connections():
    """
    Test the `process_external_connections` function.
    This test verifies that the `process_external_connections` function correctly processes
    external connections given a set of external modification connections, residue matches,
    modification graphs, and a molecular graph.
    The test uses the following fixtures:
    - `external_connections_fx`: A list of dictionaries representing external connections.
    - `external_mod_connections_fx`: A dictionary representing external modification connections.
    - `res_matches_fx`: A dictionary representing residue matches.
    - `modification_graphs_jsons_fx`: A list of dictionaries representing modification graphs in JSON format.
    - `modification_graphs_fx`: A list of NetworkX graphs converted from `modification_graphs_jsons_fx`.
    - `G_json_fx`: A dictionary representing a molecular graph in JSON format.
    - `G_fx`: A NetworkX graph converted from `G_json_fx`.
    The test asserts that the output of `process_external_connections` matches the expected
    `external_connections_fx`.
    """
    external_connections_fx = [
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
        }
    ]

    external_mod_connections_fx = {
        "1": [("Res_1_C", [(3, 6)]), ("Res_12_C", [(9, 10)])]
    }
    res_matches_fx = {
        "1": ("C", {3: "SG"}, (1, 2, 3, 4, 5, 0)),
        "12": ("C", {10: "SG"}, (12, 11, 10, 13, 14, 15)),
    }
    modification_graphs_jsons_fx = [
        {
            "nodes_tuple": [
                [6, 0, 0, 4, False, 0, 8, 0, 8],
                [6, 0, 0, 4, False, 0, 9, 0, 9],
                [6, 0, 0, 4, False, 0, 6, 0, 6],
                [7, 0, 0, 4, False, 0, 7, 0, 7],
            ],
            "nodes_columns": [
                "atomic_num",
                "chiral_tag",
                "formal_charge",
                "hybridization",
                "is_aromatic",
                "isotope",
                "node_id",
                "num_explicit_hs",
            ],
            "edges_tuple": ((1, None, 8, 7), (1, None, 8, 9), (1, None, 6, 7)),
            "edges_columns": ["bond_type", "is_peptide_bond", "bond_start", "bond_end"],
        }
    ]
    modification_graphs_fx = [mol_json_to_nx(i) for i in modification_graphs_jsons_fx]

    G_json_fx = {
        "nodes_tuple": [
            ["N", "1", "C", 7, 0, 0, 4, False, 0, 0, 1, 0],
            ["CA", "1", "C", 6, 1, 0, 4, False, 0, 1, 1, 1],
            [None, "1", "C", 6, 0, 0, 4, False, 0, 2, 0, 2],
            ["SG", "1", "C", 16, 0, 0, 4, False, 0, 3, 0, 3],
            ["CO", "1", "C", 6, 0, 0, 3, False, 0, 4, 0, 4],
            ["O", "1", "C", 8, 0, 0, 3, False, 0, 5, 0, 5],
            [None, None, None, 6, 0, 0, 4, False, 0, 6, 0, 6],
            [None, None, None, 7, 0, 0, 4, False, 0, 7, 0, 7],
            [None, None, None, 6, 0, 0, 4, False, 0, 8, 0, 8],
            [None, None, None, 6, 0, 0, 4, False, 0, 9, 0, 9],
            ["SG", "12", "C", 16, 0, 0, 4, False, 0, 10, 0, 10],
            [None, "12", "C", 6, 0, 0, 4, False, 0, 11, 0, 11],
            ["CA", "12", "C", 6, 1, 0, 4, False, 0, 12, 1, 12],
            ["CO", "12", "C", 6, 0, 0, 3, False, 0, 13, 0, 13],
            ["O", "12", "C", 8, 0, 0, 3, False, 0, 14, 0, 14],
            ["N", "12", "C", 7, 0, 0, 4, False, 0, 15, 0, 15],
        ],
        "nodes_columns": [
            "AtomName",
            "ResID",
            "ResName",
            "atomic_num",
            "chiral_tag",
            "formal_charge",
            "hybridization",
            "is_aromatic",
            "isotope",
            "node_id",
            "num_explicit_hs",
        ],
        "edges_tuple": (
            (1, None, 0, 1),
            (1, None, 1, 2),
            (1, None, 1, 4),
            (1, None, 2, 3),
            (1, None, 3, 6),
            (2, None, 4, 5),
            (1, None, 6, 7),
            (1, None, 7, 8),
            (1, None, 8, 9),
            (1, None, 9, 10),
            (1, None, 10, 11),
            (1, None, 11, 12),
            (1, None, 12, 13),
            (1, None, 12, 15),
            (2, None, 13, 14),
        ),
        "edges_columns": ["bond_type", "is_peptide_bond", "bond_start", "bond_end"],
    }

    G_fx = mol_json_to_nx(G_json_fx)

    assert external_connections_fx == process_external_connections(
        external_mod_connections_fx, res_matches_fx, modification_graphs_fx, G_fx
    )
    return


def test_full_decomposition():
    """
    Test the full decomposition of a molecule into its fragment residue names and modifications.
    This test verifies that the `full_decomposition` function correctly decomposes a given molecule
    into its constituent fragments and modifications. The test uses predefined fixtures for the
     expected fragment residue names and modifications, and compares the output of the function
     against these fixtures.
    Fixtures:
        - fragment_res_names_fixture: A dictionary mapping residue IDs to residue names.
        - fragment_modifications_fixture: A dictionary containing internal and external modifications
          of the molecule, including attachment points and SMILES strings.
    The test performs the following steps:
        1. Converts a JSON representation of a molecule to a molecule object.
        2. Retrieves a SMARTS database for chemical transformations.
        3. Calls the `full_decomposition` function with the molecule object and SMARTS database.
        4. Asserts that the fragment residue names match the expected fixture.
        5. Extracts the SMILES string from the external modifications and verifies it contains "*:".
        6. Asserts that the SMILES string matches the expected fixture.
        7. Removes the SMILES strings from the modifications and asserts that the remaining
           modifications match the expected fixture.
    """
    fragment_res_names_fixture = {"1": "C", "12": "C"}
    fragment_modifications_fixture = {
        "internal_modifications": [],
        "external_modifications": [
            {
                "smiles": "[*:1]CNCC[*:2]",
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
            }
        ],
    }
    mol0_json_fixture = {
        "nodes_tuple": [
            ["N", "1", 7, 0, 0, 4, False, 0, 0, 1, 0],
            ["CA", "1", 6, 1, 0, 4, False, 0, 1, 1, 1],
            [None, None, 6, 0, 0, 4, False, 0, 2, 0, 2],
            [None, None, 16, 0, 0, 4, False, 0, 3, 0, 3],
            ["CO", "1", 6, 0, 0, 3, False, 0, 130, 0, 130],
            ["O", "1", 8, 0, 0, 3, False, 0, 131, 0, 131],
            [None, None, 6, 0, 0, 4, False, 0, 4, 0, 4],
            [None, None, 7, 0, 0, 4, False, 0, 5, 0, 5],
            [None, None, 6, 0, 0, 4, False, 0, 6, 0, 6],
            [None, None, 6, 0, 0, 4, False, 0, 7, 0, 7],
            [None, None, 16, 0, 0, 4, False, 0, 8, 0, 8],
            [None, None, 6, 0, 0, 4, False, 0, 9, 0, 9],
            ["CA", "12", 6, 1, 0, 4, False, 0, 10, 1, 10],
            ["CO", "12", 6, 0, 0, 3, False, 0, 11, 0, 11],
            ["O", "12", 8, 0, 0, 3, False, 0, 12, 0, 12],
            ["N", "12", 7, 0, 0, 3, False, 0, 58, 0, 58],
        ],
        "nodes_columns": [
            "AtomName",
            "ResID",
            "atomic_num",
            "chiral_tag",
            "formal_charge",
            "hybridization",
            "is_aromatic",
            "isotope",
            "node_id",
            "num_explicit_hs",
        ],
        "edges_tuple": (
            (1, None, 0, 1),
            (1, None, 1, 2),
            (1, None, 1, 130),
            (1, None, 2, 3),
            (1, None, 3, 4),
            (2, None, 130, 131),
            (1, None, 4, 5),
            (1, None, 5, 6),
            (1, None, 6, 7),
            (1, None, 7, 8),
            (1, None, 8, 9),
            (1, None, 9, 10),
            (1, None, 10, 11),
            (1, None, 10, 58),
            (2, None, 11, 12),
        ),
        "edges_columns": ["bond_type", "is_peptide_bond", "bond_start", "bond_end"],
    }

    mol0_fx = mol_json_to_mol(mol0_json_fixture)
    cx_smarts_db = get_cx_smarts_db(db_json)
    n_subst_limit_fx = None

    fragment_res_names, fragment_modifications = full_decomposition(
        mol0_fx, cx_smarts_db, n_subst_limit=n_subst_limit_fx
    )
    assert fragment_res_names == fragment_res_names_fixture

    smi = fragment_modifications.get("external_modifications")[0].get("smiles")
    assert "*:" in smi
    smi_fixture = fragment_modifications_fixture.get("external_modifications")[0].get(
        "smiles"
    )

    assert smiles_are_identical(smi, smi_fixture)

    fragment_modifications.get("external_modifications")[0].pop("smiles")
    fragment_modifications_fixture.get("external_modifications")[0].pop("smiles")
    assert fragment_modifications == fragment_modifications_fixture
    return


def test_decompose_residues_internal(data_decompose_residues_internal):
    """
    Test the decompose_residues_internal function.
    This test verifies that the decompose_residues_internal function correctly decomposes
    a sequence of residues into its components, including internal and external modifications.
    Args:
        data_decompose_residues_internal (tuple): A tuple containing the following elements:
            - fragments: The fragments to be decomposed.
            - cx_smarts_db: The SMARTS database for chemical context.
            - n_subst_limit: The substitution limit.
    Expected:
        The function should return the expected sequence, internal modifications, and external modifications.
    Asserts:
        - The returned sequence matches the expected sequence.
        - The returned internal modifications match the expected internal modifications.
        - The returned external modifications match the expected external modifications.
    """
    seq_fx = "CACDAPEPsEQCGCDEF"
    internal_modifications_fx = []
    external_modifications_fx = [
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
    ]

    fragments, cx_smarts_db, n_subst_limit = data_decompose_residues_internal

    seq, internal_modifications, external_modifications = decompose_residues_internal(
        fragments, cx_smarts_db, n_subst_limit=n_subst_limit
    )
    assert seq == seq_fx
    assert internal_modifications == internal_modifications_fx
    assert external_modifications == external_modifications_fx
