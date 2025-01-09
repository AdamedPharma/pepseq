import rdkit
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import (
    mol_to_nx,
    nx_to_mol,
    nx_to_json,
    mol_json_to_nx,
    mol_json_to_mol,
)
from pepseq.tests_new.helpers import smiles_are_identical


def test_json_to_nx_to_mol():
    """
    Test the conversion from a JSON representation of a molecule to a NetworkX graph
    and then back to an RDKit molecule, ensuring the SMILES representation matches
    the expected output.
    The test performs the following steps:
    1. Converts a JSON representation of a molecule to a NetworkX graph using `mol_json_to_nx`.
    2. Converts the NetworkX graph back to an RDKit molecule using `nx_to_mol`.
    3. Converts the RDKit molecule to a SMILES string using `rdkit.Chem.MolToSmiles`.
    4. Asserts that the generated SMILES string is identical to the expected SMILES string "[*:1]CNCC[*:2]".
    The JSON representation includes:
    - `nodes_tuple`: A list of nodes with attributes such as atomic number, chiral tag, etc.
    - `nodes_columns`: A list of column names corresponding to the attributes of the nodes.
    - `edges_tuple`: A list of edges with attributes such as bond type, bond start, and bond end.
    - `edges_columns`: A list of column names corresponding to the attributes of the edges.
    """
    j2 = {
        "nodes_tuple": [
            [0, 0, "*", 0, None, False, 0, "1", 0, 0],
            [6, 0, None, 0, 4, False, 0, None, 0, 1],
            [7, 0, None, 0, 4, False, 0, None, 0, 2],
            [6, 0, None, 0, 4, False, 0, None, 0, 3],
            [6, 0, None, 0, 4, False, 0, None, 0, 4],
            [0, 0, "*", 0, None, False, 0, "2", 0, 5],
        ],
        "nodes_columns": [
            "atomic_num",
            "chiral_tag",
            "dummyLabel",
            "formal_charge",
            "hybridization",
            "is_aromatic",
            "isotope",
            "molAtomMapNumber",
            "num_explicit_hs",
            "node_id",
        ],
        "edges_tuple": (
            (1, None, 0, 1),
            (1, None, 1, 2),
            (1, None, 2, 3),
            (1, None, 3, 4),
            (1, None, 4, 5),
        ),
        "edges_columns": ["bond_type", "is_peptide_bond", "bond_start", "bond_end"],
    }

    nxg = mol_json_to_nx(j2)
    nx_mol = nx_to_mol(nxg)
    nx_smi = rdkit.Chem.MolToSmiles(nx_mol)
    assert smiles_are_identical(nx_smi, "[*:1]CNCC[*:2]")


def test_mol_to_nx_to_json():
    """
    Test the conversion of a molecule from RDKit Mol object to NetworkX graph and then to JSON format,
    and verify the correctness of the conversion.
    The test performs the following steps:
    1. Create a molecule from a SMILES string using RDKit.
    2. Convert the molecule to a NetworkX graph.
    3. Convert the NetworkX graph to a JSON representation.
    4. Assert that the JSON representation matches the expected fixture.
    5. Convert the JSON representation back to a NetworkX graph.
    6. Convert the NetworkX graph back to an RDKit Mol object.
    7. Convert the RDKit Mol object to a SMILES string.
    8. Assert that the SMILES string matches the expected SMILES string.
    The test uses two fixtures:
    - mol_json_fixture: The expected JSON representation of the molecule.
    - mol_json_fixture_2: The expected JSON representation of the molecule with a different node_id order.
    The test ensures that the conversion functions `mol_to_nx`, `nx_to_json`, and `mol_json_to_mol`
    work correctly and that the resulting molecule is identical to the original.
    """
    mol_json_fixture = {
        "nodes_tuple": [
            [0, 0, "*", 0, None, False, 0, "1", 0, 0],
            [6, 0, None, 0, 4, False, 0, None, 1, 0],
            [7, 0, None, 0, 4, False, 0, None, 2, 0],
            [6, 0, None, 0, 4, False, 0, None, 3, 0],
            [6, 0, None, 0, 4, False, 0, None, 4, 0],
            [0, 0, "*", 0, None, False, 0, "2", 5, 0],
        ],
        "nodes_columns": [
            "atomic_num",
            "chiral_tag",
            "dummyLabel",
            "formal_charge",
            "hybridization",
            "is_aromatic",
            "isotope",
            "molAtomMapNumber",
            "node_id",
            "num_explicit_hs",
        ],
        "edges_tuple": (
            (1, None, 0, 1),
            (1, None, 1, 2),
            (1, None, 2, 3),
            (1, None, 3, 4),
            (1, None, 4, 5),
        ),
        "edges_columns": ["bond_type", "is_peptide_bond", "bond_start", "bond_end"],
    }

    mol_json_fixture_2 = {
        "nodes_tuple": [
            [0, 0, "*", 0, None, False, 0, 1, 0, 0],
            [6, 0, None, 0, 4, False, 0, None, 0, 1],
            [7, 0, None, 0, 4, False, 0, None, 0, 2],
            [6, 0, None, 0, 4, False, 0, None, 0, 3],
            [6, 0, None, 0, 4, False, 0, None, 0, 4],
            [0, 0, "*", 0, None, False, 0, 2, 0, 5],
        ],
        "nodes_columns": [
            "atomic_num",
            "chiral_tag",
            "dummyLabel",
            "formal_charge",
            "hybridization",
            "is_aromatic",
            "isotope",
            "molAtomMapNumber",
            "num_explicit_hs",
            "node_id",
        ],
        "edges_tuple": (
            (1, None, 0, 1),
            (1, None, 1, 2),
            (1, None, 2, 3),
            (1, None, 3, 4),
            (1, None, 4, 5),
        ),
        "edges_columns": ["bond_type", "is_peptide_bond", "bond_start", "bond_end"],
    }

    m = rdkit.Chem.MolFromSmiles("[*:1]CNCC[*:2]")
    G = mol_to_nx(m)
    mol_json = nx_to_json(G)
    assert mol_json == mol_json_fixture_2

    mol_json_to_nx(mol_json_fixture)
    mol = mol_json_to_mol(mol_json_fixture)
    smi = rdkit.Chem.MolToSmiles(mol)
    assert smiles_are_identical(smi, "C(C[*:2])NC[*:1]")
