import json
import pkgutil
import rdkit

from pepseq.tests_new.helpers import mols_are_identical

from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import mol_to_nx, nx_to_json, \
    mol_json_to_mol, mol_json_to_nx, get_mol_json, nx_to_mol
from pepseq.Backbone import (MarkingPeptideBackbone, BreakingIntoResidueCandidateSubgraphs)

from pepseq.BuildPeptideJSONFromSMILES import \
    decompose_peptide_smiles_with_termini, decompose_peptide_smiles, get_cx_smarts_db, \
    decompose_peptide_smiles

from pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph import \
    decompose_residues_internal, get_res_matches, full_decomposition, get_matches, \
    match_molecular_graph_to_res_id, decompose, get_subgraph_tuples, get_connections, \
    split_connections_by_type, process_internal_connections, \
    process_external_connections, sequence_dict_to_string


db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)
cx_smarts_db = get_cx_smarts_db(db_json)


import importlib

def load_tests(name):
    # Load module which contains test data
    tests_module = importlib.import_module(name)
    # Tests are to be found in the variable `tests` of the module
    for test in tests_module.tests:
        yield test


# content of test_example.py
def pytest_generate_tests(metafunc):

    for fixture in metafunc.fixturenames:
        print(fixture)
        if fixture.startswith('data_'):
            # Load associated test data
            tests = load_tests(fixture)
            metafunc.parametrize(fixture, tests)


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
    "pepseq_format": "Ac~{Cys(R2)}S{Cys(R3)}ACGCK~NH2",
}

"""
correct_mol_json = {
    'nodes_tuple': (
        (6, 0, 0, 4, 0, False, 0, None, None, 0),
        (6, 0, 0, 3, 0, False, 0, None, None, 1),
        (8, 0, 0, 3, 0, False, 0, None, None, 2),
        (7, 0, 0, 3, 0, False, 0, 'N', '1', 3),
        (6, 0, 1, 4, 1, False, 0, 'CA', '1', 4),
        (6, 0, 0, 4, 0, False, 0, None, None, 5),
        (16, 0, 0, 4, 0, False, 0, None, None, 6),
        (6, 0, 0, 4, 0, False, 0, None, None, 7),
        (35, 0, 0, 4, 0, False, 0, None, None, 8),
        (6, 0, 0, 4, 0, False, 0, None, None, 9),
        (7, 0, 0, 4, 0, False, 0, None, None, 10),
        (15, 0, 0, 4, 0, False, 0, None, None, 11),
        (11, 0, 0, 1, 0, False, 0, None, None, 12),
        (16, 0, 0, 4, 0, False, 0, None, None, 13),
        (6, 0, 0, 4, 0, False, 0, None, None, 14),
        (6, 0, 1, 4, 1, False, 0, 'CA', '3', 15),
        (6, 0, 0, 3, 0, False, 0, 'CO', '3', 16),
        (8, 0, 0, 3, 0, False, 0, 'O', '3', 17),
        (7, 0, 0, 3, 0, False, 0, 'N', '4', 18),
        (6, 0, 1, 4, 1, False, 0, 'CA', '4', 19),
        (6, 0, 0, 4, 0, False, 0, None, None, 20),
        (6, 0, 0, 3, 0, False, 0, 'CO', '4', 21),
        (8, 0, 0, 3, 0, False, 0, 'O', '4', 22),
        (7, 0, 0, 3, 0, False, 0, 'N', '5', 23),
        (6, 0, 1, 4, 1, False, 0, 'CA', '5', 24),
        (6, 0, 0, 4, 0, False, 0, None, None, 25),
        (16, 0, 0, 4, 0, False, 0, None, None, 26),
        (16, 0, 0, 4, 0, False, 0, None, None, 27),
        (6, 0, 0, 4, 0, False, 0, None, None, 28),
        (6, 0, 1, 4, 1, False, 0, 'CA', '7', 29),
        (6, 0, 0, 3, 0, False, 0, 'CO', '7', 30),
        (8, 0, 0, 3, 0, False, 0, 'O', '7', 31),
        (7, 0, 0, 3, 0, False, 0, 'N', '8', 32),
        (6, 0, 1, 4, 1, False, 0, 'CA', '8', 33),
        (6, 0, 0, 4, 0, False, 0, None, None, 34),
        (6, 0, 0, 4, 0, False, 0, None, None, 35),
        (6, 0, 0, 4, 0, False, 0, None, None, 36),
        (6, 0, 0, 4, 0, False, 0, None, None, 37),
        (7, 0, 0, 4, 0, False, 0, None, None, 38),
        (6, 0, 0, 3, 0, False, 0, 'CO', '8', 39),
        (7, 0, 0, 3, 0, False, 0, None, None, 40),
        (8, 0, 0, 3, 0, False, 0, 'O', '8', 41),
        (7, 0, 0, 3, 0, False, 0, 'N', '7', 42),
        (6, 0, 0, 3, 0, False, 0, 'CO', '6', 43),
        (8, 0, 0, 3, 0, False, 0, 'O', '6', 44),
        (6, 0, 0, 4, 0, False, 0, 'CA', '6', 45),
        (7, 0, 0, 3, 0, False, 0, 'N', '6', 46),
        (6, 0, 0, 3, 0, False, 0, 'CO', '5', 47),
        (8, 0, 0, 3, 0, False, 0, 'O', '5', 48),
        (7, 0, 0, 3, 0, False, 0, 'N', '3', 49),
        (6, 0, 0, 3, 0, False, 0, 'CO', '2', 50),
        (8, 0, 0, 3, 0, False, 0, 'O', '2', 51),
        (6, 0, 2, 4, 1, False, 0, 'CA', '2', 52),
        (6, 0, 0, 4, 0, False, 0, None, None, 53),
        (8, 0, 0, 4, 0, False, 0, None, None, 54),
        (7, 0, 0, 3, 0, False, 0, 'N', '2', 55),
        (6, 0, 0, 3, 0, False, 0, 'CO', '1', 56),
        (8, 0, 0, 3, 0, False, 0, 'O', '1', 57)
    ),
    'nodes_columns': [
        'atomic_num', 'formal_charge', 'chiral_tag',
        'hybridization', 'num_explicit_hs', 'is_aromatic',
        'isotope', 'AtomName', 'ResID', 'node_id'
    ],
    'edges_tuple': (
        (1, False, 0, 1), (2, False, 1, 2), (1, False, 1, 3),
        (1, False, 3, 4), (1, False, 4, 5), (1, False, 4, 56),
        (1, False, 5, 6), (1, False, 6, 7), (1, False, 7, 8),
        (1, False, 7, 9), (1, False, 9, 10), (1, False, 10, 11),
        (1, False, 11, 12), (1, False, 11, 13), (1, False, 13, 14),
        (1, False, 14, 15), (1, False, 15, 16), (1, False, 15, 49),
        (2, False, 16, 17), (1, 'True', 16, 18), (1, False, 18, 19),
        (1, False, 19, 20), (1, False, 19, 21), (2, False, 21, 22),
        (1, 'True', 21, 23), (1, False, 23, 24), (1, False, 24, 25),
        (1, False, 24, 47), (1, False, 25, 26), (1, False, 26, 27),
        (1, False, 27, 28), (1, False, 28, 29), (1, False, 29, 30),
        (1, False, 29, 42), (2, False, 30, 31), (1, 'True', 30, 32),
        (1, False, 32, 33), (1, False, 33, 34), (1, False, 33, 39),
        (1, False, 34, 35), (1, False, 35, 36), (1, False, 36, 37),
        (1, False, 37, 38), (1, False, 39, 40), (2, False, 39, 41),
        (1, 'True', 42, 43), (2, False, 43, 44), (1, False, 43, 45),
        (1, False, 45, 46), (1, 'True', 46, 47), (2, False, 47, 48),
        (1, 'True', 49, 50), (2, False, 50, 51), (1, False, 50, 52),
        (1, False, 52, 53), (1, False, 52, 55), (1, False, 53, 54),
        (1, 'True', 55, 56), (2, False, 56, 57)
    ),
 'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']
}
"""

correct_internal_modifications = [
    {
        1: [
            {'ResID': '5', 'AtomName': 'SG', 'ResidueName': ''},
            {'ResID': '7', 'AtomName': 'SG', 'ResidueName': ''}
        ]
    }
]


correct_external_modifications = [
    {
        'smiles': '[1*]C(C)=O',
        'max_attachment_point_id': 1,
        'attachment_points_on_sequence': {
            1: {
                'attachment_point_id': 1,
                'ResID': '1',
                'AtomName': 'N',
                'ResidueName': ''
            }
        }
    },
    {
        'smiles': '[2*]C(Br)CNP([3*])[Na]',
        'max_attachment_point_id': 3,
        'attachment_points_on_sequence': {
            2: {
                'attachment_point_id': 2,
                'ResID': '1',
                'AtomName': 'SG',
                'ResidueName': ''
            },
            3: {
                'attachment_point_id': 3,
                'ResID': '3',
                'AtomName': 'SG',
                'ResidueName': ''
            }
        }
    },
    {
        'smiles': '[1*]N',
        'max_attachment_point_id': 1,
        'attachment_points_on_sequence': {
            1: {
                'attachment_point_id': 1,
                'ResID': '8',
                'AtomName': 'CO',
                'ResidueName': ''
            }
        }
    }
]




def test_decompose_peptide_smiles_db(data_mol_N_C_smiles):
    mol_N_C_smiles_val, correct_peptide_json, _ = data_mol_N_C_smiles

    peptide_json = decompose_peptide_smiles(mol_N_C_smiles_val, db_json)
    assert peptide_json == correct_peptide_json

mol_N_C_smiles_val = (
    "CC(=O)N[C@H]1CSC(Br)CNP([Na])SC[C@@H](C(=O)N[C@"
    + "@H](C)C(=O)N[C@H]2CSSC[C@@H](C(=O)N[C@@H](CCCCN)C(N)=O)NC(=O)CNC"
    + "2=O)NC(=O)[C@H](CO)NC1=O"
)


def test_decompose_peptide_smiles_db_termini():

    peptide_json_NC = decompose_peptide_smiles_with_termini(mol_N_C_smiles_val, db_json)
    assert peptide_json_NC == correct_peptide_json_NC


def test_BreakingIntoResidueCandidateSubgraphs(data_break_into_residues):
    mol_N_C_smiles_val, residues_vals = data_break_into_residues
    peptide_molecule = rdkit.Chem.MolFromSmiles(mol_N_C_smiles_val)
    peptide_molecule2 = MarkingPeptideBackbone().execute(peptide_molecule)
    residues = BreakingIntoResidueCandidateSubgraphs().execute(peptide_molecule2)


    assert json.loads(json.dumps(nx_to_json(residues[0]))) == residues_vals[0]
    assert json.loads(json.dumps(nx_to_json(residues[2]))) == residues_vals[2]
    assert json.loads(json.dumps(nx_to_json(residues[3]))) == residues_vals[3]
    return


def test_MarkingPeptideBackbone(data_mol_N_C_smiles):

    mol_N_C_smiles_val, correct_peptide_json, peptide_bonds = data_mol_N_C_smiles

    peptide_molecule = rdkit.Chem.MolFromSmiles(mol_N_C_smiles_val)
  
    peptide_molecule2 = MarkingPeptideBackbone().execute(peptide_molecule)
    G = mol_to_nx(peptide_molecule2)
    selected_edges = [(u,v) for u,v,e in G.edges(data=True) if e.get('is_peptide_bond') == 'True']
    assert sorted(selected_edges) == peptide_bonds
    #mark N, CA, CO ?


def test_decompose_residues_internal():
    #data_decompose_residues_internal
    #residues, correct_seq, correct_internal_modifications, correct_external_modifications = data_decompose_residues_internal
    residues_jsons = [
        {
            'nodes_tuple': [
                [6, 0, 0, 4, 0, False, 0, None, None, 0], [6, 0, 0, 3, 0, False, 0, None, None, 1],
                [8, 0, 0, 3, 0, False, 0, None, None, 2],[7, 0, 0, 3, 0, False, 0, 'N', '1', 3],
                [6, 0, 1, 4, 1, False, 0, 'CA', '1', 4],[6, 0, 0, 4, 0, False, 0, None, None, 5],
                [16, 0, 0, 4, 0, False, 0, None, None, 6],[6, 0, 0, 4, 0, False, 0, None, None, 7],
                [35, 0, 0, 4, 0, False, 0, None, None, 8],[6, 0, 0, 4, 0, False, 0, None, None, 9],
                [7, 0, 0, 4, 0, False, 0, None, None, 10],[15, 0, 0, 4, 0, False, 0, None, None, 11],
                [11, 0, 0, 1, 0, False, 0, None, None, 12],[16, 0, 0, 4, 0, False, 0, None, None, 13],
                [6, 0, 0, 4, 0, False, 0, None, None, 14],[6, 0, 1, 4, 1, False, 0, 'CA', '3', 15],
                [6, 0, 0, 3, 0, False, 0, 'CO', '3', 16],[8, 0, 0, 3, 0, False, 0, 'O', '3', 17],
                [7, 0, 0, 3, 0, False, 0, 'N', '3', 49],[6, 0, 0, 3, 0, False, 0, 'CO', '1', 56],
                [8, 0, 0, 3, 0, False, 0, 'O', '1', 57]
                ],
            'nodes_columns': [
                'atomic_num', 'formal_charge', 'chiral_tag', 'hybridization', 'num_explicit_hs',
                'is_aromatic','isotope','AtomName','ResID','node_id'],
            'edges_tuple': (
                (1, None, 0, 1), (2, None, 1, 2), (1, None, 1, 3),(1, None, 3, 4),
                (1, None, 4, 5),(1, None, 4, 56),(1, None, 5, 6),(1, None, 6, 7),
                (1, None, 7, 8),(1, None, 7, 9),(1, None, 9, 10),(1, None, 10, 11),
                (1, None, 11, 12),(1, None, 11, 13),(1, None, 13, 14),(1, None, 14, 15),
                (1, None, 15, 16),(1, None, 15, 49),(2, None, 16, 17),(2, None, 56, 57)),
            'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
        {'nodes_tuple': [[7, 0, 0, 3, 0, False, 0, 'N', '4', 18],
                         [6, 0, 1, 4, 1, False, 0, 'CA', '4', 19],
   [6, 0, 0, 4, 0, False, 0, None, None, 20],
   [6, 0, 0, 3, 0, False, 0, 'CO', '4', 21],
   [8, 0, 0, 3, 0, False, 0, 'O', '4', 22]],
    'nodes_columns': [
        'atomic_num', 'formal_charge', 'chiral_tag','hybridization','num_explicit_hs',
        'is_aromatic','isotope','AtomName','ResID','node_id'],
    'edges_tuple': ((1, None, 18, 19),(1, None, 19, 20),(1, None, 19, 21),
                    (2, None, 21, 22)),
  'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
  {
      'nodes_tuple': [
          [7, 0, 0, 3, 0, False, 0, 'N', '7', 42],[6, 0, 0, 3, 0, False, 0, 'CO', '5', 47],
          [8, 0, 0, 3, 0, False, 0, 'O', '5', 48],[7, 0, 0, 3, 0, False, 0, 'N', '5', 23],
          [6, 0, 1, 4, 1, False, 0, 'CA', '5', 24],[6, 0, 0, 4, 0, False, 0, None, None, 25],
          [16, 0, 0, 4, 0, False, 0, None, None, 26],[16, 0, 0, 4, 0, False, 0, None, None, 27],
          [6, 0, 0, 4, 0, False, 0, None, None, 28],[6, 0, 1, 4, 1, False, 0, 'CA', '7', 29],
          [6, 0, 0, 3, 0, False, 0, 'CO', '7', 30],[8, 0, 0, 3, 0, False, 0, 'O', '7', 31]],
        'nodes_columns': [
            'atomic_num','formal_charge','chiral_tag','hybridization','num_explicit_hs','is_aromatic',
            'isotope','AtomName','ResID','node_id'],
        'edges_tuple': (
            (1, None, 42, 29),(1, None, 47, 24),(2, None, 47, 48),(1, None, 23, 24),
            (1, None, 24, 25),(1, None, 25, 26),(1, None, 26, 27),(1, None, 27, 28),
            (1, None, 28, 29),(1, None, 29, 30),(2, None, 30, 31)),
        'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
    {
        'nodes_tuple': [
            [7, 0, 0, 3, 0, False, 0, 'N', '8', 32],[6, 0, 1, 4, 1, False, 0, 'CA', '8', 33],
            [6, 0, 0, 4, 0, False, 0, None, None, 34],[6, 0, 0, 4, 0, False, 0, None, None, 35],
            [6, 0, 0, 4, 0, False, 0, None, None, 36],[6, 0, 0, 4, 0, False, 0, None, None, 37],
            [7, 0, 0, 4, 0, False, 0, None, None, 38],[6, 0, 0, 3, 0, False, 0, 'CO', '8', 39],
            [7, 0, 0, 3, 0, False, 0, None, None, 40],[8, 0, 0, 3, 0, False, 0, 'O', '8', 41]],
        'nodes_columns': [
            'atomic_num','formal_charge','chiral_tag','hybridization','num_explicit_hs',
            'is_aromatic','isotope','AtomName','ResID','node_id'],
        'edges_tuple': (
            (1, None, 32, 33),(1, None, 33, 34),(1, None, 33, 39),(1, None, 34, 35),
            (1, None, 35, 36),(1, None, 36, 37),(1, None, 37, 38),(1, None, 39, 40),
            (2, None, 39, 41)),
        'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
    {
        'nodes_tuple': [
            [6, 0, 0, 3, 0, False, 0, 'CO', '6', 43],[8, 0, 0, 3, 0, False, 0, 'O', '6', 44],
            [6, 0, 0, 4, 0, False, 0, 'CA', '6', 45],[7, 0, 0, 3, 0, False, 0, 'N', '6', 46]],
        'nodes_columns': [
            'atomic_num','formal_charge','chiral_tag','hybridization','num_explicit_hs',
            'is_aromatic','isotope','AtomName','ResID','node_id'],
        'edges_tuple': ((2, None, 43, 44), (1, None, 43, 45), (1, None, 45, 46)),
        'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
    {
        'nodes_tuple': [
            [6, 0, 0, 3, 0, False, 0, 'CO', '2', 50],[8, 0, 0, 3, 0, False, 0, 'O', '2', 51],
            [6, 0, -1, 4, 1, False, 0, 'CA', '2', 52],[6, 0, 0, 4, 0, False, 0, None, None, 53],
            [8, 0, 0, 4, 0, False, 0, None, None, 54],[7, 0, 0, 3, 0, False, 0, 'N', '2', 55]],
        'nodes_columns': [
            'atomic_num','formal_charge','chiral_tag','hybridization','num_explicit_hs',
            'is_aromatic','isotope','AtomName','ResID','node_id'],
        'edges_tuple': (
            (2, None, 50, 51),(1, None, 50, 52),(1, None, 52, 53),(1, None, 52, 55),
            (1, None, 53, 54)),
        'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']}]
    

    residues = [mol_json_to_nx(i) for i in residues_jsons]

    seq, internal_modifications, external_modifications = decompose_residues_internal(
        residues, cx_smarts_db
    )


    assert seq == 'CSCACGCK'
    assert internal_modifications == correct_internal_modifications
    assert external_modifications == correct_external_modifications



mol_q_dict = {
    'smiles': 'N=C(N)NCCC[C@H](N)C(N)=O |atomProp:0.hybridization.SP2:0.num_explicit_hs.0:1.hybridization.SP2:1.num_explicit_hs.0:2.hybridization.SP2:2.num_explicit_hs.0:3.hybridization.SP2:3.num_explicit_hs.0:4.hybridization.SP3:4.num_explicit_hs.0:5.hybridization.SP3:5.num_explicit_hs.0:6.hybridization.SP3:6.num_explicit_hs.0:7.num_explicit_hs.1:7.ResID.14:7.hybridization.SP3:7.AtomName.CA:8.num_explicit_hs.0:8.ResID.14:8.hybridization.SP2:8.AtomName.N:9.num_explicit_hs.0:9.ResID.14:9.hybridization.SP2:9.AtomName.CO:10.hybridization.SP2:10.num_explicit_hs.0:11.num_explicit_hs.0:11.ResID.14:11.hybridization.SP2:11.AtomName.O|'
}

def test_get_res_matches():
    mol_q = rdkit.Chem.MolFromSmiles(mol_q_dict.get('smiles'))

    res_name, atom_names_dict, nums = get_res_matches(mol_q, cx_smarts_db)['14']

    assert res_name == 'R'
    assert nums == (8, 7, 6, 5, 4, 3, 1, 0, 2, 9, 11)
    #assert rdkit.Chem.MolToCXSmiles(mol) == '*[C@H](C=O)[CH2][CH2][CH2][NH]C(=*)[NH2]'
    return


def test_full_decomposition():
    residues_0_json = {
    'nodes_tuple': [
        [6, 0, 0, 4, 0, False, 0, None, None, 0], [6, 0, 0, 3, 0, False, 0, None, None, 1],
        [8, 0, 0, 3, 0, False, 0, None, None, 2], [7, 0, 0, 3, 0, False, 0, 'N', '1', 3],
        [6, 0, 1, 4, 1, False, 0, 'CA', '1', 4], [6, 0, 0, 4, 0, False, 0, None, None, 5],
        [16, 0, 0, 4, 0, False, 0, None, None, 6], [6, 0, 0, 4, 0, False, 0, None, None, 7],
        [35, 0, 0, 4, 0, False, 0, None, None, 8], [6, 0, 0, 4, 0, False, 0, None, None, 9],
        [7, 0, 0, 4, 0, False, 0, None, None, 10], [15, 0, 0, 4, 0, False, 0, None, None, 11],
        [11, 0, 0, 1, 0, False, 0, None, None, 12], [16, 0, 0, 4, 0, False, 0, None, None, 13],
        [6, 0, 0, 4, 0, False, 0, None, None, 14], [6, 0, 1, 4, 1, False, 0, 'CA', '3', 15],
        [6, 0, 0, 3, 0, False, 0, 'CO', '3', 16], [8, 0, 0, 3, 0, False, 0, 'O', '3', 17],
        [7, 0, 0, 3, 0, False, 0, 'N', '3', 49], [6, 0, 0, 3, 0, False, 0, 'CO', '1', 56],
        [8, 0, 0, 3, 0, False, 0, 'O', '1', 57]],
    'nodes_columns': [
        'atomic_num', 'formal_charge', 'chiral_tag', 'hybridization', 'num_explicit_hs',
        'is_aromatic', 'isotope', 'AtomName', 'ResID', 'node_id'],
    'edges_tuple': (
        (1, None, 0, 1), (2, None, 1, 2), (1, None, 1, 3), (1, None, 3, 4),
        (1, None, 4, 5), (1, None, 4, 56),(1, None, 5, 6), (1, None, 6, 7),
        (1, None, 7, 8), (1, None, 7, 9),(1, None, 9, 10),(1, None, 10, 11),
        (1, None, 11, 12),(1, None, 11, 13),(1, None, 13, 14), (1, None, 14, 15),
        (1, None, 15, 16),(1, None, 15, 49),(2, None, 16, 17),(2, None, 56, 57)
    ),
    'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']}


    decomposed_json = full_decomposition(  mol_json_to_mol(residues_0_json), cx_smarts_db )

    assert decomposed_json == (
    {'1': 'C', '3': 'C'},
    {
        'internal_modifications': [],
        'external_modifications': [
            {
                'smiles': '[1*]C(C)=O',
                'max_attachment_point_id': 1,
                'attachment_points_on_sequence': {
                    1: {
                        'attachment_point_id': 1,
                        'ResID': '1',
                        'AtomName': 'N',
                        'ResidueName': ''
                    }
                }
            },
            {
                'smiles': '[2*]C(Br)CNP([3*])[Na]',
                'max_attachment_point_id': 3,
                'attachment_points_on_sequence': {
                    2: {
                        'attachment_point_id': 2,
                        'ResID': '1',
                        'AtomName': 'SG',
                        'ResidueName': ''
                    },
                    3: {
                        'attachment_point_id': 3,
                        'ResID': '3',
                        'AtomName': 'SG',
                        'ResidueName': ''
                    }
                }
            }
        ]
    }
    )


    return

residues_0_json = {'nodes_tuple': [[6, 0, 0, 4, 0, False, 0, None, None, 0],
  [6, 0, 0, 3, 0, False, 0, None, None, 1],
  [8, 0, 0, 3, 0, False, 0, None, None, 2],
  [7, 0, 0, 3, 0, False, 0, 'N', '1', 3],
  [6, 0, 1, 4, 1, False, 0, 'CA', '1', 4],
  [6, 0, 0, 4, 0, False, 0, None, None, 5],
  [16, 0, 0, 4, 0, False, 0, None, None, 6],
  [6, 0, 0, 4, 0, False, 0, None, None, 7],
  [35, 0, 0, 4, 0, False, 0, None, None, 8],
  [6, 0, 0, 4, 0, False, 0, None, None, 9],
  [7, 0, 0, 4, 0, False, 0, None, None, 10],
  [15, 0, 0, 4, 0, False, 0, None, None, 11],
  [11, 0, 0, 1, 0, False, 0, None, None, 12],
  [16, 0, 0, 4, 0, False, 0, None, None, 13],
  [6, 0, 0, 4, 0, False, 0, None, None, 14],
  [6, 0, 1, 4, 1, False, 0, 'CA', '3', 15],
  [6, 0, 0, 3, 0, False, 0, 'CO', '3', 16],
  [8, 0, 0, 3, 0, False, 0, 'O', '3', 17],
  [7, 0, 0, 3, 0, False, 0, 'N', '3', 49],
  [6, 0, 0, 3, 0, False, 0, 'CO', '1', 56],
  [8, 0, 0, 3, 0, False, 0, 'O', '1', 57]],
 'nodes_columns': ['atomic_num',
  'formal_charge',
  'chiral_tag',
  'hybridization',
  'num_explicit_hs',
  'is_aromatic',
  'isotope',
  'AtomName',
  'ResID',
  'node_id'],
 'edges_tuple': ((1, None, 0, 1),
  (2, None, 1, 2),
  (1, None, 1, 3),
  (1, None, 3, 4),
  (1, None, 4, 5),
  (1, None, 4, 56),
  (1, None, 5, 6),
  (1, None, 6, 7),
  (1, None, 7, 8),
  (1, None, 7, 9),
  (1, None, 9, 10),
  (1, None, 10, 11),
  (1, None, 11, 12),
  (1, None, 11, 13),
  (1, None, 13, 14),
  (1, None, 14, 15),
  (1, None, 15, 16),
  (1, None, 15, 49),
  (2, None, 16, 17),
  (2, None, 56, 57)),
 'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']}

residues_0_mol = mol_json_to_mol(residues_0_json)


def test_get_matches_dict():

    matches_dict = get_matches(residues_0_mol, cx_smarts_db)

    C_mol_json_bkp = {
        'nodes_tuple': [
            [None, 0, 0, 0, None, False, 0, 0, 0, 0], [None, 6, 0, 0, None, False, 0, 1, 0, 1],
            [None, 6, 0, 0, None, False, 0, 2, 0, 2], ['SG', 0, 0, 0, None, False, 0, 3, 0, 3],
            [None, 6, 0, 0, None, False, 0, 4, 0, 4], [None, 8, 0, 0, None, False, 0, 5, 0, 5]],
        'nodes_columns': [
            'AtomName', 'atomic_num', 'chiral_tag', 'formal_charge', 'hybridization', 'is_aromatic',
            'isotope', 'node_id', 'num_explicit_hs'
            ],
        'edges_tuple': ((1, None, 0, 1), (1, None, 1, 2), (1, None, 1, 4), (1, None, 2, 3),
        (2, None, 4, 5)),
        'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']}

    C_mol_json = {'nodes_tuple': [[None, 0, 0, 0, None, False, 0, 0, 0],
  [None, 6, 0, 0, None, False, 0, 0, 1],
  [None, 6, 0, 0, None, False, 0, 0, 2],
  ['SG', 0, 0, 0, None, False, 0, 0, 3],
  [None, 6, 0, 0, None, False, 0, 0, 4],
  [None, 8, 0, 0, None, False, 0, 0, 5]],
 'nodes_columns': ['AtomName',
  'atomic_num',
  'chiral_tag',
  'formal_charge',
  'hybridization',
  'is_aromatic',
  'isotope',
  'num_explicit_hs',
  'node_id'],
 'edges_tuple': ((1, None, 0, 1),
  (1, None, 1, 2),
  (1, None, 1, 4),
  (1, None, 2, 3),
  (2, None, 4, 5)),
 'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']}


    assert get_mol_json(matches_dict['C'][0]) == C_mol_json
    assert matches_dict['C'][1] == ((3, 4, 5, 6, 19, 20), (18, 15, 14, 13, 16, 17))
    G = mol_to_nx(residues_0_mol)

    max_aa, max_aa_mol, max_match = match_molecular_graph_to_res_id(G, ResID='1', matches_dict=matches_dict)
    assert max_aa == 'C'
    assert max_match == (3, 4, 5, 6, 19, 20)
    #mol_json_to_mol(C_mol_json)
    #assert mols_are_identical(max_aa_mol, mol_json_to_mol(C_mol_json))
    assert get_mol_json(max_aa_mol) == C_mol_json



def test_decompose():
    mol_propagated, res_matches, modification_graphs = decompose(mol = residues_0_mol, cx_smarts_db = cx_smarts_db)
    assert [rdkit.Chem.MolToSmiles(nx_to_mol(G)) for G in modification_graphs] == ['CC=O', '[Na]PNCCBr']


    mol_propagated_json = {'nodes_tuple': [[None, None, None, 6, 0, 0, 4, False, 0, 0, 0],
  [None, None, None, 6, 0, 0, 3, False, 0, 0, 1],
  [None, None, None, 8, 0, 0, 3, False, 0, 0, 2],
  ['N', '1', 'C', 7, 0, 0, 3, False, 0, 0, 3],
  ['CA', '1', 'C', 6, 1, 0, 4, False, 0, 1, 4],
  [None, '1', 'C', 6, 0, 0, 4, False, 0, 0, 5],
  ['SG', '1', 'C', 16, 0, 0, 4, False, 0, 0, 6],
  [None, None, None, 6, 0, 0, 4, False, 0, 0, 7],
  [None, None, None, 35, 0, 0, 4, False, 0, 0, 8],
  [None, None, None, 6, 0, 0, 4, False, 0, 0, 9],
  [None, None, None, 7, 0, 0, 4, False, 0, 0, 10],
  [None, None, None, 15, 0, 0, 4, False, 0, 0, 11],
  [None, None, None, 11, 0, 0, 1, False, 0, 0, 12],
  ['SG', '3', 'C', 16, 0, 0, 4, False, 0, 0, 13],
  [None, '3', 'C', 6, 0, 0, 4, False, 0, 0, 14],
  ['CA', '3', 'C', 6, 1, 0, 4, False, 0, 1, 15],
  ['CO', '3', 'C', 6, 0, 0, 3, False, 0, 0, 16],
  ['O', '3', 'C', 8, 0, 0, 3, False, 0, 0, 17],
  ['N', '3', 'C', 7, 0, 0, 4, False, 0, 0, 18],
  ['CO', '1', 'C', 6, 0, 0, 3, False, 0, 0, 19],
  ['O', '1', 'C', 8, 0, 0, 3, False, 0, 0, 20]],
 'nodes_columns': ['AtomName',
  'ResID',
  'ResName',
  'atomic_num',
  'chiral_tag',
  'formal_charge',
  'hybridization',
  'is_aromatic',
  'isotope',
  'num_explicit_hs',
  'node_id'],
 'edges_tuple': ((1, None, 0, 1),
  (2, None, 1, 2),
  (1, None, 1, 3),
  (1, None, 3, 4),
  (1, None, 4, 5),
  (1, None, 4, 19),
  (1, None, 5, 6),
  (1, None, 6, 7),
  (1, None, 7, 8),
  (1, None, 7, 9),
  (1, None, 9, 10),
  (1, None, 10, 11),
  (1, None, 11, 12),
  (1, None, 11, 13),
  (1, None, 13, 14),
  (1, None, 14, 15),
  (1, None, 15, 16),
  (1, None, 15, 18),
  (2, None, 16, 17),
  (2, None, 19, 20)),
 'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']}

    assert nx_to_json(mol_propagated) == mol_propagated_json
    G = mol_propagated

    modification_graphs_nodes = [list(graph.nodes) for graph in modification_graphs]

    subgraph_tuples = get_subgraph_tuples(res_matches, modification_graphs_nodes, G)

    print('subgraph_tuples: ', subgraph_tuples)

    name_i, subgraph_i_nodes = subgraph_tuples[3]

    assert name_i == 'Mod_2'
    assert subgraph_i_nodes == [7,8,9,10,11,12]

    connections = [
        ('Res_1_C', 'Mod_1', [(1, 3)]),
        ('Res_1_C', 'Mod_2', [(6, 7)]),
        ('Res_3_C', 'Mod_2', [(11, 13)])]

    assert get_connections(list(G.edges), subgraph_tuples, ) == connections



    res_res_connections = []
    external_mod_connections = {
        '1': [('Res_1_C', [(1, 3)])],
        '2': [
            ('Res_1_C', [(6, 7)]),
            ('Res_3_C', [(11, 13)])]}
    
    assert split_connections_by_type(connections) == (res_res_connections, external_mod_connections)


    external_connections = [{'smiles': '[1*]C(C)=O',
  'max_attachment_point_id': 1,
  'attachment_points_on_sequence': {1: {'attachment_point_id': 1,
    'ResID': '1',
    'AtomName': 'N',
    'ResidueName': ''}}},
 {'smiles': '[2*]C(Br)CNP([3*])[Na]',
  'max_attachment_point_id': 3,
  'attachment_points_on_sequence': {2: {'attachment_point_id': 2,
    'ResID': '1',
    'AtomName': 'SG',
    'ResidueName': ''},
   3: {'attachment_point_id': 3,
    'ResID': '3',
    'AtomName': 'SG',
    'ResidueName': ''}}}]

    assert process_external_connections(
        external_mod_connections, res_matches, modification_graphs, G
    ) == external_connections

    assert {res_id: res_matches[res_id][0] for res_id in res_matches} == {'1': 'C', '3': 'C'}

    residues_2_ss_internal_connections = [{1: [{'ResID': '5', 'AtomName': 'SG', 'ResidueName': ''},
   {'ResID': '7', 'AtomName': 'SG', 'ResidueName': ''}]}]
    
    residues_2_ss_res_res_connections = [('Res_5_C', 'Res_7_C', [(6, 7)])]

    c_mol = rdkit.Chem.MolFromSmarts( '*CC(*)C=O' )
    residues_2_ss_res_matches = {'5': ('C', c_mol, (3, 4, 5, 6, 1, 2)),
        '7': ('C', c_mol, (0, 9, 8, 7, 10, 11))}
    
    residues_2_ss_json_propagated = {'nodes_tuple': [[7, 0, 0, 4, 0, False, 0, 'N', '7', 0],
        [6, 0, 0, 3, 0, False, 0, 'CO', '5', 1], [8, 0, 0, 3, 0, False, 0, 'O', '5', 2],
        [7, 0, 0, 4, 0, False, 0, 'N', '5', 3], [6, 0, 1, 4, 1, False, 0, 'CA', '5', 4],
        [6, 0, 0, 4, 0, False, 0, None, '5', 5], [16, 0, 0, 4, 0, False, 0, 'SG', '5', 6],
        [16, 0, 0, 4, 0, False, 0, 'SG', '7', 7], [6, 0, 0, 4, 0, False, 0, None, '7', 8],
        [6, 0, 1, 4, 1, False, 0, 'CA', '7', 9], [6, 0, 0, 3, 0, False, 0, 'CO', '7', 10],
        [8, 0, 0, 3, 0, False, 0, 'O', '7', 11]],
        'nodes_columns': ['atomic_num', 'formal_charge', 'chiral_tag', 'hybridization',
        'num_explicit_hs', 'is_aromatic', 'isotope', 'AtomName', 'ResID', 'node_id'],
        'edges_tuple': ((1, None, 0, 9), (1, None, 1, 4), (2, None, 1, 2), (1, None, 3, 4),
        (1, None, 4, 5), (1, None, 5, 6), (1, None, 6, 7), (1, None, 7, 8), (1, None, 8, 9),
        (1, None, 9, 10), (2, None, 10, 11)),
        'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']}
    residues_2_ss_graph_propagated = mol_json_to_nx(residues_2_ss_json_propagated)


    assert process_internal_connections(
        residues_2_ss_res_res_connections, residues_2_ss_res_matches, residues_2_ss_graph_propagated
    ) == residues_2_ss_internal_connections

    sequence_dict = {1:'A', 3:'D', 2:'C', }
    assert sequence_dict_to_string(sequence_dict) == 'ACD'
    return


