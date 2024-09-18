import rdkit
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import (mol_to_nx, nx_to_mol,
     nx_to_json, mol_json_to_nx, mol_json_to_mol)
from pepseq.tests_new.helpers import mols_are_identical, smiles_are_identical


def test_mol_to_nx_to_json():
    m = rdkit.Chem.MolFromSmiles('[*:1]CNCC[*:2]')
    G = mol_to_nx(m)
    mol_json = nx_to_json(G)

    mol_json_fixture = {'nodes_tuple': [[0, 0, '*', 0, None, False, 0, '1', 0, 0, 0],
        [6, 0, None, 0, 4, False, 0, None, 1, 0, 1],
        [7, 0, None, 0, 4, False, 0, None, 2, 0, 2],
        [6, 0, None, 0, 4, False, 0, None, 3, 0, 3],
        [6, 0, None, 0, 4, False, 0, None, 4, 0, 4],
        [0, 0, '*', 0, None, False, 0, '2', 5, 0, 5]],
        'nodes_columns': ['atomic_num',
        'chiral_tag',
        'dummyLabel',
        'formal_charge',
        'hybridization',
        'is_aromatic',
        'isotope',
        'molAtomMapNumber',
        'node_id',
        'num_explicit_hs'],
        'edges_tuple': ((1, None, 0, 1),
        (1, None, 1, 2),
        (1, None, 2, 3),
        (1, None, 3, 4),
        (1, None, 4, 5)),
        'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']}

    # we can include other columns like is_connection_to_mod

    mol_json_to_nx(mol_json_fixture)
    mol = mol_json_to_mol(mol_json_fixture)
    smi = rdkit.Chem.MolToSmiles(mol)
    assert smiles_are_identical(smi, 'C(C[*:2])NC[*:1]')


