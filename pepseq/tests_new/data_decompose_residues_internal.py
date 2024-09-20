import json
import pkgutil


from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import (mol_to_nx, nx_to_mol,
     nx_to_json, mol_json_to_nx, mol_json_to_mol)


from pepseq.BuildPeptideJSONFromSMILES import  get_cx_smarts_db


db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)

cx_smarts_db = get_cx_smarts_db(db_json)


residues_jsons_fx = [
        {
            'nodes_tuple': [
                ['N', '1', 7, 0, 0, 4, False, 0, 0, 1, 0], ['CA', '1', 6, 1, 0, 4, False, 0, 1, 1, 1],
                [None, None, 6, 0, 0, 4, False, 0, 2, 0, 2], [None, None, 16, 0, 0, 4, False, 0, 3, 0, 3],
                ['CO', '1', 6, 0, 0, 3, False, 0, 130, 0, 130], ['O', '1', 8, 0, 0, 3, False, 0, 131, 0, 131],
                [None, None, 6, 0, 0, 4, False, 0, 4, 0, 4], [None, None, 7, 0, 0, 4, False, 0, 5, 0, 5],
                [None, None, 6, 0, 0, 4, False, 0, 6, 0, 6], [None, None, 6, 0, 0, 4, False, 0, 7, 0, 7],
                [None, None, 16, 0, 0, 4, False, 0, 8, 0, 8], [None, None, 6, 0, 0, 4, False, 0, 9, 0, 9],
                ['CA', '12', 6, 1, 0, 4, False, 0, 10, 1, 10], ['CO', '12', 6, 0, 0, 3, False, 0, 11, 0, 11],
                ['O', '12', 8, 0, 0, 3, False, 0, 12, 0, 12], ['N', '12', 7, 0, 0, 3, False, 0, 58, 0, 58]
                ],
            'nodes_columns': [
                'AtomName', 'ResID', 'atomic_num', 'chiral_tag', 'formal_charge', 'hybridization',
                'is_aromatic', 'isotope', 'node_id', 'num_explicit_hs'
                ],
                'edges_tuple': (
                    (1, None, 0, 1), (1, None, 1, 2), (1, None, 1, 130), (1, None, 2, 3),
                    (1, None, 3, 4), (2, None, 130, 131), (1, None, 4, 5), (1, None, 5, 6),
                    (1, None, 6, 7), (1, None, 7, 8), (1, None, 8, 9), (1, None, 9, 10),
                    (1, None, 10, 11), (1, None, 10, 58), (2, None, 11, 12)
                    ),
                'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']
            },
        {
            'nodes_tuple': [
                ['O', '13', 8, 0, 0, 3, False, 0, 16, 0, 16], ['N', '13', 7, 0, 0, 3, False, 0, 13, 0, 13],
                ['CA', '13', 6, 0, 0, 4, False, 0, 14, 0, 14], ['CO', '13', 6, 0, 0, 3, False, 0, 15, 0, 15]
                ],
            'nodes_columns': [
                'AtomName', 'ResID', 'atomic_num', 'chiral_tag', 'formal_charge', 'hybridization',
                'is_aromatic', 'isotope', 'node_id', 'num_explicit_hs'
                ],
            'edges_tuple': ((2, None, 16, 15), (1, None, 13, 14), (1, None, 14, 15)),
            'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
            {'nodes_tuple': [
                ['N', '14', 7, 0, 0, 3, False, 0, 17, 0, 17], ['CA', '14', 6, 1, 0, 4, False, 0, 18, 1, 18],
                [None, None, 6, 0, 0, 4, False, 0, 19, 0, 19], [None, None, 16, 0, 0, 4, False, 0, 20, 0, 20],
                [None, None, 6, 0, 0, 4, False, 0, 21, 0, 21], [None, None, 7, 0, 0, 4, False, 0, 22, 0, 22],
                [None, None, 6, 0, 0, 4, False, 0, 23, 0, 23], [None, None, 6, 0, 0, 4, False, 0, 24, 0, 24],
                [None, None, 16, 0, 0, 4, False, 0, 25, 0, 25], [None, None, 15, 0, 0, 4, False, 0, 26, 0, 26],
                ['CO', '14', 6, 0, 0, 3, False, 0, 27, 0, 27], ['O', '14', 8, 0, 0, 3, False, 0, 28, 0, 28]],
                'nodes_columns': [
                    'AtomName', 'ResID', 'atomic_num', 'chiral_tag', 'formal_charge', 'hybridization',
                    'is_aromatic', 'isotope', 'node_id', 'num_explicit_hs'
                    ],
                'edges_tuple': (
                    (1, None, 17, 18), (1, None, 18, 19), (1, None, 18, 27), (1, None, 19, 20),
                    (1, None, 20, 21), (1, None, 21, 22), (1, None, 22, 23), (1, None, 23, 24),
                    (1, None, 24, 25), (1, None, 25, 26), (2, None, 27, 28)
                    ),
                'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']
                },
            {'nodes_tuple': [
                [None, None, 6, 0, 0, 3, False, 0, 32, 0, 32], [None, None, 8, 0, 0, 3, False, 0, 33, 0, 33],
                [None, None, 8, 0, 0, 3, False, 0, 34, 0, 34], ['CO', '15', 6, 0, 0, 3, False, 0, 35, 0, 35],
                ['O', '15', 8, 0, 0, 3, False, 0, 36, 0, 36], ['N', '15', 7, 0, 0, 3, False, 0, 29, 0, 29],
                ['CA', '15', 6, 1, 0, 4, False, 0, 30, 1, 30], [None, None, 6, 0, 0, 4, False, 0, 31, 0, 31]],
            'nodes_columns': [
                'AtomName',   'ResID',   'atomic_num',   'chiral_tag',   'formal_charge',
                'hybridization',   'is_aromatic',   'isotope',   'node_id',   'num_explicit_hs'],
            'edges_tuple': (
                (1, None, 32, 31),   (2, None, 32, 33),   (1, None, 32, 34),   (1, None, 35, 30),
                (2, None, 35, 36),   (1, None, 29, 30),   (1, None, 30, 31)),
            'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
            {'nodes_tuple': [
                ['N', '16', 7, 0, 0, 3, False, 0, 37, 0, 37],   ['CA', '16', 6, 1, 0, 4, False, 0, 38, 1, 38],
                [None, None, 6, 0, 0, 4, False, 0, 39, 0, 39],   [None, None, 6, 0, 0, 4, False, 0, 40, 0, 40],
                [None, None, 6, 0, 0, 3, False, 0, 41, 0, 41],   [None, None, 8, 0, 0, 3, False, 0, 42, 0, 42],
                [None, None, 8, 0, 0, 3, False, 0, 43, 0, 43], ['CO', '16', 6, 0, 0, 3, False, 0, 44, 0, 44],
                ['O', '16', 8, 0, 0, 3, False, 0, 45, 0, 45]],
            'nodes_columns': [
                'AtomName', 'ResID', 'atomic_num', 'chiral_tag', 'formal_charge',
                'hybridization', 'is_aromatic', 'isotope', 'node_id', 'num_explicit_hs'],
            'edges_tuple': (
                (1, None, 37, 38), (1, None, 38, 39), (1, None, 38, 44), (1, None, 39, 40),
                (1, None, 40, 41), (2, None, 41, 42), (1, None, 41, 43), (2, None, 44, 45)),
            'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
            {'nodes_tuple': [
                ['N', '17', 7, 0, 0, 3, False, 0, 46, 0, 46], ['CA', '17', 6, 1, 0, 4, False, 0, 47, 1, 47],
                [None, None, 6, 0, 0, 4, False, 0, 48, 0, 48], [None, None, 6, 0, 0, 3, True, 0, 49, 0, 49],
                [None, None, 6, 0, 0, 3, True, 0, 50, 0, 50], [None, None, 6, 0, 0, 3, True, 0, 51, 0, 51],
                [None, None, 6, 0, 0, 3, True, 0, 52, 0, 52], [None, None, 6, 0, 0, 3, True, 0, 53, 0, 53],
                [None, None, 6, 0, 0, 3, True, 0, 54, 0, 54], ['CO', '17', 6, 0, 0, 3, False, 0, 55, 0, 55],
                ['O', '17', 8, 0, 0, 3, False, 0, 56, 0, 56], [None, None, 8, 0, 0, 3, False, 0, 57, 0, 57]],
            'nodes_columns': [
                'AtomName', 'ResID', 'atomic_num', 'chiral_tag', 'formal_charge', 'hybridization', 'is_aromatic',
                'isotope', 'node_id', 'num_explicit_hs'],
            'edges_tuple': (
                (1, None, 46, 47), (1, None, 47, 48), (1, None, 47, 55), (1, None, 48, 49), (12, None, 49, 50),
                (12, None, 49, 54), (12, None, 50, 51), (12, None, 51, 52), (12, None, 52, 53), (12, None, 53, 54),
                (2, None, 55, 56), (1, None, 55, 57)),
            'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
            {'nodes_tuple': [
                [None, None, 6, 0, 0, 3, False, 0, 64, 0, 64], [None, None, 7, 0, 0, 3, False, 0, 65, 0, 65],
                [None, None, 8, 0, 0, 3, False, 0, 66, 0, 66], ['N', '11', 7, 0, 0, 3, False, 0, 67, 0, 67],
                ['CO', '11', 6, 0, 0, 3, False, 0, 59, 0, 59], ['O', '11', 8, 0, 0, 3, False, 0, 60, 0, 60],
                ['CA', '11', 6, -1, 0, 4, False, 0, 61, 1, 61], [None, None, 6, 0, 0, 4, False, 0, 62, 0, 62],
                [None, None, 6, 0, 0, 4, False, 0, 63, 0, 63]],
            'nodes_columns': [
                'AtomName', 'ResID', 'atomic_num', 'chiral_tag', 'formal_charge', 'hybridization',
                'is_aromatic', 'isotope', 'node_id', 'num_explicit_hs'],
            'edges_tuple': ((1, None, 64, 63), (1, None, 64, 65), (2, None, 64, 66),
            (1, None, 67, 61), (2, None, 59, 60), (1, None, 59, 61), (1, None, 61, 62),
            (1, None, 62, 63)),
            'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
            {'nodes_tuple': [
                ['CO', '10', 6, 0, 0, 3, False, 0, 68, 0, 68], ['O', '10', 8, 0, 0, 3, False, 0, 69, 0, 69],
                ['CA', '10', 6, -1, 0, 4, False, 0, 70, 1, 70], [None, None, 6, 0, 0, 4, False, 0, 71, 0, 71],
                [None, None, 6, 0, 0, 4, False, 0, 72, 0, 72], [None, None, 6, 0, 0, 3, False, 0, 73, 0, 73],
                [None, None, 8, 0, 0, 3, False, 0, 74, 0, 74], [None, None, 8, 0, 0, 3, False, 0, 75, 0, 75],
                ['N', '10', 7, 0, 0, 3, False, 0, 76, 0, 76]],
            'nodes_columns': [
                'AtomName', 'ResID', 'atomic_num', 'chiral_tag', 'formal_charge', 'hybridization', 'is_aromatic',
                'isotope', 'node_id', 'num_explicit_hs'],
            'edges_tuple': ((2, None, 68, 69), (1, None, 68, 70), (1, None, 70, 71), (1, None, 70, 76),
            (1, None, 71, 72), (1, None, 72, 73), (2, None, 73, 74), (1, None, 73, 75)),
            'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
            {'nodes_tuple': [['CO', '9', 6, 0, 0, 3, False, 0, 77, 0, 77], ['O', '9', 8, 0, 0, 3, False, 0, 78, 0, 78],
            ['CA', '9', 6, 1, 0, 4, False, 0, 79, 1, 79], [None, None, 6, 0, 0, 4, False, 0, 80, 0, 80],
            [None, None, 8, 0, 0, 4, False, 0, 81, 0, 81], ['N', '9', 7, 0, 0, 3, False, 0, 82, 0, 82]],
            'nodes_columns': ['AtomName', 'ResID', 'atomic_num', 'chiral_tag', 'formal_charge',
            'hybridization', 'is_aromatic', 'isotope', 'node_id', 'num_explicit_hs'],
            'edges_tuple': ((2, None, 77, 78), (1, None, 77, 79), (1, None, 79, 80), (1, None, 79, 82),
            (1, None, 80, 81)),
            'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
            {'nodes_tuple': [['CO', '8', 6, 0, 0, 3, False, 0, 83, 0, 83], ['O', '8', 8, 0, 0, 3, False, 0, 84, 0, 84],
            ['CA', '8', 6, -1, 0, 4, False, 0, 85, 1, 85], [None, None, 6, 0, 0, 4, False, 0, 86, 0, 86],
            [None, None, 6, 0, 0, 4, False, 0, 87, 0, 87], [None, None, 6, 0, 0, 4, False, 0, 88, 0, 88],
            ['N', '8', 7, 0, 0, 3, False, 0, 89, 0, 89]],
            'nodes_columns': ['AtomName', 'ResID', 'atomic_num', 'chiral_tag', 'formal_charge', 'hybridization',
            'is_aromatic', 'isotope', 'node_id', 'num_explicit_hs'],
            'edges_tuple': ((2, None, 83, 84), (1, None, 83, 85), (1, None, 85, 86), (1, None, 85, 89), (1, None, 86, 87),
            (1, None, 87, 88), (1, None, 88, 89)),
            'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
            {'nodes_tuple': [[None, None, 8, 0, 0, 3, False, 0, 96, 0, 96], [None, None, 8, 0, 0, 3, False, 0, 97, 0, 97],
            ['N', '7', 7, 0, 0, 3, False, 0, 98, 0, 98], ['CO', '7', 6, 0, 0, 3, False, 0, 90, 0, 90],
            ['O', '7', 8, 0, 0, 3, False, 0, 91, 0, 91], ['CA', '7', 6, -1, 0, 4, False, 0, 92, 1, 92],
            [None, None, 6, 0, 0, 4, False, 0, 93, 0, 93], [None, None, 6, 0, 0, 4, False, 0, 94, 0, 94],
            [None, None, 6, 0, 0, 3, False, 0, 95, 0, 95]],
            'nodes_columns': ['AtomName', 'ResID', 'atomic_num', 'chiral_tag', 'formal_charge', 'hybridization', 'is_aromatic',
            'isotope', 'node_id', 'num_explicit_hs'],
            'edges_tuple': ((2, None, 96, 95), (1, None, 97, 95), (1, None, 98, 92), (2, None, 90, 91), (1, None, 90, 92),
            (1, None, 92, 93), (1, None, 93, 94), (1, None, 94, 95)),
            'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
            {'nodes_tuple': [['CO', '6', 6, 0, 0, 3, False, 0, 99, 0, 99], ['O', '6', 8, 0, 0, 3, False, 0, 100, 0, 100],
            ['CA', '6', 6, -1, 0, 4, False, 0, 101, 1, 101], [None, None, 6, 0, 0, 4, False, 0, 102, 0, 102],
            [None, None, 6, 0, 0, 4, False, 0, 103, 0, 103], [None, None, 6, 0, 0, 4, False, 0, 104, 0, 104],
            ['N', '6', 7, 0, 0, 3, False, 0, 105, 0, 105]],
            'nodes_columns': ['AtomName', 'ResID', 'atomic_num', 'chiral_tag', 'formal_charge', 'hybridization', 'is_aromatic',
            'isotope', 'node_id', 'num_explicit_hs'],
  'edges_tuple': ((2, None, 99, 100),   (1, None, 99, 101),   (1, None, 101, 102),
   (1, None, 101, 105),   (1, None, 102, 103),   (1, None, 103, 104),   (1, None, 104, 105)),
  'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
 {'nodes_tuple': [['CO', '5', 6, 0, 0, 3, False, 0, 106, 0, 106],
   ['O', '5', 8, 0, 0, 3, False, 0, 107, 0, 107],   ['CA', '5', 6, -1, 0, 4, False, 0, 108, 1, 108],
   [None, None, 6, 0, 0, 4, False, 0, 109, 0, 109],   ['N', '5', 7, 0, 0, 3, False, 0, 110, 0, 110]],
  'nodes_columns': ['AtomName',   'ResID',   'atomic_num',   'chiral_tag',   'formal_charge',
   'hybridization',   'is_aromatic',   'isotope',   'node_id',   'num_explicit_hs'],
  'edges_tuple': ((2, None, 106, 107),   (1, None, 106, 108),   (1, None, 108, 109),   (1, None, 108, 110)),
  'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
 {'nodes_tuple': [['CO', '4', 6, 0, 0, 3, False, 0, 111, 0, 111], ['O', '4', 8, 0, 0, 3, False, 0, 112, 0, 112],
      ['CA', '4', 6, -1, 0, 4, False, 0, 113, 1, 113],   [None, None, 6, 0, 0, 4, False, 0, 114, 0, 114],
   [None, None, 6, 0, 0, 3, False, 0, 115, 0, 115],   [None, None, 8, 0, 0, 3, False, 0, 116, 0, 116],
   [None, None, 8, 0, 0, 3, False, 0, 117, 0, 117],   ['N', '4', 7, 0, 0, 3, False, 0, 118, 0, 118]],
  'nodes_columns': ['AtomName',   'ResID',   'atomic_num',   'chiral_tag',   'formal_charge',
   'hybridization',   'is_aromatic',   'isotope',   'node_id',   'num_explicit_hs'],
  'edges_tuple': ((2, None, 111, 112),   (1, None, 111, 113),   (1, None, 113, 114),
   (1, None, 113, 118),   (1, None, 114, 115),   (2, None, 115, 116),   (1, None, 115, 117)),
  'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
 {'nodes_tuple': [['CO', '3', 6, 0, 0, 3, False, 0, 119, 0, 119],   ['O', '3', 8, 0, 0, 3, False, 0, 120, 0, 120],
   ['CA', '3', 6, -1, 0, 4, False, 0, 121, 1, 121],   [None, None, 6, 0, 0, 4, False, 0, 122, 0, 122],
   [None, None, 16, 0, 0, 4, False, 0, 123, 0, 123],   ['N', '3', 7, 0, 0, 3, False, 0, 124, 0, 124]],
  'nodes_columns': ['AtomName',   'ResID',   'atomic_num',   'chiral_tag',   'formal_charge',
   'hybridization',   'is_aromatic',   'isotope',   'node_id',   'num_explicit_hs'],
  'edges_tuple': ((2, None, 119, 120),   (1, None, 119, 121),   (1, None, 121, 122),
   (1, None, 121, 124),   (1, None, 122, 123)),
  'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']},
 {'nodes_tuple': [[None, None, 6, 0, 0, 4, False, 0, 128, 0, 128],   ['N', '2', 7, 0, 0, 3, False, 0, 129, 0, 129],
   ['CO', '2', 6, 0, 0, 3, False, 0, 125, 0, 125],   ['O', '2', 8, 0, 0, 3, False, 0, 126, 0, 126],
   ['CA', '2', 6, -1, 0, 4, False, 0, 127, 1, 127]],
  'nodes_columns': ['AtomName',   'ResID',   'atomic_num',   'chiral_tag',   'formal_charge',
   'hybridization',   'is_aromatic',   'isotope',   'node_id',   'num_explicit_hs'],
  'edges_tuple': ((1, None, 128, 127),   (1, None, 129, 127),   (2, None, 125, 126),   (1, None, 125, 127)),
  'edges_columns': ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']}]

fragments = [mol_json_to_nx(i) for i in residues_jsons_fx]

n_subst_limit=None

ketcher=True

tests = [
    (
        (fragments, cx_smarts_db, n_subst_limit, ketcher)
    ),
]