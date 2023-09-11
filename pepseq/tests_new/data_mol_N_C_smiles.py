mol_N_C_smiles_val = (
    "CC(=O)N[C@H]1CSC(Br)CNP([Na])SC[C@@H](C(=O)N[C@"
    + "@H](C)C(=O)N[C@H]2CSSC[C@@H](C(=O)N[C@@H](CCCCN)C(N)=O)NC(=O)CNC"
    + "2=O)NC(=O)[C@H](CO)NC1=O"
)



correct_peptide_json = {
    'sequence': 'CSCACGCK',
    'internal_modifications': [
        {
            1: [
                {'ResID': '5', 'AtomName': 'SG', 'ResidueName': ''},
                {'ResID': '7', 'AtomName': 'SG', 'ResidueName': ''}
                ]
            }
        ],
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
    }


tests = [
    (mol_N_C_smiles_val, correct_peptide_json),

]
