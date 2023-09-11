import pkgutil
import json
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import mol_json_to_nx


residues_path = pkgutil.extend_path("residues.json", __name__)
with open(residues_path) as fp:
    residues_json = json.load(fp)


residues_graphs = [mol_json_to_nx(residue_json) for residue_json in residues_json]


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


sequence = "CSCACGCK"
datarow = residues_graphs, sequence, correct_internal_modifications, correct_external_modifications


tests = [
    datarow,
    ]