from typing import Any, Dict

from Peptide.models.Peptide import Peptide


def to_pepseq(peptide: Peptide, include_modifications: bool = False) -> str:
    """Encodes peptide to PepSeq string

    Args:
        peptide (Peptide): instance
        include_modifications (bool, optional):
        includes or exlude non-standard modifications. Defaults to True.

    Raises:
        Exception: _description_

    Returns:
        str: PepSeq
    """

    pepseq = peptide.peptide_json["pepseq_format"]
    return pepseq


def to_smiles(peptide: Peptide, include_modifications: bool = True) -> str:
    """Encodes peptide to canonical SMILES

    Args:
        peptide (Peptide): instance
        include_modifications (bool, optional):
        includes or exlude non-standard modifications. Defaults to True.

    Raises:
        Exception: _description_

    Returns:
        str: SMILES
    """
    return peptide.smiles


def to_json(peptide: Peptide, include_modifications: bool = True
            ) -> Dict[str, Any]:
    """Encodes peptide to JSON representation

    Args:
        peptide (Peptide): instance
        include_modifications (bool, optional):
        includes or exlude non-standard modifications. Defaults to True.

    Raises:
        Exception: _description_

    Returns:
        Dict[str, Any]: json

    Example output:
        {
            'sequence': 'CSCACGCK',
            'N_terminus': 'Ac',
            'C_terminus': 'NH2',
            'internal_modifications': {
                1: [
                    {
                        'ResID': '5',
                        'AtomName': 'SG',
                        'ResidueName': ''
                        },
                    {
                        'ResID': '7',
                        'AtomName': 'SG',
                        'ResidueName': ''
                        }
                    ]
                },
            'external_modifications': [
                {
                    'smiles': '[1*]C(Br)CNP([2*])[Na]',
                    'max_attachment_point_id': 2,
                    'attachment_points_on_sequence': {
                        1: {
                            'attachment_point_id': 1,
                            'ResID': '1',
                            'AtomName': 'SG',
                            'ResidueName': ''
                            },
                        2: {
                            'attachment_point_id': 2,
                            'ResID': '3',
                            'AtomName': 'SG',
                            'ResidueName': ''
                            }
                        }
                    }
                ]
            }

    """
    return peptide.peptide_json
