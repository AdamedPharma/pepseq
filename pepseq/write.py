from typing import Any, Dict

from pepseq.Peptide.models.Peptide import Peptide


def to_pepseq(peptide: Peptide, include_modifications: bool = False) -> str:
    """Encodes peptide to PepSeq string

    :param peptide:  Peptide instance
    :type  peptide: Peptide
    :param include_modifications: includes or exlude non-standard modifications. Defaults to True.
    :type  include_modifications: bool


    :return: PepSeq string
    :rtype: str

    Raises:
        Exception: _description_
    """

    pepseq = peptide.peptide_json["pepseq_format"]
    return pepseq


def to_smiles(peptide: Peptide, include_modifications: bool = True) -> str:
    """Encodes peptide to canonical SMILES

    :param peptide: Peptide instance
    :type  peptide: Peptide
    :param include_modifications: includes or exlude non-standard modifications. Defaults to True.
    :type  include_modifications: bool

    :return: SMILES
    :rtype:  str

    Raises:
        Exception: _description_

    """

    return peptide.smiles


def to_json(peptide: Peptide, include_modifications: bool = True) -> Dict[str, Any]:
    """Encodes peptide to JSON representation

    :param peptide: Peptide instance
    :type  peptide: Peptide
    :param include_modifications: includes or exlude non-standard modifications. Defaults to True.
    :type  include_modifications: bool


    :return: JSON representation of peptide
    :rtype: Dict[str, Any]

    Raises:
        Exception: _description_

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
