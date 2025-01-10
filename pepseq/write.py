from typing import Any, Dict

from pepseq.Peptide.models.Peptide import Peptide


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
    """
    Encodes peptide to canonical SMILES

    :param: peptide: Peptide instance
    :type: peptide: Peptide

    :param: include_modifications includes or exlude non-standard modifications. Defaults to True.
    :type: peptide: (bool, optional)

    :return: SMILES
    :rtype: str

    Raises:
        Exception: _description_
    """
    return peptide.smiles


def to_json(peptide: Peptide, include_modifications: bool = True) -> Dict[str, Any]:
    """Encodes peptide to JSON representation

    :param: peptide: Peptide instance
    :type: peptide: Peptide

    :param: include_modifications includes or exlude non-standard modifications. Defaults to True.
    :type: peptide: (bool, optional)

    :return JSON representation of peptide
    :rtype: dict
    """
    return peptide.peptide_json
