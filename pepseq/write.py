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


def to_json(peptide: Peptide, include_modifications: bool = True) -> Dict[str, Any]:
    """Encodes peptide to JSON representation

    :parameter peptide: Peptide instance
    :parameter include_modifications: includes or exlude non-standard modifications. Defaults to True.

    :return JSON representation of peptide

    """
    return peptide.peptide_json
