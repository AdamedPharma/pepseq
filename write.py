from Peptide.models.Peptide import Peptide
from typing import Dict, Any


def to_pepseq(peptide: Peptide, include_modifications: bool = False) -> str:
    """Encodes peptide to PepSeq string

    Args:
        peptide (Peptide): instance
        include_modifications (bool, optional): includes or exlude non-standard modifications. Defaults to True.

    Raises:
        Exception: _description_

    Returns:
        str: PepSeq
    """


def to_smiles(peptide: Peptide, include_modifications: bool = True) -> str:
    """Encodes peptide to canonical SMILES

    Args:
        peptide (Peptide): instance
        include_modifications (bool, optional): includes or exlude non-standard modifications. Defaults to True.

    Raises:
        Exception: _description_

    Returns:
        str: SMILES
    """
    raise Exception("TO BE IMPLEMENTED ")


def to_json(peptide: Peptide, include_modifications: bool = True) -> Dict[str, Any]:
    """Encodes peptide to JSON representation

    Args:
        peptide (Peptide): instance
        include_modifications (bool, optional): includes or exlude non-standard modifications. Defaults to True.

    Raises:
        Exception: _description_

    Returns:
        Dict[str, Any]: json

    Example output:
        {"sequence":"H{Aib}EGTFTSDVSSYLEGQAAKEFIAWLVRGRG",
        "modifications":[
            {
            "smiles": "[*1]CC(=O)NCCCC(NC(C)=O)C(=O)Nc1ccc2oc(=O)cc(CC(=O)NCCOCCOCCC(=O)NCCCCC(NC(=O)CCCCCCCCCCCCCCCCC(=O)O)C(=O)O)c2c1[*2]",
            "connecting_residues": [17, 24]
            }
        ]
    }
    where [*1] and [*2] are attachement points
    """
    raise Exception("TO BE IMPLEMENTED ")
