"""
External Modifications for Peptide are provided as list of SMILES strings

SMILES codes
"""

from typing import Union
import rdkit
from pepseq.Peptide.exceptions import InvalidSmilesError, UnattachedSmilesError


def has_attachment_point(smiles: str) -> bool:
    """
    Check if a given SMILES string has an attachment point.

    Args:
        smiles (str): The SMILES string to check.

    Returns:
        bool: True if the SMILES string has an attachment point, False otherwise.
    """
    mol = rdkit.Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            return True
    return False


def validate_attachment_points_on_smiles(smiles_codes: list[str]):
    """
    Validates the attachment points on a list of SMILES codes.

    Args:
        smiles_codes (list[str]): A list of SMILES codes to validate.

    Raises:
        UnattachedSmilesError: If any of the SMILES codes do not have an attachment point to Peptide.
    """
    invalid_ids = []
    if smiles_codes is not None:
        for i in range(len(smiles_codes)):
            smiles_code = smiles_codes[i]
            if not has_attachment_point(smiles_code):
                invalid_ids.append(i + 1)
    if invalid_ids:
        ErrorMessage = "\n".join(
            [
                "SMILES code no %d has no attachment point to Peptide" % invalid_id
                for invalid_id in invalid_ids
            ]
        )
        raise UnattachedSmilesError(ErrorMessage)


def can_create_rdkit_molecule(smiles: str) -> bool:
    """
    Check if a molecule can be created from a SMILES string.

    Args:
        smiles (str): The SMILES string to check.

    Returns:
        bool: True if a molecule can be created from the SMILES string, False otherwise.
    """
    mol = rdkit.Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    return True


def validate_structure_by_rdkit(smiles_codes: Union[list[str], None] = None):
    """
    This validates SMILES by rdkit
    Validates a list of SMILES codes.

    Args:
        smiles_codes (list[str], None): A list of SMILES codes to be validated.

    Raises:
        InvalidSmilesError: If any of the SMILES codes is invalid and cannot be constructed into a molecule.

    Returns:
        None
    """
    if smiles_codes is not None:
        for i in range(len(smiles_codes)):
            if not can_create_rdkit_molecule(smiles_codes[i]):
                raise InvalidSmilesError(
                    "SMILES code no %d was invalid. Could not construct molecule from SMILES code"
                    % (i + 1)
                )


def validate_smiles_codes(smiles_codes: Union[list[str], None] = None):
    """
    Validates a list of SMILES codes.

    Args:
        smiles_codes (list[str], None): A list of SMILES codes to be validated.

    Raises:
        InvalidSmilesError: If any of the SMILES codes is invalid and cannot be constructed into a molecule.
        UnattachedSmilesError: If any of the SMILES codes do not have an attachment point to Peptide.
    """
    validate_structure_by_rdkit(smiles_codes)
    validate_attachment_points_on_smiles(smiles_codes)
