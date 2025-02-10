from typing import Union, List

import os
import json
import rdkit

from pathlib import Path


from pepseq.get_peptide_json_from_pepseq_format import (
    get_attachment_points_on_sequence_json,
    get_pep_json,
)
from pepseq.Peptide.utils.pepseq_validation import validate_pepseq
from pepseq.Peptide.utils.smiles_validation import validate_smiles_codes

from pepseq.Peptide.exceptions import (
    AttachmentPointsMismatchError,
    AttachmentPointsNonUniqueError,
)


absolute_path = Path(__file__).parent.parent.parent
relative_db_path = "Peptide/database/db.json"
full_db_path = os.path.join(absolute_path, relative_db_path)

with open(full_db_path) as fp:
    db_json = json.load(fp)


def get_attachment_points_on_smiles(smiles_code: str) -> list:
    """
    Retrieves the attachment points on a SMILES code.

    :param smiles_code: The SMILES code.
    :type smiles_code: str

    :return: A list of attachment point IDs.
    :rtype: list
    """
    attachment_points_ids = []
    mol = rdkit.Chem.MolFromSmiles(smiles_code)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 0:
            isotope = atom.GetAtomMapNum()
            attachment_points_ids.append(isotope)
    return attachment_points_ids


def get_attachment_points_on_smiles_codes(
    smiles_codes: Union[list, None] = None
) -> set:
    """
    Retrieves the attachment points on SMILES codes.

    :param smiles_codes: List of SMILES codes.Defaults to None.
    :type smiles_codes: Union[list, None], optional

    :return: Set of attachment point IDs.
    :rtype: set

    :raises: AttachmentPointsNonUniqueError: If attachment point labels on SMILES are not unique.
    """
    attachment_points_ids = []
    if smiles_codes is not None:
        for smiles_code in smiles_codes:
            attachment_points_ids += get_attachment_points_on_smiles(smiles_code)

    unique_attachment_points_ids = set(attachment_points_ids)
    if len(attachment_points_ids) > len(unique_attachment_points_ids):
        raise AttachmentPointsNonUniqueError(
            "Attachment Points labels on SMILES are not unique."
        )

    return unique_attachment_points_ids


def validate_matching_attachment_points(pepseq: str, smiles_codes: list):
    """
    Validates if the attachment points on a peptide sequence match the attachment points on SMILES codes.

    :param pepseq: The peptide sequence.
    :type pepseq: str

    :param smiles_codes: List of SMILES codes.
    :type smiles_codes: list


    :return: None

    :raises: AttachmentPointsMismatchError: If the attachment points on the sequence
        do not match the attachment points on the SMILES codes.


    TO BE IMPROVED: the sequence needs to be parsed into residue symbols
      in order to validate the symbols agreement with SMILES codes;
    this forces us to make multiple calls to get_pep_json function just in order to validate;
    In general some features can be validated only after extracting them.

    """
    symbols = get_pep_json(pepseq)["symbols"]
    attachment_points_on_sequence = get_attachment_points_on_sequence_json(symbols)
    attachment_point_ids_on_sequence = set(attachment_points_on_sequence.keys())
    attachment_point_ids_on_smiles = get_attachment_points_on_smiles_codes(smiles_codes)

    if attachment_point_ids_on_sequence != attachment_point_ids_on_smiles:
        raise AttachmentPointsMismatchError(
            "Attachment Points on Sequence: %s do not Match Attachment Points on Smiles: %s"
            % (
                str(attachment_point_ids_on_sequence),
                str(attachment_point_ids_on_smiles),
            )
        )


def validate(pepseq: str, smiles: List[str] = [], db: dict = db_json):
    """
    Validate the system of smiles and pepseq (however it might be problematic
     I guess because we wanted to separate validating sequence in pepseq format
     from validating SMILES and then validate them together.

    :param pepseq – obligatory parameter pepseq in form like
        CSCACGCK or {CH3}-CSCACGCK-{NH2} or CS{Cys(R1)}GACG~NH2
    :type pepseq: str

    :param smiles – list of smiles codes that can be empty or full; can
    :type    smiles: List[str]

    :param db – database of monomers
    :type db: dict

    :return: None

    :raises: InvalidSymbolError: If any of the monomers are not found in the database.
    :raises: InvalidSmilesError: If any of the SMILES codes is invalid and cannot be constructed into a molecule.
    :raises: UnattachedSmilesError: If any of the SMILES codes do not have an attachment point to Peptide.
    :raises: AttachmentPointsMismatchError: If the attachment points on the
     sequence do not match the attachment points on the SMILES codes.
    :raises: AttachmentPointsNonUniqueError: If attachment point labels on SMILES are not unique.
    """
    validate_pepseq(pepseq, db)
    validate_smiles_codes(smiles)
    validate_matching_attachment_points(pepseq, smiles)
    return
