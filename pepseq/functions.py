import json
import os

from typing import Union, List

import pepseq.Peptide.utils.validation as validation

from pepseq.BuildPeptideJSONFromSMILES import \
    from_smiles_to_pepseq_and_mod_smiles_strings
from pepseq.get_peptide_json_from_pepseq_format import get_pep_json
from pepseq.read import from_json


absolute_path = os.path.dirname(__file__)
relative_db_path = "Peptide/database/db.json"
full_db_path = os.path.join(absolute_path, relative_db_path)

with open(full_db_path) as fp:
    db_json = json.load(fp)


def calculate_pepseq_and_mods(smiles: str) -> dict:
    """
    Calculate the peptide sequence and modification strings from a given SMILES string.

    Args:
        smiles (str): The SMILES string representing the peptide.

    Returns:
        dict: A dictionary containing the peptide sequence and modification strings.
    """
    return from_smiles_to_pepseq_and_mod_smiles_strings(smiles, db_json)


def validate(pepseq: str, smiles: List[str] = [], db: dict = db_json):
    """
    Validates a peptide sequence and its associated SMILES strings.

    Args:
        pepseq (str): The peptide sequence to validate.
        smiles (List[str], optional): List of SMILES strings associated with the peptide sequence. Defaults to [].
        db (dict, optional): The database to use for validation. Defaults to db_json.

    Returns:
        None
    """
    validation.validate(pepseq=pepseq, smiles=smiles)


def calculate(pepseq: str, smiles: list[str] = [], db: dict = None) -> dict:
    """
    Calculate various properties of a peptide sequence.

    Args:
        pepseq (str): The peptide sequence.
        smiles (list[str], optional): List of SMILES strings. Defaults to [].
        db (dict, optional): Database dictionary. Defaults to None.

    Returns:
        dict: A dictionary containing the calculated properties:
            - complete_smiles (str): The complete SMILES string.
            - length (int): The length of the peptide sequence.
            - mw (float): The molecular weight of the peptide.
            - sequence (str): The peptide sequence.
    """
    if smiles == []:
        smiles = None
    if db is None:
        db = db_json

    peptide_json = get_pep_json(pepseq, db, smiles)
    peptide = from_json(peptide_json, db)
    complete_smiles = peptide.complete_smiles
    sequence = peptide.sequence
    length = peptide.length
    mw = peptide.mw

    return {
        "complete_smiles": complete_smiles,
        "length": length,
        "mw": mw,
        "sequence": sequence,
    }
