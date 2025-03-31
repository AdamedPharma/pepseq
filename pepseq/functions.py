import json
import os

from typing import List, Union

import pepseq.Peptide.utils.validation as validation

from pepseq.BuildPeptideJSONFromSMILES import (
    from_smiles_to_pepseq_and_mod_smiles_strings,
)
from pepseq.get_peptide_json_from_pepseq_format import get_pep_json
from pepseq.read import from_json


absolute_path = os.path.dirname(__file__)
relative_db_path = "Peptide/database/db.json"
full_db_path = os.path.join(absolute_path, relative_db_path)

with open(full_db_path) as fp:
    db_json = json.load(fp)


def calculate_pepseq_and_mods(smiles: str) -> dict:
    """
    Calculate peptide sequence and modifications from SMILES string

    :param smiles: SMILES string
    :type smiles: str

    :return: Dictionary containing peptide sequence and modifications
    :rtype: dict
    """
    return from_smiles_to_pepseq_and_mod_smiles_strings(smiles, db_json)


def validate(pepseq: str, smiles: List[str] = [], db: dict = db_json):
    """
    Validate peptide sequence and SMILES string

    :param pepseq: Peptide sequence
    :type pepseq: str

    :param smiles: SMILES string
    :type smiles: list[str]

    :param db: Database containing information about amino acids and modifications
    :type db: dict

    :return: None
    :rtype: None
    """
    validation.validate(pepseq=pepseq, smiles=smiles)


def calculate(pepseq: str, smiles: list[str] = [], db: Union[dict, None] = None) -> dict:
    """
    Calculate properties of peptide sequence

    :param pepseq: Peptide sequence
    :type pepseq: str

    :param smiles: SMILES string
    :type smiles: list[str]

    :param db: Database containing information about amino acids and modifications
    :type db: dict

    :return: Dictionary containing properties of peptide sequence
    :rtype: dict
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
