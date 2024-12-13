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
    return from_smiles_to_pepseq_and_mod_smiles_strings(smiles, db_json)


def validate (pepseq: str, smiles: List[str] = [], db: dict = db_json):
    validation.validate(pepseq=pepseq, smiles=smiles)



def calculate(pepseq: str, smiles: list[str] = [], db: dict = None, ketcher: bool = False) -> dict:
    if smiles == []:
        smiles = None
    if db is None:
        db = db_json

    peptide_json = get_pep_json(pepseq, db, smiles, ketcher=ketcher)
    peptide = from_json(peptide_json, db, ketcher=ketcher)
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
