import json
import os

from pepseq.get_peptide_json_from_pepseq_format import get_pep_json
from pepseq.Peptide.utils.validation import (
    check_for_nested_brackets, check_parentheses,
    validate_attachment_points_on_smiles, validate_matching_attachment_points,
    validate_smiles_codes, validate_termini)
from pepseq.read import from_json

absolute_path = os.path.dirname(__file__)
relative_db_path = "Peptide/database/db.json"
full_db_path = os.path.join(absolute_path, relative_db_path)

with open(full_db_path) as fp:
    db_json = json.load(fp)


def validate_pepseq(pepseq: str):
    validate_termini(pepseq)
    check_parentheses(pepseq)
    check_for_nested_brackets(pepseq)


def calculate(pepseq: str, smiles: list[str] = []) -> dict:
    if smiles == []:
        smiles = None
    validate_pepseq(pepseq)
    validate_smiles_codes(smiles)
    validate_attachment_points_on_smiles(smiles)
    validate_matching_attachment_points(pepseq, smiles)
    peptide_json = get_pep_json(pepseq, db_json, smiles)
    peptide = from_json(peptide_json)
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
