import os
import json

from pathlib import Path
from pepseq.Peptide.exceptions import (
    ExcessTildeError,
    NestedBracketError,
    ParenthesesError,
    InvalidSymbolError,
    ValidationError,
)
from pepseq.Peptide.utils.Parser import find_termini, parse_canonical2
from pepseq.Peptide.utils.pure_parsing_functions import get_base_symbols

from pepseq.Peptide.database.db_functions import get_coding


absolute_path = Path(__file__).parent.parent.parent
relative_db_path = "Peptide/database/db.json"
full_db_path = os.path.join(absolute_path, relative_db_path)

with open(full_db_path) as fp:
    db_json = json.load(fp)


def validate_termini(s: str) -> bool:
    """
    Validates the termini of a peptide sequence.

    Args:
        s (str): The peptide sequence to be validated.

    Returns:
        bool: True if the termini are valid, False otherwise.

    Raises:
        ExcessTildeError: If the number of tildes in the sequence is greater than 2.
    """
    tilde_num = s.count("~")
    if tilde_num in [0, 1, 2]:
        return True
    elif tilde_num > 2:
        raise ExcessTildeError


def check_parentheses(s):
    """Return True if the parentheses in string s match, otherwise raise ParenthesesError.

    Args:
        s (str): The string to check for matching parentheses.

    Raises:
        ParenthesesError: If the parentheses in the string do not match.

    Returns:
        bool: True if the parentheses match, False otherwise.
    """
    j = 0
    for c in s:
        if c == "}":
            j -= 1
            if j < 0:
                raise ParenthesesError("Brackets do not match")
        elif c == "{":
            j += 1
    if j != 0:
        raise ParenthesesError("Brackets do not match")
    return True


def check_for_nested_brackets(s):
    """
    Check if the given string has nested brackets.

    Args:
        s (str): The string to be checked.

    Raises:
        NestedBracketError: If nested '{','}' brackets are found.
        ValidationError: If misplaced '}' brackets are found.
    """
    open_bracket = False

    for c in s:
        if c == "{":
            if open_bracket:
                raise NestedBracketError("Found Nested '{','}' brackets.")
            else:
                open_bracket = True
        elif c == "}":
            if open_bracket:
                open_bracket = False
            else:
                raise ValidationError("Misplaced '}' Brackets")


def get_all_available_symbols(db: dict):
    """
    Get all available symbols in the database.

    Parameters:
    db (dict): The database containing the symbols.

    Returns:
    set: A set of all available symbols.
    """
    aa_smiles_dict = db_json["smiles"].get("aa")
    coding = get_coding(db_json)
    unique_encoded_symbols = set(coding.keys())
    unique_aa = set(aa_smiles_dict.keys())
    return unique_encoded_symbols | unique_aa


def validate_monomers_in_database(pepseq_format: str, db: dict):
    """
    Validate if all of the monomers extracted from the peptide sequence are present in the database.

    Args:
        pepseq (str): The peptide sequence to validate.
        db (dict): The database containing the monomer information.

    Raises:
        InvalidSymbolError: If any of the monomers are not found in the database.

    """
    # we need to extract sequence first
    print(pepseq_format, find_termini(pepseq_format, db))
    N_terminus, C_terminus, pepseq = find_termini(pepseq_format, db)
    residue_symbols = parse_canonical2(pepseq)
    base_symbols = get_base_symbols(
        residue_symbols,
        three_to_one={"Cys": "C", "Lys": "K", "Ala": "A", "ala": "a", "Gly": "G"},
    )

    unique_residue_symbols = set(base_symbols)

    # sequence_object = Sequence(sequence, N_terminus, C_terminus)
    # residue_symbols, N_terminus_smiles, C_terminus_smiles = sequence_object.extract_residue_symbols(db_json)
    symbols_in_db = get_all_available_symbols(db)

    db_symbols_404 = unique_residue_symbols - symbols_in_db

    if db_symbols_404:
        raise InvalidSymbolError(
            "Residue Symbols: %s not found in database."
            % ", ".join(list(db_symbols_404))
        )
    return


def validate_pepseq(pepseq: str, db: dict = db_json):
    """
    Validates a peptide sequence.

    Args:
        pepseq (str): The peptide sequence to validate. format
        db (dict): The database of valid monomers (default: db_json).

    Returns:
        None
    """
    validate_termini(pepseq)
    check_parentheses(pepseq)
    check_for_nested_brackets(pepseq)
    validate_monomers_in_database(pepseq, db)
