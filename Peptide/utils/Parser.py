from collections import namedtuple
from collections.abc import Sequence
from typing import TypeVar

from Peptide.exceptions import InvalidSymbolError, ValidationError
from Peptide.models.AminoAcidInstance import AminoAcidFromSymbol, AminoAcidInstance
from Peptide.utils.chemistry.Smi2SeqObj import Smi2Seq  # , get_adict
from Peptide.utils.validation import (
    check_for_nested_brackets,
    check_parentheses,
    validate_termini,
)

AminoAcids = namedtuple("AminoAcids", "n_term amino_acids c_term")

DataBase = TypeVar("DataBase")
SeqReader = TypeVar("SeqReader")


def parentheses_locs_list(parentheses_locs: list):
    ls = []
    for k, v in parentheses_locs.items():
        l = k, v
        ls.append(l)
    ls = sorted(ls)
    return ls


def find_parentheses(s: str):
    """Find and return the location of the matching parentheses pairs in s.

    Given a string, s, return a dictionary of start: end pairs giving the
    indexes of the matching parentheses in s. Suitable exceptions are
    raised if s contains unbalanced parentheses.

    """

    # The indexes of the open parentheses are stored in a stack, implemented
    # as a list

    stack = []
    parentheses_locs = {}
    for i, c in enumerate(s):
        if c == "{":
            stack.append(i)
        elif c == "}":
            try:
                parentheses_locs[stack.pop()] = i
            except IndexError:
                raise IndexError("Too many close parentheses at index {}".format(i))
    if stack:
        raise IndexError(
            "No matching close parenthesis to open parenthesis "
            "at index {}".format(stack.pop())
        )
    return parentheses_locs_list(parentheses_locs=parentheses_locs)


def parse_seq(
    sequence_str: str,
    db_api: DataBase,
) -> list:
    """

    takes sequence string e.g. AC{Aib}E~NH2
    parses n_terminus, c_terminus and any fragments within {} brackets
    returns list of [n_terminus_code,  aa_1code, aa2_code ...  ,c_terminus_code]

     e.g. ['H' 'A', 'Aib', E, 'NH2']

    default for n_terminus is H
    default for c_terminus is OH

    """

    n_term_symbol = db_api.find_n_term(sequence_str)
    c_term_symbol = db_api.find_c_term(sequence_str)

    if n_term_symbol is None:
        n_term_symbol = "H"
        seq_start_index = 0
    else:
        seq_start_index = len(n_term_symbol) + 1

    if c_term_symbol is None:
        c_term_symbol = "OH"
        seq_end_index = len(sequence_str)
    else:
        seq_end_index = -1 * (len(c_term_symbol) + 1)

    sequence_str_wo_termini = sequence_str[seq_start_index:seq_end_index]

    s = "{%s}%s{%s}" % (n_term_symbol, sequence_str_wo_termini, c_term_symbol)

    indices_of_brackets = find_parentheses(s)

    symbols = []

    previous_close_index = 0

    for open_index, close_index in indices_of_brackets:

        one_letter_codes_fragment = s[previous_close_index:open_index]
        symbols += list(one_letter_codes_fragment)

        fragment_in_bracket = s[(open_index + 1) : close_index]

        symbols.append(fragment_in_bracket)

        previous_close_index = close_index + 1
    return symbols


def read_sequence_txt(
    sequence: str,
    db_api: DataBase = None,
) -> Sequence[AminoAcidInstance]:

    symbols = parse_seq(sequence, db_api)

    amino_acids = []

    n_terminus_code = symbols[0]
    c_terminus_code = symbols[-1]

    c_term_smiles = db_api.c_terms_smi_codes[c_terminus_code]
    n_term_smiles = db_api.n_terms_smi_codes[n_terminus_code]

    aa_symbols_list = symbols[1:-1]

    for aa_symbol in aa_symbols_list:
        amino_acid_instance = AminoAcidFromSymbol(aa_symbol, db_api)
        amino_acids.append(amino_acid_instance)

    n_term_tuple = (n_terminus_code, n_term_smiles)
    c_term_tuple = (c_terminus_code, c_term_smiles)

    peptide_tuple = (n_term_tuple, amino_acids, c_term_tuple)
    return peptide_tuple


def validate_sequence(
    sequence: str,
    db_api: DataBase = None,
) -> bool:
    if check_parentheses(sequence):
        check_for_nested_brackets(sequence)
        validate_termini(sequence)
        symbols_list = parse_seq(sequence, db_api)
        print(symbols_list)
        invalid_positions = []

        for i in range(1, len(symbols_list) - 1):
            try:
                symbol = symbols_list[i]
                aai = AminoAcidFromSymbol(symbol, db_api=db_api)
            except InvalidSymbolError:
                invalid_positions.append(i)
        if invalid_positions:
            pos_txt = ";".join(
                [
                    "'%s' (at position: %d)" % (symbols_list[pos], pos)
                    for pos in invalid_positions
                ]
            )
            raise ValidationError("Sequence Invalid! Invalid Symbols: %s" % (pos_txt))
            # position eq i because Python indexing cancels with the fact that Nterminus is the first Symbol

        return
    else:
        raise ValidationError("Sequence Invalid! (uneven number of parentheses)")


class SmilesParser(object):
    def __init__(self):
        return

    def read_smiles_txt(
        self, smiles: str, db_api: DataBase = None
    ) -> Sequence[AminoAcidInstance]:
        smi2seq_obj = Smi2Seq(smiles)
        smi2seq_obj.renumber()
        smi2seq_obj.label_CAatoms()
        seq = smi2seq_obj.get_seq(db_api.aa_smiles_dict)
        return seq
