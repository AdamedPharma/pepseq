from collections import namedtuple
from collections.abc import Sequence
from typing import Dict, TypeVar

from Peptide.exceptions import InvalidSymbolError, ValidationError
from Peptide.models.AminoAcidInstance import AminoAcidFromSymbol, AminoAcidInstance
from Peptide.utils.validation import (
    check_for_nested_brackets,
    check_parentheses,
    validate_termini,
)

AminoAcids = namedtuple("AminoAcids", "n_term amino_acids c_term")

DataBase = TypeVar("DataBase")
SeqReader = TypeVar("SeqReader")


def output_modified_residue(ResName, R_id):
    d = {"C": "Cys", "K": "Lys"}
    if ResName in d:
        ResName = d[ResName]
    s = "%s(R%s)" % (ResName, R_id)
    return s


def append_pepseq_R_info(j):
    l = parse_canonical2(j["sequence"])
    ext_mods = j["external_modifications"]
    for ext_mod in ext_mods:
        att_points = ext_mod["attachment_points_on_sequence"]
        for R_id in att_points:
            att_point = att_points[R_id]
            ResID = int(att_point["ResID"])
            ResName = l[ResID - 1]
            new_s = output_modified_residue(ResName, R_id)
            l[ResID - 1] = new_s
    new_seq = ""
    for symbol in l:
        if len(symbol) > 1:
            new_seq = new_seq + "{%s}" % symbol
        else:
            new_seq = new_seq + symbol
    return new_seq


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


def parse_canonical(canonical_sequence):
    """
    canonical_sequence = "{%s}%s{%s}" % (n_term_symbol, sequence_str_wo_termini, c_term_symbol)
    """
    indices_of_brackets = find_parentheses(canonical_sequence)
    symbols = []
    previous_close_index = 0
    for open_index, close_index in indices_of_brackets:

        one_letter_codes_fragment = canonical_sequence[previous_close_index:open_index]
        symbols += list(one_letter_codes_fragment)

        fragment_in_bracket = canonical_sequence[(open_index + 1) : close_index]

        symbols.append(fragment_in_bracket)

        previous_close_index = close_index + 1
    return symbols


def parse_canonical2(canonical_sequence):
    """
    canonical_sequence = "{%s}%s{%s}" % (n_term_symbol, sequence_str_wo_termini, c_term_symbol)
    """
    indices_of_brackets = find_parentheses(canonical_sequence)
    symbols = []
    previous_close_index = 0
    for open_index, close_index in indices_of_brackets:

        one_letter_codes_fragment = canonical_sequence[previous_close_index:open_index]
        symbols += list(one_letter_codes_fragment)

        fragment_in_bracket = canonical_sequence[(open_index + 1) : close_index]

        symbols.append(fragment_in_bracket)

        previous_close_index = close_index + 1
    symbols += canonical_sequence[previous_close_index:]
    return symbols


def find_termini(sequence_str, db_json):
    sequence_split = sequence_str.split("~")
    if len(sequence_split) == 1:
        return "H", "OH", sequence_str

    else:
        n_terms = db_json["smiles"]["n_terms"].keys()
        if sequence_split[0] in n_terms:
            n_term = sequence_split[0]
            seq_start_index = len(n_term) + 1
        else:
            n_term = "H"
            seq_start_index = 0

        c_terms = db_json["smiles"]["c_terms"].keys()

        if sequence_split[-1] in c_terms:
            c_term = sequence_split[-1]
            seq_end_index = -1 * (len(c_term) + 1)

        else:
            c_term = "OH"
            seq_end_index = len(sequence_str)
        sequence_str_wo_termini = sequence_str[seq_start_index:seq_end_index]
        return n_term, c_term, sequence_str_wo_termini


def get_canonical(sequence_str, db_json):
    n_term, c_term, sequence_str_wo_termini = find_termini(sequence_str, db_json)
    canonical_sequence = "{%s}%s{%s}" % (n_term, sequence_str_wo_termini, c_term)
    return canonical_sequence


def parse_seq(
    sequence_str: str,
    db_json: Dict,
) -> list:
    """

    takes sequence string e.g. AC{Aib}E~NH2
    parses n_terminus, c_terminus and any fragments within {} brackets
    returns list of [n_terminus_code,  aa_1code, aa2_code ...  ,c_terminus_code]

     e.g. ['H' 'A', 'Aib', E, 'NH2']

    default for n_terminus is H
    default for c_terminus is OH

    """
    canonical_sequence = get_canonical(sequence_str, db_json)
    symbols = parse_canonical(canonical_sequence)
    return symbols


def read_sequence_txt(
    sequence: str,
    db_json: Dict = None,
) -> Sequence[AminoAcidInstance]:

    symbols = parse_seq(sequence, db_json)

    amino_acids = []

    n_terminus_code = symbols[0]
    c_terminus_code = symbols[-1]

    c_term_smiles = db_json["smiles"]["c_terms"][c_terminus_code]["smiles"]
    n_term_smiles = db_json["smiles"]["n_terms"][n_terminus_code]["smiles"]

    aa_symbols_list = symbols[1:-1]

    for aa_symbol in aa_symbols_list:
        amino_acid_instance = AminoAcidFromSymbol(aa_symbol, db_json)
        amino_acids.append(amino_acid_instance)

    n_term_tuple = (n_terminus_code, n_term_smiles)
    c_term_tuple = (c_terminus_code, c_term_smiles)

    peptide_tuple = (n_term_tuple, amino_acids, c_term_tuple)
    return peptide_tuple


def validate_sequence(
    sequence: str,
    db_json: Dict = None,
) -> bool:
    if check_parentheses(sequence):
        check_for_nested_brackets(sequence)
        validate_termini(sequence)
        symbols_list = parse_seq(sequence, db_api)
        invalid_positions = []

        for i in range(1, len(symbols_list) - 1):
            try:
                symbol = symbols_list[i]
                aai = AminoAcidFromSymbol(symbol, db_json)
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
