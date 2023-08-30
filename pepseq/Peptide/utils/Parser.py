from collections import namedtuple
from typing import TypeVar

AminoAcids = namedtuple("AminoAcids", "n_term amino_acids c_term")

DataBase = TypeVar("DataBase")
SeqReader = TypeVar("SeqReader")


def output_modified_residue(ResName: str, R_id: str):
    d = {"C": "Cys", "K": "Lys"}
    if ResName in d:
        ResName = d[ResName]
    s = "%s(R%s)" % (ResName, R_id)
    return s


def parentheses_locs_list(parentheses_locs: list):
    ls = []
    for k, v in parentheses_locs.items():
        parentheses_pair = k, v
        ls.append(parentheses_pair)
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

def parse_canonical(canonical_sequence: str):
    """
    canonical_sequence = "{%s}%s{%s}" % (n_term_symbol,
      sequence_str_wo_termini, c_term_symbol)
    """
    indices_of_brackets = find_parentheses(canonical_sequence)
    symbols = []
    previous_close_index = 0
    for open_index, close_index in indices_of_brackets:

        one_letter_codes_fragment = canonical_sequence[previous_close_index:open_index]
        symbols += list(one_letter_codes_fragment)

        fragment_in_bracket = canonical_sequence[(open_index + 1): close_index]

        symbols.append(fragment_in_bracket)

        previous_close_index = close_index + 1
    return symbols


def parse_canonical2(canonical_sequence: str) -> list:
    """
    canonical_sequence = "{%s}%s{%s}" % (n_term_symbol,
      sequence_str_wo_termini, c_term_symbol)
    """
    indices_of_brackets = find_parentheses(canonical_sequence)
    symbols = []
    previous_close_index = 0
    for open_index, close_index in indices_of_brackets:

        one_letter_codes_fragment = canonical_sequence[previous_close_index:open_index]
        symbols += list(one_letter_codes_fragment)

        fragment_in_bracket = canonical_sequence[(open_index + 1): close_index]

        symbols.append(fragment_in_bracket)

        previous_close_index = close_index + 1
    symbols += canonical_sequence[previous_close_index:]
    return symbols


def find_termini(sequence_str: str, db_json) -> tuple:
    sequence_split = sequence_str.split("~")
    if len(sequence_split) == 1:
        return "H", "OH", sequence_str

    else:
        n_terms = db_json["smiles"]["n_terms"].keys()
        potential_n_term = sequence_split[0]

        if ('[' in potential_n_term) and (']' in potential_n_term):
            seq_start_index = len(potential_n_term) + 1
            n_term = potential_n_term

        elif sequence_split[0] in n_terms:
            n_term = potential_n_term
            seq_start_index = len(n_term) + 1
        else:
            n_term = "H"
            seq_start_index = 0

        c_terms = db_json["smiles"]["c_terms"].keys()
        potential_c_term = sequence_split[-1]

        if ('[' in potential_c_term) and (']' in potential_c_term):
            c_term = potential_c_term
            seq_end_index = -1 * (len(c_term) + 1)

        elif potential_c_term in c_terms:
            c_term = potential_c_term
            seq_end_index = -1 * (len(c_term) + 1)

        else:
            c_term = "OH"
            seq_end_index = len(sequence_str)
        sequence_str_wo_termini = sequence_str[seq_start_index:seq_end_index]
        return n_term, c_term, sequence_str_wo_termini


def get_canonical(sequence_str: str, db_json) -> str:
    n_term, c_term, sequence_str_wo_termini = find_termini(sequence_str, db_json)
    canonical_sequence = "{%s}%s{%s}" % (n_term, sequence_str_wo_termini, c_term)
    return canonical_sequence
