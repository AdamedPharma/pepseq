from collections import namedtuple
from typing import TypeVar

AminoAcids = namedtuple("AminoAcids", "n_term amino_acids c_term")

DataBase = TypeVar("DataBase")
SeqReader = TypeVar("SeqReader")


def output_modified_residue(ResName: str, R_id: str) -> str:
    """
    return name/symbol of radical in the form '{ResidueName}(R{radical_id})',
    e. g. Cys(R1); Lys(R2) residue names are changed to three letter for better legibility


    :parameter ResName: one letter symbol of Residue e.g. C for Cysteine; K for Lysine: str

    :parameter R_id: sequence attachment point id: str

    :return: radical_name - string in the form '{ResidueName}(R{radical_id})', e. g. Cys(R1); Lys(R2) residue names are changed to three letter for better legibility: str

    """
    d = {"C": "Cys", "K": "Lys"}
    if ResName in d:
        ResName = d[ResName]
    radical_name = "%s(R%s)" % (ResName, R_id)
    return radical_name


def parentheses_locs_list(parentheses_locs: dict) -> list:
    """
    Get ordered list of parentheses ({}) pairs locations in seuqence text

    :param parentheses_locs: locations for parentheses pairs in the form dictionary1 = {open_loc_index: close_loc_index}
    :type  parentheses_locs: dict

    :return parentheses_pairs - list of tuples in the form [(1st_pair_open_loc_index, 1st_pair_close_loc_index) ...,]
    pairs are ordered as they appear in the sequence
    :rtype list
    """
    parentheses_pairs_indices = []
    for k, v in parentheses_locs.items():
        parentheses_pair = k, v
        parentheses_pairs_indices.append(parentheses_pair)
    parentheses_pairs_indices = sorted(parentheses_pairs_indices)
    return parentheses_pairs_indices


def find_parentheses(sequence_txt: str) -> list:
    """Find and return the location of the matching parentheses pairs in sequence text.

    Given a string, s, return a dictionary of start: end pairs giving the
    indexes of the matching parentheses in s. Suitable exceptions are
    raised if s contains unbalanced parentheses.

    :param sequence_txt: sequence text including parenthese like SDC{Aib}SSED
    :type  sequence_txt: str

    :return parentheses_pairs - list of tuples in the form [(1st_pair_open_loc_index, 1st_pair_close_loc_index) ...,]
    pairs are ordered as they appear in the sequence
    :rtype list

    """

    # The indexes of the open parentheses are stored in a stack, implemented
    # as a list

    stack = []
    parentheses_locs = {}
    for i, c in enumerate(sequence_txt):
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
    parentheses_pairs_indices = parentheses_locs_list(parentheses_locs=parentheses_locs)
    return parentheses_pairs_indices


def parse_canonical(canonical_sequence: str) -> list[str]:
    """
    Parse Sequence Text in canonical format into list of amino acid residue symbols

    :param canonical_sequence: peptide sequence in form "{%s}%s{%s}" % (n_term_symbol,
      sequence_str_wo_termini, c_term_symbol)
    :type  canonical_sequence: str

    :return: list of fragments of symbols denoting residues e.g. ['CSCAH', '{Aib}', 'CDE']
    :rtype: list[str]
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
    Parse Sequence Text in canonical format into list of amino acid residue symbols

    :param canonical_sequence: peptide sequence in form "{%s}%s{%s}" % (n_term_symbol,
      sequence_str_wo_termini, c_term_symbol)
    :type  canonical_sequence: str

    :return: list of fragments of symbols denoting residues e.g. ['CSCAH', '{Aib}', 'CDE']
    :rtype: list[str]
    
    remainder is also appended to sequence
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



def find_termini(sequence_str: str, db_json: dict) -> tuple:
    """
    :param sequence_str:
    :type  sequence_str: str

    :param db_json: database for amino acid monomers/building blocks
    :type  db_json: dict

    :return output_tuple: tuple in the form of (n_term, c_term, sequence_str_wo_termini)
    """
    sequence_split = sequence_str.split("~")
    #print(sequence_split)


    n_fragments = len(sequence_split)

    #print(n_fragments)

    if n_fragments == 1:
        return "H", "OH", sequence_str
    
    elif n_fragments == 3:

        n_term, sequence_str, c_term = sequence_split
        return n_term,  c_term, sequence_str


    elif n_fragments > 3:

        return
    
    elif n_fragments == 2:
        print(sequence_split)

        n_terms = db_json["smiles"]["n_terms"].keys()
        c_terms = db_json["smiles"]["c_terms"].keys()
        potential_n_term = sequence_split[0]
        potential_c_term = sequence_split[1]
        seq_start_index = 0
        seq_end_index = len(sequence_str)
        n_term = "H"
        c_term = "OH"

        n_wrapped_in_brackets = ('[' in potential_n_term) and (']' in potential_n_term)

        if n_wrapped_in_brackets or (potential_n_term in n_terms):
            n_term = potential_n_term
            seq_start_index = len(n_term) + 1
            sequence_str_wo_termini = sequence_str[seq_start_index:seq_end_index]
            
            output_tuple = n_term, c_term, sequence_str_wo_termini
            print(output_tuple)
            return output_tuple
        else:
        
            c_wrapped_in_brackets = ('[' in potential_c_term) and (']' in potential_c_term)
            if c_wrapped_in_brackets or (potential_c_term in c_terms):
                c_term = potential_c_term
                seq_end_index = -1 * (len(c_term) + 1)
                sequence_str_wo_termini = sequence_str[seq_start_index:seq_end_index]
                output_tuple = n_term, c_term, sequence_str_wo_termini
                return output_tuple


def get_canonical(sequence_str: str, db_json: dict) -> str:
    """
    :param sequence_str: 
    :type  sequence_str: str

    :param db_json: database with info about monomers/building blocks of modified peptides
    :type  db_json: dict

    :return canonical_sequence - seq in form "{%s}%s{%s}" % (n_term, sequence_str_wo_termini, c_term)
    :rtype:  str
    """
    n_term, c_term, sequence_str_wo_termini = find_termini(sequence_str, db_json)
    canonical_sequence = "{%s}%s{%s}" % (n_term, sequence_str_wo_termini, c_term)
    return canonical_sequence
