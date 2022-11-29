from collections import namedtuple
from collections.abc import Sequence
from typing import TypeVar

from Peptide.exceptions import InvalidSymbolError, ValidationError
from Peptide.models.AminoAcidInstance import AminoAcidInstance
from Peptide.utils.chemistry.Smi2SeqObj import Smi2Seq  # , get_adict

# from Peptide.utils.chemistry.SubstructureGraph import get_ordered_nodes
from Peptide.utils.RepresentationFormat import (
    RepresentationFormat,
)  # AlignmentRepresentation,
from Peptide.utils.validation import (
    check_for_nested_brackets,
    check_parentheses,
    validate_termini,
)

Terminus = namedtuple("Terminus", "name smiles")
AminoAcids = namedtuple("AminoAcids", "n_term amino_acids c_term")

DataBase = TypeVar("DataBase")
SeqReader = TypeVar("SeqReader")


def parentheses_locs_list(parentheses_locs=None):
    ls = []
    for k, v in parentheses_locs.items():
        l = k, v
        ls.append(l)
    ls = sorted(ls)
    return ls


def find_parentheses(s):
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


def parse_seq(s, db_api: DataBase = None, representation: RepresentationFormat = None):
    if representation is None:
        representation = RepresentationFormat(name="canonical")

    n_term_symbol = "{H}"
    c_term_symbol = "{OH}"

    for smi_code in db_api.n_terms_smi_codes.keys():
        smi_code_txt = "%s~" % (smi_code)
        if s.startswith(smi_code_txt):
            n_term_symbol = "{%s}" % smi_code
            s = s[len(smi_code_txt) :]

    for smi_code in db_api.c_terms_smi_codes.keys():
        smi_code_txt = "~%s" % (smi_code)
        if s.endswith(smi_code_txt):
            c_term_symbol = "{%s}" % smi_code
            s = s[: -len(smi_code_txt)]

    s = n_term_symbol + s + c_term_symbol

    parentheses_locs = find_parentheses(s)

    symbols = []

    last_end = 0

    for bracket in parentheses_locs:

        start, end = bracket

        non_mod = s[last_end:start]
        symbols += list(non_mod)

        pos_symbol = s[(start + 1) : end]

        symbols.append(pos_symbol)

        last_end = end + 1

    non_mod = s[last_end : len(s)]
    symbols += list(non_mod)
    return symbols


"""
def reflect_dict():
    return
"""


class ParsingError(Exception):
    pass


class Parser(object):
    def __init__(self, representation: RepresentationFormat = None):
        self.representation = representation
        return

    def read_sequence_txt(
        self,
        sequence: str,
        representation: RepresentationFormat = None,
        db_api: DataBase = None,
    ) -> Sequence[AminoAcidInstance]:

        symbols_list = self.split_into_symbols_list(sequence, db_api=db_api)

        amino_acids = []

        n_term = "H"
        c_term = "OH"

        if symbols_list[-1] in db_api.c_terms_smi_codes.keys():
            c_term_smiles = db_api.c_terms_smi_codes[symbols_list[-1]]

            c_term = symbols_list[-1]
            symbols_list = symbols_list[:-1]

        if symbols_list[0] in db_api.n_terms_smi_codes.keys():
            n_term_smiles = db_api.n_terms_smi_codes[symbols_list[0]]

            n_term = symbols_list[0]
            symbols_list = symbols_list[1:]

        for symbol in symbols_list:
            aai = AminoAcidInstance.MolFromSymbol(symbol, db_api=db_api)
            amino_acids.append(aai)

        n_term_nt = Terminus(n_term, n_term_smiles)
        c_term_nt = Terminus(c_term, c_term_smiles)

        aai_namedtuple = AminoAcids(n_term_nt, amino_acids, c_term_nt)

        return aai_namedtuple

    def guess_representation(self, sequence: str) -> RepresentationFormat:
        return RepresentationFormat(name="canonical")

    def split_into_symbols_list(
        self,
        sequence: str,
        representation: RepresentationFormat = None,
        db_api: DataBase = None,
    ) -> Sequence[str]:
        return parse_seq(sequence, db_api, representation)

    def validate_sequence(
        self,
        sequence: str,
        reader: SeqReader = None,
        db_api: DataBase = None,
    ) -> bool:
        if check_parentheses(sequence):
            check_for_nested_brackets(sequence)
            validate_termini(sequence)
            symbols_list = self.split_into_symbols_list(sequence, db_api=db_api)
            invalid_positions = []

            for i in range(1, len(symbols_list) - 1):
                try:
                    symbol = symbols_list[i]
                    aai = AminoAcidInstance.MolFromSymbol(symbol, db_api=db_api)
                except InvalidSymbolError:
                    invalid_positions.append(i)
            if invalid_positions:
                pos_txt = ";".join(
                    [
                        "'%s' (at position: %d)" % (symbols_list[pos], pos)
                        for pos in invalid_positions
                    ]
                )
                raise ValidationError(
                    "Sequence Invalid! Invalid Symbols: %s" % (pos_txt)
                )
                # position eq i because Python indexing cancels with the fact that Nterminus is the first Symbol

            return
        else:
            raise ValidationError("Sequence Invalid! (uneven number of parentheses)")

        return


# take good care of exception handling
def parse_symbol(symbol):

    base_aa_symbol, substitute_radical, attachment_atom_label = parse_symbol(symbol)


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
