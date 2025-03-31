"""
This module aims to contains only pure parsing functions
This means it shall have no dependencies on other modules of the library

"""
from typing import Union, Dict
from pepseq.Peptide.exceptions import AttachmentPointsNonUniqueError


def get_attachment_point_json(
    res_id: int,
    decomposition: tuple,
    default_exit_atom_name: dict = {
        "Cys": "SG",
        "Lys": "NZ",
        "Ala": "CB",
        "Gly": "CA",
        "ala": "CB",
    },
) -> Dict:
    """
    :param res_id: Residue ID
    :type  res_id:  int

    :param decomposition: in the form of ResidueName ('Cys', 'Lys', etc.) and connecting point_id (e.g. ('Cys', 1) )
    :type  decomposition: tuple

    :param default_exit_atom_name: dictionary with default exit atom names for each residue
    :type  default_exit_atom_name: dict

    :return: att_point_json - JSON dictionary having info on Atom (Name)
      and Residue (Name and Index) to which modification is connected
    :rtype: dict
    """

    ResName, attachment_point_id = decomposition

    AtomName = default_exit_atom_name.get(ResName, "")

    att_point_json = {
        "attachment_point_id": attachment_point_id,
        "ResID": str(res_id),
        "AtomName": AtomName,
        "ResidueName": ResName,
    }
    return att_point_json


def decompose_symbol(symbol: str) -> Union[tuple, str]:
    """
    If symbol is composed e.g. Cys(R1), decompose it into (Cys, '1')

    :param s: symbol can be in normal form e.g. 'C', 'S', 'Aib', or can contain Radical with ID e.g. Cys(R1)
    :type  s: str

    :return: decomp - decomposition can be a tuple of residue name and radical id or can be an input symbol
    :rtype: tuple or str
    """
    has_round_bracket = ("(" in symbol) and (")" in symbol)

    if has_round_bracket:
        res_name, split1_2 = symbol.split("(")
        inside_bracket = split1_2.split(")")[0]

        if "R" in inside_bracket:  # e.g. Cys(R1)
            radical_id = inside_bracket[1:]
            return res_name, radical_id
    return symbol


def get_attachment_points_on_sequence_json(symbols: list) -> Dict:
    """
    Iterate over symbols of residues in amino acid sequence.
    For each symbol see if it is of form e.g. Cys(R1)
    If so: decompose it to get Cys, 1
    If not: return symbol

    :param symbols
    :type symbols: list

    :return: att_points - dictionary of attachment points on sequence
    :rtype: dict
    """
    att_points = {}

    decomposition_tuples = []

    for symbol_id in range(len(symbols)):
        symbol = symbols[symbol_id]
        decomposition = decompose_symbol(symbol)
        if type(decomposition) == tuple:
            res_name, attachment_point_id = decomposition
            res_id = symbol_id + 1


            decomposition_tuples.append((res_id, res_name, attachment_point_id))

    for res_id, res_name, attachment_point_id in decomposition_tuples:
        attachment_point_json = get_attachment_point_json(
            res_id, (res_name, attachment_point_id)
        )

        att_point_id = int(attachment_point_id)
        if att_points.get(att_point_id) is not None:
            raise AttachmentPointsNonUniqueError(
                "Attachment Points labels on sequence are not unique."
            )
        att_points[att_point_id] = attachment_point_json
    return att_points


def get_base_symbols(
    symbols: list[str],
    three_to_one: dict = {"Cys": "C", "Lys": "K", "Ala": "A", "ala": "a", "Gly": "G"},
):
    """
    Get the base symbols from a list of symbols.

    :param symbols: List of symbols.
    :type symbols: list[str]

    :param three_to_one: Dictionary mapping three-letter residue names to one-letter residue names.
                            Defaults to {"Cys": "C", "Lys": "K", "Ala": "A", "ala": "a", "Gly": "G"}.
    :type three_to_one: dict, optional

    :return: The base symbols.
    :rtype: list
    """
    base_symbols = []

    for res_id in range(len(symbols)):
        symbol = symbols[res_id]
        decomposition = decompose_symbol(symbol)
        if type(decomposition) == tuple:
            res_name, attachment_point_id = decomposition
            if res_name in three_to_one:
                res_name = three_to_one[res_name]
            symbol = res_name
        base_symbols.append(symbol)
    return base_symbols


def get_base_seq(
    symbols: list[str],
    three_to_one: dict = {"Cys": "C", "Lys": "K", "Ala": "A", "ala": "a", "Gly": "G"},
) -> str:
    """
    Get the base sequence from a list of residue symbols.

    :param symbols: List of residue symbols.
    :type symbols: list[str]

    :param three_to_one: Dictionary mapping three-letter residue names to one-letter residue names.
                         Defaults to {"Cys": "C", "Lys": "K", "Ala": "A", "ala": "a", "Gly": "G"}.
    :type three_to_one: dict, optional

    :return: The base sequence in the form of a string.
    :rtype: str
    """

    base_symbols = get_base_symbols(symbols, three_to_one=three_to_one)
    base_seq = ""

    for base_symbol in base_symbols:
        if len(base_symbol) > 1:
            base_symbol = "{%s}" % (base_symbol)
        base_seq = base_seq + base_symbol

    return base_seq
