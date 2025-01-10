import json
import os
from typing import Dict

import rdkit

from pepseq.Peptide.utils.pure_parsing_functions import (
    get_attachment_points_on_sequence_json,
    get_base_seq,
)
from pepseq.Peptide.utils.Parser import find_termini, parse_canonical2


absolute_path = os.path.dirname(__file__)
relative_db_path = "Peptide/database/db.json"
full_db_path = os.path.join(absolute_path, relative_db_path)

with open(full_db_path) as fp:
    db_json = json.load(fp)


def get_single_modification_json(
    attachment_points_on_sequence: Dict, mod_smiles: str
) -> Dict:
    """
    :param attachment_points_on_sequence: Dictionary containing attachment points on the sequence
    :type  attachment_points_on_sequence: Dict

    :param mod_smiles: SMILES code for External Modification
    :type  mod_smiles: str

    :return: JSON representation of External Modification
    :rtype: Dict
    """

    mod_mol = rdkit.Chem.MolFromSmiles(mod_smiles)
    radical_ids_str = set([])

    for atom in mod_mol.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            continue
        else:
            radical_ids_str.add(atom.GetAtomMapNum())

    min_att_points = {
        str(radical_id): attachment_points_on_sequence[int(radical_id)]
        for radical_id in radical_ids_str
    }

    max_attachment_point_id = max([int(i) for i in min_att_points.keys()])

    ext_mod = {
        "smiles": mod_smiles,
        "max_attachment_point_id": max_attachment_point_id,
        "attachment_points_on_sequence": min_att_points,
    }
    return ext_mod


def get_ext_mod_json(symbols: list[str], smiles: list) -> list:
    """
    Get the JSON representation of external modifications based on symbols and SMILES.

    :param symbols: List of symbols representing attachment points on the sequence
    :type  symbols: list[str]

    :param smiles: List of SMILES strings representing the modifications
    :type  smiles: list

    :return: List of JSON representations of external modifications
    :rtype: list
    """

    attachment_points_on_sequence = get_attachment_points_on_sequence_json(symbols)
    ext_mod_jsons = []

    if attachment_points_on_sequence.keys():
        for mod_smiles in smiles:
            ext_mod_json = get_single_modification_json(
                attachment_points_on_sequence, mod_smiles
            )
            ext_mod_jsons.append(ext_mod_json)
    return ext_mod_jsons


def get_smiles_json(symbols: str, mod_smiles_list: list[str]):
    """
    Retrieves the JSON representation of modified peptides based on their SMILES representation.

    :param mod_smiles_list: List of SMILES representations of modified peptides
    :type  mod_smiles_list: list[str]

    :param symbols: List of symbols representing attachment points on the sequence
    :type  symbols: str

    :return: List of JSON representations of modified peptides
          If no modified peptides are found, an empty list is returned.

    :rtype: list
    """
    if mod_smiles_list is not None:
        ext_mod = get_ext_mod_json(symbols, mod_smiles_list)
        if ext_mod is not None:
            return ext_mod
        else:
            return []
    else:
        return []


def get_pepseq_json(pepseq_format: str, db_json: Dict = db_json):
    """
    Convert a peptide sequence in pepseq format to a JSON representation.

    :param pepseq_format: The peptide sequence in pepseq format.
    :type  pepseq_format: str

    :param db_json: The database JSON containing the mapping of symbols to amino acids. Optional.
        Defaults to db_json.
    :type  db_json: Dict

    :return: A JSON representation of the peptide sequence.
    :rtype: Dict
    """
    N_terminus, C_terminus, pepseq = find_termini(pepseq_format, db_json)
    symbols = parse_canonical2(pepseq)
    base_seq = get_base_seq(symbols)

    all_symbols = [N_terminus] + symbols + [C_terminus]

    pep_json = {
        "length": len(symbols),
        "sequence": base_seq,
        "internal_modifications": [],
        "C_terminus": C_terminus,
        "N_terminus": N_terminus,
        "pepseq_format": pepseq_format,
        "symbols": all_symbols,
    }

    return pep_json


def get_pep_json(
    pepseq_format: str,
    db_json: Dict = db_json,
    mod_smiles_list: list = None,
) -> Dict:
    """
    Get pep_json
        peptide_json, a JSON containing info about modified peptide with
        'sequence', 'internal_modifications', 'external_modifications'
        pepseq format example is H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH

    Input:

        pepseq_string:

            str = string in pepseq format
            H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
        where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified
        amino acid; {Cys(R1)} - is amino acid
        with staple attached, {Cys(R1)} - amino acid with staple attached
        mod_smiles example is SMILES string (e.g. '[1*]C[2*]') - showing the structure of
        modification with attachment

    :param pepseq_format: string in pepseq format
    :type pepseq_format: str

    :param db_json: database JSON containing the mapping of symbols to amino acids
    :type db_json: Dict

    :param mod_smiles_list: list of SMILES strings representing the modifications
    :type mod_smiles_list: list

    :return: peptide_json: JSON representation of the peptide sequence
    :rtype: Dict
    """

    pep_json = get_pepseq_json(pepseq_format, db_json)
    pep_json["external_modifications"] = get_smiles_json(
        pep_json["symbols"][1:-1], mod_smiles_list
    )
    return pep_json
