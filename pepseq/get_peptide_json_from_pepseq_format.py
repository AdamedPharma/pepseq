import json
import os
from typing import Dict

import rdkit

from pepseq.Peptide.utils.pure_parsing_functions import get_attachment_points_on_sequence_json, get_base_seq
from pepseq.Peptide.utils.Parser import find_termini, parse_canonical2


absolute_path = os.path.dirname(__file__)
relative_db_path = "Peptide/database/db.json"
full_db_path = os.path.join(absolute_path, relative_db_path)

with open(full_db_path) as fp:
    db_json = json.load(fp)


def get_single_modification_json(attachment_points_on_sequence: Dict, mod_smiles: str) -> Dict:
    """
    :param attachment_points_on_sequence - 
    :type  attachment_points_on_sequence: Dict

    :param mod_smiles: SMILES code for External Modification

    """
    mod_mol = rdkit.Chem.MolFromSmiles(mod_smiles)
    mod_atoms = mod_mol.GetAtoms()

    radical_dummy_atoms = [atom for atom in mod_mol.GetAtoms() if (atom.GetAtomicNum() == 0)]
    #in rdkit dummy Atoms are identified by AtomicNum = 0

    radical_ids_str = set([atom.GetIsotope() for atom in radical_dummy_atoms])
    

    min_att_points = {
        radical_id: attachment_points_on_sequence[radical_id]
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

    Args:
        symbols (list[str]): List of symbols representing attachment points on the sequence.
        smiles (list): List of SMILES strings representing the modifications.

    Returns:
        list: List of JSON representations of external modifications.
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

    Args:
        mod_smiles_list (list[str]): A list of SMILES representations of modified peptides.

    Returns:
        list: A list of JSON objects representing the modified peptides. If no modified peptides are found, an empty list is returned.
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

    Args:
        pepseq_format (str): The peptide sequence in pepseq format.
        db_json (Dict, optional): The database JSON containing the mapping of symbols to amino acids. Defaults to db_json.

    Returns:
        dict: A JSON representation of the peptide sequence.

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

def get_pep_json(pepseq_format: str, db_json: Dict = db_json, mod_smiles_list: list=None) -> Dict:
    """
    Converts a peptide sequence in pepseq format to a JSON object.
    Input:
    pepseq_string:
    str = string in pepseq format
    H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
    where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified
    amino acid; {Cys(R1)} - is amino acid
    with staple attached, {Cys(R1)} - amino acid with staple attached
    mod_smiles:
    SMILES string (e.g. '[1*]C[2*]') - showing the structure of
    modification with attachment
    points:
    { Cys(R1) } <- is attached in [1*] attachment point on staple
    { Cys(R2) } <- is attached in [2*] attachment point on staple
    Output:
    peptide_json:
    JSON containing info about modified peptide with
    'sequence':
    'internal_modifications':
    'external_modifications':

    Args:
        pepseq_format (str): The peptide sequence in pepseq format.
        db_json (Dict, optional): The database JSON object. Defaults to db_json.
        mod_smiles_list (list, optional): The list of modified SMILES strings. Defaults to None.

    Returns:
        Dict: The JSON object representing the peptide sequence.
    """

    pep_json = get_pepseq_json(pepseq_format, db_json)
    pep_json["external_modifications"] = get_smiles_json(pep_json['symbols'][1:-1], mod_smiles_list)
    return pep_json


