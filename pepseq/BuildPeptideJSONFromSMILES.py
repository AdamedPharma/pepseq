import rdkit

from typing import Union
from pepseq.Backbone import (BreakingIntoResidueCandidateSubgraphs,
                             MarkingPeptideBackbone)
from pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph import \
    decompose_residues_internal, translate_external_modification
from pepseq.Peptide.utils.Parser import parse_canonical2


def get_cx_smarts_db(db_json: dict) -> dict:
    """
    Extracts CX SMARTS database from the given JSON database.

    :param db_json: The JSON database containing the SMILES information.
    :type db_json: dict

    :return: The CX SMARTS database extracted from the JSON database.
    :rtype: dict
    """
    cx_smarts_db = {}
    smiles_dict = db_json["smiles"].get("aa")

    for aa in smiles_dict:
        cx_smarts_db[aa] = smiles_dict[aa].get("smarts")

    cx_smarts_db.pop("g")

    return cx_smarts_db


def decompose_peptide_smiles(smiles: str, db_json: dict, n_subst_limit=None) -> dict:
    """
    Decomposes a peptide SMILES string into its constituent parts.
    
    :param smiles: The SMILES string representing the peptide.
    :type smiles: str

    :param db_json: The JSON database containing the SMILES/SMARTS information.
    :type db_json: dict

    :param n_subst_limit: The maximum number of substitutions allowed per residue.  Defaults to None.
    :type n_subst_limit: int or None

    :return:  A dictionary containing the decomposed peptide information, including the sequence, internal modifications, and external modifications. 
    :rtype: dict
    """

    peptide_molecule = rdkit.Chem.MolFromSmiles(smiles)
    peptide_molecule = MarkingPeptideBackbone().execute(peptide_molecule)

    residues = BreakingIntoResidueCandidateSubgraphs().execute(peptide_molecule)

    cx_smarts_db = get_cx_smarts_db(db_json)

    seq, internal_modifications, external_modifications = decompose_residues_internal(
        residues, cx_smarts_db, n_subst_limit=n_subst_limit
    )
    return {
        "sequence": seq,
        "internal_modifications": internal_modifications,
        "external_modifications": external_modifications,
    }


def get_terminal_smiles_building_block(peptide_json: dict, ResID: int, AtomName: str) -> tuple:
    """
    Get the SMILES and index of the terminal building block based on the given ResID and AtomName.

    :param peptide_json: The peptide JSON object.
    :type peptide_json: dict
    :param  ResID: The ResID of the attachment point.
    :type ResID: int
    :param AtomName: The AtomName of the attachment point.
    :type AtomName: str

    
    :return: A tuple containing the SMILES and index of the terminal building block.
    :rtype: tuple
    """

    external_modifications = peptide_json.get("external_modifications")
    if ResID == -1:
        symbols = parse_canonical2(peptide_json.get("sequence"))
        ResID = len(symbols)

    for i in range(len(external_modifications)):
        mod = external_modifications[i]
        attachment_points_on_sequence = mod.get("attachment_points_on_sequence")
        for key in attachment_points_on_sequence:
            attachment_point = attachment_points_on_sequence.get(key)
            modResID = int(attachment_point.get("ResID"))
            modAtomName = attachment_point.get("AtomName")

            if (modResID == ResID) and (modAtomName == AtomName):
                smiles = mod.get("smiles")
                return smiles, i


def get_C_terminal_smiles_building_block(peptide_json: dict) -> tuple:
    """
    Get the C-terminal SMILES building block from the given peptide JSON.

    :param peptide_json: The peptide JSON.
    :type peptide_json: dict

    :return: A tuple containing the C-terminal SMILES building block information.
    :rtype: tuple
    """

    return get_terminal_smiles_building_block(peptide_json, ResID=-1, AtomName="CO")


def get_N_terminal_smiles_building_block(peptide_json: dict) -> tuple:
    """
    Retrieves the SMILES building block for the N-terminal residue of a peptide.

    :param peptide_json: The JSON representation of the peptide.
    :type peptide_json: dict

    :return: A tuple containing the SMILES building block, the residue ID, and the atom name.
    :rtype: tuple 
    """
    return get_terminal_smiles_building_block(peptide_json, ResID=1, AtomName="N")


def smiles_are_identical(smiles1: str, smiles2: str) -> True:
    """
    Check if two SMILES strings represent identical molecules.

    :param smiles1: The first SMILES string.
    :type smiles1: str
    :param smiles2: The second SMILES string.
    :type smiles2: str

    :return: True if the molecules represented by the SMILES strings are identical, False otherwise.
    :rtype: bool 
    """

    mol1 = rdkit.Chem.MolFromSmiles(smiles1)
    mol2 = rdkit.Chem.MolFromSmiles(smiles2)
    return (mol1.HasSubstructMatch(mol2)) and (mol2.HasSubstructMatch(mol1))


def get_term_symbol(smiles: str, db_json: dict, group: str) -> Union[str, None]:
    """
    Get the term symbol for a given SMILES string and group from the database JSON.

    :param smiles: The SMILES string to match.
    :type smiles: str
    :param db_json: The database JSON containing the terms.
    :type db_json: dict
    :param group: The group to search for the term symbol.
    :type group: str

    :return: The term symbol if a match is found, None otherwise.
    :rtype: Union[str, None]
    """

    terms = db_json.get("smiles").get(group)

    for term_symbol in terms:
        term_smiles = terms.get(term_symbol).get("smiles_radical")
        if smiles_are_identical(smiles, term_smiles):
            return term_symbol


def get_c_term_from_peptide_json(peptide_json: dict, db_json: dict) -> tuple:
    """
    Retrieves the C-terminal symbol and index from a peptide JSON and a database JSON.

    :param peptide_json: The peptide JSON.
    :type peptide_json: dict

    :param db_json: The database JSON.
    :type db_json: dict

    :return: A tuple containing the C-terminal symbol and index.
    :rtype: tuple
    """

    tup = get_C_terminal_smiles_building_block(peptide_json)

    if tup is not None:
        smiles, i = tup
        term_symbol = get_term_symbol(smiles, db_json, "c_terms")
        if term_symbol is None:
            term_symbol = "SMILES:%s" % smiles
        return term_symbol, i


def get_n_term_from_peptide_json(peptide_json: dict, db_json: dict) -> tuple:
    """
    Retrieves the N-terminal symbol and index from a peptide JSON and a database JSON.

    :param peptide_json: The peptide JSON containing the information about the peptide.
    :type peptide_json: dict
    :param db_json: The database JSON containing the information about the database.
    :type db_json: dict

    :return: A tuple containing the N-terminal symbol and index.
    :rtype: tuple

    """

    tup = get_N_terminal_smiles_building_block(peptide_json)
    if tup is not None:
        smiles, i = tup
        term_symbol = get_term_symbol(smiles, db_json, "n_terms")
        if term_symbol is None:
            term_symbol = "SMILES:%s" % smiles

        return term_symbol, i


def output_modified_residue(ResName: str, R_id: str) -> str:
    """
    Returns a modified residue name with the corresponding R_id.

    :param ResName: The original residue name.
    :type ResName: str
    :param R_id: The R_id to be appended to the modified residue name.
    :type R_id: str

    :return: The modified residue name with the appended R_id.
    :rtype: str 
    """

    d = {"C": "Cys", "K": "Lys"}
    if ResName in d:
        ResName = d[ResName]
    s = "%s(R%s)" % (ResName, R_id)
    return s


def append_pepseq_R_info(j: dict) -> str:
    """
    Appends R information to the peptide sequence.

    :param j: The input dictionary containing the peptide sequence and external modifications.
    :type j: dict

    :return: The modified peptide sequence with R information appended.
    :rtype: str
    """

    seq_list = parse_canonical2(j["sequence"])
    ext_mods = j["external_modifications"]
    for ext_mod in ext_mods:
        att_points = ext_mod["attachment_points_on_sequence"]
        for R_id in att_points:
            att_point = att_points[R_id]
            ResID = int(att_point["ResID"])
            ResName = seq_list[ResID - 1]
            new_s = output_modified_residue(ResName, R_id)
            seq_list[ResID - 1] = new_s
    new_seq = ""
    for symbol in seq_list:
        if len(symbol) > 1:
            new_seq = new_seq + "{%s}" % symbol
        else:
            new_seq = new_seq + symbol
    return new_seq

    """

    """


def decompose_peptide_smiles_with_termini(smiles: str, db_json: dict, n_subst_limit=None) -> dict:
    """
    Decomposes a peptide SMILES string with termini into a JSON representation.
    Input:
    SMILES - string of peptide sequence with modified amino acids
    modification(s)
    Output:
    peptide_json:
    JSON containing info about modified peptide with
    'sequence':
    'internal_modifications':
    'external_modifications':
    mod_smiles:
    SMILES string (e.g. '[1*]C[2*]') - showing the structure of
    modification with attachment
    points:
    { Cys(R1) } <- is attached in [1*] attachment point on staple
    { Cys(R2) } <- is attached in [2*] attachment point on staple


    
    :param smiles: The peptide SMILES string.
    :type smiles: str

    :param db_json: The database JSON containing the information about the peptide.
    :type db_json: dict

    :param n_subst_limit: The maximum number of substitutions allowed. Defaults to None.
    :type n_subst_limit: int or None

    :return: The JSON representation of the decomposed peptide.
    :rtype: dict

    """

    peptide_json = decompose_peptide_smiles(smiles, db_json, n_subst_limit=n_subst_limit)
    c_terminus__c_ind = get_c_term_from_peptide_json(peptide_json, db_json)

    n_terminus__n_ind = get_n_term_from_peptide_json(peptide_json, db_json)
    external_modifications = peptide_json.get("external_modifications")

    to_del = []

    if c_terminus__c_ind is not None:
        c_terminus, c_ind = c_terminus__c_ind
        peptide_json["C_terminus"] = c_terminus
        to_del.append(c_ind)
    else:
        c_terminus = "OH"
        peptide_json["C_terminus"] = c_terminus

    if n_terminus__n_ind is not None:
        n_terminus, n_ind = n_terminus__n_ind
        peptide_json["N_terminus"] = n_terminus
        to_del.append(n_ind)
    else:
        n_terminus = "H"
        peptide_json["N_terminus"] = n_terminus

    new_ext_mods = [
        external_modifications[i]
        for i in range(len(external_modifications))
        if i not in to_del
    ]

    peptide_json["external_modifications"] = new_ext_mods

    offset = 0
    ext_mods_translated = []

    for ext_mod in peptide_json.get("external_modifications"):
        ext_mod_translated = translate_external_modification(
            ext_mod, offset=offset)
        offset += ext_mod_translated.get('max_attachment_point_id')
        ext_mods_translated.append(ext_mod_translated)

    peptide_json["external_modifications"] = ext_mods_translated

    new_seq = append_pepseq_R_info(peptide_json)

    pepseq_format = "%s~%s~%s" % (
        peptide_json["N_terminus"],
        new_seq,
        peptide_json["C_terminus"],
    )

    peptide_json["pepseq_format"] = pepseq_format

    return peptide_json


def from_smiles_to_pepseq_and_mod_smiles_strings(smiles: str, db_json: dict, n_subst_limit=None) -> tuple:
    """
    Converts a SMILES string to pepseq and mod_smiles strings.
    Input:
    SMILES - string of peptide sequence with modified amino acids
    modification(s)
    Output:
    pepseq_string:
    str = string in pepseq format
    H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
    where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified
    amino acid; {Cys(R1)} - is amino acid
    with staple attached, {Cys(R1)} - amino acid with staple attached
    modifications - external ones, with attachment points
    mod_smiles:
    SMILES string (e.g. '[1*]C[2*]') - showing the structure of
    modification with attachment
    points:
    { Cys(R1) } <- is attached in [1*] attachment point on staple
    { Cys(R2) } <- is attached in [2*] attachment point on staple
    
    :param smiles: The input SMILES string.
    :type smiles: str
    :param db_json: The database JSON containing information about the modifications.
    :type db_json: dict
    :param n_subst_limit: The maximum number of substitutions allowed. Defaults to None.
    :type n_subst_limit: int or None

    :return: A tuple containing the pepseq format string and the mod_smiles string(s).
    :rtype: tuple 
    """

    peptide_json = decompose_peptide_smiles_with_termini(smiles, db_json, n_subst_limit=n_subst_limit)
    pepseq_format = peptide_json["pepseq_format"]
    mod_smiles_list = [
        ext_mod["smiles"] for ext_mod in peptide_json["external_modifications"]
    ]
    if len(mod_smiles_list) == 1:
        return pepseq_format, mod_smiles_list[0]
    else:
        return pepseq_format, mod_smiles_list
    #return pepseq_format, mod_smiles_list

    """
    Input:

        SMILES - string of peptide sequence with modified amino acids
            modification(s)

    Output:
        pepseq_string:

            str = string in pepseq format
              H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
        where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified
          amino acid; {Cys(R1)} - is amino acid
        with staple attached, {Cys(R1)} - amino acid with staple attached

        modifications - external ones, with attachment points

        mod_smiles:

            SMILES string (e.g. '[1*]C[2*]') - showing the structure of
              modification with attachment
            points:

                { Cys(R1) } <- is attached in [1*] attachment point on staple
                { Cys(R2) } <- is attached in [2*] attachment point on staple

    """


def from_smiles_to_pepseq_and_one_mod_smiles_strings(smiles: str, db_json: dict, n_subst_limit=None):
    """
    Converts a SMILES string to a pepseq format and a list of modified SMILES strings.
    Input:
    SMILES - string of peptide sequence with modified amino acids
    modification(s)
    Output:
    pepseq_string:
    str = string in pepseq format
    H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
    where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified
    amino acid; {Cys(R1)} - is amino acid
    with staple attached, {Cys(R1)} - amino acid with staple attached
    modifications - external ones, with attachment points
    mod_smiles:
    SMILES string (e.g. '[1*]C[2*]') - showing the structure of
    modification with attachment
    points:
    { Cys(R1) } <- is attached in [1*] attachment point on staple
    { Cys(R2) } <- is attached in [2*] attachment point on staple

    Args:
        smiles (str): The input SMILES string.
        db_json (dict): The database JSON containing the mapping of modifications.
        n_subst_limit (int, optional): The maximum number of substitutions allowed. Defaults to None.

    Returns:
        tuple: A tuple containing the pepseq format and either a single modified SMILES string or a list of modified SMILES strings.
    """
    pepseq_format, mod_smiles_list = from_smiles_to_pepseq_and_mod_smiles_strings(
        smiles, db_json, n_subst_limit=n_subst_limit
    )
    if len(mod_smiles_list) == 1:
        return pepseq_format, mod_smiles_list[0]
    else:
        return pepseq_format, mod_smiles_list


def mark_external_modifications_on_seq(seq_list: list, peptide_json: dict, mod_as_X: bool = False) -> list:
    """
    Marks external modifications on a sequence list based on the provided peptide JSON.

    
    :param seq_list: The list of amino acid sequence.
    :type seq_list: list

    :param peptide_json: The peptide JSON containing information about external modifications.
    :type peptide_json: dict

    :param mod_as_X: Flag indicating whether to represent modifications as 'X'. Defaults to False.
    :type mod_as_X: bool, optional

    :return: The updated sequence list with external modifications marked.
    :rtype: list
    """

    external_modifications = peptide_json.get("external_modifications")
    for external_modification in external_modifications:
        attachment_points = external_modification.get("attachment_points_on_sequence")

        for key in attachment_points:
            attachment_point = attachment_points.get(key)
            res_id = int(attachment_point.get("ResID"))

            basic_res_name = seq_list[res_id - 1]
            if mod_as_X:
                new_res_name = "X"
            else:
                new_res_name = "{mod%s}" % (basic_res_name)

            seq_list[res_id - 1] = new_res_name
    return seq_list


def mark_internal_modifications_on_seq(seq_list: list, peptide_json: dict, mod_as_X: bool = False) -> list:
    """
    Marks internal modifications on a sequence list based on the given peptide JSON.

    :param seq_list: The list of amino acid residues in the sequence.
    :type seq_list: list

    :param peptide_json: The peptide JSON containing information about internal modifications.
    :type peptide_json: dict

    :param mod_as_X: Flag indicating whether to represent modifications as 'X'. Defaults to False.
    :type mod_as_X: bool, optional

    :return: The modified sequence list with internal modifications marked.
    :rtype: list 
    """

    internal_modifications = peptide_json.get("internal_modifications")

    for key in internal_modifications:
        mod_residues = internal_modifications[key]
        mod_residue_ids = [int(i["ResID"]) for i in mod_residues]

        for res_id in mod_residue_ids:
            basic_res_name = seq_list[res_id - 1]
            if mod_as_X:
                new_res_name = "X"
            else:
                new_res_name = "{mod%s}" % (basic_res_name)
            seq_list[res_id - 1] = new_res_name
    return seq_list


def print_sequence(peptide_json: dict, mod_as_X: bool=False) -> str:
    """
    Prints the sequence of a peptide in a formatted manner.

    :param peptide_json: A dictionary containing information about the peptide.
    :type peptide_json: dict

    :param mod_as_X: Whether to mark modifications as 'X' in the sequence. Defaults to False.
    :type mod_as_X: bool, optional

    :return: The formatted sequence of the peptide.
    :rtype: str
    """

    basic_sequence = peptide_json.get("sequence")
    aa_list = list(basic_sequence)

    N_terminus = peptide_json.get("N_terminus")
    C_terminus = peptide_json.get("C_terminus")

    aa_list = mark_external_modifications_on_seq(aa_list, peptide_json, mod_as_X)
    aa_list = mark_internal_modifications_on_seq(aa_list, peptide_json, mod_as_X)

    seq_marked = "".join(aa_list)
    if N_terminus is not None:
        seq_marked = "%s~%s" % (N_terminus, seq_marked)
    if C_terminus is not None:
        seq_marked = "%s~%s" % (seq_marked, C_terminus)

    return seq_marked
