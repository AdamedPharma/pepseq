import rdkit

from typing import Union
from pepseq.Backbone import (
    BreakingIntoResidueCandidateSubgraphs,
    MarkingPeptideBackbone,
)
from pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph import (
    decompose_residues_internal,
    translate_external_modification,
)
from pepseq.Peptide.utils.Parser import parse_canonical2


def get_cx_smarts_db(db_json: dict) -> dict:
    """
    Get dictionary of CX_SMARTS_DB from db_json

    :param db_json: JSON containing the mapping of symbols to amino acids
    :type db_json: dict

    :return: cx_smarts_db: dictionary of CX_SMARTS_DB
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
    Decompose peptide molecule given in SMILES into peptide_json representation with modification
    SMILES codes

    :param smiles: SMILES code representing modified peptide molecule
    :type smiles: str

    :param db_json: JSON containing the mapping of symbols to amino acids
    :type db_json: dict

    :param n_subst_limit: int = number of substitutions allowed
    :type n_subst_limit: int

    :return: peptide_json: JSON containing info about modified peptide with 'sequence',
    'internal_modifications', 'external_modifications':
    mod_smiles: SMILES string (e.g. '[*:1]C[*:2]') - showing the structure of

    modification with attachment points:

    { Cys(R1) } <- is attached in [*:1] attachment point on staple

    { Cys(R2) } <- is attached in [*:2] attachment point on staple

    :rtype peptide_json: dict
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


def get_terminal_smiles_building_block(
    peptide_json: dict, ResID: int, AtomName: str
) -> tuple:
    """
    Get smiles of terminal building block of peptide

    :param peptide_json: JSON containing info about modified peptide with 'sequence',
    'internal_modifications', 'external_modifications':
    mod_smiles: SMILES string (e.g. '[*:1]C[*:2]') - showing the structure of
    modification with attachment points:
    { Cys(R1) } <- is attached in [*:1] attachment point on staple
    { Cys(R2) } <- is attached in [*:2] attachment point on staple
    :type peptide_json: dict

    :param ResID: residue ID
    :type ResID: int

    :param AtomName: atom name
    :type AtomName: str

    :return: smiles, i: tuple of smiles and index
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
    Get C terminal smiles building block

    :param peptide_json: JSON containing info about modified peptide with 'sequence',
    'internal_modifications', 'external_modifications':
    mod_smiles: SMILES string (e.g. '[*:1]C[*:2]') - showing the structure of
    modification with attachment points:
    { Cys(R1) } <- is attached in [*:1] attachment point on staple
    { Cys(R2) } <- is attached in [*:2] attachment point on staple
    :type peptide_json: dict

    :return: tuple of smiles and index
    :rtype: tuple
    """
    return get_terminal_smiles_building_block(peptide_json, ResID=-1, AtomName="CO")


def get_N_terminal_smiles_building_block(peptide_json: dict) -> tuple:
    """
    Get N terminal smiles building block

    :param peptide_json: JSON containing info about modified peptide with 'sequence',
    'internal_modifications', 'external_modifications':
    mod_smiles: SMILES string (e.g. '[*:1]C[*:2]') - showing the structure of
    modification with attachment points:
    { Cys(R1) } <- is attached in [*:1] attachment point on staple
    { Cys(R2) } <- is attached in [*:2] attachment point on staple
    :type peptide_json: dict

    :return: tuple of smiles and index
    :rtype: tuple
    """
    return get_terminal_smiles_building_block(peptide_json, ResID=1, AtomName="N")


def smiles_are_identical(smiles1: str, smiles2: str) -> True:
    """
    Check if two SMILES are identical

    :param smiles1: SMILES code
    :type smiles1: str

    :param smiles2: SMILES code
    :type smiles2: str

    :return: True if SMILES are identical
    :rtype: True
    """
    mol1 = rdkit.Chem.MolFromSmiles(smiles1)
    mol2 = rdkit.Chem.MolFromSmiles(smiles2)
    return (mol1.HasSubstructMatch(mol2)) and (mol2.HasSubstructMatch(mol1))


def get_term_symbol(smiles: str, db_json: dict, group: str) -> Union[str, None]:
    """
    Get terminal symbol

    :param smiles: SMILES code
    :type smiles: str

    :param db_json: JSON containing the mapping of symbols to amino acids
    :type db_json: dict

    :param group: name of terminal group
    :type group: str

    :return: terminal symbol
    :rtype: Union[str, None]
    """
    terms = db_json.get("smiles").get(group)

    for term_symbol in terms:
        term_smiles = terms.get(term_symbol).get("smiles_radical")
        if smiles_are_identical(smiles, term_smiles):
            return term_symbol


def get_c_term_from_peptide_json(peptide_json: dict, db_json: dict) -> tuple:
    """
    Get C terminal from peptide JSON

    :param peptide_json: JSON containing info about modified peptide with 'sequence',
    'internal_modifications', 'external_modifications':
    mod_smiles: SMILES string (e.g. '[*:1]C[*:2]') - showing the structure of
    modification with attachment points:
    { Cys(R1) } <- is attached in [*:1] attachment point on staple
    { Cys(R2) } <- is attached in [*:2] attachment point on staple
    :type peptide_json: dict

    :param db_json: JSON containing the mapping of symbols to amino acids
    :type db_json: dict

    :return: tuple of terminal symbol and index
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
    Get N terminus from peptide JSON

    :param peptide_json: JSON containing info about modified peptide with 'sequence',
    'internal_modifications', 'external_modifications':
    mod_smiles: SMILES string (e.g. '[*:1]C[*:2]') - showing the structure of

    modification with attachment points:
    { Cys(R1) } <- is attached in [*:1] attachment point on staple
    { Cys(R2) } <- is attached in [*:2] attachment point on staple
    :type peptide_json: dict

    :param db_json: JSON containing the mapping of symbols to amino acids
    :type db_json: dict

    :return: tuple of terminal symbol and index
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
    Output modified residue

    :param ResName: residue name
    :type ResName: str

    :param R_id: residue ID
    :type R_id: str

    :return: modified residue
    :rtype: str
    """
    d = {"C": "Cys", "K": "Lys"}
    if ResName in d:
        ResName = d[ResName]
    s = "%s(R%s)" % (ResName, R_id)
    return s


def append_pepseq_R_info(j: dict) -> str:
    """
    Append pepseq R info

    :param j: JSON containing info about modified peptide with 'sequence',
    'internal_modifications', 'external_modifications':
    mod_smiles: SMILES string (e.g. '[*:1]C[*:2]') - showing the structure of
    modification with attachment points:
    { Cys(R1) } <- is attached in [*:1] attachment point on staple
    { Cys(R2) } <- is attached in [*:2] attachment point on staple
    :type j: dict

    :return: new_seq: new sequence
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


def decompose_peptide_smiles_with_termini(
    smiles: str, db_json: dict, n_subst_limit=None
) -> dict:
    """
    Decompose peptide molecule given in SMILES into peptide_json representation with modification
    SMILES codes
    peptide_json: JSON containing info about modified peptide with 'sequence', 'internal_modifications',
    'external_modifications':
    mod_smiles: SMILES string (e.g. '[*:1]C[*:2]') - showing the structure of
    modification with attachment points:
    { Cys(R1) } <- is attached in [*:1] attachment point on staple
    { Cys(R2) } <- is attached in [*:2] attachment point on staple

    :param smiles: SMILES code representing modified peptide molecule
    :type smiles: str

    :param db_json: JSON containing the mapping of symbols to amino acids
    :type db_json: dict

    :param n_subst_limit: number of substitutions allowed
    :type n_subst_limit: int

    :return: peptide_json, mod_smiles: JSON containing info about modified peptide with termini
    :rtype: dict
    """

    peptide_json = decompose_peptide_smiles(
        smiles, db_json, n_subst_limit=n_subst_limit
    )
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
        ext_mod_translated = translate_external_modification(ext_mod, offset=offset)
        offset += ext_mod_translated.get("max_attachment_point_id")
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


def from_smiles_to_pepseq_and_mod_smiles_strings(
    smiles: str, db_json: dict, n_subst_limit=None
) -> tuple:
    """
    Based on SMILES produce
    Output pepseq_string:
    str = string in pepseq format  H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
    where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified
    amino acid; {Cys(R1)} - is amino acid with staple attached, {Cys(R1)} - amino acid with staple attached
    modifications - external ones, with attachment points
    mod_smiles:
    SMILES string (e.g. '[*:1]C[*:2]') - showing the structure of modification with attachment
    points:
    { Cys(R1) } <- is attached in [*:1] attachment point on staple
    { Cys(R2) } <- is attached in [*:2] attachment point on staple

    :param smiles: SMILES code of the peptide
    :type smiles: str

    :param db_json: JSON containing the mapping of symbols to amino acids
    :type db_json: dict

    :param n_subst_limit: int = number of substitutions allowed
    :type n_subst_limit: int

    :return: pepseq_string, mod_smiles: string in pepseq format
    :rtype: tuple
    """

    peptide_json = decompose_peptide_smiles_with_termini(
        smiles, db_json, n_subst_limit=n_subst_limit
    )
    pepseq_format = peptide_json["pepseq_format"]
    mod_smiles_list = [
        ext_mod["smiles"] for ext_mod in peptide_json["external_modifications"]
    ]
    if len(mod_smiles_list) == 1:
        return pepseq_format, mod_smiles_list[0]
    else:
        return pepseq_format, mod_smiles_list


def from_smiles_to_pepseq_and_one_mod_smiles_strings(
    smiles: str, db_json: dict, n_subst_limit=None
) -> tuple:
    """
        Example of string in pepseq format is H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
        where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified
        amino acid; {Cys(R1)} - is amino acid
        modifications - external ones, with attachment points
        mod_smiles: SMILES string (e.g. '[*:1]C[*:2]') - showing the structure of
        modification with attachment points:
        { Cys(R1) } <- is attached in [*:1] attachment point on staple
        { Cys(R2) } <- is attached in [*:2] attachment point on staple

        :param smiles: string of peptide sequence with modified amino acids
            SMILES - string of peptide sequence with modified amino acids
                modification(s)
        :type smiles: str

        :param db_json: database JSON containing the mapping of symbols to amino acids
        :type db_json: dict

        :param n_subst_limit: number of substitutions allowed
        :type n_subst_limit: int

        :return: pepseq_string, mod_smiles_list: string in pepseq format
        :rtype: tuple
    """

    pepseq_format, mod_smiles_list = from_smiles_to_pepseq_and_mod_smiles_strings(
        smiles, db_json, n_subst_limit=n_subst_limit
    )
    if len(mod_smiles_list) == 1:
        return pepseq_format, mod_smiles_list[0]
    else:
        return pepseq_format, mod_smiles_list


def mark_external_modifications_on_seq(
    seq_list: list, peptide_json: dict, mod_as_X: bool = False
) -> list:
    """
    Mark external modifications on sequence

    :param seq_list: list of amino acids
    :type seq_list: list

    :param peptide_json: JSON containing info about modified peptide with 'sequence',
    'internal_modifications', 'external_modifications':

    :param mod_as_X: whether to mark modifications as X
    :type mod_as_X: bool

    :return: seq_list: list of amino acids
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


def mark_internal_modifications_on_seq(
    seq_list: list, peptide_json: dict, mod_as_X: bool = False
) -> list:
    """
    Mark internal modifications on sequence

    :param seq_list: list of amino acids
    :type seq_list: list

    :param peptide_json: JSON containing info about modified peptide with 'sequence',
    'internal_modifications', 'external_modifications':
    :type peptide_json: dict

    :param mod_as_X: bool = whether to mark modifications as X
    :type mod_as_X: bool

    :return: seq_list: list of amino acids
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


def print_sequence(peptide_json: dict, mod_as_X: bool = False) -> str:
    """
    Print sequence
    :param peptide_json: JSON containing info about modified peptide with 'sequence',
    'internal_modifications', 'external_modifications'
    :type peptide_json: dict

    :param mod_as_X: whether to mark modifications as X
    :type mod_as_X: bool

    :return: seq_marked: marked sequence
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
