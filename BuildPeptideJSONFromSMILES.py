import rdkit

from Backbone import BreakingIntoResidueCandidateSubgraphs, MarkingPeptideBackbone
from Peptide.utils.chemistry.ProcessResidueCandidateGraph import (
    decompose_residues_internal,
)
from Peptide.utils.Parser import parse_canonical2


def get_cx_smarts_db(db_json):

    cx_smarts_db = {}
    smiles_dict = db_json["smiles"].get("aa")

    for aa in smiles_dict:
        cx_smarts_db[aa] = smiles_dict[aa].get("smarts")

    return cx_smarts_db


def decompose_peptide_smiles(smiles, db_json):
    peptide_molecule = rdkit.Chem.MolFromSmiles(smiles)
    peptide_molecule = MarkingPeptideBackbone().execute(peptide_molecule)

    residues = BreakingIntoResidueCandidateSubgraphs().execute(peptide_molecule)

    cx_smarts_db = get_cx_smarts_db(db_json)

    seq, internal_modifications, external_modifications = decompose_residues_internal(
        residues, cx_smarts_db
    )
    return {
        "sequence": seq,
        "internal_modifications": internal_modifications,
        "external_modifications": external_modifications,
    }


def get_terminal_smiles_building_block(peptide_json, ResID, AtomName):
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


def get_C_terminal_smiles_building_block(peptide_json):
    return get_terminal_smiles_building_block(peptide_json, ResID=-1, AtomName="CO")


def get_N_terminal_smiles_building_block(peptide_json):
    return get_terminal_smiles_building_block(peptide_json, ResID=1, AtomName="N")


def smiles_are_identical(smiles1, smiles2):
    mol1 = rdkit.Chem.MolFromSmiles(smiles1)
    mol2 = rdkit.Chem.MolFromSmiles(smiles2)
    return (mol1.HasSubstructMatch(mol2)) and (mol2.HasSubstructMatch(mol1))


def get_term_symbol(smiles, db_json, group):
    terms = db_json.get("smiles").get(group)

    for term_symbol in terms:
        term_smiles = terms.get(term_symbol).get("smiles_radical")
        if smiles_are_identical(smiles, term_smiles):
            return term_symbol


def get_c_term_from_peptide_json(peptide_json, db_json):
    tup = get_C_terminal_smiles_building_block(peptide_json)

    if tup is not None:
        smiles, i = tup
        term_symbol = get_term_symbol(smiles, db_json, "c_terms")
        return term_symbol, i


def get_n_term_from_peptide_json(peptide_json, db_json):
    tup = get_N_terminal_smiles_building_block(peptide_json)
    if tup is not None:
        smiles, i = tup
        term_symbol = get_term_symbol(smiles, db_json, "n_terms")
        return term_symbol, i


def decompose_peptide_smiles_with_termini(smiles, db_json):
    peptide_json = decompose_peptide_smiles(smiles, db_json)
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

    return peptide_json


def mark_external_modifications_on_seq(l, peptide_json, mod_as_X=False):
    external_modifications = peptide_json.get("external_modifications")
    for external_modification in external_modifications:
        attachment_points = external_modification.get("attachment_points_on_sequence")

        for key in attachment_points:
            attachment_point = attachment_points.get(key)
            res_id = int(attachment_point.get("ResID"))

            basic_res_name = l[res_id - 1]
            if mod_as_X:
                new_res_name = "X"
            else:
                new_res_name = "{mod%s}" % (basic_res_name)

            l[res_id - 1] = new_res_name
    return l


def mark_internal_modifications_on_seq(l, peptide_json, mod_as_X=False):
    internal_modifications = peptide_json.get("internal_modifications")

    for key in internal_modifications:
        mod_residues = internal_modifications[key]
        mod_residue_ids = [int(i["ResID"]) for i in mod_residues]

        for res_id in mod_residue_ids:
            basic_res_name = l[res_id - 1]
            if mod_as_X:
                new_res_name = "X"
            else:
                new_res_name = "{mod%s}" % (basic_res_name)
            l[res_id - 1] = new_res_name
    return l


def print_sequence(peptide_json, mod_as_X=False):
    basic_sequence = peptide_json.get("sequence")
    l = list(basic_sequence)

    N_terminus = peptide_json.get("N_terminus")
    C_terminus = peptide_json.get("C_terminus")

    l = mark_external_modifications_on_seq(l, peptide_json, mod_as_X)
    l = mark_internal_modifications_on_seq(l, peptide_json, mod_as_X)

    seq_marked = "".join(l)
    if N_terminus is not None:
        seq_marked = "%s~%s" % (N_terminus, seq_marked)
    if C_terminus is not None:
        seq_marked = "%s~%s" % (seq_marked, C_terminus)

    return seq_marked
