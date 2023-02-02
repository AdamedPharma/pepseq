import rdkit

from Backbone import BreakingIntoResidueCandidateSubgraphs, MarkingPeptideBackbone
from Peptide.utils.chemistry.ProcessResidueCandidateGraph import (
    decompose_residues_internal,
)


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
        ResID = len(peptide_json.get("sequence"))

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
    smiles, i = get_C_terminal_smiles_building_block(peptide_json)
    if smiles is not None:
        term_symbol = get_term_symbol(smiles, db_json, "c_terms")
        return term_symbol, i


def get_n_term_from_peptide_json(peptide_json, db_json):
    smiles, i = get_N_terminal_smiles_building_block(peptide_json)
    if smiles is not None:
        term_symbol = get_term_symbol(smiles, db_json, "n_terms")
        return term_symbol, i


def decompose_peptide_smiles_with_termini(smiles, db_json):
    peptide_json = decompose_peptide_smiles(smiles, db_json)
    c_terminus, c_ind = get_c_term_from_peptide_json(peptide_json, db_json)
    n_terminus, n_ind = get_n_term_from_peptide_json(peptide_json, db_json)
    external_modifications = peptide_json.get("external_modifications")

    to_del = []

    if c_terminus is not None:
        peptide_json["C_terminus"] = c_terminus
        to_del.append(c_ind)
    else:
        c_terminus = "OH"

    if n_terminus is not None:
        peptide_json["N_terminus"] = n_terminus
        to_del.append(n_ind)

    new_ext_mods = [
        external_modifications[i]
        for i in range(len(external_modifications))
        if i not in to_del
    ]
    peptide_json["external_modifications"] = new_ext_mods

    return peptide_json
