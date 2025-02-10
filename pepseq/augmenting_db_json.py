"""
#**pepseq.augmenting_db_json**


Provide several sample math calculations.

This module allows the user to make mathematical calculations.


Examples:
    >>> from calculator import calculations
    >>> calculations.add(2, 4)
    6.0
    >>> calculations.multiply(2.0, 4.0)
    8.0
    >>> from calculator.calculations import divide
    >>> divide(4.0, 2)
    2.0

The module contains the following functions:

- `add(a, b)` - Returns the sum of two numbers.
- `subtract(a, b)` - Returns the difference of two numbers.
- `multiply(a, b)` - Returns the product of two numbers.
- `divide(a, b)` - Returns the quotient of two numbers.
"""

import copy

from typing import Union

from tqdm import tqdm
import rdkit
import networkx as nx
import pandas as pd
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import mol_to_nx, nx_to_mol


def get_R3(G: nx.Graph, label: str = "R3") -> int:
    """
    Checks 'dummyLabel' parameter for each G graph nodes and return node
     index with given label

    :param  G: networkx Graph representing (Amino Acid) Molecule
    :type   G: nx.Graph
    :param label: 'dummyLabel' to look for 'R3'
    :type  label: str

    :return i: Index of dummyAtom (radical) marked as R3
    :rtype: int
    """
    G_nodes = list(G.nodes(data=True))
    for i, node_data in G_nodes:
        if node_data.get("dummyLabel") == label:
            return i


def set_default_exit_atom(mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol:
    """
    R3 is identified in molecule and assigned extra parameter DefaultExitAtom = True

    :param mol: (Amino Acid) molecule as rdkit object
    :type  mol: rdkit.Chem.rdchem.Mol

    :return mol: (Amino Acid) molecule with R3 radical set as DefaultExitAtom for
        attaching modifications to AminoAcid (e.g. palmitoylation; stapling etc.)
    :rtype: rdkit.Chem.rdchem.Mol
    """
    G = mol_to_nx(mol)
    R3_id = get_R3(G)
    nx.set_node_attributes(G, {R3_id: True}, name="DefaultExitAtom")
    return nx_to_mol(G)


def get_basic_smiles(mol: rdkit.Chem.rdchem.Mol) -> str:
    """
    Remove nodes with 'molFileValue' = '*' from molecule

    :param mol: (Amino Acid) molecule as rdkit object
    :type  mol: rdkit.Chem.rdchem.Mol

    :return smi_basic: basic SMILES for molecule without radicals (*)
    :rtype: str
    """
    G = mol_to_nx(mol)
    nodes_to_remove = [
        i
        for i, node_data in list(G.nodes(data=True))
        if node_data.get("molFileValue") == "*"
    ]
    for node_id in nodes_to_remove:
        G.remove_node(node_id)
    smi_basic = rdkit.Chem.MolToSmiles(nx_to_mol(G))
    return smi_basic


def get_ro_json(mol: rdkit.Chem.rdchem.Mol, name: str = "Orn") -> dict:
    """
    Molecule with default exit atom (atom to attach modifications unless
    specified otherwise) is computed by set_default_exit_atom_function.
    From this CXSMILES and CXSMARTS are created together with SMILES code
    depicting basic molecule (without radicals attached).

    :param mol: (Amino Acid) molecule as rdkit object
    :type mol: rdkit.Chem.rdchem.Mol

    :param name: monomer/building block name
    :type  name: str

    :return ro_json: dictionary containing basic AminoAcid/Fragment/Group info
        necessary to include it into db.json database of Peptide building blocks
    :rtype: dict
    """

    mol_set = set_default_exit_atom(mol)
    basic_smiles = get_basic_smiles(mol_set)
    cxsmiles = rdkit.Chem.MolToCXSmiles(mol_set)
    cxsmarts = rdkit.Chem.MolToCXSmarts(mol_set)

    ro_json = {
        "smiles": basic_smiles,
        "smiles_radical": cxsmiles,
        "three_letter_code": name,
        "smarts": cxsmarts,
        "smarts_comments": "",
    }
    return ro_json


def get_ro_json_from_row(
    row: pd.core.series.Series, colname="m_abbr", mol_colname="ROMol"
) -> dict:
    """
    Process pandas.DataFrame row into dict dictionary containing basic AminoAcid/Fragment/Group info
        necessary to include it into db.json database of Peptide building blocks

    :param row: pandas DataFrame row
    :type row: pd.core.series.Series

    :param colname: name of DataFrame column containing name for the monomer/building block
    :type  colname: str

    :param mol_colname: name of DataFrame column containing monomer/building block rdkit.Mol object
    :type  mol_colname: str

    :return ro_json: dictionary containing basic AminoAcid/Fragment/Group info
        necessary to include it into db.json database of Peptide building blocks
    :rtype: dict
    """
    mol_name = getattr(row, colname)
    mol = getattr(row, mol_colname)
    ro_json = get_ro_json(mol, name=mol_name)
    return ro_json


def get_smiles_set(db_json: dict) -> set:
    """
    Return set of all SMILES codes for all monomers/building blocks found in db_json database

    :param db_json database containing info on monomers/building blocks for Modified Peptides
    :type  db_json: dict

    :return smiles_set: set of all SMILES codes for all monomers/building blocks found in db_json database
    :rtype: set
    """
    smiles_set = set()

    for aa in db_json["smiles"]["aa"]:
        smiles_set.add(db_json["smiles"]["aa"][aa]["smiles"])
    return smiles_set


def get_row_jsons(
    new_monomers_dataframe: pd.DataFrame,
    colname="m_abbr",
    mol_colname="ROMol",
    smiles_set=None,
) -> dict:
    """
    Create JSON for database of new monomers to be added to monomers database (db_json). Add only monomers
    not yet present in db_json (check for duplicate SMILES codes).

    :param df_new: DataFrame
    :type df_new: pd.DataFrame

    :param colname: name of DataFrame column containing name for the monomer/building block
    :type  colname: str

    :param mol_colname: name of DataFrame column containing monomer/building block rdkit.Mol object
    :type  mol_colname: str

    :return: new_monomers_db_json - database of new monomers
    :rtype: dict
    """
    new_monomers_db_json = {}

    rows = list(new_monomers_dataframe.iterrows())

    for row in rows:
        row_index = row[1]

        monomer_json = get_ro_json_from_row(
            row_index, colname=colname, mol_colname=mol_colname
        )
        if monomer_json["smiles"] not in smiles_set:
            code = monomer_json.get("three_letter_code")
            new_monomers_db_json[code] = copy.deepcopy(monomer_json)
    return new_monomers_db_json


def augment_db_json(
    db_json: dict,
    df_sdf: pd.DataFrame = None,
    name_column="m_abbr",
    mol_colname="ROMol",
) -> dict:
    """
    Insert new monomers to monomers database through pandas DataFrame

    :param db_json: monomers database
    :type  db_json: dict

    :param df_sdf: pandas DataFrame with rows containing info on monomers read from SDF file
    :type  df_sdf: pd.DataFrame

    :param name_column: name of DataFrame column containing name for the monomer/building block
    :type  name_column: str

    :param mol_colname: name of DataFrame column containing monomer/building block rdkit.Mol object
    :type  mol_colname: str

    :return: db_json - database of monomers inserted with new monomers read
      from SDF file through pandas DataFrame
    :rtype:  dict
    """
    aa_keys = set(db_json["smiles"]["aa"].keys())
    df_new = df_sdf[~df_sdf[name_column].isin(aa_keys)]
    smiles_set = get_smiles_set(db_json)
    dict_row_jsons = get_row_jsons(
        df_new, colname=name_column, mol_colname=mol_colname, smiles_set=smiles_set
    )
    aa_dict = db_json["smiles"].get("aa")
    aa_dict.update(dict_row_jsons)
    db_json["smiles"]["aa"] = aa_dict
    return db_json


def replace_atom(mol: rdkit.Chem.rdchem.Mol, atom_id: int, atom_smarts: str) -> str:
    """
    Update rdkit molecule pattern by replacing atom pattern and return SMARTS code

    :param mol: amino acid molecule pattern with N terminal hydrogen pattern/conditions
    to be expanded to include N terminally modified nitrogen
    :type  mol: rdkit.Chem.rdchem.Mol

    :param atom_id: id of atom (pattern) to be replaced
    :type  atom_id: int

    :param atom_smarts: SMARTS code for pattern to be assigned to atom
    :type  atom_smarts: str

    :return mol_new_smarts: SMARTS code for new updated molecule pattern
    :rtype: str
    """

    atom = mol.GetAtomWithIdx(atom_id)
    neighbors = atom.GetNeighbors()
    mol.RemoveAtom(atom_id)

    new_atom = rdkit.Chem.AtomFromSmarts(atom_smarts)
    new_atom_added = mol.AddAtom(new_atom)

    for neighbor in neighbors:
        neighbor_atom_id = neighbor.GetIdx()
        mol.AddBond(
            neighbor_atom_id, new_atom_added, rdkit.Chem.rdchem.BondType.UNSPECIFIED
        )

    mol_new_smarts = rdkit.Chem.MolToCXSmarts(mol)
    return mol_new_smarts


def N_term_mod_smarts(smarts: str) -> Union[str, None]:
    """
    Process amino acid SMARTS pattern code to match amino acid also after
    N terminal modification. We utilize the fact that N terminal N atom
    is often marked with Idx=0. Pattern is read. If N terminus is not found
    within pattern None is returned.

    takes in SMARTS for amino acid
    and returns SMARTS for said amino acid but

    fitting also N-terminal molecule with a N-terminus modification

    :param smarts: SMARTS code
    :type  smarts: str

    :return mol_new_smarts: new SMARTS code updated with N terminally modified N atom pattern
    :rtype: str
    """
    sm_start = "[$([N&X3&H2,N&X4&H3&+]),$([N&X3&H1](C)C)]"
    sm_changed = "[$([NX3H2,NX4H3+]),$([NX3H](C)(C)),$([NX3H])]"

    mol_pattern = rdkit.Chem.RWMol(rdkit.Chem.MolFromSmarts(smarts))
    atom_id = 0

    at0 = mol_pattern.GetAtomWithIdx(atom_id)
    sm0 = at0.GetSmarts()
    if sm0 != sm_start:
        return

    mol_new_smarts = replace_atom(mol_pattern, atom_id, sm_changed)

    return mol_new_smarts


def get_Nter_versions_cxsmarts_db(aa_smarts_dict: dict) -> dict:
    """
    Update amino acid SMARTS pattenr database to match N terminal modifications.

    :param aa_smarts_dict: dictionary database of SMARTS codes for amino_acids
    :type  aa_smarts_dict: dict

    :return: updated_aa_smarts_dict: updated dictionary of SMARTS codes
      for amino acids with patterns changed to include N terminal modifications
    :rtype: dict
    """
    updated_aa_smarts_dict = {}
    for aa_code in aa_smarts_dict:
        smarts = aa_smarts_dict.get(aa_code)
        new_smarts = N_term_mod_smarts(smarts)
        if new_smarts is not None:
            updated_aa_smarts_dict[aa_code] = new_smarts
    return updated_aa_smarts_dict


def get_Nter_versions(aa_smarts_dict: dict) -> dict:
    """
    Update amino acid SMARTS pattern provided in db_json format database to match N terminal modifications.

    :param aa_smarts_dict: dictionary database of SMARTS codes for amino_acids: dict

    :return: updated_aa_smarts_dict: updated dictionary of SMARTS codes
      for amino acids with patterns changed to include N terminal modifications
    :rtype: dict
    """

    updated_aa_smarts_dict = {}
    for aa_code in aa_smarts_dict:
        aa_json = aa_smarts_dict.get(aa_code).copy()
        smarts = aa_json.get("smarts")
        new_smarts = N_term_mod_smarts(smarts)
        if new_smarts is not None:
            updated_aa_smarts_dict[aa_code] = new_smarts
    return updated_aa_smarts_dict


def change_exit_atom(smiles: str) -> str:
    """

    SDF file has explicit radical (R3) connected to Exit Atom
    We change this marking by labelling a neighbouring Atom
    as DefaultExitAtom and removing the radical R3 from SMILES

    :param smiles: SMILES code for amino acid with R3 radical connected to Exit Atom
    :type  smiles: str

    :return smiles: SMILES code for amino acid with removed R3 but Exit Atom labeled as such
    :rtype: str
    """
    mol = rdkit.Chem.MolFromSmiles(smiles)
    G = mol_to_nx(mol)
    nodes_data = list(G.nodes(data=True))
    R3_id = [
        node_id
        for node_id, node_data in nodes_data
        if node_data.get("dummyLabel") == "R3"
    ][0]
    def_exit_atom = list(G.neighbors(R3_id))[0]
    nx.set_node_attributes(G, {def_exit_atom: True}, name="DefaultExitAtom")
    G.remove_node(R3_id)
    mol = nx_to_mol(G)
    smiles = rdkit.Chem.MolToCXSmiles(mol)
    return smiles


def change_exit_atoms(db_json: dict) -> dict:
    """

    SDF file has explicit radical (R3) connected to Exit Atom
    We change this marking by labelling a neighbouring Atom
    as DefaultExitAtom and removing the radical R3 from SMILES

    Done across all the monomer database

    :param db_json: database containing info on Modified Peptide monomers/building blocks
    :type  db_json: dict

    :return: db_json_copy database updated by marking exit atom and removing R3 radicals
    :rtype:  dict

    """
    db_json_copy = db_json.copy()
    smi_dict = db_json_copy.get("smiles").get("aa")
    for aa_code in smi_dict:
        smiles_radical = smi_dict.get(aa_code).get("smiles_radical")
        if "[3*]" in smiles_radical:
            smiles_radical_changed = change_exit_atom(smiles_radical)
            db_json_copy["smiles"]["aa"][aa_code][
                "smiles_radical"
            ] = smiles_radical_changed
    return db_json_copy


def get_substructure_relations(mols: list[rdkit.Chem.rdchem.Mol]) -> list:
    """
    List all cases where one amino acid is a sbustructure
    of others (e.g. glycine in serine or alanine in phenylalanine)
    (j, i )

    :param mols: list of rdkit.Chem.rdchem.Mol
    :type  mols: list[rdkit.Chem.rdchem.Mol]

    :return substructure_relations: all cases where one amino acid is a substructure
    of others (e.g. glycine in serine or alanine in phenylalanine)
    :rtype: list
    """
    substructure_relations = []
    n_mols = len(mols)
    for i in tqdm(range(n_mols)):
        mol_i = mols[i]
        for j in tqdm(range(n_mols)):
            if i != j:
                mol_j = mols[j]
                if mol_i.HasSubstructMatch(mol_j):
                    new_relation = (j, i)
                    substructure_relations.append(new_relation)
    return substructure_relations


def order_graph_nodes_from_root_to_leaves(G: nx.DiGraph, aa_codes: list[str]) -> list:
    """
    :param G: directed graph of substructure relations
    :type G:  nx.DiGraph

    :param aa_codes: list of amino acid nodes names
    :type  aa_codes: list: str

    :return new_aa_order: list of ordered amino acids from smallest to biggest
    :rtype: list

    """
    leaves = [x for x in G.nodes() if G.out_degree(x) == 0]
    leaves_to_root = []
    while leaves:
        leaves_to_root += leaves
        for leaf in leaves:
            G.remove_node(leaf)
        leaves = [x for x in G.nodes() if G.out_degree(x) == 0]
    leaves_to_root += leaves
    root_to_leaves = reversed(leaves_to_root)
    new_aa_order = [aa_codes[i] for i in root_to_leaves]
    return new_aa_order


def order_aas(db_json: dict) -> list:
    """
    Some Amino Acids can be contained within larger
    Amino Acids. ( e.g. Alanine within Serine,  Glycine within Alanine)

    When decomposing a Peptide Molecule into Monomers.
    Same sidechain would fit match few amino acids:
    e.g. Serine will also match Alanine SMARTS pattern

    In this case we select the largest Amino Acid
    that has been matched. In order to do that
    we order Amino Acids from smallest to largest.

    :param db_json: database containing info on Modified Peptide monomers/building blocks
    :type  db_json: dict

    :return new_aa_order: list of ordered amino acids from smallest to biggest
    :rtype: list

    """
    aa_smiles = db_json.get("smiles").get("aa")
    aa_codes = sorted(aa_smiles.keys())

    mols = [
        rdkit.Chem.MolFromSmiles(aa_smiles.get(aa_code).get("smiles"))
        for aa_code in aa_codes
    ]

    substructure_relations = get_substructure_relations(mols)

    G = nx.DiGraph()
    for substructure_relation in substructure_relations:
        G.add_edge(*substructure_relation)

    new_aa_order = order_graph_nodes_from_root_to_leaves(G, aa_codes)

    return new_aa_order


def remove_radicals(smarts: str) -> str:
    """

    Remove radical(s) from (CX)SMARTS code

    :param smarts: SMARTS code for amino acid
    :type  smarts: str

    :return: cx_smarts - CXSMARTS code for amino acid with removed radicals / dummyAtoms
    :rtype: str

    """
    m = rdkit.Chem.MolFromSmarts(smarts)
    atoms = list(m.GetAtoms())
    radicals = [atom.GetIdx() for atom in atoms if atom.HasProp("dummyLabel")]
    if not radicals:
        return smarts
    mw = rdkit.Chem.RWMol(m)
    for atom_id in sorted(radicals)[::-1]:
        mw.RemoveAtom(atom_id)
    cx_smarts = rdkit.Chem.MolToCXSmarts(mw)
    return cx_smarts


def remove_radicals_from_db_smarts(db_json: dict) -> dict:
    """

    SDF file SMARTS has radicals that might interfere with matching
    them as SMARTS substructure pattern. We need to remove them.

    :param db_json: database containing info on amino acids structures
    :type  db_json: dict

    :return db_json_updated: database where radicals have been removed from SMARTS codes
    :rtype: dict
    """
    db_json_updated = db_json.copy()
    aa_smiles = db_json_updated["smiles"]["aa"]
    for aa_code in aa_smiles:
        new_smarts = remove_radicals(aa_smiles.get(aa_code).get("smarts"))
        db_json_updated["smiles"]["aa"][aa_code]["smarts"] = new_smarts
    return db_json_updated
