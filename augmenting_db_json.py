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
from tqdm import tqdm
import rdkit
import networkx as nx
import pandas as pd
import rdkit
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import mol_to_nx, \
    nx_to_mol



def get_R3(G: nx.Graph, label: str = 'R3') -> int:
    """
    Input:
        G: networkx Graph representing (Amino Acid) Molecule
        label: 'dummyLabel' to look for 'R3'

    Returns:
        i: Index of dummyAtom (radical) marked as R3
    
    Action:
        checks 'dummyLabel' parameter for each G graph nodes and return node
        index with given label

    """
    G_nodes = list( G.nodes(data=True) )
    for i, node_data in G_nodes:
        if node_data.get('dummyLabel') == label:
            return i


def set_default_exit_atom(mol: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol:
    """

    Input:
        mol: (Amino Acid) molecule as rdkit object

    Output:
        mol: (Amino Acid) molecule with R3 radical set as DefaultExitAtom for
        attaching modifications to AminoAcid (e.g. palmitoylation; stapling etc.)
    
    Action:
        R3 is identified in molecule and assigned extra parameter DefaultExitAtom = True

    """
    G = mol_to_nx(mol)
    R3_id = get_R3(G)
    nx.set_node_attributes(G, {R3_id: True}, name='DefaultExitAtom')
    return nx_to_mol(G)


def get_basic_smiles(mol: rdkit.Chem.rdchem.Mol) -> str:
    """

    Input:
        mol: (Amino Acid) molecule as rdkit object
    
    Output:
        smi_basic: basic SMILES for molecule without radicals (*)
    
    Actions:
        nodes with 'molFileValue' = '*' are removed from mol

    """
    G = mol_to_nx(mol)
    nodes_to_remove = [
        i for i, node_data in list(
            G.nodes(data=True)) if node_data.get('molFileValue') == '*'
    ]
    for node_id in nodes_to_remove:
        G.remove_node(node_id)
    smi_basic = rdkit.Chem.MolToSmiles(nx_to_mol(G))
    return smi_basic


def get_ro_json(mol: rdkit.Chem.rdchem.Mol, name: str = 'Orn') -> dict:
    """

    Input:
        mol: (Amino Acid) molecule as rdkit object
    
    Output:
        ro_json: dictionary containing basic AminoAcid/Fragment/Group info
        necessary to include it into db.json database of Peptide building blocks
    
    Action:
        molecule with default exit atom (atom to attach modifications unless
          specified otherwise) is computed by set_default_exit_atom_function
        
        from this CXSMILES and CXSMARTS are created together with SMILES code
        depicting basic molecule (without radicals attached).

    """
    mol_set = set_default_exit_atom(mol)
    basic_smiles = get_basic_smiles(mol_set)
    cxsmiles = rdkit.Chem.MolToCXSmiles(mol_set)
    cxsmarts = rdkit.Chem.MolToCXSmarts(mol_set)

    ro_json = {
        'smiles': basic_smiles,
        'smiles_radical': cxsmiles,
        'three_letter_code': name,
        'smarts': cxsmarts,
        'smarts_comments': ''
    }
    return ro_json


def get_ro_json_from_row(row: pd.core.series.Series, colname='m_abbr', mol_colname='ROMol') -> dict:
    """

    """
    mol_name = getattr(row, colname)
    mol = getattr(row, mol_colname)
    ro_json = get_ro_json(
        mol, name = mol_name)
    return ro_json


def get_smiles_set(db_json):
    smiles_set = set()

    for aa in db_json['smiles']['aa']:
        smiles_set.add( db_json['smiles']['aa'][aa]['smiles'] )
    return smiles_set


def get_row_jsons(df_new, colname='m_abbr', mol_colname='ROMol', smiles_set=None):
    row_i_jsons_dict = {} 
    
    rows = list(df_new.iterrows())
    for row in rows:
        row_i = row[1]

        row_i_json = get_ro_json_from_row(
            row_i, colname=colname, mol_colname=mol_colname)
        if row_i_json['smiles'] not in smiles_set:
            code = row_i_json.get('three_letter_code')
            row_i_jsons_dict[code] = copy.deepcopy( row_i_json )
    return row_i_jsons_dict


def augment_db_json(db_json: dict, df_sdf: pd.DataFrame = None, name_column = 'm_abbr', mol_colname='ROMol'):
    """

    """
    aa_keys = set( db_json['smiles']['aa'].keys() )
    df_new = df_sdf[~df_sdf[name_column].isin(aa_keys)]
    smiles_set = get_smiles_set(db_json)
    dict_row_jsons = get_row_jsons(df_new, colname=name_column, mol_colname=mol_colname,
        smiles_set=smiles_set)
    aa_dict = db_json['smiles'].get('aa')
    aa_dict.update( dict_row_jsons )
    db_json['smiles']['aa'] = aa_dict
    return db_json


def N_term_mod_smarts(smarts: str) -> str:
    """

    takes in SMARTS for amino acid
    and returns SMARTS for said amino acid but 
    
    fitting also N-terminal molecule with a N-terminus modification
    
    return SMARTS

    """
    sm_start =   '[$([N&X3&H2,N&X4&H3&+]),$([N&X3&H1](C)C)]'
    sm_changed = '[$([NX3H2,NX4H3+]),$([NX3H](C)(C)),$([NX3H])]'
    
    pattern = rdkit.Chem.RWMol( rdkit.Chem.MolFromSmarts(smarts ))
    at0 = pattern.GetAtomWithIdx(0)
    sm0 = at0.GetSmarts()
    if (sm0 != sm_start):
        return

    neighbors = at0.GetNeighbors()

    if len(neighbors) != 1:
        return
    
    neighbor = neighbors[0]
    pattern.RemoveAtom(0)
    x = pattern.AddAtom( rdkit.Chem.AtomFromSmarts(sm_changed) )
    pattern.AddBond(neighbor.GetIdx(), x, rdkit.Chem.rdchem.BondType.UNSPECIFIED)
    return rdkit.Chem.MolToCXSmarts(pattern)


def get_Nter_versions_cxsmarts_db(d: dict) -> dict:
    """

    takes in SMARTS dictionary for amino acid
    and returns SMARTS dictionary for said amino acids but 
    
    fitting also N-terminal molecule with a N-terminus modification
    
    return SMARTS
    
    """
    d_copy = {}
    for aa_code in d:
        smarts = d.get(aa_code)
        new_smarts = N_term_mod_smarts(smarts)
        if new_smarts is not None:
            d_copy[aa_code] = new_smarts 
    return d_copy


def get_Nter_versions(d: dict) -> dict:
    """

    takes in SMARTS dictionary for amino acid
    and returns SMARTS dictionary for said amino acids but 
    
    fitting also N-terminal molecule with a N-terminus modification
    
    return SMARTS
    
    """
    d_copy = {}
    for aa_code in d:
        aa_json = d.get(aa_code).copy()
        smarts = aa_json.get('smarts')
        new_smarts = N_term_mod_smarts(smarts)
        if new_smarts is not None:
            d_copy[aa_code] = new_smarts 
    return d_copy


def change_exit_atom(smiles: str) -> str:
    """

    SDF file has explicit radical (R3) connected to Exit Atom
    We change this marking by labelling a neighbouring Atom
    as DefaultExitAtom and removing the radical R3 from SMILES

    """
    mol = rdkit.Chem.MolFromSmiles(smiles)
    G = mol_to_nx(mol)
    nodes_data = list(G.nodes(data=True))
    R3_id = [node_id for node_id, node_data in nodes_data if node_data.get(
        'dummyLabel') == 'R3'][0]
    def_exit_atom = list(G.neighbors( R3_id ))[0]
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
    
    """
    db_json_copy = db_json.copy()
    smi_dict = db_json_copy.get('smiles').get('aa')
    for aa_code in smi_dict:
        smiles_radical = smi_dict.get(aa_code).get('smiles_radical')
        if '[3*]' in smiles_radical:
            smiles_radical_changed = change_exit_atom(smiles_radical)
            db_json_copy['smiles']['aa'][aa_code][
                'smiles_radical'] = smiles_radical_changed
    return db_json_copy


def order_aas(db_json: dict) -> list:
    """
    
    Some Amino Acids can be contained within larger
    Amino Acids. (
        e.g. Alanine within Serine, 
             Glycine within Alanine
             )
    
    When decomposing a Peptide Molecule into Monomers.
    Same sidechain would fit match few amino acids:
        e.g. Serine will also match Alanine SMARTS pattern
    
    In this case we select the largest Amino Acid
    that has been matched. In order to do that 
    we order Amino Acids from smallest to largest.
    
    """
    aa_smiles = db_json.get('smiles').get('aa')
    aa_codes = sorted(aa_smiles.keys())

    mols = [rdkit.Chem.MolFromSmiles(aa_smiles.get(
        aa_code).get('smiles'))
        for aa_code in aa_codes]
    n_codes = len(aa_codes)

    edges = []

    for i in tqdm(range(n_codes)):
        for j in range(n_codes):
            if i!= j:
                if mols[i].HasSubstructMatch(mols[j]):
                    if not mols[j].HasSubstructMatch(mols[i]):
                        edges.append( (j,i) )

    G = nx.DiGraph()
    for edge in edges:
        G.add_edge(*edge)
    leaves = [x for x in G.nodes() if G.out_degree(x)==0]
    biggest = []
    while leaves:
        biggest += leaves
        for leaf in leaves: 
            G.remove_node(leaf)
        leaves = [x for x in G.nodes() if G.out_degree(x)==0]
    biggest += leaves
    new_aa_order = list(reversed([aa_codes[i] for i in biggest]))
    return new_aa_order


def remove_radicals(smarts: str) -> str:
    """

    remove radical(s) from (CX)SMARTS code
    
    """
    m = rdkit.Chem.MolFromSmarts(smarts)
    atoms = list(m.GetAtoms())
    radicals = [ atom.GetIdx() for atom in atoms if atom.HasProp('dummyLabel') ]
    if not radicals:
        return smarts
    mw = rdkit.Chem.RWMol(m)
    for atom_id in sorted(radicals)[::-1]:
        mw.RemoveAtom(atom_id)
    return rdkit.Chem.MolToCXSmarts(mw)


def remove_radicals_from_db_smarts(db_json: dict) -> dict:
    """

    SDF file SMARTS has radicals that might interfere with matching
    them as SMARTS substructure pattern. We need to remove them.
        
    """
    db_json_copy = db_json.copy()
    aa_smiles = db_json_copy['smiles']['aa']
    for aa_code in aa_smiles:
        new_smarts = remove_radicals( aa_smiles.get(aa_code).get('smarts') )
        db_json_copy['smiles']['aa'][aa_code]['smarts'] = new_smarts
    return db_json_copy