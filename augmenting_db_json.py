import copy
from tqdm import tqdm
import rdkit
import networkx as nx
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import mol_to_nx, \
    nx_to_mol



def get_R3(G_orn):
    G_orn_nodes = list( G_orn.nodes(data=True) )
    for i, node_data in G_orn_nodes:
        if node_data.get('dummyLabel') == 'R3':
            return i


def set_default_exit_atom(ro_mol):
    G_orn = mol_to_nx(ro_mol)
    R3_id = get_R3(G_orn)
    nx.set_node_attributes(G_orn, {R3_id: True}, name='DefaultExitAtom')
    return nx_to_mol(G_orn)


def get_basic_smiles(mol):
    G_set = mol_to_nx(mol)
    nodes_to_remove = [
        i for i, node_data in list(
            G_set.nodes(data=True)) if node_data.get('molFileValue') == '*'
    ]
    for node_id in nodes_to_remove:
        G_set.remove_node(node_id)
    smi_basic = rdkit.Chem.MolToSmiles(nx_to_mol(G_set))
    return smi_basic


def get_ro_json(ro_mol, name = 'Orn'):
    ro_mol_set = set_default_exit_atom(ro_mol)
    basic_smiles = get_basic_smiles(ro_mol_set)
    cxsmiles = rdkit.Chem.MolToCXSmiles(ro_mol_set)
    cxsmarts = rdkit.Chem.MolToCXSmarts(ro_mol_set)

    ro_json = {
        'smiles': basic_smiles,
        'smiles_radical': cxsmiles,
        'three_letter_code': name,
        'smarts': cxsmarts,
        'smarts_comments': ''
    }
    return ro_json


def get_ro_json_from_row(row, colname='m_abbr', mol_colname='ROMol'):
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


def augment_db_json(db_json, df_sdf=None, name_column = 'm_abbr', mol_colname='ROMol'):
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

    SDF file SMARTS has radicals that might interfere with matching
    them as SMARTS substructure pattern
    
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