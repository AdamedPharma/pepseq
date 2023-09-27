import copy
import rdkit
import networkx as nx
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import mol_to_nx, nx_to_mol



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
