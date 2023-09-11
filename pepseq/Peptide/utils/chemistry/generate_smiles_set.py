import itertools
from tqdm import tqdm
import rdkit

import networkx as nx

from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import mol_json_to_nx, \
     mol_to_nx, nx_to_mol


def find_atomLabel(G, label='_R1'):
    R = [n for n, v in G.nodes(data=True) if v.get('atomLabel') == label][0]
    return R


class G_NPS(object):
    
    def __init__(self,
                **kwargs):
        self.kwargs = kwargs
        self.R1 = kwargs.get('R1')
        self.R2 = kwargs.get('R2')
        self.R3 = kwargs.get('R3')
        self.R4 = kwargs.get('R4')
        self.R5 = kwargs.get('R5')
        self.R6 = self.kwargs.get('R6')
        self.R7 = self.kwargs.get('R7')
        self.R = self.kwargs.get('R')
        return
    
    def construct(self):
        #self.smi_NPS = '[*]N([*])C([*])([*])C([*])([*])[*][*] |$_R1;;_R2;;_R3;_R4;;_R5;_R6;_R7;_R8$,LO:9:6.10|'
        #self.smi_NPS = '[*]N([*])C([*])([*])C([*])([*])[*][*] |$_R1;;_R2;;_R3;_R4;;_R5;_R6;_R7;_R8$,LO:9:6.10|'
        self.smi_NPS = '[*]N([*])C([*])([*])C([*])([*])[*] |$_R1;;_R2;;_R3;_R4;;_R5;_R6;_R7$,LO:9:6.10|'
        self.mol_NPS = rdkit.Chem.MolFromSmiles(self.smi_NPS)
        self.G_NPS = mol_to_nx(self.mol_NPS)
        self.G = self.G_NPS

        #self.set_A()
        #self.set_R()
        
        for R_id in ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7']:
            R_param = self.__dict__.get(R_id)
            if R_param is None:
                node_id = find_atomLabel(self.G, label='_%s' %R_id)
                self.G.remove_node(node_id)
            else:
                R_node_id_param = self.__dict__.get('%s_node_id' %R_id, 0)
                self.set_Rx(G_R=R_param, R_group_node_id = R_node_id_param, R_id = R_id)
                
        return
    
    def set_Rx(self, G_R, R_group_node_id = 0, R_id = None):
        node_id = find_atomLabel(self.G, label='_%s' %R_id)
        neighbor_node_id = list(self.G.neighbors(node_id))[0]
        G_union = nx.union(self.G, G_R, rename=("", "%s_" %R_id))

        R_group_node_id = '%s_%s' %(R_id, str(R_group_node_id))
        
        G_union.add_edge(
            str(neighbor_node_id), R_group_node_id,
            bond_type=rdkit.Chem.rdchem.BondType.SINGLE)
        
        G_union.remove_node(str(node_id))
        self.G = G_union
        return


def generate_smiles_set(**kwargs):
    R1_smis = kwargs.get('R1_smis')
    R2_smis = kwargs.get('R2_smis')
    R7_smis = kwargs.get('R7_smis')
    
    smi_nps_objs = []
    
    for R1_smi, R2_smi, R7_smi in tqdm(itertools.product(R1_smis, R2_smis, R7_smis)):
        R1 = mol_to_nx(rdkit.Chem.MolFromSmiles(R1_smi))
        R2 = mol_to_nx(rdkit.Chem.MolFromSmiles(R2_smi))
        R7 = mol_to_nx(rdkit.Chem.MolFromSmiles(R7_smi))

        g_nps_obj = G_NPS(R7=R7, R1=R1, R2=R2)
        g_nps_obj.construct()
        #print(list(g_nps_obj.G.nodes(data=True)))
        smi_nps_obj = rdkit.Chem.MolToSmiles(nx_to_mol(g_nps_obj.G))
        smi_nps_objs.append( smi_nps_obj )
    return smi_nps_objs



def check_alkanes(mol_query=None, length_range=(1, 6)):
    min_length, max_length = length_range

    hexane_smi = 'C' * max_length
    m_hexane = rdkit.Chem.MolFromSmiles(hexane_smi)

    if not m_hexane.HasSubstructMatch(mol_query):
        return False
    
    if (min_length is not None) and (min_length > 0):
        min_smi = 'C' * min_length
        min_mol = rdkit.Chem.MolFromSmiles(min_smi)

        if not mol_query.HasSubstructMatch(min_mol):
            return False

    return True


def check_alkenes(mol_query=None, length_range=(0, 6)):
    min_length, max_length = length_range


    min_plus_one_alkyl_or_alkenyl_chain_smarts = '-,='.join( ['[C]'] * (min_length ) )
    mol_plus_one_alkyl_or_alkenyl_chain = rdkit.Chem.MolFromSmarts(min_plus_one_alkyl_or_alkenyl_chain_smarts)

    max_plus_one_alkyl_or_alkenyl_chain_smarts = '-,='.join( ['[C]'] * (max_length + 1) )
    mol_plus_one_alkyl_or_alkenyl_chain = rdkit.Chem.MolFromSmarts(
        max_plus_one_alkyl_or_alkenyl_chain_smarts)

    ethene_mol = rdkit.Chem.MolFromSmiles('C=C')

    if not mol_query.HasSubstructMatch(ethene_mol):
        return False

    if not mol_query.HasSubstructMatch(mol_plus_one_alkyl_or_alkenyl_chain):
        return False

    if mol_query.HasSubstructMatch(mol_plus_one_alkyl_or_alkenyl_chain):
        return False
    
    #graf no edges
    return


def check_alkynes(mol_query=None, length_range=(0, 6)):
    min_length, max_length = length_range


    min_plus_one_alkyl_or_alkenyl_chain_smarts = '-,#'.join( ['[C]'] * (min_length ) )
    mol_plus_one_alkyl_or_alkenyl_chain = rdkit.Chem.MolFromSmarts(min_plus_one_alkyl_or_alkenyl_chain_smarts)

    max_plus_one_alkyl_or_alkenyl_chain_smarts = '-,#'.join( ['[C]'] * (max_length + 1) )
    mol_plus_one_alkyl_or_alkenyl_chain = rdkit.Chem.MolFromSmarts(
        max_plus_one_alkyl_or_alkenyl_chain_smarts)

    ethyne_mol = rdkit.Chem.MolFromSmiles('C#C')

    if not mol_query.HasSubstructMatch(ethyne_mol):
        return False

    if not mol_query.HasSubstructMatch(mol_plus_one_alkyl_or_alkenyl_chain):
        return False

    if mol_query.HasSubstructMatch(mol_plus_one_alkyl_or_alkenyl_chain):
        return False
    
    #graf no edges
    return True


def check_carbon_chain(mol_query=None, length_range=(1,6) ):
    if check_alkanes( length_range=(1,6) ):
        return True
    if check_alkenes( length_range=(1,6) ):
        return True
    if check_alkynes( length_range=(1,6) ):
        return True
    return False


def get_cycloakane_mol(n=3):
    sub_smi = 'C' * (n-1)
    return rdkit.Chem.MolFromSmiles('C1%s1' %sub_smi)


def mols_are_identical(mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol) -> bool:
    are_identical = mol1.HasSubstructMatch(
        mol2, useChirality=True
    ) and mol2.HasSubstructMatch(mol1, useChirality=True)
    return are_identical


def check_cykloalkanes(mol_query=None, length_range=(3,6)):
    min_length, max_length = length_range

    for n in range(min_length, max_length+1):
        cycloalkane_mol = get_cycloakane_mol(n=n)
        if mols_are_identical(mol_query, cycloalkane_mol):
            return True
    return False


def check_benzyl(mol_query=None):
    benzene_mol = rdkit.Chem.MolFromSmiles('C1=CC=CC=C1')
    if mols_are_identical(mol_query, benzene_mol):
        return True

    return False


def check_alkilokarbonyl():
    return


def check_hydroxyl():
    return


def check_amine():
    return


def check_R1():
    if check_carbon_chain( length_range=(1,6) ):
        return True
    
    if check_cykloalkanes( length_range=(1,6) ):
        return True
    
    if check_benzyl():
        return True
    
    if check_alkilokarbonyl():
        return True

    if check_hydroxyl():
        return True
    if check_amine():
        return True
    return


def check_R1_6(G):
    if G is None:
        return True
    
    if check_carbon_chain(G, length_range=(1,6) ):
        return True
    
    if check_cykloalkanes(G, length_range=(1,6) ):
        return True
    
    if check_benzyl(G):
        return True
    
    if check_alkilokarbonyl(G):
        return True

    if check_hydroxyl(G):
        return True
    if check_amine(G):
        return True
    return


def find_Rs(mol_graph, ethylamine):
    return R1_G, R2_G, R3_G, R4_G, R5_G, R6_G


def check_has_ethylamine():
    return


def check_R7(G, templates):
    mol_query = nx_to_mol(G)
    for template in templates:
        if mols_are_identical(mol_query, template):
            return True

    return False


def check_conditions(mol_graph, ring_templates):

    ethylamine_matches = check_has_ethylamine()

    for ethylamine in ethylamine_matches:

        R1_G, R2_G, R3_G, R4_G, R5_G, R6_G = find_Rs(mol_graph, ethylamine)

        for R_G in (R1_G, R2_G, R3_G, R4_G, R5_G, R6_G):
            if not check_R1_6( R_G ):
                return False

        if not check_R7(mol_graph, ring_templates):
            return False
        #check_R8()
    return True
