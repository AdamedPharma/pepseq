import networkx as nx
import rdkit
import rdkit.Chem
import pandas as pd


def mol_to_nx_json(mol: rdkit.Chem.rdchem.Mol) -> dict:
    for atom in mol.GetAtoms():

        kwargs = {}

        for prop in atom.GetPropNames():
            prop_val = atom.GetProp(prop)
            kwargs[prop] = prop_val

    return


def mol_to_nx(mol: rdkit.Chem.rdchem.Mol) -> nx.classes.graph.Graph:
    """
    transforms rdkit.Chem.rdchem.Mol molecule
    into nx.classes.graph.Graph graph
    """

    G = nx.Graph()

    for atom in mol.GetAtoms():

        kwargs = {}

        for prop in atom.GetPropNames():
            prop_val = atom.GetProp(prop)
            kwargs[prop] = prop_val

        if kwargs.get("num_explicit_hs") is not None:
            kwargs.pop("num_explicit_hs")

        if kwargs.get("hybridization") is not None:
            kwargs.pop("hybridization")

        G.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            formal_charge=atom.GetFormalCharge(),
            chiral_tag=atom.GetChiralTag(),
            hybridization=atom.GetHybridization(),
            num_explicit_hs=atom.GetNumExplicitHs(),
            is_aromatic=atom.GetIsAromatic(),
            isotope=atom.GetIsotope(),
            **kwargs
        )

    for bond in mol.GetBonds():

        kwargs = {}
        for prop in bond.GetPropNames():
            prop_val = bond.GetProp(prop)
            kwargs[prop] = prop_val

        G.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_type=bond.GetBondType(),
            **kwargs
        )
    return G


def nx_to_mol(G: nx.classes.graph.Graph) -> rdkit.Chem.rdchem.Mol:
    """
    transforms nx.classes.graph.Graph graph
    into  rdkit.Chem.rdchem.Mol molecule
    """

    mol = rdkit.Chem.RWMol()
    atomic_nums = nx.get_node_attributes(G, "atomic_num")
    chiral_tags = nx.get_node_attributes(G, "chiral_tag")
    formal_charges = nx.get_node_attributes(G, "formal_charge")
    node_is_aromatics = nx.get_node_attributes(G, "is_aromatic")
    node_hybridizations = nx.get_node_attributes(G, "hybridization")
    num_explicit_hss = nx.get_node_attributes(G, "num_explicit_hs")
    isotope = nx.get_node_attributes(G, "isotope")

    basic_propNames = set(
        [
            "atomic_num",
            "chiral_tag",
            "formal_charge",
            "is_aromatic",
            "hybridization" "num_explicit_hs",
            "isotope",
        ]
    )

    node_to_idx = {}
    for node in G.nodes():
        a = rdkit.Chem.Atom(atomic_nums[node])

        additional_propNames = G.nodes[node].keys() - basic_propNames
        for propName in additional_propNames:
            a.SetProp(propName, str(G.nodes[node][propName]))

        a.SetChiralTag(chiral_tags[node])
        a.SetFormalCharge(formal_charges[node])
        a.SetIsAromatic(node_is_aromatics[node])
        a.SetHybridization(node_hybridizations[node])
        a.SetNumExplicitHs(num_explicit_hss[node])
        a.SetIsotope(isotope[node])
        idx = mol.AddAtom(a)
        node_to_idx[node] = idx

    basic_propNames = set(["bond_type"])

    bond_types = nx.get_edge_attributes(G, "bond_type")
    for edge in G.edges():
        first, second = edge
        ifirst = node_to_idx[first]
        isecond = node_to_idx[second]
        bond_type = bond_types[first, second]
        mol.AddBond(ifirst, isecond, bond_type)
        # find BondIdx

    rdkit.Chem.SanitizeMol(mol)
    return mol



def nx_to_json(G: nx.classes.graph.Graph) -> dict:
    nodes_list = list( G.nodes(data=True) )
    node_dicts = [ node_tuple[1] for node_tuple in nodes_list ]
    node_ids = [node_tuple[0] for node_tuple in nodes_list ]

    df = pd.DataFrame(node_dicts)
    df['node_id'] = node_ids

    rows = df.where(df.notnull(),None).values.tolist()
    nodes_tuple = tuple([tuple(row) for row in rows])
    nodes_columns = ['atomic_num', 'formal_charge', 'chiral_tag',
        'hybridization', 'num_explicit_hs', 'is_aromatic',
        'isotope', 'AtomName', 'ResID', 'node_id']

    edges_list = list(G.edges(data=True))

    df_edges = pd.DataFrame([ i[2] for i in edges_list ])

    bond_start =  [i[0] for i in edges_list]
    bond_end =  [i[1] for i in edges_list]

    df_edges['bond_start'] = bond_start
    df_edges['bond_end'] = bond_end


    df_edges = df_edges.where(df_edges.notnull(), False)

    edges_tuple = tuple([tuple(i) for i in df_edges.values.tolist()])


    edges_columns = ['bond_type', 'is_peptide_bond', 'bond_start', 'bond_end']

    nodes_columns = ['atomic_num', 'formal_charge', 'chiral_tag', 'hybridization',
                    'num_explicit_hs', 'is_aromatic', 'isotope', 'AtomName',
                    'ResID', 'node_id']


    mol_j = {
        'nodes_tuple': nodes_tuple,
        'nodes_columns': nodes_columns,
        'edges_tuple': edges_tuple,
        'edges_columns': edges_columns,
    }
    return mol_j


def get_mol_json(mol):
    G = mol_to_nx(mol)
    mol_j = nx_to_json(G)
    return mol_j

