import networkx as nx
import rdkit


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


def do_all(smiles: str, validate=False) -> nx.classes.graph.Graph:
    """
    reads SMILES string into
        nx.classes.graph.Graph graph
        with optional validation
    """
    mol = rdkit.Chem.MolFromSmiles(smiles.strip())
    can_smi = rdkit.Chem.MolToSmiles(mol)
    G = mol_to_nx(mol)
    if validate:
        mol = nx_to_mol(G)
        new_smi = rdkit.Chem.MolToSmiles(mol)
        assert new_smi == smiles
    return G
