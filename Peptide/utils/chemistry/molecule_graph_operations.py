import networkx as nx
import rdkit


def add_dummyAtom_to_G(
    G: nx.classes.graph.Graph, idx: int, isotope: int
) -> nx.classes.graph.Graph:
    name = "R%d" % (isotope)
    """
    Input:

        G - molecule graph
        idx - id of atom to connect dummyAtom to
        isotope - id of attachment point, e.g. 1 for [1*], 2 for [2*]

    Output:
        G - molecule graph with added dummyAtom

    Process:
        dummyAtom node is produced and edge is created from idx atom
        to dummyAtom
    """

    G.add_node(
        name,
        atomic_num=0,
        formal_charge=0,
        chiral_tag=rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
        hybridization=rdkit.Chem.rdchem.HybridizationType.SP3,
        num_explicit_hs=0,
        is_aromatic=False,
        isotope=isotope,
    )

    G.add_edge(idx, name, bond_type=rdkit.Chem.rdchem.BondType.SINGLE)
    return G


def merge_combo_graph(G: nx.classes.graph.Graph) -> nx.classes.graph.Graph:
    """
    Input:
        G - molecular graph containing subgraphs for molecules to be merged
            subgraphs are not connected by edge but contain dummyAtoms
            symbolizing complementary attachment points
            on both molecules e.g. [1*] and [2*]

    Output:
        G - molecular graph of new merged molecule

    Process:
        Radicals/dummyAtoms are identified as atoms with atomic_num = 0
        (just as they are identified in rdkit)

        attachment_point_IDs ([1*] or [2*] etc. ) on each of the molecules are identified from Isotope value
        (just as they are identified in rdkit)

        for each attachment_point_ID of DummyAtoms that are present on both molecules we
            identify the NEIGHBOURS of each point on each of the molecules
            and then connect the NEIGHBOURS with SingleBond

        Then the dummy atoms are removed

    Output:
        The merged Graph

    """

    radicals = []
    for node_id in G.nodes:
        node = G.nodes[node_id]
        atomic_num = node["atomic_num"]
        is_radical = atomic_num == 0
        if is_radical:
            attachment_point_id = G.nodes[node_id]["isotope"]
            tup = (
                node_id,
                attachment_point_id,
            )
            radicals.append(tup)

    #

    connections = {}

    for idx, radical_id in radicals:
        if connections.get(radical_id) is None:
            connections[radical_id] = []
        connections[radical_id].append(idx)

    to_remove = []

    for radical_id in connections:
        bond = sorted(connections[radical_id])
        to_remove += bond

        if len(bond) == 2:
            neighbour1 = list(G.neighbors(bond[0]))[0]
            neighbour2 = list(G.neighbors(bond[1]))[0]
            G.add_edge(
                neighbour1, neighbour2, bond_type=rdkit.Chem.rdchem.BondType.SINGLE
            )
    for radical in to_remove:
        G.remove_node(radical)
    return G
