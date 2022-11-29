import networkx as nx
import numpy as np
import rdkit
import rdkit.Chem


def get_molgraph(smiles_dict):
    molnames = list(smiles_dict.keys())
    num_mols = len(molnames)
    G = nx.DiGraph()

    for i in range(num_mols):

        molname_1 = molnames[i]
        smi_1 = smiles_dict[molname_1].smiles
        mol_1 = rdkit.Chem.MolFromSmiles(smi_1)
        G.add_node(molname_1)

        for j in range(num_mols):

            molname_2 = molnames[j]
            smi_2 = smiles_dict[molname_2].smiles
            mol_2 = rdkit.Chem.MolFromSmiles(smi_2)
            val = mol_1.GetSubstructMatches(mol_2, useChirality=True)
            if (i != j) and val:
                G.add_edge(molname_2, molname_1)

    bidirectional_edges = []

    for (molname_2, molname_1) in G.edges:
        if (molname_1, molname_2) in G.edges:
            bidirectional_edges.append((molname_1, molname_2))
    for edge in bidirectional_edges:
        G.remove_edge(*edge)

    return G


def remove_longshots(G):
    edges_to_remove = []

    for edge in G.edges:
        source, target = edge
        paths = [i for i in nx.all_simple_paths(G, source=source, target=target)]
        if len(paths) > 1:
            edge_to_remove = (source, target)
            edges_to_remove.append(edge_to_remove)
    for edge in edges_to_remove:
        G.remove_edge(*edge)
    return G


def order_nodes(g):
    G = g.copy()

    """
        We want to find all occurences of given AminoAcidMolecule in the PeptideMolecule.
        Therefore we define Substructure as AminoAcidMolecule and find all occurences
        of such Substructure in PeptideMolecule.
    
        The problem we might run into is that one AminoAcidMolecule can be a Substructure
        of another AminoAcidMolecule.
        For example:
            Glycine ->(is_substructure_of)-> Alanine (and other AminoAcids for that matter)
            Alanine -> aMeAla
    
        Thus:
            It is important when running this check to run (e.g. Glycine before other AAs
            and Ala before aMeAla) so that 'bigger' AminoAidMolecules occurences 'overwrite'
            occurences of smaller amino acids contained within 'bigger' Amino Acids

        G - Directed Acyclic Graph with edges going from smaller Substructure to bigger
        AminoAcidMolecule that Substructure is contained within

        E. g. 
        Glycine -> Alanine
        Alanine -> alfa-methyl-alanine
        etc.
    
        l - list of nodes such as a bigger AminoAcidMolecule always comes after its smaller Substructure
    """

    l = []
    leaves = [x for x in G.nodes() if G.out_degree(x) == 0]  # and G.in_degree(x)==1]

    while leaves:
        leaves = [
            x for x in G.nodes() if G.out_degree(x) == 0
        ]  # and G.in_degree(x)==1]

        for leaf in leaves:
            l.append(leaf)
            G.remove_node(leaf)
    for node in G.nodes:
        l.append(node)
    l.reverse()
    return l


def get_ordered_nodes(aa_smiles_dict, move_x=True):
    G = get_molgraph(aa_smiles_dict)
    ordered_nodes = order_nodes(G)

    if move_x:
        if "X" in ordered_nodes:
            i = ordered_nodes.index("X")
            ordered_nodes = ["X"] + ordered_nodes[:i] + ordered_nodes[i + 1 :]

    return ordered_nodes
