from typing import TypeVar

import networkx as nx
import rdkit
import rdkit.Chem
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import (mol_to_nx,
                                                                  nx_to_mol)

AminoAcidInstance = TypeVar("AminoAcidInstance")


def smi_to_G(smiles: str) -> nx.classes.graph.Graph:
    """
    Converts a SMILES string to a networkx graph representation.

    :param smiles: The SMILES string to convert.
    :type smiles: str

    :return: The networkx graph representation of the molecule.
    :rtype: nx.classes.graph.Graph
    """
    mol = rdkit.Chem.MolFromSmiles(smiles)
    G = mol_to_nx(mol)
    return G


def is_R(v: dict, ResID: int, r_id: int) -> bool:
    """
    Check if the given dictionary represents a monomer with the specified ResID and r_id.

    :param v: The dictionary representing the monomer.
    :type v: dict
    :param ResID: The ResID to check.
    :type ResID: int
    :param r_id: The r_id to check.
    :type r_id: int

    :return: True if the dictionary represents the specified monomer, False otherwise.
    :rtype: bool
    """
    return (v["atomic_num"] == 0) and (v["isotope"] == r_id) and (v["ResID"] == ResID)


def find_R(G: nx.classes.graph.Graph, ResID: int, r_id: int) -> dict:
    """
    Finds and returns the R node in the graph G that matches the given ResID and r_id.

    :param G: The graph to search in.
    :type G: nx.classes.graph.Graph
    :param ResID: The ResID to match.
    :type ResID: int
    :param r_id: The r_id to match.
    :type r_id: int

    :return: The R node that matches the given ResID and r_id.
    :rtype: dict
    """
    R = [n for n, v in G.nodes(data=True) if is_R(v, ResID, r_id)][0]
    return R


def find_N(G: nx.classes.graph.Graph, ResID: int) -> int:
    """
    Finds the N atom connected to the given residue ID in the graph.

    :param G: The graph representing the molecule.
    :type G: nx.classes.graph.Graph
    :param ResID: The ID of the residue.
    :type ResID: int

    :return: The ID of the N atom connected to the given residue.
    :rtype: int
    """
    R = find_R(G, ResID, r_id=1)
    N = list(G.neighbors(R))[0]
    return N


def find_CO(G: nx.classes.graph.Graph, ResID: int) -> int:
    """
    Finds the carbon-oxygen (CO) atom in a peptide graph.

    Parameters:
    :param G: The peptide graph.
    :type G: nx.classes.graph.Graph
    :param ResID: The ID of the residue.
    :type ResID: int

    :return: The ID of the carbon-oxygen (CO) atom.
    :rtype: int
    """
    R = find_R(G, ResID, r_id=2)
    N = list(G.neighbors(R))[0]
    return N


def merge_graph(G: nx.classes.graph.Graph, ResID=1) -> nx.classes.graph.Graph:
    """
    Merges two monomers in a graph by connecting the CO atom of the first monomer
    with the N atom of the second monomer.

    :param G: The graph representing the monomers.
    :type G: nx.classes.graph.Graph
    :param ResID: The ID of the first monomer. Defaults to 1.
    :type ResID: int, optional

    :return: The merged graph.
    :rtype: nx.classes.graph.Graph
    """
    ResNextID = ResID + 1

    CO = find_CO(G, ResID)
    N = find_N(G, ResNextID)

    Res1_R2 = find_R(G, ResID, r_id=2)
    Res2_R1 = find_R(G, ResNextID, r_id=1)

    G.add_edge(CO, N, bond_type=rdkit.Chem.rdchem.BondType.SINGLE)
    for node in (Res1_R2, Res2_R1):
        G.remove_node(node)
    return G


def get_residues_Gs(residue_symbols: list, smiles_building_blocks_db: dict) -> list:
    """
    Get the graph representations (Gs) of residues based on their symbols and a database of SMILES building blocks.

    :param residue_symbols: A list of residue symbols.
    :type residue_symbols: list
    :param smiles_building_blocks_db: A dictionary mapping residue symbols to SMILES building blocks.
    :type smiles_building_blocks_db: dict

    :return:  A list of graph representations (Gs) of residues.
    :rtype: list
    """
    Residue_Gs = []

    for i in range(len(residue_symbols)):
        res_symbol = residue_symbols[i]
        res_smiles = smiles_building_blocks_db[res_symbol]
        res_G = smi_to_G(res_smiles)

        nx.set_node_attributes(res_G, i + 1, name="ResID")

        Residue_Gs.append(res_G)
    return Residue_Gs


def merge_residue_graphs(graphs: list) -> nx.classes.graph.Graph:
    """
    Merges a list of residue graphs into a single peptide graph.

    :param graphs: A list of residue graphs.
    :type graphs: list

    :return: The merged peptide graph.
    :rtype: nx.classes.graph.Graph
    """

    first_residue_graph = graphs[0]
    peptide_graph = nx.union(
        first_residue_graph,
        graphs[1],
        rename=("Res%d_" % (1), "Res%d_" % (2)),
    )
    peptide_graph = merge_graph(peptide_graph, ResID=1)

    for i in range(1, len(graphs) - 1):
        peptide_graph = nx.union(
            peptide_graph, graphs[i + 1], rename=("", "Res%d_" % (i + 2))
        )
        peptide_graph = merge_graph(peptide_graph, ResID=(i + 1))
    return peptide_graph


def get_molecule_from_list_of_residue_symbols(
    residue_symbols: list, smiles_building_blocks_db
) -> rdkit.Chem.rdchem.Mol:
    """
    Converts a list of residue symbols into a molecule using the provided smiles_building_blocks_db.

    :param residue_symbols: A list of residue symbols.
    :type residue_symbols: list
    :param smiles_building_blocks_db: The database of smiles building blocks.
    :type smiles_building_blocks_db: dict

    :return: The resulting molecule.
    :rtype: rdkit.Chem.rdchem.Mol
    """
    residue_graphs = get_residues_Gs(residue_symbols, smiles_building_blocks_db)
    peptide_graph = merge_residue_graphs(residue_graphs)
    return nx_to_mol(peptide_graph)
