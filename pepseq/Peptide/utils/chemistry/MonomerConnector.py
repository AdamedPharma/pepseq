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

    Args:
        smiles (str): The SMILES string to convert.

    Returns:
        nx.classes.graph.Graph: The networkx graph representation of the molecule.
    """
    mol = rdkit.Chem.MolFromSmiles(smiles)
    G = mol_to_nx(mol)
    return G


def is_R(v: dict, ResID: int, r_id: int) -> bool:
    """
    Check if the given dictionary represents a monomer with the specified ResID and r_id.

    Args:
        v (dict): The dictionary representing the monomer.
        ResID (int): The ResID to check.
        r_id (int): The r_id to check.

    Returns:
        bool: True if the dictionary represents the specified monomer, False otherwise.
    """
    return (v["atomic_num"] == 0) and (v["isotope"] == r_id) and (v["ResID"] == ResID)


def find_R(G: nx.classes.graph.Graph, ResID: int, r_id: int) -> dict:
    """
    Finds and returns the R node in the graph G that matches the given ResID and r_id.

    Parameters:
    - G: The graph to search in.
    - ResID: The ResID to match.
    - r_id: The r_id to match.

    Returns:
    - R: The R node that matches the given ResID and r_id.
    """
    R = [n for n, v in G.nodes(data=True) if is_R(v, ResID, r_id)][0]
    return R


def find_N(G: nx.classes.graph.Graph, ResID: int) -> int:
    """
    Finds the N atom connected to the given residue ID in the graph.

    Parameters:
    - G: nx.classes.graph.Graph
        The graph representing the molecule.
    - ResID: int
        The ID of the residue.

    Returns:
    - int
        The ID of the N atom connected to the given residue.
    """
    R = find_R(G, ResID, r_id=1)
    N = list(G.neighbors(R))[0]
    return N


def find_CO(G: nx.classes.graph.Graph, ResID: int) -> int:
    """
    Finds the carbon-oxygen (CO) atom in a peptide graph.

    Parameters:
    - G (nx.classes.graph.Graph): The peptide graph.
    - ResID (int): The ID of the residue.

    Returns:
    - int: The ID of the carbon-oxygen (CO) atom.
    """
    R = find_R(G, ResID, r_id=2)
    N = list(G.neighbors(R))[0]
    return N


def merge_graph(G: nx.classes.graph.Graph, ResID=1) -> nx.classes.graph.Graph:
    """
    Merges two monomers in a graph by connecting the CO atom of the first monomer
    with the N atom of the second monomer.

    Args:
        G (nx.classes.graph.Graph): The graph representing the monomers.
        ResID (int, optional): The ID of the first monomer. Defaults to 1.

    Returns:
        nx.classes.graph.Graph: The merged graph.
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

    Args:
        residue_symbols (list): A list of residue symbols.
        smiles_building_blocks_db (dict): A dictionary mapping residue symbols to SMILES building blocks.

    Returns:
        list: A list of graph representations (Gs) of residues.
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

    Args:
        graphs (list): A list of residue graphs.

    Returns:
        nx.classes.graph.Graph: The merged peptide graph.
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

    Args:
        residue_symbols (list): A list of residue symbols.
        smiles_building_blocks_db: The database of smiles building blocks.

    Returns:
        rdkit.Chem.rdchem.Mol: The resulting molecule.
    """
    residue_graphs = get_residues_Gs(residue_symbols, smiles_building_blocks_db)
    peptide_graph = merge_residue_graphs(residue_graphs)
    return nx_to_mol(peptide_graph)
