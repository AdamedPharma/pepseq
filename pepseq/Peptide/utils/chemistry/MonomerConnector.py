from typing import TypeVar

import networkx as nx
import rdkit
import rdkit.Chem
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import mol_to_nx, nx_to_mol

AminoAcidInstance = TypeVar("AminoAcidInstance")


def smi_to_G(smiles: str) -> nx.classes.graph.Graph:
    """
    Converts a SMILES string to a nx.classes.graph.Graph object of the molecule

    :param smiles: SMILES string of a molecule
    :type smiles: str

    :return: nx.classes.graph.Graph object of the molecule
    :rtype: nx.classes.graph.Graph
    """
    mol = rdkit.Chem.MolFromSmiles(smiles)
    G = mol_to_nx(mol)
    return G


def is_R(v: dict, ResID: int, r_id: int) -> bool:
    """
    Returns True if the node is a R node with the given ResID and r_id

    :param v: Node attributes
    :type v: dict

    :param ResID: Residue ID
    :type ResID: int

    :return: True if the node is a R node with the given ResID and r_id
    :rtype: bool
    """
    return (
        (v["atomic_num"] == 0)
        and (v.get("molAtomMapNumber") == r_id)
        and (v["ResID"] == ResID)
    )


def find_R(G: nx.classes.graph.Graph, ResID: int, r_id: int) -> dict:
    """
    Returns the R node with the given ResID and r_id

    :param G: nx.classes.graph.Graph object of the molecule
    :type G: nx.classes.graph.Graph

    :param ResID: Residue ID
    :type ResID: int

    :return: The R node with the given ResID and r_id
    :rtype: dict
    """
    R = [n for n, v in G.nodes(data=True) if is_R(v, ResID, r_id)][0]
    return R


def find_N(G: nx.classes.graph.Graph, ResID: int) -> int:
    """
    Returns the N node of the residue with the given ResID

    :param G: nx.classes.graph.Graph object of the molecule
    :type G: nx.classes.graph.Graph

    :param ResID: Residue ID
    :type ResID: int

    :return: The N node of the residue with the given ResID
    :rtype: int
    """
    R = find_R(G, ResID, r_id=1)
    N = list(G.neighbors(R))[0]
    return N


def find_CO(G: nx.classes.graph.Graph, ResID: int) -> int:
    """
    Returns the CO node of the residue with the given ResID

    :param G: nx.classes.graph.Graph object of the molecule
    :type G: nx.classes.graph.Graph

    :param ResID: Residue ID
    :type ResID: int

    :return: The CO node of the residue with the given ResID
    :rtype: int
    """
    R = find_R(G, ResID, r_id=2)
    N = list(G.neighbors(R))[0]
    return N


def merge_graph(G: nx.classes.graph.Graph, ResID=1) -> nx.classes.graph.Graph:
    """
    Merges the residue with the given ResID with the next residue

    :param G: nx.classes.graph.Graph object of the molecule
    :type G: nx.classes.graph.Graph

    :param ResID: Residue ID
    :type ResID: int

    :return: nx.classes.graph.Graph object of the molecule with the residue with
     the given ResID merged with the next residue
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
    Returns a list of nx.classes.graph.Graph objects of the residues with the given residue symbols

    :param residue_symbols: List of residue symbols
    :type residue_symbols: list

    :param smiles_building_blocks_db: Dictionary of residue symbols to SMILES strings
    :type smiles_building_blocks_db: dict

    :return: List of nx.classes.graph.Graph objects of the residues with the given residue symbols
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
    Merges the list of residue graphs into a single peptide graph

    :param graphs: List of nx.classes.graph.Graph objects of the residues
    :type graphs: list

    :return: nx.classes.graph.Graph object of the peptide
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
    Returns a RDKit molecule object of the peptide with the given list of residue symbols

    :param residue_symbols: List of residue symbols
    :type residue_symbols: list

    :param smiles_building_blocks_db: Dictionary of residue symbols to SMILES strings
    :type smiles_building_blocks_db: dict

    :return: RDKit molecule object of the peptide with the given list of residue symbols
    :rtype: rdkit.Chem.rdchem.Mol
    """
    residue_graphs = get_residues_Gs(residue_symbols, smiles_building_blocks_db)
    peptide_graph = merge_residue_graphs(residue_graphs)
    return nx_to_mol(peptide_graph)
