from collections.abc import Sequence
from typing import TypeVar

import networkx as nx
import rdkit
import rdkit.Chem
from Peptide.utils.chemistry.mol_to_nx_translation import mol_to_nx, nx_to_mol

AminoAcidInstance = TypeVar("AminoAcidInstance")


def smi_to_G(smiles):
    mol = rdkit.Chem.MolFromSmiles(smiles)
    G = mol_to_nx(mol)
    return G


def is_R(v, ResID, r_id):
    return (v["atomic_num"] == 0) and (v["isotope"] == r_id) and (v["ResID"] == ResID)


def find_R(G, ResID, r_id):
    R = [n for n, v in G.nodes(data=True) if is_R(v, ResID, r_id)][0]
    return R


def find_N(G, ResID):
    R = find_R(G, ResID, r_id=1)
    N = list(G.neighbors(R))[0]
    return N


def find_CO(G, ResID):
    R = find_R(G, ResID, r_id=2)
    N = list(G.neighbors(R))[0]
    return N


def merge_graph(G, ResID=1):
    ResNextID = ResID + 1

    CO = find_CO(G, ResID)
    N = find_N(G, ResNextID)

    Res1_R2 = find_R(G, ResID, r_id=2)
    Res2_R1 = find_R(G, ResNextID, r_id=1)

    G.add_edge(CO, N, bond_type=rdkit.Chem.rdchem.BondType.SINGLE)
    for node in (Res1_R2, Res2_R1):
        G.remove_node(node)
    return G


def get_residues_Gs(residue_symbols, smiles_building_blocks_db):
    Residue_Gs = []

    for i in range(len(residue_symbols)):
        res_symbol = residue_symbols[i]
        res_smiles = smiles_building_blocks_db[res_symbol]
        res_G = smi_to_G(res_smiles)

        nx.set_node_attributes(res_G, i + 1, name="ResID")

        Residue_Gs.append(res_G)
    return Residue_Gs


def merge_residue_graphs(graphs):
    i = 0
    first_residue_graph = graphs[i]
    peptide_graph = nx.union(
        first_residue_graph,
        graphs[i + 1],
        rename=("Res%d_" % (i + 1), "Res%d_" % (i + 2)),
    )
    peptide_graph = merge_graph(peptide_graph, ResID=(i + 1))

    for i in range(1, len(graphs) - 1):
        peptide_graph = nx.union(
            peptide_graph, graphs[i + 1], rename=("", "Res%d_" % (i + 2))
        )
        peptide_graph = merge_graph(peptide_graph, ResID=(i + 1))
    return peptide_graph


def get_molecule_from_list_of_residue_symbols(
    residue_symbols, smiles_building_blocks_db
):
    residue_graphs = get_residues_Gs(residue_symbols, smiles_building_blocks_db)
    peptide_graph = merge_residue_graphs(residue_graphs)
    return nx_to_mol(peptide_graph)
