import networkx as nx
import rdkit
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import (mol_to_nx,
                                                                  nx_to_mol)
from pepseq.Peptide.utils.chemistry.MonomerConnector import find_R


def prepare_ter_G(smiles: str, ResID: int) -> nx.classes.graph.Graph:
    """
    Prepare a graph representation of a molecule with capped termini.

    Args:
        smiles (str): The SMILES string representation of the molecule.
        ResID (int): The residue ID.

    Returns:
        nx.classes.graph.Graph: The graph representation of the molecule with capped termini.
    """
    mol = rdkit.Chem.MolFromSmiles(smiles)
    G = mol_to_nx(mol)
    nx.set_node_attributes(G, ResID, "ResID")
    return G


def relabel_to_str(G: nx.classes.graph.Graph) -> nx.classes.graph.Graph:
    """
    Relabels the nodes of a graph to their string representation.

    Args:
        G (nx.classes.graph.Graph): The input graph.

    Returns:
        nx.classes.graph.Graph: The graph with relabeled nodes.
    """
    mapping = {}

    for node in G:
        if type(node) == int:
            mapping[node] = str(node)

    G = nx.relabel_nodes(G, mapping)
    return G


def find_max_ResID(G: nx.classes.graph.Graph) -> int:
    """
    Finds the maximum ResID in the given graph.

    Parameters:
    - G (nx.classes.graph.Graph): The graph to search for ResIDs.

    Returns:
    - int: The maximum ResID found in the graph.
    """
    ResIDs = nx.get_node_attributes(G, "ResID").values()
    unique_ResIDs = set(ResIDs) - set(["N_terminus", "C_terminus"])
    max_ResID = max([int(i) for i in unique_ResIDs])
    ResID = str(max_ResID)
    return ResID


def cap_terminus(mol: rdkit.Chem.rdchem.Mol, terminus: str=None, smiles_building_blocks_db: int=None,
                  TerminusResID: int=None, ResID: int=1, terminus_smiles:str=None) -> rdkit.Chem.rdchem.Mol:
    """
    Caps the terminus of a molecule with a specified building block.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The input molecule.
        terminus (str, optional): The type of terminus to cap. Defaults to None.
        smiles_building_blocks_db (int, optional): The database of building blocks. Defaults to None.
        TerminusResID (int, optional): The residue ID of the terminus. Defaults to None.
        ResID (int, optional): The residue ID of the molecule. Defaults to 1.
        terminus_smiles (str, optional): The SMILES representation of the terminus building block. Defaults to None.

    Returns:
        rdkit.Chem.rdchem.Mol: The capped molecule.
    """
    mol_G = mol_to_nx(mol)
    mol_G = relabel_to_str(mol_G)

    if ResID == -1:
        ResID = find_max_ResID(mol_G)
        mol_r_id = 2
    else:
        mol_r_id = 1

    mol_R1 = str(find_R(mol_G, ResID=str(ResID), r_id=mol_r_id))

    mol_atom = list(mol_G.neighbors(mol_R1))[0]
    if terminus_smiles is None:
        terminus_smiles = smiles_building_blocks_db[terminus]
    terminus_G = prepare_ter_G(terminus_smiles, ResID=TerminusResID)

    terminus_R1 = find_R(terminus_G, ResID=TerminusResID, r_id=1)

    terminus_atom = list(terminus_G.neighbors(terminus_R1))[0]

    union_terminus_atom = "%s_%s" % (TerminusResID, terminus_atom)
    union_terminus_R1 = "%s_%s" % (TerminusResID, terminus_R1)

    G_union = nx.union(mol_G, terminus_G, rename=("", "%s_" % TerminusResID))

    G_union.add_edge(
        mol_atom, union_terminus_atom, bond_type=rdkit.Chem.rdchem.BondType.SINGLE
    )

    for node in (mol_R1, union_terminus_R1):
        G_union.remove_node(node)
    return nx_to_mol(G_union)


def cap_N_terminus(mol: rdkit.Chem.rdchem.Mol, terminus: str=None,
                smiles_building_blocks_db: dict=None, terminus_smiles: str=None) -> rdkit.Chem.rdchem.Mol:
    """
    Caps the N-terminus of a molecule with a specified terminus.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The molecule to be modified.
        terminus (str, optional): The terminus to be added. Defaults to None.
        smiles_building_blocks_db (dict, optional): A dictionary of building blocks in SMILES format. Defaults to None.
        terminus_smiles (str, optional): The SMILES representation of the terminus. Defaults to None.

    Returns:
        rdkit.Chem.rdchem.Mol: The modified molecule with the capped N-terminus.
    """
    return cap_terminus(
        mol, terminus, smiles_building_blocks_db, TerminusResID="N_terminus", ResID=1,
        terminus_smiles=terminus_smiles
    )


def cap_C_terminus(mol: rdkit.Chem.rdchem.Mol, terminus: str=None,
                smiles_building_blocks_db=None, terminus_smiles=None) -> rdkit.Chem.rdchem.Mol:
    """
    Caps the C-terminus of a molecule with a specified terminus.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The molecule to be modified.
        terminus (str, optional): The terminus to be added. Defaults to None.
        smiles_building_blocks_db: The database of building blocks. Defaults to None.
        terminus_smiles: The SMILES representation of the terminus. Defaults to None.

    Returns:
        rdkit.Chem.rdchem.Mol: The modified molecule with the capped C-terminus.
    """
    return cap_terminus(
        mol, terminus, smiles_building_blocks_db, TerminusResID="C_terminus", ResID=-1,
        terminus_smiles=terminus_smiles
    )
