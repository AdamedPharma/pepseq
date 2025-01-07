import networkx as nx
import rdkit
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import mol_to_nx, nx_to_mol
from pepseq.Peptide.utils.chemistry.MonomerConnector import find_R


def prepare_ter_G(smiles: str, ResID: int) -> nx.classes.graph.Graph:
    mol = rdkit.Chem.MolFromSmiles(smiles)
    G = mol_to_nx(mol)
    nx.set_node_attributes(G, ResID, "ResID")
    return G


def relabel_to_str(G: nx.classes.graph.Graph) -> nx.classes.graph.Graph:
    mapping = {}

    for node in G:
        if type(node) == int:
            mapping[node] = str(node)

    G = nx.relabel_nodes(G, mapping)
    return G


def find_max_ResID(G: nx.classes.graph.Graph) -> int:
    ResIDs = nx.get_node_attributes(G, "ResID").values()
    unique_ResIDs = set(ResIDs) - set(["N_terminus", "C_terminus"])
    max_ResID = max([int(i) for i in unique_ResIDs])
    ResID = str(max_ResID)
    return ResID


def cap_terminus(
    mol: rdkit.Chem.rdchem.Mol,
    terminus: str = None,
    smiles_building_blocks_db: int = None,
    TerminusResID: int = None,
    ResID: int = 1,
    terminus_smiles: str = None,
) -> rdkit.Chem.rdchem.Mol:
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


def cap_N_terminus(
    mol: rdkit.Chem.rdchem.Mol,
    terminus: str = None,
    smiles_building_blocks_db: dict = None,
    terminus_smiles: str = None,
) -> rdkit.Chem.rdchem.Mol:
    return cap_terminus(
        mol,
        terminus,
        smiles_building_blocks_db,
        TerminusResID="N_terminus",
        ResID=1,
        terminus_smiles=terminus_smiles,
    )


def cap_C_terminus(
    mol: rdkit.Chem.rdchem.Mol,
    terminus: str = None,
    smiles_building_blocks_db=None,
    terminus_smiles=None,
) -> rdkit.Chem.rdchem.Mol:
    return cap_terminus(
        mol,
        terminus,
        smiles_building_blocks_db,
        TerminusResID="C_terminus",
        ResID=-1,
        terminus_smiles=terminus_smiles,
    )
