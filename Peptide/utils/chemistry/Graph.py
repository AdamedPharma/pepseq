import networkx as nx
import rdkit
import rdkit.Chem
from Peptide.utils.chemistry.molecule_graph_generation import mol_to_nx, nx_to_mol
from Peptide.utils.chemistry.molecule_graph_operations import (
    add_dummyAtom_to_G,
    merge_combo_graph,
)
from Peptide.utils.chemistry.smiles_to_seq import (
    break_residues,
    get_N_CA_CO,
    longest_backbone,
)


def connect_mols_on_dummyAtoms(
    mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol
) -> rdkit.Chem.rdchem.Mol:
    """
    Input:
        mol1 - first molecule (contains radical dummyAtom e.g. [1*] to be connected 'with'
                   matching dummyAtom on second molecule)
        mol2 - second molecule (contains radical dummyAtom e.g. [1*] to be connected 'with'
                   matching dummyAtom on first molecule)
    Output:
        mol_union - molecule product of merging mol1 and mol2 on matching dummyAtoms

    """
    G1 = mol_to_nx(mol1)
    G2 = mol_to_nx(mol2)

    G_union = nx.union(G1, G2, rename=("mol1_", "mol2_"))
    G_union = merge_combo_graph(G_union)
    mol_union = nx_to_mol(G_union)
    return mol_union


atomic_numbers = {"S": 16}


def get_atom_by_element(G: nx.classes.graph.Graph, atomic_num: int) -> int:
    atom_id = [i for i in G.nodes if G.nodes[i]["atomic_num"] == atomic_num][0]
    return atom_id


def find_atom_id(mol: rdkit.Chem.rdchem.Mol, res_id: int, atom_label: str) -> int:
    """ """
    G = mol_to_nx(mol)
    backbone_atoms = longest_backbone(mol)
    n_ca_co = get_N_CA_CO(G, backbone_atoms, rdkit.Chem.rdchem.BondType.DOUBLE)
    residues = break_residues(G, n_ca_co)
    res = residues[res_id - 1]
    atom_id = get_atom_by_element(res, atomic_numbers[atom_label])
    return atom_id


def add_dummyAtom(
    mol: rdkit.Chem.rdchem.Mol, atom_id: int, isotope: int
) -> rdkit.Chem.rdchem.Mol:
    G = mol_to_nx(mol)

    G = add_dummyAtom_to_G(G, atom_id, isotope)

    return nx_to_mol(G)


def connect_on_radical(
    mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol
) -> rdkit.Chem.rdchem.Mol:
    G1 = mol_to_nx(mol1)
    G2 = mol_to_nx(mol2)
    G_union = nx.union(G1, G2, rename=("mol1_", "mol2_"))
    G_union = merge_combo_graph(G_union)
    return nx_to_mol(G_union)


def reconstruct_molecule_from_json(mol_json: dict) -> rdkit.Chem.rdchem.Mol:
    """

    now to reconstruct molecule from JSON

    mol_json = {
        'attachment_points': {
            1: (2, 'C', 'S'),
            2: (4, 'C', 'S')
            },

        'sequence': 'GCGCG',

        'modification': '[1*]SCNOCN[2*]'
        }

    Output:
        mol - reconstructed molecule

    Process:
        1. rdkit Molecule is generated from Sequence directly or through FASTA
        2. For each attachment point:
            dummyAtom is Added
            2.1 To this end we need to find AtomId of AttachmentPoint

    TO DO:
        Include more sophisticated Sequences
    """
    sequence = mol_json["sequence"]

    seq_mol = rdkit.Chem.MolFromSequence(sequence)

    attachment_points = mol_json["attachment_points"]

    for dummyAtomLabel in attachment_points:
        res_id, res_name, atom_label = attachment_points[dummyAtomLabel]
        atom_id = find_atom_id(seq_mol, res_id, atom_label)
        seq_mol = add_dummyAtom(seq_mol, atom_id, dummyAtomLabel)

    mod_smi = mol_json["modification"]
    mod_mol = rdkit.Chem.MolFromSmiles(mod_smi)
    mol_reconstructed = connect_on_radical(seq_mol, mod_mol)
    return mol_reconstructed
