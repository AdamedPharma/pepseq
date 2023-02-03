import networkx as nx
import rdkit

from Peptide.utils.chemistry.cap_termini import cap_C_terminus, cap_N_terminus
from Peptide.utils.chemistry.mol_to_nx_translation import mol_to_nx, nx_to_mol
from Peptide.utils.chemistry.MonomerConnector import (
    get_molecule_from_list_of_residue_symbols,
)
from Peptide.utils.Parser import find_termini, get_canonical, parse_canonical


def add_internal_bond(G, res1_id, atom_name_1, res2_id, atom_name_2):
    Cys1_SG = [
        n
        for n, v in G.nodes(data=True)
        if (v.get("AtomName") == atom_name_1 and v.get("ResID") == str(res1_id))
    ][0]
    Cys2_SG = [
        n
        for n, v in G.nodes(data=True)
        if (v.get("AtomName") == atom_name_2 and v.get("ResID") == str(res2_id))
    ][0]
    G.add_edge(Cys1_SG, Cys2_SG, bond_type=rdkit.Chem.rdchem.BondType.SINGLE)
    return G


def add_disulfide_bond(G, res1_id, res2_id):
    return add_internal_bond(G, res1_id, "SG", res2_id, "SG")


def get_attachment_points(staple_graph):
    dummyAtoms = [
        n for n, v in staple_graph.nodes(data=True) if v.get("atomic_num") == 0
    ]
    attachment_points_on_staple_dict = {}
    for dummyAtom in dummyAtoms:

        attachment_point = list(staple_graph.neighbors(dummyAtom))[0]
        attachment_point_id = staple_graph.nodes[dummyAtom]["isotope"]

        attachment_points_on_staple_dict[attachment_point_id] = attachment_point

    for dummyAtom in dummyAtoms:
        staple_graph.remove_node(dummyAtom)
    return staple_graph, attachment_points_on_staple_dict


def find_atom(G, ResID, AtomName):
    ResIDs = nx.get_node_attributes(G, "ResID")
    AtomNames = nx.get_node_attributes(G, "AtomName")
    SGs = [i for i in AtomNames if AtomNames[i] == AtomName]
    point_list = [SG for SG in SGs if ResIDs[SG] == ResID]
    return point_list.pop()


def add_staple(
    peptide_graph, staple_graph, peptide_attachment_points, prefix="staple_1_"
):
    staple_graph, attachment_points_on_staple_dict = get_attachment_points(staple_graph)

    G_stapled_peptide_union = nx.union(peptide_graph, staple_graph, rename=("", prefix))

    for attachment_point_id in attachment_points_on_staple_dict:
        j_peptide_attachment_point = peptide_attachment_points[attachment_point_id]

        staple_attachment_point = attachment_points_on_staple_dict[attachment_point_id]

        ResID = j_peptide_attachment_point["ResID"]

        AtomName = j_peptide_attachment_point["AtomName"]

        peptide_attachment_point = find_atom(G_stapled_peptide_union, ResID, AtomName)

        G_stapled_peptide_union.add_edge(
            peptide_attachment_point,
            "%s%s" % (prefix, str(staple_attachment_point)),
            bond_type=rdkit.Chem.rdchem.BondType.SINGLE,
        )
    return G_stapled_peptide_union


def get_peptide_json_from_sequence(sequence, db_json):
    # canonical_sequence = get_canonical(sequence, db_json)
    # symbols_list_w_termini = parse_canonical(canonical_sequence)

    N_terminus, C_terminus, sequence_str_wo_termini = find_termini(sequence, db_json)

    peptide_json = {
        "sequence": sequence_str_wo_termini,
        "internal_modifications": [],
        "external_modifications": [],
        "N_terminus": N_terminus,
        "C_terminus": C_terminus,
    }

    return peptide_json


def get_molecule_from_sequence(sequence, db_json, N_terminus=None, C_terminus=None):
    canonical_sequence = get_canonical(sequence, db_json)
    symbols_list_w_termini = parse_canonical(canonical_sequence)
    if N_terminus is None:
        N_terminus = symbols_list_w_termini[0]
    if C_terminus is None:
        C_terminus = symbols_list_w_termini[-1]
    residue_symbols = symbols_list_w_termini[1:-1]

    smiles_building_blocks_db = {}

    for residue_symbol in residue_symbols:
        smiles_building_blocks_db[residue_symbol] = db_json["smiles"]["aa"].get(
            residue_symbol
        )["smiles_radical"]

    smiles_building_blocks_db[N_terminus] = db_json["smiles"]["n_terms"].get(
        N_terminus
    )["smiles_radical"]

    smiles_building_blocks_db[C_terminus] = db_json["smiles"]["c_terms"][C_terminus][
        "smiles_radical"
    ]

    mol = get_molecule_from_list_of_residue_symbols(
        residue_symbols, smiles_building_blocks_db
    )

    mol_w_n_terminus = cap_N_terminus(
        mol, terminus=N_terminus, smiles_building_blocks_db=smiles_building_blocks_db
    )

    mol_w_nc_terminus = cap_C_terminus(
        mol_w_n_terminus,
        terminus=C_terminus,
        smiles_building_blocks_db=smiles_building_blocks_db,
    )

    return mol_w_nc_terminus


def get_smiles_from_sequence(sequence, db_json):
    mol = get_molecule_from_sequence(
        sequence, db_json, N_terminus=None, C_terminus=None
    )
    smiles = rdkit.Chem.MolToSmiles(mol)
    return smiles


def get_molecule_from_json(j, db_json):
    sequence = j["sequence"]
    N_terminus = j.get("N_terminus")
    C_terminus = j.get("C_terminus")

    return get_molecule_from_sequence(sequence, db_json, N_terminus, C_terminus)


class BuildingModifiedPeptideFromPeptideJSON(object):
    """
    sequence is constructed (automatically or through connecting CXSMILES from database)
    Thus atoms serving as AttachmentPoints are marked
    For each ExternalModification edge is added between SequenceAttachmentPoint
    and ExternalModificationAttachmentPoint

    For each InternalModification edge is added between SequenceAttachmentPoint
    and matching SequenceAttachmentPoint
    """

    def __init__(self):
        return

    def execute(self, peptide_json, db_json):

        peptide_mol = get_molecule_from_json(peptide_json, db_json)

        peptide_graph = mol_to_nx(peptide_mol)

        internal_modifications = peptide_json["internal_modifications"]

        ext_mod_id = 0

        for ext_mod in peptide_json["external_modifications"]:
            ext_mod_id += 1
            staple_smi = ext_mod["smiles"]
            staple_mol = rdkit.Chem.MolFromSmiles(staple_smi)
            staple_graph = mol_to_nx(staple_mol)
            attachment_points_on_sequence = ext_mod["attachment_points_on_sequence"]

            prefix = "staple_%d_" % ext_mod_id

            peptide_graph = add_staple(
                peptide_graph,
                staple_graph.copy(),
                attachment_points_on_sequence,
                prefix,
            )

        for internal_modification_id in internal_modifications:
            internal_bond = internal_modifications[internal_modification_id]

            source, target = internal_bond

            peptide_graph = add_internal_bond(
                peptide_graph,
                source["ResID"],
                source["AtomName"],
                target["ResID"],
                target["AtomName"],
            )
        return nx_to_mol(peptide_graph)


def get_smiles_from_peptide_json(peptide_json, db_json):
    mol = BuildingModifiedPeptideFromPeptideJSON.execute(peptide_json, db_json)
    smiles = rdkit.Chem.MolToSmiles(mol)
    return smiles
