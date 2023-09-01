import networkx as nx
import rdkit
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import (mol_to_nx,
                                                                  nx_to_mol)


def match_to_res_id(mol: rdkit.Chem.rdchem.Mol, ResID: str, cx_smarts_db: dict) -> tuple:
    G = mol_to_nx(mol)
    max_cover = 0
    max_aa = None
    max_match = None
    max_aa_mol = None
    for aa in cx_smarts_db:
        aa_smarts = cx_smarts_db[aa]
        aa_mol = rdkit.Chem.MolFromSmarts(aa_smarts)
        matches = mol.GetSubstructMatches(aa_mol)

        for match in matches:
            match_subgraph = G.subgraph(match)
            match_res_ids = set(
                nx.get_node_attributes(match_subgraph, "ResID").values()
            )
            if (len(match_res_ids) == 1) and list(match_res_ids)[0] == ResID:
                cover = len(match)
                if cover > max_cover:
                    max_cover = cover
                    max_aa = aa
                    max_match = match
                    max_aa_mol = aa_mol
    return (max_aa, max_aa_mol, max_match)


def get_res_matches(mol: rdkit.Chem.rdchem.Mol, cx_smarts_db: dict) -> dict:
    res_matches = {}
    G = mol_to_nx(mol)
    ResIDs_by_atom = nx.get_node_attributes(G, "ResID")
    ResIDs = sorted(list(set(ResIDs_by_atom.values())))

    for ResID in ResIDs:
        (max_aa, max_aa_mol, max_match) = match_to_res_id(mol, ResID, cx_smarts_db)
        res_matches[ResID] = (max_aa, max_aa_mol, max_match)
    return res_matches


def propagate_matches(mol: rdkit.Chem.rdchem.Mol, res_matches: dict) -> rdkit.Chem.rdchem.Mol:
    G = mol_to_nx(mol)
    for ResID in res_matches:
        aa, aa_mol, match = res_matches[ResID]

        for i in range(len(match)):
            atom_match = aa_mol.GetAtomWithIdx(i)
            atom_id = match[i]
            G.nodes[atom_id]["ResID"] = ResID
            G.nodes[atom_id]["ResName"] = aa

            if "AtomName" in atom_match.GetPropNames():
                AtomName = atom_match.GetProp("AtomName")
                G.nodes[atom_id]["AtomName"] = AtomName
    return nx_to_mol(G)


def get_connecting_edges(union_graph: nx.classes.graph.Graph, H_graph: nx.classes.graph.Graph,
                        I_graph: nx.classes.graph.Graph) -> list:
    out = [
        edge
        for edge in union_graph.edges
        if (edge[0] in H_graph and edge[1] in I_graph)
        or (edge[0] in I_graph and edge[1] in H_graph)
    ]
    return out


def process_res_res_connection(res_atoms_1: set, res_atoms_2: set, connecting_edges: list) -> list:
    attachment_point_pairs = []
    for connecting_edge in connecting_edges:
        res_1_attachment_point = (set(connecting_edge) & res_atoms_1).pop()
        res_2_attachment_point = (set(connecting_edge) & res_atoms_2).pop()
        tup = (res_1_attachment_point, res_2_attachment_point)
        attachment_point_pairs.append(tup)
    return attachment_point_pairs


def process_internal_connections(connections, res_matches, G: nx.classes.graph.Graph) -> list:
    internal_connections = []
    connection_id = 0

    for name1, name2, connecting_edges in connections:
        list_1 = name1.split("_")
        list_2 = name2.split("_")

        g_type_1, res_id_1 = list_1[:2]
        g_type_2, res_id_2 = list_2[:2]

        res_atoms_1 = set(res_matches[res_id_1][2])
        res_atoms_2 = set(res_matches[res_id_2][2])

        attachment_point_pairs = process_res_res_connection(
            res_atoms_1, res_atoms_2, connecting_edges
        )

        for attachment_point_1, attachment_point_2 in attachment_point_pairs:
            connection_id += 1

            AtomName_1 = G.nodes[attachment_point_1].get("AtomName")
            AtomName_2 = G.nodes[attachment_point_2].get("AtomName")
            internal_connection = {
                connection_id: [
                    {
                        "ResID": res_id_1,
                        "AtomName": AtomName_1,
                        "ResidueName": "",
                    },
                    {
                        "ResID": res_id_2,
                        "AtomName": AtomName_2,
                        "ResidueName": "",
                    },
                ]
            }

            internal_connections.append(internal_connection)
    return internal_connections


def sorted_connection(connection):
    res_ids = []

    for i in range(len(connection)):
        res_name = connection[i][0]
        res_id = int(res_name.split("_")[1])
        res_ids.append((res_id, i))

    sorted_res_ids = sorted(res_ids)
    return [connection[i] for (res_id, i) in sorted_res_ids]


def process_external_connections(connections, res_matches, modification_graphs, G: nx.classes.graph.Graph):
    attachment_point_id = 0
    external_modifications = []

    for mod_id in connections:
        attachment_points_on_seq = {}
        mod_graph = modification_graphs[int(mod_id) - 1].copy()
        mod_atoms = set(mod_graph.nodes)

        connection = connections.get(mod_id)
        connection = sorted_connection(connection)

        for res_name, connecting_edges in connection:
            list_res = res_name.split("_")

            g_type_1, res_id = list_res[:2]
            res_id = list_res[1]

            res_atoms = set(res_matches[res_id][2])

            attachment_point_pairs = process_res_res_connection(
                res_atoms, mod_atoms, connecting_edges
            )

            for res_attachment_point, mod_attachment_point in attachment_point_pairs:
                attachment_point_id += 1
                AtomName = G.nodes[res_attachment_point].get("AtomName")
                attachment_points_on_seq[attachment_point_id] = {
                    "attachment_point_id": attachment_point_id,
                    "ResID": res_id,
                    "AtomName": AtomName,
                    "ResidueName": "",
                }

                mod_graph.add_node(
                    "%d*" % attachment_point_id,
                    **{
                        "atomic_num": 0,
                        "formal_charge": 0,
                        "chiral_tag": rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
                        "hybridization": rdkit.Chem.rdchem.HybridizationType.SP3,
                        "num_explicit_hs": 0,
                        "is_aromatic": False,
                        "isotope": attachment_point_id,
                    }
                )

                mod_graph.add_edge(
                    mod_attachment_point,
                    "%d*" % attachment_point_id,
                    bond_type=rdkit.Chem.rdchem.BondType.SINGLE,
                )

                mod_mol = nx_to_mol(mod_graph)
                mod_smiles = rdkit.Chem.MolToSmiles(mod_mol)

        external_modification = {
            "smiles": mod_smiles,
            "max_attachment_point_id": attachment_point_id,
            "attachment_points_on_sequence": attachment_points_on_seq,
        }
        external_modifications.append(external_modification)
    return external_modifications


def split_connections_by_type(connections):
    res_res_connections = []
    external_mod_connections = {}

    for name1, name2, connecting_edges in connections:
        if name1.split("_")[0] == name2.split("_")[0] == "Res":
            res_res_connections.append(
                (
                    name1,
                    name2,
                    connecting_edges,
                )
            )
        else:
            mod_id = name2.split("_")[1]
            if external_mod_connections.get(mod_id) is None:
                external_mod_connections[mod_id] = []
            external_mod_connections[mod_id].append(
                (
                    name1,
                    connecting_edges,
                )
            )
    return res_res_connections, external_mod_connections


def decompose(mol: rdkit.Chem.rdchem.Mol, cx_smarts_db: dict):
    res_matches = get_res_matches(mol, cx_smarts_db)
    mol = propagate_matches(mol, res_matches)
    G = mol_to_nx(mol)

    native_atom_ids = nx.get_node_attributes(G, "ResID").keys()
    external_modification_atom_ids = G.nodes - native_atom_ids

    G_external_modifications = G.subgraph(external_modification_atom_ids)
    g = (
        G_external_modifications.subgraph(c)
        for c in nx.connected_components(G_external_modifications)
    )
    modification_graphs = list(g)
    return mol, res_matches, modification_graphs


def get_subgraph_tuples(res_matches, modification_graphs, G: nx.classes.graph.Graph):
    subgraph_tuples = []

    for res_id in sorted(res_matches.keys()):
        res_name, res_mol, res_match_atoms = res_matches[res_id]
        res_subgraph = G.subgraph(res_match_atoms)
        res_id_name = "Res_%s_%s" % (res_id, res_name)
        res_tuple = (res_id_name, res_subgraph)
        subgraph_tuples.append(res_tuple)

    n_mod_graphs = len(modification_graphs)
    for i in range(n_mod_graphs):
        tup = ("Mod_%s" % (i + 1), modification_graphs[i])
        subgraph_tuples.append(tup)
    return subgraph_tuples


def get_subgraph_pair_tuples(subgraph_tuples):
    subgraph_pair_tuples = []

    n_subgraphs = len(subgraph_tuples)
    for i in range(n_subgraphs):
        name_i, subgraph_i = subgraph_tuples[i]
        for j in range(i + 1, n_subgraphs):
            name_j, subgraph_j = subgraph_tuples[j]
            subgraph_pair_tuple = (name_i, name_j, subgraph_i, subgraph_j)
            subgraph_pair_tuples.append(subgraph_pair_tuple)

    return subgraph_pair_tuples


def get_connections(subgraph_pair_tuples, G: nx.classes.graph.Graph):
    connections = []

    for res1, res2, subgraph1, subgraph2 in subgraph_pair_tuples:
        edges = get_connecting_edges(G, subgraph1, subgraph2)
        if edges:
            connections.append((res1, res2, edges))

    return connections


def full_decomposition(mol: rdkit.Chem.rdchem.Mol, cx_smarts_db: dict):

    mol, res_matches, modification_graphs = decompose(mol, cx_smarts_db)
    G = mol_to_nx(mol)

    subgraph_tuples = get_subgraph_tuples(res_matches, modification_graphs, G)
    subgraph_pair_tuples = get_subgraph_pair_tuples(subgraph_tuples)
    connections = get_connections(subgraph_pair_tuples, G)
    res_res_connections, external_mod_connections = split_connections_by_type(
        connections
    )
    internal_connections = process_internal_connections(
        res_res_connections, res_matches, G
    )
    external_connections = process_external_connections(
        external_mod_connections, res_matches, modification_graphs, G
    )
    res_return = {res_id: res_matches[res_id][0] for res_id in res_matches}
    return res_return, {
        "internal_modifications": internal_connections,
        "external_modifications": external_connections,
    }


def translate_attachment_points_on_seq(attachment_points_on_seq, offset=0):
    keys = sorted(list(attachment_points_on_seq.keys()))[::-1]
    for attachment_point_id in keys:
        attachment_points_on_seq[
            attachment_point_id + offset
        ] = attachment_points_on_seq.pop(attachment_point_id)
    return attachment_points_on_seq


def translate_mod_smiles(smiles: str, offset=0):
    mol = rdkit.Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            isotope = atom.GetIsotope()
            isotope += offset
            atom.SetIsotope(isotope)
    smiles = rdkit.Chem.MolToSmiles(mol)
    return smiles


def translate_external_modification(mod: dict, offset=0) -> dict:
    mod_smiles = mod.get("mod_smiles")
    mod_smiles = translate_mod_smiles(mod_smiles, offset=offset)

    attachment_points_on_seq = mod.get("attachment_points_on_seq")
    attachment_points_on_seq = translate_attachment_points_on_seq(
        attachment_points_on_seq, offset=offset
    )

    max_attachment_point_id = mod.get("max_attachment_point_id")
    max_attachment_point_id += offset

    mod = {
        "mod_smiles": mod_smiles,
        "max_attachment_point_id": max_attachment_point_id,
        "attachment_points_on_seq": attachment_points_on_seq,
    }
    return mod


def sequence_dict_to_string(sequence_dict: dict) -> str:
    sequence_string = ""
    for res_id in sorted(sequence_dict.keys()):
        symbol = sequence_dict.get(res_id)
        if len(symbol) > 1:
            symbol = "{%s}" % symbol
        sequence_string = sequence_string + symbol
    return sequence_string


def decompose_residues_internal(residues_internal, cx_smarts_db: dict) -> tuple:
    """

    """
    sequence_dict = {}

    internal_modifications = []
    external_modifications = []

    for residue in residues_internal:
        mol = nx_to_mol(residue)
        res_matches, modifications = full_decomposition(mol, cx_smarts_db)
        internal_modifications += modifications["internal_modifications"]

        if modifications.get("external_modifications"):
            for modification in modifications.get("external_modifications"):
                external_modifications.append(modification)

        for res_id in res_matches:
            sequence_dict[int(res_id)] = res_matches[res_id]

    sequence_string = sequence_dict_to_string(sequence_dict)
    return sequence_string, internal_modifications, external_modifications
