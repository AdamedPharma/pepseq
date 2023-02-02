import networkx as nx
import rdkit
from Peptide.utils.chemistry.mol_to_nx_translation import mol_to_nx, nx_to_mol


def match_to_res_id(mol, ResID, cx_smarts_db):
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


def get_res_matches(mol, cx_smarts_db):
    res_matches = {}
    G = mol_to_nx(mol)
    ResIDs_by_atom = nx.get_node_attributes(G, "ResID")
    ResIDs = sorted(list(set(ResIDs_by_atom.values())))

    for ResID in ResIDs:
        (max_aa, max_aa_mol, max_match) = match_to_res_id(mol, ResID, cx_smarts_db)
        res_matches[ResID] = (max_aa, max_aa_mol, max_match)
    return res_matches


def propagate_matches(mol, res_matches):
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


def get_connecting_edges(G, H, I):
    out = [e for e in G.edges if e[0] in H and e[1] in I or e[0] in I and e[1] in H]
    return out


def get_connection_info(G, G_native, G_mod, attachment_point_id=0):
    connecting_edges = get_connecting_edges(G, G_native, G_mod)
    attachment_points = {}

    G_mod_copy = G_mod.copy()

    for connecting_edge in connecting_edges:
        attachment_point_id += 1
        attachment_point_on_seq = (G_native.nodes & set(connecting_edge)).pop()
        attachment_point_on_mod = (G_mod.nodes & set(connecting_edge)).pop()

        res_id, atom_name = G.nodes[attachment_point_on_seq]["ResID"], G.nodes[
            attachment_point_on_seq
        ].get("AtomName")
        attachment_points[attachment_point_id] = res_id, atom_name

        # dodajemy dummyAtom
        G_mod_copy.add_node(
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
        G_mod_copy.add_edge(
            attachment_point_on_mod,
            "%d*" % attachment_point_id,
            **{"bond_type": rdkit.Chem.rdchem.BondType.SINGLE}
        )
        #

    mod_mol = nx_to_mol(G_mod_copy)
    mod_smiles = rdkit.Chem.MolToSmiles(mod_mol)
    return {
        "mod_smiles": mod_smiles,
        "attachment_points_on_seq": attachment_points,
        "max_attachment_point_id": attachment_point_id,
    }


def get_external_modifications_info(G, G_native, modification_graphs):
    external_modifications_info = []
    attachment_point_id = 0
    for G_mod in modification_graphs:
        external_modification_info = get_connection_info(
            G, G_native, G_mod, attachment_point_id=attachment_point_id
        )
        external_modifications_info.append(external_modification_info)
        attachment_point_id = external_modification_info["max_attachment_point_id"]
    return external_modifications_info


def process_residue_candidate_graph(mol, cx_smarts_db):
    """
    this processes case 4
    we need something to process case 3
    """
    res_matches = get_res_matches(mol, cx_smarts_db)
    mol = propagate_matches(mol, res_matches)

    #
    # now find a way to
    #

    G = mol_to_nx(mol)

    native_atom_ids = nx.get_node_attributes(G, "ResID").keys()
    external_modification_atom_ids = G.nodes - native_atom_ids

    if external_modification_atom_ids:
        # should cover case 1 and 2 and 4

        G_native = G.subgraph(native_atom_ids)

        G_copy = G.subgraph(external_modification_atom_ids)
        g = (G_copy.subgraph(c) for c in nx.connected_components(G_copy))
        modification_graphs = list(g)

        external_modifications_info = get_external_modifications_info(
            G, G_native, modification_graphs
        )
    else:
        external_modifications_info = []

    return mol, res_matches, external_modifications_info


def process_res_res_connection(res_atoms_1, res_atoms_2, connecting_edges):
    attachment_point_pairs = []
    for connecting_edge in connecting_edges:
        res_1_attachment_point = (set(connecting_edge) & res_atoms_1).pop()
        res_2_attachment_point = (set(connecting_edge) & res_atoms_2).pop()
        tup = (res_1_attachment_point, res_2_attachment_point)
        attachment_point_pairs.append(tup)
    return attachment_point_pairs


def process_internal_connections(connections, res_matches, G):
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


def process_external_connections(connections, res_matches, modification_graphs, G):
    attachment_point_id = 0
    external_modifications = []

    for mod_id in connections:
        attachment_points_on_seq = {}
        mod_graph = modification_graphs[int(mod_id) - 1].copy()
        mod_atoms = set(mod_graph.nodes)

        for res_name, connecting_edges in connections[mod_id]:
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


def decompose(mol, cx_smarts_db):
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


def get_subgraph_tuples(res_matches, modification_graphs, G):
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


def get_connections(subgraph_pair_tuples, G):
    connections = []

    for res1, res2, subgraph1, subgraph2 in subgraph_pair_tuples:
        edges = get_connecting_edges(G, subgraph1, subgraph2)
        if edges:
            connections.append((res1, res2, edges))

    return connections


def full_decomposition(mol, cx_smarts_db):
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


def translate_mod_smiles(smiles, offset=0):
    mol = rdkit.Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            isotope = atom.GetIsotope()
            isotope += offset
            atom.SetIsotope(isotope)
    smiles = rdkit.Chem.MolToSmiles(mol)
    return smiles


def translate_external_modification(mod, offset=0):
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


def decompose_residues_internal(residues_internal, cx_smarts_db):
    d_seq = {}
    # max_attachment_point_id = 0
    internal_modifications = []
    external_modifications = []

    for residue in residues_internal:
        mol = nx_to_mol(residue)
        res_matches, modifications = full_decomposition(mol, cx_smarts_db)
        internal_modifications += modifications["internal_modifications"]

        if modifications.get("external_modifications"):
            for modification in modifications.get("external_modifications"):
                # modification = translate_external_modification(
                #    modification, offset=max_attachment_point_id
                # )
                external_modifications.append(modification)
            # max_attachment_point_id = modification.get("max_attachment_point_id")

        for res_id in res_matches:
            d_seq[int(res_id)] = res_matches[res_id]

    seq = "".join(d_seq[res_id] for res_id in sorted(d_seq.keys()))
    return seq, internal_modifications, external_modifications
