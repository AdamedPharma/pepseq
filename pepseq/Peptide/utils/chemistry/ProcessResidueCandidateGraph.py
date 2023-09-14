import networkx as nx
import rdkit
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import (mol_to_nx,
                                                                  nx_to_mol)

def get_matches(mol: rdkit.Chem.rdchem.Mol, cx_smarts_db: dict) -> dict:
    """
    Input:
        FragmentMol: rdkit.Chem.rdchem.Mol

        cx_smarts_db: {'C': Cysteine_SMARTS_code, ...}

    Action:
        For each of the amino acid SMARTS codes present in database
        they are matched against FragmentMolecule
        Each matched Amino Acid specie is returned as dictionary
        together with matched substructures

    Output:
        {
            'C': (C_mol: rdkit.Chem.rdchem.Mol, (1,2,3), (4,5,6))
        }
    """
    matches_dict = {}
    for aa in cx_smarts_db:
        aa_smarts = cx_smarts_db[aa]
        aa_mol = rdkit.Chem.MolFromSmarts(aa_smarts)
        matches = mol.GetSubstructMatches(aa_mol)
        if matches:
            matches_dict[aa] = (aa_mol, matches)
    return matches_dict


def get_match_cover(G: nx.classes.graph.Graph, ResID: str):
    """
    match covers only one Residue and it is Residue studied

    Input:

    G - matched subgraph (substructure of molecular graph)
    """

    match_res_ids = set(
                nx.get_node_attributes(G, "ResID").values()
            )
    if match_res_ids != set([ResID]):
        return 0
    else:
        return len(G)


def match_molecular_graph_to_res_id(G: nx.classes.graph.Graph, ResID: str, matches_dict: dict) -> tuple:
    """
    Input:
        FragmentMol: rdkit.Chem.rdchem.Mol

        matches_dict: 
            {
                'C': (C_mol: rdkit.Chem.rdchem.Mol, (1,2,3), (4,5,6))
            }
    ResID:   for example '1' (one of the many ResIDs determined by fitting Peptide Backbone)

    Output:

    (max_aa: amino_acid specie with greatest cover over fragment portion including ResX (X=ResID)
    and only ResX,
    max_aa_mol: molecule substructure covered by amino acid specie: rdkit.Chem.rdchem.Mol
    max_match: e.g. (1,2,3) atom_ids covered by match
    )

    Action:

    for each of the substructure matches grouped by amino acid specie
    we filter the ones covering only one residue and return the match covering
    the biggest posrtion of that residue (i.e. Cysteine is preferred over Alaine)
    if Alanine and Cysteine are matched; Cysteine match will be returned

    """
    max_cover = 0
    max_aa = None
    max_match = None
    max_aa_mol = None
    
    for aa in matches_dict:
        aa_mol, matches = matches_dict[aa]

        for match in matches:
            match_subgraph = G.subgraph(match)
            match_cover = get_match_cover(match_subgraph, ResID)
            if match_cover > max_cover:
                max_cover = match_cover
                max_aa = aa
                max_match = match
                max_aa_mol = aa_mol
    return (max_aa, max_aa_mol, max_match)



def get_res_matches(mol: rdkit.Chem.rdchem.Mol, cx_smarts_db: dict) -> dict:
    """

    We use nx.classes.graph.Graph representation of modified peptide molecules

    Input:

    mol - Fragment Molecule: rdkit.Chem.rdchem.Mol
    cx_smarts_db - dictionary containing SMARTS codes for each of the Amino Acid
    species


    Output:

    {
        'ResID1':  (max_aa, max_aa_mol, max_match),
        'ResID2':  (max_aa, max_aa_mol, max_match),
        ...
        }
    
    Action:

    for each ResidueID found in fragment we find the best match by
    Amino Acid specie

    """

    matches_dict = get_matches(mol, cx_smarts_db)

    G = mol_to_nx(mol)
    res_matches = {}

    fragment_ResIDs = sorted(list(set(nx.get_node_attributes(G, "ResID").values())))


    for ResID in fragment_ResIDs:
        (max_aa, max_aa_mol, max_match) = match_molecular_graph_to_res_id(G, ResID, matches_dict)
        atom_names_dict = {}
        for i in range(len(max_match)):
            atom_match = max_aa_mol.GetAtomWithIdx(i)
            atom_id = max_match[i]
            if "AtomName" in atom_match.GetPropNames():
                atom_names_dict[atom_id] = atom_match.GetProp("AtomName")

        res_matches[ResID] = (max_aa, atom_names_dict, max_match)
    return res_matches


def propagate_matches_on_molecular_graph(G: nx.classes.graph.Graph, res_matches: dict) -> nx.classes.graph.Graph:
    """

    We use nx.classes.graph.Graph representation of modified peptide molecules

    Input:

    G - Fragment Molecular Graph: nx.classes.graph.Graph
    res_matches - dictionary containing for each sequence amino acid residue the AminoAcid specie name(symbol),
    atom_names for fitted molecule substructure; atom_ids matched by residue

    Output:

    G - Fragment Molecular Graph: nx.classes.graph.Graph with atom nodes labeled with ResID(s); ResName(s)
    and AtomNames (where present)

    """

    for ResID in res_matches:
        ResName, atom_names_dict, match = res_matches[ResID]
        nx.set_node_attributes(G, {atom_id: ResID for atom_id in match}, "ResID")
        nx.set_node_attributes(G, {atom_id: ResName for atom_id in match}, "ResName")
        nx.set_node_attributes(G, atom_names_dict, "AtomName")
    return G


def get_connecting(edges: list, nodes1: list, nodes2: list):
    return [
        edge
        for edge in edges
        if (edge[0] in nodes1 and edge[1] in nodes2)
        or (edge[0] in nodes2 and edge[1] in nodes1)
    ]


def get_connecting_edges(union_graph: nx.classes.graph.Graph, H_graph: nx.classes.graph.Graph,
                        I_graph: nx.classes.graph.Graph) -> list:
    edges = list(union_graph.edges)
    nodes1 = list(H_graph.nodes)
    nodes2 = list(I_graph.nodes)
    return get_connecting(edges, nodes1, nodes2)


def get_atom_pairs(atoms_1: set, atoms_2: set, edges: list) -> list:
    atom_pairs = []
    for edge in edges:
        atom1 = (set(edge) & atoms_1).pop()
        atom2 = (set(edge) & atoms_2).pop()
        tup = (atom1, atom2)
        atom_pairs.append(tup)
    return atom_pairs


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

        attachment_point_pairs = get_atom_pairs(
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
    """
    """
    res_ids = []

    length = len(connection)

    for i in range(length):
        res_name = connection[i][0]
        res_id = int(res_name.split("_")[1])
        res_ids.append((res_id, i))

    sorted_res_ids = sorted(res_ids)
    return [connection[i] for (res_id, i) in sorted_res_ids]





def add_connection_point_to_molecular_graph(G: nx.classes.graph.Graph, point_id: int, atom_id: int) -> nx.classes.graph.Graph:
    G.add_node(
        "%d*" % point_id,
        **{
            "atomic_num": 0,
            "formal_charge": 0,
            "chiral_tag": rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
            "hybridization": rdkit.Chem.rdchem.HybridizationType.SP3,
            "num_explicit_hs": 0,
            "is_aromatic": False,
            "isotope": point_id,
        }
    )
    G.add_edge(
        atom_id,
        "%d*" % point_id,
        bond_type=rdkit.Chem.rdchem.BondType.SINGLE,
    )

    return G


def get_residue_id_and_atoms(res_matches, res_name):
    list_res = res_name.split("_")

    res_id = list_res[1]

    res_atoms = set(res_matches[res_id][2])

    return res_id, res_atoms


def process_external_modification(G_mod: nx.classes.graph.Graph, mod_bonds, res_matches, atom_names_dict, point_id):
    """
    For each modification (like staple) we need unambigious info to 

    Input:
        G_mod - molecular Graph representing external_modification (like staple)
        bonds_to_residue - 

    """
    points_on_seq = {}
    mod_atoms = set(G_mod.nodes)

    bonds_by_residue = sorted_connection(mod_bonds)

    for res_name, bonds in bonds_by_residue:

        res_id, res_atoms = get_residue_id_and_atoms(res_matches, res_name)

        atom_pairs = get_atom_pairs(res_atoms, mod_atoms, bonds)

        for res_atom, mod_atom in atom_pairs:
            point_id += 1
            AtomName = atom_names_dict.get(res_atom)
            points_on_seq[point_id] = {
                "attachment_point_id": point_id,
                "ResID": res_id,
                "AtomName": AtomName,
                "ResidueName": "",
            }
            G_mod = add_connection_point_to_molecular_graph(G_mod,
                            point_id, mod_atom)

    mod_mol = nx_to_mol(G_mod)
    mod_smiles = rdkit.Chem.MolToSmiles(mod_mol)

    mod_json = {
        "smiles": mod_smiles,
        "max_attachment_point_id": point_id,
        "attachment_points_on_sequence": points_on_seq,
    }

    return mod_json, point_id


def process_external_connections(modifications, res_matches, modification_graphs, G: nx.classes.graph.Graph):
    attachment_point_id = 0
    external_modifications = []

    atom_names_dict = nx.get_node_attributes(G, 'AtomName')

    for mod_id in modifications:
        mod_graph = modification_graphs[int(mod_id) - 1].copy()
        mod_bonds = modifications.get(mod_id)

        external_modification, attachment_point_id = process_external_modification(
            mod_graph, mod_bonds, res_matches, atom_names_dict, attachment_point_id)
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


def get_modification_graphs_from_fragment(G: nx.classes.graph.Graph) -> list:
    """

    Input:

    G - molecular graph labeled with Sequence Amino Acid Residue ID and Residue Name

    Output:

    modification_graphs - Molecular Graphs for external modifications

    Action:

    we get atom_ids labeled with ResID (parts of amino acid sequence)

    we select nonlabeled atom_ids (they are external modification(s) e.g. staple)

    There can be more than one external modifications (they are not connected
      with each other)


    """
    native_atom_ids = nx.get_node_attributes(G, "ResID").keys()
    external_modification_atom_ids = G.nodes - native_atom_ids
    G_external_modifications = G.subgraph(external_modification_atom_ids)
    g = (
        G_external_modifications.subgraph(c)
        for c in nx.connected_components(G_external_modifications)
    )
    modification_graphs = list(g)
    return modification_graphs


def decompose(mol: rdkit.Chem.rdchem.Mol, cx_smarts_db: dict) -> tuple:
    """
    We use rdkit.Chem.rdchem.Mol representation of fragment

    """

    res_matches = get_res_matches(mol, cx_smarts_db)
    G = mol_to_nx(mol)
    
    G = propagate_matches_on_molecular_graph(G, res_matches)
    modification_graphs = get_modification_graphs_from_fragment(G)
    mol = nx_to_mol(G)
    return mol, res_matches, modification_graphs

def get_internal_connections_subgraph_tuples():
    """
    WIP
    """
    return


def get_subgraph_tuples(res_matches, modification_graphs_nodes, G: nx.classes.graph.Graph):
    subgraph_tuples = []

    for res_id in sorted(res_matches.keys()):
        res_name, res_atoms_dict, res_match_atoms = res_matches[res_id]
        res_subgraph = G.subgraph(res_match_atoms)
        res_id_name = "Res_%s_%s" % (res_id, res_name)
        res_tuple = (res_id_name, res_subgraph)
        subgraph_tuples.append(res_tuple)

    n_mod_graphs = len(modification_graphs_nodes)

    for i in range(n_mod_graphs):
        tup = ("Mod_%s" % (i + 1), modification_graphs_nodes[i])# list(modification_graphs[i].nodes) )
        subgraph_tuples.append(tup)
    return subgraph_tuples


def get_subgraph_pair_tuples(subgraph_tuples):
    subgraph_pair_tuples = []

    n_subgraphs = len(subgraph_tuples)
    for i in range(n_subgraphs):
        name_i, subgraph_i_nodes = subgraph_tuples[i]
        for j in range(i + 1, n_subgraphs):
            name_j, subgraph_j_nodes = subgraph_tuples[j]
            subgraph_pair_tuple = (name_i, name_j, list(subgraph_i_nodes), subgraph_j_nodes)
            subgraph_pair_tuples.append(subgraph_pair_tuple)

    return subgraph_pair_tuples


def get_connections(subgraph_pair_tuples, G_edges: list):
    connections = []

    for res1, res2, nodes1, nodes2 in subgraph_pair_tuples:

        edges = get_connecting(G_edges, nodes1, nodes2)
        if edges:
            connections.append((res1, res2, edges))

    return connections


def full_decomposition(mol:rdkit.Chem.rdchem.Mol, cx_smarts_db: dict):
    """
    We use rdkit.Chem.rdchem.Mol representation of residue candidate

    Input:

    mol - rdkit.Chem.rdchem.Mol

    Output:

    res_return - {
        'ResID1': 'ResName',
        ...
    }

    internal connections:

    external connections:


    Action:

    
    """

    mol, res_matches, modification_graphs = decompose(mol, cx_smarts_db)
    G = mol_to_nx(mol)

    modification_graphs_nodes = [list(graph.nodes) for graph in modification_graphs]

    subgraph_tuples = get_subgraph_tuples(res_matches, modification_graphs_nodes, G)

    subgraph_pair_tuples = get_subgraph_pair_tuples(subgraph_tuples)
    connections = get_connections(subgraph_pair_tuples, list(G.edges))


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


def decompose_residues_internal(residues_internal: list, cx_smarts_db: dict) -> tuple:
    """


    """
    sequence_dict = {}

    internal_modifications = []
    external_modifications = []

    for residue_candidate in residues_internal:
        mol = nx_to_mol(residue_candidate)
        res_matches, modifications = full_decomposition(mol, cx_smarts_db)
        internal_modifications += modifications["internal_modifications"]

        if modifications.get("external_modifications"):
            for modification in modifications.get("external_modifications"):
                external_modifications.append(modification)

        for res_id in res_matches:
            sequence_dict[int(res_id)] = res_matches[res_id]

    sequence_string = sequence_dict_to_string(sequence_dict)
    return sequence_string, internal_modifications, external_modifications
