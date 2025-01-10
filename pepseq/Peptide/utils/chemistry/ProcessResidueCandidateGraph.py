from tqdm import tqdm

import networkx as nx
import rdkit

from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import mol_to_nx, nx_to_mol
from pepseq.augmenting_db_json import get_Nter_versions_cxsmarts_db


def get_match(mol=None, aa_smarts=None, useChirality=True):
    """
    For each of the amino acid SMARTS codes present in database
    they are matched against FragmentMolecule
    Each matched Amino Acid specie is returned as dictionary
    together with matched substructures

    :param mol: peptide molecule in rdkit.Chem.rdchem.Mol format
    :type mol: rdkit.Chem.rdchem.Mol

    :param aa_smarts: SMARTS code for amino acid
    :type aa_smarts: str

    :param useChirality: whether to account for stereoisometry/ chirality
    :type useChirality: bool

    :return: tuple containing amino acid molecule and matched substructures
    :rtype: tuple
    """
    aa_mol = rdkit.Chem.MolFromSmarts(aa_smarts)
    matches = mol.GetSubstructMatches(aa_mol, useChirality=useChirality)
    return (aa_mol, matches)


def get_match_with_external_modifications():
    """
    Rationale:
        Simplest case:
        1. One Amino Acid
        One Residue
        For each fragment we need to decompose it into residue(s) and external modifications
        That is necessary because we want to limit number of external modification per Residue
        (or fragment) into one. That is why we need two subgraphs
        one matched by residue .e.g. methylserine (or X)
        and another subgraph (the rest)

    """
    return


def get_matches(
    mol: rdkit.Chem.rdchem.Mol, cx_smarts_db: dict, useChirality=True
) -> dict:
    """
    For each of the amino acid SMARTS codes present in database
    they are matched against FragmentMolecule
    Each matched Amino Acid specie is returned as dictionary
    together with matched substructures

    :param mol: peptide molecule
    :type mol: rdkit.Chem.rdchem.Mol

    :param cx_smarts_db: dictionary containing SMARTS codes for each of the Amino Acid species
    :type cx_smarts_db: dict

    :param useChirality: whether to account for stereoisometry/ chirality
    :type useChirality:

    :return dictionary showing for each amino acid specie the matched substructures
    :rtype dict
    """
    matches_dict = {}
    for aa in cx_smarts_db:
        aa_smarts = cx_smarts_db[aa]
        aa_mol, matches = get_match(
            mol=mol, aa_smarts=aa_smarts, useChirality=useChirality
        )

        if matches:
            matches_dict[aa] = (aa_mol, matches)
    return matches_dict


def get_match_cover(G: nx.classes.graph.Graph, ResID: str):
    """
    match covers only one Residue and it is Residue studied

    :param G: graph substructure of molecular graph
    :type G: nx.classes.graph.Graph

    :param ResID: residue ID
    :type ResID: str

    :return: number of nodes in the subgraph matching the Residue
    :rtype: int
    """

    match_res_ids = set(nx.get_node_attributes(G, "ResID").values())
    if match_res_ids != set([ResID]):
        return 0
    else:
        return len(G)


def match_molecular_graph_to_res_id(
    G: nx.classes.graph.Graph, ResID: str, matches_dict: dict
) -> tuple:
    """
    We match molecular graph to ResidueID
    for each of the substructure matches grouped by amino acid specie
    we filter the ones covering only one residue and return the match covering
    the biggest posrtion of that residue (i.e. Cysteine is preferred over Alanine)
    if Alanine and Cysteine are matched; Cysteine match will be returned.

    :param G: molecular graph
    :type G: nx.classes.graph.Graph

    :param: ResID: residue ID
    :type ResID: str

    :param matches_dict: dictionary containing for each amino acid specie the matched substructures
    :type matches_dict: dict

    :return: tuple
    :rtype: tuple
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


def get_mod_graphs(G, native_atom_ids, ResID):
    """
    We get molecular graphs for external modifications

    :param G: molecular graph
    :type G: nx.classes.graph.Graph

    :param native_atom_ids: set of atom ids for native atoms of the peptide
    :type native_atom_ids: set

    :param ResID: amino acid residue ID
    :type ResID: str

    :return: list of molecular graphs for external modifications
    :rtype: list
    """
    nx.set_node_attributes(G, {atom_id: ResID for atom_id in native_atom_ids}, "ResID")
    native_atom_ids = nx.get_node_attributes(G, "ResID").keys()
    external_modification_atom_ids = G.nodes - native_atom_ids
    G_external_modifications = G.subgraph(external_modification_atom_ids)
    g = (
        G_external_modifications.subgraph(c)
        for c in nx.connected_components(G_external_modifications)
    )
    modification_graphs = list(g)
    return modification_graphs


def get_n_subst_dict(G, matches_dict) -> dict:
    """
    For each of the amino acid specie matched to the fragment
    we get the number of possible substitutions

    :param G: molecular graph for fragment
    :type G: nx.classes.graph.Graph

    :param matches_dict: dictionary containing for each amino acid specie the matched substructures
    :type matches_dict: dict

    :return: dictionary containing for each amino acid specie the number of possible substitutions
    :rtype: dict
    """
    n_subst = {}
    ResID = "ResIDVal"

    for aa_name in tqdm(matches_dict):
        aa_mol, matches = matches_dict.get(aa_name)
        native_atom_ids = set([])
        for match in matches:
            native_atom_ids |= set(match)
        G_copy = G.copy()
        mod_graphs = get_mod_graphs(G_copy, native_atom_ids, ResID)
        n_subst[aa_name] = len(mod_graphs)
    return n_subst


def filter_n_subst(G, matches_dict, n_subst_limit) -> dict:
    """
    For each of the amino acid specie matched to the fragment
    we filter the ones with number of substitutions less or equal to n_subst_limit

    :param G: molecular graph for fragment nx.classes.graph.Graph
    :type G: nx.classes.graph.Graph
    :param matches_dict: dictionary containing for each amino acid specie the matched substructures dict
    :type matches_dict: dict

    :param n_subst_limit: limit on number of substitutions int
    :type n_subst_limit: int

    :return: dictionary containing for each amino acid specie the matched
     substructures with number of substitutions less or equal to n_subst_limit
    :rtype: dict
    """
    n_subst_dict = get_n_subst_dict(G, matches_dict)
    aa_names_n_subst_elt = [
        aa_name
        for aa_name in n_subst_dict
        if n_subst_dict.get(aa_name) <= n_subst_limit
    ]
    return {k: matches_dict.get(k) for k in matches_dict if k in aa_names_n_subst_elt}


def get_res_matches(
    mol: rdkit.Chem.rdchem.Mol,
    cx_smarts_db: dict,
    useChirality=True,
    n_subst_limit=None,
) -> dict:
    """

    We use nx.classes.graph.Graph representation of modified peptide molecules

    :param mol Fragment Molecule: rdkit.Chem.rdchem.Mol
    :type mol: rdkit.Chem.rdchem.Mol

    :param cx_smarts_db: dictionary containing SMARTS codes for each of the Amino Acid species: dict
    :type cx_smarts_db: dict

    :param useChirality: whether to account for stereoisometry/ chirality bool
    :type useChirality: bool

    :param n_subst_limit: limit on number of substitutions int
    :type n_subst_limit: int

    :return: dictionary containing for each sequence amino acid residue the AminoAcid
     specie name(symbol),  atom_names for fitted molecule substructure;
     atom_ids matched by residue dict
    :rtype: dict

    Output:

    {
        'ResID1':  (max_aa, atom_names_dict, max_match),
        'ResID2':  (max_aa, atom_names_dict, max_match),
        ...
        }

    Action:

    for each ResidueID found in fragment we find the best match by
    Amino Acid specie

    """
    G = mol_to_nx(mol)
    res_ids = nx.get_node_attributes(G, "ResID")

    cx_smarts_db_copy = cx_smarts_db.copy()
    nter_versions = get_Nter_versions_cxsmarts_db(cx_smarts_db_copy)
    cx_smarts_db_copy.update(nter_versions)

    if "1" in res_ids.values():
        matches_dict = get_matches(mol, cx_smarts_db_copy, useChirality=useChirality)
    else:
        matches_dict = get_matches(mol, cx_smarts_db, useChirality=useChirality)

    if n_subst_limit is not None:
        fltrd = filter_n_subst(G, matches_dict, n_subst_limit)
        if fltrd:
            matches_dict = fltrd

    res_matches = {}

    fragment_ResIDs = sorted(list(set(nx.get_node_attributes(G, "ResID").values())))

    for ResID in fragment_ResIDs:
        (max_aa, max_aa_mol, max_match) = match_molecular_graph_to_res_id(
            G, ResID, matches_dict
        )
        atom_names_dict = {}
        for i in range(len(max_match)):
            atom_match = max_aa_mol.GetAtomWithIdx(i)
            atom_id = max_match[i]
            if "AtomName" in atom_match.GetPropNames():
                atom_names_dict[atom_id] = atom_match.GetProp("AtomName")

        res_matches[ResID] = (max_aa, atom_names_dict, max_match)
    return res_matches


def propagate_matches_on_molecular_graph(
    G: nx.classes.graph.Graph, res_matches: dict
) -> nx.classes.graph.Graph:
    """

    We use nx.classes.graph.Graph representation of modified peptide molecules

    :param G - Fragment Molecular Graph: nx.classes.graph.Graph
    :type G - nx.classes.graph.Graph

    :param res_matches - dictionary containing for each sequence amino acid residue the AminoAcid specie name(symbol),
    atom_names for fitted molecule substructure; atom_ids matched by residue dict
    :type res_matches - dict

    :return G - Fragment Molecular Graph: nx.classes.graph.Graph with atom nodes labeled with ResID(s); ResName(s)
    atom_names for fitted molecule substructure; atom_ids matched by residue
    :rtype G - nx.classes.graph.Graph
    """

    for ResID in res_matches:
        ResName, atom_names_dict, match = res_matches[ResID]
        nx.set_node_attributes(G, {atom_id: ResID for atom_id in match}, "ResID")
        nx.set_node_attributes(G, {atom_id: ResName for atom_id in match}, "ResName")
        nx.set_node_attributes(G, atom_names_dict, "AtomName")
    return G


def get_connecting(edges: list, nodes1: list, nodes2: list):
    """
    Get connecting edges between two sets of nodes

    :param edges: list
    :type list

    :param nodes1: list
    :type list

    :param nodes2: list
    :type list

    :return: list of connecting edges
    :rtype: list
    """
    return [
        edge
        for edge in edges
        if (edge[0] in nodes1 and edge[1] in nodes2)
        or (edge[0] in nodes2 and edge[1] in nodes1)
    ]


def get_connecting_edges(
    union_graph: nx.classes.graph.Graph,
    H_graph: nx.classes.graph.Graph,
    I_graph: nx.classes.graph.Graph,
) -> list:
    """
    Get connecting edges between two subgraphs. Graph is representing molecule
    (peptide).

    :param union_graph: whole graph representing peptide nx.classes.graph.Graph
    :type nx.classes.graph.Graph

    :param H_graph: one subgraph representing residue nx.classes.graph.Graph
    :type nx.classes.graph.Graph

    :param I_graph: one subgraph representing rest of peptide nx.classes.graph.Graph
    :type nx.classes.graph.Graph

    :return: list of connecting edges
    :rtype: list
    """
    edges = list(union_graph.edges)
    nodes1 = list(H_graph.nodes)
    nodes2 = list(I_graph.nodes)
    return get_connecting(edges, nodes1, nodes2)


def get_atom_pairs(atoms_1: set, atoms_2: set, edges: list) -> list:
    """
    Get atom pairs from edges

    :param atoms_1: set of atoms from another molecule
    :type atoms1: set

    :param atoms_2: set of atoms from another molecule
    :type atoms2: set

    :param edges: bonds
    :type edges: list

    :return: list of atom pairs
    :rtype: list
    """
    atom_pairs = []
    for edge in edges:
        atom1 = (set(edge) & atoms_1).pop()
        atom2 = (set(edge) & atoms_2).pop()
        tup = (atom1, atom2)
        atom_pairs.append(tup)
    return atom_pairs


def process_internal_connections(
    connections, res_matches, G: nx.classes.graph.Graph
) -> list:
    """
    For each pair of Residue Subgraphs
    We find whether in Fragment MolecularGraph edges there is a connecting
    edge between Residue1 and Residue2

    :param connections: list of connections between Residue Subgraphs
    :type connections: list

    :param res_matches: dictionary containing for each sequence amino acid residue the AminoAcid specie name(symbol),
    atom_names for fitted molecule substructure; atom_ids matched by residue dict
    :type res_matches: dict

    :param G: molecular graph representing fragment
    :type G: nx.classes.graph.Graph
    """
    bond_jsons = []
    bond_id = 0

    for name1, name2, connecting_edges in connections:
        res_ids = [name.split("_")[:2][1] for name in (name1, name2)]

        res_atoms_1, res_atoms_2 = [set(res_matches[res_id][2]) for res_id in res_ids]

        bonds = get_atom_pairs(res_atoms_1, res_atoms_2, connecting_edges)

        for bond in bonds:
            bond_id += 1

            atom_names = [G.nodes[atom].get("AtomName") for atom in bond]

            bond_json = {
                bond_id: [
                    {
                        "ResID": res_ids[ind],
                        "AtomName": atom_names[ind],
                        "ResidueName": "",
                    }
                    for ind in range(2)
                ]
            }

            bond_jsons.append(bond_json)
    return bond_jsons


def sorted_connection(connection):
    """
    Sort connection atoms by ResID

    :param connection: list of atoms (pair of atoms) in connection (bond)
    :type connection: list

    :return: sorted atoms
    :rtype: list
    """
    res_ids = []

    length = len(connection)

    for i in range(length):
        res_name = connection[i][0]
        res_id = int(res_name.split("_")[1])
        res_ids.append((res_id, i))

    sorted_res_ids = sorted(res_ids)
    return [connection[i] for (res_id, i) in sorted_res_ids]


def add_connection_point_to_molecular_graph(
    G: nx.classes.graph.Graph, point_id: int, atom_id: int
) -> nx.classes.graph.Graph:
    """
    Adds a connection point to a molecular graph.
    This function adds a new node representing a connection point to the given molecular graph `G`.
    The new node is identified by `point_id` and is connected to an existing atom identified by `atom_id`.

    :param G : The molecular graph to which the connection point will be added.
    :type G : nx.classes.graph.Graph

    :param point_id : The identifier for the new connection point.
    :type point_id : int

    :param atom_id : The identifier of the existing atom to which the new connection point will be connected.
    :type atom_id : int

    :return : The updated molecular graph with the new connection point added.
    :rtype : nx.classes.graph.Graph
    """
    node_name = "*:%d" % point_id
    G.add_node(
        node_name,
        **{
            "atomic_num": 0,
            "formal_charge": 0,
            "chiral_tag": rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED,
            "hybridization": rdkit.Chem.rdchem.HybridizationType.SP3,
            "num_explicit_hs": 0,
            "is_aromatic": False,
            "isotope": 0,
            "dummyLabel": "*",
            "molAtomMapNumber": str(point_id),
        }
    )

    G.add_edge(
        atom_id,
        node_name,
        bond_type=rdkit.Chem.rdchem.BondType.SINGLE,
    )

    return G


def get_residue_id_and_atoms(res_matches, res_name):
    """
    Extracts the residue ID and atoms from the given residue matches and residue name.


    :param res_matches: A dictionary where keys are residue IDs and values are lists or tuples
    :type res_matches: dict

    :param res_name: A string representing the residue name, expected to be in the format "prefix_residueID"
    :type res_name: str

    :return: A tuple containing the residue ID and a set of atoms.
    :rtype: tuple
    """

    list_res = res_name.split("_")

    res_id = list_res[1]

    res_atoms = set(res_matches[res_id][2])

    return res_id, res_atoms


def process_external_modification(
    G_mod: nx.classes.graph.Graph, mod_bonds, res_matches, atom_names_dict, point_id
):
    """
    Processes an external modification (such as a staple) on a peptide sequence.
    This function takes a molecular graph representing an external modification and
     integrates it with a peptide sequence by identifying and adding attachment points.

    :param G_mod: Molecular graph representing the external modification.
    :type G_mod: nx.classes.graph.Graph

    :param mod_bonds: List of bonds to the residue.
    :type mod_bonds: list

    :param res_matches: Dictionary of residue matches.
    :type res_matches: dict

    :param atom_names_dict: Dictionary mapping atom names.
    :type atom_names_dict: dict

    :param point_id: Starting point ID for attachment points.
    :type point_id: int

    :return: A tuple containing:
            - mod_json (dict): A dictionary with the following keys:
                - "smiles" (str): SMILES representation of the modified molecule.
                - "max_attachment_point_id" (int): The maximum attachment point ID used.
                - "attachment_points_on_sequence" (dict): A dictionary of attachment points on the sequence.
            - point_id (int): The updated point ID after processing the modification.

    :rtype: tuple
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
            G_mod = add_connection_point_to_molecular_graph(G_mod, point_id, mod_atom)

    mod_mol = nx_to_mol(G_mod)
    mod_smiles = rdkit.Chem.MolToSmiles(mod_mol)

    mod_json = {
        "smiles": mod_smiles,
        "max_attachment_point_id": point_id,
        "attachment_points_on_sequence": points_on_seq,
    }

    return mod_json, point_id


def process_external_connections(
    modifications,
    res_matches,
    modification_graphs,
    G: nx.classes.graph.Graph,
):
    """
    Processes external connections for a given set of modifications and a residue match graph.

    :param modifications: A dictionary where keys are modification IDs and values are lists of bonds.
    :type modifications: dict

    :param res_matches: Dictionary of residue matches.
    :type res_matches: dict

    :param modification_graphs: List of NetworkX graphs representing modification structures.
    :type modification_graphs: list

    :param G: NetworkX graph representing the residue candidate graph.
    :type G: nx.classes.graph.Graph

    :return: A list of external modifications processed from the input data.
    :rtype: list
    """
    attachment_point_id = 0
    external_modifications = []

    atom_names_dict = nx.get_node_attributes(G, "AtomName")

    for mod_id in modifications:
        mod_graph = modification_graphs[int(mod_id) - 1].copy()
        mod_bonds = modifications.get(mod_id)

        external_modification, attachment_point_id = process_external_modification(
            mod_graph,
            mod_bonds,
            res_matches,
            atom_names_dict,
            attachment_point_id,
        )
        external_modifications.append(external_modification)
    return external_modifications


def split_connections_by_type(connections):
    """
    Splits a list of connections into internal bonds and external bonds.


    :param connections: A list of tuples where each tuple contains two names and a bond.
            two names (str) and a bond (any type). The names are expected to be in the
            format "Type_ID", where "Type" indicates the type of the connection (e.g., "Res")
            and "ID" is an identifier.

    :type connections: list

    :return: A tuple containing:
            - internal_bonds (list): A list of tuples representing internal bonds.
            - external_bonds_dict (dict): A dictionary of external bonds.
    :rtype: tuple
    """

    internal_bonds = []
    external_bonds_dict = {}

    for name1, name2, bonds in connections:
        types = [name.split("_")[0] for name in (name1, name2)]
        if set(types) == set(
            [
                "Res",
            ]
        ):
            internal_bonds.append(
                (
                    name1,
                    name2,
                    bonds,
                )
            )
        else:
            mod_id = name2.split("_")[1]
            if external_bonds_dict.get(mod_id) is None:
                external_bonds_dict[mod_id] = []
            external_bonds_dict[mod_id].append(
                (
                    name1,
                    bonds,
                )
            )
    return internal_bonds, external_bonds_dict


def get_modification_graphs_from_fragment(G: nx.classes.graph.Graph) -> list:
    """
    Extracts molecular graphs for external modifications from a given molecular graph.
    Args:
        G (nx.classes.graph.Graph): Molecular graph labeled with Sequence Amino Acid Residue ID and Residue Name.
    Returns:
        list: A list of molecular graphs representing external modifications.
        - Retrieves atom IDs labeled with ResID (parts of the amino acid sequence).
        - Selects non-labeled atom IDs (representing external modifications, e.g., staples).
        - Identifies and extracts subgraphs for each external modification (which are not connected to each other).

    :param G - molecular graph labeled with Sequence Amino Acid Residue ID and Residue Name
    :type G - nx.classes.graph.Graph

    :return modification_graphs - Molecular Graphs for external modifications
    :rtype modification_graphs - list

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


def decompose(
    mol: rdkit.Chem.rdchem.Mol, cx_smarts_db: dict, n_subst_limit=None
) -> tuple:
    """
    Decompose a molecule into its constituent parts based on a given SMARTS database.

    :param mol: The molecule to be decomposed.
    :type mol: rdkit.Chem.rdchem.Mol

    :param cx_smarts_db: A dictionary containing SMARTS patterns for decomposition.
    :type cx_smarts_db: dict

    :param n_subst_limit: The substitution limit for decomposition. Defaults to None.
    :type n_subst_limit: int

    :return: A tuple containing the molecular graph (networkx.Graph), residue matches (list),
     and modification graphs (list).
    :rtype: tuple
    """

    res_matches = get_res_matches(mol, cx_smarts_db, n_subst_limit=n_subst_limit)
    G = mol_to_nx(mol)

    G = propagate_matches_on_molecular_graph(G, res_matches)
    modification_graphs = get_modification_graphs_from_fragment(G)
    return G, res_matches, modification_graphs


def get_internal_connections_subgraph_tuples(
    G: nx.classes.graph.Graph, res_matches: dict
) -> list:
    """
    Get internal connections between Residue Subgraphs.
    This function takes a graph and a dictionary of residue matches, and returns a list of tuples.
    Each tuple contains a string identifier for the residue and a subgraph of the original graph
    corresponding to the atoms of that residue.

    :param G: The graph representing the molecular structure.
    :type G: nx.classes.graph.Graph

    :param res_matches
        A dictionary where keys are residue IDs and values are tuples containing residue name,
        a dictionary of residue atoms, and a list of residue atoms.
        A list of tuples, where each tuple contains a string identifier for the residue and a subgraph
        of the original graph corresponding to the atoms of that residue.
    :type res_matches: dict

    :return: A list of tuples containing the residue identifier and the corresponding subgraph.
    :rtype: list
    """
    subgraph_tuples = []

    for res_id in sorted(res_matches.keys()):
        res_name, res_atoms_dict, res_atoms = res_matches[res_id]
        res_tuple = ("Res_%s_%s" % (res_id, res_name), G.subgraph(res_atoms))
        subgraph_tuples.append(res_tuple)

    return subgraph_tuples


def get_subgraph_tuples(
    res_matches: dict, modification_graphs_nodes, G: nx.classes.graph.Graph
):
    """
    Get subgraphs by combining internal and external subgraph tuples.
    :param res_matches: A dictionary containing residue matches.
    :type res_matches: dict

    :param modification_graphs_nodes: A list of nodes in the modification graphs.
    :type modification_graphs_nodes: list

    :param G: A NetworkX graph object representing the main graph.
    :type G: nx.classes.graph.Graph

    :return: A list of tuples representing internal and external subgraphs.
    :rtype: list
    """

    internal_subgraph_tuples = get_internal_connections_subgraph_tuples(G, res_matches)

    external_subgraph_tuples = [
        ("Mod_%s" % (i + 1), modification_graphs_nodes[i])
        for i in range(len(modification_graphs_nodes))
    ]
    return internal_subgraph_tuples + external_subgraph_tuples


def get_connections(G_edges: list, subgraph_tuples: list):
    """
    Get connections between Residue Subgraphs.
    For each pair of Residue Subgraphs, this function finds whether there is a
     connecting edge between Residue1 and Residue2 in the Fragment MolecularGraph edges.

    :param G_edges: List of edges in the Fragment MolecularGraph.
    :param subgraph_tuples: List of tuples, where each tuple contains the name
                            and nodes of a Residue Subgraph.
    :return: List of tuples, where each tuple contains the names of two Residue
             Subgraphs and the edges connecting them.
    :rtype: list
    """

    bonds = []

    n_subgraphs = len(subgraph_tuples)
    for i in range(n_subgraphs):
        name_i, nodes_i = subgraph_tuples[i]
        for j in range(i + 1, n_subgraphs):
            name_j, nodes_j = subgraph_tuples[j]
            edges = get_connecting(G_edges, nodes_i, nodes_j)
            if edges:
                bonds.append((name_i, name_j, edges))

    return bonds


def full_decomposition(
    mol: rdkit.Chem.rdchem.Mol, cx_smarts_db: dict, n_subst_limit=None
):
    """
    Fully decomposes a molecule into its constituent parts and processes internal and external connections.
    This function uses the RDKit library to represent the molecule and decomposes
     it based on the provided SMARTS database.
    It then processes the internal and external connections between the decomposed parts.

    :param mol: The molecule to be decomposed.
    :type mol: rdkit.Chem.rdchem.Mol

    :param cx_smarts_db: A dictionary containing SMARTS patterns for decomposition.
    :type cx_smarts_db: dict

    :param n_subst_limit: (int, optional) A limit on the number of substitutions. Defaults to None.
    :type n_subst_limit: int

    :return: A tuple containing the residue names and the processed connections.
        - res_names (dict): A dictionary mapping residue IDs to residue names.
        - dict: A dictionary with keys "internal_modifications" and
         "external_modifications" containing the processed connections.

    :rtype: tuple
    """

    G, res_matches, modification_graphs = decompose(
        mol, cx_smarts_db, n_subst_limit=n_subst_limit
    )

    modification_graphs_nodes = [list(graph.nodes) for graph in modification_graphs]

    subgraph_tuples = get_subgraph_tuples(res_matches, modification_graphs_nodes, G)
    bonds = get_connections(list(G.edges), subgraph_tuples)

    res_res_connections, external_mod_connections = split_connections_by_type(bonds)
    internal_connections = process_internal_connections(
        res_res_connections, res_matches, G
    )
    external_connections = process_external_connections(
        external_mod_connections, res_matches, modification_graphs, G
    )
    res_names = {res_id: res_matches[res_id][0] for res_id in res_matches}
    return res_names, {
        "internal_modifications": internal_connections,
        "external_modifications": external_connections,
    }


def translate_attachment_points_on_seq(attachment_points_on_seq, offset=0):
    """
    Translates the keys of the attachment points dictionary by a given offset.
    This function takes a dictionary of attachment points on a sequence and
     shifts the keys by a specified offset. The keys are processed in reverse
     order to avoid key collisions during the translation.

    :param attachment_points_on_seq: A dictionary where keys are attachment point
                                     point identifiers and values are the
                                     corresponding attachment points.

    :type attachment_points_on_seq: dict

    :param offset: The offset by which to shift the keys. Defaults to 0.
    :type offset: int

    :return: A dictionary with the keys shifted by the specified offset.
    :rtype: dict
    """

    keys = sorted(list(attachment_points_on_seq.keys()))[::-1]
    for attachment_point_id in keys:
        attachment_points_on_seq[
            attachment_point_id + offset
        ] = attachment_points_on_seq.pop(attachment_point_id)
    return attachment_points_on_seq


def translate_mod_smiles(smiles: str, offset=0):
    """
    Translates the atom map numbers in a SMILES string by a given offset.
    This function takes a SMILES string and increments the atom map numbers
    by the specified offset. It is useful for modifying the atom map numbers
    in a molecule represented by a SMILES string.

    :param smiles: The input SMILES string representing the molecule.
    :type smiles: str

    :param offset: The offset to add to the atom map numbers. Default is 0.
    :type offset: int

    :return: The modified SMILES string with updated atom map numbers.
    :rtype: str
    """

    mol = rdkit.Chem.MolFromSmiles(smiles)

    for atom in mol.GetAtoms():
        if ("dummyLabel" in atom.GetPropNames()) and (
            "molAtomMapNumber" in atom.GetPropNames()
        ):
            molAtomMapNumber = int(atom.GetProp("molAtomMapNumber"))
            molAtomMapNumber += offset
            atom.SetProp("molAtomMapNumber", str(molAtomMapNumber))
    smiles = rdkit.Chem.MolToSmiles(mol)
    return smiles


def translate_external_modification(mod: dict, offset=0) -> dict:
    """
    Translates the external modification dictionary by applying an offset to its SMILES string,
    attachment points on the sequence, and the maximum attachment point ID.

    :param mod: A dictionary representing the external modification with keys:
           - "smiles" (str): The SMILES string of the modification.
            - "attachment_points_on_sequence" (list): A list of attachment points on the sequence.
            - "max_attachment_point_id" (int): The maximum attachment point ID.
    :type mod: dict

    :param offset: The offset to apply to the modification. Default is 0.
    :type offset: int

    :return: A new dictionary with the translated modification data, containing:
            - "smiles" (str): The translated SMILES string.
            - "attachment_points_on_sequence" (list): The translated list of attachment points on the sequence.
            - "max_attachment_point_id" (int): The translated maximum attachment point ID.
    :rtype: dict
    """
    mod_smiles = mod.get("smiles")
    mod_smiles = translate_mod_smiles(mod_smiles, offset=offset)

    attachment_points_on_seq = mod.get("attachment_points_on_sequence")
    attachment_points_on_seq = translate_attachment_points_on_seq(
        attachment_points_on_seq, offset=offset
    )

    max_attachment_point_id = mod.get("max_attachment_point_id")
    max_attachment_point_id += offset

    mod = {
        "smiles": mod_smiles,
        "max_attachment_point_id": max_attachment_point_id,
        "attachment_points_on_sequence": attachment_points_on_seq,
    }
    return mod


def sequence_dict_to_string(sequence_dict: dict) -> str:
    """
    Converts a dictionary of sequence residues to a string representation.

    :param sequence_dict: A dictionary where keys are residue IDs and values are residue symbols.
    :type sequence_dict: dict

    :return: A string representation of the sequence, where each residue symbol
     is concatenated in the order of sorted residue IDs.
             If a residue symbol has more than one character, it is enclosed in curly braces.

    :rtype: str
    """

    sequence_string = ""
    for res_id in sorted(sequence_dict.keys()):
        symbol = sequence_dict.get(res_id)
        if len(symbol) > 1:
            symbol = "{%s}" % symbol
        sequence_string = sequence_string + symbol
    return sequence_string


def decompose_residues_internal(
    fragments: list, cx_smarts_db: dict, n_subst_limit=None
) -> tuple:
    """
    Decomposes a list of molecular fragments into their constituent residues and modifications.

    :param fragments: A list of molecular fragments to be decom
    :type fragments: list

    :param cx_smarts_db: A dictionary containing SMARTS patterns for chemical substructures.
    :type cx_smarts_db: dict

    :param n_subst_limit: An optional limit on the number of substitutions. Defaults to None.
    :type n_subst_limit: int

    :return: A tuple containing:
        - sequence_string (str): A string representation of the sequence of residues.
        - internal_modifications (list): A list of internal modifications found in the fragments.
        - external_modifications (list): A list of external modifications found in the fragments.
    :rtype: tuple
    """

    sequence_dict = {}

    internal_modifications = []
    external_modifications = []

    for mol in [nx_to_mol(fragment) for fragment in fragments]:
        fragment_res_names, fragment_modifications = full_decomposition(
            mol, cx_smarts_db, n_subst_limit=n_subst_limit
        )
        internal_modifications += fragment_modifications["internal_modifications"]

        if fragment_modifications.get("external_modifications"):
            external_modifications += fragment_modifications.get(
                "external_modifications"
            )
        sequence_dict.update(
            {int(k): fragment_res_names[k] for k in fragment_res_names}
        )
    sequence_string = sequence_dict_to_string(sequence_dict)
    return sequence_string, internal_modifications, external_modifications
