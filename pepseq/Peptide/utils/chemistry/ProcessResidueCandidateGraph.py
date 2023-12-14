"""
#**pepseq.pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph**

Provide several sample math calculations.

This module allows the user to make mathematical calculations.


Examples:
    >>> from pepseq.commands import pepseq_to_smiles, calculate_json_from, read_smiles, augment_db_json_command
    >>> pepseq_to_smiles('CH3~SCAFC~NH2')
    
    >>> calculate_json_from('CH3~SC{R1}AFC~NH2', '[*1]CCC') calculations.multiply(2.0, 4.0)
    
    >>> read_smiles('mypeptide.smi', 'myppeptide_out')

    >>> augment_db_json_command('my_monomers.sdf', 'augmented_db.json')

The module contains the following functions:

- `pepseq_to_smiles(pepseq, out, db_path)` - Returns the SMILES code for pepseq.
- `calculate_json_from(sequence, mod_smiles, out, db_path)` - Returns JSON dict 
   containing info about amino acid sequence and sequence modifications
- `read_smiles(smiles_filename, out, db_path, v) - Reads SMILES from file  and
   Writes parsed Pepseq and its modifications into separate txt files.
- `augment_db_json_command(sdf_path, out)` - Reads additional (Modified) Peptide
   building blocks from SDF file. Adds them to database and outputs it to file.

"""


from tqdm import tqdm

import networkx as nx
import rdkit

from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import (mol_to_nx, nx_to_mol)
from augmenting_db_json import get_Nter_versions_cxsmarts_db

"""

"""


def get_match(mol=None, aa_smarts=None, useChirality=True):
    """
    Get the substructure matches of a molecule with a given amino acid SMARTS pattern.

    Args:
        mol (rdkit.Chem.Mol): The molecule to search for matches.
        aa_smarts (str): The SMARTS pattern of the amino acid.
        useChirality (bool, optional): Whether to consider chirality in the matching process. Defaults to True.

    Returns:
        tuple: A tuple containing the amino acid molecule and the substructure matches.

    """
    aa_mol = rdkit.Chem.MolFromSmarts(aa_smarts)
    matches = mol.GetSubstructMatches(aa_mol, useChirality=useChirality)
    return (aa_mol, matches)


def get_matches(mol: rdkit.Chem.rdchem.Mol, cx_smarts_db: dict, useChirality=True) -> dict:
    """
    Find matches of a molecule against a dictionary of SMARTS patterns. For each of the amino acid SMARTS codes present in database
    they are matched against FragmentMolecule. Each matched Amino Acid specie is returned as dictionary
    together with matched substructures


    Args:
        mol (rdkit.Chem.rdchem.Mol): The molecule to search for matches.
        cx_smarts_db (dict): The dictionary of SMARTS patterns for each amino acid.
        useChirality (bool, optional): Whether to consider chirality in the matching process. Defaults to True.

    Returns:
        dict: A dictionary containing the matched amino acids as keys and their corresponding molecule and match information as values.
            {'C': (C_mol: rdkit.Chem.rdchem.Mol, (1,2,3), (4,5,6))}

    """
    matches_dict = {}
    for aa in cx_smarts_db:
        aa_smarts = cx_smarts_db[aa]
        aa_mol, matches = get_match(mol=mol, aa_smarts=aa_smarts, useChirality=useChirality)

        if matches:
            matches_dict[aa] = (aa_mol, matches)
    return matches_dict


def get_match_cover(G: nx.classes.graph.Graph, ResID: str):
    """
    Calculates the match cover of a graph G for a given ResID.

    Parameters:
    G (nx.classes.graph.Graph): The graph to calculate the match cover for.
    ResID (str): The ResID to match.

    Returns:
    int: The match cover of the graph G for the given ResID.
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
    Matches a molecular graph to a residue ID and returns the best match. For each of the substructure matches 
    grouped by amino acid specie we filter the ones covering only one residue and return the match covering
    the biggest posrtion of that residue (i.e. Cysteine is preferred over Alaine) if Alanine and Cysteine are matched;
    Cysteine match will be returned.
 

    Args:
        G (nx.classes.graph.Graph): The molecular graph to match.
        ResID (str): The residue ID to match against. For example '1' (one of the many ResIDs determined by fitting
        Peptide Backbone)

        matches_dict (dict): A dictionary containing amino acid matches and their corresponding molecular graphs.
            {'C': (C_mol: rdkit.Chem.rdchem.Mol, (1,2,3), (4,5,6))}


    Returns:
        tuple: A tuple containing the best matching amino acid, its corresponding molecular graph, and the match itself.
            (max_aa: amino_acid specie with greatest cover over fragment portion including ResX (X=ResID)
            and only ResX,
            max_aa_mol: molecule substructure covered by amino acid specie: rdkit.Chem.rdchem.Mol
            max_match: e.g. (1,2,3) atom_ids covered by match
            )
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
    Get modification graphs from a given graph.

    Args:
        G (networkx.Graph): The input graph.
        native_atom_ids (list): List of native atom IDs.
        ResID (str): The residue ID.

    Returns:
        list: List of modification graphs.

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


def get_n_subst_dict(G, matches_dict):
    """
    Calculate the number of substitutions for each amino acid residue in the peptide graph.

    Args:
        G (networkx.Graph): The peptide graph.
        matches_dict (dict): A dictionary containing amino acid names as keys and tuples of amino acid molecule and matches as values.

    Returns:
        dict: A dictionary containing amino acid names as keys and the number of substitutions as values.
    """
    n_subst = {}
    ResID = 'ResIDVal'

    for aa_name in tqdm(matches_dict):
        aa_mol, matches = matches_dict.get(aa_name)
        native_atom_ids = set([])
        for match in matches:
            native_atom_ids |= set(match)
        G_copy = G.copy()
        mod_graphs = get_mod_graphs(G_copy, native_atom_ids, ResID)
        n_subst[aa_name] = len(mod_graphs)
    return n_subst


def filter_n_subst(G, matches_dict, n_subst_limit):
    """
    Filter the matches_dict based on the number of substitutions (n_subst_limit) allowed for each amino acid.

    Args:
        G (Graph): The graph representing the residue candidate.
        matches_dict (dict): A dictionary containing the matches for each amino acid.
        n_subst_limit (int): The maximum number of substitutions allowed for each amino acid.

    Returns:
        dict: A filtered dictionary containing the matches for amino acids that satisfy the n_subst_limit condition.
    """
    n_subst_dict = get_n_subst_dict(G, matches_dict)
    aa_names_n_subst_elt = [
        aa_name for aa_name in n_subst_dict if n_subst_dict.get(
        aa_name)  <= n_subst_limit]
    return {
        k: matches_dict.get(
            k) for k in matches_dict if k in aa_names_n_subst_elt}


def get_res_matches(mol: rdkit.Chem.rdchem.Mol, cx_smarts_db: dict, useChirality=True, n_subst_limit=None) -> dict:
    """

    We use nx.classes.graph.Graph representation of modified peptide molecules

    Input:

    mol - Fragment Molecule: rdkit.Chem.rdchem.Mol
    cx_smarts_db - dictionary containing SMARTS codes for each of the Amino Acid
    species


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

    if '1' in res_ids.values():
        
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
    """
    Returns a list of edges that connect nodes from nodes1 to nodes2 or vice versa.

    Parameters:
    edges (list): A list of edges.
    nodes1 (list): A list of nodes.
    nodes2 (list): Another list of nodes.

    Returns:
    list: A list of edges that connect nodes from nodes1 to nodes2 or vice versa.
    """
    return [
        edge
        for edge in edges
        if (edge[0] in nodes1 and edge[1] in nodes2)
        or (edge[0] in nodes2 and edge[1] in nodes1)
    ]


def get_connecting_edges(union_graph: nx.classes.graph.Graph, H_graph: nx.classes.graph.Graph,
                        I_graph: nx.classes.graph.Graph) -> list:
    """
    Returns a list of edges connecting nodes from H_graph to nodes from I_graph in the union_graph.

    Args:
        union_graph (nx.classes.graph.Graph): The union graph containing nodes and edges from both H_graph and I_graph.
        H_graph (nx.classes.graph.Graph): The first graph containing nodes.
        I_graph (nx.classes.graph.Graph): The second graph containing nodes.

    Returns:
        list: A list of edges connecting nodes from H_graph to nodes from I_graph.
    """
    edges = list(union_graph.edges)
    nodes1 = list(H_graph.nodes)
    nodes2 = list(I_graph.nodes)
    return get_connecting(edges, nodes1, nodes2)


def get_atom_pairs(atoms_1: set, atoms_2: set, edges: list) -> list:
    """
    Get atom pairs from a list of edges based on the given sets of atoms and edges.

    Args:
        atoms_1 (set): Set of atoms to consider for the first position of the pair.
        atoms_2 (set): Set of atoms to consider for the second position of the pair.
        edges (list): List of edges representing the connections between atoms.

    Returns:
        list: List of atom pairs, where each pair is represented as a tuple.

    Raises:
        KeyError: If no atom is found in the given sets for an edge.

    """
    atom_pairs = []
    for edge in edges:
        atom1 = (set(edge) & atoms_1).pop()
        atom2 = (set(edge) & atoms_2).pop()
        tup = (atom1, atom2)
        atom_pairs.append(tup)
    return atom_pairs


def process_internal_connections(connections, res_matches, G: nx.classes.graph.Graph) -> list:
    """
    Process the internal connections between residues in a peptide sequence.

    Args:
        connections (list): List of tuples representing the connections between residues.
        res_matches (dict): Dictionary containing the matches of residues.
        G (nx.classes.graph.Graph): Graph representing the peptide sequence.

    Returns:
        list: List of dictionaries representing the internal bonds between residues.
    """
    bond_jsons = []
    bond_id = 0

    for name1, name2, connecting_edges in connections:
        res_ids = [ name.split("_")[:2][1]  for name in (name1, name2) ]

        res_atoms_1, res_atoms_2 = [set(res_matches[res_id][2]) for res_id in res_ids]

        bonds = get_atom_pairs( res_atoms_1, res_atoms_2, connecting_edges)

        for bond in bonds:
            
            bond_id += 1

            atom_names = [G.nodes[atom].get("AtomName") for atom in bond]

            bond_json = {
                bond_id: [
                    {
                        "ResID": res_ids[ind],
                        "AtomName": atom_names[ind],
                        "ResidueName": "",
                    } for ind in range(2)
                ]
            }

            bond_jsons.append(bond_json)
    return bond_jsons


def sorted_connection(connection):
    """
    Sorts the connection list based on the residue IDs.

    Args:
        connection (list): The connection list to be sorted.

    Returns:
        list: The sorted connection list.
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
    """
    Adds a connection point to the molecular graph.

    Args:
        G (nx.classes.graph.Graph): The molecular graph.
        point_id (int): The ID of the connection point.
        atom_id (int): The ID of the atom to connect the connection point to.

    Returns:
        nx.classes.graph.Graph: The updated molecular graph.
    """
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
    """
    Get the residue ID and atoms for a given residue name.

    Args:
        res_matches (dict): A dictionary containing residue matches.
        res_name (str): The name of the residue.

    Returns:
        tuple: A tuple containing the residue ID and a set of atoms.

    """
    list_res = res_name.split("_")

    res_id = list_res[1]

    res_atoms = set(res_matches[res_id][2])

    return res_id, res_atoms


def process_external_modification(G_mod: nx.classes.graph.Graph, mod_bonds, res_matches, atom_names_dict, point_id):
    """
    Process the external modification by adding attachment points to the molecular graph.
    For each modification (like staple) we need unambigious info.

    Args:
        G_mod (nx.classes.graph.Graph): The molecular graph representing the modification. (like staple)
        mod_bonds: The bonds associated with the modification.
        res_matches: The matches of residues in the modification.
        atom_names_dict: A dictionary mapping residue atoms to their names.
        point_id: The starting point ID for attachment points.

    Returns:
        Tuple: A tuple containing the modified JSON representation of the modification and the updated point ID.
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
    """
    Process the external connections of a peptide.

    Args:
        modifications (dict): A dictionary of modifications.
        res_matches (list): A list of residue matches.
        modification_graphs (list): A list of modification graphs.
        G (nx.classes.graph.Graph): A graph representing the peptide.

    Returns:
        list: A list of external modifications.

    """
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
    """
    Split the connections into internal and external bonds based on the type of residues involved.

    Args:
        connections (list): A list of tuples representing the connections between residues.

    Returns:
        tuple: A tuple containing two lists. The first list contains internal bonds, and the second list contains external bonds.

    """
    internal_bonds = []
    external_bonds_dict = {}

    for name1, name2, bonds in connections:
        types = [name.split('_')[0] for name in (name1, name2)]
        if set(types) == set(['Res',]):
            internal_bonds.append( ( name1, name2, bonds,) )
        else:
            mod_id = name2.split("_")[1]
            if external_bonds_dict.get(mod_id) is None:
                external_bonds_dict[mod_id] = []
            external_bonds_dict[mod_id].append( ( name1, bonds, )
            )
    return internal_bonds, external_bonds_dict


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


def decompose(mol: rdkit.Chem.rdchem.Mol, cx_smarts_db: dict, n_subst_limit=None) -> tuple:
    """
    Decomposes a molecule into its constituent parts based on the provided SMARTS database.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The molecule to be decomposed.
        cx_smarts_db (dict): The SMARTS database used for decomposition.
        n_subst_limit (int, optional): The maximum number of substitutions allowed. Defaults to None.

    Returns:
        tuple: A tuple containing the decomposed molecular graph, the residue matches, and the modification graphs.
    """
    res_matches = get_res_matches(mol, cx_smarts_db, n_subst_limit=n_subst_limit)
    G = mol_to_nx(mol)
    
    G = propagate_matches_on_molecular_graph(G, res_matches)
    modification_graphs = get_modification_graphs_from_fragment(G)
    return G, res_matches, modification_graphs


def get_internal_connections_subgraph_tuples(G: nx.classes.graph.Graph, res_matches: dict) -> list:
    """
    Input:

    G - molecular graph representing fragment molecule
    Output:

    subgraph_tuples - [
        ( 'Res_1_C', G_residue1_molecular_subgraph ),
        ( 'Res_2_X', G_residue2_molecular_subgraph ),
        ]
    """
    subgraph_tuples = []

    for res_id in sorted(res_matches.keys()):
        res_name, res_atoms_dict, res_atoms = res_matches[res_id]
        res_tuple = ("Res_%s_%s" % (res_id, res_name),  G.subgraph(res_atoms) )
        subgraph_tuples.append(res_tuple)

    return subgraph_tuples


def get_subgraph_tuples(res_matches: dict, modification_graphs_nodes, G: nx.classes.graph.Graph):
    """
    Get subgraph tuples for internal and external connections.

    Args:
        res_matches (dict): A dictionary containing residue matches.
        modification_graphs_nodes: A list of modification graph nodes.
        G (nx.classes.graph.Graph): A graph object. molecular graph representing fragment molecule

    Returns:
        list: A list of subgraph tuples for internal and external connections.
              subgraph_tuples - [
              ( 'Res_1_C', G_residue1_molecular_subgraph ),
              ( 'Res_2_X', G_residue2_molecular_subgraph ),
              ( 'Mod_1_X', G_modification1_molecular_subgraph ),
              ]
    """
    internal_subgraph_tuples = get_internal_connections_subgraph_tuples(G, res_matches)

    external_subgraph_tuples = [("Mod_%s" % (i + 1), modification_graphs_nodes[i]
                                 ) for i in range( len(modification_graphs_nodes) )]
    return internal_subgraph_tuples + external_subgraph_tuples


def get_connections(G_edges: list, subgraph_tuples: list):
    """
    Returns a list of bonds connecting subgraphs in a graph.
    For each pair of Residue Subgraphs We find whether in Fragment MolecularGraph edges there is a connecting
    edge between Residue1 and Residue2.


    Args:
        G_edges (list): List of edges in the graph. G_edges - [
            (Res1_Node1, Res1_Node2),
            (Res1_Node1, Res2_Node2),
            ... 
            ]

        subgraph_tuples (list): List of tuples containing the name and nodes of each subgraph.
            subgraph_tuples - [(ResName_1, SubGraphRes_1),
            (ResName_2, SubGraphRes_2),
            ...
            ]


    Returns:
        list: List of tuples representing the bonds connecting subgraphs.
              Each tuple contains the names of the connected subgraphs and the edges connecting them.
              connections - [(res1, res2, [(node1_1, node2_1), (node1_2, node2_34)])]
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


def full_decomposition(mol:rdkit.Chem.rdchem.Mol, cx_smarts_db: dict, n_subst_limit=None):
    """
    Perform a full decomposition of a molecule into its constituent residues and modifications.

    Args:
        mol (rdkit.Chem.rdchem.Mol): The molecule to be decomposed.
        cx_smarts_db (dict): A dictionary containing the SMARTS patterns for each residue.
        n_subst_limit (int, optional): The maximum number of substitutions allowed. Defaults to None.

    Returns:
        tuple: A tuple containing the residue names and the internal and external modifications.
            - res_names (dict): A dictionary mapping residue IDs to residue names.
            - modifications (dict): A dictionary containing the internal and external modifications.
                - internal_modifications (dict): A dictionary containing the internal modifications.
                - external_modifications (dict): A dictionary containing the external modifications.
    """
    G, res_matches, modification_graphs = decompose(mol, cx_smarts_db, n_subst_limit=n_subst_limit)

    modification_graphs_nodes = [list(graph.nodes) for graph in modification_graphs]

    subgraph_tuples = get_subgraph_tuples(res_matches, modification_graphs_nodes, G)
    bonds = get_connections( list(G.edges),  subgraph_tuples)

    res_res_connections, external_mod_connections = split_connections_by_type(
        bonds
    )
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
    Translates the attachment points on the sequence by applying an offset.

    Args:
        attachment_points_on_seq (dict): A dictionary containing the attachment points on the sequence.
        offset (int, optional): The offset to be applied to the attachment points. Defaults to 0.

    Returns:
        dict: A dictionary with the translated attachment points.
    """
    keys = sorted(list(attachment_points_on_seq.keys()))[::-1]
    for attachment_point_id in keys:
        attachment_points_on_seq[
            attachment_point_id + offset
        ] = attachment_points_on_seq.pop(attachment_point_id)
    return attachment_points_on_seq


def translate_mod_smiles(smiles: str, offset=0):
    """
    Translates the SMILES representation of a molecule with modified atoms by applying an offset to the isotope number of each atom.

    Args:
        smiles (str): The SMILES representation of the molecule.
        offset (int, optional): The offset to be applied to the isotope number of each atom. Defaults to 0.

    Returns:
        str: The translated SMILES representation of the molecule.
    """
    mol = rdkit.Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            isotope = atom.GetIsotope()
            isotope += offset
            atom.SetIsotope(isotope)
    smiles = rdkit.Chem.MolToSmiles(mol)
    return smiles


def translate_external_modification(mod: dict, offset=0) -> dict:
    """
    Translates an external modification dictionary to a modified dictionary with updated values.

    Args:
        mod (dict): The external modification dictionary.
        offset (int, optional): The offset value to be added. Defaults to 0.

    Returns:
        dict: The modified dictionary with updated values.
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
    Converts a dictionary representing a sequence into a string.

    Args:
        sequence_dict (dict): A dictionary where the keys are residue IDs and the values are symbols.

    Returns:
        str: The sequence string.

    """
    sequence_string = ""
    for res_id in sorted(sequence_dict.keys()):
        symbol = sequence_dict.get(res_id)
        if len(symbol) > 1:
            symbol = "{%s}" % symbol
        sequence_string = sequence_string + symbol
    return sequence_string


def decompose_residues_internal(fragments: list, cx_smarts_db: dict, n_subst_limit=None) -> tuple:
    """
    Decomposes the residues in the given fragments into individual components.

    Args:
        fragments (list): A list of fragments (molecules) to decompose.
        cx_smarts_db (dict): A dictionary containing the SMARTS patterns for decomposition.
        n_subst_limit (int, optional): The maximum number of substitutions allowed. Defaults to None.

    Returns:
        tuple: A tuple containing the sequence string, internal modifications, and external modifications.
    """
    sequence_dict = {}

    internal_modifications = []
    external_modifications = []

    for mol in [nx_to_mol(fragment) for fragment in fragments]:
        fragment_res_names, fragment_modifications = full_decomposition(mol, cx_smarts_db,
                                                            n_subst_limit=n_subst_limit)
        internal_modifications += fragment_modifications["internal_modifications"]

        if fragment_modifications.get("external_modifications"):
            external_modifications += fragment_modifications.get("external_modifications")
        sequence_dict.update({int(k): fragment_res_names[k] for k in fragment_res_names} )
    sequence_string = sequence_dict_to_string(sequence_dict)
    return sequence_string, internal_modifications, external_modifications
