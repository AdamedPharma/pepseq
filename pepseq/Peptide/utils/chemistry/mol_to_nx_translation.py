import copy
import networkx as nx
import rdkit
import rdkit.Chem
import matplotlib.pyplot as plt


def get_edge_tuple(edge: tuple) -> tuple:
    """
    returns info on graph edge (representing
    molecule bond) in particular on type
    and whether the bond is a peptide
    backbone bond
    """
    bond_start, bond_end, data = edge
    bond_type = data.get("bond_type")
    if bond_type is not None:
        bond_type = int(bond_type)
    is_peptide_bond = data.get("is_peptide_bond")
    return bond_type, is_peptide_bond, bond_start, bond_end


def mol_to_nx(mol: rdkit.Chem.rdchem.Mol) -> nx.classes.graph.Graph:
    """
    transforms rdkit.Chem.rdchem.Mol molecule
    into nx.classes.graph.Graph graph
    """

    G = nx.Graph()

    for atom in mol.GetAtoms():
        kwargs = {}

        for prop in atom.GetPropNames():
            prop_val = atom.GetProp(prop)
            kwargs[prop] = prop_val

        map_num = atom.GetAtomMapNum()
        if map_num is not None:
            if map_num != 0:
                kwargs["molAtomMapNumber"] = map_num

        if kwargs.get("num_explicit_hs") is not None:
            kwargs.pop("num_explicit_hs")

        if kwargs.get("hybridization") is not None:
            kwargs.pop("hybridization")

        G.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            formal_charge=atom.GetFormalCharge(),
            chiral_tag=atom.GetChiralTag(),
            hybridization=atom.GetHybridization(),
            num_explicit_hs=atom.GetNumExplicitHs(),
            is_aromatic=atom.GetIsAromatic(),
            isotope=atom.GetIsotope(),
            **copy.deepcopy(kwargs)
        )

    for bond in mol.GetBonds():
        kwargs = {}
        for prop in bond.GetPropNames():
            prop_val = bond.GetProp(prop)
            kwargs[prop] = prop_val

        G.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_type=bond.GetBondType(),
            **kwargs
        )
    return G


def nx_to_mol(G: nx.classes.graph.Graph) -> rdkit.Chem.rdchem.Mol:
    """
    transforms nx.classes.graph.Graph graph
    into  rdkit.Chem.rdchem.Mol molecule
    """

    mol = rdkit.Chem.RWMol()
    atomic_nums = nx.get_node_attributes(G, "atomic_num")
    chiral_tags = nx.get_node_attributes(G, "chiral_tag")
    formal_charges = nx.get_node_attributes(G, "formal_charge")
    node_is_aromatics = nx.get_node_attributes(G, "is_aromatic")
    node_hybridizations = nx.get_node_attributes(G, "hybridization")
    num_explicit_hss = nx.get_node_attributes(G, "num_explicit_hs")

    basic_propNames = set(
        [
            "atomic_num",
            "chiral_tag",
            "formal_charge",
            "is_aromatic",
            "hybridization" "num_explicit_hs",
            "isotope",
        ]
    )

    node_to_idx = {}
    for node in G.nodes():
        a = rdkit.Chem.Atom(atomic_nums[node])

        additional_propNames = G.nodes[node].keys() - basic_propNames
        for propName in additional_propNames:
            a.SetProp(propName, str(G.nodes[node][propName]))

        a.SetChiralTag(chiral_tags[node])
        a.SetFormalCharge(formal_charges[node])
        a.SetIsAromatic(node_is_aromatics[node])
        if node_hybridizations[node] is not None:
            a.SetHybridization(node_hybridizations[node])
        a.SetNumExplicitHs(num_explicit_hss[node])
        a.SetIsotope(0)

        idx = mol.AddAtom(a)
        node_to_idx[node] = idx

    basic_propNames = set(["bond_type"])

    bond_types = nx.get_edge_attributes(G, "bond_type")
    for edge in G.edges():
        first, second = edge
        ifirst = node_to_idx[first]
        isecond = node_to_idx[second]
        bond_type = bond_types[first, second]
        mol.AddBond(ifirst, isecond, bond_type)

    rdkit.Chem.SanitizeMol(mol)
    return mol


def get_chiral_tag_int(chiral_tag):
    if chiral_tag == rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
        chiral_tag_int = 0
    elif chiral_tag == rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
        chiral_tag_int = 1
    else:
        chiral_tag_int = -1
    return chiral_tag_int


def get_hybridization_int(hybridization):
    if hybridization == rdkit.Chem.rdchem.HybridizationType.SP3:
        return 4
    elif hybridization == rdkit.Chem.rdchem.HybridizationType.SP2:
        return 3
    elif hybridization == rdkit.Chem.rdchem.HybridizationType.SP:
        return 2
    elif hybridization == rdkit.Chem.rdchem.HybridizationType.S:
        return 1


def get_node_tuple(node_data, keys=None):
    idx, node_dict = node_data
    if keys is None:
        keys = [
            "atomic_num",
            "formal_charge",
            "chiral_tag",
            "hybridization",
            "num_explicit_hs",
            "is_aromatic",
            "isotope",
            "AtomName",
            "ResID",
        ]
    row = []
    for key in keys:
        if key == "chiral_tag":
            val = get_chiral_tag_int(node_dict.get("chiral_tag"))
        elif key == "hybridization":
            val = get_hybridization_int(node_dict.get("hybridization"))
        else:
            val = node_dict.get(key)
        row.append(val)
    row.append(idx)
    return row


def nx_to_json(G: nx.classes.graph.Graph) -> dict:
    """
    transforms a networkx  Graph
    into JSON dict in the form

    :parameter G nx.classes.graph.Graph: networkx  Graph

    :return: JSON dict in the form

    """
    nodes_list = list(G.nodes(data=True))

    node_dicts = []
    for node_tuple in nodes_list:
        node_attrs = node_tuple[1]
        node_dicts.append(copy.deepcopy(node_attrs))

    nodekeys = set()
    for i in node_dicts:
        nodekeys |= set(i.keys())

    nodes_columns = sorted(list(nodekeys))

    node_rows = []
    for node_data in nodes_list:
        node_tuple = get_node_tuple(node_data, keys=nodes_columns)
        node_rows.append(node_tuple)

    edges_list = list(G.edges(data=True))
    edges_tuple = tuple([get_edge_tuple(edge) for edge in edges_list])

    edges_columns = ["bond_type", "is_peptide_bond", "bond_start", "bond_end"]
    nodes_columns += ["node_id"]

    mol_j = {
        "nodes_tuple": node_rows,
        "nodes_columns": nodes_columns,
        "edges_tuple": edges_tuple,
        "edges_columns": edges_columns,
    }
    return mol_j


def get_mol_json(mol: rdkit.Chem.rdchem.Mol) -> dict:
    """
    transforms rdkit.Chem.rdchem.Mol molecule object into JSON dict

    :parameter mol rdkit.Chem.rdchem.Mol: molecule object

    :return: JSON dict in the form

    """
    G = mol_to_nx(mol)
    mol_j = nx_to_json(G)
    return mol_j


def mol_json_to_nx(mol_json: dict) -> nx.classes.graph.Graph:
    """
    transforms JSON dict in the form into a networkx  Graph

    :parameter mol_json dict: JSON dict in the form

    :return: networkx  Graph

    """
    G = nx.Graph()

    nodes_tuple = mol_json["nodes_tuple"]
    nodes_columns = mol_json["nodes_columns"]

    node_dicts = []

    n_columns = len(nodes_columns)

    for node_tuple in nodes_tuple:
        node_dict = {}
        for i in range(n_columns):
            node_dict[nodes_columns[i]] = node_tuple[i]
        node_dicts.append(copy.deepcopy(node_dict))

    for node_dict in node_dicts:
        atomic_num = node_dict.get("atomic_num")
        formal_charge = node_dict.get("formal_charge")
        chiral_tag = node_dict.get("chiral_tag")
        hybridization = node_dict.get("hybridization")
        num_explicit_hs = node_dict.get("num_explicit_hs")
        is_aromatic = node_dict.get("is_aromatic")
        isotope = node_dict.get("isotope")
        AtomName = node_dict.get("AtomName")
        ResID = node_dict.get("ResID")
        node_id = node_dict.get("node_id")
        dummyLabel = node_dict.get("dummyLabel")
        molAtomMapNumber = node_dict.get("molAtomMapNumber")

        if (chiral_tag is None) or (chiral_tag == 0):
            chiral_tag = rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED
        elif chiral_tag == 1:
            chiral_tag = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
        elif chiral_tag == -1:
            chiral_tag = rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW

        if hybridization == 4:
            hybridization = rdkit.Chem.rdchem.HybridizationType.SP3
        elif hybridization == 3:
            hybridization = rdkit.Chem.rdchem.HybridizationType.SP2
        elif hybridization == 2:
            hybridization = rdkit.Chem.rdchem.HybridizationType.SP
        elif hybridization == 1:
            hybridization = rdkit.Chem.rdchem.HybridizationType.S

        kwargs = dict(
            atomic_num=atomic_num,
            formal_charge=formal_charge,
            chiral_tag=chiral_tag,
            hybridization=hybridization,
            num_explicit_hs=num_explicit_hs,
            is_aromatic=is_aromatic,
            isotope=isotope,
        )
        if dummyLabel is not None:
            kwargs["dummyLabel"] = dummyLabel
        if molAtomMapNumber is not None:
            kwargs["molAtomMapNumber"] = molAtomMapNumber

        if ResID is not None:
            kwargs["ResID"] = ResID
        if AtomName is not None:
            kwargs["AtomName"] = AtomName
        G.add_node(node_id, **kwargs)
    edges_tuple = mol_json["edges_tuple"]

    for edge_tuple in edges_tuple:
        bond_type, is_peptide_bond, bond_start, bond_end = edge_tuple

        if bond_type == 1:
            bond_type = rdkit.Chem.rdchem.BondType.SINGLE
        elif bond_type == 2:
            bond_type = rdkit.Chem.rdchem.BondType.DOUBLE
        elif bond_type == 12:
            bond_type = rdkit.Chem.rdchem.BondType.AROMATIC
            # to do include all possible bond types
        kwargs = {}
        if is_peptide_bond is not None:
            kwargs = {"is_peptide_bond": is_peptide_bond}
        G.add_edge(bond_start, bond_end, bond_type=bond_type, **kwargs)
    return G


def mol_json_to_mol(mol_json: dict) -> rdkit.Chem.rdchem.Mol:
    G = mol_json_to_nx(mol_json)
    mol = nx_to_mol(G)
    return mol


def draw_peptide_json(peptide_json: dict, out: str = "simple_path.png") -> nx.Graph:
    """
    Draws schema of a modified peptide using  networkx Graph drawing function networkx.draw.
    Nodes in the Graph represent Amino Acids in the sequence and external modifications.
    External modifications are labeled as SMILES code.
    Peptide bonds between the Amino Acids; attachment points for external modifications
    and disulfide bridges are represented as edges.
    """
    G = nx.Graph()

    sequence = peptide_json["sequence"]

    for i in range(len(sequence)):
        G.add_node(i, aa="%d%s" % (i + 1, sequence[i]), color="blue")
    for i in range(1, len(sequence)):
        G.add_edge(i - 1, i, color="black", weight=2)

    int_mods = peptide_json["internal_modifications"]

    for int_mod in int_mods:
        int_res_ids = [int(i["ResID"]) - 1 for i in list(int_mod.values())[0]]

        start_ss, end_ss = int_res_ids

        G.add_edge(start_ss, end_ss, color="orange", weight=1)

    ext_mods = peptide_json["external_modifications"]

    for i in range(len(ext_mods)):
        ext_mod = ext_mods[i]
        ext_node_name = "mod_%d" % i
        G.add_node(ext_node_name, color="red", aa=ext_mod.get("smiles"))
        att_points = [
            int(i["ResID"]) - 1
            for i in list(ext_mod["attachment_points_on_sequence"].values())
        ]
        for att_point_res_id in att_points:
            G.add_edge(att_point_res_id, ext_node_name, color="green", weight=1)

    colors = nx.get_edge_attributes(G, "color").values()
    weights = nx.get_edge_attributes(G, "weight").values()
    aas = nx.get_node_attributes(G, "aa")  # .values()
    node_colors = nx.get_node_attributes(G, "color").values()

    labeldict = aas

    nx.draw(
        G,
        edge_color=colors,
        with_labels=True,
        labels=labeldict,
        width=list(weights),
        node_color=node_colors,
    )

    plt.savefig(out)  # save as png
    return G
