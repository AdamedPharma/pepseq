import networkx as nx
import rdkit

from typing import Union
from pepseq.Peptide.exceptions import InvalidSequenceError, InvalidSymbolError
from pepseq.Peptide.utils.chemistry.cap_termini import (cap_C_terminus,
                                                        cap_N_terminus)
from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import (mol_to_nx,
                                                                  nx_to_mol)
from pepseq.Peptide.utils.chemistry.MonomerConnector import \
    get_molecule_from_list_of_residue_symbols
from pepseq.Peptide.utils.Parser import (find_termini, get_canonical,
                                         parse_canonical)

def add_internal_bond(G: nx.classes.graph.Graph, res1_id: int, atom_name_1: str,
                    res2_id: int, atom_name_2: str) -> nx.classes.graph.Graph:
    """
    Add internal bonds within peptide molecule from values defined in peptide JSON.
    :param G: Modified Peptide molecule as networkx nx.classes.graph.Graph
    :type G: nx.classes.graph.Graph
    :param res1_id: Index of first amino acid residue involved in bonding
    :type int
    :param atom_name_1: Name of atom in first amino acid residue involved in bonding
    :type str
    :param res2_id: Index of second amino acid residue involved in bonding
    :type int
    :param atom_name_2: Name of atom in second amino acid residue involved in bonding
    :type str
    :return: peptide molecule - Modified Peptide molecule as networkx nx.classes.graph.Graph
    :rtype: nx.classes.graph.Graph
    """
    
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


def add_disulfide_bond(G: nx.classes.graph.Graph, res1_id: int, res2_id: int) -> nx.classes.graph.Graph:
    """
    Add disulfide bond within peptide molecule from values defined in peptide JSON.
    Atom Name values are set to 'SG'
    :param G: Modified Peptide molecule as networkx nx.classes.graph.Graph
    :type G: nx.classes.graph.Graph
    :param res1_id: Index of first amino acid residue involved in bonding
    :type int
    :param res2_id: Index of second amino acid residue involved in bonding
    :type int
    :return: peptide molecule - Modified Peptide molecule as networkx nx.classes.graph.Graph
    :rtype: nx.classes.graph.Graph
    """

    return add_internal_bond(G, res1_id, "SG", res2_id, "SG")


def get_attachment_points(staple_graph: nx.classes.graph.Graph) -> tuple:
    """
    Get attachment points
    :param staple_graph: molecular graph representing molecular staple
    :type nx.classes.graph.Graph
    :return: tuple composed of staple_graph (molecular graph nx.classes.graph.Graph representing
    molecular staple with dummy Atoms removed) and dictionary representing atoms on staple that connect 
    to amino acid chain
    :rtype: nx.classes.graph.Graph
    """

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


def find_atom(G: nx.classes.graph.Graph, ResID, AtomName: str) -> str:
    """
    :param G: molecular graph
    :type  G: nx.classes.graph.Graph

    :param ResID: id of residue
    :type  ResID: int

    :param AtomName: name of atom
    :type  AtomName: str

    :return 

    """
    ResIDs = nx.get_node_attributes(G, "ResID")
    AtomNames = nx.get_node_attributes(G, "AtomName")
    atoms_with_right_name = [i for i in AtomNames if AtomNames[i] == AtomName]
    point_list = [atom_id for atom_id in atoms_with_right_name if ResIDs[atom_id] == ResID]
    return point_list.pop()


def add_staple(
    peptide_graph: nx.classes.graph.Graph, staple_graph: nx.classes.graph.Graph,
    peptide_attachment_points, prefix="staple_1_"
) -> nx.classes.graph.Graph:
    """
    Adds a staple to a peptide graph by connecting attachment points.

    Args:
        peptide_graph (nx.classes.graph.Graph): The graph representing the peptide.
        staple_graph (nx.classes.graph.Graph): The graph representing the staple.
        peptide_attachment_points: The attachment points on the peptide.
        prefix (str, optional): The prefix to be added to the staple attachment points. Defaults to "staple_1_".

    Returns:
        nx.classes.graph.Graph: The graph representing the peptide with the added staple.
    """
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


def get_peptide_json_from_sequence(sequence: str, db_json: dict) -> dict:
    """
    Constructs a peptide JSON object from a given peptide sequence and a database JSON.

    Args:
        sequence (str): The peptide sequence.
        db_json (dict): The database JSON containing information about the peptide.

    Returns:
        dict: The constructed peptide JSON object.
    """

    N_terminus, C_terminus, sequence_str_wo_termini = find_termini(sequence, db_json)

    peptide_json = {
        "sequence": sequence_str_wo_termini,
        "internal_modifications": [],
        "external_modifications": [],
        "N_terminus": N_terminus,
        "C_terminus": C_terminus,
    }

    return peptide_json

def get_residue_symbols_from_sequence(
    sequence: str, db_json: dict, N_terminus: Union[str, None] = None, C_terminus: Union[str, None] = None
) -> list:
    """
    Retrieves the residue symbols from a given peptide sequence.

    Args:
        sequence (str): The peptide sequence.
        db_json (dict): The database JSON containing the mapping of canonical sequences to residue symbols.
        N_terminus (str, None): The symbol representing the N-terminus of the peptide. Defaults to None.
        C_terminus (str, None): The symbol representing the C-terminus of the peptide. Defaults to None.

    Returns:
        list: The list of residue symbols.

    """
    canonical_sequence = get_canonical(sequence, db_json)
    symbols_list_w_termini = parse_canonical(canonical_sequence)
    if N_terminus is None:
        N_terminus = symbols_list_w_termini[0]
    if C_terminus is None:
        C_terminus = symbols_list_w_termini[-1]
    residue_symbols = symbols_list_w_termini[1:-1]
    return residue_symbols

class Sequence(object):
    """
    Represents a peptide sequence.

    Args:
        sequence (str): The amino acid sequence.
        N_terminus (str, optional): The symbol representing the N-terminus. Defaults to None.
        C_terminus (str, optional): The symbol representing the C-terminus. Defaults to None.
    """

    def __init__(self, sequence: str, N_terminus: Union[str, None] = None, C_terminus: Union[str, None] = None):
        self.sequence = sequence
        self.N_terminus = N_terminus
        self.C_terminus = C_terminus
        return
    
    def extract_residue_symbols(self, db_json: dict):
        """
        Extracts the residue symbols from the peptide sequence.

        Args:
            db_json (dict): The JSON containing the peptide sequence information.

        Returns:
            tuple: A tuple containing the residue symbols, N-terminus SMILES, and C-terminus SMILES.
        """
        canonical_sequence = get_canonical(self.sequence, db_json)
        symbols_list_w_termini = parse_canonical(canonical_sequence)
        if self.N_terminus is None:
            self.N_terminus = symbols_list_w_termini[0]
        if self.C_terminus is None:
            self.C_terminus = symbols_list_w_termini[-1]
        residue_symbols = symbols_list_w_termini[1:-1]
        if ('[' in self.N_terminus) and (']' in self.N_terminus):
            N_terminus_smiles = self.N_terminus[1:-1]
        else:
            N_terminus_smiles = None

        if ('[' in self.C_terminus) and (']' in self.C_terminus):
            C_terminus_smiles = self.C_terminus[1:-1]
        else:
            C_terminus_smiles = None
        if self.N_terminus == "H":
            self.N_terminus = "proton"

        return residue_symbols, N_terminus_smiles, C_terminus_smiles


def get_coding(db_json: dict) -> dict:
    """
    Retrieves the coding information from the given database JSON.

    Args:
        db_json (dict): The database JSON containing the coding information.

    Returns:
        dict: The coding information.

    """
    keys = [
        "l_proteogenic_3letter",
        "d_proteogenic_3letter",
        "d_proteogenic_2letter",
        "d_proteogenic_4letter",
        "modified_aa_codes",
    ]

    coding = db_json.get("coding").get("l_proteogenic_3letter")
    for key in keys:
        coding.update(db_json.get("coding").get(key))
    return coding
    

def get_molecule_from_sequence(sequence: str, db_json: dict, N_terminus=None,
                            C_terminus=None) -> rdkit.Chem.rdchem.Mol:
    """
    Get the molecule representation of a peptide sequence using a database JSON.

    Args:
        sequence (str): The peptide sequence.
        db_json (dict): The database JSON containing information about the peptide building blocks.
        N_terminus (str, optional): The N-terminus of the peptide. Defaults to None.
        C_terminus (str, optional): The C-terminus of the peptide. Defaults to None.

    Returns:
        rdkit.Chem.rdchem.Mol: The molecule representation of the peptide sequence.

    Raises:
        InvalidSymbolError: If any residue symbols in the sequence are not found in the database.
        InvalidSymbolError: If a residue symbol is not found in the database.
    """
    try:
        sequence_object = Sequence(sequence, N_terminus, C_terminus)
        residue_symbols, N_terminus_smiles, C_terminus_smiles = sequence_object.extract_residue_symbols(db_json)

        coding = get_coding(db_json)
        aa_smiles_dict = db_json["smiles"].get("aa")

        unique_residue_symbols = set(residue_symbols)
        unique_encoded_symbols = set(coding.keys())
        unique_aa = set(aa_smiles_dict.keys())

        db_symbols_404 = unique_residue_symbols - (unique_encoded_symbols|unique_aa)

        if db_symbols_404:
            raise InvalidSymbolError("Residue Symbols: %s not found in database." %', '.join(list(db_symbols_404)))

        for i in range(len(residue_symbols)):
            symbol = residue_symbols[i]
            if aa_smiles_dict.get(symbol) is None:
                residue_symbols[i] = coding.get(symbol)

        smiles_building_blocks_db = {}

        for residue_symbol in residue_symbols:
            residue_db_entry = aa_smiles_dict.get(residue_symbol)
            if residue_db_entry is None:
                error_msg = "Residue Symbol: %s not found in database." % residue_symbol
                raise InvalidSymbolError(error_msg)

            smiles_building_blocks_db[residue_symbol] = residue_db_entry["smiles_radical"]

        N_terminus = sequence_object.N_terminus
        C_terminus = sequence_object.C_terminus

        if N_terminus in db_json["smiles"]["n_terms"]:
            smiles_building_blocks_db[N_terminus] = db_json["smiles"]["n_terms"].get(
                N_terminus
            )["smiles_radical"]

        if C_terminus in db_json["smiles"]["c_terms"]:
            smiles_building_blocks_db[C_terminus] = db_json["smiles"]["c_terms"][C_terminus][
                "smiles_radical"
            ]

        mol = get_molecule_from_list_of_residue_symbols(
            residue_symbols, smiles_building_blocks_db
        )

        mol_w_n_terminus = cap_N_terminus(
            mol, terminus=N_terminus, smiles_building_blocks_db=smiles_building_blocks_db,
            terminus_smiles=N_terminus_smiles
        )

        mol_w_nc_terminus = cap_C_terminus(
            mol_w_n_terminus,
            terminus=C_terminus,
            smiles_building_blocks_db=smiles_building_blocks_db,
            terminus_smiles=C_terminus_smiles
        )

    except InvalidSymbolError as exc:
        raise InvalidSymbolError(exc.msg)

    """
    except Exception:

        raise InvalidSequenceError("Invalid Sequence %s" %(str(residue_symbols)))
    """
    return mol_w_nc_terminus


def get_smiles_from_sequence(sequence: str, db_json: dict) -> str:
    """
    Get the SMILES representation of a peptide sequence.

    Args:
        sequence (str): The peptide sequence.
        db_json (dict): The database JSON containing the molecule information.

    Returns:
        str: The SMILES representation of the peptide.

    """
    mol = get_molecule_from_sequence(
        sequence, db_json, N_terminus=None, C_terminus=None
    )
    smiles = rdkit.Chem.MolToSmiles(mol)
    return smiles


def get_molecule_from_json(j: dict, db_json: dict) -> rdkit.Chem.rdchem.Mol:
    """
    Get a molecule from a peptide JSON representation.

    Args:
        j (dict): The peptide JSON representation.
        db_json (dict): The database JSON representation.

    Returns:
        rdkit.Chem.rdchem.Mol: The molecule representation of the peptide.
    """
    sequence = j["sequence"]
    N_terminus = j.get("N_terminus")
    C_terminus = j.get("C_terminus")

    return get_molecule_from_sequence(sequence, db_json, N_terminus, C_terminus)


class BuildingModifiedPeptideFromPeptideJSON(object):
    """
    Class for building a modified peptide from a peptide JSON and a database JSON.
    sequence is constructed (automatically or through connecting CXSMILES
    from database)
    Thus atoms serving as AttachmentPoints are marked
    For each ExternalModification edge is added between SequenceAttachmentPoint
    and ExternalModificationAttachmentPoint
    For each InternalModification edge is added between SequenceAttachmentPoint
    and matching SequenceAttachmentPoint

    """

    def __init__(self):
        return

    def execute(self, peptide_json: dict, db_json: dict) -> rdkit.Chem.rdchem.Mol:
        """
        Builds a modified peptide molecule from a peptide JSON and a database JSON.

        Args:
            peptide_json (dict): The peptide JSON containing information about the peptide.
            db_json (dict): The database JSON containing information about the database.

        Returns:
            rdkit.Chem.rdchem.Mol: The modified peptide molecule.
        """
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


def get_smiles_from_peptide_json(peptide_json: dict, db_json: dict) -> str:
    """
    Get the SMILES representation of a peptide from a peptide JSON and a database JSON.

    Args:
        peptide_json (dict): The JSON representation of the peptide.
        db_json (dict): The JSON representation of the database.

    Returns:
        str: The SMILES representation of the peptide.
    """
    mol = BuildingModifiedPeptideFromPeptideJSON().execute(peptide_json, db_json)
    smiles = rdkit.Chem.MolToSmiles(mol)
    return smiles
