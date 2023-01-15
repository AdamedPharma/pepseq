import networkx as nx
import rdkit
from rdkit import Chem


def get_peptide_bonds(n_ca_co: list) -> list:
    """

    Input:
        n_ca_co - list of [N, Ca, Co ... N, Ca, Co]
        atoms

    Output:
        peptide_bonds

    """

    peptide_bonds = []

    for i in range(0, len(n_ca_co), 3)[:-1]:
        co_atom = i + 2
        n_atom = i + 3
        edge = (n_ca_co[co_atom], n_ca_co[n_atom])
        peptide_bonds.append(edge)

    return peptide_bonds


def sg_list(G: nx.classes.graph.Graph):
    """ """
    g = (G.subgraph(c) for c in nx.connected_components(G))
    return list(g)


def break_residues(G: nx.classes.graph.Graph, n_ca_co: list):
    """
    Input:
        G - graph of molecule
        n_ca_co list of atoms N,Ca,Co backbone atoms

    Output:
        ordered_subgraphs list of nx.classes.graph.Graph
        representing peptide residues

    Graph of peptide molecule is 'broken'
    at peptide bonds
    into subgraphs for individual residues
    subgraphs are then ordered from
    N terminus to C terminus

    """
    peptide_bonds = get_peptide_bonds(n_ca_co)
    G_copy = G.copy()
    G_copy.remove_edges_from(peptide_bonds)
    res_subgraphs = sg_list(G_copy)
    ordered_subgraphs = []
    atoms = [bond[0] for bond in peptide_bonds]
    atoms.append(peptide_bonds[-1][1])
    for atom in atoms:
        for sg in res_subgraphs:
            if atom in sg.nodes:
                ordered_subgraphs.append(sg)
    return ordered_subgraphs


def mol_to_nx(mol: rdkit.Chem.rdchem.Mol) -> nx.classes.graph.Graph:
    """
    transforms rdkit.Chem.rdchem.Mol molecule
    into nx.classes.graph.Graph graph
    """

    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            formal_charge=atom.GetFormalCharge(),
            chiral_tag=atom.GetChiralTag(),
            hybridization=atom.GetHybridization(),
            num_explicit_hs=atom.GetNumExplicitHs(),
            is_aromatic=atom.GetIsAromatic(),
        )
    for bond in mol.GetBonds():
        G.add_edge(
            bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond_type=bond.GetBondType()
        )
    return G


def nx_to_mol(G: nx.classes.graph.Graph) -> rdkit.Chem.rdchem.Mol:
    """
    transforms nx.classes.graph.Graph graph
    into  rdkit.Chem.rdchem.Mol molecule
    """

    mol = Chem.RWMol()
    atomic_nums = nx.get_node_attributes(G, "atomic_num")
    chiral_tags = nx.get_node_attributes(G, "chiral_tag")
    formal_charges = nx.get_node_attributes(G, "formal_charge")
    node_is_aromatics = nx.get_node_attributes(G, "is_aromatic")
    node_hybridizations = nx.get_node_attributes(G, "hybridization")
    num_explicit_hss = nx.get_node_attributes(G, "num_explicit_hs")
    node_to_idx = {}
    for node in G.nodes():
        a = Chem.Atom(atomic_nums[node])
        a.SetChiralTag(chiral_tags[node])
        a.SetFormalCharge(formal_charges[node])
        a.SetIsAromatic(node_is_aromatics[node])
        a.SetHybridization(node_hybridizations[node])
        a.SetNumExplicitHs(num_explicit_hss[node])
        idx = mol.AddAtom(a)
        node_to_idx[node] = idx

    bond_types = nx.get_edge_attributes(G, "bond_type")
    for edge in G.edges():
        first, second = edge
        ifirst = node_to_idx[first]
        isecond = node_to_idx[second]
        bond_type = bond_types[first, second]
        mol.AddBond(ifirst, isecond, bond_type)

    Chem.SanitizeMol(mol)
    return mol


def do_all(smiles: str, validate=False) -> nx.classes.graph.Graph:
    """
    reads SMILES string into
        nx.classes.graph.Graph graph
        with optional validation
    """
    mol = Chem.MolFromSmiles(smiles.strip())
    can_smi = Chem.MolToSmiles(mol)
    G = mol_to_nx(mol)
    if validate:
        mol = nx_to_mol(G)
        new_smi = Chem.MolToSmiles(mol)
        assert new_smi == smiles
    return G


class GraphExtended(object):
    def __init__(self, G: nx.classes.graph.Graph):
        self.G = G
        return


def MolToGraphExtended(mol: rdkit.Chem.rdchem.Mol) -> GraphExtended:
    G = mol_to_nx(mol)
    return GraphExtended(G)


def SequenceToGraphExtended(sequence: str) -> GraphExtended:
    mol = rdkit.Chem.MolFromSequence(sequence)
    G = mol_to_nx(mol)
    return GraphExtended(G)


def generate_bb_smiles(num_res: int, OH=False) -> str:
    """
    generates SMILES for polyGlycine peptide
    length (number of Gly residues) is given by num_res
        parameter
    """
    if OH:
        bb_smiles = "O" + "C(=O)CN" * num_res  # generate backbone SMILES
    else:
        bb_smiles = "C(=O)CN" * num_res  # generate backbone SMILES

    return bb_smiles


def generate_bb_mol(num_res: int, OH=False) -> rdkit.Chem.rdchem.Mol:
    """
    generates rdkit.Chem.rdchem.Mol molecule for polyGlycine peptide
    length (number of Gly residues) is given by num_res
        parameter

    """
    bbsmiles = generate_bb_smiles(num_res, OH=OH)
    bbmol = Chem.MolFromSmiles(bbsmiles)
    return bbmol


def longest_backbone(mol: rdkit.Chem.rdchem.Mol) -> list:
    """
    iteratively
        generates rdkit.Chem.Molecule(s)
        representing glycine chain (protein backbone)
        of increasing length (incrementing by 1 every iteration)
        and matches it to the self.Molecule
        when after n iterations generated GlycineChain is no longer
        a substructure of self.Molecule
        iteration loop breaks
        and  match of glycine chain generated in previous iteration
        to molecule is returned
    """

    matches = [None]
    num_res = 0
    while matches:
        bb = matches[0]
        num_res += 1
        bbmol = generate_bb_mol(num_res, OH=False)
        matches = mol.GetSubstructMatches(bbmol)
    return bb


def get_N_CA_CO(G, backbone_atoms: list, double_bond_type) -> list:
    """
    Input:
         backbone atoms: List of atom ids
         double_bond_type: rdkit.Chem.rdchem.BondType.DOUBLE
         was put as variable to avoid rdkit dependency
         which it now fails to do so i might reconsider

    Output:
        an ordered list of backbone N and C atoms
        Nterm-Ca-Co- ... -N-Ca-Coterm atoms

    """
    nc_atomic_numbers = set([6, 7])

    atomic_numbers = nx.get_node_attributes(G, "atomic_num")
    CNs = set([i for i in atomic_numbers if atomic_numbers[i] in nc_atomic_numbers])

    backbone_CNs = set(backbone_atoms) & CNs
    backbone_CNs_subgraph = G.subgraph(backbone_CNs)

    degrees = dict(backbone_CNs_subgraph.degree)
    backbone_termini = [node for node in degrees if degrees[node] == 1]

    bond_types = nx.get_edge_attributes(G, "bond_type")

    terminus0 = backbone_termini[0]
    terminus_bonds = G.edges(terminus0)

    for bond in terminus_bonds:
        if bond_types[tuple(sorted(bond))] == double_bond_type:
            backbone_termini = [backbone_termini[1], backbone_termini[0]]

    n_ca_co_path = nx.shortest_path(
        backbone_CNs_subgraph, backbone_termini[0], backbone_termini[1]
    )
    return n_ca_co_path


def get_next_item(item):
    if item is None:
        return "N"

    next_item = {"N": "CA", "CA": "CO", "CO": "N"}

    return next_item[item]


def get_labels(
    G: nx.classes.graph.Graph,
    n_ca_co: list,
) -> nx.classes.graph.Graph:
    """
    backbone_atoms are given new attributes:

    {
        'label': N/ CA/ CO,
        'res_id': id_of_residue_in_chain
        }

    """

    position = 0
    label = None
    res_id = 0
    attrs = {}

    while position < len(n_ca_co):
        label = get_next_item(label)
        if label == "N":
            res_id += 1
        atom_id = n_ca_co[position]

        attrs[atom_id] = {"label": label, "res_id": res_id}

        position += 1
    nx.set_node_attributes(G, attrs)
    return G


def append_match_to_atom(
    G: nx.classes.graph.Graph, atom_id: int, match: tuple
) -> nx.classes.graph.Graph:
    """
    map onto node/atom what amino acids have been matched on it
    """
    matched_with = nx.get_node_attributes(G, "matched_with").get(atom_id)
    if matched_with is None:
        matched_with = set()
    matched_with.add(match)
    nx.set_node_attributes(G, {atom_id: matched_with}, name="matched_with")
    return G


def mark_aa_matches(
    G: nx.classes.graph.Graph, matches: list, aa_code: str
) -> nx.classes.graph.Graph:
    """
    For each residues that have been matched to amino_acid (aa_code)
        mark that information on each atom together with
        order of residue
    """

    for i in range(len(matches)):
        atom_ids = matches[i]
        match = (aa_code, i)
        for atom_id in atom_ids:
            G = append_match_to_atom(G, atom_id, match)
    return G


def mark_all_aa_matches(
    mol: rdkit.Chem.rdchem.Mol, G: nx.classes.graph.Graph, dict_aa: dict
) -> nx.classes.graph.Graph:
    """
    for each amino acid specie in dict:
        match its SMARTS pattern to molecule
        and map the matches on the atoms in molecule graph
    """

    for aa_code in dict_aa.keys():
        smarts = dict_aa[aa_code]["smarts"]
        pattern_mol = rdkit.Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern_mol, useChirality=True)
        G = mark_aa_matches(G, matches, aa_code)
    return G


def mark_backbone(
    mol: rdkit.Chem.rdchem.Mol, G: nx.classes.graph.Graph, double_bond_type
):
    backbone_atoms = longest_backbone(mol)
    n_ca_co = get_N_CA_CO(G, backbone_atoms, double_bond_type)
    G = get_labels(G, n_ca_co)
    return G


def map_residues(
    mol: rdkit.Chem.rdchem.Mol,
    G: nx.classes.graph.Graph,
    dict_aa: dict,
    double_bond_type,
):
    G = mark_backbone(mol, G, double_bond_type)
    # now we can teoretically make subgraphs N-CAXXX-CO and then fit
    # make copy and remove CO - N edges
    # and then connected graphs
    #
    #  https://stackoverflow.com/questions/21739569/finding-separate-graphs-within-a-graph-object-in-networkx

    G = mark_all_aa_matches(mol, G, dict_aa)
    return G


def get_residue_subgraphs(G: nx.classes.graph.Graph, backbone_nc_path: list) -> list:
    """
    input:
        molecule graph G
        list of backbone_nc_path atom ids ordered from N to C terminus

    a copy of graph is made G_copy

    edges between CO - N atoms are removed
    leaving connected subgraphs containing single residues

    Output:
        Subgraphs

    Subgraphs then can be hmm, analysed for matching the amino acids
    (maybe G_copy or the molecule created based on that)
    is fitted with different amino acid species
    then amino acids that are best match (cover most atoms)
    are chosen (n_atoms)
    and then each residue is named as a residue specie (Cys, Ala, Gly etc. )

    from that we have sequence ready

    """
    return


def choose_best_matches():
    """
    in this step the biggest aa specie
        is chosen or the one that matches the
        (greatest n_atoms)

    """
    return


def set_radical(mol: rdkit.Chem.rdchem.Mol, atom_id: int, radical_id: int):
    """

    Input:
        molecule
        atom_id to be set as attachment point
        id serving to identify attachment point

    selects atom by atom_id
    atom is set as a DummyAtom (meaning for us radical/attachment point)
    attachment point is then marked with number (radical_id)
    We will be identifying Attachment Point by its radical_id
    (for purpose of attaching molecules to more than one point for example)

    Output:
        molecule with marked radical

    Example:
        m  = rdkit.Chem.MolFromSmiles('CCO')
        assert rdkit.Chem.MolToSmiles(set_radical(m, 1, 2) ) == 'C[2*]O'

    """
    atom = mol.GetAtomWithIdx(atom_id)
    atom.SetAtomicNum(0)

    atom.SetIsotope(radical_id)
    return mol


def fit_res(residue: rdkit.Chem.rdchem.Mol, sidechains: dict) -> str:
    """
    Input:
        residue:
            rdkit Molecule with the residue as it was cut from
            rest of peptide chain
            e.q. NH2-CA(-Sidechain)-CO(O) at n term and in the middle
            and  NH2-CA(-Sidechain)-CO(O)OH at c term

        sidechains:
            a dictionary mapping chiral SMARTS codes
            for amino_acid sidechains to amino acids

    Process:
        for each sidechain
        molecular pattern from SMARTS is fitted to
        residue. the biggest match is then selected
        if no matches present Glycine is selected as matched amino acid
    """

    max_cover = 0
    best_fit = "G"

    for aa_code in sidechains.keys():
        smarts = sidechains[aa_code]
        pattern_mol = rdkit.Chem.MolFromSmarts(smarts)
        matches = residue.GetSubstructMatches(pattern_mol, useChirality=True)
        if matches:
            cover = len(matches[0])
            if cover > max_cover:
                max_cover = cover
                best_fit = aa_code
    return best_fit


def fit_seq(residues: list, sidechains: list) -> str:
    """

    Input:
        residues:
            list of rdkit.Chem.rdchem.Mol
        sidechains:
            dictionary of sidechain SMARTS mapped to amino acid species

    Output:
        sequence: amino acid sequence

    Process:
        each residue in list is fitted with SMARTS for best match
        best match symbols for each residues are joined into sequence str

    """
    seq_list = []
    for res in residues:
        aa = fit_res(res, sidechains)
        seq_list.append(aa)
    sequence = "".join(seq_list)
    return sequence


def smiles_to_seq(smiles: str, sidechains_smarts: dict) -> str:
    """
    Input:
        smiles: SMILES code
        sidechains: sidechains smarts dictionary

    Output:
        seq: amino acid sequence string

    Process:
        1. backbone atoms are found
        2. peptide bonds as 'broken'
        3. each amino acid is fitted with each SMARTS with dictionary
        4. best amino acid matches for each residue are joined into
           output sequence

    """

    mol = rdkit.Chem.MolFromSmiles(smiles)
    G = mol_to_nx(mol)
    backbone_atoms = longest_backbone(mol)
    n_ca_co_ids_path = get_N_CA_CO(G, backbone_atoms, rdkit.Chem.rdchem.BondType.DOUBLE)
    res_graphs = break_residues(G, n_ca_co_ids_path)
    residues = [nx_to_mol(g) for g in res_graphs]
    seq = fit_seq(residues, sidechains_smarts)
    return seq
