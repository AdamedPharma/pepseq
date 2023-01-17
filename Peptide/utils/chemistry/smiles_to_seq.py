import networkx as nx
import rdkit
from Peptide.utils.chemistry.molecule_generation import generate_bb_mol
from Peptide.utils.chemistry.molecule_graph_generation import mol_to_nx, nx_to_mol


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


def fit_seq(residues: list, sidechains: dict) -> str:
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


def get_peptide_bonds(n_ca_co: list) -> list:
    """

    Input:
        n_ca_co - list of [N, Ca, Co ... N, Ca, Co]
        atoms

    Output:
        peptide_bonds [(1Co, 2N), (2Co, 3N), ...]

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

    TO DO:
        TEST whether it works on disulfides
            stapled molecules and such

        Probably not

        Possible Solution:
            Each Residue must contain only one backbone (N-CaCo)
            So for disulfides and staples we must choose the best
            fit for each of the NCaCO group


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


def get_N_CA_CO(
    G: nx.classes.graph.Graph, backbone_atoms: list, double_bond_type
) -> list:
    """
    Input:
         backbone atoms: List of atom ids
         double_bond_type: rdkit.Chem.rdchem.BondType.DOUBLE
         was put as variable to avoid rdkit dependency
         which it now fails to do so i might reconsider

    Output:
        an ordered list of backbone N and C atoms
        Nterm-Ca-Co- ... -N-Ca-Coterm atoms

    Process:


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


def break_molecule_into_residues(mol: rdkit.Chem.rdchem.Mol) -> list:
    """

    Process:
        1. backbone atoms are found
        2. peptide bonds are 'broken'

    """
    G = mol_to_nx(mol)
    backbone_atoms = longest_backbone(mol)
    n_ca_co_ids_path = get_N_CA_CO(G, backbone_atoms, rdkit.Chem.rdchem.BondType.DOUBLE)
    res_graphs = break_residues(G, n_ca_co_ids_path)
    residues = [nx_to_mol(g) for g in res_graphs]
    return residues


def smiles_to_seq(smiles: str, sidechains_smarts: dict) -> str:
    """
    Input:
        smiles: SMILES code
        sidechains: sidechains smarts dictionary

    Output:
        seq: amino acid sequence string

    Process:
        1. backbone atoms are found
        2. peptide bonds are 'broken'
        3. each amino acid is fitted with each SMARTS with dictionary
        4. best amino acid matches for each residue are joined into
           output sequence

    """

    mol = rdkit.Chem.MolFromSmiles(smiles)
    residues = break_molecule_into_residues(mol)
    seq = fit_seq(residues, sidechains_smarts)
    return seq
