import networkx as nx
import rdkit

from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import \
    mol_to_nx


class Functionality(object):
    """
    This class represents the functionality of the program.
    """
    pass


def generate_bb_smiles(num_res: int, OH: bool = False) -> str:
    """

    Generate SMILES code for peptide backbone (SMILES for polyGlycine peptide)

    :param num_res: number of amino_acid residues
    :type num_res: int

    :param OH: whether to generate hydroxyl group at the end
    :type OH: bool
    
    :return: bb_smiles - SMILES code for backbone
    :rtype: str

    """
    if OH:
        bb_smiles = "O" + "C(=O)CN" * num_res  # generate backbone SMILES
    else:
        bb_smiles = "NCC(=O)" * num_res  # generate backbone SMILES

    return bb_smiles


def generate_bb_mol(num_res: int, OH: bool = False) -> rdkit.Chem.rdchem.Mol:
    """

    Generate rdkit.Chem.rdchem.Mol molecule for peptide backbone (polyGlycine peptide)

    :param num_res: number of amino_acid residues
    :type num_res: int

    :param OH: whether to generate hydroxyl group at the end
    :type OH: bool
    
    :return: bb_smiles - rdkit.Chem.rdchem.Mol molecule object for backbone
    :rtype: rdkit.Chem.rdchem.Mol

    """
    bb_smiles = generate_bb_smiles(num_res, OH=OH)
    bb_mol = rdkit.Chem.MolFromSmiles(bb_smiles)
    return bb_mol


class GetLongestPolymerWithin(Functionality):
    """
    The Longest Polymer of (monomer) form -[NCC(=O) ]*n- is being
      rdkit.substructureMatch
    of rdkit.MoleculeFromSMILES ProteinBackbone Is (SMARTS * N)?
    """

    def __init__(self):
        return

    def execute(self, peptide_molecule: rdkit.Chem.rdchem.Mol) -> tuple:
        """

        Find longest peptide backbone that fits into Peptide Molecule

        :param peptide_molecule: Modified Peptide molecule created from SMILES
            and to be decomposed into amino acid chain with known sequence and its
            covalent modifications

        :type peptide_molecule: rdkit.Chem.rdchem.Mol

        :return: bb - tuple of atom indices that has been matched with backbone
            substructure (e.g.  (1,3,4, ...))
        :rtype: tuple
        
        """

        matches = [None]
        num_res = 0
        while matches:
            bb = matches[0]
            num_res += 1
            bbmol = generate_bb_mol(num_res, OH=False)
            matches = peptide_molecule.GetSubstructMatches(bbmol)
        return bb


class MarkingPeptideBackbone(Functionality):
    """
    Input:

    rdkit.Chem.rdchem.Mol object representing modified peptide molecule

    Output:

    rdkit.Chem.rdchem.Mol object representing modified peptide molecule
    with:
        - protein backbone atoms N-Ca-Co assigned new properties:
          ResidueID and AtomName (N, CA, or CO)
        
        - protein backbone peptide bonds connecting residues
          assigned new property ("is_peptide_bond" = "True")


    Purpose:

    Identifying Peptide Bonds and Protein Backbone Atoms allows as to
    identify the sequence of amino acids in peptide. In the next steps
    it allows us to separate peptide atoms from external modifications
    like staples (e.g. ornithine) and identify internal connections
    between non neighbouring residue like disulfide bridges.

    """

    def __init__(self):
        return

    def execute(self, peptide_molecule: rdkit.Chem.rdchem.Mol) -> rdkit.Chem.rdchem.Mol:
        """

        Label peptide bonds within Peptide Molecule.
        Label backbone atoms ["N", "CA", "CO", "O"] with their respective names.
        Label backbone atoms ["N", "CA", "CO", "O"] with their respective ResidueIDs
        consecutively.

        :param peptide_molecule: Modified Peptide molecule created from SMILES
            and to be decomposed into amino acid chain with known sequence and its
            covalent modifications

        :type peptide_molecule: rdkit.Chem.rdchem.Mol

        :return: peptide_molecule - rdkit.Chem.rdchem.Mol molecule object with
            peptide bonds labeled as such (added new property: "is_peptide_bond" set to "True")
          tuple of atom indices that has been matched with backbone
            substructure (e.g.  (1,3,4, ...))
        :rtype: rdkit.Chem.rdchem.Mol
        
        """

        bb_ats = ["N", "CA", "CO", "O"]

        match = GetLongestPolymerWithin().execute(peptide_molecule)
        for match_atom_id in range(len(match)):
            pep_mol_atom_id = match[match_atom_id]
            ind = int((match_atom_id % 4))
            AtomName = bb_ats[ind]
            ResID = int(match_atom_id / 4) + 1
            pep_mol_atom = peptide_molecule.GetAtomWithIdx(pep_mol_atom_id)
            pep_mol_atom.SetProp("AtomName", AtomName)
            pep_mol_atom.SetProp("ResID", str(ResID))

            if AtomName == "N" and ResID >= 2:
                prev_CO_atom_id = match[match_atom_id - 2]
                N_atom = pep_mol_atom
                N_bonds = N_atom.GetBonds()
                prev_CO_atom = peptide_molecule.GetAtomWithIdx(prev_CO_atom_id)
                prev_CO_bonds = prev_CO_atom.GetBonds()
                N_bonds_ids = set([i.GetIdx() for i in N_bonds])
                prev_CO_bond_ids = set([i.GetIdx() for i in prev_CO_bonds])
                peptide_bond_id = (N_bonds_ids & prev_CO_bond_ids).pop()
                peptide_bond = peptide_molecule.GetBondWithIdx(peptide_bond_id)
                peptide_bond.SetProp("is_peptide_bond", "True")
        return peptide_molecule


class BreakingIntoResidueCandidateSubgraphs(Functionality):
    """

    We use nx.classes.graph.Graph representation of modified peptide molecules

    Edges representing peptide bonds are deleted separating molecule object
        into several subgraphs (nx.connected_components).
    
    Unmodified Residues become separate Graphs.

    Modified Residues form same Graph with modification atoms and
    need to be separated in next steps.

    Set of Non neighboring residues connected covalently (e.g. through disulfide bonds)
    form single graph and need to be separated in next steps.

    """

    def __init__(self):
        return

    def execute(self, peptide_molecule: rdkit.Chem.rdchem.Mol) -> list[nx.classes.graph.Graph]:
        """

        Label peptide bonds within Peptide Molecule

        :param peptide_molecule: Modified Peptide molecule created from SMILES
        with peptide bonds and backbone atoms labeled with [N, CA, C or  O ] atom name
        and residue index in peptide amino acid chain.

        :type peptide_molecule: rdkit.Chem.rdchem.Mol

        :return: residue_graphs - list of fragments generated through cleavage of peptide bonds
        fragments are molecular graphs generated through translation of rdkit.Chem.rdchem.Mol
        to nx.classes.graph.Graph. Fragment can be:
        a) single residue
        b) single residue with external modification (such as palmitoylation)
        c) two or more residues bound by internal modification (such as disulfide bond, or cyclization)
        d) two or more residues bound by external modification (such as molecular staple)


          peptide_molecule - rdkit.Chem.rdchem.Mol molecule object 
        :rtype: list[nx.classes.graph.Graph]
        
        """

        peptide_molecule = MarkingPeptideBackbone().execute(peptide_molecule)
        G = mol_to_nx(peptide_molecule)
        peptide_bonds = [
            (u, v)
            for u, v, e in G.edges(data=True)
            if e.get("is_peptide_bond") == "True"
        ]
        G_copy = G.copy()
        G_copy.remove_edges_from(peptide_bonds)
        g = (G_copy.subgraph(c) for c in nx.connected_components(G_copy))
        residue_graphs = list(g)
        return residue_graphs
