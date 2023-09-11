import networkx as nx
import rdkit

from pepseq.Peptide.utils.chemistry.mol_to_nx_translation import \
    mol_to_nx


class Functionality(object):
    pass


def generate_bb_smiles(num_res: int, OH: bool = False) -> str:
    """
    generates SMILES for polyGlycine peptide
    length (number of Gly residues) is given by num_res
        parameter
    """
    if OH:
        bb_smiles = "O" + "C(=O)CN" * num_res  # generate backbone SMILES
    else:
        bb_smiles = "NCC(=O)" * num_res  # generate backbone SMILES

    return bb_smiles


def generate_bb_mol(num_res: int, OH: bool = False) -> rdkit.Chem.rdchem.Mol:
    """
    generates rdkit.Chem.rdchem.Mol molecule for polyGlycine peptide
    length (number of Gly residues) is given by num_res
        parameter

    """
    bbsmiles = generate_bb_smiles(num_res, OH=OH)
    bbmol = rdkit.Chem.MolFromSmiles(bbsmiles)
    return bbmol


class GetLongestPolymerWithin(Functionality):
    """
    The Longest Polymer of (monomer) form -[NCC(=O) ]*n- is being
      rdkit.substructureMatch
    of rdkit.MoleculeFromSMILES ProteinBackbone Is (SMARTS * N)?
    """

    def __init__(self):
        return

    def execute(self, peptide_molecule: rdkit.Chem.rdchem.Mol) -> tuple:
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

    def execute(self, peptide_molecule: rdkit.Chem.rdchem.Mol) -> list:
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
