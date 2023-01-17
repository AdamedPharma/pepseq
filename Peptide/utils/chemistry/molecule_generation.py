import rdkit


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
    bbmol = rdkit.Chem.MolFromSmiles(bbsmiles)
    return bbmol


def MolFromSequenceThroughFASTA(sequence: str) -> rdkit.Chem.rdchem.Mol:
    """
    Input:
        sequence string

    Output:
        peptide molecule

    Process:

        FASTA string is generated from sequence string

        rdkit Molecule is generated from FASTA

        The advantage over generating Molecule directly
        from Sequence string is the lowercase letters in FASTA are
        interpreted as D-amino acids, which is not the case with generating
        Molecule from Sequence str

    """

    fasta = (
        """>
    %s
    """
        % sequence
    )
    mol = rdkit.Chem.MolFromFASTA(fasta, flavor=1)
    return mol
