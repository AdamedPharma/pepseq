import rdkit
import rdkit.Chem


def mols_are_identical(
    mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol
) -> bool:
    """
    Check if two RDKit molecule objects are identical.
    This function compares two RDKit molecule objects to determine if they are
    identical by checking if each molecule is a substructure of the other,
    considering chirality.
    Args:
        mol1 (rdkit.Chem.rdchem.Mol): The first molecule to compare.
        mol2 (rdkit.Chem.rdchem.Mol): The second molecule to compare.
    Returns:
        bool: True if the molecules are identical, False otherwise.
    """
    are_identical = mol1.HasSubstructMatch(
        mol2, useChirality=True
    ) and mol2.HasSubstructMatch(mol1, useChirality=True)
    return are_identical


def smiles_are_identical(smi1: str, smi2: str) -> bool:
    """
    Check if two SMILES strings represent the same molecule.
    Args:
        smi1 (str): The first SMILES string.
        smi2 (str): The second SMILES string.
    Returns:
        bool: True if the molecules are identical, False otherwise.
    """
    mol1 = rdkit.Chem.MolFromSmiles(smi1)
    mol2 = rdkit.Chem.MolFromSmiles(smi2)
    return mols_are_identical(mol1, mol2)
