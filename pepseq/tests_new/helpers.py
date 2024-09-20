import rdkit
import rdkit.Chem


def mols_are_identical(mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol) -> bool:
    are_identical = mol1.HasSubstructMatch(
        mol2, useChirality=True
    ) and mol2.HasSubstructMatch(mol1, useChirality=True)
    return are_identical


def smiles_are_identical(smi1: str, smi2: str) -> bool:
    mol1 = rdkit.Chem.MolFromSmiles(smi1)
    mol2 = rdkit.Chem.MolFromSmiles(smi2)
    return mols_are_identical(mol1, mol2)
