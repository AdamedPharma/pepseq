import unittest

import rdkit.Chem
from reconstructing_molecules import UseCaseCombineCombo


class TestSmi2SeqObj(unittest.TestCase):
    def test1(self):
        smi = "[1*]C=CN=CC(C)=O"
        mol = rdkit.Chem.MolFromSmiles(smi)
        combo = rdkit.Chem.rdmolops.CombineMols(mol, mol)
        uc = UseCaseCombineCombo(combo)
        mol_out = uc.execute()
        assert rdkit.Chem.MolToSmiles(mol_out) == "CC(=O)C=NC=C~C=CN=CC(C)=O"

    def test2(self):
        smi = "[1*]C=CN=CC(C)O[2*]"
        mol = rdkit.Chem.MolFromSmiles(smi)
        combo = rdkit.Chem.rdmolops.CombineMols(mol, mol)
        uc = UseCaseCombineCombo(combo)
        mol_out = uc.execute()
        assert rdkit.Chem.MolToSmiles(mol_out) == "CC1C=NC=C~C=CN=CC(C)O~O1"
