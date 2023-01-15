import unittest

import rdkit
from Peptide.utils.chemistry.Graph import smiles_to_seq


class TestSmilesToSequence(unittest.TestCase):
    def test1(self):
        seq = "CACAC"
        sidechains = {"A": "[CH3X4]", "C": "[CH2X4][SX2H,SX1H0-]"}

        seq_mol = rdkit.Chem.MolFromSequence(seq)
        smiles = rdkit.Chem.MolToSmiles(seq_mol)

        assert smiles_to_seq(smiles, sidechains) == "CACAC"
