import unittest

import rdkit
from Peptide.utils.chemistry.Graph import connect_radicals, smiles_to_seq


class TestSmilesToSequence(unittest.TestCase):
    def test1(self):
        seq = "CACAC"
        sidechains = {"A": "[CH3X4]", "C": "[CH2X4][SX2H,SX1H0-]"}

        seq_mol = rdkit.Chem.MolFromSequence(seq)
        smiles = rdkit.Chem.MolToSmiles(seq_mol)

        assert smiles_to_seq(smiles, sidechains) == "CACAC"

    def test_connect_radicals(self):
        smi_1 = "[1*]NCN"
        smi_2 = "[1*]SC"

        mol1 = rdkit.Chem.MolFromSmiles(smi_1)
        mol2 = rdkit.Chem.MolFromSmiles(smi_2)

        mol_union = connect_radicals(mol1, mol2)
        smi_union = rdkit.Chem.MolToSmiles(mol_union)

        assert smi_union == "CSNCN"
