import unittest

import rdkit
import rdkit.Chem
from Peptide.utils.chemistry.Graph import (
    connect_on_radical,
    reconstruct_molecule_from_json,
)
from Peptide.utils.chemistry.smiles_to_seq import smiles_to_seq


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

        mol_union = connect_on_radical(mol1, mol2)
        smi_union = rdkit.Chem.MolToSmiles(mol_union)

        assert smi_union == "CSNCN"

    def test_reconstruct_molecule_from_json(self):
        mol_json = {
            "attachment_points": {1: (2, "C", "S"), 2: (4, "C", "S")},
            "sequence": "GCGCG",
            "modification": "[1*]SCNOCN[2*]",
        }

        mol_r = reconstruct_molecule_from_json(mol_json)
        assert (
            rdkit.Chem.MolToSmiles(mol_r)
            == "NCC(=O)N[C@H]1CSSCNOCNSC[C@@H](C(=O)NCC(=O)O)NC(=O)CNC1=O"
        )
