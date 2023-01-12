import os
import unittest
from pathlib import Path

import rdkit.Chem

from Peptide.db_api.DataBase import FileSystemDbRepo
from Peptide.utils.chemistry.Smi2SeqObj import (
    Smi2Seq,
    choose_biggest,
    generate_bb_mol,
    generate_bb_smiles,
    get_atm_id_to_res,
)
from tests.data import dermcidin_stapled_with_anchor

TEST_DATA_DB = os.path.join(
    Path(__file__).resolve().parent.parent, "Peptide/database/db.json"
)


class TestSmi2SeqObj(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.smiles = dermcidin_stapled_with_anchor
        cls.molecule = rdkit.Chem.MolFromSmiles(cls.smiles)
        cls.smi2seq_obj = Smi2Seq(molecule=cls.molecule)
        cls.fsdb_repo = FileSystemDbRepo.read_from_json(path=TEST_DATA_DB)
        return

    def test_gly_matches(self):
        gly_matches = self.smi2seq_obj.gly_matches
        assert len(gly_matches) == 49
        return

    def test_CAatoms(self):
        assert len(self.smi2seq_obj.CAatoms) == 48
        return

    def test_get_CAatoms(self):
        assert len(self.smi2seq_obj.get_CAatoms(longest_bb=True)) == 48
        return

    def test_label_CAatoms(self):
        self.smi2seq_obj.label_CAatoms(longest_bb=True)
        return

    def test_backbone(self):
        assert len(self.smi2seq_obj.backbone) == 192
        return

    def test_longest_backbone(self):
        backbone_atoms_ids = list(self.smi2seq_obj.longest_backbone)
        assert len(backbone_atoms_ids) == 192
        return

    def test_find_matches(self):
        aa = self.fsdb_repo.aa_smiles_dict.get("A")
        matches = self.smi2seq_obj.find_matches(aa_mol_code=aa.smarts, code="SMARTS")
        assert len(matches) == 3
        return

    def test_bb_list(self):
        self.smi2seq_obj.bb_list
        return

    def test_longest_bb_list(self):
        assert len(self.smi2seq_obj.longest_bb_list) == 192
        return

    def test_renumber(self):
        self.smi2seq_obj.renumber()
        return

    def test_aa_matches(self):
        aa_species_assigned_potential_atom_ids = self.smi2seq_obj.aa_matches(
            self.fsdb_repo.aa_smiles_dict, code="SMARTS"
        )
        assert len(aa_species_assigned_potential_atom_ids["D"]) == 6
        assert len(aa_species_assigned_potential_atom_ids["A"]) == 3
        return

    def test_generate_bb_smiles(self):
        assert generate_bb_smiles(3, OH=False) == "C(=O)CNC(=O)CNC(=O)CN"
        assert generate_bb_smiles(5, OH=False) == "C(=O)CNC(=O)CNC(=O)CNC(=O)CNC(=O)CN"
        return

    def test_generate_bb_mol(self):
        assert (
            rdkit.Chem.MolToSmiles(generate_bb_mol(3, OH=False))
            == "NCC(=O)NCC(=O)NCC=O"
        )
        assert (
            rdkit.Chem.MolToSmiles(generate_bb_mol(5, OH=False))
            == "NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC=O"
        )
        return

    def test_choose_biggest(self):
        aas_ordered_by_size = ["G", "A", "V", "I"]
        list_of_potential_amino_acids = ["A", "G", "I", "V"]
        assert choose_biggest(list_of_potential_amino_acids, aas_ordered_by_size) == "I"
        return

    def test_get_atm_id_to_res(self):
        """
        aa_species_assigned_potential_atom_ids = self.smi2seq_obj.aa_matches(
            self.fsdb_repo.aa_smiles_dict, code="SMARTS"
        )
        ca_atom_id_to_res_id_dict = self.uc.get_ca_atom_id_to_res_id_dict(self.smi2seq_obj)

        rows = get_rows(
            aa_species_assigned_potential_atom_ids, ca_atom_id_to_res_id_dict
        )

        atm_id_to_potential_res_ids_names_dict = (
            get_atm_id_to_potential_res_ids_names_dict(rows)
        )

        aas_ordered_by_size = get_ordered_nodes(self.repo.aa_smiles_dict)
        """
        aas_ordered_by_size = [
            "G",
            "a",
            "A",
            "s",
            "v",
            "f",
            "S",
            "V",
            "F",
            "y",
            "w",
            "t",
            "r",
            "q",
            "P",
            "n",
            "m",
            "l",
            "k",
            "i",
            "h",
            "e",
            "d",
            "c_nonmod",
            "c",
            "Y",
            "W",
            "T",
            "R",
            "Q",
            "p",
            "N",
            "M",
            "L",
            "K",
            "I",
            "H",
            "E",
            "D",
            "C_nonmod",
            "C",
            "aMeAla",
        ]

        atm_id_to_potential_res_ids_names_dict = {
            0: [(1, "H")],
            1: [(1, "H")],
            2: [(1, "H")],
            3: [(1, "H")],
            4: [(1, "H"), (2, "aMeAla")],
            5: [(2, "aMeAla")],
            6: [(2, "aMeAla")],
            7: [(2, "aMeAla")],
            8: [(3, "E"), (2, "aMeAla")],
            9: [(3, "E")],
            10: [(3, "E")],
            11: [(3, "E")],
            12: [(3, "E"), (4, "G")],
            13: [(4, "G")],
            14: [(4, "G")],
            15: [(4, "G")],
            16: [(4, "G"), (5, "T")],
            17: [(5, "T")],
            18: [(5, "T")],
            19: [(5, "T")],
            20: [(6, "F"), (5, "T")],
            21: [(6, "F")],
            22: [(6, "F")],
            23: [(6, "F")],
            24: [(6, "F"), (7, "T")],
            25: [(7, "T")],
            26: [(7, "T")],
            27: [(7, "T")],
            28: [(8, "S"), (7, "T")],
            29: [(8, "S")],
            30: [(8, "S")],
            31: [(8, "S")],
            32: [(9, "D"), (8, "S")],
            33: [(9, "D")],
            34: [(9, "D")],
            35: [(9, "D")],
            36: [(9, "D"), (10, "V")],
            37: [(10, "V")],
            38: [(10, "V")],
            39: [(10, "V")],
            40: [(11, "S"), (10, "V")],
            41: [(11, "S")],
            42: [(11, "S")],
            43: [(11, "S")],
            44: [(11, "S"), (12, "S")],
            45: [(12, "S")],
            46: [(12, "S")],
            47: [(12, "S")],
            48: [(12, "S"), (13, "Y")],
            49: [(13, "Y")],
            50: [(13, "Y")],
        }
        atm_id_to_res_dict = get_atm_id_to_res(
            atm_id_to_potential_res_ids_names_dict, aas_ordered_by_size
        )
        assert atm_id_to_res_dict == {
            0: (1, "H"),
            1: (1, "H"),
            2: (1, "H"),
            3: (1, "H"),
            4: (2, "aMeAla"),
            5: (2, "aMeAla"),
            6: (2, "aMeAla"),
            7: (2, "aMeAla"),
            8: (2, "aMeAla"),
            9: (3, "E"),
            10: (3, "E"),
            11: (3, "E"),
            12: (3, "E"),
            13: (4, "G"),
            14: (4, "G"),
            15: (4, "G"),
            16: (5, "T"),
            17: (5, "T"),
            18: (5, "T"),
            19: (5, "T"),
            20: (5, "T"),
            21: (6, "F"),
            22: (6, "F"),
            23: (6, "F"),
            24: (7, "T"),
            25: (7, "T"),
            26: (7, "T"),
            27: (7, "T"),
            28: (7, "T"),
            29: (8, "S"),
            30: (8, "S"),
            31: (8, "S"),
            32: (9, "D"),
            33: (9, "D"),
            34: (9, "D"),
            35: (9, "D"),
            36: (9, "D"),
            37: (10, "V"),
            38: (10, "V"),
            39: (10, "V"),
            40: (10, "V"),
            41: (11, "S"),
            42: (11, "S"),
            43: (11, "S"),
            44: (11, "S"),
            45: (12, "S"),
            46: (12, "S"),
            47: (12, "S"),
            48: (13, "Y"),
            49: (13, "Y"),
            50: (13, "Y"),
        }
        return
