import os
import unittest
from pathlib import Path

from Peptide.db_api.DataBase import FileSystemDbRepo
from Peptide.utils.chemistry.UseCase import read_seq_from_smiles

TEST_DATA_DB = os.path.join(
    Path(__file__).resolve().parent.parent, "Peptide/database/db.json"
)


class TestSmilesReader(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fsdb_repo = FileSystemDbRepo.read_from_json(path=TEST_DATA_DB)

        cls.our_semaglutide_smiles = "CC[C@H](C)[C@@H]1NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)C(NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@@H](NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CCC(=O)O)NC(=O)C(C)(C)NC(=O)[C@@H](N)Cc2c[nH]cn2)[C@@H](C)O)[C@@H](C)O)C(C)C)CSCC(=O)NCCCC(C(=O)Nc2ccc3oc(=O)cc(CC(=O)NCCOCCOCCC(=O)NCCCCC(NC(=O)CCCCCCCCCCCCCCCCC(=O)O)C(=O)O)c3c2)NC(=O)CSCC(C(=O)N[C@@H](Cc2c[nH]c3ccccc23)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](C(=O)N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)N[C@@H](CCCNC(=N)N)C(=O)NCC(=O)O)C(C)C)NC1=O"
        cls.our_semaglutide_seq_no_x = (
            "H{aMeAla}EGTFTSDVSSYLEG{modC}AAKEFI{modC}WLVRGRG"
        )
        cls.our_semaglutide_seq_x = "HXEGTFTSDVSSYLEGXAAKEFIXWLVRGRG"
        return

    def test_read_seq_from_smiles_mod_as_x(self):

        seq, mods_list_json = read_seq_from_smiles(
            smiles=self.our_semaglutide_smiles, repo=self.fsdb_repo, mod_as_x=True
        )

        assert seq == self.our_semaglutide_seq_x
        assert mods_list_json == [
            {
                "modification_smiles": "CC(=O)NCCCC(NC(C)=O)C(=O)Nc1ccc2oc(=O)cc(CC(=O)NCCOCCOCCC(=O)NCCCCC(NC(=O)CCCCCCCCCCCCCCCCC(=O)O)C(=O)O)c2c1",
                "connecting_residues": [(17, "C"), (24, "C")],
            }
        ]
        return

    def test_read_seq_from_smiles(self):

        seq, mods_list_json = read_seq_from_smiles(
            smiles=self.our_semaglutide_smiles, repo=self.fsdb_repo, mod_as_x=False
        )

        assert seq == "H{aMeAla}EGTFTSDVSSYLEG{modC}AAKEFI{modC}WLVRGRG"
        assert mods_list_json == [
            {
                "modification_smiles": "CC(=O)NCCCC(NC(C)=O)C(=O)Nc1ccc2oc(=O)cc(CC(=O)NCCOCCOCCC(=O)NCCCCC(NC(=O)CCCCCCCCCCCCCCCCC(=O)O)C(=O)O)c2c1",
                "connecting_residues": [(17, "C"), (24, "C")],
            }
        ]
