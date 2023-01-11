import os
import unittest
from pathlib import Path

from tests.data import (
    dermcidin_stapled_with_anchor,
    dermcidin_stapled_with_anchor_seq,
    seq_wo_mods,
)

from Peptide.db_api.DataBase import FileSystemDbRepo
from Peptide.models.AminoAcidInstance import AminoAcidInstance
from Peptide.utils.chemistry.smi2seq_utils import (
    get_atm_id_to_potential_res_ids_names_dict,
    get_rows,
)
from Peptide.utils.chemistry.Smi2SeqObj import get_atm_id_to_res
from Peptide.utils.chemistry.SubstructureGraph import get_ordered_nodes
from Peptide.utils.chemistry.UseCase import (
    UseCase,
    find_group_connected,
    get_modified_residues,
    get_rows,
    get_seq,
    read_seq_from_smiles,
)

TEST_DATA_DB = os.path.join(
    Path(__file__).resolve().parent.parent, "Peptide/database/db.json"
)


def check_peptide(peptide):
    assert (
        peptide._smiles
        == "CC(C)C[C@H](NC(=O)[C@@H](NC(=O)[C@H](CO)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCCCN)NC(=O)[C@@H](NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](Cc1c[nH]cn1)NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)CNC(=O)[C@H](CCCCN)NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](CO)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CS)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCCN)NC(=O)CNC(=O)[C@H](CC(C)C)NC(=O)CNC(=O)[C@H](CS)NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CCCCN)NC(=O)C(C)(C)NC(=O)CNC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)CNC(=O)[C@H](CCCCN)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CO)NC(=O)[C@@H](N)CO)C(C)C)C(C)C)C(C)C)C(C)C)C(C)C)C(C)C)C(C)C)C(=O)O"
    )
    assert peptide._parameters == {
        "mw": 4906.71,
        "exact_mw": 4903.6,
        "logp": -20.37,
        "hba": 72,
        "hbd": 71,
        "num_heavy_atoms": 342,
        "heteroatoms": 131,
        "rotatable_bonds": 176,
        "rings": 1,
        "formula": "C213H363*2N57O70S2",
        "atoms": 344,
        "tpsa": 2024.24,
    }

    assert peptide._n_term.name == "H"
    assert peptide._n_term.smiles == "[*][H] |$_R2;_CO$|"
    assert peptide._c_term.name == "OH"
    assert peptide._c_term.smiles == "[OH](*) |$_OH;_R$|"

    amino_acids = peptide.amino_acids

    assert amino_acids[10].symbol == "aMeAla"
    assert amino_acids[10].smiles == "CC(C)(N*)C(*)=O |$_CB;_CA;_Me;_N;_R1;_CO;_R2$|"

    assert amino_acids[16 - 1].symbol == "C"
    assert (
        amino_acids[16 - 1].smiles
        == "N(*)[C@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|"
    )

    return


class TestUseCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.smiles = dermcidin_stapled_with_anchor
        cls.seq = seq_wo_mods
        met_AASTAA_smiles = "C[C@H](N)C(=O)N[C@@H](C)C(=O)N[C@@H](CO)C(=O)N[C@@](C)([C@@H](C)O)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(O)=O"
        AfASmTAA_smiles = "C[C@H](N)C(=O)N[C@@H](CF)C(=O)N[C@@H](CO)C(=O)N[C@@](C)([C@@H](C)O)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(O)=O"
        me_smi = "C[C@H](N)C(=O)N[C@@H](C)C(=O)N[C@@](C)(CO)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(O)=O"

        cls.fsdb_repo = FileSystemDbRepo.read_from_json(path=TEST_DATA_DB)

        cls.uc = UseCase(smiles=cls.smiles, repo=cls.fsdb_repo)
        cls.smi_seq_tup = (
            (cls.smiles, cls.seq),
            (met_AASTAA_smiles, "AAS{aMeX}AA"),
            (AfASmTAA_smiles, "AXS{aMeX}AA"),
            (me_smi, "AA{aMeSer}AA"),
        )
        return

    def test_read_seq_from_smiles(self):
        seq, mods_list_json = read_seq_from_smiles(
            smiles=self.smiles, repo=self.fsdb_repo, mod_as_x=False
        )
        assert seq == dermcidin_stapled_with_anchor_seq
        assert mods_list_json == [
            {
                "modification_smiles": "CC(=O)NCCC[C@H](NC(C)=O)C(=O)NCCOCCO",
                "connecting_residues": [(16, "C"), (22, "C")],
            }
        ]
        return

    def test_get_seq(self):
        symbol_smiles_list = [
            ("H", "N(*)[C@H](C(*)=O)Cc1c[nH]cn1 |$_N;_R1;_CA;_CO;_R2;;_CB;;;$|"),
            ("aMeAla", "CC(C)(N*)C(*)=O |$_CB;_CA;_Me;_N;_R1;_CO;_R2$|"),
            ("E", "N(*)[C@H](C(*)=O)CCC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;;$|"),
            ("G", "N(*)CC(*)=O |$_N;_R1;_CA;_CO;_R2$|"),
            ("K", "N(*)[C@H](C(*)=O)CCC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;;$|"),
            ("G", "N(*)CC(*)=O |$_N;_R1;_CA;_CO;_R2$|"),
            ("C", "N(*)[C@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|"),
            ("A", "C[C@H](N*)C(*)=O |$_CB;_CA;_N;_R1;_CO;_R2;$|"),
        ]

        amino_acids = []

        for symbol, smi in symbol_smiles_list:
            aa = AminoAcidInstance()
            aa.smiles = smi
            aa.symbol = symbol

            amino_acids.append(aa)

        mod_res_nums = [(5, "K"), (7, "C")]

        seq_with_x = get_seq(amino_acids, mod_res_nums, mod_as_x=True)
        assert seq_with_x == "HXEGXGXA"

        seq_without_x = get_seq(amino_acids, mod_res_nums, mod_as_x=False)
        assert seq_without_x == "H{aMeAla}EG{modK}G{modC}A"

    def test_deduce_seq(self):
        deduced_peptide, deduced_seq, atm_id_to_res_dict = self.uc.deduce_seq()
        assert deduced_seq == self.seq

        check_peptide(deduced_peptide)

        keys = set(atm_id_to_res_dict.keys())
        mol_atoms = [
            (i, i.GetIdx(), i.GetSymbol()) for i in self.uc.molecule.GetAtoms()
        ]
        extra_atoms = [i for i in mol_atoms if (i[1] not in keys)]
        assert len(extra_atoms) == 21

        connected_atom_sets = find_group_connected([i[0] for i in extra_atoms])
        assert len(connected_atom_sets) == 2

        mod_res_nums = get_modified_residues(connected_atom_sets, atm_id_to_res_dict)

        assert mod_res_nums == [(22, "C"), (16, "C")]

        aa_species_assigned_potential_atom_ids = self.uc.smi2seq_obj.aa_matches(
            self.uc.repo.aa_smiles_dict, code="SMARTS"
        )
        assert len(aa_species_assigned_potential_atom_ids["D"]) == 6
        assert len(aa_species_assigned_potential_atom_ids["A"]) == 3
        ca_atom_id_to_res_id_dict = self.uc.get_ca_atom_id_to_res_id_dict(
            self.uc.smi2seq_obj
        )
        assert len(ca_atom_id_to_res_id_dict) == 48
        aa_species_assigned_potential_atom_ids = self.uc.smi2seq_obj.aa_matches(
            self.uc.repo.aa_smiles_dict, code="SMARTS"
        )
        no_A_residues = len(aa_species_assigned_potential_atom_ids["A"])
        assert no_A_residues == 3
        first_A_res_atoms = aa_species_assigned_potential_atom_ids["A"][0]
        assert len(first_A_res_atoms) == 6

        rows = get_rows(
            aa_species_assigned_potential_atom_ids, ca_atom_id_to_res_id_dict
        )

        atm_id_to_potential_res_ids_names_dict = get_atm_id_to_potential_res_ids_names_dict(
            rows
        )

        aas_ordered_by_size = get_ordered_nodes(self.uc.repo.aa_smiles_dict)

        atm_id_to_res_dict = get_atm_id_to_res(
            atm_id_to_potential_res_ids_names_dict, aas_ordered_by_size
        )
        symbols_list = [i[1] for i in sorted(list(set(atm_id_to_res_dict.values())))]
        for i in range(len(symbols_list)):
            if len(symbols_list[i]) > 1:
                symbols_list[i] = "{%s}" % symbols_list[i]
        deduced_seq = "".join(symbols_list)

        assert deduced_seq == self.seq
        return

    def test_values(self):
        for smi, seq in self.smi_seq_tup:
            uc = UseCase(smiles=smi, repo=self.fsdb_repo)
            deduced_peptide, deduced_seq, atm_id_to_res_dict = uc.deduce_seq()
            assert deduced_seq == seq
        return
