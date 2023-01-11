import unittest

from Peptide.db_api.DataBase import (
    AminoAcid,
    FileSystemDbRepo,
    join_dicts,
    lower_dict,
    read_aa_smiles_dict,
    read_AminoAcid_from_dict,
)


class TestDataBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.in_dict = {
            "smiles": "CC(C)(N)C=O",
            "smiles_radical": "CC(C)(N*)C(*)=O |$_CB;_CA;_Me;_N;_R1;_CO;_R2$|",
            "three_letter_code": "Aib",
            "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4]([CH3X4])([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
            "smarts_comments": "w normalnej sekwencji powinno sie outputowac jako {aMeAla}; a w zaiksowanej",
        }

        cls.aa_smiles_dict = {
            "aMeAla": {
                "smiles": "CC(C)(N)C=O",
                "smiles_radical": "CC(C)(N*)C(*)=O |$_CB;_CA;_Me;_N;_R1;_CO;_R2$|",
                "three_letter_code": "Aib",
                "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4]([CH3X4])([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
                "smarts_comments": "w normalnej sekwencji powinno sie outputowac jako {aMeAla}; a w zaiksowanej",
            },
            "A": {
                "smiles": "C[C@H](N)C=O",
                "smiles_radical": "C[C@H](N*)C(*)=O |$_CB;_CA;_N;_R1;_CO;_R2;$|",
                "three_letter_code": "Ala",
                "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
                "smarts_comments": "",
            },
            "C": {
                "smiles": "N[C@H](C=O)CS",
                "smiles_radical": "N(*)[C@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
                "three_letter_code": "Cys",
                "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CHX4]([CH2X4][SX2H,SX1H0-,SX2])[CX3](=[OX1])[OX2H,OX1-,N]",
                "smarts_comments": "Hits acid and conjugate base",
            },
        }

        return

    def test_read_AminoAcid_from_dict(self):
        amino_acid_namedtuple = read_AminoAcid_from_dict(self.in_dict)
        assert amino_acid_namedtuple.smiles == "CC(C)(N)C=O"
        assert (
            amino_acid_namedtuple.smiles_radical
            == "CC(C)(N*)C(*)=O |$_CB;_CA;_Me;_N;_R1;_CO;_R2$|"
        )
        assert amino_acid_namedtuple.three_letter_code == "Aib"
        assert (
            amino_acid_namedtuple.smarts
            == "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4]([CH3X4])([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]"
        )
        assert (
            amino_acid_namedtuple.smarts_comments
            == "w normalnej sekwencji powinno sie outputowac jako {aMeAla}; a w zaiksowanej"
        )
        return

    def test_read_aa_smiles_dict(self):
        out_dict = read_aa_smiles_dict(self.aa_smiles_dict)
        assert out_dict["aMeAla"].smiles == "CC(C)(N)C=O"
        assert (
            out_dict["aMeAla"].smiles_radical
            == "CC(C)(N*)C(*)=O |$_CB;_CA;_Me;_N;_R1;_CO;_R2$|"
        )
        assert out_dict["aMeAla"].three_letter_code == "Aib"
        assert (
            out_dict["aMeAla"].smarts
            == "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4]([CH3X4])([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]"
        )
        assert (
            out_dict["aMeAla"].smarts_comments
            == "w normalnej sekwencji powinno sie outputowac jako {aMeAla}; a w zaiksowanej"
        )
        return

    def test_lower_dict(self):
        assert lower_dict({"Ala": "A", "Cys": "C"}) == {"ala": "a", "cys": "c"}

    def test_join_dicts(self):
        dicts = [
            {"Ala": "A", "Cys": "C"},
            {"ala": "a", "cys": "c"},
            {
                "dA": "a",
                "dC": "c",
                "dD": "d",
                "dE": "e",
                "dAla": "a",
                "dCys": "c",
                "dAsp": "d",
            },
            {"aMeAla": "X", "2-Aminoisobutyric acid": "X", "2-Methylalanine": "X"},
        ]

        joined_dict = {
            "Ala": "A",
            "Cys": "C",
            "ala": "a",
            "cys": "c",
            "dA": "a",
            "dC": "c",
            "dD": "d",
            "dE": "e",
            "dAla": "a",
            "dCys": "c",
            "dAsp": "d",
            "aMeAla": "X",
            "2-Aminoisobutyric acid": "X",
            "2-Methylalanine": "X",
        }

        assert join_dicts(dicts) == joined_dict
        return


class TestFileSystemDbRepo(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fsdb_repo = FileSystemDbRepo()
        cls.fsdb_repo = fsdb_repo
        cls.db_json = {
            "l_proteogenic": ["A", "C", "D"],
            "l_proteogenic_3letter": {"Ala": "A", "Cys": "C", "Asp": "D"},
            "d_proteogenic_3letter": {"ala": "A", "cys": "C", "asp": "D"},
            "aa_codes": {"A": "ALA", "C": "CYS", "D": "ASP"},
            "d_proteogenic": ["a", "c", "d"],
            "aa_order": ["G", "A", "a", "V"],
            "aa_smiles_dict": {
                "aMeAla": {
                    "smiles": "CC(C)(N)C=O",
                    "smiles_radical": "CC(C)(N*)C(*)=O |$_CB;_CA;_Me;_N;_R1;_CO;_R2$|",
                    "three_letter_code": "Aib",
                    "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4]([CH3X4])([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
                    "smarts_comments": "w normalnej sekwencji powinno sie outputowac jako {aMeAla}; a w zaiksowanej",
                },
                "A": {
                    "smiles": "C[C@H](N)C=O",
                    "smiles_radical": "C[C@H](N*)C(*)=O |$_CB;_CA;_N;_R1;_CO;_R2;$|",
                    "three_letter_code": "Ala",
                    "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
                    "smarts_comments": "",
                },
                "C": {
                    "smiles": "N[C@H](C=O)CS",
                    "smiles_radical": "N(*)[C@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
                    "three_letter_code": "Cys",
                    "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CHX4]([CH2X4][SX2H,SX1H0-,SX2])[CX3](=[OX1])[OX2H,OX1-,N]",
                    "smarts_comments": "Hits acid and conjugate base",
                },
                "C_nonmod": {
                    "smiles": "N[C@H](C=O)CS",
                    "smiles_radical": "N(*)[C@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
                    "three_letter_code": "Cys",
                    "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][SX2H,SX1H0-])[CX3](=[OX1])[OX2H,OX1-,N]",
                    "smarts_comments": "Hits acid and conjugate base",
                },
                "D": {
                    "smiles": "N[C@H](C=O)CC(=O)O",
                    "smiles_radical": "N(*)[C@H](C(*)=O)CC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;$|",
                    "three_letter_code": "Asp",
                    "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][CX3](=[OX1])[OH0-,OH])[CX3](=[OX1])[OX2H,OX1-,N]",
                    "smarts_comments": "Hits acid and conjugate base. Also hits Glu side chain when used alone.",
                },
                "a": {
                    "smiles": "C[C@@H](N)C=O",
                    "smiles_radical": "C[C@@H](N*)C(*)=O |$_CB;_CA;_N;_R1;_CO;_R2;$|",
                    "three_letter_code": "ala",
                    "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
                    "smarts_comments": "",
                },
                "c": {
                    "smiles": "N[C@@H](C=O)CS",
                    "smiles_radical": "N(*)[C@@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
                    "three_letter_code": "cys",
                    "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][SX2H,SX1H0-,SX2])[CX3](=[OX1])[OX2H,OX1-,N]",
                    "smarts_comments": "Hits acid and conjugate base",
                },
                "c_nonmod": {
                    "smiles": "N[C@@H](C=O)CS",
                    "smiles_radical": "N(*)[C@@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
                    "three_letter_code": "cys",
                    "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][SX2H,SX1H0-])[CX3](=[OX1])[OX2H,OX1-,N]",
                    "smarts_comments": "Hits acid and conjugate base",
                },
                "d": {
                    "smiles": "N[C@@H](C=O)CC(=O)O",
                    "smiles_radical": "N(*)[C@@H](C(*)=O)CC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;$|",
                    "three_letter_code": "asp",
                    "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][CX3](=[OX1])[OH0-,OH])[CX3](=[OX1])[OX2H,OX1-,N]",
                    "smarts_comments": "Hits acid and conjugate base. Also hits Glu side chain when used alone.",
                },
            },
            "aa_smarts_dict": {
                "A": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
                "R": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([CH2X4][CH2X4][CH2X4][NHX3][CH0X3](=[NH2X3+,NHX2+0])[NH2X3])[CX3](=[OX1])[OX2H,OX1-,N]",
                "N": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([CH2X4][CX3](=[OX1])[NX3H2])[CX3](=[OX1])[OX2H,OX1-,N]",
            },
            "d_proteogenic_2letter": {"dA": "a", "dC": "c", "dD": "d"},
            "d_proteogenic_4letter": {"dAla": "a", "dCys": "c", "dAsp": "d"},
            "modified_aa_codes": {
                "aMeAla": "X",
                "2-Aminoisobutyric acid": "X",
                "2-Methylalanine": "X",
            },
            "modified_aa_codes_reverse": {"X": "aMeAla"},
            "modified_aa_smiles_dict": {
                "aMeAla": {
                    "smiles": "CC(C)(N)C=O",
                    "smiles_radical": "CC(C)(N*)C(*)=O |$_CB;_CA;_Me;_N;_R1;_CO;_R2$|",
                    "three_letter_code": "Aib",
                    "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4]([CH3X4])([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
                    "smarts_comments": "w normalnej sekwencji powinno sie outputowac jako {aMeAla}",
                }
            },
            "n_terms_order": ["H", "Ac", "Pro", "Myr", "Palm"],
            "n_terms_smi_codes": {
                "H": "[*][H] |$_R2;_CO$|",
                "Pro": "CCC(=O)(*) |$_C;_C;_CO;_O;_R$|",
                "Ac": "CC(=O)(*) |$_C;_CO;_O;_R$|",
                "Myr": "CCCCCCCCCCCCCC(=O)(*) |$C;C;C;C;C;C;C;C;C;C;C;C;C;_CO;_O;_R$|",
                "Palm": "CCCCCCCCCCCCCCCC(=O)(*) |$C;C;C;C;C;C;C;C;C;C;C;C;C;C;C;_CO;_O;_R$|",
            },
            "mod_smi_codes": {"Ac": "CC=O"},
            "c_terms_order": ["OH", "NH2", "N-Me"],
            "c_terms_smi_codes": {
                "NH2": "N(*) |$_N;_R$|",
                "OH": "[OH](*) |$_OH;_R$|",
                "N-Me": "CN(*) |$_C;_N;_R$|",
            },
            "substitute_smiles": {"CH3": "C"},
            "substitute_dict": {"aMe": "CH3"},
            "protein_letters": "ACDEFGHIKLMNPQRSTVWXY",
            "hydrophobicity_index_pH2": {"L": 100, "I": 100, "F": 92},
            "hydrophobicity_index_pH7": {"F": 100, "I": 99, "W": 97},
            "acidic_aa": {"D": "", "E": ""},
            "basic_aa": {"R": "", "K": "", "H": ""},
            "protein_weights": {"A": 89.0932, "C": 121.1582, "D": 133.1027},
        }

        return

    def test_read_json(self):
        fsdb_repo = FileSystemDbRepo()
        fsdb_repo.read_json(self.db_json)

        assert fsdb_repo.n_terms_smi_codes == {
            "H": "[*][H] |$_R2;_CO$|",
            "Pro": "CCC(=O)(*) |$_C;_C;_CO;_O;_R$|",
            "Ac": "CC(=O)(*) |$_C;_CO;_O;_R$|",
            "Myr": "CCCCCCCCCCCCCC(=O)(*) |$C;C;C;C;C;C;C;C;C;C;C;C;C;_CO;_O;_R$|",
            "Palm": "CCCCCCCCCCCCCCCC(=O)(*) |$C;C;C;C;C;C;C;C;C;C;C;C;C;C;C;_CO;_O;_R$|",
        }

        assert fsdb_repo.n_terms_order == ["H", "Ac", "Pro", "Myr", "Palm"]
        assert fsdb_repo.c_terms_smi_codes == {
            "NH2": "N(*) |$_N;_R$|",
            "OH": "[OH](*) |$_OH;_R$|",
            "N-Me": "CN(*) |$_C;_N;_R$|",
        }
        assert fsdb_repo.c_terms_order == ["OH", "NH2", "N-Me"]
        assert fsdb_repo.substitute_smiles == {"CH3": "C"}
        assert fsdb_repo.substitute_dict == {"aMe": "CH3"}
        assert fsdb_repo.aa_order == ["G", "A", "a", "V"]
        assert fsdb_repo.mod_smi_codes == {"Ac": "CC=O"}

        assert fsdb_repo.l_proteogenic_3letter == {"Ala": "A", "Cys": "C", "Asp": "D"}

        assert fsdb_repo.d_proteogenic_2letter == {"dA": "a", "dC": "c", "dD": "d"}
        assert fsdb_repo.d_proteogenic_4letter == {
            "dAla": "a",
            "dCys": "c",
            "dAsp": "d",
        }

        assert fsdb_repo.modified_aa_codes == {
            "aMeAla": "X",
            "2-Aminoisobutyric acid": "X",
            "2-Methylalanine": "X",
        }

        assert fsdb_repo.file_aa_smiles_dict == {
            "aMeAla": {
                "smiles": "CC(C)(N)C=O",
                "smiles_radical": "CC(C)(N*)C(*)=O |$_CB;_CA;_Me;_N;_R1;_CO;_R2$|",
                "three_letter_code": "Aib",
                "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4]([CH3X4])([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
                "smarts_comments": "w normalnej sekwencji powinno sie outputowac jako {aMeAla}; a w zaiksowanej",
            },
            "A": {
                "smiles": "C[C@H](N)C=O",
                "smiles_radical": "C[C@H](N*)C(*)=O |$_CB;_CA;_N;_R1;_CO;_R2;$|",
                "three_letter_code": "Ala",
                "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
                "smarts_comments": "",
            },
            "C": {
                "smiles": "N[C@H](C=O)CS",
                "smiles_radical": "N(*)[C@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
                "three_letter_code": "Cys",
                "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CHX4]([CH2X4][SX2H,SX1H0-,SX2])[CX3](=[OX1])[OX2H,OX1-,N]",
                "smarts_comments": "Hits acid and conjugate base",
            },
            "C_nonmod": {
                "smiles": "N[C@H](C=O)CS",
                "smiles_radical": "N(*)[C@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
                "three_letter_code": "Cys",
                "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][SX2H,SX1H0-])[CX3](=[OX1])[OX2H,OX1-,N]",
                "smarts_comments": "Hits acid and conjugate base",
            },
            "D": {
                "smiles": "N[C@H](C=O)CC(=O)O",
                "smiles_radical": "N(*)[C@H](C(*)=O)CC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;$|",
                "three_letter_code": "Asp",
                "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][CX3](=[OX1])[OH0-,OH])[CX3](=[OX1])[OX2H,OX1-,N]",
                "smarts_comments": "Hits acid and conjugate base. Also hits Glu side chain when used alone.",
            },
            "a": {
                "smiles": "C[C@@H](N)C=O",
                "smiles_radical": "C[C@@H](N*)C(*)=O |$_CB;_CA;_N;_R1;_CO;_R2;$|",
                "three_letter_code": "ala",
                "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
                "smarts_comments": "",
            },
            "c": {
                "smiles": "N[C@@H](C=O)CS",
                "smiles_radical": "N(*)[C@@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
                "three_letter_code": "cys",
                "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][SX2H,SX1H0-,SX2])[CX3](=[OX1])[OX2H,OX1-,N]",
                "smarts_comments": "Hits acid and conjugate base",
            },
            "c_nonmod": {
                "smiles": "N[C@@H](C=O)CS",
                "smiles_radical": "N(*)[C@@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
                "three_letter_code": "cys",
                "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][SX2H,SX1H0-])[CX3](=[OX1])[OX2H,OX1-,N]",
                "smarts_comments": "Hits acid and conjugate base",
            },
            "d": {
                "smiles": "N[C@@H](C=O)CC(=O)O",
                "smiles_radical": "N(*)[C@@H](C(*)=O)CC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;$|",
                "three_letter_code": "asp",
                "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][CX3](=[OX1])[OH0-,OH])[CX3](=[OX1])[OX2H,OX1-,N]",
                "smarts_comments": "Hits acid and conjugate base. Also hits Glu side chain when used alone.",
            },
        }

        return

    def test_process_json(self):
        fsdb_repo = FileSystemDbRepo()
        fsdb_repo.read_json(self.db_json)
        fsdb_repo.process_json()

        assert fsdb_repo.l_proteogenic == set(
            ["aMeAla", "A", "C", "C_nonmod", "D", "a", "c", "c_nonmod", "d"]
        )

        assert fsdb_repo.d_proteogenic == set(
            ["ameala", "a", "c", "c_nonmod", "d", "a", "c", "c_nonmod", "d"]
        )
        # aa_smiles_dict?

        assert fsdb_repo.aa_smiles_dict == {
            "aMeAla": AminoAcid(
                smiles="CC(C)(N)C=O",
                smiles_radical="CC(C)(N*)C(*)=O |$_CB;_CA;_Me;_N;_R1;_CO;_R2$|",
                three_letter_code="Aib",
                smarts="[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4]([CH3X4])([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
                smarts_comments="w normalnej sekwencji powinno sie outputowac jako {aMeAla}; a w zaiksowanej",
            ),
            "A": AminoAcid(
                smiles="C[C@H](N)C=O",
                smiles_radical="C[C@H](N*)C(*)=O |$_CB;_CA;_N;_R1;_CO;_R2;$|",
                three_letter_code="Ala",
                smarts="[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
                smarts_comments="",
            ),
            "C": AminoAcid(
                smiles="N[C@H](C=O)CS",
                smiles_radical="N(*)[C@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
                three_letter_code="Cys",
                smarts="[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CHX4]([CH2X4][SX2H,SX1H0-,SX2])[CX3](=[OX1])[OX2H,OX1-,N]",
                smarts_comments="Hits acid and conjugate base",
            ),
            "C_nonmod": AminoAcid(
                smiles="N[C@H](C=O)CS",
                smiles_radical="N(*)[C@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
                three_letter_code="Cys",
                smarts="[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][SX2H,SX1H0-])[CX3](=[OX1])[OX2H,OX1-,N]",
                smarts_comments="Hits acid and conjugate base",
            ),
            "D": AminoAcid(
                smiles="N[C@H](C=O)CC(=O)O",
                smiles_radical="N(*)[C@H](C(*)=O)CC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;$|",
                three_letter_code="Asp",
                smarts="[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][CX3](=[OX1])[OH0-,OH])[CX3](=[OX1])[OX2H,OX1-,N]",
                smarts_comments="Hits acid and conjugate base. Also hits Glu side chain when used alone.",
            ),
            "a": AminoAcid(
                smiles="C[C@@H](N)C=O",
                smiles_radical="C[C@@H](N*)C(*)=O |$_CB;_CA;_N;_R1;_CO;_R2;$|",
                three_letter_code="ala",
                smarts="[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
                smarts_comments="",
            ),
            "c": AminoAcid(
                smiles="N[C@@H](C=O)CS",
                smiles_radical="N(*)[C@@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
                three_letter_code="cys",
                smarts="[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][SX2H,SX1H0-,SX2])[CX3](=[OX1])[OX2H,OX1-,N]",
                smarts_comments="Hits acid and conjugate base",
            ),
            "c_nonmod": AminoAcid(
                smiles="N[C@@H](C=O)CS",
                smiles_radical="N(*)[C@@H](C(*)=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
                three_letter_code="cys",
                smarts="[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][SX2H,SX1H0-])[CX3](=[OX1])[OX2H,OX1-,N]",
                smarts_comments="Hits acid and conjugate base",
            ),
            "d": AminoAcid(
                smiles="N[C@@H](C=O)CC(=O)O",
                smiles_radical="N(*)[C@@H](C(*)=O)CC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;$|",
                three_letter_code="asp",
                smarts="[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][CX3](=[OX1])[OH0-,OH])[CX3](=[OX1])[OX2H,OX1-,N]",
                smarts_comments="Hits acid and conjugate base. Also hits Glu side chain when used alone.",
            ),
        }
        assert fsdb_repo.d_proteogenic_3letter == {"ala": "a", "cys": "c", "asp": "d"}

        assert fsdb_repo.three_letter == {
            "Ala": "A",
            "Cys": "C",
            "Asp": "D",
            "ala": "a",
            "cys": "c",
            "asp": "d",
        }

        assert fsdb_repo.symbols == {
            "dA": "a",
            "dC": "c",
            "dD": "d",
            "dAla": "a",
            "dCys": "c",
            "dAsp": "d",
            "Ala": "A",
            "Cys": "C",
            "Asp": "D",
            "ala": "a",
            "cys": "c",
            "asp": "d",
            "aMeAla": "X",
            "2-Aminoisobutyric acid": "X",
            "2-Methylalanine": "X",
        }
        return
