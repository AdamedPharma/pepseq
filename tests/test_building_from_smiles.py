from BuildPeptideJSONFromSMILES import decompose_peptide_smiles

c_smarts = "[$([N&X3&H2,N&X4&H3&+]),$([N&X3&H1](C)C)][C&H1&X4]([C&H2&X4][S&X2&H1,S&X1&H0&-,S&X2])[C](=[O&X1]) |atomProp:3.AtomName.SG|"
c_smarts = "[$([N&X3&H2,N&X4&H3&+]),$([N&X3&H1](C)C)][C&H1&X4]([C&H2&X4][S&X2&H1,S&X1&H0&-,S&X2])[C](=[O&X1]) |atomProp:3.AtomName.SG|"

K_smarts = "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0])[CX3](=[OX1])"
A_smarts = "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([CH3X4])[CX3](=[OX1])"
S_smarts = "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][OX2H])[CX3](=[OX1])"
G_smarts = "N[CX4H2][CX3](=[OX1])[O,N]"
# NH jest stapleowana

cx_smarts_db = {
    "C": c_smarts,
    "G": "N[CX4H2][CX3](=[OX1])",  # [O,N]
    "A": A_smarts,
    "S": S_smarts,
    "K": K_smarts,
}

smiles = "[H]N[C@H]1CSC(Br)CNP([Na])SC[C@@H](C(=O)N[C@@H](C)C(=O)N[C@H]2CSSC[C@@H](C(=O)N[C@@H](CCCCN)C(=O)O)NC(=O)CNC2=O)NC(=O)[C@H](CO)NC1=O"


def test_decompose_peptide_smiles():
    peptide_json = decompose_peptide_smiles(smiles, cx_smarts_db)
    assert peptide_json == {
        "sequence": "CSCACGCK",
        "internal_modifications": [
            {
                1: [
                    {"ResID": "5", "AtomName": "SG", "ResidueName": ""},
                    {"ResID": "7", "AtomName": "SG", "ResidueName": ""},
                ]
            }
        ],
        "external_modifications": [
            {
                "smiles": "[1*]C(Br)CNP([2*])[Na]",
                "max_attachment_point_id": 2,
                "attachment_points_on_sequence": {
                    1: {
                        "attachment_point_id": 1,
                        "ResID": "1",
                        "AtomName": "SG",
                        "ResidueName": "",
                    },
                    2: {
                        "attachment_point_id": 2,
                        "ResID": "3",
                        "AtomName": "SG",
                        "ResidueName": "",
                    },
                },
            },
            {
                "smiles": "[1*]O",
                "max_attachment_point_id": 1,
                "attachment_points_on_sequence": {
                    1: {
                        "attachment_point_id": 1,
                        "ResID": "8",
                        "AtomName": "CO",
                        "ResidueName": "",
                    }
                },
            },
        ],
    }
