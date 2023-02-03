dermcidin_stapled_with_anchor = "CC(C)C[C@H](NC(=O)[C@@H](NC(=O)[C@H](CO)NC(=O)[C@H](CC(O)=O)NC(=O)[C@H](CC(C)C)NC(=O)[C@@H](NC(=O)[C@H](CC(O)=O)NC(=O)[C@H](CCCCN)NC(=O)[C@@H](NC(=O)[C@H](CC(O)=O)NC(=O)[C@H](Cc1c[nH]cn1)NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)CNC(=O)[C@H](CCCCN)NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](CO)NC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(O)=O)NC(=O)[C@H](CCC(O)=O)NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)[C@H](CC(O)=O)NC(=O)[C@H](CCCCN)NC(=O)[C@H]1CSCC(=O)N[C@@H](CCCNC(=O)CSC[C@@H](NC(=O)[C@@H](NC(=O)[C@H](C)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CCCCN)NC(=O)C(C)(C)NC(=O)CNC(=O)[C@H](CC(O)=O)NC(=O)[C@H](CC(C)C)NC(=O)CNC(=O)[C@H](CCCCN)NC(=O)[C@H](CCC(O)=O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CO)NC(=O)[C@@H](N)CO)C(C)C)C(=O)NCC(=O)N[C@@H](CC(C)C)C(=O)NCC(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(C)C)C(=O)N1)C(=O)NCCOCCO)C(C)C)C(C)C)C(C)C)C(C)C)C(C)C)C(C)C)C(O)=O"
dermcidin_stapled_with_anchor_seq = (
    "SSLLEKGLDG{aMeAla}KKAV{modC}GLGKL{modC}KDAVEDLESVGKGAVHDVKDVLDSVL"
)
seq_wo_mods = "SSLLEKGLDG{aMeAla}KKAVCGLGKLCKDAVEDLESVGKGAVHDVKDVLDSVL"


physico_chemical_properties = {
    "hydrophobicity_index_pH2": {
        "C": 53,
        "A": 47,
        "G": 0,
        "S": -7,
        "K": -37,
    },
    "hydrophobicity_index_pH7": {"C": 49, "A": 41, "G": 0, "S": -5, "K": -23, "D": -55},
    "acidic_aa": {"D": "", "E": ""},
    "basic_aa": {"R": "", "K": "", "H": ""},
    "protein_weights": {
        "A": 89.0932,
        "C": 121.1582,
        "D": 133.1027,
        "K": 146.1876,
        "S": 105.0926,
    },
}


n_terms = {
    "H": {
        "smiles": "[H]",
        "smiles_radical": "[1*][H]",
        "three_letter_code": "",
        "smarts": "[H]",
        "smarts_comments": "",
    },
    "proton": {
        "smiles": "[H]",
        "smiles_radical": "[1*][H]",
        "three_letter_code": "",
        "smarts": "[H]",
        "smarts_comments": "",
    },
    "Pro": {
        "smiles": "CCC(=O)",
        "smiles_radical": "CCC(=O)([1*]) |$_C;_C;_CO;_O;_R$|",
        "three_letter_code": "Pro",
        "smarts": "CCC(=O)",
        "smarts_comments": "",
    },
    "Ac": {
        "smiles": "CC(=O)",
        "smiles_radical": "CC(=O)([1*]) |$_C;_CO;_O;_R$|",
        "three_letter_code": "Pro",
        "smarts": "CC(=O)",
        "smarts_comments": "",
    },
    "Myr": {
        "smiles": "CCCCCCCCCCCCCC(=O)",
        "smiles_radical": "CCCCCCCCCCCCCC(=O)([1*]) |$C;C;C;C;C;C;C;C;C;C;C;C;C;_CO;_O;_R$|",
        "three_letter_code": "Myr",
        "smarts": "CCCCCCCCCCCCCC(=O)",
        "smarts_comments": "",
    },
    "Palm": {
        "smiles": "CCCCCCCCCCCCCCCC(=O)",
        "smiles_radical": "CCCCCCCCCCCCCCCC(=O)([1*]) |$C;C;C;C;C;C;C;C;C;C;C;C;C;C;C;_CO;_O;_R$|",
        "three_letter_code": "Myr",
        "smarts": "CCCCCCCCCCCCCCCC(=O)",
        "smarts_comments": "",
    },
}


c_terms = {
    "NH2": {
        "smiles": "[NH2]",
        "smiles_radical": "[1*][NH2]",
        "three_letter_code": "NH2",
        "smarts": "[NH2]",
        "smarts_comments": "",
    },
    "OH": {
        "smiles": "O",
        "smiles_radical": "[1*]O",
        "three_letter_code": "OH",
        "smarts": "O",
        "smarts_comments": "",
    },
    "N-Me": {
        "smiles": "CN",
        "smiles_radical": "CN([1*])",
        "three_letter_code": "N-Me",
        "smarts": "CN",
        "smarts_comments": "",
    },
}


aa_smiles_dict = {
    "aMeAla": {
        "smiles": "CC(C)(N)C=O",
        "smiles_radical": "CC(C)(N[1*])C([2*])=O |$_CB;_CA;_Me;_N;_R1;_CO;_R2$|",
        "three_letter_code": "Aib",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4]([CH3X4])([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "w normalnej sekwencji powinno sie outputowac jako {aMeAla}; a w zaiksowanej",
    },
    "A": {
        "smiles": "C[C@H](N)C=O",
        "smiles_radical": "C[C@H](N[1*])C([2*])=O |$_CB;_CA;_N;_R1;_CO;_R2;$|",
        "three_letter_code": "Ala",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "Hits acid and conjugate base",
    },
    "C": {
        "smiles": "N[C@H](C=O)CS",
        "smiles_radical": "[1*]N[C@@H](CS)C([2*])=O |atomProp:0.dummyLabel.1*:4.AtomName.SG:4.atomLabel.SG:6.dummyLabel.2*|",
        "three_letter_code": "Cys",
        "smarts": "[$([N&X3&H2,N&X4&H3&+]),$([N&X3&H1](C)C)][C&H1&X4]([C&H2&X4][S&X2&H1,S&X1&H0&-,S&X2])[C&X3](=[O&X1])[O&X2&H1,O&X1&-,N] |atomProp:3.AtomName.SG|",
        "smarts_comments": "Hits acid and conjugate base",
    },
    "C_nonmod": {
        "smiles": "N[C@H](C=O)CS",
        "smiles_radical": "N([1*])[C@H](C([2*])=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
        "three_letter_code": "Cys",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][SX2H,SX1H0-])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "Hits acid and conjugate base",
    },
    "D": {
        "smiles": "N[C@H](C=O)CC(=O)O",
        "smiles_radical": "N([1*])[C@H](C([2*])=O)CC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;$|",
        "three_letter_code": "Asp",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][CX3](=[OX1])[OH0-,OH])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "Hits acid and conjugate base. Also hits Glu side chain when used alone.",
    },
    "G": {
        "smiles": "NCC=O",
        "smiles_radical": "N([1*])CC([2*])=O",
        "three_letter_code": "Gly",
        "smarts": "N[CX4H2][CX3](=[OX1])[O,N]",
        "smarts_comments": "",
    },
    "g": {
        "smiles": "NCC=O",
        "smiles_radical": "N([1*])CC([2*])=O |$_N;_R1;_CA;_CO;_R2$|",
        "three_letter_code": "Gly",
        "smarts": "N[CX4H2][CX3](=[OX1])[O,N]",
        "smarts_comments": "",
    },
    "K": {
        "smiles": "NCCCC[C@H](N)C=O",
        "smiles_radical": "[1*]N[C@@H](CCCCN)C([2*])=O |atomProp:0.dummyLabel.1*:7.AtomName.NZ:9.dummyLabel.2*|",
        "three_letter_code": "Lys",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "Acid and conjugate base",
    },
    "S": {
        "smiles": "N[C@H](C=O)CO",
        "smiles_radical": "N([1*])[C@H](C([2*])=O)CO |$_N;_R1;_CA;_CO;_R2;;_CB$|",
        "three_letter_code": "Ser",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][OX2H])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "",
    },
    "aMeSer": {
        "smiles": "CN[C@H](CO)C(=O)O",
        "smiles_radical": "N([1*])[C@H](C([2*])=O)CO |$_N;_R1;_CA;_CO;_R2;;_CB$|",
        "three_letter_code": "MSe",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@CX4]([CH3X4])([CH2X4][OX2H])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "",
    },
    "a": {
        "smiles": "C[C@@H](N)C=O",
        "smiles_radical": "C[C@@H](N[1*])C([2*])=O |$_CB;_CA;_N;_R1;_CO;_R2;$|",
        "three_letter_code": "ala",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH3X4])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "",
    },
    "c": {
        "smiles": "N[C@@H](C=O)CS",
        "smiles_radical": "N([1*])[C@@H](C([2*])=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
        "three_letter_code": "cys",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@@HX4]([CH2X4][SX2H,SX1H0-,SX2])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "Hits acid and conjugate base",
    },
    "c_nonmod": {
        "smiles": "N[C@@H](C=O)CS",
        "smiles_radical": "N([1*])[C@@H](C([2*])=O)CS |$_N;_R1;_CA;_CO;_R2;O;_CB;$|",
        "three_letter_code": "cys",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][SX2H,SX1H0-])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "Hits acid and conjugate base",
    },
    "d": {
        "smiles": "N[C@@H](C=O)CC(=O)O",
        "smiles_radical": "N([1*])[C@@H](C([2*])=O)CC(=O)O |$_N;_R1;_CA;_CO;_R2;;_CB;;$|",
        "three_letter_code": "asp",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][CX3](=[OX1])[OH0-,OH])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "Hits acid and conjugate base. Also hits Glu side chain when used alone.",
    },
    "k": {
        "smiles": "NCCCC[C@@H](N)C=O",
        "smiles_radical": "NCCCC[C@@H](N([1*]))C([2*])=O |$;;;;_CB;_CA;_N;_R1;_CO;_R2;$|",
        "three_letter_code": "lys",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][CH2X4][CH2X4][CH2X4][NX4+,NX3+0])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "Acid and conjugate base",
    },
    "s": {
        "smiles": "N[C@@H](C=O)CO",
        "smiles_radical": "N([1*])[C@@H](C([2*])=O)CO |$_N;_R1;_CA;_CO;_R2;;_CB$|",
        "three_letter_code": "ser",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][C@HX4]([CH2X4][OX2H])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "",
    },
    "aMeX": {
        "smiles": "CN[C@H](CO)C(=O)O",
        "smiles_radical": "N([1*])[C@H](C([2*])=O)CO |$_N;_R1;_CA;_CO;_R2;;_CB$|",
        "three_letter_code": "MSe",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4C]([CH3X4])([*])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "",
    },
    "X": {
        "smiles": "CN[C@H](CO)C(=O)O",
        "smiles_radical": "N([1*])[C@H](C([2*])=O)CO |$_N;_R1;_CA;_CO;_R2;;_CB$|",
        "three_letter_code": "MSe",
        "smarts": "[$([NX3H2,NX4H3+]),$([NX3H](C)(C))][CX4H]([*])[CX3](=[OX1])[OX2H,OX1-,N]",
        "smarts_comments": "",
    },
}


db_json = {
    "l_proteogenic": [
        "A",
        "C",
        "D",
        "K",
        "S",
        "T",
    ],
    "coding": {
        "l_proteogenic_3letter": {
            "Ala": "A",
            "Cys": "C",
            "Asp": "D",
            "Gly": "G",
            "Lys": "K",
            "Ser": "S",
        },
        "d_proteogenic_3letter": {
            "ala": "a",
            "cys": "c",
            "asp": "d",
            "gly": "g",
            "lys": "k",
            "ser": "s",
        },
        "aa_codes": {
            "A": "ALA",
            "C": "CYS",
            "D": "ASP",
            "G": "GLY",
            "K": "LYS",
            "S": "SER",
        },
        "d_proteogenic_2letter": {
            "dA": "a",
            "dC": "c",
            "dD": "d",
            "dG": "g",
            "dK": "k",
            "dS": "s",
        },
        "d_proteogenic_4letter": {
            "dAla": "a",
            "dCys": "c",
            "dAsp": "d",
            "dGly": "g",
            "dLys": "k",
            "dSer": "s",
        },
        "modified_aa_codes": {
            "aMeAla": "X",
            "2-Aminoisobutyric acid": "X",
            "2-Methylalanine": "X",
        },
        "modified_aa_codes_reverse": {"X": "aMeAla"},
    },
    "d_proteogenic": [
        "a",
        "c",
        "d",
        "g",
        "k",
        "s",
    ],
    "aa_order": [
        "G",
        "A",
        "a",
        "C",
        "c",
        "D",
        "d",
        "K",
        "k",
        "S",
        "s",
    ],
    "smiles": {
        "n_terms": n_terms,
        "c_terms": c_terms,
        "aa": aa_smiles_dict,
    },
    "n_terms_order": ["H", "Ac", "Pro", "Myr", "Palm"],
    "c_terms_order": ["OH", "NH2", "N-Me"],
    "protein_letters": "ACDEFGHIKLMNPQRSTVWXY",
    "physico_chemical_properties": physico_chemical_properties,
}
