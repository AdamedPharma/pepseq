import rdkit


s1l_vs_canonical_smiles = (
    ("ATA", rdkit.Chem.MolToSmiles(rdkit.Chem.MolFromSequence("ATA"))),
    (
        "QWERTYIPASDFGHKLCVNM",
        rdkit.Chem.MolToSmiles(rdkit.Chem.MolFromSequence("QWERTYIPASDFGHKLCVNM")),
    ),
    (
        "QWERTYIPASDFGHKLCVNM",
        "CC[C@H](C)[C@H](NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@@H](NC(=O)"
        + "[C@H](CCCNC(=N)N)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](Cc1c[nH]c2"
        + "ccccc12)NC(=O)[C@@H](N)CCC(N)=O)[C@@H](C)O)C(=O)N1CCC[C@H]1C"
        + "(=O)N[C@@H](C)C(=O)N[C@@H](CO)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@"
        + "@H](Cc1ccccc1)C(=O)NCC(=O)N[C@@H](Cc1cnc[nH]1)C(=O)N[C@@H](C"
        + "CCCN)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CS)C(=O)N[C@H](C(=O)N["
        + "C@@H](CC(N)=O)C(=O)N[C@@H](CCSC)C(=O)O)C(C)C",
    ),
    (
        "qwertyipasdfghklcvnm",
        rdkit.Chem.MolToSmiles(
            rdkit.Chem.MolFromFASTA(""">\nqwertyipasdfghklcvnm\n""", flavor=1)
        ),
    ),
    (
        "[CSPS[*:1]]~Y{aMeAla}QGTFTSDYSKYLDECAAKDFVCWLLDHHPSSGQPPPS~[CNSC[*:1]]",
        "CNSCC(=O)[C@H](CO)NC(=O)[C@@H]1CCCN1C(=O)[C@@H]1CCCN1C(=O)[C@@H]1CCCN"
        + "1C(=O)[C@H](CCC(N)=O)NC(=O)CNC(=O)[C@H](CO)NC(=O)[C@H](CO)NC(=O)[C@"
        + "@H]1CCCN1C(=O)[C@H](Cc1c[nH]cn1)NC(=O)[C@H](Cc1c[nH]cn1)NC(=O)[C@H]"
        + "(CC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](Cc1c[nH]"
        + "c2ccccc12)NC(=O)[C@H](CS)NC(=O)[C@@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)["
        + "C@H](CC(=O)O)NC(=O)[C@H](CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C"
        + "@H](CS)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CC(C)C)"
        + "NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@H](CCCCN)NC(=O)[C@H](CO)NC(=O)[C@"
        + "H](Cc1ccc(O)cc1)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](CO)NC(=O)[C@@H](NC("
        + "=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](NC(=O)CNC(=O)[C@H](CCC(N)=O)NC(=O)C"
        + "(C)(C)NC(=O)[C@H](Cc1ccc(O)cc1)NSPSC)[C@@H](C)O)[C@@H](C)O)C(C)C",
    ),
)

correct_smi = rdkit.Chem.MolToSmiles(rdkit.Chem.MolFromSequence("ADRPE"))


pepseq_vs_canonical_smiles = (
    ("H~ADRPE~OH", correct_smi),
    ("H~AD{Arg}PE~OH", correct_smi),
    ("AD{Arg}PE", correct_smi),
    ("{Ala}D{Arg}P{Glu}", correct_smi),
    ("H~{Ala}D{Arg}P{Glu}~OH", correct_smi),
    ("ADRPE", correct_smi),
)

tests = pepseq_vs_canonical_smiles + s1l_vs_canonical_smiles
