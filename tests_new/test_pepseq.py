import pytest
import rdkit

from read import from_pepseq

s1l_vs_smiles = (
    (
        "QWERTYIPASDFGHKLCVNM",
        "N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)N[C@@]([H])(CCC(=O)O)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N[C@@]([H])([C@]([H])(O)C)C(=O)N[C@@]([H])(Cc1ccc(O)cc1)C(=O)N[C@@]([H])([C@]([H])(CC)C)C(=O)N1[C@@]([H])(CCC1)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CC(=O)O)C(=O)N[C@@]([H])(Cc1ccccc1)C(=O)NCC(=O)N[C@@]([H])(CC1=CN=C-N1)C(=O)N[C@@]([H])(CCCCN)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(CS)C(=O)N[C@@]([H])(C(C)C)C(=O)N[C@@]([H])(CC(=O)N)C(=O)N[C@@]([H])(CCSC)C(=O)O",
    ),
    (
        "qwertyipasdfghklcvnm",
        "N[C@]([H])(CCC(=O)N)C(=O)N[C@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)N[C@]([H])(CCC(=O)O)C(=O)N[C@]([H])(CCCNC(=N)N)C(=O)N[C@]([H])([C@@]([H])(O)C)C(=O)N[C@]([H])(Cc1ccc(O)cc1)C(=O)N[C@]([H])([C@@]([H])(CC)C)C(=O)N1[C@]([H])(CCC1)C(=O)N[C@]([H])(C)C(=O)N[C@]([H])(CO)C(=O)N[C@]([H])(CC(=O)O)C(=O)N[C@]([H])(Cc1ccccc1)C(=O)NCC(=O)N[C@]([H])(CC1=CN=C-N1)C(=O)N[C@]([H])(CCCCN)C(=O)N[C@]([H])(CC(C)C)C(=O)N[C@]([H])(CS)C(=O)N[C@]([H])(C(C)C)C(=O)N[C@]([H])(CC(=O)N)C(=O)N[C@]([H])(CCSC)C(=O)O",
    ),
)

s1l_vs_canonical_smiles = (
    ("ATA", rdkit.Chem.MolToSmiles(rdkit.Chem.MolFromSequence("ATA"))),
    (
        "QWERTYIPASDFGHKLCVNM",
        rdkit.Chem.MolToSmiles(rdkit.Chem.MolFromSequence("QWERTYIPASDFGHKLCVNM")),
    ),
    (
        "QWERTYIPASDFGHKLCVNM",
        "CC[C@H](C)[C@H](NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@@H](NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)[C@@H](N)CCC(N)=O)[C@@H](C)O)C(=O)N1CCC[C@H]1C(=O)N[C@@H](C)C(=O)N[C@@H](CO)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](Cc1ccccc1)C(=O)NCC(=O)N[C@@H](Cc1cnc[nH]1)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CS)C(=O)N[C@H](C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H](CCSC)C(=O)O)C(C)C",
    ),
    (
        "qwertyipasdfghklcvnm",
        rdkit.Chem.MolToSmiles(
            rdkit.Chem.MolFromFASTA(""">\nqwertyipasdfghklcvnm\n""", flavor=1)
        ),
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

pepseq_vs_smiles_moded = (
    (
        "ALA{modX}",
        "N[C@@]([H])(C)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(C)C(=O)N[C@@]([H])(CCCN)C(=O)O",
    ),
    # (
    #     "qwertyipasdfghklcvnm",
    #     "N[C@]([H])(CCC(=O)N)C(=O)N[C@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)N[C@]([H])(CCC(=O)O)C(=O)N[C@]([H])(CCCNC(=N)N)C(=O)N[C@]([H])([C@@]([H])(O)C)C(=O)N[C@]([H])(Cc1ccc(O)cc1)C(=O)N[C@]([H])([C@@]([H])(CC)C)C(=O)N1[C@]([H])(CCC1)C(=O)N[C@]([H])(C)C(=O)N[C@]([H])(CO)C(=O)N[C@]([H])(CC(=O)O)C(=O)N[C@]([H])(Cc1ccccc1)C(=O)NCC(=O)N[C@]([H])(CC1=CN=C-N1)C(=O)N[C@]([H])(CCCCN)C(=O)N[C@]([H])(CC(C)C)C(=O)N[C@]([H])(CS)C(=O)N[C@]([H])(C(C)C)C(=O)N[C@]([H])(CC(=O)N)C(=O)N[C@]([H])(CCSC)C(=O)O",
    # ),
    (
        "H~ADRPE~OH",
        "C[C@H](N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N1CCC[C@@H]1C(=O)N[C@@H](CCC(=O)O)C(=O)O",
    ),
)

import rdkit


def smiles_are_identical(smi1, smi2):
    a = rdkit.Chem.MolFromSmiles(smi1)
    b = rdkit.Chem.MolFromSmiles(smi2)
    are_identical = a.HasSubstructMatch(b, useChirality=True) and b.HasSubstructMatch(
        a, useChirality=True
    )
    return are_identical


@pytest.mark.parametrize(
    "pepseq,smiles", pepseq_vs_canonical_smiles + s1l_vs_canonical_smiles
)
def test_peptide_from_pepseq(pepseq, smiles):
    peptide = from_pepseq(pepseq)
    # assert peptide.smiles == smiles
    assert smiles_are_identical(peptide.smiles, smiles)


# ValidationError
# peptide = PeptideFromPepSeq("{Ule}{Lys}", db_api=db_api)

# peptide = PeptideFromPepSeq("{Ala}{Lys}", db_api=db_api)
# peptide = PeptideFromPepSeq("H~P{dLeu}RT{Ala}{Lys}~NH2", db_api=db_api)


# print(peptide.smiles)
# print(peptide.n_term)
# print(peptide.c_term)
# print(peptide.sequence)
