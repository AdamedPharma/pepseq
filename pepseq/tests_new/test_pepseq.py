import json
import pkgutil

import pytest
import rdkit
from pepseq.BuildingModifiedPeptideFromPeptideJSON import \
    get_smiles_from_peptide_json
from pepseq.BuildPeptideJSONFromSMILES import \
    from_smiles_to_pepseq_and_one_mod_smiles_strings
from pepseq.functions import calculate
from pepseq.get_peptide_json_from_pepseq_format import get_pep_json
from pepseq.read import from_json, from_pepseq

db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


s1l_vs_smiles = (
    (
        "QWERTYIPASDFGHKLCVNM",
        "N[C@@]([H])(CCC(=O)N)C(=O)N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C"
        + "(=O)N[C@@]([H])(CCC(=O)O)C(=O)N[C@@]([H])(CCCNC(=N)N)C(=O)N["
        + "C@@]([H])([C@]([H])(O)C)C(=O)N[C@@]([H])(Cc1ccc(O)cc1)C(=O)N"
        + "[C@@]([H])([C@]([H])(CC)C)C(=O)N1[C@@]([H])(CCC1)C(=O)N[C@@]"
        + "([H])(C)C(=O)N[C@@]([H])(CO)C(=O)N[C@@]([H])(CC(=O)O)C(=O)N["
        + "C@@]([H])(Cc1ccccc1)C(=O)NCC(=O)N[C@@]([H])(CC1=CN=C-N1)C(=O"
        + ")N[C@@]([H])(CCCCN)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])("
        + "CS)C(=O)N[C@@]([H])(C(C)C)C(=O)N[C@@]([H])(CC(=O)N)C(=O)N[C@"
        + "@]([H])(CCSC)C(=O)O",
    ),
    (
        "qwertyipasdfghklcvnm",
        "N[C@]([H])(CCC(=O)N)C(=O)N[C@]([H])(CC(=CN2)C1=C2C=CC=C1)C(="
        + "O)N[C@]([H])(CCC(=O)O)C(=O)N[C@]([H])(CCCNC(=N)N)C(=O)N[C@]("
        + "[H])([C@@]([H])(O)C)C(=O)N[C@]([H])(Cc1ccc(O)cc1)C(=O)N[C@]("
        + "[H])([C@@]([H])(CC)C)C(=O)N1[C@]([H])(CCC1)C(=O)N[C@]([H])(C"
        + ")C(=O)N[C@]([H])(CO)C(=O)N[C@]([H])(CC(=O)O)C(=O)N[C@]([H])("
        + "Cc1ccccc1)C(=O)NCC(=O)N[C@]([H])(CC1=CN=C-N1)C(=O)N[C@]([H])"
        + "(CCCCN)C(=O)N[C@]([H])(CC(C)C)C(=O)N[C@]([H])(CS)C(=O)N[C@]("
        + "[H])(C(C)C)C(=O)N[C@]([H])(CC(=O)N)C(=O)N[C@]([H])(CCSC)C(=O)O",
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
        "N[C@@]([H])(C)C(=O)N[C@@]([H])(CC(C)C)C(=O)N[C@@]([H])(C)C(="
        + "O)N[C@@]([H])(CCCN)C(=O)O",
    ),
    (
        "H~ADRPE~OH",
        "C[C@H](N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N"
        + "1CCC[C@@H]1C(=O)N[C@@H](CCC(=O)O)C(=O)O",
    ),
)


correct_peptide_json = {
    "sequence": "H{aMeAla}QGTYCDAQCYS",
    "length": 13,
    "internal_modifications": [],
    "C_terminus": "NH2",
    "N_terminus": "H",
    "pepseq_format": "H~H{aMeAla}QGTY{Cys(R1)}DAQ{Cys(R2)}YS~NH2",
    "symbols": ["H", "H", "aMeAla", "Q", "G", "T", "Y", "Cys(R1)", "D", "A", "Q", "Cys(R2)", "Y", "S", "NH2"],
    "external_modifications": [
        {
            "smiles": "[1*]CC(=O)NCC[C@H](NC(=O)C[2*])C(=O)NCCC(=O)NC"
            + "COC(=O)NCC[C@H](NC(=O)CCC(=O)O)C(=O)O",
            "max_attachment_point_id": 2,
            "attachment_points_on_sequence": {
                1: {
                    "attachment_point_id": "1",
                    "ResID": "7",
                    "AtomName": "SG",
                    "ResidueName": "Cys",
                },
                2: {
                    "attachment_point_id": "2",
                    "ResID": "11",
                    "AtomName": "SG",
                    "ResidueName": "Cys",
                },
            },
        }
    ],
}


fixture_pepseq = "H~H{aMeAla}QGTY{Cys(R1)}DAQ{Cys(R2)}YS~NH2"
one_mod_smiles = (
    "[1*]CC(=O)NCC[C@H](NC(=O)C[2*])C(=O)NCCC(=O)NCCOC(="
    + "O)NCC[C@H](NC(=O)CCC(=O)O)C(=O)O"
)
correct_smiles = (
    "[H]N[C@@H](Cc1c[nH]cn1)C(=O)NC(C)(C)C(=O)N[C@@H](CC"
    + "C(N)=O)C(=O)NCC(=O)N[C@H](C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@H]1"
    + "CSCC(=O)NCC[C@@H](C(=O)NCCC(=O)NCCOC(=O)NCC[C@H](NC(=O)CCC(=O)O)"
    + "C(=O)O)NC(=O)CSC[C@@H](C(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)N[C@@H](CO"
    + ")C(N)=O)NC(=O)[C@H](CCC(N)=O)NC(=O)[C@H](C)NC(=O)[C@H](CC(=O)O)N"
    + "C1=O)[C@@H](C)O"
)


def mols_are_identical(mol1, mol2):
    are_identical = mol1.HasSubstructMatch(
        mol2, useChirality=True
    ) and mol2.HasSubstructMatch(mol1, useChirality=True)
    return are_identical


def smiles_are_identical(smi1, smi2):
    mol1 = rdkit.Chem.MolFromSmiles(smi1)
    mol2 = rdkit.Chem.MolFromSmiles(smi2)
    return mols_are_identical(mol1, mol2)


def test_calculate():
    pepseq_value = "{Cys(R1)}ACDAPEPsEQ{Cys(R2)}"
    smiles = ["[1*]CNCC[2*]"]

    r1 = calculate(pepseq_value, smiles)
    r2 = calculate("CACDAPEPsEQC", None)

    complete_smiles = (
        "[H]N[C@H]1CSCNCCSC[C@@H](C(=O)O)NC(=O)[C@H](CC"
        + "C(N)=O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](CO)NC(=O)[C@@H]2CCC"
        + "N2C(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](C)NC(=O)["
        + "C@H](CC(=O)O)NC(=O)[C@H](CS)NC(=O)[C@H](C)NC1=O"
    )

    assert r1 == {
        "complete_smiles": complete_smiles,
        "length": 12,
        "mw": 1307.45,
        "sequence": "CACDAPEPsEQC",
    }

    complete_smiles = (
        "[H]N[C@@H](CS)C(=O)N[C@@H](C)C(=O)N[C@@H](CS)C"
        + "(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](C)C(=O)N1CCC[C@H]1C(=O)N[C@@H]"
        + "(CCC(=O)O)C(=O)N1CCC[C@H]1C(=O)N[C@H](CO)C(=O)N[C@@H](CCC(=O)O)"
        + "C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CS)C(=O)O"
    )

    assert r2 == {
        "complete_smiles": complete_smiles,
        "length": 12,
        "mw": 1252.37,
        "sequence": "CACDAPEPsEQC",
    }
    pepseq_value = "{Cys(R1)}ACDAPEPsEQ{Cys(R2)}AK{Cys(R3)}"
    smiles = ["[1*]CNCC[2*]", "[3*]CNCC"]
    r3 = calculate(pepseq_value, smiles)

    assert r3 == {
        "complete_smiles": "[H]N[C@H]1CSCNCCSC[C@@H](C(=O)N[C@@H](C)C"
        + "(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CSCNCC)C(=O)O)NC(=O)[C@H](CCC"
        + "(N)=O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](CO)NC(=O)[C@@H]2CCCN"
        + "2C(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](C)NC(=O)[C"
        + "@H](CC(=O)O)NC(=O)[C@H](CS)NC(=O)[C@H](C)NC1=O",
        "length": 15,
        "mw": 1666.95,
        "sequence": "CACDAPEPsEQCAKC",
    }
    return


@pytest.mark.parametrize(
    "pepseq,smiles", pepseq_vs_canonical_smiles + s1l_vs_canonical_smiles
)
def test_peptide_from_pepseq(pepseq, smiles):
    peptide = from_pepseq(pepseq)
    assert smiles_are_identical(peptide.smiles, smiles)


def test_from_pepseq_and_one_mod_smiles_strings_to_peptide_json():

    peptide_json = get_pep_json(fixture_pepseq, db_json, [one_mod_smiles])

    assert peptide_json == correct_peptide_json
    return


def test_from_pepseq_string_and_mod_smiles_to_smiles():
    peptide_json = get_pep_json(fixture_pepseq, db_json, [one_mod_smiles])
    smiles = get_smiles_from_peptide_json(peptide_json, db_json)
    assert smiles == correct_smiles
    return


def test_from_pepseq_string_and_mod_smiles_to_peptide():
    peptide_json = get_pep_json(fixture_pepseq, db_json, [one_mod_smiles])
    peptide = from_json(peptide_json)
    peptide.sequence == correct_peptide_json["sequence"]
    peptide.length == 13
    peptide.complete_smiles == correct_smiles
    return


def test_from_smiles_to_pepseq_and_one_mod_smiles_strings():

    pepseq_format, mod_smiles = from_smiles_to_pepseq_and_one_mod_smiles_strings(
        correct_smiles, db_json
    )

    assert pepseq_format == fixture_pepseq
    assert mod_smiles == one_mod_smiles

    return
