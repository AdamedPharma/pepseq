import json
import pkgutil

import pytest
import rdkit
from BuildingModifiedPeptideFromPeptideJSON import get_smiles_from_peptide_json
from BuildPeptideJSONFromSMILES import decompose_peptide_smiles_with_termini
from get_peptide_json_from_pepseq_format import get_pep_json
from read import from_pepseq

db_path = pkgutil.extend_path("Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


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
    (
        "H~ADRPE~OH",
        "C[C@H](N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N1CCC[C@@H]1C(=O)N[C@@H](CCC(=O)O)C(=O)O",
    ),
)


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
    assert smiles_are_identical(peptide.smiles, smiles)


# or


def test_peptide_json_to_smiles():
    # this tes tests whether just based on
    return


def test_peptide_json_sequence_with_non_standard_aas_to_smiles():
    # this tes tests whether just based on
    # on PeptideJSON having non-standard amino acids we can construct
    # get peptide smiles
    return


def test_from_pepseq_and_one_mod_smiles_strings_to_peptide_json():

    """

    Input:


        pepseq_string:

            str = string in pepseq format H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
        where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified amino acid; {Cys(R1)} - is amino acid
        with staple attached, {Cys(R1)} - amino acid with staple attached


        mod_smiles:

            SMILES string (e.g. '[1*]C[2*]') - showing the structure of modification with attachment
            points:

                { Cys(R1) } <- is attached in [1*] attachment point on staple
                { Cys(R2) } <- is attached in [2*] attachment point on staple

    Output:

        peptide_json:

            JSON containing info about modified peptide with

                'sequence':

                'internal_modifications':

                'external_modifications':

    """
    pepseq_string, one_mod_smiles = "H~H{Cys(R1)}I{Cys(R2)}K~OH", "[1*]CS[2*]"
    peptide_json = get_pep_json(pepseq_string, one_mod_smiles)
    correct_peptide_json = {"placeholder": "placeholder_val"}
    assert peptide_json == correct_peptide_json


def test_from_pepseq_string_and_mod_smiles_to_smiles(pepseq_string, one_mod_smiles):
    """
    Input:


    """
    pepseq_string, one_mod_smiles = "H~H{Cys(R1)}I{Cys(R2)}K~OH", "[1*]CS[2*]"
    peptide_json = get_pep_json(pepseq_string, one_mod_smiles)
    smiles = get_smiles_from_peptide_json(peptide_json, db_json)
    correct_smiles = '"SCP |$R1;placeholder;R2$|"'
    assert smiles == correct_smiles
    return


def test_from_pepseq_and_many_mod_smiles_strings_to_peptide_json(
    pepseq_string, many_mod_smiles
):
    """

    Input:


        pepseq_string:

            str = string in pepseq format H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
        where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified amino acid; {Cys(R1)} - is amino acid
        with staple attached, {Cys(R1)} - amino acid with staple attached


        mod_smiles:

            SMILES string (e.g. '[1*]C[2*]') - showing the structure of modification with attachment
            points:

                { Cys(R1) } <- is attached in [1*] attachment point on staple
                { Cys(R2) } <- is attached in [2*] attachment point on staple

    Output:

        peptide_json:

            JSON containing info about modified peptide with

                'sequence':

                'internal_modifications':

                'external_modifications':

    """

    return


def test_from_smiles_to_peptide_json():
    """
    Input:

        SMILES - string of peptide sequence with modified amino acids
            modification(s)

    Output:
        peptide_json:

            JSON containing info about modified peptide with

                'sequence':

                'internal_modifications':

                'external_modifications':


        mod_smiles:

            SMILES string (e.g. '[1*]C[2*]') - showing the structure of modification with attachment
            points:

                { Cys(R1) } <- is attached in [1*] attachment point on staple
                { Cys(R2) } <- is attached in [2*] attachment point on staple

    """
    smiles = '"SCP |$R1;placeholder;R2$|"'

    peptide_json = decompose_peptide_smiles_with_termini(smiles, db_json)
    correct_peptide_json = {"placeholder": "placeholder_val"}
    assert peptide_json == correct_peptide_json
    return


def test_peptide_json_to_pepseq_string_and_one_mod_smiles():
    """

    Input:

        peptide_json:

            JSON containing info about modified peptide with

                'sequence':

                'internal_modifications':

                'external_modifications':


    Output:
        pepseq_string:

            str = string in pepseq format H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
        where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified amino acid; {Cys(R1)} - is amino acid
        with staple attached, {Cys(R1)} - amino acid with staple attached

        modifications - external ones, with attachment points

        mod_smiles:

            SMILES string (e.g. '[1*]C[2*]') - showing the structure of modification with attachment
            points:

                { Cys(R1) } <- is attached in [1*] attachment point on staple
                { Cys(R2) } <- is attached in [2*] attachment point on staple

    """
    smiles = '"SCP |$R1;placeholder;R2$|"'

    peptide_json = decompose_peptide_smiles_with_termini(smiles, db_json)
    correct_pepseq_string, one_mod_smiles = "H~H{Cys(R1)}I{Cys(R2)}K~OH", "[1*]CS[2*]"

    assert peptide_json["pepseq_format"] == correct_pepseq_string
    assert peptide_json["external_modifications"][0]["smiles"] == one_mod_smiles
    return


def test_peptide_json_one_mod_to_pepseq_string_many_mod_smiles():
    """
    Input:

        peptide_json:

            JSON containing info about modified peptide with

                'sequence':

                'internal_modifications':

                'external_modifications':


    Output:
        pepseq_string:

            str = string in pepseq format H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
        where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified amino acid; {Cys(R1)} - is amino acid
        with staple attached, {Cys(R1)} - amino acid with staple attached

        mod_smiles:

            SMILES string (e.g. '[1*]C[2*]') - showing the structure of modification with attachment
            points:

                { Cys(R1) } <- is attached in [1*] attachment point on staple
                { Cys(R2) } <- is attached in [2*] attachment point on staple

    """

    return


def test_from_smiles_to_pepseq_and_one_mod_smiles_strings():
    """
    Input:

        SMILES - string of peptide sequence with modified amino acids
            modification(s)

    Output:
        pepseq_string:

            str = string in pepseq format H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
        where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified amino acid; {Cys(R1)} - is amino acid
        with staple attached, {Cys(R1)} - amino acid with staple attached

        modifications - external ones, with attachment points

        mod_smiles:

            SMILES string (e.g. '[1*]C[2*]') - showing the structure of modification with attachment
            points:

                { Cys(R1) } <- is attached in [1*] attachment point on staple
                { Cys(R2) } <- is attached in [2*] attachment point on staple

    """

    return


def test_from_smiles_to_pepseq_and_many_mod_smiles_strings():
    """
    Input:

        SMILES - string of peptide sequence with modified amino acids
            modification(s)

    Output:
        pepseq_string:

            str = string in pepseq format H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
        where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified amino acid; {Cys(R1)} - is amino acid
        with staple attached, {Cys(R1)} - amino acid with staple attached

        modifications - external ones, with attachment points

        mod_smiles:

            SMILES string (e.g. '[1*]C[2*]') - showing the structure of modification with attachment
            points:

                { Cys(R1) } <- is attached in [1*] attachment point on staple
                { Cys(R2) } <- is attached in [2*] attachment point on staple

    """

    return
