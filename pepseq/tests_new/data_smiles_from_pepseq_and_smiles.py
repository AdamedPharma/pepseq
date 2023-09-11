import pkgutil
import json


db_path = pkgutil.extend_path("pepseq/Peptide/database/db.json", __name__)
with open(db_path) as fp:
    db_json = json.load(fp)


fixture_pepseq = "H~H{aMeAla}QGTY{Cys(R1)}DAQ{Cys(R2)}YS~NH2"
fixture_pepseq_2 = "CH3~Y{Gly(R1)}QGTFTSDYSKYLDECAAKDFVCWLLDHHPSSGQPPPS~NH2"
fixture_pepseq_3 = "CH3~Y{ala(R1)}QGTFTSDYSKYLDECAAKDFVCWLLDHHPSSGQPPPS~NH2"

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


correct_smiles_2 = (
    "CCCCCSC(NC(=O)[C@H](Cc1ccc(O)cc1)NC)C(=O)N[C@@H](CCC"
    + "(N)=O)C(=O)NCC(=O)N[C@H](C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@H](C(=O)"
    + "N[C@@H](CO)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N["
    + "C@@H](CO)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H"
    + "](CC(C)C)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H]("
    + "CS)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H]("
    + "CC(=O)O)C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@H](C(=O)N[C@@H](CS)C(=O)N"
    + "[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CC(C)C)"
    + "C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N[C@@H](Cc1c[n"
    + "H]cn1)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CO)C(=O)N[C@@H](CO)C(=O)NCC(=O"
    + ")N[C@@H](CCC(N)=O)C(=O)N1CCC[C@H]1C(=O)N1CCC[C@H]1C(=O)N1CCC[C@H]"
    + "1C(=O)N[C@@H](CO)C(N)=O)C(C)C)[C@@H](C)O)[C@@H](C)O"
)


correct_smiles_3 = (
    "CCCCCSC[C@@H](NC(=O)[C@H](Cc1ccc(O)cc1)NC)C(=O)N[C@@H](CCC(N)=O)C(="
    + "O)NCC(=O)N[C@H](C(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@H](C(=O)N[C@@H](C"
    + "O)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](CO)"
    + "C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](CC(C)C)"
    + "C(=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CS)C(=O)N"
    + "[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(=O)O)C"
    + "(=O)N[C@@H](Cc1ccccc1)C(=O)N[C@H](C(=O)N[C@@H](CS)C(=O)N[C@@H](Cc"
    + "1c[nH]c2ccccc12)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@"
    + "@H](CC(=O)O)C(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N[C@@H](Cc1c[nH]cn1)C(="
    + "O)N1CCC[C@H]1C(=O)N[C@@H](CO)C(=O)N[C@@H](CO)C(=O)NCC(=O)N[C@@H]("
    + "CCC(N)=O)C(=O)N1CCC[C@H]1C(=O)N1CCC[C@H]1C(=O)N1CCC[C@H]1C(=O)N[C"
    + "@@H](CO)C(N)=O)C(C)C)[C@@H](C)O)[C@@H](C)O"
)


one_mod_smiles_2 = (
    "[1*]SCCCCC"
)


tests = [
    (
        (fixture_pepseq, db_json, [one_mod_smiles]), correct_smiles
    ),
    (
        (fixture_pepseq_2, db_json, [one_mod_smiles_2]), correct_smiles_2
    ),
    (
        (fixture_pepseq_3, db_json, [one_mod_smiles_2]), correct_smiles_3
    )

]