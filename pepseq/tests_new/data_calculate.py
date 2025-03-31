from data_canonical_sequence import tests

smiles_list = ["[*:1]CNCC[*:2]"]
data_canonical_sequence_value = tests[0]

pepseq_value = "{Cys(R1)}ACDAPEPsEQ{Cys(R2)}AK{Cys(R3)}"

smiles = ["[*:1]CNCC[*:2]", "[*:3]CNCC"]

complete_smiles = (
    "[H]N[C@H]1CSCNCCSC[C@@H](C(=O)O)NC(=O)[C@H](CC"
    "C(N)=O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](CO)NC(=O)[C@@H]2CCC"
    "N2C(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](C)NC(=O)["
    "C@H](CC(=O)O)NC(=O)[C@H](CS)NC(=O)[C@H](C)NC1=O"
)

complete_smiles_no_staple = (
        "[H]N[C@@H](CS)C(=O)N[C@@H](C)C(=O)N[C@@H](CS)C(=O)N[C@@H](CC(=O)O)C(=O)"
        "N[C@@H](C)C(=O)N1CCC[C@H]1C(=O)N[C@@H](CCC(=O)O)C(=O)N1CCC[C@H]1C(=O)N"
        "[C@H](CO)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CS)"
        "C(=O)O"
)


r1_val = {
    "complete_smiles": complete_smiles,
    "length": 12,
    "mw": 1307.45,
    "sequence": "CACDAPEPsEQC",
}

no_smiles_result = {
    "complete_smiles": complete_smiles_no_staple,
    "length": 12,
    "mw": 1252.37,
    "sequence": "CACDAPEPsEQC",
}

r3_val = {
    "complete_smiles": ("[H]N[C@H]1CSCNCCSC[C@@H](C(=O)N[C@@H](C)C"
    "(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CSCNCC)C(=O)O)NC(=O)[C@H](CCC"
    "(N)=O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@@H](CO)NC(=O)[C@@H]2CCCN"
    "2C(=O)[C@H](CCC(=O)O)NC(=O)[C@@H]2CCCN2C(=O)[C@H](C)NC(=O)[C"
    "@H](CC(=O)O)NC(=O)[C@H](CS)NC(=O)[C@H](C)NC1=O"),
    "length": 15,
    "mw": 1666.95,
    "sequence": "CACDAPEPsEQCAKC",
}


tests = [
    ((data_canonical_sequence_value, smiles_list), r1_val),
    (("CACDAPEPsEQC", None), no_smiles_result),
    (("CACDAPEPsEQC", []), no_smiles_result),
    (("CACDAPEPsEQC",), no_smiles_result),
    ((pepseq_value, smiles), r3_val),
]
