
tests = (
        ("data/mypeptide.smi", "data/myppeptide_out", ["H~{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF~OH"], ["[*:1]CNCC[*:2]", "[*:3]CNCCSP"]),
        ('data/test_disulfide.smi', "data/disulfide_out", ['H~AA{Cys(R1)}AA{Cys(R1)}AA~OH'], ['NO (NON TERMINAL) SEQUENCE MODIIFICATIONS PRESENT'])
)
