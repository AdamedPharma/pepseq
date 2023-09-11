# content of conftest.py
import pytest
import rdkit
import rdkit.Chem

TERMS_LIST = [
  [[4, 3], [-2, 1]],
  [[-4, 3], [-2, 1]]
]


class Polynomial(object):
    def __init__(self, params):
        #poly = Polynomial([[4, 3], [-2, 1]])
        #exponents = [term[1] for term in poly.terms]
        self.params = params
        #self.terms = [i[0] for i in self.params]
        self.terms = self.params

        return
    
    def derivative(self):
        return self.__class__(self.params)
    


@pytest.fixture(scope='module', params=TERMS_LIST)
def poly(request):
  return Polynomial( request.param )


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
    (
        "[CSPS[1*]]~Y{aMeAla}QGTFTSDYSKYLDECAAKDFVCWLLDHHPSSGQPPPS~[CNSC[1*]]",
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
        + "(C)(C)NC(=O)[C@H](Cc1ccc(O)cc1)NSPSC)[C@@H](C)O)[C@@H](C)O)C(C)C"
    ),
)

# https://remusao.github.io/posts/pytest-param.html
