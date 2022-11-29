import rdkit
import rdkit.Chem


class Molecule(object):
    def __init__(self, Mol: rdkit.Chem.Mol = None):
        self._Mol = Mol
        return

    @classmethod
    def MolFromSequence(cls, sequence: str = None):
        Mol = rdkit.Chem.MolFromSequence(sequence)
        obj = cls(Mol)
        return obj

    @property
    def Mol(self):
        return self._Mol

    @Mol.setter
    def Mol(self, value: rdkit.Chem.Mol = None):
        self._Mol = value
        if type(value) == rdkit.Chem.Mol:
            self._smiles = rdkit.Chem.MolToSmiles(self._Mol)
        return

    @property
    def smiles(self):
        return self._smiles

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value: str = None):
        self._name = value
        return
