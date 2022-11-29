from collections import namedtuple
from collections.abc import Sequence
from typing import AnyStr, TypeVar

from rdkit.Chem import MolFromSmiles
from rdkit.Chem.Descriptors import TPSA, ExactMolWt, MolLogP, MolWt
from rdkit.Chem.Lipinski import (
    HeavyAtomCount,
    NumHAcceptors,
    NumHDonors,
    NumHeteroatoms,
    NumRotatableBonds,
    RingCount,
)
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from Peptide.models.Molecule import Molecule
from Peptide.utils.chemistry.MolExtended import MolObjExt
from Peptide.utils.chemistry.MonomerConnector import MonomerConnector

Parameters = namedtuple("Parameters", "name smiles")

AminoAcidInstance = TypeVar("AminoAcidInstance")
PeptideWriter = TypeVar("PeptideWriter")


def get_smiles_descriptors(smiles: AnyStr):
    molecule = MolFromSmiles(smiles)
    descriptors: dict = {
        "mw": round(MolWt(molecule), 2),
        "exact_mw": round(ExactMolWt(molecule), 2),
        "logp": round(MolLogP(molecule), 2),
        "hba": NumHAcceptors(molecule),
        "hbd": NumHDonors(molecule),
        "num_heavy_atoms": HeavyAtomCount(molecule),
        "heteroatoms": NumHeteroatoms(molecule),
        "rotatable_bonds": NumRotatableBonds(molecule),
        "rings": RingCount(molecule),
        "formula": CalcMolFormula(molecule),
        "atoms": molecule.GetNumAtoms(),
        "tpsa": round(TPSA(molecule), 2),
    }

    return descriptors


# notacja z prezentacji z next move software
# jak nie znmy mod to np gwiazdka
# regula 7: implicit leaving group (lys-)
# pisze testy i okreslam jakie jest oczekiwane zachowanie
# i testy na DAMBA
# moge przetestowac na bazie danych (potem)
#


class Peptide(Molecule):
    @classmethod
    def create(
        cls,
        name: str = None,
        smiles: str = None,
        amino_acids: Sequence[AminoAcidInstance] = None,
    ):
        obj = cls()
        if amino_acids is not None:
            obj.amino_acids = amino_acids
        if name is not None:
            self.name = name
        return obj

    @property
    def smiles(self):
        return self._smiles

    @smiles.setter
    def smiles(self, value: str = None):
        self._smiles = value
        return

    @property
    def writer(self):
        return self._writer

    @writer.setter
    def writer(self, value: PeptideWriter):
        self._writer = value
        return

    @property
    def sequence(self):
        if self.__dict__.get("_sequence") is None:
            if not self.amino_acids:
                return None

            core_symbols_list = []
            for i in self.amino_acids:
                if len(i.symbol) > 1:
                    symbol = "{%s}" % i.symbol
                else:
                    symbol = i.symbol
                core_symbols_list.append(symbol)
            core_seq = "".join(core_symbols_list)
            sequence = "%s~%s~%s" % (self.n_term.name, core_seq, self.c_term.name)
            self._sequence = sequence
        return self._sequence

    @property
    def canonical_sequence(self):
        if self.__dict__.get("_canonical_sequence") is None:
            core_seq = "".join([i.symbol for i in self.amino_acids])
            sequence = "%s~%s~%s" % (self.n_term.name, core_seq, self.c_term.name)
            self._canonical_sequence = sequence
        return self._canonical_sequence

    @sequence.setter
    def sequence(self, value: str):
        self._sequence = value
        return

    @property
    def amino_acids(self):
        if self.__dict__.get("_amino_acids") is None:
            self._amino_acids = []
        return self._amino_acids

    @amino_acids.setter
    def amino_acids(self, amino_acids_value: Sequence[AminoAcidInstance]):
        self._amino_acids = amino_acids_value
        self.Mol = self.calc_smiles(amino_acids_value)
        self.recalculate_parameters()
        return

    def recalculate_parameters(self):
        self._parameters = get_smiles_descriptors(self.smiles)
        return

    def calc_smiles(self, amino_acids_value: Sequence[AminoAcidInstance]):
        monomer_connector = MonomerConnector()
        return monomer_connector.connect(amino_acids_value)

    @property
    def n_term(self):
        if self.__dict__.get("_n_term") is None:
            self._n_term = None
        return self._n_term

    @n_term.setter
    def n_term(self, value):
        self._n_term = value
        mol_obj_ext = MolObjExt(self.Mol)
        mol_obj_ext.n_term = value
        self.Mol = mol_obj_ext.mol
        self.recalculate_parameters()
        return

    @property
    def c_term(self):
        if self.__dict__.get("_c_term") is None:
            self._c_term = None
        return self._c_term

    @c_term.setter
    def c_term(self, value):
        self._c_term = value
        mol_obj_ext = MolObjExt(self.Mol)
        mol_obj_ext.c_term = value
        self.Mol = mol_obj_ext.mol
        self.recalculate_parameters()
        return

    @property
    def parameters(self):
        if self.__dict__.get("_parameters") is None:
            self.recalculate_parameters()
        return self._parameters
