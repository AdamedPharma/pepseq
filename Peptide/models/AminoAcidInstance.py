from typing import TypeVar

import rdkit
import rdkit.Chem

from Peptide.exceptions import InvalidSymbolError
from Peptide.models.Molecule import Molecule

DataBase = TypeVar("DataBase")


def parse_symbol(symbol: str) -> (str, rdkit.Chem.Mol, str):
    for base_aa_symbol in base_amino_acids:
        if base_aa_symbol in symbol:
            substitute = symbol.replace(base_amino_acid, "")
            if substitute == "aMe":
                substitute_radical = rdkit.Chem.MolFromSmiles("C")
                attachment_atom_label = "CA"
            elif substitute == "aEth":
                substitute_radical = rdkit.Chem.MolFromSmiles("CC")
                attachment_atom_label = "CA"

            return base_aa_symbol, substitute_radical, attachment_atom_label


def validate_symbol(symbol: str, db_api: DataBase):
    if symbol not in db_api.valid_symbols:
        raise InvalidSymbolError("Invalid Symbol: %s" % symbol)


class AminoAcidInstance(Molecule):
    def __init__(self, smiles: str, symbol: str):
        Mol = rdkit.Chem.MolFromSmiles(smiles)
        self._Mol = Mol
        self._smiles = smiles
        self.symbol = symbol
        return

    @property
    def smiles(self):
        return self._smiles

    @smiles.setter
    def smiles(self, value: str = None):
        self._smiles = value
        self._Mol = rdkit.Chem.MolFromSmiles(value)
        return

    @property
    def smiles_with_radical(self):
        return self._smiles_with_radical

    @smiles_with_radical.setter
    def smiles_with_radical(self, value: str = None):
        self._smiles_with_radical = value
        return


def AminoAcidFromSymbol(symbol: str, db_api: DataBase):
    validate_symbol(symbol, db_api)
    smiles_radical = db_api.read_smiles_radical(symbol)
    return AminoAcidInstance(smiles_radical, symbol)
