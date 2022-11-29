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


class AminoAcidInstance(Molecule):
    @classmethod
    def MolFromSmiles(cls, smiles: str):
        Mol = rdkit.Chem.MolFromSmiles(smiles)
        obj = cls()
        obj._Mol = Mol
        obj._smiles = smiles
        return obj

    @classmethod
    def MolFromSymbol(cls, symbol: str = None, db_api: DataBase = None):
        valid_symbols_set = set(db_api.symbols.keys()) | set(
            db_api.aa_smiles_dict.keys()
        )

        if symbol not in valid_symbols_set:
            raise InvalidSymbolError("Invalid Symbol: %s" % symbol)

        amino_acid_db = db_api.aa_smiles_dict.get(symbol)
        smiles_radical = amino_acid_db.smiles_radical

        obj = cls.MolFromSmiles(smiles_radical)
        obj.symbol = symbol
        return obj

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
