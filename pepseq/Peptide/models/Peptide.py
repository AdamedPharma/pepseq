from collections import namedtuple
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

Parameters = namedtuple("Parameters", "name smiles")

AminoAcidInstance = TypeVar("AminoAcidInstance")
PeptideWriter = TypeVar("PeptideWriter")


def get_smiles_descriptors(smiles: AnyStr) -> dict:
    """
    Get the molecular descriptors for a given SMILES string.

    :param smiles: The SMILES string.
    :type smiles: str

    :return: The molecular descriptors.
    :rtype: dict
    """
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


class Peptide(object):
    """
    A class to represent a peptide.

    :param smiles: The SMILES string.
    :type   smiles: str

    :param peptide_json: The peptide JSON.
    :type   peptide_json: dict

    """
    def __init__(self, smiles: str, peptide_json: dict):
        """
        Initialize the Peptide class.

        :param smiles: The SMILES string.
        :type  smiles: str

        :param peptide_json: The peptide JSON.
        :type  peptide_json: dict

        :return: None
        :rtype:  None
        """
        self.smiles = smiles
        self.complete_smiles = smiles
        self.peptide_json = peptide_json
        self.sequence = self.peptide_json.get("sequence")
        self.length = self.peptide_json.get("length")
        descriptors = get_smiles_descriptors(self.complete_smiles)
        self.mw = descriptors["mw"]
        return
