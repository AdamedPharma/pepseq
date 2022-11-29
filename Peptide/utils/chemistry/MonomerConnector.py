import copy
from collections.abc import Sequence
from typing import TypeVar

import rdkit
import rdkit.Chem
from Peptide.utils.chemistry.MolExtended import MolObjExt

AminoAcidInstance = TypeVar("AminoAcidInstance")


def attach(mol1: rdkit.Chem.Mol, mol2: rdkit.Chem.Mol) -> rdkit.Chem.Mol:
    mol1_ext = MolObjExt(mol1)
    mol1_ext.RemoveAtomByLabel(label="_R2")

    mol2_ext = MolObjExt(mol2)
    mol2_ext.RemoveAtomByLabel(label="_R1")

    num_co_mol2 = len(mol2_ext.get_label_ids("_CO"))

    combo = rdkit.Chem.CombineMols(mol1_ext.mol, mol2_ext.mol)
    c = MolObjExt(combo)

    ind_co = -1 - num_co_mol2

    co_id = c.GetAtomIdsByLabel("_CO")[ind_co]
    n_id = c.GetAtomIdsByLabel("_N")[-1]

    edcombo = rdkit.Chem.EditableMol(combo)
    edcombo.AddBond(co_id, n_id, order=rdkit.Chem.rdchem.BondType.SINGLE)
    back = edcombo.GetMol()
    return back


class MonomerConnector(object):
    def __init__(
        self,
    ):
        return

    def connect(
        self, amino_acid_instances: Sequence[AminoAcidInstance] = None
    ) -> rdkit.Chem.Mol:
        pep_mol = amino_acid_instances[0].Mol

        for i in range(1, len(amino_acid_instances)):
            pep_mol = attach(pep_mol, amino_acid_instances[i].Mol)

        return pep_mol
