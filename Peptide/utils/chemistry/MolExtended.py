import copy
from collections.abc import Sequence

import rdkit
import rdkit.Chem


def get_mol_obj_subst(cls, smiles):
    mol = rdkit.Chem.MolFromSmiles(smiles)
    mol_obj = cls(mol)
    idx = mol_obj.GetAtomIdxByLabel("_R")

    r_atm = mol_obj.mol.GetAtomWithIdx(idx)
    binding_atm = r_atm.GetNeighbors()[0]
    binding_atm_idx = binding_atm.GetIdx()

    if binding_atm.HasProp("atomLabel"):
        displayLabel = binding_atm.GetProp("atomLabel")
        mol_obj.assign_displayLabel_to_idx(binding_atm_idx, displayLabel)

    mol_obj.assign_label_to_idx(binding_atm_idx, "_peptide_binding_point")

    mol_obj.RemoveAtomByLabel("_R")
    return mol_obj.mol


# TO DO: read n-term and c-term smiles from db, add another possible nterm smiles


class MolObjExt(object):
    def __init__(self, mol: rdkit.Chem.Mol = None, sequence: str = None):
        self.mol = mol
        self.sequence = sequence
        return

    @property
    def display_dict(self):
        return self._display_dict

    @display_dict.setter
    def display_dict(self, value):
        self._display_dict = value

        for key in self._display_dict:
            atm = self.GetAtomByLabel(key)
            atm.SetProp("_displayLabel", self._display_dict[key])
        return

    def assign_displayLabel_to_idx(self, idx: int, label: str):
        atm = self.mol.GetAtomWithIdx(idx)
        atm.SetProp("_displayLabel", label)
        return

    def RemoveAtomByLabel(
        self, label: str = "_R2", inplace: bool = True
    ) -> rdkit.Chem.Mol:

        copy_mol = copy.deepcopy(self.mol)
        ed_copy = rdkit.Chem.EditableMol(copy_mol)

        idx = self.GetAtomIdxByLabel(label)

        ed_copy.RemoveAtom(idx)
        back = ed_copy.GetMol()

        if inplace:
            self.mol = back
        else:
            return back

    @property
    def n_term(self):
        return self._n_term

    @n_term.setter
    def n_term(self, value):
        name, smiles = value
        self._n_term = name

        if name == "H":

            radical_atom = self.GetAtomByLabel("_R1")

            protonation_site = radical_atom.GetNeighbors()[0]

            self.protonate(protonation_site)
            self.RemoveAtomByLabel("_R1")

        else:
            self.substitute(smiles=smiles, terminus="N", inplace=True)
        return

    def protonate(self, atm: rdkit.Chem.Atom):
        """
        protonates given atom
        """
        nH = atm.GetNumExplicitHs()
        atm.SetNumExplicitHs(nH + 1)
        return

    def assign_label_to_idx(self, idx, label):
        atm = self.mol.GetAtomWithIdx(idx)
        atm.SetProp("atomLabel", label)
        return

    def substitute(self, smiles="N(*) |$_N;_R$|", inplace=True, terminus="C"):
        mol_subst = get_mol_obj_subst(self.__class__, smiles)

        if terminus == "C":
            subst_point_label = "_R2"
        elif terminus == "N":
            subst_point_label = "_R1"

        pep_r_idx = self.GetAtomIdxByLabel(subst_point_label)

        pep_r_atm = self.mol.GetAtomWithIdx(pep_r_idx)
        pep_r_binding_atm = pep_r_atm.GetNeighbors()[0]
        pep_r_binding_atm_idx = pep_r_binding_atm.GetIdx()

        if pep_r_binding_atm.HasProp("atomLabel"):
            displayLabel = pep_r_binding_atm.GetProp("atomLabel")
            self.assign_displayLabel_to_idx(pep_r_binding_atm_idx, displayLabel)

        self.assign_label_to_idx(pep_r_binding_atm_idx, "_subst_binding_point")

        self.RemoveAtomByLabel(subst_point_label)

        combo = rdkit.Chem.CombineMols(self.mol, mol_subst)
        mol_obj_combo = self.__class__(combo)

        subst_radical_id = mol_obj_combo.GetAtomIdxByLabel("_peptide_binding_point")
        pep_radical_id = mol_obj_combo.GetAtomIdxByLabel("_subst_binding_point")
        edcombo = rdkit.Chem.EditableMol(mol_obj_combo.mol)

        edcombo.AddBond(
            pep_radical_id, subst_radical_id, order=rdkit.Chem.rdchem.BondType.SINGLE
        )
        mol = edcombo.GetMol()
        if inplace:
            self.mol = mol
            for label in ["_subst_binding_point", "_peptide_binding_point"]:  # , ]:
                idx = self.GetAtomIdxByLabel(label)
                atm = self.mol.GetAtomWithIdx(idx)
                displayLabel = atm.GetProp("_displayLabel")
                self.assign_label_to_idx(idx, displayLabel)

        else:
            return mol

    def hydroxylate(self, label, inplace=True):
        mol_hydroxyl = rdkit.Chem.MolFromSmiles("[OH](*) |$_OH;_R$|")
        mol_obj_oh = self.__class__(mol_hydroxyl)

        mol_obj_oh.display_dict = {"_OH": "O"}
        mol_obj_oh.RemoveAtomByLabel("_R")

        mol_hydroxyl = mol_obj_oh.mol

        ids = self.GetAtomIdsByLabel("_CO")
        idx = ids[-1]

        self.assign_label_to_idx(idx, "_cterm_binding_point")
        self.assign_displayLabel_to_idx(idx, "CO")

        ids = self.GetAtomIdsByLabel("_CO")

        combo = rdkit.Chem.CombineMols(self.mol, mol_hydroxyl)
        mol_obj_combo = self.__class__(combo)

        hydroxyl_radical_id = mol_obj_combo.GetAtomIdxByLabel("_OH")
        cter_id = mol_obj_combo.GetAtomIdxByLabel("_cterm_binding_point")

        edcombo = rdkit.Chem.EditableMol(mol_obj_combo.mol)
        edcombo.AddBond(
            cter_id, hydroxyl_radical_id, order=rdkit.Chem.rdchem.BondType.SINGLE
        )
        mol = edcombo.GetMol()
        if inplace:
            self.mol = mol
        else:
            return mol

    @property
    def c_term(self):
        return self._c_term

    @c_term.setter
    def c_term(self, value):
        name, smiles = value
        self._c_term = name

        if value == "OH":
            self.RemoveAtomByLabel("_R2")
            self.hydroxylate("_CO")  # _CO_R2
        else:
            self.substitute(smiles=smiles, terminus="C", inplace=True)
        return

    def GetAtomIdxByLabel(self, label: str) -> int:
        for atm in self.mol.GetAtoms():
            if atm.HasProp("atomLabel") and atm.GetProp("atomLabel") == label:
                idx = atm.GetIdx()
                return idx

    def GetAtomIdsByLabel(self, label: str) -> Sequence[rdkit.Chem.Atom]:
        ids = []
        for atm in self.mol.GetAtoms():
            if atm.HasProp("atomLabel") and atm.GetProp("atomLabel") == label:
                idx = atm.GetIdx()
                ids.append(idx)
        return ids

    def GetAtomByLabel(self, label):
        for atm in self.mol.GetAtoms():
            if atm.HasProp("atomLabel") and atm.GetProp("atomLabel") == label:
                return atm

    def GetAtomsByLabel(self, label):
        ids = self.GetAtomIdsByLabel(label)
        return [self.mol.GetAtomWithIdx(idx) for idx in ids]

    def get_label_ids(self, label: str = "_CA") -> (rdkit.Chem.Atom, int):
        l = []
        for atm in self.mol.GetAtoms():
            if atm.HasProp("atomLabel") and atm.GetProp("atomLabel") == label:
                atm_id = atm.GetIdx()
                l.append((atm, atm_id))
        return l

    def get_CA_idx(self) -> (rdkit.Chem.Atom, int):
        return self.get_label_idx(label="_CA")

    def get_R1_idx(self) -> (rdkit.Chem.Atom, int):
        return self.get_label_idx(label="_R1")

    def get_R2_idx(self) -> (rdkit.Chem.Atom, int):
        return self.get_label_idx(label="_R2")

    def get_CO_idx(self) -> int:
        atm, atm_id = self.get_R2_idx()
        bond = atm.GetBonds()[0]
        CO = (set([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]) - set([atm_id])).pop()
        return CO

    def get_NH_idx(self) -> int:
        atm, atm_id = self.get_R1_idx()
        bond = atm.GetBonds()[0]
        NH = (set([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]) - set([atm_id])).pop()
        return NH

    def del_lab(self, label: str = "_R2", inplace: bool = True) -> rdkit.Chem.Mol:
        r2, r2_id = self.get_label_idx(label)
        ed = rdkit.Chem.EditableMol(copy.deepcopy(self.mol))
        ed.RemoveAtom(r2_id)
        back = ed.GetMol()
        if inplace:
            self.mol = back
        else:
            return back

    def ca_subsitution(
        self, substitute: rdkit.Chem.Mol = None, inplace: bool = True
    ) -> rdkit.Chem.Mol:
        aa = self.mol
        n_atoms = len(aa.GetAtoms())
        ch3 = substitute

        aa_ch3_combo = rdkit.Chem.CombineMols(aa, ch3)
        edcombo = rdkit.Chem.EditableMol(aa_ch3_combo)
        CA_idx = self.get_CA_idx()[1]

        Ch3_atm_id = (n_atoms - 1) + 1

        edcombo.AddBond(CA_idx, n_atoms, order=rdkit.Chem.rdchem.BondType.SINGLE)
        back = edcombo.GetMol()
        if inplace:
            self.mol = back
        else:
            return back

    def alfa_methylo_aa(self, inplace: bool = True) -> rdkit.Chem.Mol:
        ch3 = rdkit.Chem.MolFromSmiles("C")
        return self.ca_subsitution(substitute=ch3, inplace=inplace)
