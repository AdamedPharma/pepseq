from rdkit import Chem

gly_from_smarts = Chem.MolFromSmarts("[C:0](=[O:1])[C:2][N:3]")


def choose_biggest(list_of_potential_amino_acids, aas_ordered_by_size):
    """
    ordered_nodes <- list of atoms ordered from smallest to largest
    list_of_potential_amino_acids <- list_of_potential_amino_acids
    return the biggest potential amino acid
    """
    for aa in reversed(aas_ordered_by_size):
        if aa in list_of_potential_amino_acids:
            return aa


def get_atm_id_to_res(d, aas_ordered_by_size):
    """
    takes dictionary d of form:
            d = {
                atom_id : [ (res_id1, res_name1), (res_id2, res_name2) ],
                ...,
            }

    chooses the biggest residue out of potential residues,

    returns:
        atm_id_to_res = {
            atom_id : (res_id, res_name)
        }
    """

    atm_id_to_res = {}
    for atom_id in sorted(d.keys()):
        potential_residues = d[atom_id]  # takes potential residues
        list_of_potential_amino_acids = [i[1] for i in potential_residues]
        biggest_residue = choose_biggest(
            list_of_potential_amino_acids, aas_ordered_by_size
        )
        residue = [i for i in potential_residues if i[1] == biggest_residue][0]
        atm_id_to_res[atom_id] = residue
    return atm_id_to_res


def generate_bb_smiles(num_res, OH=False):
    if OH:
        bb_smiles = "O" + "C(=O)CN" * num_res  # generate backbone SMILES
    else:
        bb_smiles = "C(=O)CN" * num_res  # generate backbone SMILES

    return bb_smiles


def generate_bb_mol(num_res, OH=False):
    bbsmiles = generate_bb_smiles(num_res, OH=OH)
    bbmol = Chem.MolFromSmiles(bbsmiles)
    return bbmol


class Smi2Seq(object):
    # switch to stages problem is in between somehow
    # check peptide if the SMILES is alright
    # if it is really aMeAla on middle position (if in the middle)

    def __init__(self, molecule=None):
        self.m = molecule
        return

    @property
    def gly_matches(self):
        """
        find atom_ids for matched glycine residues
        """
        matches = self.m.GetSubstructMatches(gly_from_smarts)
        return matches

    @property
    def CAatoms(self):
        """
        get C alpha atoms instances (rdkit.Chem.rdchem.Atom)
        """

        if self.__dict__.get("_CAatoms") is None:
            ca_atoms = self.get_CAatoms(longest_bb=True)
            self._CAatoms = ca_atoms
        return self._CAatoms

    def get_CAatoms(self, longest_bb=False):
        """
        get C alpha atoms instances (rdkit.Chem.rdchem.Atom)
        """
        atoms = []
        ca_ids = []

        if longest_bb:
            num_residues = int(len(self.longest_backbone) / 4)
            for i in range(num_residues):
                ca_id = (i * 4) + 1
                ca_ids.append(ca_id)
        else:
            for gly_match in self.gly_matches:
                ca_id = gly_match[2]
                ca_ids.append(ca_id)

        for ca_id in ca_ids:
            atom = self.m.GetAtomWithIdx(ca_id)
            atoms.append(atom)
        return atoms

    def label_CAatoms(self, longest_bb=True):
        """
        sets name for CAatoms as ' CA '
        """
        for CAatom in self.get_CAatoms(longest_bb):
            info = Chem.AtomPDBResidueInfo()
            info.SetName(" CA ")  # spaces are important
            CAatom.SetMonomerInfo(info)
        return

    @property
    def backbone(self):
        num_res = len(self.CAatoms)
        bbsmiles = generate_bb_smiles(num_res, OH=False)
        bbmol = Chem.MolFromSmiles(bbsmiles)
        bb = self.m.GetSubstructMatches(bbmol)[0]
        return bb

    @property
    def longest_backbone(self):
        """
        iteratively
            generates rdkit.Chem.Molecule(s)
            representing glycine chain (protein backbone)
            of increasing length (incrementing by 1 every iteration)
            and matches it to the self.Molecule
            when after n iterations generated GlycineChain is no longer
            a substructure of self.Molecule
            iteration loop breaks
            and  match of glycine chain generated in previous iteration
            to molecule is returned
        """

        matches = [None]
        num_res = 0
        while matches:
            bb = matches[0]
            num_res += 1
            bbmol = generate_bb_mol(num_res, OH=False)
            matches = self.m.GetSubstructMatches(bbmol)
        return bb

    def find_matches(self, aa_mol_code=None, code="SMILES"):
        if code == "SMILES":
            aa_mol = Chem.MolFromSmiles(aa_mol_code)
        elif code == "SMARTS":
            aa_mol = Chem.MolFromSmarts(aa_mol_code)
        matches = self.m.GetSubstructMatches(aa_mol, useChirality=True)
        ordered_matches = []
        for match in matches:
            match = tuple(sorted(list(match)))
            ordered_matches.append(match)
        return ordered_matches

    @property
    def bb_list(self):
        id_list = list(self.backbone)
        id_list.reverse()
        return id_list

    @property
    def longest_bb_list(self):
        id_list = list(self.longest_backbone)
        id_list.reverse()
        return id_list

    def renumber(self, longest_bb=False):

        if longest_bb:
            id_list = self.longest_bb_list

        else:
            id_list = self.bb_list

        atom_ids = [a.GetIdx() for a in self.m.GetAtoms()]

        for idx in atom_ids:
            if idx not in id_list:
                id_list.append(idx)
        self.m = Chem.RenumberAtoms(self.m, newOrder=id_list)

        return

    def aa_matches(self, aa_smiles_dict, code="SMILES"):
        """
        For each amino_acid defined in database along with its SMILES
        finds substructures of said amino_acid in molecule

        returns dict in form:
        {
           'amino_acid1': (
                tuple_of_amino_acid1_substructure1_in_molecule_atom_ids,
                tuple_of_amino_acid1_substructure2_in_molecule_atom_ids,
                ),
            'amino_acid2': (
                tuple_of_amino_acid2_substructure1_in_molecule_atom_ids
            )
        }
        e.g.
        {
            'A': (
                (38, 25, 24, 27, 26),
                (40, 21, 20, 23, 22),
                (47, 17, 16, 19, 18),
                (52, 13, 12, 15, 14),
                (60, 9, 8, 11, 10),
                (61, 5, 4, 7, 6),
                (62, 1, 0, 3, 2),
                (67, 29, 28, 31, 30),
                (68, 33, 32, 35, 34)
                ),
            'E': (
                (16, 17, 19, 18, 47, 48, 49, 51, 50),
                ),
            'F': (
                (20, 21, 23, 22, 40, 41, 42, 43, 44, 45, 46),
                ),
            'G': (
                (0, 1, 3, 2),
                (4, 5, 7, 6),
                (8, 9, 11, 10),
                (12, 13, 15, 14),
                (16, 17, 19, 18),
                (20, 21, 23, 22),
                (24, 25, 27, 26),
                (28, 29, 31, 30),
                (32, 33, 35, 34)
                ),
                'I': (
                    (36, 37, 38, 39, 25, 24, 27, 26),
                    ),
                'K': (
                    (56, 55, 54, 53, 52, 13, 12, 15, 14),
                    ),
                'Q': (
                    (65, 64, 66, 63, 62, 1, 0, 3, 2),
                    ),
                'V': (
                    (37, 38, 39, 25, 24, 27, 26),
                    ),
                'W': (
                    (32, 33, 35, 34, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77),
                    )
                }

        """
        d_matches = {}
        aa_codes = sorted(list(aa_smiles_dict.keys()))

        for aa_code in aa_codes:
            aa = aa_smiles_dict.get(aa_code)
            if code == "SMILES":
                aa_mol_code = aa.smiles
            elif code == "SMARTS":
                aa_mol_code = aa.smarts

            matches = self.find_matches(aa_mol_code=aa_mol_code, code=code)
            if matches:
                d_matches[aa_code] = matches
        return d_matches
