import rdkit.Chem
from Peptide.db_api.DataBase import FileSystemDbRepo
from Peptide.utils.chemistry.smi2seq_utils import (
    get_atm_id_to_potential_res_ids_names_dict,
    get_rows,
)
from Peptide.utils.chemistry.Smi2SeqObj import Smi2Seq, get_atm_id_to_res
from Peptide.utils.chemistry.SubstructureGraph import get_ordered_nodes
from Peptide.utils.Parser import Parser
from Peptide.utils.PeptideReader import SeqReader
from rdkit.Chem import AllChem


def get_modified_residues(connected_atom_sets, atm_id_to_res_dict):
    modified_residues = []
    for i in connected_atom_sets:
        res = atm_id_to_res_dict.get(i)
        if res is not None:
            modified_residues.append(res)
    return modified_residues


def find_group_connected(group_atoms):
    """
    finds atoms that connect atoms in group (e.g. peptide)
     to other group (e.g. modification)
    """

    l = []
    group_atoms_ids = set([atom.GetIdx() for atom in group_atoms])

    for atom in group_atoms:

        bonds = atom.GetBonds()

        for bond in bonds:

            atom_id = bond.GetBeginAtomIdx()
            if atom_id not in group_atoms_ids:
                l.append(atom_id)
            else:
                atom_id = bond.GetEndAtomIdx()
                if atom_id not in group_atoms_ids:
                    l.append(atom_id)
    return l


def get_seq(amino_acids, mod_res_nums, mod_as_x=False):
    symbols = []

    for i in range(len(amino_acids)):
        amino_acid = amino_acids[i]
        if i in set([i[0] - 1 for i in mod_res_nums]):
            if mod_as_x:
                symbol = "X"
            else:
                symbol = "{mod%s}" % (amino_acid.symbol)
        else:
            if len(amino_acid.symbol) > 1:
                if mod_as_x:
                    symbol = "X"
                else:
                    symbol = "{%s}" % (amino_acid.symbol)
            else:
                symbol = amino_acid.symbol
        symbols.append(symbol)

    seq = "".join(symbols)
    return seq


class UseCase(object):
    def __init__(self, smiles=None, repo=None, repo_path="./Peptide/database/db.json"):
        self.seq_reader = SeqReader()
        self.parser = Parser()

        self.smiles = smiles

        if self.smiles is None:
            self.molecule = None
        else:
            self.molecule = rdkit.Chem.MolFromSmiles(self.smiles)

        if repo is None:
            repo = FileSystemDbRepo.read_from_json(path=repo_path)
        self.repo = repo
        return

    @property
    def smi2seq_obj(self):
        if self.__dict__.get("_smi2seq_obj") is None:
            self._smi2seq_obj = self.get_smi2seq_obj()
        return self._smi2seq_obj

    def get_smi2seq_obj(self):
        smi2seq_obj = Smi2Seq(molecule=self.molecule)
        smi2seq_obj.renumber(longest_bb=True)
        self.molecule = smi2seq_obj.m
        smi2seq_obj.label_CAatoms()
        return smi2seq_obj

    def deduce_seq(self):
        # USED

        aa_species_assigned_potential_atom_ids = self.smi2seq_obj.aa_matches(
            self.repo.aa_smiles_dict, code="SMARTS"
        )
        ca_atom_id_to_res_id_dict = self.get_ca_atom_id_to_res_id_dict(self.smi2seq_obj)

        rows = get_rows(
            aa_species_assigned_potential_atom_ids, ca_atom_id_to_res_id_dict
        )

        atm_id_to_potential_res_ids_names_dict = (
            get_atm_id_to_potential_res_ids_names_dict(rows)
        )

        aas_ordered_by_size = get_ordered_nodes(self.repo.aa_smiles_dict)

        atm_id_to_res_dict = get_atm_id_to_res(
            atm_id_to_potential_res_ids_names_dict, aas_ordered_by_size
        )
        symbols_list = [i[1] for i in sorted(list(set(atm_id_to_res_dict.values())))]
        for i in range(len(symbols_list)):
            if len(symbols_list[i]) > 1:
                symbols_list[i] = "{%s}" % symbols_list[i]
        deduced_seq = "".join(symbols_list)
        deduced_peptide = self.seq_reader.read(
            deduced_seq, self.parser, db_api=self.repo
        )
        return deduced_peptide, deduced_seq, atm_id_to_res_dict

    def get_ca_atom_id_to_res_id_dict(self, smi2seq_obj):
        backbone_atoms_ids = list(smi2seq_obj.longest_backbone)
        backbone_atoms_ids.reverse()
        backbone_atoms_split_by_res = [
            tuple(range(i, i + 4)) for i in range(0, len(backbone_atoms_ids), 4)
        ]
        ca_atom_id_to_res_id_dict = {
            backbone_atoms_split_by_res[i][2]: i + 1
            for i in range(len(backbone_atoms_split_by_res))
        }
        return ca_atom_id_to_res_id_dict


def find_connecting_residues(mol=None, mod_mol=None, atm_id_to_res_dict=None):
    matches = mol.GetSubstructMatches(mod_mol, useChirality=True)
    match = matches[0]
    group_atoms = [mol.GetAtomWithIdx(i) for i in match]
    connected_atom_sets = find_group_connected(group_atoms)
    mod_res_nums = get_modified_residues(connected_atom_sets, atm_id_to_res_dict)
    return mod_res_nums


def read_seq_from_smiles(smiles=None, repo=None, mod_as_x=False):
    uc = UseCase(smiles=smiles, repo=repo)
    mol = rdkit.Chem.MolFromSmiles(smiles)

    deduced_peptide, deduced_seq, atm_id_to_res_dict = uc.deduce_seq()
    keys = set(atm_id_to_res_dict.keys())
    mol_atoms = [(i, i.GetIdx(), i.GetSymbol()) for i in uc.molecule.GetAtoms()]
    extra_atoms = [i for i in mol_atoms if (i[1] not in keys)]
    connected_atom_sets = find_group_connected([i[0] for i in extra_atoms])
    mod_res_nums = get_modified_residues(connected_atom_sets, atm_id_to_res_dict)
    amino_acids = deduced_peptide.amino_acids
    seq = get_seq(amino_acids, mod_res_nums, mod_as_x=mod_as_x)

    deduced_peptide_mol = rdkit.Chem.MolFromSmiles(deduced_peptide.smiles)
    modification_mol = AllChem.DeleteSubstructs(mol, deduced_peptide_mol)
    modification_smiles = rdkit.Chem.MolToSmiles(modification_mol)

    mods_list_json = []

    if modification_smiles:

        modification_smiles_list = modification_smiles.split(".")
        # modification_mols = [ rdkit.Chem.MolFromSmiles(i) for i in modification_smiles_list ]

        for mod_smi in modification_smiles_list:
            j = {}
            j["modification_smiles"] = mod_smi
            mod_mol = rdkit.Chem.MolFromSmiles(mod_smi)
            j["connecting_residues"] = find_connecting_residues(
                mol=uc.molecule, mod_mol=mod_mol, atm_id_to_res_dict=atm_id_to_res_dict
            )
            mods_list_json.append(j)

    return seq, mods_list_json
