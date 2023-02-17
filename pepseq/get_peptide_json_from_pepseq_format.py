import rdkit

from pepseq.Peptide.utils.Parser import find_termini, parse_canonical2


def get_attachment_point_json(res_id, decomposition):
    ResName, attachment_point_id = decomposition
    d_atom_name = {"Cys": "SG", "Lys": "NZ"}

    if ResName in d_atom_name:
        AtomName = d_atom_name[ResName]
    else:
        AtomName = ""

    att_point_json = {
        "attachment_point_id": attachment_point_id,
        "ResID": str(res_id + 1),
        "AtomName": AtomName,
        "ResidueName": ResName,
    }
    return att_point_json


def decompose_symbol(s):
    if ("(" in s) and (")" in s):
        inside_bracket = s.split("(")[1].split(")")[0]
        before_bracket = s.split("(")[0]
        if "R" in inside_bracket:
            R_id = inside_bracket[1:]
            return before_bracket, R_id
    return s


def get_attachment_points_on_sequence_json(symbols):
    """ """
    att_points = {}

    for res_id in range(len(symbols)):
        symbol = symbols[res_id]
        decomposition = decompose_symbol(symbol)
        if type(decomposition) == tuple:
            res_name, attachment_point_id = decomposition
            attachment_point_json = get_attachment_point_json(res_id, decomposition)
            att_points[int(attachment_point_id)] = attachment_point_json
    return att_points


def get_base_seq(symbols):
    three_to_one = {"Cys": "C", "Lys": "K"}

    base_seq = ""

    for res_id in range(len(symbols)):
        symbol = symbols[res_id]
        decomposition = decompose_symbol(symbol)
        if type(decomposition) == tuple:
            res_name, attachment_point_id = decomposition
            if res_name in three_to_one:
                res_name = three_to_one[res_name]
            symbol = res_name

        if len(symbol) > 1:
            symbol = "{%s}" % (symbol)
        base_seq = base_seq + symbol
    return base_seq


def get_single_modification_json(attachment_points_on_sequence, mod_smiles: str):
    mod_mol = rdkit.Chem.MolFromSmiles(mod_smiles)
    mod_atoms = mod_mol.GetAtoms()

    radical_ids = set(
        [atom.GetIsotope() for atom in mod_atoms if (atom.GetAtomicNum() == 0)]
    )
    min_att_points = {
        radical_id: attachment_points_on_sequence[radical_id]
        for radical_id in radical_ids
    }

    max_attachment_point_id = max([int(i) for i in min_att_points.keys()])

    ext_mod = {
        "smiles": mod_smiles,
        "max_attachment_point_id": max_attachment_point_id,
        "attachment_points_on_sequence": min_att_points,
    }
    return ext_mod


def get_ext_mod_json(pepseq, smiles: list):
    symbols = parse_canonical2(pepseq)
    attachment_points_on_sequence = get_attachment_points_on_sequence_json(symbols)
    if attachment_points_on_sequence.keys():

        ext_mod_jsons = []
        for mod_smiles in smiles:
            ext_mod_json = get_single_modification_json(
                attachment_points_on_sequence, mod_smiles
            )
            ext_mod_jsons.append(ext_mod_json)
    return ext_mod_jsons


def get_pep_json(pepseq_format, db_json, mod_smiles_list=None):
    """

    Input:


        pepseq_string:

            str = string in pepseq format
            H~H{aMeAla}EGTFTSDVSSYLEG{Cys(R1)}AAKEFI{Cys(R2)}WLVRGRG~OH
        where H~ is N-terminus; ~OH is C_terminus, {aMeAla} is modified
        amino acid; {Cys(R1)} - is amino acid
        with staple attached, {Cys(R1)} - amino acid with staple attached


        mod_smiles:

            SMILES string (e.g. '[1*]C[2*]') - showing the structure of
                modification with attachment
            points:

                { Cys(R1) } <- is attached in [1*] attachment point on staple
                { Cys(R2) } <- is attached in [2*] attachment point on staple

    Output:

        peptide_json:

            JSON containing info about modified peptide with

                'sequence':

                'internal_modifications':

                'external_modifications':

    """

    N_terminus, C_terminus, pepseq = find_termini(pepseq_format, db_json)
    symbols = parse_canonical2(pepseq)
    base_seq = get_base_seq(symbols)

    pep_json = {
        "length": len(symbols),
        "sequence": base_seq,
        "internal_modifications": [],
        "C_terminus": C_terminus,
        "N_terminus": N_terminus,
        "pepseq_format": pepseq_format,
    }

    if mod_smiles_list is not None:
        ext_mod = get_ext_mod_json(pepseq, mod_smiles_list)
        if ext_mod is not None:
            pep_json["external_modifications"] = ext_mod
        else:
            pep_json["external_modifications"] = []
    else:
        pep_json["external_modifications"] = []

    return pep_json
