def which_residue(
    res_atoms_tuple=(2, 4, 5),
    ca_res_dict={
        2: 1,
        6: 2,
        10: 3,
    },
):
    """
    this function takes:
        res_atoms_tuple <- tuple containing atom ids for residue (C alfa amongst them)
        ca_res_dict <- dictionary assigning res_id to ca_atom_id

    finds CA atom_id in tuple and returns residue_id for given CA
    """

    for atom_id in res_atoms_tuple:
        res_num = ca_res_dict.get(atom_id)
        if res_num is not None:
            return res_num
    return


def get_rows(aa_matches, ca_res_dict):
    """ """
    rows = []
    for aa in aa_matches.keys():
        nums_tuples = aa_matches[aa]
        for nums_tup in nums_tuples:
            res = which_residue(nums_tup, ca_res_dict)
            if res is not None:
                row = (nums_tup, res, aa)
                rows.append(row)
    return rows


def get_atm_id_to_potential_res_ids_names_dict(rows):
    """
    takes list of rows in form:
        ( atom_ids, res_id, res_name )

    returns dictionary in form:
        {
            atom_id : [ (res_id1, res_name1), (res_id2, res_name2) ]
        }
    """

    d = {}
    for atom_ids, res_id, res_name in rows:
        for atom_id in atom_ids:
            if d.get(atom_id) is None:
                d[atom_id] = []
            d[atom_id].append((res_id, res_name))
    return d
