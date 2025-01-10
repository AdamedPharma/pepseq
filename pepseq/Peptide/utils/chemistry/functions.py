def get_connecting_bonds(bonds: list, native_atom_ids: list, remainder: list):
    """
    Get the connecting bonds between the native atoms and the remainder of the molecule.

    :param bonds: The bonds in the molecule.
    :type  bonds: list

    :param native_atom_ids: The native atom IDs.
    :type  native_atom_ids: list

    :param remainder: The remainder of the molecule.
    :type  remainder: list

    :return: The connecting bonds.
    :rtype:  list
    """
    connecting_bonds = []
    for bond in bonds:
        if bond[0] in native_atom_ids and bond[1] in remainder:
            connecting_bonds.append(bond)

        elif bond[0] in remainder and bond[1] in native_atom_ids:
            connecting_bonds.append((bond[1], bond[0]))
    return connecting_bonds
