def get_coding(db_json: dict) -> dict:
    """
    Get the coding dictionary from the given database JSON.

    :param db_json: The database JSON.
    :type db_json: dict

    return: The coding dictionary.
    :rtype: dict
    """
    keys = [
        "l_proteogenic_3letter",
        "d_proteogenic_3letter",
        "d_proteogenic_2letter",
        "d_proteogenic_4letter",
        "modified_aa_codes",
    ]

    coding = db_json.get("coding").get("l_proteogenic_3letter")
    for key in keys:
        coding.update(db_json.get("coding").get(key))
    return coding
