def get_coding(db_json: dict) -> dict:
    """
    Get the coding dictionary from the database JSON file.
    :param db_json: The database JSON file.

    :return: The coding dictionary.
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
