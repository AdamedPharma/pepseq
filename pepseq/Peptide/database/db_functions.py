import json
import sys
from typing import Union


def load_db_json(
    db_path: Union[str, None] = None):
    try:
        if db_path is not None:
            with open(db_path) as fp:
                db_json = json.load(fp)
            return db_json
        
    except FileNotFoundError as e:
        print(f"File {db_path} not found!", file=sys.stderr)
        return


def get_coding(db_json: dict) -> dict:
    """
    Get the coding dictionary from the database JSON file.
    :param db_json: The database JSON file.
    :type db_json: dict

    :return: The coding dictionary.
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
