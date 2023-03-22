import rdkit
from pepseq.get_peptide_json_from_pepseq_format import (
    get_attachment_points_on_sequence_json, get_pep_json)
from pepseq.Peptide.exceptions import (AttachmentPointsMismatchError,
                                       AttachmentPointsNonUniqueError,
                                       ExcessTildeError, InvalidSmilesError,
                                       NestedBracketError, ParenthesesError,
                                       UnattachedSmilesError, ValidationError)


def has_attachment_point(smiles):
    mol = rdkit.Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():
        if (atom.GetAtomicNum() == 0):
            return True
    return False


def validate_attachment_points_on_smiles(smiles_codes: list[str]):
    if smiles_codes is not None:
        for i in range(len(smiles_codes)):
            smiles_code = smiles_codes[i]
            for smiles_code in smiles_codes:
                if not has_attachment_point(smiles_code):
                    raise UnattachedSmilesError('SMILES code no %d has no attachment point to Peptide' % (i+1))


def validate_smiles_codes(smiles_codes: list[str] = None):
    if smiles_codes is not None:
        for i in range(len(smiles_codes)):
            smiles_code = smiles_codes[i]
            result = rdkit.Chem.MolFromSmiles(smiles_code)
            if result is None:
                raise InvalidSmilesError(
                    'SMILES code no %d was invalid. Could not construct molecule from SMILES code' % (i+1))


def get_attachment_points_on_smiles(smiles_code):
    attachment_points_ids = []
    mol = rdkit.Chem.MolFromSmiles(smiles_code)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 0:
            isotope = atom.GetIsotope()
            attachment_points_ids.append(isotope)
    return attachment_points_ids


def get_attachment_points_on_smiles_codes(smiles_codes):
    attachment_points_ids = []
    if smiles_codes is not None:
        for smiles_code in smiles_codes:
            attachment_points_ids += get_attachment_points_on_smiles(smiles_code)

    unique_attachment_points_ids = set(attachment_points_ids)
    if len(attachment_points_ids) > len(unique_attachment_points_ids):
        raise AttachmentPointsNonUniqueError("Attachment Points labels on SMILES are not unique.")

    return unique_attachment_points_ids


def validate_matching_attachment_points(pepseq, smiles_codes):
    symbols = get_pep_json(pepseq)['symbols']
    attachment_points_on_sequence = get_attachment_points_on_sequence_json(symbols)
    attachment_point_ids_on_sequence = set(attachment_points_on_sequence.keys())
    attachment_point_ids_on_smiles = get_attachment_points_on_smiles_codes(smiles_codes)
    if (attachment_point_ids_on_sequence != attachment_point_ids_on_smiles):
        raise AttachmentPointsMismatchError(
            'Attachment Points on Sequence: %s do not Match Attachment Points on Smiles: %s' % (
                str(attachment_point_ids_on_sequence), str(attachment_point_ids_on_smiles)))


def validate_termini(s):
    tilde_num = s.count("~")
    if tilde_num in [0, 1, 2]:
        return True
    elif tilde_num > 2:
        raise ExcessTildeError


def check_parentheses(s):
    """Return True if the parentheses in string s match, otherwise False."""
    j = 0
    for c in s:
        if c == "}":
            j -= 1
            if j < 0:
                raise ParenthesesError("Brackets do not match")
        elif c == "{":
            j += 1
    if j != 0:
        raise ParenthesesError("Brackets do not match")


def check_for_nested_brackets(s):
    open_bracket = False

    for c in s:
        if c == "{":
            if open_bracket:
                raise NestedBracketError("Found Nested '{','}' brackets.")
            else:
                open_bracket = True
        elif c == "}":
            if open_bracket:
                open_bracket = False
            else:
                raise ValidationError("Misplaced '}' Brackets")
