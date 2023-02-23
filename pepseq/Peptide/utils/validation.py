from pepseq.Peptide.exceptions import (ExcessTildeError, NestedBracketError,
                                       ParenthesesError, ValidationError)


def validate_termini(s):
    tilde_num = s.count("~")
    if tilde_num in [0, 1, 2]:
        return True
    elif tilde_num > 2:
        raise ExcessTildeError
    return


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
