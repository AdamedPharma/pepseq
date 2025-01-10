class ValidationError(Exception):
    """
    Base class for validation errors.
    """
    def __init__(self, msg):
        super().__init__()
        self.msg = msg

    def __str__(self):
        return self.msg


class AttachmentPointsMismatchError(ValidationError):
    """
    Raised when the number of attachment points on the sequence does not match the
     number of attachment points in the SMILES string.
    """
    pass


class AttachmentPointsNonUniqueError(ValidationError):
    """
    Raised when attachment points are not unique.
    """
    pass


class UnattachedSmilesError(ValidationError):
    """
    Raised when an attachment point in the SMILES string is not attached to the sequence.
    """
    pass


class InvalidSmilesError(ValidationError):
    """
    Raised when the SMILES string is invalid.
    """
    pass


class InvalidSymbolError(ValidationError):
    """
    Raised when the symbol is invalid
    """
    pass


class InvalidSequenceError(ValidationError):
    """
    Raised when the sequence is invalid.
    """
    pass


class NestedBracketError(ValidationError):
    """
    Raised when a bracket is nested within another bracket
    """
    pass


class ParenthesesError(ValidationError):
    """
    Raised when parentheses are not balanced.
    """
    pass


class TerminusError(Exception):
    """
    Raised when termini are not properly defined.
    """
    def __init__(self):
        self.message = (
            " Only one of termini has been defined with tilde (~)"
            + " separator. When defining termini, both N and C termini"
            + " must be specified."
            + " If either N or C terminus is not modified use: 'H~' for"
            + " unmodified N terminus; '~OH' for unmodified C terminus."
            + " E.g. 'H~SEQ~NH2' for aminated C terminus and N terminus"
            + " unmodified"
            + " or  'Ac~SEQ~OH  for acylated N terminus and C terminus"
            + " unmodified'."
        )
        super().__init__(self.message)

    def __str__(self) -> str:
        return f"{self.message}"


class ExcessTildeError(Exception):
    """
    Raised when there are more than 2 tildes in the sequence.
    """
    def __init__(self):
        self.message = (
            " Sequence string must contain 0 or 2"
            + " tilde (~) termini separators. Use none if both N and C terminus"
            + " aren't modified; Two if one or both termini modified."
        )

        super().__init__(self.message)

    def __str__(self) -> str:
        return f"{self.message}"
