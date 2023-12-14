class ValidationError(Exception):
    """
    Exception raised for validation errors in the Peptide module.

    Attributes:
        msg -- explanation of the error
    """

    def __init__(self, msg):
        super().__init__()
        self.msg = msg

    def __str__(self):
        return self.msg


class AttachmentPointsMismatchError(ValidationError):
    """
    Exception raised when the attachment points of a peptide are mismatched.
    """
    pass


class AttachmentPointsNonUniqueError(ValidationError):
    """
    Raised when the attachment points of a peptide are not unique.
    """
    pass


class UnattachedSmilesError(ValidationError):
    """
    Exception raised when a SMILES string is not attached to a peptide.
    """
    pass


class InvalidSmilesError(ValidationError):
    """
    Exception raised when an invalid SMILES string is encountered.
    """
    pass


class InvalidSymbolError(ValidationError):
    """Exception raised when an invalid symbol is encountered."""
    pass


class InvalidSequenceError(ValidationError):
    """Exception raised for invalid peptide sequences.

    Attributes:
        message -- explanation of the error
    """
    pass


class NestedBracketError(ValidationError):
    """
    Exception raised when there is a nested bracket error in a peptide sequence.
    """
    pass


class ParenthesesError(ValidationError):
    """Exception raised for errors related to parentheses in peptide sequences."""
    pass


class TerminusError(Exception):
    """Exception raised when only one of the termini has been defined with tilde (~) separator.
    
    When defining termini, both N and C termini must be specified.
    If either N or C terminus is not modified use: 'H~' for unmodified N terminus;
    '~OH' for unmodified C terminus.
    Example: 'H~SEQ~NH2' for aminated C terminus and N terminus unmodified
    or 'Ac~SEQ~OH' for acylated N terminus and C terminus unmodified.
    """

    def __init__(self):
        self.message = (
            "Only one of termini has been defined with tilde (~)"
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
    """Exception raised when the sequence string contains an incorrect number of tilde (~) termini separators."""

    def __init__(self):
        self.message = (
            " Sequence string must contain 0 or 2"
            + " tilde (~) termini separators. Use none if both N and C terminus"
            + " aren't modified; Two if one or both termini modified."
        )

        super().__init__(self.message)

    def __str__(self) -> str:
        return f"{self.message}"
