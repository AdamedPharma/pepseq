class ValidationError(Exception):
    def __init__(self, msg):
        super().__init__()
        self.msg = msg

    def __str__(self):
        return self.msg


class AttachmentPointsMismatchError(ValidationError):
    pass


class AttachmentPointsNonUniqueError(ValidationError):
    pass


class UnattachedSmilesError(ValidationError):
    pass


class InvalidSmilesError(ValidationError):
    pass


class InvalidSymbolError(ValidationError):
    pass


class InvalidSequenceError(ValidationError):
    pass


class NestedBracketError(ValidationError):
    pass


class ParenthesesError(ValidationError):
    pass


class TerminusError(Exception):
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
    def __init__(self):
        self.message = (
            " Sequence string must contain 0 or 2"
            + " tilde (~) termini separators. Use none if both N and C terminus"
            + " aren't modified; Two if one or both termini modified."
        )

        super().__init__(self.message)

    def __str__(self) -> str:
        return f"{self.message}"
