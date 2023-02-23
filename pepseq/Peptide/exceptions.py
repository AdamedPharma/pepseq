class ValidationError(Exception):
    pass


class InvalidSymbolError(ValidationError):
    pass


class NestedBracketError(ValidationError):
    pass


class ParenthesesError(ValidationError):
    pass


class TerminusError(Exception):
    def __init__(self):
        self.message = (
            " You have indicated only one of N and C termini, by tilde"
            + " separator. When indicating termini, both N and C termini"
            + " must be indicated."
            + " If either N or C terminus is not modified use: 'H~' for"
            + " unmodified N terminus; '~OH' for unmodified C terminus."
            + " E.g. 'H~SEQ~NH2' for aminated C terminus and N terminus"
            + " unmodified"
            + " or  'Ac~SEQ~OH  for acylated N terminus and C terminus"
            + " unmodified'."
        )
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}"


class ExcessTildeError(Exception):
    def __init__(self):
        self.message = (
            " Your sequence string contains too many more "
            + "than 2 tilde termini separators.  0 for both N and C terminus"
            + " not modified; 2 for one or both termini modified."
        )

        super().__init__(self.message)

    def __str__(self):
        return f"{self.message}"
