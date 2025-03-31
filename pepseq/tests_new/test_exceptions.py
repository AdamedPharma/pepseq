import unittest
from contextlib import contextmanager

from pepseq.functions import calculate, validate  # _pepseq
from pepseq.Peptide.exceptions import (
    AttachmentPointsMismatchError,
    AttachmentPointsNonUniqueError,
    ExcessTildeError,
    InvalidSmilesError,
    InvalidSymbolError,
    NestedBracketError,
    ParenthesesError,
    UnattachedSmilesError,
)

from pepseq.Peptide.utils.pepseq_validation import (
    check_for_nested_brackets,
    check_parentheses,
    validate_termini,
    validate_pepseq,
)


from pepseq.Peptide.utils.smiles_validation import (
    validate_attachment_points_on_smiles,
    validate_smiles_codes,
)

from pepseq.Peptide.utils.validation import (
    validate_matching_attachment_points,
)


correct_pepseqs = [
    "ACDEF",
    "H~ACDEF",
    "H~ACDEF~OH",
    "ACDEF~OH",
    "ACDEFG",
    "ACGDEF",
]
# correct_pepseqs = ['ACDE{Cys(R1)}F']

correct_pairs = [("CCC{Cys(R1)}S{Cys(R2)}", ["CCC[*:1]", "CCC[*:2]"])]


class TestExceptions(unittest.TestCase):
    """
    TestExceptions is a test case class that contains various unit tests for
     validating peptide sequences and SMILES codes.
    It includes custom assertion methods and tests for different types of errors
     that can occur during validation and calculation processes.
    Methods:
        assertNotRaises(exc_type):
            A context manager that asserts a specific exception type is not raised.
        test_invalid_symbol():
            Tests that an InvalidSymbolError is raised when an invalid residue symbol is encountered.
        test_invalid_smiles_error():
            Tests that an InvalidSmilesError is raised for invalid SMILES codes and not raised for valid ones.
        test_smiles_wo_attachment_points_error():
            Tests that an UnattachedSmilesError is raised when SMILES codes lack attachment points.
        test_validate_matching_attachment_points():
            Tests that an AttachmentPointsMismatchError is raised when attachment
             points on sequences and SMILES codes do not match.
        test_validate_unique_attachment_points():
            Tests that an AttachmentPointsNonUniqueError is raised when attachment
             points on sequences and SMILES codes are not unique.
        test_excess_tilde_error():
            Tests that an ExcessTildeError is raised for sequences with excess tildes.
        test_nested_bracket_error():
            Tests that a NestedBracketError is raised for sequences with nested brackets.
        test_parentheses_error():
            Tests that a ParenthesesError is raised for sequences with mismatched parentheses.
        test_validate_pepseq():
            Tests the validate_pepseq function for various errors including ExcessTildeError,
             NestedBracketError, and ParenthesesError.
        test_calculate():
            Tests the calculate function for various errors including ExcessTildeError,
             NestedBracketError, ParenthesesError, InvalidSmilesError, and UnattachedSmilesError.
    """
    @contextmanager
    def assertNotRaises(self, exc_type):
        """
        Asserts that a block of code does not raise a specific exception.
        This method is used as a context manager to wrap the code that is expected
        not to raise the specified exception type. If the exception is raised,
        the test will fail with an appropriate message.
        Args:
            exc_type (Exception): The type of exception that should not be raised.
        Raises:
            self.failureException: If the specified exception type is raised.
        """
        try:
            yield None
        except exc_type:
            raise self.failureException("{} raised".format(exc_type.__name__))

    def test_invalid_symbol(self):
        """
        Test that the calculate function raises an InvalidSymbolError when an invalid symbol is provided.
        This test checks that the calculate function correctly identifies and raises an error
        when it encounters a residue symbol that is not found in the database.
        Raises:
            InvalidSymbolError: If the residue symbol is not found in the database.
        """
        with self.assertRaises(InvalidSymbolError) as cm:
            calculate("CSCU")
        self.assertEqual(str(cm.exception), "Residue Symbols: U not found in database.")

        return

    def test_invalid_smiles_error(self):
        """
        Test the validation of SMILES codes and termini.
        This test case checks the following:
        - An InvalidSmilesError is raised when an invalid SMILES code is provided.
        - No exception is raised when valid termini are provided.
        The test uses the following inputs:
        - Invalid SMILES code: ["CC@%F"]
        - Valid termini: ["CCC"]
        """
        with self.assertRaises(InvalidSmilesError):
            validate_smiles_codes(["CC@%F"])

        with self.assertNotRaises(InvalidSmilesError):
            validate_termini(["CCC"])
        return

    def test_smiles_wo_attachment_points_error(self):
        """
        Test cases for the `validate_attachment_points_on_smiles` function to ensure it raises
        `UnattachedSmilesError` when SMILES strings do not contain attachment points.
        Test cases:
        - SMILES string without attachment points should raise `UnattachedSmilesError` with the
          message "SMILES code no 1 has no attachment point to Peptide".
        - A list of SMILES strings where the second string lacks attachment points should raise
          `UnattachedSmilesError` with the message "SMILES code no 2 has no attachment point to Peptide".
        - A list of SMILES strings where both strings lack attachment points should raise
          `UnattachedSmilesError` with the message "SMILES code no 1 has no
           attachment point to Peptide\nSMILES code no 2 has no attachment point to Peptide".
        - A SMILES string with an attachment point should not raise `UnattachedSmilesError`.
        """
        with self.assertRaises(UnattachedSmilesError) as cm:
            validate_attachment_points_on_smiles(["CCC"])
        self.assertEqual(
            str(cm.exception), "SMILES code no 1 has no attachment point to Peptide"
        )

        with self.assertRaises(UnattachedSmilesError) as cm:
            validate_attachment_points_on_smiles(["CCC[*:1]", "CCC"])
        self.assertEqual(
            str(cm.exception), "SMILES code no 2 has no attachment point to Peptide"
        )

        with self.assertRaises(UnattachedSmilesError) as cm:
            validate_attachment_points_on_smiles(["CCCC", "CCC"])
        self.assertEqual(
            str(cm.exception),
            "SMILES code no 1 has no attachment point to Peptide\nSMILES code no 2 has no attachment point to Peptide",
        )

        with self.assertNotRaises(UnattachedSmilesError):
            validate_attachment_points_on_smiles(["CCC[*:1]"])
        return

    def test_validate_matching_attachment_points(self):
        """
        Test the validation of matching attachment points between peptide sequences and SMILES strings.
        This test case checks both correct and incorrect pairs of peptide sequences and SMILES strings
        to ensure that the validation functions `validate_matching_attachment_points` and `calculate`
        behave as expected.
        - For correct pairs, it verifies that no `AttachmentPointsMismatchError` is raised.
        - For incorrect pairs, it verifies that an `AttachmentPointsMismatchError` is raised and
          that the error message matches the expected string.
        Correct pairs are tested to ensure that valid sequences and SMILES strings pass the validation
        without errors. Incorrect pairs are tested to ensure that mismatched attachment points are
        correctly identified and reported.
        """
        for pepseq, smiles in correct_pairs:
            with self.assertNotRaises(AttachmentPointsMismatchError):
                validate_matching_attachment_points(pepseq, smiles)
            with self.assertNotRaises(AttachmentPointsMismatchError):
                calculate(pepseq, smiles)

        incorrect_pairs = [("CCC{Cys(R1)}S", ["CCC[*:2]"])]

        for pepseq, smiles in incorrect_pairs:
            with self.assertRaises(AttachmentPointsMismatchError) as cm:
                validate_matching_attachment_points(pepseq, smiles)
            self.assertEqual(
                str(cm.exception),
                "Attachment Points on Sequence: {1} do not Match Attachment Points on Smiles: {2}",
            )

            with self.assertRaises(AttachmentPointsMismatchError):
                validate(pepseq, smiles)

    def test_validate_unique_attachment_points(self):
        """
        Test the validation of unique attachment points in peptide sequences and SMILES strings.
        This test checks that the `validate_matching_attachment_points` and `calculate` functions
        do not raise an `AttachmentPointsNonUniqueError` for correct pairs of peptide sequences
        and SMILES strings. It also verifies that the `validate_matching_attachment_points` and
         `validate` functions raise an `AttachmentPointsNonUniqueError` for incorrect pairs.
        The correct pairs are expected to have unique attachment points, while the incorrect pairs
        are expected to have non-unique attachment points.
        Raises:
            AttachmentPointsNonUniqueError: If the attachment points in the peptide sequence and
            SMILES string are not unique.
        """
        for pepseq, smiles in correct_pairs:
            with self.assertNotRaises(AttachmentPointsNonUniqueError):
                validate_matching_attachment_points(pepseq, smiles)
            with self.assertNotRaises(AttachmentPointsNonUniqueError):
                calculate(pepseq, smiles)

        incorrect_pairs = [
            ("CCC{Cys(R1)}S", ["[*:2]CCC[*:2]"]),
            ("C{Cys(R2)}CC{Cys(R1)}S{Cys(R1)}", ["CCC[*:1]", "CCC[*:2]"]),
        ]

        for pepseq, smiles in incorrect_pairs:
            with self.assertRaises(AttachmentPointsNonUniqueError):
                validate_matching_attachment_points(pepseq, smiles)
            with self.assertRaises(AttachmentPointsNonUniqueError):
                validate(pepseq, smiles)
        return

    def test_excess_tilde_error(self):
        """
        Test the validate_termini function for handling excess tilde errors.
        This test case checks if the validate_termini function raises an
        ExcessTildeError when the input string contains more tildes than expected.
        It also verifies that the function does not raise an error for valid
        peptide sequences.
        The test performs the following checks:
        1. Ensures that an ExcessTildeError is raised when the input string
           "H~CS~OH~OH" is passed to the validate_termini function.
        2. Ensures that no ExcessTildeError is raised for each peptide sequence
           in the correct_pepseqs list.
        """
        with self.assertRaises(ExcessTildeError):
            validate_termini("H~CS~OH~OH")

        for correct_pepseq in correct_pepseqs:
            with self.assertNotRaises(ExcessTildeError):
                validate_termini(correct_pepseq)

    def test_nested_bracket_error(self):
        """
        Test the check_for_nested_brackets function for handling nested brackets.
        This test verifies that the function correctly raises a NestedBracketError
        when nested brackets are present in the input string, and does not raise
        an error when the brackets are properly nested.
        - Raises NestedBracketError for input "AC{{Cys(R1)}}DEF"
        - Does not raise NestedBracketError for input "AC{Cys(R1)}DEF"
        """
        with self.assertRaises(NestedBracketError):
            check_for_nested_brackets("AC{{Cys(R1)}}DEF")

        with self.assertNotRaises(NestedBracketError):
            check_for_nested_brackets("AC{Cys(R1)}DEF")

        return

    def test_parentheses_error(self):
        """
        Test the handling of parentheses errors in peptide sequences.
        This test checks that the `check_parentheses` function raises a
         `ParenthesesError` for incorrect peptide sequences with unmatched
         parentheses. It also verifies that the `check_for_nested_brackets`
         function does not raise a `ParenthesesError` for correct peptide sequences.
        The test cases include:
        - Incorrect peptide sequences with unmatched parentheses: ["AC{DEF", "ACD}EF"]
        - Correct peptide sequences: (to be defined in the test)
        Raises:
            ParenthesesError: If the peptide sequence has unmatched parentheses.
        """
        incorrect_pepseqs = ["AC{DEF", "ACD}EF"]
        for incorrect_pepseq in incorrect_pepseqs:
            with self.assertRaises(ParenthesesError):
                check_parentheses(incorrect_pepseq)
        for correct_pepseq in correct_pepseqs:
            with self.assertNotRaises(ParenthesesError):
                check_for_nested_brackets(correct_pepseq)
        return

    def test_validate_pepseq(self):
        """
        Test the validate_pepseq function for various error conditions.
        This test case checks the following scenarios:
        - ExcessTildeError is raised when the peptide sequence contains an excess tilde.
        - NestedBracketError is raised when the peptide sequence contains nested brackets.
        - ParenthesesError is raised when the peptide sequence contains unbalanced parentheses.
        It also verifies that no errors are raised for correct peptide sequences.
        Raises:
            ExcessTildeError: If the peptide sequence contains an excess tilde.
            NestedBracketError: If the peptide sequence contains nested brackets.
            ParenthesesError: If the peptide sequence contains unbalanced parentheses.
        """
        with self.assertRaises(ExcessTildeError):
            validate_pepseq("H~CS~OH~OH")

        with self.assertRaises(NestedBracketError):
            validate_pepseq("AC{{Cys(R1)}}DEF")

        incorrect_pepseqs = ["AC{DEF", "ACD}EF"]
        for incorrect_pepseq in incorrect_pepseqs:
            with self.assertRaises(ParenthesesError):
                validate_pepseq(incorrect_pepseq)

        for correct_pepseq in correct_pepseqs:
            with self.assertNotRaises(ExcessTildeError):
                validate_pepseq(correct_pepseq)
            with self.assertNotRaises(NestedBracketError):
                validate_pepseq(correct_pepseq)
            with self.assertNotRaises(ParenthesesError):
                validate_pepseq(correct_pepseq)
        return

    def test_calculate(self):
        """
        Test the `calculate` function for various error conditions and valid cases.
        This test covers the following scenarios:
        - Excess tilde error in the peptide sequence.
        - Nested bracket error in the peptide sequence.
        - Parentheses error in the peptide sequence.
        - Valid peptide sequences that should not raise errors.
        - Valid peptide and SMILES pairs that should not raise errors.
        - Invalid SMILES strings that should raise `InvalidSmilesError`.
        - Unattached SMILES strings that should raise `UnattachedSmilesError`.
        The test uses the `validate` and `calculate` functions to check for these conditions.
        """
        with self.assertRaises(ExcessTildeError):
            validate("H~CS~OH~OH")

        with self.assertRaises(NestedBracketError):
            validate("AC{{Cys(R1)}}DEF")

        incorrect_pepseqs = ["AC{DEF", "ACD}EF"]
        for incorrect_pepseq in incorrect_pepseqs:
            with self.assertRaises(ParenthesesError):
                validate(incorrect_pepseq)

        for correct_pepseq in correct_pepseqs:
            with self.assertNotRaises(ExcessTildeError):
                calculate(correct_pepseq)
            with self.assertNotRaises(NestedBracketError):
                calculate(correct_pepseq)
            with self.assertNotRaises(ParenthesesError):
                calculate(correct_pepseq)

        correct_pairs = [
            ("ACDE{Cys(R1)}F", ["C[*:1]CC"]),
            ("AC{Cys(R1)}DEF", ["[*:1]CSC"]),
        ]
        for pepseq, smiles in correct_pairs:
            with self.assertNotRaises(InvalidSmilesError):
                calculate(pepseq, smiles)
            with self.assertNotRaises(UnattachedSmilesError):
                calculate(pepseq, smiles)

        incorrect_pairs = [
            ("ACDE{Cys(R1)}F", ["C[*:1]CCCC@%F"]),
            ("AC{Cys(R1)}DEF", ["[*:1]CSCC@%FC"]),
        ]
        for pepseq, smiles in incorrect_pairs:
            with self.assertRaises(InvalidSmilesError):
                validate(pepseq, smiles)

        unattached_pairs = [
            ("ACDE{Cys(R1)}F", ["CCC"]),
            ("ACDE{Cys(R1)}FI", ["CSC"]),
        ]

        for pepseq, smiles in unattached_pairs:
            with self.assertRaises(UnattachedSmilesError):
                validate(pepseq, smiles)

        return
