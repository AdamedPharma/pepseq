import unittest
from contextlib import contextmanager

from pepseq.functions import calculate, validate_pepseq
from pepseq.Peptide.exceptions import (ExcessTildeError, NestedBracketError,
                                       ParenthesesError)
from pepseq.Peptide.utils.validation import (check_for_nested_brackets,
                                             check_parentheses,
                                             validate_termini)


class TestExceptions(unittest.TestCase):

    @contextmanager
    def assertNotRaises(self, exc_type):
        try:
            yield None
        except exc_type:
            raise self.failureException('{} raised'.format(exc_type.__name__))

    def test_excess_tilde_error(self):
        with self.assertRaises(ExcessTildeError):
            validate_termini('H~CS~OH~OH')

        correct_pepseqs = ['ACDEF', 'H~ACDEF', 'H~ACDEF~OH', 'ACDEF~OH']

        for correct_pepseq in correct_pepseqs:
            with self.assertNotRaises(ExcessTildeError):
                validate_termini(correct_pepseq)

    def test_nested_bracket_error(self):
        with self.assertRaises(NestedBracketError):
            check_for_nested_brackets('AC{{Cys(R1)}}DEF')

        with self.assertNotRaises(NestedBracketError):
            check_for_nested_brackets('AC{Cys(R1)}DEF')

        return

    def test_parentheses_error(self):
        incorrect_pepseqs = ['AC{DEF', 'ACD}EF']
        for incorrect_pepseq in incorrect_pepseqs:
            with self.assertRaises(ParenthesesError):
                check_parentheses(incorrect_pepseq)
        correct_pepseqs = ['ACDE{Cys(R1)}F']
        for correct_pepseq in correct_pepseqs:
            with self.assertNotRaises(ParenthesesError):
                check_for_nested_brackets(correct_pepseq)
        return

    def test_validate_pepseq(self):
        with self.assertRaises(ExcessTildeError):
            validate_pepseq('H~CS~OH~OH')

        with self.assertRaises(NestedBracketError):
            validate_pepseq('AC{{Cys(R1)}}DEF')

        incorrect_pepseqs = ['AC{DEF', 'ACD}EF']
        for incorrect_pepseq in incorrect_pepseqs:
            with self.assertRaises(ParenthesesError):
                validate_pepseq(incorrect_pepseq)

        correct_pepseqs = ['ACDEF', 'H~ACDEF', 'H~ACDEF~OH', 'ACDEF~OH',
                           'ACDE{Cys(R1)}F', 'AC{Cys(R1)}DEF']

        for correct_pepseq in correct_pepseqs:
            with self.assertNotRaises(ExcessTildeError):
                validate_pepseq(correct_pepseq)
            with self.assertNotRaises(NestedBracketError):
                validate_pepseq(correct_pepseq)
            with self.assertNotRaises(ParenthesesError):
                validate_pepseq(correct_pepseq)
        return

    def test_calculate(self):
        with self.assertRaises(ExcessTildeError):
            calculate('H~CS~OH~OH')

        with self.assertRaises(NestedBracketError):
            calculate('AC{{Cys(R1)}}DEF')

        incorrect_pepseqs = ['AC{DEF', 'ACD}EF']
        for incorrect_pepseq in incorrect_pepseqs:
            with self.assertRaises(ParenthesesError):
                calculate(incorrect_pepseq)

        correct_pepseqs = ['ACDEF', 'H~ACDEF', 'H~ACDEF~OH', 'ACDEF~OH',
                           'ACDE{Cys(R1)}F', 'AC{Cys(R1)}DEF']

        for correct_pepseq in correct_pepseqs:
            with self.assertNotRaises(ExcessTildeError):
                calculate(correct_pepseq)
            with self.assertNotRaises(NestedBracketError):
                calculate(correct_pepseq)
            with self.assertNotRaises(ParenthesesError):
                calculate(correct_pepseq)
        return
