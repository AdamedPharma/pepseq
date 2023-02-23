import unittest
from contextlib import contextmanager

from pepseq.functions import calculate, validate_pepseq
from pepseq.Peptide.exceptions import (AttachmentPointsMismatchError,
                                       ExcessTildeError, InvalidSmilesError,
                                       NestedBracketError, ParenthesesError,
                                       UnattachedSmilesError)
from pepseq.Peptide.utils.validation import (
    check_for_nested_brackets, check_parentheses,
    validate_attachment_points_on_smiles, validate_matching_attachment_points,
    validate_smiles_codes, validate_termini)


class TestExceptions(unittest.TestCase):

    @contextmanager
    def assertNotRaises(self, exc_type):
        try:
            yield None
        except exc_type:
            raise self.failureException('{} raised'.format(exc_type.__name__))

    def test_invalid_smiles_error(self):
        with self.assertRaises(InvalidSmilesError):
            validate_smiles_codes(['CC@%F'])

        with self.assertNotRaises(InvalidSmilesError):
            validate_termini(['CCC'])
        return

    def test_smiles_wo_attachment_points_error(self):
        with self.assertRaises(UnattachedSmilesError):
            validate_attachment_points_on_smiles(['CCC'])

        with self.assertNotRaises(UnattachedSmilesError):
            validate_attachment_points_on_smiles(['CCC[1*]'])
        return

    def test_validate_matching_attachment_points(self):
        correct_pairs = [('CCC{Cys(R1)}S{Cys(R2)}', ['CCC[1*]', 'CCC[2*]'])]

        for pepseq, smiles in correct_pairs:
            with self.assertNotRaises(AttachmentPointsMismatchError):
                validate_matching_attachment_points(pepseq, smiles)
            with self.assertNotRaises(AttachmentPointsMismatchError):
                calculate(pepseq, smiles)

        incorrect_pairs = [('CCC{Cys(R1)}S', ['CCC[2*]'])]

        for pepseq, smiles in incorrect_pairs:
            with self.assertRaises(AttachmentPointsMismatchError):
                validate_matching_attachment_points(pepseq, smiles)
            with self.assertRaises(AttachmentPointsMismatchError):
                calculate(pepseq, smiles)

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
                           'ACDEFG', 'ACGDEF']

        for correct_pepseq in correct_pepseqs:
            with self.assertNotRaises(ExcessTildeError):
                calculate(correct_pepseq)
            with self.assertNotRaises(NestedBracketError):
                calculate(correct_pepseq)
            with self.assertNotRaises(ParenthesesError):
                calculate(correct_pepseq)

        correct_pairs = [('ACDE{Cys(R1)}F', ['C[1*]CC']), ('AC{Cys(R1)}DEF', ['[1*]CSC'])]
        for pepseq, smiles in correct_pairs:
            with self.assertNotRaises(InvalidSmilesError):
                calculate(pepseq, smiles)
            with self.assertNotRaises(UnattachedSmilesError):
                calculate(pepseq, smiles)

        incorrect_pairs = [('ACDE{Cys(R1)}F', ['C[1*]CCCC@%F']), ('AC{Cys(R1)}DEF', ['[1*]CSCC@%FC'])]
        for pepseq, smiles in incorrect_pairs:
            with self.assertRaises(InvalidSmilesError):
                calculate(pepseq, smiles)

        unattached_pairs = [('ACDE{Cys(R1)}F', ['CCC']), ('ACDE{Cys(R1)}FI', ['CSC']),]

        for pepseq, smiles in unattached_pairs:
            with self.assertRaises(UnattachedSmilesError):
                calculate(pepseq, smiles)

        return
