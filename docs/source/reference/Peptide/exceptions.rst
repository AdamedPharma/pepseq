.. _commands:

*******************************
Command Line Interface Commands
*******************************

.. currentmodule:: pepseq


Pepseq provides Command Line Interface commands.

************************************
Creating SMILES from Pepseq Format
************************************

to create SMILES code from Modified Peptide given in Pepseq Format
you can use the ``pepseq.commands.pepseq_to_smiles()`` function:

.. autofunction:: commands.pepseq_to_smiles


.. admonition:: Example

    >>> import pepseq
    >>> pepseq.commands.pepseq_to_smiles('GC')
    '[H]NCC(=O)N[C@@H](CS)C(=O)O'



.. autoclass:: pepseq.Peptide.exceptions.ValidationError

.. autoclass:: pepseq.Peptide.exceptions.AttachmentPointsMismatchError

.. autoclass:: pepseq.Peptide.exceptions.AttachmentPointsNonUniqueError

.. autoclass:: pepseq.Peptide.exceptions.UnattachedSmilesError

.. autoclass:: pepseq.Peptide.exceptions.InvalidSmilesError
   
.. autoclass:: pepseq.Peptide.exceptions.InvalidSymbolError

.. autoclass:: pepseq.Peptide.exceptions.InvalidSequenceError

.. autoclass:: pepseq.Peptide.exceptions.NestedBracketError

.. autoclass:: pepseq.Peptide.exceptions.ParenthesesError
   
.. autoclass:: pepseq.Peptide.exceptions.TerminusError
   
.. autoclass:: pepseq.Peptide.exceptions.ExcessTildeError
