.. _pepseq_validation:

pepseq_validation
^^^^^^^^^^^^^^^^^^

.. currentmodule:: pepseq


Pepseq provides Command Line Interface commands.

validate_termini
"""""""""""""""""

to create SMILES code from Modified Peptide given in Pepseq Format
you can use the ``pepseq.commands.pepseq_to_smiles()`` function:

.. autofunction:: commands.pepseq_to_smiles


.. admonition:: Example

    >>> import pepseq
    >>> pepseq.commands.pepseq_to_smiles('GC')
    '[H]NCC(=O)N[C@@H](CS)C(=O)O'


.. autofunction:: pepseq.Peptide.utils.pepseq_validation.validate_termini
   
.. autofunction:: pepseq.Peptide.utils.pepseq_validation.check_parentheses

.. autofunction:: pepseq.Peptide.utils.pepseq_validation.check_for_nested_brackets
