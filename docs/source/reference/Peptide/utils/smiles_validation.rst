.. _smiles_validation:

smiles_validation
^^^^^^^^^^^^^^^^^^

.. currentmodule:: pepseq


Pepseq provides Command Line Interface commands.

has_attachment_point
"""""""""""""""""""""

to create SMILES code from Modified Peptide given in Pepseq Format
you can use the ``pepseq.commands.pepseq_to_smiles()`` function:

.. autofunction:: commands.pepseq_to_smiles


.. admonition:: Example

    >>> import pepseq
    >>> pepseq.commands.pepseq_to_smiles('GC')
    '[H]NCC(=O)N[C@@H](CS)C(=O)O'


.. autofunction:: pepseq.Peptide.utils.smiles_validation.has_attachment_point
   
.. autofunction:: pepseq.Peptide.utils.smiles_validation.validate_attachment_points_on_smiles
