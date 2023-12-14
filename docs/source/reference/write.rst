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


.. autofunction:: pepseq.write.to_pepseq

.. autofunction:: pepseq.write.to_smiles

.. autofunction:: pepseq.write.to_json
