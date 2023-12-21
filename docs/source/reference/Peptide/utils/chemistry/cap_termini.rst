.. _cap_termini:

cap_termini
"""""""""""

.. currentmodule:: pepseq


Pepseq provides Command Line Interface commands.

************************************
prepare_ter_G
************************************

to create SMILES code from Modified Peptide given in Pepseq Format
you can use the ``pepseq.commands.pepseq_to_smiles()`` function:

.. autofunction:: commands.pepseq_to_smiles


.. admonition:: Example

    >>> import pepseq
    >>> pepseq.commands.pepseq_to_smiles('GC')
    '[H]NCC(=O)N[C@@H](CS)C(=O)O'


.. autofunction:: pepseq.Peptide.utils.chemistry.cap_termini.prepare_ter_G

.. autofunction:: pepseq.Peptide.utils.chemistry.cap_termini.relabel_to_str

.. autofunction:: pepseq.Peptide.utils.chemistry.cap_termini.find_max_ResID

.. autofunction:: pepseq.Peptide.utils.chemistry.cap_termini.cap_terminus

.. autofunction:: pepseq.Peptide.utils.chemistry.cap_termini.cap_N_terminus

.. autofunction:: pepseq.Peptide.utils.chemistry.cap_termini.cap_C_terminus
