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

.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.smi_to_G

.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.is_R

.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.find_R
   
.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.find_N
   
.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.find_CO
   
.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.merge_graph

.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.get_residues_Gs
   
.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.merge_residue_graphs

.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.get_molecule_from_list_of_residue_symbols
