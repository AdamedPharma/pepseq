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


.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.get_cx_smarts_db

.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.decompose_peptide_smiles

.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.get_terminal_smiles_building_block

.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.get_C_terminal_smiles_building_block

.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.get_N_terminal_smiles_building_block

.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.smiles_are_identical

.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.get_term_symbol

.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.get_c_term_from_peptide_json

.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.get_n_term_from_peptide_json

.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.output_modified_residue
   
.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.append_pepseq_R_info

.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.decompose_peptide_smiles_with_termini
   
.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.from_smiles_to_pepseq_and_mod_smiles_strings

.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.from_smiles_to_pepseq_and_one_mod_smiles_strings

.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.mark_external_modifications_on_seq

.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.mark_internal_modifications_on_seq
   
.. autofunction:: pepseq.BuildPeptideJSONFromSMILES.print_sequence
