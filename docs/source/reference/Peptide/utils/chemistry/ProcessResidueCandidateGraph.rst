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


*****************************************************
Reading Sequence and Modification from SMILES string
*****************************************************

.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_match

.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_matches
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_match_cover
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.match_molecular_graph_to_res_id

.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_mod_graphs
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_n_subst_dict
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.filter_n_subst
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_res_matches

.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.propagate_matches_on_molecular_graph
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_connecting
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_connecting_edges

.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_atom_pairs
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.process_internal_connections
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.sorted_connection

.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.add_connection_point_to_molecular_graph
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_residue_id_and_atoms

.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.process_external_modification
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.process_external_connections

.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.split_connections_by_type
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_modification_graphs_from_fragment
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.decompose

.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_internal_connections_subgraph_tuples
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_subgraph_tuples

.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.get_connections
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.full_decomposition
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.translate_attachment_points_on_seq

.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.translate_mod_smiles
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.translate_external_modification

.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.sequence_dict_to_string
   
.. autofunction:: pepseq.Peptide.utils.chemistry.ProcessResidueCandidateGraph.decompose_residues_internal
