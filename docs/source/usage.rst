####################################
Usage
####################################

.. _installation:

************************************
Installation
************************************

To use Pepseq, first install it using pip:

.. code-block:: console

   (.venv) $ pip install pepseq

************************************
Creating SMILES from Pepseq Format
************************************

to create SMILES code from Modified Peptide given in Pepseq Format
you can use the ``pepseq.commands.pepseq_to_smiles()`` function:

.. autofunction:: commands.pepseq_to_smiles


>>> import pepseq
>>> pepseq.commands.pepseq_to_smiles('GC')
'[H]NCC(=O)N[C@@H](CS)C(=O)O'

to calculate JSON representing Modified Peptide through peptide sequence in
Pepseq format and modification list in SMILES format
you can use the ``pepseq.commands.calculate_json_from()`` function:

.. autofunction:: commands.calculate_json_from

>>> import pepseq
>>> pepseq.commands.calculate_json_from
       ('{Cys(R1)}ACDAPEPsEQ{Cys(R2)}G{Cys(R3)}DEF',
        ['[1*]CNCC[2*]' --mod-smiles '[3*]CNCCSP']

To calculate sequence(s) in Pepseq Format and  list of modification SMILES from
Modified Peptide structure given in SMILES you can use the
``pepseq.commands.read_smiles()`` function:

.. autofunction:: commands.read_smiles


>>> import pepseq
>>> pepseq.commands.read-smiles('stuff_in.smi', db_path='augmented_db.json', out=path)


To augment monomers database with new monomers and modifications you can use
the ``pepseq.commands.augment_db_json_command()``

.. autofunction:: commands.augment_db_json_command


>>> import pepseq
>>> pepseq.commands.augment_db_json_command(sdf_path='sdf_file.sdf',
                            out='db_json_augmented.json'
                            )
>>> 

************************************
Augment Monomer Database
************************************

.. autofunction:: augmenting_db_json.get_R3

.. autofunction:: augmenting_db_json.set_default_exit_atom

.. autofunction:: augmenting_db_json.get_basic_smiles

.. autofunction:: augmenting_db_json.get_ro_json

.. autofunction:: augmenting_db_json.get_ro_json_from_row

.. autofunction:: augmenting_db_json.get_smiles_set

.. autofunction:: augmenting_db_json.get_row_jsons

.. autofunction:: augmenting_db_json.augment_db_json

.. autofunction:: augmenting_db_json.N_term_mod_smarts

.. autofunction:: augmenting_db_json.get_Nter_versions_cxsmarts_db

.. autofunction:: augmenting_db_json.get_Nter_versions

.. autofunction:: augmenting_db_json.change_exit_atom

.. autofunction:: augmenting_db_json.change_exit_atoms

.. autofunction:: augmenting_db_json.order_aas

.. autofunction:: augmenting_db_json.remove_radicals

.. autofunction:: augmenting_db_json.remove_radicals_from_db_smarts

************************************
pepseq
************************************

Backbone
########

.. autofunction:: pepseq.Backbone.generate_bb_smiles

.. autofunction:: pepseq.Backbone.generate_bb_mol

.. autoclass::  pepseq.Backbone.GetLongestPolymerWithin

.. autoclass::  pepseq.Backbone.MarkingPeptideBackbone

.. autoclass::  pepseq.Backbone.BreakingIntoResidueCandidateSubgraphs


BuildingModifiedPeptideFromPeptideJSON
########################################

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.add_internal_bond

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.add_disulfide_bond

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_attachment_points

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.find_atom

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.add_staple

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_peptide_json_from_sequence
   
.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_residue_symbols_from_sequence

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_molecule_from_sequence

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_smiles_from_sequence

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_molecule_from_json

.. autofunction:: pepseq.BuildingModifiedPeptideFromPeptideJSON.get_smiles_from_peptide_json


BuildPeptideJSONFromSMILES
################################

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


color_coding
################

Contains coding for different Modified Peptide components


draw
########

.. autofunction:: pepseq.draw.draw_pepseq

Drawer
########

.. autofunction:: pepseq.Drawer.get_start_x

.. autofunction:: pepseq.Drawer.get_rev_x_and_font_size

.. autofunction:: pepseq.Drawer.generate_kwargs_for_text_in_ellipse_balls

.. autofunction:: pepseq.Drawer.generate_kwargs_for_ellipse_balls

.. autofunction:: pepseq.Drawer.get_is_corner

.. autofunction:: pepseq.Drawer.get_direction

.. autofunction:: pepseq.Drawer.get_fragment_length

.. autofunction:: pepseq.Drawer.schema_layout_generator_from_symbols

.. autofunction:: pepseq.Drawer.get_fragment_kwargs

.. autofunction:: pepseq.Drawer.get_N_terminus_params

.. autofunction:: pepseq.Drawer.get_C_terminus_params

.. autofunction:: pepseq.Drawer.get_kwargs_from_symbols

.. autofunction:: pepseq.Drawer.draw_ellipse_ball

.. autofunction:: pepseq.Drawer.draw_ellipse_balls

.. autofunction:: pepseq.Drawer.draw_text_in_ellipse

.. autofunction:: pepseq.Drawer.draw_text_in_ellipse_balls

.. autofunction:: pepseq.Drawer.get_png_string_from_surface

.. autofunction:: pepseq.Drawer.draw_symbols


functions
################

.. autofunction:: pepseq.functions.validate_pepseq

.. autofunction:: pepseq.functions.calculate_pepseq_and_mods

.. autofunction:: pepseq.functions.calculate


get_peptide_json_from_pepseq_format
########################################

.. autofunction:: pepseq.get_peptide_json_from_pepseq_format.get_attachment_point_json

.. autofunction:: pepseq.get_peptide_json_from_pepseq_format.decompose_symbol

.. autofunction:: pepseq.get_peptide_json_from_pepseq_format.get_attachment_points_on_sequence_json

.. autofunction:: pepseq.get_peptide_json_from_pepseq_format.get_base_seq

.. autofunction:: pepseq.get_peptide_json_from_pepseq_format.get_single_modification_json

.. autofunction:: pepseq.get_peptide_json_from_pepseq_format.get_ext_mod_json
   
.. autofunction:: pepseq.get_peptide_json_from_pepseq_format.get_pep_json


read
########

.. autofunction:: pepseq.read.from_pepseq

.. autofunction:: pepseq.read.from_pepseq_and_mod_smiles

.. autofunction:: pepseq.read.from_smiles

.. autofunction:: pepseq.read.from_json


write
########

.. autofunction:: pepseq.write.to_pepseq

.. autofunction:: pepseq.write.to_smiles

.. autofunction:: pepseq.write.to_json


Peptide
########

exceptions
***********

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


models
***********


Peptide
======================


.. autofunction:: pepseq.Peptide.models.Peptide.get_smiles_descriptors

.. autoclass:: pepseq.Peptide.models.Peptide.Peptide


utils
***********


Parser 
======


.. autofunction:: pepseq.Peptide.utils.Parser.output_modified_residue
   
.. autofunction:: pepseq.Peptide.utils.Parser.parentheses_locs_list
   
.. autofunction:: pepseq.Peptide.utils.Parser.find_parentheses

.. autofunction:: pepseq.Peptide.utils.Parser.parse_canonical

.. autofunction:: pepseq.Peptide.utils.Parser.parse_canonical2
   
.. autofunction:: pepseq.Peptide.utils.Parser.find_termini
   
.. autofunction:: pepseq.Peptide.utils.Parser.get_canonical


validation 
==========


.. autofunction:: pepseq.Peptide.utils.validation.has_attachment_point
   
.. autofunction:: pepseq.Peptide.utils.validation.validate_attachment_points_on_smiles
   
.. autofunction:: pepseq.Peptide.utils.validation.validate_smiles_codes

.. autofunction:: pepseq.Peptide.utils.validation.get_attachment_points_on_smiles
   
.. autofunction:: pepseq.Peptide.utils.validation.get_attachment_points_on_smiles_codes

.. autofunction:: pepseq.Peptide.utils.validation.validate_matching_attachment_points
   
.. autofunction:: pepseq.Peptide.utils.validation.validate_termini
   
.. autofunction:: pepseq.Peptide.utils.validation.check_parentheses

.. autofunction:: pepseq.Peptide.utils.validation.check_for_nested_brackets


chemistry 
=========


cap_termini
-----------


.. autofunction:: pepseq.Peptide.utils.chemistry.cap_termini.prepare_ter_G

.. autofunction:: pepseq.Peptide.utils.chemistry.cap_termini.relabel_to_str

.. autofunction:: pepseq.Peptide.utils.chemistry.cap_termini.find_max_ResID

.. autofunction:: pepseq.Peptide.utils.chemistry.cap_termini.cap_terminus

.. autofunction:: pepseq.Peptide.utils.chemistry.cap_termini.cap_N_terminus

.. autofunction:: pepseq.Peptide.utils.chemistry.cap_termini.cap_C_terminus


functions
-----------


.. autofunction:: pepseq.Peptide.utils.chemistry.functions.get_connecting_bonds


mol_to_nx_translation
---------------------


.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.get_edge_tuple

.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.mol_to_nx
   
.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.nx_to_mol

.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.get_chiral_tag_int

.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.get_hybridization_int

.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.get_node_tuple
   
.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.nx_to_json

.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.get_mol_json
   
.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.mol_json_to_nx
   
.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.mol_json_to_mol
   
.. autofunction:: pepseq.Peptide.utils.chemistry.mol_to_nx_translation.draw_peptide_json


MonomerConnector
----------------


.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.smi_to_G

.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.is_R

.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.find_R
   
.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.find_N
   
.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.find_CO
   
.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.merge_graph

.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.get_residues_Gs
   
.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.merge_residue_graphs

.. autofunction:: pepseq.Peptide.utils.chemistry.MonomerConnector.get_molecule_from_list_of_residue_symbols


ProcessResidueCandidateGraph
----------------------------


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
