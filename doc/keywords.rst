.. _keyword_reference:

Keyword Reference
=================

Every keyword of :file:`simulation.cfg`, generated from the keyword table in the source code
(:file:`src/io/cfg_table_mod.f90`), which the configuration reader also uses: the table below is
therefore always that of the program. Each keyword links to its entry in the chapter that
describes it. A default of "required" means that the keyword must be given; keywords given
"in each phase block" belong in a phase block (see :ref:`material_behavior`). Rows of the form
"a component, then ..., one value per step" start with the component (e.g. ``33``) and give
one value per step (see :ref:`loading`).

.. include:: keywords_table.rst

.. _compatibility:

Option Compatibility
--------------------

The combinations of options that FEPX rejects (an error when the configuration file is read)
or ignores (a warning), generated from the rules in the source code
(:file:`src/io/cfg_rules_mod.f90`), which the configuration reader also uses. Each rule is also
stated at the entries of the options it concerns.

.. include:: compatibility_table.rst
