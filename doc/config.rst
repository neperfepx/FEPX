.. _config_file:

Configuration File (:file:`simulation.cfg`)
===========================================

A simulation input consists of two files with fixed names: the mesh file,
:file:`simulation.msh`, and the configuration file, :file:`simulation.cfg`. The mesh file
describes the polycrystal as a finite element mesh, where the grains are element sets, with
additional data such as phases, crystal orientations, node sets, face sets and mesh partitions.
It is normally written by Neper and not edited (see :ref:`mesh_file`). The configuration file
completes the definition of the polycrystal with the material model and its parameters, and
defines the simulation (boundary conditions, loading history, solver parameters, output). It is
written by the user. The two files must be consistent with each other (number of phases, node
and face set names, partitions).

Syntax
------

The file is a sequence of lines of the form::

    <keyword> <value1> <value2> ...

- One keyword per line, followed by its values, separated by blanks (spaces or tabs).
- ``#`` starts a comment, which runs to the end of the line; ``## Section`` headers are
  therefore comments. Blank lines are ignored. Indentation is free.
- Keywords are lower case. Word values are compared without regard to case (``fcc`` and
  ``FCC`` are the same).
- Lists are separated by blanks, e.g. ``slip_families {110}<111> {112}<111>``.
- Keywords may be given in any order, except that a per-phase keyword belongs to the block
  opened by the last ``phase <i>`` line (see :ref:`material_behavior`).
- Per-step rows give one value per step, and ``-`` for "not used in this step" (see
  :ref:`loading`).

A typical file has three parts (an elastic polycrystal pulled in two steps):

.. literalinclude:: ../tests/regression/cases/elastic_poly_fcc/simulation.cfg
   :start-at: ## Material

Every keyword, with its values, unit and default, is listed in
:ref:`keyword_reference`.

Validation
----------

FEPX checks the whole file before a simulation starts. Every problem is a fatal error that
gives the file name and line (``simulation.cfg:12: ...``) and names the keywords involved;
FEPX then stops with a non-zero exit status. In particular:

- **Syntax:** an unknown keyword; a keyword reserved for a later version (``... is not
  available in this version of FEPX``); a keyword
  given twice (except the per-component rows, one per keyword and component, and the
  repeatable ``print`` and ``set_bc`` lines); a wrong number of values; a value of the wrong
  type, out of range or not among the allowed words.

- **Structure:** a per-phase keyword outside a phase block; a phase index outside
  1..\ ``number_of_phases``, or a phase defined twice or not at all; elastic constants that do
  not match the crystal type, or that are not positive definite.

- **Steps:** a per-step row without exactly ``number_of_steps`` values; a step without its
  target or rate; a rate that moves away from its target; not exactly one of
  ``number_of_incrs``, ``dtime`` and ``dstrain`` in a step.

- **Boundary conditions:** ``set_bc`` lines with a preset; rate rows with ``general``; a
  uniaxial preset with rows on several components; a node set that the mesh does not have; a
  degree of freedom given two velocities; a ``set_mpc`` without ``vel``, a zero coefficient, a
  node set and direction given two MPCs, MPCs with ``triaxial`` or with load control, a node
  that two different MPCs make a secondary, a cycle of MPCs, velocity conditions that an MPC
  group cannot satisfy.

- **Output:** an unknown variable, or one given twice.

The commands that read the orientations, phases and per-element values from separate files, and
the restart of a simulation, are described in :ref:`special_commands`.
