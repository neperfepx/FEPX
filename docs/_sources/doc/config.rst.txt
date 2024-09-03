.. _config_file:

Configuration File (:file:`simulation.cfg`)
===========================================

A simulation input typically consists of two files with pre-defined names: the mesh file, :file:`simulation.msh`, and the configuration file, :file:`simulation.cfg`, which are read by the program at run time.  The mesh file describes the polycrystal as a finite element mesh, where the grains are defined by element sets, with additional fields to record polycrystal-related information, such as (optionally) phases, crystal orientations, or simulation-related information, such as node sets and mesh partitions.
The configuration file completes the definition of the polycrystal by the specification of the material model and the corresponding parameters, and defines the details of the simulation (boundary conditions, loading history, output, etc).
The mesh file is typically provided by Neper and does not need to be edited, and information on its format can be found in :ref:`mesh_file`.  Conversely, the configuration file is specific to the program and is to be written by the user.
Of course, the mesh file and the and configuration file are closely related and should be consistent with each other.

General Structure
-----------------

The configuration file (:file:`simulation.cfg`) is structured in several, successive blocks that define different aspects of the simulation.  Each of these blocks is headed by a line starting by :data:`#`, which provides a short description of the block.  The structure for the :file:`simulation.cfg` file is as follows::

  # Material Description

      { parameter definition }

  # Boundary Conditions / Deformation History

      { parameter definition }

  # Output

      { parameter definition }

  # Special commands

     { commands }

:data:`{ parameter definition }` represents a general block of text that follows the following syntax::

      <key_phrase1> <value1>
      <key_phrase2> <value2>
      ...


:data:`\<key_phrase\> \<value\>` is the input structure of the file, where :data:`<key_phrase>` is either a simulation parameter or a command and :data:`<value>` are one or several numerical or literal values.  A parameter takes only one value, while a command can require one or several values. The possible inputs will be defined in the following sections.

:data:`{ commands }` represents *special commands* that go beyond the definition of a standard simulation, by allowing to override mesh attributes (see :ref:`overriding`) or restarting a simulation (see :ref:`restart`).

In the following, brackets (:data:`[` / :data:`]`) will be used to denote optional inputs (they are *not* expected in the :file:`simulation.cfg` file and will generate parsing errors)::

      <required_key_phrase> <value1>
      [ <optional_key_phrase> <value> ]
     ...

Any piece of text that is preceded by a :data:`#` is considered as a comment and is ignored (which also makes block headers optional).
A single line in the file should only ever pertain to a single :data:`<key_phrase>`/:data:`<value>` pairing.  All strings are interpreted literally and should be lowercase except where otherwise stated.
Parameters within a block may generally be specified in any order; however, the material parameters and deformation history must follow specific orders, and it is recommended that the overall structure of the :ref:`config_file` follows the example structure above.  Spaces and indentations may be freely used for file formatting and do not affect parsing.

.. _consistency:

.. note::

  The following points are important to ensure the consistency between the :ref:`mesh_file`, the :ref:`config_file` and the running conditions (MPI):

    - The physical dimension of the mesh: the loading definition may be sensitive to this dimension
    - The number of phases: a material behavior must be defined for each phase
    - The definition of the node sets: the boundary conditions must be applied to existing node sets
    - The number of partitions of the mesh: it must match the number of computation units. A mesh can be "re-partitionned" by Neper at any time using:

      .. code-block:: console

        neper -M -loadmesh simulation.msh -part <num_partitions>

Special Commands
----------------

.. _overriding:

Polycrystal Overwriting
^^^^^^^^^^^^^^^^^^^^^^^

Optionally, crystallographic phase and orientations different from the ones defined in the :ref:`mesh_file` can be defined in external files by adding appropriate commands to the configuration file.

To read in orientations from an :ref:`ori_file`, the following command may be used::

    read_ori_from_file

To read in grain/phase assignments from a :ref:`phase_file`, the following command may be used::

    read_phase_from_file


Config Overwriting
^^^^^^^^^^^^^^^^^^^^^^^

To read in hardening parameters from an :ref:`opt_file`, the following command may be used::

    read_from_opt_file


.. _restart:

Simulation Restart
^^^^^^^^^^^^^^^^^^

A simulation may be restarted only if the restart files were printed as simulation output on the previous run (see :ref:`advanced_results`).  Upon restart, the restart files must be included in the simulation directory along with all simulation inputs (the :ref:`mesh_file`, the :ref:`config_file`, and any other files, such as an :ref:`ori_file` or :ref:`phase_file`).

A simulation may be restarted by adding the following line to the :file:`simulation.cfg` file::

  restart on

When restarting a simulation, the prescribed :ref:`deformation_history` should include only additional steps, as the restarted simulation will not consider steps that were completed in the previous simulation.  Consequently, step and increment indices, as well as the simulation time, are all reset to 0.

A simulation restart must be performed with the same number of cores that were used to run the original simulation.
When a simulation restarts, it will attempt to find the simulation restart files with the highest step index, :data:`<N>`. It will write output variable data to a new set of files, :file:`post.<var>.rst<N+1>.core*`, where :data:`<var>` is the requested output variable name and :data:`rst<N+1>` is the restart label applied to all new output variable files. Likewise, if restart files are again printed, their index will increase to :data:`<N+2>` (the previous restart's files will not be overwritten).
