.. _file_formats:

File and Directory Formats
==========================

.. _simulation_directory:

Simulation Directory (:file:`simulation.sim`)
---------------------------------------------

The results are written in a directory, :file:`simulation.sim`, which Neper (5.0.1 or
later) reads:

.. code-block:: none

  simulation.sim
  |-- inputs
  |   |-- simulation.cfg
  |   |-- simulation.msh
  |   `-- (simulation.ori, simulation.phase, simulation.opt, when read;
  |        simulation.tess, when present)
  |-- restart                (with print restart)
  |   |-- restart.rank0
  |   `-- ...
  `-- results
      |-- nodes
      |   |-- coo
      |   |   |-- coo.step0
      |   |   |-- coo.step1
      |   |   `-- ...
      |   `-- ...
      |-- elts
      |   |-- stress
      |   |   |-- stress.step0
      |   |   `-- ...
      |   `-- ...
      |-- forces
      |   |-- x0
      |   |-- x1
      |   `-- ...
      |-- convergence        (with print convergence)
      `-- macro              (with periodic conditions)

where

- :file:`inputs` holds copies of the input files, and of the tessellation
  :file:`simulation.tess` when it is in the working directory (FEPX does not read it, but Neper
  uses it for post-processing and visualization);
- :file:`restart` holds the restart point, one file per process (:ref:`restart`);
- :file:`results/nodes/<var>` and :file:`results/elts/<var>` hold one file per printed step,
  :file:`<var>.step<n>`, from step 0 (the initial state); each line is one node or element, in
  the order of their ids, with the components on successive columns (see
  :ref:`simulation_results`);
- :file:`results/forces/<faset>` holds the force table of each face set, one row per
  increment;
- :file:`results/convergence` and :file:`results/macro` hold the convergence table and the
  macroscopic table (:ref:`simulation_results`);
- Neper's post-processing (``neper -S``) adds its own results, e.g. :file:`results/mesh`.

.. _mesh_file:

Mesh File (:file:`simulation.msh`)
----------------------------------

The mesh file is Neper's msh format (Gmsh's msh format version 2.2, extended by Neper), as
written by Neper. FEPX reads **only mesh version 2.4** (``$MeshVersion 2.4``), written by Neper 5.0.1 and later. The sections FEPX reads:

.. code-block:: none

  $MeshFormat
  2.2 0 8
  $EndMeshFormat
  $MeshVersion
  2.4
  $EndMeshVersion
  $Nodes
  <number_of_nodes>
  <node_id> <node_x> <node_y> <node_z>
  ...
  $EndNodes
  $Elements
  <number_of_elements>
  <elt_id> <elt_type> 3 <elset_id> <elset_id> <partition> <elt_node1> ...
  ...
  $EndElements
  $NSets
  <number_of_nsets>
  <nset_label>
  <nset_node_nb>
  <nset_node1>
  ...
  $EndNSets
  $Fasets
  <number_of_fasets>
  <faset_label>
  <faset_elt_nb>
  <faset_elt_id> <faset_elt_node1> ...
  ...
  $EndFasets
  $NodePartitions
  <number_of_nodes>
  <node_id> <node_partition>
  ...
  $EndNodePartitions
  $ElsetOrientations
  <number_of_elsets> <orientation_descriptor>:<convention>
  <elset_id> <ori_param1> ...
  ...
  $EndElsetOrientations
  $ElementOrientations
  <number_of_elements> <orientation_descriptor>:<convention>
  <element_id> <ori_param1> ...
  ...
  $EndElementOrientations
  $Groups
  elset
  <number_of_elsets>
  <elset_id> <phase>
  ...
  $EndGroups

and also ``$PhysicalNames``, ``$Periodicity`` (for periodic conditions, :ref:`periodic_conditions`:
lines ``<secondary> <primary> <P_x> <P_y> <P_z>``, with :math:`\mathbf x_s = \mathbf x_p + \mathbf P
\circ \mathbf l`) and the value sections ``$ElsetCrss``, ``$ElsetCrssSat``, ``$ElsetHardCoeff``,
``$ElementCrss``, ``$ElementCrssSat``, ``$ElementHardCoeff``, ``$ElsetRho0``, ``$ElementRho0``,
``$ElsetIrrLoopDensity``, ``$ElementIrrLoopDensity``, ``$ElsetIrrLoopDiameter``,
``$ElementIrrLoopDiameter`` (see
:ref:`opt_file`). ``$Domain``, ``$Topology``,
``$ElsetCrySym`` and ``$Comments`` are accepted and not used. Any other section is skipped with
a warning. Every section is checked strictly (counts, ids in range, footer); a section given
twice is an error; the order of the sections does not matter.

- **Elements.** The 3D elements are all 10-node tetrahedra (Gmsh type 11, the default) or all
  20-node hexahedra (type 17); lower-dimensional elements (Gmsh types 15, 1, 8, 2, 3, 9, 16,
  10) are skipped. The element nodes are in Gmsh's order. The three tags are the element set
  (twice) and the partition (non-zero for a partitioned mesh).

- **Element sets** (grains) are numbered in increasing order of their ids.

- **Face sets** (``$Fasets``) give, for each boundary face, its adjacent 3D element and the
  face's nodes, in Gmsh's order. Faces are assumed to point out of the body, as Neper writes
  them (on the cut surfaces of a domain with ``-transform cut``, from Neper 5.0.1).

- **Orientations.** ``<orientation_descriptor>`` is ``rodrigues``, ``euler-bunge``,
  ``euler-kocks``, ``axis-angle`` or ``quaternion`` (see :ref:`ori_file`), and
  ``<convention>`` is ``active`` or ``passive``. Element orientations, when given, are used
  instead of element-set orientations.

- **Phases.** ``$Groups`` gives the phase of each element set (all phase 1 if absent).

- **Partitions.** The partition tags and ``$NodePartitions`` define the distribution of the
  mesh over the processes (see :ref:`running`).

Re-partition a mesh with ``neper -M -loadmesh simulation.msh -part <N> -o simulation`` (see
:ref:`running`).

.. _ori_file:

Orientation File (:file:`simulation.ori`)
-----------------------------------------

With ``read_ori_from_file``, the orientations of the mesh file are replaced by those of
:file:`simulation.ori`, which has the format of the mesh file's orientation sections, per
element set or per element::

    $ElsetOrientations
    <number_of_elsets> <orientation_descriptor>:<convention>
    <elset_id> <ori_param1> ...
    ...
    $EndElsetOrientations

    $ElementOrientations
    <number_of_elements> <orientation_descriptor>:<convention>
    <element_id> <ori_param1> ...
    ...
    $EndElementOrientations

The file may also hold the ``$ElsetCrySym`` section that Neper writes at its top (e.g.
``triclinic``); it must then be the same as the mesh file's.

where ``<number_of_elsets>`` (``<number_of_elements>``) is the number of orientations of the
section, one per element set (element), ``<orientation_descriptor>`` is the parameterization of
the orientations, ``<convention>`` is the orientation convention, ``<elset_id>``
(``<element_id>``) is the 1-indexed identifier of the element set (element), and
``<ori_param1> ...`` are the parameters of its orientation.

``<convention>`` is ``active`` or ``passive``. An *active* orientation describes the basis
transformation from the crystal frame to the sample frame (crystal-to-sample); a *passive* one
that from the sample frame to the crystal frame (sample-to-crystal).

``<orientation_descriptor>`` is one of the following, each with its per-line format:

- ``rodrigues``: each orientation is described by :math:`r_1, r_2, r_3`, with
  :math:`\mathbf{r} = \mathbf{t}\tan(\omega/2)` (:math:`\mathbf{t}` and :math:`\omega` as in
  ``axis-angle``). The per-line format is::

    <id> <r_1> <r_2> <r_3>

- ``euler-bunge``: each orientation is described by :math:`\varphi_1, \Phi, \varphi_2`, the
  rotations about the :math:`z`, :math:`x'` and :math:`z''` axes, in degrees. The per-line
  format is::

    <id> <phi_1> <Phi> <phi_2>

- ``euler-kocks``: each orientation is described by :math:`\Psi, \Theta, \phi`, the rotations
  about the :math:`z`, :math:`y'` and :math:`z''` axes, in degrees. The per-line format is::

    <id> <Psi> <Theta> <phi>

- ``axis-angle``: each orientation is described by :math:`t_1, t_2, t_3, \omega`, with
  :math:`\mathbf{t}` the unit rotation axis and :math:`\omega` the rotation angle about it, in
  degrees. The per-line format is::

    <id> <t_1> <t_2> <t_3> <omega>

- ``quaternion``: each orientation is described by :math:`q_0, q_1, q_2, q_3`, with
  :math:`q_0 = \cos(\omega/2)` and :math:`q_i = t_i \sin(\omega/2)`, :math:`i = 1, 2, 3`
  (:math:`\mathbf{t}` and :math:`\omega` as in ``axis-angle``). The per-line format is::

    <id> <q_0> <q_1> <q_2> <q_3>

where ``<id>`` is the element set's (element's) identifier. For example, two element sets with
Bunge Euler angles in the active convention::

    $ElsetOrientations
    2 euler-bunge:active
    1 30.0 45.0 60.0
    2 0.0 0.0 0.0
    $EndElsetOrientations

The output ``ori`` uses the descriptor and convention of the input.

.. _phase_file:

Phase File (:file:`simulation.phase`)
-------------------------------------

With ``read_phase_from_file``, the phases of the mesh file are replaced by those of
:file:`simulation.phase`::

    $Groups
    elset
    <number_of_elsets>
    <elset_id> <phase>
    ...
    $EndGroups

.. _opt_file:

Optional Inputs File (:file:`simulation.opt`)
---------------------------------------------

With ``read_from_opt_file``, per-element-set or per-element values of the
crystal-plasticity phases are read from :file:`simulation.opt` and replace the configuration's
values (:ref:`per_element_values`)::

    $Elset<variable>
    <number_of_elsets> <number_of_values>
    <elset_id> <value_1> ...
    ...
    $EndElset<variable>

    $Element<variable>
    <number_of_elements> <number_of_values>
    <element_id> <value_1> ...
    ...
    $EndElement<variable>

where ``<variable>`` is ``Crss`` (initial slip strengths ``tau_0``), ``CrssSat`` (saturation
strengths ``tau_s``, or :math:`\tau_{s0}` with ``saturation_evolution on``) or ``HardCoeff``
(hardening coefficients ``h_0``). Each line gives one value, one per active slip family or one
per slip system of the entity's phase (``CrssSat`` and ``HardCoeff`` one value with
``saturation_isotropic``, ``CrssSat`` one value with saturation evolution); the second header
value is not used. For each variable give either the per-element-set or the per-element
section, not both, with a line for every element set or every element of the mesh. The file
describes the whole polycrystal: values for element sets or elements of an elastic phase are
ignored, with a warning. Each element must keep ``tau_s`` above ``tau_0`` (system 1 with
``saturation_isotropic``) and ``h_0`` non-negative.

For phases with ``hardening kocks_mecking`` the only variable read is ``Rho0``,
the initial network dislocation density :math:`\rho_0` [1/mm²] (``km_rho_0``), one positive
value per element set or element (``$ElsetRho0`` or ``$ElementRho0``, not both, covering every
element set or element). ``Crss``, ``CrssSat`` and ``HardCoeff`` values are ignored for those
phases, and ``Rho0`` values for the other phases, each with a warning.

For phases with ``irradiation on`` the variables ``IrrLoopDensity`` (the loop
number density :math:`\rho_l` [1/mm³], ``irr_loop_density``, >= 0) and
``IrrLoopDiameter`` (the loop diameter :math:`d_l` [mm], ``irr_loop_diameter``, > 0)
are read, one value per element set or element each (``$ElsetIrrLoopDensity`` or
``$ElementIrrLoopDensity``, ``$ElsetIrrLoopDiameter`` or ``$ElementIrrLoopDiameter``), with
``Rho0`` and independently of it; a variable not in the file takes the configuration value for
every element. Values for elements of phases without irradiation are ignored, with a warning.
