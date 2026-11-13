.. _messages:

Error and Warning Messages
==========================

FEPX writes its messages to the terminal (standard output, and the errors also to standard
error):

.. code-block:: none

  Error  : <routine>: <message>
  Warning:   - <message>

An **error** ends the run on every process with exit status 1; a run that stops within a step
first writes the results of its last increment as the current step (:ref:`running`). A
**warning** does not stop the run. In the lists below, ``<...>`` stands for the part of a
message that names the case (a file, a line, a keyword, a value, a number). The messages of the
configuration file start with the file and the line, ``simulation.cfg:<line>:`` (or
``simulation.cfg:`` for a missing keyword). Messages that are much alike are listed once.

A message that starts with ``internal error (please report it):`` reports a failed consistency
check of the program itself, which no input should cause: please report it to the developers,
with the input files.

.. contents::
   :local:
   :depth: 1

Command Line
------------

``unknown option "<option>" (see fepx --help)``
   An argument of ``fepx`` that is not one of its options. Run ``fepx --help`` for the list.

Configuration File: Syntax and Values
-------------------------------------

These messages come from the reading of :file:`simulation.cfg` (:ref:`config_file`); the
:ref:`keyword_reference` gives the values of every keyword.

``unknown keyword "<keyword>"``
   A word at the start of a line that is not a keyword of FEPX (for example, a
   misspelling).

``<keyword> is not available in this version of FEPX``
   A keyword, or a value of ``print``, that is reserved for a later version.

``<keyword> must follow a "phase" line``
   A material keyword before the first ``phase`` line: give it inside a phase block.

``<keyword> given twice (lines <l1> and <l2>)``
   A keyword that may be given once, given twice (in the same phase block for the material
   keywords). Remove one of the lines.

``phase <i> defined twice (lines <l1> and <l2>)``
   Two ``phase`` blocks with the same number.

``<keyword>: per-step row without number_of_steps``
   A per-step row (``strain_rate``, ``dtime``, ...) and no ``number_of_steps``.

``<keyword>: <n> values for <m> steps``
   A per-step row with another number of values than ``number_of_steps`` (``-`` counts as a
   value).

``<keyword>: expected a component (11, 22 or 33) and one value per step``
   A per-component row without its component.

``<keyword>: invalid component "<c>" (expected 11, 22, 33, or with periodic conditions 23, 13, 12, 32, 31, 21)``
   A component that is not one of those listed.

``<keyword>: expected <n> value(s), found <m>``, ``<keyword>: expected <n> value(s)``, ``<keyword>: value(s) missing``, ``<keyword> takes no value``
   A line with the wrong number of values (``read_ori_from_file`` and the other flags take
   none).

``<keyword>: "-" is not allowed (every step needs a value)``
   ``-`` in a per-step row that needs a value in every step (``print_results``).

``<keyword>: invalid value "<value>" (expected one of: <values>)``
   A word that is not one of the allowed values.

``<keyword>: "<value>" is not an integer``, ``<keyword>: "<value>" is not a number``
   A value that cannot be read as an integer or a real number.

``<keyword>: value <value> out of range (expected <range>)``
   A value outside the range of the keyword (e.g. a tolerance that must be greater than 0).

``value "<value>" is longer than <n> characters``
   A value too long to be stored (usually a missing blank between two values).

Configuration File: Option Combinations
---------------------------------------

The messages of the form ``<option> is used only with <condition>``,
``<option> cannot be used with <condition>`` (errors) and ``<option> is ignored with
<condition>`` (a warning) are those of the rules listed in :ref:`compatibility` and at the
entries of the options they concern. Give the option only with its condition, or remove it.

Configuration File: Material
----------------------------

``simulation.cfg: <keyword> is required``, ``simulation.cfg: phase <i> (line <l>): <keyword> is required``
   A required keyword is missing (``number_of_phases``, ``boundary_conditions``,
   ``number_of_steps``; in a phase block ``crystal_type`` and ``model``).

``phase <i> is outside 1..<n> (number_of_phases)``, ``simulation.cfg: phase <i> is not defined (number_of_phases <n>)``
   The phase blocks must be numbered 1 to ``number_of_phases``, each once.

``simulation.cfg: phase <i> (<type>, line <l>): <constant> is required``, ``<constant> is not a constant of crystal_type <type>``
   The elastic constants of a phase are exactly those of its crystal type (:ref:`aniso_elasticity`).

``simulation.cfg: phase <i> (line <l>): the elastic constants are not positive definite``
   The stiffness of the given constants is not positive definite (e.g. :math:`C_{12} \ge
   C_{11}`, or a sign or unit error): check the constants and their units (MPa).

``simulation.cfg: phase <i> (<type>, line <l>): c_over_a is required with model crystal_plasticity``, ``<keyword> is required with model crystal_plasticity``
   A crystal-plasticity phase needs ``m``, ``gammadot_0`` and ``tau_0``, and ``c_over_a`` for
   ``hcp`` and ``bct``.

``slip_families: "<family>" is not a slip family of crystal_type <type> (expected: <families>, or all)``, ``slip_families: <family> given twice``, ``slip_families: all cannot be combined with other families``
   The slip families are written exactly as in :ref:`crystal`, each once, or ``all`` alone.

``<keyword>: <n> value(s); give one per active slip family (<f>) or one per slip system (<s>)``, ``<keyword>: <n> value(s); give one value, one per active slip family (<f>) or one per slip system (<s>)``
   A per-family parameter with another number of values than the accepted counts (see
   :ref:`crystal_plasticity`).

``<keyword> is required with hardening <law>``, ``<keyword>: hardening saturation_isotropic takes one value (one hardening increment for all systems, set by slip system 1); per family or per system: saturation_anisotropic``
   The parameters of the hardening law: ``h_0``, ``tau_s`` and ``n`` for the saturation laws,
   one value each with ``saturation_isotropic``.

``tau_s must be greater than tau_0 of slip system 1 (<values>)``, ``tau_s must be greater than tau_0 for every slip system (system <s>: tau_s <value>, tau_0 <value>)``
   The saturation strength must exceed the initial strength (the law would harden away from
   saturation).

``interaction_matrix: more than <n> rows (one per slip system)``, ``interaction_matrix: row <i> has <n> values; give <m> (one per slip system)``, ``interaction_matrix: "<value>" is not a number``, ``interaction_matrix: entry (<i>, <j>) is negative; entries must be >= 0``, ``interaction_matrix: <n> rows; give <m> (one per slip system)``
   The interaction matrix is :math:`N \times N` for the :math:`N` active slip systems, one row
   per line, entries :math:`\ge 0`.

``<keyword> is required with <modifier> on``, ``gammadot_s0 must be positive``, ``m_prime must be positive (saturation_evolution off for a constant saturation strength)``, ``tau_s: saturation_evolution on takes one value (tau_s0, one saturation strength for all systems)``
   The parameters of the saturation-strength evolution.

``cyclic_a must lie between 0 and 1 (the critical accumulated slip a (crss/tau_s)^c must stay below 1)``, ``cyclic_c must be positive``
   The parameters of cyclic hardening.

``<keyword> is required with precipitate_strengthening on``, ``simulation.cfg: precipitate_strengthening on needs the cutting parameters (a_p, b_p), the bowing parameters (c_p, b_m), or both, each set complete``, ``f_p is a volume fraction (at most 1)``, ``precipitate bowing needs a spacing r_p sqrt(pi/f_p) larger than 2 r_p (f_p < pi/4)``
   The parameters of precipitate strengthening: ``f_p``, ``r_p``, and the cutting set, the
   bowing set or both.

``simulation.cfg: tau_s must be greater than tau_0 + the precipitate strength (system 1: <values>)``, ``simulation.cfg: tau_s must be greater than tau_0 + the precipitate strength for every slip system (system <s>: <values>)``
   With precipitates, the saturation strength must exceed the initial strength plus the
   precipitate strength.

``<keyword> is required with hardening kocks_mecking``, ``<keyword> is required with irradiation on``
   The parameters of Kocks--Mecking hardening and of irradiation softening, all required.

``irr_loop_normals: "<family>" is not a loop-normal family (expected: <111>, <100>, <110>)``, ``irr_loop_normals: <family> given twice``
   The loop-normal families are ``<111>``, ``<100>`` and ``<110>``, each once.

Configuration File: Boundary Conditions and Steps
-------------------------------------------------

``simulation.cfg: boundary_conditions general needs set_bc vel lines``
   With ``general`` conditions the velocities are the ``set_bc`` lines: give at least one.

``set_bc: invalid kind "<kind>" (expected vel)``, ``set_bc vel: expected <nset> and one to three <x|y|z> <velocity> pairs``, ``set_bc vel: invalid direction "<d>" (expected x, y or z)``, ``set_bc vel: "<value>" is not a number``, ``set_bc vel: <nset> <d> given twice (lines <l1> and <l2>)``
   The syntax of ``set_bc`` (:ref:`boundary_conditions`): one velocity per node set and
   direction.

``set_mpc needs vel: set_mpc vel <nset> <x|y|z> <coefficient> <offset> [...]``, ``set_mpc1: invalid kind "<kind>" (expected vel)``, ``set_mpc1 vel: expected <nset> and one to three directions <x|y|z>``, ``set_mpc vel: expected <nset> and one to three <x|y|z> <coefficient> <offset> triples``, ``set_mpc1 vel: invalid direction "<d>" (expected x, y or z)``, ``set_mpc vel: "<value>" is not a number``, ``set_mpc vel: the coefficient must be nonzero``, ``set_mpc vel: <nset> <d> given twice (lines <l1> and <l2>)``
   The syntax of the multi-point constraints (:ref:`multi_point_constraints`).

``<keyword> <c>: <preset> loads one component, and line <l> gives <c'> (several components: boundary_conditions triaxial)``
   A uniaxial preset has its rows on one component, the loading axis.

``simulation.cfg: boundary_conditions triaxial needs control_component <11, 22 or 33> (the component that sets the step durations)``, ``control_component <c> has no strain_rate or target_strain rows``, ``control_component <c> is not the loaded component <c'> of <preset>``
   The control component (:ref:`deformation_history`): required with ``triaxial``, a component
   with rows, the loaded component of a uniaxial preset.

``simulation.cfg: boundary_conditions periodic needs control_component <ij> (an imposed component, whose target_strain or the time ends each step)``, ``control_component <c> has no strain_rate row (the control component is imposed)``, ``target_strain <c>: with boundary_conditions periodic only the control component <c'> has strain targets``, ``strain_rate <c>: step <i> has "-": with periodic conditions a component is imposed (a value in every step) or free (no row)``
   The rows of periodic conditions (:ref:`periodic_conditions`).

``simulation.cfg: step <i>, component <c>: two targets (<row>, line <l>, and <row>, line <l>); give one, "-" for the other``, ``simulation.cfg: step <i>, component <c>: no target (<rows>, target_time or target_duration)``
   The control component has exactly one target per step.

``simulation.cfg: step <i>, component <c>: no rate (boundary_conditions <preset> needs a strain_rate row, e.g. strain_rate 33 1.0e-3)``, ``simulation.cfg: step <i>, component <c>: two rates (<row>, line <l>, and <row>, line <l>); give one, "-" for the other``, ``simulation.cfg: step <i>, component <c>: no rate (<row>, line <l>, has "-")``, ``simulation.cfg: step <i>, component <c>: no rate (boundary_conditions <preset> needs a strain_rate <c> row)``
   With a preset, the control component has exactly one rate per step.

``simulation.cfg: step <i>, component <c>: <rate> 0 (line <l>) never reaches <target>: a rate 0 holds the value (a dwell); give target_time or target_duration``, ``simulation.cfg: step 1, component <c>: a dwell (<rate> 0, line <l>) cannot be the first step (the load starts at 0); load first``
   A dwell (a load or stress rate 0) holds the value: it ends at a time, and not in step 1.

``simulation.cfg: step <i>, component <c>: <row> <c> (line <l>) and <row> <c> (line <l>) both given; a component other than the control component <c'> takes one, "-" for the others``, ``simulation.cfg: step <i>, component <c>: no strain_rate, target_strain, target_load, target_stress, load_rate or stress_rate (component <c> has rows: give one of them in every step)``
   A component other than the control component has exactly one row in every step
   (:ref:`deformation_history`).

``simulation.cfg: step <i>, component <c>: <row> <c> with strain_rate <c'> 0: the load-controlled face velocities are found by perturbations of min_pert_frac times the control face velocity, which must not be 0``
   Triaxial load control needs a moving control face: give the control component a nonzero
   strain rate.

``simulation.cfg: step <i>, component <c>: target_strain <c> (line <l>) in a step that ends at <target>: the duration of the step is not known; give strain_rate <c>``, ``simulation.cfg: step <i>, component <c>: <row> <c> (line <l>) in a step that ends at <target> at a <rate>: the duration of the step is not known; give load_rate <c> or stress_rate <c>``
   In a step whose duration is not known (a load or stress target, or a strain target at a load
   or stress rate), the other components are given by rates.

``step <i>: <row> (line <l>) and <row> both give the increments; give one, "-" for the other``, ``simulation.cfg: step <i>: no increments (number_of_incrs, dtime or dstrain)``
   Each step has exactly one of ``number_of_incrs``, ``dtime`` and ``dstrain``.

``simulation.cfg: step <i>: a <target> needs dtime_min, the smallest time increment of the step``, ``step <i>: dtime_min <value> s is larger than dtime <value> s``
   A step whose duration is not known needs ``dtime_min``, at most ``dtime``.

``print <variable> given twice``
   One ``print`` line per variable.

``output_format: give ascii, hdf5, or both once each``, ``output_format hdf5: this FEPX build has no HDF5 output (configure with -DFEPX_ENABLE_HDF5=ON)``
   The output formats (:ref:`output_formats`); the HDF5 output needs a build with it
   (:ref:`installation`).

Before the Run: The Steps
-------------------------

The steps are planned when the run starts, and each step again when it starts (where its start
values are known only then). These messages start with ``simulation.cfg: step <i>...``.

``: boundary_conditions general with a target on component <c>: every set_bc velocity is zero``
   A ``target_strain`` with ``general`` conditions needs a nonzero velocity.

``: dtime_min <value> s is larger than the largest increment <value> s``
   ``dtime_min`` must not exceed the largest increment of the step.

``: the steps need <n> increments, more than max_incr <m>``, ``: <n> increments of <dt> s exceed max_incr <m>``, ``: the steps end at <t> s, after max_total_time <t'> s``, ``: step <i>, component <c>: the strain at the end of the step, <e>, exceeds max_strain <e'>``
   The run limits (:ref:`loading`) would be exceeded: raise the limit, or change the steps.

``: target_strain <c> <value> equals the strain at the start of the step``, ``: target_time <value> s is not after the start of the step (<t> s)``, ``: <target> <c> <value> is already reached at the start of the step (<value>)``
   The target of a step is its start value: the step would have no duration. Targets are
   cumulative (strain, time, load from the start of the run).

``: strain_rate <c> <rate> moves away from target_strain <c> <value> (strain at the start of the step: <value>)``, ``: strain_rate <c> <rate> moves away from <target> <c> <value> (value at the start of the step: <value>; a positive strain rate increases the load)``, ``: <rate> <c> <value> moves away from <target> <c> <value> (value at the start of the step: <value>)``
   The sign of the rate must move the value towards the target.

``: <target> <c> <value> is not more than <tol> <unit> from the initial value 0``, ``: <target> <c> <value> and <value> are not more than twice <tol> (2 x <tol>) apart (a step ends within <tol> of its target)``
   Load and stress targets must be farther apart than their tolerances (:ref:`loading`).

``: dstrain needs a loading component (a target_strain row); use number_of_incrs or dtime``, ``: dstrain with a zero strain rate; use number_of_incrs or dtime``
   ``dstrain`` is a strain increment of the control component.

``: <duration> s is not a whole number of increments of <dt> s; <consequence>.`` (warning)
   ``dtime`` or ``dstrain`` does not divide the step: the last increment is shorter.

Mesh File
---------

These messages start with the file name, :file:`simulation.msh:` (:ref:`mesh_file`).

``: no $MeshFormat section``, ``: no $Nodes section``, ``: no $Elements section``
   The file is not a Gmsh mesh, or it is incomplete.

``: no $MeshVersion section: meshes written before Neper msh file version <v> are not supported (re-save the mesh with a current Neper)``
   An old mesh: write it again with a current version of Neper.

``mesh format must be "2.2 0 8" (Gmsh 2.2, ASCII, 8-byte reals), found "<format>"``, ``expected one mesh version, found "<line>"``, ``mesh file version <v> is not supported: FEPX reads Neper msh file version <v'> only (re-save the mesh with a current Neper)``
   FEPX reads the Gmsh 2.2 ASCII format of Neper's msh files of version 2.4: write the mesh
   with a current Neper (``-format msh``).

``expected a section header ($Name), found "<line>"``, ``unexpected <text>``, ``expected $End<Section>, found "<line>" (wrong count or extra lines in $<Section>)``, ``section $<Section> given twice``
   The structure of the file is broken (sections, their counts and their end lines).

``node ids must be 1, 2, ... in order: expected <id>, found <id'>``, ``element line: expected id, type and number of tags, found "<line>"``, ``element line: expected 3 tags (elset, elset, partition) and <n> nodes, found "<line>"``, ``3D element ids must be consecutive: expected <id>, found <id'>``
   The nodes and the 3D elements are numbered from 1 in order, as Neper writes them.

``element type <type> is not supported (3D elements must be second-order tetrahedra, type 11, or hexahedra, type 17)``, ``mixed 3D element types <type> and <type'> are not supported``, ``no second-order 3D elements (type 11 or 17) in $Elements``
   FEPX uses second-order elements of one type: mesh with ``neper -M ... -order 2``.

``faset face line: expected an element id and 6 or 8 nodes, found "<line>"``, ``periodicity vector components must be -1, 0 or 1, not all 0``, ``physical name line: expected dimension, id and name, found "<line>"``, ``$Groups: the group entity must be "elset", found "<entity>"``
   An invalid line of a face set, of the periodicity relations, of the physical names or of
   the groups.

``$<Section> header: expected count and descriptor:convention, found "<line>"``, ``$<Section> header: expected count and number of values, found "<line>"``, ``negative count``, ``number of values must be positive``, ``$<Section> line: expected an id and 1 to <n> values, found "<line>"``, ``<variable> must not be negative``, ``<what>: expected one word, found "<line>"``
   An invalid section of orientations, groups or values (in the mesh file, :file:`simulation.ori`,
   :file:`simulation.phase` or :file:`simulation.opt`).

``$ElsetCrySym "<symmetry>" differs from the mesh file's "<symmetry'>"``, ``section <section> is not allowed in an auxiliary file (only orientations, groups and value sections)``
   An auxiliary file (:file:`simulation.ori`, :file:`simulation.phase`, :file:`simulation.opt`)
   holds only orientations, groups and values, consistent with the mesh file.

``<what>: expected <n> values, found <m>``, ``<what>: expected integers, found "<line>"``, ``<what>: expected numbers, found "<line>"``, ``<what>: expected an integer, found "<value>"``
   A line of an input file with too few values, or values that are not numbers.

``<size> <n> exceeds the largest integer (<max>): the mesh is too large for the 32-bit indices of FEPX``
   The mesh has more nodes, elements or degrees of freedom than FEPX's integers can count.

``Skipping unknown mesh section <section> of `<file>'`` (warning)
   A section that FEPX does not read; it is ignored.

``unsupported element type (Gmsh type code <code>, supported: 11 and 17)``
   FEPX uses 10-node tetrahedra (type 11) and 20-node hexahedra (type 17): mesh with
   ``neper -M ... -order 2``.

``: $Elements refers to a node id outside 1..<n>``, ``: node <id> belongs to no 3D element``, ``: 3D element with a non-positive elset id``
   An inconsistent mesh (nodes, elements and element sets).

``The mesh has <n> partition(s), but <m> ranks were started: partition it for <m> ranks (neper -M -loadmesh simulation.msh -part <m> -o simulation) or start as many ranks as it has partitions``, ``The mesh is not partitioned, but <m> ranks were started: ...``
   The mesh must have as many partitions as there are processes (:ref:`running`).

``: $NodePartitions given, but the 3D elements are not partitioned``, ``: some 3D elements have a partition and others not``, ``: the 3D elements are partitioned, but there is no $NodePartitions section``, ``: $NodePartitions has <n> lines for <m> nodes``, ``: $NodePartitions: node id <id> out of range``, ``: $NodePartitions: node <id> given twice``, ``: $NodePartitions: partition <p> of node <id> outside 1..<n>``
   An inconsistent partition of the mesh: partition the mesh again with Neper (``neper -M ...
   -part <n>``).

``: nset <name> has a node id outside 1..<n>``, ``: faset <name> has <n>-node faces, which do not match the <type> elements``, ``: faset <name> refers to element id <id>, which is not a 3D element``, ``: faset <name>: node <id> of face <f> is not a node of its element (id <id>)``, ``: $Periodicity refers to a node id outside 1..<n>``
   An inconsistent node set, face set or periodicity relation.

``: <section>: groups (phases) must be positive``, ``: <section> has <n> lines for <m> elsets``, ``: <section>: elset <id> does not exist``, ``: <section>: elset <id> given twice``, ``: <section> has <n> lines for <m> elements``, ``: <section>: element <id> outside 1..<n>``, ``: <section>: element <id> given twice``
   An inconsistent section of phases or orientations, per element set or per element (in the
   mesh file, :file:`simulation.ori` or :file:`simulation.phase`).

``: <section>: invalid orientation descriptor "<descriptor>" (expected descriptor:convention, e.g. rodrigues:active)``, ``: <section>: id <id> has <n> parameters, expected <m>``, ``rotmat orientation is not a rotation matrix (orthogonal, det +1) to 1e-6``, ``axis-angle orientation with a zero axis``, ``zero quaternion orientation``
   An invalid orientation (:ref:`file_formats`).

``the mesh has <n> phase(s), the configuration <m> (number_of_phases)``
   The phases of the mesh (or of :file:`simulation.phase`) and ``number_of_phases`` differ.

``the mesh has no orientations ($ElsetOrientations, $ElementOrientations or simulation.ori)``
   The orientations are in the mesh file or in :file:`simulation.ori` (``read_ori_from_file``).

``the mesh has no face set <axis>1, whose force is the load of the control component <c> (target_load, target_stress, check_necking)``, ``the mesh has no face set <axis>1, whose force is the load of the load-controlled component <c> (target_load, target_stress)``
   Loads are those of the faces ``x1``, ``y1``, ``z1`` of the mesh (Neper writes them with
   ``-faset faces``).

``node set <name> is not in the mesh``, ``line <l>: node set <name> of a multi-point constraint is not in the mesh``
   A node set of ``set_bc`` or of a multi-point constraint that the mesh does not have.

``Line <l>: node set <name> has fewer than two nodes: no multi-point constraint.`` (warning)
   A multi-point constraint on a single node has no effect.

``conflicting velocities prescribed on a node of node set <name> (component <d>)``
   Two velocity conditions on the same degree of freedom (e.g. two node sets that share a node).

``node <id> component <d> is a secondary of two different multi-point constraints (lines <l1> and <l2>)``, ``the multi-point constraints of line <l> form a cycle (a primary node is a secondary of its own secondary)``, ``conflicting velocity conditions in a multi-point constraint group``
   Multi-point constraints that contradict each other (:ref:`multi_point_constraints`).

``Node <id> component <d>: condition <value>, from its multi-point constraint group <g>.`` (warning)
   A velocity condition on one node of a constraint group fixes the whole group.

``boundary_conditions periodic needs a mesh with periodicity relations ($Periodicity; neper -T -periodicity all)``, ``the mesh is not periodic along <axis> ($Periodicity has no relation with vector +-e_<j>): periodic conditions need a mesh periodic along all three axes (neper -T -periodicity all)``, ``node <id> is the secondary of two periodicity relations``, ``the periodicity relations form a cycle (node <id>)``, ``the drive pair of axis <j> has zero length``
   Periodic conditions need a mesh periodic along the three axes (:ref:`periodic_conditions`).

``a velocity condition is on a periodic secondary dof (or on the drive dof of an imposed component), which the periodicity relations constrain: give it on an independent dof (an interior or primary node)``, ``Node <id> component <d>: a velocity condition on a dof that the periodicity relations constrain.`` (warning)
   With periodic conditions, ``set_bc`` conditions go on independent degrees of freedom.

``Periodic conditions: Fbar-dot <ij> and <ji> are both free: the rigid rotation rate in that plane is not fixed (impose one of them, or use small_strain).`` (warning)
   With ``deformation_gradient_rate``, a free shear pair leaves a rotation free.

``initialization isotropic_viscoplastic: the active slip systems of phase <i> span fewer than 5 deviatoric dimensions (at least five independent systems are needed)``
   The yield surface of the phase is not closed (e.g. ``hcp`` with basal and prismatic slip
   only): use ``initialization elastic_step``, or activate more slip families.

``midnode_linearization applies to 10-node tetrahedra only``
   ``midnode_linearization on`` with a mesh of 20-node hexahedra.

Input Files: simulation.opt
---------------------------

``simulation.opt: give $Elset<variable> or $Element<variable>, not both``
   A variable is given per element set or per element (:ref:`opt_file`).

``simulation.opt: <section> of <id> (phase <i>): <n> value(s); give <counts>``, ``simulation.opt: <section> of <id> (phase <i>): <n> value(s); give one value``
   The number of values of a strength or hardening parameter, as in the configuration.

``simulation.opt: <section>: initial strengths must be positive``, ``simulation.opt: <section>: hardening coefficients must not be negative``, ``simulation.opt: <section>: the saturation strength must exceed the initial strength``, ``simulation.opt: <section>: <variable> <value> must be positive``, ``simulation.opt: <section>: <variable> <value> must not be negative``
   A value out of range.

``Initial strengths from simulation.opt are used as given; precipitate_strengthening is ignored for them.`` (warning), ``simulation.opt gives plasticity parameters for phase <i>, which is elastic: ignored.``, ``simulation.opt gives Crss, CrssSat or HardCoeff values for phase <i>, whose hardening is kocks_mecking (only Rho0 and the irradiation loop sections are read): ignored.``, ``simulation.opt gives Rho0 values for phase <i>, whose hardening is not kocks_mecking: ignored.``, ``simulation.opt gives IrrLoopDensity or IrrLoopDiameter values for phase <i>, which has no irradiation: ignored.`` (warnings)
   Values of :file:`simulation.opt` that the phase's model does not use.

During the Run
--------------

``stopped in step <i> at increment <n> (t = <t> s): <reason>``
   The run stopped within a step (the reason follows); the results of the last increment are
   written as the current step.

``the global iteration did not converge in global_max_iters <n> iterations``, ``the global iteration is diverging``
   The increment could not be solved: shorten the increments (``number_of_incrs``, ``dtime``,
   ``dstrain``), or use ``subincrementation on`` (:ref:`solver_parameters`).

``the linear solver (CG) failed: it reached cg_max_iters, or a nodal block of the block-Jacobi preconditioner is not positive definite``
   The linear system could not be solved: raise ``cg_max_iters``, or shorten the increments.

``Material point solves did not converge (<n> at the converged velocities, at most <m> in an iteration); they keep their last iterate.`` (warning)
   Some material points did not converge: shorten the increments, or use
   ``subincrementation on``.

``The residual of the global iteration is not finite.``, ``The velocity change of the global iteration is not finite.``, ``A nodal block of the block-Jacobi preconditioner is not positive definite.``, ``Newton-Raphson diverging: back to SA.``, ``Newton-Raphson converging slowly: back to SA.`` (warnings)
   The global iteration recovers (back to successive approximation); the increment may still
   fail with one of the errors above.

``increment <n>: <reason> after <d> subdivisions (sub-increment dt/2^<d> = <dt> s); refine the time stepping (number_of_incrs, dtime or dstrain) or raise subincrement_max_depth``, ``<n> successive increments needed sub-incrementation (increment <i>; <reason>); refine the time stepping (number_of_incrs, dtime or dstrain) or raise subincrement_max_successive``
   The limits of sub-incrementation were reached.

``step <i>, increment <n>: no landing on the target <value> <unit> within the tolerance <tol> in max_bc_iter <m> attempts (bracket: <...>); the last accepted increment is written; raise max_bc_iter or the tolerance``
   A load or stress target could not be reached within its tolerance (:ref:`loading`).

``step <i>, increment <n>, boundary-velocity iteration <k>: the perturbed solves give a singular system (the face forces do not respond to the face velocities); check the load-controlled components and min_pert_frac``, ``step <i>, increment <n>: the perturbed solves of the boundary-velocity iteration give a singular system (the face forces do not respond to the face velocities); check the load-controlled components and min_pert_frac``
   The load-controlled faces do not respond to their velocities (e.g. two faces of the same
   axis): check the load-controlled components, or raise ``min_pert_frac``.

``<message>; component <c>: <value>, to reach <value>``
   The boundary-velocity iteration did not converge in ``max_bc_iter`` iterations; the message
   gives each load-controlled value and the value it must reach. Raise ``max_bc_iter``, loosen
   ``load_tol`` or ``stress_tol``, or shorten the increments.

``stopped in step <i> at increment <n> (t = <t> s): necking (check_necking on): the load on face <axis>1 fell from <f> N to <f'> N``
   The load fell (``check_necking on``, :ref:`loading`).

``initialization isotropic_viscoplastic: a nodal block of the block-Jacobi preconditioner is not positive definite``, ``initialization isotropic_viscoplastic: the conjugate gradient reached cg_max_iters``, ``initialization isotropic_viscoplastic: the velocity change is not finite``
   The isotropic viscous start could not be solved: use ``initialization elastic_step``.

``The mesh has <n> partitions, but 1 rank was started: the rank takes the whole mesh.`` (warning)
   Nothing to correct: the partitions are not needed on 1 process (:ref:`running`).

Restart
-------

``restart on: no restart file <file> (the run to continue needs print restart)``
   The run to continue wrote no restart point: it needs ``print restart`` (:ref:`restart`).

``restart on: output directory <dir> not found``, ``restart on: <dir>/<file> not found``, ``restart on: cannot open <file>``, ``restart on: cannot reopen <file>``, ``restart on: <file> is shorter than at the restart point``, ``<file>:<dataset> has fewer rows than at the restart point``
   The simulation directory of the run to continue is missing or was changed.

``restart on: the restart files of the ranks are of different restart points (increments <n> to <m>)``, ``restart on: incomplete restart file``
   The restart files were not all written by the same restart point (an interrupted run, or
   files copied from another run).

``restart on: "<line>" of the restart point's configuration is not in simulation.cfg (the material, the boundary conditions, the input files and the output format must not change)``, ``restart on: "<line>" in simulation.cfg is not in the restart point's configuration (the material, the boundary conditions, the input files and the output format must not change)``, ``Restart on: "<line>" of the restart point's configuration is not in simulation.cfg.``, ``Restart on: "<line>" in simulation.cfg is not in the restart point's configuration.`` (warnings for the solver keywords)
   The configuration of a restart keeps the material, the boundary conditions, the input files
   and the output format of the run it continues.

``restart on: the mesh (nodes, coordinates, phases or orientations) or its partition differs from the restart point's``, ``: material state <component> <n> in the restart file, <m> in this run``
   The mesh, or the material, differs from that of the run to continue.

``cannot open restart file <file>: <reason>``, ``<file> is not a FEPX restart file``, ``<file>: incomplete header``, ``<file>: restart format version <v>, this FEPX reads version <v'>``, ``<file>: written with another byte order``, ``<file>: written with other real or integer kinds``
   The restart file cannot be read by this build of FEPX.

``<file>: written on <n> rank(s); a restart needs the same number of ranks (this run: <m>)``, ``<file>: the file of rank <r>, read by rank <r'>``
   A restart runs on the same number of ranks.

``<file>: record <name> expected, end of file found``, ``<file>: record <name> expected, <name'> found``, ``<file>: record <name> has another type``, ``<file>: record <name> has an invalid rank``, ``<file>: record <name> incomplete``, ``<file>: record <name> has rank <n>, expected <m>``, ``<file>: record <name> has shape (<...>), expected (<...>)``, ``<file>: record <name>: values incomplete``
   A damaged restart file, or one of another version of FEPX.

``cannot create <file>: <reason>``, ``cannot rename <file> to <file'>``
   A restart point cannot be written (permissions, disk full).

Output and Files
----------------

``output directory <dir> already exists``
   FEPX does not overwrite a simulation directory: remove or rename :file:`simulation.sim`, or
   continue it (``restart on``).

``cannot create directory <dir>``, ``cannot create <file>``, ``cannot write <file>``, ``cannot remove <file>``, ``cannot read input file <file>``
   The simulation directory cannot be written (permissions, disk full).

``cannot open file `<file>' (status <n>)``, ``<file>: read failed after line <n>: <reason>``, ``<file>: unexpected end of file after line <n> (expected <what>)``, ``<file>:<line>: <message>``
   An input file is missing, cannot be read, or ends early.

``HDF5 file <file> not found``, ``<file>:<dataset> is not a two-dimensional dataset``, ``HDF5 error in h5get_libversion_f``, ``the HDF5 library cannot be initialized (h5open_f)``, ``<file>:<dataset>: HDF5 error in <routine>``
   The HDF5 output cannot be written or read (a damaged file, or an HDF5 installation
   problem).

System
------

``<action> failed (status <n>): <reason>``, ``<action> failed (status <n>)``
   An operation of the system (an allocation, the opening or reading of a file) failed: the
   reason is that of the system (e.g. not enough memory, a missing file).

``<MPI routine> failed (MPI error <n>): <reason>``
   An MPI call failed: see the reason, and the documentation of the MPI library.

``gathered size exceeds the largest default integer``
   An output array is too large for MPI's counts (a mesh of more than about two billion
   values): write fewer variables, or contact the developers.

``internal error (please report it): <check>``
   A failed consistency check of the program: please report it, with the input files.
