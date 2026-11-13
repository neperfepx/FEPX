:orphan:

.. _porting:

Porting FEPX v2 Inputs
======================

FEPX v3 does **not** read FEPX v2 inputs unchanged, and it has no compatibility mode: an obsolete
keyword or value is a fatal error (an unknown keyword, or an invalid value), so an old file fails
loudly instead of being silently misread. This page lists what changes and how to port an FEPX v2
simulation.

Mesh File (:file:`simulation.msh`)
----------------------------------

FEPX v3 reads only Neper msh files of **version 2.4** (the ``$MeshVersion`` section). A missing
``$MeshVersion`` or another version is a fatal error. Older files differ in ways that FEPX v3
does not try to detect: their orientation labels ``passive``/``active`` are swapped (versions
up to 2.2), versions before 2.2.1 have wrong 3D element ids in ``$Fasets``, and versions before
2.4 write the nodes of each ``$Fasets`` face in reverse order.

To port a mesh, regenerate it with a current Neper, or re-save it:

.. code-block:: console

  $ neper -M -loadmesh old.msh -o simulation

which writes :file:`simulation.msh` at version 2.4. The orientation label of a re-saved file of
version 2.2 or earlier changes (e.g. ``rodrigues:passive`` becomes ``rodrigues:active`` with the
same numbers): this is Neper correcting the old, swapped label, and is expected.

The orientation (:file:`simulation.ori`) and phase (:file:`simulation.phase`) files keep their
FEPX v2 format (see :ref:`file_formats`).

Configuration File (:file:`simulation.cfg`)
-------------------------------------------

The configuration file keeps FEPX v2's look (``keyword value(s)`` lines, ``#`` comments), but:

- **The parser is strict.** Unknown, misplaced or duplicated keywords, wrong value counts or
  types, out-of-range values and conflicting keywords are fatal errors with the file name and
  line (see :ref:`config_file`).

- **Each phase states its model** (``model elastic`` or ``model crystal_plasticity``).

- **The boundary conditions and the loading are separate.** ``boundary_conditions`` chooses the
  constrained faces (a preset) or ``general`` conditions; the loading is written as per-step
  rows whose first value is the component (see :ref:`loading`).

- **Hexagonal crystals need** ``c33``. FEPX v2 imposed :math:`C_{33} = C_{11} + C_{12} - C_{13}`;
  FEPX v3 uses the full anisotropic stiffness. To reproduce an FEPX v2 simulation, give ``c33``
  equal to ``c11 + c12 - c13``.

- **Slip strengths are** ``tau_0`` and ``tau_s`` (FEPX v2: ``g_0``, ``g_s0``), and the active slip
  families are chosen with ``slip_families`` (FEPX v2: ``g_0_bcc_112``). FEPX v2's last slip system
  never hardened under isotropic hardening (an error of FEPX v2, absent in FEPX v3),
  so results with hardening differ from those of the original FEPX v2.

- **The first increment is preceded by a small elastic step** (``initialization``, see
  :ref:`solver_parameters`), where FEPX v2 started from the velocities of an isotropic viscous solve
  and a stress estimate at a yield-surface vertex; ``initialization isotropic_viscoplastic``
  starts as FEPX v2 does, with vertices computed from the phase's slip systems and strengths (FEPX
  v2 used fixed tables), and for every kind of boundary condition (FEPX v2: not with its triaxial
  conditions), with the multi-point constraints imposed in the isotropic solve (FEPX v2 imposed them
  afterwards).

- **The conjugate gradient solver uses a nodal block-Jacobi preconditioner** without constraints
  (``cg_preconditioner``, see :ref:`solver_parameters`), where FEPX v2 used the diagonal (Jacobi)
  one; ``cg_preconditioner jacobi`` uses FEPX v2's. The solves need fewer iterations; the results
  agree to the solver tolerance. With multi-point constraints or periodic conditions FEPX v3
  uses Jacobi, as FEPX v2.

- **Per-family values are given for every family.** ``m``, ``gammadot_0`` (and later the
  strengths) take one value per active slip family or one per slip system; a single value is
  accepted only when one family is active. FEPX v2 broadcast a single ``m`` to every ``hcp`` or
  ``bct`` family: write it once per family (e.g. ``m 0.01 0.01 0.01`` for ``hcp``).
  ``gammadot_0`` may now also differ between families.

The main renamed and replaced keywords:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - FEPX v2
     - FEPX v3
   * - ``set_bc strainrate 33 <rate> grip``
     - ``boundary_conditions uniaxial_grip`` and ``strain_rate 33 <rate> ...`` (one value per
       step); ``minimal`` and ``symmetry`` likewise
   * - several ``set_bc strainrate <c> <rate> minimal`` lines
     - ``boundary_conditions triaxial``, ``control_component <c>`` and ``strain_rate <c> ...``
       rows
   * - ``target_strain33 <e1> <e2> ...``
     - ``target_strain 33 <e1> <e2> ...``
   * - ``target_load3 ...`` (``def_control_by uniaxial_load_target``)
     - ``target_load 33 ...``, with a signed ``strain_rate 33`` per step (FEPX v2 reversed the
       velocities itself) and ``dtime_min`` in every load step; or ``target_stress 33 ...``
   * - ``load_tol`` (uniaxial load target tolerance)
     - ``load_tol`` [N] (or ``stress_tol`` [MPa] with stress targets), default 0.1, greater
       than 0: FEPX v3 lands within it (an increment past the target by more is redone shorter;
       FEPX v2 ends the step at the first increment that reaches the target, default 0)
   * - ``load_tol_abs``, ``load_tol_rel`` (triaxial iteration tolerances)
     - ``load_tol`` [N] or ``stress_tol`` [MPa]: absolute tolerances only
   * - ``def_control_by triaxial_constant_strain_rate``, ``number_of_csr_load_steps``,
       ``target_csr_load <Fx> <Fy> <Fz> <dt_max> <dt_min> <print>``, ``strain_rate``,
       ``loading_direction``
     - ``boundary_conditions triaxial``, ``control_component <c>``, a signed
       ``strain_rate <c>`` row (FEPX v2 took the sign from its elastic estimate), ``target_load``
       rows for the three components (the other two load-controlled on the load path), and
       ``dtime``, ``dtime_min``, ``print_results`` rows; ``min_pert_frac``, ``max_bc_iter`` as
       in FEPX v2. FEPX v3 starts each step's load path at the values reached (FEPX v2: at the
       previous step's targets), and lands on the control target within ``load_tol`` (FEPX v2
       accepted the first increment that reached it)
   * - ``def_control_by triaxial_constant_load_rate``, ``number_of_clr_load_steps``,
       ``target_clr_load <Fx> <Fy> <Fz> <dt> <print>``, ``load_rate``, ``loading_direction``
     - ``boundary_conditions triaxial``, ``control_component <c>``, a signed
       ``load_rate <c>`` row (its sign towards each target; FEPX v2's single rate applied to the
       control component in both directions), ``target_load`` rows for the three components (the
       other two load-controlled in time), ``number_of_incrs`` or ``dtime`` and
       ``print_results`` rows; ``max_bc_iter``, ``min_pert_frac`` (now also the lower bound of
       the CLR perturbations, which FEPX v2 left unbounded). FEPX v2 cut a step into
       ``max(nint(T/dt), 1)`` equal increments: write ``number_of_incrs`` for the same
       increments (FEPX v3's ``dtime`` adds a shorter last increment when ``dt`` does not divide the
       step, whose duration starts from the load reached, not the previous target); FEPX v2
       classified loading and unloading by the signed load (loading in compression was unloading),
       FEPX v3 by its magnitude
   * - ``load_rate_jump <step> <rate>``, ``number_of_load_rate_jumps``
     - a ``load_rate`` row whose value changes at that step
   * - ``dwell_episode <step> <time> <dt> <print>``, ``number_of_dwell_episodes``,
       ``max_strain_incr``
     - a step after load step ``<step>`` with ``load_rate <c> 0``, ``target_duration <time>``,
       ``dtime <dt>`` and its ``print_results``; ``max_strain_incr`` as in FEPX v2. The dwell ends
       at its time exactly (FEPX v2 stopped when less than 0.001 ``dt`` remained), its strain bound
       uses the original length (FEPX v2: the current one)
   * - ``set_mpc <nset> <dir> <coeff> <offset> ...``
     - ``set_mpc vel <nset> <dir> <coeff> <offset> ...`` (``vel`` added, as ``set_mpc1``); the
       primary is still the first node of the node set. A node that two lines make a secondary
       is an error unless both relations are the same (FEPX v2: the later line won); offsets stay
       constant at step changes (FEPX v2 scaled them with the velocity)
   * - ``set_bc strainrate <ij> [<value>]`` (periodic, no method word),
       ``imposed_strain_rate_state``
     - ``boundary_conditions periodic`` (``periodic_kinematics small_strain``, the default),
       ``control_component <ij>``, a ``strain_rate <ij>`` row per imposed component (free: no row;
       FEPX v2: a line without value), ``target_strain`` on the control component. FEPX v2 ended a
       step at the strain of the largest rate whatever component the target names: choose the
       control component and target that give the same time
   * - ``set_bc velgrad <ij> [<value>]``
     - the same with ``periodic_kinematics deformation_gradient_rate`` (the nine components)
   * - ``pre_cond 1``
     - ``periodic_preconditioner lumped_diagonal`` (the default); ``pre_cond 2`` and ``3`` are not
       available (``exact_diagonal`` is the exact Jacobi diagonal)
   * - ``loading_direction``, ``def_control_by``
     - removed: inferred from the rows
   * - ``loading_face``, ``number_of_bcs``, ``number_of_mpc``
     - removed (FEPX v2 ignored them)
   * - ``out_format``
     - ``output_format``
   * - ``boundary_conditions user_defined``
     - ``boundary_conditions general``
   * - ``number_of_strain_steps`` (and the other step counts)
     - ``number_of_steps``
   * - ``strain_rate_jump``, ``number_of_strain_rate_jumps``
     - a ``strain_rate`` row whose value changes between steps
   * - ``g_0``, ``g_s0`` (``g_s``)
     - ``tau_0``, ``tau_s``
   * - ``g_0_bcc_112``
     - ``slip_families`` and ``tau_0``
   * - ``hard_type``
     - ``hardening``
   * - ``hardening saturation,isotropic`` (FEPX v2's default)
     - ``hardening saturation_isotropic`` (the default); ``h_0``, ``tau_s`` (``g_s0``) and
       ``n`` take one value each
   * - ``hardening anisotropic``, ``interaction_matrix_parameters``
     - ``hardening saturation_anisotropic``; ``h_0``, ``tau_s`` and ``n`` per family or per
       system; ``interaction_matrix`` (one line per row; default 1 on the diagonal, 1.4 off it).
       FEPX v2's two values ``d h`` are the default form (``d`` on the diagonal, ``h`` off it); its
       longer "coplanar" lists are block matrices, written out row by row: ``fcc`` with five
       values ``d h1 … h4``, 3 × 3 blocks for systems 1–3, …, 10–12 with ``d`` on the diagonal
       and ``h1`` … ``h4`` off it; ``bcc`` with seven values ``d h1 … h6``, 2 × 2 blocks for
       systems (1, 2), …, (11, 12); 0 between blocks
   * - ``hardening saturation_evolution`` (or ``m_prime`` and ``gammadot_s0`` alone)
     - ``saturation_evolution on`` with ``gammadot_s0``, ``m_prime``; ``tau_s`` is
       :math:`\tau_{s0}`
   * - ``hardening cyclic``
     - ``cyclic_hardening on`` (with ``saturation_isotropic``, or new with
       ``saturation_anisotropic``); ``cyclic_a``, ``cyclic_c``
   * - ``hardening precipitation`` (or the ``a_p`` ... ``c_p`` values alone)
     - ``precipitate_strengthening on`` with ``f_p``, ``r_p`` and the cutting (``a_p``,
       ``b_p``) and/or bowing (``c_p``, ``b_m``) set; ``b_p`` is not taken for bowing alone
   * - ``g_0 from_file``, ``g_s0 from_file``, ``h_0 from_file``
     - ``read_from_opt_file``, with the configuration values kept (the file replaces them)
   * - ``midnode_linearization yes|no``
     - ``midnode_linearization on|off`` (every switch is ``on``/``off``)
   * - ``nl_max_iters``, ``nl_tol_strict``, ``nl_tol_min``, ``nl_tol_loose``
     - ``global_max_iters``, ``global_tol``, ``global_tol_abs``, ``global_tol_slow``
   * - ``nr_tol_switch_ref``, ``nr_tol_conv``
     - ``global_nr_switch``, ``global_nr_slow_ratio``
   * - ``sx_tol``, ``sx_max_iters_newton``, ``sx_max_iters_state``, ``max_iter_hard_limit``
     - ``point_tol``, ``point_max_iters_newton``, ``point_max_iters_state``,
       ``point_max_iters_hardening``
   * - ``print strain-el`` (hyphenated output names)
     - ``print strain_el`` (underscores only)
   * - ``suppress <var>``
     - removed: only ``print`` selects outputs

The keywords of FEPX v3 are listed in :ref:`keyword_reference`.

Other changes of behavior:

- ``max_strain`` and ``max_eqstrain`` have no default (no limit) and apply in every loading
  mode (FEPX v2: 0.2, triaxial loading only).
- A ``dtime`` or ``dstrain`` that does not divide a step gives full increments and a shorter
  last one that ends the step exactly on its target, with a warning (FEPX v2 overshot the target).
- Rows that FEPX v2 accepted with fewer values than steps (repeating the last value) must give one
  value per step.

Example
-------

An FEPX v2 file:

.. Not run: an FEPX v2 input, which FEPX does not read.

.. code-block:: none

  number_of_phases 1
  phase 1
  crystal_type HCP
  m 0.010d0
  gammadot_0 1.0d0
  h_0 190.0d0
  g_0 1.0d10 1.0d10 1.0d10
  g_s0 530.0d0
  n 1.0d0
  c11 1.6240e5
  c12 0.9200e5
  c13 0.6900e5
  c44 0.4670e5
  c_over_a 1.587d0

  set_bc strainrate 11 -0.25e-2 minimal
  set_bc strainrate 22 -0.3e-2 minimal
  set_bc strainrate 33 1.0e-2 minimal

  number_of_steps 2
  target_strain33 0.0005 0.0010
  number_of_incrs 2 2
  print_results 1 1

  print stress

and its FEPX v3 equivalent (elastic only: the very large ``g_0`` made the FEPX v2 run elastic;
``c33`` is :math:`C_{11} + C_{12} - C_{13}`, the value FEPX v2 used):

.. literalinclude:: ../tests/regression/cases/elastic_sx_hcp_triaxial/simulation.cfg
   :start-at: ## Material

Output
------

- Results are written unrounded with 17 significant digits (FEPX v2 wrote 7, and set values below
  :math:`10^{-12}` to 0).
- The :file:`.sim` layout and the variable names are FEPX v2's (see :ref:`simulation_results`).
- ``print convergence`` writes :file:`results/convergence` inside the :file:`.sim` directory
  (FEPX v2: :file:`post.conv` in the working directory), with FEPX v2's columns; the initialization
  increment is increment 0.
- Restart (see :ref:`restart`): as in FEPX v2, ``print restart`` writes and ``restart on`` reads,
  and the second configuration lists only the new steps; FEPX v3 writes its restart point at every
  step end (FEPX v2: printed steps only) into :file:`simulation.sim/restart` (FEPX v2: ``rst<N>.*``
  in the working directory), saves the whole state (FEPX v2 lost the accumulated slip of cyclic
  hardening), plans the new steps from where the previous run ended (FEPX v2 assumed monotonic
  loading from zero), and continues :file:`simulation.sim` (FEPX v2 deleted it). FEPX v2's restart
  files cannot be read; ``restart append`` is an error.
- ``strain_pl`` keeps FEPX v2's definition (the plastic strain accumulated in the crystal frame);
  ``strain_pl_sample`` is new: the plastic strain accumulated in the sample frame.
