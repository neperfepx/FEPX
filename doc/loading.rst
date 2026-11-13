.. _loading:

Loading (:file:`simulation.cfg`)
================================

The loading is defined in two independent parts: the **boundary conditions**, which choose the
constrained nodes and degrees of freedom (``boundary_conditions``), and the **deformation
history**, a sequence of *steps* written as per-step rows. Velocities are in mm/s, strain
rates in 1/s, times in s.

.. _terminology:

Terminology
-----------

A simulation is organized in nested levels:

- the **run**: the whole simulation, a sequence of steps;
- a **step**: a part of the deformation history with its own target (a strain or a time) and
  rate, and whose results can be written (``print_results``);
- an **increment**: a time increment of a step (``number_of_incrs``, ``dtime`` or ``dstrain``);
  the velocities, stresses and orientations are solved for at the end of each increment;
- a **global iteration**: one iteration of the solution of an increment for its velocities
  (successive approximation or Newton--Raphson; keywords ``global_...``, see
  :ref:`solver_parameters`); each solves a linear system (``cg_...``);
- a **point solve**: at each global iteration, the solution for the stress at every
  quadrature point (keywords ``point_...``).

With crystal plasticity, the run starts by default with a small, unwritten elastic step that is
part of step 1, or with an isotropic viscoplastic start (``initialization``, see
:ref:`solver_parameters`).

.. _boundary_conditions:

Boundary Conditions
-------------------

::

    boundary_conditions <type>

where ``<type>`` is a preset for cubic (right-prism) domains, ``uniaxial_grip``,
``uniaxial_symmetry``, ``uniaxial_minimal`` or ``triaxial``, or ``general`` for velocity
conditions given node set by node set, or ``periodic`` for periodic conditions on a periodic
mesh (`Periodic conditions`_).

.. confval:: boundary_conditions

   The boundary conditions: a preset for cubic domains (``uniaxial_grip``,
   ``uniaxial_symmetry``, ``uniaxial_minimal``, ``triaxial``), ``general`` (velocities on node
   sets, ``set_bc``) or ``periodic`` (a periodic mesh).

   .. include:: generated/keywords/boundary_conditions.inc


Presets
~~~~~~~

The presets apply to domains whose faces are the node sets ``x0``, ``y0``, ``z0`` (minimum
coordinates) and ``x1``, ``y1``, ``z1`` (maximum coordinates), as Neper writes them; corners
are the node sets ``x<i>y<j>z<k>``. The loading axis :math:`a` of a uniaxial preset is the
component of its loading rows (``11``: :math:`x`, ``22``: :math:`y`, ``33``: :math:`z`);
:math:`b` and :math:`c` are the next axes in cyclic order (for :math:`a = z`: :math:`b = x`,
:math:`c = y`). The loading face :math:`a_1` moves along :math:`a` at the strain rate times the
original length of the domain along :math:`a`.

- ``uniaxial_grip``: face :math:`a_0` fixed (all three components zero); face :math:`a_1`
  moving along :math:`a`, with zero velocity along :math:`b` and :math:`c`.

- ``uniaxial_symmetry``: zero normal velocity on the faces :math:`a_0`, :math:`b_0` and
  :math:`c_0`; face :math:`a_1` moving along :math:`a`. The other faces are free.

- ``uniaxial_minimal``: zero velocity along :math:`a` on face :math:`a_0`; face :math:`a_1`
  moving along :math:`a`; and the rigid-body motions removed by the corner at the origin
  (zero velocity along :math:`b` and :math:`c`) and the corner at coordinate 1 along
  :math:`b` (zero velocity along :math:`c`).

- ``triaxial``: zero normal velocity on the faces ``x0``, ``y0`` and ``z0``; on each face
  ``x1``, ``y1``, ``z1`` of a component with loading rows, a normal velocity equal to that
  component's strain rate times the original length; a face whose component has no rows is
  free. Requires ``control_component`` (below). The control component may end its steps at
  load or stress targets, and the other components may be load-controlled (`Triaxial load
  control`_): their faces then move at the uniform normal velocity that gives the prescribed
  load or stress.

General Conditions
~~~~~~~~~~~~~~~~~~

With ``boundary_conditions general``, velocities are prescribed on node sets of the mesh (or on
single nodes, by node id)::

    set_bc vel <nset> <direction> <velocity> [ <direction> <velocity> ... ]
    ...

where ``<direction>`` is ``x``, ``y`` or ``z``; up to three direction/velocity pairs per line;
any number of lines. A degree of freedom given two velocities is an error. The ``set_bc``
velocities are the loading: no rate rows are given.

.. confval:: set_bc

   A velocity condition of ``boundary_conditions general``: ``vel``, a node set (or a node id),
   then up to three pairs of a direction (``x``, ``y``, ``z``) and a velocity [mm/s]. With
   ``periodic`` conditions, on independent degrees of freedom only (`Periodic conditions`_).

   .. include:: generated/keywords/set_bc.inc


.. _multi_point_constraints:

Multi-Point Constraints
~~~~~~~~~~~~~~~~~~~~~~~

A multi-point constraint (MPC) ties one velocity component of every node of a node set to that of
the set's **first node** (as the mesh file lists it), the *primary*; the other
nodes are its *secondaries*::

    set_mpc1 vel <nset> <direction> [ <direction> [ <direction> ] ]
    set_mpc  vel <nset> <direction> <coefficient> <offset> [ ... ]

``set_mpc1`` makes the component equal on the whole node set (e.g. ``set_mpc1 vel x1 x`` keeps
face ``x1`` planar). ``set_mpc`` gives every secondary :math:`s`
:math:`v_s = c\,v_p + d`, with a nonzero coefficient :math:`c` and an offset :math:`d` [mm/s],
for up to three direction/coefficient/offset triples per line; the offset is constant over the
run. Any number of lines; one relation per node set and direction. A primary may itself be a
secondary of another line (a chain): it is followed to the top primary, with the coefficients
multiplied and the offsets composed.

.. confval:: set_mpc1

   A multi-point constraint that makes velocity components equal on a node set: ``vel``, the
   node set, then one to three directions.

   .. include:: generated/keywords/set_mpc1.inc

.. confval:: set_mpc

   A multi-point constraint :math:`v_s = c\,v_p + d` on a node set: ``vel``, the node set, then
   up to three triples of a direction, a nonzero coefficient :math:`c` and an offset :math:`d`
   [mm/s].

   .. include:: generated/keywords/set_mpc.inc


- A velocity condition (``set_bc vel``, or a preset's) on any node of a group fixes the group:
  every node then moves as the relation gives. Conditions on two nodes of a group that the
  relation cannot both satisfy are an error.
- A node that two lines make a secondary is an error, unless both give the same relation (the
  same primary, coefficient and offset). In Neper's meshes two faces' node sets start with
  their common corner (``x1`` and ``y1`` both with ``x1y1z0``), so ``set_mpc1 vel x1 x`` and
  ``set_mpc1 vel y1 x`` are accepted. A cycle of relations is an error.
- MPCs are available with ``boundary_conditions general`` and the uniaxial presets at strain
  rates; they cannot be used with ``triaxial``, with load or stress targets or rates, or with
  periodic conditions (`Periodic conditions`_).

The equations are solved for the primary and unconstrained degrees of freedom only
(master–slave elimination; :ref:`solver_parameters`); the results are written for every
node, as without MPCs. The step-0 velocities already satisfy the constraints.

.. _periodic_conditions:

Periodic Conditions
~~~~~~~~~~~~~~~~~~~

With ``boundary_conditions periodic`` the domain is a periodic representative volume: the mesh
must have a ``$Periodicity`` section (Neper: ``neper -T ... -periodicity all``, then
``neper -M``), periodic along all three axes. Such meshes are usually tessellations whose grains
cross the domain boundary, so the boundary is not planar; nothing in the periodic conditions
uses faces, node sets or the bounding box. Each relation ties a *secondary* node :math:`s` to a
*primary* node :math:`p` with :math:`\mathbf x_s = \mathbf x_p + \mathbf P \circ \mathbf l`
(:math:`P_j \in \{-1, 0, 1\}`, :math:`\mathbf l` the domain lengths), and imposes

.. math::

   v_k(s) - v_k(p) = \sum_j P_j\, l_j\, \dot{\bar F}_{kj},

where :math:`\dot{\bar{\mathbf F}}` is the rate of the macroscopic deformation gradient (on the
reference periodicity vectors, so that the velocity jumps are constant at constant rate). Each
component of :math:`\dot{\bar{\mathbf F}}` is *imposed* (a ``strain_rate`` row) or *free* (no
row: its macroscopic stress component is zero). The stress that vanishes is the one conjugate to
:math:`\dot{\bar{\mathbf F}}`, the nominal (first Piola–Kirchhoff) stress :math:`\bar{\mathbf
P} = J\,\bar{\boldsymbol\sigma}\,\bar{\mathbf F}^{-\mathsf T}` of the volume-averaged Cauchy
stress written in :file:`results/macro` (with ``small_strain``, :math:`P_{ij} + P_{ji}` of a free
shear pair); the written Cauchy component of a free direction is then not exactly zero under
shear, but of the order of the shear stress times the strain::

    boundary_conditions periodic
    periodic_kinematics small_strain | deformation_gradient_rate     (default small_strain)
    control_component <ij>
    strain_rate <ij> <rate> ... <rate>        (one row per imposed component)
    target_strain <ij> <value> ... <value>    (the control component only)

- ``small_strain``: :math:`\dot{\bar{\mathbf F}}` is
  symmetric; the components are ``11 22 33 23 13 12`` and a shear row imposes :math:`ij` and
  :math:`ji` (one value, the tensor component).
- ``deformation_gradient_rate``: the nine components ``11 12 13 21 22
  23 31 32 33``. A shear pair with both :math:`ij` and :math:`ji` free leaves a rigid rotation
  rate free: allowed with a warning, but the solution is then usually not unique (the
  run may not converge); impose one of the two.
- An imposed component has a value in every step (``-`` is not allowed); the rates may change
  from step to step. ``control_component`` (required) is an imposed component; its own
  ``target_strain`` (a component of :math:`\bar{\mathbf F} - \mathbf I`), or ``target_time``
  or ``target_duration``, ends each step, with the increments of a uniaxial step (``dstrain`` on
  the control component's rate). No load or stress targets or rates; no multi-point constraints;
  no ``check_necking``.
- The macroscopic deformation gradient and its rate are those of the *drive pairs*: for each axis
  :math:`j`, the first relation of the file with vector :math:`\pm\mathbf e_j`; the
  lengths :math:`l_j` are their initial separations. A free component is carried by its drive
  pair's velocities (one unknown per free component; with ``small_strain`` one per free shear
  pair).
- Rigid translation: no node is fixed; the fluctuation :math:`\mathbf v - \dot{\bar{\mathbf
  F}}\mathbf X` has zero mean over the initial volume (with the consistent nodal weights), as in
  the formulation of F. Villette et al. (in preparation).
- ``set_bc vel`` conditions may be given on independent degrees of freedom (an interior or
  primary node, a drive pair's primary node); they are applied as given, not carried through the
  relations. A condition on a constrained degree of freedom (a secondary node, or the
  drive dof of an imposed component) is an error. A velocity component with a condition is not
  shifted to the zero mean.
- ``periodic_preconditioner lumped_diagonal | exact_diagonal`` (default ``lumped_diagonal``): the
  Jacobi diagonal of the reduced system, lumped (the element diagonals summed onto the primary
  degrees of freedom) or the exact one
  (:ref:`solver_parameters`). It changes only the number of conjugate-gradient iterations.
- Output: :file:`results/macro` (:ref:`simulation_results`); the face forces file is not
  written (``print forces`` gives a warning).

.. confval:: periodic_kinematics

   The kinematics of periodic conditions: ``small_strain`` (a symmetric :math:`\dot{\bar{\mathbf
   F}}`, components ``11 22 33 23 13 12``) or ``deformation_gradient_rate`` (the nine
   components).

   .. include:: generated/keywords/periodic_kinematics.inc

.. confval:: periodic_preconditioner

   The Jacobi diagonal of the reduced system with periodic conditions: ``lumped_diagonal`` or
   ``exact_diagonal``. It changes only the number of conjugate-gradient iterations.

   .. include:: generated/keywords/periodic_preconditioner.inc


.. _deformation_history:

Deformation History (Steps)
---------------------------

The simulation is divided into steps; each step reaches a *target* at a given *rate*, in one or
several time *increments*. Each loading row starts with its component (``11``, ``22`` or
``33``), followed by one value per step; ``-`` marks a step where the row is not used::

    number_of_steps <n>
    #                    comp   step 1   ...   step n
    strain_rate          <c>    <rate>   ...   <rate>
    [ load_rate          <c>    <rate>   ...   <rate>   ]
    [ stress_rate        <c>    <rate>   ...   <rate>   ]
    target_strain        <c>    <strain> ...   <strain>
    [ target_load        <c>    <load>   ...   <load>   ]
    [ target_stress      <c>    <stress> ...   <stress> ]
    [ target_time               <time>   ...   <time>   ]
    [ target_duration           <time>   ...   <time>   ]
    number_of_incrs | dtime | dstrain   <value> ... <value>
    [ dtime_min                 <time>   ...   <time>   ]
    [ print_results             <flag>   ...   <flag>   ]

Each step has exactly one target, of any kind (``-`` in the other target rows): the steps of a
run may mix strain, load, stress and time targets.

- ``target_strain <c>``: the engineering strain of component ``<c>`` at the end of the step
  (cumulative).
- ``target_load <c>``, ``target_stress <c>``: a load or stress target (see `Load and stress
  targets`_).
- ``target_time``: the time at the end of the step (cumulative).
- ``target_duration``: the length of the step from its actual start [s] (not with a load or
  stress target of the control component, which ends the step).
- ``strain_rate <c>``: the engineering strain rate of component ``<c>`` during the step. Its
  sign must move the strain towards the target; a rate that changes sign between steps
  unloads.
- ``load_rate <c>``, ``stress_rate <c>``: a constant load [N/s] or stress [MPa/s] rate (see
  `Load and stress rates, dwells`_).
- One of ``number_of_incrs`` (the number of increments), ``dtime`` (the time increment) or
  ``dstrain`` (the strain increment of the control component) per step. With ``dtime`` or
  ``dstrain`` that do not divide the step, the last increment is shorter and ends the step
  exactly on its target, with a warning.
- ``print_results``: ``1`` (default) writes the results at the end of the step, ``0`` does
  not. The initial state (step 0) is always written.

.. confval:: number_of_steps

   The number of steps of the run: every per-step row gives this many values.

   .. include:: generated/keywords/number_of_steps.inc

.. confval:: strain_rate

   The engineering strain rate [1/s] of a component in each step: the component, then one value
   per step (``-`` where the component has another row in the step).

   .. include:: generated/keywords/strain_rate.inc

.. confval:: target_strain

   The engineering strain of a component at the end of each step (cumulative): the component,
   then one value per step.

   .. include:: generated/keywords/target_strain.inc

.. confval:: target_time

   The time at the end of each step [s] (cumulative).

   .. include:: generated/keywords/target_time.inc

.. confval:: target_duration

   The length of each step from its actual start [s].

   .. include:: generated/keywords/target_duration.inc

.. confval:: number_of_incrs

   The number of increments of each step (one of ``number_of_incrs``, ``dtime`` and ``dstrain``
   per step).

   .. include:: generated/keywords/number_of_incrs.inc

.. confval:: dtime

   The time increment of each step [s]; in a step whose duration is not known (a load or stress
   target), the largest increment.

   .. include:: generated/keywords/dtime.inc

.. confval:: dstrain

   The strain increment of the control component in each step; in a step whose duration is not
   known, the largest increment.

   .. include:: generated/keywords/dstrain.inc

.. confval:: print_results

   ``1``: the results are written at the end of the step; ``0``: they are not. One value per
   step.

   .. include:: generated/keywords/print_results.inc


Load and Stress Targets
~~~~~~~~~~~~~~~~~~~~~~~

A step can end at a load instead of a strain. The **load** of component ``<c>`` is the
``<c>`` component of the force on the face ``x1``, ``y1`` or ``z1`` of that axis (as written in
``results/forces``; tension is positive) [N]; the **stress** is that force divided by the
face's current area [MPa]. Load and stress targets of the other components are described in
`Triaxial load control`_. The step runs at its strain rate (with ``general`` conditions: the
``set_bc`` velocities, reversed when the target is below the current value) until the first
increment whose value is within the tolerance :math:`\delta` of the target, :math:`\delta_F` =
``load_tol`` [N] for loads and :math:`\delta_\sigma` = ``stress_tol`` [MPa] for stresses (absolute,
greater than 0, default 0.1 each); its results are those of that increment,
and the log reports the value reached. The step **lands** on its target: an increment that
passes the target by more than :math:`\delta` is discarded (coordinates, velocities and state
restored) and solved again with a shorter time increment, found by false position between the
start of the increment and the attempts made (the time increment of the next attempt is where
the straight line through the bracketing attempts reaches the target), until an attempt ends
within the tolerance; only that attempt counts in the time, the increments and the output.
``max_bc_iter`` (default 10) limits the attempts of an increment; beyond it the run writes the
last accepted increment as the current step and stops with an error. The log reports each
attempt and the number of attempts that landed. Since the step's duration is not known in
advance:

- the increments are given by ``dtime`` (or ``dstrain``), the largest increment, and
  ``dtime_min``, the smallest, which is required in every load or stress target step (``-`` in
  the other steps); ``number_of_incrs`` cannot be used;
- near the target, the next increment is the time the last increment's load rate needs to reach
  the target, times ``dtime_factor`` (default 1.001, at least 1: aim slightly past the target),
  kept between ``dtime_min`` and ``dtime``; the landing attempts are not bound by
  ``dtime_min`` and ``dtime_factor``;
- the strain rate's sign must move the value towards the target (a positive strain rate is
  taken to increase the load); as the value at the start of a step is known only during the
  run (after a strain step it is the load reached; after a load step, the target within its
  tolerance), this is checked when the step starts, and before the run where possible; a
  target within its tolerance of the value at the start of its step is an error;
- before the run, a load (stress) target of step 1 must be more than :math:`\delta` from 0,
  and consecutive load (stress) targets more than :math:`2\delta` apart (each step ends within
  :math:`\delta` of its target, so the next step then starts more than :math:`\delta` from its
  own);
- under ``triaxial``, the other components with rows give a ``strain_rate``, a
  ``target_load`` or a ``target_stress`` in such a step (a ``target_strain`` would need the
  step's duration);
- the steps after a load or stress target start from the time and strains reached, and are
  planned when they start.

.. confval:: target_load

   The load [N] of a component at the end of each step: the component, then one value per step.
   Of the control component, it ends the step; of another component (``triaxial``), it is the
   end of the load path (`Triaxial load control`_).

   .. include:: generated/keywords/target_load.inc

.. confval:: target_stress

   As ``target_load``, for the stress [MPa]: the load over the face's current area.

   .. include:: generated/keywords/target_stress.inc

.. confval:: dtime_min

   The smallest time increment [s] of each step whose duration is not known (``-`` in the other
   steps); at most ``dtime``.

   .. include:: generated/keywords/dtime_min.inc

.. confval:: dtime_factor

   The factor (at least 1) on the time the last increment's rate needs to reach the target,
   which gives the next increment near a target.

   .. include:: generated/keywords/dtime_factor.inc

.. confval:: load_tol

   The tolerance :math:`\delta_F` [N] of the loads (targets and load-controlled values).

   .. include:: generated/keywords/load_tol.inc

.. confval:: stress_tol

   The tolerance :math:`\delta_\sigma` [MPa] of the stresses.

   .. include:: generated/keywords/stress_tol.inc

.. confval:: max_bc_iter

   The largest number of landing attempts of an increment, and of iterations of each loop of the
   boundary-velocity iteration.

   .. include:: generated/keywords/max_bc_iter.inc


Example: loading to 1 % strain, unloading to zero stress, reloading to 1.2 % strain:

.. literalinclude:: ../tests/regression/cases/mixed_unload_sx_fcc/simulation.cfg
   :start-at: boundary_conditions
   :end-before: global_tol 1e-10

Load and Stress Rates, Dwells
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The control component can be driven at a constant **load rate** (``load_rate <c>`` [N/s]) or
**stress rate** (``stress_rate <c>`` [MPa/s], of the load over the face's current area) instead
of a strain rate (constant load rate). Its face velocity is then found in
every increment so that its load (stress) has its value at the end of the increment, within
``load_tol`` (``stress_tol``), by an iteration (see below). With a uniaxial preset the
loaded component is driven this way; with ``triaxial`` the other components may be
strain-controlled, free or load-controlled (`Triaxial load control`_). ``general`` conditions
take no rate rows. A ``load_rate`` ends at:

- a ``target_load`` (a ``stress_rate`` at a ``target_stress``): the step lasts
  :math:`|q^* - q_\mathrm{start}|/|\dot q|` from the value at its start, and its increments
  are given as for a strain target (``number_of_incrs``, or ``dtime`` with a shorter last
  increment); the rate's sign must move the value towards the target (checked before the run
  where possible, else when the step starts), and the target must not be within the tolerance
  of the start value; before the run, step 1's target must be more than :math:`\delta` from 0
  and consecutive targets more than :math:`2\delta` apart, as for load targets;
- a ``target_time`` or ``target_duration``: the rate holds until that time;
- a ``target_strain``: the step ends at the first increment that reaches the strain of the
  control component's face (its duration is not known: ``dtime`` and ``dtime_min`` give the
  largest and smallest increments, the next increment aimed at the target as for a load target,
  ``dtime_factor``); the rate's sign must move the strain towards it (a positive rate increases
  the strain). The other components then take rates only (``strain_rate``, ``load_rate``,
  ``stress_rate``).

A rate **0** is a **dwell**: the value is held for the step's ``target_time`` or
``target_duration`` (a target cannot be reached at rate 0); step 1 cannot be a dwell (the load
starts at 0). Its increments: ``dtime`` (required; not ``number_of_incrs`` or
``dstrain``) is the largest increment; the first increment of a dwell that follows a step that
is not a dwell is 0.01 ``dtime``; then each increment is the time in which the control face, at
its last velocity, moves by ``max_strain_incr`` (default 0.001) times its original length, at
most ``dtime``; the last increment ends the dwell at its time exactly. A dwell above yield
shows creep (the load held, the strain increasing); below yield the velocities vanish.

The components other than the control component may give a ``load_rate`` or ``stress_rate``
too: their value then moves from its value at the start of the step at that rate (rate 0
holds it).

.. confval:: load_rate

   The load rate [N/s] of a component in each step: the component, then one value per step;
   ``0`` is a dwell.

   .. include:: generated/keywords/load_rate.inc

.. confval:: stress_rate

   The stress rate [MPa/s] of a component in each step, as ``load_rate``.

   .. include:: generated/keywords/stress_rate.inc

.. confval:: max_strain_incr

   The strain of the control face (times its original length) by which each increment of a dwell
   moves it, at most ``dtime``.

   .. include:: generated/keywords/max_strain_incr.inc


The iteration has two loops, each at most
``max_bc_iter`` (default 10) iterations, after which the run stops with an error naming the
loop. The first scales every load-controlled face velocity until the control value is within
:math:`\max(\delta, 0.1\,|q^* - q_\mathrm{init}|)` of its value (the factor 1.1, or 0.9 in a
dwell, then a secant through the last two solves); it stops after its first solve at the start
of a dwell or of a reload, in an unloading step, or with a zero control velocity. The second is
a finite-difference Newton iteration over the load-controlled face velocities (one solve with
each perturbed, the correction from the linear system, a solve), until every value is within
its tolerance; it makes no correction when the values already are. Its perturbations start at
1e-6 mm/s at the start of a dwell, otherwise at the larger of 0.1 times the control velocity's
change in the increment and 0.01 times the control velocity, then the size of the last
correction, and never below ``min_pert_frac`` (default 0.001) times the largest elastic
estimate of the face velocities for the run's largest targets over the increment (so that a
dwell below yield, with vanishing velocities, still gives a usable linear system). A step is
**loading** when the control value's magnitude increases (from zero, or away from zero without
crossing it), **unloading** when it decreases or crosses zero, and a dwell at rate 0; a step at
a strain rate counts by the change of its control-face load. The first increment of an episode
starts: at the run's first increment and at a reload (loading after a step that
was not loading) from an elastic estimate of the face velocities for the increment's load
change (after a solve at zero velocity, at a reload); at an unload start from a tenth of the
velocities; at a dwell start from the velocities stored at the first increment of the last dwell
that followed the same kind of step (zero the first time). Necking is not tested at a load or
stress rate. A run whose step 1 is at a load or stress rate has no initialization step
(``initialization elastic_step``).

.. confval:: min_pert_frac

   The smallest velocity perturbation of the boundary-velocity iteration, as a fraction of the
   reference velocity (above, and `Triaxial load control`_).

   .. include:: generated/keywords/min_pert_frac.inc


Example (with dwells): 33 at 2 N/s to 100 N, a 20 s dwell, to 250 N, a 10 s dwell, to 300 N;
11 and 22 load-controlled, their loads held during the dwells (``load_rate 0``):

.. literalinclude:: ../tests/regression/cases/triaxial_clr_dwell/simulation.cfg
   :start-at: boundary_conditions
   :end-before: initialization

Example: a uniaxial creep test, 100 N/s to 440 N, held for 2 s and 1 s (two dwell steps), then
unloaded at 200 N/s to 100 N:

.. literalinclude:: ../tests/verification/cases/restart_dwell/full.cfg
   :start-at: boundary_conditions
   :end-before: initialization

Triaxial Load Control
~~~~~~~~~~~~~~~~~~~~~

Under ``boundary_conditions triaxial``, a component other than the control component may give a
``target_load``, ``target_stress``, ``load_rate`` or ``stress_rate`` in a step: it is then
**load-controlled** in that step (as the control component itself at a load or stress rate,
`Load and stress rates, dwells`_, whose iteration then finds all the load-controlled face
velocities). Its
face ``x1``, ``y1`` or ``z1`` stays plane and moves at a uniform normal velocity (as a
strain-controlled face), and that velocity is found, in every increment, so that the face's load
(stress: the load over the face's current area) has the value it must have at the end of the
increment, within ``load_tol`` [N] (``stress_tol`` [MPa]):

- when the control component has a load or stress target (constant strain rate to a load), the loads
  follow the step's **load path**: each load-controlled value moves from its value at the start of
  the step to its target in proportion to the control component's progress towards its own target,
  :math:`q_j = q_{j,\mathrm{start}} + \lambda\,(q_j^* - q_{j,\mathrm{start}})`,
  :math:`\lambda = (q_p - q_{p,\mathrm{start}})/(q_p^* - q_{p,\mathrm{start}})`
  (the starting values are those reached at the start of the step; step 1 starts at 0). The
  step ends when the control component lands on its target;
- when the control component has a strain or time target (the step's duration is known), or a
  load or stress rate to a target or a time, each load-controlled value moves **linearly in
  time** from its start value to its target over the step;
- a ``load_rate`` or ``stress_rate`` row moves the value from its start value at that rate, in
  any step.

A target equal to the start value holds the value. ``target_load 22 0`` in every step is
therefore not a free face: the face ``y1`` stays plane with zero net normal force.

The face velocities are found by a boundary-velocity iteration: a solve of the increment
at the current face velocities; if a load-controlled value is off by more than its tolerance,
one more solve with each load-controlled face velocity changed by a small perturbation, and the
velocity correction from the linear system these solves give (with the path, the corrected
values lie on the path); then again, at most ``max_bc_iter`` (default 10) iterations, after
which the run stops with an error. Each solve starts from the start of the increment; the loads
are those of the end of the increment (the geometry updated). The perturbations are the size of
the last correction, and at least ``min_pert_frac`` (default 0.001) times the velocity of the
control face, which must therefore not be zero in such a step. A velocity change of a face is
spread over the body as a uniform strain rate along its axis, which makes a good initial guess
of the next solve. The face velocities are carried from increment to increment and from step to
step; a run whose step 1 has load-controlled components starts from an elastic estimate
(isotropic, with the Voigt averages of the phases' elastic constants weighted by their volume
fractions). With ``subincrementation on``, the iteration runs in every sub-increment, and a
failed solve of the iteration rejects the attempt like any failed increment. The landing on the
control component's target redoes whole increments, the iteration included. The log reports
every iteration (face velocities, values, values to reach), and the convergence file
(``print convergence``) numbers the iteration of every solve (:ref:`simulation_results`).

Example: compression along x at a constant strain rate to
-300 and -600 N, with y held at zero load and z pulled in proportion:

.. literalinclude:: ../tests/regression/cases/triaxial_csr/simulation.cfg
   :start-at: boundary_conditions
   :end-before: initialization

Example: strain control of 33 to 0.2 % and back to 0.1 %, with the stress of 11 interpolated in
time to -50 MPa and back to 0 (22 free):

.. literalinclude:: ../tests/verification/cases/fepx_iso_triaxial_time/simulation.cfg
   :start-at: boundary_conditions
   :end-before: ## Output

Control Component
~~~~~~~~~~~~~~~~~

One component controls the duration of the steps::

    control_component <c>

It is required with ``triaxial``. With a uniaxial preset it is the loaded component (inferred;
if given, it must match), and it is an error with ``general``. In every step:

- the control component has exactly one target (``target_strain <c>``, ``target_load <c>``,
  ``target_stress <c>``, ``target_time`` or ``target_duration``) and, with a preset, exactly
  one rate (``strain_rate <c>``, ``load_rate <c>`` or ``stress_rate <c>``); they set the step's
  duration;
- every other component with rows has exactly one of a ``strain_rate`` (its end strain follows
  from the duration), a ``target_strain`` (its rate is the strain change over the duration),
  or a ``target_load``, ``target_stress``, ``load_rate`` or ``stress_rate`` (it is
  load-controlled, `Triaxial load control`_).

.. confval:: control_component

   The component that sets the duration of the steps: ``11``, ``22`` or ``33``; with
   ``periodic`` conditions, an imposed component, shear components included.

   .. include:: generated/keywords/control_component.inc


Example: a triaxial run in which 33 is controlled, 11 follows a rate and 22 targets (step 2 with
``dstrain``):

.. literalinclude:: ../tests/verification/cases/fepx_hcp_sx_triaxial_elastic/simulation.cfg
   :start-at: boundary_conditions
   :end-before: ## Output

General Conditions and Steps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With ``general`` conditions, the steps have targets but no rate rows: the ``set_bc``
velocities apply. A ``target_strain <c>`` (on one component, the loading axis) is converted to
a time with the strain rate taken as the largest ``set_bc`` velocity divided by the length of
the domain along that axis; for a step towards a lower strain, the velocities are reversed.
``target_time`` steps use the velocities as given.

Run Limits
----------

FEPX stops with an error, after writing the results of the current increment as the current
step, when a limit is exceeded. The limits are checked after every increment (after the test
for the end of a step), and, where the steps allow it, before the run starts:

- ``max_incr <n>``: total number of increments (default 50000);
- ``max_total_time <t>``: simulated time (default 12000 s);
- ``max_strain <e>``: largest normal engineering strain of the domain (no default: no limit);
- ``max_eqstrain <e>``: equivalent engineering strain of the domain (no default: no limit).

.. confval:: max_incr

   The largest total number of increments.

   .. include:: generated/keywords/max_incr.inc

.. confval:: max_total_time

   The largest simulated time [s].

   .. include:: generated/keywords/max_total_time.inc

.. confval:: max_strain

   The largest normal engineering strain of the domain (no limit by default).

   .. include:: generated/keywords/max_strain.inc

.. confval:: max_eqstrain

   The largest equivalent engineering strain of the domain (no limit by default).

   .. include:: generated/keywords/max_eqstrain.inc


Necking
-------

``check_necking on`` (default ``off``) stops the run when the load on the control face (as
above) changes against the step's direction (the target's direction for a load or stress
target, the strain rate's sign otherwise), as when a specimen begins to neck. The test is made
after every increment; the results of the increment where the load first fell are written as
the current step, then the run stops with an error that names the step, the increment and the
two loads.

.. confval:: check_necking

   ``on``: the run stops when the load on the control face changes against the step's direction.

   .. include:: generated/keywords/check_necking.inc

