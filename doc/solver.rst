.. _solver_parameters:

Solver Parameters (:file:`simulation.cfg`)
==========================================

These parameters control the solution algorithms. All have defaults and
normally need no change. Each increment is solved at three nested levels, and the keywords are
named after the level they control:

1. the **global iteration** (keywords ``global_...``) finds the velocities at the end of the
   increment, by successive approximation and Newton--Raphson;
2. at each global iteration, the **material points** (keywords ``point_...``) find the stress
   at every quadrature point for the current velocities, and give the element matrices;
3. at each global iteration, the **linear solver** (keywords ``cg_...``) solves the assembled
   system for the velocity change.


Global Iteration (``global_...``)
---------------------------------

In each increment, the end-of-increment velocities are found by a global iteration:
it starts with successive approximation (SA, the secant stiffness) and switches to
Newton--Raphson (NR, the tangent stiffness) when the relative velocity change falls below
``global_nr_switch``. If NR diverges (the residual grows above 1.1 times its first value), it
returns to the last SA velocities and to SA, with a switch tolerance ten times smaller; if NR
converges slowly (its rate, measured from its second iteration, falls below ``global_nr_slow_ratio``
times the reference rate), it returns to SA for the rest of the increment.

The iteration has converged, from its second iteration on, when the norm of the velocity change
divided by the norm of the velocities is below ``global_tol``, or when the norm of the
velocity change is below ``global_tol_abs`` times the number of degrees of freedom (useful when the
velocities are near zero; with very small velocities this test can end the iteration before
``global_tol`` is reached). After a slow NR, or in SA after NR was tried, the looser
``global_tol_slow`` is accepted. It stops with an error after ``global_max_iters`` iterations, or
when the SA change grows five times in a row.

NR's tangent is the derivative of the stress with respect to the strain rate at the slip strengths
and orientations of the current iterate: their change within the increment is left out. NR
therefore converges linearly, each iteration reducing the change by a ratio that grows with the
hardening over an increment. At the default tolerances this costs a few iterations. With a very
tight ``global_tol`` (1e-10, say) and a fast-hardening material, large increments can need more
than ``global_max_iters`` iterations: use smaller increments (the ratio falls about in proportion)
or a larger ``global_max_iters``.

.. confval:: global_max_iters

   The largest number of global iterations of an increment.

   .. include:: generated/keywords/global_max_iters.inc

.. confval:: global_tol

   The relative tolerance of the global iteration: the norm of the velocity change over the norm
   of the velocities.

   .. include:: generated/keywords/global_tol.inc

.. confval:: global_tol_abs

   The absolute tolerance [mm/s] of the global iteration, times the number of degrees of
   freedom.

   .. include:: generated/keywords/global_tol_abs.inc

.. confval:: global_tol_slow

   The relative tolerance accepted after a slow NR, or in SA after NR was tried.

   .. include:: generated/keywords/global_tol_slow.inc

.. confval:: global_nr_switch

   The relative velocity change below which SA switches to NR; ``0`` keeps SA throughout.

   .. include:: generated/keywords/global_nr_switch.inc

.. confval:: global_nr_slow_ratio

   NR is slow when its rate falls below this fraction of its reference rate.

   .. include:: generated/keywords/global_nr_slow_ratio.inc


The log reports, for each increment, the number of global iterations and how many of them were
NR iterations.

Material Points (``point_...``)
-------------------------------

A *material point* (or *point*) is a quadrature point of an element: the place where the
constitutive equations are solved. Elemental results are written at the point at the element
centroid (see :ref:`simulation_results`). At every global iteration, each point is given its
deformation rate :math:`\mathbf d` and spin :math:`\mathbf w` from the current velocities,
and solves for its stress, slip-system strengths and lattice orientation at the end of the
increment. For a crystal-plasticity phase this is three nested loops (the
names are those of the theory documentation and of the program-flow chart):

1. the **state loop**, a fixed-point iteration over the state of the point: a stress solve,
   then a strength update, then an orientation update, repeated until the stress and the
   strengths settle: the change of the stress norm and of every strength between two
   iterations is below ``point_tol_state`` times :math:`\tau_0` (of the first slip system, of
   the point's element); at most ``point_max_iters_state`` iterations;
2. the **stress solve** (in each state iteration): Newton's method with a line search on the
   time-discretized stress--deformation-rate relation of the crystal, for the current
   orientation and strengths, until the norm of its residual is below ``point_tol`` [1/s]; at
   most ``point_max_iters_newton`` iterations;
3. the **strength update** (in each state iteration): backward Euler of the hardening law with
   the slip rates of the new stress, solved by Newton's method for each slip system
   :math:`\alpha` until its residual is below ``point_tol_hardening`` times
   :math:`\tau_0^\alpha`; at most ``point_max_iters_hardening`` iterations.

For a phase with irradiation softening (``irradiation on``), the state loop is damped: a
growing damage tensor lowers the strengths, which raises the slip rates, and the plain
fixed-point iteration can overshoot and oscillate. From the second state iteration on, only
a fraction :math:`\omega` of each change of the dislocation density and the damage tensor is
applied. :math:`\omega` is the Aitken factor computed from the last two changes (Küttler and
Wall, 2008, Comput. Mech. 43, 61, Eq. 44). It falls below 1 as soon as the changes show an
overshoot, never goes below 0.001 or above 1, is 1 when the changes show no overshoot, returns to 1
when the overshoot stops, and starts at 1 in every point solve. The convergence test above is
applied to the undamped change, so the converged result is that of the plain loop. Other phases are
not damped. There is no keyword.

A point whose solve does not converge keeps its last estimate, and the number of such points
is reported as a warning. For an ``elastic`` phase there is no strength update, the stress
solve is linear, and the state loop iterates between the stress and the lattice spin until
the relative change of the stress norm is below ``point_tol_state``.

.. confval:: point_tol

   The tolerance of the stress solve, on the norm of the residual of the deformation-rate
   equation [1/s].

   .. include:: generated/keywords/point_tol.inc

.. confval:: point_tol_state

   The tolerance of the state loop, relative to :math:`\tau_0` (elastic phases: to the stress
   norm).

   .. include:: generated/keywords/point_tol_state.inc

.. confval:: point_tol_hardening

   The tolerance of the strength update, relative to each system's :math:`\tau_0`.

   .. include:: generated/keywords/point_tol_hardening.inc

.. confval:: point_max_iters_newton

   The largest number of Newton iterations of the stress solve.

   .. include:: generated/keywords/point_max_iters_newton.inc

.. confval:: point_max_iters_hardening

   The largest number of Newton iterations of the strength update. Slip systems whose update
   does not converge keep their start-of-increment strength, and the point is counted in the
   warning.

   .. include:: generated/keywords/point_max_iters_hardening.inc

.. confval:: point_max_iters_state

   The largest number of state-loop iterations (also used by the elastic model).

   .. include:: generated/keywords/point_max_iters_state.inc


Smaller tolerances bring the solution closer to the exact solution of the discretized
equations, at a higher cost.

Initialization
--------------

.. confval:: initialization

   How a run with a crystal-plasticity phase starts: the velocity guess of the first increment
   and the stress its material points start from (``elastic_step`` or
   ``isotropic_viscoplastic``, below).

   .. include:: generated/keywords/initialization.inc

``elastic_step`` runs a small elastic step before the first increment, which provides the
initial estimate of the first increment: the strain that gives about 1/100 of the expected
yield stress, :math:`3 \tau_0` with the smallest ``tau_0``, at the fastest strain rate of
step 1, or half of the first increment if that is shorter. The step is part of step 1 (the
first increment is shortened by it), and its results are not written. Its velocities are the
guess of the first increment, whose material points start from the stress it reached.

``isotropic_viscoplastic`` starts in two parts. The velocity guess of the
first increment is the solution of an *isotropic viscous* problem on the initial mesh with the
conditions of step 1 (and their multi-point or periodicity constraints): a linear,
incompressible fluid of unit viscosity, its incompressibility imposed by three
augmented-Lagrangian iterations with a constant pressure per element. It is solved once, after
step 0 is written (the log shows the change of each iteration); it costs about three linear
solves. The three iterations do not make the field incompressible (uniaxially, the lateral
rates are :math:`-1/18` of the axial one, not :math:`-1/2`): it is only a starting guess. Then, in
every global iteration of the first increment, each material point starts its state loop from half
of a *viscoplastic estimate* of its stress, the stress of the rigid-viscoplastic crystal (a power
law with :math:`m = 0.4`) for its deformation rate, found by Newton's method from the vertex of
the single crystal yield surface that does the most plastic work on that rate (up to four
vertices are tried). The vertices are computed at the start of the run from each phase's
active slip systems, its ``c_over_a`` and its ``tau_0`` (Kasemer and Dawson, 2025): every
stress at which five linearly independent slip systems reach their strength and no system
exceeds it. The log gives their number and the time taken (well under a second for every
crystal type). The active slip systems must contain five independent ones (e.g. ``hcp`` with
only basal and prismatic slip does not): the run stops with an error otherwise.

Sub-Incrementation
------------------

By default, a material point whose solve does not converge keeps its last estimate (with a
warning), and a global iteration that does not converge (``global_max_iters`` reached, or
diverging) stops the run. With sub-incrementation, FEPX instead repeats such an increment as
shorter ones:

.. confval:: subincrementation

   ``on``: an increment *attempt* is rejected when a material point solve did not converge at
   the converged velocities, or when the global iteration did not converge
   (``global_max_iters``) or diverged. The attempt is undone and its time interval is solved as
   two halves, the second starting from the first; each half is divided again only if it is
   rejected itself. Each accepted part is solved like an increment (velocities, point solves,
   geometry update), but only the whole increment is written (the increments of the output and
   of the run limits are those of the step). A failure of the linear solver (``cg_max_iters``)
   still stops the run. The initialization step (``initialization elastic_step``) is not
   subdivided.

   .. include:: generated/keywords/subincrementation.inc

.. confval:: subincrement_max_depth

   An increment may be divided at most this many times :math:`n`, i.e. down to parts of
   :math:`\Delta t / 2^n`; a rejected part of that length stops the run with an error asking for
   finer increments (``number_of_incrs``, ``dtime`` or ``dstrain``) or a larger
   ``subincrement_max_depth``.

   .. include:: generated/keywords/subincrement_max_depth.inc

.. confval:: subincrement_max_successive

   At most this many successive increments may be subdivided; the next one that needs it stops
   the run with a similar error. An increment solved without subdivision resets the count.

   .. include:: generated/keywords/subincrement_max_successive.inc

The log reports each rejected
attempt (its length, the reason and the new depth) and, per subdivided increment, the number
of parts and attempts and the shortest part; the convergence file
(``print convergence``, :ref:`simulation_results`) then lists every attempt. With no rejected
attempt the results are identical, bit for bit, to those with ``subincrementation off``. The method
is the usual cutback of implicit finite element programs (M. Bäker, 2018, arXiv:1811.05753); it pays
the cost of shorter increments only where they are needed. Since a point that fails in a stable
2-cycle of its state loop may fail at every increment length, the limits keep the cost bounded.

Linear Solver (``cg_...``)
--------------------------

Each iteration solves a linear system by the conjugate gradient method, element by element (no
global matrix is assembled), with a nodal block-Jacobi preconditioner by default (Jacobi with
constraints; ``cg_preconditioner`` below). With
multi-point constraints (:ref:`multi_point_constraints`) the system is
reduced to the independent degrees of freedom, :math:`T^T K T\,\Delta\hat v = T^T r` with
:math:`v = T \hat v + g`, and preconditioned with the exact diagonal of :math:`T^T K T`; the
convergence tests use the full velocity field, as without constraints. Periodic conditions
(:ref:`periodic_conditions`) reduce the system in the same way; their preconditioner is the lumped
diagonal (``periodic_preconditioner lumped_diagonal``, the default: the element diagonals summed
onto the primary degrees of freedom, the couplings left out) or, with ``periodic_preconditioner
exact_diagonal``, the exact diagonal of :math:`T^T K T`. Without velocity conditions the periodic
system is singular by the rigid translations: its right-hand side is projected onto the range (its
translation components removed; Kaasschieter 1988, Bochev & Lehoucq 2005), so that the conjugate
gradient converges also when the residual is at round-off, and each velocity change is shifted to
the zero-mean fluctuation.

.. confval:: cg_tol

   The relative tolerance of the conjugate gradient: the iteration stops when the norm of the
   preconditioned residual :math:`M^{-1} r` falls below ``cg_tol`` times its initial value.

   .. include:: generated/keywords/cg_tol.inc

.. confval:: cg_max_iters

   The largest number of conjugate-gradient iterations of a solve.

   .. include:: generated/keywords/cg_max_iters.inc

.. confval:: cg_preconditioner

   The preconditioner :math:`M`, an approximation of the matrix that is cheap to invert; each
   iteration applies :math:`M^{-1}` to the residual. A good preconditioner lowers the number of
   iterations, and so the time of the solve, much more than it costs. The default is
   ``block_jacobi``, and ``jacobi`` with multi-point constraints or periodic conditions.

   - ``jacobi``: :math:`M` is the diagonal of the matrix. It rescales the unknowns so that they
     are comparable (the diagonal varies by orders of magnitude between corner and mid-side
     nodes, between small and large elements, and between elastic and flowing grains). It costs
     almost nothing; it is the preconditioner of runs with constraints.
   - ``block_jacobi``: :math:`M` is made of the :math:`3 \times 3` diagonal blocks of the matrix,
     one per node (the three velocity components of a node and their coupling), restricted to the
     free components. It also accounts for the coupling of a node's components, which the
     diagonal leaves out, and does not depend on the orientation of the axes. An iteration costs
     about as much as with ``jacobi`` (a :math:`3 \times 3` product per node), and the solves need
     fewer iterations: on polycrystals and single crystals, cubic and hexagonal, with 10-node
     tetrahedra and 20-node hexahedra, about 10 to 20 % fewer, which makes the run a few percent
     to 10 % faster, with the same nonlinear iterations and the same results to the solver
     tolerance. The default without constraints.
   - ``none``: no preconditioner (:math:`M` is the identity). For debugging only (to compare
     with a preconditioned solve); not intended for general use: the solves need many more
     iterations and may reach ``cg_max_iters``.

   The preconditioner applies to every solve of the run. With multi-point constraints or periodic
   conditions only ``jacobi`` is available (its diagonal is the one described above): it is the
   default there, and another value is an error.

   .. include:: generated/keywords/cg_preconditioner.inc

Mesh
----

.. confval:: midnode_linearization

   ``on``: at the start of every increment, the edge (mid-side) nodes of the 10-node tetrahedra
   are moved to the midpoints of their edges. This can postpone negative Jacobians with poor
   meshes or large deformations. A negative-Jacobian error suggests trying it.

   .. include:: generated/keywords/midnode_linearization.inc

