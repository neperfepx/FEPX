.. _running:

Convergence Parameters (:file:`simulation.cfg`)
==================================================

These options modify the tolerances and general behavior of the solution algorithms, and all have default values already defined.
They should only be modified by those who know what they are doing.

Velocity Convergence
--------------------

The velocity solver employs a hybrid successive-approximation/Newton-Raphson algorithm. Convergence of the velocity solution is based on a convergence parameter, which unless otherwise noted, is defined as the norm of the change in the velocity field, divided by the norm of the velocity field, :math:`||\Delta u||/||u||`. Other parameters are also used to assess the convergence of the velocity solution. The following parameters pertain to the convergence of the velocity solver:

- :data:`nl_max_iters \<iterations\>`: maximum allowable number of iterations of the nonlinear velocity solver (default: :data:`50`).

- :data:`nl_tol_strict \<tolerance\>`: desired tolerance on the elasto-viscoplastic velocity solution (default: :data:`5e-4`).

- :data:`nl_tol_loose \<tolerance\>`: acceptable level of convergence if the desired level of convergence cannot be reached via :data:`nl_tol_strict` (default: :data:`5e-4`).

- :data:`nl_tol_min \<tolerance\>`: tolerance on the norm of the change in velocity, divided by the number of degrees of freedom, :math:`(||u||/ \rm max(ndof))`. This parameter is useful for assessing convergence when the macroscopic velocity is near zero (default: :data:`1e-10`).

- :data:`nl_tol_switch_ref \<tolerance\>`: value of the convergence parameter at which the solution algorithm switches from successive-approximation to Newton-Raphson. To only use successive-approximations, set the value of :data:`nr_tol_switch_ref` equal to the value of :data:`nl_tol_strict` (default: :data:`1e-2`).

- :data:`nl_tol_conv \<tolerance\>`: parameter between 0 and 1 that is used to assess whether the Newton-Raphson algorithm is converging slowly (default: :data:`0.2`).

Conjugate Gradient Convergence
------------------------------

The solution of the linear system of equations :math:`[K]\{\Delta u\} = -\{R\}` is performed using a conjugate gradient solver. The following parameters pertain to the convergence of the conjugate gradient solver:

- :data:`cg_max_iters \<iterations\>`: maximum allowable number of iterations of the conjugate gradient solver (default: :data:`16000`).

- :data:`cg_tol \<tolerance\>`: desired tolerance on the conjugate gradient solver (default: :data:`1e-8`).

Material State Convergence
--------------------------

The convergence of the material stress state for both the viscoplastic and elasto-viscoplastic solutions is assessed by the following parameters:

- :data:`sx_max_iters_state \<iterations\>`: maximum number of iterations on material state (default: :data:`100`).

- :data:`sx_max_iters_newton \<iterations\>`: maximum number of iterations of the Newton algorithm used to solve for crystal stress (default: :data:`100`).

- :data:`sx_tol \<tolerance\>`: tolerance on the stress solution (default: :data:`1e-4`).

General Convergence
-------------------

These options may control standard simulation behavior or pertain to specific deformation modes.

- :data:`max_incr \<increments\>` specifies the maximum number of increments (default: :data:`50000`).

- :data:`max_total_time \<time\>` specifies the maximum deformation time (default: :data:`12000.0`).

- :data:`load_tol \<tolerance\>` is the the target load tolerance. A small positive load tolerance (e.g. 0.1 :math:`\times` control surface area) improves load control while reducing the number of small steps near target loads (default: :data:`0.0`).

- :data:`dtime_factor \<time\>` is a number greater than or equal to 1 which is used when calculating time increments near target loads (default: :data:`1.001`).

- :data:`max_bc_iter \<iterations\>` specifies the maximum number of boundary condition iterations (default: :data:`10`).

- :data:`min_pert_frac \<fraction\>` is the minimum fraction of the control velocity by which the secondary and tertiary surface velocities are perturbed during boundary condition iterations (default: :data:`0.001`).

- :data:`load_tol_abs \<tolerance\>` is the absolute tolerance on the secondary and tertiary loads. The absolute load criterion is that both loads are within the absolute load tolerance of the ideal load. Loads are considered to be within tolerance if either the absolute or relative criterion is satisfied (default: :data:`0.1`).

- :data:`load_tol_rel \<tolerance\>` is the relative load tolerance on the secondary and tertiary loads. It represents a fraction of the load in the control direction. The relative load criterion is that the difference between the load and ideal load, normalized by the load in the control direction, is less than the relative load tolerance. Loads are considered to be within tolerance if either the absolute or relative criterion is satisfied (default: :data:`0.001`).

- :data:`max_strain_incr \<increments\>` specifies the maximum strain increment for dwell episodes (default: :data:`0.001`).

- :data:`max_strain \<strain\>` specifies the maximum allowable macroscopic strain (default: :data:`0.2`).

- :data:`max_eqstrain \<strain\>` specifies the maximum allowable macroscopic equivalent strain (default: :data:`0.2`).

- :data:`max_iter_hard_limit \<iterations\>` specifies the maximum allowable iterations on the Backward Euler approximation used to update hardnesses (default: :data:`10`).

Termination
-----------

- :data:`check_necking {on,off}` specifies whether or not to terminate simulation when specimen begins to neck (default: :data:`off`).

Other Parameters
----------------

- :data:`uniform_elts \<logical\>` specifies whether orientations and critical resolved shear stresses should be uniform inside individual elements, instead of having different values at different Gauss points (default :data:`false`).
