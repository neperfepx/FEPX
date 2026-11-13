! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Boundary-velocity iteration of triaxial load control, its pure algorithm parts (step 8.2b):
!> the values the load-controlled components must reach in an increment (on the step's load
!> path, or interpolated in time), the convergence test, the signs and sizes of the velocity
!> perturbations, the correction of the face velocities from the perturbed solves (v2's linear
!> system), and the elastic estimate of the face velocities (Voigt-averaged moduli). Step 8.2c:
!> the values of rate rows, the kind of a step (loading, unloading, dwell) and the parts of
!> iteration (B) (v2 CLR: the scaling loop's factor and test, the Newton loop's perturbation
!> sizes and their lower bound). The trial solves themselves (the increment solve and the face
!> forces) are made by the caller (fepx.f90::bc_iterate, bc_iterate_rate).
!>
!> Source : port of v2 triaxial/driver_triaxcsr_mod.f90::driver_triaxcsr (the loop
!>          bc_iteration: ideal loads at the control load reached, tolerance test, perturbation
!>          signs, the 3 x 3 system of the perturbed load changes and the load path, the update
!>          of the perturbation sizes with the lower bound min_pert_frac x the control face
!>          velocity, v2-fixed C58) and triaxial/driver_triax_utilities_mod.f90::est_avg_mod
!>          (Voigt-averaged Young's modulus E_V = (A - B + 3C)(A + 2B)/(2A + 3B + C) and
!>          nu = (3K - E_V)/(6K), Hosford, "The Mechanics of Crystals and Textured
!>          Polycrystals", 1993, p. 22, as cited in v2; every phase, v2-fixed C59) and the
!>          elastic velocity estimate of driver_triaxcsr's initialization; with the semantics of
!>          dev/theory/addendum_triaxial_load_control.md §3.1, §3.4 (decisions Q-T4, Q-T8 (ii),
!>          Q-T10, Q-T12, Q-T18; ROADMAP D66). Iteration (B) (step 8.2c): port of v2
!>          triaxial/driver_triaxclr_mod.f90::driver_triaxclr (the step actions "loading",
!>          "unloading", "dwelling", by magnitude, Q-T9; bc_iteration_1: the test
!>          |F_p - F*_p| < max(load_tol, 0.1 |F*_p - F_p,init|) and the scale 1.1 (0.9 in a
!>          dwell), then the secant; bc_iteration_2: the perturbation 1e-6 at a dwell start, else
!>          max(0.1 |V_p - V_p,init|, 0.01 |V_p|), then max |Delta V|, with the lower bound of
!>          v2-fixed C60, min_pert_frac times the largest elastic estimate of the face
!>          velocities for the run's largest target magnitudes over dt, Q-T1), addendum §3.4.
!>          The linear systems (size <= 4) are solved by
!>          linalg_mod::lu_factor/lu_solve (Golub & Van Loan 2013, Algorithm 3.4.1).
!>
!> Semantics (addendum §3.4, iteration (A), strain-rate control of the control component p):
!>   - The unknowns are the normal velocities V_j of the faces j1 of the load-controlled
!>     components j in U (n = |U| <= 2); each has a controlled value q_j: its face force F_j [N]
!>     (load rows) or F_j / A_j [MPa] (stress rows, A_j the face's current area).
!>   - Path rows (the control component has a load or stress target, the step's duration is not
!>     known; v2 CSR): the values on the step's path at the control's progress,
!>     q_j^id = q_j,start + lambda (q_j* - q_j,start), lambda = (q_p - q_p,start) /
!>     (q_p* - q_p,start) (bc_path_values). Time-based rows (the control has a strain or time
!>     target, duration T known; Q-T8 (ii)): q_j^id = q_j,start + (q_j* - q_j,start)(t - t_start)
!>     / T (bc_time_values). The start values are the values measured at the step's start (Q-T4).
!>   - Converged when |q_j - q_j^id| <= tol_j for every j in U (load_tol [N] or stress_tol [MPa],
!>     absolute only, Q-T18) (bc_converged).
!>   - Otherwise each V_j is perturbed by s_j delta_j, s_j = +1 when q_j <= q_j^id, else -1
!>     (bc_pert_signs, v2), giving the column dq_j of the changes of the controlled values; the
!>     correction (bc_correction) solves, with path rows, v2's system over the rows U and p,
!>     sum_k dq_ik c_k - c_(n+1) r_i = q_i,start - q_i (r_i = (q_i* - q_i,start) /
!>     (q_p* - q_p,start), r_p = 1: the corrected values lie on the path at the control value
!>     q_p,start + c_(n+1)), and with time-based rows the n x n system
!>     sum_k dq_jk c_k = q_j^id - q_j; then Delta V_k = c_k s_k delta_k. With n = 2 path rows this
!>     is v2's 3 x 3 system.
!>   - The next perturbation sizes are delta_k = max(|Delta V_k|, min_pert_frac |V_p|)
!>     (bc_pert_update; v2-fixed C58).
!>   - Elastic estimate (bc_elastic_moduli, bc_elastic_velocities): isotropic Hooke's law with the
!>     Voigt-averaged E and nu of the phases, averaged by phase volume fraction (Q-T12; v2: by
!>     element counts): V_i = (s_i - nu sum_(j /= i) s_j) L_i / E for the normal stress changes
!>     (rates) s of the faces. The initial face velocities of a run (bc_initial_velocities): v2's
!>     estimate for step 1's targets scaled so that the control face moves at its strain rate
!>     (path rows); with time-based rows the control's stress rate from Hooke's law so that its
!>     face moves at the strain rate, the others at their targets' rates.
!>   - Rate rows (step 8.2c, addendum §3.1): a load_rate or stress_rate r gives the value
!>     q_start + r (t - t_start) (bc_mode_values, with the path and time-based values).
!>   - Step kind (Q-T9, addendum §3.4; bc_step_action): a control rate 0 is a dwell; otherwise
!>     loading when the control value's magnitude increases (|q_end| > |q_start| with q_end and
!>     q_start of the same sign, or q_start within its tolerance of 0), unloading otherwise (a
!>     zero crossing unloads); with the end value not known, by the rate's sign against
!>     q_start's.
!>   - Iteration (B) (load- or stress-rate control and dwells, v2 CLR; U includes the control
!>     component p): a scaling loop, every face velocity of U scaled by s, s = 1.1 (0.9 in a
!>     dwell) at its first iteration, then v2's secant on the control,
!>     s = (1 - V_prev/V)(q* - q_prev)/(q - q_prev) + V_prev/V (bc_rate_scale), done when
!>     |q_p - q_p*| < max(tol_p, 0.1 |q_p* - q_p,init|) (bc_rate_first_converged; the caller's
!>     other exits: a dwell or reload start, unloading, a zero control velocity); then a
!>     finite-difference Newton loop over U (the correction of bc_correction with time-based
!>     rows), with the perturbation size delta = 1e-6 mm/s at a dwell start, else
!>     max(0.1 |V_p - V_p,init|, 0.01 |V_p|) (bc_rate_pert_initial), then max_k |Delta V_k|
!>     (bc_rate_pert_update), never below min_pert_frac max_j |V_j^el|, V^el the elastic
!>     estimate for the run's largest target magnitudes of each component over dt
!>     (bc_rate_pert_floor, v2-fixed C60, Q-T1).
!>
!> Public parameters:
!>   bc_mode_path, bc_mode_time, bc_mode_rate     how a value to reach is given (path, time,
!>                                                rate)
!>   bc_action_loading, bc_action_unloading, bc_action_dwell  step kinds (bc_step_action)
!>
!> Public routines:
!>   bc_path_values         values on the step's load path at the control's progress
!>   bc_time_values         values interpolated linearly in time
!>   bc_converged           every load-controlled value within its tolerance
!>   bc_pert_signs          perturbation signs towards the values to reach
!>   bc_correction          face-velocity correction from the perturbed solves
!>   bc_pert_update         perturbation sizes for the next iteration
!>   bc_elastic_moduli      Voigt-averaged E and nu of the phases
!>   bc_elastic_velocities  face velocities of isotropic Hooke's law for normal stress rates
!>   bc_initial_velocities  initial face velocities of a run with load-controlled components
!>   bc_mode_values         values to reach of rows given on the path, in time or as rates
!>   bc_step_action         kind of a step: loading, unloading or dwell (iteration (B))
!>   bc_rate_first_converged  test of iteration (B)'s scaling loop
!>   bc_rate_scale          scale factor of iteration (B)'s scaling loop
!>   bc_rate_pert_floor     lower bound of iteration (B)'s perturbations
!>   bc_rate_pert_initial   first perturbation size of iteration (B)'s Newton loop
!>   bc_rate_pert_update    next perturbation size of iteration (B)'s Newton loop
module bc_iteration_mod

  use kinds_mod, only: rk
  use linalg_mod, only: lu_factor, lu_solve

  implicit none

  private

  public :: bc_path_values, bc_time_values, bc_converged, bc_pert_signs, bc_correction
  public :: bc_pert_update, bc_elastic_moduli, bc_elastic_velocities, bc_initial_velocities
  public :: bc_mode_values, bc_step_action, bc_rate_first_converged, bc_rate_scale
  public :: bc_rate_pert_floor, bc_rate_pert_initial, bc_rate_pert_update
  public :: bc_mode_path, bc_mode_time, bc_mode_rate
  public :: bc_action_loading, bc_action_unloading, bc_action_dwell

  ! How the value to reach of a load-controlled row is given (step 8.2c): on the step's load
  ! path (a target, the step's duration not known), linearly in time (a target, the duration
  ! known), or by a rate.
  integer, parameter :: bc_mode_path = 1, bc_mode_time = 2, bc_mode_rate = 3

  ! Step kinds of iteration (B) (v2 driver_triaxclr's curr_action "loading", "unloading",
  ! "dwelling"; by magnitude, Q-T9).
  integer, parameter :: bc_action_loading = 1, bc_action_unloading = 2, bc_action_dwell = 3

  ! v2 driver_triaxclr's heuristics (addendum §3.4, Q-T9): the first scale of the scaling loop,
  ! in a dwell and otherwise; the fraction of the increment's control change that the scaling
  ! loop accepts; the first perturbation at a dwell start [mm/s]; the fractions of the control
  ! velocity's change in the increment and of the control velocity that give the first
  ! perturbation otherwise.
  real(rk), parameter :: rate_scale_first = 1.1_rk
  real(rk), parameter :: rate_scale_first_dwell = 0.9_rk
  real(rk), parameter :: rate_first_tol_frac = 0.1_rk
  real(rk), parameter :: rate_pert_dwell = 1.0e-6_rk
  real(rk), parameter :: rate_pert_change_frac = 0.1_rk
  real(rk), parameter :: rate_pert_vel_frac = 0.01_rk

contains

  !> Values on the step's load path at the control's progress (v2 CSR: ideal_load = step_frac
  !> target_k + (1 - step_frac) target_(k-1), with the measured start values, Q-T4):
  !>   q_ideal = q_start + lambda (q_target - q_start),
  !>   lambda = (q_ctrl - q_ctrl_start)/(q_ctrl_target - q_ctrl_start).
  !>
  !> Arguments:
  !>   q_start        (in)  real(rk)(:)  start values of the step [N] or [MPa]
  !>   q_target       (in)  real(rk)(:)  targets of the step [N] or [MPa]
  !>   q_ctrl         (in)  real(rk)     control value reached [N] or [MPa]
  !>   q_ctrl_start   (in)  real(rk)     control value at the step's start [N] or [MPa]
  !>   q_ctrl_target  (in)  real(rk)     control target (/= q_ctrl_start) [N] or [MPa]
  !>   q_ideal        (out, result) real(rk)(size(q_start))  values to reach [N] or [MPa]
  !> Locals:
  !>   lambda         real(rk)     progress of the control along its step [-]
  pure function bc_path_values(q_start, q_target, q_ctrl, q_ctrl_start, q_ctrl_target) &
    & result(q_ideal)

    real(rk), intent(in) :: q_start(:), q_target(:), q_ctrl, q_ctrl_start, q_ctrl_target
    real(rk) :: q_ideal(size(q_start))

    real(rk) :: lambda

    lambda = (q_ctrl - q_ctrl_start)/(q_ctrl_target - q_ctrl_start)
    q_ideal = q_start + lambda*(q_target - q_start)

  end function bc_path_values

  !> Values interpolated linearly in time over a step of known duration (addendum §3.1, v2 CLR's
  !> targets in time): q_ideal = q_start + frac (q_target - q_start).
  !>
  !> Arguments:
  !>   q_start   (in)  real(rk)(:)  start values of the step [N] or [MPa]
  !>   q_target  (in)  real(rk)(:)  targets of the step [N] or [MPa]
  !>   frac      (in)  real(rk)     (t - t_start)/T, the fraction of the step reached [-]
  !>   q_ideal   (out, result) real(rk)(size(q_start))  values to reach [N] or [MPa]
  pure function bc_time_values(q_start, q_target, frac) result(q_ideal)

    real(rk), intent(in) :: q_start(:), q_target(:), frac
    real(rk) :: q_ideal(size(q_start))

    q_ideal = q_start + frac*(q_target - q_start)

  end function bc_time_values

  !> Whether every load-controlled value is within its tolerance of the value to reach (v2:
  !> s_load_err <= max_load_err and t_load_err <= max_load_err; absolute tolerances only, Q-T18).
  !>
  !> Arguments:
  !>   q        (in)  real(rk)(:)  values reached [N] or [MPa]
  !>   q_ideal  (in)  real(rk)(:)  values to reach [N] or [MPa]
  !>   tol      (in)  real(rk)(:)  tolerances, > 0 [N] or [MPa]
  pure logical function bc_converged(q, q_ideal, tol)

    real(rk), intent(in) :: q(:), q_ideal(:), tol(:)

    bc_converged = all(abs(q - q_ideal) <= tol)

  end function bc_converged

  !> Signs of the velocity perturbations, towards the values to reach (v2: +1 when the load is
  !> not above its ideal value, else -1; a faster face increases its force).
  !>
  !> Arguments:
  !>   q        (in)  real(rk)(:)  values reached [N] or [MPa]
  !>   q_ideal  (in)  real(rk)(:)  values to reach [N] or [MPa]
  !>   signs    (out, result) real(rk)(size(q))  +1 or -1
  pure function bc_pert_signs(q, q_ideal) result(signs)

    real(rk), intent(in) :: q(:), q_ideal(:)
    real(rk) :: signs(size(q))

    signs = merge(1.0_rk, -1.0_rk, q <= q_ideal)

  end function bc_pert_signs

  !> Correction of the face velocities from the perturbed solves (v2 CSR's linear system).
  !> Path rows (path true): the rows are the n load-controlled components then the control
  !> component (m = n + 1 values), the unknowns c_1..c_n and the control value's change
  !> c_(n+1) along the path:
  !>   sum_k dq(i, k) c_k - c_(n+1) ratio(i) = q_aim(i) - q(i),  i = 1..n+1,
  !> with q_aim the start values and ratio(i) = (q_i* - q_i,start)/(q_p* - q_p,start) (1 for the
  !> control); time-based rows (path false): the n x n system
  !>   sum_k dq(j, k) c_k = q_aim(j) - q(j),  j = 1..n,
  !> with q_aim the values to reach. The velocity changes are dvel(k) = c_k pert(k).
  !>
  !> Arguments:
  !>   dq      (in)  real(rk)(:, :)  (m, n) change of each controlled value (rows) by each
  !>                                 perturbation (columns) [N] or [MPa]
  !>   q       (in)  real(rk)(:)     (m) values at the current velocities [N] or [MPa]
  !>   q_aim   (in)  real(rk)(:)     (m) start values (path) or values to reach (time) [N] or [MPa]
  !>   ratio   (in)  real(rk)(:)     (m) path direction of each value per control change (path;
  !>                                 ignored otherwise) [-]
  !>   pert    (in)  real(rk)(:)     (n) the signed perturbations s_k delta_k [mm/s]
  !>   path    (in)  logical         path rows (else time-based)
  !>   dvel    (out) real(rk)(:)     (n) changes of the face velocities [mm/s]
  !>   status  (out) integer         0: solved; 1: singular system (dvel = 0)
  !> Locals:
  !>   n, m    integer               load-controlled components, size of the system
  !>   a, lu   real(rk)(:, :)        (m, m) the system and its LU factors
  !>   b, c    real(rk)(:)           (m) right-hand side and solution
  !>   perm    integer(:)            (m) row permutation
  !>   det     real(rk)              determinant (unused)
  pure subroutine bc_correction(dq, q, q_aim, ratio, pert, path, dvel, status)

    real(rk), intent(in) :: dq(:, :), q(:), q_aim(:), ratio(:), pert(:)
    logical, intent(in) :: path
    real(rk), intent(out) :: dvel(:)
    integer, intent(out) :: status

    integer :: n, m
    real(rk), allocatable :: a(:, :), lu(:, :), b(:), c(:)
    integer, allocatable :: perm(:)
    real(rk) :: det

    ! Size of the system: the path adds the control value's change as an unknown (v2's third
    ! column -(target_k - target_(k-1))/(load step)).
    n = size(pert)
    m = n
    if (path) m = n + 1
    allocate (a(m, m), lu(m, m), b(m), c(m), perm(m))
    a(:, 1:n) = dq(1:m, 1:n)
    if (path) a(:, m) = -ratio(1:m)
    b = q_aim(1:m) - q(1:m)

    ! Solve, then the velocity changes of the perturbed faces.
    call lu_factor(a, lu, perm, det, status)
    dvel = 0.0_rk
    if (status /= 0) return
    call lu_solve(lu, perm, b, c)
    dvel = c(1:n)*pert

  end subroutine bc_correction

  !> Perturbation sizes of the next iteration (v2: the size of the last correction, at least
  !> min_pert_frac times the control face velocity, v2-fixed C58).
  !>
  !> Arguments:
  !>   dvel           (in)  real(rk)(:)  the last velocity changes [mm/s]
  !>   vel_ctrl       (in)  real(rk)     velocity of the control face [mm/s]
  !>   min_pert_frac  (in)  real(rk)     lower bound, fraction of |vel_ctrl| [-]
  !>   delta          (out, result) real(rk)(size(dvel))  the sizes [mm/s]
  pure function bc_pert_update(dvel, vel_ctrl, min_pert_frac) result(delta)

    real(rk), intent(in) :: dvel(:), vel_ctrl, min_pert_frac
    real(rk) :: delta(size(dvel))

    delta = max(abs(dvel), min_pert_frac*abs(vel_ctrl))

  end function bc_pert_update

  !> Voigt-averaged Young's modulus and Poisson's ratio of the phases (v2 est_avg_mod, Hosford
  !> 1993 p. 22): for each phase A, B, C the means of C11, C22, C33, of C12, C13, C23 and of the
  !> (Voigt) C44, C55, C66 of its stiffness, E_V = (A - B + 3C)(A + 2B)/(2A + 3B + C), the Voigt
  !> bulk modulus K = (A + 2B)/3 and nu = (3K - E_V)/(6K); E and nu averaged with the phases'
  !> weights (volume fractions, Q-T12). Frame-independent (Voigt averages are invariants).
  !>
  !> Arguments:
  !>   stiff_c  (in)  real(rk)(:, :, :)  (6, 6, num_phases) Mandel stiffness of each phase,
  !>                                     crystal frame [MPa]
  !>   weights  (in)  real(rk)(:)        (num_phases) volume fractions, summing to 1 [-]
  !>   e        (out) real(rk)           average Young's modulus [MPa]
  !>   nu       (out) real(rk)           average Poisson's ratio [-]
  !> Locals:
  !>   i        integer                  phase
  !>   a, b, c  real(rk)                 the means A, B, C of the phase [MPa]
  !>   e_ph, k_ph  real(rk)              the phase's E_V and K [MPa]
  pure subroutine bc_elastic_moduli(stiff_c, weights, e, nu)

    real(rk), intent(in) :: stiff_c(:, :, :), weights(:)
    real(rk), intent(out) :: e, nu

    integer :: i
    real(rk) :: a, b, c, e_ph, k_ph

    e = 0.0_rk
    nu = 0.0_rk
    do i = 1, size(weights)
      ! Means of the normal, the off-diagonal normal and the shear terms (Mandel shear
      ! diagonal = 2 x Voigt).
      a = (stiff_c(1, 1, i) + stiff_c(2, 2, i) + stiff_c(3, 3, i))/3.0_rk
      b = (stiff_c(1, 2, i) + stiff_c(1, 3, i) + stiff_c(2, 3, i))/3.0_rk
      c = (stiff_c(4, 4, i) + stiff_c(5, 5, i) + stiff_c(6, 6, i))/6.0_rk

      ! The phase's moduli, weighted into the averages.
      e_ph = (a - b + 3.0_rk*c)*(a + 2.0_rk*b)/(2.0_rk*a + 3.0_rk*b + c)
      k_ph = (a + 2.0_rk*b)/3.0_rk
      e = e + weights(i)*e_ph
      nu = nu + weights(i)*(3.0_rk*k_ph - e_ph)/(6.0_rk*k_ph)
    end do

  end subroutine bc_elastic_moduli

  !> Face velocities of isotropic Hooke's law for normal stress rates of the faces (v2's
  !> elastic estimate: V_i = (F_i/A0_i - nu sum_(j /= i) F_j/A0_j) L_i/(E dt)):
  !>   V_i = (s_i - nu sum_(j /= i) s_j) L_i / E.
  !>
  !> Arguments:
  !>   stress_rate  (in)  real(rk)(3)  normal stress rates of the faces x1, y1, z1 [MPa/s]
  !>   length       (in)  real(rk)(3)  lengths along x, y, z [mm]
  !>   e            (in)  real(rk)     Young's modulus [MPa]
  !>   nu           (in)  real(rk)     Poisson's ratio [-]
  !>   vel          (out, result) real(rk)(3)  face velocities [mm/s]
  pure function bc_elastic_velocities(stress_rate, length, e, nu) result(vel)

    real(rk), intent(in) :: stress_rate(3), length(3), e, nu
    real(rk) :: vel(3)

    vel = ((1.0_rk + nu)*stress_rate - nu*sum(stress_rate))*length/e

  end function bc_elastic_velocities

  !> Initial face velocities of a run whose step 1 has load-controlled components (addendum
  !> §3.4, "Initial velocities"). Path rows (v2 CSR): the elastic estimate of the step's target
  !> stresses (load-controlled components and the control; 0 elsewhere), scaled so that the
  !> control face moves at vel_ctrl (a zero control velocity of the estimate: the others 0).
  !> Time-based rows: the targets' stress rates of the load-controlled components, the control's
  !> stress rate from Hooke's law so that its face moves at vel_ctrl,
  !> s_p = E vel_ctrl / L_p + nu sum_(j /= p) s_j, then Hooke's law.
  !>
  !> Arguments:
  !>   stress    (in)  real(rk)(3)  path: target stress changes of the step (load targets over
  !>                                the original face areas) [MPa]; time-based: their rates
  !>                                [MPa/s] (the control's entry ignored); 0 for components
  !>                                without a load or stress target
  !>   axis      (in)  integer      control component (1..3)
  !>   vel_ctrl  (in)  real(rk)     velocity of the control face (strain rate x length) [mm/s]
  !>   length    (in)  real(rk)(3)  original lengths along x, y, z [mm]
  !>   e, nu     (in)  real(rk)     elastic moduli (bc_elastic_moduli) [MPa], [-]
  !>   path      (in)  logical      path rows (else time-based)
  !>   vel       (out, result) real(rk)(3)  face velocities [mm/s]
  !> Locals:
  !>   rates     real(rk)(3)        stress rates (time-based) [MPa/s]
  pure function bc_initial_velocities(stress, axis, vel_ctrl, length, e, nu, path) result(vel)

    real(rk), intent(in) :: stress(3), vel_ctrl, length(3), e, nu
    integer, intent(in) :: axis
    logical, intent(in) :: path
    real(rk) :: vel(3)

    real(rk) :: rates(3)

    if (path) then
      ! v2: the estimate of the targets, scaled to the control strain rate.
      vel = bc_elastic_velocities(stress, length, e, nu)
      if (abs(vel(axis)) > 0.0_rk) then
        vel = vel*(vel_ctrl/vel(axis))
      else
        vel = 0.0_rk
      end if
    else
      ! The control's stress rate that moves its face at vel_ctrl.
      rates = stress
      rates(axis) = 0.0_rk
      rates(axis) = e*vel_ctrl/length(axis) + nu*sum(rates)
      vel = bc_elastic_velocities(rates, length, e, nu)
    end if
    vel(axis) = vel_ctrl

  end function bc_initial_velocities

  !> Values to reach of load-controlled rows, each given on the step's load path, linearly in
  !> time or by a rate (addendum §3.1; step 8.2c):
  !>   path:  q_start + lambda (q_target - q_start)   (bc_path_values, v2 CSR)
  !>   time:  q_start + frac (q_target - q_start)      (bc_time_values, v2 CLR)
  !>   rate:  q_start + rate elapsed                   (a load_rate or stress_rate row)
  !>
  !> Arguments:
  !>   mode      (in)  integer(:)   bc_mode_path, bc_mode_time or bc_mode_rate of each row
  !>   q_start   (in)  real(rk)(:)  start values of the step [N] or [MPa]
  !>   q_target  (in)  real(rk)(:)  targets (path and time rows) [N] or [MPa]
  !>   rate      (in)  real(rk)(:)  rates (rate rows) [N/s] or [MPa/s]
  !>   lambda    (in)  real(rk)     progress of the control along its step (path rows) [-]
  !>   frac      (in)  real(rk)     (t - t_start)/T, the fraction of the step reached (time rows) [-]
  !>   elapsed   (in)  real(rk)     t - t_start (rate rows) [s]
  !>   q_ideal   (out, result) real(rk)(size(mode))  values to reach [N] or [MPa]
  pure function bc_mode_values(mode, q_start, q_target, rate, lambda, frac, elapsed) &
    & result(q_ideal)

    integer, intent(in) :: mode(:)
    real(rk), intent(in) :: q_start(:), q_target(:), rate(:), lambda, frac, elapsed
    real(rk) :: q_ideal(size(mode))

    where (mode == bc_mode_path)
      q_ideal = q_start + lambda*(q_target - q_start)
    elsewhere (mode == bc_mode_time)
      q_ideal = q_start + frac*(q_target - q_start)
    elsewhere
      q_ideal = q_start + rate*elapsed
    end where

  end function bc_mode_values

  !> Kind of a step of iteration (B) (v2 driver_triaxclr's curr_action, by magnitude, user
  !> decision Q-T9): a rate 0 is a dwell; otherwise loading when the magnitude of the control
  !> value increases, |q_end| > |q_start| with q_end and q_start of the same sign, or q_start
  !> within the tolerance of 0 (a step from zero); unloading otherwise (a zero crossing
  !> unloads). With the end value not known (a strain target), loading when q_start is within
  !> the tolerance of 0 or the rate has q_start's sign. (v2 compared the signed values, so
  !> loading in compression was "unloading", finding F5.)
  !>
  !> Arguments:
  !>   rate       (in)  real(rk)  rate of the control value (or its change over the step) [N/s],
  !>                              [MPa/s] (or [N], [MPa])
  !>   q_start    (in)  real(rk)  control value at the step's start [N] or [MPa]
  !>   q_end      (in)  real(rk)  control value at the step's end (when end_known) [N] or [MPa]
  !>   end_known  (in)  logical   q_end is known
  !>   tol        (in)  real(rk)  tolerance of the value, > 0 [N] or [MPa]
  !>   action     (out, result) integer  bc_action_loading, bc_action_unloading or
  !>                                     bc_action_dwell
  !> Locals:
  !>   from_zero  logical   q_start within tol of 0
  !>   grows      logical   the magnitude increases without a zero crossing
  pure function bc_step_action(rate, q_start, q_end, end_known, tol) result(action)

    real(rk), intent(in) :: rate, q_start, q_end, tol
    logical, intent(in) :: end_known
    integer :: action

    logical :: from_zero, grows

    ! A rate 0 holds the value.
    if (.not. abs(rate) > 0.0_rk) then
      action = bc_action_dwell
      return
    end if

    ! Loading when the magnitude grows from the start (or from zero).
    from_zero = .not. abs(q_start) > tol
    if (end_known) then
      grows = q_start*q_end > 0.0_rk .and. abs(q_end) > abs(q_start)
    else
      grows = rate*q_start > 0.0_rk
    end if
    if (from_zero .or. grows) then
      action = bc_action_loading
    else
      action = bc_action_unloading
    end if

  end function bc_step_action

  !> Test of the scaling loop of iteration (B) (v2 driver_triaxclr, bc_iteration_1): the control
  !> value within max(tol, 0.1 |q_target - q_init|) of its value to reach, q_init the value at
  !> the start of the increment (v2: strict).
  !>
  !> Arguments:
  !>   q         (in)  real(rk)  control value reached [N] or [MPa]
  !>   q_target  (in)  real(rk)  control value to reach [N] or [MPa]
  !>   q_init    (in)  real(rk)  control value at the start of the increment [N] or [MPa]
  !>   tol       (in)  real(rk)  tolerance of the value, > 0 [N] or [MPa]
  pure logical function bc_rate_first_converged(q, q_target, q_init, tol)

    real(rk), intent(in) :: q, q_target, q_init, tol

    bc_rate_first_converged = abs(q - q_target) < max(tol, rate_first_tol_frac* &
      & abs(q_target - q_init))

  end function bc_rate_first_converged

  !> Scale factor of the face velocities in the scaling loop of iteration (B) (v2
  !> driver_triaxclr, bc_iteration_1): 1.1 at the first iteration (0.9 in a dwell), then the
  !> secant on the control component through the last two solves,
  !>   s = (1 - V_prev/V)(q_target - q_prev)/(q - q_prev) + V_prev/V,
  !> the factor that takes V to the root of the line through (V_prev, q_prev), (V, q).
  !>
  !> Arguments:
  !>   iter      (in)  integer   iteration of the loop (1: the first)
  !>   dwell     (in)  logical   the step is a dwell
  !>   vel       (in)  real(rk)  control face velocity of the last solve (/= 0) [mm/s]
  !>   vel_prev  (in)  real(rk)  that of the solve before [mm/s]
  !>   q_target  (in)  real(rk)  control value to reach [N] or [MPa]
  !>   q         (in)  real(rk)  control value of the last solve [N] or [MPa]
  !>   q_prev    (in)  real(rk)  that of the solve before [N] or [MPa]
  !>   s         (out, result) real(rk)  the factor [-]
  pure function bc_rate_scale(iter, dwell, vel, vel_prev, q_target, q, q_prev) result(s)

    integer, intent(in) :: iter
    logical, intent(in) :: dwell
    real(rk), intent(in) :: vel, vel_prev, q_target, q, q_prev
    real(rk) :: s

    if (iter == 1) then
      s = merge(rate_scale_first_dwell, rate_scale_first, dwell)
    else
      s = (1.0_rk - vel_prev/vel)*(q_target - q_prev)/(q - q_prev) + vel_prev/vel
    end if

  end function bc_rate_scale

  !> Lower bound of the perturbations of iteration (B)'s Newton loop (v2-fixed C60, user
  !> decision Q-T1): min_pert_frac times the largest elastic estimate of the face velocities
  !> (bc_elastic_velocities) for the normal stress changes of the run's largest target
  !> magnitudes of each component over the increment dt.
  !>
  !> Arguments:
  !>   stress_max     (in)  real(rk)(3)  largest target magnitude of each component, as a
  !>                                     stress (loads over the original face areas) [MPa]
  !>   length         (in)  real(rk)(3)  lengths along x, y, z [mm]
  !>   e, nu          (in)  real(rk)     elastic moduli (bc_elastic_moduli) [MPa], [-]
  !>   dt             (in)  real(rk)     the increment, > 0 [s]
  !>   min_pert_frac  (in)  real(rk)     the fraction [-]
  !>   floor          (out, result) real(rk)  the bound [mm/s]
  pure function bc_rate_pert_floor(stress_max, length, e, nu, dt, min_pert_frac) result(floor)

    real(rk), intent(in) :: stress_max(3), length(3), e, nu, dt, min_pert_frac
    real(rk) :: floor

    floor = min_pert_frac*maxval(abs(bc_elastic_velocities(stress_max/dt, length, e, nu)))

  end function bc_rate_pert_floor

  !> First perturbation size of iteration (B)'s Newton loop (v2 driver_triaxclr, before
  !> bc_iteration_2): 1e-6 mm/s at a dwell start, else max(0.1 |V - V_init|, 0.01 |V|) of the
  !> control face velocity V and its value V_init at the start of the increment; at least the
  !> lower bound (bc_rate_pert_floor).
  !>
  !> Arguments:
  !>   vel          (in)  real(rk)  control face velocity [mm/s]
  !>   vel_init     (in)  real(rk)  its value at the start of the increment [mm/s]
  !>   dwell_start  (in)  logical   the increment starts a dwell episode
  !>   floor        (in)  real(rk)  the lower bound [mm/s]
  !>   delta        (out, result) real(rk)  the size [mm/s]
  pure function bc_rate_pert_initial(vel, vel_init, dwell_start, floor) result(delta)

    real(rk), intent(in) :: vel, vel_init, floor
    logical, intent(in) :: dwell_start
    real(rk) :: delta

    if (dwell_start) then
      delta = rate_pert_dwell
    else
      delta = max(rate_pert_change_frac*abs(vel - vel_init), rate_pert_vel_frac*abs(vel))
    end if
    delta = max(delta, floor)

  end function bc_rate_pert_initial

  !> Next perturbation size of iteration (B)'s Newton loop (v2 driver_triaxclr,
  !> bc_iteration_2): the largest velocity change of the last correction, at least the lower
  !> bound (v2-fixed C60).
  !>
  !> Arguments:
  !>   dvel   (in)  real(rk)(:)  the last velocity changes [mm/s]
  !>   floor  (in)  real(rk)     the lower bound [mm/s]
  !>   delta  (out, result) real(rk)  the size [mm/s]
  pure function bc_rate_pert_update(dvel, floor) result(delta)

    real(rk), intent(in) :: dvel(:), floor
    real(rk) :: delta

    delta = max(maxval(abs(dvel)), floor)

  end function bc_rate_pert_update

end module bc_iteration_mod
