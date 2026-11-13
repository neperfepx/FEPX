! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Slip-system hardening of the crystal plasticity model (ROADMAP steps 6.4, 7.1, 7.2; D51, D59,
!> D60):
!> the laws, their modifiers, and a general time integration of the strengths over an increment,
!> per point (pure kernels).
!>
!> Theory : dev/theory/addendum_hardening.md (approved 2026-10-10) and
!>          dev/theory/addendum_hardening_modifiers.md (approved 2026-10-11):
!>   saturation_isotropic (§2, v2's base model): one increment for all systems, set by slip system
!>     1: crss_dot^a = h_0 ((tau_s - crss^1)/(tau_s - tau_0^1))^n drive for every a, zero once
!>     crss^1 >= tau_s,
!>     drive = Gammadot = sum_b |gammadot^b|; h_0, tau_s, n are the values of system 1;
!>   saturation_anisotropic (§3, v2's anisotropic hardening): every system its own rate,
!>     crss_dot^a = h_0^a ((tau_s^a - crss^a)/(tau_s^a - tau_0^a))^(n^a) drive^a, zero once
!>     crss^a >= tau_s^a,
!>     drive^a = sum_b H_ab |gammadot^b| with the interaction matrix H (default 1 on the diagonal,
!>     1.4 off it, or the user's matrix);
!>   kocks_mecking (step 7.2, dev/theory/addendum_kocks_mecking.md, approved 2026-10-12, D60):
!>     one network dislocation density km_rho per point,
!>     crss^a = tau_0^a + b mu_0 sqrt(h_n km_rho)                                           (K1)
!>     km_rho dot = Gammadot (k_1 sqrt(km_rho) - k_2 km_rho),  km_rho(0) = rho_0             (K2)
!>     k_2 = k_20 (gammadot_k0/Gammadot)^(1/n_km)                (K3; errata T8)
!>     (Barton et al. 2013, Eqs. 3 and 6, in the physical density of Beyerlein & Tome 2008,
!>     Eq. 3.6); mu_0 is Hill's average shear modulus of the phase (set by the caller). The
!>     stored state is km_rho; the update works in s = sqrt(km_rho), in which (K2) is
!>     s dot = (Gammadot/2)(k_1 - k_2 s) (addendum §3.2; Marin 2006, SAND2006-4170, Eq. 2.45),
!>     linear at fixed slip rates; Gammadot k_2 is zero when Gammadot = 0 (addendum §2);
!>   irradiation softening (step 7.3, dev/theory/addendum_irradiation.md, approved 2026-10-12,
!>     D61; only with kocks_mecking): a damage tensor H (symmetric, Mandel 6-vector, crystal
!>     frame; the field dmg_c) adds a defect term to (K1),
!>     crss^a = tau_0^a + b mu_0 (sqrt(h_n km_rho) + sqrt(h_d <N^a : H>_+))                  (I1)
!>     H dot = -eta sum_a (N^a : H) N^a |gammadot^a|                                       (I3)
!>     H_0 = 3 d_l rho_l (1/N_dirs) sum_i (I - n_i (x) n_i)                            (I4)
!>     (Barton et al. 2013, Eqs. 3-5); N^a = n^a (x) n^a with the unit slip-plane normal. (I3)
!>     is Barton's Eq. 5 as printed; it does not keep H positive semidefinite, so N^a : H can
!>     become negative and <x>_+ = max(x, 0) keeps the defect term real (errata T10, addendum
!>     §3). (I3) is linear in H at fixed slip rates: backward Euler's first Newton step is
!>     exact (6 components, full symmetric Jacobian; time_integration_mod unchanged);
!>   none: constant strengths.
!>   Modifiers (addendum_hardening_modifiers.md):
!>   - saturation_evolution (§4): one tau_s = tau_s0 (Gammadot/gammadot_s0)^m' for every system
!>     (tau_s0 = tau_s(1); tau_s0 when Gammadot = 0), from the slip rates of the current stress,
!>     so constant during the backward-Euler solve (zero derivative with respect to crss);
!>   - cyclic hardening (§3, (H-C1)-(H-C7)): only systems whose accumulated slip A^b has reached
!>     the threshold Dgamma_c drive hardening, and the drive is clamped to [1e-6, 10] d_eq
!>     (d_eq: equivalent deformation rate of the point): isotropic law, f = clamp(sum over
!>     active b of |gammadot^b|) with Dgamma_c = a (crss^1/tau_s)^c for every b (v2); anisotropic
!>     law, f^a = clamp(sum over active b of H_ab |gammadot^b|) with each system's own
!>     Dgamma_c^b = a (crss^b/tau_s^b)^c (new, approved 2026-10-11); a threshold below 2.22e-16 is
!>     zero (v2). The active set is held fixed in the derivative (v2);
!>   - precipitate strengthening (§5) acts on tau_0 at setup only (config_mod), not here.
!>   Integration (addendum_hardening.md §5): the saturation laws evolve the strengths themselves
!>   (x = crss): forward-Euler estimate, then backward Euler
!>   r(crss) = crss - crss_prev - dt crss_dot(crss) = 0 at fixed slip rates, both by the general
!>   methods of time_integration_mod (which do not know what x is); this module supplies the
!>   rates and their N x N Jacobian d(crss_dot)/d(crss).
!>   Kocks-Mecking: the law integrates s = sqrt(km_rho), not the strengths, by the same general
!>   methods (forward-Euler estimate, then backward Euler at fixed slip rates; one component,
!>   tolerance tol_hardening s_prev, addendum §3.4), and the strengths follow from (K1).
!>   Per-point values (§6, §8): tau_0 is an argument; tau_s, h_0, the accumulated slip and d_eq of
!>   the point come in a hardening_point_t (per-element values from simulation.opt; the cyclic
!>   inputs), set in place by hardening_point_set (the phase's values of hardening_params_t
!>   where none are given; then cyclic hardening sees no accumulated slip).
!>   Work arrays: every array the per-point routines need (the point's values, the rates and
!>   their Jacobian, the backward-Euler solver states) is in a hardening_work_t allocated once
!>   at the exact size by hardening_work_alloc, outside the point loops (ROADMAP D100).
!> Source : port of v2-fixed hardening_mod3.f90::isotropic_hardening (C38, C40 fixed),
!>          ::anisotropic_hardening with hardening_mod4.f90::calculate_shrate (sum_b H_ab
!>          |gammadot^b|), ::cyclic_hardening (C46 fixed; the anisotropic form is new), and of
!>          hardening_mod1.f90::update_hardening (saturation evolution, C50 fixed; the forward-
!>          Euler estimate, then backward Euler; toler_hard 1e-4 times tau_0 of each system; at most
!>          point_max_iters_hardening iterations, v2 max_iter_hard_limit 10; the solver and its
!>          deviations: time_integration_mod). A rate is evaluated only below saturation (v2
!>          multiplies a NaN by zero for a non-integer n above it). Kocks-Mecking is not in v2
!>          (addendum_kocks_mecking.md §1): Barton, Arsenlis & Marian (2013), J. Mech. Phys.
!>          Solids 61, 341, Eqs. 3, 6; Beyerlein & Tome (2008), Int. J. Plast. 24, 867, Eq. 3.6.
!>
!> Public parameters:
!>   hardening_none, hardening_isotropic, hardening_anisotropic, hardening_kocks_mecking
!>                        law codes
!>   default_latent       the off-diagonal entry of the default interaction matrix (1.4)
!>
!> Public types:
!>   hardening_params_t      parameters of a phase (components at the type)
!>   hardening_point_t       per-point values (components at the type)
!>   hardening_rates_work_t  work arrays of hardening_rates (components at the type)
!>   hardening_work_t        work arrays of the per-point routines (components at the type)
!>
!> Public routines:
!>   hardening_default_matrix  the default interaction matrix (1 diagonal, 1.4 off it)
!>   hardening_point_alloc     allocate per-point values for N slip systems
!>   hardening_point_set       per-point values from the phase's parameters, in place
!>   hardening_work_alloc      allocate the work arrays of a phase
!>   hardening_rates           hardening rates and their Jacobian d(crss_dot)/d(crss)
!>   hardening_forward         forward-Euler strengths (the estimate)
!>   hardening_backward        backward-Euler strengths (Newton with line search)
!>   hardening_km_strengths    Kocks-Mecking strengths from the density (K1)
!>   hardening_km_k2           Kocks-Mecking recovery coefficient k_2 (K3)
!>   hardening_km_rate         Kocks-Mecking density rate rho_n dot (K2)
!>   hardening_km_forward      Kocks-Mecking forward-Euler density (the estimate)
!>   hardening_km_backward     Kocks-Mecking backward-Euler density (in s = sqrt(km_rho))
module hardening_mod

  use kinds_mod, only: rk
  use tensor_mod, only: mandel_from_matrix
  use time_integration_mod, only: backward_euler_t, forward_euler, backward_euler_alloc, &
    & backward_euler_start, backward_euler_step, backward_euler_finish

  implicit none

  private

  public :: hardening_none, hardening_isotropic, hardening_anisotropic, hardening_kocks_mecking
  public :: default_latent
  public :: hardening_params_t, hardening_point_t, hardening_rates_work_t, hardening_work_t, &
    & hardening_default_matrix, hardening_point_alloc, hardening_point_set, hardening_work_alloc, &
    & hardening_rates, hardening_forward, hardening_backward, &
    & hardening_km_strengths, hardening_km_k2, hardening_km_rate, hardening_km_forward, &
    & hardening_km_backward
  public :: hardening_irr_t, hardening_irr_loop_average, hardening_irr_planes, hardening_irr_dmg_0, &
    & hardening_irr_strengths, hardening_irr_rate, hardening_irr_forward, hardening_irr_backward

  ! Law codes.
  integer, parameter :: hardening_none = 0, hardening_isotropic = 1, hardening_anisotropic = 2, &
    & hardening_kocks_mecking = 3

  ! Default latent-hardening entry of the interaction matrix (user, 2026-10-10).
  real(rk), parameter :: default_latent = 1.4_rk

  ! Cyclic hardening (v2 cyclic_hardening): thresholds below this are zero (v2 mach_eps), and
  ! the clamp of the drive relative to d_eq (v2 shr_min, shr_max).
  real(rk), parameter :: cyclic_crit_floor = 2.22e-16_rk
  real(rk), parameter :: cyclic_drive_min = 1.0e-6_rk, cyclic_drive_max = 10.0_rk

  !> Hardening parameters of a phase (one value per slip system, in the slip-system order).
  !>
  !> Components:
  !>   law          integer               hardening_none, _isotropic, _anisotropic or
  !>                                      _kocks_mecking
  !>   h_0          real(rk)(N)           hardening rate coefficient [MPa]
  !>   tau_s        real(rk)(N)           saturation strength, or tau_s0 with saturation
  !>                                      evolution (one value then) [MPa]
  !>   n            real(rk)(N)           hardening exponent [-]
  !>   interaction  real(rk)(N, N)        interaction matrix H (anisotropic) [-]
  !>   saturation_evolution  logical      tau_s = tau_s0 (Gammadot/gammadot_s0)^m_prime
  !>                                      (saturation_evolution on)
  !>   gammadot_s0  real(rk)              reference slip rate of saturation evolution [1/s]
  !>   m_prime      real(rk)              exponent of saturation evolution [-]
  !>   cyclic       logical               cyclic hardening (cyclic_hardening on)
  !>   cyclic_a     real(rk)              threshold coefficient a of cyclic hardening [-]
  !>   cyclic_c     real(rk)              threshold exponent c of cyclic hardening [-]
  !>   tol_hardening real(rk)             tolerance of the strength update (backward Euler,
  !>                                      addendum_hardening.md H-B1): |r^a| < tol_hardening
  !>                                      tau_0^a for every system (keyword point_tol_hardening;
  !>                                      v2 toler_hard) [-]
  !>   max_iters    integer               backward-Euler iteration limit
  !>   Kocks-Mecking (law hardening_kocks_mecking; addendum_kocks_mecking.md §4):
  !>   km_b         real(rk)              Burgers vector b [mm]
  !>   km_mu_0      real(rk)              shear modulus mu_0: Hill's average of the phase's
  !>                                      elastic constants (crystal_mod::shear_modulus_hill,
  !>                                      set by the caller) [MPa]
  !>   km_h_n       real(rk)              network hardening coefficient h_n [-]
  !>   km_rho_0     real(rk)              initial network dislocation density rho_0 [1/mm^2]
  !>   km_k1        real(rk)              storage coefficient k_1 [1/mm]
  !>   km_k20       real(rk)              recovery coefficient k_20 [-]
  !>   km_gammadot_k0 real(rk)            reference slip rate of k_2 [1/s]
  !>   km_n         real(rk)              exponent n_km of k_2 (power 1/n_km) [-]
  type :: hardening_params_t
    integer :: law = hardening_none
    real(rk), allocatable :: h_0(:)
    real(rk), allocatable :: tau_s(:)
    real(rk), allocatable :: n(:)
    real(rk), allocatable :: interaction(:, :)
    logical :: saturation_evolution = .false.
    real(rk) :: gammadot_s0 = 1.0_rk
    real(rk) :: m_prime = 0.0_rk
    logical :: cyclic = .false.
    real(rk) :: cyclic_a = 0.0_rk
    real(rk) :: cyclic_c = 0.0_rk
    real(rk) :: tol_hardening = 1.0e-4_rk
    integer :: max_iters = 10
    real(rk) :: km_b = 0.0_rk
    real(rk) :: km_mu_0 = 0.0_rk
    real(rk) :: km_h_n = 0.0_rk
    real(rk) :: km_rho_0 = 0.0_rk
    real(rk) :: km_k1 = 0.0_rk
    real(rk) :: km_k20 = 0.0_rk
    real(rk) :: km_gammadot_k0 = 1.0_rk
    real(rk) :: km_n = 1.0_rk
  end type hardening_params_t

  !> Per-point hardening values, every array of size N (hardening_point_alloc, then
  !> hardening_point_set).
  !>
  !> Components:
  !>   num_slip    integer                 number of slip systems N
  !>   tau_s       real(rk)(N)             saturation strength (or tau_s0) of the point's element
  !>                                       [MPa]
  !>   h_0         real(rk)(N)             hardening rate coefficient of the element [MPa]
  !>   slip_accum  real(rk)(N)             accumulated slip |gammadot| dt since the last load
  !>                                       reversal, start of the increment (cyclic) [-]
  !>   defrate_eq  real(rk)                equivalent deformation rate sqrt(2/3 d':d') of the
  !>                                       point (cyclic clamp) [1/s]
  type :: hardening_point_t
    integer :: num_slip = 0
    real(rk), allocatable :: tau_s(:)
    real(rk), allocatable :: h_0(:)
    real(rk), allocatable :: slip_accum(:)
    real(rk) :: defrate_eq = 0.0_rk
  end type hardening_point_t

  !> Work arrays of hardening_rates, every array of size N (hardening_work_alloc).
  !>
  !> Components:
  !>   tau_s       real(rk)(N)  saturation strength of each system [MPa]
  !>   active      real(rk)(N)  1 for a system that drives hardening, else 0 (cyclic) [-]
  !>   abs_rate    real(rk)(N)  |gammadot^b| of the systems that drive hardening [1/s]
  !>   drive       real(rk)(N)  driving slip rate of each system [1/s]
  !>   crss_ratio  real(rk)(N)  crss/tau_s that sets each system's cyclic threshold [-]
  !>   crit        real(rk)(N)  cyclic thresholds Dgamma_c [-]
  type :: hardening_rates_work_t
    real(rk), allocatable :: tau_s(:)
    real(rk), allocatable :: active(:)
    real(rk), allocatable :: abs_rate(:)
    real(rk), allocatable :: drive(:)
    real(rk), allocatable :: crss_ratio(:)
    real(rk), allocatable :: crit(:)
  end type hardening_rates_work_t

  !> Work arrays of the per-point hardening routines of a phase, at the exact size, allocated
  !> once outside the point loops (hardening_work_alloc; ROADMAP D100).
  !>
  !> Components:
  !>   num_slip    integer                 number of slip systems N
  !>   point       type(hardening_point_t) the point's values, set by the caller
  !>                                       (hardening_point_set) before the point's updates
  !>   rate        real(rk)(N)             hardening rates (saturation laws) [MPa/s]
  !>   jac         real(rk)(N, N)          their Jacobian d crss_dot / d crss [1/s]
  !>   toler       real(rk)(N)             backward-Euler tolerances of the strengths [MPa]
  !>   rates       type(hardening_rates_work_t)  hardening_rates' work arrays
  !>   be_crss     type(backward_euler_t)  solver state of the strengths (saturation laws; N)
  !>   be_density  type(backward_euler_t)  solver state of sqrt(km_rho) (Kocks-Mecking; 1)
  !>   be_damage   type(backward_euler_t)  solver state of the damage tensor (irradiation; 6)
  type :: hardening_work_t
    integer :: num_slip = 0
    type(hardening_point_t) :: point
    real(rk), allocatable :: rate(:)
    real(rk), allocatable :: jac(:, :)
    real(rk), allocatable :: toler(:)
    type(hardening_rates_work_t) :: rates
    type(backward_euler_t) :: be_crss
    type(backward_euler_t) :: be_density
    type(backward_euler_t) :: be_damage
  end type hardening_work_t

  !> Irradiation-softening parameters of a phase (with kocks_mecking; addendum_irradiation.md
  !> §6), a type of its own beside hardening_params_t.
  !>
  !> Components:
  !>   h_d            real(rk)            defect hardening coefficient h_d [-]
  !>   eta            real(rk)            annihilation efficiency eta [-]
  !>   loop_density   real(rk)            loop number density rho_l of the phase [1/mm^3]
  !>   loop_diameter  real(rk)            loop diameter d_l of the phase [mm]
  !>   loop_avg_c     real(rk)(6)         (1/N_dirs) sum_i (I - n_i (x) n_i) over the loop
  !>                                      normals, Mandel, crystal frame [-]
  !>   plane_c        real(rk)(6, N)      N^a = n^a (x) n^a of every slip system, Mandel,
  !>                                      crystal frame [-]
  type :: hardening_irr_t
    real(rk) :: h_d = 0.0_rk
    real(rk) :: eta = 0.0_rk
    real(rk) :: loop_density = 0.0_rk
    real(rk) :: loop_diameter = 0.0_rk
    real(rk) :: loop_avg_c(6) = 0.0_rk
    real(rk), allocatable :: plane_c(:, :)
  end type hardening_irr_t

contains

  !> The default interaction matrix: 1 on the diagonal, default_latent off it (addendum §3).
  !>
  !> Arguments:
  !>   num   (in)  integer  number of slip systems N
  !>   mat   (out, result) real(rk)(N, N)  the matrix
  !> Locals:
  !>   a     integer  system
  pure function hardening_default_matrix(num) result(mat)

    integer, intent(in) :: num
    real(rk) :: mat(num, num)

    integer :: a

    mat = default_latent
    do a = 1, num
      mat(a, a) = 1.0_rk
    end do

  end function hardening_default_matrix

  !> Allocate per-point values for N slip systems (once, outside the point loops); values
  !> already of size N are kept.
  !>
  !> Arguments:
  !>   point     (inout) type(hardening_point_t)  per-point values
  !>   num_slip  (in)    integer                  number of slip systems N
  pure subroutine hardening_point_alloc(point, num_slip)

    type(hardening_point_t), intent(inout) :: point
    integer, intent(in) :: num_slip

    ! Nothing to do when the size is already N.
    if (point%num_slip == num_slip .and. allocated(point%tau_s)) return

    ! The arrays at the new size.
    if (allocated(point%tau_s)) deallocate (point%tau_s, point%h_0, point%slip_accum)
    point%num_slip = num_slip
    allocate (point%tau_s(num_slip), point%h_0(num_slip), point%slip_accum(num_slip))

  end subroutine hardening_point_alloc

  !> Per-point values from the phase's parameters, with the given per-element values in their
  !> place (simulation.opt) and the cyclic inputs of the point, set in place (every value is
  !> set: zero where neither the phase nor the arguments give one).
  !>
  !> Arguments:
  !>   point       (inout) type(hardening_point_t)  the values (hardening_point_alloc)
  !>   hard        (in)  type(hardening_params_t)  phase parameters (law hardening_none: only
  !>                                               the cyclic inputs are set)
  !>   tau_s       (in, optional) real(rk)(N)  element's saturation strength (or tau_s0) [MPa]
  !>   h_0         (in, optional) real(rk)(N)  element's hardening rate coefficient [MPa]
  !>   slip_accum  (in, optional) real(rk)(N)  accumulated slip of the point [-]
  !>   defrate_eq  (in, optional) real(rk)     equivalent deformation rate of the point [1/s]
  !> Locals:
  !>   n           integer                     number of slip systems N
  pure subroutine hardening_point_set(point, hard, tau_s, h_0, slip_accum, defrate_eq)

    type(hardening_point_t), intent(inout) :: point
    type(hardening_params_t), intent(in) :: hard
    real(rk), intent(in), optional :: tau_s(:)
    real(rk), intent(in), optional :: h_0(:)
    real(rk), intent(in), optional :: slip_accum(:)
    real(rk), intent(in), optional :: defrate_eq

    integer :: n

    ! Zero, then the phase values, then the element's and the point's.
    n = point%num_slip
    point%tau_s(:) = 0.0_rk
    point%h_0(:) = 0.0_rk
    point%slip_accum(:) = 0.0_rk
    point%defrate_eq = 0.0_rk
    if (allocated(hard%tau_s)) then
      if (size(hard%tau_s) == n) point%tau_s(:) = hard%tau_s
    end if
    if (allocated(hard%h_0)) then
      if (size(hard%h_0) == n) point%h_0(:) = hard%h_0
    end if
    if (present(tau_s)) point%tau_s(:) = tau_s
    if (present(h_0)) point%h_0(:) = h_0
    if (present(slip_accum)) point%slip_accum(:) = slip_accum
    if (present(defrate_eq)) point%defrate_eq = defrate_eq

  end subroutine hardening_point_set

  !> Allocate the work arrays of a phase's per-point hardening routines at the exact size, once,
  !> outside the point loops (ROADMAP D100): the point's values (set to the phase's), the
  !> rates, and the backward-Euler solver state of the variable the law evolves (the strengths
  !> of the saturation laws; sqrt(km_rho) of Kocks-Mecking; the damage tensor with irradiation).
  !>
  !> Arguments:
  !>   work         (inout) type(hardening_work_t)    work arrays
  !>   hard         (in)    type(hardening_params_t)  phase parameters (law)
  !>   num_slip     (in)    integer                   number of slip systems N
  !>   irradiation  (in, optional) logical            irradiation softening (default: no)
  !> Locals:
  !>   n            integer                           N
  pure subroutine hardening_work_alloc(work, hard, num_slip, irradiation)

    type(hardening_work_t), intent(inout) :: work
    type(hardening_params_t), intent(in) :: hard
    integer, intent(in) :: num_slip
    logical, intent(in), optional :: irradiation

    integer :: n

    ! The point's values (the phase's until the caller sets the point's own).
    n = num_slip
    call hardening_point_alloc(work%point, n)
    call hardening_point_set(work%point, hard)

    ! The arrays of the rates.
    if (work%num_slip /= n .or. .not. allocated(work%rate)) then
      if (allocated(work%rate)) then
        deallocate (work%rate, work%jac, work%toler, work%rates%tau_s, work%rates%active, &
          & work%rates%abs_rate, work%rates%drive, work%rates%crss_ratio, work%rates%crit)
      end if
      allocate (work%rate(n), work%jac(n, n), work%toler(n), work%rates%tau_s(n), &
        & work%rates%active(n), work%rates%abs_rate(n), work%rates%drive(n), &
        & work%rates%crss_ratio(n), work%rates%crit(n))
    end if
    work%num_slip = n

    ! The solver state of the evolved variable.
    select case (hard%law)
    case (hardening_isotropic, hardening_anisotropic)
      call backward_euler_alloc(work%be_crss, n)
    case (hardening_kocks_mecking)
      call backward_euler_alloc(work%be_density, 1)
    case default
      continue
    end select
    if (present(irradiation)) then
      if (irradiation) call backward_euler_alloc(work%be_damage, 6)
    end if

  end subroutine hardening_work_alloc

  !> Hardening rates of every system and their Jacobian jac(a, b) = d crss_dot^a / d crss^b, at
  !> given strengths and slip rates (module header).
  !>
  !> Theory : addendum_hardening.md §2-3; addendum_hardening_modifiers.md §3-4.
  !> Source : port of v2-fixed isotropic_hardening, anisotropic_hardening, calculate_shrate,
  !>          cyclic_hardening, and update_hardening's saturation evolution.
  !>
  !> Arguments:
  !>   hard      (in)  type(hardening_params_t)  parameters
  !>   crss      (in)  real(rk)(N)     strengths [MPa]
  !>   tau_0     (in)  real(rk)(N)     initial strengths [MPa]
  !>   sliprate  (in)  real(rk)(N)     slip rates [1/s]
  !>   rate      (out) real(rk)(N)     crss_dot [MPa/s]
  !>   jac       (out) real(rk)(N, N)  d crss_dot / d crss [1/s]
  !>   point     (in)  type(hardening_point_t)  per-point values (hardening_point_set; the
  !>                                   phase's: tau_s and h_0, no accumulated slip, d_eq 0)
  !>   work      (inout) type(hardening_rates_work_t)  work arrays (hardening_work_alloc):
  !>                                   tau_s, the saturation strength of each system [MPa];
  !>                                   active, 1 for a system that drives hardening, else 0
  !>                                   (cyclic) [-]; abs_rate, |gammadot^b| of the systems that
  !>                                   drive hardening [1/s]; drive, the driving slip rate of
  !>                                   each system [1/s]: Gammadot or f (isotropic), sum_b H_ab
  !>                                   (active) |gammadot^b| (anisotropic); crss_ratio and crit,
  !>                                   cyclic_active's inputs
  !> Locals:
  !>   n         integer               number of systems
  !>   a         integer               system
  !>   shrate    real(rk)              Gammadot = sum_b |gammadot^b| [1/s]
  !>   ratio     real(rk)              (tau_s - crss)/(tau_s - tau_0) [-]
  pure subroutine hardening_rates(hard, crss, tau_0, sliprate, rate, jac, point, work)

    type(hardening_params_t), intent(in) :: hard
    real(rk), intent(in) :: crss(:)
    real(rk), intent(in) :: tau_0(:)
    real(rk), intent(in) :: sliprate(:)
    real(rk), intent(out) :: rate(:)
    real(rk), intent(out) :: jac(:, :)
    type(hardening_point_t), intent(in) :: point
    type(hardening_rates_work_t), intent(inout) :: work

    integer :: n, a
    real(rk) :: shrate, ratio

    rate = 0.0_rk
    jac = 0.0_rk
    if (hard%law == hardening_none) return
    n = size(crss)

    ! Saturation strength: the given one per system, or one evolving value (v2
    ! update_hardening: tau_s0 when no slip).
    shrate = sum(abs(sliprate))
    work%tau_s(:) = point%tau_s
    if (hard%saturation_evolution) then
      work%tau_s(:) = point%tau_s(1)
      if (shrate > 0.0_rk) work%tau_s(:) = point%tau_s(1)*(shrate/hard%gammadot_s0)**hard%m_prime
    end if

    select case (hard%law)

    case (hardening_isotropic)
      ! System 1 sets one rate for every system; it depends on crss^1 only (first column).
      if (work%tau_s(1) - crss(1) > 0.0_rk) then
        ! The drive: Gammadot, or the clamped slip rate of the systems past the threshold set
        ! by system 1 (cyclic, v2).
        if (hard%cyclic) then
          work%crss_ratio(:) = crss(1)/work%tau_s(1)
          call cyclic_active(hard, point%slip_accum, work%crss_ratio, work%crit, work%active)
          work%drive(:) = clamp_drive(sum(work%active*abs(sliprate)), point%defrate_eq)
        else
          work%drive(:) = shrate
        end if
        ratio = (work%tau_s(1) - crss(1))/(work%tau_s(1) - tau_0(1))
        rate = point%h_0(1)*ratio**hard%n(1)*work%drive(1)
        jac(:, 1) = -point%h_0(1)*hard%n(1)*ratio**(hard%n(1) - 1.0_rk)*work%drive(1)/ &
          & (work%tau_s(1) - tau_0(1))
      end if

    case (hardening_anisotropic)
      ! Each system its own rate, driven through the interaction matrix (diagonal Jacobian); with
      ! cyclic hardening only the systems past their own threshold drive it, and each drive is
      ! clamped (addendum_hardening_modifiers.md (H-C5)-(H-C7)).
      work%abs_rate(:) = abs(sliprate)
      if (hard%cyclic) then
        work%crss_ratio(:) = crss/work%tau_s
        call cyclic_active(hard, point%slip_accum, work%crss_ratio, work%crit, work%active)
        work%abs_rate(:) = work%abs_rate*work%active
      end if
      work%drive(:) = matmul(hard%interaction, work%abs_rate)
      if (hard%cyclic) then
        do a = 1, n
          work%drive(a) = clamp_drive(work%drive(a), point%defrate_eq)
        end do
      end if
      do a = 1, n
        if (.not. (work%tau_s(a) - crss(a) > 0.0_rk)) cycle
        ratio = (work%tau_s(a) - crss(a))/(work%tau_s(a) - tau_0(a))
        rate(a) = point%h_0(a)*ratio**hard%n(a)*work%drive(a)
        jac(a, a) = -point%h_0(a)*hard%n(a)*ratio**(hard%n(a) - 1.0_rk)*work%drive(a)/ &
          & (work%tau_s(a) - tau_0(a))
      end do

    case default
      continue

    end select

  end subroutine hardening_rates

  !> Cyclic hardening: which systems drive hardening (v2 cyclic_hardening): those whose
  !> accumulated slip has reached the threshold a (crss/tau_s)^c (zero below 2.22e-16).
  !>
  !> Arguments:
  !>   hard        (in)  type(hardening_params_t)  parameters (cyclic_a, cyclic_c)
  !>   slip_accum  (in)  real(rk)(N)  accumulated slip of each system [-]
  !>   crss_ratio  (in)  real(rk)(N)  crss/tau_s that sets each system's threshold (isotropic law:
  !>                                  system 1's for all, v2; anisotropic: each its own) [-]
  !>   crit        (out) real(rk)(N)  work array: the thresholds Dgamma_c [-]
  !>   active      (out) real(rk)(N)  1 for an active system, else 0
  pure subroutine cyclic_active(hard, slip_accum, crss_ratio, crit, active)

    type(hardening_params_t), intent(in) :: hard
    real(rk), intent(in) :: slip_accum(:)
    real(rk), intent(in) :: crss_ratio(:)
    real(rk), intent(out) :: crit(:)
    real(rk), intent(out) :: active(:)

    ! Thresholds (zero below the floor), then the systems that have reached theirs.
    crit = hard%cyclic_a*abs(crss_ratio)**hard%cyclic_c
    where (crit < cyclic_crit_floor) crit = 0.0_rk
    active = merge(1.0_rk, 0.0_rk, slip_accum >= crit)

  end subroutine cyclic_active

  !> Cyclic hardening: the drive clamped to [1e-6, 10] d_eq (v2 cyclic_hardening).
  !>
  !> Arguments:
  !>   drive       (in)  real(rk)  unclamped drive [1/s]
  !>   defrate_eq  (in)  real(rk)  equivalent deformation rate of the point [1/s]
  !>   clamped     (out, result) real(rk)  the clamped drive [1/s]
  pure function clamp_drive(drive, defrate_eq) result(clamped)

    real(rk), intent(in) :: drive
    real(rk), intent(in) :: defrate_eq
    real(rk) :: clamped

    clamped = min(max(drive, cyclic_drive_min*defrate_eq), cyclic_drive_max*defrate_eq)

  end function clamp_drive

  !> Forward-Euler strengths: crss = crss_prev + dt crss_dot(crss_prev) (v2 'euler_forward'; the
  !> saturation laws evolve the strengths themselves, x = crss,
  !> time_integration_mod::forward_euler).
  !>
  !> Arguments:
  !>   hard       (in)  type(hardening_params_t)  parameters
  !>   crss_prev  (in)  real(rk)(N)  strengths at the start of the increment [MPa]
  !>   tau_0      (in)  real(rk)(N)  initial strengths [MPa]
  !>   sliprate   (in)  real(rk)(N)  slip rates [1/s]
  !>   dt         (in)  real(rk)     time increment [s]
  !>   crss       (out) real(rk)(N)  estimated strengths [MPa]
  !>   work       (inout) type(hardening_work_t)  work arrays (hardening_work_alloc), with the
  !>                                  point's values in work%point; the rates and Jacobian in
  !>                                  work%rate, work%jac
  pure subroutine hardening_forward(hard, crss_prev, tau_0, sliprate, dt, crss, work)

    type(hardening_params_t), intent(in) :: hard
    real(rk), intent(in) :: crss_prev(:)
    real(rk), intent(in) :: tau_0(:)
    real(rk), intent(in) :: sliprate(:)
    real(rk), intent(in) :: dt
    real(rk), intent(out) :: crss(:)
    type(hardening_work_t), intent(inout) :: work

    ! The rates at the start strengths, then the general forward Euler.
    call hardening_rates(hard, crss_prev, tau_0, sliprate, work%rate, work%jac, work%point, &
      & work%rates)
    call forward_euler(crss_prev, work%rate, dt, crss)

  end subroutine hardening_forward

  !> Backward-Euler strengths at fixed slip rates: the saturation laws evolve the strengths
  !> themselves (x = crss), so r(crss) = crss - crss_prev - dt crss_dot(crss) = 0 is solved by
  !> the general solver time_integration_mod::backward_euler_* (Newton, v2's line search), with
  !> this law's rates and Jacobian evaluated where the solver asks; converged
  !> when |r^a| < tol_hardening tau_0^a for every system (v2 toler_hard); systems that did not
  !> converge keep crss_prev (v2).
  !>
  !> Arguments:
  !>   hard       (in)    type(hardening_params_t)  parameters
  !>   crss_prev  (in)    real(rk)(N)  strengths at the start of the increment [MPa]
  !>   tau_0      (in)    real(rk)(N)  initial strengths [MPa]
  !>   sliprate   (in)    real(rk)(N)  slip rates [1/s]
  !>   dt         (in)    real(rk)     time increment [s]
  !>   crss       (inout) real(rk)(N)  start estimate in, end-of-increment strengths out [MPa]
  !>   converged  (out)   logical      every system converged
  !>   work       (inout) type(hardening_work_t)  work arrays (hardening_work_alloc), with the
  !>                                  point's values in work%point: the solver state be_crss,
  !>                                  its tolerances toler [MPa], the rates rate [MPa/s] and
  !>                                  their Jacobian jac [1/s] at be_crss%x_eval
  !>   num_iters  (out, optional) integer  Newton iterations used
  pure subroutine hardening_backward(hard, crss_prev, tau_0, sliprate, dt, crss, converged, &
    & work, num_iters)

    type(hardening_params_t), intent(in) :: hard
    real(rk), intent(in) :: crss_prev(:)
    real(rk), intent(in) :: tau_0(:)
    real(rk), intent(in) :: sliprate(:)
    real(rk), intent(in) :: dt
    real(rk), intent(inout) :: crss(:)
    logical, intent(out) :: converged
    type(hardening_work_t), intent(inout) :: work
    integer, intent(out), optional :: num_iters

    ! No hardening: constant strengths.
    if (hard%law == hardening_none) then
      crss = crss_prev
      converged = .true.
      if (present(num_iters)) num_iters = 0
      return
    end if

    ! The solve, with this law's rates where the solver asks for them.
    work%toler(:) = hard%tol_hardening*tau_0
    call backward_euler_start(work%be_crss, crss_prev, crss, dt, work%toler, hard%max_iters)
    do while (.not. work%be_crss%done)
      call hardening_rates(hard, work%be_crss%x_eval, tau_0, sliprate, work%rate, work%jac, &
        & work%point, work%rates)
      call backward_euler_step(work%be_crss, work%rate, work%jac)
    end do
    call backward_euler_finish(work%be_crss, crss, converged, num_iters)

  end subroutine hardening_backward

  !> Kocks-Mecking strengths from the network dislocation density (K1):
  !> crss^a = tau_0^a + b mu_0 sqrt(h_n km_rho).
  !>
  !> Theory : addendum_kocks_mecking.md (K1); Barton et al. (2013) Eq. 3 without the damage term.
  !>
  !> Arguments:
  !>   hard   (in)  type(hardening_params_t)  parameters (km_b, km_mu_0, km_h_n)
  !>   tau_0  (in)  real(rk)(N)  lattice strengths of the point [MPa]
  !>   km_rho  (in)  real(rk)     network dislocation density (> 0) [1/mm^2]
  !>   crss   (out, result) real(rk)(N)  strengths [MPa]
  pure function hardening_km_strengths(hard, tau_0, km_rho) result(crss)

    type(hardening_params_t), intent(in) :: hard
    real(rk), intent(in) :: tau_0(:)
    real(rk), intent(in) :: km_rho
    real(rk) :: crss(size(tau_0))

    ! The same density term for every system.
    crss = tau_0 + hard%km_b*hard%km_mu_0*sqrt(hard%km_h_n*km_rho)

  end function hardening_km_strengths

  !> Kocks-Mecking recovery coefficient k_2 = k_20 (gammadot_k0/Gammadot)^(1/n_km) (K3); zero
  !> when Gammadot = 0 (k_2 undefined, no evolution; addendum_kocks_mecking.md §2, §7).
  !>
  !> Theory : addendum_kocks_mecking.md (K3); Barton et al. (2013) Eq. 6 with errata T8.
  !>
  !> Arguments:
  !>   hard    (in)  type(hardening_params_t)  parameters (km_k20, km_gammadot_k0, km_n)
  !>   shrate  (in)  real(rk)  Gammadot = sum_b |gammadot^b| [1/s]
  !>   k2      (out, result) real(rk)  k_2 [-]
  pure function hardening_km_k2(hard, shrate) result(k2)

    type(hardening_params_t), intent(in) :: hard
    real(rk), intent(in) :: shrate
    real(rk) :: k2

    k2 = 0.0_rk
    if (shrate > 0.0_rk) k2 = hard%km_k20*(hard%km_gammadot_k0/shrate)**(1.0_rk/hard%km_n)

  end function hardening_km_k2

  !> Rate of the network dislocation density (K2), rho_n dot = Gammadot (k_1 sqrt(rho_n) -
  !> k_2 rho_n) = 2 sqrt(rho_n) s dot, evaluated through km_sqrt_rate (Gammadot k_2 = 0 at
  !> Gammadot = 0); the output km_rhorate (step 7.3, D61).
  !>
  !> Theory : addendum_kocks_mecking.md (K2), (K2s); addendum_irradiation.md §8.
  !>
  !> Arguments:
  !>   hard    (in)  hardening_params_t  the law's parameters
  !>   shrate  (in)  real(rk)            Gammadot = sum of |slip rates| [1/s]
  !>   rho_n   (in)  real(rk)            network dislocation density [1/mm^2]
  !>   rate    (out, result) real(rk)    rho_n dot [1/(mm^2 s)]
  !>
  !> Locals:
  !>   s_rate, s_jac  real(rk)  s dot [1/(mm s)] and its derivative (unused) at s = sqrt(rho_n)
  pure function hardening_km_rate(hard, shrate, rho_n) result(rate)

    type(hardening_params_t), intent(in) :: hard
    real(rk), intent(in) :: shrate
    real(rk), intent(in) :: rho_n
    real(rk) :: rate

    real(rk) :: s_rate, s_jac

    ! Rate in s, then the chain rule back to the density.
    call km_sqrt_rate(hard, shrate, sqrt(rho_n), s_rate, s_jac)
    rate = 2.0_rk*sqrt(rho_n)*s_rate

  end function hardening_km_rate

  !> Kocks-Mecking rate in s = sqrt(km_rho) and its derivative: s dot = (Gammadot/2)(k_1 - k_2 s),
  !> d(s dot)/ds = -Gammadot k_2/2, with Gammadot k_2 = 0 when Gammadot = 0.
  !>
  !> Theory : addendum_kocks_mecking.md (K2s), §3.2 (chain rule s dot = km_rho dot/(2 sqrt(km_rho));
  !>          Marin (2006), SAND2006-4170, Eq. 2.45 with epsilon = b s).
  !>
  !> Arguments:
  !>   hard    (in)  type(hardening_params_t)  parameters
  !>   shrate  (in)  real(rk)  Gammadot [1/s]
  !>   s       (in)  real(rk)  sqrt(km_rho) [1/mm]
  !>   rate    (out) real(rk)  s dot [1/(mm s)]
  !>   jac     (out) real(rk)  d(s dot)/ds [1/s]
  pure subroutine km_sqrt_rate(hard, shrate, s, rate, jac)

    type(hardening_params_t), intent(in) :: hard
    real(rk), intent(in) :: shrate
    real(rk), intent(in) :: s
    real(rk), intent(out) :: rate
    real(rk), intent(out) :: jac

    ! Storage and recovery; Gammadot k_2 from k_2 (zero at Gammadot = 0).
    rate = 0.5_rk*shrate*(hard%km_k1 - hardening_km_k2(hard, shrate)*s)
    jac = -0.5_rk*shrate*hardening_km_k2(hard, shrate)

  end subroutine km_sqrt_rate

  !> Kocks-Mecking forward-Euler density (the estimate of the point solve): in s = sqrt(km_rho),
  !> s = s_prev + dt s dot(s_prev) (time_integration_mod::forward_euler); a non-positive
  !> estimate (slip increment above 2/k_2) is replaced by s_prev (addendum §3.3: a start value
  !> only); km_rho = s^2.
  !>
  !> Arguments:
  !>   hard      (in)  type(hardening_params_t)  parameters
  !>   rho_prev  (in)  real(rk)     density at the start of the increment (> 0) [1/mm^2]
  !>   sliprate  (in)  real(rk)(N)  slip rates [1/s]
  !>   dt        (in)  real(rk)     time increment [s]
  !>   rho       (out) real(rk)     estimated density [1/mm^2]
  !> Locals:
  !>   s_prev, s  real(rk)(1)  sqrt of the densities [1/mm]
  !>   rate, jac  real(rk)(1)  s dot [1/(mm s)] and its derivative [1/s]
  pure subroutine hardening_km_forward(hard, rho_prev, sliprate, dt, rho)

    type(hardening_params_t), intent(in) :: hard
    real(rk), intent(in) :: rho_prev
    real(rk), intent(in) :: sliprate(:)
    real(rk), intent(in) :: dt
    real(rk), intent(out) :: rho

    real(rk) :: s_prev(1), s(1), rate(1), jac(1)

    ! Forward Euler in s, then back to the density.
    s_prev(1) = sqrt(rho_prev)
    call km_sqrt_rate(hard, sum(abs(sliprate)), s_prev(1), rate(1), jac(1))
    call forward_euler(s_prev, rate, dt, s)
    if (.not. s(1) > 0.0_rk) s = s_prev
    rho = s(1)**2

  end subroutine hardening_km_forward

  !> Kocks-Mecking backward-Euler density at fixed slip rates: in s = sqrt(km_rho),
  !> r(s) = s - s_prev - dt s dot(s) = 0, solved by the general solver
  !> time_integration_mod::backward_euler_* (exact in its first Newton step, the equation being
  !> linear in s); converged when |r| < tol_hardening s_prev (addendum §3.4); a non-converged
  !> update keeps rho_prev (as the solver's other users).
  !>
  !> Arguments:
  !>   hard       (in)    type(hardening_params_t)  parameters
  !>   rho_prev   (in)    real(rk)     density at the start of the increment (> 0) [1/mm^2]
  !>   sliprate   (in)    real(rk)(N)  slip rates [1/s]
  !>   dt         (in)    real(rk)     time increment [s]
  !>   rho        (inout) real(rk)     start estimate in, end-of-increment density out [1/mm^2]
  !>   converged  (out)   logical      the update converged
  !>   work       (inout) type(hardening_work_t)  work arrays (hardening_work_alloc): the solver
  !>                                  state be_density
  !>   num_iters  (out, optional) integer  Newton iterations used
  !> Locals:
  !>   s_prev     real(rk)(1)             sqrt(rho_prev) [1/mm]
  !>   s          real(rk)(1)             the solution in s [1/mm]
  !>   shrate     real(rk)                Gammadot [1/s]
  !>   rate, jac  real(rk)(1), (1, 1)     s dot and its derivative at be%x_eval
  pure subroutine hardening_km_backward(hard, rho_prev, sliprate, dt, rho, converged, work, &
    & num_iters)

    type(hardening_params_t), intent(in) :: hard
    real(rk), intent(in) :: rho_prev
    real(rk), intent(in) :: sliprate(:)
    real(rk), intent(in) :: dt
    real(rk), intent(inout) :: rho
    logical, intent(out) :: converged
    type(hardening_work_t), intent(inout) :: work
    integer, intent(out), optional :: num_iters

    real(rk) :: s_prev(1), s(1), shrate, rate(1), jac(1, 1)

    ! The solve in s, from the start estimate, with the rate evaluated where the solver asks.
    shrate = sum(abs(sliprate))
    s_prev(1) = sqrt(rho_prev)
    call backward_euler_start(work%be_density, s_prev, [sqrt(rho)], dt, &
      & hard%tol_hardening*s_prev, hard%max_iters)
    do while (.not. work%be_density%done)
      call km_sqrt_rate(hard, shrate, work%be_density%x_eval(1), rate(1), jac(1, 1))
      call backward_euler_step(work%be_density, rate, jac)
    end do
    call backward_euler_finish(work%be_density, s, converged, num_iters)

    ! Back to the stored density.
    rho = s(1)**2

  end subroutine hardening_km_backward

  !> Average loop tensor of equally populated loop normals, (1/N_dirs) sum_i (I - n_i (x) n_i)
  !> (addendum_irradiation.md (I4); Barton et al. 2013 Eq. 4; Jerry's note Eq. 27).
  !>
  !> Arguments:
  !>   normals_c  (in)  real(rk)(3, N_dirs)  loop normals (any length; normalized here),
  !>                                         crystal frame
  !>   avg_c      (out, result) real(rk)(6)  the average, Mandel, crystal frame [-]
  !>
  !> Locals:
  !>   i, k  integer        loop-normal index; diagonal index
  !>   unit  real(rk)(3)    the normal scaled to unit length
  !>   mat   real(rk)(3,3)  I - n (x) n
  pure function hardening_irr_loop_average(normals_c) result(avg_c)

    real(rk), intent(in) :: normals_c(:, :)
    real(rk) :: avg_c(6)

    integer :: i, k
    real(rk) :: unit(3), mat(3, 3)

    ! Sum I - n (x) n over the normals, then divide by their number.
    avg_c = 0.0_rk
    do i = 1, size(normals_c, 2)
      unit = normals_c(:, i)/norm2(normals_c(:, i))
      mat = -spread(unit, 2, 3)*spread(unit, 1, 3)
      do k = 1, 3
        mat(k, k) = mat(k, k) + 1.0_rk
      end do
      avg_c = avg_c + mandel_from_matrix(mat)
    end do
    avg_c = avg_c/real(size(normals_c, 2), rk)

  end function hardening_irr_loop_average

  !> N^a = n^a (x) n^a of every slip system, Mandel, crystal frame (addendum_irradiation.md
  !> (I2)), so that N^a : H is the dot product of the Mandel vectors.
  !>
  !> Arguments:
  !>   slip_normal_c  (in)  real(rk)(3, N)  unit slip-plane normals, crystal frame
  !>   plane_c        (out, result) real(rk)(6, N)  N^a, Mandel, crystal frame [-]
  !>
  !> Locals:
  !>   a  integer  slip-system index
  pure function hardening_irr_planes(slip_normal_c) result(plane_c)

    real(rk), intent(in) :: slip_normal_c(:, :)
    real(rk) :: plane_c(6, size(slip_normal_c, 2))

    integer :: a

    ! Outer product of each normal with itself.
    do a = 1, size(slip_normal_c, 2)
      plane_c(:, a) = mandel_from_matrix(spread(slip_normal_c(:, a), 2, 3)* &
        & spread(slip_normal_c(:, a), 1, 3))
    end do

  end function hardening_irr_planes

  !> Initial damage tensor H_0 = 3 d_l rho_l (1/N_dirs) sum_i (I - n_i (x) n_i)
  !> (addendum_irradiation.md (I4), §5: Barton's (1/V) sum over loops with the loop number
  !> density; V cancels).
  !>
  !> Arguments:
  !>   irr            (in)  hardening_irr_t     the phase's irradiation parameters (loop_avg_c)
  !>   loop_density   (in)  real(rk)            rho_l of the point's element [1/mm^3]
  !>   loop_diameter  (in)  real(rk)            d_l of the point's element [mm]
  !>   dmg_c          (out, result) real(rk)(6) H_0, Mandel, crystal frame [1/mm^2]
  pure function hardening_irr_dmg_0(irr, loop_density, loop_diameter) result(dmg_c)

    type(hardening_irr_t), intent(in) :: irr
    real(rk), intent(in) :: loop_density
    real(rk), intent(in) :: loop_diameter
    real(rk) :: dmg_c(6)

    ! Loop perimeter 3 d per loop, times the number density, times the average loop tensor.
    dmg_c = loop_density*3.0_rk*loop_diameter*irr%loop_avg_c

  end function hardening_irr_dmg_0

  !> Strengths with the network and defect terms (addendum_irradiation.md (I1)):
  !> crss^a = tau_0^a + b mu_0 (sqrt(h_n km_rho) + sqrt(h_d <N^a : H>_+)).
  !>
  !> Theory : Barton et al. (2013) Eq. 3; <x>_+ = max(x, 0) (errata T10, addendum §3).
  !>
  !> Arguments:
  !>   hard    (in)  hardening_params_t  Kocks-Mecking parameters (km_*)
  !>   irr     (in)  hardening_irr_t     irradiation parameters (h_d, plane_c)
  !>   tau_0   (in)  real(rk)(N)         lattice strengths [MPa]
  !>   km_rho  (in)  real(rk)            network dislocation density [1/mm^2]
  !>   dmg_c   (in)  real(rk)(6)         damage tensor H, Mandel, crystal frame [1/mm^2]
  !>   crss    (out, result) real(rk)(N) strengths [MPa]
  !> Locals:
  !>   a        integer   slip system
  !>   network  real(rk)  b mu_0 sqrt(h_n km_rho), the Kocks-Mecking term [MPa]
  !>   contact  real(rk)  N^a : H, the system's contact density [1/mm^2]
  pure function hardening_irr_strengths(hard, irr, tau_0, km_rho, dmg_c) result(crss)

    type(hardening_params_t), intent(in) :: hard
    type(hardening_irr_t), intent(in) :: irr
    real(rk), intent(in) :: tau_0(:)
    real(rk), intent(in) :: km_rho
    real(rk), intent(in) :: dmg_c(6)
    real(rk) :: crss(size(tau_0))

    integer :: a
    real(rk) :: network, contact

    ! Kocks-Mecking strengths (hardening_km_strengths' term) plus each system's defect term from
    ! its contact density (one system at a time: no array temporary).
    network = hard%km_b*hard%km_mu_0*sqrt(hard%km_h_n*km_rho)
    do a = 1, size(tau_0)
      contact = dot_product(dmg_c, irr%plane_c(:, a))
      crss(a) = (tau_0(a) + network) + hard%km_b*hard%km_mu_0*sqrt(irr%h_d*max(contact, 0.0_rk))
    end do

  end function hardening_irr_strengths

  !> Rate of the damage tensor (addendum_irradiation.md (I3); Barton et al. 2013 Eq. 5),
  !> H dot = -eta sum_a (N^a : H) N^a |gammadot^a| = J H with the constant (at fixed slip
  !> rates) symmetric Jacobian J = -eta sum_a |gammadot^a| N^a N^a^T (Mandel).
  !>
  !> Arguments:
  !>   irr       (in)  hardening_irr_t     irradiation parameters (eta, plane_c)
  !>   sliprate  (in)  real(rk)(N)         slip rates [1/s]
  !>   dmg_c     (in)  real(rk)(6)         H, Mandel, crystal frame [1/mm^2]
  !>   rate_c    (out) real(rk)(6)         H dot, Mandel, crystal frame [1/(mm^2 s)]
  !>   jac       (out, optional) real(rk)(6, 6)  d(H dot)/dH [1/s]
  !>
  !> Locals:
  !>   a     integer         slip-system index
  !>   mat   real(rk)(6, 6)  the Jacobian
  pure subroutine hardening_irr_rate(irr, sliprate, dmg_c, rate_c, jac)

    type(hardening_irr_t), intent(in) :: irr
    real(rk), intent(in) :: sliprate(:)
    real(rk), intent(in) :: dmg_c(6)
    real(rk), intent(out) :: rate_c(6)
    real(rk), intent(out), optional :: jac(6, 6)

    integer :: a
    real(rk) :: mat(6, 6)

    ! Sum of the rank-one terms of the active systems (the outer product first, so that the
    ! matrix is exactly symmetric).
    mat = 0.0_rk
    do a = 1, size(sliprate)
      if (.not. abs(sliprate(a)) > 0.0_rk) cycle
      mat = mat - irr%eta*abs(sliprate(a))*(spread(irr%plane_c(:, a), 2, 6)* &
        & spread(irr%plane_c(:, a), 1, 6))
    end do

    ! Rate linear in H.
    rate_c = matmul(mat, dmg_c)
    if (present(jac)) jac = mat

  end subroutine hardening_irr_rate

  !> Forward-Euler estimate of the damage tensor, H = H_prev + dt H dot(H_prev), replaced by
  !> H_prev where dt eta Gammadot >= 1 (addendum_irradiation.md §4: forward Euler overshoots
  !> the contact density of an active system past zero for eta Delta Gamma > 1, which would
  !> give the state loop's first stress solve strengths far from the converged ones; a start
  !> value only, the backward-Euler result being exact in one step).
  !>
  !> Arguments:
  !>   irr       (in)  hardening_irr_t     irradiation parameters (eta, plane_c)
  !>   dmg_prev  (in)  real(rk)(6)         H at the start of the increment, Mandel, crystal
  !>                                       frame [1/mm^2]
  !>   sliprate  (in)  real(rk)(N)         slip rates [1/s]
  !>   dt        (in)  real(rk)            time increment [s]
  !>   dmg_c     (out) real(rk)(6)         the estimate [1/mm^2]
  !>
  !> Locals:
  !>   rate_c  real(rk)(6)  H dot at H_prev [1/(mm^2 s)]
  pure subroutine hardening_irr_forward(irr, dmg_prev, sliprate, dt, dmg_c)

    type(hardening_irr_t), intent(in) :: irr
    real(rk), intent(in) :: dmg_prev(6)
    real(rk), intent(in) :: sliprate(:)
    real(rk), intent(in) :: dt
    real(rk), intent(out) :: dmg_c(6)

    real(rk) :: rate_c(6)

    ! The guard, then the general forward Euler.
    if (dt*irr%eta*sum(abs(sliprate)) >= 1.0_rk) then
      dmg_c = dmg_prev
      return
    end if
    call hardening_irr_rate(irr, sliprate, dmg_prev, rate_c)
    call forward_euler(dmg_prev, rate_c, dt, dmg_c)

  end subroutine hardening_irr_forward

  !> Backward-Euler update of the damage tensor at fixed slip rates,
  !> r(H) = H - H_prev - dt H dot(H) = 0, by time_integration_mod (unchanged; full symmetric
  !> Jacobian, so the normal-equations solve and the whole-step line search; the first Newton
  !> step is exact). Tolerance |r_i| < tol_hardening max_j |H_prev,j| (addendum_irradiation.md
  !> §4); H_prev = 0 (no loops) gives H = 0 without a solve.
  !>
  !> Arguments:
  !>   hard       (in)    hardening_params_t  tolerance and iteration limit
  !>   irr        (in)    hardening_irr_t     irradiation parameters
  !>   dmg_prev   (in)    real(rk)(6)  H at the start of the increment, Mandel, crystal frame
  !>                                   [1/mm^2]
  !>   sliprate   (in)    real(rk)(N)  slip rates of the current stress [1/s]
  !>   dt         (in)    real(rk)     time increment [s]
  !>   dmg_c      (inout) real(rk)(6)  in: the start estimate; out: H at the end [1/mm^2]
  !>   converged  (out)   logical      the solve converged (otherwise dmg_c = dmg_prev)
  !>   work       (inout) type(hardening_work_t)  work arrays (hardening_work_alloc with
  !>                                   irradiation): the solver state be_damage
  !>   num_iters  (out, optional) integer  Newton iterations
  !>
  !> Locals:
  !>   scale    real(rk)          max_j |H_prev,j| [1/mm^2]
  !>   rate_c   real(rk)(6)       H dot at the solver's point [1/(mm^2 s)]
  !>   jac      real(rk)(6, 6)    its Jacobian [1/s]
  pure subroutine hardening_irr_backward(hard, irr, dmg_prev, sliprate, dt, dmg_c, converged, &
    & work, num_iters)

    type(hardening_params_t), intent(in) :: hard
    type(hardening_irr_t), intent(in) :: irr
    real(rk), intent(in) :: dmg_prev(6)
    real(rk), intent(in) :: sliprate(:)
    real(rk), intent(in) :: dt
    real(rk), intent(inout) :: dmg_c(6)
    logical, intent(out) :: converged
    type(hardening_work_t), intent(inout) :: work
    integer, intent(out), optional :: num_iters

    real(rk) :: scale, rate_c(6), jac(6, 6)

    ! No loops: nothing evolves (a zero tolerance could not be met).
    scale = maxval(abs(dmg_prev))
    if (.not. scale > 0.0_rk) then
      dmg_c = 0.0_rk
      converged = .true.
      if (present(num_iters)) num_iters = 0
      return
    end if

    ! The solve from the start estimate, with the rate evaluated where the solver asks.
    call backward_euler_start(work%be_damage, dmg_prev, dmg_c, dt, &
      & spread(hard%tol_hardening*scale, 1, 6), hard%max_iters)
    do while (.not. work%be_damage%done)
      call hardening_irr_rate(irr, sliprate, work%be_damage%x_eval, rate_c, jac)
      call backward_euler_step(work%be_damage, rate_c, jac)
    end do
    call backward_euler_finish(work%be_damage, dmg_c, converged, num_iters)

  end subroutine hardening_irr_backward

end module hardening_mod
