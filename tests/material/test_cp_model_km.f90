! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Tests of the crystal plasticity point solve with Kocks-Mecking hardening (step 7.2, D60,
!> dev/theory/addendum_kocks_mecking.md §3, §5, §9 item 1).
!>
!> Source : FEPX (ROADMAP step 7.2).
!>
!> Checks (fcc, copper-like stiffness, mu_0 = Hill's shear modulus of that stiffness; Barton's
!> Kocks-Mecking parameters in FEPX units otherwise):
!>   - state initialization: km_rho = km_rho_0 and crss = tau_0 + b mu_0 sqrt(h_n rho_0) at every
!>     point; with per-element densities (cp_model_set_element_values) each element its own;
!>   - fcc [001] under uniaxial isochoric {d} = 1e-3 [-1/2, -1/2, 1, 0, 0, 0] /s, 40 increments of
!>     0.25 s, tight tolerances: every increment converges; sqrt(km_rho) follows the backward-
!>     Euler recursion s = (s_prev + dt Gammadot k_1/2)/(1 + dt Gammadot k_2/2) with the returned
!>     total slip rate (relative 1e-9); the strengths are tau_0 + b mu_0 sqrt(h_n km_rho) (1e-13);
!>     the density does not decrease (the first, elastic increments change it by less than
!>     round-off) and has grown by more than half at the end;
!>   - without the km_rho argument the point solve starts from km_rho_0 (same result as passing
!>     it);
!>   - an annealed start (rho_0 = 1e4 /mm^2) with 2 % slip increments (dt = 8 s), where Newton
!>     on km_rho fails (addendum §3.2): every increment converges and follows the recursion.
!>
!> Locals:
!>   work   type(cp_work_t)  point-solve work arrays (cp_work_alloc before each call; D100)
!>   crys       type(crystal_t)          the fcc crystal
!>   stiff      real(rk)(6, 6)           stiffness [MPa]
!>   model      type(cp_model_t)         the model
!>   hard       type(hardening_params_t) Kocks-Mecking parameters
!>   state      type(material_state_t)   a state of 2 elements, 1 point each
!>   kk         real(rk)                 b mu_0 sqrt(h_n) [MPa mm]
!>   g0         real(rk)                 lattice strength [MPa]
program test_cp_model_km

  use kinds_mod, only: rk, ik
  use crystal_mod, only: crystal_t, crystal_init, crystal_fcc, stiffness_cubic, &
    & shear_modulus_hill
  use crystal_plasticity_model_mod, only: cp_model_t, cp_model_init, cp_point_solve, &
    & cp_model_set_element_values, cp_converged
  use crystal_plasticity_model_mod, only: cp_work_t, cp_work_alloc
  use hardening_mod, only: hardening_params_t, hardening_kocks_mecking, hardening_km_k2
  use material_state_mod, only: material_state_t, material_state_alloc
  use assert_mod, only: assert_true, assert_close, test_summary

  implicit none

  type(cp_work_t) :: work
  type(crystal_t) :: crys
  real(rk) :: stiff(6, 6), kk
  real(rk), parameter :: g0 = 100.0_rk
  type(cp_model_t) :: model
  type(hardening_params_t) :: hard
  type(material_state_t) :: state

  call crystal_init(crys, crystal_fcc)
  stiff = stiffness_cubic(168400.0_rk, 121400.0_rk, 75400.0_rk)

  ! Kocks-Mecking parameters (Barton's, FEPX units), mu_0 by Hill.
  hard%law = hardening_kocks_mecking
  hard%km_b = 2.48e-7_rk
  hard%km_mu_0 = shear_modulus_hill(stiff)
  hard%km_h_n = 0.125_rk
  hard%km_rho_0 = 2.0e7_rk
  hard%km_k1 = 4.5e5_rk
  hard%km_k20 = 14.0_rk
  hard%km_gammadot_k0 = 1.0e10_rk
  hard%km_n = 100.0_rk
  hard%tol_hardening = 1.0e-13_rk
  kk = hard%km_b*hard%km_mu_0*sqrt(hard%km_h_n)
  call make_model()

  call check_init_state()
  call check_increments(2.0e7_rk, 0.25_rk, 40, 'Barton start')
  call check_increments(1.0e4_rk, 8.0_rk, 15, 'annealed, 2 % slip increments')
  call check_default_density()

  call test_summary()

contains

  !> The model: fcc, m 0.02, gammadot_0 1 /s, tau_0 = g0, tight tolerances, hardening hard.
  subroutine make_model()

    integer :: n

    n = crys%num_slip
    call cp_model_init(model, stiff, crys%schmid_sym_c, crys%schmid_skew_c, &
      & spread(0.02_rk, 1, n), spread(1.0_rk, 1, n), spread(g0, 1, n), tol_newton=1.0e-14_rk, &
      & max_iters_state=200, tol_state=1.0e-14_rk, hard=hard)

  end subroutine make_model

  !> State initialization with the phase's and with per-element densities.
  !>
  !> Locals:
  !>   i  integer  element
  subroutine check_init_state()

    integer :: i

    call material_state_alloc(state, 1_ik, 2_ik, int(crys%num_slip, ik))
    call model%init_state([1_ik, 2_ik], state)
    call assert_close(state%km_rho(1, :), [2.0e7_rk, 2.0e7_rk], 'init: km_rho = km_rho_0', &
      & rel_tol=1.0e-15_rk)
    do i = 1, 2
      call assert_close(state%crss(:, 1, i), spread(g0 + kk*sqrt(2.0e7_rk), 1, crys%num_slip), &
        & 'init: strengths tau_0 + b mu_0 sqrt(h_n rho_0)', rel_tol=1.0e-14_rk)
    end do

    ! Per-element densities.
    call cp_model_set_element_values(model, rho_0=[1.0e6_rk, 4.0e8_rk])
    call model%init_state([1_ik, 2_ik], state)
    call assert_close(state%km_rho(1, :), [1.0e6_rk, 4.0e8_rk], 'init: per-element rho_0', &
      & rel_tol=1.0e-15_rk)
    call assert_close(state%crss(1, 1, 2), g0 + kk*sqrt(4.0e8_rk), &
      & 'init: per-element strengths', rel_tol=1.0e-14_rk)
    call make_model()

  end subroutine check_init_state

  !> Uniaxial increments on fcc [001] from a zero stress: convergence, the recursion in s, the
  !> strengths from the density, the density increasing.
  !>
  !> Arguments:
  !>   rho_start  (in) real(rk)   initial density [1/mm^2]
  !>   dt         (in) real(rk)   time increment [s]
  !>   num_incr   (in) integer    increments
  !>   label      (in) character  case name
  !> Locals:
  !>   rot, sig, e, d, w  orientation, stress, elastic strain, {d}, w
  !>   rot_n, sig_n, e_n, wp, sec, h, tan  point-solve output
  !>   crss_prev, crss, rate  strengths [MPa], slip rates [1/s]
  !>   rho, rho_prev  density [1/mm^2]
  !>   s_ref, shrate, k2  recursion value [1/mm], Gammadot [1/s], k_2 [-]
  !>   all_converged  logical
  !>   status, k, j  integer
  subroutine check_increments(rho_start, dt, num_incr, label)

    real(rk), intent(in) :: rho_start
    real(rk), intent(in) :: dt
    integer, intent(in) :: num_incr
    character(len=*), intent(in) :: label

    real(rk) :: rot(3, 3), sig(6), e(6), d(6), w(3), rot_n(3, 3), sig_n(6), e_n(6), wp(3)
    real(rk) :: sec(6, 6), h(6), tan(6, 6), rho, rho_prev, s_ref, shrate, k2
    real(rk) :: crss_prev(crys%num_slip), crss(crys%num_slip), rate(crys%num_slip)
    logical :: all_converged
    integer :: status, k, j

    rot = 0.0_rk
    do j = 1, 3
      rot(j, j) = 1.0_rk
    end do
    sig = 0.0_rk
    e = 0.0_rk
    d = 1.0e-3_rk*[-0.5_rk, -0.5_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
    w = 0.0_rk
    rho = rho_start
    crss_prev = g0 + kk*sqrt(rho)
    s_ref = sqrt(rho)
    all_converged = .true.
    do k = 1, num_incr
      rho_prev = rho
      call cp_work_alloc(work, model)
      call cp_point_solve(model, crss_prev, rot, sig, e, d, w, dt, sig_n, e_n, wp, rot_n, &
        & crss, sec, h, tan, rate, status, km_rho=rho, work=work)
      all_converged = all_converged .and. status == cp_converged

      ! The backward-Euler recursion in s with the returned total slip rate, and the strengths.
      shrate = sum(abs(rate))
      k2 = hardening_km_k2(hard, shrate)
      s_ref = (sqrt(rho_prev) + dt*shrate*hard%km_k1/2.0_rk)/(1.0_rk + dt*shrate*k2/2.0_rk)
      call assert_close(sqrt(rho), s_ref, label // ': recursion in sqrt(km_rho)', &
        & rel_tol=1.0e-9_rk)
      call assert_close(crss, spread(g0 + kk*sqrt(rho), 1, crys%num_slip), &
        & label // ': strengths from km_rho', rel_tol=1.0e-13_rk)
      call assert_true(.not. rho < rho_prev, label // ': density does not decrease')

      ! Next increment.
      crss_prev = crss
      rot = rot_n
      sig = sig_n
      e = e_n
    end do
    call assert_true(all_converged, label // ': every increment converged')
    call assert_true(rho > 1.5_rk*rho_start, label // ': density grown by plastic slip')

  end subroutine check_increments

  !> Without the km_rho argument the point solve starts from the phase's km_rho_0.
  !>
  !> Locals: as check_increments; rho, crss2 the results with the argument.
  subroutine check_default_density()

    real(rk) :: rot(3, 3), sig(6), e(6), d(6), w(3), rot_n(3, 3), sig_n(6), e_n(6), wp(3)
    real(rk) :: sec(6, 6), h(6), tan(6, 6), rho
    real(rk) :: crss_prev(crys%num_slip), crss(crys%num_slip), crss2(crys%num_slip)
    real(rk) :: rate(crys%num_slip)
    integer :: status, j

    rot = 0.0_rk
    do j = 1, 3
      rot(j, j) = 1.0_rk
    end do
    sig = 0.0_rk
    e = 0.0_rk
    d = 1.0e-3_rk*[-0.5_rk, -0.5_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
    w = 0.0_rk
    crss_prev = g0 + kk*sqrt(hard%km_rho_0)
    call cp_work_alloc(work, model)
    call cp_point_solve(model, crss_prev, rot, sig, e, d, w, 0.5_rk, sig_n, e_n, wp, rot_n, &
      & crss, sec, h, tan, rate, status, work=work)
    rho = hard%km_rho_0
    call cp_work_alloc(work, model)
    call cp_point_solve(model, crss_prev, rot, sig, e, d, w, 0.5_rk, sig_n, e_n, wp, rot_n, &
      & crss2, sec, h, tan, rate, status, km_rho=rho, work=work)
    call assert_close(crss, crss2, 'default start density km_rho_0', rel_tol=1.0e-15_rk)

  end subroutine check_default_density

end program test_cp_model_km
