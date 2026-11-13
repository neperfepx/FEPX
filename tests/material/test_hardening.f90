! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of hardening_mod (step 6.4, D51, dev/theory/addendum_hardening.md): the two
!> saturation laws and the general backward-Euler integration of the strengths.
!>
!> Source : FEPX (ROADMAP step 6.4; step 9.15: the routines with their work arrays, D100, called
!>          through the adapters rates, forward and backward below).
!>
!> Checks (fixed random seed; 12 systems):
!>   - default interaction matrix: 1 on the diagonal, 1.4 off it;
!>   - saturation_isotropic (§2): every system gets h_0 ((tau_s - crss^1)/(tau_s - tau_0^1))^n
!>     Gammadot (system 1 sets it), against the explicit formula; the Jacobian has only its first
!>     column, equal to the central finite difference in crss^1 (n = 1, 1.5, 2; relative 1e-7);
!>   - saturation_anisotropic (§3): per-system parameters and a random matrix H >= 0: the rate of
!>     every system against the explicit formula with sum_b H_ab |gammadot^b|; diagonal Jacobian
!>     equal to the finite differences;
!>   - saturation: no rate at or above tau_s (isotropic: system 1; anisotropic: per system);
!>   - forward Euler: crss_prev + dt crss_dot(crss_prev);
!>   - backward Euler, n = 1, closed forms: isotropic, every system crss_prev^a + Delta with Delta =
!>     c (tau_s - crss_prev^1)/(1 + c), c = dt h_0 Gammadot/(tau_s - tau_0^1); anisotropic, every
!>     system (crss_prev^a + c^a tau_s^a)/(1 + c^a), c^a = dt h_0^a drive^a/(tau_s^a - tau_0^a)
!>     (1e-12);
!>   - addendum §5 test 1: the anisotropic update (n = 2, per-system parameters) equals the root
!>     of each system's scalar equation found by an independent scalar Newton in this test (the
!>     update is decoupled across systems; relative 1e-10);
!>   - addendum §5 test 4: uniform parameters, one initial strength and H = all ones: the
!>     anisotropic update equals the isotropic one (1e-12);
!>   - law none: unchanged; failure (one iteration on a strongly nonlinear update): not
!>     converged and the strengths reset to crss_prev (v2).
!>   Step 7.1 (D59, addendum_hardening_modifiers.md):
!>   - per-point values (hardening_point_t): the element's tau_s and h_0 replace the phase's;
!>   - saturation evolution (§4): tau_s = tau_s0 (Gammadot/gammadot_s0)^m' in the rates of both
!>     laws (one value for all systems), against the explicit formula; Jacobian by finite
!>     differences;
!>   - cyclic hardening, isotropic law (§3.1, v2): only systems with accumulated slip A^b >= a
!>     (crss^1/tau_s)^c drive hardening, f clamped to [1e-6, 10] d_eq (both clamps reached), against
!>     the explicit formula; Jacobian column 1 by finite differences;
!>   - cyclic hardening, anisotropic law (§3.3): each system judged by its own strength,
!>     f^a = clamp(sum_b H_ab active_b |gammadot^b|), against the explicit formula; diagonal
!>     Jacobian by finite differences; with H = all ones, uniform parameters and equal strengths
!>     the anisotropic rates equal the isotropic ones.
!>
!> Locals:
!>   hard, hard_iso   type(hardening_params_t)  parameters under test
!>   tau_0, prev, crss, crss2  real(rk)(N)  initial, start and current strengths [MPa]
!>   rate             real(rk)(N)     rates [MPa/s]
!>   rp               real(rk)(N)     a finite-difference Jacobian column [1/s]
!>   jac              real(rk)(N, N)  Jacobian [1/s]
!>   rates_s          real(rk)(N)     slip rates [1/s]
!>   drive            real(rk)(N)     sum_b H_ab |gammadot^b| [1/s]
!>   dt, c, ratio, hstep, expo  scalars
!>   pt               type(hardening_point_t)  per-point values
!>   tau_s, deq, f, crit  real(rk)       saturation strength [MPa], d_eq [1/s], cyclic drive [1/s],
!>                                    threshold [-]
!>   active           real(rk)(N)     1 for systems past the cyclic threshold
!>   expected, rate2  real(rk)(N)     reference and second rates [MPa/s]
!>   conv             logical         convergence flag
!>   i, k, a          integer         indices
program test_hardening

  use kinds_mod, only: rk
  use hardening_mod, only: hardening_params_t, hardening_none, hardening_isotropic, &
    & hardening_anisotropic, hardening_default_matrix, hardening_rates, hardening_forward, &
    & hardening_backward, hardening_point_t, hardening_point_alloc, hardening_point_set, &
    & hardening_work_t, hardening_work_alloc
  use tensor_test_utils_mod, only: seed_random
  use assert_mod, only: assert_true, assert_close, test_summary

  implicit none

  integer, parameter :: num = 12
  type(hardening_params_t) :: hard, hard_iso
  real(rk) :: tau_0(num), prev(num), crss(num), crss2(num), rate(num), rp(num)
  real(rk) :: jac(num, num), rates_s(num), drive(num)
  real(rk) :: dt, c, ratio, hstep, expo(3), shrate
  logical :: conv
  integer :: i, k, a
  type(hardening_point_t) :: pt
  real(rk) :: tau_s, deq, f, crit, active(num), expected(num), rate2(num)

  call seed_random(20261011)
  call hardening_point_alloc(pt, num)
  expo = [1.0_rk, 1.5_rk, 2.0_rk]
  hstep = 1.0e-4_rk

  ! Default interaction matrix.
  jac = hardening_default_matrix(num)
  call assert_true(all([((abs(jac(a, k) - merge(1.0_rk, 1.4_rk, a == k)) < 1.0e-15_rk, &
    & a = 1, num), k = 1, num)]), 'default_matrix')

  do i = 1, 20
    ! Random state below saturation and random slip rates.
    call random_number(tau_0)
    tau_0 = 150.0_rk + 100.0_rk*tau_0
    call random_number(rates_s)
    rates_s = 1.0e-3_rk*(rates_s - 0.5_rk)
    shrate = sum(abs(rates_s))

    ! --- saturation_isotropic: one rate, set by system 1 ---
    call set_params(hard, hardening_isotropic, 200.0_rk, 330.0_rk, 1.0_rk, uniform=.true.)
    call random_number(prev)
    prev = tau_0 + (330.0_rk - tau_0)*0.8_rk*prev
    prev(1) = min(prev(1), 320.0_rk)
    do k = 1, 3
      hard%n = expo(k)
      call rates(hard, prev, tau_0, rates_s, rate, jac)
      ratio = (330.0_rk - prev(1))/(330.0_rk - tau_0(1))
      call assert_close(rate, spread(200.0_rk*ratio**expo(k)*shrate, 1, num), 'iso_rate', &
        & rel_tol=1.0e-14_rk)
      call fd_column(hard, prev, tau_0, rates_s, 1, rp)
      call assert_close(jac(:, 1), rp, 'iso_jacobian_fd', rel_tol=1.0e-7_rk)
      call assert_true(.not. any(abs(jac(:, 2:)) > 0.0_rk), 'iso_jacobian_first_column_only')
    end do

    ! --- saturation_anisotropic: per-system parameters, random H >= 0 ---
    call set_params(hard, hardening_anisotropic, 200.0_rk, 330.0_rk, 1.0_rk, uniform=.false.)
    call random_number(hard%interaction)
    hard%interaction = 2.0_rk*hard%interaction
    call random_number(prev)
    prev = tau_0 + (hard%tau_s - tau_0)*0.8_rk*prev
    drive = matmul(hard%interaction, abs(rates_s))
    do k = 1, 3
      hard%n = expo(k)
      call rates(hard, prev, tau_0, rates_s, rate, jac)
      call assert_close(rate, hard%h_0*((hard%tau_s - prev)/(hard%tau_s - tau_0))**expo(k)* &
        & drive, 'aniso_rate', rel_tol=1.0e-13_rk)
      do a = 1, num
        call fd_column(hard, prev, tau_0, rates_s, a, rp)
        call assert_close(jac(:, a), rp, 'aniso_jacobian_fd', rel_tol=1.0e-7_rk, &
          & abs_tol=1.0e-14_rk)
      end do
    end do

    ! Forward Euler.
    hard%n = 1.5_rk
    dt = 0.3_rk
    call rates(hard, prev, tau_0, rates_s, rate, jac)
    call forward(hard, prev, tau_0, rates_s, dt, crss)
    call assert_close(crss, prev + dt*rate, 'forward_euler', rel_tol=1.0e-15_rk)

    ! Backward Euler, n = 1, closed forms.
    dt = 50.0_rk
    hard%n = 1.0_rk
    crss = prev
    call backward(hard, prev, tau_0, rates_s, dt, crss, conv)
    call assert_true(conv, 'aniso_n1_converged')
    call assert_close(crss, (prev + dt*hard%h_0*drive/(hard%tau_s - tau_0)*hard%tau_s)/(1.0_rk + &
      & dt*hard%h_0*drive/(hard%tau_s - tau_0)), 'aniso_n1_closed_form', rel_tol=1.0e-12_rk)
    call set_params(hard_iso, hardening_isotropic, 200.0_rk, 330.0_rk, 1.0_rk, uniform=.true.)
    crss = prev
    call backward(hard_iso, prev, tau_0, rates_s, dt, crss, conv)
    c = dt*200.0_rk*shrate/(330.0_rk - tau_0(1))
    call assert_true(conv, 'iso_n1_converged')
    call assert_close(crss, prev + c*(330.0_rk - prev(1))/(1.0_rk + c), 'iso_n1_closed_form', &
      & rel_tol=1.0e-12_rk)

    ! Addendum §5 test 1: n = 2, the update equals the per-system scalar roots.
    hard%n = 2.0_rk
    hard%tol_hardening = 1.0e-12_rk
    hard%max_iters = 50
    crss = prev
    call backward(hard, prev, tau_0, rates_s, dt, crss, conv)
    call assert_true(conv, 'aniso_n2_converged')
    do a = 1, num
      call assert_close(crss(a), scalar_root(hard%h_0(a), hard%tau_s(a), hard%n(a), tau_0(a), &
        & prev(a), drive(a), dt), 'test1_vector_equals_scalar_roots', rel_tol=1.0e-10_rk)
    end do

    ! Addendum §5 test 4: uniform parameters, one tau_0, H = ones: anisotropic = isotropic.
    tau_0 = 200.0_rk
    prev = 260.0_rk
    call set_params(hard, hardening_anisotropic, 200.0_rk, 330.0_rk, 1.0_rk, uniform=.true.)
    hard%interaction = 1.0_rk
    hard%n = 1.5_rk
    call set_params(hard_iso, hardening_isotropic, 200.0_rk, 330.0_rk, 1.0_rk, uniform=.true.)
    hard_iso%n = 1.5_rk
    hard%tol_hardening = 1.0e-13_rk
    hard_iso%tol_hardening = 1.0e-13_rk
    hard%max_iters = 50
    hard_iso%max_iters = 50
    crss = prev
    crss2 = prev
    call backward(hard, prev, tau_0, rates_s, dt, crss, conv)
    call backward(hard_iso, prev, tau_0, rates_s, dt, crss2, conv)
    call assert_close(crss, crss2, 'test4_aniso_ones_equals_iso', rel_tol=1.0e-12_rk)
  end do

  ! Saturation: no rate at or above tau_s.
  call set_params(hard, hardening_isotropic, 200.0_rk, 330.0_rk, 1.0_rk, uniform=.true.)
  tau_0 = 200.0_rk
  prev = 300.0_rk
  prev(1) = 330.0_rk
  call rates(hard, prev, tau_0, rates_s, rate, jac)
  call assert_true(.not. any(abs(rate) > 0.0_rk), 'iso_saturated_system1')
  call set_params(hard, hardening_anisotropic, 200.0_rk, 330.0_rk, 1.0_rk, uniform=.true.)
  hard%interaction = hardening_default_matrix(num)
  call rates(hard, prev, tau_0, rates_s, rate, jac)
  call assert_true(.not. (abs(rate(1)) > 0.0_rk) .and. all(rate(2:) > 0.0_rk), &
    & 'aniso_saturated_per_system')

  ! Law none: unchanged.
  hard%law = hardening_none
  crss = 250.0_rk
  prev = 210.0_rk
  call backward(hard, prev, tau_0, rates_s, 1.0_rk, crss, conv)
  call assert_true(conv .and. all(abs(crss - prev) < 1.0e-300_rk), 'none_unchanged')

  ! Failure: one iteration on a strongly nonlinear update; reset to crss_prev.
  call set_params(hard, hardening_anisotropic, 5000.0_rk, 330.0_rk, 5.0_rk, uniform=.true.)
  hard%interaction = hardening_default_matrix(num)
  hard%max_iters = 1
  rates_s = 1.0_rk
  prev = 210.0_rk
  crss = prev
  call backward(hard, prev, tau_0, rates_s, 10.0_rk, crss, conv)
  call assert_true(.not. conv, 'failure_not_converged')
  call assert_true(all(abs(crss - prev) < 1.0e-300_rk), 'failure_reset_to_prev')

  ! --- Step 7.1: per-point values, saturation evolution, cyclic hardening ---
  call random_number(rates_s)
  rates_s = 1.0e-3_rk*(rates_s - 0.5_rk)
  shrate = sum(abs(rates_s))
  tau_0 = 200.0_rk
  call random_number(prev)
  prev = 210.0_rk + 40.0_rk*prev

  ! Per-point tau_s and h_0 replace the phase's.
  call set_params(hard, hardening_anisotropic, 200.0_rk, 330.0_rk, 1.5_rk, uniform=.true.)
  call hardening_point_set(pt, hard, tau_s=spread(400.0_rk, 1, num), h_0=spread(50.0_rk, 1, &
    & num))
  call rates(hard, prev, tau_0, rates_s, rate, jac, pt)
  call assert_close(rate, 50.0_rk*((400.0_rk - prev)/(400.0_rk - tau_0))**1.5_rk* &
    & matmul(hard%interaction, abs(rates_s)), 'point_values', rel_tol=1.0e-13_rk)

  ! Saturation evolution, both laws.
  do k = 1, 2
    call set_params(hard, merge(hardening_isotropic, hardening_anisotropic, k == 1), 200.0_rk, &
      & 330.0_rk, 1.5_rk, uniform=.true.)
    hard%saturation_evolution = .true.
    hard%gammadot_s0 = 1.0e-3_rk
    hard%m_prime = 0.05_rk
    tau_s = 330.0_rk*(shrate/1.0e-3_rk)**0.05_rk
    call rates(hard, prev, tau_0, rates_s, rate, jac)
    if (k == 1) then
      expected = 200.0_rk*((tau_s - prev(1))/(tau_s - tau_0(1)))**1.5_rk*shrate
    else
      expected = 200.0_rk*((tau_s - prev)/(tau_s - tau_0))**1.5_rk*matmul(hard%interaction, &
        & abs(rates_s))
    end if
    call assert_close(rate, expected, 'saturation_evolution_rate', rel_tol=1.0e-13_rk)
    do a = 1, num
      call fd_column(hard, prev, tau_0, rates_s, a, rp)
      call assert_close(jac(:, a), rp, 'saturation_evolution_jacobian_fd', rel_tol=1.0e-7_rk, &
        & abs_tol=1.0e-14_rk)
    end do
  end do

  ! Cyclic, isotropic law: half the systems past system 1's threshold.
  call set_params(hard, hardening_isotropic, 200.0_rk, 330.0_rk, 1.5_rk, uniform=.true.)
  hard%cyclic = .true.
  hard%cyclic_a = 0.05_rk
  hard%cyclic_c = 3.5_rk
  crit = 0.05_rk*(prev(1)/330.0_rk)**3.5_rk
  active = merge(1.0_rk, 0.0_rk, mod([(a, a = 1, num)], 2) == 0)
  deq = 1.0e-3_rk
  call hardening_point_set(pt, hard, slip_accum=merge(2.0_rk*crit, 0.5_rk*crit, &
    & active > 0.5_rk), defrate_eq=deq)
  f = min(max(sum(active*abs(rates_s)), 1.0e-6_rk*deq), 10.0_rk*deq)
  call rates(hard, prev, tau_0, rates_s, rate, jac, pt)
  call assert_close(rate, spread(200.0_rk*((330.0_rk - prev(1))/(330.0_rk - tau_0(1)))**1.5_rk* &
    & f, 1, num), 'cyclic_iso_rate', rel_tol=1.0e-13_rk)
  call fd_column(hard, prev, tau_0, rates_s, 1, rp, pt)
  call assert_close(jac(:, 1), rp, 'cyclic_iso_jacobian_fd', rel_tol=1.0e-7_rk)
  ! The clamps: no active system (lower bound), a huge slip rate (upper bound).
  pt%slip_accum = 0.0_rk
  call rates(hard, prev, tau_0, rates_s, rate, jac, pt)
  call assert_close(rate(1), 200.0_rk*((330.0_rk - prev(1))/(330.0_rk - tau_0(1)))**1.5_rk* &
    & 1.0e-6_rk*deq, 'cyclic_iso_clamp_low', rel_tol=1.0e-13_rk)
  pt%slip_accum = 1.0_rk
  call rates(hard, prev, tau_0, 1.0e3_rk*rates_s, rate, jac, pt)
  call assert_close(rate(1), 200.0_rk*((330.0_rk - prev(1))/(330.0_rk - tau_0(1)))**1.5_rk* &
    & 10.0_rk*deq, 'cyclic_iso_clamp_high', rel_tol=1.0e-13_rk)

  ! Cyclic, anisotropic law: each system judged by its own strength, clamp per system.
  call set_params(hard, hardening_anisotropic, 200.0_rk, 330.0_rk, 1.5_rk, uniform=.false.)
  hard%cyclic = .true.
  hard%cyclic_a = 0.05_rk
  hard%cyclic_c = 3.5_rk
  call hardening_point_set(pt, hard, slip_accum=merge(2.0_rk, 0.5_rk, active > 0.5_rk)* &
    & 0.05_rk*(prev/hard%tau_s)**3.5_rk, defrate_eq=deq)
  call rates(hard, prev, tau_0, rates_s, rate, jac, pt)
  do a = 1, num
    f = min(max(sum(hard%interaction(a, :)*active*abs(rates_s)), 1.0e-6_rk*deq), 10.0_rk*deq)
    expected(a) = hard%h_0(a)*((hard%tau_s(a) - prev(a))/(hard%tau_s(a) - tau_0(a)))**1.5_rk*f
  end do
  call assert_close(rate, expected, 'cyclic_aniso_rate', rel_tol=1.0e-13_rk)
  do a = 1, num
    call fd_column(hard, prev, tau_0, rates_s, a, rp, pt)
    call assert_close(jac(:, a), rp, 'cyclic_aniso_jacobian_fd', rel_tol=1.0e-7_rk, &
      & abs_tol=1.0e-14_rk)
  end do
  ! H = all ones, uniform parameters, equal strengths: anisotropic = isotropic.
  call set_params(hard, hardening_anisotropic, 200.0_rk, 330.0_rk, 1.5_rk, uniform=.true.)
  call set_params(hard_iso, hardening_isotropic, 200.0_rk, 330.0_rk, 1.5_rk, uniform=.true.)
  hard%interaction = 1.0_rk
  hard%cyclic = .true.
  hard_iso%cyclic = .true.
  hard%cyclic_a = 0.05_rk
  hard_iso%cyclic_a = 0.05_rk
  hard%cyclic_c = 3.5_rk
  hard_iso%cyclic_c = 3.5_rk
  crss = 250.0_rk
  call hardening_point_set(pt, hard, slip_accum=merge(2.0_rk, 0.5_rk, active > 0.5_rk)* &
    & 0.05_rk*(250.0_rk/330.0_rk)**3.5_rk, defrate_eq=deq)
  call rates(hard, crss, tau_0, rates_s, rate, jac, pt)
  call rates(hard_iso, crss, tau_0, rates_s, rate2, jac, pt)
  call assert_close(rate, rate2, 'cyclic_aniso_ones_equals_iso', rel_tol=1.0e-13_rk)

  call test_summary()

contains

  !> hardening_rates with the work arrays of a phase (hardening_work_alloc) and the given or
  !> the phase's per-point values.
  !>
  !> Arguments:
  !>   p, crss, t0, sr  (in)  parameters, strengths, initial strengths, slip rates
  !>   r, j   (out) real(rk)(:), (:, :)  rates and Jacobian
  !>   pp     (in, optional) type(hardening_point_t)  per-point values
  !> Locals:
  !>   work   type(hardening_work_t)  work arrays
  subroutine rates(p, crss, t0, sr, r, j, pp)

    type(hardening_params_t), intent(in) :: p
    real(rk), intent(in) :: crss(:), t0(:), sr(:)
    real(rk), intent(out) :: r(:), j(:, :)
    type(hardening_point_t), intent(in), optional :: pp

    type(hardening_work_t) :: work

    call hardening_work_alloc(work, p, size(crss))
    if (present(pp)) work%point = pp
    call hardening_rates(p, crss, t0, sr, r, j, work%point, work%rates)

  end subroutine rates

  !> hardening_forward with the work arrays of a phase and the phase's per-point values.
  !>
  !> Arguments:
  !>   p, prv, t0, sr, dtt  (in)  parameters, start and initial strengths, slip rates, dt
  !>   c      (out) real(rk)(:)    estimated strengths
  !> Locals:
  !>   work   type(hardening_work_t)  work arrays
  subroutine forward(p, prv, t0, sr, dtt, c)

    type(hardening_params_t), intent(in) :: p
    real(rk), intent(in) :: prv(:), t0(:), sr(:), dtt
    real(rk), intent(out) :: c(:)

    type(hardening_work_t) :: work

    call hardening_work_alloc(work, p, size(prv))
    call hardening_forward(p, prv, t0, sr, dtt, c, work)

  end subroutine forward

  !> hardening_backward with the work arrays of a phase and the phase's per-point values.
  !>
  !> Arguments:
  !>   p, prv, t0, sr, dtt  (in)  parameters, start and initial strengths, slip rates, dt
  !>   c      (inout) real(rk)(:)  start estimate in, strengths out
  !>   ok     (out) logical        converged
  !> Locals:
  !>   work   type(hardening_work_t)  work arrays
  subroutine backward(p, prv, t0, sr, dtt, c, ok)

    type(hardening_params_t), intent(in) :: p
    real(rk), intent(in) :: prv(:), t0(:), sr(:), dtt
    real(rk), intent(inout) :: c(:)
    logical, intent(out) :: ok

    type(hardening_work_t) :: work

    call hardening_work_alloc(work, p, size(prv))
    call hardening_backward(p, prv, t0, sr, dtt, c, ok, work)

  end subroutine backward

  !> Parameters of a law: uniform values, or per-system values spread around them.
  !>
  !> Arguments:
  !>   p        (out) type(hardening_params_t)  parameters
  !>   law      (in)  integer   law code
  !>   h0, tau_s, nn  (in) real(rk)  base values [MPa], [MPa], [-]
  !>   uniform  (in)  logical   same values for every system
  !> Locals:
  !>   noise    real(rk)(num)   random factors
  subroutine set_params(p, law, h0, tau_s, nn, uniform)

    type(hardening_params_t), intent(out) :: p
    integer, intent(in) :: law
    real(rk), intent(in) :: h0, tau_s, nn
    logical, intent(in) :: uniform

    real(rk) :: noise(num)

    p%law = law
    allocate (p%h_0(num), p%tau_s(num), p%n(num), p%interaction(num, num))
    p%h_0 = h0
    p%tau_s = tau_s
    p%n = nn
    p%interaction = hardening_default_matrix(num)
    if (.not. uniform) then
      call random_number(noise)
      p%h_0 = h0*(0.5_rk + noise)
      call random_number(noise)
      p%tau_s = tau_s + 100.0_rk*noise
    end if

  end subroutine set_params

  !> Column a of the Jacobian by central differences in crss^a.
  !>
  !> Arguments:
  !>   p, crss, t0, sr  (in)  parameters, strengths, initial strengths, slip rates
  !>   a      (in)  integer      column
  !>   col    (out) real(rk)(num)  d crss_dot / d crss^a
  !>   pp     (in, optional) type(hardening_point_t)  per-point values
  !> Locals:
  !>   crss_plus, crss_minus  real(rk)(num)  perturbed strengths
  !>   rpl, rmi  real(rk)(num)  their rates
  !>   jj        real(rk)(num, num)  scratch Jacobian
  subroutine fd_column(p, crss, t0, sr, a, col, pp)

    type(hardening_params_t), intent(in) :: p
    real(rk), intent(in) :: crss(num), t0(num), sr(num)
    integer, intent(in) :: a
    real(rk), intent(out) :: col(num)
    type(hardening_point_t), intent(in), optional :: pp

    real(rk) :: crss_plus(num), crss_minus(num), rpl(num), rmi(num), jj(num, num)

    crss_plus = crss
    crss_plus(a) = crss_plus(a) + hstep
    crss_minus = crss
    crss_minus(a) = crss_minus(a) - hstep
    call rates(p, crss_plus, t0, sr, rpl, jj, pp)
    call rates(p, crss_minus, t0, sr, rmi, jj, pp)
    col = (rpl - rmi)/(2.0_rk*hstep)

  end subroutine fd_column

  !> Root of one system's backward-Euler equation
  !> crss - crss_prev - dt h_0 ((tau_s - crss)/(tau_s - tau_0))^n drive = 0 by an independent
  !> scalar Newton (bisection safeguard), for test 1.
  !>
  !> Arguments:
  !>   h0, tau_s, nn, tau_0, crss_prev, dr, dtt  (in)  the system's parameters, start strength,
  !>                                                   drive, dt
  !> Locals:
  !>   x, lo, hi, f, df  real(rk)  iterate, bracket, residual and derivative
  !>   it    integer      iteration
  function scalar_root(h0, tau_s, nn, tau_0, crss_prev, dr, dtt) result(x)

    real(rk), intent(in) :: h0, tau_s, nn, tau_0, crss_prev, dr, dtt
    real(rk) :: x

    real(rk) :: lo, hi, f, df
    integer :: it

    ! The root lies in [crss_prev, tau_s] (the rate is positive below tau_s, zero at it).
    lo = crss_prev
    hi = tau_s
    x = 0.5_rk*(lo + hi)
    do it = 1, 200
      f = x - crss_prev - dtt*h0*((tau_s - x)/(tau_s - tau_0))**nn*dr
      df = 1.0_rk + dtt*h0*nn*((tau_s - x)/(tau_s - tau_0))**(nn - 1.0_rk)*dr/(tau_s - tau_0)
      if (abs(f) < 1.0e-13_rk*tau_s) exit
      ! Shrink the bracket, Newton step, bisection when the step leaves the bracket.
      if (f > 0.0_rk) then
        hi = x
      else
        lo = x
      end if
      x = x - f/df
      if (x <= lo .or. x >= hi) x = 0.5_rk*(lo + hi)
    end do

  end function scalar_root

end program test_hardening
