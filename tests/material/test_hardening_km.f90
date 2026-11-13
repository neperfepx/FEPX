! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of the Kocks-Mecking law of hardening_mod (step 7.2, D60,
!> dev/theory/addendum_kocks_mecking.md §9 items 1-3 and 5).
!>
!> Source : FEPX (ROADMAP step 7.2; step 9.15: the work arrays of hardening_work_alloc, D100).
!>
!> Checks (Barton et al. 2013's parameters in FEPX units, addendum §2; 12 systems with random
!> lattice strengths, fixed seed):
!>   - strengths (K1): crss^a = tau_0^a + b mu_0 sqrt(h_n km_rho), against the explicit formula;
!>     with h_n = 0 the strengths are tau_0;
!>   - the density rate (K2) (km_rhorate) against the explicit formula, zero at Gammadot = 0;
!>   - k_2 (K3): k_20 (gammadot_k0/Gammadot)^(1/n_km) against the explicit formula; zero at
!>     Gammadot = 0;
!>   - forward Euler in s = sqrt(km_rho): (s_prev + dt (Gammadot/2)(k_1 - k_2 s_prev))^2; a
!>     non-positive estimate (slip increment above 2/k_2) is replaced by rho_prev;
!>   - backward Euler in s against its closed form s = (s_prev + dt Gammadot k_1/2)/
!>     (1 + dt Gammadot k_2/2) (relative 1e-13), in one Newton iteration, from an annealed density,
!>     Barton's and one above saturation, slip increments 1e-3 to 0.1; at Gammadot = 0 the density
!>     is unchanged;
!>   - addendum §9 item 2: at constant slip rates, n increments of backward Euler converge at
!>     first order to the analytic solution s(t) = s_s + (s_0 - s_s) exp(-k_2 Gammadot t/2),
!>     s_s = k_1/k_2 (observed order 0.95-1.05 between successive halvings; finest level,
!>     slip increment 1.25e-3: relative error 2.0e-3 measured, bound 2.5e-3);
!>   - addendum §9 item 3: at constant slip rates, the Kocks-Mecking strengths equal those of
!>     saturation_isotropic with n = 1, h_0 = (K/2)(k_1 - k_2 s_0), tau_s = tau_0^1 + K k_1/k_2
!>     and the initial strengths tau_0 + K s_0 (tau_0 the Kocks-Mecking lattice strength),
!>     K = b mu_0 sqrt(h_n), increment by increment (relative 1e-12).
!>
!> Locals:
!>   hard, hard_iso     type(hardening_params_t)  Kocks-Mecking and the equivalent saturation law
!>   tau_0              real(rk)(N)   lattice strengths [MPa]
!>   crss, crss_iso, crss_prev  real(rk)(N)  strengths [MPa]
!>   rates_s            real(rk)(N)   slip rates [1/s]
!>   rho, rho_prev, rho_start  real(rk)  densities [1/mm^2]
!>   rho_list, incr_list  real(rk)(:)  start densities [1/mm^2] and slip increments [-]
!>   shrate, k2, s_prev, s_exact, kk, dt, t_end  scalars
!>   errors             real(rk)(5)   errors of the convergence study [1/mm]
!>   conv               logical       convergence flag
!>   iters              integer       Newton iterations
!>   i, j, k, num_incr  integer       indices, increments
program test_hardening_km

  use kinds_mod, only: rk
  use hardening_mod, only: hardening_params_t, hardening_kocks_mecking, hardening_isotropic, &
    & hardening_km_strengths, hardening_km_k2, hardening_km_rate, hardening_km_forward, &
    & hardening_km_backward, &
    & hardening_backward, hardening_work_t, hardening_work_alloc
  use tensor_test_utils_mod, only: seed_random
  use assert_mod, only: assert_true, assert_close, test_summary

  implicit none

  integer, parameter :: num = 12
  type(hardening_params_t) :: hard, hard_iso
  type(hardening_work_t) :: work, work_iso
  real(rk) :: tau_0(num), crss(num), crss_iso(num), crss_prev(num), rates_s(num)
  real(rk) :: rho, rho_prev, rho_start, shrate, k2, s_prev, s_exact, kk, dt, t_end
  real(rk) :: errors(5)
  real(rk), parameter :: rho_list(3) = [1.0e4_rk, 2.0e7_rk, 1.0e9_rk]
  real(rk), parameter :: incr_list(4) = [1.0e-3_rk, 1.0e-2_rk, 2.0e-2_rk, 0.1_rk]
  logical :: conv
  integer :: iters, i, j, k, num_incr

  call seed_random(20261012)

  ! Barton's parameters in FEPX units (addendum §2), mu_0 Hill's for Fe.
  hard%law = hardening_kocks_mecking
  hard%km_b = 2.48e-7_rk
  hard%km_mu_0 = 81800.0_rk
  hard%km_h_n = 0.125_rk
  hard%km_rho_0 = 2.0e7_rk
  hard%km_k1 = 4.5e5_rk
  hard%km_k20 = 14.0_rk
  hard%km_gammadot_k0 = 1.0e10_rk
  hard%km_n = 100.0_rk
  hard%tol_hardening = 1.0e-4_rk
  hard%max_iters = 10
  call hardening_work_alloc(work, hard, num)
  call random_number(tau_0)
  tau_0 = 80.0_rk + 20.0_rk*tau_0
  call random_number(rates_s)
  rates_s = rates_s - 0.5_rk
  shrate = sum(abs(rates_s))
  kk = hard%km_b*hard%km_mu_0*sqrt(hard%km_h_n)

  ! Strengths (K1); h_n = 0 gives tau_0.
  crss = hardening_km_strengths(hard, tau_0, 2.0e7_rk)
  call assert_close(crss, tau_0 + 2.48e-7_rk*81800.0_rk*sqrt(0.125_rk*2.0e7_rk), &
    & 'strengths (K1)', rel_tol=1.0e-14_rk)
  call assert_close(crss(1) - tau_0(1), 32.07_rk, 'Barton start: b mu_0 sqrt(h_n rho_0)', &
    & abs_tol=0.01_rk)
  hard_iso = hard
  hard_iso%km_h_n = 0.0_rk
  call assert_close(hardening_km_strengths(hard_iso, tau_0, 2.0e7_rk), tau_0, &
    & 'h_n = 0: strengths tau_0', rel_tol=1.0e-15_rk)

  ! k_2 (K3), zero at Gammadot = 0.
  k2 = hardening_km_k2(hard, shrate)
  call assert_close(k2, 14.0_rk*(1.0e10_rk/shrate)**(1.0_rk/100.0_rk), 'k_2 (K3)', &
    & rel_tol=1.0e-14_rk)
  call assert_close(hardening_km_k2(hard, 1.0_rk), 17.625_rk, 'k_2 at 1/s (addendum §2)', &
    & abs_tol=1.0e-3_rk)
  call assert_true(.not. abs(hardening_km_k2(hard, 0.0_rk)) > 0.0_rk, 'k_2 = 0 at Gammadot = 0')

  ! Density rate (K2) against the explicit formula; zero at Gammadot = 0 (output km_rhorate).
  call assert_close(hardening_km_rate(hard, shrate, 2.0e7_rk), shrate*(hard%km_k1* &
    & sqrt(2.0e7_rk) - k2*2.0e7_rk), 'density rate (K2)', rel_tol=1.0e-13_rk)
  call assert_true(.not. abs(hardening_km_rate(hard, 0.0_rk, 2.0e7_rk)) > 0.0_rk, &
    & 'density rate 0 at Gammadot = 0')

  ! Forward Euler in s, and the guard of a non-positive estimate.
  rho_prev = 2.0e7_rk
  dt = 0.01_rk
  s_prev = sqrt(rho_prev)
  call hardening_km_forward(hard, rho_prev, rates_s, dt, rho)
  call assert_close(rho, (s_prev + dt*0.5_rk*shrate*(hard%km_k1 - k2*s_prev))**2, &
    & 'forward Euler in s', rel_tol=1.0e-13_rk)
  call hardening_km_forward(hard, 1.0e12_rk, rates_s, 10.0_rk/shrate, rho)
  call assert_close(rho, 1.0e12_rk, 'forward Euler: non-positive estimate replaced by rho_prev', &
    & rel_tol=1.0e-15_rk)

  ! Backward Euler in s against its closed form, in one Newton iteration.
  do i = 1, size(rho_list)
    do j = 1, size(incr_list)
      rho_prev = rho_list(i)
      dt = incr_list(j)/shrate
      s_prev = sqrt(rho_prev)
      call hardening_km_forward(hard, rho_prev, rates_s, dt, rho_start)
      rho = rho_start
      call hardening_km_backward(hard, rho_prev, rates_s, dt, rho, conv, work, iters)
      s_exact = (s_prev + dt*shrate*hard%km_k1/2.0_rk)/(1.0_rk + dt*shrate*k2/2.0_rk)
      call assert_true(conv .and. iters == 1, 'backward Euler: converged in one iteration')
      call assert_close(rho, s_exact**2, 'backward Euler in s: closed form', rel_tol=1.0e-13_rk)
    end do
  end do
  rho = 3.0e6_rk
  call hardening_km_backward(hard, 3.0e6_rk, 0.0_rk*rates_s, 1.0_rk, rho, conv, work)
  call assert_true(conv, 'backward Euler at Gammadot = 0: converged')
  call assert_close(rho, 3.0e6_rk, 'backward Euler at Gammadot = 0: unchanged', &
    & rel_tol=1.0e-15_rk)

  ! Addendum §9 item 2: first-order convergence to the analytic s(t) at constant slip rates
  ! (from an annealed density 1e4 /mm^2 to a slip of 0.2).
  t_end = 0.2_rk/shrate
  do k = 1, size(errors)
    num_incr = 10*2**(k - 1)
    dt = t_end/num_incr
    rho = 1.0e4_rk
    do i = 1, num_incr
      rho_prev = rho
      call hardening_km_forward(hard, rho_prev, rates_s, dt, rho)
      call hardening_km_backward(hard, rho_prev, rates_s, dt, rho, conv, work)
    end do
    s_exact = hard%km_k1/k2 + (sqrt(1.0e4_rk) - hard%km_k1/k2)*exp(-k2*shrate*t_end/2.0_rk)
    errors(k) = abs(sqrt(rho) - s_exact)
  end do
  do k = 2, size(errors)
    call assert_true(log(errors(k - 1)/errors(k))/log(2.0_rk) > 0.95_rk .and. &
      & log(errors(k - 1)/errors(k))/log(2.0_rk) < 1.05_rk, &
      & 'first order towards the analytic solution')
  end do
  call assert_true(errors(size(errors)) < 2.5e-3_rk*s_exact, 'finest error below 2.5e-3 of s')

  ! Addendum §9 item 3: equal to saturation_isotropic with n = 1 and the mapped parameters,
  ! increment by increment (Barton's start, slip increments of 1 %).
  rho = 2.0e7_rk
  hard_iso = hard
  hard_iso%law = hardening_isotropic
  allocate (hard_iso%h_0(num), hard_iso%tau_s(num), hard_iso%n(num))
  hard_iso%n = 1.0_rk
  hard_iso%h_0 = kk/2.0_rk*(hard%km_k1 - k2*sqrt(rho))
  hard_iso%tau_s = tau_0(1) + kk*hard%km_k1/k2
  crss_iso = hardening_km_strengths(hard, tau_0, rho)
  call hardening_work_alloc(work_iso, hard_iso, num)
  dt = 0.01_rk/shrate
  do i = 1, 30
    rho_prev = rho
    crss_prev = crss_iso
    call hardening_km_backward(hard, rho_prev, rates_s, dt, rho, conv, work)
    crss = hardening_km_strengths(hard, tau_0, rho)
    call hardening_backward(hard_iso, crss_prev, hardening_km_strengths(hard, tau_0, 2.0e7_rk), &
      & rates_s, dt, crss_iso, conv, work_iso)
    call assert_close(crss, crss_iso, 'equal to saturation_isotropic, n = 1 (mapped)', &
      & rel_tol=1.0e-12_rk)
  end do

  call test_summary()

end program test_hardening_km
