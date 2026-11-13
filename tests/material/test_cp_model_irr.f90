! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Tests of the crystal plasticity point solve with irradiation softening (step 7.3, D61,
!> dev/theory/addendum_irradiation.md §4, §7, §9, §10 items 1 and 3).
!>
!> Source : FEPX (ROADMAP step 7.3).
!>
!> Checks (bcc {110}<111>, Fe stiffness, mu_0 = Hill's; Barton's Kocks-Mecking and irradiation
!> parameters in FEPX units, loops 1.63e22 m^-3 = 1.63e13 /mm^3, d = 100 b, <111> normals):
!>   - state initialization: dmg_c = H_0 = 3 d_l rho_l (2/3) I, the strengths tau_0 + b mu_0
!>     (sqrt(h_n rho_0) + sqrt(h_d N:H_0)); per-element loop densities and diameters;
!>   - bcc [001] under uniaxial isochoric {d} = 1e-3 /s, 40 increments of 1 s: every increment converges; H
!>     follows the backward-Euler recursion (I + dt eta sum |gammadot^a| N^a N^a^T) H = H_prev
!>     with the returned slip rates (relative 1e-9 of |H_prev|); the strengths are (I1) from
!>     km_rho and H (1e-13); the contact density of every plane does not increase on the first
!>     increments (every active contact density stays positive here) and the smallest has
!>     fallen below 0.7 of its start;
!>   - no loops (rho_l = 0): H stays 0 and the results equal Kocks-Mecking without
!>     irradiation exactly (the same strengths, stress, density), with the plain state loop
!>     (damp_state off; the damping, D62, changes the iterates).
!>
!> Locals:
!>   work   type(cp_work_t)  point-solve work arrays (cp_work_alloc before each call; D100)
!>   crys        type(crystal_t)           the bcc crystal
!>   stiff       real(rk)(6, 6)            stiffness [MPa]
!>   model, model_km  type(cp_model_t)     the model with and without irradiation
!>   hard        type(hardening_params_t)  Kocks-Mecking parameters
!>   irr         type(hardening_irr_t)     irradiation parameters
!>   state       type(material_state_t)    a state of 2 elements, 1 point each
!>   g0          real(rk)                  lattice strength [MPa]
program test_cp_model_irr

  use kinds_mod, only: rk, ik
  use crystal_mod, only: crystal_t, crystal_init, crystal_bcc, stiffness_cubic, &
    & shear_modulus_hill
  use crystal_plasticity_model_mod, only: cp_model_t, cp_model_init, cp_point_solve, &
    & cp_model_set_element_values, cp_converged
  use crystal_plasticity_model_mod, only: cp_work_t, cp_work_alloc
  use hardening_mod, only: hardening_params_t, hardening_irr_t, hardening_kocks_mecking, &
    & hardening_km_strengths, hardening_irr_loop_average, hardening_irr_planes, &
    & hardening_irr_dmg_0, hardening_irr_strengths
  use linalg_mod, only: normal_equations_solve
  use material_state_mod, only: material_state_t, material_state_alloc
  use assert_mod, only: assert_true, assert_close, test_summary

  implicit none

  type(cp_work_t) :: work
  type(crystal_t) :: crys
  real(rk) :: stiff(6, 6)
  real(rk), parameter :: g0 = 90.0_rk
  type(cp_model_t) :: model, model_km
  type(hardening_params_t) :: hard
  type(hardening_irr_t) :: irr
  type(material_state_t) :: state

  call crystal_init(crys, crystal_bcc)
  stiff = stiffness_cubic(231400.0_rk, 134700.0_rk, 116400.0_rk)

  ! Barton's parameters (FEPX units), mu_0 by Hill.
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
  irr%h_d = 0.675_rk
  irr%eta = 100.0_rk
  irr%loop_density = 1.63e13_rk
  irr%loop_diameter = 100.0_rk*hard%km_b
  irr%loop_avg_c = hardening_irr_loop_average(reshape([1, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1]* &
    & 1.0_rk, [3, 4]))
  irr%plane_c = hardening_irr_planes(crys%slip_normal_c)
  call make_models()

  call check_init_state()
  call check_increments()
  call check_no_loops()

  call test_summary()

contains

  !> The models: bcc, m 0.05, gammadot_0 1 /s, tau_0 = g0, tight tolerances; with and without
  !> irradiation.
  !>
  !> Locals:
  !>   n  integer  number of slip systems
  subroutine make_models()

    integer :: n

    n = crys%num_slip
    call cp_model_init(model, stiff, crys%schmid_sym_c, crys%schmid_skew_c, &
      & spread(0.05_rk, 1, n), spread(1.0_rk, 1, n), spread(g0, 1, n), tol_newton=1.0e-14_rk, &
      & max_iters_state=200, tol_state=1.0e-14_rk, hard=hard, irr=irr)
    call cp_model_init(model_km, stiff, crys%schmid_sym_c, crys%schmid_skew_c, &
      & spread(0.05_rk, 1, n), spread(1.0_rk, 1, n), spread(g0, 1, n), tol_newton=1.0e-14_rk, &
      & max_iters_state=200, tol_state=1.0e-14_rk, hard=hard)

  end subroutine make_models

  !> State initialization with the phase's and with per-element loops.
  !>
  !> Locals:
  !>   h0  real(rk)  the diagonal of H_0 [1/mm^2]
  !>   gd  real(rk)  the initial defect term [MPa]
  subroutine check_init_state()

    real(rk) :: h0, gd

    call material_state_alloc(state, 1_ik, 2_ik, int(crys%num_slip, ik))
    call model%init_state([1_ik, 2_ik], state)
    h0 = 1.63e13_rk*3.0_rk*irr%loop_diameter*2.0_rk/3.0_rk
    call assert_close(state%dmg_c(:, 1, 2), [h0, h0, h0, 0.0_rk, 0.0_rk, 0.0_rk], &
      & 'init: H_0', abs_tol=1.0e-14_rk*h0)
    gd = hard%km_b*hard%km_mu_0*sqrt(irr%h_d*h0)
    call assert_close(state%crss(:, 1, 1), spread(g0 + hard%km_b*hard%km_mu_0* &
      & sqrt(hard%km_h_n*2.0e7_rk) + gd, 1, crys%num_slip), 'init: strengths (I1)', &
      & rel_tol=1.0e-13_rk)

    ! Per-element loops: densities and diameters.
    call cp_model_set_element_values(model, loop_density=[0.0_rk, 8.15e11_rk], &
      & loop_diameter=[2.0e-5_rk, 4.0e-5_rk])
    call model%init_state([1_ik, 2_ik], state)
    call assert_true(.not. any(abs(state%dmg_c(:, 1, 1)) > 0.0_rk), 'init: no loops, H_0 = 0')
    call assert_close(state%dmg_c(1, 1, 2), 8.15e11_rk*3.0_rk*4.0e-5_rk*2.0_rk/3.0_rk, &
      & 'init: per-element loops', rel_tol=1.0e-14_rk)
    call assert_close(state%crss(:, 1, 1), hardening_km_strengths(hard, &
      & spread(g0, 1, crys%num_slip), 2.0e7_rk), 'init: no loops, Kocks-Mecking strengths', &
      & rel_tol=0.0_rk)
    call make_models()

  end subroutine check_init_state

  !> Uniaxial increments on bcc [001]: convergence, the recursion in H, the strengths.
  !>
  !> Locals:
  !>   rot, sig, e, d, w  orientation, stress, elastic strain, {d}, w
  !>   rot_n, sig_n, e_n, wp, sec, h, tan  point-solve output
  !>   crss_prev, crss, rate  strengths [MPa], slip rates [1/s]
  !>   rho, dmg, dmg_prev, dmg_ref  density [1/mm^2], H [1/mm^2]
  !>   amat  real(rk)(6, 6)  the backward-Euler matrix
  !>   amat_ta  real(rk)(6, 6)  work array of the normal equations
  !>   contact, contact_prev  real(rk)(N)  N^a : H [1/mm^2]
  !>   contact_start  real(rk)  N^1 : H_0 [1/mm^2]
  !>   all_converged  logical
  !>   status, k, j, a, ios  integer
  subroutine check_increments()

    real(rk), parameter :: dt = 1.0_rk
    real(rk) :: rot(3, 3), sig(6), e(6), d(6), w(3), rot_n(3, 3), sig_n(6), e_n(6), wp(3)
    real(rk) :: sec(6, 6), h(6), tan(6, 6), rho, dmg(6), dmg_prev(6), dmg_ref(6), amat(6, 6), amat_ta(6, 6)
    real(rk) :: crss_prev(crys%num_slip), crss(crys%num_slip), rate(crys%num_slip)
    real(rk) :: contact(crys%num_slip), contact_prev(crys%num_slip), contact_start
    logical :: all_converged
    integer :: status, k, j, a, ios

    rot = 0.0_rk
    do j = 1, 3
      rot(j, j) = 1.0_rk
    end do
    sig = 0.0_rk
    e = 0.0_rk
    d = 1.0e-3_rk*[-0.5_rk, -0.5_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
    w = 0.0_rk
    rho = hard%km_rho_0
    dmg = hardening_irr_dmg_0(irr, irr%loop_density, irr%loop_diameter)
    crss_prev = hardening_irr_strengths(hard, irr, spread(g0, 1, crys%num_slip), rho, dmg)
    all_converged = .true.
    contact_start = dot_product(dmg, irr%plane_c(:, 1))
    do k = 1, 40
      dmg_prev = dmg
      contact_prev = matmul(dmg_prev, irr%plane_c)
      call cp_work_alloc(work, model)
      call cp_point_solve(model, crss_prev, rot, sig, e, d, w, dt, sig_n, e_n, wp, rot_n, &
        & crss, sec, h, tan, rate, status, km_rho=rho, dmg_c=dmg, work=work)
      all_converged = all_converged .and. status == cp_converged

      ! The backward-Euler recursion in H with the returned slip rates, and the strengths.
      amat = 0.0_rk
      do j = 1, 6
        amat(j, j) = 1.0_rk
      end do
      do a = 1, crys%num_slip
        amat = amat + dt*irr%eta*abs(rate(a))*(spread(irr%plane_c(:, a), 2, 6)* &
          & spread(irr%plane_c(:, a), 1, 6))
      end do
      call normal_equations_solve(amat, dmg_prev, dmg_ref, amat_ta, ios)
      call assert_close(dmg, dmg_ref, 'recursion in H', abs_tol=1.0e-9_rk*maxval(abs(dmg_prev)))
      call assert_close(crss, hardening_irr_strengths(hard, irr, spread(g0, 1, &
        & crys%num_slip), rho, dmg), 'strengths (I1)', rel_tol=1.0e-13_rk)
      contact = matmul(dmg, irr%plane_c)
      call assert_true(all(contact <= contact_prev + 1.0e-12_rk*maxval(abs(contact_prev))), &
        & 'contact densities do not increase')

      ! Next increment.
      crss_prev = crss
      rot = rot_n
      sig = sig_n
      e = e_n
    end do
    call assert_true(all_converged, 'every increment converged')
    call assert_true(minval(contact) < 0.7_rk*contact_start .and. minval(contact) > 0.0_rk, &
      & 'damage evolved, contact densities positive')

  end subroutine check_increments

  !> No loops: H stays 0, and the point solve equals Kocks-Mecking without irradiation exactly.
  !>
  !> Locals: as check_increments; crss2, sig2, rho2 the results without irradiation.
  subroutine check_no_loops()

    real(rk) :: rot(3, 3), sig(6), e(6), d(6), w(3), rot_n(3, 3), sig_n(6), e_n(6), wp(3)
    real(rk) :: sec(6, 6), h(6), tan(6, 6), rho, rho2, dmg(6), sig2(6)
    real(rk) :: crss_prev(crys%num_slip), crss(crys%num_slip), crss2(crys%num_slip)
    real(rk) :: rate(crys%num_slip)
    integer :: status, j, k

    rot = 0.0_rk
    do j = 1, 3
      rot(j, j) = 1.0_rk
    end do
    sig = 0.0_rk
    e = 0.0_rk
    d = 1.0e-3_rk*[-0.5_rk, -0.5_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
    w = 0.0_rk
    rho = hard%km_rho_0
    rho2 = rho
    dmg = 0.0_rk
    crss_prev = hardening_km_strengths(hard, spread(g0, 1, crys%num_slip), rho)
    ! Bit-identity is a property of the plain state loop: the damping (D62) changes the
    ! iterates (test_cp_state_damping checks it against the plain loop to tolerance).
    model%damp_state = .false.
    do k = 1, 10
      call cp_work_alloc(work, model)
      call cp_point_solve(model, crss_prev, rot, sig, e, d, w, 0.5_rk, sig_n, e_n, wp, rot_n, &
        & crss, sec, h, tan, rate, status, km_rho=rho, dmg_c=dmg, work=work)
      call cp_work_alloc(work, model_km)
      call cp_point_solve(model_km, crss_prev, rot, sig, e, d, w, 0.5_rk, sig2, e_n, wp, &
        & rot_n, crss2, sec, h, tan, rate, status, km_rho=rho2, work=work)
      call assert_close(crss, crss2, 'no loops: strengths identical to Kocks-Mecking', &
        & abs_tol=0.0_rk, rel_tol=0.0_rk)
      call assert_close(sig_n, sig2, 'no loops: stress identical to Kocks-Mecking', &
        & abs_tol=0.0_rk, rel_tol=0.0_rk)
      call assert_close(rho, rho2, 'no loops: density identical to Kocks-Mecking', &
        & abs_tol=0.0_rk, rel_tol=0.0_rk)
      call assert_true(.not. any(abs(dmg) > 0.0_rk), 'no loops: H stays 0')
      crss_prev = crss
      rot = rot_n
      sig = sig_n
      e = e_n
    end do
    model%damp_state = .true.

  end subroutine check_no_loops

end program test_cp_model_irr
