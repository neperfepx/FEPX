! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of bc_iteration_mod (step 8.2b): the values to reach on the load path and in time,
!> the convergence test, the perturbation signs and sizes, the correction of the face
!> velocities (a linear response: one correction lands on the path, or on the time-based
!> values, exactly), the Voigt-averaged moduli (isotropic constants, v2's hcp formula, two
!> phases) and the elastic velocity estimates, against hand values. Step 8.2c: the values of
!> rows given by mode (path, time, rate), the step kinds of iteration (B) (by magnitude, a zero
!> crossing unloading, a dwell), the scaling loop's test and factor (the secant lands on a
!> linear response), and the perturbation sizes of the Newton loop with their lower bound.
!>
!> Source : FEPX (ROADMAP step 8.2b; dev/theory/addendum_triaxial_load_control.md §3.4; v2
!>          driver_triaxcsr_mod.f90::driver_triaxcsr, driver_triax_utilities_mod.f90::
!>          est_avg_mod).
!>
!> Locals:
!>   jac        real(rk)(3, 2)   a linear response: change of the values (rows s, t, p) per
!>                               unit velocity change of the faces s and t [N s/mm]
!>   q0, q_new  real(rk)(3)      values before and after the correction [N]
!>   ratio      real(rk)(3)      path direction per control change [-]
!>   pert       real(rk)(2)      signed perturbations [mm/s]
!>   dq         real(rk)(3, 2)   value changes of the perturbed solves [N]
!>   dvel       real(rk)(2)      velocity corrections [mm/s]
!>   status     integer          solve status
!>   stiff      real(rk)(6, 6, 2)  Mandel stiffnesses of two phases [MPa]
!>   e, nu      real(rk)         elastic moduli [MPa], [-]
!>   e_ref, nu_ref  real(rk)     reference values
!>   c11, c12, c13, c33, c44, c66  real(rk)  hcp constants [MPa]
!>   f, g, h, k real(rk)         v2's means and bulk modulus [MPa]
!>   vel        real(rk)(3)      face velocities [mm/s]
!>   s3         real(rk)         control stress rate (time-based estimate) [MPa/s]
!>   a, b       real(rk)         a linear response q = a + b V of the control [N], [N s/mm]
!>   v1, v2s    real(rk)         two control velocities of the scaling loop [mm/s]
!>   floor      real(rk)         lower bound of the perturbations [mm/s]
program test_bc_iteration

  use kinds_mod, only: rk
  use bc_iteration_mod, only: bc_path_values, bc_time_values, bc_converged, bc_pert_signs, &
    & bc_correction, bc_pert_update, bc_elastic_moduli, bc_elastic_velocities, &
    & bc_initial_velocities, bc_mode_values, bc_step_action, bc_rate_first_converged, &
    & bc_rate_scale, bc_rate_pert_floor, bc_rate_pert_initial, bc_rate_pert_update, &
    & bc_mode_path, bc_mode_time, bc_mode_rate, bc_action_loading, bc_action_unloading, &
    & bc_action_dwell
  use assert_mod, only: assert_equal, assert_close, assert_true, test_summary

  implicit none

  real(rk) :: jac(3, 2), q0(3), q_new(3), ratio(3), pert(2), dq(3, 2), dvel(2)
  integer :: status
  real(rk) :: stiff(6, 6, 2), e, nu, e_ref, nu_ref
  real(rk) :: c11, c12, c13, c33, c44, c66, f, g, h, k
  real(rk) :: vel(3), s3
  real(rk) :: a, b, v1, v2s, floor

  ! --- Values to reach ---
  ! Path: control from 0 to -300 N, at -150 N (lambda 1/2); start values 10 and -5 N, targets
  ! 110 and 45 N: 60 and 20 N. From a measured start: control from -290 to -600 N, at -445 N.
  call assert_close(bc_path_values([10.0_rk, -5.0_rk], [110.0_rk, 45.0_rk], -150.0_rk, 0.0_rk, &
    & -300.0_rk), [60.0_rk, 20.0_rk], 'path values, lambda 1/2', abs_tol=1.0e-12_rk)
  call assert_close(bc_path_values([290.0_rk], [600.0_rk], -445.0_rk, -290.0_rk, -600.0_rk), &
    & [445.0_rk], 'path values from a measured start', abs_tol=1.0e-12_rk)
  ! Time: a quarter of the step from 20 to 100 MPa: 40 MPa; the end: the target.
  call assert_close(bc_time_values([20.0_rk, 0.0_rk], [100.0_rk, -8.0_rk], 0.25_rk), &
    & [40.0_rk, -2.0_rk], 'time values', abs_tol=1.0e-12_rk)
  call assert_close(bc_time_values([20.0_rk], [100.0_rk], 1.0_rk), [100.0_rk], 'time values end', &
    & abs_tol=0.0_rk)

  ! --- Convergence and signs ---
  call assert_true(bc_converged([10.05_rk, -3.0_rk], [10.0_rk, -3.09_rk], [0.1_rk, 0.1_rk]), &
    & 'converged within the tolerances')
  call assert_true(.not. bc_converged([10.05_rk, -3.0_rk], [10.0_rk, -3.11_rk], &
    & [0.1_rk, 0.1_rk]), 'not converged: one value outside')
  call assert_true(.not. bc_converged([10.05_rk], [10.0_rk], [0.01_rk]), &
    & 'not converged: tolerance of the row')
  call assert_close(bc_pert_signs([1.0_rk, 5.0_rk, 2.0_rk], [2.0_rk, 4.0_rk, 2.0_rk]), &
    & [1.0_rk, -1.0_rk, 1.0_rk], 'perturbation signs (towards; +1 when equal)', abs_tol=0.0_rk)

  ! --- Correction, path rows (v2's 3 x 3 system) ---
  ! A linear response q = q0 + jac dV (rows: s, t, control p). Start values 0, targets 100, 50,
  ! -300 N: after one correction the values lie on the path, q_s = (100/-300) q_p,
  ! q_t = (50/-300) q_p.
  jac = reshape([2.0_rk, 0.3_rk, -1.0_rk, 0.5_rk, 3.0_rk, -1.0_rk], [3, 2])
  q0 = [5.0_rk, -2.0_rk, -100.0_rk]
  ratio = [100.0_rk/(-300.0_rk), 50.0_rk/(-300.0_rk), 1.0_rk]
  pert = [0.1_rk, -0.2_rk]
  dq(:, 1) = jac(:, 1)*pert(1)
  dq(:, 2) = jac(:, 2)*pert(2)
  call bc_correction(dq, q0, [0.0_rk, 0.0_rk, 0.0_rk], ratio, pert, .true., dvel, status)
  call assert_equal(status, 0, 'path correction solved')
  q_new = q0 + matmul(jac, dvel)
  call assert_close(q_new(1:2), ratio(1:2)*q_new(3), 'path correction on the path', &
    & abs_tol=1.0e-10_rk)
  ! From nonzero start values: on the path from (10, 20, -290) to (110, 70, -600).
  ratio = [100.0_rk/(-310.0_rk), 50.0_rk/(-310.0_rk), 1.0_rk]
  call bc_correction(dq, q0, [10.0_rk, 20.0_rk, -290.0_rk], ratio, pert, .true., dvel, status)
  q_new = q0 + matmul(jac, dvel)
  call assert_close(q_new(1:2) - [10.0_rk, 20.0_rk], ratio(1:2)*(q_new(3) + 290.0_rk), &
    & 'path correction from start values', abs_tol=1.0e-10_rk)
  ! One path row (n = 1, a 2 x 2 system).
  call bc_correction(dq(2:3, 2:2), q0(2:3), [0.0_rk, 0.0_rk], [0.5_rk, 1.0_rk], pert(2:2), &
    & .true., dvel(1:1), status)
  q_new(2:3) = q0(2:3) + jac(2:3, 2)*dvel(1)
  call assert_close(q_new(2), 0.5_rk*q_new(3), 'one path row', abs_tol=1.0e-10_rk)

  ! --- Correction, time-based rows ---
  ! Rows s, t only: the corrected values equal the values to reach.
  call bc_correction(dq(1:2, :), q0(1:2), [30.0_rk, -12.0_rk], [0.0_rk, 0.0_rk], pert, .false., &
    & dvel, status)
  call assert_equal(status, 0, 'time correction solved')
  q_new(1:2) = q0(1:2) + matmul(jac(1:2, :), dvel)
  call assert_close(q_new(1:2), [30.0_rk, -12.0_rk], 'time correction reaches the values', &
    & abs_tol=1.0e-10_rk)
  ! A perturbation that changes nothing: singular, no correction.
  dq = 0.0_rk
  call bc_correction(dq(1:2, :), q0(1:2), [30.0_rk, -12.0_rk], [0.0_rk, 0.0_rk], pert, .false., &
    & dvel, status)
  call assert_equal(status, 1, 'singular system reported')
  call assert_close(dvel, [0.0_rk, 0.0_rk], 'singular: no correction', abs_tol=0.0_rk)

  ! --- Perturbation sizes (v2-fixed C58: bound by min_pert_frac |V_p|, any sign of V_p) ---
  call assert_close(bc_pert_update([2.0e-6_rk, -0.2_rk], -0.01_rk, 0.001_rk), &
    & [1.0e-5_rk, 0.2_rk], 'perturbation sizes', abs_tol=1.0e-18_rk)

  ! --- Voigt-averaged moduli ---
  ! Isotropic constants (lambda 1e5, mu 5e4 MPa): E = 1.3333e5 MPa, nu = 1/3.
  stiff = 0.0_rk
  call cubic(200000.0_rk, 100000.0_rk, 50000.0_rk, stiff(:, :, 1))
  call bc_elastic_moduli(stiff(:, :, 1:1), [1.0_rk], e, nu)
  call assert_close(e, 400000.0_rk/3.0_rk, 'isotropic E', rel_tol=1.0e-14_rk)
  call assert_close(nu, 1.0_rk/3.0_rk, 'isotropic nu', rel_tol=1.0e-14_rk)
  ! hcp (LEGACY_triaxial_csr constants) with v2's implied c33 = c11 + c12 - c13: v2's
  ! est_avg_mod formula with f = (2 c11 + c33)/3, g = (c12 + 2 c13)/3, h = (2 c44 + c66)/3 and
  ! K = (c11 + c12 + c13)/3.
  c11 = 1.6966e5_rk
  c12 = 0.8866e5_rk
  c13 = 0.6166e5_rk
  c44 = 0.4250e5_rk
  c33 = c11 + c12 - c13
  c66 = 0.5_rk*(c11 - c12)
  stiff(:, :, 2) = 0.0_rk
  stiff(1, 1, 2) = c11
  stiff(2, 2, 2) = c11
  stiff(3, 3, 2) = c33
  stiff(1, 2, 2) = c12
  stiff(1, 3, 2) = c13
  stiff(2, 3, 2) = c13
  stiff(2, 1, 2) = c12
  stiff(3, 1, 2) = c13
  stiff(3, 2, 2) = c13
  stiff(4, 4, 2) = 2.0_rk*c44
  stiff(5, 5, 2) = 2.0_rk*c44
  stiff(6, 6, 2) = 2.0_rk*c66
  f = (2.0_rk*c11 + c33)/3.0_rk
  g = (c12 + 2.0_rk*c13)/3.0_rk
  h = (2.0_rk*c44 + c66)/3.0_rk
  k = (c11 + c12 + c13)/3.0_rk
  e_ref = (f - g + 3.0_rk*h)*(f + 2.0_rk*g)/(2.0_rk*f + 3.0_rk*g + h)
  nu_ref = (3.0_rk*k - e_ref)/(6.0_rk*k)
  call bc_elastic_moduli(stiff(:, :, 2:2), [1.0_rk], e, nu)
  call assert_close(e, e_ref, 'hcp E (v2 est_avg_mod)', rel_tol=1.0e-14_rk)
  call assert_close(nu, nu_ref, 'hcp nu (v2 est_avg_mod)', rel_tol=1.0e-14_rk)
  ! Two phases, weights 1/4 and 3/4: the weighted averages of E and nu.
  call bc_elastic_moduli(stiff, [0.25_rk, 0.75_rk], e, nu)
  call assert_close(e, 0.25_rk*400000.0_rk/3.0_rk + 0.75_rk*e_ref, 'two phases E', &
    & rel_tol=1.0e-14_rk)
  call assert_close(nu, 0.25_rk/3.0_rk + 0.75_rk*nu_ref, 'two phases nu', rel_tol=1.0e-14_rk)

  ! --- Elastic velocities ---
  ! Uniaxial stress rate 10 MPa/s along z, E 1e5, nu 0.3, lengths 1, 2, 4 mm: V_z = 4e-4,
  ! V_x = -3e-5, V_y = -6e-5 mm/s.
  call assert_close(bc_elastic_velocities([0.0_rk, 0.0_rk, 10.0_rk], [1.0_rk, 2.0_rk, 4.0_rk], &
    & 1.0e5_rk, 0.3_rk), [-3.0e-5_rk, -6.0e-5_rk, 4.0e-4_rk], 'uniaxial Hooke velocities', &
    & abs_tol=1.0e-18_rk)
  ! Path (LEGACY_triaxial_csr, unit cube): targets (-300, 0, 300) MPa, control x at -1e-4 mm/s:
  ! the estimate is (1 + nu) s / E, scaled: (-1e-4, 0, 1e-4).
  vel = bc_initial_velocities([-300.0_rk, 0.0_rk, 300.0_rk], 1, -1.0e-4_rk, [1.0_rk, 1.0_rk, &
    & 1.0_rk], 1.0e5_rk, 0.3_rk, .true.)
  call assert_close(vel, [-1.0e-4_rk, 0.0_rk, 1.0e-4_rk], 'initial velocities, path', &
    & abs_tol=1.0e-18_rk)
  ! Time-based: control z at 1e-3 mm/s (length 1), y at 2 MPa/s (length 2), x free: the
  ! control's stress rate s3 = E 1e-3 + nu 2, then Hooke's law.
  vel = bc_initial_velocities([0.0_rk, 2.0_rk, 999.0_rk], 3, 1.0e-3_rk, [1.0_rk, 2.0_rk, 1.0_rk], &
    & 1.0e5_rk, 0.3_rk, .false.)
  s3 = 1.0e5_rk*1.0e-3_rk + 0.3_rk*2.0_rk
  call assert_close(vel, [-0.3_rk*(2.0_rk + s3)/1.0e5_rk, (2.0_rk - 0.3_rk*s3)*2.0_rk/1.0e5_rk, &
    & 1.0e-3_rk], 'initial velocities, time-based', rel_tol=1.0e-13_rk)

  ! --- Rows by mode (step 8.2c) ---
  ! Path at lambda 1/2, time at a quarter, a rate of -4 N/s after 2.5 s, from the start values
  ! 10, 20, 30 to the targets 110, 100 (the rate row's target unused).
  call assert_close(bc_mode_values([bc_mode_path, bc_mode_time, bc_mode_rate], &
    & [10.0_rk, 20.0_rk, 30.0_rk], [110.0_rk, 100.0_rk, 999.0_rk], [0.0_rk, 0.0_rk, -4.0_rk], &
    & 0.5_rk, 0.25_rk, 2.5_rk), [60.0_rk, 40.0_rk, 20.0_rk], 'values by mode', abs_tol=1.0e-12_rk)

  ! --- Step kinds of iteration (B), by magnitude (Q-T9) ---
  call assert_equal(bc_step_action(0.0_rk, 100.0_rk, 100.0_rk, .true., 0.1_rk), &
    & bc_action_dwell, 'rate 0: dwell')
  call assert_equal(bc_step_action(2.0_rk, 0.0_rk, 100.0_rk, .true., 0.1_rk), &
    & bc_action_loading, 'from 0: loading')
  call assert_equal(bc_step_action(-2.0_rk, 0.05_rk, -100.0_rk, .true., 0.1_rk), &
    & bc_action_loading, 'from within the tolerance of 0: loading')
  call assert_equal(bc_step_action(-2.0_rk, -100.0_rk, -375.0_rk, .true., 0.1_rk), &
    & bc_action_loading, 'compression, magnitude up: loading (v2: unloading, F5)')
  call assert_equal(bc_step_action(2.0_rk, -375.0_rk, 0.0_rk, .true., 0.1_rk), &
    & bc_action_unloading, 'compression to 0: unloading (v2: loading)')
  call assert_equal(bc_step_action(-2.0_rk, 100.0_rk, -50.0_rk, .true., 0.1_rk), &
    & bc_action_unloading, 'zero crossing: unloading')
  call assert_equal(bc_step_action(-2.0_rk, 100.0_rk, 0.0_rk, .false., 0.1_rk), &
    & bc_action_unloading, 'end not known, rate against the start value: unloading')
  call assert_equal(bc_step_action(2.0_rk, 100.0_rk, 0.0_rk, .false., 0.1_rk), &
    & bc_action_loading, 'end not known, rate with the start value: loading')

  ! --- Scaling loop of iteration (B) ---
  ! Accepted within max(tol, 0.1 |q* - q_init|): from 100 to 250 N, 14.9 N off (tol 5) yes, 15.1
  ! no; strict at the bound (v2).
  call assert_true(bc_rate_first_converged(235.1_rk, 250.0_rk, 100.0_rk, 5.0_rk), &
    & 'scaling loop: within 0.1 of the change')
  call assert_true(.not. bc_rate_first_converged(234.9_rk, 250.0_rk, 100.0_rk, 5.0_rk), &
    & 'scaling loop: outside 0.1 of the change')
  call assert_true(bc_rate_first_converged(104.0_rk, 100.0_rk, 99.0_rk, 5.0_rk), &
    & 'scaling loop: within the tolerance')
  call assert_true(.not. bc_rate_first_converged(105.0_rk, 100.0_rk, 99.0_rk, 5.0_rk), &
    & 'scaling loop: strict at the tolerance')
  call assert_close(bc_rate_scale(1, .false., 2.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk), 1.1_rk, &
    & 'first scale', abs_tol=0.0_rk)
  call assert_close(bc_rate_scale(1, .true., 2.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk), 0.9_rk, &
    & 'first scale in a dwell', abs_tol=0.0_rk)
  ! The secant on a linear response q = a + b V lands: V_prev 1, V 1.1, q* 150: V' = 1.4.
  a = 10.0_rk
  b = 100.0_rk
  v1 = 1.0_rk
  v2s = 1.1_rk
  call assert_close(v2s*bc_rate_scale(2, .false., v2s, v1, 150.0_rk, a + b*v2s, a + b*v1), &
    & 1.4_rk, 'secant on a linear response', rel_tol=1.0e-13_rk)

  ! --- Perturbation sizes of iteration (B)'s Newton loop ---
  ! Lower bound: 1e-3 x the elastic estimate of (100, 0, 250) MPa over 50 s, E 1e5, nu 0.3, unit
  ! lengths: V = (2 - 0.3 x 5, -0.3 x 7, 5 - 0.3 x 2)/1e5 mm/s, largest 4.4e-5: 4.4e-8.
  floor = bc_rate_pert_floor([100.0_rk, 0.0_rk, 250.0_rk], [1.0_rk, 1.0_rk, 1.0_rk], 1.0e5_rk, &
    & 0.3_rk, 50.0_rk, 1.0e-3_rk)
  call assert_close(floor, 4.4e-8_rk, 'perturbation floor', rel_tol=1.0e-13_rk)
  call assert_close(bc_rate_pert_initial(2.0e-3_rk, 1.0e-3_rk, .false., 0.0_rk), 1.0e-4_rk, &
    & 'first perturbation: 0.1 of the change', rel_tol=1.0e-14_rk)
  call assert_close(bc_rate_pert_initial(2.0e-3_rk, 1.9e-3_rk, .false., 0.0_rk), 2.0e-5_rk, &
    & 'first perturbation: 0.01 of the velocity', rel_tol=1.0e-14_rk)
  call assert_close(bc_rate_pert_initial(2.0e-3_rk, 1.0e-3_rk, .true., 0.0_rk), 1.0e-6_rk, &
    & 'first perturbation at a dwell start', abs_tol=0.0_rk)
  call assert_close(bc_rate_pert_initial(1.0e-18_rk, 0.0_rk, .false., floor), floor, &
    & 'first perturbation: the floor (vanishing velocities, F1)', abs_tol=0.0_rk)
  call assert_close(bc_rate_pert_update([1.0e-5_rk, -3.0e-5_rk, 2.0e-5_rk], floor), 3.0e-5_rk, &
    & 'next perturbation: largest change', abs_tol=0.0_rk)
  call assert_close(bc_rate_pert_update([1.0e-12_rk, -3.0e-12_rk], floor), floor, &
    & 'next perturbation: the floor', abs_tol=0.0_rk)

  call test_summary()

contains

  !> Mandel stiffness of a cubic crystal.
  !>
  !> Arguments:
  !>   c11, c12, c44  (in)  real(rk)         constants [MPa]
  !>   stiff          (out) real(rk)(6, 6)   Mandel stiffness [MPa]
  subroutine cubic(c11, c12, c44, stiff)

    real(rk), intent(in) :: c11, c12, c44
    real(rk), intent(out) :: stiff(6, 6)

    stiff = 0.0_rk
    stiff(1:3, 1:3) = c12
    stiff(1, 1) = c11
    stiff(2, 2) = c11
    stiff(3, 3) = c11
    stiff(4, 4) = 2.0_rk*c44
    stiff(5, 5) = 2.0_rk*c44
    stiff(6, 6) = 2.0_rk*c44

  end subroutine cubic

end program test_bc_iteration
