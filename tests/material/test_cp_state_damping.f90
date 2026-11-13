! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Tests of the Aitken damping of the state loop (step 7.3, D62,
!> dev/theory/addendum_state_damping.md §7).
!>
!> Source : FEPX (ROADMAP step 7.3); Kuttler & Wall (2008), Comput. Mech. 43:61-72, Eq. 44
!>          and Remark 9 (exact in the linear scalar case).
!>
!> Checks:
!>   - cp_state_damping_factor on the damped fixed-point iteration of a linear map
!>     x -> lambda x + c (scalar, and a 7-vector with the same slope in every component): with
!>     omega_1 = 1, the second factor is 1/(1 - lambda) (lambda = -3, -1, -0.5: 1e-14) and the
!>     third iterate is the fixed point (1e-13); lambda = 0.5 gives the cap 1, lambda = -30
!>     the factor 1/31 (inside the bounds), lambda = -2000 the bound 0.001, lambda = 2 (no
!>     overshoot, estimate -1) the factor 1 (D62 revision 3); an unchanged residual keeps the
!>     previous factor;
!>   - point solve (bcc, Barton's parameters, loops 3.61e12 /mm^3, default tolerances), 160
!>     increments of a mixed {d} through the stress drop: with and without damping every
!>     increment converges, and the final density, H, strengths and stress agree (relative
!>     1e-4: both are the same fixed point to tol_state).
!>
!> Locals:
!>   work   type(cp_work_t)  point-solve work arrays (cp_work_alloc before each call; D100)
!>   lambdas  real(rk)(3)  slopes of the converging checks [-]
!>   i        integer      loop index
program test_cp_state_damping

  use kinds_mod, only: rk
  use crystal_mod, only: crystal_t, crystal_init, crystal_bcc, stiffness_cubic, &
    & shear_modulus_hill
  use crystal_plasticity_model_mod, only: cp_model_t, cp_model_init, cp_point_solve, &
    & cp_state_damping_factor, cp_converged
  use crystal_plasticity_model_mod, only: cp_work_t, cp_work_alloc
  use hardening_mod, only: hardening_params_t, hardening_irr_t, hardening_kocks_mecking, &
    & hardening_irr_loop_average, hardening_irr_planes, hardening_irr_dmg_0, &
    & hardening_irr_strengths
  use assert_mod, only: assert_true, assert_close, test_summary

  implicit none

  type(cp_work_t) :: work
  real(rk), parameter :: lambdas(3) = [-3.0_rk, -1.0_rk, -0.5_rk]
  integer :: i

  ! Linear maps: the factor and the fixed point.
  do i = 1, size(lambdas)
    call check_linear(lambdas(i), 1)
    call check_linear(lambdas(i), 7)
  end do
  call check_bounds()

  ! Point solve with and without damping.
  call check_point_solve()

  call test_summary()

contains

  !> Damped iteration of x -> lambda x + c from x_1 = 0 (c = 1, ..., n): omega_1 = 1, then
  !> the factor of iteration 2 is 1/(1 - lambda) and x_3 is the fixed point c/(1 - lambda).
  !>
  !> Arguments:
  !>   lambda  (in) real(rk)  slope of the map [-]
  !>   n       (in) integer   dimension
  !> Locals:
  !>   c                  real(rk)(n)  the map's constant
  !>   x1, x2, x3         real(rk)(n)  iterates
  !>   r1, r2             real(rk)(n)  residuals G(x) - x
  !>   omega              real(rk)     factor of iteration 2
  !>   j                  integer      component index
  !>   label              character    check label
  subroutine check_linear(lambda, n)

    real(rk), intent(in) :: lambda
    integer, intent(in) :: n

    real(rk) :: c(n), x1(n), x2(n), x3(n), r1(n), r2(n), omega
    integer :: j
    character(len=64) :: label

    write (label, '(a,f5.2,a,i0)') 'linear map lambda ', lambda, ', n ', n
    c = [(real(j, rk), j=1, n)]

    ! Iteration 1 (omega 1), iteration 2 (Aitken), the third iterate.
    x1 = 0.0_rk
    r1 = lambda*x1 + c - x1
    x2 = x1 + r1
    r2 = lambda*x2 + c - x2
    omega = cp_state_damping_factor(r1, r2, 1.0_rk)
    call assert_close(omega, 1.0_rk/(1.0_rk - lambda), trim(label)//': factor 1/(1 - lambda)', &
      & rel_tol=1.0e-14_rk)
    x3 = x2 + omega*r2
    call assert_close(x3, c/(1.0_rk - lambda), trim(label)//': fixed point at x_3', &
      & rel_tol=1.0e-13_rk)

  end subroutine check_linear

  !> The bounds: lambda = 0.5 (Aitken 2) gives 1, lambda = -30 gives 1/31,
  !> lambda = -2000 (Aitken 1/2001) the bound 0.001, lambda = 2 (Aitken -1) the factor 1 (D62
  !> revision 3), and an
  !> unchanged residual keeps the previous factor.
  !>
  !> Locals:
  !>   r1, r2  real(rk)(1)  residuals of two iterations
  subroutine check_bounds()

    real(rk) :: r1(1), r2(1)

    ! lambda = 0.5 from x_1 = 0, c = 1: r1 = 1, r2 = 0.5.
    r1 = 1.0_rk
    r2 = 0.5_rk
    call assert_close(cp_state_damping_factor(r1, r2, 1.0_rk), 1.0_rk, 'cap 1 (lambda 0.5)', &
      & abs_tol=0.0_rk)
    ! lambda = -30: r2 = (1 + (lambda - 1)) r1 = -30.
    r2 = -30.0_rk
    call assert_close(cp_state_damping_factor(r1, r2, 1.0_rk), 1.0_rk/31.0_rk, &
      & 'inside the bounds (lambda -30)', rel_tol=1.0e-14_rk)
    ! lambda = -2000: r2 = -2000.
    r2 = -2000.0_rk
    call assert_close(cp_state_damping_factor(r1, r2, 1.0_rk), 0.001_rk, &
      & 'lower bound 0.001 (lambda -2000)', abs_tol=0.0_rk)
    ! lambda = 2 (no overshoot): r2 = 2, Aitken -1, factor 1.
    r2 = 2.0_rk
    call assert_close(cp_state_damping_factor(r1, r2, 0.3_rk), 1.0_rk, &
      & 'non-positive estimate gives 1 (lambda 2)', abs_tol=0.0_rk)
    ! Unchanged residual: the previous factor.
    call assert_close(cp_state_damping_factor(r1, r1, 0.3_rk), 0.3_rk, &
      & 'unchanged residual keeps the factor', abs_tol=0.0_rk)

  end subroutine check_bounds

  !> Point solve through the stress drop with and without damping (default tolerances).
  !>
  !> Locals:
  !>   crys        type(crystal_t)           the bcc crystal
  !>   stiff       real(rk)(6, 6)            stiffness [MPa]
  !>   model       type(cp_model_t)          the model with irradiation
  !>   hard        type(hardening_params_t)  Kocks-Mecking parameters (Barton, FEPX units)
  !>   irr         type(hardening_irr_t)     irradiation parameters
  !>   n           integer                   number of slip systems
  !>   damp       integer                   0: plain loop, 1: damped
  !>   k, j        integer                   increment, index
  !>   status      integer                   point-solve status
  !>   all_converged  logical                every increment converged
  !>   rot, rot_n  real(rk)(3, 3)            orientation at the start and end of an increment
  !>   sig, sig_n, e, e_n, d, h  real(rk)(6) stress [MPa], elastic strain [-], {d} [1/s], {h}
  !>   w, wp       real(rk)(3)               spins [1/s]
  !>   sec, tan    real(rk)(6, 6)            [s], [s_tan]
  !>   rho         real(rk)                  density [1/mm^2]
  !>   dmg         real(rk)(6)               H [1/mm^2]
  !>   crss_prev, crss, rate  real(rk)(:)    strengths [MPa], slip rates [1/s]
  !>   final_rho, final_dmg, final_sig, final_crss  the results of the two runs
  subroutine check_point_solve()

    type(crystal_t) :: crys
    real(rk) :: stiff(6, 6)
    type(cp_model_t) :: model
    type(hardening_params_t) :: hard
    type(hardening_irr_t) :: irr
    integer :: n, damp, k, j, status
    logical :: all_converged
    real(rk) :: rot(3, 3), rot_n(3, 3), sig(6), sig_n(6), e(6), e_n(6), d(6), h(6)
    real(rk) :: w(3), wp(3), sec(6, 6), tan(6, 6), rho, dmg(6)
    real(rk), allocatable :: crss_prev(:), crss(:), rate(:)
    real(rk) :: final_rho(0:1), final_dmg(6, 0:1), final_sig(6, 0:1)
    real(rk), allocatable :: final_crss(:, :)

    ! Barton's parameters (FEPX units), mu_0 by Hill, the loops of the FE study.
    call crystal_init(crys, crystal_bcc)
    stiff = stiffness_cubic(231400.0_rk, 134700.0_rk, 116400.0_rk)
    hard%law = hardening_kocks_mecking
    hard%km_b = 2.48e-7_rk
    hard%km_mu_0 = shear_modulus_hill(stiff)
    hard%km_h_n = 0.125_rk
    hard%km_rho_0 = 2.0e7_rk
    hard%km_k1 = 4.5e5_rk
    hard%km_k20 = 14.0_rk
    hard%km_gammadot_k0 = 1.0e10_rk
    hard%km_n = 100.0_rk
    hard%tol_hardening = 1.0e-12_rk
    irr%h_d = 0.675_rk
    irr%eta = 100.0_rk
    irr%loop_density = 3.61e12_rk
    irr%loop_diameter = 100.0_rk*hard%km_b
    irr%loop_avg_c = hardening_irr_loop_average(reshape([1, 1, 1, 1, 1, -1, 1, -1, 1, -1, &
      & 1, 1]*1.0_rk, [3, 4]))
    irr%plane_c = hardening_irr_planes(crys%slip_normal_c)
    n = crys%num_slip
    allocate (crss_prev(n), crss(n), rate(n), final_crss(n, 0:1))
    call cp_model_init(model, stiff, crys%schmid_sym_c, crys%schmid_skew_c, &
      & spread(0.05_rk, 1, n), spread(1.0_rk, 1, n), spread(90.0_rk, 1, n), hard=hard, irr=irr)

    ! The same loading without and with damping: a mixed {d} (shear components), 160
    ! increments of 0.05 s, through the stress drop.
    do damp = 0, 1
      model%damp_state = damp == 1
      rot = 0.0_rk
      do j = 1, 3
        rot(j, j) = 1.0_rk
      end do
      sig = 0.0_rk
      e = 0.0_rk
      d = 0.05_rk*[-0.5_rk, -0.5_rk, 1.0_rk, 0.3_rk, 0.2_rk, 0.1_rk]
      w = 0.0_rk
      rho = hard%km_rho_0
      dmg = hardening_irr_dmg_0(irr, irr%loop_density, irr%loop_diameter)
      crss_prev = hardening_irr_strengths(hard, irr, spread(90.0_rk, 1, n), rho, dmg)
      all_converged = .true.
      do k = 1, 160
        call cp_work_alloc(work, model)
        call cp_point_solve(model, crss_prev, rot, sig, e, d, w, 0.05_rk, sig_n, e_n, wp, &
          & rot_n, crss, sec, h, tan, rate, status, km_rho=rho, dmg_c=dmg, work=work)
        all_converged = all_converged .and. status == cp_converged
        crss_prev = crss
        rot = rot_n
        sig = sig_n
        e = e_n
      end do
      call assert_true(all_converged, 'point solve: every increment converged')
      final_rho(damp) = rho
      final_dmg(:, damp) = dmg
      final_sig(:, damp) = sig
      final_crss(:, damp) = crss
    end do

    ! The same fixed point to tol_state.
    call assert_close(final_rho(1), final_rho(0), 'damped vs plain: density', rel_tol=1.0e-4_rk)
    call assert_close(final_dmg(:, 1), final_dmg(:, 0), 'damped vs plain: H', &
      & abs_tol=1.0e-4_rk*maxval(abs(final_dmg(:, 0))))
    call assert_close(final_crss(:, 1), final_crss(:, 0), 'damped vs plain: strengths', &
      & rel_tol=1.0e-4_rk)
    call assert_close(final_sig(:, 1), final_sig(:, 0), 'damped vs plain: stress', &
      & abs_tol=1.0e-4_rk*maxval(abs(final_sig(:, 0))))

  end subroutine check_point_solve

end program test_cp_state_damping
