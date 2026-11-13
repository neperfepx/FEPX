! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of crystal_plasticity_model_mod (step 6.2): the per-point solve.
!>
!> Source : FEPX (ROADMAP step 6.2).
!>
!> Checks (fixed random seed; fcc and hcp, Cu and Ti constants, random orientations, rates of
!> order 1e-3/s, dt 0.1-1 s):
!>   - elastic limit: with strengths of 1e8 MPa (no slip) and tight tolerances, the point solve
!>     equals elastic_model_mod::elastic_point_solve (stress, elastic strain, w^p, orientation,
!>     [s], {h}) and [s_tan] = [s] (relative 1e-9);
!>   - plastic points (strengths near the stress level): the solve converges; the residual (R1)
!>     at the returned state, recomputed by cp_residual in the final crystal frame, is below the
!>     Newton tolerance (tight tolerances: 1e-12); sigma = [s]({d} - {h}) (R10, relative 1e-9);
!>     the secant and tangent differ (plastic flow); w^p of the residual equals (R2) evaluated
!>     with explicit 3x3 algebra (relative 1e-12);
!>   - tangent: beta [s_tan] approximates the derivative of the Kirchhoff stress with respect to
!>     {d} (central finite differences; w = 0 and e^e_0 = 0 so that only the neglected w^p and
!>     beta variations and the orientation change remain; relative 2e-2 of the largest entry);
!>   - frame invariance: rotating the start orientation, stress, elastic strain and {d}, w by a
!>     rotation Q rotates the stress, strain, w^p and [s] and premultiplies the orientation
!>     (relative 1e-9);
!>   - failure: point_max_iters_newton = 1 with a tolerance of 1e-30 returns cp_newton_failed with
!>     finite outputs (the last iterate is kept, D49);
!>   - steady state: an fcc crystal along [001] (d = rate (-1/2, -1/2, 1), w = 0, the elastic
!>     strain of the steady stress as e^e_0) gives sigma_33 = sqrt(6) crss (sqrt(6) rate /
!>     (8 gammadot_0))^m / beta (8 systems with Schmid factor 1/sqrt(6)), relative 1e-6;
!>   - state-loop convergence with latent hardening (step 6.4, addendum_hardening §5 test 2;
!>     H off-diagonal 1.4 and 5, iteration counts printed);
!>   - single slip under uniaxial stress with latent hardening (step 6.4: crystal [123], H
!>     off-diagonal 1.4 and 1.0; recursion, one active system, equal inactive strengths, the
!>     integrated single-slip relation);
!>   - hardening (step 6.4): fcc [001] at 1e-3/s over 40 increments with the isotropic law
!>     (n = 1): the strengths of every increment follow the backward-Euler recursion
!>     crss_k = (crss_(k-1) + c_k tau_s)/(1 + c_k), c_k = dt h_0 Gammadot_k/(tau_s - tau_0), with
!>     Gammadot_k the total slip rate the point solve returns (1e-8 relative); all systems stay
!>     equal; the strengths increase (by 3.4 MPa over 1 % strain) and stay below tau_s;
!>     sigma_33 - sigma_11 follows the flow stress of the current strength (2e-3; measured
!>     2.6e-4: the elastic rate of a hardening crystal is not zero);
!>   - lattice rotation by slip: simple shear L = gammadot e_x (x) e_y of an fcc crystal whose
!>     system 1 has its direction along x and its plane normal along y, at its steady stress:
!>     the slip spin cancels the continuum spin and the lattice does not rotate (rotation angle
!>     below 1e-2 of gammadot dt; with the slip spin's sign reversed it is 2 x w dt =
!>     gammadot dt).
!>
!> Locals:
!>   work   type(cp_work_t)  point-solve work arrays (cp_work_alloc before each call; D100)
!>   model            type(cp_model_t)   model under test (set by each check)
!>   crss_end         real(rk)(:)        end-of-increment strengths of a point solve [MPa]
!>   crys             type(crystal_t)    slip geometry
!>   stiff            real(rk)(6, 6)     stiffness [MPa]
!>   i                integer            repetition
program test_cp_model

  use kinds_mod, only: rk
  use tensor_mod, only: matrix_from_mandel, mandel_from_matrix, matrix_from_skew
  use slip_kinetics_mod, only: resolved_shear_stress, slip_power_law
  use rotation_mod, only: transform_crystal_to_sample_mandel, transform_sample_to_crystal_mandel, &
    & transform_sample_to_crystal_skew, transform_crystal_to_sample_mat6, &
    & transform_crystal_to_sample_skew
  use crystal_mod, only: crystal_t, crystal_init, crystal_fcc, crystal_hcp, stiffness_cubic, &
    & stiffness_hexagonal
  use elastic_model_mod, only: elastic_point_solve
  use crystal_plasticity_model_mod, only: cp_model_t, cp_model_init, cp_point_solve, &
    & cp_residual, cp_converged, cp_newton_failed
  use crystal_plasticity_model_mod, only: cp_work_t, cp_work_alloc
  use hardening_mod, only: hardening_params_t, hardening_isotropic, hardening_anisotropic, &
    & hardening_default_matrix
  use linalg_mod, only: normal_equations_solve
  use tensor_test_utils_mod, only: seed_random, random_rotation, random_symmetric
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary

  implicit none

  type(cp_work_t) :: work
  type(cp_model_t) :: model
  type(crystal_t) :: crys
  real(rk), allocatable :: crss_end(:)
  real(rk) :: stiff(6, 6)
  integer :: i

  call seed_random(20261010)

  ! fcc copper and hcp titanium, several random states each.
  call crystal_init(crys, crystal_fcc)
  stiff = stiffness_cubic(168400.0_rk, 121400.0_rk, 75400.0_rk)
  do i = 1, 5
    call check_elastic_limit(stiff, 'fcc')
    call check_plastic(stiff, 100.0_rk, 'fcc')
    call check_tangent(stiff, 100.0_rk, 'fcc')
    call check_frame(stiff, 100.0_rk, 'fcc')
  end do
  call check_failure(stiff)
  call check_steady_fcc001(stiff)
  call check_single_slip_shear(stiff)
  call check_hardening_increments(stiff)
  call check_state_loop_latent(stiff)
  call check_single_slip_latent(stiff)
  call crystal_init(crys, crystal_hcp, c_over_a=1.587_rk)
  stiff = stiffness_hexagonal(162400.0_rk, 92000.0_rk, 69000.0_rk, 180700.0_rk, 46700.0_rk)
  do i = 1, 5
    call check_elastic_limit(stiff, 'hcp')
    call check_plastic(stiff, 150.0_rk, 'hcp')
    call check_tangent(stiff, 150.0_rk, 'hcp')
    call check_frame(stiff, 150.0_rk, 'hcp')
  end do

  call test_summary()

contains

  !> Model of the current crystal: per-family-like varied m (0.01-0.03), gammadot_0 1, uniform
  !> strength crss, tight tolerances unless given.
  !>
  !> Arguments:
  !>   stiff_c   (in)  real(rk)(6, 6)  stiffness [MPa]
  !>   crss      (in)  real(rk)        strength [MPa]
  !>   tol       (in)  real(rk)        Newton tolerance [1/s]
  !>   iters     (in)  integer         Newton limit
  !> Locals:
  !>   n, a      integer               systems, system
  !>   m         real(rk)(:)           rate sensitivities
  subroutine make_model(stiff_c, crss, tol, iters)

    real(rk), intent(in) :: stiff_c(6, 6)
    real(rk), intent(in) :: crss
    real(rk), intent(in) :: tol
    integer, intent(in) :: iters

    integer :: n, a
    real(rk), allocatable :: m(:)

    n = crys%num_slip
    allocate (m(n))
    do a = 1, n
      m(a) = 0.01_rk + 0.01_rk*mod(a, 3)
    end do
    call cp_model_init(model, stiff_c, crys%schmid_sym_c, crys%schmid_skew_c, m, &
      & spread(1.0_rk, 1, n), spread(crss, 1, n), tol_newton=tol, max_iters_newton=iters, &
      & max_iters_state=200, tol_state=1.0e-14_rk)
    if (allocated(crss_end)) deallocate (crss_end)
    allocate (crss_end(n))

  end subroutine make_model

  !> A random start state and loading: orientation, a random start stress scaled so that its
  !> largest resolved shear stress is crss (the systems at the strength slip), its elastic strain,
  !> {d} of order 1e-3/s, w of order 1e-3/s.
  !>
  !> Arguments:
  !>   crss     (in)  real(rk)        strength scale [MPa]
  !>   rot0     (out) real(rk)(3, 3)  start orientation
  !>   sig0_s   (out) real(rk)(6)     start Cauchy stress, sample frame [MPa]
  !>   e0_s     (out) real(rk)(6)     start elastic strain, sample frame [-]
  !>   d_s      (out) real(rk)(6)     {d} [1/s]
  !>   w_s      (out) real(rk)(3)     w [1/s]
  !> Locals:
  !>   tau_c    real(rk)(6)           start Kirchhoff stress, crystal frame [MPa]
  !>   e_c      real(rk)(6)           its elastic strain [-]
  !>   beta     real(rk)              det(I + e) [-]
  !>   sk       real(rk)(3)           random numbers
  subroutine random_state(crss, rot0, sig0_s, e0_s, d_s, w_s)

    real(rk), intent(in) :: crss
    real(rk), intent(out) :: rot0(3, 3)
    real(rk), intent(out) :: sig0_s(6)
    real(rk), intent(out) :: e0_s(6)
    real(rk), intent(out) :: d_s(6)
    real(rk), intent(out) :: w_s(3)

    real(rk) :: tau_c(6), e_c(6), beta, sk(3)

    rot0 = random_rotation()
    tau_c = mandel_from_matrix(random_symmetric())
    tau_c = crss*tau_c/maxval(abs(matmul(tau_c, model%schmid_sym_c)))
    e_c = matmul(model%compl_el_c, tau_c)
    beta = det3(matrix_from_mandel(e_c))
    sig0_s = transform_crystal_to_sample_mandel(rot0, tau_c)/beta
    e0_s = transform_crystal_to_sample_mandel(rot0, e_c)
    d_s = 2.0e-3_rk*mandel_from_matrix(random_symmetric())
    call random_number(sk)
    w_s = 2.0e-3_rk*(sk - 0.5_rk)

  end subroutine random_state

  !> Elastic limit: equal to elastic_point_solve.
  !>
  !> Arguments:
  !>   stiff_c  (in)  real(rk)(6, 6)  stiffness [MPa]
  !>   name     (in)  character       label
  !> Locals: the start state (random_state) and the outputs of both solves.
  subroutine check_elastic_limit(stiff_c, name)

    real(rk), intent(in) :: stiff_c(6, 6)
    character(len=*), intent(in) :: name

    real(rk) :: rot0(3, 3), sig0(6), e0(6), d(6), w(3), dt
    real(rk) :: sig(6), e(6), wp(3), rot(3, 3), sec(6, 6), h(6), tan(6, 6)
    real(rk) :: sig_el(6), e_el(6), wp_el(3), rot_el(3, 3), sec_el(6, 6), h_el(6)
    real(rk), allocatable :: rate(:)
    integer :: status, status_el

    call make_model(stiff_c, 1.0e8_rk, 1.0e-14_rk, 100)
    allocate (rate(model%num_slip))
    call random_state(100.0_rk, rot0, sig0, e0, d, w)
    dt = 0.5_rk
    call cp_work_alloc(work, model)
    call cp_point_solve(model, model%tau_0, rot0, sig0, e0, d, w, dt, sig, e, wp, rot, crss_end, sec, h, &
      & tan, rate, status, work=work)
    call elastic_point_solve(stiff_c, rot0, e0, d, w, dt, 1.0e-14_rk, 200, sig_el, e_el, wp_el, &
      & rot_el, sec_el, h_el, status_el)
    call assert_equal(status, cp_converged, name // '_elastic_status')
    call assert_equal(status_el, 0, name // '_elastic_status_el')
    call assert_close(sig, sig_el, name // '_elastic_stress', abs_tol=1.0e-9_rk*maxval(abs(sig_el)))
    call assert_close(e, e_el, name // '_elastic_strain', abs_tol=1.0e-9_rk*maxval(abs(e_el)))
    call assert_close(wp, wp_el, name // '_elastic_spin', abs_tol=1.0e-9_rk*maxval(abs(wp_el)))
    call assert_close(rot, rot_el, name // '_elastic_rot', abs_tol=1.0e-12_rk)
    call assert_close(sec, sec_el, name // '_elastic_secant', abs_tol=1.0e-9_rk*maxval(abs(sec_el)))
    call assert_close(h, h_el, name // '_elastic_h', abs_tol=1.0e-9_rk*maxval(abs(h_el)))
    call assert_close(tan, sec, name // '_elastic_tangent', abs_tol=1.0e-9_rk*maxval(abs(sec)))
    call assert_true(maxval(abs(rate)) < 1.0e-30_rk, name // '_elastic_no_slip')

  end subroutine check_elastic_limit

  !> Plastic point: converged residual, (R10), secant /= tangent.
  !>
  !> Arguments:
  !>   stiff_c  (in)  real(rk)(6, 6)  stiffness [MPa]
  !>   crss     (in)  real(rk)        strength [MPa]
  !>   name     (in)  character       label
  subroutine check_plastic(stiff_c, crss, name)

    real(rk), intent(in) :: stiff_c(6, 6)
    real(rk), intent(in) :: crss
    character(len=*), intent(in) :: name

    real(rk) :: rot0(3, 3), sig0(6), e0(6), d(6), w(3), dt, residual_norm, rhs(6), spin_hat(3)
    real(rk) :: sig(6), e(6), wp(3), rot(3, 3), sec(6, 6), h(6), tan(6, 6), tau_c(6), beta
    real(rk), allocatable :: rate(:), coef(:)
    integer :: status

    call make_model(stiff_c, crss, 1.0e-12_rk, 100)
    allocate (rate(model%num_slip), coef(model%num_slip))
    call random_state(crss, rot0, sig0, e0, d, w)
    dt = 0.2_rk
    call cp_work_alloc(work, model)
    call cp_point_solve(model, model%tau_0, rot0, sig0, e0, d, w, dt, sig, e, wp, rot, crss_end, sec, h, &
      & tan, rate, status, work=work)
    call assert_equal(status, cp_converged, name // '_plastic_status')
    call assert_true(maxval(abs(rate)) > 1.0e-6_rk, name // '_plastic_slips')

    ! Residual (R1) at the returned state in its crystal frame.
    beta = det3(matrix_from_mandel(e))
    tau_c = transform_sample_to_crystal_mandel(rot, beta*sig)
    call cp_residual(model, tau_c, model%tau_0, transform_sample_to_crystal_mandel(rot, d), &
      & transform_sample_to_crystal_skew(rot, w), transform_sample_to_crystal_mandel(rot, e0), &
      & dt, rhs, residual_norm, spin_hat, coef, work%residual_rss, work%residual_sliprate, &
      & work%residual_coef_sec)
    call assert_true(residual_norm < 1.0e-10_rk, name // '_plastic_residual')

    ! w^p (R2) against explicit 3x3 algebra: w - (E P - P E) + (E E0 - E0 E)/(2 dt), with
    ! P = sum gammadot p (the Newton's slip rates), in the final crystal frame.
    call check_spin_hat(tau_c, transform_sample_to_crystal_mandel(rot, e0), &
      & transform_sample_to_crystal_skew(rot, w), dt, spin_hat, name)

    ! (R10) and a tangent stiffer than the secant (plastic flow).
    call assert_close(matmul(sec, d - h), sig, name // '_plastic_R10', &
      & abs_tol=1.0e-9_rk*maxval(abs(sig)))
    call assert_true(maxval(abs(tan - sec)) > 1.0e-3_rk*maxval(abs(sec)), &
      & name // '_plastic_tangent_differs')

  end subroutine check_plastic

  !> w^p of cp_residual against explicit 3x3 algebra (R2).
  !>
  !> Arguments:
  !>   tau_c     (in)  real(rk)(6)  Kirchhoff stress, crystal frame [MPa]
  !>   e0_c      (in)  real(rk)(6)  e^e_0, crystal frame [-]
  !>   w_c       (in)  real(rk)(3)  w, crystal frame [1/s]
  !>   dt        (in)  real(rk)     time increment [s]
  !>   spin_hat  (in)  real(rk)(3)  w^p returned by cp_residual [1/s]
  !>   name      (in)  character    label
  !> Locals:
  !>   ee, ee0, pp, ww  real(rk)(3, 3)  e^e, e^e_0, d^p, w as matrices
  !>   rss, rate, cs, ct  real(rk)(:)   slip quantities
  subroutine check_spin_hat(tau_c, e0_c, w_c, dt, spin_hat, name)

    real(rk), intent(in) :: tau_c(6), e0_c(6), w_c(3), dt, spin_hat(3)
    character(len=*), intent(in) :: name

    real(rk) :: ee(3, 3), ee0(3, 3), pp(3, 3), ww(3, 3), expected(3, 3)
    real(rk), allocatable :: rss(:), rate(:), cs(:), ct(:)

    allocate (rss(model%num_slip), rate(model%num_slip), cs(model%num_slip), ct(model%num_slip))
    rss = resolved_shear_stress(model%schmid_sym_c, tau_c)
    call slip_power_law(rss, model%tau_0, model%m, model%gammadot_0, model%rss_cutoff_newton, &
      & rate, cs, ct)
    ee = matrix_from_mandel(matmul(model%compl_el_c, tau_c))
    ee0 = matrix_from_mandel(e0_c)
    pp = matrix_from_mandel(matmul(model%schmid_sym_c, rate))
    ww = matrix_from_skew(w_c)
    expected = ww - (matmul(ee, pp) - matmul(pp, ee)) + (matmul(ee, ee0) - matmul(ee0, ee))/ &
      & (2.0_rk*dt)
    call assert_close(matrix_from_skew(spin_hat), expected, name // '_spin_hat_3x3', &
      & abs_tol=1.0e-12_rk*maxval(abs(expected)))

  end subroutine check_spin_hat

  !> Tangent: beta [s_tan] against finite differences of the Kirchhoff stress in {d}, with
  !> w = 0 and e^e_0 = 0.
  !>
  !> Arguments:
  !>   stiff_c  (in)  real(rk)(6, 6)  stiffness [MPa]
  !>   crss     (in)  real(rk)        strength [MPa]
  !>   name     (in)  character       label
  subroutine check_tangent(stiff_c, crss, name)

    real(rk), intent(in) :: stiff_c(6, 6)
    real(rk), intent(in) :: crss
    character(len=*), intent(in) :: name

    real(rk) :: rot0(3, 3), sig0(6), e0(6), d(6), w(3), dt, hstep, fd(6, 6), tau_p(6), tau_m(6)
    real(rk) :: sig(6), e(6), wp(3), rot(3, 3), sec(6, 6), h(6), tan(6, 6), dd(6), beta
    real(rk), allocatable :: rate(:)
    integer :: status, k

    call make_model(stiff_c, crss, 1.0e-13_rk, 100)
    allocate (rate(model%num_slip))
    call random_state(crss, rot0, sig0, e0, d, w)
    w = 0.0_rk
    e0 = 0.0_rk
    sig0 = 0.0_rk
    d = 2.0e-3_rk*mandel_from_matrix(random_symmetric())
    dt = 1.0_rk
    call cp_work_alloc(work, model)
    call cp_point_solve(model, model%tau_0, rot0, sig0, e0, d, w, dt, sig, e, wp, rot, crss_end, sec, h, &
      & tan, rate, status, work=work)
    call assert_equal(status, cp_converged, name // '_tangent_status')
    beta = det3(matrix_from_mandel(e))
    hstep = 1.0e-6_rk*maxval(abs(d))
    do k = 1, 6
      dd = 0.0_rk
      dd(k) = hstep
      tau_p = kirchhoff(rot0, sig0, e0, d + dd, w, dt)
      tau_m = kirchhoff(rot0, sig0, e0, d - dd, w, dt)
      fd(:, k) = (tau_p - tau_m)/(2.0_rk*hstep)
    end do
    call assert_close(beta*tan, fd, name // '_tangent_fd', abs_tol=2.0e-2_rk*maxval(abs(fd)))

  end subroutine check_tangent

  !> Kirchhoff stress (sample frame) of the point solve of the current model at {d}.
  !>
  !> Arguments:
  !>   rot0, sig0, e0, dv, w, dt  (in)  start orientation, stress, elastic strain, {d}, w, dt
  !>   tau_s  (out, result) real(rk)(6)  beta sigma [MPa]
  !> Locals: the outputs of the point solve.
  function kirchhoff(rot0, sig0, e0, dv, w, dt) result(tau_s)

    real(rk), intent(in) :: rot0(3, 3), sig0(6), e0(6), dv(6), w(3), dt
    real(rk) :: tau_s(6)

    real(rk) :: s(6), ee(6), wq(3), r(3, 3), s1(6, 6), hh(6), t1(6, 6)
    real(rk), allocatable :: rate(:)
    integer :: st

    allocate (rate(model%num_slip))
    call cp_work_alloc(work, model)
    call cp_point_solve(model, model%tau_0, rot0, sig0, e0, dv, w, dt, s, ee, wq, r, crss_end, s1, hh, &
      & t1, rate, st, work=work)
    tau_s = det3(matrix_from_mandel(ee))*s

  end function kirchhoff

  !> Frame invariance under a rotation q of the sample frame.
  !>
  !> Arguments:
  !>   stiff_c  (in)  real(rk)(6, 6)  stiffness [MPa]
  !>   crss     (in)  real(rk)        strength [MPa]
  !>   name     (in)  character       label
  subroutine check_frame(stiff_c, crss, name)

    real(rk), intent(in) :: stiff_c(6, 6)
    real(rk), intent(in) :: crss
    character(len=*), intent(in) :: name

    real(rk) :: rot0(3, 3), sig0(6), e0(6), d(6), w(3), dt, q(3, 3)
    real(rk) :: sig(6), e(6), wp(3), rot(3, 3), sec(6, 6), h(6), tan(6, 6)
    real(rk) :: sig_q(6), e_q(6), wp_q(3), rot_q(3, 3), sec_q(6, 6), h_q(6), tan_q(6, 6)
    real(rk), allocatable :: rate(:)
    integer :: status

    call make_model(stiff_c, crss, 1.0e-12_rk, 100)
    allocate (rate(model%num_slip))
    call random_state(crss, rot0, sig0, e0, d, w)
    dt = 0.3_rk
    call cp_work_alloc(work, model)
    call cp_point_solve(model, model%tau_0, rot0, sig0, e0, d, w, dt, sig, e, wp, rot, crss_end, sec, h, &
      & tan, rate, status, work=work)

    ! q as a "crystal-to-sample" map takes sample-frame quantities to the rotated sample frame.
    q = random_rotation()
    call cp_work_alloc(work, model)
    call cp_point_solve(model, model%tau_0, matmul(q, rot0), &
      & transform_crystal_to_sample_mandel(q, sig0), transform_crystal_to_sample_mandel(q, e0), &
      & transform_crystal_to_sample_mandel(q, d), transform_crystal_to_sample_skew(q, w), dt, &
      & sig_q, e_q, wp_q, rot_q, crss_end, sec_q, h_q, tan_q, rate, status, work=work)
    call assert_close(sig_q, transform_crystal_to_sample_mandel(q, sig), name // '_frame_stress', &
      & abs_tol=1.0e-9_rk*maxval(abs(sig)))
    call assert_close(e_q, transform_crystal_to_sample_mandel(q, e), name // '_frame_strain', &
      & abs_tol=1.0e-9_rk*maxval(abs(e)))
    call assert_close(wp_q, transform_crystal_to_sample_skew(q, wp), name // '_frame_spin', &
      & abs_tol=1.0e-9_rk*maxval(abs(wp)))
    call assert_close(rot_q, matmul(q, rot), name // '_frame_rot', abs_tol=1.0e-12_rk)
    call assert_close(sec_q, transform_crystal_to_sample_mat6(q, sec), name // '_frame_secant', &
      & abs_tol=1.0e-9_rk*maxval(abs(sec)))
    call assert_close(h_q, transform_crystal_to_sample_mandel(q, h), name // '_frame_h', &
      & abs_tol=1.0e-9_rk*maxval(abs(h)))

  end subroutine check_frame

  !> Failure path: one Newton iteration with an unreachable tolerance.
  !>
  !> Arguments:
  !>   stiff_c  (in)  real(rk)(6, 6)  stiffness [MPa]
  subroutine check_failure(stiff_c)

    real(rk), intent(in) :: stiff_c(6, 6)

    real(rk) :: rot0(3, 3), sig0(6), e0(6), d(6), w(3)
    real(rk) :: sig(6), e(6), wp(3), rot(3, 3), sec(6, 6), h(6), tan(6, 6)
    real(rk), allocatable :: rate(:)
    integer :: status

    call make_model(stiff_c, 100.0_rk, 1.0e-30_rk, 1)
    allocate (rate(model%num_slip))
    call random_state(100.0_rk, rot0, sig0, e0, d, w)
    call cp_work_alloc(work, model)
    call cp_point_solve(model, model%tau_0, rot0, sig0, e0, d, w, 0.2_rk, sig, e, wp, rot, crss_end, sec, &
      & h, tan, rate, status, work=work)
    call assert_equal(status, cp_newton_failed, 'failure_status')
    call assert_true(all(abs(sig) < huge(1.0_rk)) .and. all(abs(sec) < huge(1.0_rk)), &
      & 'failure_finite_outputs')

  end subroutine check_failure

  !> Steady state of an fcc crystal along [001] (identity orientation).
  !>
  !> Arguments:
  !>   stiff_c  (in)  real(rk)(6, 6)  stiffness [MPa]
  subroutine check_steady_fcc001(stiff_c)

    real(rk), intent(in) :: stiff_c(6, 6)

    real(rk) :: rot0(3, 3), sig0(6), e0(6), d(6), w(3), tau_ss, crss, rate_e, beta
    real(rk) :: sig(6), e(6), wp(3), rot(3, 3), sec(6, 6), h(6), tan(6, 6)
    real(rk), allocatable :: rate(:)
    integer :: status, k

    ! Uniform m = 0.02, gammadot_0 = 1, crss = 100 MPa.
    crss = 100.0_rk
    call cp_model_init(model, stiff_c, crys%schmid_sym_c, crys%schmid_skew_c, &
      & spread(0.02_rk, 1, crys%num_slip), spread(1.0_rk, 1, crys%num_slip), &
      & spread(crss, 1, crys%num_slip), tol_newton=1.0e-14_rk, tol_state=1.0e-14_rk)
    if (allocated(crss_end)) deallocate (crss_end)
    allocate (crss_end(crys%num_slip))
    allocate (rate(model%num_slip))
    rate_e = 1.0e-3_rk
    tau_ss = sqrt(6.0_rk)*crss*(sqrt(6.0_rk)*rate_e/8.0_rk)**0.02_rk

    ! Start at the steady stress (uniaxial along 3) with its elastic strain: the elastic rate
    ! is then zero and d = d^p (isochoric, uniaxial).
    rot0 = 0.0_rk
    do k = 1, 3
      rot0(k, k) = 1.0_rk
    end do
    e0 = matmul(model%compl_el_c, [0.0_rk, 0.0_rk, tau_ss, 0.0_rk, 0.0_rk, 0.0_rk])
    beta = det3(matrix_from_mandel(e0))
    sig0 = [0.0_rk, 0.0_rk, tau_ss, 0.0_rk, 0.0_rk, 0.0_rk]/beta
    d = rate_e*[-0.5_rk, -0.5_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
    w = 0.0_rk
    call cp_work_alloc(work, model)
    call cp_point_solve(model, model%tau_0, rot0, sig0, e0, d, w, 0.5_rk, sig, e, wp, rot, crss_end, sec, &
      & h, tan, rate, status, work=work)
    call assert_equal(status, cp_converged, 'steady_status')
    call assert_close(sig(3), tau_ss/det3(matrix_from_mandel(e)), 'steady_fcc001_sigma33', &
      & rel_tol=1.0e-6_rk)
    call assert_equal(count(abs(rate) > 1.0e-6_rk), 8, 'steady_fcc001_eight_systems')

  end subroutine check_steady_fcc001

  !> Simple shear of an fcc crystal oriented for single slip on system 1: no lattice rotation.
  !>
  !> Arguments:
  !>   stiff_c  (in)  real(rk)(6, 6)  stiffness [MPa]
  !> Locals:
  !>   rot0        real(rk)(3, 3)  orientation: slip direction to x, normal to y
  !>   gdot, crss  real(rk)        shear rate [1/s], strength [MPa]
  !>   tau_c       real(rk)(6)     steady Kirchhoff stress, crystal frame [MPa]
  !>   angle       real(rk)        rotation angle of the increment [rad]
  !>   dt          real(rk)        time increment [s]
  subroutine check_single_slip_shear(stiff_c)

    real(rk), intent(in) :: stiff_c(6, 6)

    real(rk) :: rot0(3, 3), sig0(6), e0(6), d(6), w(3), gdot, crss, tau_c(6), angle, dt, beta
    real(rk) :: sig(6), e(6), wp(3), rot(3, 3), sec(6, 6), h(6), tan(6, 6), dn(3)
    real(rk), allocatable :: rate(:)
    integer :: status

    ! Uniform m = 0.02, gammadot_0 = 1, crss = 100 MPa.
    crss = 100.0_rk
    gdot = 1.0e-3_rk
    dt = 0.1_rk
    call cp_model_init(model, stiff_c, crys%schmid_sym_c, crys%schmid_skew_c, &
      & spread(0.02_rk, 1, crys%num_slip), spread(1.0_rk, 1, crys%num_slip), &
      & spread(crss, 1, crys%num_slip), tol_newton=1.0e-14_rk, tol_state=1.0e-14_rk)
    if (allocated(crss_end)) deallocate (crss_end)
    allocate (crss_end(crys%num_slip))
    allocate (rate(model%num_slip))

    ! Orientation: rows of the crystal-to-sample matrix map d -> x, n -> y, d x n -> z.
    dn = [crys%slip_dir_c(2, 1)*crys%slip_normal_c(3, 1) - crys%slip_dir_c(3, 1)* &
      & crys%slip_normal_c(2, 1), crys%slip_dir_c(3, 1)*crys%slip_normal_c(1, 1) - &
      & crys%slip_dir_c(1, 1)*crys%slip_normal_c(3, 1), crys%slip_dir_c(1, 1)* &
      & crys%slip_normal_c(2, 1) - crys%slip_dir_c(2, 1)*crys%slip_normal_c(1, 1)]
    rot0(1, :) = crys%slip_dir_c(:, 1)
    rot0(2, :) = crys%slip_normal_c(:, 1)
    rot0(3, :) = dn

    ! Steady stress of system 1 (tau = crss (gdot/gammadot_0)^m) with its elastic strain.
    tau_c = crss*(gdot)**0.02_rk*crys%schmid_sym_c(:, 1)/dot_product(crys%schmid_sym_c(:, 1), &
      & crys%schmid_sym_c(:, 1))
    e0 = transform_crystal_to_sample_mandel(rot0, matmul(model%compl_el_c, tau_c))
    beta = det3(matrix_from_mandel(e0))
    sig0 = transform_crystal_to_sample_mandel(rot0, tau_c)/beta

    ! L = gdot e_x (x) e_y: d12 = gdot/2 (Mandel sqrt2 d12), w12 = gdot/2.
    d = [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, sqrt(2.0_rk)*0.5_rk*gdot]
    w = [0.0_rk, 0.0_rk, 0.5_rk*gdot]
    call cp_work_alloc(work, model)
    call cp_point_solve(model, model%tau_0, rot0, sig0, e0, d, w, dt, sig, e, wp, rot, crss_end, sec, h, &
      & tan, rate, status, work=work)
    call assert_equal(status, cp_converged, 'shear_status')
    call assert_true(abs(rate(1)) > 0.5_rk*gdot, 'shear_system1_slips')

    ! Rotation angle of rot rot0^T.
    angle = rotation_angle(matmul(rot, transpose(rot0)))
    call assert_true(angle < 1.0e-2_rk*gdot*dt, 'shear_no_lattice_rotation')

  end subroutine check_single_slip_shear

  !> Hardening over increments: fcc [001], isotropic law (n = 1), state passed from increment to
  !> increment; strengths against the backward-Euler recursion with the returned slip rates.
  !>
  !> Arguments:
  !>   stiff_c  (in)  real(rk)(6, 6)  stiffness [MPa]
  !> Locals:
  !>   hard          type(hardening_params_t)  the law (h_0 200, tau_s 330, n 1)
  !>   rot, sig, e   orientation, stress, elastic strain carried from increment to increment
  !>   tau_0, crss_ref  real(rk)  initial strength, recursion value [MPa]
  !>   gamma_tot     real(rk)  total slip rate of the increment [1/s]
  !>   k             integer   increment
  subroutine check_hardening_increments(stiff_c)

    real(rk), intent(in) :: stiff_c(6, 6)

    type(hardening_params_t) :: hard
    real(rk) :: rot(3, 3), sig(6), e(6), d(6), w(3), dt, tau_0, crss_ref, c, gamma_tot, tau_flow
    real(rk) :: rot_n(3, 3), sig_n(6), e_n(6), wp(3), sec(6, 6), h(6), tan(6, 6), beta
    real(rk), allocatable :: rate(:), crss_prev(:)
    integer :: status, k, j
    logical :: all_converged

    tau_0 = 100.0_rk
    hard%law = hardening_isotropic
    hard%h_0 = spread(200.0_rk, 1, crys%num_slip)
    hard%tau_s = spread(330.0_rk, 1, crys%num_slip)
    hard%n = spread(1.0_rk, 1, crys%num_slip)
    hard%tol_hardening = 1.0e-13_rk
    call cp_model_init(model, stiff_c, crys%schmid_sym_c, crys%schmid_skew_c, &
      & spread(0.02_rk, 1, crys%num_slip), spread(1.0_rk, 1, crys%num_slip), &
      & spread(tau_0, 1, crys%num_slip), tol_newton=1.0e-14_rk, tol_state=1.0e-14_rk, hard=hard)
    allocate (rate(model%num_slip), crss_prev(model%num_slip))
    if (allocated(crss_end)) deallocate (crss_end)
    allocate (crss_end(model%num_slip))

    ! Identity orientation, zero start state, uniaxial (isochoric) loading along 3.
    rot = 0.0_rk
    do j = 1, 3
      rot(j, j) = 1.0_rk
    end do
    sig = 0.0_rk
    e = 0.0_rk
    d = 1.0e-3_rk*[-0.5_rk, -0.5_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
    w = 0.0_rk
    dt = 0.25_rk
    crss_prev = model%tau_0
    crss_ref = tau_0
    all_converged = .true.
    do k = 1, 40
      call cp_work_alloc(work, model)
      call cp_point_solve(model, crss_prev, rot, sig, e, d, w, dt, sig_n, e_n, wp, rot_n, &
        & crss_end, sec, h, tan, rate, status, work=work)
      all_converged = all_converged .and. status == cp_converged

      ! Recursion with the returned total slip rate.
      gamma_tot = sum(abs(rate))
      c = dt*hard%h_0(1)*gamma_tot/(hard%tau_s(1) - tau_0)
      crss_ref = (crss_ref + c*hard%tau_s(1))/(1.0_rk + c)
      call assert_close(crss_end(1), crss_ref, 'hardening_recursion', rel_tol=1.0e-8_rk)
      call assert_close(crss_end, spread(crss_end(1), 1, model%num_slip), &
        & 'hardening_all_systems_equal', rel_tol=1.0e-14_rk)

      ! Next increment from this one's end state.
      crss_prev = crss_end
      rot = rot_n
      sig = sig_n
      e = e_n
    end do
    call assert_true(all_converged, 'hardening_increments_converged')
    call assert_true(crss_end(1) > tau_0 + 1.0_rk .and. crss_end(1) < hard%tau_s(1), &
      & 'hardening_increased_below_saturation')

    ! The deviatoric stress follows the flow stress of the current strength (8 systems,
    ! Schmid factor 1/sqrt6; isochoric loading, so sigma_33 - sigma_11 carries the flow).
    beta = det3(matrix_from_mandel(e))
    tau_flow = sqrt(6.0_rk)*crss_end(1)*(sqrt(6.0_rk)*1.0e-3_rk/8.0_rk)**0.02_rk
    call assert_close(sig(3) - sig(1), tau_flow/beta, 'hardening_flow_stress', rel_tol=2.0e-3_rk)

  end subroutine check_hardening_increments

  !> Addendum_hardening §5 test 2: state-loop convergence with latent hardening. fcc, random
  !> orientations and deformation rates (random_state), 40 increments of 0.25 s from a zero
  !> stress, the anisotropic law with per-system h_0 (150-250 MPa), tau_s 330 MPa, n 1.5 and the
  !> default interaction matrix with off-diagonal 1.4 and 5, then the isotropic law for
  !> comparison; production tolerances (cp_model_init and hardening defaults). Every increment
  !> must converge; the largest number of state iterations of each law is printed (the report).
  !>
  !> Arguments:
  !>   stiff_c  (in)  real(rk)(6, 6)  stiffness [MPa]
  !> Locals:
  !>   hard         type(hardening_params_t)  the law
  !>   latent       real(rk)(3)   off-diagonal entries of H (the third unused: isotropic)
  !>   rot, sig, e  orientation, stress, elastic strain carried from increment to increment
  !>   rot0, sig0, e0  random_state output (start orientation; stress and strain unused)
  !>   d, w         deformation rate [1/s] and spin [1/s] of the sample
  !>   rot_n, sig_n, e_n, wp, sec, h, tan  point-solve output
  !>   crss_prev, rate  strengths [MPa], slip rates [1/s]
  !>   n, a, k, l, trial, it, status, max_its  counts, indices, status, state iterations
  !>   all_converged  logical  every increment converged
  subroutine check_state_loop_latent(stiff_c)

    real(rk), intent(in) :: stiff_c(6, 6)

    type(hardening_params_t) :: hard
    real(rk), parameter :: latent(3) = [1.4_rk, 5.0_rk, 1.4_rk]
    real(rk) :: rot(3, 3), sig(6), e(6), rot0(3, 3), sig0(6), e0(6), d(6), w(3)
    real(rk) :: rot_n(3, 3), sig_n(6), e_n(6), wp(3), sec(6, 6), h(6), tan(6, 6)
    real(rk), allocatable :: crss_prev(:), rate(:)
    integer :: n, a, k, l, trial, it, status, max_its
    logical :: all_converged

    n = crys%num_slip
    allocate (crss_prev(n), rate(n))
    if (allocated(crss_end)) deallocate (crss_end)
    allocate (crss_end(n))

    ! l = 1, 2: anisotropic with off-diagonal latent(l); l = 3: isotropic.
    do l = 1, 3
      hard%h_0 = [(150.0_rk + 100.0_rk*real(a - 1, rk)/real(n - 1, rk), a = 1, n)]
      hard%tau_s = spread(330.0_rk, 1, n)
      hard%n = spread(1.5_rk, 1, n)
      if (l < 3) then
        hard%law = hardening_anisotropic
        hard%interaction = spread(spread(latent(l), 1, n), 1, n)
        do a = 1, n
          hard%interaction(a, a) = 1.0_rk
        end do
      else
        hard%law = hardening_isotropic
        hard%interaction = hardening_default_matrix(n)
      end if
      call cp_model_init(model, stiff_c, crys%schmid_sym_c, crys%schmid_skew_c, &
        & spread(0.02_rk, 1, n), spread(1.0_rk, 1, n), spread(100.0_rk, 1, n), hard=hard)

      ! Random orientations and loadings, from a zero state.
      max_its = 0
      all_converged = .true.
      do trial = 1, 5
        call random_state(100.0_rk, rot0, sig0, e0, d, w)
        rot = rot0
        sig = 0.0_rk
        e = 0.0_rk
        crss_prev = model%tau_0
        do k = 1, 40
          call cp_work_alloc(work, model)
          call cp_point_solve(model, crss_prev, rot, sig, e, d, w, 0.25_rk, sig_n, e_n, wp, &
            & rot_n, crss_end, sec, h, tan, rate, status, num_state_iters=it, work=work)
          all_converged = all_converged .and. status == cp_converged
          max_its = max(max_its, it)
          crss_prev = crss_end
          rot = rot_n
          sig = sig_n
          e = e_n
        end do
      end do

      ! Report and check.
      if (l < 3) then
        write (*, '(a, f4.1, a, i0)') 'state loop, anisotropic, H off-diagonal ', latent(l), &
          & ': largest number of state iterations ', max_its
      else
        write (*, '(a, i0)') 'state loop, isotropic: largest number of state iterations ', &
          & max_its
      end if
      call assert_true(all_converged, 'state_loop_latent_converged')
    end do

  end subroutine check_state_loop_latent

  !> Single slip under uniaxial stress with latent hardening (step 6.4, user suggestion
  !> 2026-10-10; addendum_hardening §3, §5): fcc with the tensile axis along crystal [123]
  !> (largest Schmid factor 0.467 on one system a, next 0.350), loaded at d_33 = 1e-3 1/s to
  !> 3 % in 300 increments, w = 0; the other five components of {d} are found by a quasi-Newton
  !> iteration with [s_tan] so that every stress component except sigma_33 vanishes (uniaxial
  !> stress: the crystal is free to shear, so one system can carry the deformation). Anisotropic
  !> law, h_0 200, tau_s 330, tau_0 100 MPa, n = 1, m 0.02, H = 1 on the diagonal and latent (1.4,
  !> then 1.0 for comparison) off it. Checks (latent 1.4):
  !>   - uniaxial stress (off-axis components < 1e-9 sigma_33 + 1e-9 MPa);
  !>   - every strength follows the backward-Euler recursion with the returned slip rates,
  !>     crss = (crss_prev + c tau_s)/(1 + c), c = dt h_0 (H |gammadot|)/(tau_s - tau_0) (n = 1),
  !>     1e-9;
  !>   - single slip: every other system slips at most 1e-5 of system a, and that ratio falls by
  !>     at least 5x over the run (the inactive systems harden 1.4x faster, so their resolved
  !>     stress ratios fall);
  !>   - the 11 inactive strengths are equal (same drive H_ba |gammadot^a|), 1e-6;
  !>   - the integrated relation of single slip,
  !>     tau_s - crss^b = (tau_s - tau_0)((tau_s - crss^a)/(tau_s - tau_0))^1.4
  !>     (n = 1), to the time-discretization error, 3e-5 relative in tau_s - crss^b (measured
  !>     2.8e-6).
  !> With latent 1.0 every strength stays equal (2e-6), system a ends at the same strength as
  !> with 1.4 (it hardens by its own slip only, H_aa = 1; 1e-6), and the slip ratio changes
  !> little (lattice rotation only); the ratios of both runs are printed.
  !>
  !> Arguments:
  !>   stiff_c  (in)  real(rk)(6, 6)  stiffness [MPa]
  !> Locals:
  !>   hard          type(hardening_params_t)  the law
  !>   latent        real(rk)(2)    off-diagonal entries of H (1.4, 1.0)
  !>   axis, e1, e2  real(rk)(3)    tensile axis [123] and two normals to it, crystal frame
  !>   rot0          real(rk)(3, 3) start orientation (rows e1, e2, axis: the axis to sample z)
  !>   schmid        real(rk)(N)    Schmid factors of the systems for the axis
  !>   rot, sig, e, crss_prev  state carried from increment to increment
  !>   rot_n, sig_n, e_n, wp, sec, h, tan, rate  point-solve output
  !>   d, w          deformation rate and spin [1/s] (d_33 prescribed)
  !>   jff, step     real(rk)(5, 5), (5)  free block of [s_tan], quasi-Newton step
  !>   jff_tj        real(rk)(5, 5)       work array of the normal equations
  !>   crss_a_latent real(rk)      strength of system a at the end with latent 1.4 [MPa]
  !>   c, crss_ref   real(rk)(N)    recursion coefficient and strengths
  !>   leak, leak_first, leak_last  real(rk)  largest other slip rate / |gammadot^a| (first and
  !>                                last once system a slips; -1 before)
  !>   n, a, b, k, j, it, l, status, ios  counts, indices, status
  !>   free          integer(5)     the free components of {d} (Mandel 11 22 23 13 12)
  subroutine check_single_slip_latent(stiff_c)

    real(rk), intent(in) :: stiff_c(6, 6)

    type(hardening_params_t) :: hard
    real(rk), parameter :: latent(2) = [1.4_rk, 1.0_rk]
    integer, parameter :: free(5) = [1, 2, 4, 5, 6]
    real(rk) :: axis(3), e1(3), e2(3), rot0(3, 3), rot(3, 3), sig(6), e(6), d(6), w(3)
    real(rk) :: rot_n(3, 3), sig_n(6), e_n(6), wp(3), sec(6, 6), h(6), tan(6, 6)
    real(rk) :: jff(5, 5), jff_tj(5, 5), step(5), leak, leak_first, leak_last, crss_a_latent
    real(rk), allocatable :: schmid(:), crss_prev(:), rate(:), c(:), crss_ref(:)
    integer :: n, a, b, k, j, it, l, status, ios

    n = crys%num_slip
    allocate (schmid(n), crss_prev(n), rate(n), c(n), crss_ref(n))
    if (allocated(crss_end)) deallocate (crss_end)
    allocate (crss_end(n))

    ! Orientation: crystal [123] along sample z; Schmid factors P^a : (axis axis) and system a.
    axis = [1.0_rk, 2.0_rk, 3.0_rk]/sqrt(14.0_rk)
    e1 = [2.0_rk, -1.0_rk, 0.0_rk]/sqrt(5.0_rk)
    e2 = [axis(2)*e1(3) - axis(3)*e1(2), axis(3)*e1(1) - axis(1)*e1(3), &
      & axis(1)*e1(2) - axis(2)*e1(1)]
    rot0(1, :) = e1
    rot0(2, :) = e2
    rot0(3, :) = axis
    schmid = abs(matmul(mandel_from_matrix(spread(axis, 2, 3)*spread(axis, 1, 3)), &
      & crys%schmid_sym_c))
    a = maxloc(schmid, 1)
    call assert_close(schmid(a), 0.4666_rk, 'single_slip_schmid_factor', abs_tol=1.0e-4_rk)

    do l = 1, 2
      ! The anisotropic law with H = 1 on the diagonal, latent(l) off it; n = 1.
      hard%law = hardening_anisotropic
      hard%h_0 = spread(200.0_rk, 1, n)
      hard%tau_s = spread(330.0_rk, 1, n)
      hard%n = spread(1.0_rk, 1, n)
      hard%tol_hardening = 1.0e-13_rk
      hard%max_iters = 50
      hard%interaction = spread(spread(latent(l), 1, n), 1, n)
      do j = 1, n
        hard%interaction(j, j) = 1.0_rk
      end do
      call cp_model_init(model, stiff_c, crys%schmid_sym_c, crys%schmid_skew_c, &
        & spread(0.02_rk, 1, n), spread(1.0_rk, 1, n), spread(100.0_rk, 1, n), &
        & tol_newton=1.0e-14_rk, max_iters_newton=100, max_iters_state=200, &
        & tol_state=1.0e-14_rk, hard=hard)

      ! Zero start state; first guess of {d}: isochoric uniaxial.
      rot = rot0
      sig = 0.0_rk
      e = 0.0_rk
      crss_prev = model%tau_0
      d = 1.0e-3_rk*[-0.5_rk, -0.5_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]
      w = 0.0_rk
      leak_first = -1.0_rk
      leak_last = -1.0_rk
      do k = 1, 300
        ! Uniaxial stress: quasi-Newton on the free components of {d} with [s_tan].
        do it = 1, 50
          call cp_work_alloc(work, model)
          call cp_point_solve(model, crss_prev, rot, sig, e, d, w, 0.1_rk, sig_n, e_n, wp, &
            & rot_n, crss_end, sec, h, tan, rate, status, work=work)
          if (maxval(abs(sig_n(free))) <= 1.0e-9_rk*abs(sig_n(3)) + 1.0e-9_rk) exit
          jff = tan(free, free)
          call normal_equations_solve(jff, -sig_n(free), step, jff_tj, ios)
          d(free) = d(free) + step
        end do
        if (l == 1) then
          call assert_true(status == cp_converged .and. maxval(abs(sig_n(free))) <= &
            & 1.0e-9_rk*abs(sig_n(3)) + 1.0e-9_rk, 'single_slip_uniaxial_stress')

          ! Backward-Euler recursion with the returned slip rates (n = 1, below saturation).
          c = 0.1_rk*hard%h_0*matmul(hard%interaction, abs(rate))/(hard%tau_s - model%tau_0)
          crss_ref = (crss_prev + c*hard%tau_s)/(1.0_rk + c)
          call assert_close(crss_end, crss_ref, 'single_slip_recursion', rel_tol=1.0e-9_rk)
        end if

        ! Slip on the other systems relative to system a, once system a slips (above 1e-4 1/s;
        ! about 2e-3 1/s in the flow; zero below the slip cutoff in the elastic increments).
        if (abs(rate(a)) > 1.0e-4_rk) then
          leak = 0.0_rk
          do b = 1, n
            if (b /= a) leak = max(leak, abs(rate(b))/abs(rate(a)))
          end do
          if (leak_first < 0.0_rk) leak_first = leak
          leak_last = leak
          if (l == 1) call assert_true(leak < 1.0e-5_rk, 'single_slip_one_system')
        end if

        ! Next increment.
        crss_prev = crss_end
        rot = rot_n
        sig = sig_n
        e = e_n
      end do

      write (*, '(a, f3.1, a, es9.2, a, es9.2, a, 2f9.3, a)') &
        & 'single slip [123], H off-diagonal ', latent(l), ': other/primary slip rate ', &
        & leak_first, ' -> ', leak_last, &
        & '; crss^a, crss^b ', crss_end(a), crss_end(merge(2, 1, a == 1)), ' MPa'
      if (l == 1) then
        ! Latent hardening suppresses the other systems; the inactive strengths are equal.
        call assert_true(leak_last < 0.2_rk*leak_first, 'single_slip_latent_suppresses')
        b = merge(2, 1, a == 1)
        do j = 1, n
          if (j /= a) call assert_close(crss_end(j), crss_end(b), 'single_slip_inactive_equal', &
            & rel_tol=1.0e-6_rk)
        end do
        ! The integrated single-slip relation (n = 1), to the time-discretization error.
        crss_a_latent = crss_end(a)
        write (*, '(a, es9.2)') &
          & 'single slip: integrated relation, relative error in tau_s - crss^b ', &
          & abs((330.0_rk - crss_end(b))/(230.0_rk*((330.0_rk - crss_end(a))/230.0_rk)**1.4_rk) &
          & - 1.0_rk)
        call assert_close(330.0_rk - crss_end(b), 230.0_rk*((330.0_rk - crss_end(a))/ &
          & 230.0_rk)**1.4_rk, 'single_slip_integrated_relation', rel_tol=3.0e-5_rk)
      else
        ! System a hardens by its own slip only (H_aa = 1): the same as with latent 1.4.
        call assert_close(crss_end(a), crss_a_latent, 'single_slip_primary_independent', &
          & rel_tol=1.0e-6_rk)
        ! No extra latent hardening: every system hardens alike.
        call assert_close(crss_end, spread(crss_end(a), 1, n), 'single_slip_h_ones_equal', &
          & rel_tol=2.0e-6_rk)
      end if
    end do

  end subroutine check_single_slip_latent

  !> Rotation angle of a rotation matrix, acos((trace - 1)/2) [rad].
  !>
  !> Arguments:
  !>   r  (in)  real(rk)(3, 3)  rotation
  pure function rotation_angle(r) result(angle)

    real(rk), intent(in) :: r(3, 3)
    real(rk) :: angle

    angle = acos(min(1.0_rk, max(-1.0_rk, 0.5_rk*(r(1, 1) + r(2, 2) + r(3, 3) - 1.0_rk))))

  end function rotation_angle

  !> det(I + A).
  !>
  !> Arguments:
  !>   a  (in)  real(rk)(3, 3)  matrix
  !> Locals:
  !>   m  real(rk)(3, 3)        I + A
  pure function det3(a) result(det)

    real(rk), intent(in) :: a(3, 3)
    real(rk) :: det

    real(rk) :: m(3, 3)

    m = a
    m(1, 1) = m(1, 1) + 1.0_rk
    m(2, 2) = m(2, 2) + 1.0_rk
    m(3, 3) = m(3, 3) + 1.0_rk
    det = m(1, 1)*(m(2, 2)*m(3, 3) - m(2, 3)*m(3, 2)) - m(1, 2)*(m(2, 1)*m(3, 3) - &
      & m(2, 3)*m(3, 1)) + m(1, 3)*(m(2, 1)*m(3, 2) - m(2, 2)*m(3, 1))

  end function det3

end program test_cp_model
