! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of v2's increment-1 viscoplastic stress estimate
!> (crystal_plasticity_model_mod::cp_vertex_estimate, initialization isotropic_viscoplastic,
!> ROADMAP steps 6.6, 6.7, D78), started from the yield-surface vertices computed from the slip
!> systems (yield_vertices_mod::yield_vertices_compute, step 6.7, D57; stresses).
!>
!> Source : FEPX (ROADMAP step 6.6, D54 (3); step 6.7, D57; v2's tables and the comparison of
!>          the estimates from them removed in step 9.4).
!>
!> Checks:
!>   - estimate (fcc and hcp, Cu and Ti constants, m 0.01-0.03, strength crss, random rotations and
!>     deformation rates, tight Newton tolerance 1e-12): it converges (status 0) and solves the
!>     rigid-viscoplastic equations sum gammadot^a(tau / d_eq^0.4) p^a = d'/d_eq with
!>     gammadot^a = gammadot_0 |t^a|^(1/0.4 - 1) t^a, t^a = p^a : tau / (d_eq^0.4 crss) (explicit
!>     loop over the systems; every deviatoric component, so v2's three-component residual of
!>     CHANGES.md C43 fails it), to 1e-10 relative; it is deviatoric; it scales as
!>     estimate(lambda d) = lambda^0.4 estimate(d) and is odd, estimate(-d) = -estimate(d)
!>     (1e-12 relative); with the production tolerance point_tol = 1e-4, the residual of every
!>     deviatoric component is at most 1e-4 (C43: a three-component stop leaves the others
!>     above it); a zero rate gives zero; without vertices, status 2;
!>   - the point solve started from 0.5 times the estimate (tau_start_c) reaches the stress of
!>     the default start, to the state tolerance (tight tolerances: 1e-8 relative).
!>
!> Locals:
!>   work   type(cp_work_t)  point-solve work arrays (cp_work_alloc before each call; D100)
!>   crys   type(crystal_t)   slip geometry
!>   model  type(cp_model_t)  model
program test_vertex_estimate

  use kinds_mod, only: rk
  use tensor_mod, only: mandel_from_matrix, mandel_trace, mandel_norm, mandel_dev, equiv_strain
  use rotation_mod, only: transform_sample_to_crystal_mandel
  use crystal_mod, only: crystal_t, crystal_init, crystal_fcc, crystal_hcp, stiffness_cubic, &
    & stiffness_hexagonal
  use yield_vertices_mod, only: yield_vertices_compute
  use crystal_plasticity_model_mod, only: cp_model_t, cp_model_init, cp_vertex_estimate, &
    & cp_point_solve
  use crystal_plasticity_model_mod, only: cp_work_t, cp_work_alloc
  use tensor_test_utils_mod, only: seed_random, random_rotation, random_symmetric
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary

  implicit none

  type(cp_work_t) :: work
  type(crystal_t) :: crys
  type(cp_model_t) :: model

  call seed_random(20261011)

  call crystal_init(crys, crystal_fcc)
  call make_model(stiffness_cubic(168.4e3_rk, 121.4e3_rk, 75.4e3_rk), 100.0_rk)
  call check_estimate('fcc')
  call check_start('fcc')

  call crystal_init(crys, crystal_hcp, c_over_a=1.587_rk)
  call make_model(stiffness_hexagonal(162.4e3_rk, 92.0e3_rk, 69.0e3_rk, 180.7e3_rk, &
    & 46.7e3_rk), 150.0_rk)
  call check_estimate('hcp')
  call check_start('hcp')

  call test_summary()

contains

  !> Model of the current crystal with the vertices computed from its slip systems: m 0.01-0.03,
  !> gammadot_0 1, uniform strength crss, tight tolerances.
  !>
  !> Arguments:
  !>   stiff_c  (in)  real(rk)(6, 6)  stiffness [MPa]
  !>   crss     (in)  real(rk)        strength [MPa]
  !> Locals:
  !>   n, a       integer         systems, system
  !>   m          real(rk)(:)     rate sensitivities
  !>   vertices   real(rk)(:, :)  computed vertices [MPa]
  !>   num_indep  integer         independent combinations
  subroutine make_model(stiff_c, crss)

    real(rk), intent(in) :: stiff_c(6, 6)
    real(rk), intent(in) :: crss

    integer :: n, a, num_indep
    real(rk), allocatable :: m(:), vertices(:, :)

    ! Rate sensitivities varying over the systems, and the vertices at the uniform strength.
    n = crys%num_slip
    allocate (m(n))
    do a = 1, n
      m(a) = 0.01_rk + 0.01_rk*mod(a, 3)
    end do
    call yield_vertices_compute(crys%schmid_sym_c, spread(crss, 1, n), vertices, num_indep)
    call cp_model_init(model, stiff_c, crys%schmid_sym_c, crys%schmid_skew_c, m, &
      & spread(1.0_rk, 1, n), spread(crss, 1, n), tol_newton=1.0e-12_rk, &
      & max_iters_newton=200, max_iters_state=200, tol_state=1.0e-14_rk, &
      & vertices_c=vertices)

  end subroutine make_model

  !> The estimate solves its equations, scales and is odd (header).
  !>
  !> Arguments:
  !>   name  (in)  character  label
  !> Locals:
  !>   k            integer          sample
  !>   status, s2   integer          estimate statuses
  !>   d_c          real(rk)(6)      random deformation rate, crystal frame [1/s]
  !>   tau, tau2    real(rk)(6)      estimates [MPa]
  !>   crss         real(rk)(:)      strengths [MPa]
  !>   a            integer          index of the zero rate
  !>   no_table     type(cp_model_t) model without vertices
  !>   loose        type(cp_model_t) model with point_tol 1e-4
  !>   worst        real(rk)         largest full residual with point_tol 1e-4 [1/s]
  subroutine check_estimate(name)

    character(len=*), intent(in) :: name

    integer :: k, status, s2, a
    real(rk) :: d_c(6), tau(6), tau2(6), worst
    real(rk), allocatable :: crss(:)
    type(cp_model_t) :: no_table, loose

    crss = model%tau_0
    loose = model
    loose%tol_newton = 1.0e-4_rk
    worst = 0.0_rk
    do k = 1, 10
      ! A random deformation rate in the crystal frame (with a volumetric part).
      d_c = transform_sample_to_crystal_mandel(random_rotation(), &
        & 1.0e-3_rk*mandel_from_matrix(random_symmetric()))
      call cp_work_alloc(work, model)
      call cp_vertex_estimate(model, crss, d_c, tau, status, work)
      call assert_equal(status, 0, name // ' estimate converged')
      call assert_true(abs(mandel_trace(tau)) < 1.0e-10_rk*mandel_norm(tau), &
        & name // ' estimate deviatoric')

      ! The rigid-viscoplastic equations at tau / d_eq^0.4, every deviatoric component.
      call assert_true(full_residual(tau, d_c, crss) <= 1.0e-10_rk, &
        & name // ' estimate solves its equations')

      ! Production tolerance: the full residual at most point_tol.
      call cp_work_alloc(work, loose)
      call cp_vertex_estimate(loose, crss, d_c, tau2, s2, work)
      worst = max(worst, full_residual(tau2, d_c, crss))

      ! Scaling lambda^0.4 and odd symmetry.
      call cp_work_alloc(work, model)
      call cp_vertex_estimate(model, crss, 7.0_rk*d_c, tau2, s2, work)
      call assert_close(tau2, 7.0_rk**0.4_rk*tau, name // ' estimate scaling', &
        & abs_tol=1.0e-12_rk*mandel_norm(tau))
      call cp_work_alloc(work, model)
      call cp_vertex_estimate(model, crss, -d_c, tau2, s2, work)
      call assert_close(tau2, -tau, name // ' estimate odd', abs_tol=1.0e-12_rk*mandel_norm(tau))
    end do

    call assert_true(worst <= 1.0e-4_rk, name // ' full residual within point_tol 1e-4')

    ! Zero rate: zero estimate; no table: status 2.
    call cp_work_alloc(work, model)
    call cp_vertex_estimate(model, crss, [(0.0_rk, a = 1, 6)], tau, status, work)
    call assert_true(all(abs(tau) <= 0.0_rk) .and. status == 0, name // ' zero rate')
    no_table = model
    deallocate (no_table%vertices_c)
    allocate (no_table%vertices_c(6, 0))
    call cp_work_alloc(work, no_table)
    call cp_vertex_estimate(no_table, crss, d_c, tau, status, work)
    call assert_equal(status, 2, name // ' no vertices')

  end subroutine check_estimate

  !> Norm of the residual of the estimate's equations, every deviatoric component (explicit
  !> loop over the systems): sum gammadot^a p^a - d'/d_eq at tau / d_eq^0.4.
  !>
  !> Arguments:
  !>   tau   (in)  real(rk)(6)  estimate [MPa]
  !>   d_c   (in)  real(rk)(6)  deformation rate, crystal frame [1/s]
  !>   crss  (in)  real(rk)(:)  strengths [MPa]
  !>   norm  (out, result) real(rk)  residual norm [1/s]
  !> Locals:
  !>   deq, t, rate  real(rk)     equivalent rate, t^a, gammadot^a
  !>   resid         real(rk)(6)  residual
  !>   a             integer      system
  function full_residual(tau, d_c, crss) result(norm)

    real(rk), intent(in) :: tau(6)
    real(rk), intent(in) :: d_c(6)
    real(rk), intent(in) :: crss(:)
    real(rk) :: norm

    real(rk) :: deq, t, rate, resid(6)
    integer :: a

    deq = equiv_strain(d_c)
    resid = -mandel_dev(d_c)/deq
    do a = 1, crys%num_slip
      t = dot_product(crys%schmid_sym_c(:, a), tau)/(deq**0.4_rk*crss(1))
      rate = 0.0_rk
      if (abs(t) > model%rss_cutoff(a)) rate = model%gammadot_0(a)*abs(t)**(1.0_rk/0.4_rk - &
        & 1.0_rk)*t
      resid = resid + rate*crys%schmid_sym_c(:, a)
    end do
    norm = mandel_norm(resid)

  end function full_residual

  !> The point solve from 0.5 times the estimate reaches the stress of the default start
  !> (header).
  !>
  !> Arguments:
  !>   name  (in)  character  label
  !> Locals:
  !>   rot0         real(rk)(3, 3)   start orientation
  !>   d, w         real(rk)(6), (3) deformation rate [1/s], spin [1/s]
  !>   zero6        real(rk)(6)      zero stress and elastic strain
  !>   tau_est      real(rk)(6)      estimate (crystal frame) [MPa]
  !>   sig1, sig2   real(rk)(6)      final stresses (default and estimate starts) [MPa]
  !>   e, h         real(rk)(6)      elastic strain, {h} (unused)
  !>   wp           real(rk)(3)      w^p (unused)
  !>   rot          real(rk)(3, 3)   final orientation (unused)
  !>   crss, rate   real(rk)(:)      strengths and slip rates (unused)
  !>   sec, tan     real(rk)(6, 6)   point matrices (unused)
  !>   status, s1, s2  integer       statuses
  subroutine check_start(name)

    character(len=*), intent(in) :: name

    real(rk) :: rot0(3, 3), d(6), w(3), zero6(6), tau_est(6), sig1(6), sig2(6), e(6), h(6)
    real(rk) :: wp(3), rot(3, 3), sec(6, 6), tan(6, 6)
    real(rk), allocatable :: crss(:), rate(:)
    integer :: status, s1, s2

    rot0 = random_rotation()
    d = 1.0e-3_rk*mandel_dev(mandel_from_matrix(random_symmetric()))
    call random_number(w)
    w = 1.0e-4_rk*(w - 0.5_rk)
    zero6 = 0.0_rk
    allocate (crss(crys%num_slip), rate(crys%num_slip))
    call cp_work_alloc(work, model)
    call cp_vertex_estimate(model, model%tau_0, transform_sample_to_crystal_mandel(rot0, d), &
      & tau_est, status, work)
    call cp_work_alloc(work, model)
    call cp_point_solve(model, model%tau_0, rot0, zero6, zero6, d, w, 0.5_rk, sig1, e, wp, &
      & rot, crss, sec, h, tan, rate, s1, work=work)
    call cp_work_alloc(work, model)
    call cp_point_solve(model, model%tau_0, rot0, zero6, zero6, d, w, 0.5_rk, sig2, e, wp, &
      & rot, crss, sec, h, tan, rate, s2, tau_start_c=0.5_rk*tau_est, work=work)
    call assert_true(s1 == 0 .and. s2 == 0, name // ' point solves converged')
    call assert_close(sig2, sig1, name // ' same stress from the estimate start', &
      & abs_tol=1.0e-8_rk*mandel_norm(sig1))

  end subroutine check_start

end program test_vertex_estimate
