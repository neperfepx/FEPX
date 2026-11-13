! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of slip_kinetics_mod: cutoff, resolved shear stress, power law, [M*] and [M*_tan].
!>
!> Source : FEPX (ROADMAP step 6.1).
!>
!> Checks (random inputs with a fixed seed; fcc, bcc with both families, hcp, bct):
!>   - cutoff: the values printed by v2-fixed's set_tmin in the C35 runs (m = 0.010, 0.013:
!>     5.7223676593502179e-2, 5.9161131848674935e-3) and its fallback exp(-50) for
!>     1/m - 1 <= 1 (m = 0.5, 0.9), to 1e-15 relative;
!>   - resolved shear stress (R5) against explicit 3x3 algebra, d . sigma . n;
!>   - power law (R6) against the explicit formula gammadot_0 (|tau|/crss)^(1/m) sgn(tau), with m
!>     and gammadot_0 differing between systems; slip rate = secant coefficient x rss; tangent
!>     coefficient = secant coefficient / m;
!>   - v2: the slip rate and tangent coefficient against v2's power_law and compliance formulas
!>     transcribed (gammadot_0 r exp(p log|r|), gammadot/(r crss m));
!>   - cutoff behaviour: zero at and just below the cutoff, nonzero just above; the cutoff is
!>     strict (">", as v2);
!>   - [M*] (R7): [M*]{tau} = sum_a gammadot^a {p^a} (the plastic deformation rate), symmetric,
!>     positive semidefinite (eigen check through random quadratic forms), [M*]{I} = 0
!>     (deviatoric);
!>   - [M*_tan]: equals the central finite-difference Jacobian of {d^p}(tau) = sum_a
!>     gammadot^a(tau) {p^a} (relative tolerance 1e-6 of its largest entry);
!>   - secant = m x tangent when every m is equal; with m = 0.02 and 0.04 alternating,
!>     [M*] - 0.02 [M*_tan] = sum_a (c^a - 0.02 c_tan^a) p p^T, nonzero when an m = 0.04
!>     system slips (CHANGES.md C11);
!>   - frame: [M*] built in the crystal frame and transformed to the sample frame equals the one
!>     built from sample-frame Schmid tensors and stress;
!>   - systems with no slip (coefficient 0) add nothing.
!>
!> Program locals:
!>   crys        type(crystal_t)   crystal under test
!>   t_min       real(rk)          a cutoff [-]
!>   rate, coef, coef_tan  real(rk)(1)  power-law outputs of one system
module slip_kinetics_test_mod

  use kinds_mod, only: rk
  use tensor_mod, only: mandel_from_matrix, matrix_from_mandel
  use rotation_mod, only: transform_crystal_to_sample_mandel, transform_crystal_to_sample_mat6
  use crystal_mod, only: crystal_t
  use slip_kinetics_mod, only: slip_rss_cutoff, resolved_shear_stress, slip_power_law, &
    & plastic_compliance
  use tensor_test_utils_mod, only: random_symmetric, random_rotation
  use assert_mod, only: assert_true, assert_close

  implicit none

  private

  public :: check_crystal

contains

  !> All per-crystal checks with random stresses and parameters.
  !>
  !> Arguments:
  !>   crys  (in)  type(crystal_t)  crystal (slip geometry)
  !>   name  (in)  character        label
  !> Locals:
  !>   n           integer                 number of systems
  !>   a, i, j     integer                 system, repetition, random vector
  !>   stress_c    real(rk)(6)             Kirchhoff stress, crystal frame [MPa]
  !>   rss, crss, m, gdot0, cutoff, rate, coef, coef_tan  real(rk)(n)  per-system quantities
  !>   expected    real(rk)                explicit slip rate [1/s]
  !>   ratio       real(rk)                tau/crss
  !>   mstar, mtan real(rk)(6,6)           [M*], [M*_tan], crystal frame
  !>   dp          real(rk)(6)             sum_a gammadot^a p^a
  !>   jac         real(rk)(6,6)           finite-difference Jacobian
  !>   x           real(rk)(6)             random vector
  !>   rot         real(rk)(3,3)           random rotation (crystal to sample)
  !>   schmid_s    real(rk)(6,n)           Schmid tensors, sample frame
  !>   sig3        real(rk)(3,3)           stress as a matrix
  subroutine check_crystal(crys, name)

    type(crystal_t), intent(in) :: crys
    character(len=*), intent(in) :: name

    integer :: n, a, i, j
    real(rk) :: stress_c(6), expected, ratio, mstar(6, 6), mtan(6, 6), dp(6), jac(6, 6), x(6)
    real(rk) :: rot(3, 3), sig3(3, 3)
    real(rk), allocatable :: rss(:), crss(:), m(:), gdot0(:), cutoff(:), rate(:), coef(:), &
      & coef_tan(:), schmid_s(:, :)

    n = crys%num_slip
    allocate (rss(n), crss(n), m(n), gdot0(n), cutoff(n), rate(n), coef(n), coef_tan(n))
    allocate (schmid_s(6, n))

    do i = 1, 10
      ! Random parameters: crss in [100, 300) MPa, m in [0.01, 0.06), gammadot_0 in [0.5, 2) 1/s;
      ! stress scaled so that |tau/crss| reaches about 1.
      call random_number(crss)
      crss = 100.0_rk + 200.0_rk*crss
      call random_number(m)
      m = 0.01_rk + 0.05_rk*m
      call random_number(gdot0)
      gdot0 = 0.5_rk + 1.5_rk*gdot0
      cutoff = slip_rss_cutoff(m)
      sig3 = 400.0_rk*random_symmetric()
      stress_c = mandel_from_matrix(sig3)

      ! Resolved shear stress against d . sigma . n.
      rss = resolved_shear_stress(crys%schmid_sym_c, stress_c)
      do a = 1, n
        call assert_close(rss(a), dot_product(crys%slip_dir_c(:, a), &
          & matmul(sig3, crys%slip_normal_c(:, a))), name // '_rss_3x3', abs_tol=1.0e-10_rk)
      end do

      ! Power law against the explicit formula and v2's transcribed formulas.
      call slip_power_law(rss, crss, m, gdot0, cutoff, rate, coef, coef_tan)
      do a = 1, n
        ratio = rss(a)/crss(a)
        if (abs(ratio) > cutoff(a)) then
          expected = gdot0(a)*(abs(ratio))**(1.0_rk/m(a))*sign(1.0_rk, ratio)
          call assert_close(rate(a), expected, name // '_power_law', &
            & rel_tol=1.0e-12_rk, abs_tol=1.0e-300_rk)
          call assert_close(rate(a), coef(a)*rss(a), name // '_rate_is_coef_rss', &
            & rel_tol=1.0e-15_rk, abs_tol=1.0e-300_rk)
          call assert_close(coef_tan(a), coef(a)/m(a), name // '_tangent_coef', &
            & rel_tol=1.0e-15_rk, abs_tol=1.0e-300_rk)
          ! v2: power_law gdot = gammadot_0 r exp((1/m - 1) log|r|); compliance gdot/(r crss m).
          expected = gdot0(a)*ratio*exp((1.0_rk/m(a) - 1.0_rk)*log(abs(ratio)))
          call assert_close(rate(a), expected, name // '_v2_power_law', &
            & rel_tol=1.0e-14_rk, abs_tol=1.0e-300_rk)
          call assert_close(coef_tan(a), expected/(ratio*crss(a)*m(a)), name // &
            & '_v2_compliance', rel_tol=1.0e-14_rk, abs_tol=1.0e-300_rk)
        else
          call assert_true(.not. (abs(rate(a)) > 0.0_rk) .and. .not. (abs(coef(a)) > 0.0_rk) &
            & .and. .not. (abs(coef_tan(a)) > 0.0_rk), name // '_below_cutoff_zero')
        end if
      end do

      ! [M*]{tau} = sum_a gammadot^a p^a; symmetric, positive semidefinite, deviatoric.
      mstar = plastic_compliance(crys%schmid_sym_c, coef)
      dp = matmul(crys%schmid_sym_c, rate)
      call assert_close(matmul(mstar, stress_c), dp, name // '_mstar_dp', &
        & rel_tol=1.0e-12_rk, abs_tol=1.0e-12_rk*maxval(abs(dp)))
      call assert_close(mstar, transpose(mstar), name // '_mstar_symmetric', &
        & abs_tol=1.0e-15_rk*maxval(abs(mstar)))
      do j = 1, 5
        call random_number(x)
        x = x - 0.5_rk
        call assert_true(dot_product(x, matmul(mstar, x)) >= -1.0e-14_rk*maxval(abs(mstar)), &
          & name // '_mstar_psd')
      end do
      call assert_close(matmul(mstar, [1.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]), &
        & spread(0.0_rk, 1, 6), name // '_mstar_deviatoric', &
        & abs_tol=1.0e-14_rk*maxval(abs(mstar)))

      ! [M*_tan] equals the finite-difference Jacobian of {d^p}(tau).
      mtan = plastic_compliance(crys%schmid_sym_c, coef_tan)
      jac = fd_jacobian(crys, stress_c, crss, m, gdot0, cutoff)
      call assert_close(mtan, jac, name // '_mtan_fd', abs_tol=1.0e-6_rk*maxval(abs(mtan)))

      ! Frame: crystal-frame [M*] transformed equals the sample-frame build.
      rot = random_rotation()
      do a = 1, n
        schmid_s(:, a) = transform_crystal_to_sample_mandel(rot, crys%schmid_sym_c(:, a))
      end do
      call assert_close(transform_crystal_to_sample_mat6(rot, mstar), &
        & plastic_compliance(schmid_s, coef), name // '_frame', &
        & abs_tol=1.0e-13_rk*maxval(abs(mstar)))
      call assert_close(resolved_shear_stress(schmid_s, &
        & transform_crystal_to_sample_mandel(rot, stress_c)), rss, name // '_rss_frame', &
        & abs_tol=1.0e-10_rk)

      ! Equal m: secant = m x tangent. Different m on active systems: not.
      m = 0.02_rk
      cutoff = slip_rss_cutoff(m)
      call slip_power_law(rss, crss, m, gdot0, cutoff, rate, coef, coef_tan)
      mstar = plastic_compliance(crys%schmid_sym_c, coef)
      mtan = plastic_compliance(crys%schmid_sym_c, coef_tan)
      call assert_close(mstar, 0.02_rk*mtan, name // '_secant_equal_m', &
        & abs_tol=1.0e-14_rk*maxval(abs(mstar)))
      do a = 1, n
        m(a) = merge(0.02_rk, 0.04_rk, mod(a, 2) == 0)
      end do
      cutoff = slip_rss_cutoff(m)
      call slip_power_law(rss, crss, m, gdot0, cutoff, rate, coef, coef_tan)
      mstar = plastic_compliance(crys%schmid_sym_c, coef)
      mtan = plastic_compliance(crys%schmid_sym_c, coef_tan)
      call assert_close(mstar - 0.02_rk*mtan, plastic_compliance(crys%schmid_sym_c, &
        & coef - 0.02_rk*coef_tan), name // '_secant_mixed_m_difference', &
        & abs_tol=1.0e-14_rk*maxval(abs(mstar)))
      if (any(abs(coef) > 0.0_rk .and. abs(m - 0.04_rk) < 1.0e-12_rk)) then
        call assert_true(maxval(abs(mstar - 0.02_rk*mtan)) > 0.0_rk, &
          & name // '_secant_differs_mixed_m')
      end if
    end do

    ! Systems with zero coefficients add nothing.
    coef = 0.0_rk
    call assert_close(plastic_compliance(crys%schmid_sym_c, coef), &
      & reshape(spread(0.0_rk, 1, 36), [6, 6]), name // '_zero_coef', abs_tol=0.0_rk)

  end subroutine check_crystal

  !> Central finite-difference Jacobian of {d^p}(tau) = sum_a gammadot^a(tau) {p^a} with
  !> respect to the Mandel stress, at fixed strengths.
  !>
  !> Arguments:
  !>   crys      (in)  type(crystal_t)  slip geometry
  !>   stress_c  (in)  real(rk)(6)      stress [MPa]
  !>   crss, m, gdot0, cutoff  (in)  real(rk)(:)  per-system parameters
  !>   jac       (out, result) real(rk)(6,6)  Jacobian [1/(MPa s)]
  !> Locals:
  !>   k         integer      perturbed component
  !>   h         real(rk)     step [MPa]
  !>   sp, sm    real(rk)(6)  perturbed stresses
  !>   rate, coef, coef_tan, rss  real(rk)(:)  power-law outputs
  !>   dp_p, dp_m  real(rk)(6)  {d^p} at the perturbed stresses
  function fd_jacobian(crys, stress_c, crss, m, gdot0, cutoff) result(jac)

    type(crystal_t), intent(in) :: crys
    real(rk), intent(in) :: stress_c(6)
    real(rk), intent(in) :: crss(:), m(:), gdot0(:), cutoff(:)
    real(rk) :: jac(6, 6)

    integer :: k
    real(rk) :: h, sp(6), sm(6), dp_p(6), dp_m(6)
    real(rk) :: rate(size(crss)), coef(size(crss)), coef_tan(size(crss)), rss(size(crss))

    ! Step 1e-6 of the stress magnitude: truncation error O(h^2), round-off O(eps/h).
    h = 1.0e-6_rk*norm2(stress_c)
    do k = 1, 6
      sp = stress_c
      sp(k) = sp(k) + h
      sm = stress_c
      sm(k) = sm(k) - h
      rss = resolved_shear_stress(crys%schmid_sym_c, sp)
      call slip_power_law(rss, crss, m, gdot0, cutoff, rate, coef, coef_tan)
      dp_p = matmul(crys%schmid_sym_c, rate)
      rss = resolved_shear_stress(crys%schmid_sym_c, sm)
      call slip_power_law(rss, crss, m, gdot0, cutoff, rate, coef, coef_tan)
      dp_m = matmul(crys%schmid_sym_c, rate)
      jac(:, k) = (dp_p - dp_m)/(2.0_rk*h)
    end do

  end function fd_jacobian

end module slip_kinetics_test_mod

program test_slip_kinetics

  use kinds_mod, only: rk
  use crystal_mod, only: crystal_t, crystal_init, crystal_fcc, crystal_bcc, crystal_hcp, &
    & crystal_bct
  use slip_kinetics_mod, only: slip_rss_cutoff, slip_power_law
  use slip_kinetics_test_mod, only: check_crystal
  use tensor_test_utils_mod, only: seed_random
  use assert_mod, only: assert_true, assert_close, test_summary

  implicit none

  type(crystal_t) :: crys
  real(rk) :: t_min, rate(1), coef(1), coef_tan(1)

  call seed_random(20261010)

  ! Cutoff: v2-fixed's values (C35 runs) and fallback.
  call assert_close(slip_rss_cutoff(0.010_rk), 5.7223676593502179e-2_rk, 'cutoff_m0.010', &
    & rel_tol=1.0e-15_rk)
  call assert_close(slip_rss_cutoff(0.013_rk), 5.9161131848674935e-3_rk, 'cutoff_m0.013', &
    & rel_tol=1.0e-15_rk)
  call assert_close(slip_rss_cutoff(0.5_rk), 1.9287498479639178e-22_rk, 'cutoff_m0.5', &
    & rel_tol=1.0e-15_rk)
  call assert_close(slip_rss_cutoff(0.9_rk), exp(-50.0_rk), 'cutoff_m0.9', rel_tol=1.0e-15_rk)

  ! Cutoff behaviour on one system (m = 0.02): strict ">". At the cutoff, crss = 1 MPa so that
  ! |tau/crss| equals t_min exactly.
  t_min = slip_rss_cutoff(0.02_rk)
  call slip_power_law([t_min], [1.0_rk], [0.02_rk], [1.0_rk], [t_min], rate, coef, coef_tan)
  call assert_true(.not. (abs(rate(1)) > 0.0_rk), 'cutoff_at_zero')
  call slip_power_law([-0.999_rk*t_min*100.0_rk], [100.0_rk], [0.02_rk], [1.0_rk], [t_min], &
    & rate, coef, coef_tan)
  call assert_true(.not. (abs(rate(1)) > 0.0_rk) .and. .not. (abs(coef_tan(1)) > 0.0_rk), &
    & 'cutoff_below_zero')
  call slip_power_law([-1.001_rk*t_min*100.0_rk], [100.0_rk], [0.02_rk], [1.0_rk], [t_min], &
    & rate, coef, coef_tan)
  call assert_true(rate(1) < 0.0_rk .and. coef(1) > 0.0_rk .and. coef_tan(1) > 0.0_rk, &
    & 'cutoff_above_slips_negative')

  ! Per-crystal checks.
  call crystal_init(crys, crystal_fcc)
  call check_crystal(crys, 'fcc')
  call crystal_init(crys, crystal_bcc, families=[2, 1])
  call check_crystal(crys, 'bcc24')
  call crystal_init(crys, crystal_hcp, c_over_a=1.587_rk)
  call check_crystal(crys, 'hcp')
  call crystal_init(crys, crystal_bct, c_over_a=0.5456_rk)
  call check_crystal(crys, 'bct')

  call test_summary()

end program test_slip_kinetics
