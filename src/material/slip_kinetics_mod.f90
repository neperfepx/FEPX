! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Slip kinetics of crystal plasticity: resolved shear stresses (R5), the power law (R6) with a
!> rate sensitivity m and a reference slip rate gammadot_0 per slip system, the
!> resolved-shear-stress cutoff that keeps the power law within range, and the viscoplastic
!> compliances [M*] (secant, R7) and [M*_tan] (tangent, §4). Pure per-point kernels.
!>
!> Theory : recoupled_formulation.md §2 (R5)-(R7) and §4 ([M*_tan]); per-system m and
!>          gammadot_0: ROADMAP D48.
!> Source : port of v2-fixed:
!>          - solve_evp_vpstress_mod3.f90::power_law (slip rate, zero at or below the cutoff of
!>            the normalized resolved shear stress) and ::compliance (tangent coefficient
!>            gammadot/(m tau), zero at or below the cutoff);
!>          - crys_type_mod2.f90::set_tmin (the cutoff);
!>          - solveit_evp_mod4.f90::material_matrix_evp_aniso (resolved shear stress normalized
!>            by the strength, tangent compliance sum_a comp_a p_a p_a^T, secant compliance
!>            sum_a m_a comp_a p_a p_a^T as fixed in CHANGES.md C11).
!>          Deviations: m and gammadot_0 per slip system (v2: gammadot_0 per phase, m per phase
!>          or per family), hence the cutoff per slip system, from that system's m (v2: one per
!>          phase, from the phase m, which is unset with per-family m, CHANGES.md C35); 6-vector
!>          Mandel Schmid tensors (v2: deviatoric 5-vectors); the secant coefficient
!>          gammadot/tau is computed directly, as recoupled_formulation.md §2 asks (v2 divides
!>          the slip rate by the resolved shear stress).
!>
!> Notation: tau^a resolved shear stress, crss^a strength, r^a = tau^a/crss^a, m^a, gammadot_0^a,
!> p^a the Mandel Schmid tensor of system a, all in the same frame (the crystal frame here,
!> suffix _c; [M*] is transformed to the sample frame with
!> rotation_mod::transform_crystal_to_sample_mat6 by the caller).
!>   rss cutoff (set_tmin): t_min = 10^(2 - 321/(1/m - 1)) if 1/m - 1 > 1, else exp(-50); it
!>   bounds |r|^(1/m - 1) below by about 1e-321 * 100^(1/m - 1), inside the range of real64;
!>   slip rate (R6):      gammadot^a = gammadot_0^a |r^a|^(1/m^a) sgn(tau^a)  if |r^a| > t_min^a,
!>                        0 otherwise;
!>   secant coefficient:  c^a = gammadot^a/tau^a = (gammadot_0^a/crss^a) |r^a|^(1/m^a - 1)  (same
!>                        condition), so that gammadot^a = c^a tau^a;
!>   tangent coefficient: c_tan^a = d gammadot^a/d tau^a = c^a/m^a (same condition);
!>   [M*]     = sum_a c^a {p^a}{p^a}^T                (R7, the secant: {d^p} = [M*]{tau});
!>   [M*_tan] = sum_a c_tan^a {p^a}{p^a}^T            (the derivative of {d^p} with respect to
!>                                                     {tau} at fixed strengths, §4).
!> The kinetic law returns both coefficients, so that [M*] and [M*_tan] are formed without
!> knowing the law (the relation c_tan = c/m holds for the power law only). Another law (a
!> thermally activated one is an option for Stage 7, ROADMAP Q20) would return its own pair.
!>
!> Public routines:
!>   slip_rss_cutoff        cutoff t_min of the normalized resolved shear stress, from m
!>   resolved_shear_stress  tau^a = {p^a}^T {tau} of every system (R5)
!>   slip_power_law         slip rates, secant and tangent coefficients of every system (R6)
!>   plastic_compliance     sum_a coef^a {p^a}{p^a}^T: [M*] with the secant coefficients,
!>                          [M*_tan] with the tangent coefficients (R7, §4)
module slip_kinetics_mod

  use kinds_mod, only: rk

  implicit none

  private

  public :: slip_rss_cutoff, resolved_shear_stress, slip_power_law, plastic_compliance

contains

  !> Cutoff of the normalized resolved shear stress |tau/g| below which a system does not slip
  !> (prevents underflow and overflow in the power law).
  !>
  !> Source : port of v2-fixed crys_type_mod2.f90::set_tmin (same formula and constants).
  !>
  !> Arguments:
  !>   m           (in)  real(rk)  rate sensitivity (> 0) [-]
  !>   rss_cutoff  (out, result) real(rk)  cutoff t_min of |tau/g| [-]
  !> Locals:
  !>   expo        real(rk)  exponent 1/m - 1 of the secant coefficient (v2's xnn)
  elemental function slip_rss_cutoff(m) result(rss_cutoff)

    real(rk), intent(in) :: m
    real(rk) :: rss_cutoff

    real(rk) :: expo

    ! v2: t_min = 10^(-321/xnn + 2) for xnn > 1, else exp(-50).
    expo = 1.0_rk/m - 1.0_rk
    if (expo > 1.0_rk) then
      rss_cutoff = exp((-321.0_rk/expo + 2.0_rk)*log(10.0_rk))
    else
      rss_cutoff = exp(-50.0_rk)
    end if

  end function slip_rss_cutoff

  !> Resolved shear stress of every slip system, tau^a = {p^a}^T {tau} (R5).
  !>
  !> Theory : recoupled_formulation.md (R5).
  !>
  !> Arguments:
  !>   schmid_sym_c  (in)  real(rk)(6, num_slip)  Mandel Schmid tensors p^a, crystal frame
  !>   stress_c      (in)  real(rk)(6)            Kirchhoff stress, Mandel, crystal frame [MPa]
  !>   rss           (out, result) real(rk)(num_slip)  resolved shear stresses [MPa]
  pure function resolved_shear_stress(schmid_sym_c, stress_c) result(rss)

    real(rk), intent(in) :: schmid_sym_c(:, :)
    real(rk), intent(in) :: stress_c(6)
    real(rk) :: rss(size(schmid_sym_c, 2))

    ! One dot product per system (the Mandel product is the tensor contraction).
    rss = matmul(stress_c, schmid_sym_c)

  end function resolved_shear_stress

  !> Slip rates (R6), secant coefficients gammadot^a/tau^a and tangent coefficients
  !> d gammadot^a/d tau^a = (gammadot^a/tau^a)/m^a of every slip system; all zero where
  !> |tau^a/crss^a| is at or below the system's cutoff.
  !>
  !> Theory : recoupled_formulation.md (R6) and the implementation note under (R7).
  !> Source : port of v2-fixed solve_evp_vpstress_mod3.f90::power_law and ::compliance (module
  !>          header; the strict ">" of both kept).
  !>
  !> Arguments:
  !>   rss          (in)  real(rk)(num_slip)  resolved shear stresses tau^a [MPa]
  !>   crss         (in)  real(rk)(num_slip)  strengths (> 0) [MPa]
  !>   m            (in)  real(rk)(num_slip)  rate sensitivities m^a (> 0) [-]
  !>   gammadot_0   (in)  real(rk)(num_slip)  reference slip rates gammadot_0^a (> 0) [1/s]
  !>   rss_cutoff   (in)  real(rk)(num_slip)  cutoffs t_min^a of |tau^a/crss^a|
  !>                                          (slip_rss_cutoff) [-]
  !>   sliprate     (out) real(rk)(num_slip)  slip rates gammadot^a [1/s]
  !>   coef_secant  (out) real(rk)(num_slip)  gammadot^a/tau^a [1/(MPa s)]
  !>   coef_tangent (out) real(rk)(num_slip)  d gammadot^a/d tau^a [1/(MPa s)]
  !> Locals:
  !>   a            integer   slip system
  !>   ratio_abs    real(rk)  |tau^a/crss^a| [-]
  pure subroutine slip_power_law(rss, crss, m, gammadot_0, rss_cutoff, sliprate, coef_secant, &
    & coef_tangent)

    real(rk), intent(in) :: rss(:)
    real(rk), intent(in) :: crss(:)
    real(rk), intent(in) :: m(:)
    real(rk), intent(in) :: gammadot_0(:)
    real(rk), intent(in) :: rss_cutoff(:)
    real(rk), intent(out) :: sliprate(:)
    real(rk), intent(out) :: coef_secant(:)
    real(rk), intent(out) :: coef_tangent(:)

    integer :: a
    real(rk) :: ratio_abs

    do a = 1, size(rss)
      ! Normalized resolved shear stress.
      ratio_abs = abs(rss(a)/crss(a))

      ! At or below the cutoff: no slip (v2: power_law and compliance both test "> t_min").
      if (.not. (ratio_abs > rss_cutoff(a))) then
        coef_secant(a) = 0.0_rk
        coef_tangent(a) = 0.0_rk
        sliprate(a) = 0.0_rk
        cycle
      end if

      ! c^a = (gammadot_0^a/crss^a) |r^a|^(1/m^a - 1), evaluated as exp(p log|r|) like v2's
      ! power_law; then gammadot^a = c^a tau^a, and the tangent c^a/m^a (v2's compliance).
      coef_secant(a) = gammadot_0(a)/crss(a)*exp((1.0_rk/m(a) - 1.0_rk)*log(ratio_abs))
      coef_tangent(a) = coef_secant(a)/m(a)
      sliprate(a) = coef_secant(a)*rss(a)
    end do

  end subroutine slip_power_law

  !> Viscoplastic compliance sum_a coef^a {p^a}{p^a}^T: the secant [M*] (R7) with the secant
  !> coefficients gammadot^a/tau^a, the tangent [M*_tan] (§4) with the tangent coefficients
  !> d gammadot^a/d tau^a (slip_power_law).
  !>
  !> Theory : recoupled_formulation.md (R7) and §4.
  !> Source : port of the compliance accumulation of v2-fixed
  !>          solveit_evp_mod4.f90::material_matrix_evp_aniso (tan_stif_vp, sec_stif_vp_aniso;
  !>          module header).
  !>
  !> Arguments:
  !>   schmid_sym_c  (in)  real(rk)(6, num_slip)  Mandel Schmid tensors p^a, crystal frame
  !>   coef          (in)  real(rk)(num_slip)     coefficient of each system [1/(MPa s)]
  !>   compliance_c  (out, result) real(rk)(6, 6)  the compliance, crystal frame [1/(MPa s)]
  !> Locals:
  !>   a, i, j       integer  slip system, row, column
  pure function plastic_compliance(schmid_sym_c, coef) result(compliance_c)

    real(rk), intent(in) :: schmid_sym_c(:, :)
    real(rk), intent(in) :: coef(:)
    real(rk) :: compliance_c(6, 6)

    integer :: a, i, j

    compliance_c = 0.0_rk
    do a = 1, size(coef)
      ! Systems that do not slip add nothing.
      if (.not. (abs(coef(a)) > 0.0_rk)) cycle

      ! coef^a {p^a}{p^a}^T, column by column.
      do j = 1, 6
        do i = 1, 6
          compliance_c(i, j) = compliance_c(i, j) + coef(a)*schmid_sym_c(i, a)*schmid_sym_c(j, a)
        end do
      end do
    end do

  end function plastic_compliance

end module slip_kinetics_mod
