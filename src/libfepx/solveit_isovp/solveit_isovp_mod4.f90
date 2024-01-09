! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module solveit_isovp_mod4

! Material matrix for the viscoplastic solution

! Contains subroutines:
! material_matrix_vp: Material matrix for the viscoplastic solution
! isotropic: Return viscosity for the isotropic viscoplastic solution

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use kinematics_mod
  use solve_evp_vpstress_mod
  use parallel_mod
  use crys_type_mod2

  implicit none

  public

contains

  subroutine isotropic(dii, t, s_pa, visc, dsdt)

    ! Return viscosity for isotropic problem. Apparently, this routine was
    !   disabled and was replaced by a linear viscous problem with unit
    !   viscosity. mpk: Is this still true?

    ! Computes the viscosity for an element based on the effective strain rate.
    !   Uses the fit of Dawson to Hart's model for pure Aluminum.

    ! Modifications:
    ! The effective viscosity computation has been changed to:
    !   visc = (2nd invariant of sigma) / 3*(2nd invariant of strain rate)
    !   so as to be compatible with definition of strain rate invariant
    !   by ajb (9/25/87). Units are MPa, kJoule/mole-k, etc.

    !----------------------------------------------------------------------

    ! Arguments:
    ! dii: The effective strain rate for the element
    ! t: Element temperature (f)
    ! s_pa: The value of the state variable (Pascals)
    ! visc: The resulting viscosity
    ! dsdt: The derivative of the state variable

    real(rk) :: dii(elt_sub:elt_sup)
    real(rk) :: t(elt_sub:elt_sup)
    real(rk) :: s_pa(elt_sub:elt_sup)
    real(rk) :: visc(elt_sub:elt_sup)
    real(rk) :: dsdt(elt_sub:elt_sup)

    ! Locals:

    ! Constants:

    !real(rk) :: a0
    !parameter ( a0     = 9.64d52     )
    real(rk), parameter :: log_a0 = 122.0d0
    real(rk), parameter :: c0 = 6.19d-9
    real(rk), parameter :: qpr = 1.45d4
    real(rk), parameter :: ge = 24.20d3
    real(rk), parameter :: m = 7.80d0
    real(rk), parameter :: f0 = 2.12d19
    real(rk), parameter :: qr = qpr
    real(rk), parameter :: smallm = 5.0d0
    real(rk), parameter :: lambda = 0.14d0
    real(rk), parameter :: mp = 3.5d0
    real(rk), parameter :: n = 6.0d0

    real(rk) :: a(elt_sub:elt_sup)
    real(rk) :: s(elt_sub:elt_sup)
    real(rk) :: dstar(elt_sub:elt_sup)
    real(rk) :: sigma(elt_sub:elt_sup)
    real(rk) :: sigmap(elt_sub:elt_sup)
    real(rk) :: sigmav(elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    ! Place the state variable in Pascals

    visc = 1.0d0

    return

    s = s_pa/1.0d6

    ! Compute flow stress in viscous (friction) element
    !   a = a0 * exp(-qpr/t)
    a = log_a0 + (-qpr/t)
    a = exp(a)

    sigmav = ge*(dii/a)**(1.0/m)

    ! Compute flow stress in the plastic element

    dstar = f0*((s/ge)**smallm)*exp(-qr/t)
    sigmap = s*exp(-(dstar/dii)**lambda)

    ! Total flow stress

    sigma = sigmav + sigmap

    ! Compute the derivative of the state vaiable

    dsdt = c0*s*dii*((ge/s)**mp)*((sigmap/s)**n)

    ! Convert to Pascals

    !sigma = s * dii**0.10
    sigma = sigma*1.0d6

    ! Compute viscosity

    visc = sigma/(3.0d0*dii)

    print *, 'sigma', minval(sigma), maxval(sigma)
    print *, 'visc', minval(visc), maxval(visc)

    return

  end subroutine isotropic

end module solveit_isovp_mod4
