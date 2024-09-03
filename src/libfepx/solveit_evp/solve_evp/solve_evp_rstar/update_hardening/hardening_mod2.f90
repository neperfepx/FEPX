! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module hardening_mod2

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use crys_type_mod2
  use matrix_operations_mod
  use hardening_mod3
  use parallel_mod, only: par_quit
  use gather_scatter_mod
  use utils_mod

  implicit none

  public

contains

  subroutine hard_law(mesh, crys, acmslip, sliprate, func, dfunc, crss, crss_sat, shear, shrate, epseff, &
      & icode)

    ! Evaluate hardening rate or derivative of hardening rate

    !---------------------------------------------------------------------------

    ! Arguments:
    ! func: Array of computed hardening rates
    ! dfunc: Array of computed hardening rate derivatives
    ! crss: Array of critical resolved shear stress (hardness)
    ! crss_sat: Array of saturation crss
    ! shear: Shear values
    ! shrate: Shear rates
    ! epseff:
    ! icode: Hardening law

    type(mesh_type) :: mesh
    type(crys_type) :: crys (:)
    real(rk), intent(in) :: acmslip(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: sliprate(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(out) :: func(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(out) :: dfunc(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: crss(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: crss_sat(elt_sub:elt_sup)
    real(rk), intent(in) :: shrate(elt_sub:elt_sup)
    real(rk), intent(in) :: shear(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: epseff(elt_sub:elt_sup)
    integer, intent(in) :: icode

    ! Locals:
    ! mpk: Do any of these need to be defined here? They have been moved to the
    !   individual subroutines...

    real(rk) :: mysign(mesh%maxnumslip, elt_sub:elt_sup)
    integer :: i, islip
    integer :: my_phase(elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    ! rc: Initialize these to zero or else the terms that are not used in a dual
    !   phase material get assigned random values.

    func = 0.0d0
    dfunc = 0.0d0

    ! rc: In future updates this will need to be set up such that it makes a
    !   call for either a isotropic case or anisotropic hardening case

    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    mysign = 1.0d0

    do islip = 1, mesh%maxnumslip
      where ((crss_sat - crss(islip, :)) .le. 0.0d0)
        mysign(islip, :) = 0.0d0
      end where
    end do

    ! Options for which hardening model to use

    do i = 1, size(crys)
      if (crys(i)%cyclic) then
        call cyclic_hardening(mesh, crys(i), i, acmslip, func, dfunc, crss, crss_sat, shear, epseff, &
            & icode, mysign, my_phase)

      else if (.not. crys(i)%anisotropic) then
        call isotropic_hardening(mesh, crys(i), i, func, dfunc, crss, crss_sat, shrate, icode, &
            & mysign, my_phase)

      else if (crys(i)%anisotropic) then
        call anisotropic_hardening(mesh, crys(i), i, sliprate, func, dfunc, crss, crss_sat, icode, &
            & mysign, my_phase)
      else
        call par_quit('Error  :     > Invalid hardening model option input.')
      end if
    end do

  end subroutine hard_law

end module hardening_mod2
