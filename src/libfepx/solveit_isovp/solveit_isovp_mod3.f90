! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module solveit_isovp_mod3

! Material matrix for the viscoplastic solution

! Contains subroutines:
! material_matrix_vp: Material matrix for the viscoplastic solution
! isotropic: Return viscosity for the isotropic viscoplastic solution

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use solveit_isovp_mod4
  use kinematics_mod
  use solve_evp_vpstress_mod
  use parallel_mod
  use crys_type_mod2

  implicit none

  public

contains

  !> Material matrix for the viscoplastic solution
  subroutine material_matrix_vp(itype, stif, dndx, dndy, dndz, gvel, scale, &
      & epseff)

    character(len = *) :: itype
    real(rk) :: stif(5, 5, elt_sub:elt_sup)
    real(rk) :: dndx(ndim, elt_sub:elt_sup)
    real(rk) :: dndy(ndim, elt_sub:elt_sup)
    real(rk) :: dndz(ndim, elt_sub:elt_sup)
    real(rk) :: gvel(kdim, elt_sub:elt_sup)
    real(rk) :: scale(elt_sub:elt_sup)
    real(rk) :: epseff(elt_sub:elt_sup)

    integer :: i
    integer :: m_elt
    real(rk) :: d(3, 3, elt_sub:elt_sup)
    real(rk) :: cmu(elt_sub:elt_sup)
    real(rk) :: temp_k(elt_sub:elt_sup)
    real(rk) :: dsdt_iso(elt_sub:elt_sup)
    real(rk) :: state_iso(elt_sub:elt_sup)
    real(rk) :: cconst(5)

    data cconst/1.0d0, 3.0d0, 4.0d0, 4.0d0, 4.0d0/

    !---------------------------------------------------------------------------

    m_elt = elt_sup - elt_sub + 1

    call defrate(gvel, dndx, dndy, dndz, d)

    call defrate_eq(d, epseff)

    ! hr-tm
    ! we aren't going through this portion to update for two-phase compatibility
    if (itype .eq. 'anisotropic_vp') then
      call par_quit('Error  :     > anisotropic_vp is no longer implemented.')

    else if (itype .eq. 'isotropic_vp') then
      ! dbg: Note the hardwired temperature and state variable.

      temp_k = 273.0d0 + 400.0d0
      state_iso = 20.0d6

      call isotropic(epseff, temp_k, state_iso, cmu, dsdt_iso)

      stif = 0.0d0
      stif(1, 1, :) = cconst(1)*cmu
      stif(2, 2, :) = cconst(2)*cmu
      stif(3, 3, :) = cconst(3)*cmu
      stif(4, 4, :) = cconst(4)*cmu
      stif(5, 5, :) = cconst(5)*cmu
    end if

    scale = 0.0d0

    do i = 1, 5
      scale = scale + stif(i, i, :)/cconst(i)
    end do

    scale = scale/5.0d0

  end subroutine material_matrix_vp

end module solveit_isovp_mod3
