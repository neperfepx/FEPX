! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module stiffness_mod

! Module for elemental stiffness matrix for the viscoplastic problem

! Contains subroutines:
! elt_stif_vp: Form elemental stiffness matrix for viscoplastic problem
! add_to_stiffness: Add component to stiffness matrix

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use matrix_operations_mod
  use quadrature_mod
  use shape_3d_mod

  implicit none

  public

contains

  subroutine add_to_stiffness(gstiff, sij, i, j)

    ! Add component to stiffness matrix.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! gstiff:
    ! sij:
    ! i:
    ! j:

    real(rk), intent(inout) :: gstiff(kdim, kdim, elt_sub:elt_sup)
    real(rk), intent(in) :: sij(elt_sub:elt_sup)
    integer :: i
    integer :: j

    !---------------------------------------------------------------------------

    if (i .ge. j) then
      if (i .eq. j) then
        gstiff(i, j, :) = gstiff(i, j, :) + sij

      else
        gstiff(i, j, :) = gstiff(i, j, :) + sij
        gstiff(j, i, :) = gstiff(j, i, :) + sij
      end if
    end if

  end subroutine add_to_stiffness

end module stiffness_mod
