! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module solveit_evp_mod3

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use solveit_evp_mod4
  use kinematics_mod
  use solve_evp_mod
  use units_mod

  implicit none

  public

contains

  !> Material matrix for the evp solution
  subroutine material_matrix_evp(mesh, crys, exec, &
                              & results_prev, results, c, &
                              & c_tan, eforce, dndx, dndy, dndz, evel, &
                              & incr, dtime, det, ecoos)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys(:)
    type(exec_type), intent(in) :: exec
    type(results_type), intent(in) :: results_prev
    type(results_type), intent(inout) :: results
    real(rk), intent(out) :: c(5, 5, elt_sub:elt_sup, nqpt)
    real(rk), intent(out) :: c_tan(5, 5, elt_sub:elt_sup, nqpt)
    real(rk), intent(out) :: eforce(5, elt_sub:elt_sup, nqpt)
    real(rk), intent(inout) :: dndx(ndim, elt_sub:elt_sup, nqpt)
    real(rk), intent(inout) :: dndy(ndim, elt_sub:elt_sup, nqpt)
    real(rk), intent(inout) :: dndz(ndim, elt_sub:elt_sup, nqpt)
    real(rk), intent(in) :: evel(kdim, elt_sub:elt_sup)
    integer, intent(in) :: incr
    real(rk), intent(in) :: dtime
    real(rk), intent(out) :: det(elt_sub:elt_sup, nqpt)
    real(rk), intent(in) :: ecoos(kdim, elt_sub:elt_sup)

    integer :: i

    !---------------------------------------------------------------------------

    do i = 1, nqpt
      call solve_evp(mesh, crys, exec, results_prev, i, results, incr, &
                     & dtime, ecoos, evel, dndx, dndy, dndz, det)

      call material_matrix_evp_aniso(mesh, crys, exec, i, results, &
          & c(:, :, :, i), c_tan(:, :, :, i), eforce(:, :, i), &
          & dtime)
    end do

  end subroutine material_matrix_evp

end module solveit_evp_mod3
