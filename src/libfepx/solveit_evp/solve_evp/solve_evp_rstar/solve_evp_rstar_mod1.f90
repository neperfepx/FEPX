! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module solve_evp_rstar_mod

! Module to update the ori (rstar).

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use hardening_mod
  use solve_evp_vpstress_mod
  use crys_type_mod2
  use matrix_operations_mod
  use orientation_conversion_mod
  use solve_evp_rstar_mod2

  implicit none

  public

contains

  subroutine solve_evp_rstar(mesh, crys, exec, results_prev, qpt, done, results, &
                       & dtime, euler_method)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys (:)
    type(exec_type), intent(in) :: exec
    type(results_type), intent(in) :: results_prev
    logical, intent(in)  :: done(elt_sub:elt_sup)
    integer, intent(in) :: qpt
    type(results_type), intent(inout) :: results
    real(rk), intent(in) :: dtime
    character(len = *), intent(in) :: euler_method

    !---------------------------------------------------------------------------

    ! Compute slip rates
    call compute_sliprate(mesh, crys, qpt, results)

    ! Update crss
    call update_hardening(mesh, crys, exec, results_prev, qpt, done, results, &
                & dtime, euler_method)

    ! Update rstar (r*) and d_rstar(dR*)
    call update_rstar(crys, mesh, results_prev, qpt, done, dtime, results)

    ! Compute c = c_0 * r*
    call update_ori(results_prev, qpt, done, results)

  end subroutine solve_evp_rstar

end module solve_evp_rstar_mod
