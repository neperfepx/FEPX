! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module solve_evp_mod

! Module hangling elastic-viscoplastic response for polycrystals.

! Contains subroutines:
! solve_evp: evp response for polycrystals
! solve_evp: evp response for polycrystal at quad points

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use solve_stress_evp_mod, only: solve_evp_stress

  use solve_evp_mod2

  implicit none

  public

contains

  subroutine solve_evp(mesh, crys, exec, results_prev, qpt, results, &
                     & incr, dtime, ecoos, evel, dndx, dndy, dndz, det)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys(:)
    type(exec_type), intent(in) :: exec
    type(results_type), intent(in) :: results_prev
    type(results_type), intent(inout) :: results
    integer, intent(in) :: qpt
    real(rk), intent(in) :: ecoos(kdim, elt_sub:elt_sup)
    real(rk), intent(in) :: evel(kdim, elt_sub:elt_sup)
    real(rk), intent(inout) :: dndx(ndim, elt_sub:elt_sup, nqpt)
    real(rk), intent(inout) :: dndy(ndim, elt_sub:elt_sup, nqpt)
    real(rk), intent(inout) :: dndz(ndim, elt_sub:elt_sup, nqpt)
    real(rk), intent(out) :: det(elt_sub:elt_sup, nqpt)

    integer :: incr
    real(rk) :: dtime

    call solve_evp_pre(mesh, crys, results_prev, qpt, results, ecoos, evel, &
                     & dndx, dndy, dndz, det)

    call solve_evp_dev(mesh, crys, exec, results_prev, qpt, results, incr, dtime)

    call solve_evp_vol(mesh, crys, results_prev, qpt, results, dtime)

  end subroutine solve_evp

end module solve_evp_mod
