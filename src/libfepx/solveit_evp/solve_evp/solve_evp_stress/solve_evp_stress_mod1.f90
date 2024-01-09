! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module solve_stress_evp_mod

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use matrix_operations_mod
  use solve_stress_evp_mod2
  use solve_evp_vpstress_mod
  use units_mod
  use parallel_mod
  use crys_type_mod2

  implicit none

  public

contains

  subroutine solve_evp_stress(mesh, crys, exec, qpt, done, results, &
                             & dtime, converged_newton)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys (:)
    type(exec_type), intent(in) :: exec
    integer, intent(in) :: qpt
    logical, intent(in)     :: done(elt_sub:elt_sup)
    type(results_type), intent(inout) :: results
    real(rk), intent(in)    :: dtime
    logical, intent(out)  :: converged_newton

    integer  :: ier, num
    real(rk) :: sig0_avg, sig1_avg, sig2_avg, sig3_avg, sig4_avg
    logical  :: converged(elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    converged = .false.
    where (done) converged = .true.
    ! Solve nle for stresses

    ier = solve_stress_evp_newton(mesh, crys, exec, results, qpt, dtime, &
                                & converged, done)

    if (ier .ne. 0) then

      write (*, '(a)') 'Warning:       . solve_evp_stress: Some elements did not converge'
      !write(*, '(a)') 'Warning:       . These elements will receive average stress values'

      ! Adjust values to avoid divide-by-zero exceptions - deb
      num = max (count(converged), 1)

      sig0_avg = sum(results%sig_vec(1, :, qpt), mask=converged)/num
      sig1_avg = sum(results%sig_vec(2, :, qpt), mask=converged)/num
      sig2_avg = sum(results%sig_vec(3, :, qpt), mask=converged)/num
      sig3_avg = sum(results%sig_vec(4, :, qpt), mask=converged)/num
      sig4_avg = sum(results%sig_vec(5, :, qpt), mask=converged)/num

      where (.not. converged)
        results%sig_vec(1, :, qpt) = sig0_avg
        results%sig_vec(2, :, qpt) = sig1_avg
        results%sig_vec(3, :, qpt) = sig2_avg
        results%sig_vec(4, :, qpt) = sig3_avg
        results%sig_vec(5, :, qpt) = sig4_avg
      end where

      converged_newton = .false.

    end if

  end subroutine solve_evp_stress

end module solve_stress_evp_mod
