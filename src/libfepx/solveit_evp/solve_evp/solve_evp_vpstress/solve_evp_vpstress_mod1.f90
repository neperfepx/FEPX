! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module solve_evp_vpstress_mod

! Routines for solving viscoplastic crystal stress equations

! Contains subroutines:
! solve_evp_vpstress: Routine which takes care of scaling and initial guesses
! solve_newton_vp: Nonlinear solution of viscoplastic crystal stress equations
! scale_up_sigm: Rescale the stress after solution is found

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use matrix_operations_mod
  use units_mod
  use parallel_mod
  use crys_type_mod2
  use stiffness_mod
  use solve_evp_vpstress_mod2

  implicit none

  public

contains

  !> Driver routine which takes care of scaling and initial guesses
  subroutine solve_evp_vpstress(mesh, crys, exec, qpt, results, vp_log)

    type(mesh_type) :: mesh
    type(crys_type) :: crys (:)
    type(exec_type), intent(in) :: exec
    integer, intent(in) :: qpt
    type(results_type), intent(inout) :: results
    logical, intent(in) :: vp_log

    logical :: converged(elt_sub:elt_sup)
    integer :: i_edge
    integer :: ier
    integer :: ntrials
    integer :: vertex(elt_sub:elt_sup)
    real(rk) :: direction(elt_sub:elt_sup)
    integer :: maxnumvert
    real(rk), allocatable :: plwork(:,:)
    real(rk) :: sig_t(5, elt_sub:elt_sup)
    real(rk) :: crss_avg(elt_sub:elt_sup)
    real(rk) :: qr5x5(5, 5, elt_sub:elt_sup)
    real(rk) :: d_vec_lat(5, elt_sub:elt_sup)
    integer :: m_elt

    m_elt = elt_sup - elt_sub + 1

    ! results%ori(:, :, :, qpt) [3x3] --> qr5x5 [5x5]
    ! d_vec(sample coords) --> d_vec_lat(crystal coords)
    ! {d_vec_lat} = [qr5x5]'{d_vec}
    call rotmat_symm(results%ori(:, :, :, qpt), qr5x5, m_elt)
    call lattice_deform(qr5x5, results%d_vec(:, :, qpt), d_vec_lat)

    !---------------------------------------------------------------------------
    call crys_maxnumvert (crys, maxnumvert)

    allocate(plwork(maxnumvert, elt_sub:elt_sup))


    call scale_down_defr(d_vec_lat, results%defrate_eq(:, qpt))
    call compute_work(mesh, crys, plwork, d_vec_lat)

    converged = .false.

    ! This loop used to be from 1 to n_edge, in order to try all vertices as
    ! initial guesses. However, this can be very time-consuming and doesn't
    ! usually help.

    ntrials = 4

    do i_edge = 1, ntrials

      call find_vertex(mesh, crys, vertex, direction, plwork)

      call vertex_stress(mesh, crys, sig_t, vertex, direction)

      where (.not. converged)
        results%sig_vec(1, :, qpt) = sig_t(1, :)
        results%sig_vec(2, :, qpt) = sig_t(2, :)
        results%sig_vec(3, :, qpt) = sig_t(3, :)
        results%sig_vec(4, :, qpt) = sig_t(4, :)
        results%sig_vec(5, :, qpt) = sig_t(5, :)
      end where

      call solve_newton_vp(mesh, crys, exec, results%sig_vec(:, :, qpt), d_vec_lat, crss_avg, ier, &
        & exec%sx_tol, converged, vp_log)

      if (ier .eq. 0) then
        exit
      end if

    end do

    if (ier .ne. 0) then
      if ((myid .eq. 0) .and. vp_log) then
        write (*, '(a)') 'Warning:       . solve_evp_vpstress: '
        write (*, '(a)') 'Warning:       . ', count(.not. converged), '&
          & elements did not converge'
        write (*, '(a)') 'Warning:       . after ', ntrials, ' trials.'
      end if

      where (.not. converged)
        results%sig_vec(1, :, qpt) = sig_t(1, :)
        results%sig_vec(2, :, qpt) = sig_t(2, :)
        results%sig_vec(3, :, qpt) = sig_t(3, :)
        results%sig_vec(4, :, qpt) = sig_t(4, :)
        results%sig_vec(5, :, qpt) = sig_t(5, :)
      end where
    end if

    call scale_up_sigm(mesh, results%sig_vec(:, :, qpt), results%defrate_eq(:, qpt))

  end subroutine solve_evp_vpstress

end module solve_evp_vpstress_mod
