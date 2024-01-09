! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module driver_triax_utilities_mod

! This module defines common subroutines called by the various driver modules.

! Contains subroutines:
! est_avg_mod: Estimate bulk elastic moduli. Assumes a uniform texture and
!   equal valued element volumes.
! vel_iteration: Perform iteration on vel field.

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use solveit_evp_mod
  use kinematics_mod
  use matrix_operations_mod
  use solve_evp_mod
  use quadrature_mod
  use shape_3d_mod
  use surface_mod
  use units_mod
  use gather_scatter_mod
  use parallel_mod

  use finalize_res_mod

  implicit none

  private

  public :: vel_iteration
  public :: est_avg_mod

contains

  subroutine vel_iteration(mesh, crys, loading, exec, &
      & results_prev, results, printing, load, dtime, incr)

    ! Perform vel iteration

    !---------------------------------------------------------------------------

    ! Arguments:

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys (:)
    type(loading_type), intent(in) :: loading
    type(exec_type), intent(inout) :: exec
    type(results_type), intent(in) :: results_prev
    type(results_type), intent(inout) :: results
    type(printing_type), intent(in) :: printing
    real(rk), intent(out) :: load(3)
    real(rk) :: bak_sig_kk(elt_sub:elt_sup, nqpt)
    real(rk) :: dtime
    integer, intent(in) :: incr

    ! Locals:

    real(rk) :: surf_load_array(mesh%num_fasets, 3)
    real(rk) :: area(mesh%num_fasets)
    real(rk) :: bak_oris(3, 3, elt_sub:elt_sup, nqpt)
    real(rk) :: bak_rstar(3, 3, elt_sub:elt_sup, nqpt)
    real(rk) :: bak_coos(dof_sub:dof_sup)
    real(rk) :: bak_sig_vec(5, elt_sub:elt_sup, nqpt)
    real(rk) :: bak_e_elas_kk(elt_sub:elt_sup, nqpt)
    real(rk) :: bak_d_kk(elt_sub:elt_sup, nqpt)

    ! Initialization

    bak_oris = results%ori
    bak_coos = results%coo
    bak_e_elas_kk = results%e_elas_kk_bar
    bak_sig_vec = results%sig_vec
    bak_sig_kk = results%sig_kk
    bak_rstar = results%rstar
    bak_d_kk = results%d_kk

    ! Assign temporary variables for solveit_evp and temp_update_state_evp
    results%sig_vec = results_prev%sig_vec

    ! Iterate on velocity field
    call solveit_evp(mesh, crys, loading, exec, printing, &
                      & results_prev, results, dtime, incr, surf_load_array, &
                      & area, "nofinalize")

    ! Compute load
    call finalize_res_stressstrain(mesh, crys, results, surf_load_array, area)

    load(1) = surf_load_array(2, 1) ! 4,1
    load(2) = surf_load_array(4, 2) ! 6,2
    load(3) = surf_load_array(6, 3) ! 2,3

    ! Reset coordinates etc.
    results%coo = bak_coos
    results%sig_kk = bak_sig_kk
    results%e_elas_kk_bar = bak_e_elas_kk
    results%sig_vec = bak_sig_vec
    results%ori = bak_oris
    results%rstar = bak_rstar
    results%d_kk = bak_d_kk

  end subroutine vel_iteration

  subroutine est_avg_mod(mesh, crys, e_avg, nu_avg)

    ! Estimate Voigt-averaged bulk elastic moduli (Hosford p22)
    ! Assumes uniform texture and equal element volumes.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! e_avg:
    ! nu_avg:

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys
    real(rk), intent(out) :: e_avg
    real(rk), intent(out) :: nu_avg

    ! Locals:
    ! my_phase:
    ! i/iphase:
    ! part_numel_phase:
    ! numel_phase:
    ! phase_frac:
    ! k_phase:
    ! e_phase:
    ! nu_phase:
    ! c11, c12, c13, c33, c44, c66:
    ! f, g, h:

    integer :: my_phase(elt_sub:elt_sup)
    integer :: i, iphase
    real(rk) :: part_numel_phase
    real(rk) :: numel_phase
    real(rk) :: phase_frac
    real(rk) :: k_phase
    real(rk) :: e_phase
    real(rk) :: nu_phase
    real(rk) :: c11
    real(rk) :: c12
    real(rk) :: c13
    real(rk) :: c33
    real(rk) :: c44
    real(rk) :: c66
    real(rk) :: f
    real(rk) :: g
    real(rk) :: h

    !---------------------------------------------------------------------------

    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    e_avg = 0.0d0
    nu_avg = 0.0d0

    do iphase = 1, mesh%num_phases
      part_numel_phase = 0.0d0

      do i = elt_sub, elt_sup
        if (my_phase(i) .eq. iphase) then
          part_numel_phase = part_numel_phase + 1.0
        end if
      end do

      call par_sum(part_numel_phase, numel_phase)
      phase_frac = numel_phase/real(mesh%num_elts, rk)

      if ((crys%structure .eq. "fcc") .or. (crys%structure .eq. "bcc")) then
        c11 = crys%elas_coeffs(1)
        c12 = crys%elas_coeffs(2)
        c44 = crys%elas_coeffs(4)
        e_phase = (c11 - c12 + 3*c44)*(c11 + 2*c12)/ &
            & (2*c11 + 3*c12 + c44)

      else if (crys%structure .eq. "hcp") then ! Hexagonal (hcp)
        c11 = crys%elas_coeffs(1)
        c12 = crys%elas_coeffs(2)
        c13 = crys%elas_coeffs(3)
        c44 = crys%elas_coeffs(4)
        c33 = c11 + c12 - c13
        c66 = (c11 - c12)/2

        f = (2*c11 + c33)/3
        g = (c12 + 2*c13)/3
        h = (2*c44 + c66)/3
        e_phase = (f - g + 3*h)*(f + 2*g)/(2*f + 3*g + h)

      else if (crys%structure .eq. "bct") then ! Tetragonal (bct)
        c11 = crys%elas_coeffs(1)
        c12 = crys%elas_coeffs(2)
        c13 = crys%elas_coeffs(3)
        c44 = crys%elas_coeffs(4)
        c66 = crys%elas_coeffs(5)
        c33 = c11 + c12 - c13

        f = (2*c11 + c33)/3
        g = (c12 + 2*c13)/3
        h = (2*c44 + c66)/3
        e_phase = (f - g + 3*h)*(f + 2*g)/(2*f + 3*g + h)

      else
        call par_quit('Error  :     > Invalid crystal type.')
      end if

      k_phase = crys%bulk_mod
      nu_phase = (3*k_phase - e_phase)/(6*k_phase)

      e_avg = e_avg + phase_frac*e_phase
      nu_avg = nu_avg + phase_frac*nu_phase
    end do

    return

  end subroutine est_avg_mod

end module driver_triax_utilities_mod
