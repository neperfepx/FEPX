! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module boundary_conditions_mod

! Module to calculate essential boundary conditions

! Contains subroutines:
! calc_bcs: Calculate grip, symmetry, or initial triaxial boundary conditions.

  use general_mod
  use types_mod
  use, intrinsic :: iso_fortran_env, only: rk => real64

! From libfepx:

  use parallel_mod, only: par_max, par_min, par_quit, par_message
  use types_mod
  use boundary_conditions_mod2
  use boundary_conditions_mod2_legacy

  implicit none

  private

  public :: calc_bcs

contains

  !> Calculate boundary conditions from loading_options
  subroutine calc_bcs(mesh, loading_options, loading)

    type(mesh_type), intent(in) :: mesh
    type(loading_options_type), intent(inout) :: loading_options
    type(loading_type), intent(inout) :: loading

    ! Notes:
    ! The mesh face set is always defined from [1 2 3 4 5 6] as
    ! [x_min x_max y_min y_max z_min z_max]. This should be consistent
    ! with the face set provided during mesh generation. The loading
    ! directions are defined along the sample axes [x y z] as [0 1 2].

    ! This subroutine by default assumes that displacement occurs in a
    ! positive axis direction (that is, a domain loading in tension or
    ! sheared along a positive coordinate axis). In order to perform
    ! compression or negative axis direction shearing tests, the user
    ! must sign the bcs_vel in their *.loads or *.disp files.

    ! Initialize variables
    allocate (loading%bcs_vel_defined(dof_sub:dof_sup))
    allocate (loading%bcs_vel(dof_sub:dof_sup))
    loading%bcs_vel_defined = .false.
    loading%bcs_vel = 0.0d0

    loading%curr_load = 0.0d0
    loading%prev_load = 0.0d0

    loading%def_control_by = loading_options%def_control_by

    if (loading_options%num_mpcs .gt. 0) then
      loading%mpc_status = .true.
    end if

    loading%loading_face = 2 * loading_options%loading_direction

    ! triaxial related
    loading%dwell_max_strain_incr = loading_options%dwell_max_strain_incr

    loading%loading_direction = loading_options%loading_direction
    loading%loading_face = 2 * loading%loading_direction

    ! Error handling to ensure loading_options%def_control_by and boundary_conditions agree.
    if ((loading_options%def_control_by(1:8) .eq. "uniaxial") &
      & .and. (loading_options%boundary_conditions .eq. "triaxial")) then
      call par_quit('Error  :     > Parameters "def_control_by" and&
                  & "boundary_conditions" do not agree.')

    else if ((loading_options%def_control_by(1:8) .eq. "triaxial") &
      & .and. (loading_options%boundary_conditions .ne. "triaxial")) then
      call par_quit('Error  :     > Parameters "def_control_by" and&
                  & "boundary_conditions" do not agree.')
    end if

    !---------------------------------------------------------------------------

    select case (loading_options%boundary_conditions)

    ! -- legacy code that will call calc_bcs_general internally -------------
    case ("uniaxial_grip")
      call calc_bcs_uniaxial_grip(mesh, loading_options, loading)

    case ("uniaxial_symmetry")
      call calc_bcs_uniaxial_symmetry(mesh, loading_options, loading)

    case ("uniaxial_minimal")
      call calc_bcs_uniaxial_minimal(mesh, loading_options, loading)

    case ("triaxial") ! Begin Case (triaxial)
      call calc_bcs_triaxial(mesh, loading_options, loading)
    ! -- end of legacy code that will call calc_bcs_general internally ----

    case ("periodic")
      if (myid .eq. 0) then
        write (*, '(a)') 'Info   :   - Periodic Boundary Conditions (PBC):'
      end if

    case default
      call calc_bcs_general(mesh, loading_options, loading)

    end select

    select case (loading_options%def_control_by)

    ! -- legacy code -----------
    case ("triaxial_constant_strain_rate")
      call process_ctrl_data_csr(loading_options, loading)

    case ("triaxial_constant_load_rate")
      call process_ctrl_data_clr (loading_options, loading)
    ! -- end of legacy code -----------

    case default
      call calc_bcs_steps(loading_options, mesh, loading)

    end select

  end subroutine calc_bcs

end module boundary_conditions_mod
