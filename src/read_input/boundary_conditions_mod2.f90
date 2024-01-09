! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module boundary_conditions_mod2

  use general_mod
  use mesh_type_mod
  use types_mod
  use, intrinsic :: iso_fortran_env, only: rk => real64

  use parallel_mod, only: par_min, par_max, par_quit
  use types_mod
  use utils_mod
  use boundary_conditions_mod3

  implicit none

contains

  ! ============================================================================

  !> Calculate general boundary conditions
  subroutine calc_bcs_general(mesh, loading_options, loading)

    type(mesh_type), intent(in) :: mesh
    type(loading_options_type), intent(inout) :: loading_options
    type(loading_type), intent(inout) :: loading

    integer :: i, j, ind(2)
    character(len=2) :: temp

    real(rk) :: L(3,3)
    logical :: L_defined(3,3)
    integer :: L_defined_bcindex(3,3)

    ! Loop on vel bcs
    do i = 1, loading_options%num_bcs
      if (loading_options%bc_var(i) .eq. "vel") then
        call calc_bcs_general_vel(mesh, loading_options, i, loading)
      end if
    end do

    ! Loop on strainrate bcs
    L_defined = .false.
    L = 0.0d0
    do i = 1, loading_options%num_bcs
      if (loading_options%bc_var(i) .eq. "strainrate") then

        do j = 1, 2
          temp = loading_options%bc_dir(i)
          ind(j) = ichar (temp(j:j)) - ichar ('0')
        end do

        if (ind(1) .ne. ind(2)) then
          call par_quit('Error  :     > Shear components not available for strainrate (see velgrad instead).')
        end if

        L(ind(1),ind(2)) = loading_options%bc_vel(i)
        L_defined(ind(1),ind(2)) = .true.
        L_defined_bcindex(ind(1),ind(2)) = i

      else if (loading_options%bc_var(i) .eq. "velgrad") then

        do j = 1, 2
          temp = loading_options%bc_dir(i)
          ind(j) = ichar (temp(j:j)) - ichar ('0')
        end do

        L(ind(1),ind(2)) = loading_options%bc_vel(i)
        L_defined(ind(1),ind(2)) = .true.
        L_defined_bcindex(ind(1),ind(2)) = i
      end if
    end do

    do i = 1, 3
      do j = 1, 3
        if (L_defined(i,j)) then
          call calc_bcs_general_velgrad(mesh, loading_options, L, L_defined, L_defined_bcindex, loading)
        end if
      end do
    end do

  end subroutine calc_bcs_general

  ! ============================================================================

  subroutine calc_bcs_steps(loading_options, mesh, loading)

    type(loading_options_type), intent(in) :: loading_options
    type(mesh_type), intent(in) :: mesh
    type(loading_type), intent(inout) :: loading

    integer :: i
    real(rk) :: diff_stress
    real(rk), allocatable :: diff_strain(:), diff_time(:)
    real(rk) :: max_vel

    loading%num_steps = loading_options%num_steps

    if (loading%pbc_status .eqv. .false.) then
      
      allocate (loading%target_time(loading%num_steps))
      allocate (loading%step_velocity(loading%num_steps))
      allocate (loading%step_dt_min(loading%num_steps))
      allocate (loading%step_dt_max(loading%num_steps))
      allocate (loading%target_sign(loading%num_steps))
      allocate (loading%step_print(loading%num_steps))
      allocate (loading%vel_factor(loading%num_steps))
      
      call alloc_1d (diff_time, 1, loading_options%num_steps)

      ! step_velocity is computed as the largest bc velocity
      max_vel = 0.0d0
      do i = 1, loading_options%num_bcs
        if (max_vel .lt. abs (loading_options%bc_vel(i))) then
          max_vel = abs (loading_options%bc_vel(i))
        end if
      end do

      loading%step_velocity = max_vel

      ! FIXME disp
      loading%step_print = loading_options%step_print
    
    end if

    select case (loading_options%def_control_by)

    case ('uniaxial_strain_target')

      ! Allocate the arrays for the strain target sequence.

      if (loading%pbc_status .eqv. .false.) then
        allocate (loading%target_disp(loading%num_steps))
      end if

      if (loading_options%loading_direction .ne. 0) then
        if (loading%pbc_status .eqv. .false.) then
          call mesh_length(mesh, loading_options%loading_direction, loading%gage_length)
        end if
        loading%target_disp = loading_options%target_strain*loading%gage_length
      end if

      ! diff_strain ------------------------------------------------------------

      call alloc_1d (diff_strain, 1, loading_options%num_steps)

      do i = 1, loading%num_steps

        if (i .eq. 1) then
          diff_strain(i) = loading_options%target_strain(i)
        else
          diff_strain(i) = loading_options%target_strain(i) - loading_options%target_strain(i - 1)
        end if

      end do

      ! target_sign ------------------------------------------------------------

      loading%target_sign = 1

      if (loading_options%target .eq. "strain") then
        do i = 1, loading%num_steps

          if (diff_strain(i) .gt. 0.0d0) then
            loading%target_sign(i) = 1
          else
            loading%target_sign(i) = -1
          end if

        end do
      end if

      ! vel_factor -------------------------------------------------------------

      do i = 1, loading%num_steps
        ! FIXME disp
        if (i .eq. 1) then
          loading%vel_factor(i) = 1.0d0
        else
          loading%vel_factor(i) = loading%step_velocity(i) / loading%step_velocity(i - 1)
          if (loading%target_sign(i) .ne. loading%target_sign(i - 1)) then
            loading%vel_factor(i) = -loading%vel_factor(i)
          end if
        end if
      end do

      ! bcs_vel ----------------------------------------------------------------

      do i = 1, loading%num_steps

        if (i .eq. 1) then
          if (loading%target_sign(i) .eq. -1) then
            loading%bcs_vel = -loading%bcs_vel
          end if
        end if

      end do

      ! target_time ------------------------------------------------------------

      if (loading_options%target .eq. "time") then

        loading%target_time = loading_options%target_time

      else if (loading_options%target .eq. "strain") then

        if (loading%pbc_status .eqv. .false.) then
          call mesh_length(mesh, loading_options%loading_direction, loading%gage_length)
        end if

        do i = 1, loading%num_steps
          if (i .eq. 1) then
            diff_strain(i) = loading_options%target_strain(i)
          else
            diff_strain(i) = loading_options%target_strain(i) - loading_options%target_strain(i - 1)
          end if
        end do

        do i = 1, loading%num_steps
          loading%target_time(i) = abs(diff_strain(i)*loading%gage_length) / loading%step_velocity(i)
        end do
        do i = 2, loading%num_steps
          loading%target_time(i) = loading%target_time(i) + loading%target_time(i - 1)
        end do

      end if

      ! time step (step_dt_max) ------------------------------------------------

      if (allocated (loading_options%step_dt_max)) then
        loading%step_dt_max = loading_options%step_dt_max
        loading%step_dt_min = loading%step_dt_max

      else if (allocated (loading_options%step_num_incr)) then
        do i = 1, loading%num_steps
          if (i .eq. 1) then
            diff_time(i) = loading%target_time(i)
          else
            diff_time(i) = loading%target_time(i) - loading%target_time(i - 1)
          end if
        end do

        do i = 1, loading%num_steps
          loading%step_dt_min(i) = diff_time(i) / loading_options%step_num_incr(i)
          loading%step_dt_max(i) = diff_time(i) / loading_options%step_num_incr(i)
        end do

      else if (allocated (loading_options%step_dstrain_max)) then
        if (loading%pbc_status .eqv. .false.) then
          call mesh_length(mesh, loading_options%loading_direction, loading%gage_length)
        end if
        do i = 1, loading%num_steps
          loading%step_dt_max(i) = abs (loading_options%step_dstrain_max(i) * loading%gage_length / (loading%step_velocity(i)))
          loading%step_dt_min(i) = loading%step_dt_max(i)
        end do
      end if

      ! ------------------------------------------------------------------------

    case ('uniaxial_load_target')
      ! Initialize definition variables.

      diff_stress = 0.0d0

      ! Allocate the arrays for the load target sequence.

      if (loading%pbc_status .eqv. .false.) then
        allocate (loading%step_target_time_incr(loading%num_steps))
        allocate (loading%target_load(loading%num_steps, 3))
      end if

      loading%target_load = loading_options%target_load

      loading%step_dt_max = loading_options%step_dt_max
      loading%step_dt_min = loading_options%step_dt_min

      do i = 1, loading%num_steps
        ! Handle the first step alone to avoid out of bounds on arrays.

        if (i .eq. 1) then
          loading%vel_factor(i) = 1.0d0
        else
          loading%vel_factor(i) = loading%step_velocity(i)/loading%step_velocity(i - 1)
        end if

        if (i .eq. 1) then
          diff_stress = loading_options%target_load(i, 1)
        else
          diff_stress = loading_options%target_load(i, 1) &
                    & - loading_options%target_load(i - 1, 1)
        end if

        if (diff_stress .eq. 0.0d0) then
          write (*,*) "i = ", i, loading_options%target_load(i, 1)
          call par_quit('Error  :     > Load differential is zero between steps.')
        end if

        if (diff_stress .gt. 0.0d0) then
          loading%target_sign(i) = 1
        else
          loading%target_sign(i) = -1
        end if

        if (i .eq. 1) then
          if (diff_stress .lt. 0.0d0) then
            loading%bcs_vel = -loading%bcs_vel * loading%vel_factor(i)
          end if
        else
          if (loading%target_sign(i) .ne. loading%target_sign(i - 1)) then
            loading%vel_factor(i) = -loading%vel_factor(i)
          end if
        end if
      end do

    case default
      call par_quit('Error  :     > Invalid control option defined.')

    end select

    ! Define initial values for driver.

    loading%curr_step = 1
    loading%step_complete = .false.
    loading%all_steps_complete = .false.
    loading%prev_dt = loading%step_dt_max(1)

  end subroutine calc_bcs_steps

  ! ============================================================================

end module boundary_conditions_mod2
