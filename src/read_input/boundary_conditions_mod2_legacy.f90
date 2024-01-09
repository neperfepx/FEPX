! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module boundary_conditions_mod2_legacy

  use general_mod
  use mesh_type_mod
  use types_mod
  use, intrinsic :: iso_fortran_env, only: rk => real64

  use parallel_mod, only: par_min, par_max, par_quit
  use types_mod
  use utils_mod

  use boundary_conditions_mod2

  implicit none

contains

  ! ============================================================================

  !> Calculate uniaxial grip boundary conditions
  subroutine calc_bcs_uniaxial_grip(mesh, loading_options, loading)

    type(mesh_type), intent(in) :: mesh
    type(loading_options_type), intent(inout) :: loading_options
    type(loading_type), intent(inout) :: loading

    real(rk) :: lengths(3)

    call mesh_lengths (mesh, lengths)
    loading%gage_length = lengths(loading_options%loading_direction)

    loading_options%bc_var = "vel"
    loading_options%num_bcs = 6

    ! z0
    loading_options%bc_nset(1) = "z0"
    loading_options%bc_nset(2) = "z0"
    loading_options%bc_nset(3) = "z0"
    loading_options%bc_dir(1)  = "x"
    loading_options%bc_dir(2)  = "y"
    loading_options%bc_dir(3)  = "z"
    ! vel is already 0

    ! z1
    loading_options%bc_nset(4) = "z1"
    loading_options%bc_nset(5) = "z1"
    loading_options%bc_nset(6) = "z1"
    loading_options%bc_dir(4)  = "x"
    loading_options%bc_dir(5)  = "y"
    loading_options%bc_dir(6)  = "z"
    loading_options%bc_vel(6) =  loading_options%strain_rate*loading%gage_length

    call calc_bcs_general(mesh, loading_options, loading)

  end subroutine calc_bcs_uniaxial_grip

  ! ============================================================================

  !> Calculate uniaxial symmetry boundary conditions
  subroutine calc_bcs_uniaxial_symmetry(mesh, loading_options, loading)

    type(mesh_type), intent(in) :: mesh
    type(loading_options_type), intent(inout) :: loading_options
    type(loading_type), intent(inout) :: loading

    real(rk) :: lengths(3)

    call mesh_lengths (mesh, lengths)
    loading%gage_length = lengths(loading_options%loading_direction)

    loading_options%bc_var = "vel"
    loading_options%num_bcs = 4

    loading_options%bc_nset(1) = "x0"
    loading_options%bc_dir(1) = "x"
    loading_options%bc_vel(1) = 0

    loading_options%bc_nset(2) = "y0"
    loading_options%bc_dir(2) = "y"
    loading_options%bc_vel(2) = 0

    loading_options%bc_nset(3) = "z0"
    loading_options%bc_dir(3) = "z"
    loading_options%bc_vel(3) = 0

    loading_options%bc_nset(4) = "z1"
    loading_options%bc_dir(4) = "z"
    loading_options%bc_vel(4) =  loading_options%strain_rate*loading%gage_length

    call calc_bcs_general(mesh, loading_options, loading)

  end subroutine calc_bcs_uniaxial_symmetry

  ! ============================================================================

  !> Calculate uniaxial minimal boundary conditions
  ! simple tension with constraints to prevent rigid-body motions
  subroutine calc_bcs_uniaxial_minimal(mesh, loading_options, loading)

    type(mesh_type), intent(in) :: mesh
    type(loading_options_type), intent(inout) :: loading_options
    type(loading_type), intent(inout) :: loading

    real(rk) :: lengths(3)

    call mesh_lengths (mesh, lengths)
    loading%gage_length = lengths(loading_options%loading_direction)

    loading_options%bc_var = "vel"
    loading_options%num_bcs = 4

    ! vels already 0

    loading_options%bc_nset(1) = "z0"
    loading_options%bc_dir(1) = "z"

    loading_options%bc_nset(2) = "z1"
    loading_options%bc_dir(2) = "z"
    loading_options%bc_vel(2) =  loading_options%strain_rate*loading%gage_length

    loading_options%bc_nset(3) = "x0y0z0"
    loading_options%bc_dir(3) = "x"

    loading_options%bc_nset(4) = "x0y0z0"
    loading_options%bc_dir(4) = "y"

    loading_options%bc_nset(5) = "x1y0z0"
    loading_options%bc_dir(5) = "y"

    call calc_bcs_general(mesh, loading_options, loading)

  end subroutine calc_bcs_uniaxial_minimal

  ! ============================================================================

  !> Calculate triaxial boundary conditions
  subroutine calc_bcs_triaxial(mesh, loading_options, loading)

    type(mesh_type), intent(in) :: mesh
    type(loading_options_type), intent(inout) :: loading_options
    type(loading_type), intent(inout) :: loading

    loading_options%bc_var = "vel"
    loading_options%num_bcs = 6

    ! no vel?
    loading_options%bc_nset(1)= "x0"
    loading_options%bc_dir(1) = "x"
    loading_options%bc_nset(2)= "x1"
    loading_options%bc_dir(2) = "x"
    loading_options%bc_nset(3)= "y0"
    loading_options%bc_dir(3) = "y"
    loading_options%bc_nset(4)= "y1"
    loading_options%bc_dir(4) = "y"
    loading_options%bc_nset(5)= "z0"
    loading_options%bc_dir(5) = "z"
    loading_options%bc_nset(6)= "z1"
    loading_options%bc_dir(6) = "z"

    call calc_bcs_general(mesh, loading_options, loading)

  end subroutine calc_bcs_triaxial

  ! ============================================================================

  subroutine process_ctrl_data_clr (loading_options, loading)

    ! Process input data for load control at constant load rate.

    !---------------------------------------------------------------------------

    ! Locals:
    ! i/ii: Generic loop index values.
    ! iload: Index value used to access target_load sequentially by row.
    ! idwell: Index value used to access dwell_episode sequentially by row.
    ! step_load_rate: Array used to store ramp rate jump values.

    type(loading_options_type), intent(inout) :: loading_options
    type(loading_type), intent(inout) :: loading
    integer :: i, ii, iload, idwell
    real(rk), allocatable :: step_load_rate(:)

    ! Notes:
    ! User-defined input (ramp rates) should always be positive. Compression
    ! and similar signed deformation control are handled internally within the
    ! primary driver subroutine.

    !---------------------------------------------------------------------------

    ! Set loading%num_steps to account for possible dwell episodes to adjust the array

    loading%num_steps = loading_options%num_steps + &
        & loading_options%number_of_dwell_episodes

    ! Allocate the necessary arrays
    allocate (loading%target_load(loading%num_steps, 3))
    allocate (loading%control_dir(loading%num_steps))
    allocate (loading%load_rate(loading%num_steps))
    allocate (loading%dwell_time(loading%num_steps))
    allocate (loading%step_target_time_incr(loading%num_steps))
    allocate (loading%step_print(loading%num_steps))

    allocate (step_load_rate(loading_options%num_steps))

    ! Initialize arrays.

    loading%target_load = 0.0d0
    loading%control_dir = 0
    loading%control_dir = loading_options%loading_direction
    loading%load_rate = loading_options%load_rate
    loading%dwell_time = 0.0d0
    loading%step_target_time_incr = 0.0d0
    loading%step_print = .false.

    step_load_rate = 0.0d0
    step_load_rate = loading_options%load_rate

    ! Build load rate array, specifically handle load rate jumps.

    if (loading_options%number_of_load_rate_jumps .gt. 0) then
      ii = 1

      do i = 1, loading_options%num_steps
        if (loading_options%load_rate_jump(ii, 1) .eq. i) then
          step_load_rate(i) = &
              & loading_options%load_rate_jump(ii, 2)

          if (ii .eq. loading_options%number_of_load_rate_jumps) then
            ii = ii

          else
            ii = ii + 1
          end if
        end if
      end do
    end if

    ! Check that a dwell doesn't occur before the first loading step

    if (loading_options%number_of_dwell_episodes .gt. 0) then
      if (loading_options%dwell_episode(1, 1) .eq. 0) then
        call par_quit('Error  :     > a dwell episode can not be applied &
            &to an unloaded specimen on the first step.')
      end if
    end if

    ! Increment the first column (step index) of dwell episodes so it can be
    ! accessed in a logical order to build the expected arrays properly

    if (loading_options%number_of_dwell_episodes .gt. 0) then
      do i = 1, loading_options%number_of_dwell_episodes
        loading_options%dwell_episode(i, 1) = &
            & loading_options%dwell_episode(i, 1) + 1
      end do
    end if

    ! Initialize additional index values (used to access load and dwell arrays)

    iload = 1
    idwell = 1

    ! Build the input arrays for all loads (including possible dwells)

    do i = 1, loading%num_steps
      ! Check if the current index value is a dwell episode

      if (loading_options%number_of_dwell_episodes .gt. 0) then
        if (i .eq. loading_options%dwell_episode(idwell, 1)) then
          ! Assign the load values from the previous step

          loading%target_load(i, 1) = loading%target_load(i - 1, 1)
          loading%target_load(i, 2) = loading%target_load(i - 1, 2)
          loading%target_load(i, 3) = loading%target_load(i - 1, 3)

          ! Set the ramp rate to zero for dwell

          loading%load_rate(i) = 0.0d0
          loading%control_dir(i) = -1*loading%control_dir(i)

          loading%dwell_time(i) = loading_options%dwell_episode(idwell, 2)
          loading%step_target_time_incr(i) = loading_options%dwell_episode(idwell, 3)

          if (loading_options%dwell_episode(idwell, 4) .eq. 1) then
            loading%step_print(i) = .true.
          else
            loading%step_print(i) = .false.
          end if

          ! Increment the index idwell to access next available row
          ! Check if max is reached to avoid out of bounds array access

          if (idwell .eq. loading_options%number_of_dwell_episodes) then
            idwell = idwell

          else
            idwell = idwell + 1
          end if

        else ! If not a dwell episode than process a typical load step
          ! Assign the load values for this step

          loading%target_load(i, 1) = loading_options%target_load(iload, 1)
          loading%target_load(i, 2) = loading_options%target_load(iload, 2)
          loading%target_load(i, 3) = loading_options%target_load(iload, 3)

          loading%load_rate(i) = step_load_rate(iload)
          loading%step_target_time_incr(i) = loading_options%step_target_time_incr(iload)
          loading%step_print(i) = loading_options%step_print(i)

          ! Increment the index iload to access next available row
          ! Check if max is reached to avoid out of bounds array access

          if (iload .eq. loading_options%num_steps) then
            iload = iload

          else
            iload = iload + 1
          end if
        end if
      end if

      if (loading_options%number_of_dwell_episodes .eq. 0) then
        loading%target_load(i, 1) = loading_options%target_load(iload, 1)
        loading%target_load(i, 2) = loading_options%target_load(iload, 2)
        loading%target_load(i, 3) = loading_options%target_load(iload, 3)

        loading%load_rate(i) = step_load_rate(iload)
        loading%step_target_time_incr(i) = loading_options%step_target_time_incr(iload)

        loading%step_print(i) = loading_options%step_print(i)

        ! Increment the index iload to access next available row
        ! Check if max is reached to avoid out of bounds array access

        if (iload .eq. loading_options%num_steps) then
          iload = iload

        else
          iload = iload + 1
        end if
      end if
    end do

    return

  end subroutine process_ctrl_data_clr

  ! ============================================================================

  subroutine process_ctrl_data_csr(loading_options, loading)

    ! Read data for triaxial load control at constant strain rate.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! loading%bcs_vel: Global d.o.f. array storing loading%bcs_vel components.

    type(loading_options_type), intent(inout) :: loading_options
    type(loading_type), intent(inout) :: loading

    ! Locals:
    ! i/ii: Generic loop index values.
    ! primary_dir: Loading direction value used to index load arrays.
    ! diff_stress: Load target difFErential calculated between steps.

    integer :: i, primary_dir
    real(rk) :: diff_stress

    ! Notes:
    ! User-defined input (strain rates) should always be positive. Compression
    ! and similar signed deformation control are handled internally by way of
    ! the sign in the targeted load.

    !---------------------------------------------------------------------------

    ! Initialize definition variables.

    loading%num_steps = loading_options%num_steps
    primary_dir = loading_options%loading_direction
    diff_stress = 0.0d0

    ! Allocate the arrays for the load target sequence.

    allocate (loading%target_load(loading%num_steps, 3))
    allocate (loading%step_dt_min(loading%num_steps))
    allocate (loading%step_dt_max(loading%num_steps))
    allocate (loading%step_target_time_incr(loading%num_steps))
    allocate (loading%target_sign(loading%num_steps))
    allocate (loading%vel_factor(loading%num_steps))
    allocate (loading%step_print(loading%num_steps))

    allocate (loading%step_velocity(loading%num_steps))

    ! Initialize arrays

    ! FIXME disp
    loading%step_velocity = loading_options%strain_rate

    ! Process all load steps and build control arrays.

    do i = 1, loading%num_steps
      ! Handle the first step alone to avoid out of bounds on arrays.

      loading%target_load(i, 1) = loading_options%target_load(i, 1)
      loading%target_load(i, 2) = loading_options%target_load(i, 2)
      loading%target_load(i, 3) = loading_options%target_load(i, 3)
      loading%step_dt_max(i) = loading_options%step_dt_max(i)
      loading%step_dt_min(i) = loading_options%step_dt_min(i)
      loading%step_print(i) = loading%step_print(i)

      if (i .eq. 1) then
        loading%vel_factor(i) = 1.0d0
      else
        loading%vel_factor(i) = loading%step_velocity(i)/loading%step_velocity(i - 1)
      end if

      if (i .eq. 1) then
        diff_stress = loading_options%target_load(i, primary_dir)
      else
        diff_stress = loading_options%target_load(i, primary_dir) &
                  & - loading_options%target_load(i - 1, primary_dir)
      end if

      if (diff_stress .eq. 0.0d0) then
        call par_quit('Error  :     > Load differential is zero between steps.')
      end if

      if (diff_stress .gt. 0.0d0) then
        loading%target_sign(i) = 1

      else if (diff_stress .lt. 0.0d0) then
        loading%target_sign(i) = -1
      end if

      if (i .eq. 1) then
        if (diff_stress .lt. 0.0d0) then
          loading%bcs_vel = loading%bcs_vel*((-1)*(loading%vel_factor(i)))
        end if
      else
        if (loading%target_sign(i) .ne. loading%target_sign(i - 1)) then
          loading%vel_factor(i) = -loading%vel_factor(i)
        end if
      end if
    end do

    return

  end subroutine process_ctrl_data_csr

  ! ============================================================================

end module boundary_conditions_mod2_legacy
