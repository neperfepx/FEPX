! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module loading_type_mod

! Module to define types

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use parallel_mod
  use types_mod, only: mesh_type, loading_type, exec_type
  use loading_type_mod2

  implicit none

  public

  contains

  !> Return time increment according to loading step
  subroutine loading_timeincr(mesh, loading, exec, load, time, dtime)

    implicit none

    ! Arguments:
    ! dtime: time increment
    ! load: Current step load value.

    type(mesh_type), intent(in) :: mesh
    type(loading_type), intent(inout) :: loading
    type(exec_type), intent(in) :: exec
    real(rk), intent(in) :: load(mesh%num_fasets, 3)
    real(rk), intent(in) :: time
    real(rk), intent(out) :: dtime

    ! Locals:
    ! index: Generic loop index value.
    ! dload: Differential load between current and previous time steps.
    ! trial_load: Estimated load value to determine if dt needs to be changed.
    ! trial_dt: Estimated new time step value for the trial load.
    ! min_time_step: Current step minimum defined time INCRement.
    ! delta_load: DifFErence between the trial load and current load at step.

    integer :: index
    real(rk) :: dload(3), trial_load(3), trial_dt, min_time_step, delta_load

    ! Notes:
    ! Used for load control only
    ! The suggested time step is given in the input data.
    ! Near the target load, we use a linear approximation to match the load,
    ! with a minimum value to insure we get there.

    !---------------------------------------------------------------------------

    dtime = loading%step_dt_max(loading%curr_step)

    ! First, check if we are at target increment or load

    select case (loading%def_control_by)

    case ('uniaxial_strain_target')
      loading%step_complete = .false.

      ! FIXME 1e-6?
      if (time + 1e-6 .gt. loading%target_time(loading%curr_step)) then
        loading%step_complete = .true.
      end if

    case ('uniaxial_load_target')
      loading%prev_load = loading%curr_load
      loading%curr_load = load(loading%loading_face, :)

      loading%step_complete = loading_load_isinrange(loading, exec, loading%curr_load)

    case default
      call par_quit('Error  :     > Invalid control option.')

    end select

    if (loading%step_complete) then
      if (loading%curr_step .eq. loading%num_steps) then
        loading%all_steps_complete = .true.

        return
      end if

      if (.not. loading%step_complete) then
        dtime = loading%step_dt_max(loading%curr_step)
      else if (loading%curr_step .lt. loading%num_steps) then
        dtime = loading%step_dt_max(loading%curr_step + 1)
      else
        dtime = -1.0d0
      end if
    end if

    ! Now, check to see if we are close to target, and need to adjust the time step.
    if (loading%def_control_by .eq. 'uniaxial_load_target') then
      dload = (loading%curr_load - loading%prev_load)/loading%prev_dt
      trial_load = loading%curr_load + dload*dtime

      if (loading_load_isinrange(loading, exec, trial_load)) then
        index = loading%loading_direction

        !legacy note: -tm. from donald's version:
        delta_load = trial_load(index) - loading%curr_load(index)

        if (abs(delta_load) < 1.0d-16) then
          loading%prev_dt = dtime;
          return
        end if

        trial_dt = (loading%target_load(loading%curr_step, 1) - loading%curr_load(index))/&
            & (trial_load(index) - loading%curr_load(index))

        trial_dt = trial_dt*dtime
        dtime = trial_dt*exec%dtime_factor

        if (dtime > loading%step_dt_max(loading%curr_step)) then
          dtime = loading%step_dt_max(loading%curr_step)
        end if

        min_time_step = loading%step_dt_min(loading%curr_step)

        if ((dtime) < min_time_step) then
          dtime = min_time_step
        end if
      end if
    end if

    loading%prev_dt = dtime

  end subroutine loading_timeincr

  !> Check if specimen time or increment limit is tripped.
  function loading_checklimit(exec, time, incr) result(is_limit_tripped)

    implicit none

    ! Arguments:
    ! time: Current step total time.
    ! incr: Current step increment (additive).
    ! is_limit_tripped: Flag for checking if simulation controls are tripped.

    type(exec_type), intent(in) :: exec
    real(rk), intent(in) :: time
    integer, intent(in)  :: incr

    logical :: is_limit_tripped

    !---------------------------------------------------------------------------

    is_limit_tripped = .false.

    if ((time .ge. exec%max_total_time) .or. &
        & (incr .ge. exec%max_incr)) then
      is_limit_tripped = .true.

    end if

    return

    !---------------------------------------------------------------------------

  end function loading_checklimit

  function loading_isnecking(loading) result(is_necking)

    ! Check if specimen is necking (that is, the load on the target surface is
    ! decreasing from the previous step).

    !---------------------------------------------------------------------------

    implicit none

    ! Arguments:
    ! is_necking: Flag for checking if the specimen has necked.

    type(loading_type), intent(inout) :: loading
    logical :: is_necking

    ! Locals:
    ! index: Current step loading direction.

    integer :: index

    !---------------------------------------------------------------------------

    is_necking = .false.

    ! FIXME
    ! quick fix for final-state testing (or should be use curr_step - 1?)
    if (loading%curr_step .gt. loading%num_steps) return

    index = loading%loading_direction

    if ((.not. loading%step_complete) .and. &
        & ((loading%curr_load(index) - loading%prev_load(index))/loading%target_sign(loading%curr_step)&
        & .lt. 0.0)) then
      is_necking = .true.
      loading%step_complete = .true.
      loading%all_steps_complete = .true.

      ! Update step because get_print_flag checks flag from previous step.
      loading%curr_step = loading%curr_step + 1
    end if

    return

    !---------------------------------------------------------------------------

  end function loading_isnecking

end module loading_type_mod
