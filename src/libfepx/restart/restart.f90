! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module restart_mod
  
  ! Module for printing restart variables to file at the end of a step and reading
  ! said files.
  !
  ! Contains subroutines:
  !
  ! Reading restart files from the various drivers:
  ! read_restart_field: Reads field data for restarting simulation
  ! read_uniaxial_restart: Read uniaxial control restart information.
  !
  ! Writing restart files from the various drivers:
  ! write_restart_field: Writes field data for restarting a simulation.
  ! write_uniaxial_restart: Writes required information for the uniaxial restart
  ! write_triaxcsr_restart: Writes required information for the uniaxial restart
  ! write_triaxclr_restart: Writes required information for the uniaxial restart

  use general_mod
  use types_mod
  use parallel_mod

  implicit none

  public

contains

  !> Reads field data for restarting simulation
  subroutine read_restart_field(mesh, results)

    type(mesh_type), intent(out) :: mesh
    type(results_type), intent(out) :: results

    ! Locals:

    integer :: myunit
    character(len=8) :: charid ! assumes less than 10,000 processes
    integer :: rst_num
    logical :: file_exists
    character(len=8) :: rst_num_str
    character(len=64) :: filename

    intrinsic :: trim

    !----------------------------------------------------------------------

    write (charid, '(i0)') myid + 1

    rst_num = 1000
    file_exists = .false.

    ! Find max value n of rstN.control

    do while (.not. file_exists)
      write (rst_num_str, '(i0)') rst_num

      filename = 'rst'//trim(rst_num_str)//'.control'

      inquire (file=filename, exist=file_exists)
      rst_num = rst_num - 1

      if (rst_num .eq. -2) then
        call par_quit('Error  :     > Restart control file not found.')
      end if
    end do

    rst_num = rst_num + 1
    write (rst_num_str, '(i0)') rst_num

    filename = 'rst'//trim(rst_num_str)//'.field.core'//trim(charid)

    open (newunit=myunit, file=filename, form='unformatted', action='read')

    ! Velocity and coordinates.

    read (myunit) results%coo
    read (myunit) results%vel

    ! Orientations, weights and hardnesses.

    read (myunit) mesh%ori
    read (myunit) results%ori(:, :, :, cqpt)
    read (myunit) results%rstar(:, :, :, cqpt)
    read (myunit) results%crss(:, :, cqpt)

    ! Elastic Strains.

    read (myunit) results%e_elas_kk_bar
    read (myunit) results%sig_vec

    ! Equivalent Strains.

    read (myunit) results%slip

    ! Total work, plastic work, and rates.

    read (myunit) results%workrate
    read (myunit) results%workrate_pl

    ! Other integrated quantities.

    read (myunit) results%strain_pl
    read (myunit) results%strain

    close (myunit)

  end subroutine read_restart_field

  subroutine read_uniaxial_restart(loading, mesh, results, &
             & incr, time, load, area, area0)

    ! Read uniaxial control restart information.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! incr: Current step increment (additive).
    ! time: Current step total time.
    ! load: Current load on all surfaces in the mesh.
    ! loading%prev_load: Previous step load value.
    ! area/area0: Current and initial surface areas of the mesh at current step.

    type(mesh_type), intent(in) :: mesh
    type(loading_type), intent(inout) :: loading
    type(results_type) :: results
    integer, intent(out) :: incr
    real(rk), intent(out) :: time
    real(rk), intent(out) :: load(mesh%num_fasets, 3)
    real(rk), intent(out) :: area(mesh%num_fasets)
    real(rk), intent(out) :: area0(mesh%num_fasets)

    ! Locals:
    ! myunit: Current unit number to open restart file.
    ! isurf: Generic loop index to loop over mesh surfaces.

    integer :: myunit
    integer :: isurf
    integer :: rst_num
    logical :: file_exists
    character(len=8) :: rst_num_str
    character(len=64) :: filename
    real(rk) :: diff_disp
    real(rk) :: diff_stress
    integer :: i

    ! Notes:
    ! Current_load becomes the loading%prev_load inside subroutine time_increment.

    !---------------------------------------------------------------------------

    rst_num = 1000
    file_exists = .false.

    do while (.not. file_exists)
      write (rst_num_str, '(i0)') rst_num

      filename = 'rst'//trim(rst_num_str)//'.control'

      inquire (file=filename, exist=file_exists)
      rst_num = rst_num - 1

      if (rst_num .eq. -2) then
        call par_quit('Error  :     > Restart control file not found.')
      end if
    end do

    rst_num = rst_num + 1
    write (rst_num_str, '(i0)') rst_num
    filename = 'rst'//trim(rst_num_str)//'.control'

    open (newunit=myunit, file=filename, form='unformatted', action='read')

    read (myunit) loading%curr_step
    read (myunit) loading%prev_load
    read (myunit) loading%step_complete
    read (myunit) loading%prev_dt
    read (myunit) incr
    read (myunit) time

    do isurf = 1, mesh%num_fasets
      read (myunit) load(isurf, :)
    end do

    read (myunit) area
    read (myunit) area0
    read (myunit) loading%prev_disp
    read (myunit) loading%curr_disp

    ! Re-assign variables for new restart

    loading%step_complete = .true.
    loading%curr_load(1) = load(2, 1)
    loading%curr_load(2) = load(4, 2)
    loading%curr_load(3) = load(6, 3)

    ! Print statistics to console

    if (myid .eq. 0) then
      write (*, '(a)') 'Info   : Reading restart control information...'
      write (*, '(a)') 'Info   :   - Previous simulation ended with final:'
      write (*, '(a, i0)') 'Info   :     > Increments: ', incr
      write (*, '(a, i0)') 'Info   :     > Steps:      ', loading%curr_step - 1
      write (*, '(a, e14.4)') 'Info   :     > Time:  ', time
      write (*, '(a, 3(e14.4))') 'Info   :     > Normal loads:  ', &
          & loading%curr_load
    end if

    ! Reinitialize step and time if new files
    ! In the future, we will have to handle the append case, which will continue
    ! sequentially, rather than reinitializing
    loading%curr_step = 1
    loading%step_complete = .false.
    incr = 0
    time = 0.0d0

    ! Reevaluate time step for first step if strain targeting

    if (loading%def_control_by .eq. 'uniaxial_strain_target') then
      ! FIXME disp
      loading%step_dt_max(1) = 1.0d0
    end if

    ! Check for a reversal in loading direction upon restart, correct loading
    ! history if so. The loading history as assigned in process_ctrl_data
    ! assumes that we are starting from 0 load. Below rectifies any directional
    ! issues with loading due to this assumption.

    if ((loading%def_control_by .eq. 'uniaxial_load_target') .and. &
        & (((loading%curr_load(loading%loading_direction) - &
        & loading%prev_load(loading%loading_direction))* &
        & (loading%target_load(1, 1) -  &
        & loading%curr_load(loading%loading_direction))) .lt. 0.0d0)) then
      ! Change vel direction

      results%vel = -1.0d0*results%vel

      ! And check for any load reversals on subsequent steps

      if (loading%num_steps .ge. 2) then
        ! Calculate loading%target_sign for the first step

        diff_stress = loading%target_load(1, 1) - &
            & loading%curr_load(loading%loading_direction)

        if (diff_stress .gt. 0.0d0) then
          loading%target_sign(1) = 1

        else if (diff_stress .lt. 0.0d0) then
          loading%target_sign(1) = -1

        else ! The differential change is zero, thus dt is zero.
          call par_quit('Error  :     > Load differential is zero&
              & between steps.')
        end if

        ! Reset loading%vel_factor for reassignment in the below loops

        loading%vel_factor = abs(loading%vel_factor)

        ! Calculate loading%vel_factor for each step after the first

        do i = 2, loading%num_steps
          diff_stress = (loading%target_load(i, 1)) &
              & - (loading%target_load(i - 1, 1))

          if (diff_stress .gt. 0.0d0) then
            loading%target_sign(i) = 1

          else if (diff_stress .lt. 0.0d0) then
            loading%target_sign(i) = -1

          else ! The differential change is zero, thus dt is zero.
            call par_quit('Error  :     > Load differential is zero&
                & between steps.')
          end if

          ! Check if a load reversal occured since the last step.

          if (loading%target_sign(i) .ne. loading%target_sign(i - 1)) then
            loading%vel_factor(i) = loading%vel_factor(i)*(-1)

          else if (loading%target_sign(i) .eq. loading%target_sign(i - 1)) then
            loading%vel_factor(i) = loading%vel_factor(i)*(1)

          else
            call par_quit('Error  :     > Load reversal check failed.')
          end if
        end do
      end if

    else if ((loading%def_control_by .eq. 'uniaxial_strain_target') .and. &
        & (((loading%curr_disp - loading%prev_disp)* &
        & (loading%target_disp(1) - loading%curr_disp)) &
        & .lt. 0.0d0)) then
      ! Change vel direction

      results%vel = -1.0d0*results%vel

      ! And check for any load reversals on subsequent steps

      if (loading%num_steps .ge. 2) then
        ! FIXME disp
        diff_disp = loading%target_disp(1) - loading%curr_disp

        if (diff_disp .gt. 0.0d0) then
          loading%target_sign(1) = 1

        else if (diff_disp .lt. 0.0d0) then
          loading%target_sign(1) = -1

        else ! The differential change is zero, thus dt is zero.
          call par_quit('Error  :     > Strain differential is zero between steps.')
        end if

        ! Reset loading%vel_factor for reassignment in the below loops

        loading%vel_factor = abs(loading%vel_factor)

        do i = 2, loading%num_steps
        ! FIXME disp
          diff_disp = loading%target_disp(i) - loading%target_disp(i - 1)

          if (diff_disp .gt. 0.0d0) then
            loading%target_sign(i) = 1

          else if (diff_disp .lt. 0.0d0) then
            loading%target_sign(i) = -1

          else ! The differential change is zero, thus dt is zero.
            call par_quit('Error  :     > Load differential is zero&
                & between steps.')
          end if

          if (loading%target_sign(i) .ne. loading%target_sign(i - 1)) then
            loading%vel_factor(i) = loading%vel_factor(i)*(-1)

          else if (loading%target_sign(i) .eq. loading%target_sign(i - 1)) then
            loading%vel_factor(i) = loading%vel_factor(i)*(1)

          else
            call par_quit('Error  :     > Load reversal check failed.')
          end if
        end do
      end if
    end if

    close (myunit)

    return

  end subroutine read_uniaxial_restart

  subroutine write_restart_field(istep, mesh, results)

    !  Writes field data for restarting simulation

    !---------------------------------------------------------------------------

    ! Arguments:

    integer, intent(in) :: istep
    type(mesh_type), intent(inout) :: mesh
    type(results_type), intent(in) :: results

    ! Locals:

    integer :: myunit
    character(len=8) :: charid ! assumes less than 10,000 processes
    integer :: rst_num
    logical :: file_exists
    character(len=8) :: rst_num_str
    character(len=64) :: filename

    intrinsic :: trim

    !---------------------------------------------------------------------------

    write (charid, '(i0)') myid + 1

    rst_num = 1000
    file_exists = .false.

    do while (.not. file_exists)
      write (rst_num_str, '(i0)') rst_num

      filename = 'rst'//trim(rst_num_str)//'.control'

      inquire (file=filename, exist=file_exists)
      rst_num = rst_num - 1

      ! If we are on the first simulation, non-restart

      if (rst_num .eq. -2) then
        file_exists = .true.
        rst_num = -2
      end if
    end do

    if (istep .eq. 1) then
      rst_num = rst_num + 2
      write (rst_num_str, '(i0)') rst_num
    end if

    filename = 'rst'//trim(rst_num_str)//'.field.core'//trim(charid)
    open (newunit=myunit, file=filename, form='unformatted', action='write')

    ! Velocity and coordinates.

    write (myunit) results%coo
    write (myunit) results%vel

    ! Orientations, weights and hardnesses.

    write (myunit) mesh%ori
    write (myunit) results%ori(:, :, :, cqpt)
    write (myunit) results%rstar(:, :, :, cqpt)
    write (myunit) results%crss(:, :, cqpt)

    ! Elastic Strains.

    write (myunit) results%e_elas_kk_bar
    write (myunit) results%sig_vec

    ! Equivalent Strains

    write (myunit) results%slip

    ! Total work, plastic work, and rates.

    write (myunit) results%workrate
    write (myunit) results%workrate_pl

    ! Other integrated quantities

    write (myunit) results%strain_pl
    write (myunit) results%strain

    close (myunit)

  end subroutine write_restart_field

  !> write uniaxial control restart information
  subroutine write_uniaxial_restart(mesh, loading, incr, time, load, area, area0)

    ! Arguments:
    ! load: Current load on all surfaces in the mesh.
    ! area/area0: Current and initial surface areas of the mesh at current step.

    type(mesh_type), intent(in) :: mesh
    type(loading_type), intent(in) :: loading
    integer, intent(in) :: incr
    real(rk), intent(in) :: time
    real(rk), intent(in) :: load(mesh%num_fasets, 3)
    real(rk), intent(in) :: area(mesh%num_fasets)
    real(rk), intent(in) :: area0(mesh%num_fasets)

    integer :: myunit
    integer :: rst_num
    logical :: file_exists
    character(len=8) :: rst_num_str
    character(len=64) :: filename
    integer :: isurf

    !---------------------------------------------------------------------------

    if (myid .eq. 0) then

      rst_num = 1000
      file_exists = .false.

      do while (.not. file_exists)
        write (rst_num_str, '(i0)') rst_num

        filename = 'rst'//trim(rst_num_str)//'.control'

        inquire (file=filename, exist=file_exists)
        rst_num = rst_num - 1

        ! If we are on the first simulation, non-restart

        if (rst_num .eq. -2) then
          file_exists = .true.
          rst_num = -2
        end if
      end do

      if (loading%curr_step - 1 .eq. 1) then
        rst_num = rst_num + 2
        write (rst_num_str, '(i0)') rst_num
      end if

      filename = 'rst'//trim(rst_num_str)//'.control'
      open (newunit=myunit, file=filename, form='unformatted', action='write')

      write (myunit) loading%curr_step
      write (myunit) loading%prev_load
      write (myunit) loading%step_complete
      write (myunit) loading%prev_dt
      write (myunit) incr
      write (myunit) time

      do isurf = 1, mesh%num_fasets
        write (myunit) load(isurf, :)
      end do

      write (myunit) area
      write (myunit) area0
      write (myunit) loading%prev_disp
      write (myunit) loading%curr_disp

      close (myunit)
    end if

    return

  end subroutine write_uniaxial_restart

  !===========================================================================

  subroutine write_triaxcsr_restart(mesh, istep, curr_load, prev_load, step_complete, &
      & dtime, incr, time, surf_load_array, area, area0, length, length0, curr_vel, &
      & s_pert_mag, t_pert_mag)

    ! Write uniaxial control restart information.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! Needs to be defined - jc

    type(mesh_type), intent(in) :: mesh
    logical, intent(in)  :: step_complete
    integer, intent(in)  :: istep
    integer, intent(in)  :: incr
    real(rk), intent(in) :: curr_load(3), prev_load(3)
    real(rk), intent(in) :: dtime, time
    real(rk), intent(in) :: surf_load_array(mesh%num_fasets, 3)
    real(rk), intent(in) :: area(mesh%num_fasets), area0(mesh%num_fasets)
    real(rk), intent(in) :: length(3), length0(3)
    real(rk), intent(in) :: curr_vel(3)
    real(rk), intent(in) :: s_pert_mag, t_pert_mag

    ! Locals:
    ! myunit: Current unit number to open restart file.
    ! isurf: Generic loop index to loop over mesh surfaces.

    integer :: myunit
    integer :: rst_num
    logical :: file_exists
    character(len=8) :: rst_num_str
    character(len=64) :: filename
    integer :: isurf

    !---------------------------------------------------------------------------

    if (myid .eq. 0) then

      rst_num = 1000
      file_exists = .false.

      do while (.not. file_exists)
        write (rst_num_str, '(i0)') rst_num

        filename = 'rst'//trim(rst_num_str)//'.control'

        inquire (file=filename, exist=file_exists)
        rst_num = rst_num - 1

        ! If we are on the first simulation, non-restart

        if (rst_num .eq. -2) then
          file_exists = .true.
          rst_num = -2
        end if
      end do

      if (istep - 1 .eq. 1) then
        rst_num = rst_num + 2
        write (rst_num_str, '(i0)') rst_num
      end if

      filename = 'rst'//trim(rst_num_str)//'.control'
      open (newunit=myunit, file=filename, form='unformatted', action='write')

      write (myunit) istep
      write (myunit) curr_load
      write (myunit) prev_load
      write (myunit) step_complete
      write (myunit) dtime
      write (myunit) incr
      write (myunit) time

      do isurf = 1, mesh%num_fasets
        write (myunit) surf_load_array(isurf, :)
      end do

      write (myunit) area
      write (myunit) area0
      write (myunit) length
      write (myunit) length0
      write (myunit) curr_vel
      write (myunit) s_pert_mag
      write (myunit) t_pert_mag

      close (myunit)
    end if

    return

  end subroutine write_triaxcsr_restart

  !===========================================================================

  subroutine write_triaxclr_restart(mesh, istep, curr_load, prev_load, &
      & first_incr_in_step, incr, time, surf_load_array, area, area0, &
      & length, length0, curr_vel, prev_action, curr_action, &
      & initial_load_dwell_vel, initial_unload_dwell_vel)

    !  Write TriaxCLR restart information.

    !---------------------------------------------------------------------------

    ! Arguments:

    type(mesh_type), intent(in) :: mesh
    logical, intent(in) :: first_incr_in_step
    integer, intent(in) :: istep
    integer, intent(in) :: incr
    character(len=10), intent(in) :: prev_action, curr_action
    real(rk), intent(in) :: curr_load(3)
    real(rk), intent(in) :: prev_load(3)
    real(rk), intent(in) :: time
    real(rk), intent(in) :: surf_load_array(mesh%num_fasets, 3)
    real(rk), intent(in) :: area(mesh%num_fasets)
    real(rk), intent(in) :: area0(mesh%num_fasets)
    real(rk), intent(in) :: length(3), length0(3)
    real(rk), intent(in) :: curr_vel(3)
    real(rk), intent(in) :: initial_load_dwell_vel(3)
    real(rk), intent(in) :: initial_unload_dwell_vel(3)

    ! Locals:

    integer :: myunit
    integer :: rst_num
    logical :: file_exists
    character(len=8) :: rst_num_str
    character(len=64) :: filename
    integer :: isurf

    !---------------------------------------------------------------------------

    if (myid .eq. 0) then

      rst_num = 1000
      file_exists = .false.

      do while (.not. file_exists)
        write (rst_num_str, '(i0)') rst_num

        filename = 'rst'//trim(rst_num_str)//'.control'

        inquire (file=filename, exist=file_exists)
        rst_num = rst_num - 1

        ! If we are on the first simulation, non-restart

        if (rst_num .eq. -2) then
          file_exists = .true.
          rst_num = -2
        end if
      end do

      if (istep .eq. 1) then
        rst_num = rst_num + 2
        write (rst_num_str, '(i0)') rst_num
      end if

      filename = 'rst'//trim(rst_num_str)//'.control'
      open (newunit=myunit, file=filename, form='unformatted', action='write')

      write (myunit) istep
      write (myunit) curr_load
      write (myunit) prev_load
      write (myunit) first_incr_in_step
      write (myunit) incr
      write (myunit) time

      do isurf = 1, mesh%num_fasets
        write (myunit) surf_load_array(isurf, :)
      end do
      write (myunit) area
      write (myunit) area0
      write (myunit) length
      write (myunit) length0
      write (myunit) curr_vel
      write (myunit) prev_action
      write (myunit) curr_action
      write (myunit) initial_load_dwell_vel
      write (myunit) initial_unload_dwell_vel

      close (myunit)
    end if

    return

  end subroutine write_triaxclr_restart

end module restart_mod
