! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module driver_uniaxial_control_mod

! Driver for uniaxial strain control with either load or strain targets.

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use utils_mod
  use loading_type_mod
  use mesh_type_mod
  use res_init_mod
  use restart_mod
  use solveit_isovp_mod
  use driver_uniaxial_control_mod2
  use finalize_res_mod
  use solveit_evp_mod
  use kinematics_mod
  use matrix_operations_mod, only: strain_equiv_3x3, vec6_mat_symm
  use quadrature_mod, only: nqpt
  use restart_mod, only: read_restart_field
  use write_res_mod
  use gather_scatter_mod, only: trace, part_gather
  use parallel_mod, only: par_message, par_quit

  implicit none

  public

contains

  !> Primary driver for uniaxial strain control with load or strain targeting.
  subroutine driver_uniaxial_control(mesh, crys, loading, exec, results, printing)

    type(mesh_type), intent(inout) :: mesh
    type(crys_type) :: crys(:)
    type(loading_type), intent(inout) :: loading
    type(exec_type), intent(inout) :: exec
    type(results_type), intent(inout) :: results
    type(printing_type), intent(in) :: printing

    ! Locals:

    type(results_type) :: results_prev
    integer :: incr
    real(rk) :: dtime, time
    real(rk) :: load(mesh%num_fasets, 3)
    real(rk) :: area0(mesh%num_fasets), area(mesh%num_fasets)
    integer  :: ntsteps
    integer  :: restart_incr
    real(rk) :: clock_end
    real(rk) :: buff_vel(dof_sub:dof_sup)
    real(rk) :: vel_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: buff_primary_drive_vel(dof_sub:dof_sup)
    real(rk) :: primary_drive_vel_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: buff_secondary_drive_vel(dof_sub:dof_sup)
    real(rk) :: secondary_drive_vel_ebe(kdim, elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    call res_qpt_alloc (mesh, results_prev)
    call res_qpt_init_prev (results, results_prev)

    ! Initialization -----------------------------------------------------------

    ntsteps = 0
    buff_vel = 0.0d0
    vel_ebe = 0.0d0
    buff_primary_drive_vel = 0.0d0
    primary_drive_vel_ebe = 0.0d0
    buff_secondary_drive_vel = 0.0d0
    secondary_drive_vel_ebe = 0.0d0
    load = 0.0d0
    
    ! Standard simulation (no restart)
    if (.not. exec%restart) then
      ! Running isovp solution
      call write_res(0, mesh, crys, results, 0.0d0, printing)

      if (myid .eq. 0) write (*, '(a)') 'Info   :   - Initializing fields &
          &from isotropic viscoplastic solution'

      call solveit_vp(mesh, loading, exec, results)

      ! Correction of results%vel to insure right offset giving mpcs
      if ((loading%mpc_status .eqv. .true.) .and. (mesh%num_periodicity .eq. 0)) then

       call part_gather(vel_ebe, results%vel, int(loading%conn_mpc), loading%mpc_trace)
       buff_vel = 0.0d0
       call part_scatter(buff_vel, vel_ebe, mesh%elt_dofs, exec%dof_trace)
       buff_vel=buff_vel/mesh%g_ones

       results%vel = loading%coeff_ps*buff_vel+loading%offset_ps

     else if ((loading%mpc_status .eqv. .false.) .and. (mesh%num_periodicity .gt. 0)) then

       call part_gather(vel_ebe, results%vel, int(loading%conn_mpc), loading%mpc_trace)
       buff_vel = 0.0d0
       call part_scatter(buff_vel, vel_ebe, mesh%elt_dofs, exec%dof_trace)
       buff_vel=buff_vel/mesh%g_ones
       
       call part_gather(primary_drive_vel_ebe, results%vel, int(loading%primary_drive_ebe), loading%mpc_trace)
       buff_primary_drive_vel = 0.0d0
       call part_scatter(buff_primary_drive_vel, primary_drive_vel_ebe, mesh%elt_dofs, exec%dof_trace)
       buff_primary_drive_vel=buff_primary_drive_vel/mesh%g_ones
       
       call part_gather(secondary_drive_vel_ebe, results%vel, int(loading%secondary_drive_ebe), loading%mpc_trace)
       buff_secondary_drive_vel = 0.0d0
       call part_scatter(buff_secondary_drive_vel, secondary_drive_vel_ebe, mesh%elt_dofs, exec%dof_trace)
       buff_secondary_drive_vel=buff_secondary_drive_vel/mesh%g_ones

       results%vel = loading%coeff_ps*(buff_vel + loading%label_sg*(1 - loading%imposed_state)*&
                    &(buff_primary_drive_vel + buff_secondary_drive_vel))+ loading%label_sg*loading%offset_ps

      end if

      ! Compute initial area (area0)
      call mesh_surfaceareas(mesh, areas=area0)

      ! Initialize deformation control
      incr = 0
      time = 0.0d0
      dtime = loading%step_dt_max(1)

    ! Restart simulation
    else
      call read_restart_field(mesh, results)

      call read_uniaxial_restart(loading, mesh, results, incr, time, load, area, &
          & area0)

      restart_incr = incr

      ! Update initial guess of velocity field for load reversal
      ! loading%step_complete resets after call to loading_timeincr
      if (loading%step_complete) then
        results%vel = loading%vel_factor(loading%curr_step)*results%vel
      end if

      call loading_timeincr(mesh, loading, exec, load, time, dtime)
    end if

    ! Writing initial values to force and conv files ---------------------------

    if (myid .eq. 0) then

      ! Write the header to file post.force#
      if (printing%print_forces) then
        call write_force_file_headers(mesh, printing, 1)

        ! If virgin sample, write 0th step
        if (.not. exec%restart) then
          call write_force_file_data(mesh, printing, 0, 0, load, area0, 0.0d0)
        end if
      end if

      if (printing%print_conv) then
        call write_conv_file_headers(printing)
      end if
    end if

    ! Starting evp simulation -------------------------------------------------

    ! Time stepping loop: increment to find new configuration and material state
    time_stepping: do

      ! Increment
      incr = incr + 1
      ! Update the total time
      time = time + dtime

      ! If previous step just complete, starting new step
      if (loading%step_complete) then
        loading%curr_step = loading%curr_step + 1
      end if

      call res_qpt_init_prev (results, results_prev)

      ! Update initial guess of velocity field for load reversal
      if (loading%step_complete) then
        results%vel = loading%vel_factor(loading%curr_step)*results%vel
        ntsteps = 0
      end if

      if ((loading%vel_factor(loading%curr_step) .eq. -1) .and. (ntsteps .eq. 0)) then
        results%acmslip = 0.0d0
      end if

      ntsteps = ntsteps + 1

      call driver_uniaxial_control_printtoterminal (loading, printing, loading%curr_step, &
                                                  & time, dtime, restart_incr, incr)

      ! Compute: vel @(t+dt) using pcg solver
      call solveit_evp(mesh, crys, loading, exec, printing, results_prev, &
                     & results, dtime, incr, load, area)

      ! Check if the target load is reached and find new dtime
      call loading_timeincr(mesh, loading, exec, load, time, dtime)

      ! Writing results --------------------------------------------------------

      ! Write time increment results
      if (myid .eq. 0) then
        if (printing%print_forces) then
          call write_force_file_data(mesh, printing, loading%curr_step, incr, load, area, time)
        end if
      end if

      ! Write step results
      if (loading%step_complete) then

        ! write results per se
        if (loading%step_print(loading%curr_step)) then
          !
          call write_res(loading%curr_step, mesh, crys, results, dtime, printing)
          ! writing restart
          if (printing%print_restart) then
            call write_restart_field(loading%curr_step, mesh, results)

            if (loading%def_control_by .eq. 'uniaxial_strain_target') then
              if (loading%curr_step .eq. 1) then
                loading%prev_disp = 0.0d0

              else
                loading%prev_disp = loading%target_disp(loading%curr_step - 1)
              end if

              loading%curr_disp = loading%target_disp(loading%curr_step)

            else if (loading%def_control_by .eq. 'uniaxial_load_target') then
              loading%prev_disp = 0.0d0
              loading%curr_disp = 0.0d0
            end if

            call write_uniaxial_restart(mesh, loading, incr, time, load, area, area0)
          end if

        end if

        ! checking for limits
        if (loading_checklimit(exec, time, incr)) then
          loading%step_complete = .true.
          loading%all_steps_complete = .true.
          call write_dot_sim_file_complete_steps(printing, loading)
          call par_quit('Error  :     > Maximum time or maximum increments exceeded.')
        end if

        ! optionally, checking for necking
        if (exec%check_necking .and. loading_isnecking(loading)) then
          call write_dot_sim_file_complete_steps(printing, loading)
          call par_quit('Error  :     > Specimen is necking.')
        end if
      end if

      if (loading%all_steps_complete) then
        call write_dot_sim_file_complete_steps(printing, loading)
        ! Finalize clock values and print to console

        if (myid .eq. 0) then
          call cpu_time(clock_end)
          write (*, '(a, f10.3, a)') 'Info   : Elapsed time:', &
              & clock_end - exec%clock_start, ' secs.'
        end if

        call par_quit('Info   : Final step terminated. Simulation completed successfully.')
      end if

    end do time_stepping

    return

  end subroutine driver_uniaxial_control

end module driver_uniaxial_control_mod
