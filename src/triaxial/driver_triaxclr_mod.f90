! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module driver_triaxclr_mod

! Driver for evp simulation using triaxial load control at constant load rate

! Contains subroutines:
! driver_triaxclr: Driver for triaxial clr simulation
! print_headers: Print headers to output files
! read_triaxclr_restart: Read restart files for triaxial clr simulations.

  use general_mod
  use types_mod
  use mesh_type_mod
  use utils_mod
  use res_init_mod
  use, intrinsic :: iso_fortran_env, only: rk => real64
  use finalize_res_mod
  use solveit_evp_mod2
  use driver_triax_utilities_mod
  use kinematics_mod
  use matrix_operations_mod, only: solve_lin_sys_3, &
      & strain_equiv_3x3, vec6_mat_symm
  use surface_mod
  use units_mod
  use write_res_mod
  use restart_mod
  use gather_scatter_mod
  use parallel_mod

  implicit none

  public

contains

  subroutine driver_triaxclr(mesh, crys, loading, exec, &
             & results, printing)

    !---------------------------------------------------------------------------

    ! Driver for triaxial clr simulations

    ! Arguments:

    type(mesh_type), intent(inout) :: mesh
    type(crys_type) :: crys (:)
    type(loading_type) :: loading
    type(exec_type), intent(inout) :: exec
    type(results_type), intent(inout) :: results
    type(printing_type), intent(inout) :: printing

    ! fraction of loading%step_target_time_incr to use on first dwell increement
    real, parameter :: initial_time_frac = 0.01
    ! tolerance on dwell time expressed as a fraction of loading%step_target_time_incr
    real, parameter :: time_tol = 0.001
    integer :: incr
    real(rk) :: dtime, dtime_step, time
    real(rk) :: bak_vel(dof_sub:dof_sup)
    integer :: m_elt, i
    integer :: istep
    integer :: idir
    integer :: indx((dof_sup - dof_sub + 1)/3)
    integer :: indy((dof_sup - dof_sub + 1)/3)
    integer :: indz((dof_sup - dof_sub + 1)/3)
    real(rk) :: load(mesh%num_fasets, 3)
    real(rk) :: area0(mesh%num_fasets), area(mesh%num_fasets)
    real(rk) :: length(3), length0(3)
    real(rk) :: macro_eng_strain(3)
    logical :: first_incr_in_step
    logical :: start_reload, start_unload, start_dwell
    character(len=10) :: prev_action, curr_action
    integer :: nincr_step, incr_count
    integer :: pert_dir
    integer :: bc_iter_1, bc_iter_2
    real(rk) :: step_frac
    real(rk) :: initial_load(3), prev_load(3), curr_load(3)
    real(rk) :: target_load(3), trial_load(3)
    real(rk) :: i_end_load, i1_end_load
    real(rk) :: dwell_time_remaining
    real(rk) :: initial_vel(3), prev_vel(3), curr_vel(3)
    real(rk) :: initial_load_dwell_vel(3), initial_unload_dwell_vel(3)
    real(rk) :: delta_vel(3)
    real(rk) :: vel_scale_factor
    real(rk) :: pert_mag, pert_sign
    real(rk) :: a(3, 3)
    real(rk) :: b(3)
    real(rk) :: eavg, nuavg
    real(rk) :: curr_eqstrain
    character(len=16) :: field, time_string, dtime_string

    real(rk) :: clock_end

    type(results_type) :: results_prev

    call res_qpt_alloc (mesh, results_prev)
    call res_qpt_init_prev (results, results_prev)

    !---------------------------------------------------------------------------

    loading%num_steps = size(loading%target_load, 1)

    m_elt = elt_sup - elt_sub + 1

    ! Indices for x, y, and z degrees of freedom

    do i = 1, ((dof_sup - dof_sub + 1)/3)
      indx(i) = dof_sub + 3*(i - 1)
      indy(i) = dof_sub + 3*(i - 1) + 1
      indz(i) = dof_sub + 3*(i - 1) + 2
    end do

    ! Restarting options

    if (exec%restart) then
      call read_restart_field(mesh, results)

      call read_triaxclr_restart(mesh, istep, curr_load, prev_load, &
          & first_incr_in_step, incr, time, load, area, area0, &
          & length, length0, curr_vel, prev_action, curr_action, &
          & initial_load_dwell_vel, initial_unload_dwell_vel)

    else
      ! Initialize areas and load arrays

      load = 0.0d0
      curr_load = 0.0d0
      area = 0.0d0
      area0 = 0.0d0

      ! Compute initial area (area0)

      call mesh_surfaceareas(mesh, results%coo, area0)

      !if (myid .eq. 0) write(*,'(a25,6f12.6)') 'initial areas &
      !    &(area0):   ', area0

      ! Compute initial mesh dimensions

      call calc_mesh_dim(results, length, indx, indy, indz)

      length0 = length

      ! Initialize state

      results%ori = spread(mesh%ori, 4, nqpt)

      ! Initialize integrated quantities

      results%strain_pl = 0.0d0
      results%defrate = 0.0d0
      results%strain = 0.0d0
      results%slip = 0.0d0

      ! Initialize deformation control

      incr = 1
      istep = 1
      time = 0.0d0
      first_incr_in_step = .true.
      prev_action = "dwelling"
      curr_action = "dwelling"
      initial_load_dwell_vel = 0.0d0
      initial_unload_dwell_vel = 0.0d0
      curr_eqstrain = 0.0d0

      ! Print initial values
      call write_res(0, mesh, crys, results, 0.0d0, printing)
    end if

    ! Initialize load control

    start_reload = .false.
    start_unload = .false.
    start_dwell = .false.

    ! Estimate bulk elastic moduli
    ! Assume uniform texture and equal element volumes
    ! Consider first crystal
    call est_avg_mod(mesh, crys(1), eavg, nuavg)

    ! if (myid .eq. 0) write(*,'(2e16.8)') eavg, nuavg

    ! Print headers to debug files - commented out to avoid standard printing
    ! call print_headers()

    ! Debug printing
    ! write(*,*) '1', loading%target_load(1,1), loading%target_load(1,2), loading%target_load(1,3), &
    !    & loading%control_dir(1), loading%load_rate(1), loading%dwell_time(1), loading%step_target_time_incr(1), &
    !    & loading%step_print(1)

    ! Print headers to output files
    if (myid .eq. 0) then
      ! Write the header to file post.conv
      if (printing%print_conv) then
        call write_conv_file_headers(printing)
      end if

      ! Write the header to file post.force1-6
      if (printing%print_forces) then
        call write_force_file_headers(mesh, printing, 2)

        if (printing%restart .eqv. .false.) then
          ! If virgin sample, write 0th step

          call write_force_file_data(mesh, printing, 0, 0, load, area0, 0.0d0, length0)
        end if
      end if
    end if

    ! if (myid .eq. 0) write(*,'(a)') 'Info   : Running step 1...'

    ! Time stepping loop

    time_stepping: do while (.true.)

      ! Initialize gammadots

      call res_qpt_init_prev (results, results_prev)

      ! If first increment in step, calculate time increment, number of
      !   increments in the step, and initialize increment counter. Check for
      !   start of unloading, reloading, and dwell episodes.

      if (first_incr_in_step) then
        if (myid .eq. 0) then
          write (*, '(a,i0,a)') 'Info   : Running step ', istep, '...'
        end if

        ! Get end loads for steps i and i-1

        idir = abs(loading%control_dir(istep))
        i_end_load = loading%target_load(istep, idir)

        if (istep .eq. 1) then
          i1_end_load = 0.0d0

        else
          i1_end_load = loading%target_load(istep - 1, idir)
        end if

        ! Get previous and current actions

        prev_action = curr_action

        if (loading%control_dir(istep) .lt. 0) then
          curr_action = "dwelling"

        else if (i_end_load .gt. i1_end_load) then
          curr_action = "loading"

        else if (i_end_load .lt. i1_end_load) then
          curr_action = "unloading"

        else
          write (*, *) loading%control_dir(istep)
          write (*, *) i_end_load, i1_end_load
          call par_quit('Error  :     > Invalid load sequence')
        end if

        ! Check for start of loading episode

        if ((curr_action .eq. "loading") .and. (prev_action .ne. "loading") &
            & .and. (istep .ne. 1)) then
          start_reload = .true.

          if (myid .eq. 0) then
            write (*, '(a)') 'Info   :   - Starting reload episode'
          end if
        end if

        ! Check for start of unloading episode

        if ((curr_action .eq. "unloading") .and. &
            & (prev_action .ne. "unloading")) then
          start_unload = .true.

          if (myid .eq. 0) then
            write (*, '(a)') 'Info   :   - Starting unload episode'
          end if
        end if

        ! Check for start of dwell episode

        if ((curr_action .eq. "dwelling") .and. (prev_action .ne. "dwelling")) &
            & then
          start_dwell = .true.

          if (myid .eq. 0) then
            write (*, '(a)') 'Info   :   - Starting dwell episode'
          end if
        end if

        ! Calculate time step and time increment

        if (curr_action .eq. "dwelling") then
          dtime_step = loading%dwell_time(istep)
          dwell_time_remaining = dtime_step

          ! nincr_step, dtime, and incr_count are updated on each
          !   increment of the dwell episode

          nincr_step = 1

          if (start_dwell) then
            dtime = initial_time_frac*loading%step_target_time_incr(istep)

          else
            dtime = min(loading%dwell_max_strain_incr*length(idir)/ &
                & curr_vel(idir), loading%step_target_time_incr(istep), &
                & dwell_time_remaining)
          end if

          incr_count = 1

        else
          dtime_step = abs((i_end_load - i1_end_load)/loading%load_rate(istep))

          ! nincr_step and dtime are static for each load/unload step

          nincr_step = max(nint(dtime_step/loading%step_target_time_incr(istep)), 1)
          dtime = dtime_step/nincr_step
          incr_count = 1
          dwell_time_remaining = 0.0d0
        end if
      end if ! First increment in step

      ! Update time

      time = time + dtime

      ! Calculate target loads

      step_frac = real(incr_count)/real(nincr_step)

      if (istep .eq. 1) then
        target_load = loading%target_load(1, :)*step_frac

      else
        target_load = loading%target_load(istep, :)*step_frac + &
            & loading%target_load(istep - 1, :)*(1.0d0 - step_frac)
      end if

      if (myid .eq. 0) then
        write (field, '(f0.6)') time
        if (field(1:1) .eq. '.') field = '0'//field(1:15)
        time_string = field
        write (field, '(f0.6)') dtime
        if (field(1:1) .eq. '.') field = '0'//field(1:15)
        dtime_string = field

        write (*, '(a,i0,a,a,a,a,a)') 'Info   :   - &
            &Increment ', incr, ': t = ', trim(time_string), '&
            & secs, dt = ', trim(dtime_string), ' secs'
      end if

      ! Iterate on surface vel (single vel mode - fixed from
      ! previous time increment)

      ! Store initial load and update previous load

      initial_load = curr_load
      prev_load = curr_load

      ! Initialize previous vel

      prev_vel = 0

      ! For starting reload apply zero vel

      if (start_reload) then
        if (myid .eq. 0) then
          write (*, '(a)') 'Info   :   - Zero-vel iteration'
        end if

        results%vel = 0.0d0

        call vel_iteration(mesh, crys, loading, exec, &
            & results_prev, results, printing, curr_load, &
            & dtime, incr)
      end if

      ! Initialize vel for first time increment or start of reload episode

      if ((incr .eq. 1) .or. start_reload) then
        curr_vel(1) = ((target_load(1) - curr_load(1))/area0(4) &
            & - nuavg*(target_load(2) - curr_load(2))/area0(6) &
            & - nuavg*(target_load(3) - curr_load(3))/area0(2))/ &
            & eavg*length(1)/dtime
        curr_vel(2) = (-nuavg*(target_load(1) - curr_load(1))/area0(4) &
            & + (target_load(2) - curr_load(2))/area0(6) &
            & - nuavg*(target_load(3) - curr_load(3))/area0(2))/ &
            & eavg*length(2)/dtime
        curr_vel(3) = (-nuavg*(target_load(1) - curr_load(1))/area0(4) &
            & - nuavg*(target_load(2) - curr_load(2))/area0(6) &
            & + (target_load(3) - curr_load(3))/area0(2))/eavg* &
            & length(3)/dtime

        results%vel(indx) = curr_vel(1)*results%coo(indx)/length(1)
        results%vel(indy) = curr_vel(2)*results%coo(indy)/length(2)
        results%vel(indz) = curr_vel(3)*results%coo(indz)/length(3)
      end if

      ! Initialize vel for start of unload - back off yield surface

      if (start_unload) then
        curr_vel = curr_vel/10.0d0
        results%vel = results%vel/10.0d0
      end if

      ! Initialize vel for start of dwell episode

      if (start_dwell) then
        if (prev_action .eq. "loading") then
          curr_vel = initial_load_dwell_vel

        else if (prev_action .eq. "unloading") then
          curr_vel = initial_unload_dwell_vel
        end if

        results%vel(indx) = curr_vel(1)*results%coo(indx)/length(1)
        results%vel(indy) = curr_vel(2)*results%coo(indy)/length(2)
        results%vel(indz) = curr_vel(3)*results%coo(indz)/length(3)
      end if

      ! Store the vel at start of increment

      initial_vel = curr_vel

      ! Initialize iteration counter

      bc_iter_1 = 0

      bc_iteration_1: do while (.true.)

        bc_iter_1 = bc_iter_1 + 1

        if (myid .eq. 0) then
          write (*, '(a,i0)') 'Info   :     > trial_bc1: &
              &Iteration ', bc_iter_1
        end if

        call vel_iteration(mesh, crys, loading, exec, &
            & results_prev, results, printing, curr_load, &
            & dtime, incr)

        ! Print bcs_iter data

        if (myid .eq. 0) then
          ! For debug purpose only:
          ! write(printing%bcs_iter_1_u,'(3(i8),10(e15.5))') istep, incr, &
          ! & bc_iter_1, curr_vel, curr_load, target_load, dtime

          if (minval(abs(curr_vel)) .le. 0.0010) then
            write (*, '(a,3(e12.2))') 'Info   :       . &
                &Velocity:    ', curr_vel

          else
            write (*, '(a,3(f12.6))') 'Info   :       . &
                &Velocity:    ', curr_vel
          end if

          write (*, '(a,3(f12.6))') 'Info   :       . &
              &Load:        ', curr_load
          write (*, '(a,3(f12.6))') 'Info   :       . &
              &Target load: ', target_load
        end if

        ! Check if force is within range

        if ((abs(curr_load(idir) - target_load(idir)) .lt. &
            & max(exec%load_tol_abs, 0.1d0*abs(target_load(idir) - &
            & initial_load(idir)))) .or. start_dwell .or. start_reload &
            & .or. (curr_action .eq. "unloading")) then
          ! Within tolerance - exit loop

          exit

        else
          ! Not within tolerance - update guess for surface vel

          if (bc_iter_1 .eq. 1) then
            if (curr_action .eq. "dwelling") then
              vel_scale_factor = 0.9d0

            else
              vel_scale_factor = 1.1d0
            end if

          else
            vel_scale_factor = (1 - prev_vel(idir)/curr_vel(idir))* &
               & (target_load(idir) - prev_load(idir))/ &
               & (curr_load(idir) - prev_load(idir)) + &
               & prev_vel(idir)/curr_vel(idir)
          end if

          prev_load = curr_load
          prev_vel = curr_vel
          curr_vel = curr_vel*vel_scale_factor
          results%vel = results%vel*vel_scale_factor
        end if

        ! Check if maximum number of iteration are exceeded

        if (bc_iter_1 .eq. exec%max_bc_iter) then
          call par_quit('Error  :     > Maximum number of boundary&
              & condition iterations exceeded.', exec%clock_start)
        end if
      end do bc_iteration_1

      ! Iterate on surface vel (all three vel modes)

      ! Initialize iterataion counter

      bc_iter_2 = 0

      ! Calculate perturbation magnitude

      if (start_dwell) then
        pert_mag = 1.0d-6

      else
        pert_mag = max(0.1d0*abs(curr_vel(idir) - initial_vel(idir)), &
            & 0.01d0*abs(curr_vel(idir)))
      end if

      bc_iteration_2: do while (.true.)

        bc_iter_2 = bc_iter_2 + 1

        if (myid .eq. 0) then
          write (*, '(a,i0)') 'Info   :     > trial_bc2: Iteration ', &
              & bc_iter_2
        end if

        ! Perturb vel

        do pert_dir = 1, 3
          if (myid .eq. 0) then
            select case (pert_dir)

            case (1)
              write (*, '(a)') 'Info   :       . &
                  &Perturbing x-vel'

            case (2)
              write (*, '(a)') 'Info   :       . &
                  &Perturbing y-vel'

            case (3)
              write (*, '(a)') 'Info   :       . &
                  &Perturbing z-vel'

            end select
          end if

          ! Apply perturbation

          if (curr_load(pert_dir) .lt. target_load(pert_dir)) then
            pert_sign = 1.0d0

          else
            pert_sign = -1.0d0
          end if

          bak_vel = results%vel

          select case (pert_dir)

          case (1)
            results%vel(indx) = results%vel(indx) + pert_sign* &
                & pert_mag*results%coo(indx)/length(1)

          case (2)
            results%vel(indy) = results%vel(indy) + pert_sign* &
                & pert_mag*results%coo(indy)/length(2)

          case (3)
            results%vel(indz) = results%vel(indz) + pert_sign* &
                & pert_mag*results%coo(indz)/length(3)

          end select

          call vel_iteration(mesh, crys, loading, exec, &
              & results_prev, results, printing, trial_load, &
              & dtime, incr)

          results%vel = bak_vel

          if (myid .eq. 0) then
            write (*, '(a,3(f12.6))') 'Info   :     > &
                &Trial load: ', trial_load
          end if

          ! Compute matrix coefficients

          a(1, pert_dir) = (trial_load(1) - curr_load(1))/ &
              & (pert_sign*pert_mag)
          a(2, pert_dir) = (trial_load(2) - curr_load(2))/ &
              & (pert_sign*pert_mag)
          a(3, pert_dir) = (trial_load(3) - curr_load(3))/ &
              & (pert_sign*pert_mag)
        end do ! Perturbate vel

        ! Update vel field

        if (myid .eq. 0) then
          write (*, '(a)') 'Info   :     > Updating vel field'
        end if

        ! Solve for surface vel increment

        ! From right-hand side array

        b = target_load - curr_load

        ! Solve system of equations

        call solve_lin_sys_3(a, b, delta_vel)

        ! Apply new vel boundary conditions

        results%vel(indx) = results%vel(indx) + delta_vel(1)*results%coo(indx)/ &
            & length(1)
        results%vel(indy) = results%vel(indy) + delta_vel(2)*results%coo(indy)/ &
            & length(2)
        results%vel(indz) = results%vel(indz) + delta_vel(3)*results%coo(indz)/ &
            & length(3)

        curr_vel = curr_vel + delta_vel

        call vel_iteration(mesh, crys, loading, exec, &
            & results_prev, results, printing, curr_load, &
            & dtime, incr)

        ! Print bcs_iter data

        if (myid .eq. 0) then
          ! For debug purposes only
          !write(printing%bcs_iter_2_u,'(3(i8),11(e15.5))') istep, incr, &
          !    & bc_iter_2, curr_vel, curr_load, target_load, dtime, &
          !    & pert_mag

          if (minval(abs(curr_vel)) .le. 0.0010) then
            write (*, '(a,3(e12.2))') 'Info   :       . &
                &Velocity:    ', curr_vel

          else
            write (*, '(a,3(f12.6))') 'Info   :       . &
                &Velocity:    ', curr_vel
          end if

          write (*, '(a,3(f12.6))') 'Info   :       . &
                  &Load:        ', curr_load
          write (*, '(a,3(f12.6))') 'Info   :       . &
                  &Target load: ', target_load
        end if

        ! Check if force is within range

        if (maxval(abs(curr_load - target_load)) .lt. exec%load_tol_abs) then
          ! Within tolerance - exit loop

          exit
        end if

        ! Check if maximum number of iteration are exceeded

        if (bc_iter_2 .eq. exec%max_bc_iter) then
          call par_quit('Error  :     > Maximum number of boundary&
              & condition iterations exceeded', exec%clock_start)
        end if

        ! Update perturbation magnitude

        pert_mag = maxval(abs(delta_vel))
      end do bc_iteration_2

      ! Update rstar and sig_vec_n (?)

      ! Write boundary condition iteration statistics - debug only
      !if (myid .eq. 0) then
      !
      !    write(printing%bcs_iter_log_u,'(4(i12))') istep, incr, bc_iter_1, &
      !        & bc_iter_2
      !
      !end if

      ! Update state (at center of each element), geometry and tool bc's

      ! Update:
      ! - coo @(t+dt), using dt and vel @(t+dt)

      call update_state_evp(mesh, crys, exec, results_prev, results, dtime)

      call finalize_res_stressstrain(mesh, crys, results, load, area)

      curr_load(1) = load(2, 1) ! legacy format: 4,1
      curr_load(2) = load(4, 2) !                6,2
      curr_load(3) = load(6, 3) !                2,3

      ! Calculate mesh dimensions

      call calc_mesh_dim(results, length, indx, indy, indz)

      ! Calculate macroscopic engineering strain

      do i = 1, 3
        macro_eng_strain(i) = length(i)/length0(i) - 1.0d0
      end do

      ! Calculate the current increment macroscopic strain_eq

      curr_eqstrain = (2.0d0/3.0d0)* &
          & dsqrt((3.0d0/2.0d0)*((macro_eng_strain(1)**2.0d0) &
          & + (macro_eng_strain(2)**2.0d0) &
          & + (macro_eng_strain(3)**2.0d0)))

      ! Print the macroscopic strain values to console for monitoring

      if (myid .eq. 0) then
        write (*, '(a,f12.6)') 'Info   :       . &
            &Maximum eng. strain: ', maxval(abs(macro_eng_strain(:)))
        write (*, '(a,f12.6)') 'Info   :       . &
            &Current eqv. strain: ', curr_eqstrain
      end if

      ! Calculate dwell_time_remaining

      if (curr_action .eq. "dwelling") then
        dwell_time_remaining = dwell_time_remaining - dtime
      end if

      ! Write loads to post.force# files

      if (myid .eq. 0) then
        if (printing%print_forces) then
          call write_force_file_data(mesh, printing, istep, incr, load, area, time, length)
        end if
      end if

      call finalize_res (mesh, crys, dtime, results_prev,results, load, area)

      ! Store vel if at beginning of dwell episode

      if (start_dwell) then
        if (prev_action .eq. "loading") then
          initial_load_dwell_vel = curr_vel

        else if (prev_action .eq. "unloading") then
          initial_unload_dwell_vel = curr_vel
        end if
      end if

      ! Output computed quantities at end of step

      if (((curr_action .ne. "dwelling") .and. (incr_count .eq. nincr_step)) &
          & .or. ((curr_action .eq. "dwelling") .and. (dwell_time_remaining &
          & .lt. time_tol*loading%step_target_time_incr(istep)))) then

        if (loading%step_print(istep)) then

          call write_res(istep, mesh, crys, results, dtime, printing)
          
          if (printing%print_restart) then
            call write_restart_field(istep, mesh, results)

            call write_triaxclr_restart(mesh, istep, &
                & curr_load, prev_load, .true., incr + 1, time, &
                & load, area, area0, length, length0, &
                & curr_vel, prev_action, curr_action, &
                & initial_load_dwell_vel, initial_unload_dwell_vel)
          end if
        end if

        if (istep .eq. loading%num_steps) then
          ! Job finished successfully

          if (myid .eq. 0 ) call write_dot_sim_file_complete_steps(printing, loading, istep)

          call par_quit('Info   : Final step terminated. Simulation&
              & completed successfully.', exec%clock_start)
        end if

        ! Check that maximum strain has not been exceeded

        if (maxval(abs(macro_eng_strain)) .ge. exec%max_strain) then
          if (myid .eq. 0 ) call write_dot_sim_file_complete_steps(printing, loading, istep)

          call par_quit('Error  :     > Maximum eng. strain exceeded.', exec%clock_start)
        end if

        ! Check that maximum strain_eq has not been exceeded

        if (curr_eqstrain .ge. exec%max_eqstrain) then
          if (myid .eq. 0 ) call write_dot_sim_file_complete_steps(printing, loading, istep)

          call par_quit('Error  :     > Maximum eqv. strain exceeded.', exec%clock_start)
        end if

        istep = istep + 1
        first_incr_in_step = .true.
        incr_count = 1

      else
        first_incr_in_step = .false.
        incr_count = incr_count + 1

        ! Compute new time increment for dwell

        if (curr_action .eq. "dwelling") then
          dtime = min(loading%dwell_max_strain_incr*length(idir)/abs(curr_vel(idir)), &
              &loading%step_target_time_incr(istep), dwell_time_remaining)
          nincr_step = nincr_step + 1
        end if

        ! Check that maximum strain has not been exceeded

        if (maxval(abs(macro_eng_strain)) .ge. exec%max_strain) then
          ! mpk - 10/2021: Don't write that we have completed the final
          !   step. We are currently in the middle of a step.
          if (myid .eq. 0 ) call write_dot_sim_file_complete_steps(printing, loading, istep - 1)

          call par_quit('Info   :     > Maximum eng. strain exceeded.', exec%clock_start)
        end if

        ! Check that maximum strain_eq has not been exceeded

        if (curr_eqstrain .ge. exec%max_eqstrain) then
          ! mpk - 10/2021: Don't write that we have completed the final
          !   step. We are currently in the middle of a step.
          if (myid .eq. 0 ) call write_dot_sim_file_complete_steps(printing, loading, istep - 1)

          call par_quit('Info   :     > Maximum eqv. strain exceeded.', exec%clock_start)
        end if
      end if

      incr = incr + 1
      start_reload = .false.
      start_unload = .false.
      start_dwell = .false.
    end do time_stepping

    return

  end subroutine driver_triaxclr

  !===========================================================================

  subroutine print_headers(printing)

    ! Print output file headers. For debug purposes only - commented out above

    !---------------------------------------------------------------------------

    ! Locals:
    ! iostatus: Returns value about file command success

    type(printing_type), intent(in) :: printing
    integer :: iostatus

    !---------------------------------------------------------------------------

    if (myid .eq. 0) then
      ! Open post.bcs_iter_1 file

      open (printing%bcs_iter_1_u, file='post.bcs_iter_1', iostat=iostatus)

      if (iostatus .ne. 0) then
        call par_quit('Error  :     > io Failure to open &
            &post.bcs_iter_1 file.')
      end if

      ! Write the header to the post.bcs_iter_1 file

      write (printing%bcs_iter_1_u, '(a)') ' %   step    incr bc_iter_1 &
          &vel_x          vel_y          vel_z          &
          &load_x         load_y         load_z         &
          &target_load_x  target_load_y  target_load_z  &
          &dtime'

      ! Open post.bcs_iter_2 file

      open (printing%bcs_iter_2_u, file='post.bcs_iter_2', iostat=iostatus)

      if (iostatus .ne. 0) then
        call par_quit('Error  :     > io Failure to open &
            &post.BCS_iter_2 file.')
      end if

      ! Write the header to the post.BCS_iter_2 file

      write (printing%bcs_iter_2_u, '(a)') ' %   step    incr bc_iter_2 &
          &vel_x          vel_y          vel_z          &
          &load_x         load_y         load_z         &
          &target_load_x  target_load_y  target_load_z  &
          &dtime         pert_mag'

      ! Open post.bcs_iter_log file

      open (printing%bcs_iter_log_u, file='post.bcs_iter_log', &
          & iostat=iostatus)

      if (iostatus .ne. 0) then
        call par_quit('Error  :     > io Failure to open &
            &post.bcs_iter_log file.')
      end if

      ! Write the header to the post.bcs_iter_log file

      write (printing%bcs_iter_log_u, '(a)') '%       step        incr   &
          &bc_iter_1   bc_iter_2'
    end if

    return

  end subroutine print_headers

  !===========================================================================

  subroutine read_triaxclr_restart(mesh, istep, curr_load, prev_load, &
      & first_incr_in_step, incr, time, load, area, area0, &
      & length, length0, curr_vel, prev_action, curr_action, &
      & initial_load_dwell_vel, initial_unload_dwell_vel)

    ! Read TriaxCLR restart information.

    !---------------------------------------------------------------------------

    ! Arguments:

    type(mesh_type), intent(in) :: mesh
    logical, intent(out) :: first_incr_in_step
    integer, intent(out) :: istep
    integer, intent(out) :: incr
    character(len=10), intent(out) :: prev_action, curr_action
    real(rk), intent(out) :: curr_load(3)
    real(rk), intent(out) :: prev_load(3)
    real(rk), intent(out) :: time
    real(rk), intent(out) :: load(mesh%num_fasets, 3)
    real(rk), intent(out) :: area(mesh%num_fasets)
    real(rk), intent(out) :: area0(mesh%num_fasets)
    real(rk), intent(out) :: length(3), length0(3)
    real(rk), intent(out) :: curr_vel(3)
    real(rk), intent(out) :: initial_load_dwell_vel(3)
    real(rk), intent(out) :: initial_unload_dwell_vel(3)

    ! Locals:

    integer :: myunit
    integer :: rst_num
    logical :: file_exists
    character(len=8) :: rst_num_str
    character(len=64) :: filename
    integer :: isurf

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

    read (myunit) istep
    read (myunit) curr_load
    read (myunit) prev_load
    read (myunit) first_incr_in_step
    read (myunit) incr
    read (myunit) time

    do isurf = 1, mesh%num_fasets
      read (myunit) load(isurf, :)
    end do

    read (myunit) area
    read (myunit) area0
    read (myunit) length
    read (myunit) length0
    read (myunit) curr_vel
    read (myunit) prev_action
    read (myunit) curr_action
    read (myunit) initial_load_dwell_vel
    read (myunit) initial_unload_dwell_vel

    if (myid .eq. 0) then
      write (*, '(a)') 'Info   : Reading restart control information...'
      write (*, '(a)') 'Info   :   - Previous simulation ended with final:'
      write (*, '(a, i0)') 'Info   :     > Increments: ', incr - 1
      write (*, '(a, i0)') 'Info   :     > Steps:      ', istep
      write (*, '(a, e14.4)') 'Info   :     > Time:  ', time
      write (*, '(a, 3(e14.4))') 'Info   :     > Normal loads:  ', &
          & curr_load
    end if

    close (myunit)

    ! Reinitialize step and time if new files
    ! In the future, we will have to handle the append case, which will continue
    ! sequentially, rather than reinitializing
    istep = 1
    incr = 1
    time = 0.0d0

    return

  end subroutine read_triaxclr_restart

end module driver_triaxclr_mod
