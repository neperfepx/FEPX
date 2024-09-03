! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module driver_triaxcsr_mod

! Driver for triaxial load control with constant strain rate.

! Contains subroutines:
! driver_triaxcsr: Primary driver for triaxial csr control simulations.
! driver_triaxcsr_restart: Read triaxial csr restart information.

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use utils_mod
  use mesh_type_mod
  use res_init_mod
  use solveit_evp_mod2
  use finalize_res_mod
  use driver_triax_utilities_mod
  use kinematics_mod
  use matrix_operations_mod, only: solve_lin_sys_3, &
      & strain_equiv_3x3, vec6_mat_symm
  use write_res_mod
  use restart_mod
  use gather_scatter_mod, only: trace, part_gather
  use parallel_mod, only: par_message, par_quit

  implicit none

! Private

contains

  subroutine driver_triaxcsr(mesh, crys, loading, exec, &
             & results, printing)

    ! Primary driver for triaxial load control with constant strain rate.

    !---------------------------------------------------------------------------

    ! Arguments:

    type(mesh_type), intent(inout) :: mesh
    type(crys_type) :: crys (:)
    type(loading_type), intent(inout) :: loading
    type(exec_type), intent(inout) :: exec
    type(results_type), intent(inout) :: results
    type(printing_type), intent(inout) :: printing

    ! Locals:
    ! Need to be defined - jc

    logical :: step_complete
    logical :: is_necking
    logical :: is_limit_tripped
    integer :: incr
    real(rk) :: dtime, prev_dt, time
    real(rk) :: bak_vel(dof_sub:dof_sup)
    integer :: m_elt, i, ii
    integer :: istep
    real(rk) :: load(mesh%num_fasets, 3)
    real(rk) :: area0(mesh%num_fasets), area(mesh%num_fasets)
    real(rk) :: length(3), length0(3)

    ! optional parameters.
    real(rk) :: initial_vel

    ! control variables.
    integer  :: pdir, sdir, tdir
    integer  :: ind((dof_sup - dof_sub + 1)/3, 3)
    integer  :: bc_iter
    real(rk) :: eavg, nuavg
    real(rk) :: curr_vel(3)
    real(rk) :: curr_load(3), prev_load(3)
    real(rk) :: ideal_load(3), trial_load(3)
    real(rk) :: load_step
    real(rk) :: step_frac
    real(rk) :: s_pert_mag, t_pert_mag
    real(rk) :: s_pert_sign, t_pert_sign
    real(rk) :: s_delta_vel, t_delta_vel
    real(rk) :: s_load_err, t_load_err, max_load_err
    real(rk) :: a(3, 3)
    real(rk) :: b(3)
    real(rk) :: coeffs(3)
    real(rk) :: macro_eng_strain(3)
    real(rk) :: curr_eqstrain
    character(len=16) :: field, time_string, dtime_string
    real(rk) :: clock_end

    type(results_type) :: results_prev

    call res_qpt_alloc (mesh, results_prev)
    call res_qpt_init_prev (results, results_prev)

    !---------------------------------------------------------------------------

    m_elt = elt_sup - elt_sub + 1

    ! Designate primary (control), secondary, and tertiary directions

    select case (loading%loading_direction)

    case (1)
      pdir = 1
      sdir = 2
      tdir = 3

    case (2)
      pdir = 2
      sdir = 3
      tdir = 1

    case (3)
      pdir = 3
      sdir = 1
      tdir = 2

    case default
      call par_quit('Error  :     > Invalid control direction provided.')

    end select

    ! Store indices for x, y, and z degrees of freedom

    do i = 1, ((dof_sup - dof_sub + 1)/3)
      ind(i, 1) = dof_sub + 3*(i - 1)
      ind(i, 2) = dof_sub + 3*(i - 1) + 1
      ind(i, 3) = dof_sub + 3*(i - 1) + 2
    end do

    ! Recording initial vel

    initial_vel = loading%step_velocity(1)

    ! Initialization

    is_necking = .false.
    is_limit_tripped = .false.

    if (exec%restart) then
      call read_restart_field(mesh, results)

      call read_triaxcsr_restart(loading, mesh, results, &
          & printing, istep, curr_load, prev_load, step_complete, &
          & dtime, incr, time, load, area, area0, length, &
          & length0, curr_vel, s_pert_mag, t_pert_mag)

    else
      ! Initialize areas and load arrays

      load = 0.0d0
      curr_load = 0.0d0
      area = 0.0d0
      area0 = 0.0d0

      ! Compute initial area (area0)

      call mesh_surfaceareas(mesh, results%coo, area0)

      !if (myid .eq. 0) then
      !    write(*,'(a25,6f12.6)') 'initial areas (area0):   ', area0
      !end if

      ! Compute initial mesh dimensions

      call calc_mesh_dim(results, length, ind(:, 1), ind(:, 2), ind(:, 3))
      length0 = length

      ! Initialize state

      results%ori = spread(mesh%ori, 4, nqpt)

      ! Initialize integrated quantities

      results%slip = 0.0d0
      results%strain_pl = 0.0d0
      results%defrate = 0.0d0

      ! Initialize deformation control

      incr = 0
      time = 0.0d0
      istep = 1
      incr = 0
      step_complete = .false.
      dtime = loading%step_dt_max(1)
      s_pert_mag = exec%min_pert_frac*abs(initial_vel)
      t_pert_mag = exec%min_pert_frac*abs(initial_vel)
      curr_eqstrain = 0.0d0

      ! Estimate bulk elastic moduli
      ! Assume uniform texture and equal element volumes
      ! Consider first crystal
      call est_avg_mod(mesh, crys(1), eavg, nuavg)

      ! Initialize vel field using elasticity

      curr_vel(1) = (loading%target_load(1, 1)/area0(4) - &
          & nuavg*loading%target_load(1, 2)/area0(6) - &
          & nuavg*loading%target_load(1, 3)/area0(2))*length0(1)/eavg

      curr_vel(2) = (-nuavg*loading%target_load(1, 1)/area0(4) + &
          & loading%target_load(1, 2)/area0(6) - &
          & nuavg*loading%target_load(1, 3)/area0(2))*length0(2)/eavg

      curr_vel(3) = (-nuavg*loading%target_load(1, 1)/area0(4) - &
          & nuavg*loading%target_load(1, 2)/area0(6) + &
          & loading%target_load(1, 3)/area0(2))*length0(3)/eavg

      curr_vel = curr_vel*abs(initial_vel/curr_vel(pdir))

      results%vel(ind(:, 1)) = curr_vel(1)*results%coo(ind(:, 1))
      results%vel(ind(:, 2)) = curr_vel(2)*results%coo(ind(:, 2))
      results%vel(ind(:, 3)) = curr_vel(3)*results%coo(ind(:, 3))

      ! Print initial values
      call write_res(0, mesh, crys, results, 0.0d0, printing)
    end if

    if (myid .eq. 0) then
      ! Write the header to file post.force1-6

      if (printing%print_forces) then
        call write_force_file_headers(mesh, printing, 2)

        ! If virgin sample, write 0th step

        if (printing%restart .eqv. .false.) then
          call write_force_file_data(mesh, printing, 0, 0, load, area0, 0.0d0, length0)
        end if
      end if

      if (printing%print_conv) then
        call write_conv_file_headers(printing)
      end if

      ! Open post.bcs_iter file - for debug only
      ! open (printing%bcs_iter_1_u, file='post.bcs_iter_1',iostat=iostatus)

      ! if (iostatus .ne. 0) then
      ! call par_quit('Error  :     > io Failure to open &
      ! &post.bcs_iter file.')
      ! end if

      ! Write the header to file post.bcs_iter

      ! write(printing%bcs_iter_1_u,'(a)')  '%   step    incr bc_iter    &
      ! &vel_x          vel_y          vel_z          &
      ! &load_x         load_y         load_z         &
      ! &target_load_x  target_load_y  target_load_z  &
      ! &dtime          s_pert_mag     t_pert_mag'
    end if

    ! if (myid .eq. 0) write(*, '(a)') 'Info   :   - Starting t-stepping &
    ! &loop for anisotropic viscoplastic solution'
    ! if (myid .eq. 0) write(*,'(a)') 'Info   : Running step 1...'

    ! Time stepping loop

    time_stepping: do

      ! Iterate to find new configuration and material state.

      incr = incr + 1

      if (step_complete .and. (myid .eq. 0)) then
        write (*, '(a,i0,a)') 'Info   : Running step ', istep, '...'
      end if

      ! Initialize shear rates

      call res_qpt_init_prev (results, results_prev)

      ! Calculate new time increment
      ! (based on time_increment from load_disp_control_mod)

      if (incr .gt. 1) then
        prev_dt = dtime
        dtime = loading%step_dt_max(istep)

        ! Check if we are close to target and need to adjust the time step

        trial_load = curr_load + (curr_load - prev_load)*dtime/prev_dt

        if (loading%target_sign(istep)*(trial_load(pdir) - &
            & loading%target_load(istep, pdir)) .ge. 0.0) then
          if (abs(trial_load(pdir) - curr_load(pdir)) .gt. vtiny) then
            dtime = (loading%target_load(istep, pdir) - curr_load(pdir))/ &
                & (trial_load(pdir) - curr_load(pdir))*dtime* &
                & exec%dtime_factor

            if (dtime .gt. loading%step_dt_max(istep)) then
              dtime = loading%step_dt_max(istep)
            end if

            if (dtime .lt. loading%step_dt_min(istep)) then
              dtime = loading%step_dt_min(istep)
            end if
          end if
        end if
      end if

      ! Update the total time

      time = time + dtime

      if (myid .eq. 0) then
        write (field, '(f0.4)') time
        if (field(1:1) .eq. '.') field = '0'//field(1:15)
        time_string = field
        write (field, '(f0.4)') dtime
        if (field(1:1) .eq. '.') field = '0'//field(1:15)
        dtime_string = field

        write (*, '(a,i0,a,a,a,a,a)') 'Info   :   - &
            &Increment ', incr, ': t = ', trim(time_string), '&
            & secs, dt = ', trim(dtime_string), ' secs'
      end if

      ! Store previous load

      prev_load = curr_load

      ! Update initial guess of vel field for load reversal

      if (step_complete) then
        results%vel = loading%vel_factor(istep)*results%vel
        curr_vel = loading%vel_factor(istep)*curr_vel
      end if

      ! Surface vel iteration loop

      bc_iter = 0

      bc_iteration: do

        ! Check number of boundary condition iterations

        if (bc_iter .ge. exec%max_bc_iter) then
          call par_quit('Error  :     > Maximum number of &
              &boundary conditions reached.', exec%clock_start)
        end if

        bc_iter = bc_iter + 1

        if (myid .eq. 0) then
          write (*, '(a,i0)') 'Info   :     > trial_bc: Iteration ', &
              & bc_iter
        end if

        ! Begin vel iteration

        call vel_iteration(mesh, crys, loading, exec, &
            & results_prev, results, printing, curr_load, &
            & dtime, incr)

        ! Check load tolerances

        if (istep .eq. 1) then
          step_frac = curr_load(pdir)/loading%target_load(1, pdir)
          ideal_load = step_frac*(loading%target_load(1, :))

        else
          step_frac = (curr_load(pdir) - loading%target_load(istep - 1, pdir))/ &
              & (loading%target_load(istep, pdir) - loading%target_load(istep - 1, pdir))
          ideal_load = step_frac*loading%target_load(istep, :) &
              & + (1 - step_frac)*loading%target_load(istep - 1, :)
        end if

        ! Write output to console and file

        if (myid .eq. 0) then
          ! For debug only
          ! write(printing%bcs_iter_1_u,'(3(i8),12(e15.5))') istep, incr, &
          ! & bc_iter, curr_vel, curr_load, ideal_load, dtime, &
          ! & s_pert_mag, t_pert_mag

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
              &Target load: ', ideal_load
        end if

        ! Check if loads are within tolerance

        s_load_err = abs(curr_load(sdir) - ideal_load(sdir))
        t_load_err = abs(curr_load(tdir) - ideal_load(tdir))
        max_load_err = max(exec%load_tol_abs, exec%load_tol_rel* &
            & abs(curr_load(pdir)))

        if ((s_load_err .le. max_load_err) .and. &
            & (t_load_err .le. max_load_err)) then
          ! Load within tolerance

          exit
        end if

        ! Perturb secondary direction vel

        if (myid .eq. 0) then
          write (*, '(a)') 'Info   :       . Perturbing surface &
              &vel in secondary direction'
        end if

        if (curr_load(sdir) .le. ideal_load(sdir)) then
          s_pert_sign = 1.0d0
        else
          s_pert_sign = -1.0d0
        end if

        bak_vel = results%vel

        results%vel(ind(:, sdir)) = results%vel(ind(:, sdir)) &
            & + s_pert_sign*s_pert_mag*results%coo(ind(:, sdir))/ &
            & length(sdir)

        call vel_iteration(mesh, crys, loading, exec, &
            & results_prev, results, printing, trial_load, &
            & dtime, incr)

        results%vel = bak_vel

        if (myid .eq. 0) then
          write (*, '(a,3(f12.6))') 'Info   :     > Trial load: ', &
              & trial_load
        end if

        ! Compute matrix coefficients

        a(1, 1) = trial_load(1) - curr_load(1)
        a(2, 1) = trial_load(2) - curr_load(2)
        a(3, 1) = trial_load(3) - curr_load(3)

        ! Perturb tertiary direction vel

        if (myid .eq. 0) then
          write (*, '(a)') 'Info   :       . Perturbing surface &
              &vel in tertiary direction'
        end if

        if (curr_load(tdir) .le. ideal_load(tdir)) then
          t_pert_sign = 1.0d0

        else
          t_pert_sign = -1.0d0
        end if

        bak_vel = results%vel

        results%vel(ind(:, tdir)) = results%vel(ind(:, tdir)) &
            & + t_pert_sign*t_pert_mag*results%coo(ind(:, tdir))/ &
            & length(tdir)

        call vel_iteration(mesh, crys, loading, exec, &
            & results_prev, results, printing, trial_load, &
            & dtime, incr)
        results%vel = bak_vel

        if (myid .eq. 0) then
          write (*, '(a,3(f12.6))') 'Info   :     > Trial load: ', &
              & trial_load
        end if

        ! Compute matrix coefficients

        a(1, 2) = trial_load(1) - curr_load(1)
        a(2, 2) = trial_load(2) - curr_load(2)
        a(3, 2) = trial_load(3) - curr_load(3)

        ! Compute change to surface vel

        if (istep .eq. 1) then
          load_step = loading%target_load(istep, pdir)

          a(1, 3) = -loading%target_load(istep, 1)/load_step
          a(2, 3) = -loading%target_load(istep, 2)/load_step
          a(3, 3) = -loading%target_load(istep, 3)/load_step

          b(1) = -curr_load(1)
          b(2) = -curr_load(2)
          b(3) = -curr_load(3)

        else
          load_step = loading%target_load(istep, pdir) - &
              & loading%target_load(istep - 1, pdir)

          a(1, 3) = (loading%target_load(istep - 1, 1) - &
              & loading%target_load(istep, 1))/load_step
          a(2, 3) = (loading%target_load(istep - 1, 2) - &
              & loading%target_load(istep, 2))/load_step
          a(3, 3) = (loading%target_load(istep - 1, 3) - &
              & loading%target_load(istep, 3))/load_step

          b(1) = loading%target_load(istep - 1, 1) - curr_load(1)
          b(2) = loading%target_load(istep - 1, 2) - curr_load(2)
          b(3) = loading%target_load(istep - 1, 3) - curr_load(3)
        end if

        ! Solve linear system of equations

        call solve_lin_sys_3(a, b, coeffs)

        ! Extract changes to vel field

        s_delta_vel = coeffs(1)*s_pert_sign*s_pert_mag
        t_delta_vel = coeffs(2)*t_pert_sign*t_pert_mag

        ! Apply new surface vel

        results%vel(ind(:, sdir)) = results%vel(ind(:, sdir)) &
            & + s_delta_vel*results%coo(ind(:, sdir))/length(sdir)
        results%vel(ind(:, tdir)) = results%vel(ind(:, tdir)) &
            & + t_delta_vel*results%coo(ind(:, tdir))/length(tdir)

        curr_vel(sdir) = curr_vel(sdir) + s_delta_vel
        curr_vel(tdir) = curr_vel(tdir) + t_delta_vel

        ! Update perturbation magnitudes
        s_pert_mag = max(abs(s_delta_vel), exec%min_pert_frac*curr_vel(pdir))
        t_pert_mag = max(abs(t_delta_vel), exec%min_pert_frac*curr_vel(pdir))
      end do bc_iteration

      ! Update state (at center of each element), geometry

      ! Update:
      ! coo @(t+dt), using dt and vel @(t+dt)

      ! Element centroid:
      ! e_elas_kk: volumetric lattice strain @(t+dt)

      call update_state_evp(mesh, crys, exec, results_prev, results, dtime)

      call finalize_res_stressstrain(mesh, crys, results, load, area)

      curr_load(1) = load(2, 1) ! 4,1 legacy faceset values
      curr_load(2) = load(4, 2) ! 6,2
      curr_load(3) = load(6, 3) ! 2,3

      ! Compute mesh dimensions

      call calc_mesh_dim(results, length, ind(:, 1), ind(:, 2), ind(:, 3))

      ! Calculate macroscopic engineering strain (ported from TriaxCLR)

      do ii = 1, 3
        macro_eng_strain(ii) = length(ii)/length0(ii) - 1.0d0
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

      ! Write loads to post.force# files

      if (myid .eq. 0) then
        if (printing%print_forces) then
          call write_force_file_data(mesh, printing, istep, incr, load, area, time, length)
        end if
      end if

      call finalize_res (mesh, crys, dtime, results_prev,results, load, area)

      ! Check if the target load is reached

      step_complete = .false.

      !if (loading%target_sign(istep)*(curr_load(pdir)-loading%target_load(istep,pdir))&
      !    & + load_tol .ge. 0.0) then
      if (abs(curr_load(pdir) - loading%target_load(istep, pdir)) .le. &
          & exec%load_tol_abs) then
        step_complete = .true.
      end if

      ! Check for necking

      if (exec%check_necking .and. ((curr_load(pdir) - prev_load(pdir)) &
          & /curr_vel(pdir) .lt. 0.0)) then
        is_necking = .true.
      end if

      ! Check time and increment count

      if ((time .ge. exec%max_total_time) .or. &
          & (incr .ge. exec%max_incr)) then
        is_limit_tripped = .true.
      end if

      ! Output computed quantities.

      if (step_complete .or. is_necking .or. is_limit_tripped) then
        write(*,*) "step = ", loading%step_print(istep)
        if (loading%step_print(istep)) then

          call write_res(istep, mesh, crys, results, dtime, printing)
          
          if (printing%print_restart) then
            call write_restart_field(istep, mesh, results)

            call write_triaxcsr_restart(mesh, istep + 1, curr_load, prev_load, &
                & step_complete, dtime, incr, time, load, &
                & area, area0, length, length0, curr_vel, s_pert_mag, &
                & t_pert_mag)
          end if
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

        if (is_necking) then
          if (myid .eq. 0 ) call write_dot_sim_file_complete_steps(printing, loading, istep)
          call par_quit('Error  :     > Specimen is necking.', exec%clock_start)
        end if

        if (is_limit_tripped) then
          if (myid .eq. 0 ) call write_dot_sim_file_complete_steps(printing, loading, istep)
          call par_quit('Error  :     > Maximum time or maximum &
              &increments exceeded.', exec%clock_start)
        end if

        if (istep .eq. loading%num_steps) then
          if (myid .eq. 0 ) call write_dot_sim_file_complete_steps(printing, loading, istep)
          ! Finalize clock values and print to console

          if (myid .eq. 0) then
            call cpu_time(clock_end)
            write (*, '(a, f10.3, a)') 'Info   : Elapsed time:', &
                & (clock_end - exec%clock_start), ' secs.'
          end if

          call par_quit('Info   : Final step terminated. Simulation&
              & completed successfully.', exec%clock_start)
        end if

        ! Advance step count

        istep = istep + 1
      end if
    end do time_stepping

    return

  end subroutine driver_triaxcsr

  !===========================================================================

  subroutine read_triaxcsr_restart(loading, mesh, results, &
      & printing, istep, curr_load, prev_load, step_complete, &
      & dtime, incr, time, load, area, area0, length, length0, curr_vel, &
      & s_pert_mag, t_pert_mag)

    ! Read triaxial csr restart information.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! Needs to be defined - jc

    type(mesh_type), intent(in) :: mesh
    type(loading_type), intent(inout) :: loading
    type(results_type) :: results
    type(printing_type), intent(inout) :: printing
    logical, intent(out) :: step_complete
    integer, intent(out) :: istep
    integer, intent(out) :: incr
    real(rk), intent(out) :: curr_load(3)
    real(rk), intent(out) :: prev_load(3)
    real(rk), intent(out) :: dtime
    real(rk), intent(out) :: time
    real(rk), intent(out) :: load(mesh%num_fasets, 3)
    real(rk), intent(out) :: area(mesh%num_fasets)
    real(rk), intent(out) :: area0(mesh%num_fasets)
    real(rk), intent(out) :: length(3)
    real(rk), intent(out) :: length0(3)
    real(rk), intent(out) :: curr_vel(3)
    real(rk), intent(out) :: s_pert_mag
    real(rk), intent(out) :: t_pert_mag

    ! Locals:
    ! myunit: Current unit number to open restart file.
    ! isurf: Generic loop index to loop over mesh surfaces.

    integer :: myunit
    integer :: isurf
    integer :: rst_num
    logical :: file_exists
    character(len=8) :: rst_num_str
    character(len=64) :: filename
    real(rk) :: diff_stress
    integer :: i

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
    read (myunit) step_complete
    read (myunit) dtime
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
    read (myunit) s_pert_mag
    read (myunit) t_pert_mag

    if (myid .eq. 0) then
      write (*, '(a)') 'Info   : Reading restart control information...'
      write (*, '(a)') 'Info   :   - Previous simulation ended with final:'
      write (*, '(a, i0)') 'Info   :     > Increments: ', incr
      write (*, '(a, i0)') 'Info   :     > Steps:      ', istep - 1
      write (*, '(a, e14.4)') 'Info   :     > Time:  ', time
      write (*, '(a, 3(e14.4))') 'Info   :     > Normal loads:  ', &
          & curr_load
    end if

    close (myunit)

    ! Reinitialize step and time if new files
    ! In the future, we will have to handle the append case, which will continue
    ! sequentially, rather than reinitializing
    istep = 1
    step_complete = .true.
    incr = 0
    time = 0.0d0

    ! Set a temporary variable to assist in printing to .sim file

    printing%restart_initial_step = istep

    ! Check for a reversal in loading direction upon restart, correct loading
    ! history if so. The loading history as assigned in process_ctrl_data_csr
    ! assumes that we are starting from 0 load. Below rectifies any directional
    ! issues with loading due to this assumption.

    if ((curr_load(loading%loading_direction) - &
        & prev_load(loading%loading_direction))* &
        & (loading%target_load(1, loading%loading_direction) - &
        & curr_load(loading%loading_direction)) .lt. 0.0d0) then
      ! Change vel direction

      results%vel = -1.0d0*results%vel

      ! And check for any load reversals on subsequent steps

      if (loading%num_steps .ge. 2) then
        diff_stress = (loading%target_load(1, &
            & loading%loading_direction) - &
            & curr_load(loading%loading_direction))

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
          diff_stress = &
              & (loading%target_load(i, &
              & loading%loading_direction)) &
              & - (loading%target_load(i - 1, &
              & loading%loading_direction))

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
            loading%vel_factor(i) = -loading%vel_factor(i)
          end if
        end do
      end if
    end if

    return

  end subroutine read_triaxcsr_restart

end module driver_triaxcsr_mod
