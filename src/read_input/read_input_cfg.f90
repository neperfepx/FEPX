! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
module read_input_cfg_mod
!
! Module to parse configuration file.
!
! Contains subroutines/functions:
! read_config: Primary subroutine to read options from the configuration file
! initialize_options: Initialize (non crys_type) options.
! Various helper subroutines.
!
!---------------------------------------------------------------------------------

! Use modules
use, intrinsic :: iso_fortran_env, only: rk => real64
use general_mod
use types_mod
use crys_type_mod
use loading_options_type_mod
use parallel_mod, only: par_quit
use types_mod
use utils_mod

! Default for this module is private variables and subroutines
private

! Explicitly set public variables and subroutines
public :: read_config

! Contains following subroutines and functions
contains

!===============================================================================

subroutine read_config(mesh_file, ori_file, phase_file, printing, loading_options, &
                     & crys, exec)

  ! Primary subroutine to read options from the configuration file

  ! Arguments:
  ! options: General options type
  ! printing: Printing options type
  ! loading_options: Boundary conditions options type
  ! crys: Crystal parameters type
  character(len=20), intent(out) :: mesh_file
  character(len=20), intent(out) :: ori_file
  character(len=20), intent(out) :: phase_file
  type(printing_type), intent(inout) :: printing
  type(loading_options_type), intent(inout) :: loading_options
  type(crys_type), allocatable, intent(inout) :: crys (:)
  type(exec_type) :: exec
  ! Locals:
  ! line: Current line being parsed
  ! flag: Input flag as determined from line
  ! value: Input value as determined from line
  ! file_id: File ID for simulation.cfg
  ! j: Looping index (number of characters in line)
  ! j0: Index of first space (demarcation between flag and value)
  ! state: Index for determining whether line has a value or not
  ! state_begin, state_in_flag, state_in_sep: parameters for the same
  ! temp: Temporary variable for user defined boundary condition allocation
  ! message: String for message printing
  character(len=1000) :: line
  character(len=100) :: flag
  character(len=200) :: value
  integer :: file_id
  integer :: j
  integer :: j0
  integer :: state
  integer, parameter :: state_begin = 1
  integer, parameter :: state_in_flag = 2
  integer, parameter :: state_in_sep = 3
  integer :: phase, number_of_phases
  character(len=500) :: message
  character(len=2) :: version
  logical :: g_s_parsed = .false., g_s0_parsed = .false.

  !-----------------------------------------------------------------------------

  mesh_file = "simulation.msh"

  ! Initialize allocated options
  call loading_options_set_default(loading_options)

  loading_options%def_control_by = "uniaxial_strain_target"

  ! Write to terminal
  if (myid .eq. 0) then
    write (*, '(a,a,a)') "Info   :   [i] Parsing file `", "simulation.cfg", "'..."
  end if

  ! Open the configuration file
  open(newunit=file_id, file="simulation.cfg", status="old", action='read', &
    & iostat=iostatus)

  if (iostatus .ne. 0) then
    call par_quit("Error  :     > Failure to open `simulation.cfg'.")
  end if

  ! Loop over entire file once to find number of phases
  do ! loop over file

    ! Read the current line of the file
    read(file_id, "(a)", iostat=ioerr) line

    ! Error checking upon reading
    if (ioerr .lt. 0) then
      call par_quit("Error  : Cannot read from simulation.cfg")
    end if

    ! Ignore lines that are commented or empty
    if ((line(1:1) .eq. "#") .or. (len_trim(line) .eq. 0)) then
      cycle
    end if

    ! Extract input flags (flag) and value(s) (value) by checking for
    ! location of space between the input flag and the following value(s)
    state = state_begin
    do j = 1, len_trim(line) ! loop over characters of line
      if (state .eq. state_begin) then
        if (line(j:j) .ne. " ") then
          j0 = j
          state = state_in_flag
        end if
      else if (state .eq. state_in_flag) then
        if (index("= ", line(j:j)) .gt. 0) then
          flag = set_lowercase(line(j0:j-1))
          state = state_in_sep
        end if
      else if (state .eq. state_in_sep) then
        if (index(" =", line(j:j)) .eq. 0) then
          value = line(j:)
          exit
        end if
      else
        call par_quit("Error  : Cannot read line from simulation.cfg")
      end if
    end do ! loop over characters of line

    ! Extract any input flags with no following values
    if (state == state_in_flag) then
      flag = set_lowercase(line(j0:))
    elseif (state == state_begin) then
      cycle
    end if

    ! Extract number of phases
    if (flag .eq. 'number_of_phases') then
      call parse_integer(flag, value, 1, huge(0), number_of_phases)
      ! Rewind file to start for full parsing
      rewind(file_id)
      ! Exit do loop
      exit
    else
      cycle
    end if

  end do ! loop over file

  ! With number of phases, allocate and initialize crys
  allocate(crys(number_of_phases))
  do i = 1, number_of_phases
    call crys_set_default (crys(i))
  end do

  ! Loop over the entire file and parse data line-by-line
  do ! loop over file

    ! Read the current line of the file
    read(file_id, "(a)", iostat=ioerr) line

    ! Error checking upon reading, exit loop if EOF
    if (ioerr .lt. -1) then
      call par_quit("Error  : Cannot read from simulation.cfg")
    else if (ioerr .eq. -1) then
      exit
    end if

    ! Ignore lines that are commented or empty
    if ((line(1:1) .eq. "#") .or. (len_trim(line) .eq. 0)) then
      cycle
    end if

    ! Extract input flags (flag) and value(s) (value) by checking for
    ! location of space between the input flag and the following value(s)
    state = state_begin
    do j = 1, len_trim(line) ! loop over characters of line
      if (state .eq. state_begin) then
        if (line(j:j) .ne. " ") then
          j0 = j
          state = state_in_flag
        end if
      else if (state .eq. state_in_flag) then
        if (index("= ", line(j:j)) .gt. 0) then
          flag = set_lowercase(line(j0:j-1))
          state = state_in_sep
        end if
      else if (state .eq. state_in_sep) then
        if (index(" =", line(j:j)) .eq. 0) then
          value = set_lowercase(line(j:))
          exit
        end if
      else
        call par_quit("Error  : Cannot read line from simulation.cfg")
      end if
    end do ! loop over characters of line

    ! Extract any input flags with no following values
    if (state == state_in_flag) then
      flag = set_lowercase(line(j0:))
    elseif (state == state_begin) then
      cycle
    end if

    ! Extract value(s) associated with input flag

    ! Check for general optional input
    if (flag .eq. "check_necking") then
      call parse_check_necking(value, exec)
    else if (flag .eq. "def_control_by") then
      version = "1"
      call parse_def_control_by(value, loading_options)
    else if (flag .eq. "dtime_factor") then
      call parse_real(flag, value, 1.0d0, huge(0.0d0), exec%dtime_factor)
    else if (flag .eq. "load_tol") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%load_tol_rel)
    else if (flag .eq. "load_tol_abs") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%load_tol_abs)
    else if (flag .eq. "max_eqstrain") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%max_eqstrain)
    else if (flag .eq. "max_incr") then
      call parse_integer(flag, value, 1, huge(0), exec%max_incr)
    else if (flag .eq. "max_strain") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%max_strain)
    else if (flag .eq. "max_total_time") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%max_total_time)
    else if (flag .eq. "max_iter_hard_limit") then
      call parse_integer(flag, value, 1, huge(0), exec%max_iter_hard_limit)
    else if (flag .eq. "read_ori_from_file") then
      ori_file = "simulation.ori"
    else if (flag .eq. "read_phase_from_file") then
      phase_file = "simulation.phase"
    else if (flag .eq. "restart") then
      call parse_restart(value, exec, printing)

    ! Check for printing options
    else if (flag .eq. "print") then
      call parse_print(value, printing)
    else if (flag .eq. "suppress") then
      ! Option obsolete: default is suppressed output. Simply cycle.
      cycle

    ! Check for boundary conditions options
    else if (flag .eq. "boundary_conditions") then
      call parse_boundary_conditions(value, loading_options)
    else if (flag .eq. "load_rate") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), loading_options%load_rate)
    else if (flag .eq. "loading_direction") then
      call parse_loading_direction(value, loading_options)
    else if (flag .eq. "loading_face") then
      if (myid .eq. 0) then
        write (*, '(a)') 'Warning: Option loading_face obsolete, ignoring.'
      end if
      cycle
    else if (flag .eq. "strain_rate") then
      call parse_real(flag, value, -huge(0.0d0), huge(0.0d0), loading_options%strain_rate)

    ! Check for general boundary conditions options
    else if (flag .eq. "number_of_bcs") then
      write (*, '(a)') 'Warning: Option number_of_bcs obsolete, ignoring.'
    else if (flag .eq. "set_bc") then
      call parse_set_bc(value, loading_options)

! Check for user-defined mpcs
    else if (flag .eq. "number_of_mpc") then
      write (*, '(a)') 'Warning: Option number_of_mpcs obsolete, ignoring.'
    else if (flag .eq. "set_mpc1") then
      call parse_set_mpc1(value, loading_options)
    else if (flag .eq. "set_mpc") then
      call parse_set_mpc(value, loading_options)

! Check for periodic boundary conditions
    else if (flag .eq. "periodic_bc") then
      call parse_string(value, loading_options%periodicity)
    else if (flag .eq. "imposed_strain_rate_state") then
      call parse_set_strain_rate_state(value, loading_options%imposed_strain_rate_state)

! Check for crystal parameter options
    else if (flag .eq. "number_of_phases") then
      ! Already parsed this at-the-top
      cycle
    else if (flag .eq. "phase") then
      call parse_integer(flag, value, 1, huge(0), phase)
    else if (flag .eq. "crystal_type") then
      call parse_string(value, crys(phase)%structure)
    else if (flag .eq. "m") then
      call parse_m(value, crys(phase))
    else if (flag .eq. "gammadot_0") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%gammadot_0)
    else if (flag .eq. "h_0") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%h_0)
    else if (flag .eq. "g_0") then
      call parse_g_0(value, crys(phase))
    else if (flag .eq. "g_0_bcc_112") then
      call parse_g_0_bcc_112(value, crys(phase))
    else if (flag .eq. "g_s") then
      if (g_s0_parsed) then
        call par_quit ("g_s and g_s0 are mutually exclusive.")
      end if
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%g_s0)
    else if (flag .eq. "g_s0") then
      if (g_s_parsed) then
        call par_quit ("g_s and g_s0 are mutually exclusive.")
      end if
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%g_s0)
    else if (flag .eq. "n") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%n)
    else if (flag .eq. "c11") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%c11)
    else if (flag .eq. "c12") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%c12)
    else if (flag .eq. "hardening" .or. flag .eq. "hard_type") then
      call parse_hard_type(value, crys(phase))
      ! If fcc or bcc, set c13 equal to value of c12
      if ((crys(phase)%structure .eq. "fcc") .or. &
        & (crys(phase)%structure .eq. "bcc")) then
        crys(phase)%c13 = crys(phase)%c12
      end if
    else if (flag .eq. "c13") then
      ! Check crystal type, c13 input only valid for hcp or bct
      if ((crys(phase)%structure .eq. "fcc") .or. &
        & (crys(phase)%structure .eq. "bcc")) then
        call par_quit("Error  : c13 invalid for selected crys_type.")
      else
        call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%c13)
      end if
    else if (flag .eq. "c44") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%c44)
    else if (flag .eq. "c66") then
      ! Check crystal type, c66 input only valid for bct
      if ((crys(phase)%structure .eq. "fcc") .or. &
        (crys(phase)%structure .eq. "bcc") .or. &
        (crys(phase)%structure .eq. "hcp")) then
        call par_quit("Error  : c66 invalid for selected crys_type.")
      else
        call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%c66)
      end if
    else if (flag .eq. "c_over_a") then
      ! Check crystal type, c_over_a input only valid for hcp or bct
      if ((crys(phase)%structure .eq. "fcc") .or. &
        & (crys(phase)%structure .eq. "bcc")) then
        call par_quit("Error  : c_over_a invalid for selected crys_type.")
      else
        call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%c_over_a)
      end if
    else if (flag .eq. "m_prime") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%m_prime)
    else if (flag .eq. "gammadot_s0") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%gammadot_s0)
    else if (flag .eq. "cyclic_a") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%cyclic_a)
    else if (flag .eq. "cyclic_c") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%cyclic_c)
    else if (flag .eq. "interaction_matrix_parameters") then
      call parse_interaction_matrix_parameters(value, crys(phase))
    else if (flag .eq. "a_p") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%a_p)
    else if (flag .eq. "f_p") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%f_p)
    else if (flag .eq. "b_p") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%b_p)
    else if (flag .eq. "r_p") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%r_p)
    else if (flag .eq. "b_m") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%b_m)
    else if (flag .eq. "c_p") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), crys(phase)%c_p)

 ! loading / boundary conditions --------------------------------------------

    ! number of steps
    else if (flag .eq. "number_of_steps" &
        .or. flag .eq. "number_of_strain_steps" &
        .or. flag .eq. "number_of_load_steps" &
        .or. flag .eq. "number_of_csr_load_steps" &
        .or. flag .eq. "number_of_clr_load_steps") then
      call parse_integer(flag, value, 1, huge(0), loading_options%num_steps)

      ! should go elsewhere

      call alloc_1d (loading_options%target_time, 1, loading_options%num_steps)
      call alloc_1d (loading_options%target_strain, 1, loading_options%num_steps)
      call alloc_1d (loading_options%step_target_time_incr, 1, loading_options%num_steps)

      if (flag .ne. "number_of_steps") then
        write (*, '(a,a,a)') 'Warning: "', trim(flag), '" obsolete and will be removed in a'
        write (*, '(a)') '         future version.  Use "number_of_steps" instead.'
      end if

    ! targets
    else if (flag(1:7) .eq. "target_") then

      if (version .eq. "" ) then
        version = "2"
        loading_options%def_control_by = "uniaxial"
      end if

      loading_options%target = flag(8:)

      if (version .eq. "1") then
        if (loading_options%target .eq. 'time') then
          call parse_target_time_legacy(value, loading_options)
        else if (loading_options%target .eq. 'strain') then
          call parse_target_strain_legacy(value, loading_options)
        else if (loading_options%target(1:4) .eq. 'load') then
          if (.not. allocated (loading_options%target_load)) then
            call alloc_2d (loading_options%target_load, 1, loading_options%num_steps, 1, 3)
          end if
          if (.not. allocated (loading_options%step_dt_min)) then
            call alloc_1d (loading_options%step_dt_min, 1, loading_options%num_steps)
          end if
          if (.not. allocated (loading_options%step_dt_max)) then
            call alloc_1d (loading_options%step_dt_max, 1, loading_options%num_steps)
          end if
          if (.not. allocated (loading_options%step_print)) then
            call alloc_1d_logical (loading_options%step_print, 1, loading_options%num_steps)
          end if
          call parse_target_load_legacy(value, 4, loading_options)
        else if (flag .eq. "target_csr_load") then
          if (.not. allocated (loading_options%target_load)) then
            call alloc_2d (loading_options%target_load, 1, loading_options%num_steps, 1, 3)
          end if
          if (.not. allocated (loading_options%step_dt_min)) then
            call alloc_1d (loading_options%step_dt_min, 1, loading_options%num_steps)
          end if
          if (.not. allocated (loading_options%step_dt_max)) then
            call alloc_1d (loading_options%step_dt_max, 1, loading_options%num_steps)
          end if
          if (.not. allocated (loading_options%step_print)) then
            call alloc_1d_logical (loading_options%step_print, 1, loading_options%num_steps)
          end if
          call parse_target_load_legacy(value, 6, loading_options)
        else if (flag .eq. "target_clr_load") then
          if (.not. allocated (loading_options%target_load)) then
            call alloc_2d (loading_options%target_load, 1, loading_options%num_steps, 1, 3)
          end if
          if (.not. allocated (loading_options%step_dt_min)) then
            call alloc_1d (loading_options%step_dt_min, 1, loading_options%num_steps)
          end if
          if (.not. allocated (loading_options%step_dt_max)) then
            call alloc_1d (loading_options%step_dt_max, 1, loading_options%num_steps)
          end if
          if (.not. allocated (loading_options%step_print)) then
            call alloc_1d_logical (loading_options%step_print, 1, loading_options%num_steps)
          end if
          call parse_target_load_legacy(value, 5, loading_options)
        end if


      else
        read (value, '(a)', iostat=ierr) line
        if (loading_options%num_steps .eq. 0) then
          call string_numsubstrings (line, loading_options%num_steps)
        end if

        call alloc_1d (loading_options%step_dt_min, 1, loading_options%num_steps)
        call alloc_1d (loading_options%step_target_time_incr, 1, loading_options%num_steps)
        call alloc_1d (loading_options%target_time, 1, loading_options%num_steps)
        call alloc_1d (loading_options%target_strain, 1, loading_options%num_steps)
        call alloc_2d (loading_options%target_load, 1, loading_options%num_steps, 1, 3)

        ! time target
        if (loading_options%target .eq. 'time') then
          call parse_target(value, loading_options%num_steps, loading_options%target_time)

        ! strainXX target
        else if (loading_options%target(1:6) .eq. 'strain') then
          if (loading_options%target(7:) .eq. '11') then
            loading_options%loading_direction = 1
          else if (loading_options%target(7:) .eq. '22') then
            loading_options%loading_direction = 2
          else if (loading_options%target(7:) .eq. '33') then
            loading_options%loading_direction = 3
          else
            write (message, '(a,a,a)') 'Error:   Cannot process "', trim(flag), '"'
            call par_quit(message)
          end if

          loading_options%target = "strain"
          call parse_target(value, loading_options%num_steps, loading_options%target_strain)

        ! loadX target
        else if (loading_options%target(1:4) .eq. 'load') then
          loading_options%def_control_by = "uniaxial_load_target"
          if (loading_options%target(5:5) .eq. '1' .or. loading_options%target(5:5) .eq. 'x') then
            loading_options%loading_direction = 1
          else if (loading_options%target(5:5) .eq. '2' .or. loading_options%target(5:5) .eq. 'y') then
            loading_options%loading_direction = 2
          else if (loading_options%target(5:5) .eq. '3' .or. loading_options%target(5:5) .eq. 'z') then
            loading_options%loading_direction = 3
          else
            write (message, '(a,a,a)') 'Error:   Cannot process "', trim(flag), '"'
            call par_quit(message)
          end if

          loading_options%target = "load"
          call parse_target(value, loading_options%num_steps, &
                          & loading_options%target_load(:,loading_options%loading_direction))

        else
          write (message, '(a,a,a)') 'Error:   Cannot process "', trim(flag), '"'
          call par_quit(message)
        end if
      end if

    else if (flag(1:14) .eq. 'number_of_incr') then
      call alloc_1d_int (loading_options%step_num_incr, 1, loading_options%num_steps)
      call parse_target_int(value, loading_options%num_steps, loading_options%step_num_incr)

    else if (flag .eq. 'dtime') then
      call alloc_1d (loading_options%step_dt_max, 1, loading_options%num_steps)
      call parse_target(value, loading_options%num_steps, loading_options%step_dt_max)

    else if (flag .eq. 'dtime_min') then
      call alloc_1d (loading_options%step_dt_min, 1, loading_options%num_steps)
      call parse_target(value, loading_options%num_steps, loading_options%step_dt_min)

    else if (flag .eq. 'dstrain') then
      call alloc_1d (loading_options%step_dstrain_max, 1, loading_options%num_steps)
      call parse_target(value, loading_options%num_steps, loading_options%step_dstrain_max)

    else if (flag .eq. 'print_results') then
      call alloc_1d_logical (loading_options%step_print, 1, loading_options%num_steps)
      call parse_target_logical(value, loading_options%num_steps, loading_options%step_print)

    ! --------------------------------------------------------------------------

    else if (flag .eq. "number_of_strain_rate_jumps") then
      call parse_integer(flag, value, 1, huge(0), loading_options%number_of_strain_rate_jumps)
      ! Allocate necessary arrays for strain rate jumps
      allocate (loading_options%strain_rate_jump(loading_options%number_of_strain_rate_jumps, 2))
      loading_options%strain_rate_jump = -1.0d0
    else if (flag .eq. "strain_rate_jump") then
      call parse_uniaxial_strain_rate_jump(value, loading_options)
    ! Check for triaxial constant strain rate control options
    else if (flag .eq. "min_pert_frac") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%min_pert_frac)
    else if (flag .eq. "load_tol_rel") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%load_tol_rel)
























    ! Check for triaxial constant load rate control options
    else if (flag .eq. "max_strain_incr") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), loading_options%dwell_max_strain_incr)
    else if (flag .eq. "number_of_load_rate_jumps") then
      call parse_integer(flag, value, 1, huge(0), loading_options%number_of_load_rate_jumps)
      ! Allocate necessary arrays for clr load rate jumps
      allocate (loading_options%load_rate_jump(loading_options%number_of_load_rate_jumps, 2))
      loading_options%load_rate_jump = -1.0d0
    else if (flag .eq. "load_rate_jump") then
      call parse_load_rate_jump(value, loading_options)
    else if (flag .eq. "number_of_dwell_episodes") then
      call parse_integer(flag, value, 1, huge(0), loading_options%number_of_dwell_episodes)
      ! Allocate necessary arrays for clr dwell episodes
      allocate (loading_options%dwell_episode(loading_options%number_of_dwell_episodes, 4))
      loading_options%dwell_episode = -1.0d0
    else if (flag .eq. "dwell_episode") then
      call parse_dwell_episode(value, loading_options)

    ! Check for convergence parameters options
    else if (flag .eq. "cg_max_iters") then
      call parse_integer(flag, value, 1, huge(0), exec%cg_max_iters)
    else if (flag .eq. "cg_tol") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%cg_tol)
    else if (flag .eq. "nl_max_iters") then
      call parse_integer(flag, value, 1, huge(0), exec%nl_max_iters)
    else if (flag .eq. "nl_tol_strict") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%nl_tol_strict)
    else if (flag .eq. "nl_tol_loose") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%nl_tol_loose)
    else if (flag .eq. "nl_tol_min") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%nl_tol_min)
    else if (flag .eq. "nr_tol_switch_ref") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%nr_tol_switch_ref)
    else if (flag .eq. "nr_tol_conv") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%nr_tol_conv)
    else if (flag .eq. "sx_max_iters_state") then
      call parse_integer(flag, value, 1, huge(0), exec%sx_max_iters_state)
    else if (flag .eq. "sx_max_iters_newton") then
      call parse_integer(flag, value, 1, huge(0), exec%sx_max_iters_newton)
    else if (flag .eq. "sx_tol") then
      call parse_real(flag, value, 0.0d0, huge(0.0d0), exec%sx_tol)
    else if (flag .eq. "max_bc_iter") then
      call parse_integer(flag, value, 1, huge(0), exec%max_bc_iter)

    ! Option not recognized
    else
      write(message, '(a,a,a)') "Error  : Unknown option ", trim(flag), "."
      call par_quit(message)
    end if

  end do ! loop over file

  ! Close file upon completion
  close(file_id)

  ! Perform post-parsing to set other options and check errors
  call post_parsing(crys, loading_options, exec)

  ! Write values to terminal
  call print_cfg_terminal(crys, number_of_phases)

  ! Write to terminal
  if (myid .eq. 0) then
    write (*, '(a,a,a)') "Info   :   [i] Parsed file `", "simulation.cfg", "'."
  end if

end subroutine read_config

!===============================================================================

subroutine parse_integer(flag, value, min, max, var)

  ! Arguments:
  character(len=100), intent(in) :: flag
  character(len=*), intent(in) :: value
  integer, intent(in) :: min
  integer, intent(in) :: max
  integer, intent(inout) :: var
  ! Locals:
  character(len=1000) :: message

  !-----------------------------------------------------------------------------

  read(value, *, iostat=ioerr) var
  if ((ioerr .ne. 0) .or. (var .lt. min) .or. (var .gt. max)) then
    write(message, '(a,a,a)') "Error  : parse_integer failed with option ", trim(flag), "."
    call par_quit(message)
  end if

end subroutine parse_integer

!===============================================================================

subroutine parse_real(flag, value, min, max, var)

  ! Arguments:
  character(len=100), intent(in) :: flag
  character(len=*), intent(in) :: value
  real(rk), intent(in) :: min
  real(rk), intent(in) :: max
  real(rk), intent(inout) :: var
  ! Locals:
  character(len=500) :: message

  !-----------------------------------------------------------------------------

  read(value, *, iostat=ioerr) var
  if ((ioerr .ne. 0) .or. (var .lt. min) .or. (var .gt. max)) then
    write(message, '(a,a,a)') "Error  : parse_real failed with option ", trim(flag), "."
    call par_quit(message)
  end if

end subroutine parse_real

!===============================================================================

subroutine parse_string(value, string)

  ! Arguments:
  character(len=*), intent(in) :: value
  character(len=*), intent(out) :: string

  string = trim(adjustl(value))

end subroutine parse_string

!===============================================================================

subroutine parse_check_necking(value, exec)

  ! Arguments:
  character(len=*), intent(in) :: value
  type(exec_type), intent(inout) :: exec

  !-----------------------------------------------------------------------------

  select case(trim(adjustl(value)))
    case('on', 't', 'true')
      exec%check_necking = .true.
    case('off', 'f', 'false')
      exec%check_necking = .false.
    case default
      call par_quit("Error  : Supplied value invalid for option check_necking.")
  end select

end subroutine parse_check_necking

!===============================================================================

subroutine parse_def_control_by(value, loading_options)

  ! Arguments:
  character(len=*), intent(in) :: value
  type(loading_options_type), intent(inout) :: loading_options

  !-----------------------------------------------------------------------------

  loading_options%def_control_by = trim(adjustl(value))

  if (loading_options%def_control_by .ne. 'uniaxial_load_target' .and. &
    & loading_options%def_control_by .ne. 'uniaxial_strain_target' .and. &
    & loading_options%def_control_by .ne. 'triaxial_constant_strain_rate' .and. &
    & loading_options%def_control_by .ne. 'triaxial_constant_load_rate') then
      call par_quit("Error  : Supplied value invalid for option def_control_by.")
  end if

end subroutine parse_def_control_by

!===============================================================================

subroutine parse_hard_type(value, crys)

  ! Arguments:
  character(len=*), intent(in) :: value
  type(crys_type), intent(inout) :: crys

  !-----------------------------------------------------------------------------

  crys%hardening  = trim(adjustl(value))

  if (ut_list_testelt (crys%hardening, ',', 'cyclic')) then
    crys%cyclic = .true.
  end if

  if (ut_list_testelt (crys%hardening, ',', 'precipitation')) then
    crys%precipitation = .true.
  end if

  if (ut_list_testelt (crys%hardening, ',', 'anisotropic')) then
    crys%anisotropic = .true.
  end if

  if (ut_list_testelt (crys%hardening, ',', 'saturation_evolution')) then
    crys%saturation_evolution = .true.
  end if

  if (crys%anisotropic .and. crys%cyclic) then
    call par_quit('Error  :     > Anisotropic cyclic hardening is not currently supported.')
  end if

end subroutine parse_hard_type

!===============================================================================

subroutine parse_restart(value, exec, printing)

  ! Arguments:
  character(len=*), intent(in) :: value
  type(exec_type), intent(inout) :: exec
  type(printing_type), intent(inout) :: printing

  !-----------------------------------------------------------------------------

  select case(trim(adjustl(value)))
    case('append')
      call par_quit("Error  : Supplied value obsolete for option restart.")
    case('on')
      exec%restart = .true.
      printing%restart = .true.
      printing%restart_file_handling = "restart"
  end select

end subroutine parse_restart

!===============================================================================

subroutine parse_print(value, printing)

  ! Arguments:
  character(len=*), intent(in) :: value
  type(printing_type), intent(inout) :: printing

  integer :: i, num_inputs
  character(len=255) :: variable

  !-----------------------------------------------------------------------------

  call string_numsubstrings (value, num_inputs)

  do i = 1, num_inputs

    call string_substring (value, i, variable)

    select case(trim(adjustl(variable)))
      case('convergence')
        printing%print_conv = .true.
      case('coo')
        printing%print_coo = .true.
      case('crss')
        printing%print_crss = .true.
      case('defrate')
        printing%print_defrate = .true.
      case('defrate-eq', 'defrate_eq')
        printing%print_defrate_eq = .true.
      case('defrate-pl', 'defrate_pl')
        printing%print_defrate_pl = .true.
      case('defrate-pl-eq', 'defrate_pl_eq')
        printing%print_defrate_pl_eq = .true.
      case('disp')
        printing%print_disp = .true.
      case('forces')
        printing%print_forces = .true.
      case('ori')
        printing%print_ori = .true.
      case('restart')
        printing%print_restart = .true.
      case('slip')
        printing%print_slip = .true.
      case('sliprate')
        printing%print_sliprate = .true.
      case('rss')
        printing%print_rss = .true.
      case('spinrate')
        printing%print_spinrate = .true.
      case('rotrate')
        printing%print_rotrate = .true.
      case('rotrate_spin')
        printing%print_rotrate_spin = .true.
      case('rotrate_slip')
        printing%print_rotrate_slip = .true.
      case('strain')
        printing%print_strain = .true.
      case('strain-el', 'strain_el')
        printing%print_strain_el = .true.
      case('strain-el-eq', 'strain_el_eq')
        printing%print_strain_el_eq = .true.
      case('strain-eq', 'strain_eq')
        printing%print_strain_eq = .true.
      case('strain-pl', 'strain_pl')
        printing%print_strain_pl = .true.
      case('strain-pl-eq', 'strain_pl_eq')
        printing%print_strain_pl_eq = .true.
      case('stress')
        printing%print_stress = .true.
      case('stress-eq', 'stress_eq')
        printing%print_stress_eq = .true.
      case('vel')
        printing%print_vel = .true.
      case('velgrad')
        printing%print_velgrad = .true.
      case('work')
        printing%print_work = .true.
      case('work-pl', 'work_pl')
        printing%print_work_pl = .true.
      case('workrate')
        printing%print_workrate = .true.
      case('workrate-pl', 'workrate_pl')
        printing%print_workrate_pl = .true.
      case default
        write (*,*) variable
        call par_quit("Error  : Supplied value invalid for option print.")
    end select

  end do

end subroutine parse_print

!===============================================================================

subroutine parse_boundary_conditions(value, loading_options)

  ! Arguments:
  character(len=*), intent(in) :: value
  type(loading_options_type), intent(inout) :: loading_options

  !-----------------------------------------------------------------------------

  select case(trim(adjustl(value)))
    case('user_defined', 'general')
      loading_options%boundary_conditions = "general"
    case('uniaxial_grip')
      loading_options%boundary_conditions = "uniaxial_grip"
    case('uniaxial_symmetry')
      loading_options%boundary_conditions = "uniaxial_symmetry"
    case('uniaxial_minimal')
      loading_options%boundary_conditions = "uniaxial_minimal"
    case('triaxial')
      loading_options%boundary_conditions = "triaxial"
    case('periodic','PBC')
      loading_options%boundary_conditions = "periodic"
    case default
      loading_options%boundary_conditions = trim(adjustl(value))
  end select

end subroutine parse_boundary_conditions

!===============================================================================

subroutine parse_loading_direction(value, loading_options)

  ! Arguments:
  character(len=*), intent(in) :: value
  type(loading_options_type), intent(inout) :: loading_options

  !-----------------------------------------------------------------------------

  select case(trim(adjustl(value)))
    case('x', '+x', '1')
      loading_options%loading_direction = 1
    case('y', '+y', '2')
      loading_options%loading_direction = 2
    case('z', '+z', '3')
      loading_options%loading_direction = 3
    case default
      call par_quit("Error  : Supplied value invalid for option loading_direction.")
  end select

end subroutine parse_loading_direction

!===============================================================================

subroutine parse_set_bc(value, loading_options)

  ! Arguments:
  character(len=*), intent(in) :: value
  type(loading_options_type), intent(inout) :: loading_options

  integer :: num_dofs
  character(len=256) :: bc_var, bc_subvar, bc_method

  !-----------------------------------------------------------------------------

  call string_substring (value, 1, bc_var)

  if (bc_var .eq. "vel") then

    call string_numsubstrings (value, num_dofs)

    num_dofs = (num_dofs - 2) / 2

    do i = 1, num_dofs

      loading_options%num_bcs = loading_options%num_bcs + 1

      call string_substring (value, 2, loading_options%bc_nset(loading_options%num_bcs))
      loading_options%bc_var(loading_options%num_bcs) = bc_var

      call string_substring (value, 2 * i + 1, loading_options%bc_dir(loading_options%num_bcs))
      call string_substring_real (value, 2 * i + 2, loading_options%bc_vel(loading_options%num_bcs))

      ! Error message : x, y, z, vx, vy, vz not allowed names for nsets to avoid misundertanding
      if ((loading_options%bc_nset(loading_options%num_bcs) .eq. 'x') .or. &
        & (loading_options%bc_nset(loading_options%num_bcs) .eq. 'y') .or. &
        & (loading_options%bc_nset(loading_options%num_bcs) .eq. 'z') .or. &
        & (loading_options%bc_nset(loading_options%num_bcs) .eq. 'vx') .or. &
        & (loading_options%bc_nset(loading_options%num_bcs) .eq. 'vy') .or. &
      & (loading_options%bc_nset(loading_options%num_bcs) .eq. 'vz')) then
        call par_quit("Error  : set_bc: Invalid nset name.")
      end if

      ! Error message : x, y, z, vx, vy, vz not allowed names for nsets to avoid misundertanding
      if ((loading_options%bc_dir(loading_options%num_bcs) .ne. 'x') .and. &
        & (loading_options%bc_dir(loading_options%num_bcs) .ne. 'y') .and. &
        & (loading_options%bc_dir(loading_options%num_bcs) .ne. 'z') .and. &
        & (loading_options%bc_dir(loading_options%num_bcs) .ne. 'vx') .and. &
        & (loading_options%bc_dir(loading_options%num_bcs) .ne. 'vy') .and. &
        & (loading_options%bc_dir(loading_options%num_bcs) .ne. 'vz')) then
        call par_quit("Error  : set_bc: Invalid dir value.")
      end if
    end do

  else if (bc_var .eq. "strainrate" .or. bc_var .eq. "velgrad") then

    call string_numsubstrings (value, num_dofs)
    
    call string_substring (value, num_dofs, bc_subvar)

    if (bc_subvar .eq. "periodic" .or. bc_subvar .eq. "pbc") then
  
      bc_method = "periodic"
      
      if (mod(num_dofs - 2, 2) .eq. 1) then
        call string_substring (value, num_dofs, bc_method)
      end if

      num_dofs = (num_dofs - 2) / 2

      do i = 1, num_dofs
        loading_options%imposed_strain_rate_state%nb_comp = &
          & loading_options%imposed_strain_rate_state%nb_comp + 1

        call string_substring(value, 2 * i, loading_options%imposed_strain_rate_state%comp &
          & (loading_options%imposed_strain_rate_state%nb_comp))

        call string_substring_real(value, 2 * i + 1, loading_options%imposed_strain_rate_state%val &
          & (loading_options%imposed_strain_rate_state%nb_comp))

        if (any(loading_options%imposed_strain_rate_state%comp &
          & (loading_options%imposed_strain_rate_state%nb_comp) .eq. &
          & loading_options%imposed_strain_rate_state%comp &
          & (1:loading_options%imposed_strain_rate_state%nb_comp - 1))) then
          call par_quit("Error  : set_bc: Two specified values for the same component.")
        end if

        ! Error message : x, y, z, vx, vy, vz not allowed names for nsets to avoid misundertanding
        if ((loading_options%imposed_strain_rate_state%comp &
          & (loading_options%imposed_strain_rate_state%nb_comp) .ne. '11') .and. &
          & (loading_options%imposed_strain_rate_state%comp &
          & (loading_options%imposed_strain_rate_state%nb_comp) .ne. '12') .and. &
          & (loading_options%imposed_strain_rate_state%comp &
          & (loading_options%imposed_strain_rate_state%nb_comp) .ne. '13') .and. &
          & (loading_options%imposed_strain_rate_state%comp& 
          & (loading_options%imposed_strain_rate_state%nb_comp) .ne. '21') .and. &
          & (loading_options%imposed_strain_rate_state%comp &
          & (loading_options%imposed_strain_rate_state%nb_comp) .ne. '22') .and. &
          & (loading_options%imposed_strain_rate_state%comp &
          & (loading_options%imposed_strain_rate_state%nb_comp) .ne. '23') .and. &
          & (loading_options%imposed_strain_rate_state%comp &
          & (loading_options%imposed_strain_rate_state%nb_comp) .ne. '31') .and. &
          & (loading_options%imposed_strain_rate_state%comp &
          & (loading_options%imposed_strain_rate_state%nb_comp) .ne. '32') .and. &
          & (loading_options%imposed_strain_rate_state%comp &
          & (loading_options%imposed_strain_rate_state%nb_comp) .ne. '33')) then
          call par_quit("Error  : set_bc: Invalid dir value.")
        end if
        
       ! if (((loading_options%imposed_strain_rate_state%comp &
       !   & (loading_options%imposed_strain_rate_state%nb_comp) .ne. '11') .and. &
       ! (loading_options%imposed_strain_rate_state%comp &
       !   & (loading_options%imposed_strain_rate_state%nb_comp) .ne. '22') .and. &
       ! (loading_options%imposed_strain_rate_state%comp &
       !   & (loading_options%imposed_strain_rate_state%nb_comp) .ne. '33')) &
       !   & .and. (loading_options%imposed_strain_rate_state%val &
       !   & (loading_options%imposed_strain_rate_state%nb_comp) .ne. 0.0d0)) then
       !     call par_quit("Error  : set_bc: Nonzero imposed shear not allowed.")
       ! end if

      end do

    end if

    else 

    call string_numsubstrings (value, num_dofs)
    bc_method = "minimal"
    if (mod(num_dofs - 1, 2) .eq. 1) then
      call string_substring (value, num_dofs, bc_method)
    end if

    num_dofs = (num_dofs - 1) / 2

    do i = 1, num_dofs

      loading_options%num_bcs = loading_options%num_bcs + 1

      loading_options%bc_var(loading_options%num_bcs) = bc_var
      loading_options%bc_type(loading_options%num_bcs) = bc_method

      call string_substring (value, 2 * i, loading_options%bc_dir(loading_options%num_bcs))
      call string_substring_real (value, 2 * i + 1, loading_options%bc_vel(loading_options%num_bcs))

      ! Error message : x, y, z, vx, vy, vz not allowed names for nsets to avoid misundertanding
      if ((loading_options%bc_dir(loading_options%num_bcs) .ne. '11') .and. &
        & (loading_options%bc_dir(loading_options%num_bcs) .ne. '12') .and. &
        & (loading_options%bc_dir(loading_options%num_bcs) .ne. '13') .and. &
        & (loading_options%bc_dir(loading_options%num_bcs) .ne. '21') .and. &
        & (loading_options%bc_dir(loading_options%num_bcs) .ne. '22') .and. &
        & (loading_options%bc_dir(loading_options%num_bcs) .ne. '23') .and. &
        & (loading_options%bc_dir(loading_options%num_bcs) .ne. '31') .and. &
        & (loading_options%bc_dir(loading_options%num_bcs) .ne. '32') .and. &
        & (loading_options%bc_dir(loading_options%num_bcs) .ne. '33')) then
        write (*,*) loading_options%bc_dir(loading_options%num_bcs)
        call par_quit("Error  : set_bc: Invalid dir value.")
      end if

    end do

  end if
 
end subroutine parse_set_bc

!===============================================================================

subroutine parse_set_mpc(value, loading_options)

  ! Arguments
  character(len=*), intent(in) :: value
  type(loading_options_type), intent(inout) :: loading_options

  integer :: id

  loading_options%num_mpcs = loading_options%num_mpcs + 1
  id = loading_options%num_mpcs

  call string_numsubstrings (value, loading_options%general_mpc(id)%nbdof)

  loading_options%general_mpc(id)%nbdof = (loading_options%general_mpc(id)%nbdof - 1) / 3

  allocate (loading_options%general_mpc(id)%dir(loading_options%general_mpc(id)%nbdof))
  allocate (loading_options%general_mpc(id)%coeff(loading_options%general_mpc(id)%nbdof))
  allocate (loading_options%general_mpc(id)%offset(loading_options%general_mpc(id)%nbdof))

  ! Reading of arguments in the commande line set_bc
  if (loading_options%general_mpc(id)%nbdof .eq. 1) then
    read (value, *, iostat=ioerr) &
      & loading_options%general_mpc(id)%nset, &
      & loading_options%general_mpc(id)%dir(1), &
      & loading_options%general_mpc(id)%coeff(1), &
      & loading_options%general_mpc(id)%offset(1)

  else if (loading_options%general_mpc(id)%nbdof .eq. 2) then
    read (value, *, iostat=ioerr)&
      & loading_options%general_mpc(id)%nset, &
      & loading_options%general_mpc(id)%dir(1), &
      & loading_options%general_mpc(id)%coeff(1), &
      & loading_options%general_mpc(id)%offset(1), &
      & loading_options%general_mpc(id)%dir(2), &
      & loading_options%general_mpc(id)%coeff(2), &
      & loading_options%general_mpc(id)%offset(2)

  else if (loading_options%general_mpc(id)%nbdof .eq. 3) then
    read (value, *, iostat=ioerr)&
      & loading_options%general_mpc(id)%nset, &
      & loading_options%general_mpc(id)%dir(1), &
      & loading_options%general_mpc(id)%coeff(1), &
      & loading_options%general_mpc(id)%offset(1), &
      & loading_options%general_mpc(id)%dir(2), &
      & loading_options%general_mpc(id)%coeff(2), &
      & loading_options%general_mpc(id)%offset(2), &
      & loading_options%general_mpc(id)%dir(3), &
      & loading_options%general_mpc(id)%coeff(3), &
      & loading_options%general_mpc(id)%offset(3)

  else
    call par_quit("Error  : parse_set_mpc: Invalid number of degrees of freedom.")

  end if

  ! reading error checking
  if (ioerr .ne. 0) then
    call par_quit("Error  : parse_set_mpc: Error reading mpcs.")
  end if

end subroutine parse_set_mpc

!===============================================================================

subroutine parse_set_mpc1(value, loading_options)

  ! Arguments
  character(len=*), intent(in) :: value
  type(loading_options_type), intent(inout) :: loading_options

  integer :: id
  character(len=256) :: tmp

  loading_options%num_mpcs = loading_options%num_mpcs + 1
  id = loading_options%num_mpcs

  call string_numsubstrings (value, loading_options%general_mpc(id)%nbdof)

  loading_options%general_mpc(id)%nbdof = loading_options%general_mpc(id)%nbdof - 2

  allocate (loading_options%general_mpc(id)%dir(loading_options%general_mpc(id)%nbdof))
  allocate (loading_options%general_mpc(id)%coeff(loading_options%general_mpc(id)%nbdof))
  allocate (loading_options%general_mpc(id)%offset(loading_options%general_mpc(id)%nbdof))
  loading_options%general_mpc(id)%coeff = 1.0d0
  loading_options%general_mpc(id)%offset = 0.0d0

  call string_substring (value, 1, tmp)
  call string_substring (value, 2, loading_options%general_mpc(id)%nset)

  ! Reading of arguments in the commande line set_bc
  if (loading_options%general_mpc(id)%nbdof .eq. 1) then
    read (value, *, iostat=ioerr) &
      & tmp, &
      & loading_options%general_mpc(id)%nset, &
      & loading_options%general_mpc(id)%dir(1)

  else if (loading_options%general_mpc(id)%nbdof .eq. 2) then
    read (value, *, iostat=ioerr)&
      & tmp, &
      & loading_options%general_mpc(id)%nset, &
      & loading_options%general_mpc(id)%dir(1), &
      & loading_options%general_mpc(id)%dir(2)

  else if (loading_options%general_mpc(id)%nbdof .eq. 3) then
    read (value, *, iostat=ioerr)&
      & tmp, &
      & loading_options%general_mpc(id)%nset, &
      & loading_options%general_mpc(id)%dir(1), &
      & loading_options%general_mpc(id)%dir(2), &
      & loading_options%general_mpc(id)%dir(3)

  else
    call par_quit("Error  : parse_set_mpc1: Invalid number of degrees of freedom.")

  end if

  ! reading error checking
  if (ioerr .ne. 0) then
    call par_quit("Error  : parse_set_mpc1: Error reading mpcs.")
  end if

end subroutine parse_set_mpc1

subroutine parse_set_strain_rate_state(value, imposed_strain_rate_state)

  ! Arguments:
  character(len=*), intent(in) :: value
  type(imposed_strain_rate_state_type), intent(inout) :: imposed_strain_rate_state
  ! Locals:
  integer :: nvals
  integer :: i

  !-----------------------------------------------------------------------------

  ! Find number of inputs in "value" by counting spaces
  i = 1
  nvals = count((/(value(i:i), i=1, len_trim(value))/) .eq. " ") + 1
  imposed_strain_rate_state%nb_comp = floor(real(nvals) / 2)

  ! Reading of arguments in the commande line set_bc
  if (imposed_strain_rate_state%nb_comp .eq. 1) then
    read (value, *, iostat=ioerr) &
      & imposed_strain_rate_state%comp(1), &
      & imposed_strain_rate_state%val(1)

  else if (imposed_strain_rate_state%nb_comp .eq. 2) then
    read (value, *, iostat=ioerr) &
      & imposed_strain_rate_state%comp(1), &
      & imposed_strain_rate_state%val(1), &
      & imposed_strain_rate_state%comp(2), &
      & imposed_strain_rate_state%val(2)

  else if (imposed_strain_rate_state%nb_comp .eq. 3) then
    read (value, *, iostat=ioerr) &
      & imposed_strain_rate_state%comp(1), &
      & imposed_strain_rate_state%val(1), &
      & imposed_strain_rate_state%comp(2), &
      & imposed_strain_rate_state%val(2), &
      & imposed_strain_rate_state%comp(3), &
      & imposed_strain_rate_state%val(3)
  
  else if (imposed_strain_rate_state%nb_comp .eq. 4) then
    read (value, *, iostat=ioerr) &
      & imposed_strain_rate_state%comp(1), &
      & imposed_strain_rate_state%val(1), &
      & imposed_strain_rate_state%comp(2), &
      & imposed_strain_rate_state%val(2), &
      & imposed_strain_rate_state%comp(3), &
      & imposed_strain_rate_state%val(3),&
      & imposed_strain_rate_state%comp(4), &
      & imposed_strain_rate_state%val(4)
  
  else if (imposed_strain_rate_state%nb_comp .eq. 5) then
    read (value, *, iostat=ioerr) &
      & imposed_strain_rate_state%comp(1), &
      & imposed_strain_rate_state%val(1), &
      & imposed_strain_rate_state%comp(2), &
      & imposed_strain_rate_state%val(2), &
      & imposed_strain_rate_state%comp(3), &
      & imposed_strain_rate_state%val(3),&
      & imposed_strain_rate_state%comp(4), &
      & imposed_strain_rate_state%val(4),&
      & imposed_strain_rate_state%comp(5), &
      & imposed_strain_rate_state%val(5)
  
  else if (imposed_strain_rate_state%nb_comp .eq. 6) then
    read (value, *, iostat=ioerr) &
      & imposed_strain_rate_state%comp(1), &
      & imposed_strain_rate_state%val(1), &
      & imposed_strain_rate_state%comp(2), &
      & imposed_strain_rate_state%val(2), &
      & imposed_strain_rate_state%comp(3), &
      & imposed_strain_rate_state%val(3),&
      & imposed_strain_rate_state%comp(4), &
      & imposed_strain_rate_state%val(4),&
      & imposed_strain_rate_state%comp(5), &
      & imposed_strain_rate_state%val(5),&
      & imposed_strain_rate_state%comp(6), &
      & imposed_strain_rate_state%val(6)

  else
    call par_quit("Error  : parse_set_strain_rate_state: Invalid number of components.")

  end if

  ! reading error checking
  if (ioerr .ne. 0) then
    call par_quit("Error  : parse_set_strain_rate_state: Error reading imposed strain rate state.")
  end if

end subroutine parse_set_strain_rate_state

!===============================================================================

subroutine parse_m(value, my_crys)

  ! Arguments:
  character(len=*), intent(in) :: value
  type(crys_type), intent(inout) :: my_crys
  ! Locals:
  integer :: i
  integer :: nvals

  !---------------------------------------------------------------------------

  ! If fcc or bcc, read in singular m value
  if ((my_crys%structure .eq. "fcc") .or. &
      & (my_crys%structure .eq. "bcc")) then

    ! Find number of inputs in "value" by counting spaces
    i = 1
    nvals = count((/(value(i:i), i=1, len_trim(value))/) .eq. " ") + 1

    ! Based on number of values input, proceed accordingly
    select case (nvals)

      ! If one value, read in the isotropic `m'
      case (1)
        read(value, *, iostat=ioerr) my_crys%m
        if ((ioerr .ne. 0) .or. (my_crys%m .le. 0.0d0)) then
          call par_quit("Error  : parse_m failed.")
        end if

      case default
        call par_quit("Error  : parse_m: invalid number of inputs.")

    end select

  ! If hcp, try to read 1 or 3 values
  else if (my_crys%structure .eq. "hcp") then

    ! Find number of inputs in "value" by counting spaces
    i = 1
    nvals = count((/(value(i:i), i=1, len_trim(value))/) .eq. " ") + 1

    select case (nvals)

      ! If one value, read in the isotropic `m'
      case (1)
        read(value, *, iostat=ioerr) my_crys%m
        if ((ioerr .ne. 0) .or. (my_crys%m .le. 0.0d0)) then
          call par_quit("Error  : parse_m failed.")
        end if

      ! If three values, read in the anisotropic `m'
      case (3)
        my_crys%use_aniso_m = .true.
        read (value, *, iostat=ioerr) my_crys%aniso_m(1:3)
        if ((ioerr .ne. 0) .or. (any(my_crys%aniso_m(1:3) .le. 0.0d0))) then
          call par_quit("Error  : parse_m failed.")
        end if

      case default
        call par_quit("Error  : parse_m: invalid number of inputs.")

    end select

  ! If bct, try to read 1 or 10 values
  else if (my_crys%structure .eq. "bct") then

    ! Find number of inputs in "value" by counting spaces
    i = 1
    nvals = count((/(value(i:i), i=1, len_trim(value))/) .eq. " ") + 1

    select case (nvals)

      ! If one value, read in the isotropic `m'
      case (1)
        read (value, *, iostat=ioerr) my_crys%m
        if ((ioerr .ne. 0) .or. (my_crys%m .le. 0.0d0)) then
          call par_quit("Error  : parse_m failed.")
        end if

      ! If ten values, read in the anisotropic `m'
      case (10)
        my_crys%use_aniso_m = .true.
        read (value, *, iostat=ioerr) my_crys%aniso_m(1:10)
        if ((ioerr .ne. 0) .or. (any(my_crys%aniso_m(1:10) .le. 0.0d0))) then
          call par_quit("Error  : parse_m failed.")
        end if

      case default
        call par_quit("Error  : parse_m: invalid number of inputs.")

    end select

  else
    call par_quit("Error  : parse_m: invalid crystal phase.")

  end if

end subroutine parse_m

!===============================================================================

subroutine parse_g_0(value, my_crys)

  ! Arguments:
  character(len=*), intent(in) :: value
  type(crys_type), intent(inout) :: my_crys
  ! Locals:
  integer :: i
  integer :: nvals
  real(rk) :: g_0_cubic(12)
  real(rk) :: g_0_hcp(18)
  real(rk) :: g_0_bct(32)

  !-----------------------------------------------------------------------------

  ! If fcc or bcc, try to read 1 or 12 values
  if ((my_crys%structure .eq. "fcc") .or. &
      & (my_crys%structure .eq. "bcc")) then

    ! Find number of inputs in "value" by counting spaces
    i = 1
    nvals = count((/(value(i:i), i=1, len_trim(value))/) .eq. " ") + 1
    my_crys%hratio_num = nvals

    ! Based on number of values input, proceed accordingly
    select case (nvals)

      ! Isotropic g_0 (1 value input)
      case(1)
        read(value, *, iostat=ioerr) my_crys%g_0
        if ((ioerr .ne. 0) .or. (my_crys%g_0 .le. 0.0d0)) then
          call par_quit("Error  : parse_g_0 failed.")
        end if

      ! Anisotropic g_0 (12 values input)
      case(12)
        read(value, *, iostat=ioerr) g_0_cubic(1:12)
        if ((ioerr .ne. 0) .or. (any(g_0_cubic(1:12) .le. 0.0d0))) then
          call par_quit("Error  : parse_g_0 failed.")
        end if
        my_crys%g_0 = g_0_cubic(1)
        my_crys%hratio_cubic(1:12) = g_0_cubic(1:12) / my_crys%g_0

      ! Otherwise not permitted
      case default
          call par_quit("Error  : parse_g_0: invalid number of inputs.")

    end select

  ! If hcp, try to read 3 or 18 values
  else if (my_crys%structure .eq. "hcp") then

    ! Find number of inputs in "value" by counting spaces
    i = 1
    nvals = count((/(value(i:i), i=1, len_trim(value))/) .eq. " ") + 1
    my_crys%hratio_num = nvals

    ! Based on number of values input, proceed accordingly
    select case (nvals)

      ! Basal, prismatic, and pyramidal g_0 (3 value inputs)
      case(3)
        read(value, *, iostat=ioerr) g_0_hcp(1:3)
        if ((ioerr .ne. 0) .or. (any(g_0_hcp(1:3) .le. 0.0d0))) then
          call par_quit("Error  : parse_g_0 failed.")
        end if
        my_crys%g_0 = g_0_hcp(1)
        my_crys%hratio_hcp(1:3) = 1.0d0
        my_crys%hratio_hcp(4:6) = g_0_hcp(2) / g_0_hcp(1)
        my_crys%hratio_hcp(7:18) = g_0_hcp(3) / g_0_hcp(1)

      ! Fully anisotropic g_0 (18 values input)
      case(18)
        read(value, *, iostat=ioerr) g_0_hcp(1:18)
        if ((ioerr .ne. 0) .or. (any(g_0_hcp(1:18) .le. 0.0d0))) then
          call par_quit("Error  : parse_g_0 failed.")
        end if
        my_crys%g_0 = g_0_hcp(1)
        my_crys%hratio_hcp(1:18) = g_0_hcp(1:18) / my_crys%g_0

      ! Otherwise not permitted
      case default
        call par_quit("Error  : parse_g_0: invalid number of inputs.")

    end select

  ! If bct, try to read 10 or 32 values
  else if (my_crys%structure .eq. "bct") then

    ! Find number of inputs in "value" by counting spaces
    i = 1
    nvals = count((/(value(i:i), i=1, len_trim(value))/) .eq. " ") + 1
    my_crys%hratio_num = nvals

    ! Based on number of values input, proceed accordingly
    select case (nvals)

      ! Per-family g_0 (10 value inputs)
      case(10)
        read(value, *, iostat=ioerr) g_0_bct(1:10)
        if ((ioerr .ne. 0) .or. (any(g_0_bct(1:10) .le. 0.0d0))) then
          call par_quit("Error  : parse_g_0 failed.")
        end if
        my_crys%g_0 = g_0_bct(1)
        my_crys%hratio_bct(1:2) = 1.0d0
        my_crys%hratio_bct(3:4) = g_0_bct(2) / g_0_bct(1)
        my_crys%hratio_bct(5:6) = g_0_bct(3) / g_0_bct(1)
        my_crys%hratio_bct(7:10) = g_0_bct(4) / g_0_bct(1)
        my_crys%hratio_bct(11:12) = g_0_bct(5) / g_0_bct(1)
        my_crys%hratio_bct(13:16) = g_0_bct(6) / g_0_bct(1)
        my_crys%hratio_bct(17:18) = g_0_bct(7) / g_0_bct(1)
        my_crys%hratio_bct(19:20) = g_0_bct(8) / g_0_bct(1)
        my_crys%hratio_bct(21:24) = g_0_bct(9) / g_0_bct(1)
        my_crys%hratio_bct(25:32) = g_0_bct(10) / g_0_bct(1)

      ! Fully anisotropic g_0 (32 values input)
      case(32)
        read(value, *, iostat=ioerr) g_0_bct(1:32)
        if ((ioerr .ne. 0) .or. (any(g_0_bct(1:32) .le. 0.0d0))) then
          call par_quit("Error  : parse_g_0 failed.")
        end if
        my_crys%g_0 = g_0_bct(1)
        my_crys%hratio_bct(1:32) = g_0_bct(1:32) / my_crys%g_0

      ! Otherwise not permitted
      case default
        call par_quit("Error  : parse_g_0: invalid number of inputs.")

    end select

  ! Otherwise, incorrect phase assignment
  else
    call par_quit("Error  : parse_g_0: invalid crystal type.")

  end if

end subroutine parse_g_0

!===============================================================================

subroutine parse_g_0_bcc_112(value, my_crys)

  ! Arguments:
  character(len=*), intent(in) :: value
  type(crys_type), intent(inout) :: my_crys
  ! Locals:
  integer :: i
  integer :: nvals
  real(rk) :: g_0_cubic(12)

  !-----------------------------------------------------------------------------

  ! If bcc only, try to read 1 or 12 values
  if (my_crys%structure .eq. "bcc") then

    ! Find number of inputs in "value" by counting spaces
    i = 1
    nvals = count((/(value(i:i), i=1, len_trim(value))/) .eq. " ") + 1
    my_crys%hratio_num_112 = nvals

    ! Based on number of values input, proceed accordingly
    select case (nvals)

      ! Isotropic g_0 (1 value input)
      case(1)
        read(value, *, iostat=ioerr) my_crys%g_0_bcc_112
        my_crys%hratio_cubic_112 = my_crys%g_0_bcc_112 / my_crys%g_0
        if ((ioerr .ne. 0) .or. (my_crys%g_0_bcc_112 .le. 0.0d0)) then
          call par_quit("Error  : parse_g_0_112 failed.")
        end if

      ! Anisotropic g_0 (12 values input)
      case(12)
        read(value, *, iostat=ioerr) g_0_cubic(1:12)
        if ((ioerr .ne. 0) .or. (any(g_0_cubic(1:12) .le. 0.0d0))) then
          call par_quit("Error  : parse_g_0_112 failed.")
        end if
        my_crys%g_0_bcc_112 = g_0_cubic(1)
        my_crys%hratio_cubic_112(1:12) = g_0_cubic(1:12) / my_crys%g_0

      ! Otherwise not permitted
      case default
        call par_quit("Error  : parse_g_0_112: invalid number of inputs.")

    end select

  ! Otherwise, incorrect phase assignment
  else
    call par_quit("Error  : parse_g_0_112: invalid crystal type.")

  end if

end subroutine parse_g_0_bcc_112

!===============================================================================

subroutine parse_interaction_matrix_parameters(value, my_crys)

  ! Arguments:
  character(len=*), intent(in) :: value
  type(crys_type), intent(inout) :: my_crys
  ! Locals:
  integer :: i
  integer :: nvals

  !-----------------------------------------------------------------------------

  ! If fcc, try to read 2 or 5 values
  if (my_crys%structure .eq. "fcc") then

    ! Find number of inputs in "value" by counting spaces
    i = 1
    nvals = count((/(value(i:i), i=1, len_trim(value))/) .eq. " ") + 1
    my_crys%interaction_matrix_parameters_num = nvals

    select case (nvals)

      ! Fully anisotropic interaction (2 values input)
      case (2)
        read(value, *, iostat=ioerr) my_crys%interaction_matrix_parameters(1:2)
        if ((ioerr .ne. 0) .or. (any(my_crys%interaction_matrix_parameters(1:2) .le. 0.0d0))) then
          call par_quit("Error  : parse_interaction_matrix_parameters failed.")
        end if

      ! Co-planar interaction (5 values input)
      case (5)
        read(value, *, iostat=ioerr) my_crys%interaction_matrix_parameters(1:5)
        if ((ioerr .ne. 0) .or. (any(my_crys%interaction_matrix_parameters(1:5) .le. 0.0d0))) then
          call par_quit("Error  : parse_interaction_matrix_parameters failed.")
        end if

      case default
        call par_quit("Error  : parse_interaction_matrix_parameters: invalid number of inputs.")

    end select

  ! If bcc, try to read 2 or 7 values
  else if (my_crys%structure .eq. "bcc") then

    ! Find number of inputs in "value" by counting spaces
    i = 1
    nvals = count((/(value(i:i), i=1, len_trim(value))/) .eq. " ") + 1
    my_crys%interaction_matrix_parameters_num = nvals

    select case (nvals)

      ! Fully anisotropic interaction (2 values input)
      case (2)
        read(value, *, iostat=ioerr) my_crys%interaction_matrix_parameters(1:2)
        if ((ioerr .ne. 0) .or. (any(my_crys%interaction_matrix_parameters(1:2) .le. 0.0d0))) then
          call par_quit("Error  : parse_interaction_matrix_parameters failed.")
        end if

      ! Co-planar interaction (7 values input)
      case (7)
        read (value, *, iostat=ioerr) my_crys%interaction_matrix_parameters(1:7)
        if ((ioerr .ne. 0) .or. (any(my_crys%interaction_matrix_parameters(1:7) .le. 0.0d0))) then
          call par_quit("Error  : parse_interaction_matrix_parameters failed.")
        end if

      case default
        call par_quit("Error  : parse_interaction_matrix_parameters: invalid number of inputs.")

    end select

  ! If hcp, try to read 2 or 8 values
  else if (my_crys%structure .eq. "hcp") then

    ! Find number of inputs in "value" by counting spaces
    i = 1
    nvals = count((/(value(i:i), i=1, len_trim(value))/) .eq. " ") + 1
    my_crys%interaction_matrix_parameters_num = nvals

    select case (nvals)

      ! Fully anisotropic interaction (2 values input)
      case (2)
        read (value, *, iostat=ioerr) my_crys%interaction_matrix_parameters(1:2)
        if ((ioerr .ne. 0) .or. (any(my_crys%interaction_matrix_parameters(1:2) .le. 0.0d0))) then
          call par_quit("Error  : parse_interaction_matrix_parameters failed.")
        end if

      ! Co-planar interaction (8 values input)
      case (8)
        read (value, *, iostat=ioerr) my_crys%interaction_matrix_parameters(1:8)
        if ((ioerr .ne. 0) .or. (any(my_crys%interaction_matrix_parameters(1:8) .le. 0.0d0))) then
          call par_quit("Error  : parse_interaction_matrix_parameters failed.")
        end if

      case default
        call par_quit("Error  : parse_interaction_matrix_parameters: invalid number of inputs.")

    end select

  ! If bct, try to read 2 or 11 values
  else if (my_crys%structure .eq. "bct") then

    ! Find number of inputs in "value" by counting spaces
    i = 1
    nvals = count((/(value(i:i), i=1, len_trim(value))/) .eq. " ") + 1
    my_crys%interaction_matrix_parameters_num = nvals

    select case (nvals)

      ! Fully anisotropic interaction (2 values input)
      case (2)
        read (value, *, iostat=ioerr) my_crys%interaction_matrix_parameters(1:2)
        if ((ioerr .ne. 0) .or. (any(my_crys%interaction_matrix_parameters(1:2) .le. 0.0d0))) then
          call par_quit("Error  : parse_interaction_matrix_parameters failed.")
        end if

      ! Co-planar interaction (11 values input)
      case (11)
        read (value, *, iostat=ioerr) my_crys%interaction_matrix_parameters(1:11)
        if ((ioerr .ne. 0) .or. (any(my_crys%interaction_matrix_parameters(1:11) .le. 0.0d0))) then
          call par_quit("Error  : parse_interaction_matrix_parameters failed.")
        end if

      case default
        call par_quit("Error  : parse_interaction_matrix_parameters: invalid number of inputs.")

    end select

  else
    call par_quit("Error  : parse_interaction_matrix_parameters: invalid crystal type.")

  end if

end subroutine parse_interaction_matrix_parameters

!===============================================================================

subroutine parse_target_time_legacy(value, loading_options)

  character(len=*), intent(in) :: value
  type(loading_options_type), intent(inout) :: loading_options

  character(len=15) :: temp

  loading_options%target = "time"

  loading_options%curr_target = loading_options%curr_target + 1

  ! test allocation -----------
  if (.not. allocated (loading_options%step_num_incr)) then
    call alloc_1d_int (loading_options%step_num_incr, 1, loading_options%num_steps)
  end if
  if (loading_options%curr_target .gt. loading_options%num_steps) then
    call par_quit("Error  : Supplied value for num_steps invalid.")
  end if
  if (.not. allocated (loading_options%step_print)) then
    call alloc_1d_logical (loading_options%step_print, 1, loading_options%num_steps)
  end if
  ! end test allocation -------

  read (value, *, iostat=ioerr) &
    & loading_options%target_time(loading_options%curr_target), &
    & loading_options%step_num_incr(loading_options%curr_target), temp
  if (ioerr .ne. 0) then
    call par_quit("Error  : parse_target_time failed.")
  end if

  if (temp .eq. 'print_data') then
    loading_options%step_print(loading_options%curr_target) = .true.
  else if (temp .eq. 'suppress_data') then
    loading_options%step_print(loading_options%curr_target) = .false.
  else
    call par_quit("Error  : parse_target_time: invalid print control.")
  end if

end subroutine parse_target_time_legacy

!===============================================================================

subroutine parse_target_strain_legacy(value, loading_options)

  character(len=*), intent(in) :: value
  type(loading_options_type), intent(inout) :: loading_options

  character(len=15) :: temp

  loading_options%target = "strain"

  ! test allocation -----------
  if (.not. allocated (loading_options%step_num_incr)) then
    call alloc_1d_int (loading_options%step_num_incr, 1, loading_options%num_steps)
  end if
  if (loading_options%curr_target .gt. loading_options%num_steps) then
    call par_quit("Error  : Supplied value for num_steps invalid.")
  end if
  if (.not. allocated (loading_options%step_print)) then
    call alloc_1d_logical (loading_options%step_print, 1, loading_options%num_steps)
  end if
  ! end test allocation -------

  loading_options%curr_target = loading_options%curr_target + 1

  if (loading_options%curr_target .gt. loading_options%num_steps) then
    call par_quit("Error  : Supplied value for num_steps invalid.")
  end if

  read (value, *, iostat=ioerr) &
    & loading_options%target_strain(loading_options%curr_target), &
    & loading_options%step_num_incr(loading_options%curr_target), temp
  if (ioerr .ne. 0) then
    call par_quit("Error  : parse_uniaxial_target_strain failed.")
  end if

  if (temp .eq. 'print_data') then
    loading_options%step_print(loading_options%curr_target) = .true.
  else if (temp .eq. 'suppress_data') then
    loading_options%step_print(loading_options%curr_target) = .false.
  else
    call par_quit("Error  : parse_uniaxial_target_strain: invalid print control.")
  end if

end subroutine parse_target_strain_legacy

!===============================================================================

subroutine parse_target (value, size, array)

  character(len=*), intent(in) :: value
  integer, intent(in) :: size
  real(rk), intent(out) :: array(size)

  array = -1232768262 ! random value to check for initialization
  read (value, *, iostat=ioerr) array(1:size)
  do i = 2, size
    if (array(i) .eq. -1232768262) then
      array(i) = array(i - 1)
    endif
  end do

end subroutine parse_target

!===============================================================================

subroutine parse_target_int (value, size, array)

  character(len=*), intent(in) :: value
  integer, intent(in) :: size
  integer, intent(out) :: array(size)

  array = -1232768262 ! random value to check for initialization
  read (value, *, iostat=ioerr) array
  do i = 2, size
    if (array(i) .eq. -1232768262) then
      array(i) = array(i - 1)
    endif
  end do

end subroutine parse_target_int

!===============================================================================

subroutine parse_target_logical (value, size, array)

  character(len=*), intent(in) :: value
  integer, intent(in) :: size
  logical, intent(out) :: array(size)
  integer :: tmp(size)

  call parse_target_int (value, size, tmp)
  do i = 1, size
    if (tmp(i) .eq. 0) then
      array(i) = .false.
    else if (tmp(i) .eq. 1) then
      array(i) = .true.
    end if
  end do

end subroutine parse_target_logical

!===============================================================================

subroutine parse_uniaxial_strain_rate_jump(value, loading_options)

  character(len=*), intent(in) :: value
  type(loading_options_type), intent(inout) :: loading_options

  !-----------------------------------------------------------------------------

  loading_options%curr_strain_jump = loading_options%curr_strain_jump + 1

  if (loading_options%curr_strain_jump .gt. loading_options%number_of_strain_rate_jumps) then
    call par_quit("Error  : Supplied value for number_of_strain_rate_jumps invalid.")
  end if

  read (value, *, iostat=ioerr) &
    & loading_options%strain_rate_jump(loading_options%curr_strain_jump, 1:2)
  if (ioerr .ne. 0) then
    call par_quit("Error  : parse_uniaxial_strain_rate_jump failed.")
  end if

end subroutine parse_uniaxial_strain_rate_jump

!===============================================================================

subroutine parse_target_load_legacy(value, length, loading_options)

  character(len=*), intent(in) :: value
  integer, intent(in) :: length
  type(loading_options_type), intent(inout) :: loading_options

  character(len=15) :: temp

  loading_options%target = "load"

  loading_options%curr_target = loading_options%curr_target + 1

  if (loading_options%curr_target .gt. loading_options%num_steps) then
    call par_quit("Error  : Supplied value for 'num_steps' is invalid.")
  end if

  if (length .eq. 4) then
    read (value, *, iostat=ioerr) &
      & loading_options%target_load(loading_options%curr_target, 1:1), &
      & loading_options%step_dt_max(loading_options%curr_target), &
      & loading_options%step_dt_min(loading_options%curr_target), &
      & temp

  elseif (length .eq. 5) then
    read (value, *, iostat=ioerr) &
      & loading_options%target_load(loading_options%curr_target, 1:3), &
      & loading_options%step_target_time_incr(loading_options%curr_target), &
      & temp

  elseif (length .eq. 6) then
    read (value, *, iostat=ioerr) &
      & loading_options%target_load(loading_options%curr_target, 1:3), &
      & loading_options%step_dt_max(loading_options%curr_target), &
      & loading_options%step_dt_min(loading_options%curr_target), &
      & temp
  end if

  if (ioerr .ne. 0) then
    call par_quit("Error  : parse_target_load failed.")
  end if

  if (temp .eq. 'print_data') then
    loading_options%step_print(loading_options%curr_target) = .true.
  elseif (temp .eq. 'suppress_data') then
    loading_options%step_print(loading_options%curr_target) = .false.
  else
    call par_quit("Error  : parse_target_load: invalid print control.")
  end if

end subroutine parse_target_load_legacy

!===============================================================================

subroutine parse_load_rate_jump(value, loading_options)

  character(len=*), intent(in) :: value
  type(loading_options_type), intent(inout) :: loading_options

  loading_options%curr_load_rate_jump = loading_options%curr_load_rate_jump + 1

  if (loading_options%curr_load_rate_jump .gt. loading_options%number_of_load_rate_jumps) then
    call par_quit("Error  : Supplied value for 'number_of_load_rate_jumps' is invalid.")
  end if

  read (value, *, iostat=ioerr) &
    & loading_options%load_rate_jump(loading_options%curr_load_rate_jump, 1:2)
  if (ioerr .ne. 0) then
    call par_quit("Error  : parse_load_rate_jump failed.")
  end if

end subroutine parse_load_rate_jump

!===============================================================================

subroutine parse_dwell_episode(value, loading_options)

  character(len=*), intent(in) :: value
  type(loading_options_type), intent(inout) :: loading_options

  character(len=15) :: temp

  loading_options%curr_dwell_episode = loading_options%curr_dwell_episode + 1

  if (loading_options%curr_dwell_episode .gt. loading_options%number_of_dwell_episodes) then
    call par_quit("Error  : Supplied value for 'number_of_dwell_episodes' is invalid.")
  end if

  read (value, *, iostat=ioerr) &
    & loading_options%dwell_episode(loading_options%curr_dwell_episode, 1:3), temp
  if (ioerr .ne. 0) then
    call par_quit("Error  : parse_dwell_episode failed.")
  end if

  if (temp .eq. 'print_data') then
    loading_options%dwell_episode(loading_options%curr_dwell_episode, 4) = 0
  else if (temp .eq. 'suppress_data') then
    loading_options%dwell_episode(loading_options%curr_dwell_episode, 4) = 1
  else
    call par_quit("Error  : parse_dwell_episode: invalid print control.")
  end if

end subroutine parse_dwell_episode

!===============================================================================

!> Set options based on input and check for errors
subroutine post_parsing(my_crys, loading_options, exec)

  type(crys_type), intent(inout) :: my_crys (:)
  type(loading_options_type), intent(inout) :: loading_options
  type(exec_type), intent(in) :: exec

  integer :: iphase

  if (.not. allocated (loading_options%step_print)) then
    call alloc_1d_logical (loading_options%step_print, 1, loading_options%num_steps)
    loading_options%step_print = .true.
  end if

  ! Error checking of proper def_control_by for max_bc_iter
  if (exec%max_bc_iter .ne. 10) then
    if ((loading_options%def_control_by .ne. 'triaxial_constant_strain_rate') .or. &
      & (loading_options%def_control_by .ne. 'triaxial_constant_load_rate')) then
      ! Do nothing, this is proper
    else
      call par_quit("Error  : max_bc_iter invalid for selected def_control_by.")
    end if
  end if

  ! Error checking of proper def_control_by for load_tol_abs
  if (exec%load_tol_abs .ne. 0.1d0) then
    if ((loading_options%def_control_by .eq. 'triaxial_constant_strain_rate') .or. &
      & (loading_options%def_control_by .eq. 'triaxial_constant_load_rate')) then
      ! Do nothing, this is proper
    else
      call par_quit("Error  : load_tol_abs invalid for selected def_control_by.")
    end if
  end if

  ! Error checking of proper def_control_by for max_eqstrain
  if (exec%max_eqstrain .ne. 0.2d0) then
    if ((loading_options%def_control_by .ne. 'triaxial_constant_strain_rate') .or. &
      & (loading_options%def_control_by .ne. 'triaxial_constant_load_rate')) then
      ! Do nothing, this is proper
    else
      call par_quit("Error  : max_eqstrain invalid for selected def_control_by.")
    end if
  end if

  ! Error checking of proper def_control_by for max_strain
  if (exec%max_strain .ne. 0.2d0) then
    if ((loading_options%def_control_by .ne. 'triaxial_constant_strain_rate') .or. &
      & (loading_options%def_control_by .ne. 'triaxial_constant_load_rate')) then
      ! Do nothing, this is proper
    else
      call par_quit("Error  : max_strain invalid for selected def_control_by.")
    end if
  end if

  ! Error checking of my_crys%m_prime and my_crys%gammadot_s0
  do iphase = 1, size(my_crys)
    if ((my_crys(iphase)%m_prime .gt. 0.0d0) .and. &
      & (my_crys(iphase)%gammadot_s0 .gt. 0.0d0)) then
      my_crys(iphase)%saturation_evolution = .true.
    else if (((my_crys(iphase)%m_prime .gt. 0.0d0) .and. &
      & (my_crys(iphase)%gammadot_s0 .le. 0.0d0)) .or. &
      & ((my_crys(iphase)%m_prime .le. 0.0d0) .and. &
      & (my_crys(iphase)%gammadot_s0 .gt. 0.0d0))) then
      call par_quit("Error  : Insufficient saturation strength evolution parameters.")
    end if
  end do

  ! Error checking of my_crys%cyclic_a and my_crys%cyclic_c
  do iphase = 1, size(my_crys)

    if (my_crys(iphase)%cyclic .and. &
      & (my_crys(iphase)%cyclic_a .gt. 0.0d0) .and. &
      & (my_crys(iphase)%cyclic_c .le. 0.0d0)) then
      call par_quit("Error  : Insufficient cyclic hardening parameters.")
    else if (my_crys(iphase)%cyclic .and. &
      & (my_crys(iphase)%cyclic_a .le. 0.0d0) .and. &
      & (my_crys(iphase)%cyclic_c .gt. 0.0d0)) then
      call par_quit("Error  : Insufficient cyclic hardening parameters.")
    else if (.not.my_crys(iphase)%cyclic .and. &
      & ((my_crys(iphase)%cyclic_a .gt. 0.0d0) .or. &
      & (my_crys(iphase)%cyclic_c .gt. 0.0d0))) then
      call par_quit("Error  : Cyclic parameters invalid for selected hardening.")
    end if
  end do

  ! Error checking of precipitate values
  do iphase = 1, size(my_crys)
    ! Check for "correct" supply of variables and set options
    if (all((/my_crys(iphase)%a_p, my_crys(iphase)%f_p, my_crys(iphase)%b_p, &
      & my_crys(iphase)%r_p, my_crys(iphase)%b_m, my_crys(iphase)%c_p/) .gt. 0.0d0)) then
      my_crys(iphase)%precipitation = .true.
    else if (all((/my_crys(iphase)%a_p, my_crys(iphase)%f_p, my_crys(iphase)%b_p, &
      & my_crys(iphase)%r_p/) .gt. 0.0d0) .and. &
      all((/my_crys(iphase)%b_m, my_crys(iphase)%c_p/) .le. 0.0d0)) then
      my_crys(iphase)%precipitation_cutting = .true.
    else if (all((/my_crys(iphase)%f_p, my_crys(iphase)%b_p, my_crys(iphase)%r_p, &
      my_crys(iphase)%b_m, my_crys(iphase)%c_p/) .gt. 0.0d0) .and. &
      (my_crys(iphase)%a_p .le. 0.0d0)) then
      my_crys(iphase)%precipitation_bowing = .true.
    end if
    ! Next check for "incorrect" supply of input
    if (any((/my_crys(iphase)%a_p, my_crys(iphase)%f_p, my_crys(iphase)%b_p, &
      & my_crys(iphase)%r_p, my_crys(iphase)%b_m, my_crys(iphase)%c_p/) .gt. 0.0d0)) then
      ! Check for cutting-only values
      if (all((/my_crys(iphase)%b_m, my_crys(iphase)%c_p/) .le. 0.0d0) .and. &
        & any((/my_crys(iphase)%a_p, my_crys(iphase)%f_p, &
        & my_crys(iphase)%b_p, my_crys(iphase)%r_p/) .le. 0.0d0)) then
        call par_quit("Error  : Insufficient precipitate cutting parameters.")
      ! Next check for bowing-only values
      else if ((my_crys(iphase)%a_p .le. 0.0d0) .and. &
        & any((/my_crys(iphase)%f_p, my_crys(iphase)%b_p, &
        & my_crys(iphase)%r_p, my_crys(iphase)%b_m, my_crys(iphase)%c_p/) .le. 0.0d0)) then
        call par_quit("Error  : Insufficient precipitate bowing parameters.")
      ! Next check for cutting and bowing
      else if (all((/my_crys(iphase)%a_p, my_crys(iphase)%b_m, my_crys(iphase)%c_p/) .gt. 0.0d0) .and. &
        & any((/my_crys(iphase)%a_p, my_crys(iphase)%f_p, my_crys(iphase)%b_p, &
        & my_crys(iphase)%r_p, my_crys(iphase)%b_m, my_crys(iphase)%c_p/) .le. 0.0d0)) then
        call par_quit("Error  : Insufficient precipitate strengthening parameters.")
      end if
    end if
  end do

  ! Error checking of proper def_control_by for number_of_strain_rate_jumps
  if ((loading_options%def_control_by .eq. 'triaxial_constant_load_rate') .and. &
    & (loading_options%number_of_strain_rate_jumps .gt. 0)) then
    call par_quit("Error  : number_of_strain_rate_jumps invalid for selected def_control_by.")
  end if

  ! Error checking of proper def_control_by for max_strain_incr
  if ((loading_options%dwell_max_strain_incr .ne. 0.001d0) .and. &
    & (loading_options%def_control_by .ne. 'triaxial_constant_load_rate')) then
    call par_quit("Error  : max_strain_incr invalid for selected def_control_by.")
  end if

  ! Error checking of proper def_control_by for number_of_load_rate_jumps
  if ((loading_options%number_of_load_rate_jumps .gt. 0) .and. &
    & (loading_options%def_control_by .ne. 'triaxial_constant_load_rate')) then
    call par_quit("Error  : number_of_load_rate_jumps invalid for selected def_control_by.")
  end if

  ! Error checking of proper def_control_by for number_of_dwell_episodes
  if ((loading_options%number_of_dwell_episodes .gt. 0) .and. &
    & (loading_options%def_control_by .ne. 'triaxial_constant_load_rate')) then
    call par_quit("Error  : number_of_dwell_episodes invalid for selected def_control_by.")
  end if

end subroutine post_parsing

!===============================================================================

subroutine print_cfg_terminal(my_crys, number_of_phases)

  ! Print out crystal parameters to the console for user confirmation.

  ! Arguments:
  ! options: General options type
  ! my_crys: Crystal parameters type
  type(crys_type), intent(in) :: my_crys (:)
  ! Locals:
  ! iphase: Looping index over number_of_phases
  integer :: iphase

  !-----------------------------------------------------------------------------

  if (myid .eq. 0) then

    write (*, '(a)') 'Info   :   - Material parameters:'
    write (*, '(a,i0)') 'Info   :     > Number of phases: ', number_of_phases

    do iphase = 1, size(my_crys)

    write (*, '(a,i0,a,a)') 'Info   :     > phase ', &
        & iphase, ' - crystal type:   ', my_crys(iphase)%structure

    if (my_crys(iphase)%use_aniso_m .eqv. .false.) then
      write (*, '(a,e14.6)') 'Info   :     > m:             ', &
          & my_crys(iphase)%m

    else if ((my_crys(iphase)%use_aniso_m .eqv. .true.) .and. &
        & (my_crys(iphase)%structure .eq. "hcp")) then
      write (*, '(a,e14.6)') 'Info   :     > m (basal):     ', &
          & my_crys(iphase)%aniso_m(1)
      write (*, '(a,e14.6)') 'Info   :     > m (pris.):     ', &
          & my_crys(iphase)%aniso_m(2)
      write (*, '(a,e14.6)') 'Info   :     > m (pyr.):      ', &
          & my_crys(iphase)%aniso_m(3)

    else if ((my_crys(iphase)%use_aniso_m .eqv. .true.) .and. &
        & (my_crys(iphase)%structure .eq. "bct")) then
        do i = 1, 10
          write (*, '(a,i0,a,e14.6)') 'Info   :     > m (mode ', i, '):   ', &
            & my_crys(iphase)%aniso_m(i)
        end do
    end if

    write (*, '(a,e14.6)') 'Info   :     > gammadot_0:    ', &
        & my_crys(iphase)%gammadot_0
    write (*, '(a,e14.6)') 'Info   :     > h_0:           ', &
        & my_crys(iphase)%h_0

    if ((my_crys(iphase)%structure .eq. "fcc") .or. &
        & ((my_crys(iphase)%structure .eq. "bcc") .and. &
        & (my_crys(iphase)%g_0_bcc_112 .le. 0.0d0))) then
      if (my_crys(iphase)%hratio_num .eq. 1) then
        write (*, '(a,e14.6)') 'Info   :     > &
            &g_0:           ', my_crys(iphase)%g_0

      else
        do i = 1, 12
          write (*, '(a,i0,a,e14.6)') 'Info   :     > g0 (', i, '):       ', &
            & my_crys(iphase)%hratio_cubic(1)*my_crys(iphase)%g_0
        end do
      end if

    else if ((my_crys(iphase)%structure .eq. "bcc") .and. &
        & (my_crys(iphase)%g_0_bcc_112 .gt. 0.0d0)) then
      if (my_crys(iphase)%hratio_num .eq. 1) then
        write (*, '(a,e14.6)') 'Info   :     > &
            &g_0:           ', my_crys(iphase)%g_0

      else
        do i = 1, 12
          write (*, '(a,i0,a,e14.6)') 'Info   :     > g0 (', i, '):       ', &
            & my_crys(iphase)%hratio_cubic(1)*my_crys(iphase)%g_0
        end do
      end if

      if (my_crys(iphase)%hratio_num_112 .eq. 1) then
        write (*, '(a,e14.6)') 'Info   :     > &
            &g_0_112:       ', my_crys(iphase)%g_0_bcc_112

      else
        do i = 1, 12
          write (*, '(a,i0,a,e14.6)') 'Info   :     > &
              &g_0_112 (', i, '):   ', &
              & my_crys(iphase)%hratio_cubic_112(i)* &
              & my_crys(iphase)%g_0
        end do
      end if

    else if (my_crys(iphase)%structure .eq. "hcp") then
      if (my_crys(iphase)%hratio_num .eq. 3) then
        write (*, '(a,e14.6)') 'Info   :     > &
            &g_0 (Basal):   ', &
            & my_crys(iphase)%hratio_hcp(1)* &
            & my_crys(iphase)%g_0
        write (*, '(a,e14.6)') 'Info   :     > &
            &g_0 (Pris.):   ', &
            & my_crys(iphase)%hratio_hcp(4)* &
            & my_crys(iphase)%g_0
        write (*, '(a,e14.6)') 'Info   :     > &
            &g_0 (Pyr.):    ', &
            & my_crys(iphase)%hratio_hcp(7)* &
            & my_crys(iphase)%g_0

      else
        do i = 1, 18
          write (*, '(a,i0,a,e14.6)') 'Info   :     > &
            &g_0_112 (', i, '):   ', &
            & my_crys(iphase)%hratio_hcp(i)
        end do
      end if

    else if (my_crys(iphase)%structure .eq. "bct") then
      if (my_crys(iphase)%hratio_num .eq. 10) then
        do i = 1, 25, 2
          write (*, '(a,i0,a,e14.6)') 'Info   :     > &
            &g_0_112 (', i, '):   ', &
            & my_crys(iphase)%hratio_bct(i)*my_crys(iphase)%g_0
        end do

      else
        do i = 1, 32
          write (*, '(a,i0,a,e14.6)') 'Info   :     > &
            &g_0 (', i, '):   ', &
            & my_crys(iphase)%hratio_bct(i)*my_crys(iphase)%g_0
        end do
      end if
    end if

    write (*, '(a,e14.6)') 'Info   :     > g_s0:          ', &
        & my_crys(iphase)%g_s0

    if (my_crys(iphase)%saturation_evolution .eqv. .true.) then
      write (*, '(a,e14.6)') 'Info   :     > m_prime:       ', &
          & my_crys(iphase)%m_prime
      write (*, '(a,e14.6)') 'Info   :     > gammadot_s0:   ', &
          & my_crys(iphase)%gammadot_s0
    end if

    write (*, '(a,e14.6)') 'Info   :     > n:             ', &
        & my_crys(iphase)%n
    write (*, '(a,e14.6)') 'Info   :     > c11:           ', my_crys(iphase)%c11
    write (*, '(a,e14.6)') 'Info   :     > c12:           ', my_crys(iphase)%c12

    if ((my_crys(iphase)%structure .eq. "hcp") .or. &
        & (my_crys(iphase)%structure .eq. "bct")) then
      write (*, '(a,e14.6)') 'Info   :     > c13:           ', my_crys(iphase)%c13
    end if

    write (*, '(a,e14.6)') 'Info   :     > c44:           ', my_crys(iphase)%c44

    if (my_crys(iphase)%structure .eq. "bct") then
      write (*, '(a,e14.6)') 'Info   :     > c66:           ', my_crys(iphase)%c66
    end if

    ! Print out precipitate hardening parameters - if available.

    if (my_crys(iphase)%precipitation .eqv. .true.) then
      write (*, '(a,e14.6)') 'Info   :     > a_p:           ', &
          & my_crys(iphase)%a_p
      write (*, '(a,e14.6)') 'Info   :     > f_p:           ', &
          & my_crys(iphase)%f_p
      write (*, '(a,e14.6)') 'Info   :     > r_p:           ', &
          & my_crys(iphase)%r_p
      write (*, '(a,e14.6)') 'Info   :     > b_p:           ', &
          & my_crys(iphase)%b_p
      write (*, '(a,e14.6)') 'Info   :     > b_m:           ', &
          & my_crys(iphase)%b_m
      write (*, '(a,e14.6)') 'Info   :     > c_p:           ', &
          & my_crys(iphase)%c_p
    end if

    if (my_crys(iphase)%precipitation_cutting .eqv. .true.) then
      write (*, '(a,e14.6)') 'Info   :     > a_p:           ', &
          & my_crys(iphase)%a_p
      write (*, '(a,e14.6)') 'Info   :     > f_p:           ', &
          & my_crys(iphase)%f_p
      write (*, '(a,e14.6)') 'Info   :     > r_p:           ', &
          & my_crys(iphase)%r_p
      write (*, '(a,e14.6)') 'Info   :     > b_p:           ', &
          & my_crys(iphase)%b_p
    end if

    if (my_crys(iphase)%precipitation_bowing .eqv. .true.) then
      write (*, '(a,e14.6)') 'Info   :     > f_p:           ', &
          & my_crys(iphase)%f_p
      write (*, '(a,e14.6)') 'Info   :     > r_p:           ', &
          & my_crys(iphase)%r_p
      write (*, '(a,e14.6)') 'Info   :     > b_m:           ', &
          & my_crys(iphase)%b_m
      write (*, '(a,e14.6)') 'Info   :     > c_p:           ', &
          & my_crys(iphase)%c_p
    end if

    ! Print out cyclic hardening parameters - if available.

    if (ut_list_testelt (my_crys(iphase)%hardening, ',', 'cyclic')) then
      write (*, '(a,e14.6)') 'Info   :     > cyclic_a:      ', &
          & my_crys(iphase)%cyclic_a
      write (*, '(a,e14.6)') 'Info   :     > cyclic_c:      ', &
          & my_crys(iphase)%cyclic_c
    end if

    ! Print out anisotropic hardening parameters - if available.

    if (ut_list_testelt (my_crys(iphase)%hardening, ',', 'anisotropic')) then
      select case (my_crys(iphase)%structure)

        case ("fcc")
          num = 4

        case ("bcc")
          num = 6

        case ("hcp")
          num = 7

        case ("bct")
          num = 10

      end select

      write (*, '(a,e14.6)') 'Info   :     >&
          & diag:          ', my_crys(iphase)%interaction_matrix_parameters(1)
      do i = 1, num
        write (*, '(a,i0,a,e14.6)') 'Info   :     > h', i, ':   ', &
          & my_crys(iphase)%interaction_matrix_parameters(i + 1)
      end do
    end if
  end do
  end if

end subroutine print_cfg_terminal

end module read_input_cfg_mod
