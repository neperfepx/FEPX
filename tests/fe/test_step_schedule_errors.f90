! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Error cases of step_schedule_mod: each case writes a configuration to case.cfg in the working
!> directory, reads it and builds its schedule, which must end with a fatal error whose message
!> is checked by tests/core/check_run.cmake (list in tests/CMakeLists.txt).
!>
!> Source : FEPX (ROADMAP steps 5.1, 8.1, 8.2a, 8.2c).
!>
!> Cases (a bcc material, one or two steps):
!>   rate_away           strain rate -1e-3 towards the target 0.001
!>   reversal_rate_away  step 2 back from 0.002 to 0.001 with a positive rate
!>   zero_change         step 2 targets the strain step 1 ended with
!>   time_not_after      a target_time before the end time of the previous step
!>   dstrain_no_axis     general conditions with time targets only and dstrain
!>   dstrain_zero_rate   dstrain in a step with strain rate 0 (a hold to a time target)
!>   general_zero_vel    general conditions with a strain target and only zero velocities
!>   max_incr_step       one step needing more increments than max_incr
!>   max_incr_total      two steps that together need more increments than max_incr
!>   max_total_time      steps that end after max_total_time
!>   max_strain          a strain target beyond max_strain
!>   max_strain_time     a time target whose end strain exceeds max_strain
!>   dtime_min_dstrain   a load target whose dtime_min exceeds its largest increment
!>                       dstrain / |rate| (step 8.1)
!>   load_equal          two consecutive equal load targets (not more than 2 load_tol apart)
!>   load_first_close    a first load target within load_tol of 0 (step 8.2a, Q-T16)
!>   load_targets_close  consecutive load targets 0.15 N apart, load_tol 0.1 (Q-T16)
!>   stress_targets_close  consecutive stress targets 5 MPa apart, stress_tol 3 (Q-T16)
!>   load_away_previous  a load target below the previous one with a positive rate
!>   rate_load_away      a load rate -2 N/s to the load target 100 N from 0 (step 8.2c)
!>   rate_away_previous  a load rate 2 N/s to 50 N after a step to 100 N (step 8.2c)
!>   rate_targets_close  load-rate steps to 100 and 100.15 N, load_tol 0.1 (step 8.2c, Q-T16)
!>   rate_open_away      a load rate -1 N/s to the strain target 0.001 (step 8.2c, Q-T8 (i))
!>
!> Locals:
!>   case_name   character(32)          name of the case
!>   loading     character(64)(:)       the loading lines of the case
!>   rate        real(rk)               general conditions' strain rate passed to the build [1/s]
!>   config      type(config_t)         the configuration
!>   sched       type(step_schedule_t)  the schedule (never complete)
!>   material    character(64)(7)       the material lines
program test_step_schedule_errors

  use kinds_mod, only: rk
  use config_mod, only: config_t, config_read
  use step_schedule_mod, only: step_schedule_t, step_schedule_build

  implicit none

  character(len=32) :: case_name
  character(len=64), allocatable :: loading(:)
  real(rk) :: rate
  type(config_t) :: config
  type(step_schedule_t) :: sched
  character(len=64), parameter :: material(7) = [character(len=64) :: 'number_of_phases 1', &
    & 'phase 1', 'crystal_type bcc', 'model elastic', 'c11 231400', 'c12 134700', 'c44 116400']

  call get_command_argument(1, case_name)

  ! The loading lines of the case (presets unless stated) and the general strain rate.
  allocate (loading(0))
  rate = 0.0_rk
  select case (trim(case_name))
  case ('rate_away')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 1', 'strain_rate 33 -1.0e-3', 'target_strain 33 0.001', &
      & 'number_of_incrs 1']
  case ('reversal_rate_away')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 2', 'strain_rate 33 1.0e-3 1.0e-3', 'target_strain 33 0.002 0.001', &
      & 'number_of_incrs 1 1']
  case ('zero_change')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 2', 'strain_rate 33 1.0e-3 1.0e-3', 'target_strain 33 0.002 0.002', &
      & 'number_of_incrs 1 1']
  case ('time_not_after')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 2', 'strain_rate 33 1.0e-3 1.0e-3', 'target_strain 33 0.002 -', &
      & 'target_time - 1.5', 'number_of_incrs 1 1']
  case ('dstrain_no_axis')
    loading = [character(len=64) :: 'boundary_conditions general', 'set_bc vel z1 z 1.0', &
      & 'number_of_steps 1', 'target_time 1.0', 'dstrain 0.001']
  case ('dstrain_zero_rate')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 2', 'strain_rate 33 1.0e-3 0.0', 'target_strain 33 0.002 -', &
      & 'target_time - 5.0', 'dstrain 0.001 0.001']
  case ('dtime_min_dstrain')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 1', 'strain_rate 33 1.0e-3', 'target_load 33 100.0', &
      & 'dstrain 1.0e-4', 'dtime_min 0.5']
  case ('load_equal')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 2', 'strain_rate 33 1.0e-3 1.0e-3', 'target_load 33 100.0 100.0', &
      & 'dtime 0.1 0.1', 'dtime_min 0.01 0.01']
  case ('load_first_close')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 1', 'strain_rate 33 1.0e-3', 'target_load 33 0.05', &
      & 'dtime 0.1', 'dtime_min 0.01']
  case ('load_targets_close')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 2', 'strain_rate 33 1.0e-3 1.0e-3', 'target_load 33 100.0 100.15', &
      & 'dtime 0.1 0.1', 'dtime_min 0.01 0.01']
  case ('stress_targets_close')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 2', 'strain_rate 33 1.0e-3 -1.0e-3', 'target_stress 33 50.0 45.0', &
      & 'dtime 0.1 0.1', 'dtime_min 0.01 0.01', 'stress_tol 3.0']
  case ('load_away_previous')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 2', 'strain_rate 33 1.0e-3 1.0e-3', 'target_load 33 100.0 50.0', &
      & 'dtime 0.1 0.1', 'dtime_min 0.01 0.01']
  case ('rate_load_away')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 1', 'load_rate 33 -2.0', 'target_load 33 100.0', 'number_of_incrs 2']
  case ('rate_away_previous')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 2', 'load_rate 33 2.0 2.0', 'target_load 33 100.0 50.0', &
      & 'number_of_incrs 2 2']
  case ('rate_targets_close')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 2', 'load_rate 33 2.0 2.0', 'target_load 33 100.0 100.15', &
      & 'number_of_incrs 2 2']
  case ('rate_open_away')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 1', 'load_rate 33 -1.0', 'target_strain 33 0.001', 'dtime 0.1', &
      & 'dtime_min 0.01']
  case ('general_zero_vel')
    loading = [character(len=64) :: 'boundary_conditions general', 'set_bc vel z1 z 0.0', &
      & 'number_of_steps 1', 'target_strain 33 0.001', 'number_of_incrs 1']
  case ('max_incr_step')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 1', 'strain_rate 33 1.0e-3', 'target_strain 33 0.002', &
      & 'dtime 0.1', 'max_incr 10']
  case ('max_incr_total')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 2', 'strain_rate 33 1.0e-3 1.0e-3', 'target_strain 33 0.002 0.003', &
      & 'number_of_incrs 6 5', 'max_incr 10']
  case ('max_total_time')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 2', 'strain_rate 33 1.0e-3 1.0e-3', 'target_strain 33 0.002 0.003', &
      & 'number_of_incrs 1 1', 'max_total_time 2.5']
  case ('max_strain')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 2', 'strain_rate 33 1.0e-3 -1.0e-3', 'target_strain 33 0.002 -0.3', &
      & 'number_of_incrs 1 1', 'max_strain 0.2']
  case ('max_strain_time')
    loading = [character(len=64) :: 'boundary_conditions uniaxial_minimal', &
      & 'number_of_steps 1', 'strain_rate 33 1.0e-2', 'target_time 30.0', &
      & 'number_of_incrs 1', 'max_strain 0.2']
  case default
    error stop 'test_step_schedule_errors: unknown case'
  end select

  ! Write, read and schedule.
  call write_file('case.cfg', [material, loading])
  call config_read(config, 'case.cfg')
  call step_schedule_build(sched, config, rate)

contains

  !> Write lines to a file (trailing blanks removed).
  !>
  !> Arguments:
  !>   name   (in)  character        file name
  !>   text   (in)  character(*)(:)  the lines
  !> Locals:
  !>   u, i   integer                unit, line index
  subroutine write_file(name, text)

    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: text(:)

    integer :: u, i

    open (newunit=u, file=name, status='replace', action='write')
    do i = 1, size(text)
      write (u, '(a)') trim(text(i))
    end do
    close (u)

  end subroutine write_file

end program test_step_schedule_errors
