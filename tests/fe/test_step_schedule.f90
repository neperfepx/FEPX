! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of step_schedule_mod: step times, strains, strain rates, velocity factors and
!> increments built from configuration files, against values worked out by hand (module header
!> of step_schedule_mod, ROADMAP D44).
!>
!> Source : FEPX (ROADMAP step 5.1); load targets and v2's completion and next-increment rules
!>          (step 8.1, dev/theory/addendum_load_targets.md §2); the landing rule, false
!>          position and the separation of targets (step 8.2a,
!>          dev/theory/addendum_triaxial_load_control.md §3.2, §3.6).
!>
!> Locals:
!>   config     type(config_t)            a configuration
!>   sched      type(step_schedule_t)     its schedule
!>   unit, ios  integer                   log file unit, read status
!>   i          integer                   index
!>   line       character(256)            a line of the log file
!>   log_all    character(:), allocatable the log's texts joined (wrapped messages)
!>   warned     logical                   the log file holds the expected warning
!>   material   character(64)(7)          the material lines shared by every file
!>   dt_dwell   real(rk)                  a dwell increment [s]
!>   last_dwell logical                   it ends the dwell
program test_step_schedule

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk
  use log_mod, only: log_set_unit
  use config_mod, only: config_t, config_read, target_by_load, target_by_stress
  use bc_iteration_mod, only: bc_mode_time, bc_mode_rate
  use step_schedule_mod, only: step_schedule_t, step_schedule_build, step_schedule_plan, &
    & step_load_complete, step_load_next_dt, step_load_landing, step_load_false_position, &
    & step_dwell_dt, &
    & step_load_targets_apart, landing_short, landing_within, landing_passed
  use assert_mod, only: assert_equal, assert_close, assert_true, test_summary

  implicit none

  type(config_t) :: config
  type(step_schedule_t) :: sched
  integer :: unit, ios, i
  character(len=256) :: line
  character(len=:), allocatable :: log_all
  logical :: warned
  real(rk) :: dt_dwell
  logical :: last_dwell
  character(len=64), parameter :: material(7) = [character(len=64) :: 'number_of_phases 1', &
    & 'phase 1', 'crystal_type bcc', 'model elastic', 'c11 231400', 'c12 134700', 'c44 116400']

  ! --- 1: preset, two strain targets, number_of_incrs ---
  ! Step 1: 0 -> 0.001 at 1e-3/s, 1 s in 5 increments; step 2: 0.001 -> 0.002 at 2e-3/s,
  ! 0.5 s in 7 increments; step 2 not printed.
  call write_file('test_step_schedule_1.cfg', [material, [character(len=64) :: &
    & 'boundary_conditions uniaxial_minimal', 'number_of_steps 2', &
    & 'strain_rate 33 1.0e-3 2.0e-3', 'target_strain 33 0.001 0.002', &
    & 'number_of_incrs 5 7', 'print_results 1 0']])
  call config_read(config, 'test_step_schedule_1.cfg')
  call step_schedule_build(sched, config, 0.0_rk)
  call assert_equal(sched%axis, 3, 's1 axis')
  call assert_close(sched%time_start, [0.0_rk, 1.0_rk], 's1 start times')
  call assert_close(sched%time_end, [1.0_rk, 1.5_rk], 's1 end times')
  call assert_close(sched%strain_start(3, :), [0.0_rk, 0.001_rk], 's1 start strains')
  call assert_close(sched%strain_end(3, :), [0.001_rk, 0.002_rk], 's1 end strains')
  call assert_close(sched%strain_rate(3, :), [1.0e-3_rk, 2.0e-3_rk], 's1 strain rates')
  call assert_close(sched%vel_factor, [1.0_rk, 1.0_rk], 's1 velocity factors')
  call assert_equal(sched%num_incrs(1), 5, 's1 increments step 1')
  call assert_equal(sched%num_incrs(2), 7, 's1 increments step 2')
  call assert_close(sched%dt, [0.2_rk, 0.5_rk/7.0_rk], 's1 dt')
  call assert_close(sched%dt_last, sched%dt, 's1 dt_last = dt')
  call assert_equal(sched%total_incrs, 12, 's1 total increments')
  call assert_true(sched%print_results(1) .and. .not. sched%print_results(2), 's1 print flags')

  ! --- 2: preset load reversal with dtime, then a hold to a time target ---
  ! Step 1: 0 -> 0.002 at 1e-3/s, 2 s in increments of 0.5 s (4); step 2: 0.002 -> 0.001 at
  ! -1e-3/s, 1 s (2); step 3: rate 0 until t = 4 s, 1 s in 2 increments.
  call write_file('test_step_schedule_2.cfg', [material, [character(len=64) :: &
    & 'boundary_conditions uniaxial_grip', 'number_of_steps 3', &
    & 'strain_rate 11 1.0e-3 -1.0e-3 0.0', 'target_strain 11 0.002 0.001 -', &
    & 'target_time - - 4.0', 'dtime 0.5 0.5 -', 'number_of_incrs - - 2']])
  call config_read(config, 'test_step_schedule_2.cfg')
  call step_schedule_build(sched, config, 0.0_rk)
  call assert_equal(sched%axis, 1, 's2 axis')
  call assert_close(sched%time_end, [2.0_rk, 3.0_rk, 4.0_rk], 's2 end times')
  call assert_close(sched%strain_end(1, :), [0.002_rk, 0.001_rk, 0.001_rk], 's2 end strains')
  call assert_equal(sched%num_incrs(1), 4, 's2 increments step 1')
  call assert_equal(sched%num_incrs(2), 2, 's2 increments step 2')
  call assert_equal(sched%num_incrs(3), 2, 's2 increments step 3')
  call assert_close(sched%dt, [0.5_rk, 0.5_rk, 0.5_rk], 's2 dt')
  call assert_close(sched%dt_last, sched%dt, 's2 whole steps')

  ! --- 3: dtime that does not divide the step; one within round-off of dividing it ---
  ! Step 1: 0.025 s with dtime 0.01: 0.01, 0.01, 0.005 (and a warning); step 2: 0.03 s with
  ! dtime 0.01 (0.03/0.01 is not exactly 3 in floating point): 3 x 0.01, no remainder; step 3:
  ! 0.004 s with dtime 0.01: one increment of 0.004 s. The log goes to a file to check the
  ! warning. Times that come from differences of targets carry round-off (relative 1e-12).
  call write_file('test_step_schedule_3.cfg', [material, [character(len=64) :: &
    & 'boundary_conditions uniaxial_symmetry', 'number_of_steps 3', &
    & 'strain_rate 22 1.0e-3 1.0e-3 1.0e-3', 'target_time 0.025 0.055 0.059', &
    & 'dtime 0.01 0.01 0.01']])
  call config_read(config, 'test_step_schedule_3.cfg')
  open (newunit=unit, file='test_step_schedule_3.log', status='replace', action='readwrite')
  call log_set_unit(unit)
  call step_schedule_build(sched, config, 0.0_rk)
  call assert_equal(sched%num_incrs(1), 3, 's3 uneven step: increments')
  call assert_close(sched%dt(1), 0.01_rk, 's3 uneven step: full increments')
  call assert_close(sched%dt_last(1), 0.005_rk, 's3 uneven step: last increment', &
    & rel_tol=1.0e-12_rk)
  call assert_close(sched%time_start(1) + 2.0_rk*sched%dt(1) + sched%dt_last(1), 0.025_rk, &
    & 's3 uneven step ends on its target', rel_tol=1.0e-15_rk)
  call assert_equal(sched%num_incrs(2), 3, 's3 round-off step: increments')
  call assert_close(sched%dt(2), 0.01_rk, 's3 round-off step: equal increments')
  call assert_close(sched%dt_last(2), sched%dt(2), 's3 round-off step: no remainder', &
    & rel_tol=0.0_rk, abs_tol=0.0_rk)
  call assert_equal(sched%num_incrs(3), 1, 's3 short step: one increment')
  call assert_close(sched%dt_last(3), 0.004_rk, 's3 short step: its length', rel_tol=1.0e-12_rk)
  call assert_close(sched%strain_end(2, 3), 0.059e-3_rk, 's3 strain at the time target')

  ! The warnings: the uneven step (2 x 0.01 + 0.005) and the short step, not the round-off one.
  ! (The log's messages are wrapped at 72 columns, step 9.7: the lines' texts are joined.)
  rewind (unit)
  log_all = ''
  ios = 0
  do while (ios == 0)
    read (unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    log_all = log_all // ' ' // trim(adjustl(line(10:)))
  end do
  warned = index(log_all, 'step 1, component 22') > 0 .and. &
    & index(log_all, '3 increments: 2 x') > 0
  call assert_true(index(log_all, 'step 2,') == 0, 's3 no warning for the round-off step')
  call assert_true(warned, 's3 warning for the uneven step')
  call log_set_unit(output_unit)
  close (unit)

  ! --- 4: dstrain (dt = dstrain / |rate|) and max_strain met exactly ---
  ! Step 1: 0 -> -0.01 at -2e-3/s, 5 s, dstrain 0.001 -> dt 0.5 s, 10 increments.
  call write_file('test_step_schedule_4.cfg', [material, [character(len=64) :: &
    & 'boundary_conditions uniaxial_minimal', 'number_of_steps 1', &
    & 'strain_rate 33 -2.0e-3', 'target_strain 33 -0.01', 'dstrain 0.001', &
    & 'max_strain 0.01', 'max_total_time 5.0', 'max_incr 10']])
  call config_read(config, 'test_step_schedule_4.cfg')
  call step_schedule_build(sched, config, 0.0_rk)
  call assert_equal(sched%num_incrs(1), 10, 's4 dstrain increments')
  call assert_close(sched%dt(1), 0.5_rk, 's4 dstrain dt')
  call assert_close(sched%time_end(1), 5.0_rk, 's4 ends at max_total_time (allowed)')

  ! --- 5: general conditions: direction from the strain targets, time target as given ---
  ! Rate 1e-3/s (max |v| / L, given). Step 1: 0 -> 0.002, 2 s, factor +1; step 2:
  ! 0.002 -> -0.001, 3 s, factor -1; step 3: time target 6 s, factor +1, strain -0.001 + 1e-3.
  call write_file('test_step_schedule_5.cfg', [material, [character(len=64) :: &
    & 'boundary_conditions general', 'set_bc vel z0 z 0.0', 'set_bc vel z1 z 1.0e-3', &
    & 'number_of_steps 3', 'target_strain 33 0.002 -0.001 -', 'target_time - - 6.0', &
    & 'number_of_incrs 2 3 1']])
  call config_read(config, 'test_step_schedule_5.cfg')
  call step_schedule_build(sched, config, 1.0e-3_rk)
  call assert_equal(sched%axis, 3, 's5 axis from target_strain')
  call assert_close(sched%vel_factor, [1.0_rk, -1.0_rk, 1.0_rk], 's5 velocity factors')
  call assert_close(sched%strain_rate(3, :), [1.0e-3_rk, -1.0e-3_rk, 1.0e-3_rk], &
    & 's5 strain rates')
  call assert_close(sched%time_end, [2.0_rk, 5.0_rk, 6.0_rk], 's5 end times')
  call assert_close(sched%strain_end(3, :), [0.002_rk, -0.001_rk, 0.0_rk], 's5 end strains', &
    & abs_tol=1.0e-15_rk)

  ! --- 6: general conditions with time targets only (no loading component) ---
  call write_file('test_step_schedule_6.cfg', [material, [character(len=64) :: &
    & 'boundary_conditions general', 'set_bc vel z0 z 0.0', 'set_bc vel z1 z 1.0e-3', &
    & 'number_of_steps 2', 'target_time 1.0 3.0', 'dtime 0.25 0.5']])
  call config_read(config, 'test_step_schedule_6.cfg')
  call step_schedule_build(sched, config, 0.0_rk)
  call assert_equal(sched%axis, 0, 's6 no axis')
  call assert_close(reshape(sched%strain_rate, [6]), [(0.0_rk, i = 1, 6)], 's6 no strain rate')
  call assert_equal(sched%num_incrs(1), 4, 's6 increments step 1')
  call assert_equal(sched%num_incrs(2), 4, 's6 increments step 2')

  ! --- 7: triaxial strain control (D45): control 33, component 11 by rates, 22 by targets ---
  ! Step 1: 33 from 0 to 0.002 at 1e-3/s (2 s, 10 increments); 11 at -0.3e-3/s ends at -6e-4;
  ! 22 to -4e-4 at -2e-4/s. Step 2: 33 back to 0.001 at -1e-3/s (1 s, dstrain 2e-4 -> dt 0.2 s,
  ! 5 increments); 11 at 0.3e-3/s ends at -3e-4; 22 to -2e-4 at 2e-4/s.
  call write_file('test_step_schedule_7.cfg', [material, [character(len=64) :: &
    & 'boundary_conditions triaxial', 'control_component 33', 'number_of_steps 2', &
    & 'strain_rate 33 1.0e-3 -1.0e-3', 'target_strain 33 0.002 0.001', &
    & 'strain_rate 11 -0.3e-3 0.3e-3', 'target_strain 22 -0.0004 -0.0002', &
    & 'number_of_incrs 10 -', 'dstrain - 2e-4']])
  call config_read(config, 'test_step_schedule_7.cfg')
  call step_schedule_build(sched, config, 0.0_rk)
  call assert_equal(sched%axis, 3, 's7 control component')
  call assert_true(all(sched%prescribed), 's7 all three components prescribed')
  call assert_close(sched%time_end, [2.0_rk, 3.0_rk], 's7 end times')
  call assert_close(sched%strain_end(1, :), [-6.0e-4_rk, -3.0e-4_rk], 's7 11 from its rates', &
    & rel_tol=1.0e-12_rk)
  call assert_close(sched%strain_rate(2, :), [-2.0e-4_rk, 2.0e-4_rk], 's7 22 rates from targets', &
    & rel_tol=1.0e-12_rk)
  call assert_close(sched%strain_end(2, :), [-4.0e-4_rk, -2.0e-4_rk], 's7 22 end strains')
  call assert_close(sched%strain_start(:, 2), sched%strain_end(:, 1), 's7 step 2 starts at step 1')
  call assert_equal(sched%num_incrs(2), 5, 's7 dstrain with the control rate')
  call assert_close(sched%dt(2), 0.2_rk, 's7 dt = dstrain / |control rate|', rel_tol=1.0e-12_rk)

  ! --- 8: triaxial with a free component (22 has no rows) ---
  call write_file('test_step_schedule_8.cfg', [material, [character(len=64) :: &
    & 'boundary_conditions triaxial', 'control_component 11', 'number_of_steps 1', &
    & 'strain_rate 11 1.0e-3', 'target_strain 11 0.001', 'strain_rate 33 0.0', &
    & 'number_of_incrs 2']])
  call config_read(config, 'test_step_schedule_8.cfg')
  call step_schedule_build(sched, config, 0.0_rk)
  call assert_true(sched%prescribed(1) .and. .not. sched%prescribed(2) .and. &
    & sched%prescribed(3), 's8 22 free')
  call assert_close(sched%strain_end(:, 1), [0.001_rk, 0.0_rk, 0.0_rk], 's8 end strains')

  ! --- 9: load target between strain targets (step 8.1, addendum §3.5) ---
  ! Step 1: 0 -> 0.001 at 1e-3/s (1 s, 2 increments); step 2: target_load 100 N at 1e-3/s,
  ! increments 0.01..0.1 s, planned before the run without an end; step 3 (strain back to
  ! 0.0005) starts after it: not planned before the run; planned from t = 2.0 s and strain
  ! 0.0015 it lasts (0.0015 - 0.0005)/1e-3 = 1 s.
  call write_file('test_step_schedule_9.cfg', [material, [character(len=64) :: &
    & 'boundary_conditions uniaxial_minimal', 'number_of_steps 3', &
    & 'strain_rate 33 1.0e-3 1.0e-3 -1.0e-3', 'target_strain 33 0.001 - 0.0005', &
    & 'target_load 33 - 100.0 -', 'number_of_incrs 2 - 4', 'dtime - 0.1 -', &
    & 'dtime_min - 0.01 -', 'dtime_factor 1.01']])
  call config_read(config, 'test_step_schedule_9.cfg')
  call step_schedule_build(sched, config, 0.0_rk)
  call assert_equal(sched%num_planned, 2, 's9 planned before the run: up to the load target')
  call assert_true(sched%is_load(2) .and. .not. sched%is_load(1) .and. .not. sched%is_load(3), &
    & 's9 load step')
  call assert_close(sched%time_start(2), 1.0_rk, 's9 load step start')
  call assert_close(sched%dt(2), 0.1_rk, 's9 largest increment')
  call assert_close(sched%dt_min(2), 0.01_rk, 's9 smallest increment')
  call assert_close(sched%target_value(2), 100.0_rk, 's9 target')
  call assert_close(sched%dtime_factor, 1.01_rk, 's9 dtime_factor')
  call assert_close(sched%target_tol(2), 0.1_rk, 's9 load_tol default 0.1 N (Q-T18)')
  call assert_close(sched%load_sign(2), 1.0_rk, 's9 direction from the previous strain step: +1')
  call assert_equal(sched%total_incrs, 2, 's9 increments known before the run')
  call step_schedule_plan(sched, config, 3, 2.0_rk, [0.0_rk, 0.0_rk, 0.0015_rk], 0.0_rk)
  call assert_equal(sched%num_planned, 3, 's9 step 3 planned during the run')
  call assert_close(sched%time_end(3), 3.0_rk, 's9 step 3 end from the reached strain', &
    & rel_tol=1.0e-12_rk)
  call assert_equal(sched%num_incrs(3), 4, 's9 step 3 increments')

  ! --- 10: completion test and next increment (v2 loading_load_isinrange, loading_timeincr) ---
  ! Hand values: v = 80 after 60 over 0.1 s (rate 200/s), target 100.
  call assert_true(step_load_complete(100.5_rk, 100.0_rk, 1.0_rk, 0.0_rk), 's10 passed')
  call assert_true(.not. step_load_complete(99.5_rk, 100.0_rk, 1.0_rk, 0.0_rk), 's10 short')
  call assert_true(step_load_complete(99.5_rk, 100.0_rk, 1.0_rk, 1.0_rk), 's10 within tol')
  call assert_true(.not. step_load_complete(100.0_rk, 100.0_rk, 1.0_rk, 0.0_rk), &
    & 's10 exactly at the target with tol 0: not complete (v2 strict test)')
  call assert_true(step_load_complete(-100.5_rk, -100.0_rk, -1.0_rk, 0.0_rk), 's10 downwards')
  ! dt_max 0.05: extrapolation 90 does not reach 100: dt_max.
  call assert_close(step_load_next_dt(80.0_rk, 60.0_rk, 0.1_rk, 100.0_rk, 1.0_rk, 0.0_rk, &
    & 0.05_rk, 1.0e-3_rk, 1.001_rk), 0.05_rk, 's10 not reaching: dt_max')
  ! dt_max 0.2: extrapolation 120 reaches: 1.001 x 0.2 x 20/40 = 0.1001 s.
  call assert_close(step_load_next_dt(80.0_rk, 60.0_rk, 0.1_rk, 100.0_rk, 1.0_rk, 0.0_rk, &
    & 0.2_rk, 1.0e-3_rk, 1.001_rk), 0.1001_rk, 's10 aimed', rel_tol=1.0e-12_rk)
  ! dt_max 0.1: extrapolation exactly 100 (not complete with tol 0, v2): dt_max.
  call assert_close(step_load_next_dt(80.0_rk, 60.0_rk, 0.1_rk, 100.0_rk, 1.0_rk, 0.0_rk, &
    & 0.1_rk, 1.0e-3_rk, 1.001_rk), 0.1_rk, 's10 extrapolation exactly at the target')
  ! Downwards, mirrored: -80 after -60, target -100, dt_max 0.2: 0.1001 s.
  call assert_close(step_load_next_dt(-80.0_rk, -60.0_rk, 0.1_rk, -100.0_rk, -1.0_rk, 0.0_rk, &
    & 0.2_rk, 1.0e-3_rk, 1.001_rk), 0.1001_rk, 's10 aimed downwards', rel_tol=1.0e-12_rk)
  ! Close to the target: 1.001 x 0.1 x 0.01/0.99 = 1.0111e-3 s, below dt_min 0.01: dt_min.
  call assert_close(step_load_next_dt(99.99_rk, 99.0_rk, 0.1_rk, 100.0_rk, 1.0_rk, 0.0_rk, &
    & 0.1_rk, 0.01_rk, 1.001_rk), 0.01_rk, 's10 clamped to dt_min')
  ! Within the tolerance of the extrapolation: tol 25 makes 90 complete at dt_max 0.05:
  ! 1.001 x 0.05 x 20/10 = 0.1001, clamped to dt_max.
  call assert_close(step_load_next_dt(80.0_rk, 60.0_rk, 0.1_rk, 100.0_rk, 1.0_rk, 25.0_rk, &
    & 0.05_rk, 1.0e-3_rk, 1.001_rk), 0.05_rk, 's10 tolerance, clamped to dt_max')

  ! --- 11: landing rule, false position, separation of targets (step 8.2a) ---
  ! Target 100, tolerance 0.1: short, within (both sides and the bounds), past; downwards.
  call assert_equal(step_load_landing(99.8_rk, 100.0_rk, 1.0_rk, 0.1_rk), landing_short, &
    & 's11 short by 0.2')
  call assert_equal(step_load_landing(99.95_rk, 100.0_rk, 1.0_rk, 0.1_rk), landing_within, &
    & 's11 short by 0.05: within')
  call assert_equal(step_load_landing(100.05_rk, 100.0_rk, 1.0_rk, 0.1_rk), landing_within, &
    & 's11 past by 0.05: within')
  call assert_equal(step_load_landing(100.25_rk, 100.0_rk, 1.0_rk, 0.25_rk), landing_within, &
    & 's11 past by exactly the tolerance: within')
  call assert_equal(step_load_landing(99.75_rk, 100.0_rk, 1.0_rk, 0.25_rk), landing_within, &
    & 's11 short by exactly the tolerance: within')
  call assert_equal(step_load_landing(100.2_rk, 100.0_rk, 1.0_rk, 0.1_rk), landing_passed, &
    & 's11 past by 0.2')
  call assert_equal(step_load_landing(-100.2_rk, -100.0_rk, -1.0_rk, 0.1_rk), landing_passed, &
    & 's11 downwards, past')
  call assert_equal(step_load_landing(-99.8_rk, -100.0_rk, -1.0_rk, 0.1_rk), landing_short, &
    & 's11 downwards, short')
  ! False position: bracket (0 s, 80) and (0.2 s, 120), target 100: 0.1 s; (0.1 s, 98) and
  ! (0.2 s, 108), target 100: 0.1 + 0.1 x 2/10 = 0.12 s; downwards (0, -80), (0.2, -130),
  ! target -100: 0.2 x 20/50 = 0.08 s; unloading (0, 50), (0.05, -10), target 0: 0.05 x 50/60.
  call assert_close(step_load_false_position(0.0_rk, 80.0_rk, 0.2_rk, 120.0_rk, 100.0_rk), &
    & 0.1_rk, 's11 false position, linear', rel_tol=1.0e-14_rk)
  call assert_close(step_load_false_position(0.1_rk, 98.0_rk, 0.2_rk, 108.0_rk, 100.0_rk), &
    & 0.12_rk, 's11 false position, inner bracket', rel_tol=1.0e-14_rk)
  call assert_close(step_load_false_position(0.0_rk, -80.0_rk, 0.2_rk, -130.0_rk, -100.0_rk), &
    & 0.08_rk, 's11 false position, downwards', rel_tol=1.0e-14_rk)
  call assert_close(step_load_false_position(0.0_rk, 50.0_rk, 0.05_rk, -10.0_rk, 0.0_rk), &
    & 0.05_rk*50.0_rk/60.0_rk, 's11 false position, unloading to 0', rel_tol=1.0e-14_rk)
  ! Separation (Q-T16): step 1 more than one tolerance from 0, consecutive targets more than
  ! two tolerances apart.
  call assert_true(step_load_targets_apart(0.0_rk, 0.15_rk, 0.1_rk, 1), 's11 first 0.15 > 0.1')
  call assert_true(.not. step_load_targets_apart(0.0_rk, -0.1_rk, 0.1_rk, 1), &
    & 's11 first -0.1: not more than 0.1')
  call assert_true(step_load_targets_apart(100.0_rk, 100.25_rk, 0.1_rk, 2), &
    & 's11 consecutive 0.25 > 0.2')
  call assert_true(.not. step_load_targets_apart(100.0_rk, 99.85_rk, 0.1_rk, 2), &
    & 's11 consecutive 0.15: not more than 0.2')

  ! --- 12: load-controlled components other than the control component (step 8.2b) ---
  ! Control 33: step 1 to strain 0.002 at 1e-3/s (2 s, 4 increments) with 11 at -50 MPa in time
  ! (target_stress) and 22 at a strain rate; step 2 to 300 N with 11 on the path to 0 MPa and 22
  ! to 10 N (target_load). Only step 1 is planned before the run (11's strain is measured); its
  ! 11 is load-controlled (no planned strain), 22 strain-controlled.
  call write_file('test_step_schedule_12.cfg', [material, [character(len=64) :: &
    & 'boundary_conditions triaxial', 'control_component 33', 'number_of_steps 2', &
    & 'strain_rate 33 1.0e-3 1.0e-3', 'target_strain 33 0.002 -', &
    & 'target_load 33 - 300.0', 'target_stress 11 -50.0 0.0', 'strain_rate 22 -2.0e-4 -', &
    & 'target_load 22 - 10.0', 'number_of_incrs 4 -', 'dtime - 0.1', 'dtime_min - 0.01', &
    & 'load_tol 0.5', 'stress_tol 0.02', 'min_pert_frac 0.01']])
  call config_read(config, 'test_step_schedule_12.cfg')
  call assert_close(config%min_pert_frac, 0.01_rk, 's12 min_pert_frac')
  call step_schedule_build(sched, config, 0.0_rk)
  call assert_equal(sched%num_planned, 1, 's12 planned before the run: step 1 only')
  call assert_true(all(sched%bc_kind(:, 1) == [target_by_stress, 0, 0]) .and. &
    & all(sched%bc_kind(:, 2) == [target_by_stress, target_by_load, 0]), 's12 kinds')
  call assert_close(sched%bc_target(:, 1), [-50.0_rk, 0.0_rk, 0.0_rk], 's12 targets step 1')
  call assert_close(sched%bc_target(:, 2), [0.0_rk, 10.0_rk, 0.0_rk], 's12 targets step 2')
  call assert_close(sched%bc_tol(:, 2), [0.02_rk, 0.5_rk, 0.0_rk], 's12 tolerances step 2')
  call assert_close(sched%strain_rate(:, 1), [0.0_rk, -2.0e-4_rk, 1.0e-3_rk], &
    & 's12 strain rates: none planned for 11')
  call assert_close(sched%strain_end(:, 1), [0.0_rk, -4.0e-4_rk, 0.002_rk], &
    & 's12 end strains: 11 not planned (its start)')
  call assert_equal(sched%num_incrs(1), 4, 's12 increments step 1')

  ! --- 13: load rates of the control component (step 8.2c, iteration (B)) ---
  ! Control 33 at 2 N/s to 100 N (50 s: dtime 15 gives 3 x 15 s + 5 s), 11 to -37.5 N in time, 22
  ! by its load rate 0.5 N/s; a 20 s dwell (dtime 5); at -2 N/s back to 0 N in 4 increments. Only
  ! step 1 is planned before the run (the others from the values measured at their starts): the
  ! dwell from 99.97 N at 50 s, step 3 from 100.02 N at 70 s (50.01 s).
  call write_file('test_step_schedule_13.cfg', [material, [character(len=64) :: &
    & 'boundary_conditions triaxial', 'control_component 33', 'number_of_steps 3', &
    & 'load_rate 33 2.0 0.0 -2.0', 'target_load 33 100.0 - 0.0', &
    & 'target_duration - 20.0 -', 'target_load 11 -37.5 -37.5 0.0', &
    & 'load_rate 22 0.5 0.0 -0.5', 'dtime 15.0 5.0 -', 'number_of_incrs - - 4', &
    & 'load_tol 0.5', 'max_strain_incr 0.002']])
  call config_read(config, 'test_step_schedule_13.cfg')
  call step_schedule_build(sched, config, 0.0_rk)
  call assert_equal(sched%num_planned, 1, 's13 planned before the run: step 1 only')
  call assert_true(all(sched%is_rate) .and. all(sched%is_dwell .eqv. [.false., .true., &
    & .false.]) .and. .not. any(sched%is_load .or. sched%is_open), 's13 step kinds')
  call assert_true(all(sched%bc_kind(:, 1) == target_by_load) .and. &
    & all(sched%bc_mode(:, 1) == [bc_mode_time, bc_mode_rate, bc_mode_time]), &
    & 's13 step 1: all three load-controlled, 22 by its rate, 33 and 11 in time')
  call assert_close(sched%bc_rate(:, 1), [0.0_rk, 0.5_rk, 2.0_rk], 's13 rates step 1')
  call assert_close(sched%bc_target(:, 1), [-37.5_rk, 0.0_rk, 100.0_rk], 's13 targets step 1')
  call assert_close(sched%time_end(1), 50.0_rk, 's13 step 1 duration from 0 N', rel_tol=1.0e-15_rk)
  call assert_equal(sched%num_incrs(1), 4, 's13 step 1 increments (D44)')
  call assert_close([sched%dt(1), sched%dt_last(1)], [15.0_rk, 5.0_rk], 's13 step 1 dt', &
    & rel_tol=1.0e-14_rk)
  call assert_true(all(sched%bc_mode(:, 2) == bc_mode_rate .or. sched%bc_mode(:, 2) == &
    & bc_mode_time) .and. sched%bc_mode(3, 2) == bc_mode_rate, 's13 dwell: the control held')
  call step_schedule_plan(sched, config, 2, 50.0_rk, [0.0_rk, 0.0_rk, 0.0_rk], 0.0_rk, 99.97_rk)
  call assert_close(sched%time_end(2), 70.0_rk, 's13 dwell end', rel_tol=1.0e-15_rk)
  call assert_equal(sched%num_incrs(2), 0, 's13 dwell increments chosen during the run')
  call assert_close(sched%dt(2), 5.0_rk, 's13 dwell dtime')
  call step_schedule_plan(sched, config, 3, 70.0_rk, [0.0_rk, 0.0_rk, 0.0_rk], 0.0_rk, 100.02_rk)
  call assert_close(sched%time_end(3), 120.01_rk, 's13 step 3 from the measured value', &
    & rel_tol=1.0e-14_rk)
  call assert_close(sched%dt(3), 12.5025_rk, 's13 step 3 equal increments', rel_tol=1.0e-14_rk)
  call assert_close(sched%load_sign(3), -1.0_rk, 's13 step 3 downwards')

  ! --- Dwell increments (step 8.2c, v2 driver_triaxclr) ---
  ! Episode start: 0.01 dtime; then the strain bound 0.002 x 1 mm / 1e-3 mm/s = 2 s below dtime
  ! 5 s; V = 0: dtime; at most the remaining time, all of it when less than 0.001 dtime would
  ! remain (5 s of 5.004 s left: 5.004).
  call step_dwell_dt(.true., 1.0e-3_rk, 1.0_rk, 0.002_rk, 5.0_rk, 20.0_rk, dt_dwell, last_dwell)
  call assert_close(dt_dwell, 0.05_rk, 'dwell start: 0.01 dtime', rel_tol=1.0e-15_rk)
  call assert_true(.not. last_dwell, 'dwell start: not the last')
  call step_dwell_dt(.false., -1.0e-3_rk, 1.0_rk, 0.002_rk, 5.0_rk, 19.95_rk, dt_dwell, &
    & last_dwell)
  call assert_close(dt_dwell, 2.0_rk, 'dwell: strain bound', rel_tol=1.0e-15_rk)
  call step_dwell_dt(.false., 0.0_rk, 1.0_rk, 0.002_rk, 5.0_rk, 19.95_rk, dt_dwell, last_dwell)
  call assert_close(dt_dwell, 5.0_rk, 'dwell: zero velocity, dtime', rel_tol=1.0e-15_rk)
  call step_dwell_dt(.false., 1.0e-12_rk, 1.0_rk, 0.002_rk, 5.0_rk, 4.95_rk, dt_dwell, last_dwell)
  call assert_close(dt_dwell, 4.95_rk, 'dwell: the remaining time', abs_tol=0.0_rk)
  call assert_true(last_dwell, 'dwell: the remaining time ends the step')
  call step_dwell_dt(.false., 1.0e-12_rk, 1.0_rk, 0.002_rk, 5.0_rk, 5.004_rk, dt_dwell, &
    & last_dwell)
  call assert_close(dt_dwell, 5.004_rk, 'dwell: a small remainder taken in', abs_tol=0.0_rk)
  call assert_true(last_dwell, 'dwell: the small remainder ends the step')
  call step_dwell_dt(.false., 1.0e-12_rk, 1.0_rk, 0.002_rk, 5.0_rk, 5.006_rk, dt_dwell, &
    & last_dwell)
  call assert_close(dt_dwell, 5.0_rk, 'dwell: a larger remainder left', abs_tol=0.0_rk)
  call assert_true(.not. last_dwell, 'dwell: the larger remainder does not end the step')

  call test_summary()

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

end program test_step_schedule
