! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The loading steps as a time schedule, planned step by step: for each step its start and end
!> times, the strains of the controlled normal components at both ends, their strain rates,
!> and its increments (count and sizes), with the checks that need the start-of-step values;
!> for load and stress targets (step 8.1) the step completion test and the next increment.
!>
!> Source : port of v2's step setup for strain and time targets
!>          (boundary_conditions_mod2.f90::calc_bcs_steps, case uniaxial_strain_target, and
!>          driver_uniaxial_control_mod1.f90), v2-fixed C29-C32, with the semantics decided in
!>          ROADMAP steps 5.1 (D44) and 5.3 (D45: control component, triaxial strain control);
!>          load and stress targets: port of v2's uniaxial load target
!>          (general/loading_type_mod1.f90::loading_timeincr, loading_type_mod2.f90::
!>          loading_load_isinrange, calc_bcs_steps case uniaxial_load_target; v2-fixed C52),
!>          dev/theory/addendum_load_targets.md §2-§3 (step 8.1); landing on a load or stress
!>          target (step 8.2a, dev/theory/addendum_triaxial_load_control.md §3.2, §3.6,
!>          decisions Q-T13, Q-T14, Q-T16): the false position step on the increment length,
!>          Press, Teukolsky, Vetterling & Flannery (2007), Numerical Recipes, 3rd ed., §9.2
!>          (step_load_false_position). Otherwise time-step bookkeeping and v2's linear
!>          extrapolation of the load.
!>
!> Semantics (v2 unless noted):
!>   - A step ends at the target of the control component (D45). A strain target
!>     (engineering strain) is reached at constant strain rate, so the step lasts
!>     |de|/|rate|, de the change of the target from the start of the step (v2: |de| L / v
!>     with v = rate L).
!>   - The other components with rows (triaxial, D45) follow the step's duration T: a given
!>     strain_rate r gives the end strain e0 + r T; a given target_strain e gives the rate
!>     (e - e0)/T. Components without rows are free (strain 0 here: not prescribed).
!>   - Presets: the rate is the step's strain_rate; its sign must move towards a strain target
!>     (an error otherwise; v2 silently used the magnitude, dev/cfg_keywords.md §1).
!>   - General: the strain rate of a strain-target step is max |set_bc velocity| / gage length
!>     (v2: step_velocity / gage_length), signed towards the target: the set_bc velocities are
!>     multiplied by vel_factor = +1 or -1 (v2's target_sign / vel_factor). Time-target steps
!>     use the set_bc velocities as given (vel_factor +1, as v2).
!>   - A time target is cumulative and must be later than the start of the step.
!>   - Increments: number_of_incrs n gives n equal increments. dtime (or dstrain, as
!>     dtime = dstrain / |rate| with the control component's rate) gives full increments of
!>     that size and a shorter last one so that the step ends exactly on its target, with a
!>     warning (D44; v2 overshot the target with a full last increment). A step within
!>     round-off of a whole number of increments (relative 1e-9) is cut into that number of
!>     equal increments, without a warning.
!>   - Load and stress targets (step 8.1, addendum §2-§3): the step's duration is not known.
!>     Its largest increment is dtime (or dstrain / |rate|), its smallest dtime_min; the step
!>     ends after the first increment whose value v (force on the face x1/y1/z1 of the control
!>     component [N], or that force over the face's current area [MPa]) is within the
!>     tolerance tol (load_tol [N] or stress_tol [MPa], > 0) of the target v*
!>     (step_load_landing; s = sign(v* - v_start)). An increment that passes the target by more
!>     than tol is discarded by the caller and redone with the length given by false position
!>     on the bracket of the increment lengths and values reached (step_load_false_position;
!>     landing rule, step 8.2a, addendum_triaxial_load_control.md §3.6). The next increment is
!>     dt_max, or, when the linear extrapolation of the last increment's rate would reach the
!>     target within dt_max, clamp(f (v* - v)/vdot, dt_min, dt_max), f = dtime_factor
!>     (step_load_next_dt). Presets: the strain rate's sign must move v towards its target (a
!>     positive rate is taken to increase the load; user decision Q-L6); general: the set_bc
!>     velocities are reversed (vel_factor -1) when the target is below the start value (v2).
!>     The direction is checked at the start of the step from the actual value
!>     (step_schedule_load_start), and before the run for step 1 (from 0) and for a load or
!>     stress target that follows a target of the same kind (from that target). Before the run
!>     (Q-T16, step_load_targets_apart): a load or stress target of step 1 must differ from 0
!>     by more than tol, and consecutive targets of the same kind by more than 2 tol (each
!>     step lands within tol of its target, so the next one starts more than tol from its
!>     own).
!>   - Planning (addendum §3.5, §3.7): steps are planned from their actual start (time and
!>     strains). Before the run (step_schedule_build) the steps up to and including the first
!>     load- or stress-target step, or the first step with load-controlled components (8.2b),
!>     are planned; every later step is planned at its start
!>     during the run (step_schedule_plan), from the time and strains reached.
!>   - Load-controlled components (step 8.2b, dev/theory/addendum_triaxial_load_control.md
!>     §3.1): a component other than the control component with a target_load or
!>     target_stress row in a step has its face velocity found by the boundary-velocity
!>     iteration (bc_iteration_mod); the schedule records its kind, target and tolerance
!>     (bc_kind, bc_target, bc_tol); its strain is not planned (rate 0 here, end strain
!>     measured during the run).
!>   - Load and stress rates of the control component (step 8.2c,
!>     dev/theory/addendum_triaxial_load_control.md §3.1-§3.4; iteration (B)): the control
!>     component is load-controlled (bc_kind, its value to reach in time); a step to a load or
!>     stress target lasts T = |v* - v_start|/|rate| from the value v_start measured at its
!>     start (Q-T4; planned at its start, with the rate's direction towards the target and the
!>     target not already reached checked there, and before the run for step 1 and between
!>     consecutive targets of the same kind), its increments by D44's rule (Q-T5); a step to a
!>     time or duration target (target_duration: its length from its actual start, Q-T7) lasts
!>     that time; a rate 0 is a dwell, whose increments are chosen during the run
!>     (step_dwell_dt; num_incrs set when it ends); a strain target at a load or stress rate
!>     (Q-T8 (i)) has no known duration: dtime is its largest increment, dtime_min its smallest,
!>     the next increment by v2's extrapolation of the strain (step_load_next_dt with tolerance
!>     0: the step ends at the first increment that reaches the strain), the rate's sign must
!>     move the strain towards its target (a positive rate increases the strain).
!>   - Rows of load-controlled components (step 8.2c): bc_mode gives how each value to reach
!>     is defined: on the step's load path (a target in a step of unknown duration at a strain
!>     rate, v2 CSR), in time (a target, the duration known), or by a rate (bc_rate).
!>   - Run limits (v2-fixed C32): a run whose planned increments or end time exceed max_incr
!>     or max_total_time is an error before it starts (ending exactly at a limit is allowed;
!>     only the steps planned before the run count); with max_strain, a step whose end strain
!>     of a prescribed component is known and exceeds it is an error when it is planned (D44:
!>     the limit applies to every control mode, no default).
!>
!> Restart (step 8.4, dev/theory/addendum_restart.md §3 R1): a restarted run's configuration
!> lists the steps after the restart point; the last step of the run that wrote the restart
!> files is the previous step of its step 1 (step_record_t, given to step_schedule_build), and
!> the steps are planned from that step's end time and strains, with the checks of step 1 made
!> against it as against any previous step.
!>
!> Public types:
!>   step_schedule_t  the schedule (components at the type)
!>   step_record_t    one step of a schedule (the restart point's step)
!> Public routines:
!>   step_schedule_build       allocate the schedule and plan the steps known before the run
!>   step_schedule_record      one step of a schedule as a step_record_t
!>   step_schedule_plan        plan one step from its actual start (during the run)
!>   step_schedule_load_start  direction and checks of a load- or stress-target step at its
!>                             start
!>   step_load_complete        v2's completion test of a load or stress target (used for the
!>                             extrapolated value of step_load_next_dt)
!>   step_load_next_dt         next increment of a load- or stress-target step (v2)
!>   step_load_landing         short of, within the tolerance of, or past a target (8.2a)
!>   step_load_false_position  increment length of the next landing attempt (8.2a)
!>   step_load_targets_apart   targets separated by more than a number of tolerances (8.2a)
!>   step_dwell_dt             next increment of a dwell (8.2c)
module step_schedule_mod

  use kinds_mod, only: rk
  use error_mod, only: fatal_error
  use log_mod, only: log_warning
  use config_mod, only: config_t, incr_by_count, incr_by_dtime, target_by_load, &
    & target_by_stress, target_by_strain, rate_by_strain, rate_by_load, rate_by_stress
  use bc_iteration_mod, only: bc_mode_path, bc_mode_time, bc_mode_rate

  implicit none

  private

  public :: step_schedule_t, step_record_t, step_schedule_build, step_schedule_record
  public :: step_schedule_plan, step_schedule_load_start
  public :: step_load_complete, step_load_next_dt, step_load_landing, step_load_false_position
  public :: step_load_targets_apart, step_dwell_dt
  public :: landing_short, landing_within, landing_passed

  ! Results of step_load_landing: the value is short of the target by more than the tolerance,
  ! within the tolerance of it, or past it by more than the tolerance.
  integer, parameter :: landing_short = -1, landing_within = 0, landing_passed = 1

  ! Relative tolerance under which a step counts as a whole number of increments (round-off
  ! of T/dt, about 1e-16 relative, kept far below any intended remainder).
  real(rk), parameter :: tol_whole = 1.0e-9_rk

  ! Smallest extrapolated load change for which v2 adjusts the increment (loading_timeincr).
  real(rk), parameter :: load_change_min = 1.0e-16_rk

  ! Dwells (v2 driver_triaxclr, step 8.2c): the first increment of a dwell episode as a fraction
  ! of dtime (initial_time_frac), and the remaining time, as a fraction of dtime, below which
  ! the step ends (time_tol; FEPX: taken into the last increment, which ends the step at its
  ! time exactly).
  real(rk), parameter :: dwell_first_frac = 0.01_rk
  real(rk), parameter :: dwell_time_tol = 0.001_rk

  !> One step of a schedule (the components of step_schedule_t at one step; see there), and
  !> the step's per_rate row of the configuration (periodic conditions): the previous step of
  !> a restarted run's step 1 (step 8.4).
  type :: step_record_t
    integer :: target_kind = 0
    logical :: is_load = .false.
    real(rk) :: target_value = 0.0_rk
    real(rk) :: target_tol = 0.0_rk
    real(rk) :: load_sign = 1.0_rk
    real(rk) :: dt_min = 0.0_rk
    real(rk) :: strain_rate(3) = 0.0_rk
    real(rk) :: vel_factor = 1.0_rk
    real(rk) :: time_start = 0.0_rk
    real(rk) :: time_end = 0.0_rk
    real(rk) :: strain_start(3) = 0.0_rk
    real(rk) :: strain_end(3) = 0.0_rk
    integer :: num_incrs = 0
    real(rk) :: dt = 0.0_rk
    real(rk) :: dt_last = 0.0_rk
    logical :: print_results = .false.
    integer :: bc_kind(3) = 0
    real(rk) :: bc_target(3) = 0.0_rk
    real(rk) :: bc_tol(3) = 0.0_rk
    integer :: bc_mode(3) = 0
    real(rk) :: bc_rate(3) = 0.0_rk
    integer :: rate_kind = 0
    logical :: is_rate = .false.
    logical :: is_dwell = .false.
    logical :: is_open = .false.
    real(rk) :: ctrl_rate = 0.0_rk
    real(rk) :: per_rate(9) = 0.0_rk
  end type step_record_t

  !> The schedule.
  !>
  !> Components:
  !>   num_steps     integer                  number of steps
  !>   num_planned   integer                  steps planned so far (1..num_planned)
  !>   axis          integer                  control component (1..3), 0 if none (general
  !>                                          with time targets only)
  !>   comp_label    character(2)(3)          the components as text, for messages (11, 22,
  !>                                          33; periodic conditions: the control component
  !>                                          ij as component 1, step 8.3b)
  !>   prescribed    logical(3)               normal components whose strain is prescribed
  !>                                          (control component and components with rows)
  !>   target_kind   integer(num_steps)       target of each step (config_mod target_by_*)
  !>   is_load       logical(num_steps)       the step ends at a load or stress target
  !>   target_value  real(rk)(num_steps)      load [N] or stress [MPa] target (load steps)
  !>   target_tol    real(rk)(num_steps)      its tolerance [N] or [MPa] (load steps)
  !>   load_sign     real(rk)(num_steps)      s = sign(target - value at the step start), +1
  !>                                          or -1, set by step_schedule_load_start (load
  !>                                          steps)
  !>   dt_min        real(rk)(num_steps)      smallest increment of a load step [s]
  !>   dtime_factor  real(rk)                 factor on the increment aimed at a target [-]
  !>   strain_rate   real(rk)(3, num_steps)   engineering strain rate of each prescribed
  !>                                          normal component during the step (0 otherwise)
  !>                                          [1/s]
  !>   vel_factor    real(rk)(num_steps)      general: factor (+1 or -1) on the set_bc
  !>                                          velocities; presets: 1 [-]
  !>   time_start    real(rk)(num_steps)      time at the start of the step [s]
  !>   time_end      real(rk)(num_steps)      time at the end of the step (the target; load
  !>                                          steps: set by the caller when reached) [s]
  !>   strain_start  real(rk)(3, num_steps)   engineering strain of each prescribed normal
  !>                                          component at the start of the step (0
  !>                                          otherwise) [-]
  !>   strain_end    real(rk)(3, num_steps)   the same at the end of the step (load steps: set
  !>                                          by the caller when reached) [-]
  !>   num_incrs     integer(num_steps)       increments of the step (load steps: 0 when
  !>                                          planned, the count when completed)
  !>   dt            real(rk)(num_steps)      size of the increments but the last (load steps:
  !>                                          the largest increment) [s]
  !>   dt_last       real(rk)(num_steps)      size of the last increment (= dt when the step
  !>                                          is a whole number of increments) [s]
  !>   print_results logical(num_steps)       write the results at the end of the step
  !>   total_incrs   integer                  increments of the steps planned before the run
  !>   bc_kind       integer(3, num_steps)    load-controlled components other than the
  !>                                          control component: target_by_load or
  !>                                          target_by_stress (0: not load-controlled)
  !>                                          (step 8.2b)
  !>   bc_target     real(rk)(3, num_steps)   their load [N] or stress [MPa] targets
  !>   bc_tol        real(rk)(3, num_steps)   their tolerances (load_tol [N], stress_tol [MPa])
  !>   bc_mode       integer(3, num_steps)    how their values to reach are given:
  !>                                          bc_mode_path, bc_mode_time or bc_mode_rate
  !>                                          (step 8.2c; 0: not load-controlled)
  !>   bc_rate       real(rk)(3, num_steps)   their rates (bc_mode_rate) [N/s] or [MPa/s]
  !>   rate_kind     integer(num_steps)       rate of the control component (config_mod
  !>                                          rate_by_*; step 8.2c)
  !>   is_rate       logical(num_steps)       the control component has a load or stress rate
  !>                                          (iteration (B); it is then load-controlled)
  !>   is_dwell      logical(num_steps)       that rate is 0: a dwell
  !>   is_open       logical(num_steps)       a strain target at a load or stress rate (duration
  !>                                          not known); its target in target_value [-],
  !>                                          tolerance 0
  !>   ctrl_rate     real(rk)(num_steps)      the control component's load or stress rate (0
  !>                                          otherwise) [N/s] or [MPa/s]
  !>   max_strain_incr  real(rk)              largest strain increment of a dwell [-]
  !>   has_prev      logical                  step 1 has a previous step (a restarted run,
  !>                                          step 8.4): prev
  !>   prev          type(step_record_t)      that step: the last step of the run that wrote
  !>                                          the restart files
  type :: step_schedule_t
    integer :: num_steps = 0
    integer :: num_planned = 0
    integer :: axis = 0
    character(len=2) :: comp_label(3) = ['11', '22', '33']
    logical :: prescribed(3) = .false.
    integer, allocatable :: target_kind(:)
    logical, allocatable :: is_load(:)
    real(rk), allocatable :: target_value(:)
    real(rk), allocatable :: target_tol(:)
    real(rk), allocatable :: load_sign(:)
    real(rk), allocatable :: dt_min(:)
    real(rk) :: dtime_factor = 1.0_rk
    real(rk), allocatable :: strain_rate(:, :)
    real(rk), allocatable :: vel_factor(:)
    real(rk), allocatable :: time_start(:)
    real(rk), allocatable :: time_end(:)
    real(rk), allocatable :: strain_start(:, :)
    real(rk), allocatable :: strain_end(:, :)
    integer, allocatable :: num_incrs(:)
    real(rk), allocatable :: dt(:)
    real(rk), allocatable :: dt_last(:)
    logical, allocatable :: print_results(:)
    integer :: total_incrs = 0
    integer, allocatable :: bc_kind(:, :)
    real(rk), allocatable :: bc_target(:, :)
    real(rk), allocatable :: bc_tol(:, :)
    integer, allocatable :: bc_mode(:, :)
    real(rk), allocatable :: bc_rate(:, :)
    integer, allocatable :: rate_kind(:)
    logical, allocatable :: is_rate(:)
    logical, allocatable :: is_dwell(:)
    logical, allocatable :: is_open(:)
    real(rk), allocatable :: ctrl_rate(:)
    real(rk) :: max_strain_incr = 0.0_rk
    logical :: has_prev = .false.
    type(step_record_t) :: prev
  end type step_schedule_t

contains

  !> Allocate the schedule of a configuration and plan, from time 0 and strain 0, every step up
  !> to and including the first load- or stress-target step (the later ones start at values
  !> known only during the run), with the checks listed in the module header that can be made
  !> before the run: the direction of step 1 and of a load or stress target following a target
  !> of the same kind, dtime_min not larger than the largest increment of every load step, and
  !> the run limits of the planned steps (fatal errors naming the step and the rows).
  !>
  !> Arguments:
  !>   sched         (out) type(step_schedule_t)  the schedule
  !>   config        (in)  type(config_t)         the configuration
  !>   rate_general  (in)  real(rk)               general conditions: max |set_bc velocity| /
  !>                                              gage length of the loading axis (> 0 when
  !>                                              axis > 0); ignored with a preset [1/s]
  !>   prev          (in, optional) type(step_record_t)  a restarted run (step 8.4): the last
  !>                                              step of the run that wrote the restart files,
  !>                                              the previous step of step 1; the steps are
  !>                                              planned from its end, and a step 1 at a load
  !>                                              or stress rate during the run
  !>   incrs_done    (in, optional) integer       increments done before (a restarted run; for
  !>                                              max_incr; default 0)
  !> Locals:
  !>   i_step        integer                      step
  !>   axis          integer                      control component (0: none)
  !>   comp          character(:), allocatable    the loading component as text
  !>   step          character(:), allocatable    "<file>: step <i>, component <c>"
  !>   t0            real(rk)                     time at the start of the next step [s]
  !>   e0            real(rk)(3)                  strains at the start of the next step [-]
  !>   total         real(rk)                     planned increments, as a real (overflow-safe
  !>                                              sum)
  !>   rate          real(rk)                     control strain rate of a load step [1/s]
  !>   dt_max        real(rk)                     largest increment of a load step [s]
  !>   value_prev    real(rk)                     start value of a load step known before the
  !>                                              run (0 for step 1, the previous target) [N]
  !>                                              or [MPa]
  !>   c             integer                      a normal component
  subroutine step_schedule_build(sched, config, rate_general, prev, incrs_done)

    type(step_schedule_t), intent(out) :: sched
    type(config_t), intent(in) :: config
    real(rk), intent(in) :: rate_general
    type(step_record_t), intent(in), optional :: prev
    integer, intent(in), optional :: incrs_done

    integer :: i_step, axis, c
    character(len=:), allocatable :: comp, step
    real(rk) :: t0, e0(3), total, rate, dt_max, value_prev

    ! A restarted run: the previous step of step 1 (step 8.4).
    if (present(prev)) then
      sched%has_prev = .true.
      sched%prev = prev
    end if

    ! Allocate the per-step arrays.
    sched%num_steps = config%num_steps
    allocate (sched%strain_rate(3, sched%num_steps), sched%vel_factor(sched%num_steps))
    allocate (sched%time_start(sched%num_steps), sched%time_end(sched%num_steps))
    allocate (sched%strain_start(3, sched%num_steps), sched%strain_end(3, sched%num_steps))
    allocate (sched%num_incrs(sched%num_steps), sched%dt(sched%num_steps))
    allocate (sched%dt_last(sched%num_steps))
    allocate (sched%target_kind(sched%num_steps), sched%is_load(sched%num_steps))
    allocate (sched%target_value(sched%num_steps), sched%target_tol(sched%num_steps))
    allocate (sched%load_sign(sched%num_steps), sched%dt_min(sched%num_steps))
    allocate (sched%rate_kind(sched%num_steps), sched%is_rate(sched%num_steps))
    allocate (sched%is_dwell(sched%num_steps), sched%is_open(sched%num_steps))
    allocate (sched%ctrl_rate(sched%num_steps))
    sched%strain_rate = 0.0_rk
    sched%vel_factor = 1.0_rk
    sched%time_start = 0.0_rk
    sched%time_end = 0.0_rk
    sched%strain_start = 0.0_rk
    sched%strain_end = 0.0_rk
    sched%num_incrs = 0
    sched%dt = 0.0_rk
    sched%dt_last = 0.0_rk
    sched%print_results = config%print_results
    axis = config%control_axis
    sched%axis = axis
    sched%comp_label = config%comp_label
    sched%prescribed = config%comp_rows
    if (axis > 0) sched%prescribed(axis) = .true.
    comp = ''
    if (axis > 0) comp = trim(sched%comp_label(axis))

    ! Targets of the steps: kind, and for load and stress targets the value and tolerance.
    sched%target_kind = config%target_kind
    sched%rate_kind = config%rate_kind
    sched%is_rate = config%rate_kind /= rate_by_strain
    sched%is_load = (config%target_kind == target_by_load .or. &
      & config%target_kind == target_by_stress) .and. .not. sched%is_rate
    sched%is_open = sched%is_rate .and. config%target_kind == target_by_strain
    sched%ctrl_rate = 0.0_rk
    if (axis > 0) then
      where (config%rate_kind == rate_by_load) sched%ctrl_rate = config%load_rate(axis, :)
      where (config%rate_kind == rate_by_stress) sched%ctrl_rate = config%stress_rate(axis, :)
    end if
    sched%is_dwell = sched%is_rate .and. .not. abs(sched%ctrl_rate) > 0.0_rk
    sched%max_strain_incr = config%max_strain_incr
    sched%target_value = 0.0_rk
    sched%target_tol = 0.0_rk
    sched%load_sign = 1.0_rk
    sched%dt_min = config%dtime_min
    sched%dtime_factor = config%dtime_factor
    do i_step = 1, sched%num_steps
      if (sched%target_kind(i_step) == target_by_load) then
        sched%target_value(i_step) = config%target_load(axis, i_step)
        sched%target_tol(i_step) = config%load_tol
      else if (sched%target_kind(i_step) == target_by_stress) then
        sched%target_value(i_step) = config%target_stress(axis, i_step)
        sched%target_tol(i_step) = config%stress_tol
      else if (sched%is_open(i_step)) then
        sched%target_value(i_step) = config%target_strain(axis, i_step)
        sched%target_tol(i_step) = 0.0_rk
      end if
    end do

    ! Load-controlled components (step 8.2b): kind, target, tolerance, and how the value to
    ! reach is given (step 8.2c): the other components with load or stress rows (targets on the
    ! load path of a strain-rate step of unknown duration, otherwise in time; rates), and the
    ! control component at a load or stress rate (its target in time, otherwise its rate).
    allocate (sched%bc_kind(3, sched%num_steps), sched%bc_target(3, sched%num_steps))
    allocate (sched%bc_tol(3, sched%num_steps), sched%bc_mode(3, sched%num_steps))
    allocate (sched%bc_rate(3, sched%num_steps))
    sched%bc_kind = 0
    sched%bc_target = 0.0_rk
    sched%bc_tol = 0.0_rk
    sched%bc_mode = 0
    sched%bc_rate = 0.0_rk
    do i_step = 1, sched%num_steps
      do c = 1, 3
        if (c == axis) cycle
        if (config%has_target_load(c, i_step)) then
          sched%bc_kind(c, i_step) = target_by_load
          sched%bc_target(c, i_step) = config%target_load(c, i_step)
          sched%bc_tol(c, i_step) = config%load_tol
        else if (config%has_target_stress(c, i_step)) then
          sched%bc_kind(c, i_step) = target_by_stress
          sched%bc_target(c, i_step) = config%target_stress(c, i_step)
          sched%bc_tol(c, i_step) = config%stress_tol
        else if (config%has_load_rate(c, i_step)) then
          sched%bc_kind(c, i_step) = target_by_load
          sched%bc_rate(c, i_step) = config%load_rate(c, i_step)
          sched%bc_tol(c, i_step) = config%load_tol
          sched%bc_mode(c, i_step) = bc_mode_rate
        else if (config%has_stress_rate(c, i_step)) then
          sched%bc_kind(c, i_step) = target_by_stress
          sched%bc_rate(c, i_step) = config%stress_rate(c, i_step)
          sched%bc_tol(c, i_step) = config%stress_tol
          sched%bc_mode(c, i_step) = bc_mode_rate
        end if
        if (sched%bc_kind(c, i_step) /= 0 .and. sched%bc_mode(c, i_step) == 0) &
          & sched%bc_mode(c, i_step) = merge(bc_mode_path, bc_mode_time, sched%is_load(i_step))
      end do
      if (sched%is_rate(i_step)) then
        if (sched%rate_kind(i_step) == rate_by_load) then
          sched%bc_kind(axis, i_step) = target_by_load
          sched%bc_tol(axis, i_step) = config%load_tol
        else
          sched%bc_kind(axis, i_step) = target_by_stress
          sched%bc_tol(axis, i_step) = config%stress_tol
        end if
        sched%bc_rate(axis, i_step) = sched%ctrl_rate(i_step)
        sched%bc_mode(axis, i_step) = bc_mode_rate
        if (sched%target_kind(i_step) == target_by_load .or. &
          & sched%target_kind(i_step) == target_by_stress) then
          sched%bc_target(axis, i_step) = sched%target_value(i_step)
          sched%bc_mode(axis, i_step) = bc_mode_time
        end if
      end if
    end do

    ! General conditions with a loading component need moving set_bc velocities.
    if (config%bc_type == 'general' .and. axis > 0 .and. .not. rate_general > 0.0_rk) &
      & call fatal_error('cfg', config%file_name // ': boundary_conditions general with ' // &
      & 'a target on component ' // comp // ': every set_bc velocity is zero')

    ! Load and stress steps, before the run: dtime_min within the largest increment; step 1's
    ! target more than one tolerance from 0 and consecutive targets of the same kind more than
    ! two tolerances apart (Q-T16); the direction where the start value is known (step 1 from
    ! 0; a target following a target of the same kind, from that target).
    do i_step = 1, sched%num_steps
      if (.not. (sched%target_kind(i_step) == target_by_load .or. &
        & sched%target_kind(i_step) == target_by_stress)) cycle
      step = config%file_name // ': step ' // int_text(i_step) // ', component ' // comp
      if (sched%is_load(i_step)) then
        rate = control_rate(config, i_step, rate_general)
        dt_max = largest_dt(config, i_step, rate, step)
        if (sched%dt_min(i_step) > dt_max) call fatal_error('cfg', step // ': dtime_min ' // &
          & real_text(sched%dt_min(i_step)) // ' s is larger than the largest increment ' // &
          & real_text(dt_max) // ' s')
      end if
      if (i_step == 1 .and. sched%has_prev) then
        ! A restarted run (step 8.4): from the previous run's last step, when it had a target
        ! of the same kind (as between any two steps); otherwise checked at the step's start.
        if (sched%prev%target_kind /= sched%target_kind(1)) cycle
        value_prev = sched%prev%target_value
        call check_apart(sched, i_step, config%file_name // ': the last step of the ' // &
          & 'restart point and step 1, component ' // comp, .false., value_prev)
      else if (i_step == 1) then
        value_prev = 0.0_rk
        call check_apart(sched, i_step, step, .true.)
      else if (sched%target_kind(i_step - 1) == sched%target_kind(i_step)) then
        value_prev = sched%target_value(i_step - 1)
        call check_apart(sched, i_step, config%file_name // ': steps ' // &
          & int_text(i_step - 1) // ' and ' // int_text(i_step) // ', component ' // comp, &
          & .false.)
      else
        cycle
      end if
      if (sched%is_rate(i_step)) then
        ! A load or stress rate (step 8.2c): its sign towards the target.
        call check_rate_direction(sched, i_step, sched%target_value(i_step) - value_prev, &
          & value_prev, step)
      else
        call check_direction(sched, config, i_step, value_prev, rate_general, step)
      end if
    end do

    ! Plan from time 0 and zero strains, up to and including the first load step, or the first
    ! step with load-controlled components (step 8.2b: their strains are measured); a step at a
    ! load or stress rate after step 1 starts from a value measured during the run (step 8.2c:
    ! planned then).
    ! (A restarted run, step 8.4: from the previous run's last step; a step 1 at a load or
    ! stress rate is planned at its start, from the value measured then.)
    t0 = 0.0_rk
    e0 = 0.0_rk
    total = 0.0_rk
    if (sched%has_prev) then
      t0 = sched%prev%time_end
      e0 = sched%prev%strain_end
    end if
    if (present(incrs_done)) total = real(incrs_done, rk)
    do i_step = 1, sched%num_steps
      if (sched%is_rate(i_step) .and. (i_step > 1 .or. sched%has_prev)) exit
      call step_schedule_plan(sched, config, i_step, t0, e0, rate_general, 0.0_rk)
      total = total + real(sched%num_incrs(i_step), rk)
      if (sched%is_load(i_step) .or. any(sched%bc_kind(:, i_step) /= 0)) exit
      t0 = sched%time_end(i_step)
      e0 = sched%strain_end(:, i_step)
    end do

    ! Run limits of the planned steps (v2-fixed C32: ending exactly at a limit is allowed).
    if (total > real(config%max_incr, rk)) call fatal_error('cfg', config%file_name // &
      & ': the steps need ' // real_text(total) // ' increments, more than max_incr ' // &
      & int_text(config%max_incr))
    sched%total_incrs = nint(total)
    if (t0 > config%max_total_time*(1.0_rk + tol_whole)) call fatal_error('cfg', &
      & config%file_name // ': the steps end at ' // real_text(t0) // ' s, after ' // &
      & 'max_total_time ' // real_text(config%max_total_time) // ' s')

  end subroutine step_schedule_build

  !> One step of a schedule as a step_record_t (step 8.4: the restart point's step), with the
  !> configuration's per_rate row of the step (zero without periodic conditions).
  !>
  !> Arguments:
  !>   sched   (in)  type(step_schedule_t)  the schedule
  !>   config  (in)  type(config_t)         its configuration
  !>   i_step  (in)  integer                the step (1..num_steps)
  !>   rec     (out, result) type(step_record_t)  the step
  function step_schedule_record(sched, config, i_step) result(rec)

    type(step_schedule_t), intent(in) :: sched
    type(config_t), intent(in) :: config
    integer, intent(in) :: i_step
    type(step_record_t) :: rec

    ! Every per-step component at the step.
    rec%target_kind = sched%target_kind(i_step)
    rec%is_load = sched%is_load(i_step)
    rec%target_value = sched%target_value(i_step)
    rec%target_tol = sched%target_tol(i_step)
    rec%load_sign = sched%load_sign(i_step)
    rec%dt_min = sched%dt_min(i_step)
    rec%strain_rate = sched%strain_rate(:, i_step)
    rec%vel_factor = sched%vel_factor(i_step)
    rec%time_start = sched%time_start(i_step)
    rec%time_end = sched%time_end(i_step)
    rec%strain_start = sched%strain_start(:, i_step)
    rec%strain_end = sched%strain_end(:, i_step)
    rec%num_incrs = sched%num_incrs(i_step)
    rec%dt = sched%dt(i_step)
    rec%dt_last = sched%dt_last(i_step)
    rec%print_results = sched%print_results(i_step)
    rec%bc_kind = sched%bc_kind(:, i_step)
    rec%bc_target = sched%bc_target(:, i_step)
    rec%bc_tol = sched%bc_tol(:, i_step)
    rec%bc_mode = sched%bc_mode(:, i_step)
    rec%bc_rate = sched%bc_rate(:, i_step)
    rec%rate_kind = sched%rate_kind(i_step)
    rec%is_rate = sched%is_rate(i_step)
    rec%is_dwell = sched%is_dwell(i_step)
    rec%is_open = sched%is_open(i_step)
    rec%ctrl_rate = sched%ctrl_rate(i_step)

    ! The periodic rates of the step.
    rec%per_rate = 0.0_rk
    if (allocated(config%per_rate)) rec%per_rate = config%per_rate(:, i_step)

  end function step_schedule_record

  !> Plan one step from its start (time and strains): its strain rates, end time and strains,
  !> and increments, with the checks of the module header (time target after the start,
  !> strain target different from the start strain and reached by the rate's sign, max_strain
  !> on the known end strains, max_incr on the step's increments). A load- or stress-target
  !> step gets its rates and its largest increment only (its end is reached during the run).
  !> A step at a load or stress rate of the control component (step 8.2c) is planned from the
  !> control value at its start: to a load or stress target its duration |v* - v_start|/|rate|
  !> (the rate towards the target, the target not already reached), to a time or a duration
  !> that time; a dwell gets its largest increment dtime (its increments are chosen during the
  !> run); a strain target at a load or stress rate its rates and largest increment only, with
  !> the rate's sign checked towards the strain target.
  !> Called by step_schedule_build for the steps known before the run, and during the run at
  !> the start of every later step.
  !>
  !> Arguments:
  !>   sched         (inout) type(step_schedule_t)  the schedule (allocated)
  !>   config        (in)    type(config_t)         the configuration
  !>   i_step        (in)    integer                the step
  !>   t0            (in)    real(rk)               time at the start of the step [s]
  !>   e0            (in)    real(rk)(3)            engineering strains of the prescribed
  !>                                                components at the start of the step [-]
  !>   rate_general  (in)    real(rk)               general: max |set_bc velocity| / gage
  !>                                                length [1/s]
  !>   value_start   (in, optional) real(rk)        the control value at the step's start (a
  !>                                                load or stress rate step; default 0) [N] or
  !>                                                [MPa]
  !> Locals:
  !>   v0            real(rk)                       value_start, or 0 [N] or [MPa]
  !>   is_preset     logical                        a preset (uniaxial or triaxial)
  !>   axis          integer                        control component (0: none)
  !>   c             integer                        a normal component
  !>   step          character(:), allocatable      "<file>: step <i>[, component <c>]"
  !>   comp          character(:), allocatable      the loading component as text
  !>   rate          real(rk)                       control strain rate of the step [1/s]
  !>   de            real(rk)                       strain change to a strain target [-]
  !>   duration      real(rk)                       length of the step [s]
  !>   dt_given      real(rk)                       increment size given (dtime, or dstrain
  !>                                                over |rate|) [s]
  !>   incr_text     character(:), allocatable      the increments of a step, for the warning
  subroutine step_schedule_plan(sched, config, i_step, t0, e0, rate_general, value_start)

    type(step_schedule_t), intent(inout) :: sched
    type(config_t), intent(in) :: config
    integer, intent(in) :: i_step
    real(rk), intent(in) :: t0, e0(3), rate_general
    real(rk), intent(in), optional :: value_start

    integer :: axis, c
    logical :: is_preset
    character(len=:), allocatable :: step, comp, incr_text
    real(rk) :: rate, de, duration, dt_given, v0

    is_preset = config%bc_type /= 'general'
    axis = sched%axis
    comp = ''
    if (axis > 0) comp = trim(sched%comp_label(axis))
    step = config%file_name // ': step ' // int_text(i_step)
    if (axis > 0) step = step // ', component ' // comp
    sched%num_planned = max(sched%num_planned, i_step)
    sched%time_start(i_step) = t0
    sched%strain_start(:, i_step) = e0

    ! Strain rate of the step: the strain_rate row with a preset; with general, the rate of
    ! the set_bc velocities, signed towards a strain target (v2's target_sign), positive for
    ! a time target (v2); a load target's sign is set at the step start (vel_factor, kept).
    if (.not. sched%is_load(i_step)) sched%vel_factor(i_step) = 1.0_rk
    rate = control_rate(config, i_step, rate_general)
    if (.not. is_preset .and. axis > 0) then
      rate = sched%vel_factor(i_step)*rate
      if (config%has_target_strain(axis, i_step)) then
        if (config%target_strain(axis, i_step) < e0(axis)) then
          sched%vel_factor(i_step) = -1.0_rk
          rate = -rate_general
        end if
      end if
    end if

    v0 = 0.0_rk
    if (present(value_start)) v0 = value_start

    ! A load or stress target, or a strain target at a load or stress rate (step 8.2c, its
    ! rate's sign towards the strain target): the rates and the largest increment; the end is
    ! reached during the run.
    if (sched%is_open(i_step)) then
      sched%load_sign(i_step) = sign(1.0_rk, sched%target_value(i_step) - e0(axis))
      if (.not. abs(sched%target_value(i_step) - e0(axis)) > 0.0_rk) call fatal_error('cfg', &
        & step // ': target_strain ' // comp // ' ' // real_text(sched%target_value(i_step)) &
        & // ' equals the strain at the start of the step')
      call check_rate_direction(sched, i_step, sched%target_value(i_step) - e0(axis), e0(axis), &
        & step)
    end if
    if (sched%is_load(i_step) .or. sched%is_open(i_step)) then
      sched%strain_rate(axis, i_step) = rate
      do c = 1, 3
        if (c == axis .or. .not. config%comp_rows(c)) cycle
        sched%strain_rate(c, i_step) = config%strain_rate(c, i_step)
      end do
      sched%time_end(i_step) = t0
      sched%strain_end(:, i_step) = e0
      sched%num_incrs(i_step) = 0
      sched%dt(i_step) = largest_dt(config, i_step, rate, step)
      sched%dt_last(i_step) = sched%dt(i_step)
      return
    end if

    ! Length of the step from its target.
    if (config%has_target_time(i_step)) then
      ! Time target (cumulative): after the start of the step; the strain follows the rate.
      duration = config%target_time(i_step) - t0
      if (.not. duration > 0.0_rk) call fatal_error('cfg', step // ': target_time ' // &
        & real_text(config%target_time(i_step)) // ' s is not after the start of the ' // &
        & 'step (' // real_text(t0) // ' s)')
      de = rate*duration
    else if (config%has_target_duration(i_step)) then
      ! Duration target (step 8.2c, Q-T7): the step's length from its actual start.
      duration = config%target_duration(i_step)
      de = rate*duration
    else if (sched%is_rate(i_step)) then
      ! A load or stress target at a load or stress rate (step 8.2c, addendum §3.3): from the
      ! value at the start, not already reached, the rate towards the target.
      call check_reached(sched, i_step, v0, step)
      call check_rate_direction(sched, i_step, sched%target_value(i_step) - v0, v0, step)
      sched%load_sign(i_step) = sign(1.0_rk, sched%target_value(i_step) - v0)
      duration = abs(sched%target_value(i_step) - v0)/abs(sched%ctrl_rate(i_step))
      de = 0.0_rk
    else
      ! Strain target: a change from the start value, reached at the step's rate.
      de = config%target_strain(axis, i_step) - e0(axis)
      if (.not. abs(de) > 0.0_rk) call fatal_error('cfg', step // ': target_strain ' // &
        & comp // ' ' // real_text(config%target_strain(axis, i_step)) // &
        & ' equals the strain at the start of the step')
      if (.not. rate*de > 0.0_rk) call fatal_error('cfg', step // ': strain_rate ' // &
        & comp // ' ' // real_text(rate) // ' moves away from target_strain ' // comp // &
        & ' ' // real_text(config%target_strain(axis, i_step)) // ' (strain at the ' // &
        & 'start of the step: ' // real_text(e0(axis)) // ')')
      duration = de/rate
    end if
    if (config%has_target_time(i_step)) then
      sched%time_end(i_step) = config%target_time(i_step)
    else
      sched%time_end(i_step) = t0 + duration
    end if

    ! The control component: its rate, and its end strain (the target, or start + de; at a
    ! load or stress rate not planned, measured during the run).
    sched%strain_end(:, i_step) = 0.0_rk
    if (axis > 0) then
      sched%strain_rate(axis, i_step) = rate
      if (sched%target_kind(i_step) == target_by_strain) then
        sched%strain_end(axis, i_step) = config%target_strain(axis, i_step)
      else
        sched%strain_end(axis, i_step) = e0(axis) + de
      end if
    end if

    ! The other components with rows (triaxial, D45): a rate gives the end strain, a target
    ! gives the rate over the step's duration; a load-controlled component's strain is not
    ! planned (step 8.2b: measured during the run).
    do c = 1, 3
      if (c == axis .or. .not. config%comp_rows(c)) cycle
      if (sched%bc_kind(c, i_step) /= 0) then
        sched%strain_rate(c, i_step) = 0.0_rk
        sched%strain_end(c, i_step) = e0(c)
      else if (config%has_strain_rate(c, i_step)) then
        sched%strain_rate(c, i_step) = config%strain_rate(c, i_step)
        sched%strain_end(c, i_step) = e0(c) + config%strain_rate(c, i_step)*duration
      else
        sched%strain_end(c, i_step) = config%target_strain(c, i_step)
        sched%strain_rate(c, i_step) = (config%target_strain(c, i_step) - e0(c))/duration
      end if
    end do

    ! max_strain: an error when an end strain is known (presets: every prescribed
    ! component; general: strain targets; general time-target strains and free components
    ! are checked during the run, step 5.2).
    if (config%has_max_strain .and. axis > 0) then
      do c = 1, 3
        if (.not. sched%prescribed(c) .or. sched%bc_kind(c, i_step) /= 0) cycle
        if (.not. (is_preset .or. .not. config%has_target_time(i_step))) cycle
        if (abs(sched%strain_end(c, i_step)) > config%max_strain) call fatal_error('cfg', &
          & config%file_name // ': step ' // int_text(i_step) // ', component ' // &
          & trim(sched%comp_label(c)) // ': the strain at the end of the step, ' // &
          & real_text(sched%strain_end(c, i_step)) // ', exceeds max_strain ' // &
          & real_text(config%max_strain))
      end do
    end if

    ! Increments of the step (a dwell: chosen during the run, at most dtime, step 8.2c).
    if (sched%is_dwell(i_step)) then
      sched%num_incrs(i_step) = 0
      sched%dt(i_step) = config%dtime(i_step)
      sched%dt_last(i_step) = config%dtime(i_step)
      return
    end if
    select case (config%incr_kind(i_step))
    case (incr_by_count)
      ! number_of_incrs: equal increments.
      sched%num_incrs(i_step) = config%num_incrs(i_step)
      sched%dt(i_step) = duration/real(config%num_incrs(i_step), rk)
      sched%dt_last(i_step) = sched%dt(i_step)
    case default
      ! dtime, or dstrain converted to a time increment with the step's strain rate.
      dt_given = largest_dt(config, i_step, rate, step)

      ! Too many increments for max_incr (also guards the integer conversion below).
      if (duration/dt_given > real(config%max_incr, rk)) call fatal_error('cfg', step // &
        & ': ' // real_text(duration/dt_given) // ' increments of ' // real_text(dt_given) &
        & // ' s exceed max_incr ' // int_text(config%max_incr))
      call split_step(duration, dt_given, sched%num_incrs(i_step), sched%dt(i_step), &
        & sched%dt_last(i_step))

      ! A shortened last increment is reported (D44): "<n-1> x <dt> s + <dt_last> s", or
      ! the single increment of a step shorter than dt_given.
      if (sched%dt_last(i_step) < sched%dt(i_step)) then
        if (sched%num_incrs(i_step) == 1) then
          incr_text = '1 increment of ' // real_text(sched%dt_last(i_step)) // ' s'
        else
          incr_text = int_text(sched%num_incrs(i_step)) // ' increments: ' // &
            & int_text(sched%num_incrs(i_step) - 1) // ' x ' // &
            & real_text(sched%dt(i_step)) // ' s + ' // &
            & real_text(sched%dt_last(i_step)) // ' s'
        end if
        call log_warning(step // ': ' // real_text(duration) // ' s is not a whole ' // &
          & 'number of increments of ' // real_text(dt_given) // ' s; ' // incr_text // '.', 2)
      end if
    end select

  end subroutine step_schedule_plan

  !> Start of a load- or stress-target step during the run (addendum §3.3): its direction
  !> s = sign(target - value), checked against the strain rate with a preset (an error when
  !> the rate moves away from the target), or setting the set_bc velocity factor with general
  !> (v2: reversed for a target below the start value); a target already reached within its
  !> tolerance is an error.
  !>
  !> Arguments:
  !>   sched         (inout) type(step_schedule_t)  the schedule (step planned)
  !>   config        (in)    type(config_t)         the configuration
  !>   i_step        (in)    integer                the step (a load or stress target)
  !>   value         (in)    real(rk)               the current value: force on the control
  !>                                                face [N] or that force over its area
  !>                                                [MPa]
  !>   rate_general  (in)    real(rk)               general: max |set_bc velocity| / gage
  !>                                                length [1/s]
  !> Locals:
  !>   step          character(:), allocatable      "<file>: step <i>, component <c>"
  subroutine step_schedule_load_start(sched, config, i_step, value, rate_general)

    type(step_schedule_t), intent(inout) :: sched
    type(config_t), intent(in) :: config
    integer, intent(in) :: i_step
    real(rk), intent(in) :: value, rate_general

    character(len=:), allocatable :: step

    step = config%file_name // ': step ' // int_text(i_step) // ', component ' // &
      & trim(sched%comp_label(sched%axis))
    call check_direction(sched, config, i_step, value, rate_general, step)

  end subroutine step_schedule_load_start

  !> Direction of a load- or stress-target step from a start value (addendum §3.3): sets
  !> load_sign; presets: the strain rate must have that sign (Q-L6); general: vel_factor and
  !> the control strain rate take it (v2); a target within its tolerance of the start value is
  !> an error (v2: "Load differential is zero between steps").
  !>
  !> Arguments:
  !>   sched         (inout) type(step_schedule_t)  the schedule
  !>   config        (in)    type(config_t)         the configuration
  !>   i_step        (in)    integer                the step (a load or stress target)
  !>   value         (in)    real(rk)               the start value [N] or [MPa]
  !>   rate_general  (in)    real(rk)               general: max |set_bc velocity| / gage
  !>                                                length [1/s]
  !>   step          (in)    character              "<file>: step <i>, component <c>"
  !> Locals:
  !>   name          character(:), allocatable      the target row ("target_load 33")
  !>   unit          character(:), allocatable      its unit
  !>   rate          real(rk)                       the preset's strain rate [1/s]
  subroutine check_direction(sched, config, i_step, value, rate_general, step)

    type(step_schedule_t), intent(inout) :: sched
    type(config_t), intent(in) :: config
    integer, intent(in) :: i_step
    real(rk), intent(in) :: value, rate_general
    character(len=*), intent(in) :: step

    character(len=:), allocatable :: name, unit
    real(rk) :: rate

    if (sched%target_kind(i_step) == target_by_load) then
      name = 'target_load ' // trim(sched%comp_label(sched%axis))
      unit = ' N'
    else
      name = 'target_stress ' // trim(sched%comp_label(sched%axis))
      unit = ' MPa'
    end if

    ! A target already reached.
    if (.not. abs(sched%target_value(i_step) - value) > sched%target_tol(i_step)) &
      & call fatal_error('cfg', step // ': ' // name // ' ' // &
      & real_text(sched%target_value(i_step)) // unit // ' is already reached at the start ' &
      & // 'of the step (' // real_text(value) // unit // ')')
    sched%load_sign(i_step) = sign(1.0_rk, sched%target_value(i_step) - value)

    if (config%bc_type == 'general') then
      ! The set_bc velocities, reversed for a target below the start value (v2).
      sched%vel_factor(i_step) = sched%load_sign(i_step)
      sched%strain_rate(sched%axis, i_step) = sched%load_sign(i_step)*rate_general
    else
      ! The strain rate must move the value towards the target (a positive rate increases
      ! the load).
      rate = config%strain_rate(sched%axis, i_step)
      if (.not. rate*sched%load_sign(i_step) > 0.0_rk) call fatal_error('cfg', step // &
        & ': strain_rate ' // trim(sched%comp_label(sched%axis)) // ' ' // real_text(rate) // &
        & ' moves away from ' // name // ' ' // real_text(sched%target_value(i_step)) // &
        & unit // ' (value at the start of the step: ' // real_text(value) // unit // &
        & '; a positive strain rate increases the load)')
    end if

  end subroutine check_direction

  !> A load or stress target already reached at the start of a step at a load or stress rate
  !> (step 8.2c; the same rule as check_direction's): an error.
  !>
  !> Arguments:
  !>   sched   (in)  type(step_schedule_t)  the schedule
  !>   i_step  (in)  integer                the step (a load or stress target)
  !>   value   (in)  real(rk)               the value at the start [N] or [MPa]
  !>   step    (in)  character              "<file>: step <i>, component <c>"
  !> Locals:
  !>   name, unit  character(:), allocatable  the target row and its unit
  subroutine check_reached(sched, i_step, value, step)

    type(step_schedule_t), intent(in) :: sched
    integer, intent(in) :: i_step
    real(rk), intent(in) :: value
    character(len=*), intent(in) :: step

    character(len=:), allocatable :: name, unit

    call target_name(sched, i_step, name, unit)
    if (.not. abs(sched%target_value(i_step) - value) > sched%target_tol(i_step)) &
      & call fatal_error('cfg', step // ': ' // name // ' ' // &
      & real_text(sched%target_value(i_step)) // unit // ' is already reached at the start ' &
      & // 'of the step (' // real_text(value) // unit // ')')

  end subroutine check_reached

  !> Direction of a load or stress rate of the control component (step 8.2c, addendum §3.2,
  !> the sign rule Q-L6): the rate must have the sign of the change to its target (a load or
  !> stress target, or a strain target: a positive rate increases the strain); an error naming
  !> the rate, the target and the start value otherwise.
  !>
  !> Arguments:
  !>   sched   (in)  type(step_schedule_t)  the schedule
  !>   i_step  (in)  integer                the step (at a load or stress rate)
  !>   change  (in)  real(rk)               target minus start value [N], [MPa] or [-]
  !>   value   (in)  real(rk)               the start value [N], [MPa] or [-]
  !>   step    (in)  character              "<file>: step <i>, component <c>"
  !> Locals:
  !>   name, unit  character(:), allocatable  the target row and its unit
  !>   rate_name   character(:), allocatable  the rate row
  subroutine check_rate_direction(sched, i_step, change, value, step)

    type(step_schedule_t), intent(in) :: sched
    integer, intent(in) :: i_step
    real(rk), intent(in) :: change, value
    character(len=*), intent(in) :: step

    character(len=:), allocatable :: name, unit, rate_name

    call target_name(sched, i_step, name, unit)
    rate_name = 'load_rate '
    if (sched%rate_kind(i_step) == rate_by_stress) rate_name = 'stress_rate '
    rate_name = rate_name // trim(sched%comp_label(sched%axis))
    if (.not. sched%ctrl_rate(i_step)*change > 0.0_rk) call fatal_error('cfg', step // ': ' // &
      & rate_name // ' ' // real_text(sched%ctrl_rate(i_step)) // ' moves away from ' // name &
      & // ' ' // real_text(sched%target_value(i_step)) // unit // ' (value at the start of ' &
      & // 'the step: ' // real_text(value) // unit // ')')

  end subroutine check_rate_direction

  !> The target row of a step's control component and its unit, for messages.
  !>
  !> Arguments:
  !>   sched   (in)  type(step_schedule_t)       the schedule
  !>   i_step  (in)  integer                     the step
  !>   name    (out) character(:), allocatable   "target_load 33", "target_stress 33",
  !>                                             "target_strain 33"
  !>   unit    (out) character(:), allocatable   " N", " MPa", ""
  subroutine target_name(sched, i_step, name, unit)

    type(step_schedule_t), intent(in) :: sched
    integer, intent(in) :: i_step
    character(len=:), allocatable, intent(out) :: name, unit

    select case (sched%target_kind(i_step))
    case (target_by_load)
      name = 'target_load '
      unit = ' N'
    case (target_by_stress)
      name = 'target_stress '
      unit = ' MPa'
    case default
      name = 'target_strain '
      unit = ''
    end select
    name = name // trim(sched%comp_label(sched%axis))

  end subroutine target_name

  !> Before the run (Q-T16, addendum_triaxial_load_control.md §3.2): step 1's load or stress
  !> target more than one tolerance from the initial value 0, or a target more than two
  !> tolerances from the previous step's target of the same kind; an error naming the step(s)
  !> otherwise.
  !>
  !> Arguments:
  !>   sched       (in)  type(step_schedule_t)  the schedule (targets set)
  !>   i_step      (in)  integer                the step (a load or stress target; > 1 when
  !>                                            not from_start)
  !>   step        (in)  character              the step(s), for the message
  !>   from_start  (in)  logical                step 1 from 0 (else from the previous target)
  !>   prev_value  (in, optional) real(rk)      the previous target (default: that of step
  !>                                            i_step - 1; a restarted run's step 1: the
  !>                                            previous run's last target, step 8.4)
  !> Locals:
  !>   name        character(:), allocatable    the target row ("target_load 33")
  !>   unit        character(:), allocatable    its unit
  !>   tol_name    character(:), allocatable    its tolerance keyword
  subroutine check_apart(sched, i_step, step, from_start, prev_value)

    type(step_schedule_t), intent(in) :: sched
    integer, intent(in) :: i_step
    character(len=*), intent(in) :: step
    logical, intent(in) :: from_start
    real(rk), intent(in), optional :: prev_value

    character(len=:), allocatable :: name, unit, tol_name
    real(rk) :: v_prev

    if (sched%target_kind(i_step) == target_by_load) then
      name = 'target_load ' // trim(sched%comp_label(sched%axis))
      unit = ' N'
      tol_name = 'load_tol'
    else
      name = 'target_stress ' // trim(sched%comp_label(sched%axis))
      unit = ' MPa'
      tol_name = 'stress_tol'
    end if

    if (from_start) then
      ! Step 1 from the initial value 0.
      if (.not. step_load_targets_apart(0.0_rk, sched%target_value(i_step), &
        & sched%target_tol(i_step), 1)) call fatal_error('cfg', step // ': ' // name // ' ' &
        & // real_text(sched%target_value(i_step)) // unit // ' is not more than ' // &
        & tol_name // ' ' // real_text(sched%target_tol(i_step)) // unit // ' from the ' // &
        & 'initial value 0')
    else
      ! Consecutive targets: the first step ends within one tolerance of its target.
      if (present(prev_value)) then
        v_prev = prev_value
      else
        v_prev = sched%target_value(i_step - 1)
      end if
      if (.not. step_load_targets_apart(v_prev, &
        & sched%target_value(i_step), sched%target_tol(i_step), 2)) call fatal_error('cfg', &
        & step // ': ' // name // ' ' // real_text(v_prev) // &
        & ' and ' // real_text(sched%target_value(i_step)) // unit // ' are not more than ' &
        & // 'twice ' // tol_name // ' (2 x ' // real_text(sched%target_tol(i_step)) // unit // &
        & ') apart (a step ends within ' // tol_name // ' of its target)')
    end if

  end subroutine check_apart

  !> Completion test of a load or stress target (v2 loading_load_isinrange): the value has
  !> reached the target, or is within the tolerance short of it, in the step's direction. Since
  !> step 8.2a it tests only the extrapolated value of step_load_next_dt; a step's completion
  !> is step_load_landing (the landing rule).
  !>
  !> Arguments:
  !>   value      (in)  real(rk)  current value [N] or [MPa]
  !>   target     (in)  real(rk)  target [N] or [MPa]
  !>   load_sign  (in)  real(rk)  direction of the step, +1 or -1
  !>   tol        (in)  real(rk)  tolerance, >= 0 [N] or [MPa]
  pure logical function step_load_complete(value, target, load_sign, tol)

    real(rk), intent(in) :: value, target, load_sign, tol

    step_load_complete = (value - target)*load_sign > -tol

  end function step_load_complete

  !> Next increment of a load- or stress-target step that is not complete (v2
  !> loading_timeincr; addendum §2 item 2): dt_max, unless the value extrapolated linearly
  !> over dt_max with the last increment's rate would complete the step; then
  !> f dt_max (target - value)/(trial - value) = f (target - value)/vdot, clamped to
  !> [dt_min, dt_max] (an extrapolated change below 1e-16 keeps dt_max, v2).
  !>
  !> Arguments:
  !>   value       (in)  real(rk)  value after the last increment [N] or [MPa]
  !>   value_prev  (in)  real(rk)  value before it [N] or [MPa]
  !>   dt_done     (in)  real(rk)  length of the last increment, > 0 [s]
  !>   target      (in)  real(rk)  the step's target [N] or [MPa]
  !>   load_sign   (in)  real(rk)  direction of the step, +1 or -1
  !>   tol         (in)  real(rk)  tolerance of the target [N] or [MPa]
  !>   dt_max      (in)  real(rk)  largest increment of the step [s]
  !>   dt_min      (in)  real(rk)  smallest increment of the step [s]
  !>   factor      (in)  real(rk)  dtime_factor, >= 1 [-]
  !>   dt          (out, result) real(rk)  the next increment [s]
  !> Locals:
  !>   rate        real(rk)  value rate of the last increment [N/s] or [MPa/s]
  !>   trial       real(rk)  value extrapolated over dt_max [N] or [MPa]
  pure function step_load_next_dt(value, value_prev, dt_done, target, load_sign, tol, &
    & dt_max, dt_min, factor) result(dt)

    real(rk), intent(in) :: value, value_prev, dt_done, target, load_sign, tol, dt_max, &
      & dt_min, factor
    real(rk) :: dt

    real(rk) :: rate, trial

    ! The largest increment, unless the extrapolation reaches the target.
    dt = dt_max
    rate = (value - value_prev)/dt_done
    trial = value + rate*dt_max
    if (.not. step_load_complete(trial, target, load_sign, tol)) return
    if (abs(trial - value) < load_change_min) return

    ! Aim slightly past the target (factor >= 1), within [dt_min, dt_max].
    dt = factor*dt_max*(target - value)/(trial - value)
    dt = max(dt_min, min(dt_max, dt))

  end function step_load_next_dt

  !> Where a value lies with respect to a load or stress target (landing rule, step 8.2a,
  !> addendum_triaxial_load_control.md §3.6): short of it by more than the tolerance
  !> (s (v - v*) < -tol: the step goes on), within the tolerance (|v - v*| <= tol: the step is
  !> complete), or past it by more than the tolerance (s (v - v*) > tol: the increment is
  !> redone shorter).
  !>
  !> Arguments:
  !>   value      (in)  real(rk)  value after the increment [N] or [MPa]
  !>   target     (in)  real(rk)  the step's target [N] or [MPa]
  !>   load_sign  (in)  real(rk)  direction of the step, +1 or -1
  !>   tol        (in)  real(rk)  tolerance, > 0 [N] or [MPa]
  !>   position   (out, result) integer  landing_short, landing_within or landing_passed
  !> Locals:
  !>   ahead      real(rk)  s (v - v*): the signed distance past the target [N] or [MPa]
  pure function step_load_landing(value, target, load_sign, tol) result(position)

    real(rk), intent(in) :: value, target, load_sign, tol
    integer :: position

    real(rk) :: ahead

    ahead = (value - target)*load_sign
    if (ahead < -tol) then
      position = landing_short
    else if (ahead > tol) then
      position = landing_passed
    else
      position = landing_within
    end if

  end function step_load_landing

  !> Length of the next landing attempt by false position (Press et al. 2007, Numerical
  !> Recipes 3rd ed., §9.2; addendum_triaxial_load_control.md §3.6): the root of the chord
  !> through the bracket (dt_lo, v_lo), (dt_hi, v_hi) of increment lengths and the values they
  !> reached, v_lo short of the target and v_hi past it (each by more than the tolerance, so
  !> v_hi /= v_lo):
  !>   dt = dt_lo + (dt_hi - dt_lo) (v* - v_lo)/(v_hi - v_lo),  dt_lo < dt < dt_hi.
  !>
  !> Arguments:
  !>   dt_lo     (in)  real(rk)  length of the bracket's short end (0: the start state) [s]
  !>   value_lo  (in)  real(rk)  value reached with it [N] or [MPa]
  !>   dt_hi     (in)  real(rk)  length of the bracket's passed end [s]
  !>   value_hi  (in)  real(rk)  value reached with it [N] or [MPa]
  !>   target    (in)  real(rk)  the step's target [N] or [MPa]
  !>   dt        (out, result) real(rk)  length of the next attempt [s]
  pure function step_load_false_position(dt_lo, value_lo, dt_hi, value_hi, target) result(dt)

    real(rk), intent(in) :: dt_lo, value_lo, dt_hi, value_hi, target
    real(rk) :: dt

    dt = dt_lo + (dt_hi - dt_lo)*(target - value_lo)/(value_hi - value_lo)

  end function step_load_false_position

  !> Whether two load or stress values are more than a number of tolerances apart (Q-T16,
  !> addendum_triaxial_load_control.md §3.2): |target - target_prev| > num_tol tol, with
  !> num_tol 1 for step 1's target from the initial value 0, 2 for consecutive targets.
  !>
  !> Arguments:
  !>   target_prev  (in)  real(rk)  previous target, or the initial value [N] or [MPa]
  !>   target       (in)  real(rk)  the target [N] or [MPa]
  !>   tol          (in)  real(rk)  tolerance, > 0 [N] or [MPa]
  !>   num_tol      (in)  integer   tolerances needed between the two (1 or 2)
  pure logical function step_load_targets_apart(target_prev, target, tol, num_tol)

    real(rk), intent(in) :: target_prev, target, tol
    integer, intent(in) :: num_tol

    step_load_targets_apart = abs(target - target_prev) > real(num_tol, rk)*tol

  end function step_load_targets_apart

  !> Next increment of a dwell (step 8.2c; dev/theory/addendum_triaxial_load_control.md §3.3;
  !> v2 driver_triaxclr): 0.01 dtime at the start of a dwell episode (a dwell that follows a
  !> step that is not a dwell), otherwise min(max_strain_incr L0 / |V|, dtime) (dtime when
  !> V = 0), V the control face velocity and L0 the original length along the control
  !> component (the strain increment of the control); at most the remaining time, and the
  !> whole remaining time when less than 0.001 dtime would remain (v2 ended the step there, short
  !> of its time; FEPX ends the step at its time exactly, D44).
  !>
  !> Arguments:
  !>   episode_start    (in)  logical   the increment starts a dwell episode
  !>   vel_ctrl         (in)  real(rk)  control face velocity of the last increment [mm/s]
  !>   length0          (in)  real(rk)  original length along the control component [mm]
  !>   max_strain_incr  (in)  real(rk)  largest strain increment of the control, > 0 [-]
  !>   dtime            (in)  real(rk)  the step's dtime, > 0 [s]
  !>   remaining        (in)  real(rk)  time left in the step, > 0 [s]
  !>   dt               (out) real(rk)  the increment [s]
  !>   last             (out) logical   the increment ends the step
  pure subroutine step_dwell_dt(episode_start, vel_ctrl, length0, max_strain_incr, dtime, &
    & remaining, dt, last)

    logical, intent(in) :: episode_start
    real(rk), intent(in) :: vel_ctrl, length0, max_strain_incr, dtime, remaining
    real(rk), intent(out) :: dt
    logical, intent(out) :: last

    ! v2's increment: a small first one, then bounded by the strain increment and dtime.
    if (episode_start) then
      dt = dwell_first_frac*dtime
    else if (abs(vel_ctrl) > 0.0_rk) then
      dt = min(max_strain_incr*length0/abs(vel_ctrl), dtime)
    else
      dt = dtime
    end if

    ! The remaining time, all of it when little would be left.
    last = .not. remaining - dt > dwell_time_tol*dtime
    if (last) dt = remaining

  end subroutine step_dwell_dt

  !> The control component's strain rate of a step: the strain_rate row with a preset; with
  !> general the set_bc rate (positive; its sign follows the target).
  !>
  !> Arguments:
  !>   config        (in)  type(config_t)  the configuration
  !>   i_step        (in)  integer         the step
  !>   rate_general  (in)  real(rk)        general: max |set_bc velocity| / gage length [1/s]
  !>   rate          (out, result) real(rk)  the rate [1/s]
  pure function control_rate(config, i_step, rate_general) result(rate)

    type(config_t), intent(in) :: config
    integer, intent(in) :: i_step
    real(rk), intent(in) :: rate_general
    real(rk) :: rate

    if (config%bc_type /= 'general') then
      rate = config%strain_rate(config%control_axis, i_step)
    else if (config%control_axis > 0) then
      rate = rate_general
    else
      rate = 0.0_rk
    end if

  end function control_rate

  !> The increment size a step asks for: dtime, or dstrain over |rate| (the control
  !> component's rate); for a load or stress target the largest increment.
  !>
  !> Arguments:
  !>   config  (in)  type(config_t)  the configuration
  !>   i_step  (in)  integer         the step (dtime or dstrain given)
  !>   rate    (in)  real(rk)        control strain rate of the step [1/s]
  !>   step    (in)  character       the step, for messages
  !>   dt      (out, result) real(rk)  the increment [s]
  function largest_dt(config, i_step, rate, step) result(dt)

    type(config_t), intent(in) :: config
    integer, intent(in) :: i_step
    real(rk), intent(in) :: rate
    character(len=*), intent(in) :: step
    real(rk) :: dt

    if (config%incr_kind(i_step) == incr_by_dtime) then
      dt = config%dtime(i_step)
    else
      if (config%control_axis == 0) call fatal_error('cfg', step // ': dstrain needs a ' // &
        & 'loading component (a target_strain row); use number_of_incrs or dtime')
      if (.not. abs(rate) > 0.0_rk) call fatal_error('cfg', step // ': dstrain with ' // &
        & 'a zero strain rate; use number_of_incrs or dtime')
      dt = config%dstrain(i_step)/abs(rate)
    end if

  end function largest_dt

  !> Cut a step of a given length into increments of a given size: full increments and a
  !> shorter last one that ends the step exactly (D44). A step within round-off of a whole
  !> number n of increments (|T/dt - n| <= tol_whole T/dt) gets n equal increments T/n.
  !>
  !> Arguments:
  !>   duration   (in)  real(rk)  length of the step, > 0 [s]
  !>   dt_given   (in)  real(rk)  increment size asked for, > 0 [s]
  !>   num_incrs  (out) integer   number of increments (>= 1)
  !>   dt         (out) real(rk)  size of the increments but the last [s]
  !>   dt_last    (out) real(rk)  size of the last increment (<= dt; < dt when the step is not
  !>                              a whole number of increments) [s]
  !> Locals:
  !>   ratio      real(rk)        duration over dt_given (increments, as a real)
  !>   n_near     integer         nearest whole number of increments
  pure subroutine split_step(duration, dt_given, num_incrs, dt, dt_last)

    real(rk), intent(in) :: duration, dt_given
    integer, intent(out) :: num_incrs
    real(rk), intent(out) :: dt, dt_last

    real(rk) :: ratio
    integer :: n_near

    ! Increments needed, as a real, and the nearest whole number (at least one).
    ratio = duration/dt_given
    n_near = max(1, nint(ratio))

    if (abs(ratio - real(n_near, rk)) <= tol_whole*ratio) then
      ! A whole number of increments: equal increments ending exactly on the target.
      num_incrs = n_near
      dt = duration/real(num_incrs, rk)
      dt_last = dt
    else
      ! Full increments, then the remainder (between 0 and dt_given, exclusive) as the last
      ! (a step shorter than dt_given is one increment of its own length).
      num_incrs = ceiling(ratio)
      dt = dt_given
      dt_last = duration - real(num_incrs - 1, rk)*dt_given
    end if

  end subroutine split_step

  !> Decimal text of an integer.
  !>
  !> Arguments:
  !>   value  (in)  integer  the integer
  !> Locals:
  !>   buf    character(24)  formatting buffer
  pure function int_text(value) result(text)

    integer, intent(in) :: value
    character(len=:), allocatable :: text

    character(len=24) :: buf

    write (buf, '(i0)') value
    text = trim(buf)

  end function int_text

  !> Short text of a real for messages (scientific, 5 significant digits, no blanks).
  !>
  !> Arguments:
  !>   value  (in)  real(rk)  the real
  !> Locals:
  !>   buf    character(32)   formatting buffer
  pure function real_text(value) result(text)

    real(rk), intent(in) :: value
    character(len=:), allocatable :: text

    character(len=32) :: buf

    write (buf, '(es12.4)') value
    text = trim(adjustl(buf))

  end function real_text

end module step_schedule_mod
