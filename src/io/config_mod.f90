! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The typed configuration: built from a validated cfg store (cfg_store_mod) with the checks that
!> involve several keywords (dev/cfg_keywords.md §1a, for the keywords active in step 4.3).
!>
!> Source : FEPX (ROADMAP steps 4.3, 5.1, 5.3 and 6.1, D42, D45, D48; dev/cfg_keywords.md,
!>          io_format.md
!>          §2-3). The positive
!>          definiteness test of the elastic constants uses the LDL^T factorization of
!>          linalg_mod::ldlt_factor (a symmetric matrix is positive definite if and only if every
!>          pivot of its LDL^T factorization is positive; Golub & Van Loan, Matrix Computations,
!>          4th ed., 2013, §4.1-4.2). Contains no other numerical method.
!>
!> Checks (fatal, "<file>:<line>: <message>" or "<file>: <message>" for a missing keyword):
!>   - required keywords: number_of_phases, boundary_conditions, number_of_steps; per phase
!>     crystal_type and model;
!>   - phases: every phase 1..number_of_phases defined, none outside that range;
!>   - elastic constants: exactly those of the crystal type (cubic c11 c12 c44; hcp adds c13
!>     c33; bct adds c13 c33 c66), positive definite; c_over_a only for hcp and bct;
!>   - slip (steps 6.1, 6.2, dev/cfg_keywords.md §2.1, ROADMAP D48, D49): slip_families, m,
!>     gammadot_0 and tau_0 only with model crystal_plasticity, which requires m, gammadot_0
!>     and tau_0 (and c_over_a for
!>     hcp and bct); slip_families tokens matched exactly against the crystal type's families,
!>     or "all" alone, each family once; m, gammadot_0, tau_0 with one value per active family (in
!>     the slip_families order) or one per slip system (in the slip-system order), tau_0
!>     likewise;
!>   - boundary conditions and control rows (§1 table): presets (uniaxial, triaxial) take no
!>     set_bc lines; a uniaxial preset has all control rows on one normal component (the
!>     loading axis); triaxial may have rows on several; general needs set_bc vel lines, takes
!>     no strain_rate row, and its target_strain rows on one normal component (the loading
!>     axis, v2's loading direction); set_bc syntax (vel <nset> <x|y|z> <v>, up to three
!>     pairs), one velocity per node set and direction;
!>   - multi-point constraints (step 8.3a, D73, dev/theory/addendum_periodic_mpc.md §3.2):
!>     set_mpc1 vel <nset> <x|y|z> [<x|y|z> [<x|y|z>]] and set_mpc vel <nset> <x|y|z> <coeff>
!>     <offset> (up to three triples); v2's set_mpc without vel is an error naming the form;
!>     a nonzero coefficient (v2); one relation per node set and direction; not with
!>     boundary_conditions triaxial, nor with a load or stress target or rate (Q-P8). The
!>     node sets, overlaps and chains are checked against the mesh by constraint_mod;
!>   - periodic conditions (step 8.3b, D71-D73, dev/theory/addendum_periodic_mpc.md §3.3, §4):
!>     periodic_kinematics and periodic_preconditioner only with boundary_conditions periodic;
!>     shear components only with it (small_strain: 23 13 12; deformation_gradient_rate: the
!>     nine); a strain_rate row (a value in every step) per imposed component; control_component
!>     required, imposed, the only component with target_strain rows; one target per step
!>     (its target_strain, target_time or target_duration); no load or stress rows, no
!>     multi-point constraints (v2-fixed C68), no check_necking; print forces ignored with a
!>     warning (Q-P11); set_bc vel lines allowed (Q-P6);
!>   - control component (D45): control_component required with triaxial (and it must have
!>     rows), checked against the loaded component of a uniaxial preset, an error with general;
!>   - per step (§1a rule 3 as amended by D45): the control component has exactly one target
!>     (target_strain, target_load or target_stress on it, target_time, or target_duration,
!>     step 8.2c) and, with a preset, exactly one rate (strain_rate, load_rate or stress_rate on
!>     it, step 8.2c); every other component with rows has exactly one of strain_rate,
!>     target_strain, target_load, target_stress, load_rate and stress_rate (load and stress
!>     rows on a component other than the control component: triaxial load control, steps
!>     8.2b, 8.2c, dev/theory/addendum_triaxial_load_control.md §3.1), and not a target_strain
!>     in a step whose duration is not known (a load or stress target at a strain rate, step
!>     8.1, dev/theory/addendum_load_targets.md §3.2; a strain target at a load or stress
!>     rate, where a target_load or target_stress is not allowed either, step 8.2c);
!>   - load and stress rates of the control component (step 8.2c, addendum §3.1-§3.3): not
!>     with general; a load_rate ends at a target_load, a stress_rate at a target_stress (or
!>     either at a target_strain, target_time or target_duration); a rate 0 (a dwell) only
!>     with target_time or target_duration, not in step 1, and with dtime (not number_of_incrs
!>     or dstrain); dstrain needs a strain_rate on the control component; max_strain_incr
!>     (> 0, default 0.001, v2) only with a dwell;
!>   - per step (§1a rule 6, step 5.1): exactly one of number_of_incrs, dtime, dstrain; in a
!>     load- or stress-target step dtime or dstrain (the largest increment), not
!>     number_of_incrs, and dtime_min (the smallest increment, required, <= dtime; "-" or
!>     absent in the other steps) (step 8.1, addendum §3.4, user decision Q-L4);
!>   - load and stress targets (step 8.1): dtime_factor only with a load or stress target;
!>     check_necking needs a loaded component; load_tol (> 0) only with a target_load row,
!>     stress_tol (> 0) only with a target_stress row, max_bc_iter only with a load or stress
!>     target (step 8.2a, dev/theory/addendum_triaxial_load_control.md §3.6, §4, Q-T18), on
!>     any component (step 8.2b), load and stress rates included (step 8.2c); min_pert_frac
!>     (> 0) only with a load-controlled component (a load or stress row on a component other
!>     than the control component, step 8.2b, or a load or stress rate of the control
!>     component, step 8.2c; addendum §3.4, §4);
!>   - sub-incrementation (Q35, D63): subincrement_max_depth and subincrement_max_successive
!>     only with subincrementation on;
!>   - print: a variable given twice.
!> The checks of option combinations (an option used only with, or not with, another) take
!> their messages from the rules table cfg_rules_mod (step 9.9, D104), from which the manual's
!> notes and compatibility table are generated.
!> The checks that need the start-of-step values (the sign of the rate towards the target, a
!> target equal to the start value, time targets after the start of the step) and the run
!> limits are made where the steps are scheduled (step_schedule_mod). Node sets of set_bc lines
!> are checked against the mesh by bc_mod::bc_resolve.
!>
!> Public types:
!>   cfg_phase_t   material of one phase (components at the type)
!>   cfg_vel_bc_t  one velocity condition of a set_bc vel line (components at the type)
!>   cfg_mpc_t     one multi-point constraint of a set_mpc1 or set_mpc line (components at the
!>                 type)
!>   config_t      the configuration (components at the type)
!>
!> Public routines:
!>   config_read   read, validate and build the configuration of a file
!>   config_build  build the configuration from a validated store
!>   config_log_material  log the material part of the input summary (step 9.7)
!>   config_log_phase     log one phase's parameters (step 9.7)
!>   config_log_options   log the options other than the defaults and the output (step 9.7)
module config_mod

  use kinds_mod, only: rk, ik
  use error_mod, only: fatal_error
  use log_mod, only: log_warning, log_info, log_real_text
  use linalg_mod, only: ldlt_factor
  use hardening_mod, only: hardening_default_matrix
  use crystal_mod, only: stiffness_cubic, stiffness_hexagonal, stiffness_tetragonal, &
    & shear_modulus_hill, crystal_structure_from_name, crystal_num_families, crystal_family_token, &
    & crystal_family_from_token, crystal_family_size, crystal_default_families, &
    & crystal_loop_normals
  use cfg_table_mod, only: cfg_keywords, cfg_find_keyword, cfg_num_material, cfg_per_step, &
    & cfg_comp_per_step, cfg_keyword_t, cfg_active, cfg_global, cfg_phase, cfg_single, cfg_list, &
    & cfg_flag, cfg_int, cfg_real, cfg_word
  use hdf5_file_mod, only: hdf5_available
  use cfg_rules_mod
  use cfg_store_mod, only: cfg_store_t, cfg_store_read, cfg_find, cfg_error, cfg_get_real, &
    & cfg_get_int, cfg_get_word, cfg_is_dash

  implicit none

  private

  public :: cfg_phase_t, cfg_vel_bc_t, cfg_mpc_t, config_t, config_read, config_build
  public :: config_log_material, config_log_phase, config_log_options
  public :: incr_by_count, incr_by_dtime, incr_by_dstrain
  public :: target_by_strain, target_by_time, target_by_load, target_by_stress
  public :: target_by_duration, rate_by_strain, rate_by_load, rate_by_stress

  ! How the increments of a step are given (config_t%incr_kind): number_of_incrs, dtime or
  ! dstrain.
  integer, parameter :: incr_by_count = 1, incr_by_dtime = 2, incr_by_dstrain = 3
  ! Kind of target of a step (config_t%target_kind), in the order of the target rows.
  integer, parameter :: target_by_strain = 1, target_by_load = 2, target_by_stress = 3, &
    & target_by_time = 4, target_by_duration = 5
  ! Rate of the control component of a step (config_t%rate_kind; step 8.2c): strain_rate,
  ! load_rate or stress_rate (the set_bc velocities of general count as a strain rate).
  integer, parameter :: rate_by_strain = 1, rate_by_load = 2, rate_by_stress = 3

  !> Material of one phase.
  !>
  !> Components:
  !>   crystal_type  character(8)    fcc, bcc, hcp or bct
  !>   model         character(24)   constitutive model (elastic, crystal_plasticity)
  !>   c11 ... c66   real(rk)        Voigt elastic constants (0 when not used by the type) [MPa]
  !>   c_over_a      real(rk)        c/a (hcp, bct; 0 when not given) [-]
  !>   stiff_el_c    real(rk)(6, 6)  elastic stiffness, Mandel, crystal frame [MPa]
  !>   families      integer(:)      active slip families (crystal_mod canonical families), in
  !>                                 the slip_families order (empty with model elastic)
  !>   num_slip      integer         slip systems of the active families (0 with model elastic)
  !>   m             real(rk)(num_slip)  rate sensitivity per slip system, in the slip-system
  !>                                 order (family after family as in families) [-]
  !>   gammadot_0    real(rk)(num_slip)  reference slip rate per slip system [1/s]
  !>   tau_0         real(rk)(num_slip)  initial slip strength per slip system [MPa]
  !>   hardening     character(24)   hardening law (saturation_isotropic, saturation_anisotropic,
  !>                                 kocks_mecking, none; '' with model elastic)
  !>   h_0           real(rk)(:)     hardening rate coefficient per slip system (saturation laws;
  !>                                 empty otherwise) [MPa]
  !>   tau_s         real(rk)(:)     saturation strength per slip system [MPa]
  !>   n_hard        real(rk)(:)     hardening exponent n per slip system [-]
  !>   interaction   real(rk)(:, :)  interaction matrix (saturation_anisotropic; empty
  !>                                 otherwise) [-]
  !>   saturation_evolution  logical  saturation_evolution on (step 7.1, D59)
  !>   gammadot_s0   real(rk)        reference slip rate of saturation evolution [1/s]
  !>   m_prime       real(rk)        exponent of saturation evolution [-]
  !>   cyclic        logical         cyclic_hardening on (step 7.1, D59)
  !>   cyclic_a      real(rk)        threshold coefficient a, 0 < a < 1 [-]
  !>   cyclic_c      real(rk)        threshold exponent c [-]
  !>   precipitation logical         precipitate_strengthening on (step 7.1, D59)
  !>   tau_p         real(rk)        precipitate strength added to tau_0 of every system
  !>                                 (0 when off; not added to strengths from simulation.opt)
  !>                                 [MPa]
  !>   precip_regime character(8)    cutting or bowing (the smaller strength), '' when off
  !>   km_b, km_h_n, km_rho_0, km_k1, km_k20, km_gammadot_k0, km_n  real(rk)  Kocks-Mecking
  !>                                 parameters (step 7.2, D60): Burgers vector [mm], network
  !>                                 hardening coefficient [-], initial density [1/mm^2], storage
  !>                                 coefficient [1/mm], recovery coefficient [-], reference
  !>                                 slip rate of k_2 [1/s], its exponent n [-] (0 otherwise)
  !>   irradiation   logical         irradiation on (step 7.3, D61; with kocks_mecking)
  !>   irr_h_d, irr_eta, irr_loop_density, irr_loop_diameter  real(rk)  defect hardening
  !>                                 coefficient [-], annihilation efficiency [-], loop number
  !>                                 density [1/mm^3], loop diameter [mm] (0 otherwise)
  !>   irr_loop_normals_c real(rk)(3, :)  unit loop normals of the irr_loop_normals families,
  !>                                 crystal frame (size 0 otherwise)
  type :: cfg_phase_t
    character(len=8) :: crystal_type = ''
    character(len=24) :: model = ''
    real(rk) :: c11 = 0.0_rk, c12 = 0.0_rk, c13 = 0.0_rk, c33 = 0.0_rk, c44 = 0.0_rk, &
      & c66 = 0.0_rk
    real(rk) :: c_over_a = 0.0_rk
    real(rk) :: stiff_el_c(6, 6) = 0.0_rk
    integer, allocatable :: families(:)
    integer :: num_slip = 0
    real(rk), allocatable :: m(:)
    real(rk), allocatable :: gammadot_0(:)
    real(rk), allocatable :: tau_0(:)
    character(len=24) :: hardening = ''
    real(rk), allocatable :: h_0(:)
    real(rk), allocatable :: tau_s(:)
    real(rk), allocatable :: n_hard(:)
    real(rk), allocatable :: interaction(:, :)
    logical :: saturation_evolution = .false.
    real(rk) :: gammadot_s0 = 1.0_rk
    real(rk) :: m_prime = 0.0_rk
    logical :: cyclic = .false.
    real(rk) :: cyclic_a = 0.0_rk
    real(rk) :: cyclic_c = 0.0_rk
    logical :: precipitation = .false.
    real(rk) :: tau_p = 0.0_rk
    character(len=8) :: precip_regime = ''
    real(rk) :: km_b = 0.0_rk, km_h_n = 0.0_rk, km_rho_0 = 0.0_rk, km_k1 = 0.0_rk, &
      & km_k20 = 0.0_rk, km_gammadot_k0 = 0.0_rk, km_n = 0.0_rk
    logical :: irradiation = .false.
    real(rk) :: irr_h_d = 0.0_rk, irr_eta = 0.0_rk, irr_loop_density = 0.0_rk, &
      & irr_loop_diameter = 0.0_rk
    real(rk), allocatable :: irr_loop_normals_c(:, :)
  end type cfg_phase_t

  !> One velocity condition (a node set or node id, a direction and a value).
  !>
  !> Components:
  !>   nset   character(:), allocatable  node-set label or node id
  !>   comp   integer                    direction (1, 2, 3 = x, y, z)
  !>   value  real(rk)                   velocity [mm/s]
  !>   line   integer(ik)                line in the configuration file
  type :: cfg_vel_bc_t
    character(len=:), allocatable :: nset
    integer :: comp = 0
    real(rk) :: value = 0.0_rk
    integer(ik) :: line = 0_ik
  end type cfg_vel_bc_t

  !> One multi-point constraint (step 8.3a): every node s of the node set but the first (p)
  !> gets v_s = coeff v_p + offset in the direction.
  !>
  !> Components:
  !>   nset    character(:), allocatable  node-set label
  !>   comp    integer                    direction (1, 2, 3 = x, y, z)
  !>   coeff   real(rk)                   coefficient (nonzero; 1 for set_mpc1) [-]
  !>   offset  real(rk)                   offset (0 for set_mpc1) [mm/s]
  !>   line    integer(ik)                line in the configuration file
  type :: cfg_mpc_t
    character(len=:), allocatable :: nset
    integer :: comp = 0
    real(rk) :: coeff = 1.0_rk
    real(rk) :: offset = 0.0_rk
    integer(ik) :: line = 0_ik
  end type cfg_mpc_t

  !> The configuration.
  !>
  !> Components:
  !>   file_name          character(:), allocatable  configuration file
  !>   num_phases         integer                    number of phases
  !>   phases             type(cfg_phase_t)(:)       (num_phases) materials
  !>   bc_type            character(24)              boundary_conditions value
  !>   control_axis       integer                    control component (1..3, D45):
  !>                                                 control_component with triaxial; the
  !>                                                 loaded component of a uniaxial preset, or
  !>                                                 of general's target rows; 0 if none
  !>   comp_rows          logical(3)                 normal component with strain_rate,
  !>                                                 target_strain, target_load,
  !>                                                 target_stress, load_rate or stress_rate
  !>                                                 rows
  !>   comp_label         character(2)(3)            the components of the control rows as
  !>                                                 text, for messages: 11, 22, 33; periodic:
  !>                                                 the control component ij first (step 8.3b)
  !>   periodic_kinematics  character(32)            periodic: small_strain or
  !>                                                 deformation_gradient_rate (step 8.3b)
  !>   periodic_preconditioner  character(16)        periodic: lumped_diagonal or exact_diagonal
  !>   per_imposed        logical(9)                 periodic: imposed components of Fbar-dot,
  !>                                                 index 3 (i - 1) + j for ij (small strain: ij
  !>                                                 and ji together)
  !>   per_rate           real(rk)(9, num_steps)     periodic: imposed rates of Fbar-dot per step
  !>                                                 (0 where free) [1/s]
  !>   per_control        integer                    periodic: index of the control component
  !>                                                 (1..9; also control_axis 1 for the steps)
  !>   vel_bcs            type(cfg_vel_bc_t)(:)      set_bc vel conditions (general, periodic)
  !>   mpcs               type(cfg_mpc_t)(:)         multi-point constraints, in file order
  !>   num_steps          integer                    number of steps
  !>   strain_rate        real(rk)(3, num_steps)     strain rate per normal component and step
  !>                                                 (0 where not given) [1/s]
  !>   has_strain_rate    logical(3, num_steps)      strain rate given
  !>   target_strain      real(rk)(3, num_steps)     strain target (0 where not given) [-]
  !>   has_target_strain  logical(3, num_steps)      strain target given
  !>   target_load        real(rk)(3, num_steps)     load target: force on the face x1/y1/z1,
  !>                                                 that component (0 where not given) [N]
  !>   has_target_load    logical(3, num_steps)      load target given
  !>   target_stress      real(rk)(3, num_steps)     stress target: that force over the current
  !>                                                 face area (0 where not given) [MPa]
  !>   has_target_stress  logical(3, num_steps)      stress target given
  !>   load_rate          real(rk)(3, num_steps)     load rate: of the force on the face
  !>                                                 x1/y1/z1 (0 where not given) [N/s]
  !>   has_load_rate      logical(3, num_steps)      load rate given
  !>   stress_rate        real(rk)(3, num_steps)     stress rate: of that force over the current
  !>                                                 face area (0 where not given) [MPa/s]
  !>   has_stress_rate    logical(3, num_steps)      stress rate given
  !>   target_time        real(rk)(num_steps)        time target (0 where not given) [s]
  !>   has_target_time    logical(num_steps)         time target given
  !>   target_duration    real(rk)(num_steps)        length of the step from its actual start
  !>                                                 (0 where not given) [s]
  !>   has_target_duration logical(num_steps)        target_duration given
  !>   target_kind        integer(num_steps)         target of the step: target_by_strain,
  !>                                                 target_by_load, target_by_stress,
  !>                                                 target_by_time or target_by_duration
  !>   rate_kind          integer(num_steps)         rate of the control component:
  !>                                                 rate_by_strain, rate_by_load or
  !>                                                 rate_by_stress (step 8.2c)
  !>   incr_kind          integer(num_steps)         increments given by: incr_by_count,
  !>                                                 incr_by_dtime or incr_by_dstrain
  !>   num_incrs          integer(num_steps)         increments per step (number_of_incrs; 0
  !>                                                 where not given)
  !>   dtime              real(rk)(num_steps)        time increment (dtime; 0 where not given) [s]
  !>   dstrain            real(rk)(num_steps)        strain increment on the loading component
  !>                                                 (dstrain; 0 where not given) [-]
  !>   dtime_min          real(rk)(num_steps)        smallest time increment of a load- or
  !>                                                 stress-target step (0 elsewhere) [s]
  !>   dtime_factor       real(rk)                   factor on the increment aimed at a load or
  !>                                                 stress target (v2, default 1.001) [-]
  !>   load_tol           real(rk)                   absolute tolerance of load rows: landing
  !>                                                 on a load target (step 8.2a) [N]
  !>   stress_tol         real(rk)                   absolute tolerance of stress rows [MPa]
  !>   max_bc_iter        integer                    landing attempts of an increment, and
  !>                                                 iterations of the boundary-velocity
  !>                                                 iteration (step 8.2b)
  !>   min_pert_frac      real(rk)                   lower bound of the velocity perturbations,
  !>                                                 fraction of the control face velocity
  !>                                                 (step 8.2b, v2) or of the elastic velocity
  !>                                                 estimate (step 8.2c, v2-fixed C60) [-]
  !>   max_strain_incr    real(rk)                   largest strain increment of the control
  !>                                                 component in a dwell (step 8.2c, v2) [-]
  !>   check_necking      logical                    stop when the control-face load falls
  !>   print_results      logical(num_steps)         write the results of the step
  !>   max_incr           integer                    run limit: increments
  !>   max_total_time     real(rk)                   run limit: time [s]
  !>   max_strain         real(rk)                   run limit: largest normal engineering
  !>                                                 strain (when has_max_strain) [-]
  !>   has_max_strain     logical                    max_strain given (default: no limit)
  !>   max_eqstrain       real(rk)                   run limit: equivalent strain (when
  !>                                                 has_max_eqstrain) [-]
  !>   has_max_eqstrain   logical                    max_eqstrain given (default: no limit)
  !>   cg_tol             real(rk)                   CG relative tolerance [-]
  !>   cg_max_iters       integer                    CG iteration limit
  !>   cg_preconditioner  character(16)              jacobi, block_jacobi or none (step 9.2)
  !>   global_max_iters       integer                    global iteration limit
  !>   global_tol      real(rk)                   relative velocity-change tolerance [-]
  !>   global_tol_abs         real(rk)                   absolute velocity-change tolerance [mm/s]
  !>   global_tol_slow       real(rk)                   tolerance accepted when NR or SA is slow [-]
  !>   global_nr_switch  real(rk)                   SA-to-NR switch tolerance (0: never) [-]
  !>   global_nr_slow_ratio        real(rk)                   NR slow-convergence fraction [-]
  !>   midnode_linearization  logical                reset tet10 mid-nodes (D26)
  !>   initialization     character(32)              start of crystal-plasticity runs
  !>                                                 (elastic_step or isotropic_viscoplastic;
  !>                                                 D49, D54, D78)
  !>   point_tol             real(rk)                   point Newton tolerance [1/s]
  !>   point_max_iters_newton  integer                  point Newton iteration limit
  !>   point_max_iters_state integer                    point state-loop iteration limit
  !>   point_max_iters_hardening  integer               point strength-update iteration limit
  !>   point_tol_state       real(rk)                   point state-loop tolerance, relative to
  !>                                                    tau_0 (crystal plasticity) or to the
  !>                                                    stress norm (elastic) [-] (D55)
  !>   point_tol_hardening   real(rk)                   point strength-update tolerance,
  !>                                                    relative to tau_0 of each system [-] (D55)
  !>   subincrementation     logical                    repeat a rejected increment as
  !>                                                    halves (Q35, D63)
  !>   subincrement_max_depth  integer                  subdivisions of one increment (L1)
  !>   subincrement_max_successive  integer             successive subdivided increments (L2)
  !>   read_ori_from_file, read_phase_from_file, read_from_opt_file  logical  input files
  !>   print_vars         character(32)(:)           output variables, in input order (restart
  !>                                                    not among them: print_restart)
  !>   print_restart      logical                    write the restart files at every step end
  !>   output_ascii       logical                    write the .sim text files (output_format
  !>                                                    ascii, the default; step 9.1, D93)
  !>   output_hdf5        logical                    write the HDF5 and XDMF files (output_format
  !>                                                    hdf5)
  !>                                                    (print restart; step 8.4)
  !>   restart            logical                    continue simulation.sim from its restart
  !>                                                    files (restart on; step 8.4)
  !>   fingerprint        character(512)(:)          every record but the per-step rows,
  !>                                                    number_of_steps and restart, as text
  !>                                                    "<keyword>[ phase <p>][ <comp>]:
  !>                                                    <values>", in file order (the restart
  !>                                                    check, addendum_restart.md §3 R1)
  !>   fingerprint_strict logical(:)                 the record must be unchanged on a restart
  !>                                                    (material, boundary conditions and
  !>                                                    constraints, input files); others may
  !>                                                    change, with a warning
  !>   store              type(cfg_store_t)          the validated records (set by config_read;
  !>                                                    the input summary, step 9.7)
  type :: config_t
    character(len=:), allocatable :: file_name
    integer :: num_phases = 0
    type(cfg_phase_t), allocatable :: phases(:)
    character(len=24) :: bc_type = ''
    integer :: control_axis = 0
    logical :: comp_rows(3) = .false.
    character(len=2) :: comp_label(3) = ['11', '22', '33']
    character(len=32) :: periodic_kinematics = ''
    character(len=16) :: periodic_preconditioner = ''
    logical :: per_imposed(9) = .false.
    real(rk), allocatable :: per_rate(:, :)
    integer :: per_control = 0
    type(cfg_vel_bc_t), allocatable :: vel_bcs(:)
    type(cfg_mpc_t), allocatable :: mpcs(:)
    integer :: num_steps = 0
    real(rk), allocatable :: strain_rate(:, :)
    logical, allocatable :: has_strain_rate(:, :)
    real(rk), allocatable :: target_strain(:, :)
    logical, allocatable :: has_target_strain(:, :)
    real(rk), allocatable :: target_load(:, :)
    logical, allocatable :: has_target_load(:, :)
    real(rk), allocatable :: target_stress(:, :)
    logical, allocatable :: has_target_stress(:, :)
    real(rk), allocatable :: load_rate(:, :)
    logical, allocatable :: has_load_rate(:, :)
    real(rk), allocatable :: stress_rate(:, :)
    logical, allocatable :: has_stress_rate(:, :)
    real(rk), allocatable :: target_time(:)
    logical, allocatable :: has_target_time(:)
    real(rk), allocatable :: target_duration(:)
    logical, allocatable :: has_target_duration(:)
    integer, allocatable :: target_kind(:)
    integer, allocatable :: rate_kind(:)
    integer, allocatable :: incr_kind(:)
    integer, allocatable :: num_incrs(:)
    real(rk), allocatable :: dtime(:)
    real(rk), allocatable :: dstrain(:)
    real(rk), allocatable :: dtime_min(:)
    real(rk) :: dtime_factor = 0.0_rk
    real(rk) :: load_tol = 0.0_rk
    real(rk) :: stress_tol = 0.0_rk
    integer :: max_bc_iter = 0
    real(rk) :: min_pert_frac = 0.0_rk
    real(rk) :: max_strain_incr = 0.0_rk
    logical :: check_necking = .false.
    logical, allocatable :: print_results(:)
    integer :: max_incr = 0
    real(rk) :: max_total_time = 0.0_rk
    real(rk) :: max_strain = 0.0_rk
    logical :: has_max_strain = .false.
    real(rk) :: max_eqstrain = 0.0_rk
    logical :: has_max_eqstrain = .false.
    real(rk) :: cg_tol = 0.0_rk
    integer :: cg_max_iters = 0
    character(len=16) :: cg_preconditioner = ''
    integer :: global_max_iters = 0
    real(rk) :: global_tol = 0.0_rk
    real(rk) :: global_tol_abs = 0.0_rk
    real(rk) :: global_tol_slow = 0.0_rk
    real(rk) :: global_nr_switch = 0.0_rk
    real(rk) :: global_nr_slow_ratio = 0.0_rk
    logical :: midnode_linearization = .false.
    character(len=32) :: initialization = ''
    real(rk) :: point_tol = 0.0_rk
    integer :: point_max_iters_newton = 0
    integer :: point_max_iters_state = 0
    integer :: point_max_iters_hardening = 0
    real(rk) :: point_tol_state = 0.0_rk
    real(rk) :: point_tol_hardening = 0.0_rk
    logical :: subincrementation = .false.
    integer :: subincrement_max_depth = 0
    integer :: subincrement_max_successive = 0
    logical :: read_ori_from_file = .false.
    logical :: read_phase_from_file = .false.
    logical :: read_from_opt_file = .false.
    character(len=32), allocatable :: print_vars(:)
    logical :: print_restart = .false.
    logical :: output_ascii = .true.
    logical :: output_hdf5 = .false.
    logical :: restart = .false.
    character(len=512), allocatable :: fingerprint(:)
    logical, allocatable :: fingerprint_strict(:)
    type(cfg_store_t) :: store
  end type config_t

  ! Axis letters of set_bc directions.
  character(len=1), parameter :: axis_letter(3) = ['x', 'y', 'z']

contains

  !> Read a configuration file: validate every line (cfg_store_read), then build the typed
  !> configuration with the cross-keyword checks (config_build).
  !>
  !> Arguments:
  !>   config     (out) type(config_t)  the configuration
  !>   file_name  (in)  character       configuration file
  !> Locals:
  !>   store      type(cfg_store_t)     validated records
  subroutine config_read(config, file_name)

    type(config_t), intent(out) :: config
    character(len=*), intent(in) :: file_name

    type(cfg_store_t) :: store

    call cfg_store_read(store, file_name)
    call config_build(config, store)
    config%store = store

  end subroutine config_read

  !> Build the configuration from a validated store, with the checks listed in the module
  !> header.
  !>
  !> Arguments:
  !>   config  (out) type(config_t)     the configuration
  !>   store   (in)  type(cfg_store_t)  validated records
  !> Locals:
  !>   i_rec   integer                  record index
  subroutine config_build(config, store)

    type(config_t), intent(out) :: config
    type(cfg_store_t), intent(in) :: store

    integer :: i_rec

    config%file_name = store%file_name

    ! Restart (step 8.4): known first (step 1 of a restarted run follows the restart point).
    config%restart = word_or_default(store, 'restart') == 'on'

    ! Material, boundary conditions and loading, solver and output.
    call build_phases(config, store)
    call build_loading(config, store)
    config%cg_tol = real_or_default(store, 'cg_tol')
    config%cg_max_iters = nint(real_or_default(store, 'cg_max_iters'))
    ! The preconditioner (step 9.2, addendum_block_jacobi.md §2, §3, Q-B1, Q-B2): default
    ! block_jacobi; with constraints (MPCs, periodicity) only jacobi, the diagonal chosen by
    ! periodic_preconditioner: the default there is jacobi, another value given is an error.
    config%cg_preconditioner = word_or_default(store, 'cg_preconditioner')
    i_rec = cfg_find(store, 'cg_preconditioner', 0, 0)
    if (i_rec == 0 .and. (config%bc_type == 'periodic' .or. size(config%mpcs) > 0)) &
      & config%cg_preconditioner = 'jacobi'
    if (i_rec > 0 .and. config%cg_preconditioner /= 'jacobi') then
      if (config%bc_type == 'periodic') call cfg_error(store, i_rec, &
        & cfg_rule_message(rule_preconditioner_periodic, 'cg_preconditioner ' // &
        & trim(config%cg_preconditioner)))
      if (size(config%mpcs) > 0) call cfg_error(store, i_rec, &
        & cfg_rule_message(rule_preconditioner_mpc, 'cg_preconditioner ' // &
        & trim(config%cg_preconditioner)))
    end if
    config%global_max_iters = nint(real_or_default(store, 'global_max_iters'))
    config%global_tol = real_or_default(store, 'global_tol')
    config%global_tol_abs = real_or_default(store, 'global_tol_abs')
    config%global_tol_slow = real_or_default(store, 'global_tol_slow')
    config%global_nr_switch = real_or_default(store, 'global_nr_switch')
    config%global_nr_slow_ratio = real_or_default(store, 'global_nr_slow_ratio')
    config%midnode_linearization = word_or_default(store, 'midnode_linearization') == 'on'
    config%initialization = word_or_default(store, 'initialization')
    config%point_tol = real_or_default(store, 'point_tol')
    config%point_max_iters_newton = nint(real_or_default(store, 'point_max_iters_newton'))
    config%point_max_iters_state = nint(real_or_default(store, 'point_max_iters_state'))
    config%point_max_iters_hardening = nint(real_or_default(store, 'point_max_iters_hardening'))
    config%point_tol_state = real_or_default(store, 'point_tol_state')
    config%point_tol_hardening = real_or_default(store, 'point_tol_hardening')
    ! Sub-incrementation (Q35, D63, addendum_subincrementation.md §4): its limits only with the
    ! switch on.
    config%subincrementation = word_or_default(store, 'subincrementation') == 'on'
    config%subincrement_max_depth = nint(real_or_default(store, 'subincrement_max_depth'))
    config%subincrement_max_successive = nint(real_or_default(store, &
      & 'subincrement_max_successive'))
    do i_rec = 1, store%num_records
      associate (name => cfg_keywords(store%records(i_rec)%key)%name)
        if ((name == 'subincrement_max_depth' .or. name == 'subincrement_max_successive') &
          & .and. .not. config%subincrementation) call cfg_error(store, i_rec, &
          & cfg_rule_message(rule_subincrementation, trim(name)))
      end associate
    end do
    config%read_ori_from_file = cfg_find(store, 'read_ori_from_file', 0, 0) > 0
    config%read_phase_from_file = cfg_find(store, 'read_phase_from_file', 0, 0) > 0
    config%read_from_opt_file = cfg_find(store, 'read_from_opt_file', 0, 0) > 0

    ! Output variables in input order; a variable given twice is an error.
    allocate (config%print_vars(0))
    do i_rec = 1, store%num_records
      if (cfg_keywords(store%records(i_rec)%key)%name /= 'print') cycle
      if (any(config%print_vars == cfg_get_word(store, i_rec, 1))) call cfg_error(store, i_rec, &
        & 'print ' // cfg_get_word(store, i_rec, 1) // ' given twice')
      ! Periodic conditions (step 8.3b, Q-P11): no face forces file (the face loads are not
      ! meaningful); the macroscopic quantities are written instead.
      if (config%bc_type == 'periodic' .and. cfg_get_word(store, i_rec, 1) == 'forces') then
        call log_warning(store%file_name // ':' // int_text(int(store%records(i_rec)%line)) // &
          & ': ' // cfg_rule_message(rule_print_forces_periodic, 'print forces') // '.', 1)
        cycle
      end if
      ! Restart files (step 8.4): not a result variable of the writer.
      if (cfg_get_word(store, i_rec, 1) == 'restart') then
        config%print_restart = .true.
        cycle
      end if
      config%print_vars = [character(len=32) :: config%print_vars, cfg_get_word(store, i_rec, 1)]
      ! The Kocks-Mecking outputs need a Kocks-Mecking phase (step 7.2, D60).
      if (any(cfg_get_word(store, i_rec, 1) == [character(len=10) :: 'km_rho', 'km_k2', &
        & 'km_rhorate']) .and. .not. any(config%phases(:)%hardening == 'kocks_mecking')) &
        & call cfg_error(store, i_rec, cfg_rule_message(rule_print_kocks_mecking, 'print ' // &
        & cfg_get_word(store, i_rec, 1)))
      ! The irradiation outputs need a phase with irradiation on (step 7.3, D61).
      if (any(cfg_get_word(store, i_rec, 1) == [character(len=10) :: 'dmg', 'dmgrate']) .and. &
        & .not. any(config%phases(:)%irradiation)) call cfg_error(store, i_rec, &
        & cfg_rule_message(rule_print_irradiation, 'print ' // cfg_get_word(store, i_rec, 1)))
    end do

    ! Output formats (step 9.1, D93): ascii by default; a format given twice is an error; hdf5
    ! needs a build with the HDF5 output.
    i_rec = cfg_find(store, 'output_format', 0, 0)
    if (i_rec > 0) then
      associate (words => store%records(i_rec)%values)
        if (size(words) > 2) call cfg_error(store, i_rec, 'output_format: give ascii, hdf5, ' // &
          & 'or both once each')
        if (size(words) == 2) then
          if (words(1) == words(2)) call cfg_error(store, i_rec, 'output_format: give ascii, ' // &
            & 'hdf5, or both once each')
        end if
        config%output_ascii = any(words == 'ascii')
        config%output_hdf5 = any(words == 'hdf5')
      end associate
      if (config%output_hdf5 .and. .not. hdf5_available) call cfg_error(store, i_rec, &
        & 'output_format hdf5: this FEPX build has no HDF5 output (configure with ' // &
        & '-DFEPX_ENABLE_HDF5=ON)')
    end if

    ! Restart (step 8.4): the fingerprint of the records.
    call build_fingerprint(config, store)

  end subroutine config_build

  !> The restart fingerprint (step 8.4, addendum_restart.md §3 R1, §5): every record of the
  !> file but the per-step rows (the new steps of a restarted run), number_of_steps and restart,
  !> as one line of text each, in file order; strict (unchanged on a restart) for the material
  !> keywords, the boundary-condition and constraint definitions, the control component, the
  !> periodic kinematics, the input-file switches and the output formats (step 9.1). The output
  !> formats are one line built from the formats written, "output_format: ascii", "hdf5" or
  !> "ascii hdf5", whether or not the keyword is given and in whatever order its values are: the
  !> default and an explicit output_format ascii are the same output, and are not a change.
  !>
  !> Arguments:
  !>   config  (inout) type(config_t)     the configuration (fingerprint set)
  !>   store   (in)    type(cfg_store_t)  validated records
  !> Locals:
  !>   i_rec, j  integer                  record and value indices
  !>   n         integer                  lines so far
  !>   line      character(:), allocatable  the line of a record
  !>   strict_names  character(*)         the strict keywords outside the material group
  !>   formats   character(:), allocatable  the output formats written, in a fixed order
  subroutine build_fingerprint(config, store)

    type(config_t), intent(inout) :: config
    type(cfg_store_t), intent(in) :: store

    integer :: i_rec, j, n
    character(len=:), allocatable :: line, formats
    character(len=*), parameter :: strict_names = ' boundary_conditions set_bc set_mpc1 ' // &
      & 'set_mpc control_component periodic_kinematics read_ori_from_file ' // &
      & 'read_phase_from_file read_from_opt_file '

    ! One line per record and one for the output formats (allocated for every record, the first n
    ! used).
    allocate (config%fingerprint(store%num_records + 1), &
      & config%fingerprint_strict(store%num_records + 1))
    n = 0
    do i_rec = 1, store%num_records
      associate (rec => store%records(i_rec), key => cfg_keywords(store%records(i_rec)%key))
        if (key%form == cfg_per_step .or. key%form == cfg_comp_per_step) cycle
        if (key%name == 'number_of_steps' .or. key%name == 'restart') cycle
        if (key%name == 'output_format') cycle
        ! "<keyword>[ phase <p>][ <comp>]: <values>"
        line = trim(key%name)
        if (rec%phase > 0) line = line // ' phase ' // int_text(rec%phase)
        if (rec%comp > 0) line = line // ' ' // int_text(rec%comp)
        line = line // ':'
        do j = 1, size(rec%values)
          line = line // ' ' // trim(rec%values(j))
        end do
        n = n + 1
        config%fingerprint(n) = line
        config%fingerprint_strict(n) = rec%key <= cfg_num_material .or. &
          & index(strict_names, ' ' // trim(key%name) // ' ') > 0
      end associate
    end do

    ! The output formats written (strict), the same line for the default and for the keyword.
    formats = ''
    if (config%output_ascii) formats = formats // ' ascii'
    if (config%output_hdf5) formats = formats // ' hdf5'
    n = n + 1
    config%fingerprint(n) = 'output_format:' // formats
    config%fingerprint_strict(n) = .true.
    config%fingerprint = config%fingerprint(1:n)
    config%fingerprint_strict = config%fingerprint_strict(1:n)

  end subroutine build_fingerprint

  !> The phases: number_of_phases, one block per phase with its crystal type, model and
  !> elastic constants (checked against the crystal type and for positive definiteness).
  !>
  !> Arguments:
  !>   config  (inout) type(config_t)     the configuration
  !>   store   (in)    type(cfg_store_t)  validated records
  !> Locals:
  !>   i_rec, i_phase  integer            record and phase indices
  !>   num             integer            number_of_phases record
  !>   const           character(3)(6)    names of the elastic constants
  !>   needed          logical(6)         constants of the crystal type
  !>   values          real(rk)(6)        constants given (0 where absent)
  !>   i_const         integer            constant index
  !>   status          integer            ldlt_factor status
  !>   fact            real(rk)(6, 6)     the stiffness, then its LDL^T factors (unused)
  !>   phase_rec       integer            the record of "phase i_phase"
  subroutine build_phases(config, store)

    type(config_t), intent(inout) :: config
    type(cfg_store_t), intent(in) :: store

    integer :: i_rec, i_phase, num, i_const, status, phase_rec
    character(len=3), parameter :: const(6) = ['c11', 'c12', 'c13', 'c33', 'c44', 'c66']
    logical :: needed(6)
    real(rk) :: values(6), fact(6, 6)

    ! Number of phases (required); every "phase i" within 1..number_of_phases.
    num = require(store, 'number_of_phases')
    config%num_phases = cfg_get_int(store, num, 1)
    allocate (config%phases(config%num_phases))
    do i_rec = 1, store%num_records
      if (cfg_keywords(store%records(i_rec)%key)%name /= 'phase') cycle
      if (cfg_get_int(store, i_rec, 1) > config%num_phases) call cfg_error(store, i_rec, &
        & 'phase ' // cfg_get_word(store, i_rec, 1) // ' is outside 1..' // &
        & int_text(config%num_phases) // ' (number_of_phases)')
    end do

    do i_phase = 1, config%num_phases
      ! The block must exist.
      phase_rec = 0
      do i_rec = 1, store%num_records
        if (cfg_keywords(store%records(i_rec)%key)%name == 'phase') then
          if (cfg_get_int(store, i_rec, 1) == i_phase) phase_rec = i_rec
        end if
      end do
      if (phase_rec == 0) call fatal_error('cfg', store%file_name // ': phase ' // &
        & int_text(i_phase) // ' is not defined (number_of_phases ' // &
        & int_text(config%num_phases) // ')')

      associate (ph => config%phases(i_phase))
        ! Crystal type and model (required).
        ph%crystal_type = cfg_get_word(store, require_phase(store, 'crystal_type', i_phase, &
          & phase_rec), 1)
        ph%model = cfg_get_word(store, require_phase(store, 'model', i_phase, phase_rec), 1)

        ! The constants of the crystal type are required, the others are errors.
        select case (ph%crystal_type)
        case ('fcc', 'bcc')
          needed = [.true., .true., .false., .false., .true., .false.]
        case ('hcp')
          needed = [.true., .true., .true., .true., .true., .false.]
        case default
          needed = .true.
        end select
        do i_const = 1, 6
          i_rec = cfg_find(store, const(i_const), i_phase, 0)
          if (needed(i_const) .and. i_rec == 0) call fatal_error('cfg', store%file_name // &
            & ': phase ' // int_text(i_phase) // ' (' // trim(ph%crystal_type) // &
            & ', line ' // line_text(store, phase_rec) // '): ' // const(i_const) // &
            & ' is required')
          if (.not. needed(i_const) .and. i_rec > 0) call cfg_error(store, i_rec, &
            & const(i_const) // ' is not a constant of crystal_type ' // trim(ph%crystal_type))
          values(i_const) = 0.0_rk
          if (i_rec > 0) values(i_const) = cfg_get_real(store, i_rec, 1)
        end do
        ph%c11 = values(1)
        ph%c12 = values(2)
        ph%c13 = values(3)
        ph%c33 = values(4)
        ph%c44 = values(5)
        ph%c66 = values(6)

        ! c/a: hcp and bct only.
        i_rec = cfg_find(store, 'c_over_a', i_phase, 0)
        if (i_rec > 0) then
          if (ph%crystal_type == 'fcc' .or. ph%crystal_type == 'bcc') call cfg_error(store, &
            & i_rec, cfg_rule_message(rule_c_over_a, 'c_over_a'))
          ph%c_over_a = cfg_get_real(store, i_rec, 1)
        end if

        ! Stiffness, positive definite (all LDL^T pivots positive).
        select case (ph%crystal_type)
        case ('fcc', 'bcc')
          ph%stiff_el_c = stiffness_cubic(ph%c11, ph%c12, ph%c44)
        case ('hcp')
          ph%stiff_el_c = stiffness_hexagonal(ph%c11, ph%c12, ph%c13, ph%c33, ph%c44)
        case default
          ph%stiff_el_c = stiffness_tetragonal(ph%c11, ph%c12, ph%c13, ph%c33, ph%c44, ph%c66)
        end select
        fact = ph%stiff_el_c
        call ldlt_factor(fact, status)
        if (status /= 0) call fatal_error('cfg', store%file_name // ': phase ' // &
          & int_text(i_phase) // ' (line ' // line_text(store, phase_rec) // &
          & '): the elastic constants are not positive definite')

        ! Slip families and kinetics, then hardening (crystal_plasticity).
        call build_slip(ph, store, i_phase, phase_rec)
        call build_hardening(ph, store, i_phase, phase_rec)
        call build_modifiers(ph, store, i_phase, phase_rec)
        call build_kocks_mecking(ph, store, i_phase, phase_rec)
        call build_irradiation(ph, store, i_phase, phase_rec)
      end associate
    end do

  end subroutine build_phases

  !> Slip families, slip kinetics and strengths of a phase (dev/cfg_keywords.md §2.1, ROADMAP
  !> D48, D49): slip_families, m, gammadot_0 and tau_0 are errors with model elastic; with
  !> crystal_plasticity, the active families (slip_families, matched exactly, or the crystal
  !> type's default), and m, gammadot_0 and tau_0 (required) expanded to one value per slip
  !> system with the counting rule. hcp and bct need c_over_a.
  !>
  !> Arguments:
  !>   ph         (inout) type(cfg_phase_t)  the phase (crystal type and model read)
  !>   store      (in)    type(cfg_store_t)  validated records
  !>   i_phase    (in)    integer            phase
  !>   phase_rec  (in)    integer            the record of "phase i_phase"
  !> Locals:
  !>   structure  integer                    crystal_mod structure code
  !>   names      character(13)(4)           the four slip keywords
  !>   recs       integer(4)                 their records (0: absent)
  !>   i, j       integer                    keyword or value index, earlier value index
  !>   family     integer                    family of a token
  !>   nvals      integer                    number of slip_families values
  !>   word       character(:), allocatable  a slip_families value
  !>   where      character(:), allocatable  "<file>: phase <i> (<type>, line <l>)" for messages
  subroutine build_slip(ph, store, i_phase, phase_rec)

    type(cfg_phase_t), intent(inout) :: ph
    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_phase, phase_rec

    character(len=13), parameter :: names(4) = [character(len=13) :: 'slip_families', 'm', &
      & 'gammadot_0', 'tau_0']
    integer :: structure, recs(4), i, j, family, nvals
    character(len=:), allocatable :: word, where

    structure = crystal_structure_from_name(ph%crystal_type)
    recs = [(cfg_find(store, trim(names(i)), i_phase, 0), i = 1, 4)]
    where = store%file_name // ': phase ' // int_text(i_phase) // ' (' // &
      & trim(ph%crystal_type) // ', line ' // line_text(store, phase_rec) // ')'

    ! Elastic: no slip; the slip keywords are errors.
    if (ph%model /= 'crystal_plasticity') then
      do i = 1, 4
        if (recs(i) > 0) call cfg_error(store, recs(i), cfg_rule_message(rule_slip_model, &
          & trim(names(i)), detail='phase ' // int_text(i_phase) // ' is ' // trim(ph%model)))
      end do
      allocate (ph%families(0), ph%m(0), ph%gammadot_0(0), ph%tau_0(0))
      ph%num_slip = 0
      return
    end if

    ! hcp and bct slip systems depend on c/a.
    if ((ph%crystal_type == 'hcp' .or. ph%crystal_type == 'bct') .and. &
      & cfg_find(store, 'c_over_a', i_phase, 0) == 0) call fatal_error('cfg', where // &
      & ': c_over_a is required with model crystal_plasticity')

    ! Active families: the default, "all" (alone), or tokens matched exactly, each once.
    if (recs(1) == 0) then
      allocate (ph%families(size(crystal_default_families(structure))))
      ph%families(:) = crystal_default_families(structure)
    else
      nvals = size(store%records(recs(1))%values)
      allocate (ph%families(nvals))
      do i = 1, nvals
        word = cfg_get_word(store, recs(1), i)
        if (word == 'all') then
          if (nvals > 1) call cfg_error(store, recs(1), 'slip_families: all cannot be ' // &
            & 'combined with other families')
          deallocate (ph%families)
          allocate (ph%families(crystal_num_families(structure)))
          ph%families(:) = [(j, j = 1, crystal_num_families(structure))]
          exit
        end if
        family = crystal_family_from_token(structure, word)
        if (family == 0) call cfg_error(store, recs(1), 'slip_families: "' // word // &
          & '" is not a slip family of crystal_type ' // trim(ph%crystal_type) // &
          & ' (expected: ' // family_list(structure) // ', or all)')
        do j = 1, i - 1
          if (ph%families(j) == family) call cfg_error(store, recs(1), 'slip_families: ' // &
            & word // ' given twice')
        end do
        ph%families(i) = family
      end do
    end if
    ph%num_slip = sum([(crystal_family_size(structure, ph%families(i)), &
      & i = 1, size(ph%families))])

    ! m, gammadot_0 and tau_0: required, one value per active family or per slip system.
    do i = 2, 4
      if (recs(i) == 0) call fatal_error('cfg', where // ': ' // trim(names(i)) // &
        & ' is required with model crystal_plasticity')
    end do
    allocate (ph%m(ph%num_slip), ph%gammadot_0(ph%num_slip), ph%tau_0(ph%num_slip))
    ph%m(:) = per_system_values(store, recs(2), 'm', structure, ph%families, ph%num_slip)
    ph%gammadot_0(:) = per_system_values(store, recs(3), 'gammadot_0', structure, &
      & ph%families, ph%num_slip)
    ph%tau_0(:) = per_system_values(store, recs(4), 'tau_0', structure, ph%families, &
      & ph%num_slip)

  end subroutine build_slip

  !> Hardening law and parameters of a phase (step 6.4, ROADMAP D51, addendum_hardening.md):
  !> hardening, h_0, tau_s, n and interaction_matrix are errors with model elastic and with
  !> hardening none. With crystal_plasticity the law is hardening (default saturation_isotropic):
  !>   - saturation_isotropic: h_0, tau_s and n required, one value each (one increment for all
  !>     systems, set by slip system 1, §2); interaction_matrix is an error; tau_s > tau_0 of
  !>     system 1;
  !>   - saturation_anisotropic: h_0, tau_s and n required, one value or one per active family
  !>     or per system (§1); tau_s > tau_0 for every system; interaction_matrix optional: exactly
  !>     num_slip lines of num_slip values >= 0, rows in the slip-system order (default: 1 on the
  !>     diagonal, 1.4 off it, §3);
  !>   - kocks_mecking (step 7.2, D60): h_0, tau_s, n and interaction_matrix are errors (its
  !>     parameters: build_kocks_mecking).
  !> The parameters are stored per system (one value repeated).
  !>
  !> Arguments:
  !>   ph         (inout) type(cfg_phase_t)  the phase (model, families and tau_0 read)
  !>   store      (in)    type(cfg_store_t)  validated records
  !>   i_phase    (in)    integer            phase
  !>   phase_rec  (in)    integer            the record of "phase i_phase"
  !> Locals:
  !>   names      character(18)(5)           the hardening keywords
  !>   recs       integer(5)                 their records (0: absent; the first matrix row)
  !>   i, a       integer                    keyword index, slip system
  !>   structure  integer                    crystal_mod structure code
  !>   where      character(:), allocatable  "<file>: phase <i> (<type>, line <l>)" for messages
  subroutine build_hardening(ph, store, i_phase, phase_rec)

    type(cfg_phase_t), intent(inout) :: ph
    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_phase, phase_rec

    character(len=18), parameter :: names(5) = [character(len=18) :: 'hardening', 'h_0', &
      & 'tau_s', 'n', 'interaction_matrix']
    integer :: recs(5), i, a, structure
    character(len=:), allocatable :: where

    recs = [(cfg_find(store, trim(names(i)), i_phase, 0), i = 1, 5)]
    where = store%file_name // ': phase ' // int_text(i_phase) // ' (' // &
      & trim(ph%crystal_type) // ', line ' // line_text(store, phase_rec) // ')'
    allocate (ph%h_0(0), ph%tau_s(0), ph%n_hard(0), ph%interaction(0, 0))

    ! Elastic: no hardening keywords.
    if (ph%model /= 'crystal_plasticity') then
      do i = 1, 5
        if (recs(i) > 0) call cfg_error(store, recs(i), cfg_rule_message(rule_hardening_model, &
          & trim(names(i)), detail='phase ' // int_text(i_phase) // ' is ' // trim(ph%model)))
      end do
      ph%hardening = ''
      return
    end if

    ! The law (default saturation_isotropic).
    ph%hardening = 'saturation_isotropic'
    if (recs(1) > 0) ph%hardening = cfg_get_word(store, recs(1), 1)

    ! none: constant strengths, no parameters.
    if (ph%hardening == 'none') then
      do i = 2, 5
        if (recs(i) > 0) call cfg_error(store, recs(i), cfg_rule_message(rule_saturation_law, &
          & trim(names(i)), detail='phase ' // int_text(i_phase) // ' has hardening none'))
      end do
      return
    end if

    ! Kocks-Mecking: none of the saturation-law parameters (its own: build_kocks_mecking).
    if (ph%hardening == 'kocks_mecking') then
      do i = 2, 5
        if (recs(i) > 0) call cfg_error(store, recs(i), cfg_rule_message(rule_saturation_law, &
          & trim(names(i)), detail='phase ' // int_text(i_phase) // ' has hardening kocks_mecking'))
      end do
      return
    end if

    ! Both saturation laws: h_0, tau_s and n required.
    do i = 2, 4
      if (recs(i) == 0) call fatal_error('cfg', where // ': ' // trim(names(i)) // &
        & ' is required with hardening ' // trim(ph%hardening))
    end do
    deallocate (ph%h_0, ph%tau_s, ph%n_hard, ph%interaction)
    allocate (ph%h_0(ph%num_slip), ph%tau_s(ph%num_slip), ph%n_hard(ph%num_slip))
    structure = crystal_structure_from_name(ph%crystal_type)

    if (ph%hardening == 'saturation_isotropic') then
      ! One value each (one increment for all systems, set by system 1); no matrix.
      do i = 2, 4
        if (size(store%records(recs(i))%values) /= 1) call cfg_error(store, recs(i), &
          & trim(names(i)) // ': hardening saturation_isotropic takes one value (one ' // &
          & 'hardening increment for all systems, set by slip system 1); per family or per ' // &
          & 'system: saturation_anisotropic')
      end do
      if (recs(5) > 0) call cfg_error(store, recs(5), &
        & cfg_rule_message(rule_interaction_matrix, 'interaction_matrix'))
      ph%h_0 = cfg_get_real(store, recs(2), 1)
      ph%tau_s = cfg_get_real(store, recs(3), 1)
      ph%n_hard = cfg_get_real(store, recs(4), 1)
      if (.not. (ph%tau_s(1) > ph%tau_0(1))) call cfg_error(store, recs(3), 'tau_s must be ' // &
        & 'greater than tau_0 of slip system 1 (' // trim(real_text(ph%tau_0(1))) // ')')
      allocate (ph%interaction(0, 0))
      return
    end if

    ! saturation_anisotropic: one value, per family or per system; tau_s > tau_0 everywhere.
    ph%h_0(:) = per_system_values(store, recs(2), 'h_0', structure, ph%families, ph%num_slip, &
      & one_for_all=.true.)
    ph%tau_s(:) = per_system_values(store, recs(3), 'tau_s', structure, ph%families, &
      & ph%num_slip, one_for_all=.true.)
    ph%n_hard(:) = per_system_values(store, recs(4), 'n', structure, ph%families, &
      & ph%num_slip, one_for_all=.true.)
    do a = 1, ph%num_slip
      if (.not. (ph%tau_s(a) > ph%tau_0(a))) call cfg_error(store, recs(3), 'tau_s must be ' // &
        & 'greater than tau_0 for every slip system (system ' // int_text(a) // ': tau_s ' // &
        & trim(real_text(ph%tau_s(a))) // ', tau_0 ' // trim(real_text(ph%tau_0(a))) // ')')
    end do

    ! The interaction matrix: the default, or the user's N rows of N values >= 0.
    allocate (ph%interaction(ph%num_slip, ph%num_slip))
    if (recs(5) == 0) then
      ph%interaction(:, :) = hardening_default_matrix(ph%num_slip)
    else
      call read_interaction_matrix(store, i_phase, ph%num_slip, ph%interaction)
    end if

  end subroutine build_hardening

  !> Hardening modifiers of a phase (step 7.1, ROADMAP D59, addendum_hardening_modifiers.md
  !> §3-5): saturation_evolution, cyclic_hardening and precipitate_strengthening (on|off,
  !> default off), each only with a saturation law, its parameters required when on and errors
  !> when off:
  !>   - saturation_evolution on: gammadot_s0 > 0, m_prime > 0 (one value each); tau_s one
  !>     value (tau_s0, one tau_s for all systems);
  !>   - cyclic_hardening on: cyclic_a in (0, 1), cyclic_c > 0 (one value each); a < 1 keeps
  !>     the threshold a (crss/tau_s)^c below 1 wherever a rate is evaluated (crss < tau_s), the
  !>     bound v2 checks at run time;
  !> With hardening kocks_mecking every switch is an error, even off (no modifiers, D60).
  !>   - precipitate_strengthening on: f_p in (0, 1] and r_p, with the cutting set (a_p, b_p)
  !>     and/or the bowing set (c_p, b_m), each complete; tau_p = a_p (f_p r_p/b_p)^(1/2)
  !>     (cutting), c_p G_m b_m/(L - 2 r_p) with L = r_p sqrt(pi/f_p) and G_m the Hill shear
  !>     modulus of the phase (bowing, L > 2 r_p), the smaller of both when both sets are given
  !>     (v2); tau_s > tau_0 + tau_p (system 1 for saturation_isotropic, every system otherwise).
  !>
  !> Arguments:
  !>   ph         (inout) type(cfg_phase_t)  the phase (law, tau_0 and tau_s read)
  !>   store      (in)    type(cfg_store_t)  validated records
  !>   i_phase    (in)    integer            phase
  !>   phase_rec  (in)    integer            the record of "phase i_phase"
  !> Locals:
  !>   switches   character(25)(3)           the three switches
  !>   params     character(11)(10)          their parameters (2 saturation evolution, 2 cyclic,
  !>                                         6 precipitation)
  !>   owner      integer(10)                switch of each parameter
  !>   recs       integer(10)                parameter records (0: absent)
  !>   srecs      integer(3)                 switch records (0: absent)
  !>   on         logical(3)                 switch values
  !>   i, a       integer                    index, slip system
  !>   cutting, bowing  logical               the complete cutting / bowing sets are given
  !>   tau_p_cut, tau_p_bow, spacing  real(rk)  cutting and bowing strengths [MPa], L [mm]
  !>   where      character(:), allocatable  "<file>: phase <i> (<type>, line <l>)" for messages
  subroutine build_modifiers(ph, store, i_phase, phase_rec)

    type(cfg_phase_t), intent(inout) :: ph
    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_phase, phase_rec

    character(len=25), parameter :: switches(3) = [character(len=25) :: &
      & 'saturation_evolution', 'cyclic_hardening', 'precipitate_strengthening']
    character(len=11), parameter :: params(10) = [character(len=11) :: 'gammadot_s0', &
      & 'm_prime', 'cyclic_a', 'cyclic_c', 'a_p', 'b_p', 'c_p', 'b_m', 'f_p', 'r_p']
    integer, parameter :: owner(10) = [1, 1, 2, 2, 3, 3, 3, 3, 3, 3]
    integer, parameter :: switch_rules(3) = [rule_saturation_evolution, rule_cyclic_hardening, &
      & rule_precipitate_strengthening]
    integer :: recs(10), srecs(3), i, a
    logical :: on(3), cutting, bowing
    real(rk) :: tau_p_cut, tau_p_bow, spacing
    character(len=:), allocatable :: where

    srecs = [(cfg_find(store, trim(switches(i)), i_phase, 0), i = 1, 3)]
    recs = [(cfg_find(store, trim(params(i)), i_phase, 0), i = 1, 10)]
    where = store%file_name // ': phase ' // int_text(i_phase) // ' (' // &
      & trim(ph%crystal_type) // ', line ' // line_text(store, phase_rec) // ')'

    ! Switches: only with a saturation law (Kocks-Mecking takes no modifiers, D60).
    on = .false.
    do i = 1, 3
      if (srecs(i) == 0) cycle
      on(i) = cfg_get_word(store, srecs(i), 1) == 'on'
      if (ph%model /= 'crystal_plasticity' .or. ph%hardening == 'none') call cfg_error(store, &
        & srecs(i), cfg_rule_message(rule_saturation_law, trim(switches(i)), &
        & detail='phase ' // int_text(i_phase)))
      if (ph%hardening == 'kocks_mecking') call cfg_error(store, srecs(i), &
        & cfg_rule_message(rule_saturation_law, trim(switches(i)), detail='phase ' // &
        & int_text(i_phase) // ' has hardening kocks_mecking, which takes no modifiers'))
    end do

    ! Parameters: required by their switch when on, errors when off.
    do i = 1, 10
      if (recs(i) > 0 .and. .not. on(owner(i))) call cfg_error(store, recs(i), &
        & cfg_rule_message(switch_rules(owner(i)), trim(params(i))))
      if (owner(i) /= 3 .and. on(owner(i)) .and. recs(i) == 0) call fatal_error('cfg', &
        & where // ': ' // trim(params(i)) // ' is required with ' // &
        & trim(switches(owner(i))) // ' on')
    end do

    ! Saturation evolution: one tau_s0.
    if (on(1)) then
      ph%saturation_evolution = .true.
      ph%gammadot_s0 = cfg_get_real(store, recs(1), 1)
      ph%m_prime = cfg_get_real(store, recs(2), 1)
      if (.not. (ph%gammadot_s0 > 0.0_rk)) call cfg_error(store, recs(1), &
        & 'gammadot_s0 must be positive')
      if (.not. (ph%m_prime > 0.0_rk)) call cfg_error(store, recs(2), 'm_prime must be ' // &
        & 'positive (saturation_evolution off for a constant saturation strength)')
      if (any(abs(ph%tau_s - ph%tau_s(1)) > 0.0_rk)) call cfg_error(store, &
        & cfg_find(store, 'tau_s', i_phase, 0), 'tau_s: saturation_evolution on takes one ' // &
        & 'value (tau_s0, one saturation strength for all systems)')
    end if

    ! Cyclic hardening: 0 < a < 1, c > 0.
    if (on(2)) then
      ph%cyclic = .true.
      ph%cyclic_a = cfg_get_real(store, recs(3), 1)
      ph%cyclic_c = cfg_get_real(store, recs(4), 1)
      if (.not. (ph%cyclic_a > 0.0_rk .and. ph%cyclic_a < 1.0_rk)) call cfg_error(store, &
        & recs(3), 'cyclic_a must lie between 0 and 1 (the critical accumulated slip a ' // &
        & '(crss/tau_s)^c must stay below 1)')
      if (.not. (ph%cyclic_c > 0.0_rk)) call cfg_error(store, recs(4), &
        & 'cyclic_c must be positive')
    end if

    ! Precipitate strengthening: complete sets, tau_p, and tau_s above the raised tau_0.
    if (on(3)) then
      ph%precipitation = .true.
      do i = 9, 10
        if (recs(i) == 0) call fatal_error('cfg', where // ': ' // trim(params(i)) // &
          & ' is required with precipitate_strengthening on')
      end do
      cutting = recs(5) > 0 .and. recs(6) > 0
      bowing = recs(7) > 0 .and. recs(8) > 0
      if ((recs(5) > 0 .neqv. recs(6) > 0) .or. (recs(7) > 0 .neqv. recs(8) > 0) .or. &
        & .not. (cutting .or. bowing)) call fatal_error('cfg', where // ': ' // &
        & 'precipitate_strengthening on needs the cutting parameters (a_p, b_p), the ' // &
        & 'bowing parameters (c_p, b_m), or both, each set complete')
      if (.not. (cfg_get_real(store, recs(9), 1) <= 1.0_rk)) call cfg_error(store, recs(9), &
        & 'f_p is a volume fraction (at most 1)')
      tau_p_cut = huge(1.0_rk)
      tau_p_bow = huge(1.0_rk)
      if (cutting) tau_p_cut = cfg_get_real(store, recs(5), 1)*sqrt(cfg_get_real(store, recs(9), &
        & 1)*cfg_get_real(store, recs(10), 1)/cfg_get_real(store, recs(6), 1))
      if (bowing) then
        spacing = cfg_get_real(store, recs(10), 1)*sqrt(4.0_rk*atan(1.0_rk)/ &
          & cfg_get_real(store, recs(9), 1))
        if (.not. (spacing > 2.0_rk*cfg_get_real(store, recs(10), 1))) call cfg_error(store, &
          & recs(10), 'precipitate bowing needs a spacing r_p sqrt(pi/f_p) larger than ' // &
          & '2 r_p (f_p < pi/4)')
        tau_p_bow = cfg_get_real(store, recs(7), 1)*shear_modulus_hill(ph%stiff_el_c)* &
          & cfg_get_real(store, recs(8), 1)/(spacing - 2.0_rk*cfg_get_real(store, recs(10), 1))
      end if
      ph%tau_p = min(tau_p_cut, tau_p_bow)
      ph%precip_regime = merge('cutting', 'bowing ', tau_p_cut <= tau_p_bow)
      if (ph%hardening == 'saturation_isotropic') then
        if (.not. (ph%tau_s(1) > ph%tau_0(1) + ph%tau_p)) call fatal_error('cfg', where // &
          & ': tau_s must be greater than tau_0 + the precipitate strength (system 1: ' // &
          & trim(real_text(ph%tau_0(1) + ph%tau_p)) // ')')
      else
        do a = 1, ph%num_slip
          if (.not. (ph%tau_s(a) > ph%tau_0(a) + ph%tau_p)) call fatal_error('cfg', where // &
            & ': tau_s must be greater than tau_0 + the precipitate strength for every ' // &
            & 'slip system (system ' // int_text(a) // ': ' // &
            & trim(real_text(ph%tau_0(a) + ph%tau_p)) // ')')
        end do
      end if
    end if

  end subroutine build_modifiers

  !> Kocks-Mecking parameters of a phase (step 7.2, ROADMAP D60, addendum_kocks_mecking.md §4):
  !> with hardening kocks_mecking, km_b, km_h_n, km_rho_0, km_k1, km_k20, km_gammadot_k0 and
  !> km_n are required (one value each; their ranges are checked by the keyword table:
  !> km_b, km_rho_0, km_gammadot_k0, km_n > 0; km_h_n, km_k1, km_k20 >= 0); with any other law
  !> or model they are errors. mu_0 is not an input (Hill's average, set by the driver).
  !>
  !> Arguments:
  !>   ph         (inout) type(cfg_phase_t)  the phase (model and law read)
  !>   store      (in)    type(cfg_store_t)  validated records
  !>   i_phase    (in)    integer            phase
  !>   phase_rec  (in)    integer            the record of "phase i_phase"
  !> Locals:
  !>   names      character(14)(7)           the keywords
  !>   recs       integer(7)                 their records (0: absent)
  !>   values     real(rk)(7)                their values
  !>   i          integer                    keyword index
  !>   where      character(:), allocatable  "<file>: phase <i> (<type>, line <l>)" for messages
  subroutine build_kocks_mecking(ph, store, i_phase, phase_rec)

    type(cfg_phase_t), intent(inout) :: ph
    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_phase, phase_rec

    character(len=14), parameter :: names(7) = [character(len=14) :: 'km_b', 'km_h_n', &
      & 'km_rho_0', 'km_k1', 'km_k20', 'km_gammadot_k0', 'km_n']
    integer :: recs(7), i
    real(rk) :: values(7)
    character(len=:), allocatable :: where

    recs = [(cfg_find(store, trim(names(i)), i_phase, 0), i = 1, 7)]
    where = store%file_name // ': phase ' // int_text(i_phase) // ' (' // &
      & trim(ph%crystal_type) // ', line ' // line_text(store, phase_rec) // ')'

    ! Other laws and models: no Kocks-Mecking keywords.
    if (ph%hardening /= 'kocks_mecking') then
      do i = 1, 7
        if (recs(i) > 0) call cfg_error(store, recs(i), cfg_rule_message(rule_kocks_mecking, &
          & trim(names(i)), detail='phase ' // int_text(i_phase)))
      end do
      return
    end if

    ! Kocks-Mecking: every parameter required.
    do i = 1, 7
      if (recs(i) == 0) call fatal_error('cfg', where // ': ' // trim(names(i)) // &
        & ' is required with hardening kocks_mecking')
      values(i) = cfg_get_real(store, recs(i), 1)
    end do
    ph%km_b = values(1)
    ph%km_h_n = values(2)
    ph%km_rho_0 = values(3)
    ph%km_k1 = values(4)
    ph%km_k20 = values(5)
    ph%km_gammadot_k0 = values(6)
    ph%km_n = values(7)

  end subroutine build_kocks_mecking

  !> Irradiation-softening parameters of a phase (step 7.3, ROADMAP D61,
  !> addendum_irradiation.md §6): the switch irradiation (on|off) is used only with hardening
  !> kocks_mecking, and only for fcc and bcc phases (hcp, bct: loop normals not defined yet);
  !> with irradiation on, irr_h_d, irr_eta, irr_loop_density, irr_loop_diameter (one value each;
  !> ranges by the keyword table) and irr_loop_normals (one or more family tokens of
  !> crystal_loop_normals, each once) are required; otherwise they are errors.
  !>
  !> Arguments:
  !>   ph         (inout) type(cfg_phase_t)  the phase (model, law and crystal type read)
  !>   store      (in)    type(cfg_store_t)  validated records
  !>   i_phase    (in)    integer            phase
  !>   phase_rec  (in)    integer            the record of "phase i_phase"
  !> Locals:
  !>   names      character(17)(5)           the parameters
  !>   recs       integer(5)                 their records (0: absent)
  !>   srec       integer                    the switch record (0: absent)
  !>   i, j       integer                    indices
  !>   word       character(:), allocatable  a family token
  !>   normals    real(rk)(3, :), allocatable  the normals of one family
  !>   tokens     character(8)(:), allocatable  the tokens given so far
  !>   where      character(:), allocatable  "<file>: phase <i> (<type>, line <l>)" for messages
  subroutine build_irradiation(ph, store, i_phase, phase_rec)

    type(cfg_phase_t), intent(inout) :: ph
    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_phase, phase_rec

    character(len=17), parameter :: names(5) = [character(len=17) :: 'irr_h_d', 'irr_eta', &
      & 'irr_loop_density', 'irr_loop_diameter', 'irr_loop_normals']
    integer :: recs(5), srec, i, j
    character(len=:), allocatable :: word, where
    real(rk), allocatable :: normals(:, :)
    character(len=8), allocatable :: tokens(:)

    srec = cfg_find(store, 'irradiation', i_phase, 0)
    recs = [(cfg_find(store, trim(names(i)), i_phase, 0), i = 1, 5)]
    where = store%file_name // ': phase ' // int_text(i_phase) // ' (' // &
      & trim(ph%crystal_type) // ', line ' // line_text(store, phase_rec) // ')'
    allocate (ph%irr_loop_normals_c(3, 0))

    ! The switch: only with Kocks-Mecking, and (on) only for cubic crystals.
    if (srec > 0) then
      if (ph%hardening /= 'kocks_mecking') call cfg_error(store, srec, &
        & cfg_rule_message(rule_irradiation_law, 'irradiation', detail='phase ' // &
        & int_text(i_phase)))
      ph%irradiation = cfg_get_word(store, srec, 1) == 'on'
      if (ph%irradiation .and. ph%crystal_type /= 'fcc' .and. ph%crystal_type /= 'bcc') &
        & call cfg_error(store, srec, cfg_rule_message(rule_irradiation_crystal, &
        & 'irradiation on', detail='phase ' // int_text(i_phase) // ' has crystal_type ' // &
        & trim(ph%crystal_type)))
    end if

    ! Off: no parameters.
    if (.not. ph%irradiation) then
      do i = 1, 5
        if (recs(i) > 0) call cfg_error(store, recs(i), &
          & cfg_rule_message(rule_irradiation_parameters, trim(names(i)), detail='phase ' // &
          & int_text(i_phase)))
      end do
      return
    end if

    ! On: every parameter required.
    do i = 1, 5
      if (recs(i) == 0) call fatal_error('cfg', where // ': ' // trim(names(i)) // &
        & ' is required with irradiation on')
    end do
    ph%irr_h_d = cfg_get_real(store, recs(1), 1)
    ph%irr_eta = cfg_get_real(store, recs(2), 1)
    ph%irr_loop_density = cfg_get_real(store, recs(3), 1)
    ph%irr_loop_diameter = cfg_get_real(store, recs(4), 1)

    ! Loop-normal families: tokens of the crystal type, each once; their normals together.
    allocate (tokens(0))
    do i = 1, size(store%records(recs(5))%values)
      word = cfg_get_word(store, recs(5), i)
      call crystal_loop_normals(crystal_structure_from_name(trim(ph%crystal_type)), word, &
        & normals)
      if (size(normals, 2) == 0) call cfg_error(store, recs(5), 'irr_loop_normals: "' // &
        & word // '" is not a loop-normal family (expected: <111>, <100>, <110>)')
      do j = 1, size(tokens)
        if (tokens(j) == word) call cfg_error(store, recs(5), 'irr_loop_normals: ' // word // &
          & ' given twice')
      end do
      tokens = [character(len=8) :: tokens, word]
      ph%irr_loop_normals_c = reshape([ph%irr_loop_normals_c, normals], &
        & [3, size(ph%irr_loop_normals_c, 2) + size(normals, 2)])
    end do

  end subroutine build_irradiation

  !> The user's interaction matrix of a phase: its interaction_matrix lines, in file order, are
  !> the rows 1..N in the slip-system order; exactly N lines of N numbers >= 0.
  !>
  !> Arguments:
  !>   store     (in)  type(cfg_store_t)  validated records
  !>   i_phase   (in)  integer            phase
  !>   num       (in)  integer            number of slip systems N
  !>   mat       (out) real(rk)(N, N)     the matrix
  !> Locals:
  !>   i_rec     integer                  record index
  !>   row, col  integer                  row (line count) and column
  !>   ios       integer                  read status
  !>   x         real(rk)                 an entry
  !>   word      character(:), allocatable  a value as written
  !>   last_rec  integer                  the last interaction_matrix record (for messages)
  subroutine read_interaction_matrix(store, i_phase, num, mat)

    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_phase, num
    real(rk), intent(out) :: mat(num, num)

    integer :: i_rec, row, col, ios, last_rec
    real(rk) :: x
    character(len=:), allocatable :: word

    row = 0
    last_rec = 0
    do i_rec = 1, store%num_records
      if (cfg_keywords(store%records(i_rec)%key)%name /= 'interaction_matrix') cycle
      if (store%records(i_rec)%phase /= i_phase) cycle
      last_rec = i_rec
      row = row + 1
      if (row > num) call cfg_error(store, i_rec, 'interaction_matrix: more than ' // &
        & int_text(num) // ' rows (one per slip system)')
      if (size(store%records(i_rec)%values) /= num) call cfg_error(store, i_rec, &
        & 'interaction_matrix: row ' // int_text(row) // ' has ' // &
        & int_text(size(store%records(i_rec)%values)) // ' values; give ' // int_text(num) // &
        & ' (one per slip system)')
      do col = 1, num
        word = cfg_get_word(store, i_rec, col)
        read (word, *, iostat=ios) x
        if (ios /= 0 .or. verify(word, '+-.0123456789eEdD') /= 0) call cfg_error(store, i_rec, &
          & 'interaction_matrix: "' // word // '" is not a number')
        if (x < 0.0_rk) call cfg_error(store, i_rec, 'interaction_matrix: entry (' // &
          & int_text(row) // ', ' // int_text(col) // ') is negative; entries must be >= 0')
        mat(row, col) = x
      end do
    end do
    if (row < num) call cfg_error(store, last_rec, 'interaction_matrix: ' // int_text(row) // &
      & ' rows; give ' // int_text(num) // ' (one per slip system)')

  end subroutine read_interaction_matrix

  !> Short text of a real, for messages.
  !>
  !> Arguments:
  !>   x     (in)  real(rk)  the value
  !> Locals:
  !>   buf   character(32)   formatting buffer
  function real_text(x) result(text)

    real(rk), intent(in) :: x
    character(len=:), allocatable :: text

    character(len=32) :: buf

    write (buf, '(g0.6)') x
    text = trim(buf)

  end function real_text

  !> Values of a per-family keyword expanded to one per slip system (dev/cfg_keywords.md §2.1
  !> counting rule): one value per active family (in the families order; each family's
  !> systems get its value) or one per slip system (taken as given); any other count is an
  !> error naming the accepted counts.
  !>
  !> Arguments:
  !>   store      (in)  type(cfg_store_t)  validated records
  !>   i_rec      (in)  integer            the keyword's record
  !>   name       (in)  character          keyword, for messages
  !>   structure  (in)  integer            crystal_mod structure code
  !>   families   (in)  integer(:)         active families, in order
  !>   num_slip   (in)  integer            slip systems of the active families
  !>   one_for_all (in, optional) logical    also accept one value for all systems (hardening
  !>                                         parameters; default false)
  !>   values     (out, result) real(rk)(num_slip)  one value per slip system
  !> Locals:
  !>   nvals      integer                  number of values given
  !>   i, k, a    integer                  family, system counter, system of the family
  !>   single     logical                  one value for all systems is accepted
  function per_system_values(store, i_rec, name, structure, families, num_slip, one_for_all) &
    & result(values)

    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_rec
    character(len=*), intent(in) :: name
    integer, intent(in) :: structure
    integer, intent(in) :: families(:)
    integer, intent(in) :: num_slip
    logical, intent(in), optional :: one_for_all
    real(rk) :: values(num_slip)

    integer :: nvals, i, k, a
    logical :: single

    single = .false.
    if (present(one_for_all)) single = one_for_all
    nvals = size(store%records(i_rec)%values)
    if (single .and. nvals == 1) then
      ! One value for all systems (hardening parameters, addendum_hardening.md §1).
      values = cfg_get_real(store, i_rec, 1)
    else if (nvals == size(families)) then
      ! One per family: repeated over the family's systems, family after family.
      k = 0
      do i = 1, size(families)
        do a = 1, crystal_family_size(structure, families(i))
          k = k + 1
          values(k) = cfg_get_real(store, i_rec, i)
        end do
      end do
    else if (nvals == num_slip) then
      ! One per system, as given.
      values = [(cfg_get_real(store, i_rec, k), k = 1, num_slip)]
    else
      if (single) then
        call cfg_error(store, i_rec, name // ': ' // int_text(nvals) // ' value(s); give one ' &
          & // 'value, one per active slip family (' // int_text(size(families)) // ') or ' // &
          & 'one per slip system (' // int_text(num_slip) // ')')
      else
        call cfg_error(store, i_rec, name // ': ' // int_text(nvals) // ' value(s); give one ' &
          & // 'per active slip family (' // int_text(size(families)) // ') or one per slip ' // &
          & 'system (' // int_text(num_slip) // ')')
      end if
    end if

  end function per_system_values

  !> The slip family tokens of a structure, comma separated, for messages.
  !>
  !> Arguments:
  !>   structure  (in)  integer  crystal_mod structure code
  !>   text       (out, result) character(:), allocatable  the list
  !> Locals:
  !>   i          integer  family
  function family_list(structure) result(text)

    integer, intent(in) :: structure
    character(len=:), allocatable :: text

    integer :: i

    text = ''
    do i = 1, crystal_num_families(structure)
      if (i > 1) text = text // ', '
      text = text // crystal_family_token(structure, i)
    end do

  end function family_list

  !> Boundary conditions and loading rows, with the checks of the module header.
  !>
  !> Arguments:
  !>   config  (inout) type(config_t)     the configuration
  !>   store   (in)    type(cfg_store_t)  validated records
  !> Locals:
  !>   i_rec, i_step, comp_idx  integer   record, step and component indices
  !>   first_row  integer                 first row setting the loading axis (one-component modes)
  !>   name       character(:), allocatable  keyword of a record
  !>   is_preset  logical                 a preset (uniaxial or triaxial), not general
  !>   is_triax   logical                 triaxial preset
  !>   is_periodic  logical               periodic conditions (step 8.3b)
  subroutine build_loading(config, store)

    type(config_t), intent(inout) :: config
    type(cfg_store_t), intent(in) :: store

    integer :: i_rec, i_step, comp_idx, first_row
    character(len=:), allocatable :: name
    logical :: is_preset, is_triax, is_periodic

    ! Boundary-condition type and steps (required).
    config%bc_type = cfg_get_word(store, require(store, 'boundary_conditions'), 1)
    is_preset = config%bc_type /= 'general' .and. config%bc_type /= 'periodic'
    is_triax = config%bc_type == 'triaxial'
    is_periodic = config%bc_type == 'periodic'
    config%num_steps = cfg_get_int(store, require(store, 'number_of_steps'), 1)
    allocate (config%strain_rate(3, config%num_steps), config%target_strain(3, config%num_steps))
    allocate (config%has_strain_rate(3, config%num_steps))
    allocate (config%has_target_strain(3, config%num_steps))
    allocate (config%target_time(config%num_steps), config%has_target_time(config%num_steps))
    allocate (config%target_load(3, config%num_steps), config%target_stress(3, config%num_steps))
    allocate (config%has_target_load(3, config%num_steps))
    allocate (config%has_target_stress(3, config%num_steps))
    allocate (config%target_kind(config%num_steps), config%rate_kind(config%num_steps))
    allocate (config%load_rate(3, config%num_steps), config%stress_rate(3, config%num_steps))
    allocate (config%has_load_rate(3, config%num_steps))
    allocate (config%has_stress_rate(3, config%num_steps))
    allocate (config%target_duration(config%num_steps))
    allocate (config%has_target_duration(config%num_steps))
    config%load_rate = 0.0_rk
    config%stress_rate = 0.0_rk
    config%has_load_rate = .false.
    config%has_stress_rate = .false.
    config%target_duration = 0.0_rk
    config%has_target_duration = .false.
    config%rate_kind = rate_by_strain
    config%strain_rate = 0.0_rk
    config%target_strain = 0.0_rk
    config%target_load = 0.0_rk
    config%target_stress = 0.0_rk
    config%target_time = 0.0_rk
    config%has_strain_rate = .false.
    config%has_target_strain = .false.
    config%has_target_load = .false.
    config%has_target_stress = .false.
    config%has_target_time = .false.
    config%target_kind = 0

    ! The keywords of periodic conditions only with them (step 8.3b).
    if (.not. is_periodic) then
      do i_rec = 1, store%num_records
        name = trim(cfg_keywords(store%records(i_rec)%key)%name)
        if (name == 'periodic_kinematics' .or. name == 'periodic_preconditioner') &
          & call cfg_error(store, i_rec, cfg_rule_message(rule_periodic_keywords, name))
      end do
    end if

    ! Component rows: values per component and step. All rows of a uniaxial preset, and the
    ! target rows of general, are on one component (the loading axis); triaxial may use
    ! several; general takes no rate row (strain, load or stress rate, step 8.2c); shear
    ! components only with periodic conditions (step 8.3b).
    first_row = 0
    do i_rec = 1, store%num_records
      if (is_periodic) exit
      name = trim(cfg_keywords(store%records(i_rec)%key)%name)
      if (name /= 'strain_rate' .and. name /= 'target_strain' .and. name /= 'target_load' &
        & .and. name /= 'target_stress' .and. name /= 'load_rate' .and. &
        & name /= 'stress_rate') cycle
      if (mod(store%records(i_rec)%comp, 11) /= 0) call cfg_error(store, i_rec, &
        & cfg_rule_message(rule_shear_components, 'shear components', prefix=name // ' ' // &
        & int_text(store%records(i_rec)%comp), plural=.true.))
      comp_idx = store%records(i_rec)%comp/11
      if ((name == 'strain_rate' .or. name == 'load_rate' .or. name == 'stress_rate') .and. &
        & config%bc_type == 'general') call cfg_error(store, i_rec, &
        & cfg_rule_message(rule_rates_general, name))
      config%comp_rows(comp_idx) = .true.
      if (is_triax) then
        ! Several components allowed; the control component is named below.
        continue
      else if (first_row == 0) then
        first_row = i_rec
        config%control_axis = comp_idx
      else if (comp_idx /= config%control_axis) then
        call cfg_error(store, i_rec, name // ' ' // int_text(11*comp_idx) // ': ' // &
          & trim(config%bc_type) // ' loads one component, and line ' // &
          & line_text(store, first_row) // ' gives ' // int_text(11*config%control_axis) // &
          & ' (several components: boundary_conditions triaxial)')
      end if
      do i_step = 1, config%num_steps
        if (cfg_is_dash(store, i_rec, i_step)) cycle
        select case (name)
        case ('strain_rate')
          config%strain_rate(comp_idx, i_step) = cfg_get_real(store, i_rec, i_step)
          config%has_strain_rate(comp_idx, i_step) = .true.
        case ('target_strain')
          config%target_strain(comp_idx, i_step) = cfg_get_real(store, i_rec, i_step)
          config%has_target_strain(comp_idx, i_step) = .true.
        case ('target_load')
          config%target_load(comp_idx, i_step) = cfg_get_real(store, i_rec, i_step)
          config%has_target_load(comp_idx, i_step) = .true.
        case ('load_rate')
          config%load_rate(comp_idx, i_step) = cfg_get_real(store, i_rec, i_step)
          config%has_load_rate(comp_idx, i_step) = .true.
        case ('stress_rate')
          config%stress_rate(comp_idx, i_step) = cfg_get_real(store, i_rec, i_step)
          config%has_stress_rate(comp_idx, i_step) = .true.
        case default
          config%target_stress(comp_idx, i_step) = cfg_get_real(store, i_rec, i_step)
          config%has_target_stress(comp_idx, i_step) = .true.
        end select
      end do
    end do

    ! Time targets.
    i_rec = cfg_find(store, 'target_time', 0, 0)
    if (i_rec > 0) then
      do i_step = 1, config%num_steps
        if (cfg_is_dash(store, i_rec, i_step)) cycle
        config%target_time(i_step) = cfg_get_real(store, i_rec, i_step)
        config%has_target_time(i_step) = .true.
      end do
    end if

    ! Step lengths (step 8.2c, Q-T7).
    i_rec = cfg_find(store, 'target_duration', 0, 0)
    if (i_rec > 0) then
      do i_step = 1, config%num_steps
        if (cfg_is_dash(store, i_rec, i_step)) cycle
        config%target_duration(i_step) = cfg_get_real(store, i_rec, i_step)
        config%has_target_duration(i_step) = .true.
      end do
    end if

    if (is_periodic) then
      ! Periodic conditions (step 8.3b): their rows, control component and step checks.
      call build_periodic(config, store)
    else
      ! The control component (D45).
      call build_control_component(config, store, is_triax)

      ! Each step: the control component's target and rate, the other components' rate or
      ! target (§1a rule 3, D45).
      call check_step_controls(config, store, is_preset)
    end if

    ! Increments of each step: exactly one of number_of_incrs, dtime, dstrain (§1a rule 6);
    ! dtime_min in load- and stress-target steps.
    call build_increments(config, store)

    ! Load and stress targets: time-step factor, necking (step 8.1), tolerances and landing
    ! attempts (step 8.2a).
    call build_load_options(config, store)

    ! Print flags per step (default: every step) and run limits (defaults from the table;
    ! max_strain and max_eqstrain: no limit unless given).
    allocate (config%print_results(config%num_steps))
    config%print_results = .true.
    i_rec = cfg_find(store, 'print_results', 0, 0)
    if (i_rec > 0) config%print_results = [(cfg_get_int(store, i_rec, i_step) == 1, &
      & i_step = 1, config%num_steps)]
    config%max_incr = nint(real_or_default(store, 'max_incr'))
    config%max_total_time = real_or_default(store, 'max_total_time')
    i_rec = cfg_find(store, 'max_strain', 0, 0)
    config%has_max_strain = i_rec > 0
    if (i_rec > 0) config%max_strain = cfg_get_real(store, i_rec, 1)
    i_rec = cfg_find(store, 'max_eqstrain', 0, 0)
    config%has_max_eqstrain = i_rec > 0
    if (i_rec > 0) config%max_eqstrain = cfg_get_real(store, i_rec, 1)

    ! Velocity conditions: general only, and required there.
    call build_vel_bcs(config, store)
    if (is_preset .and. size(config%vel_bcs) > 0) call fatal_error('cfg', store%file_name // &
      & ':' // int_text(int(config%vel_bcs(1)%line)) // ': ' // &
      & cfg_rule_message(rule_set_bc_preset, 'set_bc', detail='given: boundary_conditions ' // &
      & trim(config%bc_type)))
    if (config%bc_type == 'general' .and. size(config%vel_bcs) == 0) call fatal_error('cfg', &
      & store%file_name // ': boundary_conditions general needs set_bc vel lines')

    ! Multi-point constraints (step 8.3a), and their interactions (Q-P8: errors; with periodic
    ! conditions as v2-fixed C68).
    call build_mpcs(config, store)
    if (size(config%mpcs) > 0) then
      if (is_periodic) call fatal_error('cfg', store%file_name // ':' // &
        & int_text(int(config%mpcs(1)%line)) // ': ' // &
        & cfg_rule_message(rule_mpc_periodic, 'multi-point constraints', plural=.true.))
      if (is_triax) call fatal_error('cfg', store%file_name // ':' // &
        & int_text(int(config%mpcs(1)%line)) // ': ' // &
        & cfg_rule_message(rule_mpc_triaxial, 'multi-point constraints', plural=.true.))
      if (any(config%target_kind == target_by_load .or. config%target_kind == &
        & target_by_stress .or. config%rate_kind == rate_by_load .or. config%rate_kind == &
        & rate_by_stress)) call fatal_error('cfg', store%file_name // ':' // &
        & int_text(int(config%mpcs(1)%line)) // ': ' // &
        & cfg_rule_message(rule_mpc_load_control, 'multi-point constraints', plural=.true.))
    end if

  end subroutine build_loading

  !> The multi-point constraint lines (step 8.3a, D73): "set_mpc1 vel <nset> <x|y|z> [<x|y|z>
  !> [<x|y|z>]]" (coefficient 1, offset 0) and "set_mpc vel <nset> <x|y|z> <coeff> <offset>
  !> [...]" (up to three triples; v2 read_input_cfg.f90::parse_set_mpc1, parse_set_mpc, with
  !> "vel" added to set_mpc, Q-P7). A zero coefficient is an error (v2 read_mpcs); one relation
  !> per node set and direction.
  !>
  !> Arguments:
  !>   config  (inout) type(config_t)     the configuration
  !>   store   (in)    type(cfg_store_t)  validated records
  !> Locals:
  !>   i_rec, i_dir, i, ios  integer      record, direction and relation indices, read status
  !>   nvals   integer                    values of the line
  !>   width   integer                    values per direction (1 set_mpc1, 3 set_mpc)
  !>   name    character(:), allocatable  keyword of the record
  !>   mpc     type(cfg_mpc_t)            one relation
  !>   word    character(:), allocatable  a value
  !>   grown   type(cfg_mpc_t)(:)         larger copy of the list
  subroutine build_mpcs(config, store)

    type(config_t), intent(inout) :: config
    type(cfg_store_t), intent(in) :: store

    integer :: i_rec, i_dir, i, ios, nvals, width
    character(len=:), allocatable :: name
    type(cfg_mpc_t) :: mpc
    character(len=:), allocatable :: word
    type(cfg_mpc_t), allocatable :: grown(:)

    allocate (config%mpcs(0))
    do i_rec = 1, store%num_records
      name = trim(cfg_keywords(store%records(i_rec)%key)%name)
      if (name /= 'set_mpc1' .and. name /= 'set_mpc') cycle
      width = merge(1, 3, name == 'set_mpc1')

      ! The kind (vel) and the number of values.
      nvals = size(store%records(i_rec)%values)
      word = cfg_get_word(store, i_rec, 1)
      if (word /= 'vel') then
        if (name == 'set_mpc') call cfg_error(store, i_rec, 'set_mpc needs vel: set_mpc ' // &
          & 'vel <nset> <x|y|z> <coefficient> <offset> [...]')
        call cfg_error(store, i_rec, 'set_mpc1: invalid kind "' // word // '" (expected vel)')
      end if
      if (nvals < 2 + width .or. nvals > 2 + 3*width .or. mod(nvals - 2, width) /= 0) then
        if (name == 'set_mpc1') call cfg_error(store, i_rec, 'set_mpc1 vel: expected ' // &
          & '<nset> and one to three directions <x|y|z>')
        call cfg_error(store, i_rec, 'set_mpc vel: expected <nset> and one to three ' // &
          & '<x|y|z> <coefficient> <offset> triples')
      end if

      ! Each direction (with its coefficient and offset for set_mpc).
      mpc%nset = cfg_get_word(store, i_rec, 2)
      mpc%line = store%records(i_rec)%line
      do i_dir = 1, (nvals - 2)/width
        word = cfg_get_word(store, i_rec, 3 + width*(i_dir - 1))
        mpc%comp = findloc(axis_letter, word, dim=1)
        if (mpc%comp == 0) call cfg_error(store, i_rec, name // ' vel: invalid direction "' // &
          & word // '" (expected x, y or z)')
        mpc%coeff = 1.0_rk
        mpc%offset = 0.0_rk
        if (width == 3) then
          word = cfg_get_word(store, i_rec, 4 + width*(i_dir - 1))
          read (word, *, iostat=ios) mpc%coeff
          if (ios /= 0 .or. verify(word, '+-.0123456789eEdD') /= 0) call cfg_error(store, &
            & i_rec, 'set_mpc vel: "' // word // '" is not a number')
          if (.not. (abs(mpc%coeff) > 0.0_rk)) call cfg_error(store, i_rec, 'set_mpc vel: ' // &
            & 'the coefficient must be nonzero')
          word = cfg_get_word(store, i_rec, 5 + width*(i_dir - 1))
          read (word, *, iostat=ios) mpc%offset
          if (ios /= 0 .or. verify(word, '+-.0123456789eEdD') /= 0) call cfg_error(store, &
            & i_rec, 'set_mpc vel: "' // word // '" is not a number')
        end if
        do i = 1, size(config%mpcs)
          if (config%mpcs(i)%nset == mpc%nset .and. config%mpcs(i)%comp == mpc%comp) then
            call cfg_error(store, i_rec, name // ' vel: ' // mpc%nset // ' ' // &
              & axis_letter(mpc%comp) // ' given twice (lines ' // &
              & int_text(int(config%mpcs(i)%line)) // ' and ' // int_text(int(mpc%line)) // ')')
          end if
        end do

        ! Append (element by element, as build_vel_bcs).
        allocate (grown(size(config%mpcs) + 1))
        do i = 1, size(config%mpcs)
          grown(i) = config%mpcs(i)
        end do
        grown(size(grown)) = mpc
        call move_alloc(grown, config%mpcs)
      end do
    end do

  end subroutine build_mpcs

  !> Periodic conditions (step 8.3b, D71-D73; dev/theory/addendum_periodic_mpc.md §3.3, §4):
  !> the kinematics and preconditioner words; a strain_rate row per imposed component of the
  !> rate of the macroscopic deformation gradient Fbar-dot, with a value in every step (a
  !> component is imposed or free for the whole run; free: no row; Q-P5); small_strain takes
  !> 11 22 33 23 13 12 (a shear row imposes ij and ji, v2's strainrate), deformation_gradient_rate
  !> the nine ij (v2's velgrad); control_component (required) is an imposed component, the only
  !> one with target_strain rows (a component of Fbar - I); each step ends at exactly one of its
  !> target_strain, target_time or target_duration; no load or stress rows (Q-P10). The steps
  !> are then planned as a one-component problem: control_axis 1 with the control component's
  !> rate and targets (comp_label(1) its name).
  !>
  !> Arguments:
  !>   config  (inout) type(config_t)     the configuration (steps and time targets read)
  !>   store   (in)    type(cfg_store_t)  validated records
  !> Locals:
  !>   i_rec, i_step  integer             record, step
  !>   name    character(:), allocatable  keyword of a record
  !>   comp    integer                    a component ij as 10 i + j
  !>   c, ct   integer                    its index 3 (i - 1) + j, and that of ji
  !>   small   logical                    small_strain kinematics
  !>   ctrl_rec, target_rec  integer      control_component and its target_strain records
  !>   label   character(:), allocatable  the control component as text
  !>   step    character(:), allocatable  "step <i>" with the component, for messages
  !>   given   logical(3)                 the step's target_strain, target_time,
  !>                                      target_duration
  !>   names   character(32)(3)           those rows as text
  !>   recs    integer(3)                 their records
  !>   k, first  integer                  a target, the first found
  !>   kinds   integer(3)                 the target_by_* code of each of the three targets
  subroutine build_periodic(config, store)

    type(config_t), intent(inout) :: config
    type(cfg_store_t), intent(in) :: store

    integer :: i_rec, i_step, comp, c, ct, ctrl_rec, target_rec, recs(3), k, first
    integer, parameter :: kinds(3) = [target_by_strain, target_by_time, target_by_duration]
    character(len=:), allocatable :: name, label, step
    logical :: small, given(3)
    character(len=32) :: names(3)

    ! Kinematics and preconditioner.
    config%periodic_kinematics = word_or_default(store, 'periodic_kinematics')
    config%periodic_preconditioner = word_or_default(store, 'periodic_preconditioner')
    small = config%periodic_kinematics == 'small_strain'
    allocate (config%per_rate(9, config%num_steps))
    config%per_rate = 0.0_rk
    config%per_imposed = .false.

    ! Imposed components: strain_rate rows, a value in every step; no load or stress rows.
    do i_rec = 1, store%num_records
      name = trim(cfg_keywords(store%records(i_rec)%key)%name)
      if (name == 'target_load' .or. name == 'target_stress' .or. name == 'load_rate' .or. &
        & name == 'stress_rate') call cfg_error(store, i_rec, &
        & cfg_rule_message(rule_load_rows_periodic, name))
      if (name /= 'strain_rate' .and. name /= 'target_strain') cycle
      comp = store%records(i_rec)%comp
      if (small .and. (comp == 32 .or. comp == 31 .or. comp == 21)) call cfg_error(store, &
        & i_rec, cfg_rule_message(rule_small_strain_components, name // ' ' // int_text(comp)))
      if (name /= 'strain_rate') cycle
      c = 3*(comp/10 - 1) + mod(comp, 10)
      do i_step = 1, config%num_steps
        if (cfg_is_dash(store, i_rec, i_step)) call cfg_error(store, i_rec, 'strain_rate ' // &
          & int_text(comp) // ': step ' // int_text(i_step) // ' has "-": with periodic ' // &
          & 'conditions a component is imposed (a value in every step) or free (no row)')
        config%per_rate(c, i_step) = cfg_get_real(store, i_rec, i_step)
      end do
      config%per_imposed(c) = .true.
      ! Small strain: a shear row imposes ij and ji (v2's strainrate).
      ct = 3*(mod(comp, 10) - 1) + comp/10
      if (small) then
        config%per_rate(ct, :) = config%per_rate(c, :)
        config%per_imposed(ct) = .true.
      end if
    end do

    ! The control component: required, imposed; the only one with target_strain rows.
    ctrl_rec = cfg_find(store, 'control_component', 0, 0)
    if (ctrl_rec == 0) call fatal_error('cfg', store%file_name // ': boundary_conditions ' // &
      & 'periodic needs control_component <ij> (an imposed component, whose target_strain ' // &
      & 'or the time ends each step)')
    label = cfg_get_word(store, ctrl_rec, 1)
    comp = cfg_get_int(store, ctrl_rec, 1)
    if (small .and. (comp == 32 .or. comp == 31 .or. comp == 21)) call cfg_error(store, &
      & ctrl_rec, cfg_rule_message(rule_small_strain_components, 'control_component ' // label))
    config%per_control = 3*(comp/10 - 1) + mod(comp, 10)
    if (cfg_find(store, 'strain_rate', 0, comp) == 0) call cfg_error(store, ctrl_rec, &
      & 'control_component ' // label // ' has no strain_rate row (the control component ' // &
      & 'is imposed)')
    do i_rec = 1, store%num_records
      if (cfg_keywords(store%records(i_rec)%key)%name /= 'target_strain') cycle
      if (store%records(i_rec)%comp /= comp) call cfg_error(store, i_rec, 'target_strain ' // &
        & int_text(store%records(i_rec)%comp) // ': with boundary_conditions periodic only ' // &
        & 'the control component ' // label // ' has strain targets')
    end do

    ! The steps as a one-component problem: the control component's rate and targets.
    config%control_axis = 1
    config%comp_rows = [.true., .false., .false.]
    config%comp_label(1) = label
    config%strain_rate(1, :) = config%per_rate(config%per_control, :)
    config%has_strain_rate(1, :) = .true.
    target_rec = cfg_find(store, 'target_strain', 0, comp)
    if (target_rec > 0) then
      do i_step = 1, config%num_steps
        if (cfg_is_dash(store, target_rec, i_step)) cycle
        config%target_strain(1, i_step) = cfg_get_real(store, target_rec, i_step)
        config%has_target_strain(1, i_step) = .true.
      end do
    end if

    ! Each step: exactly one target (the control component's strain, a time or a duration).
    recs = [target_rec, cfg_find(store, 'target_time', 0, 0), &
      & cfg_find(store, 'target_duration', 0, 0)]
    names = [character(len=32) :: 'target_strain ' // label, 'target_time', 'target_duration']
    do i_step = 1, config%num_steps
      step = store%file_name // ': step ' // int_text(i_step) // ', component ' // label
      given = [config%has_target_strain(1, i_step), config%has_target_time(i_step), &
        & config%has_target_duration(i_step)]
      first = 0
      do k = 1, 3
        if (.not. given(k)) cycle
        if (first > 0) call fatal_error('cfg', step // ': two targets (' // &
          & trim(names(first)) // ', line ' // line_text(store, recs(first)) // ', and ' // &
          & trim(names(k)) // ', line ' // line_text(store, recs(k)) // '); give one, "-" ' // &
          & 'for the other')
        first = k
      end do
      if (first == 0) call fatal_error('cfg', step // ': no target (' // trim(names(1)) // &
        & ', target_time or target_duration)')
      config%target_kind(i_step) = kinds(first)
      config%rate_kind(i_step) = rate_by_strain
    end do

  end subroutine build_periodic

  !> The control component (D45): control_component, required with triaxial (with rows on
  !> it); with a uniaxial preset the loaded component, which a given control_component must
  !> match; an error with general (whose loading axis is that of its target_strain rows).
  !>
  !> Arguments:
  !>   config    (inout) type(config_t)     the configuration (component rows read)
  !>   store     (in)    type(cfg_store_t)  validated records
  !>   is_triax  (in)    logical            triaxial preset
  !> Locals:
  !>   i_rec     integer                    control_component record (0: absent)
  !>   comp_idx  integer                    its component index (1..3)
  subroutine build_control_component(config, store, is_triax)

    type(config_t), intent(inout) :: config
    type(cfg_store_t), intent(in) :: store
    logical, intent(in) :: is_triax

    integer :: i_rec, comp_idx

    i_rec = cfg_find(store, 'control_component', 0, 0)
    comp_idx = 0
    if (i_rec > 0) then
      if (mod(cfg_get_int(store, i_rec, 1), 11) /= 0) call cfg_error(store, i_rec, &
        & cfg_rule_message(rule_shear_components, 'shear components', prefix='control_component ' &
        & // cfg_get_word(store, i_rec, 1), plural=.true.))
      comp_idx = cfg_get_int(store, i_rec, 1)/11
    end if

    if (is_triax) then
      ! Required, and on a component with rows.
      if (i_rec == 0) call fatal_error('cfg', store%file_name // ': boundary_conditions ' // &
        & 'triaxial needs control_component <11, 22 or 33> (the component that sets the ' // &
        & 'step durations)')
      if (.not. config%comp_rows(comp_idx)) call cfg_error(store, i_rec, 'control_component ' &
        & // int_text(11*comp_idx) // ' has no strain_rate or target_strain rows')
      config%control_axis = comp_idx
    else if (config%bc_type == 'general') then
      ! Not used: the steps end at their target_strain or target_time.
      if (i_rec > 0) call cfg_error(store, i_rec, &
        & cfg_rule_message(rule_control_general, 'control_component'))
    else
      ! Uniaxial: the loaded component; a given value must match it.
      if (i_rec > 0 .and. config%control_axis > 0 .and. comp_idx /= config%control_axis) &
        & call cfg_error(store, i_rec, 'control_component ' // int_text(11*comp_idx) // &
        & ' is not the loaded component ' // int_text(11*config%control_axis) // ' of ' // &
        & trim(config%bc_type))
    end if

  end subroutine build_control_component

  !> Targets and rates of each step (dev/cfg_keywords.md §1a rule 3, D45): exactly one target,
  !> target_strain, target_load or target_stress on the control component, target_time or
  !> target_duration (step 8.1: any kind in any step; target_duration step 8.2c); with a preset,
  !> exactly one rate, strain_rate, load_rate or stress_rate on the control component (general
  !> takes no rate row: its velocities are the set_bc lines); every other component with rows
  !> (triaxial), exactly one of strain_rate, target_strain, target_load, target_stress,
  !> load_rate and stress_rate (load and stress rows: triaxial load control, steps 8.2b, 8.2c,
  !> addendum_triaxial_load_control.md §3.1), not a target_strain in a step whose duration is
  !> not known (a load or stress target at a strain rate, addendum_load_targets.md §3.2; a strain
  !> target at a load or stress rate, where target_load and target_stress are not allowed
  !> either). Load and stress rates of the control component (step 8.2c, addendum §3.1, §3.2):
  !> a load_rate does not end at a target_stress, nor a stress_rate at a target_load; a rate 0
  !> (a dwell) only with target_time or target_duration, not in step 1 (v2: a dwell before the
  !> first load step). Sets target_kind and rate_kind. The messages name the step, the
  !> component and the rows with their lines.
  !>
  !> Arguments:
  !>   config     (inout) type(config_t)     the configuration (control rows read)
  !>   store      (in)    type(cfg_store_t)  validated records
  !>   is_preset  (in)    logical            a preset (uniaxial or triaxial)
  !> Locals:
  !>   i_step     integer                    step
  !>   axis       integer                    control component (0: no component rows)
  !>   other      integer                    another component with rows
  !>   other_txt  character(:), allocatable  its component as text
  !>   comp       character(:), allocatable  the loading component as text ("33"), or ""
  !>   rate_recs  integer(3)                 records of strain_rate, load_rate, stress_rate on the
  !>                                         loading component (0: none), in rate_by_* order
  !>   target_recs  integer(5)               records of target_strain, target_load,
  !>                                         target_stress (on the loading component),
  !>                                         target_time and target_duration (0: none), in
  !>                                         target_by_* order
  !>   given      logical(5)                 the step gives each of the five targets
  !>   given_rate logical(3)                 the step gives each of the three rates
  !>   names      character(:), allocatable  the five target rows as text, for messages
  !>   rate_names character(32)(3)           the three rate rows as text, for messages
  !>   k, first   integer                    target kind, the first target found
  !>   kr, first_rate  integer               rate kind, the first rate found
  !>   rate       real(rk)                   the control component's load or stress rate
  !>   open_end   logical                    the step's duration is not known (a load or stress
  !>                                         target at a strain rate, a strain target at a load
  !>                                         or stress rate)
  !>   step       character(:), allocatable  "step <i>" with the component, for messages
  !>   rows       character(13)(6)           the rows of another component (strain_rate,
  !>                                         target_strain, target_load, target_stress,
  !>                                         load_rate, stress_rate)
  !>   given_other  logical(6)               the step gives each of them on that component
  !>   k1, k2     integer                    the first two of them given
  subroutine check_step_controls(config, store, is_preset)

    type(config_t), intent(inout) :: config
    type(cfg_store_t), intent(in) :: store
    logical, intent(in) :: is_preset

    integer :: i_step, axis, rate_recs(3), target_recs(5), other, k, first, k1, k2, kr, &
      & first_rate
    logical :: given(5), given_rate(3), given_other(6), open_end
    real(rk) :: rate
    character(len=:), allocatable :: comp, step, other_txt
    character(len=32) :: names(5), rate_names(3)
    character(len=13), parameter :: rows(6) = [character(len=13) :: 'strain_rate', &
      & 'target_strain', 'target_load', 'target_stress', 'load_rate', 'stress_rate']

    ! The rows on the loading component (none if no component row was given).
    axis = config%control_axis
    comp = ''
    rate_recs = 0
    target_recs = 0
    if (axis > 0) then
      comp = int_text(11*axis)
      rate_recs(rate_by_strain) = cfg_find(store, 'strain_rate', 0, 11*axis)
      rate_recs(rate_by_load) = cfg_find(store, 'load_rate', 0, 11*axis)
      rate_recs(rate_by_stress) = cfg_find(store, 'stress_rate', 0, 11*axis)
      target_recs(target_by_strain) = cfg_find(store, 'target_strain', 0, 11*axis)
      target_recs(target_by_load) = cfg_find(store, 'target_load', 0, 11*axis)
      target_recs(target_by_stress) = cfg_find(store, 'target_stress', 0, 11*axis)
      names(target_by_strain) = 'target_strain ' // comp
      names(target_by_load) = 'target_load ' // comp
      names(target_by_stress) = 'target_stress ' // comp
      rate_names = [character(len=32) :: 'strain_rate ' // comp, 'load_rate ' // comp, &
        & 'stress_rate ' // comp]
    else
      names(1:3) = [character(len=32) :: 'target_strain <component>', &
        & 'target_load <component>', 'target_stress <component>']
      rate_names = [character(len=32) :: 'strain_rate', 'load_rate', 'stress_rate']
    end if
    target_recs(target_by_time) = cfg_find(store, 'target_time', 0, 0)
    names(target_by_time) = 'target_time'
    target_recs(target_by_duration) = cfg_find(store, 'target_duration', 0, 0)
    names(target_by_duration) = 'target_duration'

    do i_step = 1, config%num_steps
      step = store%file_name // ': step ' // int_text(i_step)
      if (axis > 0) step = step // ', component ' // comp

      ! Exactly one target (step 8.1: any of them in any step; target_duration step 8.2c).
      given = .false.
      if (axis > 0) given(1:3) = [config%has_target_strain(axis, i_step), &
        & config%has_target_load(axis, i_step), config%has_target_stress(axis, i_step)]
      given(target_by_time) = config%has_target_time(i_step)
      given(target_by_duration) = config%has_target_duration(i_step)
      first = 0
      do k = 1, 5
        if (.not. given(k)) cycle
        if (first > 0) call fatal_error('cfg', step // ': two targets (' // &
          & trim(names(first)) // ', line ' // line_text(store, target_recs(first)) // &
          & ', and ' // trim(names(k)) // ', line ' // line_text(store, target_recs(k)) // &
          & '); give one, "-" for the other')
        first = k
      end do
      if (first == 0) call fatal_error('cfg', step // ': no target (' // trim(names(1)) // &
        & ', ' // trim(names(2)) // ', ' // trim(names(3)) // ', target_time or ' // &
        & 'target_duration)')
      config%target_kind(i_step) = first

      ! Exactly one rate with a preset (one row per kind and component; step 8.2c: strain, load
      ! or stress rate).
      config%rate_kind(i_step) = rate_by_strain
      if (is_preset) then
        if (axis == 0) call fatal_error('cfg', step // ': no rate (boundary_conditions ' // &
          & trim(config%bc_type) // ' needs a strain_rate row, e.g. strain_rate 33 1.0e-3)')
        given_rate = [config%has_strain_rate(axis, i_step), config%has_load_rate(axis, i_step), &
          & config%has_stress_rate(axis, i_step)]
        first_rate = 0
        do kr = 1, 3
          if (.not. given_rate(kr)) cycle
          if (first_rate > 0) call fatal_error('cfg', step // ': two rates (' // &
            & trim(rate_names(first_rate)) // ', line ' // &
            & line_text(store, rate_recs(first_rate)) // ', and ' // trim(rate_names(kr)) // &
            & ', line ' // line_text(store, rate_recs(kr)) // '); give one, "-" for the other')
          first_rate = kr
        end do
        if (first_rate == 0) then
          if (any(rate_recs > 0)) then
            kr = findloc(rate_recs > 0, .true., dim=1)
            call fatal_error('cfg', step // ': no rate (' // trim(rate_names(kr)) // &
              & ', line ' // line_text(store, rate_recs(kr)) // ', has "-")')
          else
            call fatal_error('cfg', step // ': no rate (boundary_conditions ' // &
              & trim(config%bc_type) // ' needs a strain_rate ' // comp // ' row)')
          end if
        end if
        config%rate_kind(i_step) = first_rate
      end if

      ! A load or stress rate of the control component (step 8.2c, addendum §3.1, §3.2): it
      ! ends at a target of its kind, a strain target or a time; a rate 0 holds the value (a
      ! dwell) and needs a time, and cannot start the run.
      if (config%rate_kind(i_step) /= rate_by_strain) then
        kr = config%rate_kind(i_step)
        if (kr == rate_by_load) then
          rate = config%load_rate(axis, i_step)
        else
          rate = config%stress_rate(axis, i_step)
        end if
        if (kr == rate_by_load .and. first == target_by_stress) call fatal_error('cfg', &
          & cfg_rule_message(rule_load_rate_target, trim(rate_names(kr)) // ' (line ' // &
          & line_text(store, rate_recs(kr)) // ')', prefix=step, detail=trim(names(first)) // &
          & ', line ' // line_text(store, target_recs(first))))
        if (kr == rate_by_stress .and. first == target_by_load) call fatal_error('cfg', &
          & cfg_rule_message(rule_stress_rate_target, trim(rate_names(kr)) // ' (line ' // &
          & line_text(store, rate_recs(kr)) // ')', prefix=step, detail=trim(names(first)) // &
          & ', line ' // line_text(store, target_recs(first))))
        if (.not. abs(rate) > 0.0_rk) then
          if (first /= target_by_time .and. first /= target_by_duration) call fatal_error( &
            & 'cfg', step // ': ' // trim(rate_names(kr)) // ' 0 (line ' // &
            & line_text(store, rate_recs(kr)) // ') never reaches ' // trim(names(first)) // &
            & ': a rate 0 holds the value (a dwell); give target_time or target_duration')
          ! (After a restart point, step 1 follows the restart point's step: allowed, step 8.4.)
          if (i_step == 1 .and. .not. config%restart) call fatal_error('cfg', step // &
            & ': a dwell (' // &
            & trim(rate_names(kr)) // ' 0, line ' // line_text(store, rate_recs(kr)) // &
            & ') cannot be the first step (the load starts at 0); load first')
        end if
      end if

      ! The step's duration is not known: a load or stress target at a strain rate (step 8.1),
      ! or a strain target at a load or stress rate (step 8.2c).
      open_end = is_load_step(config, i_step)

      ! Every other component with rows: one rate or target (D45; a load or stress target or
      ! rate: triaxial load control, steps 8.2b, 8.2c).
      do other = 1, 3
        if (other == axis .or. .not. config%comp_rows(other)) cycle
        other_txt = int_text(11*other)
        step = store%file_name // ': step ' // int_text(i_step) // ', component ' // other_txt
        given_other = [config%has_strain_rate(other, i_step), &
          & config%has_target_strain(other, i_step), config%has_target_load(other, i_step), &
          & config%has_target_stress(other, i_step), config%has_load_rate(other, i_step), &
          & config%has_stress_rate(other, i_step)]
        if (count(given_other) > 1) then
          k1 = findloc(given_other, .true., dim=1)
          k2 = k1 + findloc(given_other(k1 + 1:), .true., dim=1)
          call fatal_error('cfg', step // ': ' // trim(rows(k1)) // ' ' // other_txt // &
            & ' (line ' // line_text(store, cfg_find(store, trim(rows(k1)), 0, 11*other)) // &
            & ') and ' // trim(rows(k2)) // ' ' // other_txt // ' (line ' // &
            & line_text(store, cfg_find(store, trim(rows(k2)), 0, 11*other)) // ') both ' // &
            & 'given; a component other than the control component ' // comp // ' takes one, ' &
            & // '"-" for the others')
        end if
        if (.not. any(given_other)) call fatal_error('cfg', step // ': no strain_rate, ' // &
          & 'target_strain, target_load, target_stress, load_rate or stress_rate (component ' &
          & // other_txt // ' has rows: give one of them in every step)')
        ! The boundary-velocity iteration (A) perturbs the load-controlled face velocities by
        ! min_pert_frac times the control face velocity (v2): the control face must move.
        if (any(given_other(3:6)) .and. is_preset) then
          if (config%has_strain_rate(axis, i_step) .and. &
            & .not. abs(config%strain_rate(axis, i_step)) > 0.0_rk) then
            k = findloc(given_other(3:6), .true., dim=1) + 2
            call fatal_error('cfg', step // ': ' // trim(rows(k)) // ' ' // other_txt // &
              & ' with strain_rate ' // comp // ' 0: the load-controlled face velocities ' // &
              & 'are found by perturbations of min_pert_frac times the control face ' // &
              & 'velocity, which must not be 0')
          end if
        end if
        ! A target over a duration that is not known.
        if (open_end .and. config%has_target_strain(other, i_step)) call fatal_error('cfg', &
          & step // ': target_strain ' // other_txt // ' (line ' // line_text(store, &
          & cfg_find(store, 'target_strain', 0, 11*other)) // ') in a step that ends at ' // &
          & trim(names(first)) // ': the duration of the step is not known; give ' // &
          & 'strain_rate ' // other_txt)
        if (open_end .and. config%rate_kind(i_step) /= rate_by_strain .and. &
          & any(given_other(3:4))) then
          k = findloc(given_other(3:4), .true., dim=1) + 2
          call fatal_error('cfg', step // ': ' // trim(rows(k)) // ' ' // other_txt // &
            & ' (line ' // line_text(store, cfg_find(store, trim(rows(k)), 0, 11*other)) // &
            & ') in a step that ends at ' // trim(names(first)) // ' at a ' // &
            & trim(rate_names(config%rate_kind(i_step))) // ': the duration of the step is ' // &
            & 'not known; give load_rate ' // other_txt // ' or stress_rate ' // other_txt)
        end if
      end do
    end do

  end subroutine check_step_controls

  !> Increments of each step (dev/cfg_keywords.md §1a rule 6): exactly one of
  !> number_of_incrs, dtime and dstrain per step ("-" where another one is used).
  !>
  !> Arguments:
  !>   config  (inout) type(config_t)     the configuration
  !>   store   (in)    type(cfg_store_t)  validated records
  !> Locals:
  !>   rows    character(15)(3)           the three keywords
  !>   recs    integer(3)                 their records (0: absent)
  !>   i_step, i_row, first  integer      step, row, the row already found for the step
  !>   rec_min  integer                   dtime_min record (0: absent)
  !>   has_min  logical                   dtime_min gives a value for the step
  !>   text     character(:), allocatable message part
  subroutine build_increments(config, store)

    type(config_t), intent(inout) :: config
    type(cfg_store_t), intent(in) :: store

    character(len=15), parameter :: rows(3) = [character(len=15) :: 'number_of_incrs', &
      & 'dtime', 'dstrain']
    integer :: recs(3), i_step, i_row, first, rec_min
    logical :: has_min
    character(len=:), allocatable :: text

    allocate (config%incr_kind(config%num_steps), config%num_incrs(config%num_steps))
    allocate (config%dtime(config%num_steps), config%dstrain(config%num_steps))
    config%incr_kind = 0
    config%num_incrs = 0
    config%dtime = 0.0_rk
    config%dstrain = 0.0_rk
    recs = [(cfg_find(store, trim(rows(i_row)), 0, 0), i_row = 1, 3)]

    do i_step = 1, config%num_steps
      first = 0
      do i_row = 1, 3
        ! The row is absent or has "-" for this step.
        if (recs(i_row) == 0) cycle
        if (cfg_is_dash(store, recs(i_row), i_step)) cycle

        ! A second row for the same step.
        if (first > 0) call cfg_error(store, recs(i_row), 'step ' // int_text(i_step) // &
          & ': ' // trim(rows(first)) // ' (line ' // line_text(store, recs(first)) // &
          & ') and ' // trim(rows(i_row)) // ' both give the increments; give one, "-" ' // &
          & 'for the other')
        first = i_row

        ! The value (incr_by_count, incr_by_dtime, incr_by_dstrain follow the row order).
        config%incr_kind(i_step) = i_row
        select case (i_row)
        case (incr_by_count)
          config%num_incrs(i_step) = cfg_get_int(store, recs(i_row), i_step)
        case (incr_by_dtime)
          config%dtime(i_step) = cfg_get_real(store, recs(i_row), i_step)
        case default
          config%dstrain(i_step) = cfg_get_real(store, recs(i_row), i_step)
        end select
      end do
      if (first == 0) call fatal_error('cfg', store%file_name // ': step ' // &
        & int_text(i_step) // ': no increments (number_of_incrs, dtime or dstrain)')

      ! A load- or stress-target step has no known duration: its increments are bounded by
      ! dtime (or dstrain) and dtime_min (step 8.1, addendum §3.4; a strain target at a load or
      ! stress rate, step 8.2c).
      if (is_load_step(config, i_step) .and. config%incr_kind(i_step) == incr_by_count) &
        & call cfg_error(store, recs(incr_by_count), &
        & cfg_rule_message(rule_number_of_incrs_open, 'number_of_incrs', prefix='step ' // &
        & int_text(i_step), detail='this step: a ' // open_target_text(config, i_step)))

      ! At a load or stress rate of the control component (step 8.2c): dstrain needs a strain
      ! rate; a dwell (rate 0) takes dtime (v2's target time increment, addendum §3.3).
      if (config%rate_kind(i_step) /= rate_by_strain) then
        if (config%incr_kind(i_step) == incr_by_dstrain) call cfg_error(store, &
          & recs(incr_by_dstrain), cfg_rule_message(rule_dstrain_rate, 'dstrain', &
          & prefix='step ' // int_text(i_step)))
        if (is_dwell_step(config, i_step) .and. config%incr_kind(i_step) /= incr_by_dtime) &
          & call cfg_error(store, recs(config%incr_kind(i_step)), &
          & cfg_rule_message(rule_number_of_incrs_dwell, trim(rows(config%incr_kind(i_step))), &
          & prefix='step ' // int_text(i_step)))
      end if
    end do

    ! dtime_min: required in every load- or stress-target step (user decision Q-L4), "-"
    ! elsewhere; not larger than dtime.
    allocate (config%dtime_min(config%num_steps))
    config%dtime_min = 0.0_rk
    rec_min = cfg_find(store, 'dtime_min', 0, 0)
    do i_step = 1, config%num_steps
      has_min = .false.
      if (rec_min > 0) has_min = .not. cfg_is_dash(store, rec_min, i_step)
      if (is_load_step(config, i_step)) then
        if (.not. has_min) then
          text = ''
          if (rec_min > 0) text = ' (dtime_min, line ' // line_text(store, rec_min) // &
            & ', has "-")'
          call fatal_error('cfg', store%file_name // ': step ' // int_text(i_step) // &
            & ': a ' // open_target_text(config, i_step) // ' needs ' // &
            & 'dtime_min, the smallest time increment of the step' // text)
        end if
        config%dtime_min(i_step) = cfg_get_real(store, rec_min, i_step)
        if (config%incr_kind(i_step) == incr_by_dtime) then
          if (config%dtime_min(i_step) > config%dtime(i_step)) call cfg_error(store, &
            & rec_min, 'step ' // int_text(i_step) // ': dtime_min ' // &
            & real_text(config%dtime_min(i_step)) // ' s is larger than dtime ' // &
            & real_text(config%dtime(i_step)) // ' s')
        end if
      else if (has_min) then
        call cfg_error(store, rec_min, cfg_rule_message(rule_dtime_min, 'dtime_min', &
          & prefix='step ' // int_text(i_step)))
      end if
    end do

  end subroutine build_increments

  !> Options of load and stress targets (step 8.1, addendum_load_targets.md §3.4, §3.6):
  !> dtime_factor (default 1.001, v2) only with a load- or stress-target step; check_necking
  !> on|off (default off), which needs a loaded component. Step 8.2a
  !> (addendum_triaxial_load_control.md §3.6, §4, decisions Q-T13, Q-T18): load_tol [N] only
  !> with a target_load row, stress_tol [MPa] only with a target_stress row (absolute, > 0,
  !> default 0.1), max_bc_iter (>= 1, default 10, v2) only with a load or stress target. Step
  !> 8.2b (addendum §3.4, §4): the load and stress rows of every component count;
  !> min_pert_frac (> 0, default 0.001, v2) only with a load or stress row on a component other
  !> than the control component. Step 8.2c: load and stress rates count as load and stress
  !> rows (also for min_pert_frac, iteration (B)'s lower bound, Q-T1); max_strain_incr (> 0,
  !> default 0.001, v2) only with a dwell.
  !>
  !> Arguments:
  !>   config  (inout) type(config_t)     the configuration (steps and targets read)
  !>   store   (in)    type(cfg_store_t)  validated records
  !> Locals:
  !>   i_rec   integer                    a record (0: absent)
  !>   any_load, any_stress  logical      some step has a load, a stress target or row
  !>   any_other  logical                 some step has a load or stress row on a component
  !>                                      other than the control component, or a load or
  !>                                      stress rate of the control component
  !>   any_dwell  logical                 some step is a dwell
  !>   other, i_step  integer             a component, a step
  subroutine build_load_options(config, store)

    type(config_t), intent(inout) :: config
    type(cfg_store_t), intent(in) :: store

    integer :: i_rec, other, i_step
    logical :: any_load, any_stress, any_other, any_dwell

    any_load = any(config%target_kind == target_by_load) .or. any(config%has_target_load) .or. &
      & any(config%has_load_rate)
    any_stress = any(config%target_kind == target_by_stress) .or. &
      & any(config%has_target_stress) .or. any(config%has_stress_rate)
    any_other = any(config%rate_kind /= rate_by_strain)
    do other = 1, 3
      if (other == config%control_axis) cycle
      any_other = any_other .or. any(config%has_target_load(other, :)) .or. &
        & any(config%has_target_stress(other, :)) .or. any(config%has_load_rate(other, :)) .or. &
        & any(config%has_stress_rate(other, :))
    end do
    any_dwell = .false.
    do i_step = 1, config%num_steps
      any_dwell = any_dwell .or. is_dwell_step(config, i_step)
    end do

    ! Time-step factor near a target.
    config%dtime_factor = real_or_default(store, 'dtime_factor')
    i_rec = cfg_find(store, 'dtime_factor', 0, 0)
    if (i_rec > 0 .and. .not. (any_load .or. any_stress)) call cfg_error(store, i_rec, &
      & cfg_rule_message(rule_dtime_factor, 'dtime_factor'))

    ! Tolerances of the load and stress rows (Q-T18), and the landing attempts (Q-T13).
    config%load_tol = real_or_default(store, 'load_tol')
    i_rec = cfg_find(store, 'load_tol', 0, 0)
    if (i_rec > 0 .and. .not. any_load) call cfg_error(store, i_rec, &
      & cfg_rule_message(rule_load_tol, 'load_tol'))
    config%stress_tol = real_or_default(store, 'stress_tol')
    i_rec = cfg_find(store, 'stress_tol', 0, 0)
    if (i_rec > 0 .and. .not. any_stress) call cfg_error(store, i_rec, &
      & cfg_rule_message(rule_stress_tol, 'stress_tol'))
    config%max_bc_iter = nint(real_or_default(store, 'max_bc_iter'))
    i_rec = cfg_find(store, 'max_bc_iter', 0, 0)
    if (i_rec > 0 .and. .not. (any_load .or. any_stress)) call cfg_error(store, i_rec, &
      & cfg_rule_message(rule_max_bc_iter, 'max_bc_iter'))

    ! Lower bound of the velocity perturbations of the boundary-velocity iteration (step 8.2b).
    config%min_pert_frac = real_or_default(store, 'min_pert_frac')
    i_rec = cfg_find(store, 'min_pert_frac', 0, 0)
    if (i_rec > 0 .and. .not. any_other) call cfg_error(store, i_rec, &
      & cfg_rule_message(rule_min_pert_frac, 'min_pert_frac'))

    ! Dwell increments (step 8.2c, v2).
    config%max_strain_incr = real_or_default(store, 'max_strain_incr')
    i_rec = cfg_find(store, 'max_strain_incr', 0, 0)
    if (i_rec > 0 .and. .not. any_dwell) call cfg_error(store, i_rec, &
      & cfg_rule_message(rule_max_strain_incr, 'max_strain_incr'))

    ! Necking: the load of the control component's face x1/y1/z1.
    i_rec = cfg_find(store, 'check_necking', 0, 0)
    if (i_rec > 0) then
      config%check_necking = cfg_get_word(store, i_rec, 1) == 'on'
      if (config%check_necking .and. config%bc_type == 'periodic') call cfg_error(store, &
        & i_rec, cfg_rule_message(rule_necking_periodic, 'check_necking on'))
      if (config%check_necking .and. config%control_axis == 0) call cfg_error(store, i_rec, &
        & cfg_rule_message(rule_necking_loaded, 'check_necking on'))
    end if

  end subroutine build_load_options

  !> Whether a step's duration is not known: it ends at a load or stress target at a strain
  !> rate (step 8.1), or at a strain target at a load or stress rate (step 8.2c, Q-T8 (i)). Such
  !> a step takes dtime or dstrain (the largest increment) and dtime_min.
  !>
  !> Arguments:
  !>   config  (in)  type(config_t)  the configuration (target_kind and rate_kind set)
  !>   i_step  (in)  integer         step
  pure logical function is_load_step(config, i_step)

    type(config_t), intent(in) :: config
    integer, intent(in) :: i_step

    if (config%rate_kind(i_step) == rate_by_strain) then
      is_load_step = config%target_kind(i_step) == target_by_load .or. &
        & config%target_kind(i_step) == target_by_stress
    else
      is_load_step = config%target_kind(i_step) == target_by_strain
    end if

  end function is_load_step

  !> Whether a step is a dwell: a load or stress rate 0 of the control component (step 8.2c).
  !>
  !> Arguments:
  !>   config  (in)  type(config_t)  the configuration (rates and rate_kind set)
  !>   i_step  (in)  integer         step
  pure logical function is_dwell_step(config, i_step)

    type(config_t), intent(in) :: config
    integer, intent(in) :: i_step

    select case (config%rate_kind(i_step))
    case (rate_by_load)
      is_dwell_step = .not. abs(config%load_rate(config%control_axis, i_step)) > 0.0_rk
    case (rate_by_stress)
      is_dwell_step = .not. abs(config%stress_rate(config%control_axis, i_step)) > 0.0_rk
    case default
      is_dwell_step = .false.
    end select

  end function is_dwell_step

  !> The target of a step whose duration is not known, for messages: "load target", "stress
  !> target", or "strain target at a load rate" ("... stress rate").
  !>
  !> Arguments:
  !>   config  (in)  type(config_t)  the configuration (target_kind and rate_kind set)
  !>   i_step  (in)  integer         step
  !>   text    (out, result) character(:), allocatable
  pure function open_target_text(config, i_step) result(text)

    type(config_t), intent(in) :: config
    integer, intent(in) :: i_step
    character(len=:), allocatable :: text

    text = trim(target_word(config%target_kind(i_step))) // ' target'
    if (config%rate_kind(i_step) == rate_by_load) text = text // ' at a load rate'
    if (config%rate_kind(i_step) == rate_by_stress) text = text // ' at a stress rate'

  end function open_target_text

  !> The word of a target kind, for messages ("strain", "load", "stress", "time").
  !>
  !> Arguments:
  !>   kind  (in)  integer  target_by_* code
  pure function target_word(kind) result(word)

    integer, intent(in) :: kind
    character(len=8) :: word

    select case (kind)
    case (target_by_strain)
      word = 'strain'
    case (target_by_load)
      word = 'load'
    case (target_by_stress)
      word = 'stress'
    case (target_by_duration)
      word = 'duration'
    case default
      word = 'time'
    end select

  end function target_word

  !> The set_bc lines: "set_bc vel <nset> <x|y|z> <v> [<x|y|z> <v> [<x|y|z> <v>]]"; one
  !> velocity per node set and direction.
  !>
  !> Arguments:
  !>   config  (inout) type(config_t)     the configuration
  !>   store   (in)    type(cfg_store_t)  validated records
  !> Locals:
  !>   i_rec, i_pair, i, ios  integer     record, pair and condition indices, read status
  !>   nvals   integer                    values of the line
  !>   bc      type(cfg_vel_bc_t)         one condition
  !>   word    character(:), allocatable  a value
  !>   grown   type(cfg_vel_bc_t)(:)      larger copy of the list
  subroutine build_vel_bcs(config, store)

    type(config_t), intent(inout) :: config
    type(cfg_store_t), intent(in) :: store

    integer :: i_rec, i_pair, i, ios, nvals
    type(cfg_vel_bc_t) :: bc
    character(len=:), allocatable :: word
    type(cfg_vel_bc_t), allocatable :: grown(:)

    allocate (config%vel_bcs(0))
    do i_rec = 1, store%num_records
      if (cfg_keywords(store%records(i_rec)%key)%name /= 'set_bc') cycle

      ! Kind of condition.
      nvals = size(store%records(i_rec)%values)
      word = cfg_get_word(store, i_rec, 1)
      select case (word)
      case ('vel')
      case default
        call cfg_error(store, i_rec, 'set_bc: invalid kind "' // word // '" (expected vel)')
      end select
      if (nvals < 4 .or. nvals > 8 .or. mod(nvals, 2) /= 0) call cfg_error(store, i_rec, &
        & 'set_bc vel: expected <nset> and one to three <x|y|z> <velocity> pairs')

      ! Each direction-value pair.
      bc%nset = cfg_get_word(store, i_rec, 2)
      bc%line = store%records(i_rec)%line
      do i_pair = 1, (nvals - 2)/2
        word = cfg_get_word(store, i_rec, 1 + 2*i_pair)
        bc%comp = findloc(axis_letter, word, dim=1)
        if (bc%comp == 0) call cfg_error(store, i_rec, 'set_bc vel: invalid direction "' // &
          & word // '" (expected x, y or z)')
        word = cfg_get_word(store, i_rec, 2 + 2*i_pair)
        read (word, *, iostat=ios) bc%value
        if (ios /= 0 .or. verify(word, '+-.0123456789eEdD') /= 0) call cfg_error(store, i_rec, &
          & 'set_bc vel: "' // word // '" is not a number')
        do i = 1, size(config%vel_bcs)
          if (config%vel_bcs(i)%nset == bc%nset .and. config%vel_bcs(i)%comp == bc%comp) then
            call cfg_error(store, i_rec, 'set_bc vel: ' // bc%nset // ' ' // &
              & axis_letter(bc%comp) // ' given twice (lines ' // &
              & int_text(int(config%vel_bcs(i)%line)) // ' and ' // &
              & int_text(int(bc%line)) // ')')
          end if
        end do

        ! Append (element by element: derived-type array constructors with allocatable
        ! components lost data with gfortran 13 -O3, ROADMAP section G).
        allocate (grown(size(config%vel_bcs) + 1))
        do i = 1, size(config%vel_bcs)
          grown(i) = config%vel_bcs(i)
        end do
        grown(size(grown)) = bc
        call move_alloc(grown, config%vel_bcs)
      end do
    end do

  end subroutine build_vel_bcs

  !> Record of a required global keyword (a missing keyword is a fatal error).
  !>
  !> Arguments:
  !>   store  (in)  type(cfg_store_t)  validated records
  !>   name   (in)  character          keyword
  function require(store, name) result(i_rec)

    type(cfg_store_t), intent(in) :: store
    character(len=*), intent(in) :: name
    integer :: i_rec

    i_rec = cfg_find(store, name, 0, 0)
    if (i_rec == 0) call fatal_error('cfg', store%file_name // ': ' // name // ' is required')

  end function require

  !> Record of a required per-phase keyword (a missing keyword is a fatal error).
  !>
  !> Arguments:
  !>   store      (in)  type(cfg_store_t)  validated records
  !>   name       (in)  character          keyword
  !>   i_phase    (in)  integer            phase
  !>   phase_rec  (in)  integer            the record of "phase i_phase"
  function require_phase(store, name, i_phase, phase_rec) result(i_rec)

    type(cfg_store_t), intent(in) :: store
    character(len=*), intent(in) :: name
    integer, intent(in) :: i_phase, phase_rec
    integer :: i_rec

    i_rec = cfg_find(store, name, i_phase, 0)
    if (i_rec == 0) call fatal_error('cfg', store%file_name // ': phase ' // int_text(i_phase) &
      & // ' (line ' // line_text(store, phase_rec) // '): ' // name // ' is required')

  end function require_phase

  !> Value of an optional numeric keyword, or its default from the keyword table.
  !>
  !> Arguments:
  !>   store  (in)  type(cfg_store_t)  validated records
  !>   name   (in)  character          keyword
  !> Locals:
  !>   i_rec  integer                  record index
  function real_or_default(store, name) result(x)

    type(cfg_store_t), intent(in) :: store
    character(len=*), intent(in) :: name
    real(rk) :: x

    integer :: i_rec

    i_rec = cfg_find(store, name, 0, 0)
    if (i_rec > 0) then
      x = cfg_get_real(store, i_rec, 1)
    else
      read (cfg_keywords(cfg_find_keyword(name))%default_text, *) x
    end if

  end function real_or_default

  !> Value of an optional word keyword, or its default from the keyword table.
  !>
  !> Arguments:
  !>   store  (in)  type(cfg_store_t)  validated records
  !>   name   (in)  character          keyword
  !> Locals:
  !>   i_rec  integer                  record index
  function word_or_default(store, name) result(word)

    type(cfg_store_t), intent(in) :: store
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: word

    integer :: i_rec

    i_rec = cfg_find(store, name, 0, 0)
    if (i_rec > 0) then
      word = cfg_get_word(store, i_rec, 1)
    else
      word = trim(cfg_keywords(cfg_find_keyword(name))%default_text)
    end if

  end function word_or_default

  !> Line number of a record, as text.
  !>
  !> Arguments:
  !>   store  (in)  type(cfg_store_t)  validated records
  !>   i_rec  (in)  integer            record index
  function line_text(store, i_rec) result(text)

    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_rec
    character(len=:), allocatable :: text

    text = int_text(int(store%records(i_rec)%line))

  end function line_text

  !> Log the material part of the input summary (step 9.7, ROADMAP D94,
  !> dev/plan/terminal_output.md §3.2; Q61): "- Material parameters:" and "> Number of
  !> phases: <n>"; the phases follow with config_log_phase. Rank 0 writes (log_mod).
  !>
  !> Arguments:
  !>   config  (in)  type(config_t)  the configuration
  subroutine config_log_material(config)

    type(config_t), intent(in) :: config

    ! The heading and the number of phases.
    call log_info('Material parameters:', 1)
    call log_info('Number of phases: ' // int_text(config%num_phases), 2)

  end subroutine config_log_material

  !> Log one phase of the input summary (step 9.7; Q61): "> Phase <i>: <crystal_type>,
  !> <model>", then one ". <keyword>: <values> <unit>" line per material keyword, in the
  !> keyword table's order but with hardening before tau_0 (the hardening parameters under the
  !> hardening law, user 2026-10-26); a keyword not given is shown with its default and
  !> "(def.)" when it has a default value and applies to the phase (keyword_applies).
  !> Values are written as given, reals in ES format with 4 digits. Rank 0 writes.
  !>
  !> Arguments:
  !>   config   (in)  type(config_t)  the configuration (with its store, config_read)
  !>   i_phase  (in)  integer         the phase
  !> Locals:
  !>   k        integer               keyword entry
  subroutine config_log_phase(config, i_phase)

    type(config_t), intent(in) :: config
    integer, intent(in) :: i_phase

    integer :: k

    associate (ph => config%phases(i_phase))
      call log_info('Phase ' // int_text(i_phase) // ': ' // trim(ph%crystal_type) // ', ' // &
        & trim(ph%model), 2)
      ! The phase keywords in the table's order, hardening moved before tau_0.
      do k = 1, cfg_num_material
        if (cfg_keywords(k)%scope /= cfg_phase) cycle
        select case (trim(cfg_keywords(k)%name))
        case ('crystal_type', 'model', 'hardening')
          cycle
        case ('tau_0')
          call log_parameter(config%store, i_phase, ph, cfg_find_keyword('hardening'))
        end select
        call log_parameter(config%store, i_phase, ph, k)
      end do
    end associate

  end subroutine config_log_phase

  !> Log the options and the output of the input summary (step 9.7; Q61): "- Options other
  !> than the defaults:" with one "> <keyword> <values>" line per solver or run keyword given
  !> with a value other than its default (option_differs), or "(none)"; then "- Output: <print
  !> variables> (format <output_format>)". Rank 0 writes.
  !>
  !> Arguments:
  !>   config  (in)  type(config_t)  the configuration (with its store, config_read)
  !> Locals:
  !>   k, j     integer                    keyword entry, value index
  !>   i_rec    integer                    record index
  !>   num_opts integer                    options listed
  !>   vars, formats  character(:), allocatable  the print variables and the output formats
  !>   value    character(:), allocatable  an option's values as text
  subroutine config_log_options(config)

    type(config_t), intent(in) :: config

    integer :: k, j, i_rec, num_opts
    character(len=:), allocatable :: vars, formats, value

    ! --- Options other than the defaults (solver and run keywords) ---
    num_opts = 0
    do k = cfg_num_material + 1, size(cfg_keywords)
      associate (kw => cfg_keywords(k))
        if (kw%status /= cfg_active .or. kw%scope /= cfg_global) cycle
        if (kw%form /= cfg_single .and. kw%form /= cfg_flag) cycle
        ! (The loading and the output have their own lines.)
        select case (trim(kw%name))
        case ('boundary_conditions', 'control_component', 'number_of_steps')
          cycle
        end select
        i_rec = cfg_find(config%store, trim(kw%name), 0, 0)
        if (i_rec == 0) cycle
        if (.not. option_differs(config%store, i_rec, kw)) cycle
        num_opts = num_opts + 1
        if (num_opts == 1) call log_info('Options other than the defaults:', 1)
        value = record_values(config%store, i_rec, kw)
        if (len(value) > 0) then
          call log_info(trim(kw%name) // ' ' // value, 2)
        else
          call log_info(trim(kw%name), 2)
        end if
      end associate
    end do
    if (num_opts == 0) call log_info('Options other than the defaults: (none)', 1)

    ! --- Output: the print variables, in file order, and the formats ---
    vars = ''
    do i_rec = 1, config%store%num_records
      if (config%store%records(i_rec)%key /= cfg_find_keyword('print')) cycle
      do j = 1, size(config%store%records(i_rec)%values)
        vars = vars // ' ' // trim(config%store%records(i_rec)%values(j))
      end do
    end do
    if (len(vars) == 0) vars = ' (none)'
    i_rec = cfg_find(config%store, 'output_format', 0, 0)
    if (i_rec > 0) then
      formats = record_values(config%store, i_rec, cfg_keywords(cfg_find_keyword('output_format')))
    else
      formats = trim(cfg_keywords(cfg_find_keyword('output_format'))%default_text)
    end if
    call log_info('Output:' // vars // ' (format ' // formats // ')', 1)

  end subroutine config_log_options

  !> Log one material keyword of a phase (config_log_phase): ". <keyword>: <values> <unit>"
  !> for each record given (one per row for interaction_matrix), else its default with
  !> "(def.)" when it has a default value and applies to the phase; nothing otherwise. A
  !> line longer than the log's width continues under the values (hang).
  !>
  !> Arguments:
  !>   store    (in)  type(cfg_store_t)  validated records
  !>   i_phase  (in)  integer            phase
  !>   ph       (in)  type(cfg_phase_t)  the phase (keyword_applies)
  !>   k        (in)  integer            keyword entry
  !> Locals:
  !>   i_rec    integer                    record index
  !>   found    logical                    the keyword is given in the phase
  !>   label    character(27)              "<keyword>:" padded to align the values
  !>   unit     character(:), allocatable  " <unit>", or empty for "" and "-"
  subroutine log_parameter(store, i_phase, ph, k)

    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_phase
    type(cfg_phase_t), intent(in) :: ph
    integer, intent(in) :: k

    integer :: i_rec
    logical :: found
    character(len=27) :: label
    character(len=:), allocatable :: unit

    associate (kw => cfg_keywords(k))
      if (kw%status /= cfg_active) return
      label = trim(kw%name) // ':'
      unit = ''
      if (len_trim(kw%unit) > 0 .and. kw%unit /= '-') unit = ' ' // trim(kw%unit)
      ! The records given in the phase (several for a repeatable line).
      found = .false.
      do i_rec = 1, store%num_records
        if (store%records(i_rec)%key /= k .or. store%records(i_rec)%phase /= i_phase) cycle
        found = .true.
        call log_info(label // record_values(store, i_rec, kw) // unit, 3, hang=len(label))
      end do
      if (found) return
      ! Not given: the default, when there is one and it applies.
      if (.not. has_default_value(kw)) return
      if (.not. keyword_applies(trim(kw%name), ph)) return
      call log_info(label // trim(kw%default_text) // ' (def.)', 3, hang=len(label))
    end associate

  end subroutine log_parameter

  !> Whether a material keyword with a default applies to a phase (config_log_phase): the
  !> hardening law to crystal plasticity, interaction_matrix to saturation_anisotropic, the
  !> modifiers (saturation_evolution, cyclic_hardening, precipitate_strengthening) to the
  !> saturation laws (build_modifiers), irradiation to kocks_mecking (D61).
  !>
  !> Arguments:
  !>   name  (in)  character          the keyword
  !>   ph    (in)  type(cfg_phase_t)  the phase
  !> Locals:
  !>   saturation  logical  the phase has a saturation hardening law
  logical function keyword_applies(name, ph) result(applies)

    character(len=*), intent(in) :: name
    type(cfg_phase_t), intent(in) :: ph

    logical :: saturation

    ! Crystal plasticity first; then by the hardening law.
    applies = .false.
    if (ph%model /= 'crystal_plasticity') return
    saturation = ph%hardening == 'saturation_isotropic' .or. &
      & ph%hardening == 'saturation_anisotropic'
    select case (name)
    case ('hardening')
      applies = .true.
    case ('interaction_matrix')
      applies = ph%hardening == 'saturation_anisotropic'
    case ('saturation_evolution', 'cyclic_hardening', 'precipitate_strengthening')
      applies = saturation
    case ('irradiation')
      applies = ph%hardening == 'kocks_mecking'
    end select

  end function keyword_applies

  !> Whether a keyword's default text is a value (not empty, not "required ...", not a rule
  !> such as "per crystal_type" or "one of 3 per step").
  !>
  !> Arguments:
  !>   kw  (in)  type(cfg_keyword_t)  the keyword entry
  logical function has_default_value(kw) result(has)

    type(cfg_keyword_t), intent(in) :: kw

    ! A value: not empty and not one of the rule texts.
    has = len_trim(kw%default_text) > 0
    if (index(kw%default_text, 'required') == 1) has = .false.
    if (index(kw%default_text, 'per ') == 1) has = .false.
    if (index(kw%default_text, 'one of') == 1) has = .false.

  end function has_default_value

  !> Whether a given option differs from its default (config_log_options): numbers compared as
  !> values (relative difference above 1e-12), words as words; a keyword without a default value
  !> (a flag, a run limit not set by default) always differs.
  !>
  !> Arguments:
  !>   store  (in)  type(cfg_store_t)    validated records
  !>   i_rec  (in)  integer              the option's record
  !>   kw     (in)  type(cfg_keyword_t)  its keyword entry
  !> Locals:
  !>   given, default  real(rk)  the values (numbers)
  logical function option_differs(store, i_rec, kw) result(differs)

    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_rec
    type(cfg_keyword_t), intent(in) :: kw

    real(rk) :: given, default

    ! No default value: listed whenever given.
    differs = .true.
    if (.not. has_default_value(kw) .or. kw%form == cfg_flag) return
    select case (kw%kind)
    case (cfg_int, cfg_real)
      given = cfg_get_real(store, i_rec, 1)
      read (kw%default_text, *) default
      differs = abs(given - default) > 1.0e-12_rk*max(abs(given), abs(default))
    case (cfg_word)
      differs = cfg_get_word(store, i_rec, 1) /= trim(kw%default_text)
    end select

  end function option_differs

  !> The values of a record as text, separated by blanks: reals in ES format with 4 digits,
  !> the others as given (words in lower case, cfg_store_read).
  !>
  !> Arguments:
  !>   store  (in)  type(cfg_store_t)    validated records
  !>   i_rec  (in)  integer              record index
  !>   kw     (in)  type(cfg_keyword_t)  its keyword entry
  !> Result:
  !>   text   character(:), allocatable  the values
  !> Locals:
  !>   i      integer  value index
  function record_values(store, i_rec, kw) result(text)

    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_rec
    type(cfg_keyword_t), intent(in) :: kw
    character(len=:), allocatable :: text

    integer :: i

    text = ''
    if (.not. allocated(store%records(i_rec)%values)) return
    do i = 1, size(store%records(i_rec)%values)
      if (i > 1) text = text // ' '
      if (kw%kind == cfg_real .and. .not. cfg_is_dash(store, i_rec, i)) then
        text = text // log_real_text(cfg_get_real(store, i_rec, i), 4)
      else
        text = text // trim(store%records(i_rec)%values(i))
      end if
    end do

  end function record_values

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

end module config_mod
