! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The configuration keyword table (dev/cfg_keywords.md, ROADMAP D42): one entry per keyword,
!> with its status, scope, form, value kind, range or allowed words, unit, default and
!> description. The parser (cfg_store_mod) validates every line against its entry; the builder
!> (config_mod) uses the validated values. Adding or renaming a keyword is one entry.
!>
!> Source : FEPX (ROADMAP step 4.3); keyword names and meanings from v2-fixed
!>          read_input_cfg.f90 and its manual (doc/*.rst) as listed in dev/cfg_keywords.md.
!>          Contains no numerical method.
!>
!> Entry fields (cfg_keyword_t):
!>   name     keyword
!>   status   cfg_active (implemented), cfg_planned (recognized, implemented in the given
!>            stage: an error until then)
!>   scope    cfg_global, or cfg_phase (inside a phase block)
!>   form     cfg_single (fixed number of values, nvals), cfg_per_step (one value per step,
!>            "-" allowed when step_dash), cfg_comp_per_step (a component, then one value per
!>            step), cfg_repeat (a line that may appear several times; values checked by the
!>            builder), cfg_flag (no value), cfg_list (one or more values, the count checked
!>            by the builder, e.g. one per slip family or per slip system)
!>   kind     cfg_int, cfg_real, cfg_word (one of the allowed words), cfg_text (any token)
!>   nvals    number of values (cfg_single)
!>   lo, hi, lo_open   allowed range of numeric values (lo_open: lo excluded)
!>   words    allowed words (cfg_word), separated by blanks
!>   unit, default_text, info   documentation (default_text: "required", "required (<when>)",
!>            empty, or the default value)
!>   note     stage of a planned keyword
!>   later_words, later_note  planned values of an active cfg_word keyword (an error naming
!>            the stage in later_note until then)
!>   values_text  documentation of the values of a cfg_text or cfg_repeat keyword whose values
!>            the builder checks (their syntax, reStructuredText for the manual); empty when the
!>            kind, range and form describe them (cfg_manual_mod)
!>
!> Public parameters: the status, scope, form and kind codes; cfg_keywords (the table);
!>   cfg_num_material (its leading material entries).
!> Public routines:
!>   cfg_find_keyword    index of a keyword in the table (0 if unknown)
!>   cfg_table_markdown  write the table as a Markdown document (documentation generation; the
!>                       manual's pages are written by cfg_manual_mod)
module cfg_table_mod

  use kinds_mod, only: rk

  implicit none

  private

  public :: cfg_keyword_t, cfg_keywords, cfg_find_keyword, cfg_table_markdown
  public :: cfg_num_material
  public :: cfg_active, cfg_planned
  public :: cfg_global, cfg_phase
  public :: cfg_single, cfg_per_step, cfg_comp_per_step, cfg_repeat, cfg_flag, cfg_list
  public :: cfg_int, cfg_real, cfg_word, cfg_text

  ! Status, scope, form and kind codes.
  integer, parameter :: cfg_active = 1, cfg_planned = 2
  integer, parameter :: cfg_global = 1, cfg_phase = 2
  integer, parameter :: cfg_single = 1, cfg_per_step = 2, cfg_comp_per_step = 3, &
    & cfg_repeat = 4, cfg_flag = 5, cfg_list = 6
  integer, parameter :: cfg_int = 1, cfg_real = 2, cfg_word = 3, cfg_text = 4

  ! Bounds used in the table.
  real(rk), parameter :: big = huge(1.0_rk)

  !> One keyword (fields documented in the module header).
  type :: cfg_keyword_t
    character(len=32) :: name = ''
    integer :: status = cfg_active
    integer :: scope = cfg_global
    integer :: form = cfg_single
    integer :: kind = cfg_real
    integer :: nvals = 1
    real(rk) :: lo = -big
    real(rk) :: hi = big
    logical :: lo_open = .false.
    logical :: step_dash = .false.
    character(len=512) :: words = ''
    character(len=16) :: unit = ''
    character(len=48) :: default_text = ''
    character(len=96) :: info = ''
    character(len=96) :: note = ''
    character(len=512) :: later_words = ''
    character(len=64) :: later_note = ''
    character(len=240) :: values_text = ''
  end type cfg_keyword_t

  ! The table. Grouped as dev/cfg_keywords.md: material, boundary conditions and loading,
  ! solver, input files and output; planned keywords last.
  type(cfg_keyword_t), parameter :: cfg_keywords_material(*) = [ &
    ! --- Material (active) ---
    & cfg_keyword_t('number_of_phases', cfg_active, cfg_global, cfg_single, cfg_int, 1, 1.0_rk, &
    &   big, .false., .false., '', '', 'required', 'number of phases', ''), &
    & cfg_keyword_t('phase', cfg_active, cfg_global, cfg_single, cfg_int, 1, 1.0_rk, big, &
    &   .false., .false., '', '', '', 'start of the block of phase i', ''), &
    & cfg_keyword_t('crystal_type', cfg_active, cfg_phase, cfg_single, cfg_word, 1, -big, big, &
    &   .false., .false., 'fcc bcc hcp bct', '', 'required', 'crystal structure', ''), &
    & cfg_keyword_t('model', cfg_active, cfg_phase, cfg_single, cfg_word, 1, -big, big, &
    &   .false., .false., 'elastic crystal_plasticity', '', 'required', 'constitutive model', &
    &   ''), &
    & cfg_keyword_t('c11', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', 'MPa', 'required', 'elastic constant C11 (Voigt)', ''), &
    & cfg_keyword_t('c12', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', 'MPa', 'required', 'elastic constant C12 (Voigt)', ''), &
    & cfg_keyword_t('c13', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', 'MPa', 'required (hcp, bct)', &
    &   'elastic constant C13 (Voigt)', ''), &
    & cfg_keyword_t('c33', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', 'MPa', 'required (hcp, bct)', &
    &   'elastic constant C33 (Voigt; independent)', ''), &
    & cfg_keyword_t('c44', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', 'MPa', 'required', 'elastic constant C44 (Voigt)', ''), &
    & cfg_keyword_t('c66', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', 'MPa', 'required (bct)', 'elastic constant C66 (Voigt)', ''), &
    & cfg_keyword_t('c_over_a', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', '-', 'required (hcp, bct plasticity)', &
    &   'lattice ratio c/a (hcp, bct)', ''), &
    & cfg_keyword_t('slip_families', cfg_active, cfg_phase, cfg_list, cfg_text, 0, -big, big, &
    &   .false., .false., '', '', 'per crystal_type', &
    &   'active slip families, {plane}<direction>, or all (crystal_plasticity)', ''), &
    & cfg_keyword_t('m', cfg_active, cfg_phase, cfg_list, cfg_real, 0, 0.0_rk, big, .true., &
    &   .false., '', '-', 'required (crystal_plasticity)', &
    &   'rate sensitivity, one per active slip family or per slip system', ''), &
    & cfg_keyword_t('gammadot_0', cfg_active, cfg_phase, cfg_list, cfg_real, 0, 0.0_rk, big, &
    &   .true., .false., '', '1/s', 'required (crystal_plasticity)', &
    &   'reference slip rate, one per active slip family or per slip system', ''), &
    & cfg_keyword_t('tau_0', cfg_active, cfg_phase, cfg_list, cfg_real, 0, 0.0_rk, big, &
    &   .true., .false., '', 'MPa', 'required (crystal_plasticity)', &
    &   'initial slip strength, one per active slip family or per slip system', ''), &
    & cfg_keyword_t('hardening', cfg_active, cfg_phase, cfg_single, cfg_word, 1, -big, big, &
    &   .false., .false., 'saturation_isotropic saturation_anisotropic kocks_mecking none', &
    &   '', &
    &   'saturation_isotropic', 'slip-system hardening law (crystal_plasticity)', ''), &
    & cfg_keyword_t('h_0', cfg_active, cfg_phase, cfg_list, cfg_real, 0, 0.0_rk, big, .false., &
    &   .false., '', 'MPa', 'required (saturation laws)', &
    &   'hardening rate coefficient: one value, or per family or system (anisotropic)', ''), &
    & cfg_keyword_t('tau_s', cfg_active, cfg_phase, cfg_list, cfg_real, 0, 0.0_rk, big, .true., &
    &   .false., '', 'MPa', 'required (saturation laws)', &
    &   'saturation slip strength: one value, or per family or system (anisotropic)', ''), &
    & cfg_keyword_t('n', cfg_active, cfg_phase, cfg_list, cfg_real, 0, 0.0_rk, big, .false., &
    &   .false., '', '-', 'required (saturation laws)', &
    &   'hardening exponent: one value, or per family or system (anisotropic)', ''), &
    & cfg_keyword_t('interaction_matrix', cfg_active, cfg_phase, cfg_repeat, cfg_text, 0, -big, &
    &   big, .false., .false., '', '-', '1 diag., 1.4 off', &
    &   'one row of the interaction matrix (N lines of N values; saturation_anisotropic)', '', &
    &   values_text='*N* reals >= 0, one row of the matrix per line, *N* lines in the ' // &
    &   'slip-system order (*N* the number of active slip systems)'), &
    & cfg_keyword_t('saturation_evolution', cfg_active, cfg_phase, cfg_single, cfg_word, 1, &
    &   -big, big, .false., .false., 'on off', '', 'off', &
    &   'saturation strength tau_s0 (Gammadot/gammadot_s0)^m_prime (saturation laws)', ''), &
    & cfg_keyword_t('gammadot_s0', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', '1/s', 'required (saturation_evolution)', &
    &   'reference slip rate of saturation evolution', ''), &
    & cfg_keyword_t('m_prime', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', '-', 'required (saturation_evolution)', &
    &   'exponent of saturation evolution', ''), &
    & cfg_keyword_t('cyclic_hardening', cfg_active, cfg_phase, cfg_single, cfg_word, 1, -big, &
    &   big, .false., .false., 'on off', '', 'off', &
    &   'cyclic hardening: only slip past a critical accumulated slip hardens', ''), &
    & cfg_keyword_t('cyclic_a', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', '-', 'required (cyclic_hardening on)', &
    &   'critical accumulated slip a (crss/tau_s)^c: coefficient a (0 < a < 1)', ''), &
    & cfg_keyword_t('cyclic_c', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', '-', 'required (cyclic_hardening on)', &
    &   'critical accumulated slip a (crss/tau_s)^c: exponent c', ''), &
    & cfg_keyword_t('precipitate_strengthening', cfg_active, cfg_phase, cfg_single, cfg_word, &
    &   1, -big, big, .false., .false., 'on off', '', 'off', &
    &   'precipitate strength tau_p added to tau_0 (cutting or bowing)', ''), &
    & cfg_keyword_t('a_p', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, .true., &
    &   .false., '', 'MPa', '', 'precipitate cutting coefficient (cutting set: a_p, b_p)', ''), &
    & cfg_keyword_t('b_p', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, .true., &
    &   .false., '', 'mm', '', 'Burgers vector of the precipitate (cutting set)', ''), &
    & cfg_keyword_t('c_p', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, .true., &
    &   .false., '', '-', '', 'precipitate bowing coefficient (bowing set: c_p, b_m)', ''), &
    & cfg_keyword_t('b_m', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, .true., &
    &   .false., '', 'mm', '', 'Burgers vector of the matrix (bowing set)', ''), &
    & cfg_keyword_t('f_p', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, 1.0_rk, &
    &   .true., .false., '', '-', 'required (precipitate_strengthening on)', &
    &   'precipitate volume fraction', ''), &
    & cfg_keyword_t('r_p', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, .true., &
    &   .false., '', 'mm', 'required (precipitate_strengthening on)', &
    &   'mean precipitate radius', ''), &
    ! Kocks-Mecking (step 7.2, D60, addendum_kocks_mecking.md §4); mu_0 is Hill's average.
    & cfg_keyword_t('km_b', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, .true., &
    &   .false., '', 'mm', 'required (kocks_mecking)', 'Burgers vector b (Kocks-Mecking)', ''), &
    & cfg_keyword_t('km_h_n', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .false., .false., '', '-', 'required (kocks_mecking)', &
    &   'network hardening coefficient h_n (Kocks-Mecking)', ''), &
    & cfg_keyword_t('km_rho_0', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', '1/mm^2', 'required (kocks_mecking)', &
    &   'initial network dislocation density rho_0 (Kocks-Mecking)', ''), &
    & cfg_keyword_t('km_k1', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .false., .false., '', '1/mm', 'required (kocks_mecking)', &
    &   'dislocation storage coefficient k_1 (Kocks-Mecking)', ''), &
    & cfg_keyword_t('km_k20', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .false., .false., '', '-', 'required (kocks_mecking)', &
    &   'dynamic recovery coefficient k_20 (Kocks-Mecking)', ''), &
    & cfg_keyword_t('km_gammadot_k0', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, &
    &   big, .true., .false., '', '1/s', 'required (kocks_mecking)', &
    &   'reference slip rate of k_2 = k_20 (gammadot_k0/Gammadot)^(1/km_n)', ''), &
    & cfg_keyword_t('km_n', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', '-', 'required (kocks_mecking)', &
    &   'rate exponent n of k_2 = k_20 (gammadot_k0/Gammadot)^(1/km_n)', ''), &
    ! Irradiation softening (step 7.3, D61, addendum_irradiation.md §6), with kocks_mecking.
    & cfg_keyword_t('irradiation', cfg_active, cfg_phase, cfg_single, cfg_word, 1, -big, big, &
    &   .false., .false., 'on off', '', 'off', &
    &   'irradiation softening: damage tensor from dislocation loops (kocks_mecking)', ''), &
    & cfg_keyword_t('irr_h_d', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .false., .false., '', '-', 'required (irradiation)', &
    &   'defect hardening coefficient h_d (irradiation)', ''), &
    & cfg_keyword_t('irr_eta', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .false., .false., '', '-', 'required (irradiation)', &
    &   'loop annihilation efficiency eta (irradiation)', ''), &
    & cfg_keyword_t('irr_loop_density', cfg_active, cfg_phase, cfg_single, cfg_real, 1, 0.0_rk, &
    &   big, .false., .false., '', '1/mm^3', 'required (irradiation)', &
    &   'dislocation loop number density rho_l (irradiation)', ''), &
    & cfg_keyword_t('irr_loop_diameter', cfg_active, cfg_phase, cfg_single, cfg_real, 1, &
    &   0.0_rk, big, .true., .false., '', 'mm', 'required (irradiation)', &
    &   'dislocation loop diameter d_l (irradiation)', ''), &
    & cfg_keyword_t('irr_loop_normals', cfg_active, cfg_phase, cfg_list, cfg_text, 0, -big, big, &
    &   .false., .false., '', '', 'required (irradiation)', &
    &   'loop-normal families, <111> <100> <110> (irradiation; fcc, bcc)', '', &
    &   values_text='one or more of ``<111>``, ``<100>``, ``<110>``, each once')]

  ! The rest of the table (split from the material keywords: gfortran's limit of 255
  ! continuation lines per statement).
  type(cfg_keyword_t), parameter :: cfg_keywords_rest(*) = [ &
    ! --- Boundary conditions and loading (active) ---
    & cfg_keyword_t('boundary_conditions', cfg_active, cfg_global, cfg_single, cfg_word, 1, &
    &   -big, big, .false., .false., &
    &   'uniaxial_grip uniaxial_symmetry uniaxial_minimal triaxial general periodic', '', &
    &   'required', 'boundary-condition preset (constrained nodes and dofs)', ''), &
    & cfg_keyword_t('set_bc', cfg_active, cfg_global, cfg_repeat, cfg_text, 0, -big, big, &
    &   .false., .false., '', 'mm/s', '', &
    &   'set_bc vel <nset> <x, y or z> <v> [...] (general)', '', &
    &   values_text='``vel <nset> <d> <v> [<d> <v> [<d> <v>]]``: a node set (or a node id), ' // &
    &   'then one to three pairs of a direction ``x``, ``y`` or ``z`` and a velocity; the ' // &
    &   'line may be repeated'), &
    ! Multi-point constraints (step 8.3a, D73; dev/theory/addendum_periodic_mpc.md §3.2, §4).
    & cfg_keyword_t('set_mpc1', cfg_active, cfg_global, cfg_repeat, cfg_text, 0, -big, big, &
    &   .false., .false., '', '', '', &
    &   'set_mpc1 vel <nset> <x, y or z> [...]: equal velocity components on a node set', '', &
    &   values_text='``vel <nset> <d> [<d> [<d>]]``: a node set, then one to three ' // &
    &   'directions ``x``, ``y`` or ``z``; the line may be repeated'), &
    & cfg_keyword_t('set_mpc', cfg_active, cfg_global, cfg_repeat, cfg_text, 0, -big, big, &
    &   .false., .false., '', 'mm/s', '', &
    &   'set_mpc vel <nset> <x, y or z> <c> <d> [...]: v = c v(first node) + d', '', &
    &   values_text='``vel <nset> <d> <c> <o> [...]``: a node set, then one to three ' // &
    &   'triples of a direction ``x``, ``y`` or ``z``, a nonzero coefficient and an offset; ' // &
    &   'the line may be repeated'), &
    & cfg_keyword_t('control_component', cfg_active, cfg_global, cfg_single, cfg_word, 1, -big, &
    &   big, .false., .false., '11 22 33 23 13 12 32 31 21', '', 'required (triaxial, periodic)', &
    &   'component that sets the step durations (shear: periodic)', ''), &
    ! Periodic conditions (step 8.3b, D71-D73; dev/theory/addendum_periodic_mpc.md §3.3, §4).
    & cfg_keyword_t('periodic_kinematics', cfg_active, cfg_global, cfg_single, cfg_word, 1, &
    &   -big, big, .false., .false., 'small_strain deformation_gradient_rate', '', &
    &   'small_strain', 'periodic: imposed rates of the symmetric or the full Fbar-dot', ''), &
    & cfg_keyword_t('periodic_preconditioner', cfg_active, cfg_global, cfg_single, cfg_word, 1, &
    &   -big, big, .false., .false., 'lumped_diagonal exact_diagonal', '', 'lumped_diagonal', &
    &   'periodic: Jacobi diagonal of the reduced system', ''), &
    & cfg_keyword_t('number_of_steps', cfg_active, cfg_global, cfg_single, cfg_int, 1, 1.0_rk, &
    &   big, .false., .false., '', '', 'required', 'number of loading steps', ''), &
    & cfg_keyword_t('strain_rate', cfg_active, cfg_global, cfg_comp_per_step, cfg_real, 0, &
    &   -big, big, .false., .true., '', '1/s', '', &
    &   'strain rate of a component (periodic: of Fbar), per step', ''), &
    & cfg_keyword_t('target_strain', cfg_active, cfg_global, cfg_comp_per_step, cfg_real, 0, &
    &   -big, big, .false., .true., '', '-', '', &
    &   'strain target of a component (periodic: of Fbar - I), per step', ''), &
    & cfg_keyword_t('target_load', cfg_active, cfg_global, cfg_comp_per_step, cfg_real, 0, &
    &   -big, big, .false., .true., '', 'N', '', &
    &   'load target: force on the face x1/y1/z1 of the control component, per step', ''), &
    & cfg_keyword_t('target_stress', cfg_active, cfg_global, cfg_comp_per_step, cfg_real, 0, &
    &   -big, big, .false., .true., '', 'MPa', '', &
    &   'stress target: that force over the current face area, per step', ''), &
    & cfg_keyword_t('target_time', cfg_active, cfg_global, cfg_per_step, cfg_real, 0, 0.0_rk, &
    &   big, .true., .true., '', 's', '', 'time target (cumulative), per step', ''), &
    ! Load and stress rates, step lengths and dwells (step 8.2c, D66, Q-T7, Q-T8;
    ! dev/theory/addendum_triaxial_load_control.md §3.1-§3.4, §4).
    & cfg_keyword_t('load_rate', cfg_active, cfg_global, cfg_comp_per_step, cfg_real, 0, -big, &
    &   big, .false., .true., '', 'N/s', '', &
    &   'load rate of a component (force on its face x1/y1/z1), per step; 0 holds it', ''), &
    & cfg_keyword_t('stress_rate', cfg_active, cfg_global, cfg_comp_per_step, cfg_real, 0, &
    &   -big, big, .false., .true., '', 'MPa/s', '', &
    &   'stress rate of a component (that force over the current face area), per step', ''), &
    & cfg_keyword_t('target_duration', cfg_active, cfg_global, cfg_per_step, cfg_real, 0, &
    &   0.0_rk, big, .true., .true., '', 's', '', &
    &   'length of the step from its actual start, per step (or a target)', ''), &
    & cfg_keyword_t('number_of_incrs', cfg_active, cfg_global, cfg_per_step, cfg_int, 0, &
    &   1.0_rk, big, .false., .true., '', '', 'required (or dtime, dstrain)', &
    &   'increments per step (one of number_of_incrs, dtime, dstrain per step)', ''), &
    & cfg_keyword_t('dtime', cfg_active, cfg_global, cfg_per_step, cfg_real, 0, 0.0_rk, big, &
    &   .true., .true., '', 's', 'required (or number_of_incrs, dstrain)', &
    &   'time increment per step (or number_of_incrs, dstrain)', ''), &
    & cfg_keyword_t('dstrain', cfg_active, cfg_global, cfg_per_step, cfg_real, 0, 0.0_rk, big, &
    &   .true., .true., '', '-', 'required (or number_of_incrs, dtime)', &
    &   'strain increment per step, loading component (or number_of_incrs, dtime)', ''), &
    & cfg_keyword_t('dtime_min', cfg_active, cfg_global, cfg_per_step, cfg_real, 0, 0.0_rk, &
    &   big, .true., .true., '', 's', 'required (unknown step duration)', &
    &   'smallest time increment of a load or stress target step', ''), &
    & cfg_keyword_t('dtime_factor', cfg_active, cfg_global, cfg_single, cfg_real, 1, 1.0_rk, &
    &   big, .false., .false., '', '-', '1.001', &
    &   'factor on the increment aimed at a load or stress target', ''), &
    ! Tolerances of load and stress rows and the landing attempts (step 8.2a, D66, Q-T13,
    ! Q-T18; dev/theory/addendum_triaxial_load_control.md §3.6, §4).
    & cfg_keyword_t('load_tol', cfg_active, cfg_global, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', 'N', '0.1', &
    &   'absolute tolerance of load rows: a load target is landed on within it', ''), &
    & cfg_keyword_t('stress_tol', cfg_active, cfg_global, cfg_single, cfg_real, 1, 0.0_rk, &
    &   big, .true., .false., '', 'MPa', '0.1', &
    &   'absolute tolerance of stress rows: a stress target is landed on within it', ''), &
    & cfg_keyword_t('max_bc_iter', cfg_active, cfg_global, cfg_single, cfg_int, 1, 1.0_rk, big, &
    &   .false., .false., '', '', '10', &
    &   'landing attempts of an increment; iterations of each boundary-velocity loop', ''), &
    ! Lower bound of the velocity perturbations of triaxial load control (steps 8.2b, 8.2c, D66).
    & cfg_keyword_t('min_pert_frac', cfg_active, cfg_global, cfg_single, cfg_real, 1, 0.0_rk, &
    &   big, .true., .false., '', '-', '0.001', &
    &   'smallest velocity perturbation: a fraction of the control velocity or elastic estimate', ''), &
    ! Dwell increments (step 8.2c, D66; v2's max_strain_incr).
    & cfg_keyword_t('max_strain_incr', cfg_active, cfg_global, cfg_single, cfg_real, 1, 0.0_rk, &
    &   big, .true., .false., '', '-', '0.001', &
    &   'largest strain increment of the control component in a dwell', ''), &
    & cfg_keyword_t('check_necking', cfg_active, cfg_global, cfg_single, cfg_word, 1, -big, &
    &   big, .false., .false., 'on off', '', 'off', &
    &   'stop when the control-face load falls (checked every increment)', ''), &
    & cfg_keyword_t('print_results', cfg_active, cfg_global, cfg_per_step, cfg_int, 0, 0.0_rk, &
    &   1.0_rk, .false., .false., '', '', '1', 'write the step results: 1 yes, 0 no', ''), &
    & cfg_keyword_t('max_incr', cfg_active, cfg_global, cfg_single, cfg_int, 1, 1.0_rk, big, &
    &   .false., .false., '', '', '50000', 'run limit: increments', ''), &
    & cfg_keyword_t('max_total_time', cfg_active, cfg_global, cfg_single, cfg_real, 1, 0.0_rk, &
    &   big, .true., .false., '', 's', '12000', 'run limit: time', ''), &
    & cfg_keyword_t('max_strain', cfg_active, cfg_global, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', '-', '', &
    &   'run limit: largest normal engineering strain (no limit if absent)', ''), &
    & cfg_keyword_t('max_eqstrain', cfg_active, cfg_global, cfg_single, cfg_real, 1, 0.0_rk, &
    &   big, .true., .false., '', '-', '', &
    &   'run limit: equivalent engineering strain (no limit if absent)', ''), &
    ! --- Solver (active) ---
    & cfg_keyword_t('cg_tol', cfg_active, cfg_global, cfg_single, cfg_real, 1, 0.0_rk, 1.0_rk, &
    &   .true., .false., '', '-', '1e-8', 'CG relative tolerance', ''), &
    & cfg_keyword_t('cg_max_iters', cfg_active, cfg_global, cfg_single, cfg_int, 1, 1.0_rk, &
    &   big, .false., .false., '', '', '16000', 'CG iteration limit', ''), &
    ! Preconditioner (step 9.2; dev/theory/addendum_block_jacobi.md §3; default block_jacobi
    ! after the study, Q-B1: jacobi with constraints, config_mod).
    & cfg_keyword_t('cg_preconditioner', cfg_active, cfg_global, cfg_single, cfg_word, 1, &
    &   -big, big, .false., .false., 'jacobi block_jacobi none', '', 'block_jacobi', &
    &   'CG preconditioner (jacobi with MPCs or periodicity; none: for debugging only)', ''), &
    & cfg_keyword_t('global_max_iters', cfg_active, cfg_global, cfg_single, cfg_int, 1, 1.0_rk, &
    &   big, .false., .false., '', '', '50', 'global iteration limit', ''), &
    & cfg_keyword_t('global_tol', cfg_active, cfg_global, cfg_single, cfg_real, 1, 0.0_rk, &
    &   1.0_rk, .true., .false., '', '-', '5e-4', 'relative velocity-change tolerance', &
    &   ''), &
    & cfg_keyword_t('global_tol_abs', cfg_active, cfg_global, cfg_single, cfg_real, 1, 0.0_rk, &
    &   big, .true., .false., '', 'mm/s', '1e-10', 'absolute velocity-change tolerance', &
    &   ''), &
    & cfg_keyword_t('global_tol_slow', cfg_active, cfg_global, cfg_single, cfg_real, 1, 0.0_rk, &
    &   1.0_rk, .true., .false., '', '-', '5e-4', &
    &   'relative tolerance accepted after a slow NR, or in SA after NR', ''), &
    & cfg_keyword_t('global_nr_switch', cfg_active, cfg_global, cfg_single, cfg_real, 1, &
    &   0.0_rk, big, .false., .false., '', '-', '1e-2', &
    &   'relative velocity change below which SA switches to NR', ''), &
    & cfg_keyword_t('global_nr_slow_ratio', cfg_active, cfg_global, cfg_single, cfg_real, 1, &
    &   0.0_rk, &
    &   big, .true., .false., '', '-', '0.2', 'NR is slow below this fraction of its rate', ''), &
    & cfg_keyword_t('midnode_linearization', cfg_active, cfg_global, cfg_single, cfg_word, 1, &
    &   -big, big, .false., .false., 'on off', '', 'off', &
    &   'move the mid-nodes of tet10 elements to the edge midpoints every increment', ''), &
    & cfg_keyword_t('initialization', cfg_active, cfg_global, cfg_single, cfg_word, 1, -big, &
    &   big, .false., .false., 'elastic_step isotropic_viscoplastic', '', 'elastic_step', &
    &   'start of crystal-plasticity runs (velocity guess and stress start of increment 1)', ''), &
    & cfg_keyword_t('point_tol', cfg_active, cfg_global, cfg_single, cfg_real, 1, 0.0_rk, big, &
    &   .true., .false., '', '1/s', '1e-4', 'material point Newton tolerance (residual norm)', &
    &   ''), &
    & cfg_keyword_t('point_tol_state', cfg_active, cfg_global, cfg_single, cfg_real, 1, 0.0_rk, &
    &   big, .true., .false., '', '-', '1e-5', &
    &   'material point state-loop tolerance (relative to tau_0)', ''), &
    & cfg_keyword_t('point_tol_hardening', cfg_active, cfg_global, cfg_single, cfg_real, 1, &
    &   0.0_rk, big, .true., .false., '', '-', '1e-4', &
    &   'material point strength-update tolerance (relative to tau_0)', ''), &
    & cfg_keyword_t('point_max_iters_newton', cfg_active, cfg_global, cfg_single, cfg_int, 1, &
    &   1.0_rk, big, .false., .false., '', '', '100', 'material point Newton iteration limit', &
    &   ''), &
    & cfg_keyword_t('point_max_iters_hardening', cfg_active, cfg_global, cfg_single, cfg_int, &
    &   1, 1.0_rk, big, .false., .false., '', '', '10', &
    &   'material point strength-update iteration limit', ''), &
    & cfg_keyword_t('point_max_iters_state', cfg_active, cfg_global, cfg_single, cfg_int, 1, &
    &   1.0_rk, big, .false., .false., '', '', '100', &
    &   'material point state-loop iteration limit', ''), &
    & cfg_keyword_t('subincrementation', cfg_active, cfg_global, cfg_single, cfg_word, 1, -big, &
    &   big, .false., .false., 'on off', '', 'off', &
    &   'repeat a rejected increment as halves', ''), &
    & cfg_keyword_t('subincrement_max_depth', cfg_active, cfg_global, cfg_single, cfg_int, 1, &
    &   1.0_rk, big, .false., .false., '', '', '4', &
    &   'subdivisions of one increment, down to dt/2^value (subincrementation on)', ''), &
    & cfg_keyword_t('subincrement_max_successive', cfg_active, cfg_global, cfg_single, cfg_int, &
    &   1, 1.0_rk, big, .false., .false., '', '', '10', &
    &   'successive increments that may be subdivided (subincrementation on)', ''), &
    ! --- Input files and output (active) ---
    & cfg_keyword_t('read_ori_from_file', cfg_active, cfg_global, cfg_flag, cfg_text, 0, -big, &
    &   big, .false., .false., '', '', '', 'orientations from simulation.ori', ''), &
    & cfg_keyword_t('read_phase_from_file', cfg_active, cfg_global, cfg_flag, cfg_text, 0, -big, &
    &   big, .false., .false., '', '', '', 'phases from simulation.phase', ''), &
    & cfg_keyword_t('read_from_opt_file', cfg_active, cfg_global, cfg_flag, cfg_text, 0, -big, &
    &   big, .false., .false., '', '', '', 'values from simulation.opt', ''), &
    ! Restart (step 8.4, dev/theory/addendum_restart.md §6).
    & cfg_keyword_t('restart', cfg_active, cfg_global, cfg_single, cfg_word, 1, -big, big, &
    &   .false., .false., 'on off', '', 'off', &
    &   'continue simulation.sim from its restart files', ''), &
    & cfg_keyword_t('print', cfg_active, cfg_global, cfg_repeat, cfg_word, 1, -big, big, &
    &   .false., .false., 'coo disp vel stress strain_el ori forces crss slip sliprate rss ' // &
    &   'strain_pl strain_pl_eq defrate_pl spinrate rotrate work strain strain_eq ' // &
    &   'strain_el_eq stress_eq defrate defrate_eq defrate_pl_eq velgrad rotrate_spin ' // &
    &   'rotrate_slip work_pl workrate workrate_pl strain_pl_sample strain_pl_sample_eq ' // &
    &   'convergence km_rho km_k2 km_rhorate dmg dmgrate restart', '', '', &
    &   'output variable (one per line); restart: restart files at every step end', '', &
    &   later_words='strain_th temp twin_frac', &
    &   later_note='Stages 7-11, dev/cfg_keywords.md section 6'), &
    ! Output formats (step 9.1, D93): the .sim text files, the HDF5 file, or both.
    & cfg_keyword_t('output_format', cfg_active, cfg_global, cfg_list, cfg_word, 0, -big, big, &
    &   .false., .false., 'ascii hdf5', '', 'ascii', &
    &   'output formats: ascii (.sim text files), hdf5 (simulation.h5 and .xdmf), or both', ''), &
    ! --- Planned (recognized; an error until their stage) ---
    & cfg_keyword_t(name='temperature', status=cfg_planned, note='Stage 10'), &
    & cfg_keyword_t(name='alpha_11', status=cfg_planned, note='Stage 10'), &
    & cfg_keyword_t(name='alpha_22', status=cfg_planned, note='Stage 10'), &
    & cfg_keyword_t(name='alpha_33', status=cfg_planned, note='Stage 10'), &
    & cfg_keyword_t(name='reference_temperature', status=cfg_planned, note='Stage 10')]

  ! The whole table.
  type(cfg_keyword_t), parameter :: cfg_keywords(*) = [cfg_keywords_material, cfg_keywords_rest]

  !> Number of material keywords: the entries 1..cfg_num_material of cfg_keywords (the restart
  !> fingerprint, step 8.4, requires them unchanged).
  integer, parameter :: cfg_num_material = size(cfg_keywords_material)

contains

  !> Index of a keyword in the table, 0 if unknown.
  !>
  !> Arguments:
  !>   name  (in)  character  keyword
  !> Locals:
  !>   i     integer          table index
  pure function cfg_find_keyword(name) result(index)

    character(len=*), intent(in) :: name
    integer :: index

    integer :: i

    index = 0
    do i = 1, size(cfg_keywords)
      if (trim(cfg_keywords(i)%name) == name) then
        index = i
        return
      end if
    end do

  end function cfg_find_keyword

  !> Write the table as a Markdown document (keyword, status, scope, values, unit, default,
  !> description or note): the generated keyword reference (io_format.md §2).
  !>
  !> Arguments:
  !>   unit  (in)  integer  output unit
  !> Locals:
  !>   i       integer         table index
  !>   k       type(cfg_keyword_t)  the entry
  !>   status, scope, values  character(:), allocatable  column texts
  subroutine cfg_table_markdown(unit)

    integer, intent(in) :: unit

    integer :: i
    type(cfg_keyword_t) :: k
    character(len=:), allocatable :: status, scope, values

    write (unit, '(a)') '# FEPX configuration keywords (generated from cfg_table_mod)'
    write (unit, '(a)') ''
    write (unit, '(a)') '| keyword | status | scope | values | unit | default | description |'
    write (unit, '(a)') '|---|---|---|---|---|---|---|'
    do i = 1, size(cfg_keywords)
      k = cfg_keywords(i)
      select case (k%status)
      case (cfg_active)
        status = 'active'
      case default
        status = 'planned (' // trim(k%note) // ')'
      end select
      scope = merge('phase ', 'global', k%scope == cfg_phase)
      if (k%status /= cfg_active) scope = ''
      values = form_text(k)
      write (unit, '(a)') '| `' // trim(k%name) // '` | ' // status // ' | ' // trim(scope) &
        & // ' | ' // values // ' | ' // trim(k%unit) // ' | ' // trim(k%default_text) // &
        & ' | ' // trim(k%info) // ' |'
    end do

  end subroutine cfg_table_markdown

  !> Description of an entry's values for the documentation.
  !>
  !> Arguments:
  !>   k  (in)  type(cfg_keyword_t)  entry
  !> Locals:
  !>   kind_text  character(:), allocatable  value kind
  function form_text(k) result(text)

    type(cfg_keyword_t), intent(in) :: k
    character(len=:), allocatable :: text

    character(len=:), allocatable :: kind_text

    if (k%status /= cfg_active) then
      text = ''
      return
    end if
    select case (k%kind)
    case (cfg_int)
      kind_text = 'integer'
    case (cfg_real)
      kind_text = 'real'
    case (cfg_word)
      kind_text = 'one of: ' // trim(k%words)
      if (k%later_words /= '') kind_text = kind_text // ' (later: ' // trim(k%later_words) // &
        & ', ' // trim(k%later_note) // ')'
    case default
      kind_text = 'text'
    end select
    select case (k%form)
    case (cfg_per_step)
      text = kind_text // ', one per step'
    case (cfg_comp_per_step)
      text = 'component, then ' // kind_text // ', one per step'
    case (cfg_repeat)
      text = kind_text // ' (repeatable line)'
    case (cfg_flag)
      text = '(no value)'
    case (cfg_list)
      text = kind_text // ', one or more'
    case default
      text = kind_text
    end select

  end function form_text

end module cfg_table_mod
