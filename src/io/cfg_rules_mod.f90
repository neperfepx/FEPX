! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The option-compatibility rules of the configuration (ROADMAP D85, D104): every combination of
!> keywords (or of keyword values) that the configuration builder (config_mod) rejects or
!> ignores, in one table. config_mod tests the condition of each rule and takes the message
!> from the table (cfg_rule_message); the manual's notes at each option and its compatibility
!> table are generated from the same table (cfg_manual_mod, fepx_keywords), so that the
!> program and the manual state the same rules and no note is missing on either side.
!>
!> Source : FEPX (ROADMAP step 9.9, D85, D104); the rules are those of config_mod (steps 4.3 to
!>          9.2, cited there). Contains no numerical method.
!>
!> A rule says that its subjects (one or more keywords, or values of a keyword) are used only
!> with a condition (rule_only_with), cannot be used with it (rule_not_with), or are ignored
!> with it (rule_ignored_with, a warning). Texts may mark keywords and values with backquotes:
!> the message drops them, the manual turns them into literals.
!>
!> Entry fields (cfg_rule_t):
!>   name                short identifier, equal to the name of its index parameter without
!>                       the "rule_" prefix (checked by the unit test)
!>   subjects            the keywords the rule constrains, separated by blanks: the manual
!>                       gives the rule at each of them
!>   subject_doc         the subjects as the manual states them, when not simply the keywords
!>                       (a value, e.g. "`cg_preconditioner block_jacobi` or `none`"); empty
!>                       otherwise
!>   plural              subject_doc names several options ("are used only with")
!>   relation            rule_only_with, rule_not_with or rule_ignored_with
!>   condition           the other option(s), e.g. "`boundary_conditions periodic`"
!>   condition_keywords  the keywords of the condition, separated by blanks: the manual gives
!>                       the rule at each of them too (both directions, D85)
!>   reason              explanation, in parentheses after the condition; empty if none
!>
!> Public parameters: the relation codes; the rule indices rule_<name>; cfg_rules (the table).
!> Public routines:
!>   cfg_rule_message   the message of a rule for one option as given
!>   cfg_rule_sentence  the sentence of a rule for the manual (reStructuredText)
!>   cfg_rule_involves  whether a rule concerns a keyword (as a subject or in its condition)
module cfg_rules_mod
  implicit none

  private

  public :: cfg_rule_t, cfg_rules, cfg_rule_message, cfg_rule_sentence, cfg_rule_involves
  public :: rule_only_with, rule_not_with, rule_ignored_with

  ! Relation codes.
  integer, parameter :: rule_only_with = 1, rule_not_with = 2, rule_ignored_with = 3

  !> One rule (fields documented in the module header).
  type :: cfg_rule_t
    character(len=32) :: name = ''
    character(len=160) :: subjects = ''
    character(len=200) :: subject_doc = ''
    logical :: plural = .false.
    integer :: relation = rule_only_with
    character(len=200) :: condition = ''
    character(len=160) :: condition_keywords = ''
    character(len=240) :: reason = ''
  end type cfg_rule_t

  ! The rule indices, in table order (material, loading, solver, output).
  public :: rule_c_over_a, rule_slip_model, rule_hardening_model, rule_saturation_law
  public :: rule_interaction_matrix, rule_saturation_evolution, rule_cyclic_hardening
  public :: rule_precipitate_strengthening, rule_kocks_mecking, rule_irradiation_law
  public :: rule_irradiation_crystal, rule_irradiation_parameters
  public :: rule_periodic_keywords, rule_shear_components, rule_small_strain_components
  public :: rule_rates_general, rule_set_bc_preset, rule_mpc_periodic, rule_mpc_triaxial
  public :: rule_mpc_load_control, rule_load_rows_periodic, rule_control_general
  public :: rule_load_rate_target, rule_stress_rate_target, rule_dtime_factor, rule_load_tol
  public :: rule_stress_tol, rule_max_bc_iter, rule_min_pert_frac, rule_max_strain_incr
  public :: rule_dtime_min, rule_number_of_incrs_open, rule_number_of_incrs_dwell
  public :: rule_dstrain_rate, rule_necking_periodic, rule_necking_loaded
  public :: rule_subincrementation, rule_preconditioner_periodic, rule_preconditioner_mpc
  public :: rule_print_kocks_mecking, rule_print_irradiation, rule_print_forces_periodic
  integer, parameter :: rule_c_over_a = 1, rule_slip_model = 2, rule_hardening_model = 3, &
    & rule_saturation_law = 4, rule_interaction_matrix = 5, rule_saturation_evolution = 6, &
    & rule_cyclic_hardening = 7, rule_precipitate_strengthening = 8, rule_kocks_mecking = 9, &
    & rule_irradiation_law = 10, rule_irradiation_crystal = 11, &
    & rule_irradiation_parameters = 12, rule_periodic_keywords = 13, &
    & rule_shear_components = 14, rule_small_strain_components = 15, rule_rates_general = 16, &
    & rule_set_bc_preset = 17, rule_mpc_periodic = 18, rule_mpc_triaxial = 19, &
    & rule_mpc_load_control = 20, rule_load_rows_periodic = 21, rule_control_general = 22, &
    & rule_load_rate_target = 23, rule_stress_rate_target = 24, rule_dtime_factor = 25, &
    & rule_load_tol = 26, rule_stress_tol = 27, rule_max_bc_iter = 28, &
    & rule_min_pert_frac = 29, rule_max_strain_incr = 30, rule_dtime_min = 31, &
    & rule_number_of_incrs_open = 32, rule_number_of_incrs_dwell = 33, &
    & rule_dstrain_rate = 34, rule_necking_periodic = 35, rule_necking_loaded = 36, &
    & rule_subincrementation = 37, rule_preconditioner_periodic = 38, &
    & rule_preconditioner_mpc = 39, rule_print_kocks_mecking = 40, &
    & rule_print_irradiation = 41, rule_print_forces_periodic = 42

  ! The table.
  type(cfg_rule_t), parameter :: cfg_rules(*) = [ &
    ! --- Material ---
    ! c/a: the hexagonal and tetragonal crystals only (step 4.3).
    & cfg_rule_t(name='c_over_a', subjects='c_over_a', &
    &   condition='`crystal_type hcp` and `bct`', condition_keywords='crystal_type'), &
    ! Slip and hardening: crystal plasticity only (steps 6.1, 6.2, D48, D49).
    & cfg_rule_t(name='slip_model', subjects='slip_families m gammadot_0 tau_0', &
    &   condition='`model crystal_plasticity`', condition_keywords='model'), &
    & cfg_rule_t(name='hardening_model', subjects='hardening h_0 tau_s n interaction_matrix', &
    &   condition='`model crystal_plasticity`', condition_keywords='model'), &
    ! The saturation laws' parameters and modifiers (steps 6.2, 7.1, 7.2; D60: Kocks-Mecking
    ! takes no modifiers).
    & cfg_rule_t(name='saturation_law', subjects='h_0 tau_s n interaction_matrix ' // &
    &   'saturation_evolution cyclic_hardening precipitate_strengthening', &
    &   condition='a saturation hardening law', condition_keywords='hardening', &
    &   reason='`hardening saturation_isotropic` or `saturation_anisotropic`'), &
    & cfg_rule_t(name='interaction_matrix', subjects='interaction_matrix', &
    &   condition='`hardening saturation_anisotropic`', condition_keywords='hardening'), &
    & cfg_rule_t(name='saturation_evolution', subjects='gammadot_s0 m_prime', &
    &   condition='`saturation_evolution on`', condition_keywords='saturation_evolution'), &
    & cfg_rule_t(name='cyclic_hardening', subjects='cyclic_a cyclic_c', &
    &   condition='`cyclic_hardening on`', condition_keywords='cyclic_hardening'), &
    & cfg_rule_t(name='precipitate_strengthening', subjects='a_p b_p c_p b_m f_p r_p', &
    &   condition='`precipitate_strengthening on`', &
    &   condition_keywords='precipitate_strengthening'), &
    ! Kocks-Mecking and irradiation (steps 7.2, 7.3, D60, D61).
    & cfg_rule_t(name='kocks_mecking', &
    &   subjects='km_b km_h_n km_rho_0 km_k1 km_k20 km_gammadot_k0 km_n', &
    &   condition='`hardening kocks_mecking`', condition_keywords='hardening'), &
    & cfg_rule_t(name='irradiation_law', subjects='irradiation', &
    &   condition='`hardening kocks_mecking`', condition_keywords='hardening'), &
    & cfg_rule_t(name='irradiation_crystal', subjects='irradiation', &
    &   subject_doc='`irradiation on`', condition='`crystal_type fcc` and `bcc`', &
    &   condition_keywords='crystal_type', &
    &   reason='the dislocation loop normals are defined only for fcc and bcc'), &
    & cfg_rule_t(name='irradiation_parameters', &
    &   subjects='irr_h_d irr_eta irr_loop_density irr_loop_diameter irr_loop_normals', &
    &   condition='`irradiation on`', condition_keywords='irradiation'), &
    ! --- Boundary conditions and loading ---
    ! Periodic conditions (step 8.3b, D71-D73, addendum_periodic_mpc.md §3.3, §4).
    & cfg_rule_t(name='periodic_keywords', subjects='periodic_kinematics ' // &
    &   'periodic_preconditioner', condition='`boundary_conditions periodic`', &
    &   condition_keywords='boundary_conditions'), &
    & cfg_rule_t(name='shear_components', subjects='strain_rate target_strain target_load ' // &
    &   'target_stress load_rate stress_rate control_component', &
    &   subject_doc='The shear components (`23`, `13`, `12`, `32`, `31`, `21`) of ' // &
    &   '`strain_rate`, `target_strain`, `target_load`, `target_stress`, `load_rate`, ' // &
    &   '`stress_rate` and `control_component`', plural=.true., &
    &   condition='`boundary_conditions periodic`', condition_keywords='boundary_conditions'), &
    & cfg_rule_t(name='small_strain_components', subjects='strain_rate target_strain ' // &
    &   'control_component', subject_doc='The components `32`, `31` and `21` of ' // &
    &   '`strain_rate`, `target_strain` and `control_component`', plural=.true., &
    &   relation=rule_not_with, condition='`periodic_kinematics small_strain`', &
    &   condition_keywords='periodic_kinematics', reason='it takes the shear components ' // &
    &   '23, 13 and 12, one per pair; `deformation_gradient_rate` takes the nine components'), &
    ! Rates and velocity conditions (steps 4.3, 8.2c).
    & cfg_rule_t(name='rates_general', subjects='strain_rate load_rate stress_rate', &
    &   relation=rule_not_with, condition='`boundary_conditions general`', &
    &   condition_keywords='boundary_conditions', &
    &   reason='the velocities are given by the `set_bc` lines'), &
    & cfg_rule_t(name='set_bc_preset', subjects='set_bc', relation=rule_not_with, &
    &   condition='`boundary_conditions uniaxial_grip`, `uniaxial_symmetry`, ' // &
    &   '`uniaxial_minimal` or `triaxial`', condition_keywords='boundary_conditions', &
    &   reason='the preset sets the velocities'), &
    ! Multi-point constraints (step 8.3a, D73, Q-P8; with periodic conditions v2-fixed C68).
    & cfg_rule_t(name='mpc_periodic', subjects='set_mpc set_mpc1', &
    &   subject_doc='Multi-point constraints (`set_mpc`, `set_mpc1`)', plural=.true., &
    &   relation=rule_not_with, condition='`boundary_conditions periodic`', &
    &   condition_keywords='boundary_conditions'), &
    & cfg_rule_t(name='mpc_triaxial', subjects='set_mpc set_mpc1', &
    &   subject_doc='Multi-point constraints (`set_mpc`, `set_mpc1`)', plural=.true., &
    &   relation=rule_not_with, condition='`boundary_conditions triaxial`', &
    &   condition_keywords='boundary_conditions'), &
    & cfg_rule_t(name='mpc_load_control', subjects='set_mpc set_mpc1', &
    &   subject_doc='Multi-point constraints (`set_mpc`, `set_mpc1`)', plural=.true., &
    &   relation=rule_not_with, condition='load or stress targets or rates', &
    &   condition_keywords='target_load target_stress load_rate stress_rate', &
    &   reason='`target_load`, `target_stress`, `load_rate`, `stress_rate`: load control'), &
    & cfg_rule_t(name='load_rows_periodic', subjects='target_load target_stress load_rate ' // &
    &   'stress_rate', relation=rule_not_with, condition='`boundary_conditions periodic`', &
    &   condition_keywords='boundary_conditions', &
    &   reason='no macroscopic load or stress control'), &
    ! The control component (D45).
    & cfg_rule_t(name='control_general', subjects='control_component', &
    &   relation=rule_not_with, condition='`boundary_conditions general`', &
    &   condition_keywords='boundary_conditions', &
    &   reason='the velocities are given by the `set_bc` lines'), &
    ! Load and stress rates of the control component (step 8.2c, addendum_triaxial_load_
    ! control.md §3.1-§3.3).
    & cfg_rule_t(name='load_rate_target', subjects='load_rate', relation=rule_not_with, &
    &   condition='a `target_stress` of the same component in the same step', &
    &   condition_keywords='target_stress', reason='a `load_rate` ends at a ' // &
    &   '`target_load`, `target_strain`, `target_time` or `target_duration`'), &
    & cfg_rule_t(name='stress_rate_target', subjects='stress_rate', relation=rule_not_with, &
    &   condition='a `target_load` of the same component in the same step', &
    &   condition_keywords='target_load', reason='a `stress_rate` ends at a ' // &
    &   '`target_stress`, `target_strain`, `target_time` or `target_duration`'), &
    ! The options of load and stress control (steps 8.1, 8.2a-8.2c; addendum_load_targets.md
    ! §3.4, addendum_triaxial_load_control.md §3.4, §3.6, §4).
    & cfg_rule_t(name='dtime_factor', subjects='dtime_factor', &
    &   condition='a `target_load` or `target_stress` (or a `load_rate` or `stress_rate` row)', &
    &   condition_keywords='target_load target_stress load_rate stress_rate'), &
    & cfg_rule_t(name='load_tol', subjects='load_tol', &
    &   condition='a `target_load` or `load_rate` row', &
    &   condition_keywords='target_load load_rate'), &
    & cfg_rule_t(name='stress_tol', subjects='stress_tol', &
    &   condition='a `target_stress` or `stress_rate` row', &
    &   condition_keywords='target_stress stress_rate'), &
    & cfg_rule_t(name='max_bc_iter', subjects='max_bc_iter', &
    &   condition='a `target_load` or `target_stress` (or a `load_rate` or `stress_rate` row)', &
    &   condition_keywords='target_load target_stress load_rate stress_rate'), &
    & cfg_rule_t(name='min_pert_frac', subjects='min_pert_frac', &
    &   condition='a `target_load` or `target_stress` row on a component other than the ' // &
    &   'control component, or a `load_rate` or `stress_rate` row', &
    &   condition_keywords='target_load target_stress load_rate stress_rate'), &
    & cfg_rule_t(name='max_strain_incr', subjects='max_strain_incr', condition='a dwell', &
    &   condition_keywords='load_rate stress_rate', reason='a step with `load_rate 0` or ' // &
    &   '`stress_rate 0` on the control component'), &
    ! The increments of a step whose duration is not known (steps 8.1, 8.2c, Q-L4).
    & cfg_rule_t(name='dtime_min', subjects='dtime_min', condition='a load or stress ' // &
    &   'target at a strain rate, or a strain target at a load or stress rate', &
    &   condition_keywords='target_load target_stress load_rate stress_rate', &
    &   reason='give "-" for the other steps'), &
    & cfg_rule_t(name='number_of_incrs_open', subjects='number_of_incrs', &
    &   relation=rule_not_with, condition='a load or stress target at a strain rate, or a ' // &
    &   'strain target at a load or stress rate', &
    &   condition_keywords='target_load target_stress load_rate stress_rate', &
    &   reason='the duration of the step is not known: give `dtime` or `dstrain` (the ' // &
    &   'largest increment) and `dtime_min`'), &
    & cfg_rule_t(name='number_of_incrs_dwell', subjects='number_of_incrs', &
    &   relation=rule_not_with, condition='a dwell', &
    &   condition_keywords='load_rate stress_rate', reason='a step with `load_rate 0` or ' // &
    &   '`stress_rate 0` on the control component: give `dtime`, its largest increment; ' // &
    &   'its increments follow `max_strain_incr`'), &
    & cfg_rule_t(name='dstrain_rate', subjects='dstrain', &
    &   condition='a `strain_rate` on the control component', &
    &   condition_keywords='strain_rate', reason='not in a step at a `load_rate` or ' // &
    &   '`stress_rate`: give `number_of_incrs` or `dtime`'), &
    ! Necking (step 8.1, step 8.3b).
    & cfg_rule_t(name='necking_periodic', subjects='check_necking', &
    &   subject_doc='`check_necking on`', relation=rule_not_with, &
    &   condition='`boundary_conditions periodic`', condition_keywords='boundary_conditions', &
    &   reason='no face loads'), &
    & cfg_rule_t(name='necking_loaded', subjects='check_necking', &
    &   subject_doc='`check_necking on`', condition='a loaded component', &
    &   condition_keywords='target_strain target_load target_stress', &
    &   reason='a `target_strain`, `target_load` or `target_stress` row; the load of its ' // &
    &   'face x1, y1 or z1 is checked'), &
    ! --- Solver ---
    ! Sub-incrementation (Q35, D63, addendum_subincrementation.md §4).
    & cfg_rule_t(name='subincrementation', &
    &   subjects='subincrement_max_depth subincrement_max_successive', &
    &   condition='`subincrementation on`', condition_keywords='subincrementation'), &
    ! The preconditioner with constraints (step 9.2, addendum_block_jacobi.md §2, §3, Q-B1).
    & cfg_rule_t(name='preconditioner_periodic', subjects='cg_preconditioner', &
    &   subject_doc='`cg_preconditioner block_jacobi` and `none`', plural=.true., &
    &   relation=rule_not_with, condition='`boundary_conditions periodic`', &
    &   condition_keywords='boundary_conditions', reason='only `jacobi`'), &
    & cfg_rule_t(name='preconditioner_mpc', subjects='cg_preconditioner', &
    &   subject_doc='`cg_preconditioner block_jacobi` and `none`', plural=.true., &
    &   relation=rule_not_with, condition='multi-point constraints', &
    &   condition_keywords='set_mpc set_mpc1', reason='`set_mpc`, `set_mpc1`: only `jacobi`'), &
    ! --- Output ---
    ! Outputs of a model (steps 7.2, 7.3, D60, D61) and of the face forces (step 8.3b, Q-P11).
    & cfg_rule_t(name='print_kocks_mecking', subjects='print', &
    &   subject_doc='`print km_rho`, `print km_k2` and `print km_rhorate`', plural=.true., &
    &   condition='a phase with `hardening kocks_mecking`', condition_keywords='hardening'), &
    & cfg_rule_t(name='print_irradiation', subjects='print', &
    &   subject_doc='`print dmg` and `print dmgrate`', plural=.true., &
    &   condition='a phase with `irradiation on`', condition_keywords='irradiation'), &
    & cfg_rule_t(name='print_forces_periodic', subjects='print', &
    &   subject_doc='`print forces`', relation=rule_ignored_with, &
    &   condition='`boundary_conditions periodic`', condition_keywords='boundary_conditions', &
    &   reason='no face forces file; the macroscopic stress is written to results/macro')]

contains

  !> The message of a rule for one option as given: "[<prefix>: ]<subject> <relation>
  !> <condition>[ (<reason>[; <detail>])]", without the backquotes of the table.
  !>
  !> Arguments:
  !>   id       (in)  integer    rule index
  !>   subject  (in)  character  the option as given, e.g. "km_b" or "cg_preconditioner none"
  !>   prefix   (in, optional)  character  context before the message, e.g. "step 2"
  !>   detail   (in, optional)  character  context in the parentheses, e.g. "phase 1"
  !>   plural   (in, optional)  logical    the subject names several options (default false)
  !> Locals:
  !>   rule     type(cfg_rule_t)           the entry
  !>   inside   character(:), allocatable  the text in the parentheses
  !>   many     logical                    plural subject
  function cfg_rule_message(id, subject, prefix, detail, plural) result(text)

    integer, intent(in) :: id
    character(len=*), intent(in) :: subject
    character(len=*), intent(in), optional :: prefix, detail
    logical, intent(in), optional :: plural
    character(len=:), allocatable :: text

    type(cfg_rule_t) :: rule
    character(len=:), allocatable :: inside
    logical :: many

    rule = cfg_rules(id)
    many = .false.
    if (present(plural)) many = plural

    ! The subject, the relation and the condition.
    text = subject // ' ' // relation_text(rule%relation, many) // ' ' // &
      & plain_text(trim(rule%condition))
    if (present(prefix)) text = prefix // ': ' // text

    ! The reason and the detail, in one pair of parentheses.
    inside = plain_text(trim(rule%reason))
    if (present(detail)) then
      if (len(inside) > 0 .and. len(detail) > 0) then
        inside = inside // '; ' // detail
      else
        inside = inside // detail
      end if
    end if
    if (len(inside) > 0) text = text // ' (' // inside // ')'

  end function cfg_rule_message

  !> The sentence of a rule for the manual (reStructuredText): the subjects (subject_doc, or
  !> the keywords as literals), the relation, the condition and the reason, with literals for
  !> the backquoted texts, ending with a period.
  !>
  !> Arguments:
  !>   id  (in)  integer  rule index
  !> Locals:
  !>   rule      type(cfg_rule_t)           the entry
  !>   subject   character(:), allocatable  the subjects as written
  !>   many      logical                    plural subject
  function cfg_rule_sentence(id) result(text)

    integer, intent(in) :: id
    character(len=:), allocatable :: text

    type(cfg_rule_t) :: rule
    character(len=:), allocatable :: subject
    logical :: many

    rule = cfg_rules(id)

    ! The subjects: as stated by the table, or the keywords listed as literals.
    if (len_trim(rule%subject_doc) > 0) then
      subject = rst_text(trim(rule%subject_doc))
      many = rule%plural
    else
      subject = literal_list(trim(rule%subjects))
      many = index(trim(rule%subjects), ' ') > 0
    end if

    ! The relation, the condition and the reason.
    text = subject // ' ' // relation_text(rule%relation, many) // ' ' // &
      & rst_text(trim(rule%condition))
    if (len_trim(rule%reason) > 0) text = text // ' (' // rst_text(trim(rule%reason)) // ')'
    text = text // '.'

    ! A sentence starts with a capital letter (a literal subject keeps its case).
    if (text(1:1) >= 'a' .and. text(1:1) <= 'z') text(1:1) = achar(iachar(text(1:1)) - 32)

  end function cfg_rule_sentence

  !> Whether a rule concerns a keyword: one of its subjects or of its condition keywords.
  !>
  !> Arguments:
  !>   id       (in)  integer    rule index
  !>   keyword  (in)  character  keyword
  pure function cfg_rule_involves(id, keyword) result(involved)

    integer, intent(in) :: id
    character(len=*), intent(in) :: keyword
    logical :: involved

    involved = index(' ' // trim(cfg_rules(id)%subjects) // ' ', ' ' // trim(keyword) // ' ') &
      & > 0 .or. index(' ' // trim(cfg_rules(id)%condition_keywords) // ' ', ' ' // &
      & trim(keyword) // ' ') > 0

  end function cfg_rule_involves

  !> The words of a relation.
  !>
  !> Arguments:
  !>   relation  (in)  integer  relation code
  !>   plural    (in)  logical  plural subject
  pure function relation_text(relation, plural) result(text)

    integer, intent(in) :: relation
    logical, intent(in) :: plural
    character(len=:), allocatable :: text

    select case (relation)
    case (rule_only_with)
      text = trim(merge('are', 'is ', plural)) // ' used only with'
    case (rule_not_with)
      text = 'cannot be used with'
    case default
      text = trim(merge('are', 'is ', plural)) // ' ignored with'
    end select

  end function relation_text

  !> A text of the table without its backquotes (for messages).
  !>
  !> Arguments:
  !>   text  (in)  character  table text
  !> Locals:
  !>   i     integer          character index
  pure function plain_text(text) result(plain)

    character(len=*), intent(in) :: text
    character(len=:), allocatable :: plain

    integer :: i

    plain = ''
    do i = 1, len(text)
      if (text(i:i) /= '`') plain = plain // text(i:i)
    end do

  end function plain_text

  !> A text of the table as reStructuredText: each backquote doubled (a literal).
  !>
  !> Arguments:
  !>   text  (in)  character  table text
  !> Locals:
  !>   i     integer          character index
  pure function rst_text(text) result(rst)

    character(len=*), intent(in) :: text
    character(len=:), allocatable :: rst

    integer :: i

    rst = ''
    do i = 1, len(text)
      if (text(i:i) == '`') then
        rst = rst // '``'
      else
        rst = rst // text(i:i)
      end if
    end do

  end function rst_text

  !> A blank-separated list of keywords as reStructuredText literals: "``a``", "``a`` and
  !> ``b``", "``a``, ``b`` and ``c``".
  !>
  !> Arguments:
  !>   words  (in)  character  the keywords, separated by single blanks
  !> Locals:
  !>   rest   character(:), allocatable  the words not yet written
  !>   word   character(:), allocatable  the current word
  !>   blank  integer                    position of the next blank in rest
  pure function literal_list(words) result(list)

    character(len=*), intent(in) :: words
    character(len=:), allocatable :: list

    character(len=:), allocatable :: rest, word
    integer :: blank

    list = ''
    rest = trim(adjustl(words))
    do while (len(rest) > 0)
      ! The next word, and the separator before it.
      blank = index(rest, ' ')
      if (blank == 0) then
        word = rest
        rest = ''
      else
        word = rest(1:blank - 1)
        rest = trim(adjustl(rest(blank + 1:)))
      end if
      if (len(list) > 0 .and. len(rest) == 0) then
        list = list // ' and '
      else if (len(list) > 0) then
        list = list // ', '
      end if
      list = list // '``' // word // '``'
    end do

  end function literal_list

end module cfg_rules_mod
