! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of the option-compatibility rules (cfg_rules_mod): every index parameter names its
!> entry; every keyword of a rule is an implemented keyword of the keyword table; every text has
!> balanced backquotes; the messages (with and without prefix, reason, detail, plural subject)
!> and the manual's sentences (subject lists of one, two and three keywords, subject_doc, a
!> capital first letter) are as specified; cfg_rule_involves finds subjects and condition
!> keywords, and not substrings of them.
!>
!> Source : FEPX (ROADMAP step 9.9, D104). The use of each rule by the configuration builder is
!>          covered by test_config_errors and checked by the docs test docs_rules_used.
!>
!> Locals:
!>   names     character(32)(42)          the rule names in index order (the index parameters)
!>   ids       integer(42)                the index parameters in the same order
!>   i         integer                    rule index
!>   words     character(:), allocatable  the keywords of a rule (subjects and condition)
!>   word      character(:), allocatable  one keyword
!>   blank     integer                    position of the next blank in words
!>   k         integer                    keyword-table index of word
program test_cfg_rules

  use cfg_table_mod, only: cfg_keywords, cfg_find_keyword, cfg_active
  use cfg_rules_mod
  use assert_mod, only: assert_equal, assert_true, test_summary

  implicit none

  character(len=32), parameter :: names(42) = [character(len=32) :: 'c_over_a', 'slip_model', &
    & 'hardening_model', 'saturation_law', 'interaction_matrix', 'saturation_evolution', &
    & 'cyclic_hardening', 'precipitate_strengthening', 'kocks_mecking', 'irradiation_law', &
    & 'irradiation_crystal', 'irradiation_parameters', 'periodic_keywords', &
    & 'shear_components', 'small_strain_components', 'rates_general', 'set_bc_preset', &
    & 'mpc_periodic', 'mpc_triaxial', 'mpc_load_control', 'load_rows_periodic', &
    & 'control_general', 'load_rate_target', 'stress_rate_target', 'dtime_factor', 'load_tol', &
    & 'stress_tol', 'max_bc_iter', 'min_pert_frac', 'max_strain_incr', 'dtime_min', &
    & 'number_of_incrs_open', 'number_of_incrs_dwell', 'dstrain_rate', 'necking_periodic', &
    & 'necking_loaded', 'subincrementation', 'preconditioner_periodic', 'preconditioner_mpc', &
    & 'print_kocks_mecking', 'print_irradiation', 'print_forces_periodic']
  integer, parameter :: ids(42) = [rule_c_over_a, rule_slip_model, rule_hardening_model, &
    & rule_saturation_law, rule_interaction_matrix, rule_saturation_evolution, &
    & rule_cyclic_hardening, rule_precipitate_strengthening, rule_kocks_mecking, &
    & rule_irradiation_law, rule_irradiation_crystal, rule_irradiation_parameters, &
    & rule_periodic_keywords, rule_shear_components, rule_small_strain_components, &
    & rule_rates_general, rule_set_bc_preset, rule_mpc_periodic, rule_mpc_triaxial, &
    & rule_mpc_load_control, rule_load_rows_periodic, rule_control_general, &
    & rule_load_rate_target, rule_stress_rate_target, rule_dtime_factor, rule_load_tol, &
    & rule_stress_tol, rule_max_bc_iter, rule_min_pert_frac, rule_max_strain_incr, &
    & rule_dtime_min, rule_number_of_incrs_open, rule_number_of_incrs_dwell, &
    & rule_dstrain_rate, rule_necking_periodic, rule_necking_loaded, rule_subincrementation, &
    & rule_preconditioner_periodic, rule_preconditioner_mpc, rule_print_kocks_mecking, &
    & rule_print_irradiation, rule_print_forces_periodic]
  integer :: i, blank, k
  character(len=:), allocatable :: words, word

  ! The table: one entry per index parameter, in order, each named as its parameter.
  call assert_equal(size(cfg_rules), size(names), 'number of rules')
  do i = 1, size(names)
    call assert_equal(ids(i), i, 'index of rule ' // trim(names(i)))
    call assert_equal(trim(cfg_rules(ids(i))%name), trim(names(i)), 'name of rule ' // &
      & trim(names(i)))
  end do

  ! Every rule: subjects, a relation, a condition with its keywords; every keyword implemented;
  ! backquotes in pairs.
  do i = 1, size(cfg_rules)
    associate (rule => cfg_rules(i))
      call assert_true(len_trim(rule%subjects) > 0, 'subjects of ' // trim(rule%name))
      call assert_true(len_trim(rule%condition) > 0, 'condition of ' // trim(rule%name))
      call assert_true(len_trim(rule%condition_keywords) > 0, 'condition keywords of ' // &
        & trim(rule%name))
      call assert_true(any(rule%relation == [rule_only_with, rule_not_with, &
        & rule_ignored_with]), 'relation of ' // trim(rule%name))
      call assert_true(mod(count_char(rule%subject_doc), 2) == 0 .and. &
        & mod(count_char(rule%condition), 2) == 0 .and. &
        & mod(count_char(rule%reason), 2) == 0, 'backquotes in pairs in ' // trim(rule%name))
      call assert_true(len_trim(rule%subject_doc) > 0 .or. .not. rule%plural, &
        & 'plural only with subject_doc in ' // trim(rule%name))
      words = trim(rule%subjects) // ' ' // trim(rule%condition_keywords)
      do while (len(words) > 0)
        ! The next keyword, and the rest of the list.
        blank = index(words, ' ')
        if (blank == 0) then
          word = words
          words = ''
        else
          word = words(1:blank - 1)
          words = trim(adjustl(words(blank + 1:)))
        end if
        k = cfg_find_keyword(word)
        call assert_true(k > 0, 'keyword ' // word // ' of rule ' // trim(rule%name) // &
          & ' in the keyword table')
        if (k > 0) call assert_true(cfg_keywords(k)%status == cfg_active, 'keyword ' // &
          & word // ' of rule ' // trim(rule%name) // ' implemented')
      end do
    end associate
  end do

  ! Messages: subject, relation and condition without backquotes; prefix; reason and detail in
  ! one pair of parentheses; plural subject.
  call assert_equal(cfg_rule_message(rule_load_tol, 'load_tol'), &
    & 'load_tol is used only with a target_load or load_rate row', 'message, no reason')
  call assert_equal(cfg_rule_message(rule_kocks_mecking, 'km_b', detail='phase 2'), &
    & 'km_b is used only with hardening kocks_mecking (phase 2)', 'message, detail only')
  call assert_equal(cfg_rule_message(rule_preconditioner_periodic, &
    & 'cg_preconditioner none'), 'cg_preconditioner none cannot be used with ' // &
    & 'boundary_conditions periodic (only jacobi)', 'message, reason only')
  call assert_equal(cfg_rule_message(rule_saturation_law, 'h_0', detail='phase 1 has ' // &
    & 'hardening none'), 'h_0 is used only with a saturation hardening law (hardening ' // &
    & 'saturation_isotropic or saturation_anisotropic; phase 1 has hardening none)', &
    & 'message, reason and detail')
  call assert_equal(cfg_rule_message(rule_dstrain_rate, 'dstrain', prefix='step 3'), &
    & 'step 3: dstrain is used only with a strain_rate on the control component (not in a ' // &
    & 'step at a load_rate or stress_rate: give number_of_incrs or dtime)', 'message, prefix')
  call assert_equal(cfg_rule_message(rule_shear_components, 'shear components', &
    & prefix='strain_rate 12', plural=.true.), 'strain_rate 12: shear components are used ' // &
    & 'only with boundary_conditions periodic', 'message, plural subject')
  call assert_equal(cfg_rule_message(rule_print_forces_periodic, 'print forces'), &
    & 'print forces is ignored with boundary_conditions periodic (no face forces file; the ' // &
    & 'macroscopic stress is written to results/macro)', 'message, ignored')

  ! Sentences: literal subject lists, subject_doc, condition literals, reason, period.
  call assert_equal(cfg_rule_sentence(rule_c_over_a), '``c_over_a`` is used only with ' // &
    & '``crystal_type hcp`` and ``bct``.', 'sentence, one subject')
  call assert_equal(cfg_rule_sentence(rule_saturation_evolution), '``gammadot_s0`` and ' // &
    & '``m_prime`` are used only with ``saturation_evolution on``.', 'sentence, two subjects')
  call assert_equal(cfg_rule_sentence(rule_subincrementation), &
    & '``subincrement_max_depth`` and ``subincrement_max_successive`` are used only with ' // &
    & '``subincrementation on``.', 'sentence, two subjects (long)')
  call assert_equal(cfg_rule_sentence(rule_cyclic_hardening), '``cyclic_a`` and ' // &
    & '``cyclic_c`` are used only with ``cyclic_hardening on``.', 'sentence, cyclic')
  call assert_equal(cfg_rule_sentence(rule_slip_model), '``slip_families``, ``m``, ' // &
    & '``gammadot_0`` and ``tau_0`` are used only with ``model crystal_plasticity``.', &
    & 'sentence, four subjects')
  call assert_equal(cfg_rule_sentence(rule_preconditioner_periodic), &
    & '``cg_preconditioner block_jacobi`` and ``none`` cannot be used with ' // &
    & '``boundary_conditions periodic`` (only ``jacobi``).', 'sentence, subject_doc, reason')
  call assert_equal(cfg_rule_sentence(rule_max_strain_incr), '``max_strain_incr`` is used ' // &
    & 'only with a dwell (a step with ``load_rate 0`` or ``stress_rate 0`` on the control ' // &
    & 'component).', 'sentence, plain condition')
  call assert_equal(cfg_rule_sentence(rule_mpc_triaxial), 'Multi-point constraints ' // &
    & '(``set_mpc``, ``set_mpc1``) cannot be used with ``boundary_conditions triaxial``.', &
    & 'sentence, capital letter')

  ! Involvement: subjects and condition keywords, whole words only.
  call assert_true(cfg_rule_involves(rule_kocks_mecking, 'km_b'), 'involves a subject')
  call assert_true(cfg_rule_involves(rule_kocks_mecking, 'hardening'), &
    & 'involves a condition keyword')
  call assert_true(.not. cfg_rule_involves(rule_kocks_mecking, 'km'), 'not a substring')
  call assert_true(.not. cfg_rule_involves(rule_kocks_mecking, 'model'), 'not another keyword')

  call test_summary()

contains

  !> Number of backquotes in a text.
  !>
  !> Arguments:
  !>   text  (in)  character  the text
  !> Locals:
  !>   j     integer          character index
  pure function count_char(text) result(n)

    character(len=*), intent(in) :: text
    integer :: n

    integer :: j

    n = 0
    do j = 1, len(text)
      if (text(j:j) == '`') n = n + 1
    end do

  end function count_char

end program test_cfg_rules
