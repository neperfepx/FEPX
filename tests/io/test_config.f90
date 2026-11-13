! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of the configuration reader (cfg_table_mod, cfg_store_mod, config_mod): the keyword
!> table's consistency, and the typed configuration built from valid files (comments, blank
!> lines, tabs, word case, defaults, "-" in per-step rows, phases of several crystal types,
!> presets and general velocity conditions, output variables; step 6.1: slip families in the
!> order given, per-family and per-system m and gammadot_0, defaults and "all", D48).
!>
!> Source : FEPX (ROADMAP steps 4.3 and 6.1). The error paths are covered by
!>          test_config_errors.
!>
!> Locals:
!>   config    type(config_t)            configuration read back
!>   k         type(cfg_keyword_t)       a table entry
!>   i, j      integer                   indices
!>   unit      integer                   output unit
!>   x         real(rk)                  a parsed default
!>   ios       integer                   read status
!>   lines     integer                   lines of the generated Markdown table
!>   buffer    character(512)            a line read back
!>   text      character(256)            a generated cfg line
!>   rows      character(256)(15)        generated interaction_matrix lines
program test_config

  use kinds_mod, only: rk
  use crystal_mod, only: stiffness_hexagonal, stiffness_cubic, stiffness_tetragonal
  use cfg_table_mod, only: cfg_keyword_t, cfg_keywords, cfg_find_keyword, cfg_table_markdown, &
    & cfg_active, cfg_planned, cfg_real, cfg_int, cfg_word, cfg_single
  use config_mod, only: config_t, config_read, incr_by_count, incr_by_dtime, incr_by_dstrain, &
    & target_by_load, target_by_stress, target_by_duration, rate_by_strain, rate_by_load, &
    & rate_by_stress
  use assert_mod, only: assert_equal, assert_close, assert_true, test_summary

  implicit none

  type(config_t) :: config
  type(cfg_keyword_t) :: k
  integer :: i, j, unit, ios, lines
  real(rk) :: x
  character(len=512) :: buffer
  character(len=256) :: text, rows(15)

  ! --- The keyword table ---

  ! Names unique, lower case, found by cfg_find_keyword; planned entries name their stage; a
  ! default is "required...", empty or a value: numeric defaults parse and lie in the range,
  ! word defaults are allowed words.
  do i = 1, size(cfg_keywords)
    k = cfg_keywords(i)
    call assert_equal(cfg_find_keyword(trim(k%name)), i, 'find ' // trim(k%name))
    call assert_true(verify(trim(k%name), 'abcdefghijklmnopqrstuvwxyz0123456789_') == 0, &
      & 'lower case ' // trim(k%name))
    do j = i + 1, size(cfg_keywords)
      call assert_true(k%name /= cfg_keywords(j)%name, 'unique ' // trim(k%name))
    end do
    if (k%status == cfg_planned) then
      call assert_true(len_trim(k%note) > 0, 'note ' // trim(k%name))
    end if
    if (k%status == cfg_active .and. k%form == cfg_single .and. len_trim(k%default_text) > 0 &
      & .and. index(k%default_text, 'required') /= 1) then
      select case (k%kind)
      case (cfg_real, cfg_int)
        read (k%default_text, *, iostat=ios) x
        call assert_equal(ios, 0, 'default parses ' // trim(k%name))
        call assert_true(x >= k%lo .and. x <= k%hi, 'default in range ' // trim(k%name))
      case (cfg_word)
        call assert_true(index(' ' // trim(k%words) // ' ', ' ' // trim(k%default_text) // &
          & ' ') > 0, 'default allowed ' // trim(k%name))
      end select
    end if
  end do
  call assert_equal(cfg_find_keyword('no_such_keyword'), 0, 'unknown keyword')

  ! Generated documentation: header, separator and one row per entry.
  open (newunit=unit, file='test_config_table.md', status='replace', action='write')
  call cfg_table_markdown(unit)
  close (unit)
  open (newunit=unit, file='test_config_table.md', status='old', action='read')
  lines = 0
  do
    read (unit, '(a)', iostat=ios) buffer
    if (ios /= 0) exit
    lines = lines + 1
    if (lines == 6) call assert_true(index(buffer, '| `phase` | active | global |') == 1, &
      & 'markdown row of phase')
  end do
  close (unit)
  call assert_equal(lines, 4 + size(cfg_keywords), 'markdown lines')

  ! --- File 1: hcp, uniaxial_minimal on 33, defaults ---
  call write_file('test_config_1.cfg', [character(len=64) :: &
    & '## Material  # a section header is a comment', &
    & 'number_of_phases 1', &
    & '', &
    & 'phase 1', &
    & '  crystal_type HCP', &
    & '  model elastic', &
    & '  c11 162400.0', &
    & achar(9) // 'c12' // achar(9) // '92000.0', &
    & '  c13 69000.0', &
    & '  c33 180700.0   # D3', &
    & '  c44 46700.0', &
    & '  c_over_a 1.587', &
    & 'boundary_conditions uniaxial_minimal', &
    & 'number_of_steps 2', &
    & '#              comp step1  step2', &
    & 'strain_rate    33   1.0e-3 2.0e-3', &
    & 'target_strain  33   0.001  0.002', &
    & 'number_of_incrs     5      7', &
    & 'print coo', &
    & 'print stress'])
  call config_read(config, 'test_config_1.cfg')
  call assert_equal(config%num_phases, 1, 'f1 phases')
  call assert_equal(config%phases(1)%crystal_type, 'hcp', 'f1 crystal type (lower case)')
  call assert_equal(config%phases(1)%model, 'elastic', 'f1 model')
  call assert_close(config%phases(1)%c12, 92000.0_rk, 'f1 c12 (tabs)')
  call assert_close(config%phases(1)%c33, 180700.0_rk, 'f1 c33')
  call assert_close(config%phases(1)%c66, 0.0_rk, 'f1 c66 unused')
  call assert_close(config%phases(1)%c_over_a, 1.587_rk, 'f1 c/a')
  call assert_close(config%phases(1)%stiff_el_c, stiffness_hexagonal(162400.0_rk, 92000.0_rk, &
    & 69000.0_rk, 180700.0_rk, 46700.0_rk), 'f1 stiffness', rel_tol=0.0_rk, abs_tol=0.0_rk)
  call assert_equal(config%bc_type, 'uniaxial_minimal', 'f1 bc type')
  call assert_equal(config%control_axis, 3, 'f1 load axis')
  call assert_equal(size(config%vel_bcs), 0, 'f1 no set_bc')
  call assert_equal(config%num_steps, 2, 'f1 steps')
  call assert_close(config%strain_rate(3, :), [1.0e-3_rk, 2.0e-3_rk], 'f1 strain rate 33')
  call assert_close(config%target_strain(3, :), [0.001_rk, 0.002_rk], 'f1 target 33')
  call assert_true(all(config%has_strain_rate(3, :)) .and. &
    & .not. any(config%has_strain_rate(1:2, :)), 'f1 strain rate given on 33 only')
  call assert_true(all(config%has_target_strain(3, :)) .and. &
    & .not. any(config%has_target_strain(1:2, :)), 'f1 target given on 33 only')
  call assert_true(.not. any(config%has_target_time), 'f1 no time target')
  call assert_equal(config%num_incrs(1), 5, 'f1 incrs step 1')
  call assert_equal(config%num_incrs(2), 7, 'f1 incrs step 2')
  call assert_true(all(config%incr_kind == incr_by_count), 'f1 increments by count')
  call assert_true(all(config%print_results), 'f1 default print_results')
  call assert_equal(config%max_incr, 50000, 'f1 default max_incr')
  call assert_close(config%max_total_time, 12000.0_rk, 'f1 default max_total_time')
  call assert_true(.not. (config%has_max_strain .or. config%has_max_eqstrain), &
    & 'f1 no strain limits by default')
  call assert_close(config%cg_tol, 1.0e-8_rk, 'f1 default cg_tol')
  call assert_equal(config%cg_max_iters, 16000, 'f1 default cg_max_iters')
  call assert_equal(config%cg_preconditioner, 'block_jacobi', &
    & 'f1 default cg_preconditioner (no constraints: block_jacobi)')
  call assert_close(config%point_tol_state, 1.0e-5_rk, 'f1 default point_tol_state')
  call assert_close(config%point_tol_hardening, 1.0e-4_rk, 'f1 default point_tol_hardening')
  call assert_equal(config%global_max_iters, 50, 'f1 default global_max_iters')
  call assert_close(config%global_tol, 5.0e-4_rk, 'f1 default global_tol')
  call assert_close(config%global_tol_abs, 1.0e-10_rk, 'f1 default global_tol_abs')
  call assert_equal(config%initialization, 'elastic_step', 'f1 default initialization')
  call assert_close(config%point_tol, 1.0e-4_rk, 'f1 default point_tol')
  call assert_close(config%global_tol_slow, 5.0e-4_rk, 'f1 default global_tol_slow')
  call assert_close(config%global_nr_switch, 1.0e-2_rk, 'f1 default global_nr_switch')
  call assert_close(config%global_nr_slow_ratio, 0.2_rk, 'f1 default global_nr_slow_ratio')
  call assert_equal(config%point_max_iters_newton, 100, 'f1 default point_max_iters_newton')
  call assert_equal(config%point_max_iters_state, 100, 'f1 default point_max_iters_state')
  call assert_equal(config%point_max_iters_hardening, 10, 'f1 default point_max_iters_hardening')
  call assert_true(.not. config%midnode_linearization, 'f1 default midnode_linearization')
  call assert_true(.not. config%subincrementation, 'f1 default subincrementation')
  call assert_equal(config%subincrement_max_depth, 4, 'f1 default subincrement_max_depth')
  call assert_equal(config%subincrement_max_successive, 10, &
    & 'f1 default subincrement_max_successive')
  call assert_true(.not. (config%read_ori_from_file .or. config%read_phase_from_file .or. &
    & config%read_from_opt_file), 'f1 no input-file flags')
  call assert_equal(size(config%print_vars), 2, 'f1 print count')
  call assert_equal(config%print_vars(1), 'coo', 'f1 print 1')
  call assert_equal(config%print_vars(2), 'stress', 'f1 print 2')

  ! --- File 2: fcc and bct phases, general conditions, "-" entries, every option set ---
  call write_file('test_config_2.cfg', [character(len=64) :: &
    & 'number_of_phases 2', &
    & 'phase 2', &
    & '  crystal_type bct', &
    & '  model elastic', &
    & '  c11 72300', &
    & '  c12 59400', &
    & '  c13 35800', &
    & '  c33 88400', &
    & '  c44 22000', &
    & '  c66 24000', &
    & '  c_over_a 0.546', &
    & 'phase 1', &
    & '  crystal_type fcc', &
    & '  model elastic', &
    & '  c44 75400', &
    & '  c12 151400', &
    & '  c11 245000', &
    & 'boundary_conditions general', &
    & 'set_bc vel z0 z 0.0', &
    & 'set_bc vel x0 x 0.0 y 0.0', &
    & 'set_bc vel 1 y 0.0', &
    & 'set_bc vel z1 z 1.5e-3 x 0 y -2', &
    & 'set_mpc1 vel x1 x y', &
    & 'set_mpc vel y1 z -0.5 2e-3 x 1.0 0', &
    & 'number_of_steps 3', &
    & 'target_strain 33 0.001 - -', &
    & 'target_time - 2.0 3.5', &
    & 'number_of_incrs 1 - -', &
    & 'dtime - 0.5 -', &
    & 'dstrain - - 2e-4', &
    & 'print_results 1 0 1', &
    & 'max_incr 100', &
    & 'max_total_time 10.0', &
    & 'max_strain 0.1', &
    & 'max_eqstrain 0.15', &
    & 'cg_tol 1.0e-10', &
    & 'cg_max_iters 500', &
    & 'global_max_iters 7', &
    & 'global_tol 1e-6', &
    & 'global_tol_abs 1e-12', &
    & 'midnode_linearization on', &
    & 'subincrementation on', &
    & 'subincrement_max_depth 6', &
    & 'subincrement_max_successive 3', &
    & 'read_ori_from_file', &
    & 'read_phase_from_file', &
    & 'read_from_opt_file', &
    & 'print forces', &
    & 'print strain_el', &
    & 'print disp'])
  call config_read(config, 'test_config_2.cfg')
  call assert_equal(config%num_phases, 2, 'f2 phases')
  call assert_equal(config%phases(1)%crystal_type, 'fcc', 'f2 phase 1 type')
  call assert_equal(config%phases(2)%crystal_type, 'bct', 'f2 phase 2 type')
  call assert_close(config%phases(1)%stiff_el_c, stiffness_cubic(245000.0_rk, 151400.0_rk, &
    & 75400.0_rk), 'f2 fcc stiffness', rel_tol=0.0_rk, abs_tol=0.0_rk)
  call assert_close(config%phases(2)%stiff_el_c, stiffness_tetragonal(72300.0_rk, 59400.0_rk, &
    & 35800.0_rk, 88400.0_rk, 22000.0_rk, 24000.0_rk), 'f2 bct stiffness', rel_tol=0.0_rk, &
    & abs_tol=0.0_rk)
  call assert_close(config%phases(1)%c_over_a, 0.0_rk, 'f2 fcc c/a unset')
  call assert_close(config%phases(2)%c_over_a, 0.546_rk, 'f2 bct c/a')
  call assert_equal(config%bc_type, 'general', 'f2 bc type')
  call assert_equal(config%control_axis, 3, 'f2 load axis from the target_strain row')
  call assert_equal(size(config%vel_bcs), 7, 'f2 velocity conditions')
  call check_bc(2, 'x0', 1, 0.0_rk, 20)
  call check_bc(3, 'x0', 2, 0.0_rk, 20)
  call check_bc(4, '1', 2, 0.0_rk, 21)
  call check_bc(5, 'z1', 3, 1.5e-3_rk, 22)
  call check_bc(7, 'z1', 2, -2.0_rk, 22)
  ! Multi-point constraints (step 8.3a): one relation per node set and direction.
  call assert_equal(size(config%mpcs), 4, 'f2 multi-point constraints')
  call assert_equal(config%mpcs(1)%nset, 'x1', 'f2 mpc 1 nset')
  call assert_equal(config%mpcs(2)%comp, 2, 'f2 mpc 2 comp (set_mpc1 second direction)')
  call assert_close(config%mpcs(2)%coeff, 1.0_rk, 'f2 mpc 2 coefficient 1')
  call assert_equal(int(config%mpcs(2)%line), 23, 'f2 mpc 2 line')
  call assert_equal(config%mpcs(3)%nset, 'y1', 'f2 mpc 3 nset')
  call assert_equal(config%mpcs(3)%comp, 3, 'f2 mpc 3 comp')
  call assert_close(config%mpcs(3)%coeff, -0.5_rk, 'f2 mpc 3 coefficient')
  call assert_close(config%mpcs(3)%offset, 2.0e-3_rk, 'f2 mpc 3 offset')
  call assert_equal(config%mpcs(4)%comp, 1, 'f2 mpc 4 comp')
  call assert_close(config%mpcs(4)%offset, 0.0_rk, 'f2 mpc 4 offset')
  call assert_true(config%has_target_strain(3, 1) .and. &
    & .not. any(config%has_target_strain(3, 2:3)), 'f2 strain target in step 1 only')
  call assert_true(.not. any(config%has_strain_rate), 'f2 no strain rate')
  call assert_true(.not. config%has_target_time(1) .and. all(config%has_target_time(2:3)), &
    & 'f2 time targets in steps 2 and 3')
  call assert_close(config%target_time, [0.0_rk, 2.0_rk, 3.5_rk], 'f2 time targets')
  call assert_true(all(config%incr_kind == [incr_by_count, incr_by_dtime, incr_by_dstrain]), &
    & 'f2 increments by count, dtime, dstrain')
  call assert_equal(config%num_incrs(1), 1, 'f2 number_of_incrs step 1')
  call assert_close(config%dtime, [0.0_rk, 0.5_rk, 0.0_rk], 'f2 dtime')
  call assert_close(config%dstrain, [0.0_rk, 0.0_rk, 2.0e-4_rk], 'f2 dstrain')
  call assert_true(config%print_results(1) .and. .not. config%print_results(2) .and. &
    & config%print_results(3), 'f2 print_results')
  call assert_equal(config%max_incr, 100, 'f2 max_incr')
  call assert_close(config%max_total_time, 10.0_rk, 'f2 max_total_time')
  call assert_true(config%has_max_strain .and. config%has_max_eqstrain, 'f2 strain limits given')
  call assert_close(config%max_strain, 0.1_rk, 'f2 max_strain')
  call assert_close(config%max_eqstrain, 0.15_rk, 'f2 max_eqstrain')
  call assert_close(config%cg_tol, 1.0e-10_rk, 'f2 cg_tol')
  call assert_equal(config%cg_max_iters, 500, 'f2 cg_max_iters')
  call assert_equal(config%global_max_iters, 7, 'f2 global_max_iters')
  call assert_close(config%global_tol, 1.0e-6_rk, 'f2 global_tol')
  call assert_close(config%global_tol_abs, 1.0e-12_rk, 'f2 global_tol_abs')
  call assert_true(config%midnode_linearization, 'f2 midnode_linearization')
  call assert_true(config%subincrementation, 'f2 subincrementation')
  call assert_equal(config%subincrement_max_depth, 6, 'f2 subincrement_max_depth')
  call assert_equal(config%subincrement_max_successive, 3, 'f2 subincrement_max_successive')
  call assert_true(config%read_ori_from_file .and. config%read_phase_from_file .and. &
    & config%read_from_opt_file, 'f2 input-file flags')
  call assert_equal(config%print_vars(1), 'forces', 'f2 print order 1')
  call assert_equal(config%print_vars(3), 'disp', 'f2 print order 3')

  ! --- File 3: uniaxial_grip on 11, then on 22 (loading axis from the component) ---
  do i = 1, 2
    call write_file('test_config_3.cfg', [character(len=64) :: &
      & 'number_of_phases 1', 'phase 1', 'crystal_type bcc', 'model elastic', &
      & 'c11 231400', 'c12 134700', 'c44 116400', &
      & 'boundary_conditions uniaxial_grip', 'number_of_steps 1', &
      & 'strain_rate ' // merge('11', '22', i == 1) // ' -1e-3', &
      & 'target_strain ' // merge('11', '22', i == 1) // ' -0.001', &
      & 'number_of_incrs 4'])
    call config_read(config, 'test_config_3.cfg')
    call assert_equal(config%control_axis, i, 'f3 load axis from the component')
    call assert_close(config%strain_rate(i, 1), -1.0e-3_rk, 'f3 negative rate')
    call assert_equal(size(config%print_vars), 0, 'f3 no print')
  end do

  ! --- File 4 (step 6.1): crystal_plasticity phases (D48) and an elastic one ---
  ! Phase 1 bcc, {112} first, per-family values; phase 2 hcp, default families, m per system;
  ! phase 3 fcc, one value each; phase 4 bct, all, per-family m and per-system gammadot_0;
  ! phase 5 hcp, pyramidal then basal; phase 6 elastic.
  write (text, '(a,18(1x,f5.3))') 'm', [(0.010_rk + 0.001_rk*i, i = 1, 18)]
  ! Phase 5's user interaction matrix (15 systems): H(a, b) = 1 on the diagonal, 0.5 + a/100 +
  ! b/1000 off it.
  do i = 1, 15
    write (rows(i), '(a,15(1x,f6.3))') 'interaction_matrix', [(merge(1.0_rk, 0.5_rk + 0.01_rk*i &
      & + 0.001_rk*j, i == j), j = 1, 15)]
  end do
  call write_file('test_config_4.cfg', [character(len=256) :: &
    & 'number_of_phases 6', &
    & 'phase 1', 'crystal_type bcc', 'model crystal_plasticity', &
    & 'c11 231400', 'c12 134700', 'c44 116400', &
    & 'slip_families {112}<111> {110}<111>', 'm 0.02 0.05', 'gammadot_0 1.0 2.0', &
    & 'tau_0 240 200', 'hardening saturation_anisotropic', 'h_0 200', 'tau_s 330 360', &
    & 'n 1.5', &
    & 'phase 2', 'crystal_type hcp', 'model crystal_plasticity', 'c11 162400', 'c12 92000', &
    & 'c13 69000', 'c33 180700', 'c44 46700', 'c_over_a 1.587', trim(text), &
    & 'gammadot_0 1.0 2.0 3.0', 'tau_0 100 120 300', 'hardening none', &
    & 'phase 3', 'crystal_type fcc', 'model crystal_plasticity', &
    & 'c11 168400', 'c12 121400', 'c44 75400', 'm 0.05', 'gammadot_0 1.0', 'tau_0 50', &
    & 'h_0 100', 'tau_s 150', 'n 1', &
    & 'phase 4', 'crystal_type bct', 'model crystal_plasticity', 'c11 72300', 'c12 59400', &
    & 'c13 35800', 'c33 88400', 'c44 22000', 'c66 24000', 'c_over_a 0.5456', &
    & 'slip_families all', 'm 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10', &
    & 'gammadot_0' // repeat(' 1.5', 32), 'tau_0' // repeat(' 20', 10), &
    & 'hardening kocks_mecking', 'km_b 2.48e-7', 'km_h_n 0.125', 'km_rho_0 2e7', 'km_k1 4.5e5', &
    & 'km_k20 14', 'km_gammadot_k0 1e10', 'km_n 100', &
    & 'phase 5', 'crystal_type hcp', 'model crystal_plasticity', 'c11 162400', 'c12 92000', &
    & 'c13 69000', 'c33 180700', 'c44 46700', 'c_over_a 1.587', &
    & 'slip_families {10-11}<11-23> {0001}<11-20>', 'm 0.03 0.01', 'gammadot_0 1 1', &
    & 'tau_0 300 100', 'hardening saturation_anisotropic', 'h_0 50 60', 'tau_s 400 200', &
    & 'n 1', rows, 'initialization elastic_step', 'point_tol 1e-6', &
    & 'point_max_iters_hardening 7', 'point_tol_state 1e-9', 'point_tol_hardening 2e-8', &
    & 'point_max_iters_newton 50', 'point_max_iters_state 30', 'global_tol_slow 1e-3', &
    & 'global_nr_switch 0', 'global_nr_slow_ratio 0.5', &
    & 'phase 6', 'crystal_type fcc', 'model elastic', &
    & 'c11 168400', 'c12 121400', 'c44 75400', &
    & 'boundary_conditions uniaxial_minimal', 'number_of_steps 1', 'strain_rate 33 1e-3', &
    & 'target_strain 33 0.001', 'number_of_incrs 1'])
  call config_read(config, 'test_config_4.cfg')
  associate (ph => config%phases)
    ! bcc, {112} first: systems 1-12 are {112} (0.02, 1.0), 13-24 {110} (0.05, 2.0).
    call assert_equal(ph(1)%model, 'crystal_plasticity', 'f4 bcc model')
    call assert_true(all(ph(1)%families == [2, 1]), 'f4 bcc families in the given order')
    call assert_equal(ph(1)%num_slip, 24, 'f4 bcc num_slip')
    call assert_close(ph(1)%m, [(0.02_rk, i = 1, 12), (0.05_rk, i = 1, 12)], 'f4 bcc m')
    call assert_close(ph(1)%gammadot_0, [(1.0_rk, i = 1, 12), (2.0_rk, i = 1, 12)], &
      & 'f4 bcc gammadot_0')
    call assert_close(ph(1)%tau_0, [(240.0_rk, i = 1, 12), (200.0_rk, i = 1, 12)], &
      & 'f4 bcc tau_0 in the slip_families order')
    call assert_equal(ph(1)%hardening, 'saturation_anisotropic', 'f4 bcc saturation_anisotropic')
    call assert_close(ph(1)%h_0, [(200.0_rk, i = 1, 24)], 'f4 bcc h_0 one value for all')
    call assert_close(ph(1)%tau_s, [(330.0_rk, i = 1, 12), (360.0_rk, i = 1, 12)], &
      & 'f4 bcc tau_s per family, slip_families order')
    call assert_close(ph(1)%n_hard, [(1.5_rk, i = 1, 24)], 'f4 bcc n')
    call assert_true(all(shape(ph(1)%interaction) == [24, 24]), 'f4 bcc default matrix shape')
    call assert_true(all([((abs(ph(1)%interaction(i, j) - merge(1.0_rk, 1.4_rk, i == j)) &
      & < 1.0e-15_rk, i = 1, 24), j = 1, 24)]), 'f4 bcc default matrix 1 / 1.4')
    call assert_equal(ph(2)%hardening, 'none', 'f4 hcp hardening none')
    call assert_equal(ph(3)%hardening, 'saturation_isotropic', 'f4 fcc default hardening')
    call assert_close([ph(3)%h_0(12), ph(3)%tau_s(12), ph(3)%n_hard(12)], [100.0_rk, 150.0_rk, &
      & 1.0_rk], 'f4 fcc isotropic values for every system')
    call assert_true(size(ph(3)%interaction) == 0, 'f4 fcc isotropic: no matrix')
    call assert_close(ph(5)%h_0, [(50.0_rk, i = 1, 12), (60.0_rk, i = 1, 3)], 'f4 hcp h_0')
    call assert_true(all([((abs(ph(5)%interaction(i, j) - merge(1.0_rk, 0.5_rk + 0.01_rk*i + &
      & 0.001_rk*j, i == j)) < 1.0e-12_rk, i = 1, 15), j = 1, 15)]), &
      & 'f4 hcp user matrix, rows in file order')
    call assert_equal(ph(6)%hardening, '', 'f4 elastic: no hardening')
    ! hcp, default families: m per system as given, gammadot_0 per family.
    call assert_true(all(ph(2)%families == [1, 2, 3]), 'f4 hcp default families')
    call assert_equal(ph(2)%num_slip, 18, 'f4 hcp num_slip')
    call assert_close(ph(2)%m, [(0.010_rk + 0.001_rk*i, i = 1, 18)], 'f4 hcp m per system', &
      & abs_tol=1.0e-12_rk)
    call assert_close(ph(2)%gammadot_0, [(1.0_rk, i = 1, 3), (2.0_rk, i = 1, 3), &
      & (3.0_rk, i = 1, 12)], 'f4 hcp gammadot_0 per family')
    ! fcc, one family: a single value each.
    call assert_true(all(ph(3)%families == [1]), 'f4 fcc families')
    call assert_close(ph(3)%m, [(0.05_rk, i = 1, 12)], 'f4 fcc m')
    ! bct, all: per-family m over the v2 family ranges, per-system gammadot_0.
    call assert_true(all(ph(4)%families == [(i, i = 1, 10)]), 'f4 bct all families')
    call assert_equal(ph(4)%num_slip, 32, 'f4 bct num_slip')
    call assert_close(ph(4)%m, 0.01_rk*[1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, &
      & 7, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10], 'f4 bct m per family', &
      & abs_tol=1.0e-15_rk)
    call assert_close(ph(4)%gammadot_0, [(1.5_rk, i = 1, 32)], 'f4 bct gammadot_0')
    ! bct, Kocks-Mecking (step 7.2): its seven parameters; no saturation-law values.
    call assert_equal(ph(4)%hardening, 'kocks_mecking', 'f4 bct kocks_mecking')
    call assert_close([ph(4)%km_b, ph(4)%km_h_n, ph(4)%km_rho_0, ph(4)%km_k1, ph(4)%km_k20, &
      & ph(4)%km_gammadot_k0, ph(4)%km_n], [2.48e-7_rk, 0.125_rk, 2.0e7_rk, 4.5e5_rk, 14.0_rk, &
      & 1.0e10_rk, 100.0_rk], 'f4 bct Kocks-Mecking parameters', rel_tol=1.0e-15_rk)
    call assert_true(size(ph(4)%h_0) == 0 .and. size(ph(4)%tau_s) == 0, &
      & 'f4 bct Kocks-Mecking: no saturation parameters')
    ! hcp, pyramidal then basal: 12 + 3 systems.
    call assert_true(all(ph(5)%families == [3, 1]), 'f4 hcp subset order')
    call assert_equal(ph(5)%num_slip, 15, 'f4 hcp subset num_slip')
    call assert_close(ph(5)%m, [(0.03_rk, i = 1, 12), (0.01_rk, i = 1, 3)], 'f4 hcp subset m')
    call assert_close(ph(5)%tau_0, [(300.0_rk, i = 1, 12), (100.0_rk, i = 1, 3)], &
      & 'f4 hcp subset tau_0')
    call assert_close(ph(2)%tau_0, [(100.0_rk, i = 1, 3), (120.0_rk, i = 1, 3), &
      & (300.0_rk, i = 1, 12)], 'f4 hcp tau_0 per family')
    ! elastic: no slip.
    call assert_equal(ph(6)%num_slip, 0, 'f4 elastic num_slip')
    call assert_true(size(ph(6)%families) == 0 .and. size(ph(6)%m) == 0 .and. &
      & size(ph(6)%gammadot_0) == 0 .and. size(ph(6)%tau_0) == 0, 'f4 elastic empty slip arrays')
  end associate
  call assert_equal(config%initialization, 'elastic_step', 'f4 initialization')
  call assert_close(config%point_tol, 1.0e-6_rk, 'f4 point_tol')
  call assert_equal(config%point_max_iters_newton, 50, 'f4 point_max_iters_newton')
  call assert_equal(config%point_max_iters_state, 30, 'f4 point_max_iters_state')
  call assert_equal(config%point_max_iters_hardening, 7, 'f4 point_max_iters_hardening')
  call assert_close(config%point_tol_state, 1.0e-9_rk, 'f4 point_tol_state')
  call assert_close(config%point_tol_hardening, 2.0e-8_rk, 'f4 point_tol_hardening')
  call assert_close(config%global_tol_slow, 1.0e-3_rk, 'f4 global_tol_slow')
  call assert_close(config%global_nr_switch, 0.0_rk, 'f4 global_nr_switch 0 (SA only)')
  call assert_close(config%global_nr_slow_ratio, 0.5_rk, 'f4 global_nr_slow_ratio')

  ! --- File 5 (step 8.2b): triaxial load control, load and stress rows on other components ---
  ! Control 11 at a strain rate to a load target; 33 load-controlled (target_load), 22 by a
  ! stress (target_stress); load_tol and stress_tol by the rows of other components; defaults of
  ! min_pert_frac and max_bc_iter.
  call write_file('test_config_5.cfg', [character(len=64) :: &
    & 'number_of_phases 1', 'phase 1', '  crystal_type fcc', '  model elastic', &
    & '  c11 200000', '  c12 100000', '  c44 50000', &
    & 'boundary_conditions triaxial', 'control_component 11', 'number_of_steps 2', &
    & 'strain_rate 11 -1e-4 -1e-4', 'target_load 11 -300 -600', &
    & 'target_stress 22 0 -10', 'target_load 33 300 600', 'dtime 5 10', &
    & 'dtime_min 0.25 0.25', 'stress_tol 0.05'])
  call config_read(config, 'test_config_5.cfg')
  call assert_equal(config%control_axis, 1, 'f5 control 11')
  call assert_true(all(config%comp_rows), 'f5 rows on every component')
  call assert_true(all(config%has_target_load(3, :)) .and. all(config%has_target_stress(2, :)) &
    & .and. .not. any(config%has_strain_rate(2:3, :)), 'f5 load and stress rows of 22, 33')
  call assert_close(config%target_load(3, :), [300.0_rk, 600.0_rk], 'f5 target_load 33')
  call assert_close(config%target_stress(2, :), [0.0_rk, -10.0_rk], 'f5 target_stress 22')
  call assert_close(config%stress_tol, 0.05_rk, 'f5 stress_tol')
  call assert_close(config%load_tol, 0.1_rk, 'f5 default load_tol')
  call assert_close(config%min_pert_frac, 0.001_rk, 'f5 default min_pert_frac (v2)')
  call assert_equal(config%max_bc_iter, 10, 'f5 default max_bc_iter (v2)')

  ! --- File 6 (step 8.2c): load and stress rates, durations, dwells ---
  ! Control 33 at a load rate to 100 N, a dwell (load_rate 0, target_duration), at a stress rate
  ! to 0 MPa; 11 by its stress rate, 22 by a load target; max_strain_incr and min_pert_frac with
  ! a dwell and the load-controlled control component.
  call write_file('test_config_6.cfg', [character(len=64) :: &
    & 'number_of_phases 1', 'phase 1', '  crystal_type fcc', '  model elastic', &
    & '  c11 200000', '  c12 100000', '  c44 50000', &
    & 'boundary_conditions triaxial', 'control_component 33', 'number_of_steps 3', &
    & 'load_rate 33 2.0 0.0 -', 'stress_rate 33 - - -5.0', 'target_load 33 100.0 - -', &
    & 'target_stress 33 - - 0.0', 'target_duration - 20.0 -', 'stress_rate 11 -1.0 0.0 1.0', &
    & 'target_load 22 0.0 0.0 0.0', 'number_of_incrs 2 - 2', 'dtime - 5.0 -', &
    & 'max_strain_incr 0.002', 'min_pert_frac 0.01'])
  call config_read(config, 'test_config_6.cfg')
  call assert_true(all(config%rate_kind == [rate_by_load, rate_by_load, rate_by_stress]), &
    & 'f6 control rates: load, load (a dwell), stress')
  call assert_true(all(config%target_kind == [target_by_load, target_by_duration, &
    & target_by_stress]), 'f6 targets: load, duration, stress')
  call assert_close(config%load_rate(3, 1:2), [2.0_rk, 0.0_rk], 'f6 load_rate 33')
  call assert_close(config%stress_rate(3, 3), -5.0_rk, 'f6 stress_rate 33')
  call assert_close(config%target_duration(2), 20.0_rk, 'f6 target_duration')
  call assert_true(all(config%has_stress_rate(1, :)) .and. all(config%has_target_load(2, :)), &
    & 'f6 rows of 11 and 22')
  call assert_close(config%max_strain_incr, 0.002_rk, 'f6 max_strain_incr')
  call assert_close(config%min_pert_frac, 0.01_rk, 'f6 min_pert_frac')
  call assert_true(.not. any(rate_by_strain == config%rate_kind) .and. &
    & config%target_kind(1) /= target_by_stress, 'f6 no strain rate on 33')

  ! --- Files 7 and 8 (step 9.2, Q-B1, Q-B2): the default preconditioner with constraints ---
  ! Periodic conditions and an MPC, cg_preconditioner absent: jacobi (block_jacobi, the default
  ! without constraints, is not available with them).
  call write_file('test_config_7.cfg', [character(len=64) :: &
    & 'number_of_phases 1', 'phase 1', '  crystal_type fcc', '  model elastic', &
    & '  c11 200000', '  c12 100000', '  c44 50000', &
    & 'boundary_conditions periodic', 'control_component 33', 'number_of_steps 1', &
    & 'strain_rate 33 1.0e-3', 'target_strain 33 0.001', 'number_of_incrs 2'])
  call config_read(config, 'test_config_7.cfg')
  call assert_equal(config%cg_preconditioner, 'jacobi', &
    & 'f7 default cg_preconditioner with periodic conditions: jacobi')
  call write_file('test_config_8.cfg', [character(len=64) :: &
    & 'number_of_phases 1', 'phase 1', '  crystal_type fcc', '  model elastic', &
    & '  c11 200000', '  c12 100000', '  c44 50000', &
    & 'boundary_conditions uniaxial_minimal', 'number_of_steps 1', &
    & 'strain_rate 33 1.0e-3', 'target_strain 33 0.001', 'number_of_incrs 2', &
    & 'set_mpc1 vel x1 x'])
  call config_read(config, 'test_config_8.cfg')
  call assert_equal(config%cg_preconditioner, 'jacobi', &
    & 'f8 default cg_preconditioner with an MPC: jacobi')

  call test_summary()

contains

  !> Write lines to a file (trailing blanks removed).
  !>
  !> Arguments:
  !>   name   (in)  character  file name
  !>   text   (in)  character(*)(:)  the lines
  !> Locals:
  !>   u, i   integer          unit, line index
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

  !> Check velocity condition i of the configuration.
  !>
  !> Arguments:
  !>   i_bc   (in)  integer    condition index
  !>   nset   (in)  character  node set
  !>   comp   (in)  integer    direction
  !>   value  (in)  real(rk)   velocity [mm/s]
  !>   line   (in)  integer    line in the file
  !> Locals:
  !>   label  character(8)     check label prefix
  subroutine check_bc(i_bc, nset, comp, value, line)

    integer, intent(in) :: i_bc
    character(len=*), intent(in) :: nset
    integer, intent(in) :: comp
    real(rk), intent(in) :: value
    integer, intent(in) :: line

    character(len=8) :: label

    write (label, '(a,i0)') 'f2 bc ', i_bc
    call assert_equal(config%vel_bcs(i_bc)%nset, nset, trim(label) // ' nset')
    call assert_equal(config%vel_bcs(i_bc)%comp, comp, trim(label) // ' comp')
    call assert_close(config%vel_bcs(i_bc)%value, value, trim(label) // ' value')
    call assert_equal(int(config%vel_bcs(i_bc)%line), line, trim(label) // ' line')

  end subroutine check_bc

end program test_config
