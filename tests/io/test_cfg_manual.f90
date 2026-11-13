! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of the manual's generated keyword pages (cfg_manual_mod): the values of keywords of
!> every kind and form, the table's values_text, the per-family parameters (and that they are
!> list keywords), the slip families of every crystal type; the entry of a keyword with several rules (a bulleted note), with one
!> rule (a plain note), in a phase block, and without a unit; the keyword reference (one row per
!> implemented keyword, no planned keyword) and the compatibility table (one row per rule).
!>
!> Source : FEPX (ROADMAP step 9.9, D104). The committed pages are checked against the generator
!>          by the docs test docs_keywords_current.
!>
!> Locals:
!>   unit      integer          file unit
!>   ios       integer          read status
!>   i         integer          keyword-table index
!>   n_rows    integer          rows of a generated table (lines starting a list-table row)
!>   n_active  integer          implemented keywords
!>   n_bullets integer          bulleted lines of an entry's note
!>   found     logical(5)       expected lines found in an entry
!>   planned   logical          a planned keyword found in the reference
!>   buffer    character(1024)  a line read back
program test_cfg_manual

  use cfg_table_mod, only: cfg_keywords, cfg_find_keyword, cfg_active, cfg_planned, cfg_list, &
    & cfg_real, cfg_phase
  use cfg_rules_mod, only: cfg_rules
  use cfg_manual_mod, only: cfg_manual_keyword, cfg_manual_reference, cfg_manual_rules, &
    & cfg_manual_values, per_family_names, saturation_names
  use assert_mod, only: assert_equal, assert_true, test_summary

  implicit none

  integer :: unit, ios, i, n_rows, n_active, n_bullets
  logical :: found(5), planned
  character(len=1024) :: buffer

  ! Values: kinds, ranges and forms.
  call assert_equal(cfg_manual_values(cfg_keywords(cfg_find_keyword('cg_preconditioner'))), &
    & 'one of ``jacobi``, ``block_jacobi``, ``none``', 'values, words')
  call assert_equal(cfg_manual_values(cfg_keywords(cfg_find_keyword('number_of_phases'))), &
    & 'integer >= 1', 'values, integer with a lower bound')
  call assert_equal(cfg_manual_values(cfg_keywords(cfg_find_keyword('c_over_a'))), &
    & 'real > 0', 'values, real with an open lower bound')
  call assert_equal(cfg_manual_values(cfg_keywords(cfg_find_keyword('print_results'))), &
    & 'integer from 0 to 1, one value per step', 'values, bounded, per step')
  call assert_equal(cfg_manual_values(cfg_keywords(cfg_find_keyword('target_time'))), &
    & 'real > 0, one value per step (``-`` for a step without a value)', &
    & 'values, per step with dashes')
  call assert_equal(cfg_manual_values(cfg_keywords(cfg_find_keyword('strain_rate'))), &
    & 'a component, then real, one value per step (``-`` for a step without a value)', &
    & 'values, component per step')
  call assert_equal(cfg_manual_values(cfg_keywords(cfg_find_keyword('read_ori_from_file'))), &
    & 'none (the keyword alone)', 'values, flag')
  call assert_equal(cfg_manual_values(cfg_keywords(cfg_find_keyword('output_format'))), &
    & 'one of ``ascii``, ``hdf5``, one or more values', 'values, list of words')
  call assert_true(index(cfg_manual_values(cfg_keywords(cfg_find_keyword('set_bc'))), &
    & '``vel <nset> <d> <v> [<d> <v> [<d> <v>]]``') == 1, 'values, the table''s values_text')
  call assert_true(index(cfg_manual_values(cfg_keywords(cfg_find_keyword('tau_0'))), &
    & 'real > 0: one value per active slip family, in the order of ``slip_families``') == 1, &
    & 'values, per family')
  call assert_true(index(cfg_manual_values(cfg_keywords(cfg_find_keyword('h_0'))), &
    & 'real >= 0: one value; with ``hardening saturation_anisotropic``') == 1, &
    & 'values, per family with saturation_anisotropic')

  ! The per-family keywords are list keywords of real values in a phase block.
  call check_names(per_family_names)
  call check_names(saturation_names)

  ! The slip families: the families of every crystal type, and the defaults (crystal_mod).
  call write_entry('slip_families')
  call read_entry([character(len=64) :: '   - ``fcc``: ``{111}<110>`` (12 systems)', &
    & '   - ``bcc``: ``{110}<111>`` (12 systems), ``{112}<111>``', &
    & '   - ``hcp``: ``{0001}<11-20>`` (3 systems)', &
    & '   :Default: ``fcc``: all families; ``bcc``: ``{110}<111>``;', &
    & '   - ``bct``: ``{100}<001>`` (2 systems)'])
  call assert_true(all(found), 'slip_families: families of every crystal type, defaults')

  ! A run-wide keyword with two rules: its facts and a bulleted note.
  call write_entry('cg_preconditioner')
  call read_entry([character(len=64) :: '   :Values: one of ``jacobi``, ``block_jacobi``, ``none``', &
    & '   :Default: ``block_jacobi``                               ', &
    & '   :Given: once, for the whole run                          ', &
    & '   .. note::                                                ', &
    & '   :Unit:                                                   '])
  call assert_true(all(found(1:4)), 'cg_preconditioner: values, default, given, note')
  call assert_true(.not. found(5), 'cg_preconditioner: no unit')
  call assert_equal(n_bullets, 2, 'cg_preconditioner: two rules')

  ! A phase keyword with one rule and no unit ("-"): a plain note.
  call write_entry('c_over_a')
  call read_entry([character(len=64) :: '   :Values: real > 0                                       ', &
    & '   :Default: required (hcp, bct plasticity)                 ', &
    & '   :Given: in each phase block (``phase``)                  ', &
    & '   ``c_over_a`` is used only with ``crystal_type hcp`` and  ', &
    & '   :Unit:                                                   '])
  call assert_true(all(found(1:4)), 'c_over_a: values, default, given, the rule')
  call assert_true(.not. found(5), 'c_over_a: no unit')
  call assert_equal(n_bullets, 0, 'c_over_a: one rule, not bulleted')

  ! A keyword with an open and a closed bound and no rule: no note.
  call write_entry('cg_tol')
  call read_entry([character(len=64) :: '   :Given: once, for the whole run                          ', &
    & '   .. note::                                                ', &
    & '   :Unit:                                                   ', &
    & '   :Values: real > 0 and <= 1                               ', &
    & '   :Default: ``1e-8``                                       '])
  call assert_true(found(1) .and. found(4) .and. found(5), 'cg_tol: given, values, default')
  call assert_true(.not. found(2), 'cg_tol: no note')

  ! The keyword reference: one row per implemented keyword, the first number_of_phases, no
  ! planned keyword.
  n_active = count(cfg_keywords%status == cfg_active)
  open (newunit=unit, file='test_cfg_manual_reference.rst', status='replace', action='write')
  call cfg_manual_reference(unit)
  close (unit)
  open (newunit=unit, file='test_cfg_manual_reference.rst', status='old', action='read')
  n_rows = 0
  planned = .false.
  do
    read (unit, '(a)', iostat=ios) buffer
    if (ios /= 0) exit
    if (index(buffer, '   * - :confval:`') == 1) then
      n_rows = n_rows + 1
      if (n_rows == 1) call assert_true(index(buffer, '`number_of_phases`') > 0, &
        & 'reference: first row')
    end if
    do i = 1, size(cfg_keywords)
      if (cfg_keywords(i)%status == cfg_planned .and. index(buffer, ':confval:`' // &
        & trim(cfg_keywords(i)%name) // '`') > 0) planned = .true.
    end do
  end do
  close (unit)
  call assert_equal(n_rows, n_active, 'reference: one row per implemented keyword')
  call assert_true(.not. planned, 'reference: no planned keyword')

  ! The compatibility table: one row per rule.
  open (newunit=unit, file='test_cfg_manual_rules.rst', status='replace', action='write')
  call cfg_manual_rules(unit)
  close (unit)
  open (newunit=unit, file='test_cfg_manual_rules.rst', status='old', action='read')
  n_rows = 0
  do
    read (unit, '(a)', iostat=ios) buffer
    if (ios /= 0) exit
    if (index(buffer, '   * - :confval:`') == 1) n_rows = n_rows + 1
  end do
  close (unit)
  call assert_equal(n_rows, size(cfg_rules), 'compatibility table: one row per rule')

  call test_summary()

contains

  !> Check that every keyword of a blank-separated list is a list keyword of real values given
  !> in a phase block (the per-family parameters).
  !>
  !> Arguments:
  !>   names  (in)  character  the keywords
  !> Locals:
  !>   rest   character(:), allocatable  the keywords not yet checked
  !>   word   character(:), allocatable  the current keyword
  !>   blank  integer                    position of the next blank in rest
  !>   k      integer                    keyword-table index
  subroutine check_names(names)

    character(len=*), intent(in) :: names

    character(len=:), allocatable :: rest, word
    integer :: blank, k

    rest = trim(adjustl(names))
    do while (len(rest) > 0)
      ! The next keyword, and the rest.
      blank = index(rest, ' ')
      if (blank == 0) then
        word = rest
        rest = ''
      else
        word = rest(1:blank - 1)
        rest = trim(adjustl(rest(blank + 1:)))
      end if
      k = cfg_find_keyword(word)
      call assert_true(k > 0, 'per-family keyword ' // word // ' in the table')
      if (k > 0) call assert_true(cfg_keywords(k)%form == cfg_list .and. &
        & cfg_keywords(k)%kind == cfg_real .and. cfg_keywords(k)%scope == cfg_phase, &
        & 'per-family keyword ' // word // ': a list of reals in a phase block')
    end do

  end subroutine check_names

  !> Write the entry of a keyword to test_cfg_manual_entry.rst.
  !>
  !> Arguments:
  !>   name  (in)  character  keyword
  !> Locals:
  !>   out   integer          file unit
  subroutine write_entry(name)

    character(len=*), intent(in) :: name

    integer :: out

    open (newunit=out, file='test_cfg_manual_entry.rst', status='replace', action='write')
    call cfg_manual_keyword(out, cfg_find_keyword(name))
    close (out)

  end subroutine write_entry

  !> Read test_cfg_manual_entry.rst back: found(j) when a line starts with the text of
  !> expected(j) (blank-trimmed, the indentation of three blanks removed for the field lines
  !> and the note directive, kept for the note's body); n_bullets counts the lines "   - ".
  !>
  !> Arguments:
  !>   expected  (in)  character(*)(5)  the line starts looked for
  !> Locals:
  !>   in        integer                file unit
  !>   j         integer                expected-line index
  !>   status    integer                read status
  !>   line      character(1024)        a line read back
  !>   start     character(:), allocatable  expected(j) as it starts a line
  subroutine read_entry(expected)

    character(len=*), intent(in) :: expected(5)

    integer :: in, j, status
    character(len=1024) :: line
    character(len=:), allocatable :: start

    found = .false.
    n_bullets = 0
    open (newunit=in, file='test_cfg_manual_entry.rst', status='old', action='read')
    do
      read (in, '(a)', iostat=status) line
      if (status /= 0) exit
      if (index(line, '   - ') == 1) n_bullets = n_bullets + 1
      do j = 1, size(expected)
        ! Field lines and the directive start in column 1; the note's body is indented.
        start = trim(expected(j))
        if (start(4:4) == ':' .or. start(4:5) == '..') start = start(4:)
        if (index(line, start) == 1) found(j) = .true.
      end do
    end do
    close (in)

  end subroutine read_entry

end program test_cfg_manual
