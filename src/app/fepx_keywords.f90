! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> FEPX keyword documentation generator: writes the manual's generated keyword pages from the
!> keyword table (cfg_table_mod) and the option-compatibility rules (cfg_rules_mod), through
!> cfg_manual_mod: the keyword reference (keywords_table.rst), the compatibility table
!> (compatibility_table.rst), and the facts and notes of every implemented keyword's entry
!> (generated/keywords/<name>.inc); or the list of the rules, for the docs test of their use.
!>
!> Usage: fepx_keywords <directory>       write the pages into <directory> (the manual,
!>                                        doc; <directory>/generated/keywords must
!>                                        exist); existing files are replaced
!>        fepx_keywords --rules <file>    write one line per rule to <file>: "<name>|<message
!>                                        for the subject <subject>>" (cfg_rule_message)
!>
!> Source : FEPX (ROADMAP task M.1, step 9.9, D104). Contains no numerical method. Serial (no
!>          MPI).
!>
!> Locals:
!>   arg     character(1024)  first command argument
!>   path    character(1024)  output directory or file
!>   nargs   integer          number of command arguments
!>   unit    integer          output unit
!>   i       integer          keyword or rule index
program fepx_keywords

  use error_mod, only: fatal_error
  use cfg_table_mod, only: cfg_keywords, cfg_active
  use cfg_rules_mod, only: cfg_rules, cfg_rule_message
  use cfg_manual_mod, only: cfg_manual_keyword, cfg_manual_reference, cfg_manual_rules

  implicit none

  character(len=1024) :: arg, path
  integer :: nargs, unit, i

  ! One argument (the manual directory), or --rules and the output file.
  nargs = command_argument_count()
  if (nargs < 1) call fatal_error('fepx_keywords', 'usage: fepx_keywords <directory> | ' // &
    & '--rules <file>')
  call get_command_argument(1, arg)

  if (trim(arg) == '--rules') then
    ! The rules: name and message, the subject being the first keyword of the rule.
    if (nargs /= 2) call fatal_error('fepx_keywords', 'usage: fepx_keywords --rules <file>')
    call get_command_argument(2, path)
    call open_output(trim(path), unit)
    do i = 1, size(cfg_rules)
      write (unit, '(a)') trim(cfg_rules(i)%name) // '|' // cfg_rule_message(i, 'subject')
    end do
    close (unit)
  else
    ! The manual's pages.
    if (nargs /= 1) call fatal_error('fepx_keywords', 'usage: fepx_keywords <directory>')
    path = arg
    call open_output(trim(path) // '/keywords_table.rst', unit)
    call cfg_manual_reference(unit)
    close (unit)
    call open_output(trim(path) // '/compatibility_table.rst', unit)
    call cfg_manual_rules(unit)
    close (unit)
    do i = 1, size(cfg_keywords)
      if (cfg_keywords(i)%status /= cfg_active) cycle
      call open_output(trim(path) // '/generated/keywords/' // trim(cfg_keywords(i)%name) // &
        & '.inc', unit)
      call cfg_manual_keyword(unit, i)
      close (unit)
    end do
  end if

contains

  !> Open an output file (replaced if it exists); a fatal error if it cannot be opened.
  !>
  !> Arguments:
  !>   file_name  (in)   character  the file
  !>   unit       (out)  integer    its unit
  !> Locals:
  !>   ios        integer           I/O status of the open
  subroutine open_output(file_name, unit)

    character(len=*), intent(in) :: file_name
    integer, intent(out) :: unit

    integer :: ios

    open (newunit=unit, file=file_name, status='replace', action='write', iostat=ios)
    if (ios /= 0) call fatal_error('fepx_keywords', 'cannot open ' // file_name)

  end subroutine open_output

end program fepx_keywords
