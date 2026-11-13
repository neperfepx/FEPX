! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit test of log_mod: line layout (prefixes, depths, file markers), wrapping at 72 columns,
!> rank selection, the all-ranks buffer, the hold mode, the warning counter and the text helpers.
!>
!> Source : FEPX (ROADMAP steps 2.1, 9.7, D94). Expected lines follow
!>          dev/plan/terminal_output.md §3 (v2's and Neper's prefixes, depth markers "- ",
!>          "> ", ". ", file markers "[i]", "[o]", 72 columns, continuation lines aligned under
!>          the text, two columns in at depth 0).
!>
!> The log is redirected to a scratch file (log_set_unit), which is then read back line by line.
!> Checks: depths 0 to 3 and the file markers; a long message wrapped at blanks, each line at
!> most 72 columns, the continuation aligned under the text (depth 2) or two columns in (depth
!> 0); a word longer than a line kept whole; another rank (3) writes nothing, its all_ranks
!> lines are buffered with its rank tag and handed over by log_buffer_take; the hold mode keeps
!> rank 0's lines and hands them over; the warnings counted; the text helpers; nothing else is
!> written.
!>
!> Locals:
!>   unit_log   integer           scratch unit the log is written to
!>   ios        integer           iostat of the reads
!>   line       character(256)    one line read back
!>   long_text  character(:), allocatable  a message longer than a line
!>   taken      character(:), allocatable  buffered or held lines handed over
!>   expected   character(100)(15)  the expected lines, in order
!>   i_line     integer           line index
program test_log

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: ik, rk
  use log_mod, only: log_info, log_warning, log_file, log_text, log_wrap, log_buffer_take, &
    & log_lines_write, log_warning_count, log_set_rank, log_get_rank, log_set_unit, &
    & log_real_text, log_int_text, log_count_text, log_hold_start, log_hold_take, log_width
  use assert_mod, only: assert_equal, assert_true, test_summary

  implicit none

  integer :: unit_log
  integer :: ios
  character(len=256) :: line
  character(len=:), allocatable :: long_text, taken
  character(len=100) :: expected(15)
  integer :: i_line

  ! Redirect the log to a scratch file.
  open (newunit=unit_log, status='scratch', action='readwrite', form='formatted', iostat=ios)
  call assert_equal(ios, 0, 'open scratch file')
  call log_set_unit(unit_log)

  ! Rank 0: depths, file markers, a rule, a warning (counted).
  call log_info('Loading simulation...')
  call log_file('i', "Parsing file `simulation.cfg'...")
  call log_info('Material parameters:', 1)
  call log_info('Phase 1: fcc, crystal_plasticity', 2)
  call log_info('c11: 2.4500E+05 MPa', 3)
  call log_warning('Something to note.', 2)
  call log_text(repeat('=', 72))

  ! A long message at depth 2 and at depth 0 (wrapped), and a word longer than a line.
  long_text = 'The mesh has 8 partitions, but 4 ranks were started: the elements and ' // &
    & 'nodes are split in even blocks instead (partitioning even); performance may be lower.'
  call log_warning(long_text, 2)
  call log_info('A depth-0 message long enough to need a second line of the log, to show ' // &
    & 'its continuation two columns in.')
  call log_info(repeat('x', 80))

  ! Expected lines, in order.
  expected(1) = 'Info   : Loading simulation...'
  expected(2) = "Info   :   [i] Parsing file `simulation.cfg'..."
  expected(3) = 'Info   :   - Material parameters:'
  expected(4) = 'Info   :     > Phase 1: fcc, crystal_plasticity'
  expected(5) = 'Info   :       . c11: 2.4500E+05 MPa'
  expected(6) = 'Warning:     > Something to note.'
  expected(7) = repeat('=', 72)
  expected(8) = 'Warning:     > The mesh has 8 partitions, but 4 ranks were started: the'
  expected(9) = 'Warning:       elements and nodes are split in even blocks instead'
  expected(10) = 'Warning:       (partitioning even); performance may be lower.'
  expected(11) = 'Info   : A depth-0 message long enough to need a second line of the log,'
  expected(12) = 'Info   :   to show its continuation two columns in.'
  expected(13) = 'Info   : ' // repeat('x', 80)

  ! Read the log back and compare line by line; every wrapped line within the width.
  rewind (unit_log)
  do i_line = 1, 13
    read (unit_log, '(a)', iostat=ios) line
    call assert_equal(ios, 0, 'read line')
    call assert_equal(trim(line), trim(expected(i_line)), 'log line')
    if (i_line /= 13) call assert_true(len_trim(line) <= log_width, 'line within 72 columns')
  end do
  read (unit_log, '(a)', iostat=ios) line
  call assert_true(ios /= 0, 'no extra line')
  call assert_equal(log_warning_count(), 2, 'warnings counted on rank 0')

  ! Rank 3: nothing written; its all_ranks lines buffered with the rank, then handed over.
  rewind (unit_log)
  call log_set_rank(3_ik)
  call assert_true(log_get_rank() == 3_ik, 'rank set')
  call log_info('not written by rank 3')
  call log_warning('not written by rank 3 either')
  call log_info('Rank-specific information.', 1, all_ranks=.true.)
  call log_warning('Rank-specific warning.', 2, all_ranks=.true.)
  call assert_equal(log_warning_count(), 3, 'all_ranks warning counted on its rank')
  call log_lines_write('not written by rank 3 (lines)' // new_line('a'))
  call log_buffer_take(taken)
  call assert_equal(taken, 'Info   :   - [Rank 3] Rank-specific information.' // new_line('a') &
    & // 'Warning:     > [Rank 3] Rank-specific warning.' // new_line('a'), 'buffered lines')
  call log_buffer_take(taken)
  call assert_equal(len(taken), 0, 'buffer emptied')
  call log_set_rank(0_ik)
  endfile (unit_log)
  rewind (unit_log)
  read (unit_log, '(a)', iostat=ios) line
  call assert_true(ios /= 0, 'nothing written by rank 3')

  ! Hold mode: rank 0's lines kept, then handed over; written after the hold.
  rewind (unit_log)
  call log_hold_start()
  call log_info('Held line.', 2)
  call log_hold_take(taken)
  call assert_equal(taken, 'Info   :     > Held line.' // new_line('a'), 'held lines')
  call log_info('After the hold.')
  endfile (unit_log)
  rewind (unit_log)
  read (unit_log, '(a)', iostat=ios) line
  call assert_equal(trim(line), 'Info   : After the hold.', 'written after the hold')
  read (unit_log, '(a)', iostat=ios) line
  call assert_true(ios /= 0, 'held line not written')

  ! The wrapping alone, and the text helpers.
  call assert_equal(log_wrap('Error  : ', '', 'short'), 'Error  : short' // new_line('a'), &
    & 'log_wrap, one line')
  call assert_equal(log_real_text(1.0_rk, 4), '1.0000E+00', 'log_real_text')
  call assert_equal(log_real_text(-2.91e-4_rk, 4), '-2.9100E-04', 'log_real_text, negative')
  call assert_equal(log_real_text(1.200003e2_rk, 6), '1.200003E+02', 'log_real_text, 6 digits')
  call assert_equal(log_int_text(613_ik), '613', 'log_int_text')
  call assert_equal(log_count_text(1_ik, 'increment'), '1 increment', 'log_count_text, one')
  call assert_equal(log_count_text(4_ik, 'increment'), '4 increments', 'log_count_text, four')

  ! Back to standard output, then report.
  close (unit_log)
  call log_set_unit(output_unit)
  call test_summary()

end program test_log
