! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of text_file_mod: line reading (any length, CRLF, trailing blanks), line numbers,
!> end of file, token counting and extraction, strict integer and real parsing.
!>
!> Source : FEPX (ROADMAP step 2.2). Error paths (wrong counts, non-numbers, missing file,
!>          unexpected end of file) are covered by test_mesh_file_errors.
!>
!> Locals:
!>   name      character(*)              scratch file name (working directory)
!>   unit      integer                   unit used to write it
!>   file      text_file_t               the file being read
!>   line      character(:), allocatable  a line read back
!>   eof       logical                   end-of-file flag
!>   long      character(:), allocatable  a 3000-character line of tokens
!>   ints      integer(ik)(3)            parsed integers
!>   reals     real(rk)(2)               parsed reals
!>   i         integer                   token index
program test_text_file

  use kinds_mod, only: rk, ik
  use text_file_mod, only: text_file_t, text_open, text_close, text_next_line, &
    & text_require_line, count_tokens, get_token, parse_ints, parse_reals, parse_int_token
  use assert_mod, only: assert_equal, assert_close, assert_true, test_summary

  implicit none

  character(len=*), parameter :: name = 'test_text_file.txt'
  integer :: unit
  type(text_file_t) :: file
  character(len=:), allocatable :: line
  logical :: eof
  character(len=:), allocatable :: long
  integer(ik) :: ints(3)
  real(rk) :: reals(2)
  integer :: i

  ! A line longer than the read chunk: 750 tokens "1234" separated by blanks (3749 chars).
  long = '1234'
  do i = 2, 750
    long = long // ' 1234'
  end do

  ! Write the test file: long line, CRLF line, trailing blanks, tabs, blank line, numbers.
  open (newunit=unit, file=name, status='replace', action='write')
  write (unit, '(a)') long
  write (unit, '(a)') 'crlf line' // achar(13)
  write (unit, '(a)') 'trailing   '
  write (unit, '(a)') achar(9) // 'a' // achar(9) // achar(9) // 'bb   ccc' // achar(9)
  write (unit, '(a)') ''
  write (unit, '(a)') ' -12  +7 3 '
  write (unit, '(a)') '1.5e3 -2.25d0'
  close (unit)

  ! Read it back.
  call text_open(file, name)
  call assert_equal(int(file%line_num), 0, 'no line read yet')

  ! Long line, whole.
  call text_next_line(file, line, eof)
  call assert_true(.not. eof, 'line 1 not eof')
  call assert_equal(len(line), len(long), 'long line length')
  call assert_true(line == long, 'long line content')
  call assert_equal(count_tokens(line), 750, 'long line tokens')
  call assert_equal(get_token(line, 750), '1234', 'long line last token')

  ! Carriage return and trailing blanks removed.
  call text_next_line(file, line, eof)
  call assert_equal(line, 'crlf line', 'crlf stripped')
  call text_require_line(file, line, 'line 3')
  call assert_equal(line, 'trailing', 'trailing blanks removed')
  call assert_equal(int(file%line_num), 3, 'line number 3')

  ! Tabs separate tokens.
  call text_next_line(file, line, eof)
  call assert_equal(count_tokens(line), 3, 'tab tokens')
  call assert_equal(get_token(line, 1), 'a', 'token 1')
  call assert_equal(get_token(line, 2), 'bb', 'token 2')
  call assert_equal(get_token(line, 3), 'ccc', 'token 3')
  call assert_equal(get_token(line, 4), '', 'token 4 missing')

  ! Blank line.
  call text_next_line(file, line, eof)
  call assert_true(.not. eof, 'blank line not eof')
  call assert_equal(len(line), 0, 'blank line empty')
  call assert_equal(count_tokens(line), 0, 'blank line no tokens')

  ! Signed integers.
  call text_next_line(file, line, eof)
  call parse_ints(file, line, ints, 'integer line')
  call assert_equal(ints(1), -12_ik, 'int 1')
  call assert_equal(ints(2), 7_ik, 'int 2')
  call assert_equal(ints(3), 3_ik, 'int 3')
  call assert_equal(parse_int_token(file, get_token(line, 1), 'token'), -12_ik, 'int token')

  ! Reals in E and D exponent forms.
  call text_next_line(file, line, eof)
  call parse_reals(file, line, reals, 'real line')
  call assert_close(reals, [1500.0_rk, -2.25_rk], 'reals')

  ! End of file: eof set, empty line, line number unchanged.
  call text_next_line(file, line, eof)
  call assert_true(eof, 'eof')
  call assert_equal(len(line), 0, 'eof line empty')
  call assert_equal(int(file%line_num), 7, 'line number at eof')
  call text_close(file)
  call assert_equal(file%unit, 0, 'closed')

  ! Remove the scratch file.
  open (newunit=unit, file=name, status='old')
  close (unit, status='delete')

  call test_summary()

end program test_text_file
