! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Line-oriented reading of text input files with strict checks and file:line error messages.
!>
!> Source : FEPX (ROADMAP step 2.2). Contains no numerical method. Used by the mesh reader
!>          (mesh_file_mod) and, later, the cfg parser (io_format.md §2: errors give file:line).
!>
!> Policy:
!>   - A file is read one whole line (record) at a time, of any length; trailing blanks are
!>     removed. Tokens are separated by blanks or tabs.
!>   - Values are parsed from a line only after its token count has been checked, so that a line
!>     with too few or too many values is an error (Fortran list-directed input alone would read
!>     on into the next line, or ignore extra values).
!>   - Every error is fatal (error_mod) with the message "<file>:<line>: <what>".
!>
!> Public entities:
!>   text_file_t        an open text file being read: unit, name, current line number
!>   text_open          open a file for reading
!>   text_close         close it
!>   text_next_line     read the next line (end of file reported through an argument)
!>   text_require_line  read the next line; end of file is a fatal error
!>   text_error         fatal error at the current line
!>   count_tokens       number of tokens of a line
!>   get_token          the i-th token of a line
!>   parse_ints         read exactly n integers from a line
!>   parse_reals        read exactly n reals from a line
!>   parse_int_token    read one integer token
module text_file_mod

  use kinds_mod, only: rk, ik
  use error_mod, only: fatal_error, check_status

  implicit none

  private

  public :: text_file_t, text_open, text_close, text_next_line, text_require_line, text_error
  public :: count_tokens, get_token, parse_ints, parse_reals, parse_int_token

  !> Length of the read buffer for one chunk of a line (lines may be longer: they are read in
  !> chunks).
  integer, parameter :: chunk_len = 1024

  !> An open text file being read.
  !>
  !> Components:
  !>   unit       integer                    Fortran unit (0 when closed)
  !>   name       character(:), allocatable  file name, for messages
  !>   line_num   integer(ik)                number of the line read last (0 before the first)
  type :: text_file_t
    integer :: unit = 0
    character(len=:), allocatable :: name
    integer(ik) :: line_num = 0_ik
  end type text_file_t

contains

  !> Open a file for reading; a missing or unreadable file is a fatal error, whose message is
  !> the run-time library's text (e.g. "Cannot open file 'simulation.cfg': No such file or
  !> directory"; step 9.7, dev/plan/terminal_output.md §3.5), or, without one,
  !> "cannot open file `<name>' (status <ios>)".
  !>
  !> Arguments:
  !>   file  (out)  type(text_file_t)  the opened file
  !>   name  (in)   character          file name
  !> Locals:
  !>   ios      integer         iostat of the open
  !>   message  character(256)  iomsg of the open
  !>   status_text  character(12)  ios as text
  subroutine text_open(file, name)

    type(text_file_t), intent(out) :: file
    character(len=*), intent(in) :: name

    integer :: ios
    character(len=256) :: message
    character(len=12) :: status_text

    ! Open the file for sequential formatted reading.
    message = ''
    open (newunit=file%unit, file=name, status='old', action='read', form='formatted', &
      & access='sequential', iostat=ios, iomsg=message)
    if (ios /= 0) then
      ! The run-time library's message is the whole error (it names the file and the reason).
      if (len_trim(message) > 0) call fatal_error('', trim(message))
      write (status_text, '(i0)') ios
      call fatal_error('', 'cannot open file `' // trim(name) // "' (status " // &
        & trim(status_text) // ')')
    end if

    ! Record the name; no line has been read yet.
    file%name = trim(name)
    file%line_num = 0_ik

  end subroutine text_open

  !> Close a file opened with text_open.
  !>
  !> Arguments:
  !>   file  (inout)  type(text_file_t)  the file; its unit is 0 afterwards
  subroutine text_close(file)

    type(text_file_t), intent(inout) :: file

    ! Close the unit, if open.
    if (file%unit /= 0) close (file%unit)
    file%unit = 0

  end subroutine text_close

  !> Read the next line of a file, of any length, without trailing blanks.
  !>
  !> Arguments:
  !>   file  (inout)  type(text_file_t)          the file; its line number is advanced
  !>   line  (out)    character(:), allocatable  the line (empty at end of file)
  !>   eof   (out)    logical                    .true. if the end of the file was reached
  !> Locals:
  !>   chunk    character(chunk_len)  one chunk of the line
  !>   num_read integer               number of characters read into chunk
  !>   ios      integer               iostat of the read
  !>   message  character(256)        iomsg of the read
  subroutine text_next_line(file, line, eof)

    type(text_file_t), intent(inout) :: file
    character(len=:), allocatable, intent(out) :: line
    logical, intent(out) :: eof

    character(len=chunk_len) :: chunk
    integer :: num_read
    integer :: ios
    character(len=256) :: message

    line = ''
    eof = .false.

    ! Read the record in chunks until its end (non-advancing input; iostat_eor ends a record).
    do
      message = ''
      read (file%unit, '(a)', advance='no', size=num_read, iostat=ios, iomsg=message) chunk
      if (is_iostat_end(ios)) then
        ! End of file: only valid at the start of a record.
        eof = .true.
        if (len(line) > 0) exit
        return
      end if
      line = line // chunk(1:num_read)
      if (is_iostat_eor(ios)) exit
      if (ios /= 0) then
        call fatal_error('text_next_line', file%name // ': read failed after line ' // &
          & int_text(file%line_num) // ': ' // trim(message))
      end if
    end do

    ! One more line read; drop trailing blanks (and a carriage return of a CRLF file).
    eof = .false.
    file%line_num = file%line_num + 1_ik
    line = trim(strip_cr(line))

  end subroutine text_next_line

  !> Read the next line; the end of the file is a fatal error naming what was expected.
  !>
  !> Arguments:
  !>   file      (inout)  type(text_file_t)          the file
  !>   line      (out)    character(:), allocatable  the line
  !>   expected  (in)     character                  what the line should contain (message)
  !> Locals:
  !>   eof  logical  end-of-file flag
  subroutine text_require_line(file, line, expected)

    type(text_file_t), intent(inout) :: file
    character(len=:), allocatable, intent(out) :: line
    character(len=*), intent(in) :: expected

    logical :: eof

    ! Read the line; end of file is an error.
    call text_next_line(file, line, eof)
    if (eof) then
      call fatal_error('text_require_line', file%name // ': unexpected end of file after line ' &
        & // int_text(file%line_num) // ' (expected ' // expected // ')')
    end if

  end subroutine text_require_line

  !> Fatal error at the current line of a file: "<file>:<line>: <message>".
  !>
  !> Arguments:
  !>   file     (in)  type(text_file_t)  the file
  !>   message  (in)  character          what is wrong
  subroutine text_error(file, message)

    type(text_file_t), intent(in) :: file
    character(len=*), intent(in) :: message

    ! Report with the position.
    call fatal_error('text_file', file%name // ':' // int_text(file%line_num) // ': ' // message)

  end subroutine text_error

  !> Number of tokens of a line (maximal runs of characters other than blank and tab).
  !>
  !> Arguments:
  !>   line  (in)  character  the line
  !>   num   (out, result) integer  number of tokens
  !> Locals:
  !>   i         integer  character index
  !>   in_token  logical  whether the previous character belongs to a token
  pure function count_tokens(line) result(num)

    character(len=*), intent(in) :: line
    integer :: num

    integer :: i
    logical :: in_token

    ! Count the starts of tokens.
    num = 0
    in_token = .false.
    do i = 1, len(line)
      if (is_separator(line(i:i))) then
        in_token = .false.
      else
        if (.not. in_token) num = num + 1
        in_token = .true.
      end if
    end do

  end function count_tokens

  !> The i-th token of a line (empty if the line has fewer tokens).
  !>
  !> Arguments:
  !>   line   (in)  character  the line
  !>   index  (in)  integer    token number (1-based)
  !>   token  (out, result) character(:), allocatable  the token
  !> Locals:
  !>   i        integer  character index
  !>   start    integer  first character of the current token
  !>   num      integer  tokens seen so far
  pure function get_token(line, index) result(token)

    character(len=*), intent(in) :: line
    integer, intent(in) :: index
    character(len=:), allocatable :: token

    integer :: i
    integer :: start
    integer :: num

    token = ''
    num = 0
    i = 1

    ! Walk through the tokens until the requested one.
    do while (i <= len(line))
      if (is_separator(line(i:i))) then
        i = i + 1
        cycle
      end if
      start = i
      do while (i <= len(line))
        if (is_separator(line(i:i))) exit
        i = i + 1
      end do
      num = num + 1
      if (num == index) then
        token = line(start:i - 1)
        return
      end if
    end do

  end function get_token

  !> Read exactly size(values) integers from a line; another token count, or a token that is not
  !> an integer, is a fatal error at the current line of the file.
  !>
  !> Arguments:
  !>   file    (in)   type(text_file_t)  the file (for messages)
  !>   line    (in)   character          the line
  !>   values  (out)  integer(ik)(:)     the values
  !>   what    (in)   character          what the line is, for messages (e.g. "node line")
  !> Locals:
  !>   num_tokens  integer         tokens on the line
  !>   ios         integer         iostat of the internal read
  subroutine parse_ints(file, line, values, what)

    type(text_file_t), intent(in) :: file
    character(len=*), intent(in) :: line
    integer(ik), intent(out) :: values(:)
    character(len=*), intent(in) :: what

    integer :: num_tokens
    integer :: ios

    ! The line must have exactly the expected number of values.
    num_tokens = count_tokens(line)
    if (num_tokens /= size(values)) then
      call text_error(file, what // ': expected ' // int_text(int(size(values), ik)) // &
        & ' values, found ' // int_text(int(num_tokens, ik)))
    end if

    ! Read them; a non-integer token is an error.
    values = 0_ik
    if (size(values) == 0) return
    read (line, *, iostat=ios) values
    if (ios /= 0 .or. .not. all_integer_tokens(line)) then
      call text_error(file, what // ': expected integers, found "' // line // '"')
    end if

  end subroutine parse_ints

  !> Read exactly size(values) reals from a line (as parse_ints).
  !>
  !> Arguments:
  !>   file    (in)   type(text_file_t)  the file (for messages)
  !>   line    (in)   character          the line
  !>   values  (out)  real(rk)(:)        the values
  !>   what    (in)   character          what the line is, for messages
  !> Locals:
  !>   num_tokens  integer  tokens on the line
  !>   ios         integer  iostat of the internal read
  subroutine parse_reals(file, line, values, what)

    type(text_file_t), intent(in) :: file
    character(len=*), intent(in) :: line
    real(rk), intent(out) :: values(:)
    character(len=*), intent(in) :: what

    integer :: num_tokens
    integer :: ios

    ! The line must have exactly the expected number of values.
    num_tokens = count_tokens(line)
    if (num_tokens /= size(values)) then
      call text_error(file, what // ': expected ' // int_text(int(size(values), ik)) // &
        & ' values, found ' // int_text(int(num_tokens, ik)))
    end if

    ! Read them; a token that is not a number is an error.
    values = 0.0_rk
    if (size(values) == 0) return
    read (line, *, iostat=ios) values
    if (ios /= 0) then
      call text_error(file, what // ': expected numbers, found "' // line // '"')
    end if

  end subroutine parse_reals

  !> Read one integer token (e.g. the token of a mixed line); not an integer is a fatal error.
  !>
  !> Arguments:
  !>   file   (in)  type(text_file_t)  the file (for messages)
  !>   token  (in)  character          the token
  !>   what   (in)  character          what the token is, for messages
  !>   value  (out, result) integer(ik)  the value
  !> Locals:
  !>   ios  integer  iostat of the internal read
  function parse_int_token(file, token, what) result(value)

    type(text_file_t), intent(in) :: file
    character(len=*), intent(in) :: token
    character(len=*), intent(in) :: what
    integer(ik) :: value

    integer :: ios

    ! Read the token; it must be an integer.
    value = 0_ik
    ios = 1
    if (len(token) > 0 .and. all_integer_tokens(token)) read (token, *, iostat=ios) value
    if (ios /= 0) call text_error(file, what // ': expected an integer, found "' // token // '"')

  end function parse_int_token

  !> Whether a character separates tokens (blank or tab).
  !>
  !> Arguments:
  !>   c    (in)  character(1)  the character
  !>   sep  (out, result) logical
  pure function is_separator(c) result(sep)

    character(len=1), intent(in) :: c
    logical :: sep

    sep = (c == ' ' .or. c == achar(9))

  end function is_separator

  !> Whether every token of a line is an optionally signed string of digits (list-directed input
  !> would accept e.g. "3.7" for an integer with some compilers, or stop at a slash).
  !>
  !> Arguments:
  !>   line  (in)  character  the line
  !>   ok    (out, result) logical
  !> Locals:
  !>   i           integer  character index
  !>   token_pos   integer  position within the current token (0 between tokens)
  !>   num_digits  integer  digits of the current token
  pure function all_integer_tokens(line) result(ok)

    character(len=*), intent(in) :: line
    logical :: ok

    integer :: i
    integer :: token_pos
    integer :: num_digits

    ok = .true.
    token_pos = 0
    num_digits = 0
    do i = 1, len(line) + 1
      ! End of a token (or of the line): it must have had a digit.
      if (i > len(line)) then
        if (token_pos > 0 .and. num_digits == 0) ok = .false.
        exit
      end if
      if (is_separator(line(i:i))) then
        if (token_pos > 0 .and. num_digits == 0) ok = .false.
        token_pos = 0
        num_digits = 0
        cycle
      end if
      ! Inside a token: a sign only first, otherwise digits.
      token_pos = token_pos + 1
      if (line(i:i) >= '0' .and. line(i:i) <= '9') then
        num_digits = num_digits + 1
      else if (.not. (token_pos == 1 .and. (line(i:i) == '-' .or. line(i:i) == '+'))) then
        ok = .false.
      end if
      if (.not. ok) exit
    end do

  end function all_integer_tokens

  !> A line without a trailing carriage return (files written on Windows).
  !>
  !> Arguments:
  !>   line      (in)  character  the line
  !>   stripped  (out, result) character(:), allocatable
  pure function strip_cr(line) result(stripped)

    character(len=*), intent(in) :: line
    character(len=:), allocatable :: stripped

    stripped = line
    if (len(stripped) > 0) then
      if (stripped(len(stripped):len(stripped)) == achar(13)) then
        stripped = stripped(1:len(stripped) - 1)
      end if
    end if

  end function strip_cr

  !> An integer as text, without blanks.
  !>
  !> Arguments:
  !>   value  (in)  integer(ik)  the integer
  !>   text   (out, result) character(:), allocatable
  !> Locals:
  !>   buffer  character(24)  formatted value
  pure function int_text(value) result(text)

    integer(ik), intent(in) :: value
    character(len=:), allocatable :: text

    character(len=24) :: buffer

    write (buffer, '(i0)') value
    text = trim(buffer)

  end function int_text

end module text_file_mod
