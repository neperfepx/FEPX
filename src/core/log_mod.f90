! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Log messages (information and warnings) for FEPX.
!>
!> Source : FEPX (ROADMAP steps 2.1, 9.7, D94; dev/plan/terminal_output.md §3). Contains no
!>          numerical method. Message format and rank selection follow v2-fixed
!>          parallel_mod.f90::par_message and v2's and Neper's message layout (prefixes
!>          "Info   : ", "Warning: ", "Error  : "; depth markers "- ", "> ", ". "; file markers
!>          "[i]", "[o]"; e.g. v2 fepx2.f90, read_input_cfg.f90; Neper's ut_print).
!>
!> Policy:
!>   - By default a message is written by rank 0 only (the other ranks return without writing).
!>     Callers on all ranks may therefore call the routines unconditionally.
!>   - Layout of a line: the prefix (9 columns), then two blanks per depth and the depth's marker
!>     (depth 1 "- ", depth 2 "> ", depth 3 ". "; depth 0: none), then the text. A file line
!>     (log_file) has "  [i] " or "  [o] " in place of the depth and marker.
!>   - No line is wider than log_width (72) columns: a longer message continues on further lines
!>     that start with the prefix and are aligned under the text, without the marker; a message
!>     of depth 0 continues two columns in (as depth 1 without its marker), so that a
!>     continuation line never looks like a new message. A message may ask for its continuation
!>     lines to be further in (log_info's hang: the input summary's parameter lines continue
!>     under their values). Lines are broken at blanks, but not
!>     where the continuation would start with a marker ("- ", "> ", ". ", "[i] ", "[o] "), so
!>     that the lines of a message can be joined again unambiguously (the tests); a word longer
!>     than the room left stays whole (the line is then wider).
!>   - With all_ranks = .true. every rank keeps its line, tagged "[Rank <r>] " after the marker,
!>     in a buffer; rank 0 writes the lines of every rank, in rank order, at the next collective
!>     point of the program (parallel_mod::par_log_flush, which uses log_buffer_take and
!>     log_lines_write). Use this for rank-specific information (e.g. a warning detected on one
!>     rank only).
!>   - Hold mode (log_hold_start, log_hold_take): rank 0 keeps its lines instead of writing
!>     them, for the caller to write them later at their place (the setup summary of fepx,
!>     which is written after the setup, dev/plan/terminal_output.md §3.2 revision 7).
!>     log_hold_release writes the kept lines at once (error_mod, before a fatal error, so that
!>     the lines leading to it are not lost).
!>   - Every warning is counted on the rank that writes or buffers it (log_warning_count); the
!>     run's total is the sum over the ranks (the closing summary).
!>   - Lines go to standard output (as in v2) and the unit is flushed after each message, so
!>     that the log of a batch job is current. Fatal errors are written by error_mod, which uses
!>     log_wrap and log_get_rank.
!>   - This module does not use MPI (MPI belongs in src/parallel/, conventions §6). The rank is
!>     set by parallel_mod::par_init through log_set_rank; serial programs keep the default
!>     rank 0 and therefore log.
!>
!> Public entities:
!>   log_info           write an information line (rank 0, or every rank through the buffer)
!>   log_warning        write a warning line (likewise), counted
!>   log_file           write a file line, "[i]" (file read) or "[o]" (file written)
!>   log_text           write a line as given, without prefix (rank 0; the banner's rules)
!>   log_wrap           the lines of a message, wrapped at log_width (used by error_mod)
!>   log_buffer_take    move this rank's buffered lines out of the buffer (par_log_flush)
!>   log_lines_write    write lines as given (rank 0; par_log_flush)
!>   log_warning_count  warnings counted on this rank
!>   log_hold_start     start keeping rank 0's lines (hold mode)
!>   log_hold_take      end the hold mode and hand over the kept lines
!>   log_hold_release   end the hold mode and write the kept lines (error_mod)
!>   log_set_rank       set the rank of this process (called by parallel_mod::par_init)
!>   log_get_rank       the rank of this process
!>   log_set_unit       set the output unit (default: standard output; used by the tests)
!>   log_real_text      a real as text in ES format without blanks, e.g. "1.0000E+00"
!>   log_int_text       an integer as text, e.g. "613"
!>   log_count_text     a count and its noun, plural unless 1, e.g. "4 increments"
!>   log_width          integer  largest width of a line (72, Neper's)
!>   log_info_prefix, log_warning_prefix, log_error_prefix  character(9)  the prefixes
!>
!> Module variables:
!>   log_rank      integer(ik)                rank of this process (0-based); 0 until
!>                                            log_set_rank is called
!>   log_unit      integer                    Fortran unit the lines are written to
!>   num_warnings  integer                    warnings written or buffered on this rank
!>   buffer        character(:), allocatable  buffered lines of this rank, each ended by a
!>                                            new-line character (unallocated: none)
!>   holding       logical                    hold mode on
!>   held          character(:), allocatable  lines kept in hold mode (rank 0), each ended
!>                                            by a new-line character
module log_mod

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik

  implicit none

  private

  public :: log_info, log_warning, log_file, log_text, log_wrap, log_buffer_take, &
    & log_lines_write, log_warning_count, log_set_rank, log_get_rank, log_set_unit, &
    & log_real_text, log_int_text, log_count_text, log_hold_start, log_hold_take, &
    & log_hold_release
  public :: log_width, log_info_prefix, log_warning_prefix, log_error_prefix

  !> Largest width of a line [columns] (Neper's).
  integer, parameter :: log_width = 72

  !> Line prefixes (9 columns, v2's and Neper's).
  character(len=9), parameter :: log_info_prefix = 'Info   : '
  character(len=9), parameter :: log_warning_prefix = 'Warning: '
  character(len=9), parameter :: log_error_prefix = 'Error  : '

  integer(ik), save :: log_rank = 0_ik
  integer, save :: log_unit = output_unit
  integer, save :: num_warnings = 0
  character(len=:), allocatable, save :: buffer
  logical, save :: holding = .false.
  character(len=:), allocatable, save :: held

contains

  !> Write an information line, "Info   : <depth marker><message>".
  !>
  !> Arguments:
  !>   message    (in)            character  text of the line
  !>   depth      (in, optional)  integer    depth 0 to 3 (default 0; see the module policy)
  !>   all_ranks  (in, optional)  logical    .true.: every rank's line, tagged with its rank,
  !>                                         buffered for rank 0; .false. (default): rank 0 only
  !>   hang       (in, optional)  integer    columns the continuation lines are further in than
  !>                                         the text (default 0)
  subroutine log_info(message, depth, all_ranks, hang)

    character(len=*), intent(in) :: message
    integer, intent(in), optional :: depth
    logical, intent(in), optional :: all_ranks
    integer, intent(in), optional :: hang

    ! Write (or buffer) the message with the information prefix.
    call put_message(log_info_prefix, depth_lead(depth), message, all_ranks, hang)

  end subroutine log_info

  !> Write a warning line, "Warning: <depth marker><message>", and count it.
  !>
  !> Arguments:
  !>   message    (in)            character  text of the line
  !>   depth      (in, optional)  integer    depth 0 to 3 (default 0)
  !>   all_ranks  (in, optional)  logical    as in log_info
  !> Locals:
  !>   every_rank  logical  whether every rank keeps the line (all_ranks, default .false.)
  subroutine log_warning(message, depth, all_ranks)

    character(len=*), intent(in) :: message
    integer, intent(in), optional :: depth
    logical, intent(in), optional :: all_ranks

    logical :: every_rank

    ! Count the warning on the rank that writes it (rank 0) or buffers it (every rank).
    every_rank = .false.
    if (present(all_ranks)) every_rank = all_ranks
    if (every_rank .or. log_rank == 0_ik) num_warnings = num_warnings + 1

    ! Write (or buffer) the message with the warning prefix.
    call put_message(log_warning_prefix, depth_lead(depth), message, all_ranks)

  end subroutine log_warning

  !> Write a file line, "Info   :   [i] <message>" (a file read) or "[o]" (a file written),
  !> on rank 0.
  !>
  !> Arguments:
  !>   marker   (in)  character  'i' (input) or 'o' (output)
  !>   message  (in)  character  text of the line, e.g. "Parsing file `simulation.cfg'..."
  subroutine log_file(marker, message)

    character(len=1), intent(in) :: marker
    character(len=*), intent(in) :: message

    ! Write the message after the file marker, at the indentation of depth 1.
    call put_message(log_info_prefix, '  [' // marker // '] ', message)

  end subroutine log_file

  !> Write one line as given, without prefix or wrapping (rank 0): the banner's rules.
  !>
  !> Arguments:
  !>   line  (in)  character  the line
  subroutine log_text(line)

    character(len=*), intent(in) :: line

    ! Rank 0 only.
    if (log_rank /= 0_ik) return
    write (log_unit, '(a)') line
    flush (log_unit)

  end subroutine log_text

  !> The lines of a message: the prefix and the lead (indentation and marker) on the first line,
  !> the prefix and as many blanks as the lead (at least two) on the following ones, the text
  !> broken at blanks so that no line is wider than log_width (a word longer than the room
  !> stays whole).
  !>
  !> Arguments:
  !>   prefix   (in)  character  line prefix, e.g. "Info   : "
  !>   lead     (in)  character  indentation and marker of the first line, e.g. "    > "
  !>   message  (in)  character  the text (trailing blanks ignored)
  !>   hang     (in, optional)  integer  further indentation of the continuation lines
  !>                            (default 0)
  !> Result:
  !>   lines    character(:), allocatable  the lines, each ended by a new-line character
  !> Locals:
  !>   text     character(:), allocatable  the text not yet placed
  !>   head     character(:), allocatable  prefix and lead of the current line
  !>   room     integer                    columns left for the text on a line
  !>   cut      integer                    length of the text placed on the current line
  !>   i        integer                    position in the text
  function log_wrap(prefix, lead, message, hang) result(lines)

    character(len=*), intent(in) :: prefix
    character(len=*), intent(in) :: lead
    character(len=*), intent(in) :: message
    integer, intent(in), optional :: hang
    character(len=:), allocatable :: lines

    character(len=:), allocatable :: text, head
    integer :: room, cut, i

    lines = ''
    text = trim(message)
    head = prefix // lead
    room = max(1, log_width - len(prefix) - len(lead))

    ! Place the text line by line.
    do
      if (len(text) <= room) then
        ! The rest fits on this line.
        lines = lines // head // text // new_line('a')
        exit
      end if
      ! The last blank within the room (a break at column room + 1 keeps the line full) after
      ! which the text does not start with a marker.
      cut = 0
      do i = room + 1, 2, -1
        if (text(i:i) /= ' ') cycle
        if (starts_with_marker(adjustl(text(i:)))) cycle
        cut = i - 1
        exit
      end do
      ! No blank in the room: break after the first word (the line is wider than the limit).
      if (cut == 0) then
        cut = index(text, ' ') - 1
        if (cut < 0) cut = len(text)
      end if
      ! Write the part, then continue with the rest without its leading blanks.
      lines = lines // head // trim(text(1:cut)) // new_line('a')
      text = adjustl(text(cut + 1:))
      text = trim(text)
      if (len(text) == 0) exit
      ! The continuation lines: the prefix, then blanks as wide as the lead (at least two),
      ! and the hang.
      head = prefix // repeat(' ', max(2, len(lead)))
      if (present(hang)) head = head // repeat(' ', max(0, hang))
      room = max(1, log_width - len(head))
    end do

  end function log_wrap

  !> Whether a text starts with a depth or file marker: "- ", "> ", ". ", "[i] " or "[o] ".
  !>
  !> Arguments:
  !>   text  (in)  character  the text
  logical function starts_with_marker(text) result(marked)

    character(len=*), intent(in) :: text

    ! The five markers.
    marked = .false.
    if (len(text) >= 2) marked = text(1:2) == '- ' .or. text(1:2) == '> ' .or. &
      & text(1:2) == '. '
    if (len(text) >= 4) marked = marked .or. text(1:4) == '[i] ' .or. text(1:4) == '[o] '

  end function starts_with_marker

  !> Move this rank's buffered lines (all_ranks messages) out of the buffer, leaving it empty.
  !>
  !> Arguments:
  !>   text  (out)  character(:), allocatable  the lines, each ended by a new-line character
  !>                                           (empty when there are none)
  subroutine log_buffer_take(text)

    character(len=:), allocatable, intent(out) :: text

    ! Hand the buffer over and empty it.
    if (allocated(buffer)) then
      call move_alloc(buffer, text)
    else
      text = ''
    end if

  end subroutine log_buffer_take

  !> Write lines as given (rank 0), each ended by a new-line character in the text: the
  !> buffered lines of every rank (par_log_flush).
  !>
  !> Arguments:
  !>   text  (in)  character  the lines
  !> Locals:
  !>   start, eol  integer  start of the current line and position of its new-line character
  subroutine log_lines_write(text)

    character(len=*), intent(in) :: text

    integer :: start, eol

    ! Rank 0 only.
    if (log_rank /= 0_ik) return

    ! One record per line.
    start = 1
    do while (start <= len(text))
      eol = index(text(start:), new_line('a'))
      if (eol == 0) then
        write (log_unit, '(a)') text(start:)
        exit
      end if
      write (log_unit, '(a)') text(start:start + eol - 2)
      start = start + eol
    end do
    flush (log_unit)

  end subroutine log_lines_write

  !> Start the hold mode: rank 0's lines (log_info, log_warning, log_file) are kept, not
  !> written, until log_hold_take or log_hold_release.
  subroutine log_hold_start()

    ! Keep from now on, starting empty.
    holding = .true.
    held = ''

  end subroutine log_hold_start

  !> End the hold mode and hand over the lines kept (empty on the other ranks).
  !>
  !> Arguments:
  !>   text  (out)  character(:), allocatable  the kept lines, each ended by a new-line
  !>                                           character
  subroutine log_hold_take(text)

    character(len=:), allocatable, intent(out) :: text

    ! Hand over the kept lines and stop keeping.
    holding = .false.
    if (allocated(held)) then
      call move_alloc(held, text)
    else
      text = ''
    end if

  end subroutine log_hold_take

  !> End the hold mode and write the lines kept so far (rank 0), e.g. before a fatal error.
  !>
  !> Locals:
  !>   text  character(:), allocatable  the kept lines
  subroutine log_hold_release()

    character(len=:), allocatable :: text

    ! Nothing kept unless holding.
    if (.not. holding) return
    call log_hold_take(text)
    call log_lines_write(text)

  end subroutine log_hold_release

  !> Warnings written or buffered on this rank so far.
  !>
  !> Result:
  !>   count  integer  the number of warnings
  integer function log_warning_count() result(count)

    ! The module's counter.
    count = num_warnings

  end function log_warning_count

  !> Set the rank of this process (0-based). Called once by parallel_mod::par_init.
  !>
  !> Arguments:
  !>   rank  (in)  integer(ik)  rank of this process in the run (0-based)
  subroutine log_set_rank(rank)

    integer(ik), intent(in) :: rank

    ! Store the rank used to select the writing process.
    log_rank = rank

  end subroutine log_set_rank

  !> The rank of this process (0-based; 0 in serial programs).
  !>
  !> Result:
  !>   rank  integer(ik)  the rank set by log_set_rank
  integer(ik) function log_get_rank() result(rank)

    ! The stored rank.
    rank = log_rank

  end function log_get_rank

  !> Set the Fortran unit the log lines are written to (default: standard output).
  !>
  !> Arguments:
  !>   unit  (in)  integer  an open, formatted, sequential unit
  subroutine log_set_unit(unit)

    integer, intent(in) :: unit

    ! Store the unit.
    log_unit = unit

  end subroutine log_set_unit

  !> A real as text in ES format with the given digits after the decimal point, without blanks:
  !> "1.0000E+00" (4 digits), "-2.9100E-04" (the terminal output's numbers, step 9.7).
  !>
  !> Arguments:
  !>   value   (in)  real(rk)  the value
  !>   digits  (in)  integer   digits after the decimal point (1 to 15)
  !> Result:
  !>   text    character(:), allocatable  the value as text
  !> Locals:
  !>   buffer  character(32)  the formatted value
  !>   form    character(16)  the edit descriptor
  function log_real_text(value, digits) result(text)

    real(rk), intent(in) :: value
    integer, intent(in) :: digits
    character(len=:), allocatable :: text

    character(len=32) :: buffer
    character(len=16) :: form

    ! ES with a width that fits any exponent, then without the leading blanks.
    write (form, '(a,i0,a,i0,a)') '(es', max(1, min(15, digits)) + 9, '.', &
      & max(1, min(15, digits)), ')'
    write (buffer, form) value
    text = trim(adjustl(buffer))

  end function log_real_text

  !> An integer as text, without blanks.
  !>
  !> Arguments:
  !>   value  (in)  integer(ik)  the value
  !> Result:
  !>   text   character(:), allocatable  the value as text
  !> Locals:
  !>   buffer  character(24)  the formatted value
  function log_int_text(value) result(text)

    integer(ik), intent(in) :: value
    character(len=:), allocatable :: text

    character(len=24) :: buffer

    ! Minimal width.
    write (buffer, '(i0)') value
    text = trim(buffer)

  end function log_int_text

  !> A count and its noun, the noun in the plural (with "s") unless the count is 1:
  !> "1 increment", "4 increments".
  !>
  !> Arguments:
  !>   count  (in)  integer(ik)  the count
  !>   noun   (in)  character    the noun in the singular
  !> Result:
  !>   text   character(:), allocatable  the count and the noun
  function log_count_text(count, noun) result(text)

    integer(ik), intent(in) :: count
    character(len=*), intent(in) :: noun
    character(len=:), allocatable :: text

    ! The count, then the noun (plural unless 1).
    text = log_int_text(count) // ' ' // trim(noun)
    if (count /= 1_ik) text = text // 's'

  end function log_count_text

  !> The indentation and marker of a depth: "" (0), "  - " (1), "    > " (2), "      . " (3).
  !>
  !> Arguments:
  !>   depth  (in, optional)  integer  depth (default 0; values outside 0 to 3 are clamped)
  !> Result:
  !>   lead   character(:), allocatable  the indentation and the marker
  !> Locals:
  !>   d      integer  the clamped depth
  function depth_lead(depth) result(lead)

    integer, intent(in), optional :: depth
    character(len=:), allocatable :: lead

    integer :: d

    ! Resolve and clamp the depth.
    d = 0
    if (present(depth)) d = max(0, min(3, depth))

    ! Two blanks per depth, then the depth's marker.
    select case (d)
    case (1)
      lead = '  - '
    case (2)
      lead = '    > '
    case (3)
      lead = '      . '
    case default
      lead = ''
    end select

  end function depth_lead

  !> Write a message on rank 0, or buffer it on every rank tagged with the rank (all_ranks).
  !>
  !> Arguments:
  !>   prefix     (in)            character  "Info   : " or "Warning: "
  !>   lead       (in)            character  indentation and marker (depth_lead, or a file
  !>                                         marker)
  !>   message    (in)            character  text of the message
  !>   all_ranks  (in, optional)  logical    as in log_info
  !>   hang       (in, optional)  integer    as in log_info
  !> Locals:
  !>   every_rank  logical        whether every rank keeps the message (default .false.)
  !>   rank_text   character(12)  log_rank as text
  subroutine put_message(prefix, lead, message, all_ranks, hang)

    character(len=*), intent(in) :: prefix
    character(len=*), intent(in) :: lead
    character(len=*), intent(in) :: message
    logical, intent(in), optional :: all_ranks
    integer, intent(in), optional :: hang

    logical :: every_rank
    character(len=12) :: rank_text

    ! Resolve the optional argument.
    every_rank = .false.
    if (present(all_ranks)) every_rank = all_ranks

    if (every_rank) then
      ! Every rank keeps its lines, tagged with its rank (v2 par_message, allwrite = .true.),
      ! for rank 0 to write at the next par_log_flush.
      write (rank_text, '(i0)') log_rank
      if (.not. allocated(buffer)) buffer = ''
      buffer = buffer // log_wrap(prefix, lead, '[Rank ' // trim(rank_text) // '] ' // message, &
        & hang)
    else
      ! Rank 0 only: kept in hold mode, otherwise written.
      if (log_rank /= 0_ik) return
      if (holding) then
        held = held // log_wrap(prefix, lead, message, hang)
      else
        call log_lines_write(log_wrap(prefix, lead, message, hang))
      end if
    end if

  end subroutine put_message

end module log_mod
