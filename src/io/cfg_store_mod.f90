! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Configuration file reading, stages 1 and 2: the file is read into records (keyword, values,
!> phase block, component, line), and every record is validated against its keyword's entry in
!> cfg_table_mod. The validated store is then turned into the typed configuration by config_mod
!> (stage 3, cross-keyword checks).
!>
!> Source : FEPX (ROADMAP step 4.3, D42; dev/cfg_keywords.md, io_format.md §2). Contains no
!>          numerical method.
!>
!> Syntax: one keyword per line, followed by its values; tokens are separated by blanks; "#"
!> starts a comment (so "## Section" headers are comments); blank lines are ignored. A per-phase
!> keyword belongs to the block of the last "phase i" line. Keywords are lower case; word values
!> are compared without regard to case (v2 files write e.g. "crystal_type FCC").
!>
!> Errors (fatal, "<file>:<line>: <message>"): unknown keyword; planned keyword (names its stage);
!> a per-phase keyword outside a phase block; a keyword given
!> twice (same phase and component), except repeatable lines; a wrong number of values (a list,
!> cfg_list, needs at least one; its count is checked by config_mod); a value of
!> the wrong type, out of range, or not among the allowed words; "-" where a value is required; a
!> per-step row whose length differs from number_of_steps; a component other than 11, 22, 33
!> (shear components: Stage 8).
!>
!> Public types:
!>   cfg_record_t   one validated line (components at the type)
!>   cfg_store_t    all records of a file (components at the type)
!>
!> Public routines:
!>   cfg_store_read  read and validate a file
!>   cfg_find        index of the record of a keyword (phase, component), 0 if absent
!>   cfg_error       fatal error at a record's line
!>   cfg_get_real, cfg_get_int, cfg_get_word  value i of a record, converted
!>   cfg_is_dash     whether value i of a record is "-"
module cfg_store_mod

  use kinds_mod, only: rk, ik
  use error_mod, only: fatal_error
  use text_file_mod, only: text_file_t, text_open, text_close, text_next_line, count_tokens, &
    & get_token
  use cfg_table_mod, only: cfg_keyword_t, cfg_keywords, cfg_find_keyword, cfg_active, &
    & cfg_planned, cfg_phase, cfg_single, cfg_per_step, cfg_comp_per_step, &
    & cfg_repeat, cfg_flag, cfg_list, cfg_int, cfg_real, cfg_word

  implicit none

  private

  public :: cfg_record_t, cfg_store_t, cfg_store_read, cfg_find, cfg_error
  public :: cfg_get_real, cfg_get_int, cfg_get_word, cfg_is_dash

  ! Length of a value token.
  integer, parameter :: token_len = 64

  !> One validated line.
  !>
  !> Components:
  !>   key     integer                 index of the keyword in cfg_keywords
  !>   phase   integer                 phase block (0: global keyword)
  !>   comp    integer                 component (11, 22, 33; shear 23, 13, 12, 32, 31, 21) of a
  !>                                   component row, else 0
  !>   line    integer(ik)             line number in the file
  !>   values  character(64)(:)        the values (without the component), as written
  type :: cfg_record_t
    integer :: key = 0
    integer :: phase = 0
    integer :: comp = 0
    integer(ik) :: line = 0_ik
    character(len=token_len), allocatable :: values(:)
  end type cfg_record_t

  !> All records of a file.
  !>
  !> Components:
  !>   file_name    character(:), allocatable  file name, for messages
  !>   num_records  integer                    records stored
  !>   records      type(cfg_record_t)(:)      the records (first num_records used)
  !>   num_steps    integer                    number_of_steps (0 if not given)
  type :: cfg_store_t
    character(len=:), allocatable :: file_name
    integer :: num_records = 0
    type(cfg_record_t), allocatable :: records(:)
    integer :: num_steps = 0
  end type cfg_store_t

contains

  !> Read a configuration file and validate every line against the keyword table.
  !>
  !> Arguments:
  !>   store      (out) type(cfg_store_t)  validated records
  !>   file_name  (in)  character          configuration file
  !> Locals:
  !>   file       type(text_file_t)        the open file
  !>   line       character(:), allocatable  current line, comment removed
  !>   eof        logical                  end of file
  !>   phase      integer                  current phase block (0 before the first)
  !>   rec        type(cfg_record_t)       the record of the current line
  !>   name       character(:), allocatable  keyword
  !>   ntok, i    integer                  token count, index
  !>   ios        integer                  read status
  subroutine cfg_store_read(store, file_name)

    type(cfg_store_t), intent(out) :: store
    character(len=*), intent(in) :: file_name

    type(text_file_t) :: file
    character(len=:), allocatable :: line, name
    logical :: eof
    integer :: phase, ntok, i, ios
    type(cfg_record_t) :: rec

    store%file_name = file_name
    allocate (store%records(16))
    phase = 0
    call text_open(file, file_name)

    do
      ! Next line without its comment; blank lines skipped.
      call text_next_line(file, line, eof)
      if (eof) exit
      i = index(line, '#')
      if (i > 0) line = line(1:i - 1)
      ntok = count_tokens(line)
      if (ntok == 0) cycle

      ! Keyword and its table entry (unknown and planned keywords are errors).
      name = get_token(line, 1)
      rec%key = lookup(file_name, file%line_num, name)
      rec%line = file%line_num
      rec%comp = 0

      ! Values (a component row's first value is its component).
      associate (k => rec%key)
        if (cfg_keywords(k)%form == cfg_comp_per_step) then
          if (ntok < 2) call line_error(file_name, rec%line, trim(name) // ': expected a ' // &
            & 'component (11, 22 or 33) and one value per step')
          rec%comp = component(file_name, rec%line, name, get_token(line, 2))
          rec%values = [(token(file_name, rec%line, get_token(line, i)), i = 3, ntok)]
        else
          rec%values = [(token(file_name, rec%line, get_token(line, i)), i = 2, ntok)]
        end if
      end associate

      ! Phase block: "phase i" opens block i; per-phase keywords need one.
      if (name == 'phase') then
        call check_values(store, rec, name)
        read (rec%values(1), *, iostat=ios) phase
      end if
      if (cfg_keywords(rec%key)%scope == cfg_phase) then
        if (phase == 0) call line_error(file_name, rec%line, trim(name) // ' must follow a ' // &
          & '"phase" line')
        rec%phase = phase
      else
        rec%phase = 0
      end if

      ! Values against the entry; a non-repeatable keyword given twice is an error.
      call check_values(store, rec, name)
      if (cfg_keywords(rec%key)%form /= cfg_repeat .and. name /= 'phase') then
        i = cfg_find(store, name, rec%phase, rec%comp)
        if (i > 0) then
          if (rec%comp > 0) name = name // ' ' // int_text(int(rec%comp, ik))
          call line_error(file_name, rec%line, name // ' given twice (lines ' // &
            & int_text(store%records(i)%line) // ' and ' // int_text(rec%line) // ')')
        end if
      end if
      if (name == 'phase') then
        do i = 1, store%num_records
          if (cfg_keywords(store%records(i)%key)%name == 'phase' .and. &
            & store%records(i)%values(1) == rec%values(1)) then
            call line_error(file_name, rec%line, 'phase ' // trim(rec%values(1)) // &
              & ' defined twice (lines ' // int_text(store%records(i)%line) // ' and ' // &
              & int_text(rec%line) // ')')
          end if
        end do
      end if
      if (name == 'number_of_steps') read (rec%values(1), *, iostat=ios) store%num_steps

      ! Store the record.
      call append(store, rec)
    end do
    call text_close(file)

    ! Per-step rows: exactly number_of_steps values.
    do i = 1, store%num_records
      associate (r => store%records(i))
        if (cfg_keywords(r%key)%form == cfg_per_step .or. &
          & cfg_keywords(r%key)%form == cfg_comp_per_step) then
          if (store%num_steps == 0) call line_error(file_name, r%line, &
            & trim(cfg_keywords(r%key)%name) // ': per-step row without number_of_steps')
          if (size(r%values) /= store%num_steps) call line_error(file_name, r%line, &
            & trim(cfg_keywords(r%key)%name) // ': ' // int_text(int(size(r%values), ik)) // &
            & ' values for ' // int_text(int(store%num_steps, ik)) // ' steps')
        end if
      end associate
    end do

  end subroutine cfg_store_read

  !> Table index of a keyword; unknown and planned keywords are errors.
  !>
  !> Arguments:
  !>   file_name  (in)  character    file, for messages
  !>   line       (in)  integer(ik)  line, for messages
  !>   name       (in)  character    keyword
  !>   key        (out, result) integer  table index
  function lookup(file_name, line, name) result(key)

    character(len=*), intent(in) :: file_name
    integer(ik), intent(in) :: line
    character(len=*), intent(in) :: name
    integer :: key

    ! Table entry.
    key = cfg_find_keyword(name)
    if (key == 0) call line_error(file_name, line, 'unknown keyword "' // name // '"')
    select case (cfg_keywords(key)%status)
    case (cfg_planned)
      ! (The stage of the note is a development record, not shown to the user, D104.)
      call line_error(file_name, line, name // ' is not available in this version of FEPX')
    end select

  end function lookup

  !> Component of a component row: 11, 22, 33, or a shear component 23, 13, 12, 32, 31, 21
  !> (periodic conditions only, step 8.3b: config_mod checks it).
  !>
  !> Arguments:
  !>   file_name, line, name  (in)  file, line and keyword, for messages
  !>   text   (in)  character   the component token
  !>   comp   (out, result) integer  the component ij as the integer 10 i + j
  function component(file_name, line, name, text) result(comp)

    character(len=*), intent(in) :: file_name
    integer(ik), intent(in) :: line
    character(len=*), intent(in) :: name, text
    integer :: comp

    select case (text)
    case ('11', '22', '33', '12', '13', '23', '21', '31', '32')
      read (text, *) comp
    case default
      call line_error(file_name, line, name // ': invalid component "' // text // &
        & '" (expected 11, 22, 33, or with periodic conditions 23, 13, 12, 32, 31, 21)')
    end select

  end function component

  !> Check a record's values against its table entry (count, type, range, allowed words, "-").
  !>
  !> Arguments:
  !>   store  (in)    type(cfg_store_t)   the store (file name)
  !>   rec    (inout) type(cfg_record_t)  the record (word values converted to lower case)
  !>   name   (in)    character           keyword, for messages
  !> Locals:
  !>   k      type(cfg_keyword_t)         the entry
  !>   i      integer                     value index
  subroutine check_values(store, rec, name)

    type(cfg_store_t), intent(in) :: store
    type(cfg_record_t), intent(inout) :: rec
    character(len=*), intent(in) :: name

    type(cfg_keyword_t) :: k
    integer :: i

    k = cfg_keywords(rec%key)

    ! Count (per-step rows are checked once number_of_steps is known).
    select case (k%form)
    case (cfg_single)
      if (size(rec%values) /= k%nvals) call line_error(store%file_name, rec%line, name // &
        & ': expected ' // int_text(int(k%nvals, ik)) // ' value(s), found ' // &
        & int_text(int(size(rec%values), ik)))
    case (cfg_flag)
      if (size(rec%values) /= 0) call line_error(store%file_name, rec%line, name // &
        & ' takes no value')
    case (cfg_list)
      if (size(rec%values) == 0) call line_error(store%file_name, rec%line, name // &
        & ': value(s) missing')
    case (cfg_repeat)
      if (size(rec%values) == 0) call line_error(store%file_name, rec%line, name // &
        & ': value(s) missing')
      if (k%nvals > 0 .and. size(rec%values) /= k%nvals) call line_error(store%file_name, &
        & rec%line, name // ': expected ' // int_text(int(k%nvals, ik)) // ' value(s)')
    end select

    ! Each value: "-" where allowed, then type, range and allowed words.
    if (k%form == cfg_repeat .and. k%nvals == 0) return
    do i = 1, size(rec%values)
      if (trim(rec%values(i)) == '-') then
        if (.not. k%step_dash) call line_error(store%file_name, rec%line, name // &
          & ': "-" is not allowed (every step needs a value)')
        cycle
      end if
      select case (k%kind)
      case (cfg_int)
        call check_number(store, rec, name, i, .true.)
      case (cfg_real)
        call check_number(store, rec, name, i, .false.)
      case (cfg_word)
        rec%values(i) = lower(rec%values(i))
        if (in_list(k%later_words, rec%values(i))) then
          call line_error(store%file_name, rec%line, name // ' ' // trim(rec%values(i)) // &
            & ' is not available in this version of FEPX')
        else if (.not. in_list(k%words, rec%values(i))) then
          call line_error(store%file_name, rec%line, name // ': invalid value "' // &
            & trim(rec%values(i)) // '" (expected one of: ' // trim(k%words) // ')')
        end if
      end select
    end do

  end subroutine check_values

  !> Check that value i of a record is a number of the right type within the entry's range.
  !>
  !> Arguments:
  !>   store, rec, name  (in)  store, record and keyword
  !>   i          (in)  integer  value index
  !>   is_int     (in)  logical  integer expected
  !> Locals:
  !>   ival, ios  integer   integer value, read status
  !>   x          real(rk)  value
  !>   k          type(cfg_keyword_t)  the entry
  !>   ok         logical   within the range
  subroutine check_number(store, rec, name, i, is_int)

    type(cfg_store_t), intent(in) :: store
    type(cfg_record_t), intent(in) :: rec
    character(len=*), intent(in) :: name
    integer, intent(in) :: i
    logical, intent(in) :: is_int

    integer :: ival, ios
    real(rk) :: x
    type(cfg_keyword_t) :: k
    logical :: ok

    k = cfg_keywords(rec%key)
    if (is_int) then
      read (rec%values(i), '(i64)', iostat=ios) ival
      x = real(ival, rk)
      if (ios /= 0 .or. verify(trim(rec%values(i)), '+-0123456789') /= 0) then
        call line_error(store%file_name, rec%line, name // ': "' // trim(rec%values(i)) // &
          & '" is not an integer')
      end if
    else
      read (rec%values(i), *, iostat=ios) x
      if (ios /= 0 .or. verify(trim(rec%values(i)), '+-.0123456789eEdD') /= 0) then
        call line_error(store%file_name, rec%line, name // ': "' // trim(rec%values(i)) // &
          & '" is not a number')
      end if
    end if
    ok = x <= k%hi .and. (x > k%lo .or. (.not. k%lo_open .and. .not. (x < k%lo)))
    if (.not. ok) call line_error(store%file_name, rec%line, name // ': value ' // &
      & trim(rec%values(i)) // ' out of range' // range_text(k))

  end subroutine check_number

  !> Index of the record of a keyword with a phase and a component, 0 if absent (the first,
  !> for repeatable lines).
  !>
  !> Arguments:
  !>   store  (in)  type(cfg_store_t)  the store
  !>   name   (in)  character          keyword
  !>   phase  (in)  integer            phase block (0: global)
  !>   comp   (in)  integer            component (0: none)
  !> Locals:
  !>   i      integer                  record index
  function cfg_find(store, name, phase, comp) result(index_rec)

    type(cfg_store_t), intent(in) :: store
    character(len=*), intent(in) :: name
    integer, intent(in) :: phase, comp
    integer :: index_rec

    integer :: i

    index_rec = 0
    do i = 1, store%num_records
      if (trim(cfg_keywords(store%records(i)%key)%name) == name .and. &
        & store%records(i)%phase == phase .and. store%records(i)%comp == comp) then
        index_rec = i
        return
      end if
    end do

  end function cfg_find

  !> Fatal error at the line of a record.
  !>
  !> Arguments:
  !>   store    (in)  type(cfg_store_t)  the store
  !>   i_rec    (in)  integer            record index
  !>   message  (in)  character          what is wrong
  subroutine cfg_error(store, i_rec, message)

    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_rec
    character(len=*), intent(in) :: message

    call line_error(store%file_name, store%records(i_rec)%line, message)

  end subroutine cfg_error

  !> Value i of record i_rec as a real (validated by cfg_store_read).
  !>
  !> Arguments:
  !>   store, i_rec, i  (in)  store, record and value index
  function cfg_get_real(store, i_rec, i) result(x)

    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_rec, i
    real(rk) :: x

    read (store%records(i_rec)%values(i), *) x

  end function cfg_get_real

  !> Value i of record i_rec as an integer (validated by cfg_store_read).
  !>
  !> Arguments:
  !>   store, i_rec, i  (in)  store, record and value index
  function cfg_get_int(store, i_rec, i) result(n)

    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_rec, i
    integer :: n

    read (store%records(i_rec)%values(i), *) n

  end function cfg_get_int

  !> Value i of record i_rec as a word (lower case, validated by cfg_store_read).
  !>
  !> Arguments:
  !>   store, i_rec, i  (in)  store, record and value index
  function cfg_get_word(store, i_rec, i) result(word)

    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_rec, i
    character(len=:), allocatable :: word

    word = trim(store%records(i_rec)%values(i))

  end function cfg_get_word

  !> Whether value i of record i_rec is "-".
  !>
  !> Arguments:
  !>   store, i_rec, i  (in)  store, record and value index
  function cfg_is_dash(store, i_rec, i) result(dash)

    type(cfg_store_t), intent(in) :: store
    integer, intent(in) :: i_rec, i
    logical :: dash

    dash = trim(store%records(i_rec)%values(i)) == '-'

  end function cfg_is_dash

  !> Append a record to the store (the record array grows by doubling).
  !>
  !> Arguments:
  !>   store  (inout) type(cfg_store_t)   the store
  !>   rec    (in)    type(cfg_record_t)  the record
  !> Locals:
  !>   grown  type(cfg_record_t)(:)       larger copy
  !>   i      integer                     index
  subroutine append(store, rec)

    type(cfg_store_t), intent(inout) :: store
    type(cfg_record_t), intent(in) :: rec

    type(cfg_record_t), allocatable :: grown(:)
    integer :: i

    ! Grow by doubling; records copied one by one (derived-type array constructors with
    ! allocatable components lost data with gfortran 13 -O3, ROADMAP section G).
    if (store%num_records == size(store%records)) then
      allocate (grown(2*size(store%records)))
      do i = 1, store%num_records
        grown(i) = store%records(i)
      end do
      call move_alloc(grown, store%records)
    end if
    store%num_records = store%num_records + 1
    store%records(store%num_records) = rec

  end subroutine append

  !> A value token, checked to fit the token length.
  !>
  !> Arguments:
  !>   file_name, line  (in)  file and line, for messages
  !>   text  (in)  character  token
  function token(file_name, line, text) result(tok)

    character(len=*), intent(in) :: file_name
    integer(ik), intent(in) :: line
    character(len=*), intent(in) :: text
    character(len=token_len) :: tok

    if (len(text) > token_len) call line_error(file_name, line, 'value "' // text // &
      & '" is longer than ' // int_text(int(token_len, ik)) // ' characters')
    tok = text

  end function token

  !> Text of an entry's range, for messages.
  !>
  !> Arguments:
  !>   k  (in)  type(cfg_keyword_t)  the entry
  function range_text(k) result(text)

    type(cfg_keyword_t), intent(in) :: k
    character(len=:), allocatable :: text

    text = ' (expected'
    if (k%lo > -huge(1.0_rk)) text = text // ' ' // trim(merge('> ', '>=', k%lo_open)) // ' ' &
      & // number_text(k%lo)
    if (k%lo > -huge(1.0_rk) .and. k%hi < huge(1.0_rk)) text = text // ' and'
    if (k%hi < huge(1.0_rk)) text = text // ' <= ' // number_text(k%hi)
    text = text // ')'

  end function range_text

  !> Short text of a range bound: an integer when the value is whole, else g0.
  !>
  !> Arguments:
  !>   x  (in)  real(rk)  the bound
  !> Locals:
  !>   buf  character(32)  formatting buffer
  function number_text(x) result(text)

    real(rk), intent(in) :: x
    character(len=:), allocatable :: text

    character(len=32) :: buf

    if (abs(x) < 1.0e15_rk .and. .not. abs(x - aint(x)) > 0.0_rk) then
      write (buf, '(i0)') nint(x, ik)
    else
      write (buf, '(g0)') x
    end if
    text = trim(buf)

  end function number_text

  !> Whether a word is in a blank-separated list.
  !>
  !> Arguments:
  !>   list  (in)  character  blank-separated words
  !>   word  (in)  character  the word
  pure function in_list(list, word) result(found)

    character(len=*), intent(in) :: list, word
    logical :: found

    found = index(' ' // trim(list) // ' ', ' ' // trim(word) // ' ') > 0

  end function in_list

  !> Fatal error "<file>:<line>: <message>".
  !>
  !> Arguments:
  !>   file_name, line, message  (in)  position and message
  subroutine line_error(file_name, line, message)

    character(len=*), intent(in) :: file_name
    integer(ik), intent(in) :: line
    character(len=*), intent(in) :: message

    call fatal_error('cfg', file_name // ':' // int_text(line) // ': ' // message)

  end subroutine line_error

  !> Lower-case copy of a token.
  !>
  !> Arguments:
  !>   text  (in)  character  token
  !> Locals:
  !>   i, c  integer  index, character code
  pure function lower(text) result(low)

    character(len=*), intent(in) :: text
    character(len=len(text)) :: low

    integer :: i, c

    low = text
    do i = 1, len(text)
      c = iachar(text(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) low(i:i) = achar(c + 32)
    end do

  end function lower

  !> Decimal text of an integer.
  !>
  !> Arguments:
  !>   value  (in)  integer(ik)  the integer
  !> Locals:
  !>   buf    character(24)      formatting buffer
  pure function int_text(value) result(text)

    integer(ik), intent(in) :: value
    character(len=:), allocatable :: text

    character(len=24) :: buf

    write (buf, '(i0)') value
    text = trim(buf)

  end function int_text

end module cfg_store_mod
