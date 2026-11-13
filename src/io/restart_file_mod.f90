! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Restart files (step 8.4, dev/theory/addendum_restart.md §5): binary stream files of named,
!> self-describing records, one file per rank, written and read in the same order.
!>
!> Source : FEPX (ROADMAP step 8.4; addendum_restart.md §5). Contains no numerical method. The
!>          concept (the state saved at a step end, one file per rank) is v2's
!>          (src/libfepx/restart/restart.f90); the format is FEPX's.
!>
!> Format (Fortran access='stream', form='unformatted': raw bytes, no record markers):
!>   header   magic (character(16) 'FEPX restart'), format version (int32), byte-order markers
!>            (int32 1, real64 1.0), the kinds rk and ik (int32 each), the number of ranks, the
!>            rank
!>   records  each: name (character(32)), type code (int32: 1 real(rk), 2 integer (default
!>            kind, = ik), 4 logical (stored as int32 0 or 1), 5 character), rank (int32),
!>            shape (int64 per dimension; character: its length), then the values in array
!>            element order
!> The reader asks for the records in the order they were written; a record of another name,
!> type or shape is an error naming the file, the record and both values (no silently skipped
!> or misread data).
!>
!> Public entities:
!>   restart_format_version  the format version written (int32)
!>   restart_file_t          an open restart file (unit, name, writing or reading)
!>   restart_create          create a file and write its header
!>   restart_open            open a file and check its header
!>   restart_close           close it
!>   restart_put             write a record (generic: real(rk) ranks 0-4, integer (= ik) ranks
!>                           0-1, logical ranks 0-1, character)
!>   restart_get             read the next record into an existing variable of the same shape
!>                           (generic, as restart_put); character into an allocatable string
!>   restart_put_state, restart_get_state  every component of a material_state_t
!>   restart_put_step, restart_get_step    every component of a step_record_t
!>   restart_rename          rename a file, replacing the target (the C library's rename,
!>                           through iso_c_binding: an atomic replacement on POSIX systems)
module restart_file_mod

  use, intrinsic :: iso_fortran_env, only: int32, int64, real64
  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
  use kinds_mod, only: rk, ik
  use error_mod, only: fatal_error, internal_error
  use material_state_mod, only: material_state_t
  use step_schedule_mod, only: step_record_t

  implicit none

  private

  public :: restart_format_version, restart_file_t, restart_create, restart_open, restart_close
  public :: restart_put, restart_get, restart_put_state, restart_get_state
  public :: restart_put_step, restart_get_step, restart_rename

  !> Format version written; another version is an error at reading.
  integer(int32), parameter :: restart_format_version = 2_int32

  !> Magic string at the start of every file.
  character(len=16), parameter :: magic = 'FEPX restart'

  !> Record type codes.
  integer(int32), parameter :: type_real = 1_int32, type_int = 2_int32, type_logical = 4_int32, &
    & type_char = 5_int32

  !> An open restart file.
  !>
  !> Components:
  !>   unit     integer                    Fortran unit (0 when closed)
  !>   name     character(:), allocatable  file name, for messages
  !>   writing  logical                    opened for writing (else reading)
  type :: restart_file_t
    integer :: unit = 0
    character(len=:), allocatable :: name
    logical :: writing = .false.
  end type restart_file_t

  !> Write a record.
  interface restart_put
    module procedure put_real0, put_real1, put_real2, put_real3, put_real4, put_int0, put_int1, &
      & put_logical0, put_logical1, put_char
  end interface restart_put

  !> The C library's rename (stdio.h): 0 on success.
  interface
    function c_rename(old, new) bind(c, name='rename') result(status)
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: old(*), new(*)
      integer(c_int) :: status
    end function c_rename
  end interface

  !> Read the next record.
  interface restart_get
    module procedure get_real0, get_real1, get_real2, get_real3, get_real4, get_int0, get_int1, &
      & get_logical0, get_logical1, get_char
  end interface restart_get

contains

  !> Create a restart file (replacing one of the same name) and write its header.
  !>
  !> Arguments:
  !>   file          (out) type(restart_file_t)  the file
  !>   name          (in)  character             file name
  !>   num_ranks     (in)  integer               number of ranks of the run
  !>   rank          (in)  integer               this rank (0-based)
  !> Locals:
  !>   ios      integer         iostat
  !>   message  character(256)  iomsg
  subroutine restart_create(file, name, num_ranks, rank)

    type(restart_file_t), intent(out) :: file
    character(len=*), intent(in) :: name
    integer, intent(in) :: num_ranks, rank

    integer :: ios
    character(len=256) :: message

    ! Open for stream writing.
    message = ''
    open (newunit=file%unit, file=name, access='stream', form='unformatted', &
      & status='replace', action='write', iostat=ios, iomsg=message)
    if (ios /= 0) call fatal_error('restart_create', 'cannot create ' // name // ': ' // &
      & trim(message))
    file%name = name
    file%writing = .true.

    ! Header: magic, version, byte-order markers, kinds, ranks, rank.
    write (file%unit) magic, restart_format_version, 1_int32, 1.0_real64, &
      & int(rk, int32), int(ik, int32), int(num_ranks, int32), int(rank, int32)

  end subroutine restart_create

  !> Open a restart file for reading and check its header: magic, format version, byte order,
  !> kinds, number of ranks and rank must be those given (errors naming the
  !> file and both values).
  !>
  !> Arguments:
  !>   file          (out) type(restart_file_t)  the file
  !>   name          (in)  character             file name
  !>   num_ranks     (in)  integer               number of ranks of this run
  !>   rank          (in)  integer               this rank (0-based)
  !> Locals:
  !>   ios      integer         iostat
  !>   message  character(256)  iomsg
  !>   magic_r  character(16)   magic read
  !>   version, one_int, rk_r, ik_r, ranks_r, rank_r  integer(int32)  header values read
  !>   one_real real(real64)    byte-order marker read
  subroutine restart_open(file, name, num_ranks, rank)

    type(restart_file_t), intent(out) :: file
    character(len=*), intent(in) :: name
    integer, intent(in) :: num_ranks, rank

    integer :: ios
    character(len=256) :: message
    character(len=16) :: magic_r
    integer(int32) :: version, one_int, rk_r, ik_r, ranks_r, rank_r
    real(real64) :: one_real

    ! Open for stream reading.
    message = ''
    open (newunit=file%unit, file=name, access='stream', form='unformatted', status='old', &
      & action='read', iostat=ios, iomsg=message)
    if (ios /= 0) call fatal_error('restart_open', 'cannot open restart file ' // name // &
      & ': ' // trim(message))
    file%name = name
    file%writing = .false.

    ! Magic and version first (a file of another kind or version stops here).
    read (file%unit, iostat=ios) magic_r
    if (ios /= 0 .or. magic_r /= magic) call fatal_error('restart_open', name // &
      & ' is not a FEPX restart file')
    read (file%unit, iostat=ios) version, one_int, one_real, rk_r, ik_r, ranks_r, rank_r
    if (ios /= 0) call fatal_error('restart_open', name // ': incomplete header')
    if (version /= restart_format_version) call fatal_error('restart_open', name // &
      & ': restart format version ' // itext(int(version)) // ', this FEPX reads version ' // &
      & itext(int(restart_format_version)))
    if (one_int /= 1_int32 .or. .not. (abs(one_real - 1.0_real64) < tiny(1.0_real64))) &
      & call fatal_error('restart_open', name // ': written with another byte order')
    if (rk_r /= rk .or. ik_r /= ik) call fatal_error('restart_open', name // &
      & ': written with other real or integer kinds')

    ! The run's layout: same number of ranks and rank (addendum §7).
    if (ranks_r /= num_ranks) call fatal_error('restart_open', name // ': written on ' // &
      & itext(int(ranks_r)) // ' rank(s); a restart needs the same number of ranks (this ' // &
      & 'run: ' // itext(num_ranks) // ')')
    if (rank_r /= rank) call fatal_error('restart_open', name // ': the file of rank ' // &
      & itext(int(rank_r)) // ', read by rank ' // itext(rank))

  end subroutine restart_open

  !> Close a restart file.
  !>
  !> Arguments:
  !>   file  (inout) type(restart_file_t)  the file (unit 0 afterwards)
  subroutine restart_close(file)

    type(restart_file_t), intent(inout) :: file

    if (file%unit /= 0) close (file%unit)
    file%unit = 0

  end subroutine restart_close

  ! --- Record headers ------------------------------------------------------------------------

  !> Write a record header: name, type code, rank and shape.
  !>
  !> Arguments:
  !>   file   (in)  type(restart_file_t)  the file (writing)
  !>   name   (in)  character             record name (at most 32 characters)
  !>   code   (in)  integer(int32)        type code
  !>   shp    (in)  integer(int64)(:)     shape (size = rank)
  !> Locals:
  !>   name32 character(32)               the name, fixed length
  subroutine put_header(file, name, code, shp)

    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    integer(int32), intent(in) :: code
    integer(int64), intent(in) :: shp(:)

    character(len=32) :: name32

    if (.not. file%writing) call internal_error('restart_put', file%name // &
      & ': not open for writing')
    if (len_trim(name) > 32) call internal_error('restart_put', 'record name ' // name // &
      & ' longer than 32 characters')
    name32 = name
    write (file%unit) name32, code, int(size(shp), int32), shp

  end subroutine put_header

  !> Read the next record header and check it against the expected name, type and shape
  !> (character: the shape is not checked, its length is returned).
  !>
  !> Arguments:
  !>   file   (in)  type(restart_file_t)  the file (reading)
  !>   name   (in)  character             expected record name
  !>   code   (in)  integer(int32)        expected type code
  !>   shp    (in)  integer(int64)(:)     expected shape (size = rank; ignored for character)
  !>   length (out, optional) integer(int64)  the length of a character record
  !> Locals:
  !>   name32  character(32)              name read
  !>   code_r, rank_r  integer(int32)     type code and rank read
  !>   shp_r   integer(int64)(:)          shape read
  !>   ios     integer                    iostat
  subroutine get_header(file, name, code, shp, length)

    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    integer(int32), intent(in) :: code
    integer(int64), intent(in) :: shp(:)
    integer(int64), intent(out), optional :: length

    character(len=32) :: name32
    integer(int32) :: code_r, rank_r
    integer(int64), allocatable :: shp_r(:)
    integer :: ios

    ! Name, type and rank.
    read (file%unit, iostat=ios) name32, code_r, rank_r
    if (ios /= 0) call fatal_error('restart_get', file%name // ': record ' // name // &
      & ' expected, end of file found')
    if (name32 /= name) call fatal_error('restart_get', file%name // ': record ' // name // &
      & ' expected, ' // trim(name32) // ' found')
    if (code_r /= code) call fatal_error('restart_get', file%name // ': record ' // name // &
      & ' has another type')
    if (rank_r < 0 .or. rank_r > 7) call fatal_error('restart_get', file%name // ': record ' &
      & // name // ' has an invalid rank')

    ! Shape: equal to the expected one (character: its length).
    allocate (shp_r(rank_r))
    read (file%unit, iostat=ios) shp_r
    if (ios /= 0) call fatal_error('restart_get', file%name // ': record ' // name // &
      & ' incomplete')
    if (code == type_char) then
      if (present(length)) length = shp_r(1)
      return
    end if
    if (size(shp_r) /= size(shp)) call fatal_error('restart_get', file%name // ': record ' // &
      & name // ' has rank ' // itext(int(rank_r)) // ', expected ' // itext(size(shp)))
    if (any(shp_r /= shp)) call fatal_error('restart_get', file%name // ': record ' // name // &
      & ' has shape (' // shape_text(shp_r) // '), expected (' // shape_text(shp) // ')')

  end subroutine get_header

  !> Read the values of a record; a short file is an error naming the record.
  !>
  !> Arguments:
  !>   file  (in)  type(restart_file_t)  the file
  !>   name  (in)  character             record name (for the message)
  !>   ios   (in)  integer               iostat of the read
  subroutine check_read(file, name, ios)

    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    integer, intent(in) :: ios

    if (ios /= 0) call fatal_error('restart_get', file%name // ': record ' // name // &
      & ': values incomplete')

  end subroutine check_read

  ! --- restart_put ---------------------------------------------------------------------------
  ! Each: (in) file type(restart_file_t), name character, value(s) of the given type and rank.

  subroutine put_real0(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    real(rk), intent(in) :: x
    call put_header(file, name, type_real, [integer(int64) ::])
    write (file%unit) x
  end subroutine put_real0

  subroutine put_real1(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    real(rk), intent(in) :: x(:)
    call put_header(file, name, type_real, int(shape(x), int64))
    write (file%unit) x
  end subroutine put_real1

  subroutine put_real2(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    real(rk), intent(in) :: x(:, :)
    call put_header(file, name, type_real, int(shape(x), int64))
    write (file%unit) x
  end subroutine put_real2

  subroutine put_real3(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    real(rk), intent(in) :: x(:, :, :)
    call put_header(file, name, type_real, int(shape(x), int64))
    write (file%unit) x
  end subroutine put_real3

  subroutine put_real4(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    real(rk), intent(in) :: x(:, :, :, :)
    call put_header(file, name, type_real, int(shape(x), int64))
    write (file%unit) x
  end subroutine put_real4

  subroutine put_int0(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    integer, intent(in) :: x
    call put_header(file, name, type_int, [integer(int64) ::])
    write (file%unit) x
  end subroutine put_int0

  subroutine put_int1(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    integer, intent(in) :: x(:)
    call put_header(file, name, type_int, int(shape(x), int64))
    write (file%unit) x
  end subroutine put_int1



  subroutine put_logical0(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    logical, intent(in) :: x
    call put_header(file, name, type_logical, [integer(int64) ::])
    write (file%unit) merge(1_int32, 0_int32, x)
  end subroutine put_logical0

  subroutine put_logical1(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    logical, intent(in) :: x(:)
    call put_header(file, name, type_logical, int(shape(x), int64))
    write (file%unit) merge(1_int32, 0_int32, x)
  end subroutine put_logical1

  subroutine put_char(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: x
    call put_header(file, name, type_char, [int(len(x), int64)])
    write (file%unit) x
  end subroutine put_char

  ! --- restart_get ---------------------------------------------------------------------------
  ! Each: (in) file, name; (inout) x, of the record's type and shape (allocated by the caller).
  ! Locals: ios integer (iostat); for logical values v integer(int32) (the stored 0 or 1).

  subroutine get_real0(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    real(rk), intent(inout) :: x
    integer :: ios
    call get_header(file, name, type_real, [integer(int64) ::])
    read (file%unit, iostat=ios) x
    call check_read(file, name, ios)
  end subroutine get_real0

  subroutine get_real1(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    real(rk), intent(inout) :: x(:)
    integer :: ios
    call get_header(file, name, type_real, int(shape(x), int64))
    read (file%unit, iostat=ios) x
    call check_read(file, name, ios)
  end subroutine get_real1

  subroutine get_real2(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    real(rk), intent(inout) :: x(:, :)
    integer :: ios
    call get_header(file, name, type_real, int(shape(x), int64))
    read (file%unit, iostat=ios) x
    call check_read(file, name, ios)
  end subroutine get_real2

  subroutine get_real3(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    real(rk), intent(inout) :: x(:, :, :)
    integer :: ios
    call get_header(file, name, type_real, int(shape(x), int64))
    read (file%unit, iostat=ios) x
    call check_read(file, name, ios)
  end subroutine get_real3

  subroutine get_real4(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    real(rk), intent(inout) :: x(:, :, :, :)
    integer :: ios
    call get_header(file, name, type_real, int(shape(x), int64))
    read (file%unit, iostat=ios) x
    call check_read(file, name, ios)
  end subroutine get_real4

  subroutine get_int0(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    integer, intent(inout) :: x
    integer :: ios
    call get_header(file, name, type_int, [integer(int64) ::])
    read (file%unit, iostat=ios) x
    call check_read(file, name, ios)
  end subroutine get_int0

  subroutine get_int1(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    integer, intent(inout) :: x(:)
    integer :: ios
    call get_header(file, name, type_int, int(shape(x), int64))
    read (file%unit, iostat=ios) x
    call check_read(file, name, ios)
  end subroutine get_int1



  subroutine get_logical0(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    logical, intent(inout) :: x
    integer :: ios
    integer(int32) :: v
    call get_header(file, name, type_logical, [integer(int64) ::])
    read (file%unit, iostat=ios) v
    call check_read(file, name, ios)
    x = v /= 0_int32
  end subroutine get_logical0

  subroutine get_logical1(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    logical, intent(inout) :: x(:)
    integer :: ios
    integer(int32), allocatable :: v(:)
    call get_header(file, name, type_logical, int(shape(x), int64))
    allocate (v(size(x)))
    read (file%unit, iostat=ios) v
    call check_read(file, name, ios)
    x = v /= 0_int32
  end subroutine get_logical1

  !> Character record: x (inout, allocatable) gets the stored string. Locals: length
  !> integer(int64) its length; buf character(:) the read buffer.
  subroutine get_char(file, name, x)
    type(restart_file_t), intent(in) :: file
    character(len=*), intent(in) :: name
    character(len=:), allocatable, intent(inout) :: x
    integer :: ios
    integer(int64) :: length
    character(len=:), allocatable :: buf
    call get_header(file, name, type_char, [integer(int64) ::], length)
    allocate (character(len=length) :: buf)
    read (file%unit, iostat=ios) buf
    call check_read(file, name, ios)
    x = buf
  end subroutine get_char

  ! --- Material state ------------------------------------------------------------------------

  !> Write every component of a material state (addendum_restart.md §4.1: the whole state, so
  !> that no field read by a later increment can be missed). A component added to
  !> material_state_t must be added here and in restart_get_state (test_restart_file checks the
  !> round trip of every component).
  !>
  !> Arguments:
  !>   file   (in)  type(restart_file_t)    the file (writing)
  !>   state  (in)  type(material_state_t)  the state (allocated)
  subroutine restart_put_state(file, state)

    type(restart_file_t), intent(in) :: file
    type(material_state_t), intent(in) :: state

    ! Sizes, then the arrays in the order of the type.
    call restart_put(file, 'state%num_qpts', state%num_qpts)
    call restart_put(file, 'state%num_elts', state%num_elts)
    call restart_put(file, 'state%num_slip', state%num_slip)
    call restart_put(file, 'state%start_from_vertex', state%start_from_vertex)
    call restart_put(file, 'state%stress_s', state%stress_s)
    call restart_put(file, 'state%strain_el_s', state%strain_el_s)
    call restart_put(file, 'state%rot_c2s', state%rot_c2s)
    call restart_put(file, 'state%spin_lat_s', state%spin_lat_s)
    call restart_put(file, 'state%crss', state%crss)
    call restart_put(file, 'state%sliprate', state%sliprate)
    call restart_put(file, 'state%rss', state%rss)
    call restart_put(file, 'state%slip', state%slip)
    call restart_put(file, 'state%slip_accum', state%slip_accum)
    call restart_put(file, 'state%strain_pl_c', state%strain_pl_c)
    call restart_put(file, 'state%defrate_pl_s', state%defrate_pl_s)
    call restart_put(file, 'state%spin_pl_s', state%spin_pl_s)
    call restart_put(file, 'state%rotrate_s', state%rotrate_s)
    call restart_put(file, 'state%workrate', state%workrate)
    call restart_put(file, 'state%work', state%work)
    call restart_put(file, 'state%defrate_s', state%defrate_s)
    call restart_put(file, 'state%spin_s', state%spin_s)
    call restart_put(file, 'state%strain_s', state%strain_s)
    call restart_put(file, 'state%strain_pl_s', state%strain_pl_s)
    call restart_put(file, 'state%spin_hat_update_s', state%spin_hat_update_s)
    call restart_put(file, 'state%spin_slip_update_s', state%spin_slip_update_s)
    call restart_put(file, 'state%workrate_pl', state%workrate_pl)
    call restart_put(file, 'state%work_pl', state%work_pl)
    call restart_put(file, 'state%km_rho', state%km_rho)
    call restart_put(file, 'state%km_k2', state%km_k2)
    call restart_put(file, 'state%km_rhorate', state%km_rhorate)
    call restart_put(file, 'state%dmg_c', state%dmg_c)
    call restart_put(file, 'state%dmgrate_c', state%dmgrate_c)

  end subroutine restart_put_state

  !> Read every component of a material state written by restart_put_state into a state
  !> allocated with the same sizes (material_state_alloc); other sizes are an error.
  !>
  !> Arguments:
  !>   file   (in)    type(restart_file_t)    the file (reading)
  !>   state  (inout) type(material_state_t)  the state (allocated)
  !> Locals:
  !>   n  integer(ik)  a size read
  subroutine restart_get_state(file, state)

    type(restart_file_t), intent(in) :: file
    type(material_state_t), intent(inout) :: state

    integer(ik) :: n

    ! Sizes: those of this run's state.
    n = -1_ik
    call restart_get(file, 'state%num_qpts', n)
    if (n /= state%num_qpts) call size_error('num_qpts', n, state%num_qpts)
    call restart_get(file, 'state%num_elts', n)
    if (n /= state%num_elts) call size_error('num_elts', n, state%num_elts)
    call restart_get(file, 'state%num_slip', n)
    if (n /= state%num_slip) call size_error('num_slip', n, state%num_slip)

    ! The arrays in the order of restart_put_state.
    call restart_get(file, 'state%start_from_vertex', state%start_from_vertex)
    call restart_get(file, 'state%stress_s', state%stress_s)
    call restart_get(file, 'state%strain_el_s', state%strain_el_s)
    call restart_get(file, 'state%rot_c2s', state%rot_c2s)
    call restart_get(file, 'state%spin_lat_s', state%spin_lat_s)
    call restart_get(file, 'state%crss', state%crss)
    call restart_get(file, 'state%sliprate', state%sliprate)
    call restart_get(file, 'state%rss', state%rss)
    call restart_get(file, 'state%slip', state%slip)
    call restart_get(file, 'state%slip_accum', state%slip_accum)
    call restart_get(file, 'state%strain_pl_c', state%strain_pl_c)
    call restart_get(file, 'state%defrate_pl_s', state%defrate_pl_s)
    call restart_get(file, 'state%spin_pl_s', state%spin_pl_s)
    call restart_get(file, 'state%rotrate_s', state%rotrate_s)
    call restart_get(file, 'state%workrate', state%workrate)
    call restart_get(file, 'state%work', state%work)
    call restart_get(file, 'state%defrate_s', state%defrate_s)
    call restart_get(file, 'state%spin_s', state%spin_s)
    call restart_get(file, 'state%strain_s', state%strain_s)
    call restart_get(file, 'state%strain_pl_s', state%strain_pl_s)
    call restart_get(file, 'state%spin_hat_update_s', state%spin_hat_update_s)
    call restart_get(file, 'state%spin_slip_update_s', state%spin_slip_update_s)
    call restart_get(file, 'state%workrate_pl', state%workrate_pl)
    call restart_get(file, 'state%work_pl', state%work_pl)
    call restart_get(file, 'state%km_rho', state%km_rho)
    call restart_get(file, 'state%km_k2', state%km_k2)
    call restart_get(file, 'state%km_rhorate', state%km_rhorate)
    call restart_get(file, 'state%dmg_c', state%dmg_c)
    call restart_get(file, 'state%dmgrate_c', state%dmgrate_c)

  contains

    !> A state size of the file differs from this run's. Arguments: what (in) character the
    !> size; saved, now (in) integer(ik) the file's and this run's values.
    subroutine size_error(what, saved, now)
      character(len=*), intent(in) :: what
      integer(ik), intent(in) :: saved, now
      call fatal_error('restart_get_state', file%name // ': material state ' // what // ' ' // &
        & itext(int(saved)) // ' in the restart file, ' // itext(int(now)) // ' in this run')
    end subroutine size_error

  end subroutine restart_get_state

  ! --- Schedule step -------------------------------------------------------------------------

  !> Write every component of a schedule step (step_record_t), as records 'step%<component>'.
  !>
  !> Arguments:
  !>   file  (in)  type(restart_file_t)  the file (writing)
  !>   rec   (in)  type(step_record_t)   the step
  subroutine restart_put_step(file, rec)

    type(restart_file_t), intent(in) :: file
    type(step_record_t), intent(in) :: rec

    ! The components in the order of the type.
    call restart_put(file, 'step%target_kind', rec%target_kind)
    call restart_put(file, 'step%is_load', rec%is_load)
    call restart_put(file, 'step%target_value', rec%target_value)
    call restart_put(file, 'step%target_tol', rec%target_tol)
    call restart_put(file, 'step%load_sign', rec%load_sign)
    call restart_put(file, 'step%dt_min', rec%dt_min)
    call restart_put(file, 'step%strain_rate', rec%strain_rate)
    call restart_put(file, 'step%vel_factor', rec%vel_factor)
    call restart_put(file, 'step%time_start', rec%time_start)
    call restart_put(file, 'step%time_end', rec%time_end)
    call restart_put(file, 'step%strain_start', rec%strain_start)
    call restart_put(file, 'step%strain_end', rec%strain_end)
    call restart_put(file, 'step%num_incrs', rec%num_incrs)
    call restart_put(file, 'step%dt', rec%dt)
    call restart_put(file, 'step%dt_last', rec%dt_last)
    call restart_put(file, 'step%print_results', rec%print_results)
    call restart_put(file, 'step%bc_kind', rec%bc_kind)
    call restart_put(file, 'step%bc_target', rec%bc_target)
    call restart_put(file, 'step%bc_tol', rec%bc_tol)
    call restart_put(file, 'step%bc_mode', rec%bc_mode)
    call restart_put(file, 'step%bc_rate', rec%bc_rate)
    call restart_put(file, 'step%rate_kind', rec%rate_kind)
    call restart_put(file, 'step%is_rate', rec%is_rate)
    call restart_put(file, 'step%is_dwell', rec%is_dwell)
    call restart_put(file, 'step%is_open', rec%is_open)
    call restart_put(file, 'step%ctrl_rate', rec%ctrl_rate)
    call restart_put(file, 'step%per_rate', rec%per_rate)

  end subroutine restart_put_step

  !> Read a schedule step written by restart_put_step.
  !>
  !> Arguments:
  !>   file  (in)    type(restart_file_t)  the file (reading)
  !>   rec   (inout) type(step_record_t)   the step
  subroutine restart_get_step(file, rec)

    type(restart_file_t), intent(in) :: file
    type(step_record_t), intent(inout) :: rec

    ! The components in the order of restart_put_step.
    call restart_get(file, 'step%target_kind', rec%target_kind)
    call restart_get(file, 'step%is_load', rec%is_load)
    call restart_get(file, 'step%target_value', rec%target_value)
    call restart_get(file, 'step%target_tol', rec%target_tol)
    call restart_get(file, 'step%load_sign', rec%load_sign)
    call restart_get(file, 'step%dt_min', rec%dt_min)
    call restart_get(file, 'step%strain_rate', rec%strain_rate)
    call restart_get(file, 'step%vel_factor', rec%vel_factor)
    call restart_get(file, 'step%time_start', rec%time_start)
    call restart_get(file, 'step%time_end', rec%time_end)
    call restart_get(file, 'step%strain_start', rec%strain_start)
    call restart_get(file, 'step%strain_end', rec%strain_end)
    call restart_get(file, 'step%num_incrs', rec%num_incrs)
    call restart_get(file, 'step%dt', rec%dt)
    call restart_get(file, 'step%dt_last', rec%dt_last)
    call restart_get(file, 'step%print_results', rec%print_results)
    call restart_get(file, 'step%bc_kind', rec%bc_kind)
    call restart_get(file, 'step%bc_target', rec%bc_target)
    call restart_get(file, 'step%bc_tol', rec%bc_tol)
    call restart_get(file, 'step%bc_mode', rec%bc_mode)
    call restart_get(file, 'step%bc_rate', rec%bc_rate)
    call restart_get(file, 'step%rate_kind', rec%rate_kind)
    call restart_get(file, 'step%is_rate', rec%is_rate)
    call restart_get(file, 'step%is_dwell', rec%is_dwell)
    call restart_get(file, 'step%is_open', rec%is_open)
    call restart_get(file, 'step%ctrl_rate', rec%ctrl_rate)
    call restart_get(file, 'step%per_rate', rec%per_rate)

  end subroutine restart_get_step

  ! --- Files ---------------------------------------------------------------------------------

  !> Rename a file, replacing an existing target (POSIX rename: atomic, so that a reader sees
  !> the old file or the new one, never a part; addendum_restart.md §5).
  !>
  !> Arguments:
  !>   old_name, new_name  (in)  character  file names
  !> Locals:
  !>   status  integer(c_int)  the C call's result
  subroutine restart_rename(old_name, new_name)

    character(len=*), intent(in) :: old_name, new_name

    integer(c_int) :: status

    status = c_rename(old_name // c_null_char, new_name // c_null_char)
    if (status /= 0_c_int) call fatal_error('restart_rename', 'cannot rename ' // old_name // &
      & ' to ' // new_name)

  end subroutine restart_rename

  ! --- Text helpers --------------------------------------------------------------------------

  !> An integer as text. Arguments: i (in) integer; t (out, result) character(:). Locals: buf
  !> character(24).
  pure function itext(i) result(t)
    integer, intent(in) :: i
    character(len=:), allocatable :: t
    character(len=24) :: buf
    write (buf, '(i0)') i
    t = trim(buf)
  end function itext

  !> A shape as "n1 n2 ...". Arguments: shp (in) integer(int64)(:); t (out, result)
  !> character(:). Locals: buf character(24); j integer.
  pure function shape_text(shp) result(t)
    integer(int64), intent(in) :: shp(:)
    character(len=:), allocatable :: t
    character(len=24) :: buf
    integer :: j
    t = ''
    do j = 1, size(shp)
      write (buf, '(i0)') shp(j)
      t = t // trim(buf)
      if (j < size(shp)) t = t // ' '
    end do
  end function shape_text

end module restart_file_mod
