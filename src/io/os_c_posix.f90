! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The C library's file-system calls, on POSIX systems (Linux, macOS).
!>
!> Source : FEPX. Contains no numerical method. The build compiles this file, or os_c_windows.f90
!>          on Windows (the same module), for dir_mod and restart_file_mod.
!>
!> Public entities:
!>   c_make_dir      create one directory (its parent must exist)
!>   c_replace_file  rename a file, replacing an existing target
!>   c_is_dir        whether a path is an existing directory
!>   c_list_dir      the entries of a directory (ls -1A through popen, fgets and pclose; one per
!>                   line, sorted by ls)
module os_c_mod

  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char, c_ptr, c_associated

  implicit none

  private

  public :: c_make_dir, c_replace_file, c_is_dir, c_list_dir

  !> Interface to the C library's mkdir (POSIX: int mkdir(const char *path, mode_t mode)).
  !>
  !> Arguments:
  !>   name  (in) character(kind=c_char)(*)  path, null-terminated
  !>   mode  (in) integer(c_int)             permissions (before the umask)
  !> Result: 0 on success, -1 on failure (e.g. the directory exists)
  interface
    function c_mkdir(name, mode) bind(C, name='mkdir')
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int), value :: mode
      integer(c_int) :: c_mkdir
    end function c_mkdir
  end interface

  !> Interface to the C library's rename (stdio.h: int rename(const char *old, const char *new);
  !> POSIX: replaces an existing new, atomically).
  !>
  !> Arguments:
  !>   old, new  (in) character(kind=c_char)(*)  file names, null-terminated
  !> Result: 0 on success
  interface
    function c_rename(old, new) bind(C, name='rename')
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: old(*), new(*)
      integer(c_int) :: c_rename
    end function c_rename
  end interface

  !> The C library's popen, fgets and pclose (POSIX), for c_list_dir.
  interface
    function c_popen(command, mode) bind(C, name='popen') result(stream)
      import :: c_ptr, c_char
      character(kind=c_char), intent(in) :: command(*), mode(*)
      type(c_ptr) :: stream
    end function c_popen
    function c_fgets(buffer, size, stream) bind(C, name='fgets') result(line)
      import :: c_ptr, c_char, c_int
      character(kind=c_char), intent(inout) :: buffer(*)
      integer(c_int), value :: size
      type(c_ptr), value :: stream
      type(c_ptr) :: line
    end function c_fgets
    function c_pclose(stream) bind(C, name='pclose') result(status)
      import :: c_ptr, c_int
      type(c_ptr), value :: stream
      integer(c_int) :: status
    end function c_pclose
  end interface

contains

  !> Create one directory, readable and writable by all (before the umask), as mkdir does.
  !>
  !> Arguments:
  !>   path    (in)  character  the directory
  !>   status  (out, result) integer(c_int)  0 on success, -1 on failure (e.g. it exists)
  function c_make_dir(path) result(status)

    character(len=*), intent(in) :: path
    integer(c_int) :: status

    status = c_mkdir(path // c_null_char, int(o'777', c_int))

  end function c_make_dir

  !> Rename a file, replacing an existing target atomically (a reader sees the old file or the
  !> new one, never a part).
  !>
  !> Arguments:
  !>   old_name, new_name  (in)  character  file names
  !>   status  (out, result) integer(c_int)  0 on success
  function c_replace_file(old_name, new_name) result(status)

    character(len=*), intent(in) :: old_name, new_name
    integer(c_int) :: status

    status = c_rename(old_name // c_null_char, new_name // c_null_char)

  end function c_replace_file

  !> Whether a path is an existing directory (a directory's '.' entry exists; the Fortran
  !> runtimes of POSIX systems answer INQUIRE on it).
  !>
  !> Arguments:
  !>   path    (in)  character  the path
  !>   is_dir  (out, result) logical  .true. if it is an existing directory
  function c_is_dir(path) result(is_dir)

    character(len=*), intent(in) :: path
    logical :: is_dir

    inquire (file=trim(path) // '/.', exist=is_dir)

  end function c_is_dir

  !> The entries of an existing directory (ls -1A, through popen; one per line, sorted by ls),
  !> without "." and "..".
  !>
  !> Arguments:
  !>   path   (in)  character            the directory
  !>   names  (out) character(1024)(:), allocatable  its entries
  !>   ok     (out) logical              the directory was listed (ls succeeded)
  !> Locals:
  !>   stream  type(c_ptr)               the pipe
  !>   line    type(c_ptr)               fgets's result (null at the end)
  !>   buffer  character(c_char)(1025)   one line read
  !>   n       integer                   its length
  !>   status  integer(c_int)            pclose's status (that of ls)
  !>   entry   character(1024)           one entry
  subroutine c_list_dir(path, names, ok)

    character(len=*), intent(in) :: path
    character(len=1024), allocatable, intent(out) :: names(:)
    logical, intent(out) :: ok

    type(c_ptr) :: stream, line
    character(kind=c_char) :: buffer(1025)
    integer :: n
    integer(c_int) :: status
    character(len=1024) :: entry

    allocate (names(0))
    ok = .false.

    ! ls through a pipe, one entry per line.
    stream = c_popen('ls -1A "' // path // '" 2>/dev/null' // c_null_char, 'r' // c_null_char)
    if (.not. c_associated(stream)) return
    do
      buffer = c_null_char
      line = c_fgets(buffer, int(size(buffer), c_int), stream)
      if (.not. c_associated(line)) exit
      ! The line up to its null character, without the new-line character.
      n = 0
      do while (n < size(buffer) - 1)
        if (buffer(n + 1) == c_null_char) exit
        n = n + 1
      end do
      if (n > 0) then
        if (buffer(n) == achar(10)) n = n - 1
      end if
      if (n == 0) cycle
      entry = transfer(buffer(1:n), repeat(' ', n))
      if (entry == '.' .or. entry == '..') cycle
      names = [character(len=1024) :: names, entry]
    end do
    status = c_pclose(stream)
    ok = status == 0_c_int

  end subroutine c_list_dir

end module os_c_mod
