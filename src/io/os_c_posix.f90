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
module os_c_mod

  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char

  implicit none

  private

  public :: c_make_dir, c_replace_file, c_is_dir

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

end module os_c_mod
