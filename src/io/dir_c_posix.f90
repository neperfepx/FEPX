! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The C library's directory creation, on POSIX systems (Linux, macOS).
!>
!> Source : FEPX. Contains no numerical method. The build compiles this file, or dir_c_windows.f90
!>          on Windows (the same module), for dir_mod.
!>
!> Public entities:
!>   c_make_dir  create one directory (its parent must exist)
module dir_c_mod

  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char

  implicit none

  private

  public :: c_make_dir

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

end module dir_c_mod
