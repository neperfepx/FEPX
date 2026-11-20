! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The C library's directory creation, on Windows (its C runtime, e.g. with MSYS2's gfortran).
!>
!> Source : FEPX. Contains no numerical method. The build compiles this file on Windows instead
!>          of dir_c_posix.f90 (the same module), for dir_mod: Windows' C runtime has
!>          _mkdir(path), without the permissions of POSIX's mkdir.
!>
!> Public entities:
!>   c_make_dir  create one directory (its parent must exist)
module dir_c_mod

  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char

  implicit none

  private

  public :: c_make_dir

  !> Interface to the C runtime's _mkdir (Windows: int _mkdir(const char *path)).
  !>
  !> Arguments:
  !>   name  (in) character(kind=c_char)(*)  path, null-terminated
  !> Result: 0 on success, -1 on failure (e.g. the directory exists)
  interface
    function c_mkdir(name) bind(C, name='_mkdir')
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int) :: c_mkdir
    end function c_mkdir
  end interface

contains

  !> Create one directory.
  !>
  !> Arguments:
  !>   path    (in)  character  the directory ('/' or '\' as separators)
  !>   status  (out, result) integer(c_int)  0 on success, -1 on failure (e.g. it exists)
  function c_make_dir(path) result(status)

    character(len=*), intent(in) :: path
    integer(c_int) :: status

    status = c_mkdir(path // c_null_char)

  end function c_make_dir

end module dir_c_mod
