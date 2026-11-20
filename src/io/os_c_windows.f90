! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The file-system calls on Windows (64-bit; its C runtime and kernel32, e.g. with MSYS2's
!> gfortran).
!>
!> Source : FEPX. Contains no numerical method. The build compiles this file on Windows instead
!>          of os_c_posix.f90 (the same module), for dir_mod and restart_file_mod: Windows' C
!>          runtime has _mkdir(path), without the permissions of POSIX's mkdir, and its rename
!>          fails when the target exists, so a file is replaced with kernel32's MoveFileExA
!>          (MOVEFILE_REPLACE_EXISTING; Microsoft's documentation of MoveFileExA).
!>
!> Public entities:
!>   c_make_dir      create one directory (its parent must exist)
!>   c_replace_file  rename a file, replacing an existing target
module os_c_mod

  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char, c_long

  implicit none

  private

  public :: c_make_dir, c_replace_file

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

  !> MOVEFILE_REPLACE_EXISTING (winbase.h): MoveFileExA replaces an existing target.
  integer(c_long), parameter :: movefile_replace_existing = 1_c_long

  !> Interface to kernel32's MoveFileExA (BOOL MoveFileExA(LPCSTR existing, LPCSTR new,
  !> DWORD flags); 64-bit Windows has one calling convention, that of C).
  !>
  !> Arguments:
  !>   old, new  (in) character(kind=c_char)(*)  file names, null-terminated
  !>   flags     (in) integer(c_long)            MOVEFILE_* flags (DWORD)
  !> Result: non-zero on success
  interface
    function c_move_file_ex(old, new, flags) bind(C, name='MoveFileExA')
      import :: c_int, c_char, c_long
      character(kind=c_char), intent(in) :: old(*), new(*)
      integer(c_long), value :: flags
      integer(c_int) :: c_move_file_ex
    end function c_move_file_ex
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

  !> Rename a file, replacing an existing target (MoveFileExA with MOVEFILE_REPLACE_EXISTING: a
  !> rename on the same volume).
  !>
  !> Arguments:
  !>   old_name, new_name  (in)  character  file names
  !>   status  (out, result) integer(c_int)  0 on success, -1 on failure (as rename)
  function c_replace_file(old_name, new_name) result(status)

    character(len=*), intent(in) :: old_name, new_name
    integer(c_int) :: status

    if (c_move_file_ex(old_name // c_null_char, new_name // c_null_char, &
      & movefile_replace_existing) /= 0_c_int) then
      status = 0_c_int
    else
      status = -1_c_int
    end if

  end function c_replace_file

end module os_c_mod
