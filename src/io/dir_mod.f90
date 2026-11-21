! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Creation of directories (Fortran has no directory statement).
!>
!> Source : FEPX. Contains no numerical method. Used by the simulation-directory writer
!>          (sim_writer_mod).
!>
!> Policy:
!>   - A directory is created with the C library's mkdir (POSIX; _mkdir on Windows: os_c_mod,
!>     os_c_posix.f90 or os_c_windows.f90), called through iso_c_binding,
!>     not with execute_command_line('mkdir -p ...'): that one forks a shell from the MPI
!>     process, which fails on clusters whose MPI uses InfiniBand (Open MPI's openib transport:
!>     "cannot create directory" at the start of multi-node runs).
!>   - As mkdir -p, the parent directories are created first, and an existing directory is not
!>     an error; failing to get the directory (e.g. no write permission) is a fatal error.
!>
!> Public entities:
!>   make_dir    create a directory and its parents
!>   dir_exists  whether a path is an existing directory (os_c_mod)
module dir_mod

  use, intrinsic :: iso_c_binding, only: c_int
  use error_mod, only: fatal_error
  use os_c_mod, only: c_make_dir, c_is_dir

  implicit none

  private

  public :: make_dir, dir_exists

contains

  !> Create a directory and its parents (as mkdir -p); a fatal error if the directory does not
  !> exist afterwards.
  !>
  !> Arguments:
  !>   path    (in) character  directory
  !>   caller  (in) character  name of the calling routine, for the error message
  !> Locals:
  !>   i       integer         position in path
  !>   ierr    integer(c_int)  result of mkdir (not used: an existing parent fails, and the
  !>                           directory is checked at the end)
  !>   exists  logical         whether the directory exists at the end
  subroutine make_dir(path, caller)

    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: caller

    integer :: i
    integer(c_int) :: ierr
    logical :: exists

    ! the parents, at each '/' (a leading '/' is the root), then the directory itself
    do i = 2, len_trim(path)
      if (path(i:i) == '/') ierr = c_make_dir(path(1:i - 1))
    end do
    ierr = c_make_dir(trim(path))

    ! the directory must exist now
    exists = dir_exists(path)
    if (.not. exists) call fatal_error(caller, 'cannot create directory ' // path)

  end subroutine make_dir

  !> Whether a path is an existing directory (os_c_mod: INQUIRE on its '.' entry on POSIX
  !> systems, GetFileAttributesA on Windows).
  !>
  !> Arguments:
  !>   path    (in)  character  the path
  !>   exists  (out, result) logical  .true. if it is an existing directory
  function dir_exists(path) result(exists)

    character(len=*), intent(in) :: path
    logical :: exists

    exists = c_is_dir(path)

  end function dir_exists

end module dir_mod
