! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit test of the file-system calls (os_c_mod, os_c_posix.f90 or os_c_windows.f90) and of
!> dir_mod, which uses them.
!>
!> Source : FEPX (PR #382). Checks, in the working directory:
!>   - make_dir creates a directory and its parents; dir_exists (c_is_dir) is .true. for it, and
!>     .false. for a missing path and for a file;
!>   - c_list_dir lists a directory's entries (files and directories, without "." and ".."),
!>     sorted (lowercase names only: ls sorts by the locale on POSIX systems), and fails on a
!>     missing directory;
!>   - c_replace_file renames a file over an existing one (the target then has the new content,
!>     and the old name is gone), as the restart files are written.
!> The directories are left in place (no rmdir); every file is written again at each run, so that
!> a second run sees the same entries.
!>
!> Locals:
!>   base       character(*)              the directory of the listing checks
!>   repl       character(*)              the directory of the replacement checks
!>   names      character(1024)(:), allocatable  entries listed
!>   ok         logical                   c_list_dir's status
!>   exists     logical                   whether a file exists
!>   unit, ios  integer                   a file's unit and I/O status
!>   status     integer(c_int)            c_replace_file's result
!>   line       character(16)             a line read back
program test_os_c

  use, intrinsic :: iso_c_binding, only: c_int
  use os_c_mod, only: c_list_dir, c_replace_file, c_is_dir
  use dir_mod, only: make_dir, dir_exists
  use assert_mod, only: assert_equal, assert_true, test_summary

  implicit none

  character(len=*), parameter :: base = 'test_os_c_dir'
  character(len=*), parameter :: repl = 'test_os_c_replace'
  character(len=1024), allocatable :: names(:)
  logical :: ok
  logical :: exists
  integer :: unit, ios
  integer(c_int) :: status
  character(len=16) :: line

  ! --- Directories ---------------------------------------------------------------------------

  ! A directory and its parent; a missing path; a file.
  call make_dir(base // '/sub', 'test_os_c')
  call assert_true(dir_exists(base), 'dir_exists: parent created')
  call assert_true(dir_exists(base // '/sub'), 'dir_exists: directory created')
  call assert_true(.not. dir_exists(base // '/missing'), 'dir_exists: missing path')
  call write_file(base // '/a.txt', 'a')
  call write_file(base // '/b.txt', 'b')
  call assert_true(.not. c_is_dir(base // '/a.txt'), 'c_is_dir: a file is not a directory')

  ! --- Listing -------------------------------------------------------------------------------

  ! The two files and the directory, sorted, without "." and "..".
  call c_list_dir(base, names, ok)
  call assert_true(ok, 'c_list_dir: ok')
  call assert_equal(size(names), 3, 'c_list_dir: 3 entries')
  if (size(names) == 3) then
    call assert_equal(trim(names(1)), 'a.txt', 'c_list_dir: entry 1')
    call assert_equal(trim(names(2)), 'b.txt', 'c_list_dir: entry 2')
    call assert_equal(trim(names(3)), 'sub', 'c_list_dir: entry 3')
  end if

  ! An empty directory, and a missing one.
  call c_list_dir(base // '/sub', names, ok)
  call assert_true(ok, 'c_list_dir: empty directory ok')
  call assert_equal(size(names), 0, 'c_list_dir: empty directory, no entry')
  call c_list_dir(base // '/missing', names, ok)
  call assert_true(.not. ok, 'c_list_dir: missing directory')
  call assert_equal(size(names), 0, 'c_list_dir: missing directory, no entry')

  ! --- Replacement ---------------------------------------------------------------------------

  ! A new file renamed over an existing one.
  call make_dir(repl, 'test_os_c')
  call write_file(repl // '/target', 'old')
  call write_file(repl // '/target.new', 'new')
  status = c_replace_file(repl // '/target.new', repl // '/target')
  call assert_equal(int(status), 0, 'c_replace_file: status')
  inquire (file=repl // '/target.new', exist=exists)
  call assert_true(.not. exists, 'c_replace_file: old name gone')
  open (newunit=unit, file=repl // '/target', status='old', action='read', iostat=ios)
  call assert_equal(ios, 0, 'c_replace_file: target opens')
  if (ios == 0) then
    read (unit, '(a)', iostat=ios) line
    close (unit)
    call assert_equal(trim(line), 'new', 'c_replace_file: target has the new content')
  end if

  ! A missing source: an error status.
  status = c_replace_file(repl // '/missing', repl // '/target')
  call assert_true(status /= 0_c_int, 'c_replace_file: missing source fails')

  call test_summary()

contains

  !> Write a file of one line (replacing it if it exists).
  !>
  !> Arguments:
  !>   name  (in)  character  the file
  !>   text  (in)  character  its line
  !> Locals:
  !>   u     integer          its unit
  subroutine write_file(name, text)

    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: text

    integer :: u

    open (newunit=u, file=name, status='replace', action='write')
    write (u, '(a)') text
    close (u)

  end subroutine write_file

end program test_os_c
