! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Error cases of restart_file_mod (step 8.4, dev/theory/addendum_restart.md §5): each case,
!> named by the first command-line argument, writes a restart file in the working directory and
!> reads it in a way that must end with a fatal error, whose message is checked by
!> tests/core/check_run.cmake (list in tests/CMakeLists.txt).
!>
!> Source : FEPX (ROADMAP step 8.4).
!>
!> Cases:
!>   not_restart   a file that is not a restart file
!>   version       a restart file of another format version
!>   ranks         written on 2 ranks, read by a run on 3
!>   rank          the file of rank 1 read by rank 0
!>   missing       a file that does not exist
!>   wrong_name    record a written, record b expected
!>   wrong_shape   a real(3) record read as real(4)
!>   wrong_type    an integer record read as a real
!>   end_of_file   a record read after the last one
!>   state_size    a state of 2 slip systems read into one of 3
!>
!> Locals:
!>   case_name  character(32)          name of the case
!>   name       character(*)           the file name
!>   file       type(restart_file_t)   the file
!>   unit       integer                a unit (hand-written files)
!>   x3, x4     real(rk)(3), (4)       record values
!>   n          integer                an integer record
!>   state, other  type(material_state_t)  states of 2 and 3 slip systems
program test_restart_file_errors

  use, intrinsic :: iso_fortran_env, only: int32, real64
  use kinds_mod, only: rk, ik
  use restart_file_mod, only: restart_file_t, restart_create, restart_open, restart_close, &
    & restart_put, restart_get, restart_put_state, restart_get_state
  use material_state_mod, only: material_state_t, material_state_alloc

  implicit none

  character(len=32) :: case_name
  character(len=*), parameter :: name = 'errors.rst'
  type(restart_file_t) :: file
  integer :: unit
  real(rk) :: x3(3), x4(4)
  integer :: n
  type(material_state_t) :: state, other

  call get_command_argument(1, case_name)
  x3 = 1.0_rk
  x4 = 1.0_rk
  n = 1

  select case (trim(case_name))
  case ('not_restart')
    ! A text file.
    open (newunit=unit, file=name, status='replace', action='write')
    write (unit, '(a)') 'not a restart file at all'
    close (unit)
    call restart_open(file, name, 1, 0)
  case ('version')
    ! The header by hand, with version 99.
    open (newunit=unit, file=name, access='stream', form='unformatted', status='replace', &
      & action='write')
    write (unit) 'FEPX restart    ', 99_int32, 1_int32, 1.0_real64, int(rk, int32), &
      & int(ik, int32), 1_int32, 0_int32, 'even    '
    close (unit)
    call restart_open(file, name, 1, 0)
  case ('ranks')
    call restart_create(file, name, 2, 0)
    call restart_close(file)
    call restart_open(file, name, 3, 0)
  case ('rank')
    call restart_create(file, name, 2, 1)
    call restart_close(file)
    call restart_open(file, name, 2, 0)
  case ('missing')
    call restart_open(file, 'no_such_file.rst', 1, 0)
  case ('wrong_name')
    call restart_create(file, name, 1, 0)
    call restart_put(file, 'a', x3)
    call restart_close(file)
    call restart_open(file, name, 1, 0)
    call restart_get(file, 'b', x3)
  case ('wrong_shape')
    call restart_create(file, name, 1, 0)
    call restart_put(file, 'a', x3)
    call restart_close(file)
    call restart_open(file, name, 1, 0)
    call restart_get(file, 'a', x4)
  case ('wrong_type')
    call restart_create(file, name, 1, 0)
    call restart_put(file, 'a', n)
    call restart_close(file)
    call restart_open(file, name, 1, 0)
    call restart_get(file, 'a', x3(1))
  case ('end_of_file')
    call restart_create(file, name, 1, 0)
    call restart_put(file, 'a', n)
    call restart_close(file)
    call restart_open(file, name, 1, 0)
    call restart_get(file, 'a', n)
    call restart_get(file, 'a', n)
  case ('state_size')
    call material_state_alloc(state, 2_ik, 1_ik, 2_ik)
    call material_state_alloc(other, 2_ik, 1_ik, 3_ik)
    call restart_create(file, name, 1, 0)
    call restart_put_state(file, state)
    call restart_close(file)
    call restart_open(file, name, 1, 0)
    call restart_get_state(file, other)
  case default
    write (*, '(a)') 'unknown case ' // trim(case_name)
    stop 2
  end select

  ! Not reached when the case fails as expected.
  write (*, '(a)') 'case ' // trim(case_name) // ': no error'
  stop 3

end program test_restart_file_errors
