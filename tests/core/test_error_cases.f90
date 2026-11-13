! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Abort handlers used by test_error_cases. They are module procedures, not internal procedures
!> of the program: passing an internal procedure as an argument needs an executable stack
!> (trampoline) with gfortran.
!>
!> Source : FEPX (ROADMAP step 1.1).
!>
!> Public routines:
!>   stopping_handler   reports its call and ends the run with status 7
!>   returning_handler  reports its call and returns
module test_error_handlers_mod

  use, intrinsic :: iso_fortran_env, only: error_unit

  implicit none

  private

  public :: stopping_handler, returning_handler

contains

  !> Abort handler that reports its call (on standard error, the unit of fatal_error's message,
  !> so that the order of the two lines is preserved) and ends the run with status 7.
  !>
  !> Arguments:
  !>   exit_status  (in)  integer  exit status passed by fatal_error
  subroutine stopping_handler(exit_status)

    integer, intent(in) :: exit_status

    ! Report the call and the status received, then end the run with a distinct status.
    write (error_unit, '(a,i0)') 'stopping handler called with status ', exit_status
    flush (error_unit)
    error stop 7

  end subroutine stopping_handler

  !> Abort handler that (wrongly) returns after reporting its call on standard error.
  !>
  !> Arguments:
  !>   exit_status  (in)  integer  exit status passed by fatal_error
  subroutine returning_handler(exit_status)

    integer, intent(in) :: exit_status

    ! Report the call and the status received, then return.
    write (error_unit, '(a,i0)') 'returning handler called with status ', exit_status
    flush (error_unit)

  end subroutine returning_handler

end module test_error_handlers_mod

!> Test driver for error_mod: runs one case, named by the first command-line argument.
!>
!> Source : FEPX (ROADMAP step 1.1).
!>
!> Each case ends the program in a known way. The exit status, standard output and standard
!> error are checked by tests/core/check_run.cmake, as set up in tests/CMakeLists.txt; every
!> error message must appear on both units. Cases:
!>   fatal_default       fatal_error without a handler: message line, exit status 1
!>   fatal_handler       registered handler is called with status 1 after the message line; it
!>                       ends the run with status 7
!>   internal            internal_error: the message with its opening, exit status 1
!>   fatal_returning     handler that returns: fatal_error still ends the run, status 1
!>   fatal_reset         handler registered then removed: `error stop`, status 1
!>   status_ok           check_status with status 0 returns; the program ends normally, status 0
!>   status_errmsg       check_status with a non-zero status and a message, status 1
!>   status_no_errmsg    check_status with a non-zero status and no message, status 1
!>   status_blank_errmsg check_status with a non-zero status and a blank message, status 1
!>   status_iostat       real use: opening a missing file with iostat/iomsg, status 1
!>
!> Locals:
!>   case_name  character(64)   name of the case to run
!>   iostat     integer         iostat of the open statement in case status_iostat
!>   iomsg      character(256)  iomsg of the open statement in case status_iostat
!>   unit       integer         unit number of the open statement in case status_iostat
program test_error_cases

  use, intrinsic :: iso_fortran_env, only: output_unit
  use error_mod, only: fatal_error, internal_error, check_status, set_abort_handler, &
    & reset_abort_handler
  use test_error_handlers_mod, only: stopping_handler, returning_handler

  implicit none

  character(len=64) :: case_name
  integer :: iostat
  character(len=256) :: iomsg
  integer :: unit

  ! Read the case name.
  call get_command_argument(1, case_name)

  ! Run the case.
  select case (trim(case_name))

  case ('fatal_default')
    call fatal_error('test_routine', 'expected failure')

  case ('internal')
    call internal_error('test_routine', 'a table out of order')

  case ('fatal_handler')
    call set_abort_handler(stopping_handler)
    call fatal_error('test_routine', 'failure through the handler')

  case ('fatal_returning')
    call set_abort_handler(returning_handler)
    call fatal_error('test_routine', 'failure with a returning handler')

  case ('fatal_reset')
    call set_abort_handler(stopping_handler)
    call reset_abort_handler()
    call fatal_error('test_routine', 'failure after reset')

  case ('status_ok')
    call check_status(0, 'test_routine', 'allocate nothing', 'unused message')
    write (output_unit, '(a)') 'check_status returned'

  case ('status_errmsg')
    call check_status(5, 'test_routine', 'open file missing.cfg', 'No such file')

  case ('status_no_errmsg')
    call check_status(-3, 'test_routine', 'read line 4')

  case ('status_blank_errmsg')
    call check_status(2, 'test_routine', 'read line 5', '   ')

  case ('status_iostat')
    iomsg = ''
    open (newunit=unit, file='this_file_does_not_exist.cfg', status='old', action='read', &
      & iostat=iostat, iomsg=iomsg)
    call check_status(iostat, 'test_routine', 'open file this_file_does_not_exist.cfg', iomsg)
    write (output_unit, '(a)') 'open succeeded unexpectedly'

  case default
    write (output_unit, '(a)') 'unknown case: ' // trim(case_name)
    error stop 99

  end select

end program test_error_cases
