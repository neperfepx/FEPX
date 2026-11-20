! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The launcher of the tests in a build without MPI (FEPX_ENABLE_MPI=OFF, ROADMAP D148), in place
!> of mpiexec: `fepx_serial_launch -n <ranks> <program> [<argument>...]` runs the program with
!> its arguments when <ranks> is 1, and fails otherwise (a build without MPI runs on one rank).
!>
!> Source : FEPX (test infrastructure). Contains no numerical method.
!>
!> Policy:
!>   - The test scripts start FEPX as `${MPIEXEC} ${MPIEXEC_NUMPROC} <ranks> <program> ...`;
!>     in a build without MPI, MPIEXEC is this program and MPIEXEC_NUMPROC is -n (FEPXTesting.cmake),
!>     so that the scripts run unchanged on one rank.
!>   - A number of ranks other than 1 is an error (exit status 2) with a message: the tests on
!>     several ranks are left out of a build without MPI (fepx_ranks, FEPXTesting.cmake), and a
!>     test that still asks for them is to be guarded too.
!>   - The program is run through execute_command_line, its arguments quoted (single quotes on
!>     POSIX systems, double quotes on Windows); the exit status is the program's.
!>
!> Locals:
!>   num_args    integer                    number of command-line arguments
!>   i           integer                    argument index
!>   length      integer                    length of an argument
!>   arg         character(:), allocatable  an argument
!>   command     character(:), allocatable  the command line run
!>   ranks       integer                    number of ranks asked for
!>   ios         integer                    status of the conversion of <ranks>
!>   exit_stat   integer                    exit status of the program
!>   cmd_stat    integer                    status of execute_command_line
!>   cmd_msg     character(256)             its message
!>   windows     logical                    on Windows (environment variable OS = Windows_NT)
!>   os_name     character(32)              value of the environment variable OS
program serial_launch

  use, intrinsic :: iso_fortran_env, only: error_unit

  implicit none

  integer :: num_args
  integer :: i
  integer :: length
  character(len=:), allocatable :: arg
  character(len=:), allocatable :: command
  integer :: ranks
  integer :: ios
  integer :: exit_stat
  integer :: cmd_stat
  character(len=256) :: cmd_msg
  logical :: windows
  character(len=32) :: os_name

  ! -n <ranks> <program>, at least.
  num_args = command_argument_count()
  if (num_args < 3) call usage_error('expected -n <ranks> <program> [<argument>...]')
  call get_arg(1, arg)
  if (arg /= '-n') call usage_error('expected -n <ranks> as the first arguments')
  call get_arg(2, arg)
  read (arg, *, iostat=ios) ranks
  if (ios /= 0) call usage_error('the number of ranks "' // arg // '" is not an integer')

  ! One rank only.
  if (ranks /= 1) then
    write (error_unit, '(a,i0,a)') 'fepx_serial_launch: ', ranks, ' ranks asked for, but this ' // &
      & 'FEPX is built without MPI (FEPX_ENABLE_MPI=OFF) and runs on one rank'
    stop 2, quiet=.true.
  end if

  ! The command line, each argument quoted.
  call get_environment_variable('OS', os_name)
  windows = os_name == 'Windows_NT'
  command = ''
  do i = 3, num_args
    call get_arg(i, arg)
    if (i > 3) command = command // ' '
    command = command // quoted(arg)
  end do

  ! Run it; its exit status is this program's.
  cmd_msg = ''
  call execute_command_line(command, exitstat=exit_stat, cmdstat=cmd_stat, cmdmsg=cmd_msg)
  if (cmd_stat /= 0) then
    write (error_unit, '(a)') 'fepx_serial_launch: cannot run ' // command // ': ' // trim(cmd_msg)
    stop 2, quiet=.true.
  end if
  if (exit_stat /= 0) error stop exit_stat, quiet=.true.

contains

  !> A command-line argument.
  !>
  !> Arguments:
  !>   index  (in)   integer                    its position
  !>   value  (out)  character(:), allocatable  the argument
  subroutine get_arg(index, value)

    integer, intent(in) :: index
    character(len=:), allocatable, intent(out) :: value

    call get_command_argument(index, length=length)
    allocate (character(len=length) :: value)
    call get_command_argument(index, value)

  end subroutine get_arg

  !> An argument quoted for the shell of execute_command_line: in single quotes on POSIX
  !> systems (a single quote as '\''), in double quotes on Windows.
  !>
  !> Arguments:
  !>   text    (in)  character  the argument
  !>   result  (out, result) character(:), allocatable  the quoted argument
  !> Locals:
  !>   j  integer  character index
  function quoted(text) result(res)

    character(len=*), intent(in) :: text
    character(len=:), allocatable :: res

    integer :: j

    if (windows) then
      res = '"' // text // '"'
      return
    end if
    res = "'"
    do j = 1, len(text)
      if (text(j:j) == "'") then
        res = res // "'\''"
      else
        res = res // text(j:j)
      end if
    end do
    res = res // "'"

  end function quoted

  !> Report a wrong command line and stop with exit status 2.
  !>
  !> Arguments:
  !>   message  (in)  character  what is wrong
  subroutine usage_error(message)

    character(len=*), intent(in) :: message

    write (error_unit, '(a)') 'fepx_serial_launch: ' // message
    stop 2, quiet=.true.

  end subroutine usage_error

end program serial_launch
