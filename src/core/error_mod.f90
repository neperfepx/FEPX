! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Fatal-error handling for FEPX.
!>
!> Source : FEPX (ROADMAP steps 1.1, 9.7, D94; dev/plan/terminal_output.md §3.5). Contains no
!>          numerical method. Behaviour follows v2-fixed parallel_mod.f90::par_quit after
!>          CHANGES.md C6/C10: an error detected on any one rank prints its message and ends the
!>          whole run with a non-zero exit status.
!>
!> Policy:
!>   - A fatal error is reported with fatal_error(routine, message). The message,
!>     "Error  : <routine>: <message>" ("Error  : <message>" when routine is empty), is wrapped
!>     at 72 columns as the log lines (log_mod::log_wrap). Rank 0 writes it to standard output
!>     (so that it is in the log when the output is redirected, as v2 does) and to standard
!>     error (so that it is in a batch job's error file); both units are flushed, and the program
!>     ends with exit status error_exit_status (1). fatal_error never returns. In an interactive
!>     terminal the message therefore appears twice (user decision 2026-10-08).
!>   - Another rank (k > 0) tags its message "Error  : [Rank k] ..." and writes it to standard
!>     error at once (a line is not lost when the first MPI_Abort ends the other ranks), then,
!>     when an abort handler is registered (an MPI run), waits error_stdout_delay seconds before
!>     writing it to standard output: when rank 0 detected the error too (e.g. an input error,
!>     which every rank detects), rank 0's abort ends this rank first and standard output has
!>     the message once; when only this rank detected it, the message reaches the log (Q62).
!>   - Ending the program is delegated to an abort handler. Without one (serial programs, unit
!>     tests) the program ends with `error stop`. This module does not use MPI (MPI belongs in
!>     src/parallel/, conventions §6): the MPI layer (step 2.1) registers a handler with
!>     set_abort_handler that calls MPI_Abort, so that the other ranks do not hang. If a
!>     handler returns, which it must not, fatal_error still ends the program with `error stop`.
!>   - Pure procedures (material kernels, conventions §6) cannot call fatal_error; they report
!>     failure through an integer status argument (0 = success), and the impure caller decides.
!>   - An internal consistency check (a condition that no input can cause: a program error) is
!>     reported with internal_error(routine, message), a fatal error whose message starts with
!>     "internal error (please report it): ". The manual lists the other messages, those caused
!>     by the input or the system (doc/messages.rst, ROADMAP D104).
!>   - check_status turns the stat/iostat value of an allocate, open, read, ... into a fatal
!>     error with the accompanying errmsg/iomsg text.
!>
!> Public entities:
!>   fatal_error          report a fatal error and end the run (never returns)
!>   internal_error       report a failed internal consistency check (never returns)
!>   check_status         call fatal_error if a stat/iostat value is non-zero
!>   set_abort_handler    register the procedure that ends the run after a fatal error
!>   reset_abort_handler  remove the registered handler (back to `error stop`)
!>   abort_handler_iface  abstract interface of an abort handler
!>   error_exit_status    integer  exit status of a run ended by a fatal error (1)
!>
!> Module variables:
!>   abort_handler  procedure pointer (abort_handler_iface)  registered handler, or null
module error_mod

  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit, int64
  use kinds_mod, only: ik
  use log_mod, only: log_wrap, log_get_rank, log_error_prefix, log_hold_release

  implicit none

  private

  public :: fatal_error, internal_error, check_status, set_abort_handler, reset_abort_handler
  public :: abort_handler_iface, error_exit_status

  !> Exit status of a run ended by a fatal error.
  integer, parameter :: error_exit_status = 1

  !> Wait of a rank other than 0 before it writes its message to standard output [s]: longer
  !> than rank 0's wait before MPI_Abort (parallel_mod abort_drain_seconds, 0.5 s).
  real, parameter :: error_stdout_delay = 2.0

  abstract interface
    !> Abort handler: ends the run (on every MPI rank) with the given exit status; must not
    !> return.
    !>
    !> Arguments:
    !>   exit_status  (in)  integer  exit status of the run
    subroutine abort_handler_iface(exit_status)
      implicit none
      integer, intent(in) :: exit_status
    end subroutine abort_handler_iface
  end interface

  procedure(abort_handler_iface), pointer, save :: abort_handler => null()

contains

  !> Report a fatal error and end the run. Never returns.
  !>
  !> Arguments:
  !>   routine  (in)  character  name of the routine that detected the error (empty: none shown)
  !>   message  (in)  character  description of the error (what was expected, what was found)
  !> Locals:
  !>   text     character(:), allocatable  the message with the routine and the rank tag
  !>   lines    character(:), allocatable  its wrapped lines, each ended by a new-line character
  !>   rank     integer(ik)                rank of this process (log_mod)
  !>   rank_text  character(12)            rank as text
  !>   clock_start, clock_now, clock_rate  integer(int64)  wall clock of the wait
  subroutine fatal_error(routine, message)

    character(len=*), intent(in) :: routine
    character(len=*), intent(in) :: message

    character(len=:), allocatable :: text, lines
    integer(ik) :: rank
    character(len=12) :: rank_text
    integer(int64) :: clock_start, clock_now, clock_rate

    ! Build the message: the routine (if any), tagged with the rank on ranks other than 0.
    if (len_trim(routine) > 0) then
      text = trim(routine) // ': ' // trim(message)
    else
      text = trim(message)
    end if
    rank = log_get_rank()
    if (rank /= 0_ik) then
      write (rank_text, '(i0)') rank
      text = '[Rank ' // trim(rank_text) // '] ' // text
    end if
    lines = log_wrap(log_error_prefix, '', text)

    ! Lines kept by the log's hold mode (the setup summary) first, so that they are not lost.
    call log_hold_release()

    ! Standard error at once, flushed so that the message is not lost when the run is ended.
    call write_lines(error_unit, lines)

    ! Standard output: rank 0 at once; another rank in an MPI run after the wait (see the
    ! module policy), unless rank 0's abort ends it first.
    if (rank /= 0_ik .and. associated(abort_handler)) then
      call system_clock(clock_start, clock_rate)
      do
        call system_clock(clock_now)
        if (real(clock_now - clock_start)/real(max(clock_rate, 1_int64)) >= &
          & error_stdout_delay) exit
      end do
    end if
    call write_lines(output_unit, lines)

    ! End the run through the registered handler (e.g. MPI_Abort), if any.
    if (associated(abort_handler)) call abort_handler(error_exit_status)

    ! No handler, or a handler that returned: end this process.
    error stop error_exit_status

  end subroutine fatal_error

  !> Report a failed internal consistency check (a program error, not caused by the input) as a
  !> fatal error: "internal error (please report it): <message>". Never returns.
  !>
  !> Arguments:
  !>   routine  (in)  character  name of the routine that detected the error
  !>   message  (in)  character  the check that failed
  subroutine internal_error(routine, message)

    character(len=*), intent(in) :: routine
    character(len=*), intent(in) :: message

    call fatal_error(routine, 'internal error (please report it): ' // message)

  end subroutine internal_error

  !> Call fatal_error if a stat/iostat value reports a failure (non-zero value).
  !>
  !> Arguments:
  !>   status   (in)  integer              stat= or iostat= value; 0 means success
  !>   routine  (in)  character            name of the calling routine
  !>   action   (in)  character            what was attempted, e.g. "allocate coo_s" or
  !>                                       "open file simulation.cfg"
  !>   errmsg   (in, optional) character   errmsg= or iomsg= text returned with the status
  !> Locals:
  !>   status_text  character(24)  status value as text
  subroutine check_status(status, routine, action, errmsg)

    integer, intent(in) :: status
    character(len=*), intent(in) :: routine
    character(len=*), intent(in) :: action
    character(len=*), intent(in), optional :: errmsg

    character(len=24) :: status_text

    ! Nothing to do on success.
    if (status == 0) return

    ! Report the failure, with the returned message text when there is one.
    write (status_text, '(i0)') status
    if (present(errmsg)) then
      if (len_trim(errmsg) > 0) then
        call fatal_error(routine, trim(action) // ' failed (status ' // trim(status_text) // &
          & '): ' // trim(errmsg))
      end if
    end if
    call fatal_error(routine, trim(action) // ' failed (status ' // trim(status_text) // ')')

  end subroutine check_status

  !> Write lines, each ended by a new-line character, to a unit with one write statement (no
  !> record advance: the new-line characters end the lines), then flush it: the message leaves
  !> the process in one piece, so that the lines of two ranks' messages do not interleave on a
  !> shared standard error.
  !>
  !> Arguments:
  !>   unit   (in)  integer    the unit
  !>   lines  (in)  character  the lines
  subroutine write_lines(unit, lines)

    integer, intent(in) :: unit
    character(len=*), intent(in) :: lines

    ! The whole message at once.
    write (unit, '(a)', advance='no') lines
    flush (unit)

  end subroutine write_lines

  !> Register the procedure that ends the run after a fatal error (replaces any previous one).
  !>
  !> Arguments:
  !>   handler  (in)  procedure(abort_handler_iface)  abort handler; must not return
  subroutine set_abort_handler(handler)

    procedure(abort_handler_iface) :: handler

    ! Point the module's handler at the given procedure.
    abort_handler => handler

  end subroutine set_abort_handler

  !> Remove the registered abort handler; fatal errors then end the process with `error stop`.
  subroutine reset_abort_handler()

    ! Disassociate the module's handler.
    abort_handler => null()

  end subroutine reset_abort_handler

end module error_mod
