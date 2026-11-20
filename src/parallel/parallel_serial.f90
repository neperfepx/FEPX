! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The parallel environment of a build without MPI (FEPX_ENABLE_MPI=OFF): one rank.
!>
!> Source : FEPX (ROADMAP D148). Contains no numerical method. The same module as
!>          parallel_mod.f90 (module parallel_mod, the same public entities and arguments), which
!>          the build compiles instead of it, as io/hdf5_file_stub.f90 replaces
!>          io/hdf5_file_mod.f90 (step 9.1).
!>
!> Policy:
!>   - One rank: num_ranks = 1, my_rank = root_rank = 0. Every collective operation of
!>     parallel_mod.f90 is what MPI does on one rank: a reduction or a broadcast returns its
!>     argument, a scatter, gather or all-to-all copies its one part. The results of a run are
!>     therefore those of the MPI build on one rank, to the bit.
!>   - The checks of parallel_mod.f90 are kept (the routine is called between par_init and
!>     par_finalize, the sizes of the arguments), so that a call that is wrong in the MPI build
!>     is wrong here too.
!>   - par_neighbor_exchange: one rank has no neighbor; a call with neighbors is an internal
!>     error. The broadcasts and par_neighbor_exchange leave some arguments unused (the build
!>     compiles this file without the unused-dummy-argument warning).
!>   - par_init registers an abort handler with error_mod that ends the run with error_mod's
!>     exit status without the compiler's "ERROR STOP" line (error stop ... quiet), so that a
!>     fatal error shows the same output as in the MPI build.
!>
!> Public entities: those of parallel_mod.f90 (par_init, par_finalize, par_barrier, par_wtime,
!> par_library_version, par_sum, par_max, par_min, par_all, par_any, par_bcast, par_scatterv,
!> par_gatherv, par_alltoall, par_alltoallv, par_log_flush, par_neighbor_exchange, num_ranks,
!> my_rank, root_rank).
!>
!> Module variables:
!>   is_started  logical  .true. between par_init and par_finalize
module parallel_mod

  use, intrinsic :: iso_fortran_env, only: int64
  use kinds_mod, only: rk, ik
  use error_mod, only: internal_error, set_abort_handler, reset_abort_handler
  use log_mod, only: log_set_rank, log_buffer_take, log_lines_write

  implicit none

  private

  public :: par_init, par_finalize, par_barrier, par_wtime, par_library_version
  public :: par_sum, par_max, par_min, par_all, par_any, par_bcast, par_scatterv, par_gatherv
  public :: par_alltoall, par_alltoallv, par_log_flush, par_neighbor_exchange
  public :: num_ranks, my_rank, root_rank

  !> Rank that reads input, logs and writes output.
  integer(ik), parameter :: root_rank = 0_ik

  integer(ik), protected, save :: num_ranks = 0_ik
  integer(ik), protected, save :: my_rank = 0_ik

  logical, save :: is_started = .false.

  !> Sum over all ranks.
  interface par_sum
    module procedure par_sum_real, par_sum_real_rank1, par_sum_int, par_sum_int_rank1
  end interface par_sum

  !> Maximum over all ranks.
  interface par_max
    module procedure par_max_real, par_max_int
  end interface par_max

  !> Minimum over all ranks.
  interface par_min
    module procedure par_min_real, par_min_int
  end interface par_min

  !> Scatter variable-size parts of a root array to every rank.
  interface par_scatterv
    module procedure par_scatterv_real, par_scatterv_int
  end interface par_scatterv

  !> Gather the ranks' parts of an array on root_rank, concatenated in rank order.
  interface par_gatherv
    module procedure par_gatherv_real, par_gatherv_int
  end interface par_gatherv

  !> Broadcast from root_rank to every rank (in place).
  interface par_bcast
    module procedure par_bcast_real, par_bcast_real_rank1, par_bcast_int, par_bcast_int_rank1, &
      & par_bcast_logical, par_bcast_char
  end interface par_bcast

contains

  !> Set up the parallel environment of one rank: num_ranks = 1, my_rank = 0, the rank passed
  !> to log_mod, and the abort handler registered with error_mod. Must be called once, before
  !> any other routine of this module.
  subroutine par_init()

    ! Once per run.
    if (is_started) call internal_error('par_init', 'the parallel environment is already started')

    ! A fatal error ends the run through the handler (no "ERROR STOP" line).
    call set_abort_handler(par_abort_handler)

    ! One rank, which logs.
    num_ranks = 1_ik
    my_rank = 0_ik
    call log_set_rank(my_rank)

    is_started = .true.

  end subroutine par_init

  !> End the parallel environment at the normal end of the run. Fatal errors after this call end
  !> the process with `error stop` (no handler).
  subroutine par_finalize()

    ! It must have been started.
    call require_started('par_finalize')

    ! Remove the abort handler.
    call reset_abort_handler()
    is_started = .false.

  end subroutine par_finalize

  !> Wait until every rank has reached this call: nothing to wait for on one rank.
  subroutine par_barrier()

    ! It must have been started.
    call require_started('par_barrier')

  end subroutine par_barrier

  !> Wall-clock time [s], from an arbitrary origin (system_clock); differences of two calls
  !> give elapsed times.
  !>
  !> Returns:
  !>   seconds  real(rk)  wall-clock time [s]
  !> Locals:
  !>   count, rate  integer(int64)  clock count and counts per second
  function par_wtime() result(seconds)

    real(rk) :: seconds

    integer(int64) :: count
    integer(int64) :: rate

    ! It must have been started.
    call require_started('par_wtime')

    ! Read the clock (64-bit counts: a fine resolution that does not wrap during a run).
    call system_clock(count, rate)
    seconds = real(count, rk)/real(max(rate, 1_int64), rk)

  end function par_wtime

  !> The MPI library of the build, for the banner's "Built with" line: none.
  !>
  !> Returns:
  !>   text  character(:), allocatable  "no MPI"
  function par_library_version() result(text)

    character(len=:), allocatable :: text

    text = 'no MPI'

  end function par_library_version

  !> Sum of a real scalar over all ranks: the value itself.
  !>
  !> Arguments:
  !>   part   (in)   real(rk)  this rank's contribution
  !>   whole  (out)  real(rk)  the sum (part)
  subroutine par_sum_real(part, whole)

    real(rk), intent(in) :: part
    real(rk), intent(out) :: whole

    call require_started('par_sum')
    whole = part

  end subroutine par_sum_real

  !> Componentwise sum of a real array over all ranks: the array itself.
  !>
  !> Arguments:
  !>   part   (in)   real(rk)(n)  this rank's contributions
  !>   whole  (out)  real(rk)(n)  the sums (part)
  subroutine par_sum_real_rank1(part, whole)

    real(rk), intent(in) :: part(:)
    real(rk), intent(out) :: whole(:)

    call require_started('par_sum')
    call check_sizes(size(part), size(whole), 'par_sum')
    whole = part

  end subroutine par_sum_real_rank1

  !> Sum of an integer scalar over all ranks: the value itself.
  !>
  !> Arguments:
  !>   part   (in)   integer(ik)  this rank's contribution
  !>   whole  (out)  integer(ik)  the sum (part)
  subroutine par_sum_int(part, whole)

    integer(ik), intent(in) :: part
    integer(ik), intent(out) :: whole

    call require_started('par_sum')
    whole = part

  end subroutine par_sum_int

  !> Componentwise sum of an integer array over all ranks: the array itself.
  !>
  !> Arguments:
  !>   part   (in)   integer(ik)(n)  this rank's contributions
  !>   whole  (out)  integer(ik)(n)  the sums (part)
  subroutine par_sum_int_rank1(part, whole)

    integer(ik), intent(in) :: part(:)
    integer(ik), intent(out) :: whole(:)

    call require_started('par_sum')
    call check_sizes(size(part), size(whole), 'par_sum')
    whole = part

  end subroutine par_sum_int_rank1

  !> Maximum of a real scalar over all ranks: the value itself.
  !>
  !> Arguments:
  !>   part     (in)   real(rk)  this rank's value
  !>   maxpart  (out)  real(rk)  the maximum (part)
  subroutine par_max_real(part, maxpart)

    real(rk), intent(in) :: part
    real(rk), intent(out) :: maxpart

    call require_started('par_max')
    maxpart = part

  end subroutine par_max_real

  !> Maximum of an integer scalar over all ranks: the value itself.
  !>
  !> Arguments:
  !>   part     (in)   integer(ik)  this rank's value
  !>   maxpart  (out)  integer(ik)  the maximum (part)
  subroutine par_max_int(part, maxpart)

    integer(ik), intent(in) :: part
    integer(ik), intent(out) :: maxpart

    call require_started('par_max')
    maxpart = part

  end subroutine par_max_int

  !> Minimum of a real scalar over all ranks: the value itself.
  !>
  !> Arguments:
  !>   part     (in)   real(rk)  this rank's value
  !>   minpart  (out)  real(rk)  the minimum (part)
  subroutine par_min_real(part, minpart)

    real(rk), intent(in) :: part
    real(rk), intent(out) :: minpart

    call require_started('par_min')
    minpart = part

  end subroutine par_min_real

  !> Minimum of an integer scalar over all ranks: the value itself.
  !>
  !> Arguments:
  !>   part     (in)   integer(ik)  this rank's value
  !>   minpart  (out)  integer(ik)  the minimum (part)
  subroutine par_min_int(part, minpart)

    integer(ik), intent(in) :: part
    integer(ik), intent(out) :: minpart

    call require_started('par_min')
    minpart = part

  end subroutine par_min_int

  !> Logical AND over all ranks: the value itself.
  !>
  !> Arguments:
  !>   part   (in)   logical  this rank's value
  !>   whole  (out)  logical  the AND (part)
  subroutine par_all(part, whole)

    logical, intent(in) :: part
    logical, intent(out) :: whole

    call require_started('par_all')
    whole = part

  end subroutine par_all

  !> Logical OR over all ranks: the value itself.
  !>
  !> Arguments:
  !>   part   (in)   logical  this rank's value
  !>   whole  (out)  logical  the OR (part)
  subroutine par_any(part, whole)

    logical, intent(in) :: part
    logical, intent(out) :: whole

    call require_started('par_any')
    whole = part

  end subroutine par_any

  !> Broadcast a real scalar from root_rank: unchanged on one rank.
  !>
  !> Arguments:
  !>   value  (inout)  real(rk)  the value (unchanged)
  subroutine par_bcast_real(value)

    real(rk), intent(inout) :: value

    call require_started('par_bcast')

  end subroutine par_bcast_real

  !> Broadcast a real array from root_rank: unchanged on one rank.
  !>
  !> Arguments:
  !>   values  (inout)  real(rk)(n)  the values (unchanged)
  subroutine par_bcast_real_rank1(values)

    real(rk), intent(inout) :: values(:)

    call require_started('par_bcast')

  end subroutine par_bcast_real_rank1

  !> Broadcast an integer scalar from root_rank: unchanged on one rank.
  !>
  !> Arguments:
  !>   value  (inout)  integer(ik)  the value (unchanged)
  subroutine par_bcast_int(value)

    integer(ik), intent(inout) :: value

    call require_started('par_bcast')

  end subroutine par_bcast_int

  !> Broadcast an integer array from root_rank: unchanged on one rank.
  !>
  !> Arguments:
  !>   values  (inout)  integer(ik)(n)  the values (unchanged)
  subroutine par_bcast_int_rank1(values)

    integer(ik), intent(inout) :: values(:)

    call require_started('par_bcast')

  end subroutine par_bcast_int_rank1

  !> Broadcast a logical scalar from root_rank: unchanged on one rank.
  !>
  !> Arguments:
  !>   value  (inout)  logical  the value (unchanged)
  subroutine par_bcast_logical(value)

    logical, intent(inout) :: value

    call require_started('par_bcast')

  end subroutine par_bcast_logical

  !> Broadcast a character string from root_rank: unchanged on one rank.
  !>
  !> Arguments:
  !>   value  (inout)  character(*)  the string (unchanged)
  subroutine par_bcast_char(value)

    character(len=*), intent(inout) :: value

    call require_started('par_bcast')

  end subroutine par_bcast_char

  !> Send each rank its part of a real array held by root_rank: on one rank, the whole array.
  !>
  !> Arguments:
  !>   send    (in)   real(rk)(:)               the parts of all ranks (one part)
  !>   counts  (in)   integer(ik)(:)            size of each rank's part (one count, size(send))
  !>   recv    (out)  real(rk)(:), allocatable  this rank's part (a copy of send)
  subroutine par_scatterv_real(send, counts, recv)

    real(rk), intent(in) :: send(:)
    integer(ik), intent(in) :: counts(:)
    real(rk), allocatable, intent(out) :: recv(:)

    call require_started('par_scatterv')
    call check_counts(counts, size(send), 'par_scatterv')
    recv = send

  end subroutine par_scatterv_real

  !> Send each rank its part of an integer array held by root_rank: on one rank, the whole
  !> array.
  !>
  !> Arguments:
  !>   send    (in)   integer(ik)(:)               the parts of all ranks (one part)
  !>   counts  (in)   integer(ik)(:)               size of each rank's part (one count)
  !>   recv    (out)  integer(ik)(:), allocatable  this rank's part (a copy of send)
  subroutine par_scatterv_int(send, counts, recv)

    integer(ik), intent(in) :: send(:)
    integer(ik), intent(in) :: counts(:)
    integer(ik), allocatable, intent(out) :: recv(:)

    call require_started('par_scatterv')
    call check_counts(counts, size(send), 'par_scatterv')
    recv = send

  end subroutine par_scatterv_int

  !> Gather each rank's part of a real array on root_rank: on one rank, the part itself.
  !>
  !> Arguments:
  !>   send  (in)   real(rk)(:)               this rank's part
  !>   recv  (out)  real(rk)(:), allocatable  the parts of all ranks (a copy of send)
  subroutine par_gatherv_real(send, recv)

    real(rk), intent(in) :: send(:)
    real(rk), allocatable, intent(out) :: recv(:)

    call require_started('par_gatherv')
    recv = send

  end subroutine par_gatherv_real

  !> Gather each rank's part of an integer array on root_rank: on one rank, the part itself.
  !>
  !> Arguments:
  !>   send  (in)   integer(ik)(:)               this rank's part
  !>   recv  (out)  integer(ik)(:), allocatable  the parts of all ranks (a copy of send)
  subroutine par_gatherv_int(send, recv)

    integer(ik), intent(in) :: send(:)
    integer(ik), allocatable, intent(out) :: recv(:)

    call require_started('par_gatherv')
    recv = send

  end subroutine par_gatherv_int

  !> Every rank sends one integer to every rank: on one rank, the value to itself.
  !>
  !> Arguments:
  !>   send  (in)   integer(ik)(1)  value for each rank
  !>   recv  (out)  integer(ik)(1)  value from each rank (send)
  subroutine par_alltoall(send, recv)

    integer(ik), intent(in) :: send(:)
    integer(ik), intent(out) :: recv(:)

    call require_started('par_alltoall')
    if (size(send) /= num_ranks .or. size(recv) /= num_ranks) then
      call internal_error('par_alltoall', 'send and recv must have one entry per rank')
    end if
    recv = send

  end subroutine par_alltoall

  !> Every rank sends a variable number of integers to every rank: on one rank, its values to
  !> itself.
  !>
  !> Arguments:
  !>   send         (in)   integer(ik)(:)               values for all ranks
  !>   send_counts  (in)   integer(ik)(1)               number of values for each rank
  !>   recv         (out)  integer(ik)(:), allocatable  values from all ranks (send)
  !>   recv_counts  (in)   integer(ik)(1)               number of values from each rank (as
  !>                                                    sent)
  subroutine par_alltoallv(send, send_counts, recv, recv_counts)

    integer(ik), intent(in) :: send(:)
    integer(ik), intent(in) :: send_counts(:)
    integer(ik), allocatable, intent(out) :: recv(:)
    integer(ik), intent(in) :: recv_counts(:)

    call require_started('par_alltoallv')
    if (size(send_counts) /= num_ranks .or. size(recv_counts) /= num_ranks) then
      call internal_error('par_alltoallv', 'expected one send and one receive count per rank')
    end if
    if (any(send_counts < 0_ik) .or. any(recv_counts < 0_ik) .or. &
      & sum(int(send_counts)) /= size(send) .or. any(recv_counts /= send_counts)) then
      call internal_error('par_alltoallv', 'invalid counts')
    end if
    recv = send

  end subroutine par_alltoallv

  !> Write the log lines buffered for all ranks (log_mod, all_ranks messages) and empty the
  !> buffer: on one rank, its own lines.
  !>
  !> Locals:
  !>   text  character(:), allocatable  the buffered lines
  subroutine par_log_flush()

    character(len=:), allocatable :: text

    call require_started('par_log_flush')
    call log_buffer_take(text)
    if (len(text) > 0) call log_lines_write(text)

  end subroutine par_log_flush

  !> Send real values to, and receive real values from, a list of neighbor ranks: one rank has
  !> no neighbor, and a call with neighbors is an internal error (the arguments are those of
  !> parallel_mod.f90).
  !>
  !> Arguments:
  !>   neighbor     (in)     integer(ik)(:)  neighbor ranks (none)
  !>   block        (in)     integer(ik)     values per entry (not used)
  !>   send_buf     (in)     real(rk)(:)     the values to send (not used)
  !>   send_offset  (in)     integer(ik)(:)  the parts of send_buf (one offset)
  !>   recv_buf     (inout)  real(rk)(:)     the values received (unchanged)
  !>   recv_offset  (in)     integer(ik)(:)  the parts of recv_buf (one offset)
  !>   tag          (in)     integer         message tag (not used)
  !>   routine      (in)     character       the calling exchange, for the error messages
  subroutine par_neighbor_exchange(neighbor, block, send_buf, send_offset, recv_buf, &
    & recv_offset, tag, routine)

    integer(ik), intent(in) :: neighbor(:)
    integer(ik), intent(in) :: block
    real(rk), intent(in), asynchronous, contiguous :: send_buf(:)
    integer(ik), intent(in) :: send_offset(:)
    real(rk), intent(inout), asynchronous, contiguous :: recv_buf(:)
    integer(ik), intent(in) :: recv_offset(:)
    integer, intent(in) :: tag
    character(len=*), intent(in) :: routine

    ! Started; one part per neighbor in each direction (as parallel_mod.f90); no neighbor.
    call require_started(routine)
    if (size(send_offset) /= size(neighbor) + 1 .or. size(recv_offset) /= size(neighbor) + 1) &
      & call internal_error(routine, 'expected one offset per neighbor, plus one')
    if (size(neighbor) > 0) call internal_error(routine, 'a build without MPI has no neighbor rank')

  end subroutine par_neighbor_exchange

  !> Abort handler registered with error_mod by par_init: ends the run with the given exit
  !> status, without the compiler's "ERROR STOP" line (error_mod has written the message).
  !> Does not return.
  !>
  !> Arguments:
  !>   exit_status  (in)  integer  exit status of the run (error_mod's error_exit_status)
  subroutine par_abort_handler(exit_status)

    integer, intent(in) :: exit_status

    error stop exit_status, quiet=.true.

  end subroutine par_abort_handler

  !> Fatal error if the parallel environment has not been started by par_init (or has been
  !> ended).
  !>
  !> Arguments:
  !>   routine  (in)  character  name of the calling public routine
  subroutine require_started(routine)

    character(len=*), intent(in) :: routine

    if (.not. is_started) call internal_error(routine, &
      & 'the parallel environment is not started (call par_init first)')

  end subroutine require_started

  !> Fatal error if the result array of a reduction does not have the size of the contribution.
  !>
  !> Arguments:
  !>   size_part   (in)  integer    size of this rank's contribution
  !>   size_whole  (in)  integer    size of the result array
  !>   routine     (in)  character  name of the calling public routine
  !> Locals:
  !>   text        character(64)  the two sizes as text
  subroutine check_sizes(size_part, size_whole, routine)

    integer, intent(in) :: size_part
    integer, intent(in) :: size_whole
    character(len=*), intent(in) :: routine

    character(len=64) :: text

    ! Nothing to do when the sizes agree.
    if (size_part == size_whole) return

    ! Report both sizes.
    write (text, '(a,i0,a,i0)') 'contribution size ', size_part, ', result size ', size_whole
    call internal_error(routine, 'array sizes differ: ' // trim(text))

  end subroutine check_sizes

  !> Fatal error if the counts of a scatter are not one non-negative count per rank, adding up
  !> to the send size (as parallel_mod.f90's scatterv_layout).
  !>
  !> Arguments:
  !>   counts     (in)  integer(ik)(:)  size of each rank's part
  !>   send_size  (in)  integer         size of the send array
  !>   routine    (in)  character       name of the calling public routine
  subroutine check_counts(counts, send_size, routine)

    integer(ik), intent(in) :: counts(:)
    integer, intent(in) :: send_size
    character(len=*), intent(in) :: routine

    if (size(counts) /= num_ranks .or. any(counts < 0_ik)) then
      call internal_error(routine, 'expected one non-negative count per rank')
    end if
    if (sum(int(counts)) /= send_size) then
      call internal_error(routine, 'the counts do not add up to the send size')
    end if

  end subroutine check_counts

end module parallel_mod
