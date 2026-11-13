! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> MPI environment and collective operations for FEPX: the only module that uses MPI.
!>
!> Source : port of v2-fixed parallel_mod.f90 (par_init, par_sum, par_max, par_min, par_barrier,
!>          par_quit's abort path after CHANGES.md C6/C10), rewritten with `use mpi_f08`
!>          (ROADMAP D5, conventions §6). Contains no numerical method. v2's par_partition /
!>          mre_decomp1d (equal contiguous blocks) are not used (D108); its gather routines became par_gatherv (step 4.3).
!>
!> Policy:
!>   - All communication is on MPI_COMM_WORLD; the number of ranks and the rank of this process
!>     are the public, read-only variables num_ranks and my_rank (0-based); root_rank (0) is the
!>     rank that reads input, logs and writes output.
!>   - par_init starts MPI, determines the MPI datatypes that match the FEPX kinds rk and ik
!>     (MPI_Type_match_size, so that no kind is assumed to be MPI's default), tells log_mod the
!>     rank, and registers an abort handler with error_mod: from then on, a fatal_error on any
!>     one rank ends the whole run through MPI_Abort with error_mod's exit status, so that the
!>     other ranks do not hang (as v2-fixed after C6/C10). par_finalize removes the handler.
!>   - Every MPI call's ierror is checked; a failure is a fatal error with MPI's error string.
!>     Every routine other than par_init checks that MPI has been started by par_init.
!>   - The reductions are collective: every rank must call them, in the same order. Arrays passed
!>     to a reduction or broadcast must have the same size on every rank. MPI implementations
!>     return the same allreduce result on every rank (MPI standard, advice to implementors);
!>     decisions that must agree on all ranks (e.g. convergence flags, CHANGES.md C16) should
!>     nevertheless be taken with par_all / par_any.
!>
!> Public entities:
!>   par_init       start MPI, set num_ranks/my_rank, register the abort handler
!>   par_finalize   end MPI (normal end of the run)
!>   par_barrier    wait until every rank has reached this call
!>   par_wtime      wall-clock time [s] (MPI_Wtime)
!>   par_sum        sum over all ranks (real(rk) or integer(ik); scalar or rank-1 array)
!>   par_max        maximum over all ranks (real(rk) or integer(ik) scalar)
!>   par_min        minimum over all ranks (real(rk) or integer(ik) scalar)
!>   par_all        logical AND over all ranks
!>   par_any        logical OR over all ranks
!>   par_bcast      broadcast from root_rank (real(rk), integer(ik): scalar or rank-1 array;
!>                  logical and character scalars)
!>   par_scatterv   send each rank its part of a rank-1 array held by root_rank (real(rk) or
!>                  integer(ik); the parts may have different sizes; step 2.2)
!>   par_gatherv    gather the ranks' parts of a rank-1 array on root_rank, in rank order
!>                  (real(rk) or integer(ik); step 4.3)
!>   par_alltoall   every rank sends one integer(ik) to every rank (step 2.3)
!>   par_alltoallv  every rank sends a variable number of integer(ik) to every rank (step 2.3)
!>   par_log_flush  rank 0 writes the log lines every rank buffered (log_mod all_ranks; step 9.7)
!>   num_ranks      integer(ik)  number of ranks (0 before par_init)
!>   my_rank        integer(ik)  rank of this process, 0-based
!>   root_rank      integer(ik)  rank that reads input, logs and writes output (0)
!>   mpi_type_rk    type(MPI_Datatype)  MPI datatype matching real(rk) (read-only; for the other
!>                  modules of src/parallel/, e.g. halo_mod)
!>   par_check_mpi  fatal error if an MPI call failed (for the other modules of src/parallel/)
!>
!> Module variables:
!>   mpi_type_rk  type(MPI_Datatype)  MPI datatype matching real(rk) (public, protected)
!>   mpi_type_ik  type(MPI_Datatype)  MPI datatype matching integer(ik)
!>   is_started   logical             .true. between par_init and par_finalize
module parallel_mod

  use mpi_f08, only: MPI_Init, MPI_Finalize, MPI_Initialized, MPI_Comm_size, MPI_Comm_rank, &
    & MPI_Abort, MPI_Barrier, MPI_Wtime, MPI_Allreduce, MPI_Bcast, MPI_Scatter, MPI_Scatterv, &
    & MPI_Alltoall, MPI_Alltoallv, MPI_Gather, MPI_Gatherv, &
    & MPI_Type_match_size, MPI_INTEGER, &
    & MPI_Error_string, MPI_Datatype, MPI_Op, MPI_COMM_WORLD, MPI_SUCCESS, MPI_SUM, MPI_MAX, &
    & MPI_MIN, MPI_LAND, MPI_LOR, MPI_LOGICAL, MPI_CHARACTER, MPI_TYPECLASS_REAL, &
    & MPI_TYPECLASS_INTEGER, MPI_MAX_ERROR_STRING, MPI_Get_library_version, &
    & MPI_MAX_LIBRARY_VERSION_STRING
  use kinds_mod, only: rk, ik
  use error_mod, only: fatal_error, internal_error, set_abort_handler, reset_abort_handler
  use log_mod, only: log_set_rank, log_buffer_take, log_lines_write

  implicit none

  private

  public :: par_init, par_finalize, par_barrier, par_wtime, par_library_version
  public :: par_sum, par_max, par_min, par_all, par_any, par_bcast, par_scatterv, par_gatherv
  public :: par_alltoall, par_alltoallv, par_log_flush
  public :: num_ranks, my_rank, root_rank
  public :: mpi_type_rk, par_check_mpi

  !> Rank that reads input, logs and writes output.
  integer(ik), parameter :: root_rank = 0_ik

  ! Wait before MPI_Abort, so that the launcher forwards the error message (par_abort_handler).
  real(rk), parameter :: abort_drain_seconds = 0.5_rk

  integer(ik), protected, save :: num_ranks = 0_ik
  integer(ik), protected, save :: my_rank = 0_ik

  type(MPI_Datatype), protected, save :: mpi_type_rk
  type(MPI_Datatype), save :: mpi_type_ik
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

  !> Broadcast from root_rank to every rank (in place).
  !> Gather the ranks' parts of an array on root_rank, concatenated in rank order (real or
  !> integer(ik), rank 1).
  interface par_gatherv
    module procedure par_gatherv_real, par_gatherv_int
  end interface par_gatherv

  interface par_bcast
    module procedure par_bcast_real, par_bcast_real_rank1, par_bcast_int, par_bcast_int_rank1, &
      & par_bcast_logical, par_bcast_char
  end interface par_bcast

contains

  !> Start MPI and set up the FEPX parallel environment (v2 par_init): set num_ranks and
  !> my_rank, find the MPI datatypes of rk and ik, pass the rank to log_mod, and register the
  !> MPI abort handler with error_mod. Must be called once, by every rank, before any other
  !> routine of this module.
  !>
  !> Locals:
  !>   ierror       integer      MPI error code
  !>   already      logical      whether MPI was already initialized
  !>   size_ranks   integer      number of ranks, as returned by MPI
  !>   rank_mpi     integer      rank of this process, as returned by MPI
  subroutine par_init()

    integer :: ierror
    logical :: already
    integer :: size_ranks
    integer :: rank_mpi

    ! MPI may be initialized only once per run.
    call MPI_Initialized(already, ierror)
    call check_mpi(ierror, 'par_init', 'MPI_Initialized')
    if (already) call internal_error('par_init', 'MPI is already initialized')

    ! Start MPI.
    call MPI_Init(ierror)
    call check_mpi(ierror, 'par_init', 'MPI_Init')

    ! From now on, a fatal error on any rank ends the whole run (C6/C10).
    call set_abort_handler(par_abort_handler)

    ! Number of ranks and rank of this process.
    call MPI_Comm_size(MPI_COMM_WORLD, size_ranks, ierror)
    call check_mpi(ierror, 'par_init', 'MPI_Comm_size')
    call MPI_Comm_rank(MPI_COMM_WORLD, rank_mpi, ierror)
    call check_mpi(ierror, 'par_init', 'MPI_Comm_rank')
    num_ranks = int(size_ranks, ik)
    my_rank = int(rank_mpi, ik)

    ! MPI datatypes matching the FEPX kinds (sizes in bytes).
    call MPI_Type_match_size(MPI_TYPECLASS_REAL, storage_size(1.0_rk)/8, mpi_type_rk, ierror)
    call check_mpi(ierror, 'par_init', 'MPI_Type_match_size (real)')
    call MPI_Type_match_size(MPI_TYPECLASS_INTEGER, storage_size(1_ik)/8, mpi_type_ik, ierror)
    call check_mpi(ierror, 'par_init', 'MPI_Type_match_size (integer)')

    ! Log from the root rank only.
    call log_set_rank(my_rank)

    is_started = .true.

  end subroutine par_init

  !> End MPI at the normal end of the run (v2 par_quit, "finalize" path). Every rank must call it.
  !> Fatal errors after this call end the process with `error stop` (no handler).
  !>
  !> Locals:
  !>   ierror  integer  MPI error code
  subroutine par_finalize()

    integer :: ierror

    ! MPI must have been started.
    call require_started('par_finalize')

    ! Remove the abort handler (MPI_Abort may not be called after MPI_Finalize), then end MPI.
    call reset_abort_handler()
    is_started = .false.
    call MPI_Finalize(ierror)
    call check_mpi(ierror, 'par_finalize', 'MPI_Finalize')

  end subroutine par_finalize

  !> Wait until every rank has reached this call (v2 par_barrier).
  !>
  !> Locals:
  !>   ierror  integer  MPI error code
  subroutine par_barrier()

    integer :: ierror

    ! MPI must have been started.
    call require_started('par_barrier')

    ! Synchronize all ranks.
    call MPI_Barrier(MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_barrier', 'MPI_Barrier')

  end subroutine par_barrier

  !> Wall-clock time of this rank [s], from an arbitrary origin (MPI_Wtime); differences of two
  !> calls on the same rank give elapsed times.
  !>
  !> Returns:
  !>   seconds  real(rk)  wall-clock time [s]
  function par_wtime() result(seconds)

    real(rk) :: seconds

    ! MPI must have been started.
    call require_started('par_wtime')

    ! Read the MPI wall clock.
    seconds = real(MPI_Wtime(), rk)

  end function par_wtime

  !> The MPI library's name and version, the first part of its version string (up to the
  !> first comma; MPI_Get_library_version), as name and number: "Open MPI v5.0.10" gives
  !> "Open MPI 5.0.10", "MPICH Version: 4.3.2" gives "MPICH 4.3.2" (the banner's "Built
  !> with" line, steps 9.1, 9.7).
  !>
  !> Returns:
  !>   text     character(:), allocatable  the name and version
  !> Locals:
  !>   buffer   character(MPI_MAX_LIBRARY_VERSION_STRING)  the version string
  !>   length, ierror  integer                its length, MPI status
  !>   comma    integer                       position of its first comma (0: none)
  !>   i        integer                       position in the text
  function par_library_version() result(text)

    character(len=:), allocatable :: text

    character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: buffer
    integer :: length, ierror, comma, i

    ! The version string (callable before MPI_Init), its first line up to a comma.
    call MPI_Get_library_version(buffer, length, ierror)
    call check_mpi(ierror, 'par_library_version', 'MPI_Get_library_version')
    text = buffer(1:length)
    if (index(text, new_line('a')) > 0) text = text(1:index(text, new_line('a')) - 1)
    comma = index(text, ',')
    if (comma > 0) text = text(1:comma - 1)
    ! Name and number only: tabs as blanks, without "Version:" and the "v" before the number,
    ! single blanks.
    do i = 1, len(text)
      if (text(i:i) == achar(9)) text(i:i) = ' '
    end do
    i = index(text, 'Version:')
    if (i > 0) text = text(1:i - 1) // text(i + 8:)
    i = index(text, ' v')
    if (i > 0 .and. i + 2 <= len(text)) then
      if (verify(text(i + 2:i + 2), '0123456789') == 0) text = text(1:i) // text(i + 2:)
    end if
    i = index(text, '  ')
    do while (i > 0)
      text = text(1:i) // text(i + 2:)
      i = index(text, '  ')
    end do
    text = trim(adjustl(text))

  end function par_library_version

  !> Sum of a real scalar over all ranks (v2 par_sum).
  !>
  !> Arguments:
  !>   part   (in)   real(rk)  this rank's contribution
  !>   whole  (out)  real(rk)  sum of the contributions of all ranks (same on every rank)
  subroutine par_sum_real(part, whole)

    real(rk), intent(in) :: part
    real(rk), intent(out) :: whole

    ! Reduce.
    call allreduce_real(part, whole, 1, MPI_SUM, 'par_sum')

  end subroutine par_sum_real

  !> Componentwise sum of a real array over all ranks.
  !>
  !> Arguments:
  !>   part   (in)   real(rk)(n)  this rank's contributions (same n on every rank)
  !>   whole  (out)  real(rk)(n)  componentwise sums over all ranks
  subroutine par_sum_real_rank1(part, whole)

    real(rk), intent(in) :: part(:)
    real(rk), intent(out) :: whole(:)

    ! The result must have the size of the contribution.
    call check_sizes(size(part), size(whole), 'par_sum')

    ! Reduce.
    call allreduce_real(part, whole, size(part), MPI_SUM, 'par_sum')

  end subroutine par_sum_real_rank1

  !> Sum of an integer scalar over all ranks.
  !>
  !> Arguments:
  !>   part   (in)   integer(ik)  this rank's contribution
  !>   whole  (out)  integer(ik)  sum over all ranks (same on every rank)
  subroutine par_sum_int(part, whole)

    integer(ik), intent(in) :: part
    integer(ik), intent(out) :: whole

    ! Reduce.
    call allreduce_int(part, whole, 1, MPI_SUM, 'par_sum')

  end subroutine par_sum_int

  !> Componentwise sum of an integer array over all ranks.
  !>
  !> Arguments:
  !>   part   (in)   integer(ik)(n)  this rank's contributions (same n on every rank)
  !>   whole  (out)  integer(ik)(n)  componentwise sums over all ranks
  subroutine par_sum_int_rank1(part, whole)

    integer(ik), intent(in) :: part(:)
    integer(ik), intent(out) :: whole(:)

    ! The result must have the size of the contribution.
    call check_sizes(size(part), size(whole), 'par_sum')

    ! Reduce.
    call allreduce_int(part, whole, size(part), MPI_SUM, 'par_sum')

  end subroutine par_sum_int_rank1

  !> Maximum of a real scalar over all ranks (v2 par_max).
  !>
  !> Arguments:
  !>   part     (in)   real(rk)  this rank's value
  !>   maxpart  (out)  real(rk)  maximum over all ranks
  subroutine par_max_real(part, maxpart)

    real(rk), intent(in) :: part
    real(rk), intent(out) :: maxpart

    ! Reduce.
    call allreduce_real(part, maxpart, 1, MPI_MAX, 'par_max')

  end subroutine par_max_real

  !> Maximum of an integer scalar over all ranks.
  !>
  !> Arguments:
  !>   part     (in)   integer(ik)  this rank's value
  !>   maxpart  (out)  integer(ik)  maximum over all ranks
  subroutine par_max_int(part, maxpart)

    integer(ik), intent(in) :: part
    integer(ik), intent(out) :: maxpart

    ! Reduce.
    call allreduce_int(part, maxpart, 1, MPI_MAX, 'par_max')

  end subroutine par_max_int

  !> Minimum of a real scalar over all ranks (v2 par_min).
  !>
  !> Arguments:
  !>   part     (in)   real(rk)  this rank's value
  !>   minpart  (out)  real(rk)  minimum over all ranks
  subroutine par_min_real(part, minpart)

    real(rk), intent(in) :: part
    real(rk), intent(out) :: minpart

    ! Reduce.
    call allreduce_real(part, minpart, 1, MPI_MIN, 'par_min')

  end subroutine par_min_real

  !> Minimum of an integer scalar over all ranks.
  !>
  !> Arguments:
  !>   part     (in)   integer(ik)  this rank's value
  !>   minpart  (out)  integer(ik)  minimum over all ranks
  subroutine par_min_int(part, minpart)

    integer(ik), intent(in) :: part
    integer(ik), intent(out) :: minpart

    ! Reduce.
    call allreduce_int(part, minpart, 1, MPI_MIN, 'par_min')

  end subroutine par_min_int

  !> Logical AND over all ranks: .true. only if part is .true. on every rank.
  !>
  !> Arguments:
  !>   part   (in)   logical  this rank's value
  !>   whole  (out)  logical  AND over all ranks (same on every rank)
  !> Locals:
  !>   ierror  integer  MPI error code
  subroutine par_all(part, whole)

    logical, intent(in) :: part
    logical, intent(out) :: whole

    integer :: ierror

    ! MPI must have been started.
    call require_started('par_all')

    ! Reduce with logical AND.
    call MPI_Allreduce(part, whole, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_all', 'MPI_Allreduce')

  end subroutine par_all

  !> Logical OR over all ranks: .true. if part is .true. on at least one rank.
  !>
  !> Arguments:
  !>   part   (in)   logical  this rank's value
  !>   whole  (out)  logical  OR over all ranks (same on every rank)
  !> Locals:
  !>   ierror  integer  MPI error code
  subroutine par_any(part, whole)

    logical, intent(in) :: part
    logical, intent(out) :: whole

    integer :: ierror

    ! MPI must have been started.
    call require_started('par_any')

    ! Reduce with logical OR.
    call MPI_Allreduce(part, whole, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_any', 'MPI_Allreduce')

  end subroutine par_any

  !> Broadcast a real scalar from root_rank.
  !>
  !> Arguments:
  !>   value  (inout)  real(rk)  in: the value on root_rank (ignored elsewhere);
  !>                             out: root_rank's value, on every rank
  !> Locals:
  !>   ierror  integer  MPI error code
  subroutine par_bcast_real(value)

    real(rk), intent(inout) :: value

    integer :: ierror

    ! MPI must have been started.
    call require_started('par_bcast')

    ! Broadcast.
    call MPI_Bcast(value, 1, mpi_type_rk, int(root_rank), MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_bcast', 'MPI_Bcast')

  end subroutine par_bcast_real

  !> Broadcast a real array from root_rank.
  !>
  !> Arguments:
  !>   values  (inout)  real(rk)(n)  as par_bcast_real, componentwise (same n on every rank)
  !> Locals:
  !>   ierror  integer  MPI error code
  subroutine par_bcast_real_rank1(values)

    real(rk), intent(inout) :: values(:)

    integer :: ierror

    ! MPI must have been started.
    call require_started('par_bcast')

    ! Broadcast.
    call MPI_Bcast(values, size(values), mpi_type_rk, int(root_rank), MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_bcast', 'MPI_Bcast')

  end subroutine par_bcast_real_rank1

  !> Broadcast an integer scalar from root_rank.
  !>
  !> Arguments:
  !>   value  (inout)  integer(ik)  as par_bcast_real
  !> Locals:
  !>   ierror  integer  MPI error code
  subroutine par_bcast_int(value)

    integer(ik), intent(inout) :: value

    integer :: ierror

    ! MPI must have been started.
    call require_started('par_bcast')

    ! Broadcast.
    call MPI_Bcast(value, 1, mpi_type_ik, int(root_rank), MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_bcast', 'MPI_Bcast')

  end subroutine par_bcast_int

  !> Broadcast an integer array from root_rank.
  !>
  !> Arguments:
  !>   values  (inout)  integer(ik)(n)  as par_bcast_real, componentwise (same n on every rank)
  !> Locals:
  !>   ierror  integer  MPI error code
  subroutine par_bcast_int_rank1(values)

    integer(ik), intent(inout) :: values(:)

    integer :: ierror

    ! MPI must have been started.
    call require_started('par_bcast')

    ! Broadcast.
    call MPI_Bcast(values, size(values), mpi_type_ik, int(root_rank), MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_bcast', 'MPI_Bcast')

  end subroutine par_bcast_int_rank1

  !> Broadcast a logical scalar from root_rank.
  !>
  !> Arguments:
  !>   value  (inout)  logical  as par_bcast_real
  !> Locals:
  !>   ierror  integer  MPI error code
  subroutine par_bcast_logical(value)

    logical, intent(inout) :: value

    integer :: ierror

    ! MPI must have been started.
    call require_started('par_bcast')

    ! Broadcast.
    call MPI_Bcast(value, 1, MPI_LOGICAL, int(root_rank), MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_bcast', 'MPI_Bcast')

  end subroutine par_bcast_logical

  !> Broadcast a character string from root_rank (all len(value) characters, trailing blanks
  !> included).
  !>
  !> Arguments:
  !>   value  (inout)  character(*)  as par_bcast_real (same length on every rank)
  !> Locals:
  !>   ierror  integer  MPI error code
  subroutine par_bcast_char(value)

    character(len=*), intent(inout) :: value

    integer :: ierror

    ! MPI must have been started.
    call require_started('par_bcast')

    ! Broadcast.
    call MPI_Bcast(value, len(value), MPI_CHARACTER, int(root_rank), MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_bcast', 'MPI_Bcast')

  end subroutine par_bcast_char

  !> Send each rank its part of a real array held by root_rank (MPI_Scatter of the counts, then
  !> MPI_Scatterv). The parts are consecutive in send, in rank order.
  !>
  !> Arguments:
  !>   send    (in)   real(rk)(:)               on root_rank: the parts of all ranks, concatenated
  !>                                            in rank order (ignored elsewhere)
  !>   counts  (in)   integer(ik)(:)            on root_rank: size of each rank's part,
  !>                                            num_ranks entries, summing to size(send)
  !>                                            (ignored elsewhere)
  !>   recv    (out)  real(rk)(:), allocatable  this rank's part (allocated to its size)
  !> Locals:
  !>   counts_mpi  integer(:)   counts as MPI integers (root)
  !>   displs      integer(:)   offset of each part in send (root)
  !>   my_count    integer      size of this rank's part
  !>   ierror      integer      MPI error code
  subroutine par_scatterv_real(send, counts, recv)

    real(rk), intent(in) :: send(:)
    integer(ik), intent(in) :: counts(:)
    real(rk), allocatable, intent(out) :: recv(:)

    integer, allocatable :: counts_mpi(:)
    integer, allocatable :: displs(:)
    integer :: my_count
    integer :: ierror

    ! MPI must have been started; on the root the counts must describe send.
    call require_started('par_scatterv')
    call scatterv_layout(size(send), counts, counts_mpi, displs)

    ! Each rank's count, then the parts.
    call MPI_Scatter(counts_mpi, 1, MPI_INTEGER, my_count, 1, MPI_INTEGER, int(root_rank), &
      & MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_scatterv', 'MPI_Scatter')
    allocate (recv(my_count))
    call MPI_Scatterv(send, counts_mpi, displs, mpi_type_rk, recv, my_count, mpi_type_rk, &
      & int(root_rank), MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_scatterv', 'MPI_Scatterv')

  end subroutine par_scatterv_real

  !> Send each rank its part of an integer array held by root_rank (as par_scatterv_real).
  !>
  !> Arguments:
  !>   send    (in)   integer(ik)(:)               on root_rank: the parts, concatenated
  !>   counts  (in)   integer(ik)(:)               on root_rank: size of each rank's part
  !>   recv    (out)  integer(ik)(:), allocatable  this rank's part
  !> Locals:
  !>   counts_mpi  integer(:)   counts as MPI integers (root)
  !>   displs      integer(:)   offset of each part in send (root)
  !>   my_count    integer      size of this rank's part
  !>   ierror      integer      MPI error code
  subroutine par_scatterv_int(send, counts, recv)

    integer(ik), intent(in) :: send(:)
    integer(ik), intent(in) :: counts(:)
    integer(ik), allocatable, intent(out) :: recv(:)

    integer, allocatable :: counts_mpi(:)
    integer, allocatable :: displs(:)
    integer :: my_count
    integer :: ierror

    ! MPI must have been started; on the root the counts must describe send.
    call require_started('par_scatterv')
    call scatterv_layout(size(send), counts, counts_mpi, displs)

    ! Each rank's count, then the parts.
    call MPI_Scatter(counts_mpi, 1, MPI_INTEGER, my_count, 1, MPI_INTEGER, int(root_rank), &
      & MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_scatterv', 'MPI_Scatter')
    allocate (recv(my_count))
    call MPI_Scatterv(send, counts_mpi, displs, mpi_type_ik, recv, my_count, mpi_type_ik, &
      & int(root_rank), MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_scatterv', 'MPI_Scatterv')

  end subroutine par_scatterv_int

  !> Gather each rank's part of a real array on root_rank, concatenated in rank order
  !> (MPI_Gather of the counts, then MPI_Gatherv); the inverse of par_scatterv_real.
  !>
  !> Arguments:
  !>   send    (in)   real(rk)(:)               this rank's part (any size, including 0)
  !>   recv    (out)  real(rk)(:), allocatable  on root_rank: the parts of all ranks in rank
  !>                                            order; size 0 elsewhere
  !> Locals:
  !>   counts_mpi  integer(:)   size of each rank's part (root)
  !>   displs      integer(:)   offset of each part in recv (root)
  !>   ierror      integer      MPI error code
  subroutine par_gatherv_real(send, recv)

    real(rk), intent(in) :: send(:)
    real(rk), allocatable, intent(out) :: recv(:)

    integer, allocatable :: counts_mpi(:)
    integer, allocatable :: displs(:)
    integer :: ierror

    ! MPI must have been started; the root learns every rank's count.
    call require_started('par_gatherv')
    call gatherv_layout(size(send), counts_mpi, displs)

    ! The parts, on the root.
    allocate (recv(sum(int(counts_mpi, ik))))
    call MPI_Gatherv(send, size(send), mpi_type_rk, recv, counts_mpi, displs, mpi_type_rk, &
      & int(root_rank), MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_gatherv', 'MPI_Gatherv')

  end subroutine par_gatherv_real

  !> Gather each rank's part of an integer array on root_rank (as par_gatherv_real).
  !>
  !> Arguments:
  !>   send    (in)   integer(ik)(:)               this rank's part
  !>   recv    (out)  integer(ik)(:), allocatable  on root_rank: the parts in rank order; size 0
  !>                                               elsewhere
  !> Locals:
  !>   counts_mpi  integer(:)   size of each rank's part (root)
  !>   displs      integer(:)   offset of each part in recv (root)
  !>   ierror      integer      MPI error code
  subroutine par_gatherv_int(send, recv)

    integer(ik), intent(in) :: send(:)
    integer(ik), allocatable, intent(out) :: recv(:)

    integer, allocatable :: counts_mpi(:)
    integer, allocatable :: displs(:)
    integer :: ierror

    ! MPI must have been started; the root learns every rank's count.
    call require_started('par_gatherv')
    call gatherv_layout(size(send), counts_mpi, displs)

    ! The parts, on the root.
    allocate (recv(sum(int(counts_mpi, ik))))
    call MPI_Gatherv(send, size(send), mpi_type_ik, recv, counts_mpi, displs, mpi_type_ik, &
      & int(root_rank), MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_gatherv', 'MPI_Gatherv')

  end subroutine par_gatherv_int

  !> Counts and displacements of a gather: every rank's count, gathered on root_rank, and their
  !> running sums (zero-size placeholders elsewhere). A total beyond the default integer range
  !> (the MPI displacement type) is a fatal error.
  !>
  !> Arguments:
  !>   my_count    (in)   integer                  size of this rank's part
  !>   counts_mpi  (out)  integer(:), allocatable  size of each rank's part (root)
  !>   displs      (out)  integer(:), allocatable  offset of each part (root)
  !> Locals:
  !>   i       integer  rank index
  !>   ierror  integer  MPI error code
  subroutine gatherv_layout(my_count, counts_mpi, displs)

    integer, intent(in) :: my_count
    integer, allocatable, intent(out) :: counts_mpi(:)
    integer, allocatable, intent(out) :: displs(:)

    integer :: i
    integer :: ierror

    ! Every rank's count on the root.
    if (my_rank == root_rank) then
      allocate (counts_mpi(num_ranks), displs(num_ranks))
    else
      allocate (counts_mpi(0), displs(0))
    end if
    call MPI_Gather(my_count, 1, MPI_INTEGER, counts_mpi, 1, MPI_INTEGER, int(root_rank), &
      & MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_gatherv', 'MPI_Gather')

    ! Offsets (root), within the default integer range.
    if (my_rank == root_rank) then
      if (sum(int(counts_mpi, ik)) > int(huge(0), ik)) call fatal_error('par_gatherv', &
        & 'gathered size exceeds the largest default integer')
      displs(1) = 0
      do i = 2, int(num_ranks)
        displs(i) = displs(i - 1) + counts_mpi(i - 1)
      end do
    end if

  end subroutine gatherv_layout

  !> Counts and displacements of a scatter, on root_rank (checked against the send size);
  !> zero-size placeholders elsewhere.
  !>
  !> Arguments:
  !>   send_size   (in)   integer         size of the send array
  !>   counts      (in)   integer(ik)(:)  size of each rank's part (root)
  !>   counts_mpi  (out)  integer(:), allocatable  counts as MPI integers
  !>   displs      (out)  integer(:), allocatable  offset of each part
  !> Locals:
  !>   i_rank  integer  rank index (1-based)
  subroutine scatterv_layout(send_size, counts, counts_mpi, displs)

    integer, intent(in) :: send_size
    integer(ik), intent(in) :: counts(:)
    integer, allocatable, intent(out) :: counts_mpi(:)
    integer, allocatable, intent(out) :: displs(:)

    integer :: i_rank

    ! Only the root's arguments are significant.
    if (my_rank /= root_rank) then
      allocate (counts_mpi(1), displs(1))
      counts_mpi = 0
      displs = 0
      return
    end if

    ! One non-negative count per rank, summing to the send size.
    if (size(counts) /= num_ranks .or. any(counts < 0_ik)) then
      call internal_error('par_scatterv', 'expected one non-negative count per rank')
    end if
    if (sum(int(counts)) /= send_size) then
      call internal_error('par_scatterv', 'the counts do not add up to the send size')
    end if
    counts_mpi = int(counts)
    allocate (displs(size(counts)))
    displs(1) = 0
    do i_rank = 2, size(counts)
      displs(i_rank) = displs(i_rank - 1) + counts_mpi(i_rank - 1)
    end do

  end subroutine scatterv_layout

  !> Every rank sends one integer to every rank (MPI_Alltoall): recv(q + 1) on rank p is
  !> send(p + 1) of rank q.
  !>
  !> Arguments:
  !>   send  (in)   integer(ik)(num_ranks)  value for each rank (rank q at q + 1)
  !>   recv  (out)  integer(ik)(num_ranks)  value from each rank
  !> Locals:
  !>   ierror  integer  MPI error code
  subroutine par_alltoall(send, recv)

    integer(ik), intent(in) :: send(:)
    integer(ik), intent(out) :: recv(:)

    integer :: ierror

    ! MPI must have been started; one value per rank.
    call require_started('par_alltoall')
    if (size(send) /= num_ranks .or. size(recv) /= num_ranks) then
      call internal_error('par_alltoall', 'send and recv must have one entry per rank')
    end if

    ! Exchange.
    call MPI_Alltoall(send, 1, mpi_type_ik, recv, 1, mpi_type_ik, MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_alltoall', 'MPI_Alltoall')

  end subroutine par_alltoall

  !> Every rank sends a variable number of integers to every rank (MPI_Alltoallv). The parts are
  !> consecutive in send and in recv, in rank order.
  !>
  !> Arguments:
  !>   send         (in)   integer(ik)(:)               values for all ranks, concatenated
  !>   send_counts  (in)   integer(ik)(num_ranks)       number of values for each rank
  !>   recv         (out)  integer(ik)(:), allocatable  values from all ranks, concatenated
  !>                                                    (allocated to sum(recv_counts))
  !>   recv_counts  (in)   integer(ik)(num_ranks)       number of values from each rank (as
  !>                                                    sent; e.g. from par_alltoall)
  !> Locals:
  !>   scounts, rcounts  integer(:)  counts as MPI integers
  !>   sdispls, rdispls  integer(:)  offsets of the parts
  !>   i_rank            integer     rank index (1-based)
  !>   ierror            integer     MPI error code
  subroutine par_alltoallv(send, send_counts, recv, recv_counts)

    integer(ik), intent(in) :: send(:)
    integer(ik), intent(in) :: send_counts(:)
    integer(ik), allocatable, intent(out) :: recv(:)
    integer(ik), intent(in) :: recv_counts(:)

    integer, allocatable :: scounts(:)
    integer, allocatable :: rcounts(:)
    integer, allocatable :: sdispls(:)
    integer, allocatable :: rdispls(:)
    integer :: i_rank
    integer :: ierror

    ! MPI must have been started; counts per rank, consistent with send.
    call require_started('par_alltoallv')
    if (size(send_counts) /= num_ranks .or. size(recv_counts) /= num_ranks) then
      call internal_error('par_alltoallv', 'expected one send and one receive count per rank')
    end if
    if (any(send_counts < 0_ik) .or. any(recv_counts < 0_ik) .or. &
      & sum(int(send_counts)) /= size(send)) then
      call internal_error('par_alltoallv', 'invalid counts')
    end if

    ! Offsets of the parts.
    scounts = int(send_counts)
    rcounts = int(recv_counts)
    allocate (sdispls(num_ranks), rdispls(num_ranks))
    sdispls(1) = 0
    rdispls(1) = 0
    do i_rank = 2, int(num_ranks)
      sdispls(i_rank) = sdispls(i_rank - 1) + scounts(i_rank - 1)
      rdispls(i_rank) = rdispls(i_rank - 1) + rcounts(i_rank - 1)
    end do

    ! Exchange.
    allocate (recv(sum(rcounts)))
    call MPI_Alltoallv(send, scounts, sdispls, mpi_type_ik, recv, rcounts, rdispls, &
      & mpi_type_ik, MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_alltoallv', 'MPI_Alltoallv')

  end subroutine par_alltoallv

  !> Write the log lines every rank buffered (log_mod, all_ranks messages) on rank 0, in rank
  !> order, and empty the buffers (step 9.7, dev/plan/terminal_output.md §3.7). Collective:
  !> called at the collective points of the program (the end of the setup, of each increment,
  !> of the run) and before a collective fatal error that follows all_ranks warnings.
  !>
  !> Locals:
  !>   text        character(:), allocatable  this rank's buffered lines
  !>   all_text    character(:), allocatable  every rank's lines, in rank order (root)
  !>   counts_mpi  integer(:), allocatable    length of each rank's lines (root)
  !>   displs      integer(:), allocatable    offset of each rank's lines (root)
  !>   total       integer(ik)                total length, over all ranks
  !>   ierror      integer                    MPI error code
  subroutine par_log_flush()

    character(len=:), allocatable :: text, all_text
    integer, allocatable :: counts_mpi(:), displs(:)
    integer(ik) :: total
    integer :: ierror

    ! MPI must have been started; nothing to do when no rank has buffered lines.
    call require_started('par_log_flush')
    call log_buffer_take(text)
    call par_sum(int(len(text), ik), total)
    if (total == 0_ik) return

    ! Every rank's lines on the root, in rank order.
    call gatherv_layout(len(text), counts_mpi, displs)
    if (my_rank == root_rank) then
      allocate (character(len=sum(counts_mpi)) :: all_text)
    else
      allocate (character(len=0) :: all_text)
    end if
    call MPI_Gatherv(text, len(text), MPI_CHARACTER, all_text, counts_mpi, displs, &
      & MPI_CHARACTER, int(root_rank), MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, 'par_log_flush', 'MPI_Gatherv')

    ! Written by the root (log_lines_write writes on rank 0 only).
    call log_lines_write(all_text)

  end subroutine par_log_flush

  !> Abort handler registered with error_mod by par_init: ends every rank of the run with the
  !> given exit status (v2-fixed par_quit, error path, C6/C10). Does not return.
  !> It writes nothing itself (step 9.7, D94: the launchers report MPI_Abort; the former line
  !> "MPI_Abort: ending all ranks" is removed). The handler waits abort_drain_seconds before MPI_Abort: MPICH's Hydra launcher kills
  !> the job on MPI_Abort and can drop output that the processes had written (and flushed) but
  !> it had not yet forwarded, so the error message was sometimes lost (CI, step 4.3).
  !>
  !> Arguments:
  !>   exit_status  (in)  integer  exit status of the run (error_mod's error_exit_status)
  !> Locals:
  !>   ierror  integer   MPI error code (not checked: the run is ending)
  !>   start   real(rk)  MPI_Wtime at the start of the wait [s]
  subroutine par_abort_handler(exit_status)

    integer, intent(in) :: exit_status

    integer :: ierror
    real(rk) :: start

    ! Give the launcher time to forward the output written so far.
    start = MPI_Wtime()
    do while (MPI_Wtime() - start < abort_drain_seconds)
    end do

    ! End all ranks. MPI_Abort does not return; if it did, error_mod ends this process.
    call MPI_Abort(MPI_COMM_WORLD, exit_status, ierror)

  end subroutine par_abort_handler

  !> Real allreduce on MPI_COMM_WORLD with error checking (common part of the real reductions).
  !>
  !> Arguments:
  !>   part     (in)   real(rk)(..)  this rank's values (scalar or array)
  !>   whole    (out)  real(rk)(..)  reduced values
  !>   count    (in)   integer       number of values
  !>   op       (in)   type(MPI_Op)  reduction operation
  !>   routine  (in)   character     name of the public routine, for error messages
  !> Locals:
  !>   ierror   integer  MPI error code
  subroutine allreduce_real(part, whole, count, op, routine)

    real(rk), intent(in) :: part(..)
    real(rk), intent(out) :: whole(..)
    integer, intent(in) :: count
    type(MPI_Op), intent(in) :: op
    character(len=*), intent(in) :: routine

    integer :: ierror

    ! MPI must have been started.
    call require_started(routine)

    ! Reduce.
    call MPI_Allreduce(part, whole, count, mpi_type_rk, op, MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, routine, 'MPI_Allreduce')

  end subroutine allreduce_real

  !> Integer allreduce on MPI_COMM_WORLD with error checking (as allreduce_real).
  !>
  !> Arguments:
  !>   part     (in)   integer(ik)(..)  this rank's values (scalar or array)
  !>   whole    (out)  integer(ik)(..)  reduced values
  !>   count    (in)   integer          number of values
  !>   op       (in)   type(MPI_Op)     reduction operation
  !>   routine  (in)   character        name of the public routine, for error messages
  !> Locals:
  !>   ierror   integer  MPI error code
  subroutine allreduce_int(part, whole, count, op, routine)

    integer(ik), intent(in) :: part(..)
    integer(ik), intent(out) :: whole(..)
    integer, intent(in) :: count
    type(MPI_Op), intent(in) :: op
    character(len=*), intent(in) :: routine

    integer :: ierror

    ! MPI must have been started.
    call require_started(routine)

    ! Reduce.
    call MPI_Allreduce(part, whole, count, mpi_type_ik, op, MPI_COMM_WORLD, ierror)
    call check_mpi(ierror, routine, 'MPI_Allreduce')

  end subroutine allreduce_int

  !> Fatal error if MPI has not been started by par_init (or has been finalized).
  !>
  !> Arguments:
  !>   routine  (in)  character  name of the calling public routine
  subroutine require_started(routine)

    character(len=*), intent(in) :: routine

    ! Collective operations need MPI.
    if (.not. is_started) call internal_error(routine, 'MPI is not started (call par_init first)')

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

  !> Fatal error if an MPI call failed (public form of check_mpi, for the other modules of
  !> src/parallel/ that call MPI directly).
  !>
  !> Arguments:
  !>   ierror   (in)  integer    error code returned by the MPI call
  !>   routine  (in)  character  name of the calling routine
  !>   action   (in)  character  the MPI call
  subroutine par_check_mpi(ierror, routine, action)

    integer, intent(in) :: ierror
    character(len=*), intent(in) :: routine
    character(len=*), intent(in) :: action

    call check_mpi(ierror, routine, action)

  end subroutine par_check_mpi

  !> Fatal error if an MPI call failed, with MPI's description of the error code.
  !>
  !> Arguments:
  !>   ierror   (in)  integer    error code returned by the MPI call
  !>   routine  (in)  character  name of the calling routine
  !>   action   (in)  character  the MPI call, e.g. "MPI_Allreduce"
  !> Locals:
  !>   error_text  character(MPI_MAX_ERROR_STRING)  MPI's description of ierror
  !>   text_len    integer                          length of error_text
  !>   ierror2     integer                          error code of MPI_Error_string (not checked)
  !>   code_text   character(24)                    ierror as text
  subroutine check_mpi(ierror, routine, action)

    integer, intent(in) :: ierror
    character(len=*), intent(in) :: routine
    character(len=*), intent(in) :: action

    character(len=MPI_MAX_ERROR_STRING) :: error_text
    integer :: text_len
    integer :: ierror2
    character(len=24) :: code_text

    ! Nothing to do on success.
    if (ierror == MPI_SUCCESS) return

    ! Report the failure with MPI's error string.
    error_text = ''
    text_len = 0
    call MPI_Error_string(ierror, error_text, text_len, ierror2)
    write (code_text, '(i0)') ierror
    call fatal_error(routine, action // ' failed (MPI error ' // trim(code_text) // '): ' // &
      & error_text(1:max(text_len, 0)))

  end subroutine check_mpi

end module parallel_mod
