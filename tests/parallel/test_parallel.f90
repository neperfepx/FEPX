! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> MPI test of parallel_mod (and of log_mod under MPI): environment, collectives, broadcast.
!>
!> Source : FEPX (ROADMAP step 2.1). Replaces the step 0.1 smoke test, which used mpi_f08
!>          directly; this test uses only the src/parallel/ wrappers.
!>
!> Checks, on every rank, for any number of ranks n (run on 1, 2 and 4 ranks):
!>   - 0 <= my_rank < num_ranks, root_rank = 0; the ranks are distinct (sum of ranks n(n-1)/2,
!>     sum of squares n(n-1)(2n-1)/6);
!>   - par_sum (real and integer, scalars and arrays), par_max, par_min against closed forms;
!>     the real contributions are halves, whose sums are exact in binary for any summation order;
!>   - the launcher started one job of the expected number of ranks (environment variable
!>     FEPX_EXPECTED_RANKS, set by fepx_add_mpi_test), not that many singletons;
!>   - par_all / par_any on values true on all, some or no ranks;
!>   - par_bcast of each type: the root's value arrives on every rank, unchanged bit for bit;
!>   - par_gatherv of each type: parts of different sizes (rank 0's empty) arrive on the root
!>     concatenated in rank order, bit for bit; other ranks receive nothing;
!>   - par_wtime is non-decreasing across a barrier;
!>   - log_mod writes on root_rank only, and on every rank with all_ranks (each rank logs to its
!>     own scratch file, read back).
!>
!> Locals:
!>   n             integer(ik)        number of ranks (num_ranks)
!>   r             integer(ik)        rank of this process (my_rank)
!>   sum_int       integer(ik)        integer reduction result
!>   max_int       integer(ik)        integer maximum
!>   min_int       integer(ik)        integer minimum
!>   part_int      integer(ik)(3)     integer array contribution
!>   whole_int     integer(ik)(3)     integer array sum
!>   sum_real      real(rk)           real reduction result
!>   max_real      real(rk)           real maximum
!>   min_real      real(rk)           real minimum
!>   part_real     real(rk)(4)        real array contribution
!>   whole_real    real(rk)(4)        real array sum
!>   flag          logical            logical reduction result
!>   bcast_real    real(rk)           broadcast real scalar
!>   bcast_reals   real(rk)(3)        broadcast real array
!>   bcast_int     integer(ik)        broadcast integer scalar
!>   bcast_ints    integer(ik)(2)     broadcast integer array
!>   bcast_flag    logical            broadcast logical
!>   bcast_text    character(16)      broadcast string
!>   time_before   real(rk)           par_wtime before the barrier [s]
!>   time_after    real(rk)           par_wtime after the barrier [s]
!>   unit_log      integer            this rank's scratch log unit
!>   ios           integer            iostat of the log reads
!>   line          character(128)     one log line read back
!>   rank_text     character(12)      r as text
!>   nf            real(rk)           n as a real
!>   gath_real     real(rk)(:)        par_gatherv result (real)
!>   gath_int      integer(ik)(:)     par_gatherv result (integer)
!>   expect_real   real(rk)(:)        expected gathered reals (root)
!>   expect_int    integer(ik)(:)     expected gathered integers (root)
!>   i, k          integer(ik)        loop indices
!>   expected      character(8)       FEPX_EXPECTED_RANKS
!>   env_len, env_stat  integer       its length and retrieval status
program test_parallel

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use parallel_mod, only: par_init, par_finalize, par_barrier, par_wtime, par_sum, par_max, &
    & par_min, par_all, par_any, par_bcast, par_gatherv, num_ranks, my_rank, root_rank, &
    & par_log_flush
  use log_mod, only: log_info, log_warning, log_set_unit
  use assert_mod, only: assert_equal, assert_close, assert_true, test_summary

  implicit none

  integer(ik) :: n
  integer(ik) :: r
  integer(ik) :: sum_int
  integer(ik) :: max_int
  integer(ik) :: min_int
  integer(ik) :: part_int(3)
  integer(ik) :: whole_int(3)
  real(rk) :: sum_real
  real(rk) :: max_real
  real(rk) :: min_real
  real(rk) :: part_real(4)
  real(rk) :: whole_real(4)
  logical :: flag
  real(rk) :: bcast_real
  real(rk) :: bcast_reals(3)
  integer(ik) :: bcast_int
  integer(ik) :: bcast_ints(2)
  logical :: bcast_flag
  character(len=16) :: bcast_text
  real(rk), allocatable :: gath_real(:), expect_real(:)
  integer(ik), allocatable :: gath_int(:), expect_int(:)
  integer(ik) :: i, k
  character(len=8) :: expected
  integer :: env_len, env_stat
  real(rk) :: time_before
  real(rk) :: time_after
  integer :: unit_log
  integer :: ios
  character(len=128) :: line
  character(len=12) :: rank_text
  real(rk) :: nf

  ! Start MPI through the wrappers.
  call par_init()
  n = num_ranks
  r = my_rank

  ! --- Launch: one job of the expected number of ranks ---------------------------------------

  call get_environment_variable('FEPX_EXPECTED_RANKS', expected, env_len, env_stat)
  call assert_equal(env_stat, 0, 'FEPX_EXPECTED_RANKS is set')
  write (rank_text, '(i0)') n
  call assert_equal(trim(rank_text), trim(expected), 'number of ranks as launched')
  nf = real(n, rk)

  ! --- Environment ---------------------------------------------------------------------------

  call assert_true(n >= 1_ik, 'num_ranks >= 1')
  call assert_true(r >= 0_ik .and. r < n, '0 <= my_rank < num_ranks')
  call assert_equal(root_rank, 0_ik, 'root_rank')

  ! The ranks are 0, ..., n - 1, each once.
  call par_sum(r, sum_int)
  call assert_equal(sum_int, n*(n - 1_ik)/2_ik, 'sum of ranks')
  call par_sum(r*r, sum_int)
  call assert_equal(sum_int, n*(n - 1_ik)*(2_ik*n - 1_ik)/6_ik, 'sum of squared ranks')

  ! --- Integer reductions --------------------------------------------------------------------

  call par_sum(1_ik, sum_int)
  call assert_equal(sum_int, n, 'par_sum integer scalar')
  call par_max(r, max_int)
  call assert_equal(max_int, n - 1_ik, 'par_max integer')
  call par_min(-r, min_int)
  call assert_equal(min_int, 1_ik - n, 'par_min integer')
  part_int = [1_ik, r, -2_ik*r]
  call par_sum(part_int, whole_int)
  call assert_equal(whole_int(1), n, 'par_sum integer array 1')
  call assert_equal(whole_int(2), n*(n - 1_ik)/2_ik, 'par_sum integer array 2')
  call assert_equal(whole_int(3), -n*(n - 1_ik), 'par_sum integer array 3')

  ! --- Real reductions (halves: exact sums) --------------------------------------------------

  call par_sum(0.5_rk*real(r + 1_ik, rk), sum_real)
  call assert_close(sum_real, nf*(nf + 1.0_rk)/4.0_rk, 'par_sum real scalar')
  call par_max(-0.5_rk*real(r, rk), max_real)
  call assert_close(max_real, 0.0_rk, 'par_max real')
  call par_min(-0.5_rk*real(r, rk), min_real)
  call assert_close(min_real, -0.5_rk*(nf - 1.0_rk), 'par_min real')
  part_real = [0.5_rk, -1.5_rk*real(r, rk), 0.0_rk, 2.0_rk**(-20)]
  call par_sum(part_real, whole_real)
  call assert_close(whole_real, [0.5_rk*nf, -0.75_rk*nf*(nf - 1.0_rk), 0.0_rk, &
    & nf*2.0_rk**(-20)], 'par_sum real array')

  ! --- Logical reductions --------------------------------------------------------------------

  call par_all(.true., flag)
  call assert_true(flag, 'par_all true everywhere')
  call par_all(r /= n - 1_ik, flag)
  call assert_true(.not. flag, 'par_all false on the last rank')
  call par_any(.false., flag)
  call assert_true(.not. flag, 'par_any false everywhere')
  call par_any(r == n - 1_ik, flag)
  call assert_true(flag, 'par_any true on the last rank only')
  call par_all(r == 0_ik, flag)
  call assert_true(flag .eqv. (n == 1_ik), 'par_all true on rank 0 only')

  ! --- Broadcast from the root ---------------------------------------------------------------

  ! Non-root ranks start from different values, which must be overwritten.
  if (r == root_rank) then
    bcast_real = 1.0_rk/3.0_rk
    bcast_reals = [-1.0e300_rk, tiny(1.0_rk), huge(1.0_rk)]
    bcast_int = -huge(1_ik)
    bcast_ints = [7_ik, huge(1_ik)]
    bcast_flag = .true.
    bcast_text = 'euler-bunge'
  else
    bcast_real = -1.0_rk
    bcast_reals = 0.0_rk
    bcast_int = r
    bcast_ints = -r
    bcast_flag = .false.
    bcast_text = 'not the root'
  end if
  call par_bcast(bcast_real)
  call par_bcast(bcast_reals)
  call par_bcast(bcast_int)
  call par_bcast(bcast_ints)
  call par_bcast(bcast_flag)
  call par_bcast(bcast_text)
  call assert_close(bcast_real, 1.0_rk/3.0_rk, 'par_bcast real', abs_tol=0.0_rk, rel_tol=0.0_rk)
  call assert_close(bcast_reals, [-1.0e300_rk, tiny(1.0_rk), huge(1.0_rk)], &
    & 'par_bcast real array', abs_tol=0.0_rk, rel_tol=0.0_rk)
  call assert_equal(bcast_int, -huge(1_ik), 'par_bcast integer')
  call assert_equal(bcast_ints(1), 7_ik, 'par_bcast integer array 1')
  call assert_equal(bcast_ints(2), huge(1_ik), 'par_bcast integer array 2')
  call assert_true(bcast_flag, 'par_bcast logical')
  call assert_equal(bcast_text, 'euler-bunge     ', 'par_bcast character')

  ! --- Gather: rank k sends k values 1000 k + i + 1/3 (rank 0 sends none) -------------------

  call par_gatherv([(1000.0_rk*r + i + 1.0_rk/3.0_rk, i = 1_ik, r)], gath_real)
  call par_gatherv([(1000_ik*r + i, i = 1_ik, r)], gath_int)
  allocate (expect_real(0), expect_int(0))
  if (r == root_rank) then
    expect_real = [((1000.0_rk*k + i + 1.0_rk/3.0_rk, i = 1_ik, k), k = 0_ik, n - 1_ik)]
    expect_int = [((1000_ik*k + i, i = 1_ik, k), k = 0_ik, n - 1_ik)]
  end if
  call assert_equal(size(gath_real), size(expect_real), 'par_gatherv real size')
  call assert_equal(size(gath_int), size(expect_int), 'par_gatherv integer size')
  if (size(gath_real) == size(expect_real)) call assert_close(gath_real, expect_real, &
    & 'par_gatherv real values', abs_tol=0.0_rk, rel_tol=0.0_rk)
  if (size(gath_int) == size(expect_int)) call assert_true(all(gath_int == expect_int), &
    & 'par_gatherv integer values')

  ! --- Wall clock ----------------------------------------------------------------------------

  time_before = par_wtime()
  call par_barrier()
  time_after = par_wtime()
  call assert_true(time_after >= time_before, 'par_wtime non-decreasing')

  ! --- Logging under MPI: rank 0 only, unless all_ranks --------------------------------------

  ! Each rank logs to its own scratch file.
  open (newunit=unit_log, status='scratch', action='readwrite', form='formatted', iostat=ios)
  call assert_equal(ios, 0, 'open scratch log')
  call log_set_unit(unit_log)
  call log_info('root only')
  call log_warning('root only')
  call log_info('every rank', all_ranks=.true.)
  call par_log_flush()
  call log_set_unit(output_unit)

  ! Read back: rank 0 has the two plain lines, then every rank's tagged line in rank order
  ! (buffered, written by rank 0 at par_log_flush, step 9.7); the other ranks have nothing.
  rewind (unit_log)
  if (r == root_rank) then
    read (unit_log, '(a)', iostat=ios) line
    call assert_equal(trim(line), 'Info   : root only', 'root info line')
    read (unit_log, '(a)', iostat=ios) line
    call assert_equal(trim(line), 'Warning: root only', 'root warning line')
    do i = 0_ik, num_ranks - 1_ik
      write (rank_text, '(i0)') i
      read (unit_log, '(a)', iostat=ios) line
      call assert_equal(ios, 0, 'read tagged line')
      call assert_equal(trim(line), 'Info   : [Rank ' // trim(rank_text) // '] every rank', &
        & 'tagged line')
    end do
  end if
  read (unit_log, '(a)', iostat=ios) line
  call assert_true(ios /= 0, 'no extra log line')
  close (unit_log)

  ! Report (summary line on the root only; a failing rank stops the job), then end MPI.
  call test_summary(verbose=(r == root_rank))
  call par_finalize()

end program test_parallel
