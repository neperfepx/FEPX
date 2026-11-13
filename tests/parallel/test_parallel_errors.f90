! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Error and abort cases of parallel_mod: runs one case, named by the first command-line argument.
!>
!> Source : FEPX (ROADMAP step 2.1). Each case ends the run in a known way; the exit status and
!>          the message on standard output and standard error are checked by
!>          tests/core/check_run.cmake (tests/CMakeLists.txt), with a CTest timeout so that a hang
!>          fails. The behaviour required is v2-fixed's after CHANGES.md C6/C10: an error detected
!>          on one rank only ends every rank with a non-zero exit status.
!>
!> Cases (run through the MPI launcher unless stated):
!>   one_rank       the last rank calls fatal_error while the others wait in a barrier
!>   all_ranks      every rank calls fatal_error
!>   collective     the last rank calls fatal_error while the others wait in a reduction
!>   size_mismatch  par_sum with a result array of the wrong size (every rank)
!>   double_init    par_init called twice
!>   not_started    par_sum before par_init (run without the launcher: MPI is never started)
!>   after_finalize fatal_error after par_finalize: no handler, ends with `error stop`
!>   halo_size      halo_forward with a vector of the wrong size (halo_mod, step 2.4)
!>   dist_dot_size  dist_dot with more owned entries than the vector has
!>
!> Locals:
!>   case_name   character(32)   name of the case
!>   sum_real    real(rk)        reduction result (unused: the cases end the run)
!>   part        real(rk)(3)     array contribution
!>   whole       real(rk)(2)     result array of the wrong size
!>   mesh        mesh_t          a 2-node mesh without neighbors (halo cases)
!>   x5          real(rk)(5)     a vector of 5 entries
!>   dot         real(rk)        dot product (unused: the case ends the run)
program test_parallel_errors

  use kinds_mod, only: rk
  use error_mod, only: fatal_error
  use parallel_mod, only: par_init, par_finalize, par_barrier, par_sum, num_ranks, my_rank
  use mesh_mod, only: mesh_t
  use halo_mod, only: halo_forward, dist_dot

  implicit none

  character(len=32) :: case_name
  real(rk) :: sum_real
  real(rk) :: part(3)
  real(rk) :: whole(2)
  type(mesh_t) :: mesh
  real(rk) :: x5(5)
  real(rk) :: dot

  ! Read the case name.
  call get_command_argument(1, case_name)

  select case (trim(case_name))

  case ('one_rank')
    ! Error on the last rank only; the others must not hang in the barrier.
    call par_init()
    if (my_rank == num_ranks - 1) call fatal_error('test_routine', 'error on one rank')
    call par_barrier()

  case ('all_ranks')
    ! Error on every rank.
    call par_init()
    call fatal_error('test_routine', 'error on every rank')

  case ('collective')
    ! Error on the last rank only while the others wait in a reduction.
    call par_init()
    if (my_rank == num_ranks - 1) call fatal_error('test_routine', 'error during a reduction')
    call par_sum(1.0_rk, sum_real)

  case ('size_mismatch')
    ! Result array of the wrong size.
    call par_init()
    part = 1.0_rk
    call par_sum(part, whole)

  case ('double_init')
    ! MPI may be started only once.
    call par_init()
    call par_init()

  case ('not_started')
    ! Collective before par_init.
    call par_sum(1.0_rk, sum_real)

  case ('after_finalize')
    ! After par_finalize the handler is removed: error stop with error_mod's status.
    call par_init()
    call par_finalize()
    call fatal_error('test_routine', 'error after finalize')

  case ('halo_size')
    ! 5 entries for 3 values per node on 2 local nodes.
    call par_init()
    mesh%num_nodes_loc = 2
    x5 = 0.0_rk
    call halo_forward(mesh, x5, 3)

  case ('dist_dot_size')
    ! 6 owned entries claimed for vectors of 5.
    call par_init()
    x5 = 1.0_rk
    dot = dist_dot(x5, x5, 6)
    write (*, '(es12.4)') dot

  case default
    call fatal_error('test_parallel_errors', 'unknown case ' // trim(case_name))

  end select

  ! Reaching this point is a test failure (every case must end the run with an error).
  write (*, '(a)') 'case returned normally'
  error stop 99

end program test_parallel_errors
