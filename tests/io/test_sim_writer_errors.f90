! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Error cases of sim_writer_mod: each case must end every rank with a fatal error and its
!> message (checked by tests/core/check_run.cmake through mpiexec).
!>
!> Source : FEPX (ROADMAP step 4.3).
!>
!> Usage: mpiexec -n N test_sim_writer_errors <case> <io fixture dir> (tests/io/fixtures)
!> Cases (on tiny.msh; output directory err_<case>_np<N>.sim, removed first):
!>   dir_exists     the output directory exists already
!>   unknown_var    print of a variable the writer does not write yet (strain_th)
!>   missing_input  an input file that does not exist
!>   node_shape     node values with 2 rows instead of 3
!>   elt_shape      element values with 3 rows instead of 6
!>   forces_shape   one force column fewer than the face sets
!>   slip_shape     slip-system counts for one element fewer than the local elements
!>   conv_shape     convergence record with five norms per iteration instead of six
!>
!> Locals:
!>   case_name, fixtures  character(512)      case and fixture directory
!>   mesh                 type(mesh_t)        rank-local mesh
!>   sim                  type(sim_writer_t)  the writer
!>   dir                  character(:), allocatable  output directory
!>   inputs               character(512)(1)   input files
!>   num_text             character(8)        number of ranks as text
!>   unit_log             integer             scratch unit receiving the log
!>   n_fasets             integer             number of face sets
program test_sim_writer_errors

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_init, par_finalize, par_barrier, num_ranks, my_rank, root_rank
  use mesh_mod, only: mesh_t
  use mesh_blocks_mod, only: mesh_load_blocks
  use sim_writer_mod, only: sim_writer_t, sim_open, sim_write_nodes, sim_write_elts, &
    & sim_write_forces, sim_write_elts_slip, sim_write_convergence

  implicit none

  character(len=512) :: case_name, fixtures
  type(mesh_t) :: mesh
  type(sim_writer_t) :: sim
  character(len=:), allocatable :: dir
  character(len=512) :: inputs(1)
  character(len=8) :: num_text
  integer :: unit_log, n_fasets

  call par_init()
  call get_command_argument(1, case_name)
  call get_command_argument(2, fixtures)
  open (newunit=unit_log, status='scratch', action='readwrite', form='formatted')
  call log_set_unit(unit_log)
  call mesh_load_blocks(mesh, trim(fixtures) // '/tiny.msh')
  call log_set_unit(output_unit)
  close (unit_log)

  ! A fresh output directory for this case and rank count.
  write (num_text, '(i0)') num_ranks
  dir = 'err_' // trim(case_name) // '_np' // trim(num_text) // '.sim'
  if (my_rank == root_rank) call execute_command_line('rm -rf ' // dir)
  call par_barrier()
  inputs(1) = trim(fixtures) // '/tiny.msh'
  n_fasets = size(mesh%fasets)

  select case (trim(case_name))
  case ('dir_exists')
    if (my_rank == root_rank) call execute_command_line('mkdir ' // dir)
    call par_barrier()
    call sim_open(sim, dir, [character(len=32) :: 'coo'], mesh, inputs)
  case ('unknown_var')
    call sim_open(sim, dir, [character(len=32) :: 'coo', 'strain_th'], mesh, inputs)
  case ('missing_input')
    inputs(1) = 'no_such_file.cfg'
    call sim_open(sim, dir, [character(len=32) :: 'coo'], mesh, inputs)
  case ('node_shape')
    call sim_open(sim, dir, [character(len=32) :: 'coo'], mesh, inputs)
    call sim_write_nodes(sim, mesh, 'coo', mesh%coo_s(1:2, :), 0)
  case ('elt_shape')
    call sim_open(sim, dir, [character(len=32) :: 'stress'], mesh, inputs)
    call sim_write_elts(sim, mesh, 'stress', mesh%coo_s(:, 1:mesh%num_elts), 0)
  case ('forces_shape')
    call sim_open(sim, dir, [character(len=32) :: 'forces'], mesh, inputs)
    call sim_write_forces(sim, 0, 0, spread([0.0_rk, 0.0_rk, 0.0_rk], 2, n_fasets - 1), &
      & spread(1.0_rk, 1, n_fasets), 0.0_rk)
  case ('slip_shape')
    call sim_open(sim, dir, [character(len=32) :: 'crss'], mesh, inputs)
    call sim_write_elts_slip(sim, mesh, 'crss', mesh%coo_s(:, 1:mesh%num_elts), &
      & spread(3_ik, 1, mesh%num_elts - 1_ik), 0)
  case ('conv_shape')
    call sim_open(sim, dir, [character(len=32) :: 'convergence'], mesh, inputs)
    call sim_write_convergence(sim, 1, [character(len=2) :: 'SA'], &
      & spread(spread(0.0_rk, 1, 5), 2, 1), [3])
  case default
    write (*, '(a)') 'unknown case ' // trim(case_name)
    stop 2
  end select
  write (*, '(a)') 'no error raised'
  call par_finalize()

end program test_sim_writer_errors
