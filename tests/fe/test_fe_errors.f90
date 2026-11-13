! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Error cases of bc_mod and constraint_mod: each case must end every rank with a fatal error and
!> its message
!> (checked by tests/core/check_run.cmake through mpiexec).
!>
!> Source : FEPX (ROADMAP step 4.2; MPC cases step 8.3a).
!>
!> Usage: mpiexec -n N test_fe_errors <case> <io fixture dir> (tests/io/fixtures)
!> Cases (on tiny.msh, node sets "face" and "apex"):
!>   unknown_nset    a condition on node set "x9", which the mesh does not have
!>   conflict        node set "face" given 0 and then 1 in component x
!>   bad_axis        a preset along axis 4
!>   bad_component   a condition with component 0
!>   mpc_unknown_nset  an MPC on node set "x9"
!>   mpc_two_lines     node 3 a secondary of node 2 (set A = 2 3) and of node 4 (set B = 4 3)
!>   mpc_cycle         node 3 a secondary of 2 (A = 2 3) and 2 a secondary of 3 (B = 3 2)
!>   mpc_conflict      set_mpc1 x on "face" (2 3 4) with x conditions 1 on node 3 and 2 on
!>                     node 4 (C66: the group cannot satisfy both)
!>
!> Locals:
!>   case_name, fixtures  character(512)        case and fixture directory
!>   mesh                 type(mesh_t)          rank-local mesh
!>   specs                type(bc_spec_t)(:)    conditions
!>   bc_mask              logical(:)            prescribed dofs
!>   bc_vel               real(rk)(:)           prescribed values
!>   unit_log             integer               scratch unit receiving the log
!>   mpcs                 type(mpc_spec_t)(:)   MPC relations
!>   cons                 type(constraint_t)    constraints
program test_fe_errors

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_init, par_finalize
  use mesh_mod, only: mesh_t
  use mesh_blocks_mod, only: mesh_load_blocks
  use bc_mod, only: bc_spec_t, bc_preset_specs, bc_resolve, bc_preset_minimal
  use mesh_file_mod, only: mesh_nset_t
  use constraint_mod, only: mpc_spec_t, constraint_t, constraint_build_mpc, &
    & constraint_propagate_bcs

  implicit none

  character(len=512) :: case_name, fixtures
  type(mesh_t) :: mesh
  type(bc_spec_t), allocatable :: specs(:)
  logical, allocatable :: bc_mask(:)
  real(rk), allocatable :: bc_vel(:)
  integer :: unit_log
  type(mpc_spec_t), allocatable :: mpcs(:)
  type(constraint_t) :: cons

  call par_init()
  call get_command_argument(1, case_name)
  call get_command_argument(2, fixtures)
  open (newunit=unit_log, status='scratch', action='readwrite', form='formatted')
  call log_set_unit(unit_log)
  call mesh_load_blocks(mesh, trim(fixtures) // '/tiny.msh')
  call log_set_unit(output_unit)
  close (unit_log)
  allocate (bc_mask(3*mesh%num_nodes_loc), bc_vel(3*mesh%num_nodes_loc))

  select case (trim(case_name))
  case ('unknown_nset')
    allocate (specs(1))
    specs(1)%nset = 'x9'
    specs(1)%comp = 1
    call bc_resolve(mesh, specs, bc_mask, bc_vel)
  case ('conflict')
    allocate (specs(2))
    specs(1)%nset = 'face'
    specs(1)%comp = 1
    specs(1)%value = 0.0_rk
    specs(2)%nset = 'face'
    specs(2)%comp = 1
    specs(2)%value = 1.0_rk
    call bc_resolve(mesh, specs, bc_mask, bc_vel)
  case ('bad_axis')
    call bc_preset_specs(bc_preset_minimal, 4, 1.0_rk, specs)
  case ('bad_component')
    allocate (specs(1))
    specs(1)%nset = 'apex'
    specs(1)%comp = 0
    call bc_resolve(mesh, specs, bc_mask, bc_vel)
  case ('mpc_unknown_nset')
    mpcs = [mpc_spec_t('x9', 1, 1.0_rk, 0.0_rk, 1_ik)]
    call constraint_build_mpc(cons, mesh, mpcs)
  case ('mpc_two_lines')
    call add_sets([2_ik, 3_ik], [4_ik, 3_ik])
    mpcs = [mpc_spec_t('A', 1, 1.0_rk, 0.0_rk, 1_ik), mpc_spec_t('B', 1, 1.0_rk, 0.0_rk, 2_ik)]
    call constraint_build_mpc(cons, mesh, mpcs)
  case ('mpc_cycle')
    call add_sets([2_ik, 3_ik], [3_ik, 2_ik])
    mpcs = [mpc_spec_t('A', 1, 1.0_rk, 0.0_rk, 1_ik), mpc_spec_t('B', 1, 1.0_rk, 0.0_rk, 2_ik)]
    call constraint_build_mpc(cons, mesh, mpcs)
  case ('mpc_conflict')
    mpcs = [mpc_spec_t('face', 1, 1.0_rk, 0.0_rk, 1_ik)]
    call constraint_build_mpc(cons, mesh, mpcs)
    allocate (specs(2))
    specs(1)%nset = '3'
    specs(1)%comp = 1
    specs(1)%value = 1.0_rk
    specs(2)%nset = '4'
    specs(2)%comp = 1
    specs(2)%value = 2.0_rk
    call bc_resolve(mesh, specs, bc_mask, bc_vel)
    call constraint_propagate_bcs(cons, mesh, bc_mask, bc_vel)
  case default
    write (*, '(a)') 'unknown case ' // trim(case_name)
    stop 2
  end select
  write (*, '(a)') 'no error raised'
  call par_finalize()

contains

  !> Two node sets "A" and "B" appended to the whole-mesh node sets (MPC cases).
  !>
  !> Arguments:
  !>   nodes_a, nodes_b  (in)  integer(ik)(:)  global node ids, in order
  !> Locals:
  !>   grown  type(mesh_nset_t)(:)  the node sets with A and B
  !>   i      integer               node-set index
  subroutine add_sets(nodes_a, nodes_b)

    integer(ik), intent(in) :: nodes_a(:), nodes_b(:)

    type(mesh_nset_t), allocatable :: grown(:)
    integer :: i

    allocate (grown(size(mesh%nsets_glob) + 2))
    do i = 1, size(mesh%nsets_glob)
      grown(i) = mesh%nsets_glob(i)
    end do
    grown(size(grown) - 1)%label = 'A'
    grown(size(grown) - 1)%nodes = nodes_a
    grown(size(grown))%label = 'B'
    grown(size(grown))%nodes = nodes_b
    call move_alloc(grown, mesh%nsets_glob)

  end subroutine add_sets

end program test_fe_errors
