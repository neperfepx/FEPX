! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> MPI tests of the velocity constraints v = T v_hat + g (constraint_mod, constraint_exchange_mod;
!> ROADMAP step 8.3a; dev/theory/addendum_periodic_mpc.md §3.1, §3.2, §5), on 1, 2 and 4 ranks.
!>
!> Source : FEPX (ROADMAP step 8.3a).
!>
!> Usage: mpiexec -n N test_constraint <regression mesh dir> (tests/regression/cases/meshes).
!>
!> Mesh: poly_part1.msh (unit cube, 656 tet10; node sets of Neper's faces, edges and corners),
!> distributed with partitioning "even", so that MPC groups cross the partitions; material:
!> isotropic [s] = [C] (E = 200 GPa, nu = 0.3) at every point. Checks:
!>   - the table (constraint_table_mpc, synthetic node sets): rows of a chain (coefficients
!>     multiplied, offsets composed), a repeated identical relation kept once, a node set with
!>     fewer than two nodes giving no row;
!>   - the table on the mesh: one row per secondary dof, the first node of each node set its
!>     primary; local rows and columns consistent on every rank;
!>   - T and T^T adjoint: (T x) . y = x . (T^T y) over the owned dofs, x reduced, y arbitrary
!>     (relative 1e-12);
!>   - the reduced operator (constraint_matvec) equal to T^T K T x of a serial assembly of the
!>     whole mesh with the explicit table (relative 1e-12);
!>   - the exact diagonal (constraint_diagonal) equal to (T e_q)^T K (T e_q) of the serial
!>     assembly, for every primary dof q and a sample of the other independent dofs
!>     (relative 1e-12), and zero at the constrained dofs;
!>   - velocity conditions (constraint_propagate_bcs, C66): a condition on a secondary node
!>     prescribes the whole group, values given on two nodes of a group averaged when equal;
!>   - solves (V4): uniaxial stress along z (symmetry conditions) with set_mpc1 x on x1 (the
!>     constraint inactive): the exact homogeneous field (1e-10 e, cg_tol 1e-12); with a linear
!>     MPC (c = -1 on y1, y) and set_mpc1 x on x1: the relations hold to round-off (1e-14 e;
!>     v = v0 + dv adds the correction on both sides separately, as v2) and the reduced
!>     residual T^T (K v) vanishes at the free independent dofs (1e-8 of the force scale).
!>
!> Locals:
!>   ework  type(element_work_t)  element work arrays (element_work_alloc; D100)
!>   young, nu, rate  real(rk)        Young's modulus [MPa], Poisson's ratio, strain rate [1/s]
!>   meshes     character(512)        regression mesh directory
!>   mesh       type(mesh_t)          rank-local mesh
!>   ref        type(mesh_file_t)     whole mesh (serial references)
!>   elt        type(element_t)       element descriptor
!>   mat_s      real(rk)(:,:,:)       [s] at the quadrature points of an element [MPa]
!>   k_elt      real(rk)(:,:,:)       local element matrices
!>   k_ref_all  real(rk)(:,:,:)       element matrices of the whole mesh (serial references)
!>   i_elt      integer(ik)           element
!>   unit_log   integer               scratch unit receiving the log
program test_constraint

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_init, par_finalize, par_sum, par_max
  use mesh_file_mod, only: mesh_file_t, mesh_file_load, mesh_nset_t
  use mesh_mod, only: mesh_t, find_sorted
  use mesh_blocks_mod, only: mesh_load_blocks
  use element_mod, only: element_t, element_init
  use element_matrix_mod, only: element_stiffness
  use halo_mod, only: halo_forward, dist_dot
  use ebe_mod, only: ebe_matvec
  use pcg_mod, only: pcg_solve
  use bc_mod, only: bc_spec_t, bc_preset_specs, bc_resolve, bc_preset_symmetry
  use constraint_mod, only: mpc_spec_t, constraint_t, constraint_table_mpc, &
    & constraint_build_mpc, constraint_expand, constraint_reduce, constraint_matvec, &
    & constraint_diagonal, constraint_propagate_bcs, constraint_project
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary
  use element_matrix_mod, only: element_work_t, element_work_alloc

  implicit none

  type(element_work_t) :: ework
  real(rk), parameter :: young = 200.0e3_rk
  real(rk), parameter :: nu = 0.3_rk
  real(rk), parameter :: rate = 1.0e-3_rk

  character(len=512) :: meshes
  type(mesh_t) :: mesh
  type(mesh_file_t) :: ref
  type(element_t) :: elt
  real(rk), allocatable :: mat_s(:, :, :), k_elt(:, :, :), k_ref_all(:, :, :)
  integer :: unit_log
  integer(ik) :: i_elt

  call par_init()
  call get_command_argument(1, meshes)

  ! Mesh (log to a scratch unit) and the whole mesh for the serial references.
  open (newunit=unit_log, status='scratch', action='readwrite', form='formatted')
  call log_set_unit(unit_log)
  call mesh_load_blocks(mesh, trim(meshes) // '/poly_part1.msh')
  call mesh_file_load(ref, trim(meshes) // '/poly_part1.msh')

  ! Isotropic [s] at every point and the element matrices.
  call element_init(elt, mesh%elt_type)
  allocate (mat_s(6, 6, elt%num_qpts))
  mat_s = spread(isotropic_mandel(young, nu), 3, int(elt%num_qpts))
  allocate (k_elt(elt%num_dofs, elt%num_dofs, mesh%num_elts))
  do i_elt = 1_ik, mesh%num_elts
    call element_work_alloc(ework, elt)
    call element_stiffness(elt, mesh%coo_s(:, mesh%elt_nodes_loc(:, i_elt)), mat_s, &
      & k_elt(:, :, i_elt), ework)
  end do
  allocate (k_ref_all(elt%num_dofs, elt%num_dofs, ref%num_elts))
  do i_elt = 1_ik, ref%num_elts
    call element_work_alloc(ework, elt)
    call element_stiffness(elt, ref%coo(:, ref%elt_nodes(:, i_elt)), mat_s, &
      & k_ref_all(:, :, i_elt), ework)
  end do

  call check_table_synthetic()
  call check_operators()
  call check_propagation()
  call check_solves()
  call log_set_unit(output_unit)
  close (unit_log)

  call test_summary()
  call par_finalize()

contains

  !> The table of synthetic node sets: chain A = (10, 20, 30) x with c = 2, d = 1 and
  !> B = (40, 10) x with c = 3, d = 0.5 gives v10 = 3 v40 + 0.5 and v20 = v30 = 6 v40 + 2; a
  !> repeated identical relation (C = (10, 20) x, c = 2, d = 1) is kept once; a one-node set
  !> gives no row.
  !>
  !> Locals:
  !>   sets   type(mesh_nset_t)(4)    synthetic node sets
  !>   specs  type(mpc_spec_t)(4)     relations
  !>   cons   type(constraint_t)      the table
  subroutine check_table_synthetic()

    type(mesh_nset_t) :: sets(4)
    type(mpc_spec_t) :: specs(4)
    type(constraint_t) :: cons

    sets(1)%label = 'A'
    sets(1)%nodes = [10_ik, 20_ik, 30_ik]
    sets(2)%label = 'B'
    sets(2)%nodes = [40_ik, 10_ik]
    sets(3)%label = 'C'
    sets(3)%nodes = [10_ik, 20_ik]
    sets(4)%label = 'D'
    sets(4)%nodes = [50_ik]
    specs(1) = mpc_spec_t('A', 1, 2.0_rk, 1.0_rk, 1_ik)
    specs(2) = mpc_spec_t('B', 1, 3.0_rk, 0.5_rk, 2_ik)
    specs(3) = mpc_spec_t('C', 1, 2.0_rk, 1.0_rk, 3_ik)
    specs(4) = mpc_spec_t('D', 2, 1.0_rk, 0.0_rk, 4_ik)
    call constraint_table_mpc(specs, sets, cons)

    call assert_equal(int(cons%num_rows_glob), 3, 'table: three rows')
    call assert_true(all(cons%glob_dof == [28_ik, 58_ik, 88_ik]), 'table: rows of nodes 10, ' // &
      & '20, 30 (x), increasing')
    call assert_true(all(cons%glob_col == 118_ik), 'table: every row on node 40 (x)')
    call assert_close(cons%glob_coef, [3.0_rk, 6.0_rk, 6.0_rk], 'table: composed coefficients', &
      & abs_tol=1e-15_rk)
    call assert_close(cons%glob_offset, [0.5_rk, 2.0_rk, 2.0_rk], 'table: composed offsets', &
      & abs_tol=1e-15_rk)
    call assert_true(all(cons%all_cols == [118_ik]), 'table: one distinct column')

  end subroutine check_table_synthetic

  !> MPCs of the mesh used by the operator checks: x1 (x, c = 1), y1 (y, c = -0.5,
  !> d = 2e-3), z0 (z, c = 1.5) and x0y0 (x, c = 2): distinct components or disjoint sets.
  !>
  !> Arguments:
  !>   specs  (out)  type(mpc_spec_t)(:), allocatable  the relations
  subroutine mesh_specs(specs)

    type(mpc_spec_t), allocatable, intent(out) :: specs(:)

    allocate (specs(4))
    specs(1) = mpc_spec_t('x1', 1, 1.0_rk, 0.0_rk, 1_ik)
    specs(2) = mpc_spec_t('y1', 2, -0.5_rk, 2.0e-3_rk, 2_ik)
    specs(3) = mpc_spec_t('z0', 3, 1.5_rk, 0.0_rk, 3_ik)
    specs(4) = mpc_spec_t('x0y0', 1, 2.0_rk, 0.0_rk, 4_ik)

  end subroutine mesh_specs

  !> The table on the mesh, T / T^T adjointness, the reduced operator and the exact diagonal
  !> against the serial assembly with the explicit table.
  !>
  !> Locals:
  !>   specs      type(mpc_spec_t)(:)  relations
  !>   cons       type(constraint_t)   constraints
  !>   x, y, tx, tty, ky, d  real(rk)(:)  local vectors
  !>   x_ref, ttkx   real(rk)(:)  (3 num_nodes) whole-mesh vectors
  !>   mask       logical(:)           no prescribed dofs
  !>   e_q, te_q, kte_q      real(rk)(:)  unit vector of dof q, T e_q, K T e_q (whole mesh)
  !>   i, k, g, q     integer(ik)      indices
  !>   num_owned_rows integer(ik)      local rows of owned dofs
  !>   expected_rows integer(ik)       secondary dofs of the relations
  !>   num_rows   integer(ik)          local rows of the owned dofs, summed over the ranks
  !>   a, b       real(rk)             the two sides of the adjoint identity
  !>   scale, err real(rk)             reference scale and largest error
  !>   dq         real(rk)             (T e_q)^T K (T e_q)
  !>   ok         logical              all entries agree
  subroutine check_operators()

    type(mpc_spec_t), allocatable :: specs(:)
    type(constraint_t) :: cons
    real(rk), allocatable :: x(:), y(:), tx(:), tty(:), ky(:), d(:)
    real(rk), allocatable :: x_ref(:), ttkx(:), e_q(:), te_q(:), kte_q(:)
    logical, allocatable :: mask(:)
    integer(ik) :: i, k, g, q, expected_rows, num_rows, num_owned_rows
    real(rk) :: a, b, scale, err, dq
    logical :: ok

    ! (Allocated empty first: gfortran 15 warns of an unset descriptor otherwise.)
    allocate (specs(0))
    call mesh_specs(specs)
    call constraint_build_mpc(cons, mesh, specs)

    ! Rows: every node of each set but the first, for the set's component.
    expected_rows = 0_ik
    do i = 1_ik, int(size(specs), ik)
      do k = 1_ik, int(size(ref%nsets), ik)
        if (ref%nsets(k)%label == specs(i)%nset) expected_rows = expected_rows + &
          & int(size(ref%nsets(k)%nodes), ik) - 1_ik
      end do
    end do
    call assert_true(cons%active, 'mesh table: active')
    call assert_equal(int(cons%num_rows_glob), int(expected_rows), 'mesh table: one row per ' // &
      & 'secondary dof')
    num_owned_rows = int(count(cons%row_dof <= 3_ik*mesh%num_nodes_owned), ik)
    call par_sum(num_owned_rows, num_rows)
    call assert_equal(int(num_rows), int(expected_rows), 'local rows: each constrained dof ' // &
      & 'owned once')
    call assert_true(all(cons%glob_col == 3_ik*(11_ik - 1_ik) + 1_ik .or. &
      & cons%glob_col == 3_ik*(11_ik - 1_ik) + 2_ik .or. &
      & cons%glob_col == 3_ik*(ref%nsets(find_set('z0'))%nodes(1) - 1_ik) + 3_ik .or. &
      & cons%glob_col == 3_ik*(ref%nsets(find_set('x0y0'))%nodes(1) - 1_ik) + 1_ik), &
      & 'mesh table: the primaries are the first nodes of the node sets')

    ! Vectors from the global ids: x reduced (zero at the constrained dofs), y arbitrary.
    allocate (x(3*mesh%num_nodes_loc), y(3*mesh%num_nodes_loc), mask(3*mesh%num_nodes_loc))
    do i = 1_ik, mesh%num_nodes_loc
      do k = 1_ik, 3_ik
        x(3*(i - 1) + k) = sin(0.1_rk*mesh%node_id_glob(i) + k)
        y(3*(i - 1) + k) = cos(0.07_rk*mesh%node_id_glob(i) + 2*k)
      end do
    end do
    where (cons%is_constrained) x = 0.0_rk
    mask = .false.

    ! Adjointness over the owned dofs.
    tx = x
    call constraint_expand(cons, tx, .false.)
    tty = y
    call constraint_reduce(cons, 3_ik*mesh%num_nodes_owned, tty)
    a = dist_dot(tx, y, 3_ik*mesh%num_nodes_owned)
    b = dist_dot(x, tty, 3_ik*mesh%num_nodes_owned)
    call assert_close(a, b, 'T and T^T adjoint', rel_tol=1e-12_rk)
    call assert_true(all(abs(tty(cons%row_dof)) < tiny(1.0_rk)), 'T^T: zero at the ' // &
      & 'constrained dofs')

    ! The reduced operator against the serial T^T K T x.
    allocate (x_ref(3*ref%num_nodes))
    do g = 1_ik, ref%num_nodes
      x_ref(3*g - 2:3*g) = [(sin(0.1_rk*g + k), k = 1_ik, 3_ik)]
    end do
    x_ref(cons%glob_dof) = 0.0_rk
    ttkx = serial_ttkt(cons, x_ref)
    call halo_forward(mesh, x, 3_ik)
    allocate (ky(size(x)))
    call constraint_matvec(cons, mesh, k_elt, mask, x, ky)
    scale = maxval(abs(ttkx))
    ok = .true.
    do i = 1_ik, mesh%num_nodes_owned
      g = mesh%node_id_glob(i)
      if (any(abs(ky(3*i - 2:3*i) - ttkx(3*g - 2:3*g)) > 1e-12_rk*scale)) ok = .false.
    end do
    call assert_true(ok, 'reduced operator = serial T^T K T')

    ! Exact diagonal: every primary dof and every 37th dof, against (T e_q)^T K (T e_q).
    allocate (d(3*mesh%num_nodes_loc))
    call constraint_diagonal(cons, mesh, k_elt, d)
    call assert_true(all(abs(d(cons%row_dof)) < tiny(1.0_rk)), 'diagonal: zero at the ' // &
      & 'constrained dofs')
    allocate (e_q(3*ref%num_nodes))
    err = 0.0_rk
    scale = 0.0_rk
    do q = 1_ik, 3_ik*ref%num_nodes
      if (find_sorted(cons%glob_dof, q) > 0_ik) cycle
      if (find_sorted(cons%all_cols, q) == 0_ik .and. mod(q, 37_ik) /= 0_ik) cycle
      e_q = 0.0_rk
      e_q(q) = 1.0_rk
      te_q = serial_expand(cons, e_q)
      kte_q = serial_k(te_q)
      dq = dot_product(te_q, kte_q)
      scale = max(scale, abs(dq))
      i = node_local_owned((q - 1_ik)/3_ik + 1_ik)
      if (i > 0_ik) err = max(err, abs(d(3*(i - 1_ik) + mod(q - 1_ik, 3_ik) + 1_ik) - dq))
    end do
    call par_max(err, a)
    call assert_true(a <= 1e-12_rk*scale, 'exact diagonal = (T e_q)^T K (T e_q)')

  end subroutine check_operators

  !> Velocity conditions through the groups (C66): set_mpc1 x on x1 with a condition on its
  !> corner x1y1z1 (node 41, a secondary: x1's first node is 11) prescribes every x1 node at
  !> that value; the same value given on two nodes is averaged to itself.
  !>
  !> Locals:
  !>   specs    type(mpc_spec_t)(1)  set_mpc1 x on x1
  !>   cons     type(constraint_t)   constraints (save: gfortran 15 -O2 with -fcheck=all, the
  !>                                 Checked build, otherwise warns of uninitialized
  !>                                 descriptors, ROADMAP section G)
  !>   bcs      type(bc_spec_t)(:)   conditions
  !>   bc_mask  logical(:)           prescribed dofs
  !>   bc_vel   real(rk)(:)          prescribed values
  !>   i, n     integer(ik)          node indices
  !>   ok       logical              every x1 node prescribed at the value
  subroutine check_propagation()

    type(mpc_spec_t) :: specs(1)
    type(constraint_t), save :: cons
    type(bc_spec_t), allocatable :: bcs(:)
    logical, allocatable :: bc_mask(:)
    real(rk), allocatable :: bc_vel(:)
    integer(ik) :: i, n
    integer :: j
    logical :: ok

    specs(1) = mpc_spec_t('x1', 1, 1.0_rk, 0.0_rk, 1_ik)
    call constraint_build_mpc(cons, mesh, specs)
    allocate (bc_mask(3*mesh%num_nodes_loc), bc_vel(3*mesh%num_nodes_loc))
    do j = 1, 2
      if (allocated(bcs)) deallocate (bcs)
      allocate (bcs(j))
      bcs(1)%nset = 'x1y1z1'
      bcs(1)%comp = 1
      bcs(1)%value = 5.0e-3_rk
      if (j == 2) then
        bcs(2)%nset = 'x1y0z0'
        bcs(2)%comp = 1
        bcs(2)%value = 5.0e-3_rk
      end if
      call bc_resolve(mesh, bcs, bc_mask, bc_vel)
      call constraint_propagate_bcs(cons, mesh, bc_mask, bc_vel)
      ok = .true.
      do n = 1_ik, int(size(ref%nsets(find_set('x1'))%nodes), ik)
        i = node_local(ref%nsets(find_set('x1'))%nodes(n))
        if (i == 0_ik) cycle
        if (.not. bc_mask(3*i - 2) .or. abs(bc_vel(3*i - 2) - 5.0e-3_rk) > 1e-15_rk) ok = .false.
        if (bc_mask(3*i - 1) .or. bc_mask(3*i)) ok = .false.
      end do
      call assert_true(ok, 'propagation: the whole x1 group prescribed at the secondary''s ' // &
        & 'value')
      call assert_equal(int(count(bc_mask(1:3*mesh%num_nodes_owned))), &
        & int(count(cons%row_dof <= 3*mesh%num_nodes_owned)) + &
        & merge(1, 0, node_local_owned(11_ik) > 0_ik), 'propagation: only the group prescribed')
    end do

  end subroutine check_propagation

  !> Solves (V4) with the constrained PCG.
  !>
  !> Locals:
  !>   specs    type(mpc_spec_t)(:)  relations
  !>   cons     type(constraint_t)   constraints (save: gfortran 15 -O3 otherwise warns of
  !>                                 uninitialized descriptors, ROADMAP section G)
  !>   bcs      type(bc_spec_t)(:)   symmetry conditions along z
  !>   bc_mask  logical(:)           prescribed dofs
  !>   bc_vel, vel, rhs, delta, diag, r  real(rk)(:)  vectors
  !>   exact    real(rk)(3)          exact velocity of a node
  !>   err, scale, err_all  real(rk) errors and scales
  !>   num_iters, status  integer    PCG results
  !>   i, n, p, s  integer(ik)       indices
  !>   ok       logical              relations exact
  subroutine check_solves()

    type(mpc_spec_t), allocatable :: specs(:)
    type(constraint_t), save :: cons
    type(bc_spec_t), allocatable :: bcs(:)
    logical, allocatable :: bc_mask(:)
    real(rk), allocatable :: bc_vel(:), vel(:), rhs(:), delta(:), diag(:), r(:)
    real(rk) :: exact(3), err, scale, err_all, vp
    integer :: num_iters, status, case_i
    integer(ik) :: i, n, p
    logical :: ok

    allocate (bc_mask(3*mesh%num_nodes_loc), bc_vel(3*mesh%num_nodes_loc), &
      & vel(3*mesh%num_nodes_loc), rhs(3*mesh%num_nodes_loc), delta(3*mesh%num_nodes_loc), &
      & diag(3*mesh%num_nodes_loc), r(3*mesh%num_nodes_loc))
    allocate (bcs(0))
    call bc_preset_specs(bc_preset_symmetry, 3, rate, bcs)

    do case_i = 1, 2
      ! Case 1: set_mpc1 x on x1 (inactive for the homogeneous field); case 2: also y1 (y,
      ! c = -1).
      if (allocated(specs)) deallocate (specs)
      allocate (specs(case_i))
      specs(1) = mpc_spec_t('x1', 1, 1.0_rk, 0.0_rk, 1_ik)
      if (case_i == 2) specs(2) = mpc_spec_t('y1', 2, -1.0_rk, 0.0_rk, 2_ik)
      call constraint_build_mpc(cons, mesh, specs)
      call bc_resolve(mesh, bcs, bc_mask, bc_vel)
      call constraint_propagate_bcs(cons, mesh, bc_mask, bc_vel)

      ! Start: the prescribed values, projected; [T^T K T] dv = T^T (-K v0).
      vel = 0.0_rk
      where (bc_mask) vel = bc_vel
      call constraint_project(cons, vel)
      call ebe_matvec(mesh, k_elt, bc_mask, vel, r)
      rhs = -r
      delta = 0.0_rk
      call constraint_diagonal(cons, mesh, k_elt, diag)
      call pcg_solve(mesh, k_elt, bc_mask, diag, rhs, delta, 1e-12_rk, 16000, num_iters, &
        & status, cons)
      call assert_equal(status, 0, 'constrained solve converged')
      vel = vel + delta

      if (case_i == 1) then
        ! The exact homogeneous uniaxial-stress field.
        err = 0.0_rk
        do i = 1_ik, mesh%num_nodes_owned
          exact = -nu*rate*mesh%coo_s(:, i)
          exact(3) = rate*mesh%coo_s(3, i)
          err = max(err, maxval(abs(vel(3*i - 2:3*i) - exact)))
        end do
        call par_max(err, err_all)
        call assert_close(err_all, 0.0_rk, 'MPC solve: exact uniaxial field', abs_tol=1e-10_rk* &
          & rate)
      else
        ! The relations to round-off: vy(y1 node) = -vy(11) (both sets start
        ! with node 11), on the local nodes, where node 11 is local too.
        p = node_local(11_ik)
        ok = .true.
        do n = 1_ik, int(size(ref%nsets(find_set('y1'))%nodes), ik)
          i = node_local(ref%nsets(find_set('y1'))%nodes(n))
          if (i == 0_ik .or. i == p .or. p == 0_ik) cycle
          vp = proj_value(cons, vel, 3_ik*(i - 1_ik) + 2_ik)
          if (abs(vel(3*i - 1) + vp) > 1e-14_rk*rate) ok = .false.
        end do
        call assert_true(ok, 'MPC solve: vy(y1 node) = -vy(primary) to round-off')

        ! The reduced residual T^T (K v) at the free independent dofs.
        call halo_forward(mesh, vel, 3_ik)
        call ebe_matvec(mesh, k_elt, bc_mask, vel, r)
        scale = maxval(abs(r(1:3*mesh%num_nodes_owned)))
        call par_max(scale, err_all)
        scale = err_all
        call constraint_reduce(cons, 3_ik*mesh%num_nodes_owned, r)
        err = maxval(abs(r(1:3*mesh%num_nodes_owned)))
        call par_max(err, err_all)
        call assert_true(err_all <= 1e-8_rk*scale, 'MPC solve: T^T K v = 0 at the free dofs')
      end if
    end do

  end subroutine check_solves

  !> The value a constrained dof's relation gives from its column (one column: MPC), read from
  !> the local copy of the column (owned or ghost; the velocity is valid on both here).
  !>
  !> Arguments:
  !>   cons  (in)  type(constraint_t)  constraints
  !>   vel   (in)  real(rk)(:)         velocities (owned and ghost valid)
  !>   dof   (in)  integer(ik)         a local constrained dof
  !>   value (out, result) real(rk)    the column's velocity (0 if the column is not local)
  !> Locals:
  !>   r  integer(ik)  local row
  function proj_value(cons, vel, dof) result(value)

    type(constraint_t), intent(in) :: cons
    real(rk), intent(in) :: vel(:)
    integer(ik), intent(in) :: dof
    real(rk) :: value

    integer(ik) :: r

    r = cons%row_of_dof(dof)
    value = 0.0_rk
    if (r == 0_ik) return
    if (cons%col_local(cons%row_start(r)) == 0_ik) return
    value = vel(cons%col_local(cons%row_start(r)))

  end function proj_value

  !> Serial T x (with no offsets) on the whole mesh, from the replicated table.
  !>
  !> Arguments:
  !>   cons  (in)  type(constraint_t)  table
  !>   x     (in)  real(rk)(:)         (3 num_nodes) independent values (zero elsewhere)
  !>   tx    (out, result) real(rk)(:)
  !> Locals:
  !>   r, e  integer(ik)  row, column entry
  function serial_expand(cons, x) result(tx)

    type(constraint_t), intent(in) :: cons
    real(rk), intent(in) :: x(:)
    real(rk), allocatable :: tx(:)

    integer(ik) :: r, e

    tx = x
    do r = 1_ik, cons%num_rows_glob
      tx(cons%glob_dof(r)) = 0.0_rk
      do e = cons%glob_start(r), cons%glob_start(r + 1_ik) - 1_ik
        tx(cons%glob_dof(r)) = tx(cons%glob_dof(r)) + cons%glob_coef(e)*x(cons%glob_col(e))
      end do
    end do

  end function serial_expand

  !> Serial K x on the whole mesh (element by element, k_ref_all).
  !>
  !> Arguments:
  !>   x   (in)  real(rk)(:)  (3 num_nodes) vector
  !>   kx  (out, result) real(rk)(:)
  !> Locals:
  !>   e      integer          element
  !>   dofs   integer(ik)(:)   element dofs
  !>   a      integer          element node
  function serial_k(x) result(kx)

    real(rk), intent(in) :: x(:)
    real(rk), allocatable :: kx(:)

    integer :: e, a
    integer(ik), allocatable :: dofs(:)

    allocate (kx(size(x)), dofs(elt%num_dofs))
    kx = 0.0_rk
    do e = 1, int(ref%num_elts)
      do a = 1, int(elt%num_nodes)
        dofs(3*a - 2:3*a) = 3_ik*(ref%elt_nodes(a, e) - 1_ik) + [1_ik, 2_ik, 3_ik]
      end do
      kx(dofs) = kx(dofs) + matmul(k_ref_all(:, :, e), x(dofs))
    end do

  end function serial_k

  !> Serial T^T K T x on the whole mesh.
  !>
  !> Arguments:
  !>   cons  (in)  type(constraint_t)  table
  !>   x     (in)  real(rk)(:)         (3 num_nodes) independent values
  !>   y     (out, result) real(rk)(:)
  !> Locals:
  !>   r, e  integer(ik)  row, column entry
  function serial_ttkt(cons, x) result(y)

    type(constraint_t), intent(in) :: cons
    real(rk), intent(in) :: x(:)
    real(rk), allocatable :: y(:)

    integer(ik) :: r, e

    y = serial_k(serial_expand(cons, x))
    do r = 1_ik, cons%num_rows_glob
      do e = cons%glob_start(r), cons%glob_start(r + 1_ik) - 1_ik
        y(cons%glob_col(e)) = y(cons%glob_col(e)) + cons%glob_coef(e)*y(cons%glob_dof(r))
      end do
      y(cons%glob_dof(r)) = 0.0_rk
    end do

  end function serial_ttkt

  !> Index of a node set of the whole mesh by label.
  !>
  !> Arguments:
  !>   label  (in)  character  node-set label
  !>   idx    (out, result) integer
  !> Locals:
  !>   k  integer  node-set index
  function find_set(label) result(idx)

    character(len=*), intent(in) :: label
    integer :: idx

    integer :: k

    idx = 0
    do k = 1, size(ref%nsets)
      if (ref%nsets(k)%label == label) idx = k
    end do

  end function find_set

  !> Local number of a global node (owned or ghost), 0 if not local.
  !>
  !> Arguments:
  !>   g  (in)  integer(ik)  global node id
  !>   i  (out, result) integer(ik)
  !> Locals:
  !>   k  integer(ik)  local node
  function node_local(g) result(i)

    integer(ik), intent(in) :: g
    integer(ik) :: i

    integer(ik) :: k

    i = 0_ik
    do k = 1_ik, mesh%num_nodes_loc
      if (mesh%node_id_glob(k) == g) i = k
    end do

  end function node_local

  !> Local number of a global node owned by this rank, 0 otherwise.
  !>
  !> Arguments:
  !>   g  (in)  integer(ik)  global node id
  !>   i  (out, result) integer(ik)
  function node_local_owned(g) result(i)

    integer(ik), intent(in) :: g
    integer(ik) :: i

    i = node_local(g)
    if (i > mesh%num_nodes_owned) i = 0_ik

  end function node_local_owned

  !> Isotropic Mandel stiffness lambda {I}{I}^T + 2 mu [I].
  !>
  !> Arguments:
  !>   e_mod, poisson  (in)  real(rk)  Young's modulus [MPa], Poisson's ratio [-]
  !> Locals:
  !>   lam, mu   real(rk)  Lame constants [MPa]
  !>   i, j      integer   indices
  pure function isotropic_mandel(e_mod, poisson) result(stiff)

    real(rk), intent(in) :: e_mod, poisson
    real(rk) :: stiff(6, 6)

    real(rk) :: lam, mu
    integer :: i, j

    lam = e_mod*poisson/((1.0_rk + poisson)*(1.0_rk - 2.0_rk*poisson))
    mu = e_mod/(2.0_rk*(1.0_rk + poisson))
    stiff = 0.0_rk
    do i = 1, 3
      do j = 1, 3
        stiff(i, j) = lam
      end do
    end do
    do i = 1, 6
      stiff(i, i) = stiff(i, i) + 2.0_rk*mu
    end do

  end function isotropic_mandel

end program test_constraint
