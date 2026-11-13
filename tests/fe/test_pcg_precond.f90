! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> MPI tests of the CG preconditioners (pcg_mod: Jacobi, nodal block-Jacobi, none) and of the
!> nodal 3 x 3 diagonal blocks (ebe_mod::ebe_block_diagonal), on 1, 2 and 4 ranks.
!>
!> Source : FEPX (ROADMAP step 9.2; dev/theory/addendum_block_jacobi.md §6, tests T1-T5).
!>
!> Usage: mpiexec -n N test_pcg_precond <regression mesh dir> <verification mesh dir>
!>        (tests/regression/cases/meshes, tests/verification/cases/meshes).
!>
!> Meshes: poly_part1.msh (unit cube, 656 tet10) and annulus_hex20_n4.msh (hex20), distributed
!> with partitioning "even"; material: an anisotropic [s] (isotropic, E = 200 GPa, nu = 0.3,
!> plus a rank-one term coupling every component) at every point, so that the nodal blocks are
!> full. Checks:
!>   - T1 blocks: ebe_block_diagonal equals the 3 x 3 diagonal blocks of a serial assembly of
!>     the whole mesh read on every rank (relative 1e-12), tet10 and hex20, on the owned and the
!>     ghost nodes;
!>   - T2 preconditioner step: one PCG iteration gives x = alpha M^-1 r, so x is parallel to the
!>     serial M^-1 r (relative 1e-11) for block-Jacobi (blocks restricted to the free components:
!>     identity rows and columns at the prescribed ones; z zero there), Jacobi (r / diag) and none
!>     (r); with a mask that prescribes one, two or three components of different nodes;
!>   - T3 invariance: for nodal rotations Q = blockdiag(R_a) (R_a from the global node id), PCG
!>     with block-Jacobi and with none on Q^T K Q, Q^T r gives Q^T times the iterates of K, r
!>     after 8 iterations (relative 1e-10); Jacobi does not (difference above 1e-6: the test is
!>     sensitive);
!>   - T4 solves: block-Jacobi and none reach the exact homogeneous uniaxial-stress field (isotropic
!>     material, minimal preset along x, y, z) to 1e-10 e with cg_tol 1e-12; on the anisotropic
!>     material the block-Jacobi and the Jacobi solutions (cg_tol 1e-12) agree to 1e-9 relative;
!>   - status 2: element matrices negated (every block negative definite): block-Jacobi returns
!>     status 2 with no iteration, on every rank.
!>   T5 (rank independence): the same checks on 1, 2 and 4 ranks.
!>
!> Locals:
!>   ework  type(element_work_t)  element work arrays (element_work_alloc; D100)
!>   young, nu, rate   real(rk)        Young's modulus [MPa], Poisson's ratio, strain rate [1/s]
!>   meshes_reg, meshes_ver  character(512)  mesh directories
!>   mesh              type(mesh_t)       rank-local mesh (the current one)
!>   ref               type(mesh_file_t)  whole mesh (serial references)
!>   elt               type(element_t)    element descriptor
!>   mat_iso, mat_aniso  real(rk)(6,6)    isotropic and anisotropic [s] [MPa s]
!>   k_elt             real(rk)(:,:,:)    local element matrices (anisotropic)
!>   diag              real(rk)(:)        assembled diagonal
!>   unit_log          integer            scratch unit receiving the log
program test_pcg_precond

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_init, par_finalize, par_sum, par_max, par_all
  use mesh_file_mod, only: mesh_file_t, mesh_file_load
  use mesh_mod, only: mesh_t
  use mesh_blocks_mod, only: mesh_load_blocks
  use element_mod, only: element_t, element_init
  use element_matrix_mod, only: element_stiffness
  use ebe_mod, only: ebe_matvec, ebe_diagonal, ebe_block_diagonal
  use pcg_mod, only: pcg_solve, pcg_precond_none, pcg_precond_jacobi, pcg_precond_block_jacobi
  use linalg_mod, only: ldlt_factor, ldlt_substitute
  use bc_mod, only: bc_spec_t, bc_preset_specs, bc_resolve, bc_gage_length, bc_preset_minimal
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary
  use element_matrix_mod, only: element_work_t, element_work_alloc

  implicit none

  type(element_work_t) :: ework
  real(rk), parameter :: young = 200.0e3_rk
  real(rk), parameter :: nu = 0.3_rk
  real(rk), parameter :: rate = 1.0e-3_rk

  character(len=512) :: meshes_reg, meshes_ver
  type(mesh_t) :: mesh
  type(mesh_file_t) :: ref
  type(element_t) :: elt
  real(rk) :: mat_iso(6, 6), mat_aniso(6, 6)
  real(rk), allocatable :: k_elt(:, :, :), diag(:)
  integer :: unit_log

  call par_init()
  call get_command_argument(1, meshes_reg)
  call get_command_argument(2, meshes_ver)
  mat_iso = isotropic_mandel(young, nu)
  mat_aniso = anisotropic_mandel()

  ! The tet10 cube: every check.
  call load(trim(meshes_reg) // '/poly_part1.msh')
  call build_matrices(mat_aniso)
  call check_blocks('tet10')
  call check_precond_step()
  call check_invariance()
  call check_solves()
  call check_not_positive_definite()

  ! The hex20 annulus: the blocks and the preconditioner step.
  call load(trim(meshes_ver) // '/annulus_hex20_n4.msh')
  call build_matrices(mat_aniso)
  call check_blocks('hex20')
  call check_precond_step()

  call test_summary()
  call par_finalize()

contains

  !> Load a mesh distributed (log to a scratch unit) and whole (serial references).
  !>
  !> Arguments:
  !>   path  (in)  character  mesh file
  subroutine load(path)

    character(len=*), intent(in) :: path

    open (newunit=unit_log, status='scratch', action='readwrite', form='formatted')
    call log_set_unit(unit_log)
    call mesh_load_blocks(mesh, path)
    call mesh_file_load(ref, path)
    call log_set_unit(output_unit)
    close (unit_log)
    call element_init(elt, mesh%elt_type)

  end subroutine load

  !> Element matrices of the local elements for a material, and the assembled diagonal.
  !>
  !> Arguments:
  !>   mat     (in)  real(rk)(6,6)  [s] at every point [MPa s]
  !> Locals:
  !>   mat_s   real(rk)(:,:,:)      [s] at the quadrature points
  !>   i_elt   integer(ik)          local element
  subroutine build_matrices(mat)

    real(rk), intent(in) :: mat(6, 6)

    real(rk), allocatable :: mat_s(:, :, :)
    integer(ik) :: i_elt

    if (allocated(k_elt)) deallocate (k_elt, diag)
    allocate (k_elt(elt%num_dofs, elt%num_dofs, mesh%num_elts), diag(3*mesh%num_nodes_loc))
    mat_s = spread(mat, 3, int(elt%num_qpts))
    do i_elt = 1_ik, mesh%num_elts
      call element_work_alloc(ework, elt)
      call element_stiffness(elt, mesh%coo_s(:, mesh%elt_nodes_loc(:, i_elt)), mat_s, &
        & k_elt(:, :, i_elt), ework)
    end do
    call ebe_diagonal(mesh, k_elt, diag)

  end subroutine build_matrices

  !> The serial nodal blocks of the whole mesh (anisotropic material).
  !>
  !> Arguments:
  !>   blocks_ref  (out) real(rk)(:,:,:)  (3, 3, num_nodes) blocks
  !> Locals:
  !>   mat_s   real(rk)(:,:,:)  [s] at the quadrature points
  !>   k_ref   real(rk)(:,:)    element matrix
  !>   e, i    integer          element, element node
  !>   g       integer(ik)      global node
  subroutine serial_blocks(blocks_ref)

    real(rk), allocatable, intent(out) :: blocks_ref(:, :, :)

    real(rk), allocatable :: mat_s(:, :, :), k_ref(:, :)
    integer :: e, i
    integer(ik) :: g

    allocate (blocks_ref(3, 3, ref%num_nodes), k_ref(elt%num_dofs, elt%num_dofs))
    mat_s = spread(mat_aniso, 3, int(elt%num_qpts))
    blocks_ref = 0.0_rk
    do e = 1, int(ref%num_elts)
      call element_work_alloc(ework, elt)
      call element_stiffness(elt, ref%coo(:, ref%elt_nodes(:, e)), mat_s, k_ref, ework)
      do i = 1, int(elt%num_nodes)
        g = ref%elt_nodes(i, e)
        blocks_ref(:, :, g) = blocks_ref(:, :, g) + k_ref(3*i - 2:3*i, 3*i - 2:3*i)
      end do
    end do

  end subroutine serial_blocks

  !> T1: the distributed blocks against the serial assembly, on the owned and ghost nodes.
  !>
  !> Arguments:
  !>   label       (in)  character        element type, for the messages
  !> Locals:
  !>   blocks      real(rk)(:)            (9 num_nodes_loc) distributed blocks
  !>   blocks_ref  real(rk)(:,:,:)        serial blocks
  !>   i           integer                local node
  !>   scale       real(rk)               largest serial entry
  !>   ok, all_ok  logical                agreement on this rank, on every rank
  subroutine check_blocks(label)

    character(len=*), intent(in) :: label

    real(rk), allocatable :: blocks(:), blocks_ref(:, :, :)
    integer :: i
    real(rk) :: scale
    logical :: ok, all_ok

    allocate (blocks(9*mesh%num_nodes_loc))
    call ebe_block_diagonal(mesh, k_elt, blocks)
    call serial_blocks(blocks_ref)
    scale = maxval(abs(blocks_ref))
    ok = .true.
    do i = 1, int(mesh%num_nodes_loc)
      if (any(abs(reshape(blocks(9*i - 8:9*i), [3, 3]) - &
        & blocks_ref(:, :, mesh%node_id_glob(i))) > 1e-12_rk*scale)) ok = .false.
    end do
    call par_all(ok, all_ok)
    call assert_true(all_ok, 'blocks (' // label // '): = serial assembly, owned and ghosts')
    ! (The blocks are full: an off-diagonal entry of the order of the diagonal.)
    call assert_true(maxval(abs(blocks_ref(1, 2, :))) > 1e-3_rk*maxval(abs(blocks_ref(1, 1, :))), &
      & 'blocks (' // label // '): coupled components')

  end subroutine check_blocks

  !> A mask prescribing one, two or three components of different nodes (by global id, the
  !> same on every rank count): node id mod 7 = 0: component 1; = 1: components 1, 2; = 2: all.
  !>
  !> Arguments:
  !>   g     (in)  integer(ik)  global node id
  !>   mask  (out) logical(3)   prescribed components
  pure subroutine node_mask(g, mask)

    integer(ik), intent(in) :: g
    logical, intent(out) :: mask(3)

    mask = .false.
    select case (mod(g, 7_ik))
    case (0_ik)
      mask(1) = .true.
    case (1_ik)
      mask(1:2) = .true.
    case (2_ik)
      mask = .true.
    end select

  end subroutine node_mask

  !> T2: one PCG iteration from x = 0 gives x = alpha M^-1 r; compared with the serial M^-1 r
  !> for each preconditioner.
  !>
  !> Locals:
  !>   bc_mask      logical(:)        prescribed dofs
  !>   rhs, x       real(rk)(:)       right-hand side and the iterate
  !>   z_ref        real(rk)(:)       serial M^-1 r on the local dofs
  !>   blocks_ref   real(rk)(:,:,:)   serial blocks
  !>   block        real(rk)(3,3)     a restricted block
  !>   m            logical(3)        a node's prescribed components
  !>   r_a, z_a     real(rk)(3)       a node's residual and preconditioned residual
  !>   i, k, p      integer           local node, component, preconditioner
  !>   num_iters, status  integer     PCG results
  !>   ldlt_status  integer           ldlt_factor status
  !>   dots, sums   real(rk)(2)       x.z_ref and z_ref.z_ref on this rank, over the ranks
  !>   c, err, err_max, scale  real(rk)  factor alpha, largest deviation, its maximum, scale
  !>   names        character(12)(3)  the preconditioners' names
  subroutine check_precond_step()

    logical, allocatable :: bc_mask(:)
    real(rk), allocatable :: rhs(:), x(:), z_ref(:), blocks_ref(:, :, :)
    real(rk) :: block(3, 3), r_a(3), z_a(3), dots(2), sums(2), c, err, err_max, scale
    logical :: m(3)
    integer :: i, k, p, num_iters, status, ldlt_status
    integer, parameter :: codes(3) = [pcg_precond_block_jacobi, pcg_precond_jacobi, &
      & pcg_precond_none]
    character(len=12), parameter :: names(3) = [character(len=12) :: 'block_jacobi', 'jacobi', &
      & 'none']

    allocate (bc_mask(3*mesh%num_nodes_loc), rhs(3*mesh%num_nodes_loc), &
      & x(3*mesh%num_nodes_loc), z_ref(3*mesh%num_nodes_loc))
    call serial_blocks(blocks_ref)

    ! Mask and right-hand side from the global ids.
    do i = 1, int(mesh%num_nodes_loc)
      call node_mask(mesh%node_id_glob(i), m)
      bc_mask(3*i - 2:3*i) = m
      do k = 1, 3
        rhs(3*(i - 1) + k) = sin(0.37_rk*mesh%node_id_glob(i) + 1.3_rk*k)
      end do
    end do

    do p = 1, 3
      ! Serial M^-1 r at the owned nodes (r = rhs, zero at the prescribed dofs).
      z_ref = 0.0_rk
      do i = 1, int(mesh%num_nodes_owned)
        m = bc_mask(3*i - 2:3*i)
        r_a = merge(0.0_rk, rhs(3*i - 2:3*i), m)
        select case (codes(p))
        case (pcg_precond_block_jacobi)
          block = blocks_ref(:, :, mesh%node_id_glob(i))
          do k = 1, 3
            if (m(k)) then
              block(k, :) = 0.0_rk
              block(:, k) = 0.0_rk
              block(k, k) = 1.0_rk
            end if
          end do
          call ldlt_factor(block, ldlt_status)
          z_a = r_a
          call ldlt_substitute(block, z_a)
        case (pcg_precond_jacobi)
          z_a = r_a/[(blocks_ref(k, k, mesh%node_id_glob(i)), k = 1, 3)]
        case default
          z_a = r_a
        end select
        z_ref(3*i - 2:3*i) = z_a
      end do

      ! One iteration (limit 1: status 1).
      x = 0.0_rk
      call pcg_solve(mesh, k_elt, bc_mask, diag, rhs, x, 1e-14_rk, 1, num_iters, status, &
        & precond=codes(p))
      call assert_equal(status, 1, 'precond step (' // trim(names(p)) // '): one iteration')

      ! x = c z_ref on the owned dofs.
      dots = [sum(x(1:3*mesh%num_nodes_owned)*z_ref(1:3*mesh%num_nodes_owned)), &
        & sum(z_ref(1:3*mesh%num_nodes_owned)**2)]
      call par_sum(dots, sums)
      c = sums(1)/sums(2)
      call par_max(maxval(abs(x(1:3*mesh%num_nodes_owned) - c*z_ref(1:3*mesh%num_nodes_owned))), &
        & err_max)
      call par_max(maxval(abs(c*z_ref)), scale)
      err = err_max/scale
      call assert_close(err, 0.0_rk, 'precond step (' // trim(names(p)) // '): x = alpha M^-1 r', &
        & abs_tol=1e-11_rk)
      call assert_true(all(abs(x(1:3*mesh%num_nodes_owned)) <= 0.0_rk .or. &
        & .not. bc_mask(1:3*mesh%num_nodes_owned)), 'precond step (' // trim(names(p)) // &
        & '): zero at the prescribed dofs')
    end do

  end subroutine check_precond_step

  !> Rotation of node g (global id): angle 0.3 + 0.01 g about the axis (sin g, cos 2g, 1),
  !> normalized (Rodrigues' formula).
  !>
  !> Arguments:
  !>   g    (in)  integer(ik)    global node id
  !>   rot  (out, result) real(rk)(3,3)  rotation matrix
  !> Locals:
  !>   axis   real(rk)(3)    unit axis
  !>   angle  real(rk)       angle [rad]
  !>   w      real(rk)(3,3)  skew matrix of the axis
  !>   k      integer        index
  pure function node_rotation(g) result(rot)

    integer(ik), intent(in) :: g
    real(rk) :: rot(3, 3)

    real(rk) :: axis(3), angle, w(3, 3)
    integer :: k

    axis = [sin(real(g, rk)), cos(2.0_rk*g), 1.0_rk]
    axis = axis/norm2(axis)
    angle = 0.3_rk + 0.01_rk*g
    w = reshape([0.0_rk, axis(3), -axis(2), -axis(3), 0.0_rk, axis(1), axis(2), -axis(1), &
      & 0.0_rk], [3, 3])
    rot = sin(angle)*w + (1.0_rk - cos(angle))*matmul(w, w)
    do k = 1, 3
      rot(k, k) = rot(k, k) + 1.0_rk
    end do

  end function node_rotation

  !> T3: PCG on Q^T K Q, Q^T r after 8 iterations equals Q^T times the iterate on K, r for
  !> block-Jacobi and none; not for Jacobi. Mask: all three components of the nodes on z = 0
  !> (unchanged by the rotations, and K positive definite on the free dofs).
  !>
  !> Locals:
  !>   k_rot        real(rk)(:,:,:)   rotated element matrices
  !>   q_e          real(rk)(:,:)     an element's blockdiag(R_a)
  !>   rot          real(rk)(:,:,:)   (3, 3, num_nodes_loc) node rotations
  !>   diag_rot     real(rk)(:)       diagonal of the rotated matrix
  !>   bc_mask      logical(:)        prescribed dofs
  !>   rhs, rhs_rot, x, x_rot  real(rk)(:)  right-hand sides and iterates
  !>   i, e, p      integer           local node, element node, preconditioner
  !>   i_elt        integer(ik)       local element
  !>   num_iters, status  integer     PCG results
  !>   err, scale   real(rk)          largest |x_rot - Q^T x| and |x| over the ranks
  !>   err_all      real(rk)          err's maximum over the ranks
  subroutine check_invariance()

    real(rk), allocatable :: k_rot(:, :, :), q_e(:, :), rot(:, :, :), diag_rot(:), rhs(:), &
      & rhs_rot(:), x(:), x_rot(:)
    logical, allocatable :: bc_mask(:)
    integer :: i, e, p, num_iters, status
    integer(ik) :: i_elt
    real(rk) :: err, err_all, scale
    integer, parameter :: codes(3) = [pcg_precond_block_jacobi, pcg_precond_none, &
      & pcg_precond_jacobi]
    character(len=12), parameter :: names(3) = [character(len=12) :: 'block_jacobi', 'none', &
      & 'jacobi']

    allocate (k_rot, mold=k_elt)
    allocate (q_e(elt%num_dofs, elt%num_dofs), rot(3, 3, mesh%num_nodes_loc), &
      & diag_rot(3*mesh%num_nodes_loc), bc_mask(3*mesh%num_nodes_loc), &
      & rhs(3*mesh%num_nodes_loc), rhs_rot(3*mesh%num_nodes_loc), x(3*mesh%num_nodes_loc), &
      & x_rot(3*mesh%num_nodes_loc))

    ! Rotations, mask, right-hand side and its rotation Q^T r.
    do i = 1, int(mesh%num_nodes_loc)
      rot(:, :, i) = node_rotation(mesh%node_id_glob(i))
      bc_mask(3*i - 2:3*i) = abs(mesh%coo_s(3, i)) < 1e-9_rk
      rhs(3*i - 2:3*i) = [sin(0.37_rk*mesh%node_id_glob(i)), cos(0.11_rk*mesh%node_id_glob(i)), &
        & sin(0.05_rk*mesh%node_id_glob(i) + 1.0_rk)]
      rhs_rot(3*i - 2:3*i) = matmul(transpose(rot(:, :, i)), rhs(3*i - 2:3*i))
    end do

    ! Rotated element matrices Q_e^T k_e Q_e.
    do i_elt = 1_ik, mesh%num_elts
      q_e = 0.0_rk
      do e = 1, int(elt%num_nodes)
        q_e(3*e - 2:3*e, 3*e - 2:3*e) = rot(:, :, mesh%elt_nodes_loc(e, i_elt))
      end do
      k_rot(:, :, i_elt) = matmul(transpose(q_e), matmul(k_elt(:, :, i_elt), q_e))
    end do
    call ebe_diagonal(mesh, k_rot, diag_rot)

    do p = 1, 3
      x = 0.0_rk
      call pcg_solve(mesh, k_elt, bc_mask, diag, rhs, x, 1e-14_rk, 8, num_iters, status, &
        & precond=codes(p))
      x_rot = 0.0_rk
      call pcg_solve(mesh, k_rot, bc_mask, diag_rot, rhs_rot, x_rot, 1e-14_rk, 8, num_iters, &
        & status, precond=codes(p))
      err = 0.0_rk
      do i = 1, int(mesh%num_nodes_owned)
        err = max(err, maxval(abs(x_rot(3*i - 2:3*i) - &
          & matmul(transpose(rot(:, :, i)), x(3*i - 2:3*i)))))
      end do
      call par_max(err, err_all)
      err = err_all
      call par_max(maxval(abs(x(1:3*mesh%num_nodes_owned))), scale)
      if (codes(p) == pcg_precond_jacobi) then
        call assert_true(err > 1e-6_rk*scale, 'invariance (' // trim(names(p)) // &
          & '): not invariant (the test is sensitive)')
      else
        call assert_close(err/scale, 0.0_rk, 'invariance (' // trim(names(p)) // &
          & '): Q^T K Q iterates = Q^T iterates', abs_tol=1e-10_rk)
      end if
    end do

  end subroutine check_invariance

  !> T4: exact uniaxial fields (isotropic material) with block-Jacobi and none; block-Jacobi
  !> and Jacobi solutions agree on the anisotropic material. The anisotropic matrices are
  !> rebuilt at the end.
  !>
  !> Locals:
  !>   specs      type(bc_spec_t)(:)  conditions
  !>   bc_mask    logical(:)          prescribed dofs
  !>   bc_vel, rhs, delta, delta_j, y  real(rk)(:)  vectors
  !>   exact      real(rk)(3)         exact velocity of a node
  !>   axis, i, p integer             loading axis, local node, preconditioner
  !>   num_iters, status  integer     PCG results
  !>   err, scale real(rk)            largest error over the ranks, scale
  !>   err_all    real(rk)            err's maximum over the ranks
  subroutine check_solves()

    type(bc_spec_t), allocatable :: specs(:)
    logical, allocatable :: bc_mask(:)
    real(rk), allocatable :: bc_vel(:), rhs(:), delta(:), delta_j(:), y(:)
    real(rk) :: exact(3), err, err_all, scale
    integer :: axis, i, p, num_iters, status
    integer, parameter :: codes(2) = [pcg_precond_block_jacobi, pcg_precond_none]
    character(len=12), parameter :: names(2) = [character(len=12) :: 'block_jacobi', 'none']

    allocate (bc_mask(3*mesh%num_nodes_loc), bc_vel(3*mesh%num_nodes_loc), &
      & rhs(3*mesh%num_nodes_loc), delta(3*mesh%num_nodes_loc), &
      & delta_j(3*mesh%num_nodes_loc), y(3*mesh%num_nodes_loc))
    ! (Allocated empty first: gfortran 15 -O3 warns of an unset descriptor otherwise.)
    allocate (specs(0))

    ! Isotropic material: the exact homogeneous uniaxial-stress field.
    call build_matrices(mat_iso)
    do p = 1, 2
      do axis = 1, 3
        call bc_preset_specs(bc_preset_minimal, axis, rate*bc_gage_length(mesh, axis), specs)
        call bc_resolve(mesh, specs, bc_mask, bc_vel)
        call ebe_matvec(mesh, k_elt, bc_mask, bc_vel, y)
        rhs = -y
        delta = 0.0_rk
        call pcg_solve(mesh, k_elt, bc_mask, diag, rhs, delta, 1e-12_rk, 16000, num_iters, &
          & status, precond=codes(p))
        call assert_equal(status, 0, 'solve (' // trim(names(p)) // '): converged')
        err = 0.0_rk
        do i = 1, int(mesh%num_nodes_owned)
          exact = -nu*rate*mesh%coo_s(:, i)
          exact(axis) = rate*mesh%coo_s(axis, i)
          err = max(err, maxval(abs(bc_vel(3*i - 2:3*i) + delta(3*i - 2:3*i) - exact)))
        end do
        call par_max(err, err_all)
        err = err_all
        call assert_close(err, 0.0_rk, 'solve (' // trim(names(p)) // '): exact uniaxial field', &
          & abs_tol=1e-10_rk*rate)
      end do
    end do

    ! Anisotropic material: block-Jacobi and Jacobi solutions agree.
    call build_matrices(mat_aniso)
    call bc_preset_specs(bc_preset_minimal, 3, rate*bc_gage_length(mesh, 3), specs)
    call bc_resolve(mesh, specs, bc_mask, bc_vel)
    call ebe_matvec(mesh, k_elt, bc_mask, bc_vel, y)
    rhs = -y
    delta = 0.0_rk
    call pcg_solve(mesh, k_elt, bc_mask, diag, rhs, delta, 1e-12_rk, 16000, num_iters, status, &
      & precond=pcg_precond_block_jacobi)
    delta_j = 0.0_rk
    call pcg_solve(mesh, k_elt, bc_mask, diag, rhs, delta_j, 1e-12_rk, 16000, num_iters, &
      & status, precond=pcg_precond_jacobi)
    call par_max(maxval(abs(delta(1:3*mesh%num_nodes_owned) - &
      & delta_j(1:3*mesh%num_nodes_owned))), err)
    call par_max(maxval(abs(delta_j(1:3*mesh%num_nodes_owned))), scale)
    call assert_close(err/scale, 0.0_rk, 'solve (anisotropic): block_jacobi = jacobi', &
      & abs_tol=1e-9_rk)

  end subroutine check_solves

  !> Status 2: with every element matrix negated, every block is negative definite.
  !>
  !> Locals:
  !>   bc_mask   logical(:)    no prescribed dofs
  !>   rhs, x    real(rk)(:)   vectors
  !>   num_iters, status  integer  PCG results
  subroutine check_not_positive_definite()

    logical, allocatable :: bc_mask(:)
    real(rk), allocatable :: rhs(:), x(:)
    integer :: num_iters, status

    allocate (bc_mask(3*mesh%num_nodes_loc), rhs(3*mesh%num_nodes_loc), &
      & x(3*mesh%num_nodes_loc))
    bc_mask = .false.
    rhs = 1.0_rk
    x = 0.0_rk
    call pcg_solve(mesh, -k_elt, bc_mask, diag, rhs, x, 1e-8_rk, 16000, num_iters, status, &
      & precond=pcg_precond_block_jacobi)
    call assert_equal(status, 2, 'block not positive definite: status 2')
    call assert_equal(num_iters, 0, 'block not positive definite: no iteration')

  end subroutine check_not_positive_definite

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

  !> Anisotropic Mandel stiffness: the isotropic one plus 1e5 v v^T, v coupling every
  !> component (symmetric positive definite).
  !>
  !> Locals:
  !>   v  real(rk)(6)  coupling vector
  !>   i  integer      index
  function anisotropic_mandel() result(stiff)

    real(rk) :: stiff(6, 6)

    real(rk) :: v(6)
    integer :: i

    v = [1.0_rk, 0.5_rk, -0.3_rk, 0.2_rk, 0.7_rk, -0.4_rk]
    stiff = isotropic_mandel(young, nu)
    do i = 1, 6
      stiff(:, i) = stiff(:, i) + 1.0e5_rk*v*v(i)
    end do

  end function anisotropic_mandel

end program test_pcg_precond
