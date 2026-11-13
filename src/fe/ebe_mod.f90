! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Element-by-element (EBE) global operations on the rank-local mesh: matrix-vector product with
!> the stored element matrices, assembly of the global diagonal and of element vectors.
!>
!> Distributed vectors have the local dof layout of mesh_mod (dof 3 (i - 1) + k of local node i;
!> owned nodes first, then ghosts). A vector is complete on its owned dofs after these routines;
!> its ghost dofs are made consistent by halo_forward when needed (before a matrix-vector
!> product, whose input must be valid on the ghosts). Prescribed (Dirichlet) dofs are given by a
!> mask: the matrix-vector product returns zero there, so that, with a right-hand side that is
!> zero there too, CG works on the free dofs only (v2's treatment).
!>
!> Source : port of v2 matrix_operations_mod.f90::sparse_matvec_ebe (gather to the elements,
!>          element matrix times element vector, scatter-add, zero at the prescribed dofs) and
!>          conjugate_gradient_mod1.f90::assemble_diagonals ("general": scatter-add of the element
!>          diagonals); v2's part_gather/part_scatter are replaced by the local connectivity and
!>          halo_mod (step 2.4). The nodal 3 x 3 diagonal blocks (step 9.2, block-Jacobi
!>          preconditioner): dev/theory/addendum_block_jacobi.md §2 (Barrett et al. 1994,
!>          Templates, §3.2.1, blocks formed by grouping the equations per node), the same
!>          scatter-add with the nine entries of each node's block.
!>
!> Public routines:
!>   ebe_matvec          y = [K] x, element by element, zero at the prescribed dofs
!>   ebe_diagonal        assembled diagonal of [K]
!>   ebe_block_diagonal  assembled 3 x 3 diagonal blocks of [K], one per node
!>   ebe_assemble        assembled vector of element vectors
!>   elt_dof_list        local dofs of an element (3 per node, node by node)
module ebe_mod

  use kinds_mod, only: rk, ik
  use mesh_mod, only: mesh_t
  use halo_mod, only: halo_forward, halo_reverse_add

  implicit none

  private

  public :: ebe_matvec, ebe_diagonal, ebe_block_diagonal, ebe_assemble, elt_dof_list

contains

  !> Local dofs of element i_elt: dof 3 (a_loc - 1) + k for node a of the element (local number
  !> a_loc) and component k, in the element's node order (mesh_t%elt_dofs_loc, built by
  !> mesh_localize).
  !>
  !> Arguments:
  !>   mesh   (in)  type(mesh_t)  rank-local mesh
  !>   i_elt  (in)  integer(ik)   local element
  !>   dofs   (out, result) integer(ik)(3 num_elt_nodes)  local dofs
  pure function elt_dof_list(mesh, i_elt) result(dofs)

    type(mesh_t), intent(in) :: mesh
    integer(ik), intent(in) :: i_elt
    integer(ik) :: dofs(3*mesh%num_elt_nodes)

    dofs = mesh%elt_dofs_loc(:, i_elt)

  end function elt_dof_list

  !> y = [K] x, element by element (v2 sparse_matvec_ebe): x must be valid on the ghost dofs;
  !> y is complete on the owned dofs and zero at the prescribed dofs (its ghost dofs hold
  !> partial sums and are not meaningful).
  !>
  !> Arguments:
  !>   mesh      (in)  type(mesh_t)        rank-local mesh
  !>   k_elt     (in)  real(rk)(:, :, :)   (num_dofs_elt, num_dofs_elt, num_elts) element matrices
  !>   bc_mask   (in)  logical(:)          (3 num_nodes_loc) prescribed dofs
  !>   x         (in)  real(rk)(:)         (3 num_nodes_loc) input vector
  !>   y         (out) real(rk)(:)         (3 num_nodes_loc) result
  !> Locals:
  !>   i_elt     integer(ik)               element
  !>   i, j      integer(ik)               element dofs (row, column)
  !>   n         integer(ik)               dofs per element
  !>   xe, ye    real(rk)(n)               element parts of x and y
  !>
  !> Loop as v2 gen_matrix_vector_mult (element by element, the product accumulated column by
  !> column) on the stored dof table, without array temporaries (step 4.4 performance gate).
  subroutine ebe_matvec(mesh, k_elt, bc_mask, x, y)

    type(mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: k_elt(:, :, :)
    logical, intent(in) :: bc_mask(:)
    real(rk), intent(in) :: x(:)
    real(rk), intent(out) :: y(:)

    integer(ik) :: i_elt, i, j, n
    real(rk) :: xe(size(k_elt, 1)), ye(size(k_elt, 1))

    ! Gather, multiply (column by column), scatter-add.
    n = size(k_elt, 1, kind=ik)
    y = 0.0_rk
    do i_elt = 1_ik, mesh%num_elts
      do j = 1_ik, n
        xe(j) = x(mesh%elt_dofs_loc(j, i_elt))
      end do
      ye = 0.0_rk
      do j = 1_ik, n
        ye = ye + k_elt(:, j, i_elt)*xe(j)
      end do
      do i = 1_ik, n
        y(mesh%elt_dofs_loc(i, i_elt)) = y(mesh%elt_dofs_loc(i, i_elt)) + ye(i)
      end do
    end do

    ! Ghost contributions to their owners; zero at the prescribed dofs.
    call halo_reverse_add(mesh, y, 3_ik)
    where (bc_mask) y = 0.0_rk

  end subroutine ebe_matvec

  !> Assembled diagonal of [K] (v2 assemble_diagonals, "general"), complete on the owned and,
  !> after a forward exchange, on the ghost dofs.
  !>
  !> Arguments:
  !>   mesh   (in)  type(mesh_t)       rank-local mesh
  !>   k_elt  (in)  real(rk)(:, :, :)  element matrices
  !>   diag   (out) real(rk)(:)        (3 num_nodes_loc) assembled diagonal
  !> Locals:
  !>   i_elt  integer(ik)             element
  !>   i      integer(ik)             element dof
  !>   dofs   integer(ik)(:)          local dofs of the element
  subroutine ebe_diagonal(mesh, k_elt, diag)

    type(mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: k_elt(:, :, :)
    real(rk), intent(out) :: diag(:)

    integer(ik) :: i_elt, i
    integer(ik) :: dofs(3*mesh%num_elt_nodes)

    ! Scatter-add of the element diagonals, then complete owners and ghosts.
    diag = 0.0_rk
    do i_elt = 1_ik, mesh%num_elts
      dofs = elt_dof_list(mesh, i_elt)
      do i = 1_ik, 3_ik*mesh%num_elt_nodes
        diag(dofs(i)) = diag(dofs(i)) + k_elt(i, i, i_elt)
      end do
    end do
    call halo_reverse_add(mesh, diag, 3_ik)
    call halo_forward(mesh, diag, 3_ik)

  end subroutine ebe_diagonal

  !> Assembled 3 x 3 diagonal blocks of [K], one per node (addendum_block_jacobi.md §2): the
  !> block of node a holds the rows and columns of its three dofs, (K_aa)_kl = sum over the
  !> elements of k_elt(3 (i - 1) + k, 3 (i - 1) + l) for the element's node i equal to a;
  !> complete on the owned and, after a forward exchange, on the ghost nodes.
  !>
  !> Arguments:
  !>   mesh    (in)  type(mesh_t)       rank-local mesh
  !>   k_elt   (in)  real(rk)(:, :, :)  element matrices
  !>   blocks  (out) real(rk)(:)        (9 num_nodes_loc) the blocks, node by node, each by
  !>                                    columns: entry 9 (a - 1) + 3 (l - 1) + k is (K_aa)_kl
  !> Locals:
  !>   i_elt   integer(ik)             element
  !>   i       integer(ik)             element node
  !>   k, l    integer(ik)             components (row, column)
  !>   base_e  integer(ik)             element dof before node i's, 3 (i - 1)
  !>   base_n  integer(ik)             entry before node a's block, 9 (a - 1)
  subroutine ebe_block_diagonal(mesh, k_elt, blocks)

    type(mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: k_elt(:, :, :)
    real(rk), intent(out) :: blocks(:)

    integer(ik) :: i_elt, i, k, l, base_e, base_n

    ! Scatter-add of the element's nodal blocks, then complete owners and ghosts.
    blocks = 0.0_rk
    do i_elt = 1_ik, mesh%num_elts
      do i = 1_ik, mesh%num_elt_nodes
        base_e = 3_ik*(i - 1_ik)
        base_n = 9_ik*(mesh%elt_nodes_loc(i, i_elt) - 1_ik)
        do l = 1_ik, 3_ik
          do k = 1_ik, 3_ik
            blocks(base_n + 3_ik*(l - 1_ik) + k) = blocks(base_n + 3_ik*(l - 1_ik) + k) + &
              & k_elt(base_e + k, base_e + l, i_elt)
          end do
        end do
      end do
    end do
    call halo_reverse_add(mesh, blocks, 9_ik)
    call halo_forward(mesh, blocks, 9_ik)

  end subroutine ebe_block_diagonal

  !> Assembled vector of element vectors, complete on the owned and the ghost dofs.
  !>
  !> Arguments:
  !>   mesh   (in)  type(mesh_t)     rank-local mesh
  !>   f_elt  (in)  real(rk)(:, :)   (num_dofs_elt, num_elts) element vectors
  !>   f      (out) real(rk)(:)      (3 num_nodes_loc) assembled vector
  !> Locals:
  !>   i_elt  integer(ik)            element
  !>   dofs   integer(ik)(:)         local dofs of the element
  subroutine ebe_assemble(mesh, f_elt, f)

    type(mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: f_elt(:, :)
    real(rk), intent(out) :: f(:)

    integer(ik) :: i_elt
    integer(ik) :: dofs(3*mesh%num_elt_nodes)

    ! Scatter-add, then complete owners and ghosts.
    f = 0.0_rk
    do i_elt = 1_ik, mesh%num_elts
      dofs = elt_dof_list(mesh, i_elt)
      f(dofs) = f(dofs) + f_elt(:, i_elt)
    end do
    call halo_reverse_add(mesh, f, 3_ik)
    call halo_forward(mesh, f, 3_ik)

  end subroutine ebe_assemble

end module ebe_mod
