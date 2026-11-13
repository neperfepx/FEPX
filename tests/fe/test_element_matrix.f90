! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of element_matrix_mod: quadrature-point kinematics, element matrix and element
!> vector (R14).
!>
!> Source : FEPX (ROADMAP step 4.2).
!>
!> For tet10 and hex20, 20 random curved elements each, [s] at every point a random symmetric
!> positive definite 6 x 6 matrix (anisotropic, different at each point):
!>   - kinematics: {d} = [B]{V} and w = skw(grad v) with grad v(i, j) = sum_a V_a(i) dN_a/dx_j
!>     (explicit 3x3 algebra, 1e-13); a rigid motion v = t + W x gives {d} = 0 and w = (W23,
!>     W13, W12) (1e-12);
!>   - [k_e]: symmetric (1e-12 relative); V^T [k_e] V = sum_q w det J {d}^T [s] {d} (the
!>     element's power from the kinematics, 1e-12); the six rigid-body motions are in its null
!>     space (|[k_e] V_rigid| <= 1e-11 |[k_e]| |V_rigid|); V^T [k_e] V > 0 for random V;
!>   - {f_h,e}: equals [k_e] V when {h} is the {d} of the velocities V (consistency of the two
!>     sums of (R14), 1e-11).
!>
!> Locals:
!>   ework  type(element_work_t)  element work arrays (element_work_alloc; D100)
!>   tet, hex  type(element_t)  descriptors
program test_element_matrix

  use kinds_mod, only: rk, ik
  use element_mod, only: element_t, element_init, element_geometry_qpt, mandel_b_matrix
  use element_matrix_mod, only: element_kinematics_qpt, element_stiffness, element_force_h
  use tet10_mod, only: tet10_gmsh_type, tet10_node_coo_ref
  use hex20_mod, only: hex20_gmsh_type, hex20_node_coo_ref
  use tensor_mod, only: mandel_from_matrix, skew_from_matrix
  use tensor_test_utils_mod, only: seed_random, random_matrix, random_skew, determinant3
  use assert_mod, only: assert_true, assert_close, test_summary
  use element_matrix_mod, only: element_work_t, element_work_alloc

  implicit none

  type(element_work_t) :: ework
  type(element_t) :: tet, hex

  call seed_random(20261009)
  call element_init(tet, tet10_gmsh_type)
  call element_init(hex, hex20_gmsh_type)
  call check_element(tet, 'tet10')
  call check_element(hex, 'hex20')
  call test_summary()

contains

  !> All checks for one element type.
  !>
  !> Arguments:
  !>   elt    (in)  type(element_t)  descriptor
  !>   label  (in)  character        check label prefix
  !> Locals:
  !>   coo       real(rk)(:,:)       element coordinates [mm]
  !>   mat_s     real(rk)(:,:,:)     [s] at the points
  !>   k_e       real(rk)(:,:)       element matrix
  !>   f_e       real(rk)(:)         element vector
  !>   vel       real(rk)(:,:)       velocities
  !>   vec_h     real(rk)(:,:)       {h} at the points
  !>   dshape_s  real(rk)(:,:)       dN/dx
  !>   defrate, spin  real(rk)(6), (3)  kinematics
  !>   power     real(rk)            sum_q w det J d^T s d
  !>   i_elt, q, a, m  integer       loop indices
  subroutine check_element(elt, label)

    type(element_t), intent(in) :: elt
    character(len=*), intent(in) :: label

    real(rk), allocatable :: coo(:, :), mat_s(:, :, :), k_e(:, :), f_e(:), vel(:, :), &
      & vec_h(:, :), dshape_s(:, :), v_rigid(:, :), bmat(:, :)
    real(rk) :: defrate(6), spin(3), velgrad(3, 3), power, det_jac, wmat(3, 3), shift(3)
    real(rk) :: k_norm
    integer :: i_elt, q, a, m

    allocate (mat_s(6, 6, elt%num_qpts), k_e(elt%num_dofs, elt%num_dofs), f_e(elt%num_dofs), &
      & vel(3, elt%num_nodes), vec_h(6, elt%num_qpts), dshape_s(3, elt%num_nodes), &
      & v_rigid(3, elt%num_nodes), bmat(6, elt%num_dofs))

    do i_elt = 1, 20
      coo = random_element(elt)
      do q = 1, int(elt%num_qpts)
        mat_s(:, :, q) = random_spd6()
      end do
      call random_number(vel)
      vel = 2.0_rk*vel - 1.0_rk

      ! Kinematics against explicit algebra.
      power = 0.0_rk
      do q = 1, int(elt%num_qpts)
        call element_work_alloc(ework, elt)
        call element_kinematics_qpt(elt, coo, vel, int(q, ik), defrate, spin, ework)
        call element_geometry_qpt(elt, coo, int(q, ik), dshape_s, det_jac)
        velgrad = matmul(vel, transpose(dshape_s))
        call mandel_b_matrix(dshape_s, bmat)
        call assert_close(defrate, matmul(bmat, reshape(vel, &
          & [elt%num_dofs])), label // ' {d} = [B]{V}', abs_tol=1e-13_rk*maxval(abs(velgrad)))
        call assert_close(defrate, mandel_from_matrix(velgrad), label // ' {d} = {sym grad v}', &
          & abs_tol=1e-13_rk*maxval(abs(velgrad)))
        call assert_close(spin, skew_from_matrix(velgrad), label // ' w = skw grad v', &
          & abs_tol=1e-13_rk*maxval(abs(velgrad)))
        power = power + elt%qpt_weight(q)*det_jac*dot_product(defrate, matmul(mat_s(:, :, q), &
          & defrate))
      end do

      ! Element matrix: symmetry, power, positivity.
      call element_work_alloc(ework, elt)
      call element_stiffness(elt, coo, mat_s, k_e, ework)
      k_norm = maxval(abs(k_e))
      call assert_close(reshape(k_e, [elt%num_dofs**2]), reshape(transpose(k_e), &
        & [elt%num_dofs**2]), label // ' [k_e] symmetric', abs_tol=1e-12_rk*k_norm)
      call assert_close(dot_product(reshape(vel, [elt%num_dofs]), matmul(k_e, reshape(vel, &
        & [elt%num_dofs]))), power, label // ' V^T [k_e] V = power', rel_tol=1e-12_rk)
      call assert_true(power > 0.0_rk, label // ' V^T [k_e] V > 0')

      ! Rigid-body motions (three translations, three rotations) in the null space; their {d}
      ! is zero and their w the rotation.
      do m = 1, 6
        shift = merge(1.0_rk, 0.0_rk, [1, 2, 3] == m)
        wmat = 0.0_rk
        if (m > 3) wmat = random_skew()
        do a = 1, int(elt%num_nodes)
          v_rigid(:, a) = shift + matmul(wmat, coo(:, a))
        end do
        call assert_true(maxval(abs(matmul(k_e, reshape(v_rigid, [elt%num_dofs])))) <= &
          & 1e-11_rk*k_norm*maxval(abs(v_rigid)), label // ' rigid motion in the null space')
        do q = 1, int(elt%num_qpts)
          call element_work_alloc(ework, elt)
          call element_kinematics_qpt(elt, coo, v_rigid, int(q, ik), defrate, spin, ework)
          call assert_close(defrate, [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk], label // &
            & ' rigid motion: {d} = 0', abs_tol=1e-12_rk)
          call assert_close(spin, skew_from_matrix(wmat), label // ' rigid motion: w', &
            & abs_tol=1e-12_rk)
        end do
      end do

      ! Element vector: {h} = {d} of V gives {f_h,e} = [k_e] V.
      do q = 1, int(elt%num_qpts)
        call element_work_alloc(ework, elt)
        call element_kinematics_qpt(elt, coo, vel, int(q, ik), vec_h(:, q), spin, ework)
      end do
      call element_work_alloc(ework, elt)
      call element_force_h(elt, coo, mat_s, vec_h, f_e, ework)
      call assert_close(f_e, matmul(k_e, reshape(vel, [elt%num_dofs])), label // &
        & ' {f_h} = [k_e]{V} for {h} = {d}', abs_tol=1e-11_rk*maxval(abs(f_e)))
    end do

  end subroutine check_element

  !> A random symmetric positive definite 6 x 6 matrix (A A^T + 6 I, entries of A in [-1, 1)).
  !>
  !> Locals:
  !>   amat  real(rk)(6,6)  random factor
  !>   i     integer        diagonal index
  function random_spd6() result(mat)

    real(rk) :: mat(6, 6)

    real(rk) :: amat(6, 6)
    integer :: i

    call random_number(amat)
    amat = 2.0_rk*amat - 1.0_rk
    mat = matmul(amat, transpose(amat))
    do i = 1, 6
      mat(i, i) = mat(i, i) + 6.0_rk
    end do

  end function random_spd6

  !> A random curved element: an affine image of the reference element close to the identity,
  !> with the edge nodes moved by up to 3 % of the reference size.
  !>
  !> Arguments:
  !>   elt  (in)  type(element_t)  descriptor
  !> Locals:
  !>   amat   real(rk)(3,3)  affine map
  !>   move   real(rk)(3)    edge-node displacement
  !>   a, first  integer     node, first edge node
  !>   scale  real(rk)       reference size
  function random_element(elt) result(coo)

    type(element_t), intent(in) :: elt
    real(rk), allocatable :: coo(:, :)

    real(rk) :: amat(3, 3), move(3), scale
    integer :: a, first

    do
      amat = 0.3_rk*random_matrix()
      amat(1, 1) = amat(1, 1) + 1.0_rk
      amat(2, 2) = amat(2, 2) + 1.0_rk
      amat(3, 3) = amat(3, 3) + 1.0_rk
      if (determinant3(amat) > 0.5_rk) exit
    end do
    if (elt%gmsh_type == tet10_gmsh_type) then
      coo = matmul(amat, tet10_node_coo_ref)
      first = 5
      scale = 1.0_rk
    else
      coo = matmul(amat, hex20_node_coo_ref)
      first = 9
      scale = 2.0_rk
    end if
    do a = first, int(elt%num_nodes)
      call random_number(move)
      coo(:, a) = coo(:, a) + 0.03_rk*scale*(2.0_rk*move - 1.0_rk)
    end do

  end function random_element

end program test_element_matrix
