! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Helpers for the property tests of the tensor algebra: reproducible random inputs and explicit
!> 3x3 / fourth-order tensor algebra written independently of src/math/tensor_mod.f90.
!>
!> Source : FEPX test infrastructure (ROADMAP step 1.2). Contains no numerical method of the
!>          solver. The rotation from a unit quaternion and the exponential of a skew matrix
!>          (Rodrigues' formula) are textbook formulas (e.g. Shuster, M.D. (1993), "A survey of
!>          attitude representations", J. Astronaut. Sci. 41(4):439-517); the tests check that
!>          their results are proper rotations. The eigen-decomposition of a symmetric matrix is
!>          the cyclic Jacobi method (Golub & Van Loan, Matrix Computations, 4th ed., 2013,
!>          §8.5.2-8.5.3, Algorithms 8.5.1 and 8.5.3), used only to build exact solutions in tests.
!>
!> Public routines:
!>   seed_random            set the intrinsic generator to a fixed, reproducible state
!>   random_matrix          3x3 matrix with entries uniform in [-1, 1)
!>   random_symmetric       symmetric 3x3 matrix with entries uniform in [-1, 1)
!>   random_skew            skew 3x3 matrix with entries uniform in [-1, 1)
!>   random_rotation        rotation matrix from a random unit quaternion
!>   rotation_from_quaternion  rotation matrix of a unit quaternion (q0, q1, q2, q3)
!>   exp_skew               exp(W) of a skew 3x3 matrix W (Rodrigues' formula)
!>   mandel_pair            Mandel component (1..6) of the tensor index pair (i, j)
!>   tensor4_from_mandel    fourth-order tensor C(3,3,3,3) of a 6x6 Mandel matrix
!>   mandel_from_tensor4    6x6 Mandel matrix of a fourth-order tensor with minor symmetries
!>   rotate_tensor4         C'_ijkl = R_ip R_jq R_kr R_ls C_pqrs
!>   determinant3           determinant of a 3x3 matrix
!>   symmetric_eigen        eigenvalues and orthonormal eigenvectors of a symmetric n x n matrix
module tensor_test_utils_mod

  use kinds_mod, only: rk

  implicit none

  private

  public :: seed_random, random_matrix, random_symmetric, random_skew, random_rotation
  public :: rotation_from_quaternion, exp_skew
  public :: mandel_pair, tensor4_from_mandel, mandel_from_tensor4, rotate_tensor4
  public :: determinant3, symmetric_eigen

contains

  !> Set the intrinsic random number generator to a fixed state, so that every run of a test
  !> uses the same inputs.
  !>
  !> Arguments:
  !>   seed  (in)  integer  base value of the seed
  !> Locals:
  !>   seed_size   integer               number of integers in the generator state
  !>   seed_array  integer(:), allocatable  generator state
  !>   i           integer               loop index
  subroutine seed_random(seed)

    integer, intent(in) :: seed

    integer :: seed_size
    integer, allocatable :: seed_array(:)
    integer :: i

    ! Fill the state with distinct values derived from the seed.
    call random_seed(size=seed_size)
    allocate (seed_array(seed_size))
    do i = 1, seed_size
      seed_array(i) = seed + 7919*i
    end do
    call random_seed(put=seed_array)

  end subroutine seed_random

  !> 3x3 matrix with entries uniform in [-1, 1).
  !>
  !> Arguments:
  !>   mat  (out, result) real(rk)(3,3)  random matrix
  function random_matrix() result(mat)

    real(rk) :: mat(3, 3)

    ! Map [0, 1) to [-1, 1).
    call random_number(mat)
    mat = 2.0_rk*mat - 1.0_rk

  end function random_matrix

  !> Symmetric 3x3 matrix with entries uniform in [-1, 1).
  !>
  !> Arguments:
  !>   mat  (out, result) real(rk)(3,3)  random symmetric matrix
  !> Locals:
  !>   i, j  integer  row and column
  function random_symmetric() result(mat)

    real(rk) :: mat(3, 3)

    integer :: i, j

    ! Draw a full matrix and copy its upper triangle to the lower triangle.
    mat = random_matrix()
    do j = 1, 3
      do i = j + 1, 3
        mat(i, j) = mat(j, i)
      end do
    end do

  end function random_symmetric

  !> Skew 3x3 matrix with off-diagonal entries uniform in [-1, 1).
  !>
  !> Arguments:
  !>   mat  (out, result) real(rk)(3,3)  random skew matrix
  !> Locals:
  !>   i, j  integer  row and column
  function random_skew() result(mat)

    real(rk) :: mat(3, 3)

    integer :: i, j

    ! Draw a full matrix, zero the diagonal, negate the upper triangle into the lower one.
    mat = random_matrix()
    do j = 1, 3
      mat(j, j) = 0.0_rk
      do i = j + 1, 3
        mat(i, j) = -mat(j, i)
      end do
    end do

  end function random_skew

  !> Rotation matrix from a random unit quaternion (normalized from 4 entries uniform in
  !> [-1, 1), redrawn when the norm is below 0.1 to avoid a loss of precision).
  !>
  !> Arguments:
  !>   rot  (out, result) real(rk)(3,3)  rotation matrix
  !> Locals:
  !>   quat   real(rk)(4)  quaternion
  !>   norm   real(rk)     norm of quat before normalization
  function random_rotation() result(rot)

    real(rk) :: rot(3, 3)

    real(rk) :: quat(4)
    real(rk) :: norm

    ! Draw until the quaternion is not too short, then normalize it.
    norm = 0.0_rk
    do while (norm < 0.1_rk)
      call random_number(quat)
      quat = 2.0_rk*quat - 1.0_rk
      norm = sqrt(sum(quat**2))
    end do
    rot = rotation_from_quaternion(quat/norm)

  end function random_rotation

  !> Rotation matrix of a unit quaternion q = (q0, q1, q2, q3), q0 the scalar part: the
  !> matrix of the map v -> q v q*.
  !>
  !> Arguments:
  !>   quat  (in)  real(rk)(4)    unit quaternion
  !>   rot   (out, result) real(rk)(3,3)  rotation matrix
  !> Locals:
  !>   q0, q1, q2, q3  real(rk)  quaternion components
  pure function rotation_from_quaternion(quat) result(rot)

    real(rk), intent(in) :: quat(4)
    real(rk) :: rot(3, 3)

    real(rk) :: q0, q1, q2, q3

    q0 = quat(1)
    q1 = quat(2)
    q2 = quat(3)
    q3 = quat(4)

    ! Rows of the rotation matrix.
    rot(1, :) = [1.0_rk - 2.0_rk*(q2*q2 + q3*q3), 2.0_rk*(q1*q2 - q0*q3), &
      & 2.0_rk*(q1*q3 + q0*q2)]
    rot(2, :) = [2.0_rk*(q1*q2 + q0*q3), 1.0_rk - 2.0_rk*(q1*q1 + q3*q3), &
      & 2.0_rk*(q2*q3 - q0*q1)]
    rot(3, :) = [2.0_rk*(q1*q3 - q0*q2), 2.0_rk*(q2*q3 + q0*q1), &
      & 1.0_rk - 2.0_rk*(q1*q1 + q2*q2)]

  end function rotation_from_quaternion

  !> Exponential of a skew 3x3 matrix W (Rodrigues' formula):
  !>   exp(W) = I + (sin t / t) W + ((1 - cos t) / t^2) W W,  t = sqrt(W:W / 2).
  !>
  !> Arguments:
  !>   skw_mat  (in)  real(rk)(3,3)  skew matrix W
  !>   rot      (out, result) real(rk)(3,3)  rotation exp(W)
  !> Locals:
  !>   angle  real(rk)  rotation angle t [rad]
  !>   i      integer   loop index over the diagonal
  pure function exp_skew(skw_mat) result(rot)

    real(rk), intent(in) :: skw_mat(3, 3)
    real(rk) :: rot(3, 3)

    real(rk) :: angle
    integer :: i

    ! Angle of the rotation.
    angle = sqrt(0.5_rk*sum(skw_mat*skw_mat))

    ! Identity; the formula's terms only when the angle is not zero.
    rot = 0.0_rk
    do i = 1, 3
      rot(i, i) = 1.0_rk
    end do
    if (angle > 0.0_rk) then
      rot = rot + (sin(angle)/angle)*skw_mat &
        & + ((1.0_rk - cos(angle))/angle**2)*matmul(skw_mat, skw_mat)
    end if

  end function exp_skew

  !> Mandel component (1..6, order 11 22 33 23 13 12) of the tensor index pair (i, j).
  !>
  !> Arguments:
  !>   i, j  (in)  integer  tensor indices (1..3)
  !>   p     (out, result) integer  Mandel component
  pure function mandel_pair(i, j) result(p)

    integer, intent(in) :: i, j
    integer :: p

    ! Diagonal pairs map to 1..3; off-diagonal pairs by the missing index k = 6 - i - j:
    ! (2,3) -> 4, (1,3) -> 5, (1,2) -> 6, i.e. p = 3 + k.
    if (i == j) then
      p = i
    else
      p = 3 + (6 - i - j)
    end if

  end function mandel_pair

  !> Fourth-order tensor C(3,3,3,3) of a 6x6 Mandel matrix: C_ijkl = M(p,q) / (s_p s_q), with
  !> p = (i,j), q = (k,l) and s = 1 (normal) or sqrt(2) (shear) (Mehrabadi & Cowin 1990).
  !>
  !> Arguments:
  !>   mat6     (in)  real(rk)(6,6)        Mandel matrix
  !>   tensor4  (out, result) real(rk)(3,3,3,3)  fourth-order tensor with minor symmetries
  !> Locals:
  !>   i, j, k, l  integer   tensor indices
  !>   p, q        integer   Mandel components
  !>   scale       real(rk)(6)  s_p
  pure function tensor4_from_mandel(mat6) result(tensor4)

    real(rk), intent(in) :: mat6(6, 6)
    real(rk) :: tensor4(3, 3, 3, 3)

    integer :: i, j, k, l
    integer :: p, q
    real(rk) :: scale(6)

    scale = [1.0_rk, 1.0_rk, 1.0_rk, sqrt(2.0_rk), sqrt(2.0_rk), sqrt(2.0_rk)]
    do l = 1, 3
      do k = 1, 3
        do j = 1, 3
          do i = 1, 3
            p = mandel_pair(i, j)
            q = mandel_pair(k, l)
            tensor4(i, j, k, l) = mat6(p, q)/(scale(p)*scale(q))
          end do
        end do
      end do
    end do

  end function tensor4_from_mandel

  !> 6x6 Mandel matrix of a fourth-order tensor with minor symmetries: M(p,q) = s_p s_q C_ijkl.
  !>
  !> Arguments:
  !>   tensor4  (in)  real(rk)(3,3,3,3)  fourth-order tensor
  !>   mat6     (out, result) real(rk)(6,6)  Mandel matrix
  !> Locals:
  !>   ii, jj  integer(6)  tensor indices of each Mandel component (11 22 33 23 13 12)
  !>   p, q    integer     Mandel components
  !>   scale   real(rk)(6) s_p
  pure function mandel_from_tensor4(tensor4) result(mat6)

    real(rk), intent(in) :: tensor4(3, 3, 3, 3)
    real(rk) :: mat6(6, 6)

    integer, parameter :: ii(6) = [1, 2, 3, 2, 1, 1]
    integer, parameter :: jj(6) = [1, 2, 3, 3, 3, 2]
    integer :: p, q
    real(rk) :: scale(6)

    scale = [1.0_rk, 1.0_rk, 1.0_rk, sqrt(2.0_rk), sqrt(2.0_rk), sqrt(2.0_rk)]
    do q = 1, 6
      do p = 1, 6
        mat6(p, q) = scale(p)*scale(q)*tensor4(ii(p), jj(p), ii(q), jj(q))
      end do
    end do

  end function mandel_from_tensor4

  !> Rotation of a fourth-order tensor, C'_ijkl = R_ip R_jq R_kr R_ls C_pqrs.
  !>
  !> Arguments:
  !>   rot      (in)  real(rk)(3,3)      rotation R
  !>   tensor4  (in)  real(rk)(3,3,3,3)  tensor C
  !>   rotated  (out, result) real(rk)(3,3,3,3)  rotated tensor C'
  !> Locals:
  !>   i, j, k, l, p, q, r, s  integer  tensor indices
  pure function rotate_tensor4(rot, tensor4) result(rotated)

    real(rk), intent(in) :: rot(3, 3)
    real(rk), intent(in) :: tensor4(3, 3, 3, 3)
    real(rk) :: rotated(3, 3, 3, 3)

    integer :: i, j, k, l, p, q, r, s

    rotated = 0.0_rk
    do l = 1, 3
      do k = 1, 3
        do j = 1, 3
          do i = 1, 3
            do s = 1, 3
              do r = 1, 3
                do q = 1, 3
                  do p = 1, 3
                    rotated(i, j, k, l) = rotated(i, j, k, l) &
                      & + rot(i, p)*rot(j, q)*rot(k, r)*rot(l, s)*tensor4(p, q, r, s)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end function rotate_tensor4

  !> Determinant of a 3x3 matrix (cofactor expansion along the first row).
  !>
  !> Arguments:
  !>   mat  (in)  real(rk)(3,3)  matrix
  !>   det  (out, result) real(rk)  determinant
  pure function determinant3(mat) result(det)

    real(rk), intent(in) :: mat(3, 3)
    real(rk) :: det

    det = mat(1, 1)*(mat(2, 2)*mat(3, 3) - mat(2, 3)*mat(3, 2)) &
      & - mat(1, 2)*(mat(2, 1)*mat(3, 3) - mat(2, 3)*mat(3, 1)) &
      & + mat(1, 3)*(mat(2, 1)*mat(3, 2) - mat(2, 2)*mat(3, 1))

  end function determinant3

  !> Eigenvalues and orthonormal eigenvectors of a symmetric n x n matrix A: A V = V diag(values),
  !> by the cyclic Jacobi method (Golub & Van Loan 2013, Algorithm 8.5.3, with the 2x2 symmetric
  !> Schur decomposition of Algorithm 8.5.1), sweeps until the off-diagonal norm is below
  !> 1e-15 of the Frobenius norm (at most 50 sweeps).
  !>
  !> Arguments:
  !>   mat      (in)  real(rk)(:, :)  symmetric matrix A (n x n)
  !>   values   (out) real(rk)(:)     (n) eigenvalues (unordered)
  !>   vectors  (out) real(rk)(:, :)  (n, n) eigenvectors, one per column
  !> Locals:
  !>   a        real(rk)(:, :)        the matrix being diagonalized, J^T A J
  !>   col_p, col_q  real(rk)(:)      columns (or rows) p and q before a rotation
  !>   norm_f   real(rk)              Frobenius norm of A
  !>   off      real(rk)              off-diagonal norm
  !>   tau, t, c, s  real(rk)         rotation parameters (Algorithm 8.5.1)
  !>   n, p, q, k, sweep  integer     size, rotation plane, index, sweep
  subroutine symmetric_eigen(mat, values, vectors)

    real(rk), intent(in) :: mat(:, :)
    real(rk), intent(out) :: values(:)
    real(rk), intent(out) :: vectors(:, :)

    real(rk), allocatable :: a(:, :), col_p(:), col_q(:)
    real(rk) :: norm_f, off, tau, t, c, s
    integer :: n, p, q, k, sweep

    ! Start: A itself and V = I.
    n = size(mat, 1)
    a = mat
    vectors = 0.0_rk
    do k = 1, n
      vectors(k, k) = 1.0_rk
    end do
    norm_f = sqrt(sum(a**2))

    do sweep = 1, 50
      ! Off-diagonal norm: stop when negligible.
      off = sqrt(max(sum(a**2) - sum([(a(k, k)**2, k = 1, n)]), 0.0_rk))
      if (off <= 1.0e-15_rk*norm_f) exit
      ! One sweep over every plane (p, q).
      do p = 1, n - 1
        do q = p + 1, n
          if (abs(a(p, q)) > 0.0_rk) then
            ! 2x2 symmetric Schur decomposition (Algorithm 8.5.1).
            tau = (a(q, q) - a(p, p))/(2.0_rk*a(p, q))
            if (tau >= 0.0_rk) then
              t = 1.0_rk/(tau + sqrt(1.0_rk + tau**2))
            else
              t = -1.0_rk/(-tau + sqrt(1.0_rk + tau**2))
            end if
            c = 1.0_rk/sqrt(1.0_rk + t**2)
            s = t*c
            ! A <- A J (columns p, q), then J^T A (rows p, q); V <- V J.
            col_p = a(:, p)
            col_q = a(:, q)
            a(:, p) = c*col_p - s*col_q
            a(:, q) = s*col_p + c*col_q
            col_p = a(p, :)
            col_q = a(q, :)
            a(p, :) = c*col_p - s*col_q
            a(q, :) = s*col_p + c*col_q
            col_p = vectors(:, p)
            col_q = vectors(:, q)
            vectors(:, p) = c*col_p - s*col_q
            vectors(:, q) = s*col_p + c*col_q
          end if
        end do
      end do
    end do

    ! The eigenvalues are the diagonal.
    values = [(a(k, k), k = 1, n)]

  end subroutine symmetric_eigen

end module tensor_test_utils_mod
