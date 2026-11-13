! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of linalg_mod (LDL^T and normal-equations solves of small dense systems).
!>
!> Source : FEPX (ROADMAP steps 4.1 and 9.15: ldlt_factor and ldlt_substitute in place of
!>          ldlt_solve, D100).
!>
!> Checks (sizes 1 to 8, 200 random systems each):
!>   - ldlt_factor + ldlt_substitute (ldlt_solve_copy below) on symmetric positive definite
!>     A = B B^T + n I: recovers a known solution (relative 1e-12), reads only the lower triangle
!>     (garbage in the upper triangle gives the same x bit for bit), status 0; ldlt_factor leaves
!>     the strict upper triangle unchanged; one factorization solves a second right-hand side
!>     (relative 1e-12); a matrix with a nonpositive pivot (negative diagonal; zero matrix;
!>     indefinite [[1, 2], [2, 1]]) gives status 1 (and x = 0 from ldlt_solve_copy);
!>   - normal_equations_solve on nonsingular nonsymmetric A = n I + B (B random, entries in
!>     [-1, 1)) and on (1/dt) I + [W^] (the elastic point-solve matrix, antisymmetric [W^]):
!>     recovers a known solution (relative 1e-10: the normal equations square the condition
!>     number); an exactly singular A ([[1, 1], [1, 1]]: zero pivot) gives status 1 and x = 0.
!>     (A nearly singular A can give a tiny positive pivot and is not detected, as in v2.)
!>   - lu_factor / lu_solve (step 6.7) on general A (entries in [-1, 1), no diagonal dominance):
!>     backward error |A x - b| <= 1e-13 |A| |x| (partial pivoting is backward stable);
!>     det(A B) = det A det B (relative 1e-11); det of a 3 x 3 against the explicit formula;
!>     a zero leading entry needs the row swap ([[0, 1], [1, 0]]: x and det = -1 exact);
!>     several right-hand sides from one factorization; exactly singular A (equal rows; zero
!>     matrix) gives status 1 and det = 0.
!>
!> Locals:
!>   n, i_trial, status  integer        size, trial, solve status
!>   a, bmat, a_garbage  real(rk)(:,:)  matrices
!>   x_true, b, x, x2    real(rk)(:)    solution, right-hand side, computed solutions
!>   fact, ata           real(rk)(:,:)  LDL^T factors; A^T A work array of the normal equations
program test_linalg

  use kinds_mod, only: rk
  use linalg_mod, only: ldlt_factor, ldlt_substitute, normal_equations_solve, lu_factor, lu_solve
  use tensor_mod, only: spin_commutator
  use tensor_test_utils_mod, only: seed_random, determinant3
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary

  implicit none

  integer :: n, i_trial, status, i, j
  real(rk), allocatable :: a(:, :), bmat(:, :), a_garbage(:, :), fact(:, :), ata(:, :)
  real(rk), allocatable :: x_true(:), b(:), x(:), x2(:)
  real(rk) :: w(3)
  real(rk), allocatable :: lu(:, :), lu2(:, :)
  integer, allocatable :: perm(:), perm2(:)
  real(rk) :: det, det2, det_prod

  call seed_random(20261009)

  do n = 1, 8
    allocate (a(n, n), bmat(n, n), a_garbage(n, n), x_true(n), b(n), x(n), x2(n), fact(n, n), &
      & ata(n, n))
    do i_trial = 1, 200
      ! Symmetric positive definite system with a known solution.
      call random_number(bmat)
      bmat = 2.0_rk*bmat - 1.0_rk
      a = matmul(bmat, transpose(bmat))
      do i = 1, n
        a(i, i) = a(i, i) + real(n, rk)
      end do
      call random_number(x_true)
      x_true = 2.0_rk*x_true - 1.0_rk
      b = matmul(a, x_true)
      call ldlt_solve_copy(a, b, x, status)
      call assert_equal(status, 0, 'ldlt status 0')
      call assert_close(x, x_true, 'ldlt solution', abs_tol=1e-12_rk*maxval(abs(x_true)))

      ! Only the lower triangle is read.
      a_garbage = a
      do j = 2, n
        do i = 1, j - 1
          a_garbage(i, j) = 1.0e30_rk
        end do
      end do
      call ldlt_solve_copy(a_garbage, b, x2, status)
      call assert_close(x2, x, 'ldlt reads the lower triangle only')

      ! The factorization leaves the strict upper triangle; it solves a second right-hand side.
      fact = a_garbage
      call ldlt_factor(fact, status)
      do j = 2, n
        call assert_close(fact(1:j - 1, j), a_garbage(1:j - 1, j), &
          & 'ldlt_factor leaves the upper triangle', abs_tol=0.0_rk)
      end do
      call random_number(x2)
      b = matmul(a, x2)
      call ldlt_substitute(fact, b)
      call assert_close(b, x2, 'ldlt second right-hand side', abs_tol=1e-12_rk*maxval(abs(x2)))
      b = matmul(a, x_true)

      ! Nonsymmetric nonsingular system through the normal equations.
      call random_number(bmat)
      a = 2.0_rk*bmat - 1.0_rk
      do i = 1, n
        a(i, i) = a(i, i) + real(n, rk)
      end do
      b = matmul(a, x_true)
      call normal_equations_solve(a, b, x, ata, status)
      call assert_equal(status, 0, 'normal equations status 0')
      call assert_close(x, x_true, 'normal equations solution', &
        & abs_tol=1e-10_rk*maxval(abs(x_true)))
    end do

    deallocate (a, bmat, a_garbage, x_true, b, x, x2, fact, ata)
  end do

  ! The elastic point-solve matrix (1/dt) I + [W^] (dt = 1e-2, |w| ~ 1).
  allocate (a(6, 6), x_true(6), b(6), x(6), ata(6, 6))
  do i_trial = 1, 200
    call random_number(w)
    a = spin_commutator(2.0_rk*w - 1.0_rk)
    do i = 1, 6
      a(i, i) = a(i, i) + 100.0_rk
    end do
    call random_number(x_true)
    b = matmul(a, x_true)
    call normal_equations_solve(a, b, x, ata, status)
    call assert_equal(status, 0, 'elastic matrix status 0')
    call assert_close(x, x_true, 'elastic matrix solution', abs_tol=1e-12_rk)
  end do
  deallocate (a, x_true, b, x, ata)

  ! Exactly singular nonsymmetric system (equal rows of ones: A^T A has an exactly zero pivot).
  allocate (a(2, 2), b(2))
  a = 1.0_rk
  b = [1.0_rk, 1.0_rk]
  allocate (x(2), ata(2, 2))
  call normal_equations_solve(a, b, x, ata, status)
  call assert_equal(status, 1, 'normal equations singular status 1')
  call assert_true(all(abs(x) < tiny(1.0_rk)), 'normal equations failure returns x = 0')
  deallocate (a, b, x, ata)

  ! Nonpositive pivots: negative diagonal, zero matrix, indefinite.
  allocate (x(2))
  call ldlt_solve_copy(reshape( &
    & [-1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk], [2, 2]), [1.0_rk, 1.0_rk], x, &
    & status)
  call assert_equal(status, 1, 'ldlt negative pivot status 1')
  call assert_true(all(abs(x) < tiny(1.0_rk)), 'ldlt failure returns x = 0')
  call ldlt_solve_copy(reshape( &
    & [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk], [2, 2]), [1.0_rk, 1.0_rk], x, &
    & status)
  call assert_equal(status, 1, 'ldlt zero matrix status 1')
  call ldlt_solve_copy(reshape( &
    & [1.0_rk, 2.0_rk, 2.0_rk, 1.0_rk], [2, 2]), [1.0_rk, 1.0_rk], x, &
    & status)
  call assert_equal(status, 1, 'ldlt indefinite status 1')
  deallocate (x)

  ! LU with partial pivoting on general matrices.
  do n = 1, 8
    allocate (a(n, n), bmat(n, n), lu(n, n), lu2(n, n), perm(n), perm2(n), x_true(n), b(n), &
      & x(n))
    do i_trial = 1, 200
      call random_number(a)
      a = 2.0_rk*a - 1.0_rk
      call random_number(bmat)
      bmat = 2.0_rk*bmat - 1.0_rk
      call lu_factor(a, lu, perm, det, status)
      call assert_equal(status, 0, 'lu status 0')
      ! Two right-hand sides from one factorization; backward error.
      do j = 1, 2
        call random_number(x_true)
        x_true = 2.0_rk*x_true - 1.0_rk
        b = matmul(a, x_true)
        call lu_solve(lu, perm, b, x)
        call assert_true(norm2(matmul(a, x) - b) <= 1e-13_rk*norm2(a)*norm2(x), &
          & 'lu backward error')
      end do
      ! det(A B) = det A det B.
      call lu_factor(bmat, lu2, perm2, det2, status)
      det_prod = det*det2
      call lu_factor(matmul(a, bmat), lu2, perm2, det2, status)
      call assert_close(det2, det_prod, 'lu det of a product', &
        & abs_tol=1e-11_rk*max(abs(det_prod), 1e-3_rk))
      if (n == 3) call assert_close(det, determinant3(a), 'lu det 3x3', abs_tol=1e-14_rk)
    end do
    deallocate (a, bmat, lu, lu2, perm, perm2, x_true, b, x)
  end do

  ! A zero leading entry needs the row swap; exact results.
  allocate (lu(2, 2), perm(2), x(2))
  call lu_factor(reshape([0.0_rk, 1.0_rk, 1.0_rk, 0.0_rk], [2, 2]), lu, perm, det, status)
  call assert_equal(status, 0, 'lu swap status 0')
  call assert_close(det, -1.0_rk, 'lu swap det', abs_tol=0.0_rk)
  call lu_solve(lu, perm, [2.0_rk, 3.0_rk], x)
  call assert_close(x, [3.0_rk, 2.0_rk], 'lu swap solution', abs_tol=0.0_rk)

  ! Exactly singular: equal rows, zero matrix.
  call lu_factor(reshape([1.0_rk, 1.0_rk, 2.0_rk, 2.0_rk], [2, 2]), lu, perm, det, status)
  call assert_equal(status, 1, 'lu equal rows status 1')
  call assert_true(.not. (abs(det) > 0.0_rk), 'lu singular det 0')
  call lu_factor(reshape([0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk], [2, 2]), lu, perm, det, status)
  call assert_equal(status, 1, 'lu zero matrix status 1')
  deallocate (lu, perm, x)

  call test_summary()

contains

  !> A x = b by ldlt_factor and ldlt_substitute on a copy of A; x = 0 when the factorization
  !> fails (the behavior of the former ldlt_solve, which the checks above were written for).
  !>
  !> Arguments:
  !>   a_in    (in)  real(rk)(:, :)  matrix (lower triangle read)
  !>   b_in    (in)  real(rk)(:)     right-hand side
  !>   x_out   (out) real(rk)(:)     solution
  !>   ios     (out) integer         ldlt_factor status
  !> Locals:
  !>   work    real(rk)(:, :)        the copy, then its factors
  subroutine ldlt_solve_copy(a_in, b_in, x_out, ios)

    real(rk), intent(in) :: a_in(:, :)
    real(rk), intent(in) :: b_in(:)
    real(rk), intent(out) :: x_out(:)
    integer, intent(out) :: ios

    real(rk), allocatable :: work(:, :)

    ! Factorize the copy; solve in place, or return zero.
    work = a_in
    call ldlt_factor(work, ios)
    x_out = 0.0_rk
    if (ios /= 0) return
    x_out = b_in
    call ldlt_substitute(work, x_out)

  end subroutine ldlt_solve_copy

end program test_linalg
