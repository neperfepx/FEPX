! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Small dense linear systems (one material point at a time): LDL^T solution of a symmetric
!> positive definite system, and the solution of a general nonsingular system through its normal
!> equations, as v2 solves its per-point Newton systems (no work arrays of their own: the caller
!> provides them at their exact size, ROADMAP D100); LU factorization with partial pivoting
!> (determinant and several right-hand sides of one general matrix; the yield-surface vertices,
!> step 6.7).
!>
!> Source : port of v2 solve_evp_vpstress_mod3.f90::solvit (LDL^T factorization without
!>          pivoting of a symmetric positive definite matrix, then L y = b, D z = y,
!>          L^T x = z; v2 unrolls it for 5 x 5, here for any size) and
!>          solve_evp_stress_mod3.f90::symmetrize_jac (normal equations A^T A x = A^T b of a
!>          nonsymmetric system, which v2 then solves with solvit). LU: Golub, G.H. and
!>          Van Loan, C.F., Matrix Computations, 4th ed. (2013), Algorithm 3.4.1 (Gaussian
!>          elimination with partial pivoting) and §3.1 (forward and back substitution)
!>          (ROADMAP D57).
!>
!> Public routines:
!>   ldlt_factor             A = L D L^T in place, A symmetric positive definite (lower triangle)
!>   ldlt_substitute         solve A x = b in place from the factors of ldlt_factor
!>   normal_equations_solve  solve A x = b, A general nonsingular, through A^T A x = A^T b
!>                           (A^T A in the caller's work array)
!>   lu_factor               P A = L U with partial pivoting, and det A
!>   lu_solve                solve A x = b from the factors of lu_factor
module linalg_mod

  use kinds_mod, only: rk

  implicit none

  private

  public :: ldlt_factor, ldlt_substitute, normal_equations_solve, lu_factor, lu_solve

contains

  !> LDL^T factorization of a symmetric positive definite A in place (no pivoting; only the lower
  !> triangle of A is read and overwritten: L strictly below the diagonal, D on it; the strict
  !> upper triangle is left unchanged). A pivot d_j <= 0 (A not positive definite, or exactly
  !> singular) stops the factorization with status 1 (the factor is then incomplete). A nearly
  !> singular A can give a tiny positive pivot and is not detected (as in v2, whose
  !> check_diagonals_evp tests only for a zero diagonal); the callers' iterations detect the
  !> resulting non-convergence. The caller keeps a copy of A if it needs it (no work array here,
  !> ROADMAP D100; the split into factorization and substitution is that of LAPACK dpotrf and
  !> dpotrs).
  !>
  !> Source : port of v2 solve_evp_vpstress_mod3.f90::solvit (the same recurrences, for n x n):
  !>            d_j = a_jj - sum_{i<j} l_ji^2 d_i,
  !>            l_kj = (a_kj - sum_{i<j} l_ki l_ji d_i) / d_j   (k > j);
  !>          the in-place storage of L and D: Golub & Van Loan (2013), §4.1.2 (Algorithm 4.1.2,
  !>          LDL^T).
  !>
  !> Arguments:
  !>   a       (inout) real(rk)(:, :)  n x n symmetric positive definite matrix in (lower triangle
  !>                                   read); out: L strictly below the diagonal, D on it
  !>   status  (out)   integer         0: factorized; 1: nonpositive pivot
  !> Locals:
  !>   n       integer                 size
  !>   i, j, k integer                 indices
  !>   lji_di  real(rk)                l_ji d_i for the current column j and row i < j
  !>   pivot   real(rk)                sum_{i<j} l_ji^2 d_i, accumulated for the pivot d_j
  pure subroutine ldlt_factor(a, status)

    real(rk), intent(inout) :: a(:, :)
    integer, intent(out) :: status

    integer :: n, i, j, k
    real(rk) :: lji_di, pivot

    n = size(a, 1)
    status = 0

    ! Factorization, column by column (v2's "j = 1 ... 5" steps). Row j already holds l_ji
    ! (i < j, set by the earlier columns); l_ji d_i is formed where it is used (no work vector),
    ! the products and the order of the sums being those of v2's recurrences.
    do j = 1, n
      ! The pivot d_j = a_jj - sum_{i<j} (l_ji d_i) l_ji.
      pivot = 0.0_rk
      do i = 1, j - 1
        lji_di = a(j, i)*a(i, i)
        pivot = pivot + a(j, i)*lji_di
      end do
      a(j, j) = a(j, j) - pivot
      if (.not. (a(j, j) > 0.0_rk)) then
        status = 1
        return
      end if
      ! Column j of L below the diagonal: l_kj = (a_kj - sum_{i<j} l_ki (l_ji d_i)) / d_j.
      do k = j + 1, n
        pivot = 0.0_rk
        do i = 1, j - 1
          lji_di = a(j, i)*a(i, i)
          pivot = pivot + a(k, i)*lji_di
        end do
        a(k, j) = (a(k, j) - pivot)/a(j, j)
      end do
    end do

  end subroutine ldlt_factor

  !> Solve A x = b in place from the LDL^T factors of ldlt_factor: L y = b (unit lower), D z = y,
  !> L^T x = z.
  !>
  !> Source : port of v2 solve_evp_vpstress_mod3.f90::solvit (its substitutions, for n x n).
  !>
  !> Arguments:
  !>   fact  (in)    real(rk)(:, :)  n x n factors of ldlt_factor (status 0)
  !>   x     (inout) real(rk)(:)     right-hand side b in, solution x out (n)
  !> Locals:
  !>   n     integer                 size
  !>   k     integer                 row
  pure subroutine ldlt_substitute(fact, x)

    real(rk), intent(in) :: fact(:, :)
    real(rk), intent(inout) :: x(:)

    integer :: n, k

    n = size(fact, 1)

    ! L y = b (unit lower triangular).
    do k = 2, n
      x(k) = x(k) - dot_product(fact(k, 1:k - 1), x(1:k - 1))
    end do

    ! D z = y.
    do k = 1, n
      x(k) = x(k)/fact(k, k)
    end do

    ! L^T x = z.
    do k = n - 1, 1, -1
      x(k) = x(k) - dot_product(fact(k + 1:n, k), x(k + 1:n))
    end do

  end subroutine ldlt_substitute

  !> Solve A x = b for a general nonsingular A through the normal equations
  !> A^T A x = A^T b, whose matrix is symmetric positive definite (ldlt_factor, ldlt_substitute).
  !> A^T A is formed in the caller's array ata (no work array here, ROADMAP D100).
  !>
  !> Source : port of v2 solve_evp_stress_mod3.f90::symmetrize_jac (A^T A and A^T b), followed
  !>          by solvit, as in v2 solve_evp_stress_mod2.f90.
  !>
  !> Arguments:
  !>   a       (in)  real(rk)(:, :)   n x n nonsingular matrix
  !>   b       (in)  real(rk)(:)      right-hand side (n)
  !>   x       (out) real(rk)(:)      solution (n); 0 when status is 1
  !>   ata     (out) real(rk)(:, :)   n x n work array: A^T A, then its LDL^T factors
  !>   status  (out) integer          0: solved; 1: A^T A not positive definite (A singular)
  pure subroutine normal_equations_solve(a, b, x, ata, status)

    real(rk), intent(in) :: a(:, :)
    real(rk), intent(in) :: b(:)
    real(rk), intent(out) :: x(:)
    real(rk), intent(out) :: ata(:, :)
    integer, intent(out) :: status

    ! A^T A and A^T b.
    ata = matmul(transpose(a), a)
    x = matmul(transpose(a), b)

    ! A^T A x = A^T b by LDL^T; x = 0 when A^T A is not positive definite.
    call ldlt_factor(ata, status)
    if (status /= 0) then
      x = 0.0_rk
      return
    end if
    call ldlt_substitute(ata, x)

  end subroutine normal_equations_solve

  !> LU factorization with partial pivoting, P A = L U (L unit lower triangular, U upper
  !> triangular, both stored in lu), and the determinant det A = sign(P) prod_k u_kk. An exactly
  !> zero pivot column (A singular) stops the factorization with status 1; det is then 0. A
  !> nearly singular A is not detected here: callers test det (e.g. against Hadamard's bound).
  !>
  !> Source : Golub & Van Loan (2013), Algorithm 3.4.1 (Gaussian elimination with partial
  !>          pivoting): for k = 1..n-1, the row p >= k with the largest |a_pk| is swapped with
  !>          row k (perm records it), l_ik = a_ik / a_kk and a_ij -= l_ik a_kj for i, j > k.
  !>
  !> Arguments:
  !>   a       (in)  real(rk)(:, :)   n x n matrix
  !>   lu      (out) real(rk)(:, :)   n x n factors: L below the diagonal (unit diagonal not
  !>                                  stored), U on and above it
  !>   perm    (out) integer(:)       row k of P A is row perm(k) of A (n)
  !>   det     (out) real(rk)         det A
  !>   status  (out) integer          0: factorized; 1: an exactly zero pivot (A singular)
  !> Locals:
  !>   n       integer                size
  !>   k, i, p integer                elimination step, row, pivot row
  !>   row     real(rk)(n)            a row being swapped
  !>   tmp     integer                a permutation entry being swapped
  pure subroutine lu_factor(a, lu, perm, det, status)

    real(rk), intent(in) :: a(:, :)
    real(rk), intent(out) :: lu(:, :)
    integer, intent(out) :: perm(:)
    real(rk), intent(out) :: det
    integer, intent(out) :: status

    integer :: n, k, i, p, tmp
    real(rk) :: row(size(a, 1))

    n = size(a, 1)
    lu = a
    status = 0
    det = 1.0_rk
    do k = 1, n
      perm(k) = k
    end do

    do k = 1, n
      ! Pivot: the row at or below k with the largest |a_ik| (first on ties).
      p = k - 1 + maxloc(abs(lu(k:n, k)), dim=1)
      if (.not. (abs(lu(p, k)) > 0.0_rk)) then
        status = 1
        det = 0.0_rk
        return
      end if

      ! Swap rows k and p (each swap changes the sign of the determinant).
      if (p /= k) then
        row = lu(k, :)
        lu(k, :) = lu(p, :)
        lu(p, :) = row
        tmp = perm(k)
        perm(k) = perm(p)
        perm(p) = tmp
        det = -det
      end if
      det = det*lu(k, k)

      ! Multipliers l_ik and the update of the trailing submatrix.
      do i = k + 1, n
        lu(i, k) = lu(i, k)/lu(k, k)
        lu(i, k + 1:n) = lu(i, k + 1:n) - lu(i, k)*lu(k, k + 1:n)
      end do
    end do

  end subroutine lu_factor

  !> Solve A x = b from the factors of lu_factor: L y = P b (forward substitution), U x = y
  !> (back substitution).
  !>
  !> Source : Golub & Van Loan (2013), §3.1 (Algorithms 3.1.1 and 3.1.2) applied to P b
  !>          (§3.4.4).
  !>
  !> Arguments:
  !>   lu      (in)  real(rk)(:, :)   factors of lu_factor (status 0)
  !>   perm    (in)  integer(:)       row permutation of lu_factor
  !>   b       (in)  real(rk)(:)      right-hand side (n)
  !>   x       (out) real(rk)(:)      solution (n)
  !> Locals:
  !>   n       integer                size
  !>   k       integer                row
  pure subroutine lu_solve(lu, perm, b, x)

    real(rk), intent(in) :: lu(:, :)
    integer, intent(in) :: perm(:)
    real(rk), intent(in) :: b(:)
    real(rk), intent(out) :: x(:)

    integer :: n, k

    n = size(lu, 1)

    ! L y = P b (unit lower triangular).
    do k = 1, n
      x(k) = b(perm(k)) - dot_product(lu(k, 1:k - 1), x(1:k - 1))
    end do

    ! U x = y.
    do k = n, 1, -1
      x(k) = (x(k) - dot_product(lu(k, k + 1:n), x(k + 1:n)))/lu(k, k)
    end do

  end subroutine lu_solve

end module linalg_mod
