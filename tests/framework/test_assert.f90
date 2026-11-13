! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Self-test of assert_mod: checks that must pass.
!>
!> Source : FEPX test infrastructure (ROADMAP step 0.1).
!>
!> Every check below is constructed to pass; the test fails if any of them is reported as a
!> failure. The complementary test_assert_failures checks that failures are detected.
!>
!> Locals:
!>   vec_a, vec_b    real(rk)(6)    rank-1 arrays that agree within a relative tolerance
!>   mat_a, mat_b    real(rk)(3,3)  rank-2 arrays that agree within an absolute tolerance
!>   pos_inf         real(rk)       +infinity
!>   i               integer        loop index over matrix columns
program test_assert

  use, intrinsic :: iso_fortran_env, only: rk => real64, i8 => int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary

  implicit none

  real(rk) :: vec_a(6)
  real(rk) :: vec_b(6)
  real(rk) :: mat_a(3, 3)
  real(rk) :: mat_b(3, 3)
  real(rk) :: pos_inf
  integer :: i

  ! Logical and exact comparisons.
  call assert_true(.true., 'true_is_true')
  call assert_equal(42, 42, 'int32_equal')
  call assert_equal(huge(0_i8), huge(0_i8), 'int64_equal')
  call assert_equal('fcc', 'fcc', 'char_equal')
  call assert_equal('fcc', 'fcc   ', 'char_equal_trailing_blanks')

  ! Scalar comparisons: exact (default tolerances), absolute, relative, both at the boundary.
  call assert_close(1.5_rk, 1.5_rk, 'scalar_exact')
  call assert_close(1.0_rk, 1.0_rk + 1.0e-9_rk, 'scalar_abs_tol', abs_tol=1.0e-8_rk)
  call assert_close(1.0e6_rk, 1.0e6_rk*(1.0_rk + 1.0e-13_rk), 'scalar_rel_tol', &
    & rel_tol=1.0e-12_rk)
  call assert_close(0.0_rk, 1.0e-20_rk, 'scalar_near_zero_abs_tol', abs_tol=1.0e-15_rk, &
    & rel_tol=1.0e-12_rk)
  call assert_close(2.0_rk, 2.5_rk, 'scalar_abs_tol_boundary', abs_tol=0.5_rk)

  ! Infinities of the same sign agree.
  pos_inf = ieee_value(pos_inf, ieee_positive_inf)
  call assert_close(pos_inf, pos_inf, 'scalar_same_infinity')

  ! Rank-1 arrays within a relative tolerance.
  vec_a = [1.0_rk, -2.0_rk, 3.0_rk, 1.0e3_rk, -1.0e-3_rk, 0.5_rk]
  vec_b = vec_a*(1.0_rk + 1.0e-14_rk)
  call assert_close(vec_a, vec_b, 'rank1_rel_tol', rel_tol=1.0e-12_rk)

  ! Rank-2 arrays within an absolute tolerance.
  do i = 1, 3
    mat_a(:, i) = real(i, rk)*[1.0_rk, 2.0_rk, 3.0_rk]
  end do
  mat_b = mat_a + 1.0e-12_rk
  call assert_close(mat_a, mat_b, 'rank2_abs_tol', abs_tol=1.0e-10_rk)

  ! Report; error stop on any failure.
  call test_summary()

end program test_assert
