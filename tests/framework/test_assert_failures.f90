! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Self-test of assert_mod: checks that must fail.
!>
!> Source : FEPX test infrastructure (ROADMAP step 0.1).
!>
!> Run with one command-line argument naming a case. Each case performs one check that passes
!> and one check that must fail, then calls test_summary. CTest (tests/CMakeLists.txt) runs
!> every case and requires the output "test summary: 2 checks, 1 failed", which proves that the
!> failure was detected and counted; one case is also run with WILL_FAIL to check the non-zero
!> exit status. The case "no_checks" performs no check and must be reported as a failure.
!>
!> Locals:
!>   case_name   character(64)  case selected on the command line
!>   nan         real(rk)       quiet NaN
!>   pos_inf     real(rk)       +infinity
!>   neg_inf     real(rk)       -infinity
!>   vec         real(rk)(3)    rank-1 reference array
!>   mat         real(rk)(2,2)  rank-2 reference array
program test_assert_failures

  use, intrinsic :: iso_fortran_env, only: rk => real64, i8 => int64, error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, &
    & ieee_negative_inf
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary

  implicit none

  character(len=64) :: case_name
  real(rk) :: nan
  real(rk) :: pos_inf
  real(rk) :: neg_inf
  real(rk) :: vec(3)
  real(rk) :: mat(2, 2)

  ! Special values and reference arrays used by the cases.
  nan = ieee_value(nan, ieee_quiet_nan)
  pos_inf = ieee_value(pos_inf, ieee_positive_inf)
  neg_inf = ieee_value(neg_inf, ieee_negative_inf)
  vec = [1.0_rk, 2.0_rk, 3.0_rk]
  mat = reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [2, 2])

  ! Read the case name.
  if (command_argument_count() /= 1) then
    write (error_unit, '(a)') 'usage: test_assert_failures <case>'
    error stop 2
  end if
  call get_command_argument(1, case_name)

  ! The case "no_checks" calls test_summary without any check.
  if (case_name == 'no_checks') then
    call test_summary()
    stop
  end if

  ! One check that passes, common to all other cases.
  call assert_true(.true., 'passing_check')

  ! One check that must fail.
  select case (trim(case_name))
  case ('true')
    call assert_true(.false., 'case_true')
  case ('int32')
    call assert_equal(1, 2, 'case_int32')
  case ('int64')
    call assert_equal(1_i8, -1_i8, 'case_int64')
  case ('char')
    call assert_equal('fcc', 'bcc', 'case_char')
  case ('scalar_exact')
    call assert_close(1.0_rk, 1.0_rk + epsilon(1.0_rk), 'case_scalar_exact')
  case ('scalar_abs')
    call assert_close(1.0_rk, 1.1_rk, 'case_scalar_abs', abs_tol=0.05_rk)
  case ('scalar_rel')
    call assert_close(1.0e6_rk, 1.0e6_rk + 2.0_rk, 'case_scalar_rel', rel_tol=1.0e-6_rk)
  case ('scalar_nan')
    call assert_close(nan, nan, 'case_scalar_nan', abs_tol=1.0_rk)
  case ('scalar_inf_sign')
    call assert_close(pos_inf, neg_inf, 'case_scalar_inf_sign', abs_tol=1.0_rk)
  case ('scalar_inf_finite')
    call assert_close(pos_inf, 1.0_rk, 'case_scalar_inf_finite', rel_tol=1.0_rk)
  case ('rank1_value')
    call assert_close(vec, [1.0_rk, 2.5_rk, 3.0_rk], 'case_rank1_value', abs_tol=0.1_rk)
  case ('rank1_nan')
    call assert_close(vec, [1.0_rk, nan, 3.0_rk], 'case_rank1_nan', abs_tol=1.0_rk)
  case ('rank1_size')
    call assert_close(vec, [1.0_rk, 2.0_rk], 'case_rank1_size', abs_tol=1.0_rk)
  case ('rank2_value')
    call assert_close(mat, reshape([1.0_rk, 2.0_rk, 3.0_rk, 5.0_rk], [2, 2]), &
      & 'case_rank2_value', rel_tol=0.1_rk)
  case ('rank2_shape')
    call assert_close(mat, reshape([1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk], [4, 1]), &
      & 'case_rank2_shape', abs_tol=1.0_rk)
  case default
    write (error_unit, '(a)') 'unknown case: ' // trim(case_name)
    error stop 2
  end select

  ! Report; must print "2 checks, 1 failed" and error stop.
  call test_summary()

end program test_assert_failures
