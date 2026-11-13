! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit test of kinds_mod: the kind parameters have the sizes required by dev/conventions.md §6.
!>
!> Source : FEPX (ROADMAP step 1.1).
!>
!> Checks:
!>   - rk, ik, i8 are the iso_fortran_env kinds real64, int32, int64;
!>   - real(rk) is IEEE binary64: 64 bits, 53-bit significand, decimal exponent range >= 307;
!>   - integer(ik) has 32 bits and integer(i8) 64 bits, so i8 holds counts above 2^31 - 1.
!>
!> Locals: none.
program test_kinds

  use, intrinsic :: iso_fortran_env, only: real64, int32, int64
  use, intrinsic :: ieee_arithmetic, only: ieee_support_datatype
  use kinds_mod, only: rk, ik, i8
  use assert_mod, only: assert_true, assert_equal, test_summary

  implicit none

  ! The kinds are the ones named in the conventions.
  call assert_equal(rk, real64, 'rk_is_real64')
  call assert_equal(ik, int32, 'ik_is_int32')
  call assert_equal(i8, int64, 'i8_is_int64')

  ! real(rk) is IEEE double precision.
  call assert_equal(storage_size(1.0_rk), 64, 'rk_storage_bits')
  call assert_equal(digits(1.0_rk), 53, 'rk_significand_bits')
  call assert_true(range(1.0_rk) >= 307, 'rk_exponent_range')
  call assert_true(ieee_support_datatype(1.0_rk), 'rk_is_ieee')

  ! Integer sizes; i8 holds values above the int32 range (global counts).
  call assert_equal(storage_size(1_ik), 32, 'ik_storage_bits')
  call assert_equal(storage_size(1_i8), 64, 'i8_storage_bits')
  call assert_true(huge(1_i8) > int(huge(1_ik), i8), 'i8_exceeds_ik')

  ! Report; error stop on any failure.
  call test_summary()

end program test_kinds
