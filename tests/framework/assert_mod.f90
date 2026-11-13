! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Minimal in-house assertion module for FEPX unit, verification and MPI tests.
!>
!> Source : FEPX test infrastructure (ROADMAP step 0.1). Contains no numerical method.
!>
!> Usage: a test program calls the assert_* routines any number of times, then calls
!> test_summary once at the end. Each failed check prints one line starting with "FAIL" to
!> standard error. test_summary prints the counts to standard output and ends the program with
!> a non-zero exit status (error stop 1) if any check failed, or if no check was performed, so
!> that CTest marks the test as failed. The module does not use MPI: in an MPI test every rank
!> keeps its own counts, and a failing rank stops the whole job through error stop.
!>
!> Tolerance rule of assert_close (applied element by element for arrays):
!>   |actual - expected| <= max(abs_tol, rel_tol * max(|actual|, |expected|))
!> Both tolerances are optional and default to zero, i.e. exact equality. A NaN in either value
!> always fails. Infinite values pass only if both values are the same infinity.
!>
!> Public routines:
!>   assert_true    check that a logical condition holds
!>   assert_equal   exact comparison of integers (int32, int64) or character strings
!>   assert_close   tolerance comparison of real64 scalars, rank-1 and rank-2 arrays
!>   test_summary   print the counts; stop with a non-zero exit status on failure
!>
!> Module variables:
!>   num_checks    integer(i8)   number of checks performed so far by this process
!>   num_failures  integer(i8)   number of failed checks so far in this process
!> Module parameters:
!>   real_fmt      character     edit descriptor used to print real values (full precision)
!>   detail_len    integer       length of the buffer holding a failure description
module assert_mod

  use, intrinsic :: iso_fortran_env, only: rk => real64, ik => int32, i8 => int64, &
    & output_unit, error_unit
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_class, &
    & ieee_class_type, operator(==)

  implicit none

  private

  public :: assert_true, assert_equal, assert_close, test_summary

  !> Exact comparison of integers or character strings.
  interface assert_equal
    module procedure assert_equal_int32, assert_equal_int64, assert_equal_char
  end interface assert_equal

  !> Tolerance comparison of real64 scalars, rank-1 and rank-2 arrays.
  interface assert_close
    module procedure assert_close_scalar, assert_close_rank1, assert_close_rank2
  end interface assert_close

  integer(i8) :: num_checks = 0_i8
  integer(i8) :: num_failures = 0_i8

  character(len=*), parameter :: real_fmt = 'es25.17'
  integer, parameter :: detail_len = 512

contains

  !> Check that a logical condition holds.
  !>
  !> Arguments:
  !>   condition  (in)  logical    condition that must be true
  !>   label      (in)  character  name of the check, printed on failure
  subroutine assert_true(condition, label)

    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    ! Record the result; the description is only printed on failure.
    call record_check(condition, label, 'condition is false')

  end subroutine assert_true

  !> Check that two default-kind (int32) integers are equal.
  !>
  !> Arguments:
  !>   actual    (in)  integer(ik)  value produced by the code under test
  !>   expected  (in)  integer(ik)  reference value
  !>   label     (in)  character    name of the check, printed on failure
  !> Locals:
  !>   detail    character(detail_len)  failure description
  subroutine assert_equal_int32(actual, expected, label)

    integer(ik), intent(in) :: actual
    integer(ik), intent(in) :: expected
    character(len=*), intent(in) :: label

    character(len=detail_len) :: detail

    ! Build the failure description (cheap; only printed on failure).
    write (detail, '(a,i0,a,i0)') 'actual ', actual, ', expected ', expected

    ! Record the result of the exact comparison.
    call record_check(actual == expected, label, trim(detail))

  end subroutine assert_equal_int32

  !> Check that two int64 integers are equal.
  !>
  !> Arguments:
  !>   actual    (in)  integer(i8)  value produced by the code under test
  !>   expected  (in)  integer(i8)  reference value
  !>   label     (in)  character    name of the check, printed on failure
  !> Locals:
  !>   detail    character(detail_len)  failure description
  subroutine assert_equal_int64(actual, expected, label)

    integer(i8), intent(in) :: actual
    integer(i8), intent(in) :: expected
    character(len=*), intent(in) :: label

    character(len=detail_len) :: detail

    ! Build the failure description.
    write (detail, '(a,i0,a,i0)') 'actual ', actual, ', expected ', expected

    ! Record the result of the exact comparison.
    call record_check(actual == expected, label, trim(detail))

  end subroutine assert_equal_int64

  !> Check that two character strings are equal. Fortran comparison rules apply: the shorter
  !> string is padded with blanks, so trailing blanks are not significant.
  !>
  !> Arguments:
  !>   actual    (in)  character  value produced by the code under test
  !>   expected  (in)  character  reference value
  !>   label     (in)  character  name of the check, printed on failure
  subroutine assert_equal_char(actual, expected, label)

    character(len=*), intent(in) :: actual
    character(len=*), intent(in) :: expected
    character(len=*), intent(in) :: label

    ! Record the result of the comparison, quoting both strings in the description.
    call record_check(actual == expected, label, &
      & 'actual "' // actual // '", expected "' // expected // '"')

  end subroutine assert_equal_char

  !> Check that two real64 scalars agree within the tolerances (module header for the rule).
  !>
  !> Arguments:
  !>   actual    (in)  real(rk)            value produced by the code under test
  !>   expected  (in)  real(rk)            reference value
  !>   label     (in)  character           name of the check, printed on failure
  !>   abs_tol   (in, optional) real(rk)   absolute tolerance, >= 0 (default 0)
  !>   rel_tol   (in, optional) real(rk)   relative tolerance, >= 0 (default 0)
  !> Locals:
  !>   tol_abs   real(rk)                  absolute tolerance actually used
  !>   tol_rel   real(rk)                  relative tolerance actually used
  !>   passed    logical                   result of the comparison
  !>   detail    character(detail_len)     failure description
  subroutine assert_close_scalar(actual, expected, label, abs_tol, rel_tol)

    real(rk), intent(in) :: actual
    real(rk), intent(in) :: expected
    character(len=*), intent(in) :: label
    real(rk), intent(in), optional :: abs_tol
    real(rk), intent(in), optional :: rel_tol

    real(rk) :: tol_abs
    real(rk) :: tol_rel
    logical :: passed
    character(len=detail_len) :: detail

    ! Resolve the optional tolerances.
    call resolve_tolerances(abs_tol, rel_tol, tol_abs, tol_rel)

    ! Compare.
    passed = is_close(actual, expected, tol_abs, tol_rel)

    ! On failure, describe the values and the tolerances.
    detail = ''
    if (.not. passed) then
      write (detail, '(a,' // real_fmt // ',a,' // real_fmt // ',a,es10.3,a,es10.3)') &
        & 'actual ', actual, ', expected ', expected, ', abs_tol ', tol_abs, ', rel_tol ', tol_rel
    end if

    ! Record the result.
    call record_check(passed, label, trim(detail))

  end subroutine assert_close_scalar

  !> Check that two real64 rank-1 arrays have the same size and agree element by element
  !> within the tolerances. On failure the first failing element and the number of failing
  !> elements are reported.
  !>
  !> Arguments:
  !>   actual    (in)  real(rk)(:)         values produced by the code under test
  !>   expected  (in)  real(rk)(:)         reference values
  !>   label     (in)  character           name of the check, printed on failure
  !>   abs_tol   (in, optional) real(rk)   absolute tolerance, >= 0 (default 0)
  !>   rel_tol   (in, optional) real(rk)   relative tolerance, >= 0 (default 0)
  !> Locals:
  !>   tol_abs   real(rk)                  absolute tolerance actually used
  !>   tol_rel   real(rk)                  relative tolerance actually used
  !>   mask      logical(:)                element-wise comparison result
  !>   i_first   integer                   index of the first failing element
  !>   detail    character(detail_len)     failure description
  subroutine assert_close_rank1(actual, expected, label, abs_tol, rel_tol)

    real(rk), intent(in) :: actual(:)
    real(rk), intent(in) :: expected(:)
    character(len=*), intent(in) :: label
    real(rk), intent(in), optional :: abs_tol
    real(rk), intent(in), optional :: rel_tol

    real(rk) :: tol_abs
    real(rk) :: tol_rel
    logical, allocatable :: mask(:)
    integer :: i_first
    character(len=detail_len) :: detail

    ! A size mismatch is a failure; the values are not compared.
    if (size(actual) /= size(expected)) then
      write (detail, '(a,i0,a,i0)') 'size mismatch: actual ', size(actual), &
        & ', expected ', size(expected)
      call record_check(.false., label, trim(detail))
      return
    end if

    ! Resolve the optional tolerances and compare element by element.
    call resolve_tolerances(abs_tol, rel_tol, tol_abs, tol_rel)
    mask = is_close(actual, expected, tol_abs, tol_rel)

    ! On failure, describe the first failing element.
    detail = ''
    if (.not. all(mask)) then
      i_first = findloc(mask, .false., dim=1)
      write (detail, '(i0,a,i0,a,i0,a,' // real_fmt // ',a,' // real_fmt // &
        & ',a,es10.3,a,es10.3)') count(.not. mask), ' of ', size(mask), &
        & ' elements differ; first (', i_first, '): actual ', actual(i_first), &
        & ', expected ', expected(i_first), ', abs_tol ', tol_abs, ', rel_tol ', tol_rel
    end if

    ! Record the result.
    call record_check(all(mask), label, trim(detail))

  end subroutine assert_close_rank1

  !> Check that two real64 rank-2 arrays have the same shape and agree element by element
  !> within the tolerances. On failure the first failing element (column-major order) and
  !> the number of failing elements are reported.
  !>
  !> Arguments:
  !>   actual    (in)  real(rk)(:,:)       values produced by the code under test
  !>   expected  (in)  real(rk)(:,:)       reference values
  !>   label     (in)  character           name of the check, printed on failure
  !>   abs_tol   (in, optional) real(rk)   absolute tolerance, >= 0 (default 0)
  !>   rel_tol   (in, optional) real(rk)   relative tolerance, >= 0 (default 0)
  !> Locals:
  !>   tol_abs   real(rk)                  absolute tolerance actually used
  !>   tol_rel   real(rk)                  relative tolerance actually used
  !>   mask      logical(:,:)              element-wise comparison result
  !>   i_first   integer(2)                (row, column) of the first failing element
  !>   detail    character(detail_len)     failure description
  subroutine assert_close_rank2(actual, expected, label, abs_tol, rel_tol)

    real(rk), intent(in) :: actual(:, :)
    real(rk), intent(in) :: expected(:, :)
    character(len=*), intent(in) :: label
    real(rk), intent(in), optional :: abs_tol
    real(rk), intent(in), optional :: rel_tol

    real(rk) :: tol_abs
    real(rk) :: tol_rel
    logical, allocatable :: mask(:, :)
    integer :: i_first(2)
    character(len=detail_len) :: detail

    ! A shape mismatch is a failure; the values are not compared.
    if (any(shape(actual) /= shape(expected))) then
      write (detail, '(a,i0,a,i0,a,i0,a,i0,a)') 'shape mismatch: actual (', &
        & size(actual, 1), ',', size(actual, 2), '), expected (', &
        & size(expected, 1), ',', size(expected, 2), ')'
      call record_check(.false., label, trim(detail))
      return
    end if

    ! Resolve the optional tolerances and compare element by element.
    call resolve_tolerances(abs_tol, rel_tol, tol_abs, tol_rel)
    mask = is_close(actual, expected, tol_abs, tol_rel)

    ! On failure, describe the first failing element.
    detail = ''
    if (.not. all(mask)) then
      i_first = findloc(mask, .false.)
      write (detail, '(i0,a,i0,a,i0,a,i0,a,' // real_fmt // ',a,' // real_fmt // &
        & ',a,es10.3,a,es10.3)') count(.not. mask), ' of ', size(mask), &
        & ' elements differ; first (', i_first(1), ',', i_first(2), '): actual ', &
        & actual(i_first(1), i_first(2)), ', expected ', expected(i_first(1), i_first(2)), &
        & ', abs_tol ', tol_abs, ', rel_tol ', tol_rel
    end if

    ! Record the result.
    call record_check(all(mask), label, trim(detail))

  end subroutine assert_close_rank2

  !> Print the number of checks and failures, and stop with a non-zero exit status if any
  !> check failed or if no check was performed (a test that checks nothing is a mistake).
  !> When it returns, all checks passed; an MPI test then finalizes MPI normally.
  !>
  !> Arguments:
  !>   verbose  (in, optional) logical  print the summary line when all checks passed
  !>                                    (default .true.; MPI tests pass it only on rank 0).
  !>                                    The summary is always printed on failure.
  !> Locals:
  !>   print_summary  logical           resolved value of verbose
  subroutine test_summary(verbose)

    logical, intent(in), optional :: verbose

    logical :: print_summary

    ! Resolve the optional argument.
    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    ! Print the counts (always on failure).
    if (print_summary .or. num_failures > 0_i8 .or. num_checks == 0_i8) then
      write (output_unit, '(a,i0,a,i0,a)') 'test summary: ', num_checks, ' checks, ', &
        & num_failures, ' failed'
      flush (output_unit)
    end if

    ! A test that performed no check fails.
    if (num_checks == 0_i8) then
      write (error_unit, '(a)') 'FAIL  test_summary : no checks were performed'
      flush (error_unit)
      error stop 1
    end if

    ! Any failed check makes the test fail.
    if (num_failures > 0_i8) error stop 1

  end subroutine test_summary

  !> Count one check, and count and report it if it failed.
  !>
  !> Arguments:
  !>   passed  (in)  logical    result of the check
  !>   label   (in)  character  name of the check
  !>   detail  (in)  character  description printed after the label on failure
  subroutine record_check(passed, label, detail)

    logical, intent(in) :: passed
    character(len=*), intent(in) :: label
    character(len=*), intent(in) :: detail

    ! Count the check.
    num_checks = num_checks + 1_i8

    ! Count and report a failure.
    if (.not. passed) then
      num_failures = num_failures + 1_i8
      write (error_unit, '(a)') 'FAIL  ' // label // ' : ' // detail
      flush (error_unit)
    end if

  end subroutine record_check

  !> Replace absent optional tolerances by their default (zero).
  !>
  !> Arguments:
  !>   abs_tol  (in, optional) real(rk)  absolute tolerance given by the caller
  !>   rel_tol  (in, optional) real(rk)  relative tolerance given by the caller
  !>   tol_abs  (out) real(rk)           absolute tolerance to use
  !>   tol_rel  (out) real(rk)           relative tolerance to use
  pure subroutine resolve_tolerances(abs_tol, rel_tol, tol_abs, tol_rel)

    real(rk), intent(in), optional :: abs_tol
    real(rk), intent(in), optional :: rel_tol
    real(rk), intent(out) :: tol_abs
    real(rk), intent(out) :: tol_rel

    ! Default: exact comparison.
    tol_abs = 0.0_rk
    tol_rel = 0.0_rk

    ! Use the caller's values when given.
    if (present(abs_tol)) tol_abs = abs_tol
    if (present(rel_tol)) tol_rel = rel_tol

  end subroutine resolve_tolerances

  !> Tolerance comparison of two real64 values (rule in the module header). Non-finite values
  !> are classified before any arithmetic, so that no floating-point exception is raised
  !> (Debug builds trap on invalid operations).
  !>
  !> Arguments:
  !>   actual    (in)  real(rk)  value produced by the code under test
  !>   expected  (in)  real(rk)  reference value
  !>   tol_abs   (in)  real(rk)  absolute tolerance, >= 0
  !>   tol_rel   (in)  real(rk)  relative tolerance, >= 0
  !>   (result)        logical   .true. if the values agree
  !> Locals:
  !>   class_actual    ieee_class_type  IEEE class of actual (used for infinities)
  !>   class_expected  ieee_class_type  IEEE class of expected
  elemental function is_close(actual, expected, tol_abs, tol_rel) result(close_enough)

    real(rk), intent(in) :: actual
    real(rk), intent(in) :: expected
    real(rk), intent(in) :: tol_abs
    real(rk), intent(in) :: tol_rel
    logical :: close_enough

    type(ieee_class_type) :: class_actual
    type(ieee_class_type) :: class_expected

    if (ieee_is_nan(actual) .or. ieee_is_nan(expected)) then
      ! A NaN never agrees with anything.
      close_enough = .false.
    else if (.not. (ieee_is_finite(actual) .and. ieee_is_finite(expected))) then
      ! At least one infinity: agree only if both are the same infinity.
      class_actual = ieee_class(actual)
      class_expected = ieee_class(expected)
      close_enough = class_actual == class_expected
    else
      ! Finite values: absolute or relative tolerance.
      close_enough = abs(actual - expected) &
        & <= max(tol_abs, tol_rel*max(abs(actual), abs(expected)))
    end if

  end function is_close

end module assert_mod
