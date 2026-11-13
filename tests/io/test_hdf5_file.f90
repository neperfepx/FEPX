! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of hdf5_file_mod (step 9.1, ROADMAP D93), built with FEPX_ENABLE_HDF5 only: the
!> round trip of real and integer tables bit for bit (thirds, negative, tiny, huge, a NaN), their
!> dimensions as rows and columns, a dataset replaced, text and integer attributes written,
!> replaced and read back (on the root and on a dataset), tables of rows appended to, truncated
!> and appended to again, paths that exist or not (a missing group is not an error), a path
!> deleted, and the file reopened read-only.
!>
!> Source : FEPX (ROADMAP step 9.1).
!>
!> Locals:
!>   name      character(*)            scratch file name (working directory)
!>   file      type(h5_file_t)         the file
!>   reals, back   real(rk)(:, :)      a real table written, read back
!>   ints, iback   integer(ik)(:, :)   an integer table written, read back
!>   rows      real(rk)(:, :)          rows appended
!>   text      character(:), allocatable  a text attribute read back
!>   ivals     integer(ik)(:), allocatable  an integer attribute read back
!>   num_rows, num_cols  integer       dimensions of a dataset
!>   unit, k   integer                 unit (to remove an old file), index
!>   exists    logical                 whether the scratch file exists
program test_hdf5_file

  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  use kinds_mod, only: rk, ik, i8
  use hdf5_file_mod, only: hdf5_available, h5_file_t, h5_create, h5_open, h5_close, h5_flush, &
    & h5_write_real, h5_write_int, h5_append_rows, h5_truncate_rows, h5_read_real, h5_read_int, &
    & h5_dims, h5_exists, h5_delete, h5_write_attr_text, h5_write_attr_int, h5_read_attr_text, &
    & h5_read_attr_int
  use assert_mod, only: assert_true, assert_equal, test_summary

  implicit none

  character(len=*), parameter :: name = 'test_hdf5_file.h5'
  type(h5_file_t) :: file
  real(rk), allocatable :: reals(:, :), back(:, :), rows(:, :)
  integer(ik), allocatable :: ints(:, :), iback(:, :), ivals(:)
  character(len=:), allocatable :: text
  integer :: num_rows, num_cols, unit, k
  logical :: exists

  call assert_true(hdf5_available, 'hdf5_available')

  ! A previous run's file removed.
  inquire (file=name, exist=exists)
  if (exists) then
    open (newunit=unit, file=name, status='old')
    close (unit, status='delete')
  end if

  ! Values with every bit used, and a NaN (written as is: the padding of slip variables).
  allocate (reals(3, 4))
  reals = reshape([(real(k, rk)/3.0_rk - 1.0_rk, k = 1, 12)], [3, 4])
  reals(1, 1) = tiny(1.0_rk)
  reals(2, 1) = -huge(1.0_rk)
  reals(3, 4) = ieee_value(1.0_rk, ieee_quiet_nan)
  ints = reshape([(k*(-1)**k, k = 1, 10)], [2, 5])
  ints(1, 1) = huge(1_ik)

  ! Tables, attributes and rows written.
  call h5_create(file, name)
  call h5_write_real(file, '/results/elts/x/x.step0', reals)
  call h5_write_int(file, '/mesh/elts', ints)
  call h5_write_attr_text(file, '/', 'format', 'first')
  call h5_write_attr_text(file, '/', 'format', 'fepx-hdf5 1')
  call h5_write_attr_text(file, '/', 'empty', '')
  call h5_write_attr_int(file, '/', 'sizes', [7_ik, -8_ik, 9_ik])
  call h5_write_attr_text(file, '/results/elts/x/x.step0', 'columns', 'a b c')
  rows = reshape([(real(k, rk), k = 1, 6)], [3, 2])
  call h5_append_rows(file, '/results/table', rows)
  call h5_append_rows(file, '/results/table', rows(:, 1:1) + 10.0_rk)
  call h5_flush(file)

  ! Dimensions: rows are the entities, columns their values.
  call h5_dims(file, '/results/elts/x/x.step0', num_rows, num_cols)
  call assert_equal(num_rows, 4, 'real table rows')
  call assert_equal(num_cols, 3, 'real table columns')
  call h5_dims(file, '/results/table', num_rows, num_cols)
  call assert_equal(num_rows, 3, 'rows appended')
  call assert_equal(num_cols, 3, 'row length')

  ! A dataset replaced by one of another shape.
  call h5_write_real(file, '/results/elts/y/y.step0', reals)
  call h5_write_real(file, '/results/elts/y/y.step0', reals(1:2, 1:3))
  call h5_dims(file, '/results/elts/y/y.step0', num_rows, num_cols)
  call assert_true(num_rows == 3 .and. num_cols == 2, 'dataset replaced')

  ! Paths: a dataset, a group, a missing dataset, a missing group on the path.
  call assert_true(h5_exists(file, '/'), 'root exists')
  call assert_true(h5_exists(file, '/results/elts/x/x.step0'), 'dataset exists')
  call assert_true(h5_exists(file, '/results/elts'), 'group exists')
  call assert_true(.not. h5_exists(file, '/results/elts/x/x.step1'), 'missing dataset')
  call assert_true(.not. h5_exists(file, '/nothing/here/x'), 'missing group')

  ! Deleted, then deleting again does nothing.
  call h5_delete(file, '/results/elts/y/y.step0')
  call assert_true(.not. h5_exists(file, '/results/elts/y/y.step0'), 'deleted')
  call h5_delete(file, '/results/elts/y/y.step0')

  ! Truncated, then appended to again.
  call h5_truncate_rows(file, '/results/table', 1)
  call h5_append_rows(file, '/results/table', rows(:, 2:2) + 20.0_rk)
  call h5_close(file)

  ! Reopened read-only: everything read back bit for bit.
  call h5_open(file, name, read_only=.true.)
  call h5_read_real(file, '/results/elts/x/x.step0', back)
  call assert_true(all(shape(back) == shape(reals)), 'real table shape')
  call assert_true(same_bits(back(:, 1:3), reals(:, 1:3)), 'real table bit for bit')
  call assert_true(same_bits(back(1:2, 4:4), reals(1:2, 4:4)) .and. ieee_is_nan(back(3, 4)), &
    & 'real table last column, NaN kept')
  call h5_read_int(file, '/mesh/elts', iback)
  call assert_true(all(shape(iback) == shape(ints)), 'integer table shape')
  call assert_true(all(iback == ints), 'integer table')
  call h5_read_attr_text(file, '/', 'format', text)
  call assert_equal(text, 'fepx-hdf5 1', 'text attribute replaced')
  call h5_read_attr_text(file, '/', 'empty', text)
  call assert_equal(text, ' ', 'empty text attribute stored as a blank')
  call h5_read_attr_int(file, '/', 'sizes', ivals, 3)
  call assert_true(all(ivals == [7_ik, -8_ik, 9_ik]), 'integer attribute')
  call h5_read_attr_text(file, '/results/elts/x/x.step0', 'columns', text)
  call assert_equal(text, 'a b c', 'attribute of a dataset')
  call h5_read_real(file, '/results/table', back)
  call assert_true(all(shape(back) == [3, 2]), 'table of rows truncated and appended to')
  call assert_true(same_bits(back(:, 1:1), rows(:, 1:1)) .and. same_bits(back(:, 2:2), &
    & rows(:, 2:2) + 20.0_rk), 'table of rows values')
  call h5_close(file)

  call test_summary()

contains

  !> Whether two real tables of the same shape have the same bits.
  !>
  !> Arguments:
  !>   a, b  (in) real(rk)(:, :)  the tables
  !>   same  (out, result) logical
  function same_bits(a, b) result(same)

    real(rk), intent(in) :: a(:, :), b(:, :)
    logical :: same

    same = all(shape(a) == shape(b))
    if (same) same = all(transfer(a, 1_i8, size(a)) == transfer(b, 1_i8, size(b)))

  end function same_bits

end program test_hdf5_file
