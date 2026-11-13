! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> HDF5 files, stub for builds without HDF5 (FEPX_ENABLE_HDF5=OFF; step 9.1, ROADMAP D93): the
!> module of hdf5_file_mod.f90 with the same name and interface, hdf5_available = .false., and
!> every routine a fatal error. The configuration rejects output_format hdf5 in such a build
!> (config_mod), so these routines are not reached by a run.
!>
!> Source : FEPX (ROADMAP step 9.1). Contains no numerical method. The arguments are those of
!>          hdf5_file_mod.f90 (documented there) and are not used here (the file is compiled with
!>          -Wno-unused-dummy-argument, src/CMakeLists.txt).
module hdf5_file_mod

  use kinds_mod, only: rk, ik
  use error_mod, only: internal_error

  implicit none

  private

  public :: hdf5_available, h5_file_t
  public :: h5_create, h5_open, h5_close, h5_flush
  public :: h5_write_real, h5_write_int, h5_append_rows, h5_truncate_rows
  public :: h5_read_real, h5_read_int, h5_dims, h5_exists, h5_delete
  public :: h5_write_attr_text, h5_write_attr_int, h5_read_attr_text, h5_read_attr_int
  public :: h5_library_version

  !> This build does not write HDF5.
  logical, parameter :: hdf5_available = .false.

  !> The message of every routine.
  character(len=*), parameter :: no_hdf5 = 'this FEPX build has no HDF5 output (configure ' // &
    & 'with -DFEPX_ENABLE_HDF5=ON)'

  !> An open file (never opened here).
  !>
  !> Components:
  !>   id    integer                    -1
  !>   name  character(:), allocatable  file name
  type :: h5_file_t
    integer :: id = -1
    character(len=:), allocatable :: name
  end type h5_file_t

contains

  !> Stub: a fatal error.
  subroutine h5_create(file, name)
    type(h5_file_t), intent(out) :: file
    character(len=*), intent(in) :: name
    call internal_error('h5_create', no_hdf5)
  end subroutine h5_create

  !> Stub: a fatal error.
  subroutine h5_open(file, name, read_only)
    type(h5_file_t), intent(out) :: file
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: read_only
    call internal_error('h5_open', no_hdf5)
  end subroutine h5_open

  !> Stub: nothing (no file is ever open).
  subroutine h5_close(file)
    type(h5_file_t), intent(inout) :: file
  end subroutine h5_close

  !> Stub: a fatal error.
  subroutine h5_flush(file)
    type(h5_file_t), intent(in) :: file
    call internal_error('h5_flush', no_hdf5)
  end subroutine h5_flush

  !> Stub: a fatal error.
  subroutine h5_write_real(file, path, values)
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    real(rk), intent(in) :: values(:, :)
    call internal_error('h5_write_real', no_hdf5)
  end subroutine h5_write_real

  !> Stub: a fatal error.
  subroutine h5_write_int(file, path, values)
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(ik), intent(in) :: values(:, :)
    call internal_error('h5_write_int', no_hdf5)
  end subroutine h5_write_int

  !> Stub: a fatal error.
  subroutine h5_append_rows(file, path, rows)
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    real(rk), intent(in) :: rows(:, :)
    call internal_error('h5_append_rows', no_hdf5)
  end subroutine h5_append_rows

  !> Stub: a fatal error.
  subroutine h5_truncate_rows(file, path, num_rows)
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer, intent(in) :: num_rows
    call internal_error('h5_truncate_rows', no_hdf5)
  end subroutine h5_truncate_rows

  !> Stub: a fatal error.
  subroutine h5_read_real(file, path, values)
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    real(rk), allocatable, target, intent(out) :: values(:, :)
    call internal_error('h5_read_real', no_hdf5)
  end subroutine h5_read_real

  !> Stub: a fatal error.
  subroutine h5_read_int(file, path, values)
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(ik), allocatable, target, intent(out) :: values(:, :)
    call internal_error('h5_read_int', no_hdf5)
  end subroutine h5_read_int

  !> Stub: a fatal error.
  subroutine h5_dims(file, path, num_rows, num_cols)
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer, intent(out) :: num_rows, num_cols
    num_rows = 0
    num_cols = 0
    call internal_error('h5_dims', no_hdf5)
  end subroutine h5_dims

  !> Stub: a fatal error.
  function h5_exists(file, path) result(exists)
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    logical :: exists
    exists = .false.
    call internal_error('h5_exists', no_hdf5)
  end function h5_exists

  !> Stub: a fatal error.
  subroutine h5_delete(file, path)
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    call internal_error('h5_delete', no_hdf5)
  end subroutine h5_delete

  !> Stub: a fatal error.
  subroutine h5_write_attr_text(file, path, name, value)
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path, name, value
    call internal_error('h5_write_attr_text', no_hdf5)
  end subroutine h5_write_attr_text

  !> Stub: a fatal error.
  subroutine h5_write_attr_int(file, path, name, values)
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path, name
    integer(ik), intent(in) :: values(:)
    call internal_error('h5_write_attr_int', no_hdf5)
  end subroutine h5_write_attr_int

  !> Stub: a fatal error.
  subroutine h5_read_attr_text(file, path, name, value)
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path, name
    character(len=:), allocatable, intent(out) :: value
    call internal_error('h5_read_attr_text', no_hdf5)
  end subroutine h5_read_attr_text

  !> Stub: a fatal error.
  subroutine h5_read_attr_int(file, path, name, values, num)
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path, name
    integer(ik), allocatable, target, intent(out) :: values(:)
    integer, intent(in) :: num
    call internal_error('h5_read_attr_int', no_hdf5)
  end subroutine h5_read_attr_int

  !> Stub: no library (an empty text).
  function h5_library_version() result(text)
    character(len=:), allocatable :: text
    text = ''
  end function h5_library_version

end module hdf5_file_mod
