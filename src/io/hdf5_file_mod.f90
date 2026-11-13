! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> HDF5 files (step 9.1, ROADMAP D10, D93): the small set of operations the output writer
!> needs, on a file opened by one process (serial HDF5; FEPX writes it from the root rank only).
!> This is the only module that uses the HDF5 library; builds without HDF5
!> (FEPX_ENABLE_HDF5=OFF) compile hdf5_file_stub.f90 instead, a module of the same name and
!> interface whose routines are fatal errors.
!>
!> Source : FEPX (ROADMAP step 9.1, dev/io_format.md §5). Contains no numerical method. The
!>          library calls follow the HDF5 Fortran reference manual (HDF5 1.14, the
!>          "Fortran 2003" interfaces with type(c_ptr) buffers for numeric data). v2 wrote HDF5
!>          through parallel MPI-IO (write_res_hdf5.f90, in the tag v2-final-comparison); FEPX
!>          writes serially (D10).
!>
!> Conventions:
!>   - Every dataset is two-dimensional. A Fortran array table(k, n) (k values of each of n
!>     entities, column-major) is stored with the HDF5 dimensions (n, k) (row-major, as C,
!>     h5py and h5dump show them): row j holds the k values of entity j, as one line of the
!>     .sim text files.
!>   - Reals are stored as 64-bit IEEE floats (H5T_NATIVE_DOUBLE), integers as 32-bit
!>     integers (the kind ik), text attributes as fixed-length strings, null-padded.
!>   - Paths are absolute ("/results/elts/stress/stress.step1"); missing groups on the path
!>     are created.
!>   - Tables of rows (h5_append_rows) are extendable in the number of rows (chunked), so that
!>     a restarted run can cut them back to their rows at the restart point (h5_truncate_rows).
!>   - Every error of the library is a fatal error naming the file, the path and the call.
!>
!> Public parameters:
!>   hdf5_available   .true. here (.false. in the stub)
!> Public types:
!>   h5_file_t        an open file (identifier and name)
!> Public routines:
!>   h5_create, h5_open, h5_close, h5_flush   create (an existing file is an error), open,
!>                    close, flush a file
!>   h5_write_real, h5_write_int   write a new two-dimensional dataset (an existing one is
!>                    replaced)
!>   h5_append_rows   append rows to a table of rows (created with its row length when absent)
!>   h5_truncate_rows keep the first rows of a table of rows
!>   h5_read_real, h5_read_int     read a dataset
!>   h5_dims          the dimensions (rows, columns) of a dataset
!>   h5_exists        whether a path (group or dataset) exists
!>   h5_delete        remove a path if it exists
!>   h5_write_attr_text, h5_write_attr_int   write (or replace) an attribute of a group or
!>                    dataset
!>   h5_read_attr_text, h5_read_attr_int     read one
!>   h5_library_version  "HDF5 <major>.<minor>.<release>" of the library linked
module hdf5_file_mod

  use, intrinsic :: iso_c_binding, only: c_loc, c_ptr
  use kinds_mod, only: rk, ik
  use error_mod, only: fatal_error, internal_error
  use hdf5, only: hid_t, hsize_t, size_t, h5open_f, h5fcreate_f, h5fopen_f, h5fclose_f, &
    & h5fflush_f, h5screate_simple_f, h5screate_f, h5sclose_f, h5sselect_hyperslab_f, &
    & h5sget_simple_extent_dims_f, h5pcreate_f, h5pclose_f, h5pset_create_inter_group_f, &
    & h5pset_chunk_f, h5dcreate_f, h5dopen_f, h5dclose_f, h5dwrite_f, h5dread_f, &
    & h5dget_space_f, h5dset_extent_f, h5lexists_f, h5ldelete_f, h5oopen_f, h5oclose_f, &
    & h5acreate_f, h5aopen_f, h5aclose_f, h5awrite_f, h5aread_f, h5aexists_f, h5adelete_f, &
    & h5aget_type_f, h5get_libversion_f, h5tcopy_f, h5tset_size_f, h5tset_strpad_f, h5tget_size_f, h5tclose_f, &
    & h5kind_to_type, H5F_ACC_EXCL_F, H5F_ACC_RDWR_F, H5F_ACC_RDONLY_F, H5F_SCOPE_GLOBAL_F, &
    & H5S_SCALAR_F, H5S_SELECT_SET_F, H5S_UNLIMITED_F, H5P_LINK_CREATE_F, &
    & H5P_DATASET_CREATE_F, H5T_FORTRAN_S1, H5T_STR_NULLPAD_F, H5_REAL_KIND, H5_INTEGER_KIND

  implicit none

  private

  public :: hdf5_available, h5_file_t
  public :: h5_create, h5_open, h5_close, h5_flush
  public :: h5_write_real, h5_write_int, h5_append_rows, h5_truncate_rows
  public :: h5_read_real, h5_read_int, h5_dims, h5_exists, h5_delete
  public :: h5_write_attr_text, h5_write_attr_int, h5_read_attr_text, h5_read_attr_int
  public :: h5_library_version

  !> This build writes HDF5.
  logical, parameter :: hdf5_available = .true.

  !> Rows per chunk of a table of rows.
  integer(hsize_t), parameter :: chunk_rows = 64_hsize_t

  !> An open file.
  !>
  !> Components:
  !>   id    integer(hid_t)             file identifier (-1: not open)
  !>   name  character(:), allocatable  file name, for messages
  type :: h5_file_t
    integer(hid_t) :: id = -1_hid_t
    character(len=:), allocatable :: name
  end type h5_file_t

  !> The library's Fortran interface has been initialized (h5open_f, once per run).
  logical, save :: library_open = .false.

contains

  !> Create a new file (an existing file of that name is an error).
  !>
  !> Arguments:
  !>   file  (out) type(h5_file_t)  the open file
  !>   name  (in)  character        file name
  !> Locals:
  !>   err   integer                library status
  subroutine h5_create(file, name)

    type(h5_file_t), intent(out) :: file
    character(len=*), intent(in) :: name

    integer :: err

    ! The library, then the file (H5F_ACC_EXCL_F: fails if it exists).
    call init_library()
    file%name = name
    call h5fcreate_f(name, H5F_ACC_EXCL_F, file%id, err)
    call check(err, file, '/', 'h5fcreate_f (does the file exist already?)')

  end subroutine h5_create

  !> Open an existing file for reading and writing, or for reading only.
  !>
  !> Arguments:
  !>   file       (out) type(h5_file_t)  the open file
  !>   name       (in)  character        file name
  !>   read_only  (in, optional) logical open for reading only (default false)
  !> Locals:
  !>   err        integer                library status
  !>   flags      integer                access flags
  !>   exists     logical                whether the file exists
  subroutine h5_open(file, name, read_only)

    type(h5_file_t), intent(out) :: file
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: read_only

    integer :: err, flags
    logical :: exists

    ! A missing file is an error with its name (before the library's own message).
    inquire (file=name, exist=exists)
    if (.not. exists) call fatal_error('h5_open', 'HDF5 file ' // name // ' not found')
    call init_library()
    flags = H5F_ACC_RDWR_F
    if (present(read_only)) then
      if (read_only) flags = H5F_ACC_RDONLY_F
    end if
    file%name = name
    call h5fopen_f(name, flags, file%id, err)
    call check(err, file, '/', 'h5fopen_f')

  end subroutine h5_open

  !> Close a file (nothing if it is not open).
  !>
  !> Arguments:
  !>   file  (inout) type(h5_file_t)  the file
  !> Locals:
  !>   err   integer                  library status
  subroutine h5_close(file)

    type(h5_file_t), intent(inout) :: file

    integer :: err

    if (file%id < 0_hid_t) return
    call h5fclose_f(file%id, err)
    call check(err, file, '/', 'h5fclose_f')
    file%id = -1_hid_t

  end subroutine h5_close

  !> Write the file's buffers to disk (so that a run that stops keeps what was written).
  !>
  !> Arguments:
  !>   file  (in) type(h5_file_t)  the open file
  !> Locals:
  !>   err   integer               library status
  subroutine h5_flush(file)

    type(h5_file_t), intent(in) :: file

    integer :: err

    call h5fflush_f(file%id, H5F_SCOPE_GLOBAL_F, err)
    call check(err, file, '/', 'h5fflush_f')

  end subroutine h5_flush

  !> Write a real table as a new dataset (an existing dataset of that path is replaced).
  !>
  !> Arguments:
  !>   file    (in) type(h5_file_t)    the open file
  !>   path    (in) character          dataset path
  !>   values  (in) real(rk)(:, :)     (k, n) values, stored as (n, k)
  !> Locals:
  !>   buf     real(rk)(:, :), target  contiguous copy (for c_loc)
  subroutine h5_write_real(file, path, values)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    real(rk), intent(in) :: values(:, :)

    real(rk), allocatable, target :: buf(:, :)

    buf = values
    call write_new(file, path, h5kind_to_type(rk, H5_REAL_KIND), shape(buf, kind=hsize_t), &
      & c_loc(buf))

  end subroutine h5_write_real

  !> Write an integer table as a new dataset (an existing dataset of that path is replaced).
  !>
  !> Arguments:
  !>   file    (in) type(h5_file_t)       the open file
  !>   path    (in) character             dataset path
  !>   values  (in) integer(ik)(:, :)     (k, n) values, stored as (n, k)
  !> Locals:
  !>   buf     integer(ik)(:, :), target  contiguous copy (for c_loc)
  subroutine h5_write_int(file, path, values)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(ik), intent(in) :: values(:, :)

    integer(ik), allocatable, target :: buf(:, :)

    buf = values
    call write_new(file, path, h5kind_to_type(ik, H5_INTEGER_KIND), shape(buf, kind=hsize_t), &
      & c_loc(buf))

  end subroutine h5_write_int

  !> Create a dataset of the given type and dimensions (Fortran order (k, n)) and write it,
  !> replacing an existing one; the groups on the path are created.
  !>
  !> Arguments:
  !>   file    (in) type(h5_file_t)    the open file
  !>   path    (in) character          dataset path
  !>   dtype   (in) integer(hid_t)     memory and file type
  !>   dims    (in) integer(hsize_t)(2) dimensions (k, n)
  !>   ptr     (in) type(c_ptr)        the values
  !> Locals:
  !>   space, lcpl, dset  integer(hid_t)  dataspace, link-creation list, dataset
  !>   err     integer                 library status
  subroutine write_new(file, path, dtype, dims, ptr)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: dtype
    integer(hsize_t), intent(in) :: dims(2)
    type(c_ptr), intent(in) :: ptr

    integer(hid_t) :: space, lcpl, dset
    integer :: err

    ! A dataset replaced (restart: a step written again).
    call h5_delete(file, path)

    ! Dataspace, the link list that creates the missing groups, the dataset.
    call h5screate_simple_f(2, dims, space, err)
    call check(err, file, path, 'h5screate_simple_f')
    lcpl = link_list(file, path)
    call h5dcreate_f(file%id, path, dtype, space, dset, err, lcpl_id=lcpl)
    call check(err, file, path, 'h5dcreate_f')

    ! The values, then the identifiers closed.
    call h5dwrite_f(dset, dtype, ptr, err)
    call check(err, file, path, 'h5dwrite_f')
    call h5dclose_f(dset, err)
    call check(err, file, path, 'h5dclose_f')
    call h5pclose_f(lcpl, err)
    call check(err, file, path, 'h5pclose_f')
    call h5sclose_f(space, err)
    call check(err, file, path, 'h5sclose_f')

  end subroutine write_new

  !> Append rows to a real table of rows, created (extendable, chunked) with the rows' length
  !> when absent; rows of another length are an error.
  !>
  !> Arguments:
  !>   file    (in) type(h5_file_t)    the open file
  !>   path    (in) character          dataset path
  !>   rows    (in) real(rk)(:, :)     (k, m) m rows of k values
  !> Locals:
  !>   buf     real(rk)(:, :), target  contiguous copy (for c_loc)
  !>   dims    integer(hsize_t)(2)     (k, rows) of the table; then its new dimensions
  !>   maxdims integer(hsize_t)(2)     (k, unlimited)
  !>   counts  integer(hsize_t)(2)     (k, m)
  !>   offsets integer(hsize_t)(2)     (0, rows before)
  !>   space, mem, dcpl, lcpl, dset  integer(hid_t)  identifiers
  !>   err     integer                 library status
  subroutine h5_append_rows(file, path, rows)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    real(rk), intent(in) :: rows(:, :)

    real(rk), allocatable, target :: buf(:, :)
    integer(hsize_t) :: dims(2), maxdims(2), counts(2), offsets(2)
    integer(hid_t) :: space, mem, dcpl, lcpl, dset
    integer :: err

    buf = rows
    counts = shape(buf, kind=hsize_t)
    if (h5_exists(file, path)) then
      ! The existing table: its row length must match.
      call h5dopen_f(file%id, path, dset, err)
      call check(err, file, path, 'h5dopen_f')
      dims = table_dims(file, path, dset)
      if (dims(1) /= counts(1)) call internal_error('h5_append_rows', file%name // ':' // path // &
        & ': rows of another length than the table''s')
    else
      ! A new, empty table extendable in its rows (chunks of chunk_rows rows).
      dims = [counts(1), 0_hsize_t]
      maxdims = [counts(1), H5S_UNLIMITED_F]
      call h5screate_simple_f(2, dims, space, err, maxdims)
      call check(err, file, path, 'h5screate_simple_f')
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl, err)
      call check(err, file, path, 'h5pcreate_f')
      call h5pset_chunk_f(dcpl, 2, [counts(1), chunk_rows], err)
      call check(err, file, path, 'h5pset_chunk_f')
      lcpl = link_list(file, path)
      call h5dcreate_f(file%id, path, h5kind_to_type(rk, H5_REAL_KIND), space, dset, err, &
        & dcpl_id=dcpl, lcpl_id=lcpl)
      call check(err, file, path, 'h5dcreate_f')
      call h5pclose_f(lcpl, err)
      call check(err, file, path, 'h5pclose_f')
      call h5pclose_f(dcpl, err)
      call check(err, file, path, 'h5pclose_f')
      call h5sclose_f(space, err)
      call check(err, file, path, 'h5sclose_f')
    end if

    ! The table extended by the new rows, which are written into the selected block.
    offsets = [0_hsize_t, dims(2)]
    dims(2) = dims(2) + counts(2)
    call h5dset_extent_f(dset, dims, err)
    call check(err, file, path, 'h5dset_extent_f')
    call h5dget_space_f(dset, space, err)
    call check(err, file, path, 'h5dget_space_f')
    call h5sselect_hyperslab_f(space, H5S_SELECT_SET_F, offsets, counts, err)
    call check(err, file, path, 'h5sselect_hyperslab_f')
    call h5screate_simple_f(2, counts, mem, err)
    call check(err, file, path, 'h5screate_simple_f')
    call h5dwrite_f(dset, h5kind_to_type(rk, H5_REAL_KIND), c_loc(buf), err, mem_space_id=mem, &
      & file_space_id=space)
    call check(err, file, path, 'h5dwrite_f')

    ! Identifiers closed.
    call h5sclose_f(mem, err)
    call check(err, file, path, 'h5sclose_f')
    call h5sclose_f(space, err)
    call check(err, file, path, 'h5sclose_f')
    call h5dclose_f(dset, err)
    call check(err, file, path, 'h5dclose_f')

  end subroutine h5_append_rows

  !> Keep the first num_rows rows of a table of rows (step 8.4: a restarted run continues
  !> from the rows at its restart point); more rows than the table has is an error.
  !>
  !> Arguments:
  !>   file      (in) type(h5_file_t)  the open file
  !>   path      (in) character        dataset path
  !>   num_rows  (in) integer          rows to keep
  !> Locals:
  !>   dims      integer(hsize_t)(2)   (k, rows)
  !>   dset      integer(hid_t)        dataset
  !>   err       integer               library status
  subroutine h5_truncate_rows(file, path, num_rows)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer, intent(in) :: num_rows

    integer(hsize_t) :: dims(2)
    integer(hid_t) :: dset
    integer :: err

    call h5dopen_f(file%id, path, dset, err)
    call check(err, file, path, 'h5dopen_f')
    dims = table_dims(file, path, dset)
    if (int(num_rows, hsize_t) > dims(2)) call fatal_error('h5_truncate_rows', file%name // &
      & ':' // path // ' has fewer rows than at the restart point')
    dims(2) = int(num_rows, hsize_t)
    call h5dset_extent_f(dset, dims, err)
    call check(err, file, path, 'h5dset_extent_f')
    call h5dclose_f(dset, err)
    call check(err, file, path, 'h5dclose_f')

  end subroutine h5_truncate_rows

  !> Read a real dataset.
  !>
  !> Arguments:
  !>   file    (in)  type(h5_file_t)          the open file
  !>   path    (in)  character                dataset path
  !>   values  (out) real(rk)(:, :), allocatable, target  (k, n) values
  !> Locals:
  !>   dset    integer(hid_t)                 dataset
  !>   dims    integer(hsize_t)(2)            (k, n)
  !>   ptr     type(c_ptr)                    address of the values (the library's buffer)
  !>   err     integer                        library status
  subroutine h5_read_real(file, path, values)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    real(rk), allocatable, target, intent(out) :: values(:, :)

    integer(hid_t) :: dset
    integer(hsize_t) :: dims(2)
    type(c_ptr) :: ptr
    integer :: err

    call h5dopen_f(file%id, path, dset, err)
    call check(err, file, path, 'h5dopen_f')
    dims = table_dims(file, path, dset)
    allocate (values(dims(1), dims(2)))
    if (size(values) > 0) then
      ptr = c_loc(values)
      call h5dread_f(dset, h5kind_to_type(rk, H5_REAL_KIND), ptr, err)
      call check(err, file, path, 'h5dread_f')
    end if
    call h5dclose_f(dset, err)
    call check(err, file, path, 'h5dclose_f')

  end subroutine h5_read_real

  !> Read an integer dataset.
  !>
  !> Arguments:
  !>   file    (in)  type(h5_file_t)          the open file
  !>   path    (in)  character                dataset path
  !>   values  (out) integer(ik)(:, :), allocatable, target  (k, n) values
  !> Locals:
  !>   dset    integer(hid_t)                 dataset
  !>   dims    integer(hsize_t)(2)            (k, n)
  !>   ptr     type(c_ptr)                    address of the values (the library's buffer)
  !>   err     integer                        library status
  subroutine h5_read_int(file, path, values)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(ik), allocatable, target, intent(out) :: values(:, :)

    integer(hid_t) :: dset
    integer(hsize_t) :: dims(2)
    type(c_ptr) :: ptr
    integer :: err

    call h5dopen_f(file%id, path, dset, err)
    call check(err, file, path, 'h5dopen_f')
    dims = table_dims(file, path, dset)
    allocate (values(dims(1), dims(2)))
    if (size(values) > 0) then
      ptr = c_loc(values)
      call h5dread_f(dset, h5kind_to_type(ik, H5_INTEGER_KIND), ptr, err)
      call check(err, file, path, 'h5dread_f')
    end if
    call h5dclose_f(dset, err)
    call check(err, file, path, 'h5dclose_f')

  end subroutine h5_read_int

  !> The dimensions of a dataset as the HDF5 tools show them: rows (entities), columns (values
  !> per entity).
  !>
  !> Arguments:
  !>   file     (in)  type(h5_file_t)  the open file
  !>   path     (in)  character        dataset path
  !>   num_rows, num_cols  (out) integer  its dimensions
  !> Locals:
  !>   dset     integer(hid_t)         dataset
  !>   dims     integer(hsize_t)(2)    (k, n)
  !>   err      integer                library status
  subroutine h5_dims(file, path, num_rows, num_cols)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer, intent(out) :: num_rows, num_cols

    integer(hid_t) :: dset
    integer(hsize_t) :: dims(2)
    integer :: err

    call h5dopen_f(file%id, path, dset, err)
    call check(err, file, path, 'h5dopen_f')
    dims = table_dims(file, path, dset)
    num_rows = int(dims(2))
    num_cols = int(dims(1))
    call h5dclose_f(dset, err)
    call check(err, file, path, 'h5dclose_f')

  end subroutine h5_dims

  !> The Fortran-order dimensions (k, n) of an open two-dimensional dataset.
  !>
  !> Arguments:
  !>   file, path  (in) type(h5_file_t), character  for messages
  !>   dset        (in) integer(hid_t)   the open dataset
  !>   dims        (out, result) integer(hsize_t)(2)  (k, n)
  !> Locals:
  !>   space       integer(hid_t)        its dataspace
  !>   maxdims     integer(hsize_t)(2)   its largest dimensions (unused)
  !>   err         integer               library status (the rank on success)
  function table_dims(file, path, dset) result(dims)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: dset
    integer(hsize_t) :: dims(2)

    integer(hid_t) :: space
    integer(hsize_t) :: maxdims(2)
    integer :: err

    call h5dget_space_f(dset, space, err)
    call check(err, file, path, 'h5dget_space_f')
    call h5sget_simple_extent_dims_f(space, dims, maxdims, err)
    if (err /= 2) call fatal_error('h5_file', file%name // ':' // path // &
      & ' is not a two-dimensional dataset')
    call h5sclose_f(space, err)
    call check(err, file, path, 'h5sclose_f')

  end function table_dims

  !> Whether a path exists: every link on it is checked in turn (a missing group is not an
  !> error).
  !>
  !> Arguments:
  !>   file    (in) type(h5_file_t)  the open file
  !>   path    (in) character        absolute path
  !>   exists  (out, result) logical
  !> Locals:
  !>   pos     integer               end of the current prefix
  !>   err     integer               library status
  function h5_exists(file, path) result(exists)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    logical :: exists

    integer :: pos, err

    ! The root always exists; then each prefix "/a", "/a/b", ..., then the whole path.
    exists = .true.
    if (path == '/') return
    do pos = 2, len(path)
      if (path(pos:pos) /= '/') cycle
      call h5lexists_f(file%id, path(1:pos - 1), exists, err)
      call check(err, file, path, 'h5lexists_f')
      if (.not. exists) return
    end do
    call h5lexists_f(file%id, path, exists, err)
    call check(err, file, path, 'h5lexists_f')

  end function h5_exists

  !> Remove a path (its link; the space is not reclaimed) if it exists.
  !>
  !> Arguments:
  !>   file  (in) type(h5_file_t)  the open file
  !>   path  (in) character        absolute path
  !> Locals:
  !>   err   integer               library status
  subroutine h5_delete(file, path)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path

    integer :: err

    if (.not. h5_exists(file, path)) return
    call h5ldelete_f(file%id, path, err)
    call check(err, file, path, 'h5ldelete_f')

  end subroutine h5_delete

  !> Write a text attribute of a group or dataset (an existing attribute of that name is
  !> replaced).
  !>
  !> Arguments:
  !>   file   (in) type(h5_file_t)  the open file
  !>   path   (in) character        the object (e.g. "/", "/results/convergence")
  !>   name   (in) character        attribute name
  !>   value  (in) character        its text (at least one character is stored)
  !> Locals:
  !>   obj, atype, space, attr  integer(hid_t)  identifiers
  !>   text   character(:), allocatable  the stored text (a blank for an empty value)
  !>   err    integer               library status
  subroutine h5_write_attr_text(file, path, name, value)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path, name, value

    integer(hid_t) :: obj, atype, space, attr
    character(len=:), allocatable :: text
    integer :: err

    ! The object, the previous attribute removed.
    text = value
    if (len(text) == 0) text = ' '
    obj = open_object(file, path, name)

    ! A fixed-length string type of the text's length, null-padded; a scalar attribute.
    call h5tcopy_f(H5T_FORTRAN_S1, atype, err)
    call check(err, file, path, 'h5tcopy_f')
    call h5tset_size_f(atype, int(len(text), size_t), err)
    call check(err, file, path, 'h5tset_size_f')
    call h5tset_strpad_f(atype, H5T_STR_NULLPAD_F, err)
    call check(err, file, path, 'h5tset_strpad_f')
    call h5screate_f(H5S_SCALAR_F, space, err)
    call check(err, file, path, 'h5screate_f')
    call h5acreate_f(obj, name, atype, space, attr, err)
    call check(err, file, path, 'h5acreate_f ' // name)
    call h5awrite_f(attr, atype, text, [1_hsize_t], err)
    call check(err, file, path, 'h5awrite_f ' // name)

    ! Identifiers closed.
    call h5aclose_f(attr, err)
    call check(err, file, path, 'h5aclose_f')
    call h5sclose_f(space, err)
    call check(err, file, path, 'h5sclose_f')
    call h5tclose_f(atype, err)
    call check(err, file, path, 'h5tclose_f')
    call h5oclose_f(obj, err)
    call check(err, file, path, 'h5oclose_f')

  end subroutine h5_write_attr_text

  !> Write an integer attribute (one or more values) of a group or dataset (an existing
  !> attribute of that name is replaced).
  !>
  !> Arguments:
  !>   file    (in) type(h5_file_t)    the open file
  !>   path    (in) character          the object
  !>   name    (in) character          attribute name
  !>   values  (in) integer(ik)(:)     its values
  !> Locals:
  !>   buf     integer(ik)(:), target  contiguous copy (for c_loc)
  !>   obj, space, attr  integer(hid_t)  identifiers
  !>   err     integer                 library status
  subroutine h5_write_attr_int(file, path, name, values)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path, name
    integer(ik), intent(in) :: values(:)

    integer(ik), allocatable, target :: buf(:)
    integer(hid_t) :: obj, space, attr
    integer :: err

    buf = values
    obj = open_object(file, path, name)
    call h5screate_simple_f(1, shape(buf, kind=hsize_t), space, err)
    call check(err, file, path, 'h5screate_simple_f')
    call h5acreate_f(obj, name, h5kind_to_type(ik, H5_INTEGER_KIND), space, attr, err)
    call check(err, file, path, 'h5acreate_f ' // name)
    call h5awrite_f(attr, h5kind_to_type(ik, H5_INTEGER_KIND), c_loc(buf), err)
    call check(err, file, path, 'h5awrite_f ' // name)
    call h5aclose_f(attr, err)
    call check(err, file, path, 'h5aclose_f')
    call h5sclose_f(space, err)
    call check(err, file, path, 'h5sclose_f')
    call h5oclose_f(obj, err)
    call check(err, file, path, 'h5oclose_f')

  end subroutine h5_write_attr_int

  !> Read a text attribute.
  !>
  !> Arguments:
  !>   file   (in)  type(h5_file_t)  the open file
  !>   path   (in)  character        the object
  !>   name   (in)  character        attribute name
  !>   value  (out) character(:), allocatable  its text
  !> Locals:
  !>   obj, attr, atype  integer(hid_t)  identifiers
  !>   length integer(size_t)        the string type's length
  !>   err    integer                library status
  subroutine h5_read_attr_text(file, path, name, value)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path, name
    character(len=:), allocatable, intent(out) :: value

    integer(hid_t) :: obj, attr, atype
    integer(size_t) :: length
    integer :: err

    ! The attribute and its string length.
    call h5oopen_f(file%id, path, obj, err)
    call check(err, file, path, 'h5oopen_f')
    call h5aopen_f(obj, name, attr, err)
    call check(err, file, path, 'h5aopen_f ' // name)
    call h5aget_type_f(attr, atype, err)
    call check(err, file, path, 'h5aget_type_f')
    call h5tget_size_f(atype, length, err)
    call check(err, file, path, 'h5tget_size_f')

    ! The text, then the identifiers closed.
    allocate (character(len=length) :: value)
    call h5aread_f(attr, atype, value, [1_hsize_t], err)
    call check(err, file, path, 'h5aread_f ' // name)
    call h5tclose_f(atype, err)
    call check(err, file, path, 'h5tclose_f')
    call h5aclose_f(attr, err)
    call check(err, file, path, 'h5aclose_f')
    call h5oclose_f(obj, err)
    call check(err, file, path, 'h5oclose_f')

  end subroutine h5_read_attr_text

  !> Read an integer attribute of num values.
  !>
  !> Arguments:
  !>   file    (in)  type(h5_file_t)  the open file
  !>   path    (in)  character        the object
  !>   name    (in)  character        attribute name
  !>   values  (out) integer(ik)(:), allocatable, target  its values
  !>   num     (in)  integer          number of values expected
  !> Locals:
  !>   obj, attr  integer(hid_t)      identifiers
  !>   ptr     type(c_ptr)            address of the values (the library's buffer argument)
  !>   err     integer                library status
  subroutine h5_read_attr_int(file, path, name, values, num)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path, name
    integer(ik), allocatable, target, intent(out) :: values(:)
    integer, intent(in) :: num

    integer(hid_t) :: obj, attr
    type(c_ptr) :: ptr
    integer :: err

    allocate (values(num))
    ptr = c_loc(values)
    call h5oopen_f(file%id, path, obj, err)
    call check(err, file, path, 'h5oopen_f')
    call h5aopen_f(obj, name, attr, err)
    call check(err, file, path, 'h5aopen_f ' // name)
    call h5aread_f(attr, h5kind_to_type(ik, H5_INTEGER_KIND), ptr, err)
    call check(err, file, path, 'h5aread_f ' // name)
    call h5aclose_f(attr, err)
    call check(err, file, path, 'h5aclose_f')
    call h5oclose_f(obj, err)
    call check(err, file, path, 'h5oclose_f')

  end subroutine h5_read_attr_int

  !> Open a group or dataset to write an attribute, the existing attribute of that name
  !> removed.
  !>
  !> Arguments:
  !>   file  (in) type(h5_file_t)  the open file
  !>   path  (in) character        the object
  !>   name  (in) character        attribute name
  !>   obj   (out, result) integer(hid_t)  the open object
  !> Locals:
  !>   exists  logical             the attribute exists
  !>   err     integer             library status
  function open_object(file, path, name) result(obj)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path, name
    integer(hid_t) :: obj

    logical :: exists
    integer :: err

    call h5oopen_f(file%id, path, obj, err)
    call check(err, file, path, 'h5oopen_f')
    call h5aexists_f(obj, name, exists, err)
    call check(err, file, path, 'h5aexists_f ' // name)
    if (exists) then
      call h5adelete_f(obj, name, err)
      call check(err, file, path, 'h5adelete_f ' // name)
    end if

  end function open_object

  !> A link-creation property list that creates the missing groups of a path.
  !>
  !> Arguments:
  !>   file, path  (in) type(h5_file_t), character  for messages
  !>   lcpl        (out, result) integer(hid_t)     the list (closed by the caller)
  !> Locals:
  !>   err         integer                          library status
  function link_list(file, path) result(lcpl)

    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(hid_t) :: lcpl

    integer :: err

    call h5pcreate_f(H5P_LINK_CREATE_F, lcpl, err)
    call check(err, file, path, 'h5pcreate_f')
    call h5pset_create_inter_group_f(lcpl, 1, err)
    call check(err, file, path, 'h5pset_create_inter_group_f')

  end function link_list

  !> The HDF5 library linked, "HDF5 <major>.<minor>.<release>" (the build information line).
  !>
  !> Arguments:
  !>   text  (out, result) character(:), allocatable
  !> Locals:
  !>   major, minor, release  integer  the version numbers
  !>   err   integer                   library status
  !>   buf   character(48)             formatting buffer
  function h5_library_version() result(text)

    character(len=:), allocatable :: text

    integer :: major, minor, release, err
    character(len=48) :: buf

    call init_library()
    call h5get_libversion_f(major, minor, release, err)
    if (err < 0) call fatal_error('h5_file', 'HDF5 error in h5get_libversion_f')
    write (buf, '(a,i0,a,i0,a,i0)') 'HDF5 ', major, '.', minor, '.', release
    text = trim(buf)

  end function h5_library_version

  !> Initialize the library's Fortran interface once (h5open_f).
  !>
  !> Locals:
  !>   err  integer  library status
  subroutine init_library()

    integer :: err

    if (library_open) return
    call h5open_f(err)
    if (err < 0) call fatal_error('h5_file', 'the HDF5 library cannot be initialized (h5open_f)')
    library_open = .true.

  end subroutine init_library

  !> A negative library status is a fatal error naming the file, the path and the call.
  !>
  !> Arguments:
  !>   err   (in) integer          library status
  !>   file  (in) type(h5_file_t)  the file
  !>   path  (in) character        the path
  !>   call_name (in) character    the library call
  subroutine check(err, file, path, call_name)

    integer, intent(in) :: err
    type(h5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path, call_name

    if (err < 0) call fatal_error('h5_file', file%name // ':' // path // ': HDF5 error in ' // &
      & call_name)

  end subroutine check

end module hdf5_file_mod
