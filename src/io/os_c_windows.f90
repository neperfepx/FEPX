! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The file-system calls on Windows (64-bit; its C runtime and kernel32, e.g. with MSYS2's
!> gfortran).
!>
!> Source : FEPX. Contains no numerical method. The build compiles this file on Windows instead
!>          of os_c_posix.f90 (the same module), for dir_mod and restart_file_mod: Windows' C
!>          runtime has _mkdir(path), without the permissions of POSIX's mkdir, and its rename
!>          fails when the target exists, so a file is replaced with kernel32's MoveFileExA
!>          (MOVEFILE_REPLACE_EXISTING; Microsoft's documentation of MoveFileExA).
!>
!> Public entities:
!>   c_make_dir      create one directory (its parent must exist)
!>   c_replace_file  rename a file, replacing an existing target
!>   c_is_dir        whether a path is an existing directory (kernel32's GetFileAttributesA:
!>                   the Fortran runtimes do not all answer INQUIRE on a directory there)
!>   c_list_dir      the entries of a directory, sorted (kernel32's FindFirstFileA, FindNextFileA
!>                   and FindClose; no shell: the C library's popen runs cmd.exe, without ls)
module os_c_mod

  use, intrinsic :: iso_c_binding, only: c_int, c_char, c_null_char, c_long, c_ptr, &
    & c_intptr_t

  implicit none

  private

  public :: c_make_dir, c_replace_file, c_is_dir, c_list_dir

  !> Interface to the C runtime's _mkdir (Windows: int _mkdir(const char *path)).
  !>
  !> Arguments:
  !>   name  (in) character(kind=c_char)(*)  path, null-terminated
  !> Result: 0 on success, -1 on failure (e.g. the directory exists)
  interface
    function c_mkdir(name) bind(C, name='_mkdir')
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: name(*)
      integer(c_int) :: c_mkdir
    end function c_mkdir
  end interface

  !> INVALID_FILE_ATTRIBUTES and FILE_ATTRIBUTE_DIRECTORY (winnt.h, fileapi.h): the result of
  !> GetFileAttributesA for a path that does not exist, and the bit of a directory.
  integer(c_long), parameter :: invalid_file_attributes = -1_c_long
  integer(c_long), parameter :: file_attribute_directory = 16_c_long

  !> Interface to kernel32's GetFileAttributesA (DWORD GetFileAttributesA(LPCSTR name)).
  !>
  !> Arguments:
  !>   name  (in) character(kind=c_char)(*)  path, null-terminated
  !> Result: the attribute bits, or INVALID_FILE_ATTRIBUTES (all bits set)
  interface
    function c_get_file_attributes(name) bind(C, name='GetFileAttributesA')
      import :: c_char, c_long
      character(kind=c_char), intent(in) :: name(*)
      integer(c_long) :: c_get_file_attributes
    end function c_get_file_attributes
  end interface

  !> The size of WIN32_FIND_DATAA (minwinbase.h: the attributes, three times, the size, two
  !> reserved words, then cFileName, 260 characters at byte offset 44, and cAlternateFileName,
  !> 14 characters; 320 bytes with its alignment), and the offset and length of cFileName.
  integer, parameter :: find_data_size = 320
  integer, parameter :: find_data_name_offset = 44
  integer, parameter :: find_data_name_length = 260

  !> Interfaces to kernel32's FindFirstFileA (HANDLE FindFirstFileA(LPCSTR pattern,
  !> LPWIN32_FIND_DATAA data); INVALID_HANDLE_VALUE, -1, when nothing matches), FindNextFileA
  !> (BOOL FindNextFileA(HANDLE, LPWIN32_FIND_DATAA); zero after the last entry) and FindClose
  !> (BOOL FindClose(HANDLE)). The structure is passed as its bytes.
  interface
    function c_find_first_file(pattern, data) bind(C, name='FindFirstFileA')
      import :: c_char, c_ptr
      character(kind=c_char), intent(in) :: pattern(*)
      character(kind=c_char), intent(inout) :: data(*)
      type(c_ptr) :: c_find_first_file
    end function c_find_first_file
    function c_find_next_file(handle, data) bind(C, name='FindNextFileA')
      import :: c_char, c_ptr, c_int
      type(c_ptr), value :: handle
      character(kind=c_char), intent(inout) :: data(*)
      integer(c_int) :: c_find_next_file
    end function c_find_next_file
    function c_find_close(handle) bind(C, name='FindClose')
      import :: c_ptr, c_int
      type(c_ptr), value :: handle
      integer(c_int) :: c_find_close
    end function c_find_close
  end interface

  !> MOVEFILE_REPLACE_EXISTING (winbase.h): MoveFileExA replaces an existing target.
  integer(c_long), parameter :: movefile_replace_existing = 1_c_long

  !> Interface to kernel32's MoveFileExA (BOOL MoveFileExA(LPCSTR existing, LPCSTR new,
  !> DWORD flags); 64-bit Windows has one calling convention, that of C).
  !>
  !> Arguments:
  !>   old, new  (in) character(kind=c_char)(*)  file names, null-terminated
  !>   flags     (in) integer(c_long)            MOVEFILE_* flags (DWORD)
  !> Result: non-zero on success
  interface
    function c_move_file_ex(old, new, flags) bind(C, name='MoveFileExA')
      import :: c_int, c_char, c_long
      character(kind=c_char), intent(in) :: old(*), new(*)
      integer(c_long), value :: flags
      integer(c_int) :: c_move_file_ex
    end function c_move_file_ex
  end interface

contains

  !> Create one directory.
  !>
  !> Arguments:
  !>   path    (in)  character  the directory ('/' or '\' as separators)
  !>   status  (out, result) integer(c_int)  0 on success, -1 on failure (e.g. it exists)
  function c_make_dir(path) result(status)

    character(len=*), intent(in) :: path
    integer(c_int) :: status

    status = c_mkdir(path // c_null_char)

  end function c_make_dir

  !> Rename a file, replacing an existing target (MoveFileExA with MOVEFILE_REPLACE_EXISTING: a
  !> rename on the same volume).
  !>
  !> Arguments:
  !>   old_name, new_name  (in)  character  file names
  !>   status  (out, result) integer(c_int)  0 on success, -1 on failure (as rename)
  function c_replace_file(old_name, new_name) result(status)

    character(len=*), intent(in) :: old_name, new_name
    integer(c_int) :: status

    if (c_move_file_ex(old_name // c_null_char, new_name // c_null_char, &
      & movefile_replace_existing) /= 0_c_int) then
      status = 0_c_int
    else
      status = -1_c_int
    end if

  end function c_replace_file

  !> Whether a path is an existing directory (GetFileAttributesA: not
  !> INVALID_FILE_ATTRIBUTES, and FILE_ATTRIBUTE_DIRECTORY set).
  !>
  !> Arguments:
  !>   path    (in)  character  the path
  !>   is_dir  (out, result) logical  .true. if it is an existing directory
  !> Locals:
  !>   attrs   integer(c_long)  the attribute bits
  function c_is_dir(path) result(is_dir)

    character(len=*), intent(in) :: path
    logical :: is_dir

    integer(c_long) :: attrs

    attrs = c_get_file_attributes(trim(path) // c_null_char)
    is_dir = attrs /= invalid_file_attributes .and. iand(attrs, file_attribute_directory) /= 0

  end function c_is_dir

  !> The entries of an existing directory, without "." and "..", sorted by byte value (as ls in
  !> the C locale; Windows returns them in the order of its file system).
  !>
  !> Arguments:
  !>   path   (in)  character            the directory
  !>   names  (out) character(1024)(:), allocatable  its entries
  !>   ok     (out) logical              the directory was listed
  !> Locals:
  !>   data    character(c_char)(find_data_size)  one WIN32_FIND_DATAA
  !>   handle  type(c_ptr)                        the search
  !>   status  integer(c_int)                     FindNextFileA's and FindClose's results
  !>   entry   character(1024)                    one entry
  !>   swap    character(1024)                    an entry being sorted
  !>   n, i, j integer                            name length; indices of the sort
  subroutine c_list_dir(path, names, ok)

    character(len=*), intent(in) :: path
    character(len=1024), allocatable, intent(out) :: names(:)
    logical, intent(out) :: ok

    character(kind=c_char) :: data(find_data_size)
    type(c_ptr) :: handle
    integer(c_int) :: status
    character(len=1024) :: entry, swap
    integer :: n, i, j

    allocate (names(0))
    ok = .false.

    ! Every entry of the directory ("." and ".." included, for an existing directory).
    data = c_null_char
    handle = c_find_first_file(trim(path) // '/*' // c_null_char, data)
    if (transfer(handle, 0_c_intptr_t) == -1_c_intptr_t) return
    do
      ! The name, up to its null character.
      n = 0
      do while (n < find_data_name_length)
        if (data(find_data_name_offset + n + 1) == c_null_char) exit
        n = n + 1
      end do
      if (n > 0) then
        entry = transfer(data(find_data_name_offset + 1:find_data_name_offset + n), &
          & repeat(' ', n))
        if (entry /= '.' .and. entry /= '..') names = [character(len=1024) :: names, entry]
      end if
      ! The next one.
      data = c_null_char
      status = c_find_next_file(handle, data)
      if (status == 0_c_int) exit
    end do
    status = c_find_close(handle)
    ok = .true.

    ! Sorted by byte value (insertion sort; a directory has few entries).
    do i = 2, size(names)
      swap = names(i)
      j = i - 1
      do while (j >= 1)
        if (.not. llt(swap, names(j))) exit
        names(j + 1) = names(j)
        j = j - 1
      end do
      names(j + 1) = swap
    end do

  end subroutine c_list_dir

end module os_c_mod
