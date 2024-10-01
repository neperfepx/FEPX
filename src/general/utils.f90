! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module utils_mod

  use general_mod
  use, intrinsic :: iso_fortran_env, only: rk => real64
  implicit none

  public

contains

  !> Change an input string to all lowercase letters
! https://stackoverflow.com/questions/10759375/how-can-i-write-a-to-upper-or-to-lower-function-in-f90
  function set_lowercase(instring) result(outstring)

    ! Change an input string to all lowercase letters

    ! Arguments:
    ! instring: Input string
    ! outstring: Output string (all lowercase)
    character(*), intent(in) :: instring
    character(len(instring)) :: outstring
    ! Locals:
    ! i: Looping index (over length of instring)
    ! icap: Index for capital letter to replace (if found)
    ! capital: Capital letters
    ! lowercase: Lowercase Letters
    integer :: icap, i
    character(26), parameter :: capital = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: lowercase = 'abcdefghijklmnopqrstuvwxyz'

    !-----------------------------------------------------------------------------

    outstring = instring
    do i = 1, len_trim(instring)
      icap = index(capital, instring(i:i))
      if (icap .gt. 0) outstring(i:i) = lowercase(icap:icap)
    end do

  end function

  !> Change an input string to all uppercase letters
! https://stackoverflow.com/questions/10759375/how-can-i-write-a-to-upper-or-to-lower-function-in-f90
  function set_uppercase(instring) result(outstring)

    ! Change an input string to all uppercase letters

    ! Arguments:
    ! instring: Input string
    ! outstring: Output string (all uppercase)
    character(*), intent(in) :: instring
    character(len(instring)) :: outstring
    ! Locals:
    ! i: Looping index (over length of instring)
    ! ilow: Index for capital letter to replace (if found)
    ! capital: Capital letters
    ! uppercase: Lowercase Letters
    integer :: ilow, i
    character(26), parameter :: capital = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(26), parameter :: lowercase = 'abcdefghijklmnopqrstuvwxyz'

    !-----------------------------------------------------------------------------

    outstring = instring
    do i = 1, len_trim(instring)
      ilow = index(lowercase, instring(i:i))
      if (ilow .gt. 0) outstring(i:i) = capital(ilow:ilow)
    end do

  end function

  subroutine alloc_1d (array, beg1, end1, method)

    real(rk), allocatable, intent(inout) :: array (:)
    integer :: beg1, end1
    character(len=*), optional :: method

    if (present(method) .and. method .eq. "soft" .and. allocated(array)) then
      return
    end if

    if (allocated(array)) then
      deallocate(array)
    end if

    allocate (array(beg1:end1))

    array = 0.0d0

  end subroutine alloc_1d

  subroutine alloc_2d (array, beg1, end1, beg2, end2, method)

    real(rk), allocatable, intent(inout) :: array (:, :)
    integer :: beg1, end1, beg2, end2
    character(len=*), optional :: method

    if (present(method) .and. method .eq. "soft" .and. allocated(array)) then
      return
    end if

    if (allocated(array)) then
      deallocate(array)
    end if

    allocate (array(beg1:end1, beg2:end2))

    array = 0.0d0

  end subroutine alloc_2d

  subroutine alloc_3d (array, beg1, end1, beg2, end2, beg3, end3, method)

    real(rk), allocatable, intent(inout) :: array (:, :, :)
    integer :: beg1, end1, beg2, end2, beg3, end3
    character(len=*), optional :: method

    if (present(method) .and. method .eq. "soft" .and. allocated(array)) then
      return
    end if

    if (allocated(array)) then
      deallocate(array)
    end if

    allocate (array(beg1:end1, beg2:end2, beg3:end3))

    array = 0.0d0

  end subroutine alloc_3d

  subroutine alloc_4d (array, beg1, end1, beg2, end2, beg3, end3, beg4, end4, method)

    real(rk), allocatable, intent(inout) :: array (:, :, :, :)
    integer :: beg1, end1, beg2, end2, beg3, end3, beg4, end4
    character(len=*), optional :: method

    if (present(method) .and. method .eq. "soft" .and. allocated(array)) then
      return
    end if

    if (allocated(array)) then
      deallocate(array)
    end if

    allocate (array(beg1:end1, beg2:end2, beg3:end3, beg4:end4))

    array = 0.0d0

  end subroutine alloc_4d

  subroutine alloc_1d_int (array, beg1, end1, method)

    integer, allocatable, intent(inout) :: array (:)
    integer :: beg1, end1
    character(len=*), optional :: method

    if (present(method) .and. method .eq. "soft" .and. allocated(array)) then
      return
    end if

    if (allocated(array)) then
      deallocate(array)
    end if

    allocate (array(beg1:end1))

    array = 0

  end subroutine alloc_1d_int

  subroutine alloc_1d_logical (array, beg1, end1, method)

    logical, allocatable, intent(inout) :: array (:)
    integer :: beg1, end1
    character(len=*), optional :: method

    if (present(method) .and. method .eq. "soft" .and. allocated(array)) then
      return
    end if

    if (allocated(array)) then
      deallocate(array)
    end if

    allocate (array(beg1:end1))

    array = .false.

  end subroutine alloc_1d_logical

subroutine string_numsubstrings (myString, count)

  character(len=*), intent(in) :: myString
  integer, intent(out) :: count

  character(len=10000) :: myString2
  character(len=10000) :: delimiter
  integer :: i, length

  ! Set the delimiter to a space
  delimiter = ' '

  ! Initialize count to zero
  count = 0

  myString2 = trim (myString)

  ! Loop to find and count non-whitespace substrings
  length = len_trim(myString2)

  do i = 1, length
    if (i .eq. 1) then
      count = 1
    else if (myString2(i-1:i-1) .eq. ' ' .and. myString2(i:i) .ne. ' ') then
      count = count + 1
    end if
  end do

end subroutine string_numsubstrings

subroutine string_substring (myString, pos, Out)

  character(len=*), intent(in) :: myString
  integer, intent(in) :: pos
  character(len=*), intent(inout) :: Out

  character(len=255) :: myString2
  integer :: i, length, count, tmp

  Out = ""
  tmp = 1

  ! Initialize count to zero
  count = 0

  myString2 = trim (myString)

  ! Loop to find and count non-whitespace substrings
  length = len_trim(myString2)

  do i = 1, length
    if (i .eq. 1) then
      count = 1
    else if (myString2(i-1:i-1) .eq. ' ' .and. myString2(i:i) .ne. ' ') then
      count = count + 1
    end if

    if (count .eq. pos .and. myString2(i:i) .ne. ' ') then
      Out(tmp:tmp) = myString2(i:i)
      tmp = tmp + 1
    end if
  end do

end subroutine string_substring

subroutine string_substring_real (myString, pos, Out)

  character(len=*), intent(in) :: myString
  integer, intent(in) :: pos
  real(rk), intent(inout) :: Out

  character(len=255) :: tmp
  integer :: ioerr

  call string_substring (myString, pos, tmp)
  read (tmp, *, iostat=ioerr) Out

end subroutine string_substring_real

function ut_list_testelt (string, c, part) result(found)

  character(len=*), intent(in) :: string
  character, intent(in) :: c
  character(len=*), intent(in) :: part

  integer :: i, length1, length2
  logical :: found

  found = .false.

  ! Loop to find and count non-whitespace substrings
  length1 = len_trim(string)
  length2 = len_trim(part)
  do i = 1, length1
    if (i .eq. 1 .or. string(i:i) .eq. c) then
      if (string(i:i + length2 - 1) .eq. part) then
        if (i + length2 - 1 .eq. length1 .or. string(i+length2:i+length2) .eq. c) then
          found = .true.
          exit
        else
          found = .false.
        end if
      end if
    end if
  end do

end function ut_list_testelt

subroutine ut_dir_remove (directory)

    character(len=*), intent(in) :: directory

    character(len=:), allocatable :: command

    allocate(character(len=1000) :: command)
    write(command, '(A,A)') 'rm -rf ', trim(directory)

    call execute_command_line(command)

    deallocate(command)

end subroutine ut_dir_remove

function ut_file_exists (filename) result(file_exists)

    character(len=*), intent(in) :: filename

    logical :: file_exists

    inquire (file=filename, exist=file_exists)

end function ut_file_exists

!> This writes the file names of the user-defined output files from the
! simulation.cfg.
subroutine ut_string_addtoarray(input_string, string_array)

  !---------------------------------------------------------------------------

  ! Arguments:
  ! input_string: Character array denoting file name to print.

  character(len=*), intent(in) :: input_string
  character(len=16), intent(inout), allocatable :: string_array(:)

  character(len=16), allocatable :: temp_array(:)
  integer :: arraysize, i

  !---------------------------------------------------------------------------

  ! Append the input_string to the present array and reallocate.
  if (allocated(string_array)) then
    arraysize = size(string_array)
    allocate (temp_array(arraysize + 1))

    do i = 1, arraysize
      temp_array(i) = string_array(i)
    end do

    temp_array(arraysize + 1) = input_string
    deallocate (string_array)
    call move_alloc(temp_array, string_array)
  else
    allocate (string_array(1))
    string_array(1) = input_string
  end if

end subroutine ut_string_addtoarray

subroutine ut_string_version(version, major, minor, patch)

  !---------------------------------------------------------------------------

  character(len=10), intent(in) :: version
  integer, intent(out) :: major, minor, patch

  character(len=10) :: temp

  minor = 0
  patch = 0

  read(version, '(I1, A)') major, temp
  if (scan(temp, ".") /= 0) then
    read(temp(2:), '(I1, A)') minor, temp
    if (scan(temp, ".") /= 0) then
      read(temp(2:), '(I1)') patch
    end if
  end if

end subroutine ut_string_version

end module utils_mod
