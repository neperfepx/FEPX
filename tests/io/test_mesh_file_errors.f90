! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Line-editing helpers for test_mesh_file_errors (module procedures, see tests/io).
!>
!> Public routines: read_lines, write_lines, write_text, replace_line, delete_lines,
!> insert_after.
!>
!> The program below: error cases of mesh_file_mod and text_file_mod; runs one case, named by
!> the first command-line argument; the second argument is the fixture directory
!> (tests/io/fixtures).
!>
!> Source : FEPX (ROADMAP step 2.2). Each case writes a copy of tiny.msh with one modification
!>          to case.msh in the working directory (or an auxiliary file case.aux) and loads it; the
!>          load must end with a fatal error, whose message is checked by
!>          tests/core/check_run.cmake (list in tests/CMakeLists.txt).
!>
!> Locals:
!>   case_name  character(32)     name of the case
!>   fixtures   character(512)    fixture directory
!>   lines      line_t(:)         the lines of tiny.msh, modified
!>   mesh       mesh_file_t       the parsed mesh (never complete)
!>   aux        character(16)(1)  auxiliary file name list
module mesh_case_lines_mod

  implicit none

  private

  public :: line_t, read_lines, write_lines, write_text, replace_line, delete_lines, insert_after

  !> One line of a text file.
  !>
  !> Components:
  !>   text  character(:), allocatable  the line
  type :: line_t
    character(len=:), allocatable :: text
  end type line_t

contains

  !> Read all lines of a file (lines up to 1024 characters, trailing blanks removed).
  !>
  !> Arguments:
  !>   name   (in)   character     file name
  !>   lines  (out)  line_t(:)     its lines
  !> Locals:
  !>   unit    integer         file unit
  !>   ios     integer         iostat of the reads
  !>   buffer  character(1024) one line
  !>   num     integer         number of lines
  !>   i       integer         line index
  subroutine read_lines(name, lines)

    character(len=*), intent(in) :: name
    type(line_t), allocatable, intent(out) :: lines(:)

    integer :: unit
    integer :: ios
    character(len=1024) :: buffer
    integer :: num
    integer :: i

    ! Count the lines, then read them into a list of that size (element-by-element assignment:
    ! appending with array constructors of this type lost the text with gfortran 13 -O3).
    open (newunit=unit, file=name, status='old', action='read')
    num = 0
    do
      read (unit, '(a)', iostat=ios) buffer
      if (ios /= 0) exit
      num = num + 1
    end do
    allocate (lines(num))
    rewind (unit)
    do i = 1, num
      read (unit, '(a)') buffer
      lines(i)%text = trim(buffer)
    end do
    close (unit)

  end subroutine read_lines

  !> Write lines to a file (replacing it).
  !>
  !> Arguments:
  !>   name   (in)  character  file name
  !>   lines  (in)  line_t(:)  the lines
  !> Locals:
  !>   unit  integer  file unit
  !>   i     integer  line index
  subroutine write_lines(name, lines)

    character(len=*), intent(in) :: name
    type(line_t), intent(in) :: lines(:)

    integer :: unit
    integer :: i

    open (newunit=unit, file=name, status='replace', action='write')
    do i = 1, size(lines)
      write (unit, '(a)') lines(i)%text
    end do
    close (unit)

  end subroutine write_lines

  !> Write a text (lines separated by new_line) to a file, replacing it.
  !>
  !> Arguments:
  !>   name  (in)  character  file name
  !>   text  (in)  character  the text
  !> Locals:
  !>   unit  integer  file unit
  subroutine write_text(name, text)

    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: text

    integer :: unit

    open (newunit=unit, file=name, status='replace', action='write')
    write (unit, '(a)') text
    close (unit)

  end subroutine write_text

  !> Replace the n-th line equal to old (default the first) with new; stop if not found.
  !>
  !> Arguments:
  !>   lines  (inout)         line_t(:)  the lines
  !>   old    (in)            character  line to replace (whole line)
  !>   new    (in)            character  replacement
  !>   nth    (in, optional)  integer    which occurrence (default 1)
  !> Locals:
  !>   i       integer  line index
  !>   seen    integer  occurrences found so far
  !>   target  integer  resolved nth
  subroutine replace_line(lines, old, new, nth)

    type(line_t), intent(inout) :: lines(:)
    character(len=*), intent(in) :: old
    character(len=*), intent(in) :: new
    integer, intent(in), optional :: nth

    integer :: i
    integer :: seen
    integer :: target

    target = 1
    if (present(nth)) target = nth
    seen = 0
    do i = 1, size(lines)
      if (lines(i)%text == old) then
        seen = seen + 1
        if (seen == target) then
          lines(i)%text = new
          return
        end if
      end if
    end do
    error stop 'replace_line: line not found (test setup error)'

  end subroutine replace_line

  !> Delete the lines from the first line equal to first to the next line equal to last
  !> (inclusive); stop if not found.
  !>
  !> Arguments:
  !>   lines  (inout)  line_t(:), allocatable  the lines
  !>   first  (in)     character               first line of the range
  !>   last   (in)     character               last line of the range
  !> Locals:
  !>   i        integer  line index
  !>   i_first  integer  index of the first line (0: not found yet)
  subroutine delete_lines(lines, first, last)

    type(line_t), allocatable, intent(inout) :: lines(:)
    character(len=*), intent(in) :: first
    character(len=*), intent(in) :: last

    integer :: i
    integer :: i_first

    i_first = 0
    do i = 1, size(lines)
      if (i_first == 0 .and. lines(i)%text == first) i_first = i
      if (i_first > 0 .and. lines(i)%text == last) then
        call remove_range(lines, i_first, i)
        return
      end if
    end do
    error stop 'delete_lines: range not found (test setup error)'

  end subroutine delete_lines

  !> Insert a line after the first line equal to after; stop if not found.
  !>
  !> Arguments:
  !>   lines  (inout)  line_t(:), allocatable  the lines
  !>   after  (in)     character               line after which to insert
  !>   new    (in)     character               inserted line
  !> Locals:
  !>   i  integer  line index
  subroutine insert_after(lines, after, new)

    type(line_t), allocatable, intent(inout) :: lines(:)
    character(len=*), intent(in) :: after
    character(len=*), intent(in) :: new

    integer :: i

    do i = 1, size(lines)
      if (lines(i)%text == after) then
        call insert_at(lines, i + 1, new)
        return
      end if
    end do
    error stop 'insert_after: line not found (test setup error)'

  end subroutine insert_after

  !> Remove lines first..last from a list (copies element by element).
  !>
  !> Arguments:
  !>   lines  (inout)  line_t(:), allocatable  the lines
  !>   first  (in)     integer                 first line removed
  !>   last   (in)     integer                 last line removed
  !> Locals:
  !>   kept  line_t(:), allocatable  the remaining lines
  !>   i     integer                 line index
  !>   j     integer                 index in kept
  subroutine remove_range(lines, first, last)

    type(line_t), allocatable, intent(inout) :: lines(:)
    integer, intent(in) :: first
    integer, intent(in) :: last

    type(line_t), allocatable :: kept(:)
    integer :: i
    integer :: j

    allocate (kept(size(lines) - (last - first + 1)))
    j = 0
    do i = 1, size(lines)
      if (i >= first .and. i <= last) cycle
      j = j + 1
      kept(j)%text = lines(i)%text
    end do
    call move_alloc(kept, lines)

  end subroutine remove_range

  !> Insert a line at a position of a list (copies element by element).
  !>
  !> Arguments:
  !>   lines  (inout)  line_t(:), allocatable  the lines
  !>   pos    (in)     integer                 position of the new line
  !>   new    (in)     character               the new line
  !> Locals:
  !>   grown  line_t(:), allocatable  the lines with the new one
  !>   i      integer                 line index
  subroutine insert_at(lines, pos, new)

    type(line_t), allocatable, intent(inout) :: lines(:)
    integer, intent(in) :: pos
    character(len=*), intent(in) :: new

    type(line_t), allocatable :: grown(:)
    integer :: i

    allocate (grown(size(lines) + 1))
    do i = 1, pos - 1
      grown(i)%text = lines(i)%text
    end do
    grown(pos)%text = new
    do i = pos, size(lines)
      grown(i + 1)%text = lines(i)%text
    end do
    call move_alloc(grown, lines)

  end subroutine insert_at

end module mesh_case_lines_mod

program test_mesh_file_errors

  use mesh_file_mod, only: mesh_file_t, mesh_file_load
  use mesh_case_lines_mod, only: line_t, read_lines, write_lines, write_text, replace_line, &
    & delete_lines, insert_after

  implicit none

  character(len=32) :: case_name
  character(len=512) :: fixtures
  type(line_t), allocatable :: lines(:)
  type(mesh_file_t) :: mesh
  character(len=16) :: aux(1)

  call get_command_argument(1, case_name)
  call get_command_argument(2, fixtures)
  call read_lines(trim(fixtures) // '/tiny.msh', lines)
  aux(1) = 'case.aux'

  ! Apply the case's modification.
  select case (trim(case_name))
  case ('missing_file')
    call mesh_file_load(mesh, 'nonexistent.msh')
  case ('missing_version')
    call delete_lines(lines, '$MeshVersion', '$EndMeshVersion')
  case ('version_2_2_3')
    call replace_line(lines, '2.4', '2.2.3')
  case ('version_2_3')
    call replace_line(lines, '2.4', '2.3')
  case ('format_binary')
    call replace_line(lines, '2.2 0 8', '2.2 1 8')
  case ('missing_format')
    call delete_lines(lines, '$MeshFormat', '$EndMeshFormat')
  case ('missing_elements')
    call delete_lines(lines, '$Elements', '$EndElements')
  case ('not_a_header')
    call insert_after(lines, '$EndDomain', 'stray text')
  case ('unexpected_end')
    call insert_after(lines, '$EndDomain', '$EndNodes')
  case ('node_order')
    call replace_line(lines, '2 1.0 0.0 0.0', '3 1.0 0.0 0.0')
  case ('node_values')
    call replace_line(lines, '2 1.0 0.0 0.0', '2 1.0 0.0')
  case ('node_not_number')
    call replace_line(lines, '2 1.0 0.0 0.0', '2 1.0 abc 0.0')
  case ('count_not_integer')
    call replace_line(lines, '14', '14.0')
  case ('node_count_footer')
    call replace_line(lines, '14', '13')
  case ('overflow_dofs')
    ! 3 x 800,000,000 degrees of freedom exceed huge(ik) (D28): stopped before allocating.
    call replace_line(lines, '14', '800000000')
  case ('unused_node')
    call replace_line(lines, '14', '15')
    call insert_after(lines, '14 0.5 0.5 1.0', '15 2.0 2.0 2.0')
  case ('elt_first_order')
    call replace_line(lines, '4 11 3 7 7 2 2 3 4 5 7 10 11 12 14 13', '4 4 3 7 7 2 2 3 4 5')
  case ('elt_tags')
    call replace_line(lines, '4 11 3 7 7 2 2 3 4 5 7 10 11 12 14 13', &
      & '4 11 2 7 7 2 3 4 5 7 10 11 12 14 13')
  case ('elt_id_gap')
    call replace_line(lines, '4 11 3 7 7 2 2 3 4 5 7 10 11 12 14 13', &
      & '5 11 3 7 7 2 2 3 4 5 7 10 11 12 14 13')
  case ('elt_node_range')
    call replace_line(lines, '4 11 3 7 7 2 2 3 4 5 7 10 11 12 14 13', &
      & '4 11 3 7 7 2 2 3 4 5 7 10 11 12 14 15')
  case ('elt_none')
    call replace_line(lines, '4', '2')
    call delete_lines(lines, '3 11 3 3 3 1 1 2 3 4 6 7 8 9 10 11', &
      & '4 11 3 7 7 2 2 3 4 5 7 10 11 12 14 13')
  case ('duplicate_section')
    ! A second (empty) $Periodicity right after the first (inserted in reverse order).
    call insert_after(lines, '$EndPeriodicity', '$EndPeriodicity')
    call insert_after(lines, '$EndPeriodicity', '0')
    call insert_after(lines, '$EndPeriodicity', '$Periodicity')
  case ('partial_partition')
    call replace_line(lines, '4 11 3 7 7 2 2 3 4 5 7 10 11 12 14 13', &
      & '4 11 3 7 7 0 2 3 4 5 7 10 11 12 14 13')
  case ('missing_node_parts')
    call delete_lines(lines, '$NodePartitions', '$EndNodePartitions')
  case ('node_part_range')
    call replace_line(lines, '5 2', '5 3')
  case ('nset_range')
    call replace_line(lines, '5', '15')
  case ('faset_wrong_elt')
    call replace_line(lines, '3 1 3 2 8 7 6', '4 1 3 2 8 7 6')
  case ('faset_not_3d')
    call replace_line(lines, '3 1 3 2 8 7 6', '2 1 3 2 8 7 6')
  case ('periodicity_vector')
    call replace_line(lines, '5 1 1 1 1', '5 1 0 0 0')
  case ('ori_unknown_elset')
    call replace_line(lines, '7 10.0 20.0 30.0', '5 10.0 20.0 30.0')
  case ('ori_twice')
    call replace_line(lines, '7 10.0 20.0 30.0', '3 10.0 20.0 30.0')
  case ('ori_count')
    call replace_line(lines, '2 euler-bunge:active', '1 euler-bunge:active')
    call delete_lines(lines, '3 40.0 50.0 60.0', '3 40.0 50.0 60.0')
  case ('ori_descriptor')
    call replace_line(lines, '2 euler-bunge:active', '2 euler:active')
  case ('ori_params')
    call replace_line(lines, '7 10.0 20.0 30.0', '7 10.0 20.0')
  case ('groups_entity')
    call replace_line(lines, 'elset', 'elt')
  case ('truncated')
    call delete_lines(lines, '$EndGroups', '$EndGroups')
  case ('groups_phase')
    call replace_line(lines, '7 2', '7 0')
  case ('aux_forbidden')
    call write_lines('case.msh', lines)
    call write_text('case.aux', '$Nodes' // new_line('a') // '0' // new_line('a') // '$EndNodes')
    call mesh_file_load(mesh, 'case.msh', aux)
  case ('aux_crysym')
    call write_lines('case.msh', lines)
    call write_text('case.aux', '$ElsetCrySym' // new_line('a') // 'cubic' // new_line('a') // &
      & '$EndElsetCrySym')
    call mesh_file_load(mesh, 'case.msh', aux)
  case ('aux_elt_ids')
    call write_lines('case.msh', lines)
    call write_text('case.aux', '$ElementCrss' // new_line('a') // '2 12' // new_line('a') // &
      & '1 100.0' // new_line('a') // '3 100.0' // new_line('a') // '$EndElementCrss')
    call mesh_file_load(mesh, 'case.msh', aux)
  case default
    error stop 'unknown case'
  end select

  ! Load the modified mesh: must end with a fatal error.
  call write_lines('case.msh', lines)
  call mesh_file_load(mesh, 'case.msh')
  write (*, '(a)') 'case loaded without error'
  error stop 99

end program test_mesh_file_errors
