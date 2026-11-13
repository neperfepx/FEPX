! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Error cases of element_mod: each case must end the program with a fatal error (status 1)
!> and its message (checked by tests/core/check_run.cmake).
!>
!> Source : FEPX (ROADMAP step 3.4).
!>
!> Usage: test_element_errors <case>
!> Cases:
!>   unknown_type   element_init with Gmsh type 4 (a first-order tetrahedron)
!>
!> Locals:
!>   case_name  character(64)    case to run
!>   elt        type(element_t)  descriptor
program test_element_errors

  use kinds_mod, only: ik
  use element_mod, only: element_t, element_init

  implicit none

  character(len=64) :: case_name
  type(element_t) :: elt

  call get_command_argument(1, case_name)
  select case (trim(case_name))
  case ('unknown_type')
    call element_init(elt, 4_ik)
  case default
    write (*, '(a)') 'unknown case ' // trim(case_name)
    stop 2
  end select
  write (*, '(a)') 'no error raised'

end program test_element_errors
