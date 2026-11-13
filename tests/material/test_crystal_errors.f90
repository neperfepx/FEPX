! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Test driver for the input validation of crystal_mod::crystal_init: runs one case, named by
!> the first command-line argument; each case must end with a fatal error (exit status 1) and a
!> message on standard output and standard error (tests/core/check_run.cmake).
!>
!> Source : FEPX (ROADMAP steps 1.4 and 6.1).
!>
!> Cases:
!>   hcp_no_c_over_a   hcp without c/a
!>   bct_zero_c_over_a bct with c/a = 0
!>   bad_structure     structure code 5
!>   no_family         an empty family list
!>   bad_family        fcc family 2
!>   family_twice      hcp families 1, 3, 1
!>
!> Locals:
!>   case_name  character(64)    case to run
!>   crys       type(crystal_t)  crystal (unused result)
program test_crystal_errors

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk
  use crystal_mod, only: crystal_t, crystal_init, crystal_fcc, crystal_hcp, crystal_bct

  implicit none

  character(len=64) :: case_name
  type(crystal_t) :: crys

  call get_command_argument(1, case_name)

  select case (trim(case_name))
  case ('hcp_no_c_over_a')
    call crystal_init(crys, crystal_hcp)
  case ('bct_zero_c_over_a')
    call crystal_init(crys, crystal_bct, c_over_a=0.0_rk)
  case ('bad_structure')
    call crystal_init(crys, 5)
  case ('no_family')
    call crystal_init(crys, crystal_fcc, families=[integer ::])
  case ('bad_family')
    call crystal_init(crys, crystal_fcc, families=[2])
  case ('family_twice')
    call crystal_init(crys, crystal_hcp, c_over_a=1.587_rk, families=[1, 3, 1])
  case default
    write (output_unit, '(a)') 'unknown case: ' // trim(case_name)
    error stop 99
  end select

  ! Reaching this point means no error was detected.
  write (output_unit, '(a)') 'no error detected'

end program test_crystal_errors
