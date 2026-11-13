! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Test driver for the input validation of orientation_mod: runs one case, named by the first
!> command-line argument; each case must end with a fatal error (exit status 1) and a message
!> on standard output and standard error, checked by tests/core/check_run.cmake.
!>
!> Source : FEPX (ROADMAP step 1.3).
!>
!> Cases:
!>   bad_descriptor      descriptor code 8
!>   bad_convention      convention code 3
!>   too_few_params      euler-bunge with 2 parameters
!>   bad_rotmat          rotmat input that is not orthogonal (a 1e-3 perturbation)
!>   reflection_rotmat   rotmat input with det -1
!>   zero_axis           axis-angle with a zero axis
!>   zero_quaternion     quaternion (0, 0, 0, 0)
!>   bad_output          descriptor_from_rot_c2s with descriptor code 0
!>
!> Locals:
!>   case_name  character(64)  case to run
!>   rot        real(rk)(3,3)  rotation (result, unused)
!>   params     real(rk)(9)    descriptor parameters (result, unused)
program test_orientation_errors

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk
  use constants_mod, only: identity_3
  use orientation_mod, only: rot_c2s_from_descriptor, descriptor_from_rot_c2s, &
    & ori_des_euler_bunge, ori_des_rotmat, ori_des_axis_angle, ori_des_quaternion, &
    & ori_conv_passive

  implicit none

  character(len=64) :: case_name
  real(rk) :: rot(3, 3)
  real(rk) :: params(9)

  call get_command_argument(1, case_name)

  select case (trim(case_name))
  case ('bad_descriptor')
    rot = rot_c2s_from_descriptor(8, ori_conv_passive, [0.0_rk, 0.0_rk, 0.0_rk])
  case ('bad_convention')
    rot = rot_c2s_from_descriptor(ori_des_euler_bunge, 3, [0.0_rk, 0.0_rk, 0.0_rk])
  case ('too_few_params')
    rot = rot_c2s_from_descriptor(ori_des_euler_bunge, ori_conv_passive, [0.0_rk, 0.0_rk])
  case ('bad_rotmat')
    rot = rot_c2s_from_descriptor(ori_des_rotmat, ori_conv_passive, &
      & [1.001_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk])
  case ('reflection_rotmat')
    rot = rot_c2s_from_descriptor(ori_des_rotmat, ori_conv_passive, &
      & [-1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk])
  case ('zero_axis')
    rot = rot_c2s_from_descriptor(ori_des_axis_angle, ori_conv_passive, &
      & [0.0_rk, 0.0_rk, 0.0_rk, 30.0_rk])
  case ('zero_quaternion')
    rot = rot_c2s_from_descriptor(ori_des_quaternion, ori_conv_passive, &
      & [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk])
  case ('bad_output')
    params = descriptor_from_rot_c2s(0, ori_conv_passive, identity_3)
  case default
    write (output_unit, '(a)') 'unknown case: ' // trim(case_name)
    error stop 99
  end select

  ! Reaching this point means no error was detected.
  write (output_unit, '(a)') 'no error detected'

end program test_orientation_errors
