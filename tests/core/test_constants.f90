! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit test of constants_mod: every constant against an independent value.
!>
!> Source : FEPX (ROADMAP step 1.1).
!>
!> Each real constant is compared with its value written as a 37-digit decimal literal (computed
!> independently with Python's decimal module at 45 digits; the literal is rounded to rk by the
!> compiler), with a relative tolerance of 2 epsilon(1.0_rk), i.e. at most about one unit in the
!> last place for a constant computed with one or two roundings. The algebraic relations used by
!> the Mandel and deviatoric forms are checked too, and the identity tensors element by element.
!>
!> Locals:
!>   tol            real(rk)       relative tolerance, 2 epsilon(1.0_rk)
!>   identity_ref   real(rk)(3,3)  3x3 identity built element by element
!>   i              integer        loop index over rows and columns
program test_constants

  use kinds_mod, only: rk
  use constants_mod, only: pi, deg_to_rad, rad_to_deg, sqrt2, inv_sqrt2, sqrt3, inv_sqrt3, &
    & sqrt6, inv_sqrt6, sqrt_3_2, sqrt_2_3, one_third, two_thirds, identity_3, identity_mandel
  use assert_mod, only: assert_close, test_summary

  implicit none

  real(rk) :: tol
  real(rk) :: identity_ref(3, 3)
  integer :: i

  tol = 2.0_rk*epsilon(1.0_rk)

  ! Constants against high-precision decimal values.
  call assert_close(pi, 3.141592653589793238462643383279502884_rk, 'pi', rel_tol=tol)
  call assert_close(deg_to_rad, 1.745329251994329576923690768488612713e-2_rk, 'deg_to_rad', &
    & rel_tol=tol)
  call assert_close(rad_to_deg, 5.729577951308232087679815481410517033e+1_rk, 'rad_to_deg', &
    & rel_tol=tol)
  call assert_close(sqrt2, 1.414213562373095048801688724209698079_rk, 'sqrt2', rel_tol=tol)
  call assert_close(inv_sqrt2, 7.071067811865475244008443621048490393e-1_rk, 'inv_sqrt2', &
    & rel_tol=tol)
  call assert_close(sqrt3, 1.732050807568877293527446341505872367_rk, 'sqrt3', rel_tol=tol)
  call assert_close(inv_sqrt3, 5.773502691896257645091487805019574556e-1_rk, 'inv_sqrt3', &
    & rel_tol=tol)
  call assert_close(sqrt6, 2.449489742783178098197284074705891392_rk, 'sqrt6', rel_tol=tol)
  call assert_close(inv_sqrt6, 4.082482904638630163662140124509818987e-1_rk, 'inv_sqrt6', &
    & rel_tol=tol)
  call assert_close(sqrt_3_2, 1.224744871391589049098642037352945696_rk, 'sqrt_3_2', &
    & rel_tol=tol)
  call assert_close(sqrt_2_3, 8.164965809277260327324280249019637973e-1_rk, 'sqrt_2_3', &
    & rel_tol=tol)
  call assert_close(one_third, 3.333333333333333333333333333333333333e-1_rk, 'one_third', &
    & rel_tol=tol)
  call assert_close(two_thirds, 6.666666666666666666666666666666666667e-1_rk, 'two_thirds', &
    & rel_tol=tol)

  ! Algebraic relations between the constants.
  call assert_close(sqrt2*sqrt2, 2.0_rk, 'sqrt2_squared', rel_tol=tol)
  call assert_close(sqrt2*inv_sqrt2, 1.0_rk, 'sqrt2_times_inverse', rel_tol=tol)
  call assert_close(sqrt3*inv_sqrt3, 1.0_rk, 'sqrt3_times_inverse', rel_tol=tol)
  call assert_close(sqrt6*inv_sqrt6, 1.0_rk, 'sqrt6_times_inverse', rel_tol=tol)
  call assert_close(sqrt2*sqrt3, sqrt6, 'sqrt2_times_sqrt3', rel_tol=tol)
  call assert_close(sqrt_3_2*sqrt_2_3, 1.0_rk, 'sqrt_3_2_times_sqrt_2_3', rel_tol=tol)
  call assert_close(deg_to_rad*rad_to_deg, 1.0_rk, 'deg_rad_round_trip', rel_tol=tol)
  call assert_close(180.0_rk*deg_to_rad, pi, 'half_turn_in_rad', rel_tol=tol)

  ! Identity tensors: 3x3 matrix element by element, Mandel vector {1,1,1,0,0,0}.
  identity_ref = 0.0_rk
  do i = 1, 3
    identity_ref(i, i) = 1.0_rk
  end do
  call assert_close(identity_3, identity_ref, 'identity_3')
  call assert_close(identity_mandel, [1.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk], &
    & 'identity_mandel')

  ! Mandel inner product preserves the tensor inner product: I:I = 3.
  call assert_close(dot_product(identity_mandel, identity_mandel), 3.0_rk, &
    & 'identity_mandel_inner_product')

  ! Report; error stop on any failure.
  call test_summary()

end program test_constants
