! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Mathematical constants used throughout FEPX.
!>
!> Source : FEPX (ROADMAP step 1.1). Contains no numerical method. pi = 4 atan(1) and the
!>          degree conversion pi/180 are those of v2 (orientation_conversion_mod.f90, `pi`,
!>          `pi_over_180`); the square-root factors are those of the Mandel and deviatoric
!>          representations of dev/conventions.md §1.1 and §1.3.
!>
!> All values are constant expressions evaluated by the compiler at real kind rk. Only
!> dimensionless mathematical constants belong here: material parameters, physical constants
!> and solver tolerances are kept with the models and methods that use them.
!>
!> Public parameters (all real(rk), dimensionless unless stated):
!>   pi               pi = 4 atan(1)
!>   deg_to_rad       pi/180, converts degrees to radians [rad/deg]
!>   rad_to_deg       180/pi, converts radians to degrees [deg/rad]
!>   sqrt2            sqrt(2), Mandel factor of shear components (conventions §1.1)
!>   inv_sqrt2        1/sqrt(2)
!>   sqrt3            sqrt(3)
!>   inv_sqrt3        1/sqrt(3)
!>   sqrt6            sqrt(6)
!>   inv_sqrt6        1/sqrt(6)
!>   sqrt_3_2         sqrt(3/2): equivalent stress sqrt(3/2)|s'|; 5-vector component 2
!>                    (conventions §1.3)
!>   sqrt_2_3         sqrt(2/3): equivalent strain (rate) sqrt(2/3)|e'|
!>   one_third        1/3
!>   two_thirds       2/3
!>   identity_3       real(rk)(3,3)  3x3 identity matrix
!>   identity_mandel  real(rk)(6)    Mandel 6-vector of the 3x3 identity, {1,1,1,0,0,0}
!>                                   (order 11 22 33 23 13 12)
module constants_mod

  use kinds_mod, only: rk

  implicit none

  private

  public :: pi, deg_to_rad, rad_to_deg
  public :: sqrt2, inv_sqrt2, sqrt3, inv_sqrt3, sqrt6, inv_sqrt6, sqrt_3_2, sqrt_2_3
  public :: one_third, two_thirds
  public :: identity_3, identity_mandel

  ! Angles.
  real(rk), parameter :: pi = 4.0_rk*atan(1.0_rk)
  real(rk), parameter :: deg_to_rad = pi/180.0_rk
  real(rk), parameter :: rad_to_deg = 180.0_rk/pi

  ! Square-root factors of the Mandel (6-vector) and deviatoric (5-vector) forms.
  real(rk), parameter :: sqrt2 = sqrt(2.0_rk)
  real(rk), parameter :: inv_sqrt2 = 1.0_rk/sqrt(2.0_rk)
  real(rk), parameter :: sqrt3 = sqrt(3.0_rk)
  real(rk), parameter :: inv_sqrt3 = 1.0_rk/sqrt(3.0_rk)
  real(rk), parameter :: sqrt6 = sqrt(6.0_rk)
  real(rk), parameter :: inv_sqrt6 = 1.0_rk/sqrt(6.0_rk)
  real(rk), parameter :: sqrt_3_2 = sqrt(1.5_rk)
  real(rk), parameter :: sqrt_2_3 = sqrt(2.0_rk/3.0_rk)

  ! Fractions.
  real(rk), parameter :: one_third = 1.0_rk/3.0_rk
  real(rk), parameter :: two_thirds = 2.0_rk/3.0_rk

  ! Identity tensor as a 3x3 matrix (filled column by column) and as a Mandel 6-vector.
  real(rk), parameter :: identity_3(3, 3) = reshape( &
    & [1.0_rk, 0.0_rk, 0.0_rk, &
    &  0.0_rk, 1.0_rk, 0.0_rk, &
    &  0.0_rk, 0.0_rk, 1.0_rk], [3, 3])
  real(rk), parameter :: identity_mandel(6) = &
    & [1.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]

end module constants_mod
