! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module quadrature_mod

! Module contains element quadrature rules. Considers 10 node tetrahedral
!   elements containing 15 integration points, and 6 node triangular elements
!   (for surface calculations) containing 7 integration points.

! Rule from p. Keast, "Moderate-Degree Tetrahedral Quadrature Formulas",
!   Computer Meth. Appl. Mechanics Eng. 55 (1986), pp 339-348.

! Integrates 5th order polynomial exactly on tetrahedron

! Contains subroutines:
! quadrature_init: Sets up quadrature rules

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod

  implicit none

  public

contains

  subroutine quadrature_init()

    ! Set up quadrature rule.

    !---------------------------------------------------------------------------

    implicit none

    ! Arguments:
    ! None

    ! Locals:
    ! None

    !---------------------------------------------------------------------------

    ! 3d quadrature point locations

    qploc(1, 1) = 0.333333333333333333d0
    qploc(1, 2) = 0.333333333333333333d0
    qploc(1, 3) = 0.333333333333333333d0
    qploc(1, 4) = 0.0d0
    qploc(1, 5) = 0.25d0
    qploc(1, 6) = 0.909090909090909091d-1
    qploc(1, 7) = 0.909090909090909091d-1
    qploc(1, 8) = 0.909090909090909091d-1
    qploc(1, 9) = 0.727272727272727273d0
    qploc(1, 10) = 0.665501535736642813d-1
    qploc(1, 11) = 0.665501535736642813d-1
    qploc(1, 12) = 0.665501535736642813d-1
    qploc(1, 13) = 0.433449846426335728d0
    qploc(1, 14) = 0.433449846426335728d0
    qploc(1, 15) = 0.433449846426335728d0

    qploc(2, 1) = 0.333333333333333333d0
    qploc(2, 2) = 0.333333333333333333d0
    qploc(2, 3) = 0.0d0
    qploc(2, 4) = 0.333333333333333333d0
    qploc(2, 5) = 0.25d0
    qploc(2, 6) = 0.909090909090909091d-1
    qploc(2, 7) = 0.909090909090909091d-1
    qploc(2, 8) = 0.727272727272727273d0
    qploc(2, 9) = 0.909090909090909091d-1
    qploc(2, 10) = 0.665501535736642813d-1
    qploc(2, 11) = 0.433449846426335728d0
    qploc(2, 12) = 0.433449846426335728d0
    qploc(2, 13) = 0.665501535736642813d-1
    qploc(2, 14) = 0.665501535736642813d-1
    qploc(2, 15) = 0.433449846426335728d0

    qploc(3, 1) = 0.333333333333333333d0
    qploc(3, 2) = 0.0d0
    qploc(3, 3) = 0.333333333333333333d0
    qploc(3, 4) = 0.333333333333333333d0
    qploc(3, 5) = 0.25d0
    qploc(3, 6) = 0.909090909090909091d-1
    qploc(3, 7) = 0.727272727272727273d0
    qploc(3, 8) = 0.909090909090909091d-1
    qploc(3, 9) = 0.909090909090909091d-1
    qploc(3, 10) = 0.433449846426335728d0
    qploc(3, 11) = 0.665501535736642813d-1
    qploc(3, 12) = 0.433449846426335728d0
    qploc(3, 13) = 0.665501535736642813d-1
    qploc(3, 14) = 0.433449846426335728d0
    qploc(3, 15) = 0.665501535736642813d-1

    ! 3d quadrature point weights

    wtqp(1, 1:4) = 0.602678571428571597d-2
    wtqp(1, 5) = 0.302836780970891856d-1
    wtqp(1, 6:9) = 0.116452490860289742d-1
    wtqp(1, 10:15) = 0.109491415613864534d-1

    ! 6 node triangular element (2d)

    ! 2d quadrature point locations

    qp2d(1, 1) = 0.33333333333333d0
    qp2d(1, 2) = 0.05971587178977d0
    qp2d(1, 3) = 0.47014206410512d0
    qp2d(1, 4) = 0.47014206410512d0
    qp2d(1, 5) = 0.79742698535309d0
    qp2d(1, 6) = 0.10128650732346d0
    qp2d(1, 7) = 0.10128650732346d0

    qp2d(2, 1) = 0.33333333333333d0
    qp2d(2, 2) = 0.47014206410512d0
    qp2d(2, 3) = 0.05971587178977d0
    qp2d(2, 4) = 0.47014206410512d0
    qp2d(2, 5) = 0.10128650732346d0
    qp2d(2, 6) = 0.79742698535309d0
    qp2d(2, 7) = 0.10128650732346d0

    ! 2d quadrature point weights

    wt2d(1) = 0.1125d0
    wt2d(2) = 0.06619707639425d0
    wt2d(3) = 0.06619707639425d0
    wt2d(4) = 0.06619707639425d0
    wt2d(5) = 0.06296959027241d0
    wt2d(6) = 0.06296959027241d0
    wt2d(7) = 0.06296959027241d0

  end subroutine quadrature_init

end module quadrature_mod
