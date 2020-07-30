! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE quadrature_mod
  !
  ! Quadrature rules.
  !
  USE IntrinsicTypesModule, RK=>REAL_KIND
  !
  IMPLICIT NONE
  !
  ! 10-node tetrahedraon 15 integration points
  ! 6-node triangle 7 integration points
  INTEGER :: iqpt_1, iqpt_2, iqpt_3
  INTEGER, PARAMETER :: nqpt = 15
  INTEGER, PARAMETER :: nqpt1 = nqpt - 1
  INTEGER, PARAMETER, PRIVATE :: MAX_Q=nqpt1
  INTEGER, PARAMETER :: MAXQP3D=nqpt, MAXQP2D=7, NDIM3=3, NDIM2=2
  INTEGER :: nqp3d=MAXQP3D, nqp2d=MAXQP2D
  REAL(RK) :: qploc(0:2, 0:MAX_Q)
  REAL(RK) :: wtqp (0:2, 0:MAX_Q)
!
  REAL(RK) :: qp3d(NDIM3, MAXQP3D), qp2d(NDIM2, MAXQP2D)
  REAL(RK) :: wt3d(MAXQP3D), wt2d(MAXQP2D)

CONTAINS
  !
  SUBROUTINE initialize()
    !
    !     Set up quadrature rule.
    !
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !     Arguments:
    !
    !     *NONE* 
    !
    !     Locals:
    !
    INTEGER i, j, k
    !
    !----------------------------------------------------------------------
    !

! ** 10-noded tetrahedral element **
    ! excerpt from Nathan's code
    ! ---
    ! quadrature point
    !
    ! Rule from P. Keast, "Moderate-Degree Tetrahedral
    ! Quadrature Formulas", Computer Meth. Appl. Mechanics Eng. 55
    ! (1986), pp 339-348.
    !
    ! original coding -- D.Mika
    !
    ! Integrates 5th order polynomial exactly on tetrahedron
    !
    !

    qploc(0,0) =  0.333333333333333333d0 
    qploc(0,1) =  0.333333333333333333d0 
    qploc(0,2) =  0.333333333333333333d0 
    qploc(0,3) =  0.0d0 
    qploc(0,4) =  0.25d0 
    qploc(0,5) =  0.909090909090909091d-1 
    qploc(0,6) =  0.909090909090909091d-1 
    qploc(0,7) =  0.909090909090909091d-1 
    qploc(0,8) =  0.727272727272727273d0 
    qploc(0,9) =  0.665501535736642813d-1 
    qploc(0,10) = 0.665501535736642813d-1 
    qploc(0,11) = 0.665501535736642813d-1 
    qploc(0,12) = 0.433449846426335728d0 
    qploc(0,13) = 0.433449846426335728d0 
    qploc(0,14) = 0.433449846426335728d0 

    qploc(1,0) =  0.333333333333333333d0      
    qploc(1,1) =  0.333333333333333333d0 
    qploc(1,2) =  0.0d0 
    qploc(1,3) =  0.333333333333333333d0 
    qploc(1,4) =  0.25d0 
    qploc(1,5) =  0.909090909090909091d-1 
    qploc(1,6) =  0.909090909090909091d-1 
    qploc(1,7) =  0.727272727272727273d0 
    qploc(1,8) =  0.909090909090909091d-1 
    qploc(1,9) =  0.665501535736642813d-1 
    qploc(1,10) = 0.433449846426335728d0 
    qploc(1,11) = 0.433449846426335728d0 
    qploc(1,12) = 0.665501535736642813d-1 
    qploc(1,13) = 0.665501535736642813d-1 
    qploc(1,14) = 0.433449846426335728d0 

    qploc(2,0) =  0.333333333333333333d0
    qploc(2,1) =  0.0d0
    qploc(2,2) =  0.333333333333333333d0
    qploc(2,3) =  0.333333333333333333d0
    qploc(2,4) =  0.25d0
    qploc(2,5) =  0.909090909090909091d-1
    qploc(2,6) =  0.727272727272727273d0
    qploc(2,7) =  0.909090909090909091d-1
    qploc(2,8) =  0.909090909090909091d-1
    qploc(2,9) =  0.433449846426335728d0 
    qploc(2,10) = 0.665501535736642813d-1
    qploc(2,11) = 0.433449846426335728d0
    qploc(2,12) = 0.665501535736642813d-1                    
    qploc(2,13) = 0.433449846426335728d0                    
    qploc(2,14) = 0.665501535736642813d-1                    
    
    ! weights

    wtqp(0,0:3)   = 0.602678571428571597d-2
    wtqp(0,4)     = 0.302836780970891856d-1
    wtqp(0,5:8)   = 0.116452490860289742d-1
    wtqp(0,9:14)  = 0.109491415613864534d-1

! ** 6-noded triangular element **

    ! quadrature points

    qp2d(1,1)=0.33333333333333d0            
    qp2d(1,2)=0.05971587178977d0            
    qp2d(1,3)=0.47014206410512d0            
    qp2d(1,4)=0.47014206410512d0            
    qp2d(1,5)=0.79742698535309d0            
    qp2d(1,6)=0.10128650732346d0            
    qp2d(1,7)=0.10128650732346d0            

    qp2d(2,1)=0.33333333333333d0
    qp2d(2,2)=0.47014206410512d0
    qp2d(2,3)=0.05971587178977d0
    qp2d(2,4)=0.47014206410512d0
    qp2d(2,5)=0.10128650732346d0
    qp2d(2,6)=0.79742698535309d0
    qp2d(2,7)=0.10128650732346d0

    ! weight

    wt2d(1) = 0.1125d0
    wt2d(2) = 0.06619707639425d0
    wt2d(3) = 0.06619707639425d0
    wt2d(4) = 0.06619707639425d0
    wt2d(5) = 0.06296959027241d0
    wt2d(6) = 0.06296959027241d0
    wt2d(7) = 0.06296959027241d0

  END SUBROUTINE initialize
  !
  !**********************************************************************
  !
END MODULE quadrature_mod
!

