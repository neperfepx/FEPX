! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE DimsModule
  !
  !  *** Program Unit:  Module
  !  ***    Unit Name:  DimsModule
  !  ***  Description:
  !
  !  This module contains hardwired dimensions.
  !
  !  *** Use Statements:
  !
  USE IntrinsicTypesModule, &
       &  RK=>REAL_KIND, IK=>INTEGER_KIND, LK=>LOGICAL_KIND
  !
  !  *** End:
  !
  IMPLICIT NONE
  !
  !--------------* Parameters and hardwired dimensions.
  !
  !    DIMS    : number of dimensions
  !    TVEC    : number of dimensions of deviatoric space
  !    N_PARM  : Number of crystal PARAMETERs.
  !    MAXSLIP : Maximum number of slip systems.
  !    MAX_VERT: max number of vertices of rate independent yield surface
  !
  INTEGER, PARAMETER :: N_PARM = 12,&
       &   DIMS=3,        DIMS1=2,&
       &   TVEC=5,        TVEC1=4
  REAL(RK), PARAMETER :: VSMALL = 1.0e-8
  REAL(RK), PARAMETER :: DTINY = 1.0e-15
  REAL(RK), PARAMETER :: VTINY = 1.0e-16
  REAL(RK), PARAMETER :: TOLER_STATE = 1.d-5
   
  !   
  INTEGER, PARAMETER :: ISOTROPIC_VP = 0
  INTEGER, PARAMETER :: ANISOTROPIC_VP = 1
  INTEGER, PARAMETER :: ANISOTROPIC_EVPS = 2   
  !
  INTEGER, PARAMETER :: ngrain = 1 
  INTEGER, PARAMETER :: ngrain1 = 0 
  !
  INTEGER ::     MAXSLIP, MAXSLIP1, MAX_VERT, MAX_VERT1

CONTAINS

    SUBROUTINE set_maxslip(nslip,nvert)

        INTEGER, INTENT(IN) :: nslip, nvert

        MAXSLIP=nslip
        MAXSLIP1=MAXSLIP-1
        MAX_VERT=nvert
        MAX_VERT1=MAX_VERT-1

    END SUBROUTINE set_maxslip

END MODULE DimsModule
!
