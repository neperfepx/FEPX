! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE DIMENSIONS_MOD
!
! This module contains hardwired dimensions.
!
! Contains subroutines:
! SET_MAXSLIP: Sets the maximum number of slip systems and SCYS vertices
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND, IK=>INTEGER_KIND, &
    & LK=>LOGICAL_KIND
!
IMPLICIT NONE
!
! Parameters and hardwired dimensions:
! DIMS: Number of dimensions
! TVEC: Number of dimensions of deviatoric space
! N_PARM: Number of crystal PARAMETERs.
! MAXSLIP: Max number of slip systems.
! MAX_VERT: Max number of vertices of rate independent yield surface
! NGRAIN: Number of grains per element (legacy, always 1)
! Various numerical tolerances / limits
! Various solution flags
!
INTEGER, PARAMETER :: DIMS = 3
INTEGER, PARAMETER :: DIMS1 = 2
!
INTEGER, PARAMETER :: TVEC = 5
INTEGER, PARAMETER :: TVEC1 = 4
!
INTEGER, PARAMETER :: N_PARM = 12
!
INTEGER :: MAXSLIP
INTEGER :: MAXSLIP1
!
INTEGER :: MAX_VERT
INTEGER :: MAX_VERT1
!
INTEGER, PARAMETER :: NGRAIN = 1
INTEGER, PARAMETER :: NGRAIN1 = 0
!
REAL(RK), PARAMETER :: VSMALL = 1.0D-8
REAL(RK), PARAMETER :: DTINY = 1.0D-15
REAL(RK), PARAMETER :: VTINY = 1.0D-16
REAL(RK), PARAMETER :: TOLER_STATE = 1.0D-5
!
INTEGER, PARAMETER :: ISOTROPIC_VP = 0
INTEGER, PARAMETER :: ANISOTROPIC_VP = 1
INTEGER, PARAMETER :: ANISOTROPIC_EVPS = 2
!
CONTAINS
    !
    SUBROUTINE SET_MAXSLIP(NSLIP, NVERT)
    !
    ! Sets the maximum number of slip systems and SCYS vertices
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NSLIP: Number of slip systems
    ! NVERT: Number of SCYS vertices
    !
    INTEGER, INTENT(IN) :: NSLIP
    INTEGER, INTENT(IN) :: NVERT
    !
    !---------------------------------------------------------------------------
    !
    !
    MAXSLIP = NSLIP
    MAXSLIP1 = MAXSLIP - 1
    MAX_VERT = NVERT
    MAX_VERT1 = MAX_VERT - 1
    !
    END SUBROUTINE SET_MAXSLIP
    !
END MODULE DIMENSIONS_MOD
