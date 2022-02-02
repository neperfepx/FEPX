! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE CRYSTAL_TYPE_MOD
!
! This module handles the basic CrystalType object.
!
! To do:
!   - Add VERTICES to object structure
!   - Add P*P^T to "get list"
!
! Contains subroutines:
! SCHMIDTENSORS: Create Schmid tensors for specified crystal type
! CRYSTALTYPEGET: Return deviatoric parts of Schmid tensors
!
! Contains functions:
! CRYSTALTYPECREATE: Create a crystal type object
!
! From libf95:
!
USE FILES_MOD
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND, IK=>INTEGER_KIND, &
    & LK=>LOGICAL_KIND
!
! From libfepx:
!
USE MATRIX_OPERATIONS_MOD
!
IMPLICIT NONE
!
PRIVATE
!
INTEGER, PARAMETER :: CLASS_FCC = 1
INTEGER, PARAMETER :: CLASS_BCC = 2
INTEGER, PARAMETER :: CLASS_HCP = 3
INTEGER, PARAMETER :: CLASS_BCT = 4
INTEGER :: DECOMP_DFLT = DECOMP_MPSIM
!
TYPE CRYSTALTYPETYPE
    !
    !PRIVATE  ! allowing access to components
    !
    !CHARACTER, POINTER :: name(:)
    !
    INTEGER :: CLASS   ! FCC, BCC, HCP
    INTEGER :: DECOMP  ! DECOMPosition convention
    INTEGER :: NUMSLIP
    INTEGER :: NUMVERTICES
    !
    REAL(RK), POINTER :: SCHMID_3X3(:, :, :)
    REAL(RK), POINTER :: DEV(:, :)
    REAL(RK), POINTER :: SKW(:, :)
    REAL(RK), POINTER :: PPTRANS(:, :, :)
    REAL(RK), POINTER :: VERTICES(:, :)
    REAL(RK), POINTER :: VERTICES3X3(:, :, :)
    !
END TYPE CRYSTALTYPETYPE
!
! Public
!
PUBLIC :: CRYSTALTYPETYPE
PUBLIC :: CLASS_FCC
PUBLIC :: CLASS_BCC
PUBLIC :: CLASS_HCP
PUBLIC :: CLASS_BCT
PUBLIC :: DECOMP_MPSIM
PUBLIC :: DECOMP_FEMEVPS
!
PUBLIC :: CRYSTALTYPECREATE
PUBLIC :: CRYSTALTYPEGET
!
CONTAINS
    !
    FUNCTION CRYSTALTYPECREATE(CTYPE, C_OVER_A, HRATIO_HCP_PYR, &
        & HRATIO_HCP_PRISM, HRATIO_BCT_A, HRATIO_BCT_B, HRATIO_BCT_C, &
        & HRATIO_BCT_D, HRATIO_BCT_E, HRATIO_BCT_F, HRATIO_BCT_G, &
        & HRATIO_BCT_H, HRATIO_BCT_I, DECOMP) RESULT(SELF)
    !
    ! Create a crystal type object
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! CTYPE: Crystal type
    ! C_OVER_A: Hexagonal c/a ratio
    ! HRATIO_HCP: Ratio of pyramidal strengths to basal/prismatic
    ! DECOMP: Decomposition convention indicator
    ! SELF: (Result)
    !
    INTEGER, INTENT(IN) :: CTYPE
    REAL(RK), INTENT(IN), OPTIONAL :: C_OVER_A
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_HCP_PYR
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_HCP_PRISM
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_A
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_B
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_C
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_D
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_E
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_F
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_G
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_H
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_I
    INTEGER, INTENT(IN), OPTIONAL :: DECOMP
    TYPE(CRYSTALTYPETYPE), POINTER :: SELF
    !
    ! Locals:
    !
    INTEGER(IK) :: MYSTAT
    !
    REAL(RK), PARAMETER :: DFLT_HRATIO_HCP_PYR = 1.0D0
    REAL(RK), PARAMETER :: DFLT_HRATIO_HCP_PRISM = 1.0D0
    REAL(RK), PARAMETER :: DFLT_CRATIO = 1.0D0
    !
    REAL(RK), PARAMETER :: DFLT_HRATIO_BCT_A = 1.0D0
    REAL(RK), PARAMETER :: DFLT_HRATIO_BCT_B = 1.0D0
    REAL(RK), PARAMETER :: DFLT_HRATIO_BCT_C = 1.0D0
    REAL(RK), PARAMETER :: DFLT_HRATIO_BCT_D = 1.0D0
    REAL(RK), PARAMETER :: DFLT_HRATIO_BCT_E = 1.0D0
    REAL(RK), PARAMETER :: DFLT_HRATIO_BCT_F = 1.0D0
    REAL(RK), PARAMETER :: DFLT_HRATIO_BCT_G = 1.0D0
    REAL(RK), PARAMETER :: DFLT_HRATIO_BCT_H = 1.0D0
    REAL(RK), PARAMETER :: DFLT_HRATIO_BCT_I = 1.0D0
    !
    REAL(RK) :: HRHCP_PYR
    REAL(RK) :: HRHCP_PRISM
    REAL(RK) :: CRATIO
    !
    REAL(RK) :: HRBCT_A
    REAL(RK) :: HRBCT_B
    REAL(RK) :: HRBCT_C
    REAL(RK) :: HRBCT_D
    REAL(RK) :: HRBCT_E
    REAL(RK) :: HRBCT_F
    REAL(RK) :: HRBCT_G
    REAL(RK) :: HRBCT_H
    REAL(RK) :: HRBCT_I
    !
    !---------------------------------------------------------------------------
    !
    ALLOCATE(SELF, STAT = MYSTAT)
    !
    NULLIFY(SELF%SCHMID_3X3, SELF%DEV, SELF%SKW, SELF%PPTRANS, SELF%VERTICES, &
        & SELF%VERTICES3X3)
    !
    SELF%CLASS = CTYPE
    !
    ! Expression is a scalar of type integer, character, or logical.
    !
    SELECT CASE(CTYPE)
    !
    CASE (CLASS_FCC, CLASS_BCC)
        !
        CALL SCHMIDTENSORS(CTYPE, SELF%SCHMID_3X3)
        !
    CASE (CLASS_HCP)
        !
        IF (PRESENT(C_OVER_A) ) THEN
            !
            CRATIO = C_OVER_A
            !
        ELSE
            !
            CRATIO = DFLT_CRATIO
            !
        END IF
        !
        IF (PRESENT(HRATIO_HCP_PYR) ) THEN
            !
            HRHCP_PYR = HRATIO_HCP_PYR
            !
        ELSE
            !
            HRHCP_PYR = DFLT_HRATIO_HCP_PYR
            !
        END IF
        !
        IF (PRESENT(HRATIO_HCP_PRISM) ) THEN
            !
            HRHCP_PRISM = HRATIO_HCP_PRISM
            !
        ELSE
            !
            HRHCP_PRISM = DFLT_HRATIO_HCP_PRISM
            !
        END IF
        !
        CALL SCHMIDTENSORS(CTYPE, SELF%SCHMID_3X3, C_OVER_A = CRATIO, &
            & HRATIO_HCP_PYR = HRHCP_PYR, HRATIO_HCP_PRISM = HRHCP_PRISM)
        !
    CASE (CLASS_BCT)
        !
        IF (PRESENT(C_OVER_A) ) THEN
            !
            CRATIO = C_OVER_A
            !
        ELSE
            !
            CRATIO = DFLT_CRATIO
            !
        END IF
        !
        IF (PRESENT(HRATIO_BCT_A) ) THEN
            !
            HRBCT_A = HRATIO_BCT_A
            !
        ELSE
            !
            HRBCT_A = DFLT_HRATIO_BCT_A
            !
        END IF
        !
        IF (PRESENT(HRATIO_BCT_B) ) THEN
            !
            HRBCT_B = HRATIO_BCT_B
            !
        ELSE
            !
            HRBCT_B = DFLT_HRATIO_BCT_B
            !
        END IF
        !
        IF (PRESENT(HRATIO_BCT_C) ) THEN
            !
            HRBCT_C = HRATIO_BCT_C
            !
        ELSE
            !
            HRBCT_C = DFLT_HRATIO_BCT_C
            !
        END IF
        !
        IF (PRESENT(HRATIO_BCT_D) ) THEN
            !
            HRBCT_D = HRATIO_BCT_D
            !
        ELSE
            !
            HRBCT_D = DFLT_HRATIO_BCT_D
            !
        END IF
        !
        IF (PRESENT(HRATIO_BCT_E) ) THEN
            !
            HRBCT_E = HRATIO_BCT_E
            !
        ELSE
            !
            HRBCT_E = DFLT_HRATIO_BCT_E
            !
        END IF
        !
        IF (PRESENT(HRATIO_BCT_F) ) THEN
            !
            HRBCT_F = HRATIO_BCT_F
            !
        ELSE
            !
            HRBCT_F = DFLT_HRATIO_BCT_F
            !
        END IF
        !
        IF (PRESENT(HRATIO_BCT_G) ) THEN
            !
            HRBCT_G = HRATIO_BCT_G
            !
        ELSE
            !
            HRBCT_G = DFLT_HRATIO_BCT_G
            !
        END IF
        !
        IF (PRESENT(HRATIO_BCT_H) ) THEN
            !
            HRBCT_H = HRATIO_BCT_H
            !
        ELSE
            !
            HRBCT_H = DFLT_HRATIO_BCT_H
            !
        END IF
        !
        IF (PRESENT(HRATIO_BCT_I) ) THEN
            !
            HRBCT_I = HRATIO_BCT_I
            !
        ELSE
            !
            HRBCT_I = DFLT_HRATIO_BCT_I
            !
        END IF
        !
        CALL SCHMIDTENSORS(CTYPE, SELF%SCHMID_3X3, C_OVER_A = CRATIO, &
            & HRATIO_BCT_A = HRBCT_A, HRATIO_BCT_B = HRBCT_B, &
            & HRATIO_BCT_C = HRBCT_C, HRATIO_BCT_D = HRBCT_D, &
            & HRATIO_BCT_E = HRBCT_E, HRATIO_BCT_F = HRBCT_F, &
            & HRATIO_BCT_G = HRBCT_G, HRATIO_BCT_H = HRBCT_H, &
            & HRATIO_BCT_I = HRBCT_I)
        !
    CASE DEFAULT
    !
    END SELECT
    !
    SELF%NUMSLIP = SIZE(SELF%SCHMID_3X3, DIM = 3)
    !
    IF (PRESENT(DECOMP)) THEN
        !
        SELF%DECOMP = DECOMP
        !
    ELSE
        !
        SELF%DECOMP = DECOMP_DFLT
        !
    END IF
    !
    CALL CRYSTALTYPEGET(SELF, DEV = SELF%DEV, SKW = SELF%SKW, &
        & PPTRANS = SELF%PPTRANS)
    !
    ! Retrieve the single-crystal yield surface VERTICES for crystal type.
    ! Note: FCC/BCC should be the same. These values are presently hard-coded
    ! from the values in the legacy VERTICES files from DPLab. - JC
    !
    IF (CTYPE .EQ. 1) THEN ! CLASS_FCC
        !
        ! Initialize SELF values for FCC.
        SELF%NUMVERTICES = 28
        ALLOCATE(SELF%VERTICES3X3(3, 3, SELF%NUMVERTICES))
        !
        ! Return the VERTICES 3X3s from data storage subroutine.
        CALL GET_VERTICES(1, SELF%VERTICES3X3)
        !
        ! Allocate array for tensor DECOMPosition.
        ALLOCATE(SELF%VERTICES(5, SELF%NUMVERTICES))
        CALL TENSOR3DDECOMPOSE(SELF%VERTICES3X3, DEV = SELF%VERTICES, &
            & DECOMP = SELF%DECOMP)
        !
    ELSE IF (CTYPE .EQ. 2) THEN ! CLASS_BCC
        !
        ! Initialize SELF values for BCC.
        SELF%NUMVERTICES = 28
        ALLOCATE(SELF%VERTICES3X3(3, 3, SELF%NUMVERTICES))
        !
        ! Return the VERTICES 3X3s from data storage subroutine.
        CALL GET_VERTICES(2, SELF%VERTICES3X3)
        !
        ! Allocate array for tensor DECOMPosition.
        ALLOCATE(SELF%VERTICES(5, SELF%NUMVERTICES))
        CALL TENSOR3DDECOMPOSE(SELF%VERTICES3X3, DEV = SELF%VERTICES, &
            & DECOMP = SELF%DECOMP)
        !
    ELSE IF ((CTYPE .EQ. 3) .OR. (CTYPE .EQ. 4)) THEN ! CLASS_HCP, CLASS_BCT
        !
        ! Initialize SELF values for HCP.
        SELF%NUMVERTICES = 240
        ALLOCATE(SELF%VERTICES3X3(3, 3, SELF%NUMVERTICES))
        !
        ! Return the VERTICES 3X3s from data storage subroutine.
        CALL GET_VERTICES(3, SELF%VERTICES3X3)
        !
        ! Allocate array for tensor DECOMPosition.
        ALLOCATE(SELF%VERTICES(5, SELF%NUMVERTICES))
        CALL TENSOR3DDECOMPOSE(SELF%VERTICES3X3, DEV = SELF%VERTICES, &
            & DECOMP = SELF%DECOMP)
        !
    ELSE
        !
        SELF%NUMVERTICES = 0
        !
    END IF
    !
    END FUNCTION CRYSTALTYPECREATE
    !
    !===========================================================================
    !
    SUBROUTINE SCHMIDTENSORS(CTYPE, SCHMID, C_OVER_A, HRATIO_HCP_PYR, &
        & HRATIO_HCP_PRISM, HRATIO_BCT_A, HRATIO_BCT_B, HRATIO_BCT_C, &
        & HRATIO_BCT_D, HRATIO_BCT_E, HRATIO_BCT_F, HRATIO_BCT_G, &
        & HRATIO_BCT_H, HRATIO_BCT_I)
    !
    ! Create Schmid tensors for specified crystal type
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! CTYPE: Crystal type
    ! SCHMID: Schmid tensors
    ! C_OVER_A: Hexagonal c/a ratio
    ! HRATIO_HCP: Ratio of pyramidal strengths to basal/prismatic
    ! HRATIO_BCT: Ratios of BCT slip strengths
    !
    INTEGER, INTENT(IN) :: CTYPE
    REAL(RK), POINTER :: SCHMID(:, :, :)
    REAL(RK), INTENT(IN), OPTIONAL :: C_OVER_A
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_HCP_PYR
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_HCP_PRISM
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_A
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_B
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_C
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_D
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_E
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_F
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_G
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_H
    REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_BCT_I
    !
    ! Locals:
    !
    INTEGER :: I
    !
    ! Cubic data
    !
    REAL(RK) :: CUB_110_C(3, 12) = 0.0D0
    REAL(RK) :: CUB_111_C(3, 12) = 0.0D0
    !
    ! 110, Miller notation
    !
    REAL(RK), PARAMETER, DIMENSION(3, 12) :: CUB_110 = RESHAPE(SOURCE = (/ &
        & 0.0D0,  1.0D0, -1.0D0, &
        & 1.0D0,  0.0D0, -1.0D0, &
        & 1.0D0, -1.0D0,  0.0D0, &
        & 0.0D0,  1.0D0,  1.0D0, &
        & 1.0D0,  0.0D0,  1.0D0, &
        & 1.0D0, -1.0D0,  0.0D0, &
        & 0.0D0,  1.0D0,  1.0D0, &
        & 1.0D0,  0.0D0, -1.0D0, &
        & 1.0D0,  1.0D0,  0.0D0, &
        & 0.0D0,  1.0D0, -1.0D0, &
        & 1.0D0,  0.0D0,  1.0D0, &
        & 1.0D0,  1.0D0,  0.0D0  &
        & /), SHAPE = (/3, 12/))
    !
    ! 111, Miller notation
    !
    REAL(RK), PARAMETER, DIMENSION(3, 12) :: CUB_111 = RESHAPE(SOURCE = (/ &
        & 1.0D0,  1.0D0,  1.0D0, &
        & 1.0D0,  1.0D0,  1.0D0, &
        & 1.0D0,  1.0D0,  1.0D0, &
        & 1.0D0,  1.0D0, -1.0D0, &
        & 1.0D0,  1.0D0, -1.0D0, &
        & 1.0D0,  1.0D0, -1.0D0, &
        & 1.0D0, -1.0D0,  1.0D0, &
        & 1.0D0, -1.0D0,  1.0D0, &
        & 1.0D0, -1.0D0,  1.0D0, &
        & 1.0D0, -1.0D0, -1.0D0, &
        & 1.0D0, -1.0D0, -1.0D0, &
        & 1.0D0, -1.0D0, -1.0D0  &
        & /), SHAPE = (/3, 12/))
    !
    ! 211, Miller notation (higher order, not currently utilized)
    !
    REAL(RK), PARAMETER, DIMENSION(3, 12) :: CUB_211 = RESHAPE(SOURCE = (/ &
        & -1.0D0, -1.0D0,  2.0D0, &
        & -1.0D0,  2.0D0, -1.0D0, &
        &  2.0D0, -1.0D0, -1.0D0, &
        & -1.0D0, -1.0D0, -2.0D0, &
        & -1.0D0,  2.0D0,  1.0D0, &
        &  2.0D0, -1.0D0,  1.0D0, &
        & -1.0D0,  1.0D0,  2.0D0, &
        & -1.0D0, -2.0D0, -1.0D0, &
        &  2.0D0,  1.0D0, -1.0D0, &
        & -1.0D0,  1.0D0, -2.0D0, &
        & -1.0D0, -2.0D0,  1.0D0, &
        &  2.0D0,  1.0D0,  1.0D0  &
        & /), SHAPE = (/3, 12/))
    !
    ! 123 (A), Miller notation (higher order, not currently utilized)
    !
    REAL(RK), PARAMETER, DIMENSION(3, 12) :: CUB_123A = RESHAPE(SOURCE = (/ &
        & -1.0D0, -2.0D0,  3.0D0, &
        & -1.0D0,  3.0D0, -2.0D0, &
        &  3.0D0, -1.0D0, -2.0D0, &
        & -1.0D0, -2.0D0, -3.0D0, &
        & -1.0D0,  3.0D0,  2.0D0, &
        &  3.0D0, -1.0D0,  2.0D0, &
        & -1.0D0,  2.0D0,  3.0D0, &
        & -1.0D0, -3.0D0, -2.0D0, &
        &  3.0D0,  1.0D0, -2.0D0, &
        & -1.0D0,  2.0D0, -3.0D0, &
        & -1.0D0, -3.0D0,  2.0D0, &
        &  3.0D0,  1.0D0,  2.0D0  &
        & /), SHAPE = (/3, 12/))
    !
    ! 123 (B), Miller notation (higher order, not currently utilized)
    !
    REAL(RK), PARAMETER, DIMENSION(3, 12) :: CUB_123B = RESHAPE(SOURCE = (/ &
        & -2.0D0, -1.0D0,  3.0D0, &
        & -2.0D0,  3.0D0, -1.0D0, &
        &  3.0D0, -2.0D0, -1.0D0, &
        & -2.0D0, -1.0D0, -3.0D0, &
        & -2.0D0,  3.0D0,  1.0D0, &
        &  3.0D0, -2.0D0,  1.0D0, &
        & -2.0D0,  1.0D0,  3.0D0, &
        & -2.0D0, -3.0D0, -1.0D0, &
        &  3.0D0,  2.0D0, -1.0D0, &
        & -2.0D0,  1.0D0, -3.0D0, &
        & -2.0D0, -3.0D0,  1.0D0, &
        &  3.0D0,  2.0D0,  1.0D0  &
        & /), SHAPE = (/3, 12/))
    !
    ! HCP data
    !
    REAL(RK) :: HCP_N_C(3, 18) = 0.0D0
    REAL(RK) :: HCP_D_C(3, 18) = 0.0D0
    REAL(RK) :: SCALE = 0.0D0
    !
    ! Plane normals, Miller-Bravais notation
    !
    REAL(RK), PARAMETER, DIMENSION(4, 18) :: HCP_N = RESHAPE(SOURCE = (/ &
        &  0.0D0,  0.0D0,  0.0D0,  1.0D0, &
        &  0.0D0,  0.0D0,  0.0D0,  1.0D0, &
        &  0.0D0,  0.0D0,  0.0D0,  1.0D0, &
        &  0.0D0,  1.0D0, -1.0D0,  0.0D0, &
        & -1.0D0,  0.0D0,  1.0D0,  0.0D0, &
        &  1.0D0, -1.0D0,  0.0D0,  0.0D0, &
        &  1.0D0,  0.0D0, -1.0D0,  1.0D0, &
        &  1.0D0,  0.0D0, -1.0D0,  1.0D0, &
        &  0.0D0,  1.0D0, -1.0D0,  1.0D0, &
        &  0.0D0,  1.0D0, -1.0D0,  1.0D0, &
        & -1.0D0,  1.0D0,  0.0D0,  1.0D0, &
        & -1.0D0,  1.0D0,  0.0D0,  1.0D0, &
        & -1.0D0,  0.0D0,  1.0D0,  1.0D0, &
        & -1.0D0,  0.0D0,  1.0D0,  1.0D0, &
        &  0.0D0, -1.0D0,  1.0D0,  1.0D0, &
        &  0.0D0, -1.0D0,  1.0D0,  1.0D0, &
        &  1.0D0, -1.0D0,  0.0D0,  1.0D0, &
        &  1.0D0, -1.0D0,  0.0D0,  1.0D0  &
        & /), SHAPE = (/4, 18/))
    !
    ! Slip directions, Miller-Bravais notation
    !
    REAL(RK), PARAMETER, DIMENSION(4, 18) :: HCP_D = RESHAPE(SOURCE = (/ &
        &  2.0D0, -1.0D0, -1.0D0,  0.0D0, &
        & -1.0D0,  2.0D0, -1.0D0,  0.0D0, &
        & -1.0D0, -1.0D0,  2.0D0,  0.0D0, &
        &  2.0D0, -1.0D0, -1.0D0,  0.0D0, &
        & -1.0D0,  2.0D0, -1.0D0,  0.0D0, &
        & -1.0D0, -1.0D0,  2.0D0,  0.0D0, &
        & -2.0D0,  1.0D0,  1.0D0,  3.0D0, &
        & -1.0D0, -1.0D0,  2.0D0,  3.0D0, &
        & -1.0D0, -1.0D0,  2.0D0,  3.0D0, &
        &  1.0D0, -2.0D0,  1.0D0,  3.0D0, &
        &  1.0D0, -2.0D0,  1.0D0,  3.0D0, &
        &  2.0D0, -1.0D0, -1.0D0,  3.0D0, &
        &  2.0D0, -1.0D0, -1.0D0,  3.0D0, &
        &  1.0D0,  1.0D0, -2.0D0,  3.0D0, &
        &  1.0D0,  1.0D0, -2.0D0,  3.0D0, &
        & -1.0D0,  2.0D0, -1.0D0,  3.0D0, &
        & -1.0D0,  2.0D0, -1.0D0,  3.0D0, &
        & -2.0D0,  1.0D0,  1.0D0,  3.0D0  &
        & /),SHAPE = (/4, 18/))
    !
    ! BCT data
    !
    REAL(RK) :: BCT_N_C(3, 32) = 0.0D0
    REAL(RK) :: BCT_D_C(3, 32) = 0.0D0
    !
    ! Plane normals, Miller notation
    !
    REAL(RK), PARAMETER, DIMENSION(3, 32) :: BCT_SN = RESHAPE(SOURCE = (/ &
        &  1.0D0,  0.0D0,  0.0D0, &
        &  0.0D0,  1.0D0,  0.0D0, &
        &  1.0D0,  1.0D0,  0.0D0, &
        &  1.0D0, -1.0D0,  0.0D0, &
        &  1.0D0,  0.0D0,  0.0D0, &
        &  0.0D0,  1.0D0,  0.0D0, &
        &  1.0D0,  1.0D0,  0.0D0, &
        &  1.0D0,  1.0D0,  0.0D0, &
        &  1.0D0, -1.0D0,  0.0D0, &
        &  1.0D0, -1.0D0,  0.0D0, &
        &  1.0D0,  1.0D0,  0.0D0, &
        &  1.0D0, -1.0D0,  0.0D0, &
        &  1.0D0,  0.0D0,  0.0D0, &
        &  1.0D0,  0.0D0,  0.0D0, &
        &  0.0D0,  1.0D0,  0.0D0, &
        &  0.0D0,  1.0D0,  0.0D0, &
        &  0.0D0,  0.0D0,  1.0D0, &
        &  0.0D0,  0.0D0,  1.0D0, &
        &  0.0D0,  0.0D0,  1.0D0, &
        &  0.0D0,  0.0D0,  1.0D0, &
        &  1.0D0,  0.0D0,  1.0D0, &
        &  1.0D0,  0.0D0, -1.0D0, &
        &  0.0D0,  1.0D0,  1.0D0, &
        &  0.0D0,  1.0D0, -1.0D0, &
        &  1.0D0,  2.0D0,  1.0D0, &
        & -1.0D0,  2.0D0,  1.0D0, &
        & -1.0D0, -2.0D0,  1.0D0, &
        &  1.0D0, -2.0D0,  1.0D0, &
        &  2.0D0,  1.0D0,  1.0D0, &
        & -2.0D0,  1.0D0,  1.0D0, &
        & -2.0D0, -1.0D0,  1.0D0, &
        &  2.0D0, -1.0D0,  1.0D0 &
        & /), SHAPE = (/3,  32/))
    !
    ! Slip directions, Miller notation
    !
    REAL(RK), PARAMETER, DIMENSION(3, 32) :: BCT_SD = RESHAPE(SOURCE = (/ &
        &  0.0D0,  0.0D0,  1.0D0, &
        &  0.0D0,  0.0D0,  1.0D0, &
        &  0.0D0,  0.0D0,  1.0D0, &
        &  0.0D0,  0.0D0,  1.0D0, &
        &  0.0D0,  1.0D0,  0.0D0, &
        &  1.0D0,  0.0D0,  0.0D0, &
        &  1.0D0, -1.0D0,  1.0D0, &
        & -1.0D0,  1.0D0,  1.0D0, &
        &  1.0D0,  1.0D0,  1.0D0, &
        & -1.0D0, -1.0D0,  1.0D0, &
        & -1.0D0,  1.0D0,  0.0D0, &
        &  1.0D0,  1.0D0,  0.0D0, &
        &  0.0D0,  1.0D0,  1.0D0, &
        &  0.0D0,  1.0D0, -1.0D0, &
        &  1.0D0,  0.0D0,  1.0D0, &
        &  1.0D0,  0.0D0, -1.0D0, &
        &  1.0D0,  0.0D0,  0.0D0, &
        &  0.0D0,  1.0D0,  0.0D0, &
        &  1.0D0,  1.0D0,  0.0D0, &
        &  1.0D0, -1.0D0,  0.0D0, &
        &  1.0D0,  0.0D0, -1.0D0, &
        &  1.0D0,  0.0D0,  1.0D0, &
        &  0.0D0,  1.0D0, -1.0D0, &
        &  0.0D0,  1.0D0,  1.0D0, &
        & -1.0D0,  0.0D0,  1.0D0, &
        &  1.0D0,  0.0D0,  1.0D0, &
        &  1.0D0,  0.0D0,  1.0D0, &
        & -1.0D0,  0.0D0,  1.0D0, &
        &  0.0D0, -1.0D0,  1.0D0, &
        &  0.0D0, -1.0D0,  1.0D0, &
        &  0.0D0,  1.0D0,  1.0D0, &
        &  0.0D0,  1.0D0,  1.0D0 &
        & /), SHAPE = (/3,  32/))
    !
    ! Construct Schmid tensors
    !
    SELECT CASE(CTYPE)
    !
    CASE (CLASS_FCC)
        !
        ALLOCATE(SCHMID(3, 3, 12))
        SCHMID = 0.0D0
        !
        ! Normalize Cartesian vectors
        !
        DO I = 1, 12
            !
            CUB_110_C(:, I) = CUB_110(:, I) / NORM2(CUB_110(:, I))
            CUB_111_C(:, I) = CUB_111(:, I) / NORM2(CUB_111(:, I))
            !
            SCHMID(:, :, I) = MATMUL(RESHAPE(CUB_110_C(:, I), (/3, 1/)), &
                & RESHAPE(CUB_111_C(:, I), (/1, 3/)))
            !
        END DO
        !
    CASE (CLASS_BCC)
        !
        ALLOCATE(SCHMID(3, 3, 12))
        SCHMID = 0.0D0
        !
        ! Normalize Cartesian vectors
        !
        DO I = 1, 12
            !
            CUB_110_C(:, I) = CUB_110(:, I) / NORM2(CUB_110(:, I))
            CUB_111_C(:, I) = CUB_111(:, I) / NORM2(CUB_111(:, I))
            !
            SCHMID(:, :, I) = MATMUL(RESHAPE(CUB_111_C(:, I), (/3, 1/)), &
                & RESHAPE(CUB_110_C(:, I), (/1, 3/)))
            !
        END DO
        !
    CASE (CLASS_HCP)
        !
        ALLOCATE(SCHMID(3, 3, 18))
        SCHMID = 0.0D0
        !
        ! Convert Miller-Bravais notation to Cartesian (orthonormal)
        !
        HCP_N_C(1, :) = HCP_N(1, :)
        HCP_N_C(2, :) = (2.0D0 * HCP_N(2, :) + HCP_N(1, :)) / DSQRT(3.0D0)
        HCP_N_C(3, :) = HCP_N(4, :) / C_OVER_A
        !
        HCP_D_C(1, :) = 1.5D0 * HCP_D(1, :)
        HCP_D_C(2, :) = (HCP_D(2, :) + 0.5D0 * HCP_D(1, :)) * DSQRT(3.0D0)
        HCP_D_C(3, :) = HCP_D(4, :) * C_OVER_A
        !
        ! Normalize Cartesian vectors
        !
        DO I = 1, 18
            !
            HCP_N_C(:, I) = HCP_N_C(:, I) / NORM2(HCP_N_C(:, I))
            HCP_D_C(:, I) = HCP_D_C(:, I) / NORM2(HCP_D_C(:, I))
            !
            SCHMID(:, :, I) = MATMUL(RESHAPE(HCP_D_C(:, I), (/3, 1/)), &
                & RESHAPE(HCP_N_C(:, I), (/1, 3/)))
            !
            IF ((I .GE. 4) .AND. (I .LE. 6)) THEN
                !
                SCALE  = 1.0D0 / HRATIO_HCP_PRISM
                SCHMID(:, :, I) = SCALE * SCHMID(:, :, I)
                !
            ELSE IF ((I .GE. 7) .AND. (I .LE. 18)) THEN
                !
                SCALE  = 1.0D0 / HRATIO_HCP_PYR
                SCHMID(:, :, I) = SCALE * SCHMID(:, :, I)
                !
            END IF
            !
        END DO
        !
    CASE (CLASS_BCT)
        !
        ALLOCATE(SCHMID(3, 3, 32))
        SCHMID = 0.0D0
        !
        ! Convert Miller notation to Cartesian
        !
        BCT_N_C(1,:) = BCT_SN(1,:)
        BCT_N_C(2,:) = BCT_SN(2,:)
        BCT_N_C(3,:) = BCT_SN(3,:) / C_OVER_A
        !
        BCT_D_C(1,:) = BCT_SD(1,:)
        BCT_D_C(2,:) = BCT_SD(2,:)
        BCT_D_C(3,:) = BCT_SD(3,:) * C_OVER_A
        !
        ! Normalize Cartesian vectors, construct Schmid tensors
        !
        DO I = 1, 32
            !
            ! Normalize
            !
            BCT_N_C(:, I) = BCT_N_C(:, I) / NORM2(BCT_N_C(:, I))
            BCT_D_C(:, I) = BCT_D_C(:, I) / NORM2(BCT_D_C(:, I))
            !
            ! D*N^T
            !
            SCHMID(:, :, I) = MATMUL(RESHAPE(BCT_D_C(:, I), (/3, 1/)), &
                & RESHAPE(BCT_N_C(:, I), (/1, 3/)))
            !
            ! Scale by relative strengths
            !
            IF ((I .GE. 3) .AND. (I .LE. 4)) THEN
                !
                SCALE  = 1.0D0 / HRATIO_BCT_A
                SCHMID(:, :, I) = SCALE * SCHMID(:, :, I)
                !
            ELSE IF ((I .GE. 5) .AND. (I .LE. 6)) THEN
                !
                SCALE  = 1.0D0 / HRATIO_BCT_B
                SCHMID(:, :, I) = SCALE * SCHMID(:, :, I)
                !
            ELSE IF ((I .GE. 7) .AND. (I .LE. 10)) THEN
                !
                SCALE  = 1.0D0 / HRATIO_BCT_C
                SCHMID(:, :, I) = SCALE * SCHMID(:, :, I)
                !
            ELSE IF ((I .GE. 11) .AND. (I .LE. 12)) THEN
                !
                SCALE  = 1.0D0 / HRATIO_BCT_D
                SCHMID(:, :, I) = SCALE * SCHMID(:, :, I)
                !
            ELSE IF ((I .GE. 13) .AND. (I .LE. 16)) THEN
                !
                SCALE  = 1.0D0 / HRATIO_BCT_E
                SCHMID(:, :, I) = SCALE * SCHMID(:, :, I)
                !
            ELSE IF ((I .GE. 17) .AND. (I .LE. 18)) THEN
                !
                SCALE  = 1.0D0 / HRATIO_BCT_F
                SCHMID(:, :, I) = SCALE * SCHMID(:, :, I)
                !
            ELSE IF ((I .GE. 19) .AND. (I .LE. 20)) THEN
                !
                SCALE  = 1.0D0 / HRATIO_BCT_G
                SCHMID(:, :, I) = SCALE * SCHMID(:, :, I)
                !
            ELSE IF ((I .GE. 21) .AND. (I .LE. 24)) THEN
                !
                SCALE  = 1.0D0 / HRATIO_BCT_H
                SCHMID(:, :, I) = SCALE * SCHMID(:, :, I)
                !
            ELSE IF ((I .GE. 25) .AND. (I .LE. 32)) THEN
                !
                SCALE  = 1.0D0 / HRATIO_BCT_I
                SCHMID(:, :, I) = SCALE * SCHMID(:, :, I)
                !
            END IF
            !
        END DO
        !
    CASE DEFAULT
        !
    END SELECT
    !
    END SUBROUTINE SCHMIDTENSORS
    !
    !===========================================================================
    !
    SUBROUTINE CRYSTALTYPEGET(SELF, DEV, SKW, PPTRANS, VERTICES, NUMSLIP, &
        & NUMVERTICES)
    !
    ! Return deviatoric parts of Schmid tensors
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    ! SELF: The CrystalType object
    ! DEV: Deviatoric part of Schmid tensors
    ! SKW: Skew part of Schmid tensors
    ! PPTRANS: Matrices of Schmid diads
    ! VERTICES: Vertices in 5-vector form
    ! NUMSLIP: Number of slip systems
    ! NUMVERTICES: Number of vertices
    !
    TYPE(CRYSTALTYPETYPE) :: SELF
    REAL(RK), POINTER, OPTIONAL :: DEV(:, :)
    REAL(RK), POINTER, OPTIONAL :: SKW(:, :)
    REAL(RK), POINTER, OPTIONAL :: PPTRANS(:, :, :)
    REAL(RK), POINTER, OPTIONAL :: VERTICES(:, :)
    INTEGER, OPTIONAL :: NUMSLIP
    INTEGER, OPTIONAL :: NUMVERTICES
    !
    ! Locals:
    !
    INTEGER :: DSHAPE(2)
    INTEGER :: SSHAPE(2)
    INTEGER :: ARGPPTSHAPE(3)
    INTEGER :: MYPPTSHAPE(3)
    INTEGER :: I
    !
    REAL(RK), POINTER :: MYDEV(:, :)
    !
    !---------------------------------------------------------------------------
    !
    IF (PRESENT(DEV)) THEN
        !
        SSHAPE = (/5, SELF%NUMSLIP/)
        !SHAPE(SELF%schmid_sym)
        !
        IF (ASSOCIATED(DEV)) THEN
            !
            ! Could check dimensions here ...
            !
            DSHAPE = SHAPE(DEV)
            !
            IF ( (DSHAPE(1) /= SSHAPE(1)) .OR. (DSHAPE(2) /= SSHAPE(2))) THEN
                !
                DEALLOCATE(DEV)
                ALLOCATE(DEV(SSHAPE(1), SSHAPE(2)))
                !
            END IF
            !
        ELSE
            !
            ALLOCATE(DEV(SSHAPE(1), SSHAPE(2)))
            !
        END IF
        !
        CALL TENSOR3DDECOMPOSE(SELF%SCHMID_3X3, DEV = DEV, DECOMP = SELF%DECOMP)
    !
    END IF
    !
    IF (PRESENT(SKW)) THEN
        !
        SSHAPE = (/3, SELF%NUMSLIP/)
        !
        IF (ASSOCIATED(SKW)) THEN
            !
            ! Could check dimensions here ...
            !
            DSHAPE = SHAPE(SKW)
            !
            IF ( (DSHAPE(1) /= SSHAPE(1)) .OR. (DSHAPE(2) /= SSHAPE(2))) THEN
                !
                DEALLOCATE(SKW)
                ALLOCATE(SKW(SSHAPE(1), SSHAPE(2)))
                !
            END IF
            !
        ELSE
            !
            ALLOCATE(SKW(SSHAPE(1), SSHAPE(2)))
            !
        END IF
        !
        CALL TENSOR3DDECOMPOSE(SELF%SCHMID_3X3, SKW = SKW, DECOMP = SELF%DECOMP)
        !
    END IF
    !
    IF (PRESENT(PPTRANS)) THEN
        !
        MYPPTSHAPE = (/5, 5, SELF%NUMSLIP/)
        !
        IF (ASSOCIATED(PPTRANS)) THEN
            !
            !  Could check dimensions here ...
            !
            ARGPPTSHAPE = SHAPE(PPTRANS)
            !
            IF ((ARGPPTSHAPE(1) /= MYPPTSHAPE(1)) .OR. &
                & (ARGPPTSHAPE(2) /= MYPPTSHAPE(2)) .OR. &
                & (ARGPPTSHAPE(3) /= MYPPTSHAPE(3))  ) THEN
                !
                DEALLOCATE(PPTRANS)
                ALLOCATE(PPTRANS(MYPPTSHAPE(1), MYPPTSHAPE(2), MYPPTSHAPE(3)))
                !
            END IF
            !
        ELSE
            !
            ALLOCATE(PPTRANS(MYPPTSHAPE(1), MYPPTSHAPE(2), MYPPTSHAPE(3)))
            !
        END IF
        !
        ALLOCATE(MYDEV(5, SELF%NUMSLIP))
        !
        CALL TENSOR3DDECOMPOSE(SELF%SCHMID_3X3, DEV = MYDEV, &
            & DECOMP = SELF%DECOMP)
        !
        DO I = 1, SELF%NUMSLIP
            !
            PPTRANS(:, :, I) = MATMUL(&
                & RESHAPE(MYDEV(:, I), SHAPE = (/5, 1/)), &
                & RESHAPE(MYDEV(:, I), SHAPE=  (/1, 5/)) )
            !
        END DO
        !
        DEALLOCATE(MYDEV)
        !
    END IF
    !
    IF (PRESENT(VERTICES)) THEN
        !
        SSHAPE = (/5, SELF%NUMVERTICES/)
        !
        IF (ASSOCIATED(VERTICES)) THEN
            !
            !  Could check dimensions here ...
            !
            DSHAPE = SHAPE(VERTICES)
            !
            IF ((DSHAPE(1) /= SSHAPE(1)) .OR. (DSHAPE(2) /= SSHAPE(2))) THEN
                !
                DEALLOCATE(VERTICES)
                ALLOCATE(VERTICES(SSHAPE(1), SSHAPE(2)))
                !
            END IF
            !
        ELSE
            !
            ALLOCATE(VERTICES(SSHAPE(1), SSHAPE(2)))
            !
        END IF
        !
        VERTICES = SELF%VERTICES
        !
    END IF
    !
    IF (PRESENT(NUMSLIP)) THEN
        !
        NUMSLIP = SELF%NUMSLIP
        !
    END IF
    !
    IF (PRESENT(NUMVERTICES)) THEN
        !
        NUMVERTICES = SELF%NUMVERTICES
        !
    END IF
    !
    ! Need call to get vertices as 3X3.
    !
    END SUBROUTINE CRYSTALTYPEGET
    !
    !===========================================================================
    !
    SUBROUTINE GET_VERTICES(CTYPE, VERTICES)
    !
    ! Return 3X3 vertex tensors in place of legacy VERTICES files.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! CTYPE: Crystal CLASS type (FCC, BCC, HCP)
    ! VERTICES: 3D array containing VERTICES for a specific crystal type.
    !
    INTEGER, INTENT(IN) :: CTYPE
    REAL(RK), POINTER, INTENT(OUT) :: VERTICES(:, :, :)
    !
    ! Locals:
    ! VERTICES_FCC_DAT: 1D array containing vertex file data for FCC/BCC.
    ! VERT_FCC: 3D array reshape of VERTICES_FCC_DAT.
    ! VERTICES_HCP_DAT: 1D array containing vertex file data for HCP.
    ! VERT_HCP: 3D array reshape of VERTICES_HCP_DAT.
    !
    REAL(RK), PARAMETER, DIMENSION(252) :: VERTICES_FCC_DAT = (/&
    &  1.6329855197502110,      0.0000000000000000,      0.0000000000000000,&
    &  0.0000000000000000,     -0.8164889388224847,      0.0000000000000000,&
    &  0.0000000000000000,      0.0000000000000000,     -0.8164965809277259,&
    & -0.8164889388224847,      0.0000000000000000,      0.0000000000000000,&
    &  0.0000000000000000,      1.6329855197502110,      0.0000000000000000,&
    &  0.0000000000000000,      0.0000000000000000,     -0.8164965809277259,&
    & -0.8164965809277259,      0.0000000000000000,      0.0000000000000000,&
    &  0.0000000000000000,     -0.8164965809277259,      0.0000000000000000,&
    &  0.0000000000000000,      0.0000000000000000,      1.6329931618554521,&
    &  0.0000000000000000,      1.2247372292863481,      1.2247372292863481,&
    &  1.2247372292863481,      0.0000000000000000,      1.2247372292863481,&
    &  1.2247372292863481,      1.2247372292863481,      0.0000000000000000,&
    &  0.0000000000000000,      1.2247372292863481,      1.2247372292863481,&
    &  1.2247372292863481,      0.0000000000000000,     -1.2247372292863481,&
    &  1.2247372292863481,     -1.2247372292863481,      0.0000000000000000,&
    &  0.0000000000000000,      1.2247372292863481,     -1.2247372292863481,&
    &  1.2247372292863481,      0.0000000000000000,      1.2247372292863481,&
    & -1.2247372292863481,      1.2247372292863481,      0.0000000000000000,&
    &  0.0000000000000000,     -1.2247372292863481,      1.2247372292863481,&
    & -1.2247372292863481,      0.0000000000000000,      1.2247372292863481,&
    &  1.2247372292863481,      1.2247372292863481,      0.0000000000000000,&
    &  0.0000000000000000,      0.0000000000000000,      0.0000000000000000,&
    &  0.0000000000000000,      0.0000000000000000,      2.4494886007083192,&
    &  0.0000000000000000,      2.4494886007083192,      0.0000000000000000,&
    &  0.0000000000000000,      0.0000000000000000,      2.4494886007083192,&
    &  0.0000000000000000,      0.0000000000000000,      0.0000000000000000,&
    &  2.4494886007083192,      0.0000000000000000,      0.0000000000000000,&
    &  0.0000000000000000,      2.4494886007083192,      0.0000000000000000,&
    &  2.4494886007083192,      0.0000000000000000,      0.0000000000000000,&
    &  0.0000000000000000,      0.0000000000000000,      0.0000000000000000,&
    &  0.8164927598751053,      1.2247372292863481,      1.2247372292863481,&
    &  1.2247372292863481,     -0.4082444694112423,      0.0000000000000000,&
    &  1.2247372292863481,      0.0000000000000000,     -0.4082482904638630,&
    &  0.8164927598751053,     -1.2247372292863481,     -1.2247372292863481,&
    & -1.2247372292863481,     -0.4082444694112423,      0.0000000000000000,&
    & -1.2247372292863481,      0.0000000000000000,     -0.4082482904638630,&
    &  0.8164927598751053,     -1.2247372292863481,      1.2247372292863481,&
    & -1.2247372292863481,     -0.4082444694112423,      0.0000000000000000,&
    &  1.2247372292863481,      0.0000000000000000,     -0.4082482904638630,&
    &  0.8164927598751053,      1.2247372292863481,     -1.2247372292863481,&
    &  1.2247372292863481,     -0.4082444694112423,      0.0000000000000000,&
    & -1.2247372292863481,      0.0000000000000000,     -0.4082482904638630,&
    & -0.4082444694112423,      1.2247372292863481,      0.0000000000000000,&
    &  1.2247372292863481,      0.8164927598751053,      1.2247372292863481,&
    &  0.0000000000000000,      1.2247372292863481,     -0.4082482904638630,&
    & -0.4082444694112423,     -1.2247372292863481,      0.0000000000000000,&
    & -1.2247372292863481,      0.8164927598751053,     -1.2247372292863481,&
    &  0.0000000000000000,     -1.2247372292863481,     -0.4082482904638630,&
    & -0.4082444694112423,      1.2247372292863481,      0.0000000000000000,&
    &  1.2247372292863481,      0.8164927598751053,     -1.2247372292863481,&
    &  0.0000000000000000,     -1.2247372292863481,     -0.4082482904638630,&
    & -0.4082444694112423,     -1.2247372292863481,      0.0000000000000000,&
    & -1.2247372292863481,      0.8164927598751053,      1.2247372292863481,&
    &  0.0000000000000000,      1.2247372292863481,     -0.4082482904638630,&
    & -0.4082482904638630,      0.0000000000000000,      1.2247372292863481,&
    &  0.0000000000000000,     -0.4082482904638630,      1.2247372292863481,&
    &  1.2247372292863481,      1.2247372292863481,      0.8164965809277259,&
    & -0.4082482904638630,      0.0000000000000000,     -1.2247372292863481,&
    &  0.0000000000000000,     -0.4082482904638630,     -1.2247372292863481,&
    & -1.2247372292863481,     -1.2247372292863481,      0.8164965809277259,&
    & -0.4082482904638630,      0.0000000000000000,     -1.2247372292863481,&
    &  0.0000000000000000,     -0.4082482904638630,      1.2247372292863481,&
    & -1.2247372292863481,      1.2247372292863481,      0.8164965809277259,&
    & -0.4082482904638630,      0.0000000000000000,      1.2247372292863481,&
    &  0.0000000000000000,     -0.4082482904638630,     -1.2247372292863481,&
    &  1.2247372292863481,     -1.2247372292863481,      0.8164965809277259,&
    &  0.0000038210526208,      0.0000000000000000,      0.0000000000000000,&
    &  0.0000000000000000,      1.2247410503389680,      1.2247372292863481,&
    &  0.0000000000000000,      1.2247372292863481,     -1.2247448713915889,&
    &  0.0000038210526208,      0.0000000000000000,      0.0000000000000000,&
    &  0.0000000000000000,      1.2247410503389680,     -1.2247372292863481,&
    &  0.0000000000000000,     -1.2247372292863481,     -1.2247448713915889,&
    & -1.2247410503389680,      0.0000000000000000,      1.2247372292863481,&
    &  0.0000000000000000,     -0.0000038210526208,      0.0000000000000000,&
    &  1.2247372292863481,      0.0000000000000000,      1.2247448713915889,&
    & -1.2247410503389680,      0.0000000000000000,     -1.2247372292863481,&
    &  0.0000000000000000,     -0.0000038210526208,      0.0000000000000000,&
    & -1.2247372292863481,      0.0000000000000000,      1.2247448713915889,&
    &  1.2247372292863470,      1.2247372292863481,      0.0000000000000000,&
    &  1.2247372292863481,     -1.2247372292863470,      0.0000000000000000,&
    &  0.0000000000000000,      0.0000000000000000,      0.0000000000000000,&
    &  1.2247372292863470,     -1.2247372292863481,      0.0000000000000000,&
    & -1.2247372292863481,     -1.2247372292863470,      0.0000000000000000,&
    &  0.0000000000000000,      0.0000000000000000,      0.0000000000000000&
    &/)
    !
    REAL(RK), PARAMETER, DIMENSION(3, 3, 28) :: &
        &   VERT_FCC = RESHAPE(SOURCE=VERTICES_FCC_DAT, SHAPE=(/3, 3, 28/))
    !
    REAL(RK), PARAMETER, DIMENSION(2160) :: VERTICES_HCP_DAT = (/&
    &  3.2380406000000002,      0.0000000000000000,     -0.7116558600000000, &
    &  0.0000000000000000,      0.9286395100000000,     -0.7438258300000000, &
    & -0.7116558600000000,     -0.7438258300000000,     -4.1666800999999998, &
    &  3.2380406000000002,      0.0000000000000000,      0.7116558600000000, &
    &  0.0000000000000000,      0.9286395100000000,      0.7438258300000000, &
    &  0.7116558600000000,      0.7438258300000000,     -4.1666800999999998, &
    &  3.0950593000000000,      0.0000000000000000,     -1.0000000000000000, &
    &  0.0000000000000000,      0.7856582600000001,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,     -3.8807176000000001, &
    &  3.0950593000000000,      0.0000000000000000,      1.0000000000000000, &
    &  0.0000000000000000,      0.7856582600000001,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,     -3.8807176000000001, &
    &  2.8851517000000002,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,      0.5757506200000000,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,     -3.4609022999999999, &
    &  2.8851517000000002,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,      0.5757506200000000,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,     -3.4609022999999999, &
    &  2.7421704000000000,      0.0000000000000000,      1.0000000000000000, &
    &  0.0000000000000000,      0.4327693700000000,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,     -3.1749398000000002, &
    &  2.7421704000000000,      0.0000000000000000,     -1.0000000000000000, &
    &  0.0000000000000000,      0.4327693700000000,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,     -3.1749398000000002, &
    & -0.4327693800000000,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,     -2.7421704000000000,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,      3.1749398000000002, &
    & -0.4327693800000000,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,     -2.7421704000000000,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,      3.1749398000000002, &
    & -0.5757506200000000,      0.0000000000000000,     -1.0000000000000000, &
    &  0.0000000000000000,     -2.8851517000000002,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,      3.4609022999999999, &
    & -0.5757506200000000,      0.0000000000000000,      1.0000000000000000, &
    &  0.0000000000000000,     -2.8851517000000002,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,      3.4609022999999999, &
    & -0.5757506200000000,      0.0000000000000000,     -1.0000000000000000, &
    &  0.0000000000000000,     -2.8851517000000002,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,      3.4609022999999999, &
    & -0.7856582700000000,      0.0000000000000000,      1.0000000000000000, &
    &  0.0000000000000000,     -3.0950593000000000,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,      3.8807176000000001, &
    & -0.8236856900000000,      0.0000000000000000,     -1.0000000000000000, &
    &  0.0000000000000000,     -3.1330868000000001,      0.0000000000000000, &
    & -1.0000000000000000,      0.0000000000000000,      3.9567724000000002, &
    & -0.9286395100000000,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,     -3.2380406000000002,      0.0000000000000000, &
    &  0.0000000000000000,      0.0000000000000000,      4.1666800999999998, &
    & -0.9286395100000000,      0.0000000000000000,      1.0000000000000000, &
    &  0.0000000000000000,     -3.2380406000000002,      0.2443991400000000, &
    &  1.0000000000000000,      0.2443991400000000,      4.1666800999999998, &
    & -1.1502695999999999,      0.3789152600000000,      1.0000000000000000, &
    &  0.3789152600000000,     -3.0221371000000001,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,      4.1724066999999998, &
    & -0.9978762799999999,      0.4386858500000000,      1.0000000000000000, &
    &  0.4386858500000000,     -2.8007266000000000,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,      3.7986027999999998, &
    & -0.9978762799999999,      0.4386858500000000,     -1.0000000000000000, &
    &  0.4386858500000000,     -2.8007266000000000,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,      3.7986027999999998, &
    & -1.0065039000000000,      0.4685726900000000,      0.8592050800000000, &
    &  0.4685726900000000,     -2.7748438000000002,     -0.6586382500000000, &
    &  0.8592050800000000,     -0.6586382500000000,      3.7813477000000000, &
    & -1.2104325000000000,      0.5023988600000000,     -1.0000000000000000, &
    &  0.5023988600000000,     -2.9397133000000002,      0.2251480400000000, &
    & -1.0000000000000000,      0.2251480400000000,      4.1501459000000001, &
    & -0.9357977900000000,      0.5227624700000000,      0.0000000000000000, &
    &  0.5227624700000000,     -2.6415647999999998,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,      3.5773625999999998, &
    & -1.1005197000000000,     -0.5453560500000000,     -1.0000000000000000, &
    & -0.5453560500000000,     -2.7801979000000001,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,      3.8807176000000001, &
    & -1.2435010000000000,     -0.5453560500000000,     -1.0000000000000000, &
    & -0.5453560500000000,     -2.9231791000000000,     -0.2443991400000000, &
    & -1.0000000000000000,     -0.2443991400000000,      4.1666800999999998, &
    & -1.3537442000000000,      0.5453560500000000,      1.0000000000000000, &
    &  0.5453560500000000,     -3.0334222999999998,      0.0123171000000000, &
    &  1.0000000000000000,      0.0123171000000000,      4.3871665000000002, &
    & -1.4534085999999999,      0.5453560500000000,      1.0000000000000000, &
    &  0.5453560500000000,     -3.1330868000000001,      0.2443991400000000, &
    &  1.0000000000000000,      0.2443991400000000,      4.5864953999999996, &
    & -1.4547916999999999,      0.5471527700000000,      1.0000000000000000, &
    &  0.5471527700000000,     -3.1323951999999999,      0.2452043300000000, &
    &  1.0000000000000000,      0.2452043300000000,      4.5871868999999998, &
    & -1.4618940000000000,      0.5453560500000000,      1.0000000000000000, &
    &  0.5453560500000000,     -3.1288440999999998,      0.2443991400000000, &
    &  1.0000000000000000,      0.2443991400000000,      4.5907380000000000, &
    & -1.0508058000000000,      0.6220392500000000,     -0.1362340200000000, &
    &  0.6220392500000000,     -2.6419378000000000,      1.0760457999999999, &
    & -0.1362340200000000,      1.0760457999999999,      3.6927436999999999, &
    & -1.5016307000000000,      0.6380901900000000,      1.0000000000000000, &
    &  0.6380901900000000,     -3.0742286999999999,      0.2320154900000000, &
    &  1.0000000000000000,      0.2320154900000000,      4.5758593999999997, &
    &  3.1067741999999998,      0.6820803100000000,     -0.7116558600000000, &
    &  0.6820803100000000,      1.5849716000000000,     -0.7438258300000000, &
    & -0.7116558600000000,     -0.7438258300000000,     -4.6917457999999996, &
    &  3.1067741999999998,      0.6820803100000000,      0.7116558600000000, &
    &  0.6820803100000000,      1.5849716000000000,      0.7438258300000000, &
    &  0.7116558600000000,      0.7438258300000000,     -4.6917457999999996, &
    &  2.7284468999999998,      0.8142619200000000,      0.0000000000000000, &
    &  0.8142619200000000,      1.3592744999999999,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,     -4.0877214000000004, &
    &  2.7284468999999998,      0.8142619200000000,      0.0000000000000000, &
    &  0.8142619200000000,      1.3592744999999999,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,     -4.0877214000000004, &
    & -1.6591488000000001,      0.5453560500000000,      1.0000000000000000, &
    &  0.5453560500000000,     -3.0133475999999999,      0.2182111200000000, &
    &  1.0000000000000000,      0.2182111200000000,      4.6724962999999997, &
    & -1.6497605000000000,      0.5919432300000000,      1.0000000000000000, &
    &  0.5919432300000000,     -3.0030958000000001,      0.2158866000000000, &
    &  1.0000000000000000,      0.2158866000000000,      4.6528562999999998, &
    & -1.5956191000000000,      0.8438615600000000,      0.1655457000000000, &
    &  0.8438615600000000,     -2.9306128000000000,      0.3781731900000000, &
    &  0.1655457000000000,      0.3781731900000000,      4.5262319000000000, &
    & -1.6438908999999999,      0.5453560500000000,      1.0000000000000000, &
    &  0.5453560500000000,     -2.9669104000000002,      0.1342776300000000, &
    &  1.0000000000000000,      0.1342776300000000,      4.6108013000000003, &
    & -1.7194469999999999,      0.5077001799999999,      0.9309517700000000, &
    &  0.5077001799999999,     -2.9891972999999998,      0.2275238100000000, &
    &  0.9309517700000000,      0.2275238100000000,      4.7086442000000002, &
    & -1.7144132000000001,      0.4578808100000000,      1.0000000000000000, &
    &  0.4578808100000000,     -2.9773325000000002,      0.2443991400000000, &
    &  1.0000000000000000,      0.2443991400000000,      4.6917457999999996, &
    &  2.6844529000000001,     -0.9588419900000000,     -0.7116558600000000, &
    & -0.9588419900000000,      1.4822272000000001,     -0.7438258300000000, &
    & -0.7116558600000000,     -0.7438258300000000,     -4.1666800999999998, &
    &  2.9370479000000000,      0.3881058300000000,      0.7116558600000000, &
    &  0.3881058300000000,      1.7546978000000000,      0.7438258300000000, &
    &  0.7116558600000000,      0.7438258300000000,     -4.6917457999999996, &
    &  2.6745150999999998,     -0.9760548000000000,      0.7116558600000000, &
    & -0.9760548000000000,      1.4921650000000000,      0.7438258300000000, &
    &  0.7116558600000000,      0.7438258300000000,     -4.1666800999999998, &
    &  3.0455904999999999,      1.0000000000000000,      0.0000000000000001, &
    &  1.0000000000000000,      1.8908900000000000,      0.0000000000000002, &
    &  0.0000000000000001,      0.0000000000000002,     -4.9364803999999998, &
    &  2.9026092000000001,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,      1.7479087000000000,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,     -4.6505178999999996, &
    &  2.9026092000000001,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,      1.7479087000000000,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,     -4.6505178999999996, &
    &  2.9026092000000001,      1.0000000000000000,     -0.9999999899999999, &
    &  1.0000000000000000,      1.7479087000000000,     -0.5773502700000001, &
    & -0.9999999899999999,     -0.5773502700000001,     -4.6505178999999996, &
    &  2.9026092000000001,      1.0000000000000000,      0.9999999899999999, &
    &  1.0000000000000000,      1.7479087000000000,      0.5773502700000001, &
    &  0.9999999899999999,      0.5773502700000001,     -4.6505178999999996, &
    &  2.9026092000000001,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,      1.7479087000000000,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,     -4.6505178999999996, &
    &  2.9026092000000001,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,      1.7479087000000000,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,     -4.6505178999999996, &
    &  2.6606903000000002,     -1.0000000000000000,     -0.6961355000000000, &
    & -1.0000000000000000,      1.5059898000000000,     -0.7276038900000000, &
    & -0.6961355000000000,     -0.7276038900000000,     -4.1666800999999998, &
    &  2.6606903000000002,     -1.0000000000000000,      0.6950136200000000, &
    & -1.0000000000000000,      1.5059898000000000,      0.7264312900000000, &
    &  0.6950136200000000,      0.7264312900000000,     -4.1666800999999998, &
    &  2.6535365000000000,     -1.0000000000000000,      0.7260825500000000, &
    & -1.0000000000000000,      1.4988360000000001,      0.7354965800000000, &
    &  0.7260825500000000,      0.7354965800000000,     -4.1523725999999996, &
    &  2.6519420000000000,     -1.0000000000000000,     -0.7292981400000000, &
    & -1.0000000000000000,      1.4972414999999999,     -0.7336400600000000, &
    & -0.7292981400000000,     -0.7336400600000000,     -4.1491835000000004, &
    &  2.6477395000000001,     -1.0000000000000000,      0.6855385400000000, &
    & -1.0000000000000000,      1.4930390000000000,      0.7589046800000000, &
    &  0.6855385400000000,      0.7589046800000000,     -4.1407784999999997, &
    &  2.6303709000000000,     -1.0000000000000000,     -0.6505120300000000, &
    & -1.0000000000000000,      1.4756704000000000,     -0.7791272400000000, &
    & -0.6505120300000000,     -0.7791272400000000,     -4.1060413000000002, &
    &  2.6212110000000002,      1.0000000000000000,      0.4999999800000000, &
    &  1.0000000000000000,      1.4665104000000000,     -0.8660254100000000, &
    &  0.4999999800000000,     -0.8660254100000000,     -4.0877214000000004, &
    &  2.6212110000000002,      1.0000000000000000,     -0.4999999800000000, &
    &  1.0000000000000000,      1.4665104000000000,      0.8660254100000000, &
    & -0.4999999800000000,      0.8660254100000000,     -4.0877214000000004, &
    &  2.5497204000000000,      1.0000000000000000,      0.0000000000000000, &
    &  1.0000000000000000,      1.3950198000000000,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,     -3.9447402000000000, &
    &  2.5497204000000000,      1.0000000000000000,      0.0000000000000000, &
    &  1.0000000000000000,      1.3950198000000000,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,     -3.9447402000000000, &
    &  2.5497203000000002,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,      1.3950198000000000,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,     -3.9447402000000000, &
    &  2.5497203000000002,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,      1.3950198000000000,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,     -3.9447402000000000, &
    &  2.5177090999999998,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,      1.3630085000000001,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,     -3.8807176000000001, &
    &  2.5177090999999998,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,      1.3630085000000001,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,     -3.8807176000000001, &
    &  2.3078013999999998,     -1.0000000000000000,      0.0000000000000000, &
    & -1.0000000000000000,      1.1531009000000001,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,     -3.4609022999999999, &
    &  2.3078013999999998,     -1.0000000000000000,      0.0000000000000000, &
    & -1.0000000000000000,      1.1531009000000001,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,     -3.4609022999999999, &
    &  2.1648201999999999,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,      1.0101195999999999,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,     -3.1749398000000002, &
    &  2.1648201999999999,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,      1.0101195999999999,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,     -3.1749398000000002, &
    & -1.0101195999999999,     -1.0000000000000000,      0.0000000000000000, &
    & -1.0000000000000000,     -2.1648201999999999,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,      3.1749398000000002, &
    & -1.0101195999999999,      1.0000000000000000,      0.0000000000000000, &
    &  1.0000000000000000,     -2.1648201999999999,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,      3.1749398000000002, &
    & -1.0101195999999999,     -1.0000000000000000,      0.0000000000000000, &
    & -1.0000000000000000,     -2.1648201999999999,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,      3.1749398000000002, &
    & -1.1326072000000000,      1.0000000000000000,      0.0000000000000000, &
    &  1.0000000000000000,     -2.2873077999999998,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,      3.4199150000000000, &
    & -1.1531009000000001,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,     -2.3078013999999998,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,      3.4609022999999999, &
    & -1.1531009000000001,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,     -2.3078013999999998,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,      3.4609022999999999, &
    & -1.1531009000000001,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,     -2.3078013999999998,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,      3.4609022999999999, &
    & -1.1599136999999999,      1.0000000000000000,      0.4244793700000000, &
    &  1.0000000000000000,     -2.3146141999999998,     -0.9096272500000000, &
    &  0.4244793700000000,     -0.9096272500000000,      3.4745279000000000, &
    & -1.1599136999999999,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,     -2.3146141999999998,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,      3.4745279000000000, &
    & -1.1599136999999999,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,     -2.3146141999999998,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,      3.4745279000000000, &
    & -1.1599136999999999,      1.0000000000000000,     -0.0502257340000000, &
    &  1.0000000000000000,     -2.3146141999999998,      1.1257026999999999, &
    & -0.0502257340000000,      1.1257026999999999,      3.4745279000000000, &
    & -1.3630085000000001,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,     -2.5177090999999998,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,      3.8807176000000001, &
    & -1.3630085000000001,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,     -2.5177090999999998,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,      3.8807176000000001, &
    & -1.4064007000000001,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,     -2.5611012000000000,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,      3.9675018999999998, &
    & -1.4348236999999999,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,     -2.5895242000000001,      0.0786792310000000, &
    & -1.0000000000000000,      0.0786792310000000,      4.0243479000000004, &
    & -1.5059898000000000,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,     -2.6606903000000002,     -0.2443991400000000, &
    & -1.0000000000000000,     -0.2443991400000000,      4.1666800999999998, &
    & -1.5059898000000000,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,     -2.6606903000000002,      0.2443991400000000, &
    &  1.0000000000000000,      0.2443991400000000,      4.1666800999999998, &
    & -1.6257071999999999,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,     -2.7804077000000000,      0.2675726300000000, &
    &  1.0000000000000000,      0.2675726300000000,      4.4061149000000004, &
    & -1.6768235000000000,      1.0000000000000000,      0.2199620200000000, &
    &  1.0000000000000000,     -2.8315239999999999,      0.3440507400000000, &
    &  0.2199620200000000,      0.3440507400000000,      4.5083473999999999, &
    & -1.6844881000000000,      1.0000000000000000,      0.7530396000000000, &
    &  1.0000000000000000,     -2.8391886999999998,      0.2316152600000000, &
    &  0.7530396000000000,      0.2316152600000000,      4.5236767999999996, &
    & -1.6286773000000001,      0.4639611400000000,     -1.0000000000000000, &
    &  0.4639611400000000,     -2.7305909000000002,      0.2079223400000000, &
    & -1.0000000000000000,      0.2079223400000000,      4.3592683000000001, &
    & -1.6845969000000001,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,     -2.7414003000000000,      0.2634078000000000, &
    &  1.0000000000000000,      0.2634078000000000,      4.4259972000000003, &
    & -1.7677229000000001,      1.0000000000000000,      0.7377517000000000, &
    &  1.0000000000000000,     -2.7901045999999998,      0.2237601200000000, &
    &  0.7377517000000000,      0.2237601200000000,      4.5578275000000001, &
    &  2.8157193000000000,      0.1779585300000000,     -0.7116558600000000, &
    &  0.1779585300000000,      1.8760264000000000,     -0.7438258300000000, &
    & -0.7116558600000000,     -0.7438258300000000,     -4.6917457999999996, &
    & -1.6122532000000001,      0.3760898600000000,      1.0000000000000000, &
    &  0.3760898600000000,     -2.4754681999999999,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,      4.0877214000000004, &
    & -1.6122532000000001,      0.3760898600000000,     -1.0000000000000000, &
    &  0.3760898600000000,     -2.4754681999999999,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,      4.0877214000000004, &
    & -1.6764927999999999,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,     -2.4584896999999999,      0.0628446440000000, &
    & -1.0000000000000000,      0.0628446440000000,      4.1349824999999996, &
    &  2.4067390999999998,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,      1.6809822999999999,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,     -4.0877214000000004, &
    &  2.4067390999999998,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,      1.6809822999999999,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,     -4.0877214000000004, &
    & -1.9962358000000000,      0.7781293400000000,      0.0761197940000000, &
    &  0.7781293400000000,     -2.7162259999999998,      0.3487155600000000, &
    &  0.0761197940000000,      0.3487155600000000,      4.7124617999999998, &
    & -1.6867365000000001,      0.4190927700000000,      0.7684755100000000, &
    &  0.4190927700000000,     -2.4009849000000001,     -0.7110209900000000, &
    &  0.7684755100000000,     -0.7110209900000000,      4.0877214000000004, &
    &  2.6401203999999998,      0.5453560500000000,      1.0000000000000000, &
    &  0.5453560500000000,      2.0103974999999998,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,     -4.6505178999999996, &
    &  2.3427164999999999,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,      1.7129935999999999,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,     -4.0557100999999998, &
    & -2.0437029000000000,      1.0000000000000000,      0.1638357500000000, &
    &  1.0000000000000000,     -2.6154294000000000,      0.3070738300000000, &
    &  0.1638357500000000,      0.3070738300000000,      4.6591322999999996, &
    &  2.1055524000000001,     -1.0000000000000000,      0.0000000000000000, &
    & -1.0000000000000000,      1.5575988999999999,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,     -3.6631513000000000, &
    & -2.1893172999999999,      0.6488078400000000,      0.2781854400000000, &
    &  0.6488078400000000,     -2.6514966000000002,      0.2907606400000000, &
    &  0.2781854400000000,      0.2907606400000000,      4.8408138999999997, &
    & -2.1301480000000002,      1.0000000000000000,      0.2629858300000000, &
    &  1.0000000000000000,     -2.5670744999999999,      0.2748739500000000, &
    &  0.2629858300000000,      0.2748739500000000,      4.6972224000000002, &
    &  2.2583326000000001,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,      1.8293889000000001,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,     -4.0877214000000004, &
    &  2.0658824999999998,     -1.0000000000000000,      0.0000000000000000, &
    & -1.0000000000000000,      1.6369388000000000,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,     -3.7028211999999998, &
    &  2.5561210999999999,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,      2.1356245999999999,     -0.2443991400000000, &
    & -1.0000000000000000,     -0.2443991400000000,     -4.6917457999999996, &
    &  2.5561210999999999,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,      2.1356245999999999,      0.2443991400000000, &
    &  1.0000000000000000,      0.2443991400000000,     -4.6917457999999996, &
    & -1.8400506999999999,      0.5076087800000000,     -0.2919118400000000, &
    &  0.5076087800000000,     -2.2476707000000000,      0.9861651500000000, &
    & -0.2919118400000000,      0.9861651500000000,      4.0877214000000004, &
    & -1.6637124999999999,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,     -2.0123350000000002,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,      3.6760475000000001, &
    & -1.6637124999999999,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,     -2.0123350000000002,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,      3.6760475000000001, &
    &  2.2116950000000002,     -0.0807787190000000,      0.0000000000000000, &
    & -0.0807787190000000,      1.8760265000000000,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,     -4.0877214000000004, &
    &  2.2260198999999998,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,      1.9463866999999999,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,     -4.1724066999999998, &
    & -2.2936323000000001,     -0.5453560500000000,     -1.0000000000000000, &
    & -0.5453560500000000,     -2.3981135000000000,     -0.2443991400000000, &
    & -1.0000000000000000,     -0.2443991400000000,      4.6917457999999996, &
    &  1.9180891000000000,     -1.0000000000000000,      0.1697576400000000, &
    & -1.0000000000000000,      1.8597090000000001,     -1.0566909000000000, &
    &  0.1697576400000000,     -1.0566909000000000,     -3.7777981000000000, &
    &  2.3455477999999998,      0.0351412330000000,     -1.0000000000000000, &
    &  0.0351412330000000,      2.3049702000000001,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,     -4.6505178999999996, &
    &  1.8818313000000000,     -1.0000000000000000,     -0.2880969500000000, &
    & -1.0000000000000000,      1.8814636000000000,      0.9883676800000000, &
    & -0.2880969500000000,      0.9883676800000000,     -3.7632949999999998, &
    &  2.4682401999999999,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,      2.4682401999999999,      0.0000000000000000, &
    &  0.0000000000000000,      0.0000000000000000,     -4.9364803999999998, &
    & -2.4682401999999999,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,     -2.4682401999999999,      0.0000000000000000, &
    &  0.0000000000000000,      0.0000000000000000,      4.9364803999999998, &
    & -1.8818313000000000,      1.0000000000000000,      0.2880969500000000, &
    &  1.0000000000000000,     -1.8814636000000000,     -0.9883676800000000, &
    &  0.2880969500000000,     -0.9883676800000000,      3.7632949999999998, &
    & -2.3455477999999998,     -0.0351412330000000,      1.0000000000000000, &
    & -0.0351412330000000,     -2.3049702000000001,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,      4.6505178999999996, &
    & -1.9180891000000000,      1.0000000000000000,     -0.1697576400000000, &
    &  1.0000000000000000,     -1.8597090000000001,      1.0566909000000000, &
    & -0.1697576400000000,      1.0566909000000000,      3.7777981000000000, &
    &  2.2936323000000001,      0.5453560500000000,      1.0000000000000000, &
    &  0.5453560500000000,      2.3981135000000000,      0.2443991400000000, &
    &  1.0000000000000000,      0.2443991400000000,     -4.6917457999999996, &
    & -2.2260198999999998,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,     -1.9463866999999999,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,      4.1724066999999998, &
    & -2.2116950000000002,      0.0807787190000000,      0.0000000000000000, &
    &  0.0807787190000000,     -1.8760265000000000,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,      4.0877214000000004, &
    &  1.6637124999999999,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,      2.0123350000000002,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,     -3.6760475000000001, &
    &  1.6637124999999999,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,      2.0123350000000002,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,     -3.6760475000000001, &
    &  1.8400506999999999,     -0.5076087800000000,      0.2919118400000000, &
    & -0.5076087800000000,      2.2476707000000000,     -0.9861651500000000, &
    &  0.2919118400000000,     -0.9861651500000000,     -4.0877214000000004, &
    & -2.5561210999999999,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,     -2.1356245999999999,     -0.2443991400000000, &
    & -1.0000000000000000,     -0.2443991400000000,      4.6917457999999996, &
    & -2.5561210999999999,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,     -2.1356245999999999,      0.2443991400000000, &
    &  1.0000000000000000,      0.2443991400000000,      4.6917457999999996, &
    & -2.0658824999999998,      1.0000000000000000,      0.0000000000000000, &
    &  1.0000000000000000,     -1.6369388000000000,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,      3.7028211999999998, &
    & -2.2583326000000001,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,     -1.8293889000000001,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,      4.0877214000000004, &
    &  2.1301480000000002,     -1.0000000000000000,     -0.2629858300000000, &
    & -1.0000000000000000,      2.5670744999999999,     -0.2748739500000000, &
    & -0.2629858300000000,     -0.2748739500000000,     -4.6972224000000002, &
    &  2.1893172999999999,     -0.6488078400000000,     -0.2781854400000000, &
    & -0.6488078400000000,      2.6514966000000002,     -0.2907606400000000, &
    & -0.2781854400000000,     -0.2907606400000000,     -4.8408138999999997, &
    & -2.1055524000000001,      1.0000000000000000,      0.0000000000000000, &
    &  1.0000000000000000,     -1.5575988999999999,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,      3.6631513000000000, &
    &  2.0437029000000000,     -1.0000000000000000,     -0.1638357500000000, &
    & -1.0000000000000000,      2.6154294000000000,     -0.3070738300000000, &
    & -0.1638357500000000,     -0.3070738300000000,     -4.6591322999999996, &
    & -2.3427164999999999,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,     -1.7129935999999999,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,      4.0557100999999998, &
    & -2.6401203999999998,     -0.5453560500000000,     -1.0000000000000000, &
    & -0.5453560500000000,     -2.0103974999999998,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,      4.6505178999999996, &
    &  1.6867365000000001,     -0.4190927700000000,     -0.7684755100000000, &
    & -0.4190927700000000,      2.4009849000000001,      0.7110209900000000, &
    & -0.7684755100000000,      0.7110209900000000,     -4.0877214000000004, &
    &  1.9962358000000000,     -0.7781293400000000,     -0.0761197940000000, &
    & -0.7781293400000000,      2.7162259999999998,     -0.3487155600000000, &
    & -0.0761197940000000,     -0.3487155600000000,     -4.7124617999999998, &
    & -2.4067390999999998,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,     -1.6809822999999999,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,      4.0877214000000004, &
    & -2.4067390999999998,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,     -1.6809822999999999,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,      4.0877214000000004, &
    &  1.6764927999999999,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,      2.4584896999999999,     -0.0628446440000000, &
    &  1.0000000000000000,     -0.0628446440000000,     -4.1349824999999996, &
    &  1.6122532000000001,     -0.3760898600000000,      1.0000000000000000, &
    & -0.3760898600000000,      2.4754681999999999,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,     -4.0877214000000004, &
    &  1.6122532000000001,     -0.3760898600000000,     -1.0000000000000000, &
    & -0.3760898600000000,      2.4754681999999999,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,     -4.0877214000000004, &
    & -2.8157193000000000,     -0.1779585300000000,      0.7116558600000000, &
    & -0.1779585300000000,     -1.8760264000000000,      0.7438258300000000, &
    &  0.7116558600000000,      0.7438258300000000,      4.6917457999999996, &
    &  1.7677229000000001,     -1.0000000000000000,     -0.7377517000000000, &
    & -1.0000000000000000,      2.7901045999999998,     -0.2237601200000000, &
    & -0.7377517000000000,     -0.2237601200000000,     -4.5578275000000001, &
    &  1.6845969000000001,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,      2.7414003000000000,     -0.2634078000000000, &
    & -1.0000000000000000,     -0.2634078000000000,     -4.4259972000000003, &
    &  1.6286773000000001,     -0.4639611400000000,      1.0000000000000000, &
    & -0.4639611400000000,      2.7305909000000002,     -0.2079223400000000, &
    &  1.0000000000000000,     -0.2079223400000000,     -4.3592683000000001, &
    &  1.6844881000000000,     -1.0000000000000000,     -0.7530396000000000, &
    & -1.0000000000000000,      2.8391886999999998,     -0.2316152600000000, &
    & -0.7530396000000000,     -0.2316152600000000,     -4.5236767999999996, &
    &  1.6768235000000000,     -1.0000000000000000,     -0.2199620200000000, &
    & -1.0000000000000000,      2.8315239999999999,     -0.3440507400000000, &
    & -0.2199620200000000,     -0.3440507400000000,     -4.5083473999999999, &
    &  1.6257071999999999,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,      2.7804077000000000,     -0.2675726300000000, &
    & -1.0000000000000000,     -0.2675726300000000,     -4.4061149000000004, &
    &  1.5059898000000000,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,      2.6606903000000002,     -0.2443991400000000, &
    & -1.0000000000000000,     -0.2443991400000000,     -4.1666800999999998, &
    &  1.5059898000000000,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,      2.6606903000000002,      0.2443991400000000, &
    &  1.0000000000000000,      0.2443991400000000,     -4.1666800999999998, &
    &  1.4348236999999999,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,      2.5895242000000001,     -0.0786792310000000, &
    &  1.0000000000000000,     -0.0786792310000000,     -4.0243479000000004, &
    &  1.4064007000000001,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,      2.5611012000000000,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,     -3.9675018999999998, &
    &  1.3630085000000001,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,      2.5177090999999998,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,     -3.8807176000000001, &
    &  1.3630085000000001,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,      2.5177090999999998,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,     -3.8807176000000001, &
    &  1.1599136999999999,     -1.0000000000000000,      0.0502257340000000, &
    & -1.0000000000000000,      2.3146141999999998,     -1.1257026999999999, &
    &  0.0502257340000000,     -1.1257026999999999,     -3.4745279000000000, &
    &  1.1599136999999999,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,      2.3146141999999998,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,     -3.4745279000000000, &
    &  1.1599136999999999,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,      2.3146141999999998,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,     -3.4745279000000000, &
    &  1.1599136999999999,     -1.0000000000000000,     -0.4244793700000000, &
    & -1.0000000000000000,      2.3146141999999998,      0.9096272500000000, &
    & -0.4244793700000000,      0.9096272500000000,     -3.4745279000000000, &
    &  1.1531009000000001,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,      2.3078013999999998,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,     -3.4609022999999999, &
    &  1.1531009000000001,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,      2.3078013999999998,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,     -3.4609022999999999, &
    &  1.1531009000000001,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,      2.3078013999999998,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,     -3.4609022999999999, &
    &  1.1326072000000000,     -1.0000000000000000,      0.0000000000000000, &
    & -1.0000000000000000,      2.2873077999999998,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,     -3.4199150000000000, &
    &  1.0101195999999999,      1.0000000000000000,      0.0000000000000000, &
    &  1.0000000000000000,      2.1648201999999999,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,     -3.1749398000000002, &
    &  1.0101195999999999,     -1.0000000000000000,      0.0000000000000000, &
    & -1.0000000000000000,      2.1648201999999999,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,     -3.1749398000000002, &
    &  1.0101195999999999,      1.0000000000000000,      0.0000000000000000, &
    &  1.0000000000000000,      2.1648201999999999,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,     -3.1749398000000002, &
    & -2.1648201999999999,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,     -1.0101195999999999,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,      3.1749398000000002, &
    & -2.1648201999999999,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,     -1.0101195999999999,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,      3.1749398000000002, &
    & -2.3078013999999998,      1.0000000000000000,      0.0000000000000000, &
    &  1.0000000000000000,     -1.1531009000000001,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,      3.4609022999999999, &
    & -2.3078013999999998,      1.0000000000000000,      0.0000000000000000, &
    &  1.0000000000000000,     -1.1531009000000001,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,      3.4609022999999999, &
    & -2.5177090999999998,      1.0000000000000000,     -1.0000000000000000, &
    &  1.0000000000000000,     -1.3630085000000001,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,      3.8807176000000001, &
    & -2.5177090999999998,      1.0000000000000000,      1.0000000000000000, &
    &  1.0000000000000000,     -1.3630085000000001,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,      3.8807176000000001, &
    & -2.5497203000000002,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,     -1.3950198000000000,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,      3.9447402000000000, &
    & -2.5497203000000002,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,     -1.3950198000000000,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,      3.9447402000000000, &
    & -2.5497204000000000,     -1.0000000000000000,      0.0000000000000000, &
    & -1.0000000000000000,     -1.3950198000000000,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,      3.9447402000000000, &
    & -2.5497204000000000,     -1.0000000000000000,      0.0000000000000000, &
    & -1.0000000000000000,     -1.3950198000000000,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,      3.9447402000000000, &
    & -2.6212110000000002,     -1.0000000000000000,      0.4999999800000000, &
    & -1.0000000000000000,     -1.4665104000000000,     -0.8660254100000000, &
    &  0.4999999800000000,     -0.8660254100000000,      4.0877214000000004, &
    & -2.6212110000000002,     -1.0000000000000000,     -0.4999999800000000, &
    & -1.0000000000000000,     -1.4665104000000000,      0.8660254100000000, &
    & -0.4999999800000000,      0.8660254100000000,      4.0877214000000004, &
    & -2.6303709000000000,      1.0000000000000000,      0.6505120300000000, &
    &  1.0000000000000000,     -1.4756704000000000,      0.7791272400000000, &
    &  0.6505120300000000,      0.7791272400000000,      4.1060413000000002, &
    & -2.6477395000000001,      1.0000000000000000,     -0.6855385400000000, &
    &  1.0000000000000000,     -1.4930390000000000,     -0.7589046800000000, &
    & -0.6855385400000000,     -0.7589046800000000,      4.1407784999999997, &
    & -2.6519420000000000,      1.0000000000000000,      0.7292981400000000, &
    &  1.0000000000000000,     -1.4972414999999999,      0.7336400600000000, &
    &  0.7292981400000000,      0.7336400600000000,      4.1491835000000004, &
    & -2.6535365000000000,      1.0000000000000000,     -0.7260825500000000, &
    &  1.0000000000000000,     -1.4988360000000001,     -0.7354965800000000, &
    & -0.7260825500000000,     -0.7354965800000000,      4.1523725999999996, &
    & -2.6606903000000002,      1.0000000000000000,     -0.6950136200000000, &
    &  1.0000000000000000,     -1.5059898000000000,     -0.7264312900000000, &
    & -0.6950136200000000,     -0.7264312900000000,      4.1666800999999998, &
    & -2.6606903000000002,      1.0000000000000000,      0.6961355000000000, &
    &  1.0000000000000000,     -1.5059898000000000,      0.7276038900000000, &
    &  0.6961355000000000,      0.7276038900000000,      4.1666800999999998, &
    & -2.9026092000000001,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,     -1.7479087000000000,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,      4.6505178999999996, &
    & -2.9026092000000001,     -1.0000000000000000,     -1.0000000000000000, &
    & -1.0000000000000000,     -1.7479087000000000,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,      4.6505178999999996, &
    & -2.9026092000000001,     -1.0000000000000000,     -0.9999999899999999, &
    & -1.0000000000000000,     -1.7479087000000000,     -0.5773502700000001, &
    & -0.9999999899999999,     -0.5773502700000001,      4.6505178999999996, &
    & -2.9026092000000001,     -1.0000000000000000,      0.9999999899999999, &
    & -1.0000000000000000,     -1.7479087000000000,      0.5773502700000001, &
    &  0.9999999899999999,      0.5773502700000001,      4.6505178999999996, &
    & -2.9026092000000001,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,     -1.7479087000000000,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,      4.6505178999999996, &
    & -2.9026092000000001,     -1.0000000000000000,      1.0000000000000000, &
    & -1.0000000000000000,     -1.7479087000000000,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,      4.6505178999999996, &
    & -3.0455904999999999,     -1.0000000000000000,      0.0000000000000000, &
    & -1.0000000000000000,     -1.8908900000000000,      0.0000000000000000, &
    &  0.0000000000000000,      0.0000000000000000,      4.9364803999999998, &
    & -2.6745150999999998,      0.9760548000000000,     -0.7116558600000000, &
    &  0.9760548000000000,     -1.4921650000000000,     -0.7438258300000000, &
    & -0.7116558600000000,     -0.7438258300000000,      4.1666800999999998, &
    & -2.9370479000000000,     -0.3881058300000000,     -0.7116558600000000, &
    & -0.3881058300000000,     -1.7546978000000000,     -0.7438258300000000, &
    & -0.7116558600000000,     -0.7438258300000000,      4.6917457999999996, &
    & -2.6844529000000001,      0.9588419900000000,      0.7116558600000000, &
    &  0.9588419900000000,     -1.4822272000000001,      0.7438258300000000, &
    &  0.7116558600000000,      0.7438258300000000,      4.1666800999999998, &
    &  1.7144132000000001,     -0.4578808100000000,     -1.0000000000000000, &
    & -0.4578808100000000,      2.9773325000000002,     -0.2443991400000000, &
    & -1.0000000000000000,     -0.2443991400000000,     -4.6917457999999996, &
    &  1.7194469999999999,     -0.5077001799999999,     -0.9309517700000000, &
    & -0.5077001799999999,      2.9891972999999998,     -0.2275238100000000, &
    & -0.9309517700000000,     -0.2275238100000000,     -4.7086442000000002, &
    &  1.6438908999999999,     -0.5453560500000000,     -1.0000000000000000, &
    & -0.5453560500000000,      2.9669104000000002,     -0.1342776300000000, &
    & -1.0000000000000000,     -0.1342776300000000,     -4.6108013000000003, &
    &  1.5956191000000000,     -0.8438615600000000,     -0.1655457000000000, &
    & -0.8438615600000000,      2.9306128000000000,     -0.3781731900000000, &
    & -0.1655457000000000,     -0.3781731900000000,     -4.5262319000000000, &
    &  1.6497605000000000,     -0.5919432300000000,     -1.0000000000000000, &
    & -0.5919432300000000,      3.0030958000000001,     -0.2158866000000000, &
    & -1.0000000000000000,     -0.2158866000000000,     -4.6528562999999998, &
    &  1.6591488000000001,     -0.5453560500000000,     -1.0000000000000000, &
    & -0.5453560500000000,      3.0133475999999999,     -0.2182111200000000, &
    & -1.0000000000000000,     -0.2182111200000000,     -4.6724962999999997, &
    & -2.7284468999999998,     -0.8142619200000000,      0.0000000000000000, &
    & -0.8142619200000000,     -1.3592744999999999,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,      4.0877214000000004, &
    & -2.7284468999999998,     -0.8142619200000000,      0.0000000000000000, &
    & -0.8142619200000000,     -1.3592744999999999,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,      4.0877214000000004, &
    & -3.1067741999999998,     -0.6820803100000000,     -0.7116558600000000, &
    & -0.6820803100000000,     -1.5849716000000000,     -0.7438258300000000, &
    & -0.7116558600000000,     -0.7438258300000000,      4.6917457999999996, &
    & -3.1067741999999998,     -0.6820803100000000,      0.7116558600000000, &
    & -0.6820803100000000,     -1.5849716000000000,      0.7438258300000000, &
    &  0.7116558600000000,      0.7438258300000000,      4.6917457999999996, &
    &  1.5016307000000000,     -0.6380901900000000,     -1.0000000000000000, &
    & -0.6380901900000000,      3.0742286999999999,     -0.2320154900000000, &
    & -1.0000000000000000,     -0.2320154900000000,     -4.5758593999999997, &
    &  1.0508058000000000,     -0.6220392500000000,      0.1362340200000000, &
    & -0.6220392500000000,      2.6419378000000000,     -1.0760457999999999, &
    &  0.1362340200000000,     -1.0760457999999999,     -3.6927436999999999, &
    &  1.4618940000000000,     -0.5453560500000000,     -1.0000000000000000, &
    & -0.5453560500000000,      3.1288440999999998,     -0.2443991400000000, &
    & -1.0000000000000000,     -0.2443991400000000,     -4.5907380000000000, &
    &  1.4547916999999999,     -0.5471527700000000,     -1.0000000000000000, &
    & -0.5471527700000000,      3.1323951999999999,     -0.2452043300000000, &
    & -1.0000000000000000,     -0.2452043300000000,     -4.5871868999999998, &
    &  1.4534085999999999,     -0.5453560500000000,     -1.0000000000000000, &
    & -0.5453560500000000,      3.1330868000000001,     -0.2443991400000000, &
    & -1.0000000000000000,     -0.2443991400000000,     -4.5864953999999996, &
    &  1.3537442000000000,     -0.5453560500000000,     -1.0000000000000000, &
    & -0.5453560500000000,      3.0334222999999998,     -0.0123171000000000, &
    & -1.0000000000000000,     -0.0123171000000000,     -4.3871665000000002, &
    &  1.2435010000000000,      0.5453560500000000,      1.0000000000000000, &
    &  0.5453560500000000,      2.9231791000000000,      0.2443991400000000, &
    &  1.0000000000000000,      0.2443991400000000,     -4.1666800999999998, &
    &  1.1005197000000000,      0.5453560500000000,      1.0000000000000000, &
    &  0.5453560500000000,      2.7801979000000001,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,     -3.8807176000000001, &
    &  0.9357977900000000,     -0.5227624700000000,      0.0000000000000000, &
    & -0.5227624700000000,      2.6415647999999998,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,     -3.5773625999999998, &
    &  1.2104325000000000,     -0.5023988600000000,      1.0000000000000000, &
    & -0.5023988600000000,      2.9397133000000002,     -0.2251480400000000, &
    &  1.0000000000000000,     -0.2251480400000000,     -4.1501459000000001, &
    &  1.0065039000000000,     -0.4685726900000000,     -0.8592050800000000, &
    & -0.4685726900000000,      2.7748438000000002,      0.6586382500000000, &
    & -0.8592050800000000,      0.6586382500000000,     -3.7813477000000000, &
    &  0.9978762799999999,     -0.4386858500000000,      1.0000000000000000, &
    & -0.4386858500000000,      2.8007266000000000,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,     -3.7986027999999998, &
    &  0.9978762799999999,     -0.4386858500000000,     -1.0000000000000000, &
    & -0.4386858500000000,      2.8007266000000000,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,     -3.7986027999999998, &
    &  1.1502695999999999,     -0.3789152600000000,     -1.0000000000000000, &
    & -0.3789152600000000,      3.0221371000000001,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,     -4.1724066999999998, &
    &  0.9286395100000000,      0.0000000000000000,     -1.0000000000000000, &
    &  0.0000000000000000,      3.2380406000000002,     -0.2443991400000000, &
    & -1.0000000000000000,     -0.2443991400000000,     -4.1666800999999998, &
    &  0.9286395100000000,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,      3.2380406000000002,      0.0000000000000000, &
    &  0.0000000000000000,      0.0000000000000000,     -4.1666800999999998, &
    &  0.8236856900000000,      0.0000000000000000,      1.0000000000000000, &
    &  0.0000000000000000,      3.1330868000000001,      0.0000000000000000, &
    &  1.0000000000000000,      0.0000000000000000,     -3.9567724000000002, &
    &  0.7856582700000000,      0.0000000000000000,     -1.0000000000000000, &
    &  0.0000000000000000,      3.0950593000000000,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,     -3.8807176000000001, &
    &  0.5757506200000000,      0.0000000000000000,      1.0000000000000000, &
    &  0.0000000000000000,      2.8851517000000002,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,     -3.4609022999999999, &
    &  0.5757506200000000,      0.0000000000000000,     -1.0000000000000000, &
    &  0.0000000000000000,      2.8851517000000002,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,     -3.4609022999999999, &
    &  0.5757506200000000,      0.0000000000000000,      1.0000000000000000, &
    &  0.0000000000000000,      2.8851517000000002,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,     -3.4609022999999999, &
    &  0.4327693800000000,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,      2.7421704000000000,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,     -3.1749398000000002, &
    &  0.4327693800000000,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,      2.7421704000000000,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,     -3.1749398000000002, &
    & -2.7421704000000000,      0.0000000000000000,      1.0000000000000000, &
    &  0.0000000000000000,     -0.4327693700000000,     -0.5773502700000001, &
    &  1.0000000000000000,     -0.5773502700000001,      3.1749398000000002, &
    & -2.7421704000000000,      0.0000000000000000,     -1.0000000000000000, &
    &  0.0000000000000000,     -0.4327693700000000,      0.5773502700000001, &
    & -1.0000000000000000,      0.5773502700000001,      3.1749398000000002, &
    & -2.8851517000000002,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,     -0.5757506200000000,     -1.1547004999999999, &
    &  0.0000000000000000,     -1.1547004999999999,      3.4609022999999999, &
    & -2.8851517000000002,      0.0000000000000000,      0.0000000000000000, &
    &  0.0000000000000000,     -0.5757506200000000,      1.1547004999999999, &
    &  0.0000000000000000,      1.1547004999999999,      3.4609022999999999, &
    & -3.0950593000000000,      0.0000000000000000,     -1.0000000000000000, &
    &  0.0000000000000000,     -0.7856582600000001,     -0.5773502700000001, &
    & -1.0000000000000000,     -0.5773502700000001,      3.8807176000000001, &
    & -3.0950593000000000,      0.0000000000000000,      1.0000000000000000, &
    &  0.0000000000000000,     -0.7856582600000001,      0.5773502700000001, &
    &  1.0000000000000000,      0.5773502700000001,      3.8807176000000001, &
    & -3.2380406000000002,      0.0000000000000000,     -0.7116558600000000, &
    &  0.0000000000000000,     -0.9286395100000000,     -0.7438258300000000, &
    & -0.7116558600000000,     -0.7438258300000000,      4.1666800999999998, &
    & -3.2380406000000002,      0.0000000000000000,      0.7116558600000000, &
    &  0.0000000000000000,     -0.9286395100000000,      0.7438258300000000, &
    &  0.7116558600000000,      0.7438258300000000,      4.1666800999999998  &
    &/)
    !
    REAL(RK), PARAMETER, DIMENSION(3, 3, 240) :: &
        &   VERT_HCP = RESHAPE(SOURCE=VERTICES_HCP_DAT, SHAPE=(/3, 3, 240/))
    !
    !---------------------------------------------------------------------------
    !
    ! Check CLASS type and return appropriate tensors.
    IF (CTYPE .EQ. 1) THEN ! CLASS_FCC
        !
        VERTICES = VERT_FCC
        !
    ELSE IF (CTYPE .EQ. 2) THEN ! CLASS_BCC
        !
        VERTICES = VERT_FCC
        !
    ELSE IF (CTYPE .EQ. 3) THEN ! CLASS_HCP
        !
        VERTICES = VERT_HCP
        !
    ELSE IF (CTYPE .EQ. 4) THEN ! CLASS_BCT
        !
        VERTICES = VERT_HCP
        !
    END IF
    !
    !---------------------------------------------------------------------------
    !
    END SUBROUTINE GET_VERTICES
    !
END MODULE CRYSTAL_TYPE_MOD
