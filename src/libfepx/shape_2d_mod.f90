! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE SHAPE_2D_MOD
!
! Shape function library for "2D" elements (faces of 3D elements). Contains full
!   deescriptions of 3-node (linear) and 6-node (quadratic) triangular
!   elements, and 4-node (linear), 8-node (quadratic) and 9-node quadrilateral
!   elements.
!
! Contains subroutines:
! SF2D: Evaluates 2D shape functions at a list of points
! SF2DG: Evaluates 2D shape function gradients at a list of points
! SF2D03: 3-node triangular shape function values
! SF2DG03: 3-node triangular shape function gradients
! SF2D04: 4-node quadrilateral shape function values
! SF2DG04: 4-node quadrilateral shape function gradients
! SF2D06: 6-node triangular shape function values
! SF2DG06: 6-node triangular shape function gradients
! SF2D08: 8-node quadrilateral shape function values
! SF2DG08: 8-node quadrilateral shape function gradients
! SF2D09: 9-node quadrilateral shape function values
! SF2DG09: 9-node quadrilateral shape function gradients
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
!
IMPLICIT NONE
!
CONTAINS
    !
    SUBROUTINE SF2D(ITYPE, NPTS, PTS, VALS, MV1, ISTAT)
    !
    ! Evaluate shape functions at a list of points.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! ITYPE: Element type
    ! NPTS: Number of points
    ! PTS: Points
    ! VALS: Point values
    ! MV1:
    ! ISTAT: Error flagging
    !
    INTEGER :: ITYPE
    INTEGER :: NPTS
    REAL(RK) :: PTS(2, *)
    REAL(RK) :: VALS(MV1, *)
    INTEGER :: MV1
    INTEGER :: ISTAT
    !
    !---------------------------------------------------------------------------
    !
    ISTAT = 0
    !
    IF (ITYPE .EQ. 3) THEN
        !
        CALL SF2D03(NPTS, PTS, VALS, MV1)
        !
    ELSE IF (ITYPE .EQ. 4) THEN
        !
        CALL SF2D04(NPTS, PTS, VALS, MV1)
        !
    ELSE IF (ITYPE .EQ. 6) THEN
        !
        CALL SF2D06(NPTS, PTS, VALS, MV1)
        !
    ELSE IF (ITYPE .EQ. 8) THEN
        !
        CALL SF2D08(NPTS, PTS, VALS, MV1)
        !
    ELSE IF (ITYPE .EQ. 9) THEN
        !
        CALL SF2D09(NPTS, PTS, VALS, MV1)
        !
    ELSE
        !
        ! Unrecognized element type.
        !
        ISTAT = 1
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE SF2D
    !
    !===========================================================================
    !
    SUBROUTINE SF2DG(ITYPE, NPTS, PTS, GRADS, MG1, ISTAT)
    !
    ! Evaluate shape function gradients at a list of points.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! ITYPE: Element type
    ! NPTS: Number of points
    ! PTS: Points
    ! GRADS: Point gradients
    ! MG1:
    ! ISTAT: Error flagging
    !
    INTEGER :: ITYPE
    INTEGER :: NPTS
    REAL(RK) :: PTS(2, *)
    REAL(RK) :: GRADS(2, MG1, *)
    INTEGER :: MG1
    INTEGER :: ISTAT
    !
    !---------------------------------------------------------------------------
    !
    ISTAT = 0
    !
    IF (ITYPE .EQ. 3) THEN
        !
        CALL SF2DG03(NPTS, PTS, GRADS, MG1)
        !
    ELSE IF (ITYPE .EQ. 4) THEN
        !
        CALL SF2DG04(NPTS, PTS, GRADS, MG1)
        !
    ELSE IF (ITYPE .EQ. 6) THEN
        !
        CALL SF2DG06(NPTS, PTS, GRADS, MG1)
        !
    ELSE IF (ITYPE .EQ. 8) THEN
        !
        CALL SF2DG08(NPTS, PTS, GRADS, MG1)
        !
    ELSE IF (ITYPE .EQ. 9) THEN
        !
        CALL SF2DG09(NPTS, PTS, GRADS, MG1)
        !
    ELSE
        !
        ! Unrecognized element type.
        !
        ISTAT = 1
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE SF2DG
    !
    !===========================================================================
    !
    SUBROUTINE SF2D03(NPTS, PNTS, VALUE, MV1)
    !
    ! 3-node triangular shape function values
    !
    !---------------------------------------------------------------------------
    !     
    ! Arguments:
    ! NPTS: Number of points
    ! PNTS: Points
    ! VALUE: Point values
    ! MV1:
    !
    INTEGER :: NPTS
    REAL(RK) :: PNTS(2, *)
    REAL(RK) :: VALUE(MV1, *)
    INTEGER :: MV1
    !
    ! Locals:
    !
    INTEGER :: I
    REAL(RK) :: XI
    REAL(RK) :: ETA
    REAL(RK) :: ZETA
    !
    !---------------------------------------------------------------------------
    !
    DO I = 1, NPTS
        !
        XI = PNTS(1, I)
        ETA = PNTS(2, I)
        ZETA = 1.0D0 - XI - ETA
        !
        ! Nodal locations:
        !
        ! 2
        ! 31
        !
        VALUE(1, I) = XI
        VALUE(2, I) = ETA
        VALUE(3, I) = ZETA
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SF2D03
    !
    !===========================================================================
    !
    SUBROUTINE SF2DG03(NPTS, PNTS, GRAD, MG1)
    !
    ! 3-node triangular shape function gradients
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NPTS: Number of points
    ! PNTS: Points
    ! GRAD: Point gradients
    ! MG1:
    !
    INTEGER :: NPTS
    REAL(RK) :: PNTS(2, *)
    REAL(RK) :: GRAD(2, MG1, *)
    INTEGER :: MG1
    !
    ! Locals:
    !
    INTEGER :: I
    REAL(RK) :: XI
    REAL(RK) :: ETA
    REAL(RK) :: ZETA
    !
    !---------------------------------------------------------------------------
    !
    DO I = 1, NPTS
        !
        XI = PNTS(1, I)
        ETA = PNTS(2, I)
        ZETA = 1.0D0 - XI - ETA
        !
        GRAD(1, 1, I) = 1.0D0
        GRAD(1, 2, I) = 0.0D0
        GRAD(1, 3, I) = -1.0D0
        !
        GRAD(2, 1, I) = 0.0D0
        GRAD(2, 2, I) = 1.0D0
        GRAD(2, 3, I) = -1.0D0
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SF2DG03
    !
    !===========================================================================
    !
    SUBROUTINE SF2D04(NPTS, PNTS, VALUE, MV1)
    !
    ! 4-node quadrilateral shape function values
    !     
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NPTS: Number of points
    ! PNTS: Points
    ! VALUE: Point values
    ! MV1:
    !
    INTEGER :: NPTS
    REAL(RK) :: PNTS(2, *)
    REAL(RK) :: VALUE(MV1, *)
    INTEGER :: MV1
    !
    ! Locals:
    !
    INTEGER :: I
    REAL(RK) :: XI
    REAL(RK) :: ETA
    !
    !---------------------------------------------------------------------------
    !
    DO I = 1, NPTS
        !
        XI  = PNTS(1, I)
        ETA = PNTS(2, I)
        !
        VALUE(1, I) = (1.0D0 - XI) * (1.0D0 - ETA)
        VALUE(2, I) = XI * (1.0D0 - ETA)
        VALUE(3, I) = XI * ETA
        VALUE(4, I) = (1.0D0 - XI) * ETA
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SF2D04
    !
    !===========================================================================
    !
    SUBROUTINE SF2DG04(NPTS, PNTS, GRAD, MG1)
    !
    ! 4-node quadrilateral shape function gradients
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NPTS: Number of points
    ! PNTS: Points
    ! GRAD: Point gradients
    ! MG1:
    !
    INTEGER :: NPTS
    REAL(RK) :: PNTS(2,*)
    REAL(RK) :: GRAD(2,MG1,*)
    INTEGER :: MG1
    !
    ! Locals:
    !
    INTEGER :: I
    REAL(RK) :: XI
    REAL(RK) :: ETA
    !
    !---------------------------------------------------------------------------
    !
    DO I = 1, NPTS
        !
        XI  = PNTS(1, I)
        ETA = PNTS(2, I)
        !
        GRAD(1,1, I) = ETA - 1.0D0
        GRAD(1,2, I) = 1.0D0 - ETA
        GRAD(1,3, I) = ETA
        GRAD(1,4, I) = -ETA
        !
        GRAD(2,1, I) = XI - 1.0D0
        GRAD(2,2, I) = -XI
        GRAD(2,3, I) = XI
        GRAD(2,4, I) = 1.0D0 - XI
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SF2DG04
    !
    !===========================================================================
    !
    SUBROUTINE SF2D06(NPTS, PNTS, VALUE, MV1)
    !
    ! 6-node triangular shape function values
    !     
    !---------------------------------------------------------------------------
    !
    ! Arguements:
    ! NPTS: Number of points
    ! PNTS: Points
    ! VALUE: Point values
    ! MV1:
    !
    INTEGER :: NPTS
    REAL(RK) :: PNTS(2, *)
    REAL(RK) :: VALUE(MV1, *)
    INTEGER :: MV1
    !
    ! Locals:
    !
    INTEGER :: I
    REAL(RK) :: XI
    REAL(RK) :: ETA
    REAL(RK) :: ZETA
    !
    !---------------------------------------------------------------------------
    !
    DO I = 1, NPTS
        !
        XI = PNTS(1, I)
        ETA = PNTS(2, I)
        ZETA = 1.0D0 - XI - ETA
        !
        ! Nodal locations:
        ! 3
        ! 42
        ! 561
        !
        VALUE(1, I) = (2.0D0 * XI - 1.0D0) * XI
        VALUE(2, I) = 4.0D0 * ETA * XI
        VALUE(3, I) = (2.0D0 * ETA - 1.0D0) * ETA
        VALUE(4, I) = 4.0D0 * ETA * ZETA
        VALUE(5, I) = (2.0D0 * ZETA - 1.0D0) * ZETA
        VALUE(6, I) = 4.0D0 * XI * ZETA
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SF2D06
    !
    !===========================================================================
    !
    SUBROUTINE SF2DG06(NPTS, PNTS, GRAD, MG1)
    !
    ! 6-node triangular shape function gradients
    !     
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NPTS: Number of points
    ! PNTS: Points
    ! GRAD: Point gradients
    ! MG1:
    !
    INTEGER :: NPTS
    REAL(RK) :: PNTS(2, *)
    REAL(RK) :: GRAD(2, MG1, *)
    INTEGER :: MG1
    !
    ! Locals:
    !
    INTEGER :: I
    REAL(RK) :: XI
    REAL(RK) :: ETA
    REAL(RK) :: ZETA
    !
    !---------------------------------------------------------------------------
    !
    DO I = 1, NPTS
        !
        XI = PNTS(1, I)
        ETA = PNTS(2, I)
        ZETA = 1.0D0 - XI - ETA
        !
        GRAD(1, 1, I) = 4.0D0 * XI - 1.0D0
        GRAD(1, 2, I) = 4.0D0 * ETA
        GRAD(1, 3, I) = 0.0D0
        GRAD(1, 4, I) = -4.0D0 * ETA
        GRAD(1, 5, I) = -4.0D0 * ZETA + 1.0D0
        GRAD(1, 6, I) = 4.0D0 * ZETA - 4.0D0 * XI
        !
        GRAD(2, 1, I) = 0.0D0
        GRAD(2, 2, I) = 4.0D0 * XI
        GRAD(2, 3, I) = 4.0D0 * ETA - 1.0D0
        GRAD(2, 4, I) = 4.0D0 * ZETA - 4.0D0 * ETA
        GRAD(2, 5, I) = -4.0D0 * ZETA + 1.0D0
        GRAD(2, 6, I) = -4.0D0 * XI
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SF2DG06
    !
    !===========================================================================
    !
    SUBROUTINE SF2D08(NPTS, PNTS, VALUE, MV1)
    !
    ! 8-node quadrilateral shape function values
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NPTS: Number of points
    ! PNTS: Points
    ! VALUE: Point values
    ! MV1:
    !
    INTEGER :: NPTS
    REAL(RK) :: PNTS(2, *)
    REAL(RK) :: VALUE(MV1, *)
    INTEGER :: MV1
    !
    ! Locals:
    !
    INTEGER :: I
    REAL(RK) :: XI
    REAL(RK) :: ETA
    !
    !---------------------------------------------------------------------------
    !
    DO I = 1, NPTS
        !
        XI = PNTS(1, I)
        ETA = PNTS(2, I)
        !
        VALUE(1, I) = (1.0D0 - XI) * (1.0D0 - ETA) * (-2.0D0 * XI - 2.0D0 &
            & * ETA + 1.0D0)
        VALUE(2, I) = 4.0D0 * XI * (1.0D0 - XI) * (1.0D0 - ETA)
        VALUE(3, I) = XI * (1.0D0 - ETA) * (2.0D0 * XI - 2.0D0 * ETA - &
            & 1.0D0)
        VALUE(4, I) = 4.0D0 * XI * ETA * (1.0D0 - ETA)
        VALUE(5, I) = XI * ETA * (2.0D0 * XI + 2.0D0 * ETA - 3.0D0)
        VALUE(6, I) = 4.0D0 * XI * ETA * (1.0D0 - XI)
        VALUE(7, I) = ETA * (1.0D0 - XI) * (-2.0D0 * XI + 2.0D0 * ETA - &
            & 1.0D0)
        VALUE(8, I) = 4.0D0 * ETA * (1.0D0 - XI) * (1.0D0 - ETA)
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SF2D08
    !
    !===========================================================================
    !
    SUBROUTINE SF2DG08(NPTS, PNTS, GRAD, MG1)
    !
    ! 8-node quadrilateral shape function gradients
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NPTS: Number of points
    ! PNTS: Points
    ! GRAD: Point gradients
    ! MG1:
    !
    INTEGER :: NPTS
    REAL(RK) :: PNTS(2, *)
    REAL(RK) :: GRAD(2, MG1, *)
    INTEGER :: MG1
    !
    ! Locals:
    !
    INTEGER :: I
    REAL(RK) :: XI
    REAL(RK) :: ETA
    !
    !---------------------------------------------------------------------------
    !
    DO I = 1, NPTS
        !
        XI = PNTS(1, I)
        ETA = PNTS(2, I)
        !
        GRAD(1, 1, I) = -1.0D0 * (1.0D0 - ETA) * (3.0D0 - 4.0D0 * XI - &
            & 2.0D0 * ETA)
        GRAD(1, 2, I) = 4.0D0 * (1.0D0 - ETA) * (1.0D0 - 2.0D0 * XI)
        GRAD(1, 3, I) = (1.0D0 - ETA) * (4.0D0 * XI - 2.0D0 * ETA - 1.0D0)
        GRAD(1, 4, I) = 4.0D0 * ETA * (1.0D0 - ETA)
        GRAD(1, 5, I) = ETA * (4.0D0 * XI + 2.0D0 * ETA - 3.0D0)
        GRAD(1, 6, I) = 4.0D0 * ETA * (1.0D0 - 2.0D0 * XI)
        GRAD(1, 7, I) = -ETA * (-4.0D0 * XI + 2.0D0 * ETA + 1.0D0)
        GRAD(1, 8, I) = -4.0D0 * ETA * (1.0D0 - ETA)
        !
        GRAD(2, 1, I) = -1.0D0 * (1.0D0 - XI) * (3.0D0 - 4.0D0 * ETA - &
            & 2.0D0 * XI)
        GRAD(2, 2, I) = -4.0D0 * XI * (1.0D0 - XI)
        GRAD(2, 3, I) = -XI * (2.0D0 * XI - 4.0D0 * ETA + 1.0D0)
        GRAD(2, 4, I) = 4.0D0 * XI * (1.0D0 - 2.0D0 * ETA)
        GRAD(2, 5, I) = XI * (2.0D0 * XI + 4.0D0 * ETA - 3.0D0)
        GRAD(2, 6, I) = 4.0D0 * XI * (1.0D0 - XI)
        GRAD(2, 7, I) = (1.0D0 - XI) * (-2.0D0 * XI + 4.0D0 * ETA - 1.0D0)
        GRAD(2, 8, I) = 4.0D0 * (1.0D0 - XI) * (1.0D0 - 2.0D0 * ETA)
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SF2DG08
    !
    !===========================================================================
    !
    SUBROUTINE SF2D09(NPTS, PNTS, VALUE, MV1)
    !
    ! 9-node quadrilateral shape function values
    ! 
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NPTS: Number of points
    ! PNTS: Points
    ! VALUE: Point values
    ! MV1:
    !
    INTEGER :: NPTS
    REAL(RK) :: PNTS(2, *)
    REAL(RK) :: VALUE(MV1, *)
    INTEGER :: MV1
    !
    ! Locals:
    !
    INTEGER :: I
    REAL(RK) :: XI
    REAL(RK) :: ETA
    !
    !---------------------------------------------------------------------------
    !
    DO I=1, NPTS
        !
        XI  = PNTS(1, I)
        ETA = PNTS(2, I)
        !
        VALUE(1, I) = (2.0D0 * XI - 1.0D0) * (2.0D0 * ETA - 1.0D0) * (XI - &
            & 1.0D0) * (ETA - 1.0D0)
        VALUE(2, I) = 4.0D0 * XI * (1.0D0 - XI) * (ETA - 1.0D0) * (2.0D0 * &
            & ETA - 1.0D0)
        VALUE(3, I) = XI * (2.0D0 * XI - 1.0D0) * (2.0D0 * ETA - 1.0D0) * &
            & (ETA - 1.0D0)
        VALUE(4, I) = 4.0D0 * XI * ETA * (2.0D0 * XI - 1.0D0) * (1.0D0 - &
            & ETA)
        VALUE(5, I) = XI * ETA * (2.0D0 * XI - 1.0D0) * (2.0D0 * ETA - &
            & 1.0D0)
        VALUE(6, I) = 4.0D0 * XI * ETA * (1.0D0 - XI) * (2.0D0 * ETA - &
            & 1.0D0)
        VALUE(7, I) = (XI - 1.0D0) * (2.0D0 * XI - 1.0D0) * (2.0D0 * ETA - &
            & 1.0D0) * ETA
        VALUE(8, I) = 4.0D0 * (2.0D0 * XI - 1.0D0) * (XI - 1.0D0) * &
            & (1.0D0 - ETA) * ETA
        VALUE(9, I) = 16.0D0 * XI * ETA * (1.0D0 - XI) * (1.0D0 - ETA)
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SF2D09
    !
    !===========================================================================
    !
    SUBROUTINE SF2DG09(NPTS, PNTS, GRAD, MG1)
    !
    ! 9-node quadrilateral shape function gradients
    ! 
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NPTS: Number of points
    ! PNTS: Points
    ! GRAD: Point gradients
    ! MG1:
    !
    INTEGER :: NPTS
    REAL(RK) :: PNTS(2, *)
    REAL(RK) :: GRAD(2, MG1, *)
    INTEGER :: MG1
    !
    ! Locals:
    !
    INTEGER :: I
    REAL(RK) :: XI
    REAL(RK) :: ETA
    !
    !---------------------------------------------------------------------------
    !
    DO I = 1, NPTS
        !
        XI = PNTS(1, I)
        ETA = PNTS(2, I)
        !
        GRAD(1, 1, I) = (4.0D0 * XI - 3.0D0) * (2.0D0 * ETA - 1.0D0) * &
            & (ETA - 1.0D0)
        GRAD(1, 2, I) = 4.0D0 * (1.0D0 - 2.0D0 * XI) * (ETA - 1.0D0) * &
            & (2.0D0 * ETA - 1.0D0)
        GRAD(1, 3, I) = (4.0D0 * XI - 1.0D0) * (2.0D0 * ETA - 1.0D0) * &
            & (ETA - 1.0D0)
        GRAD(1, 4, I) = 4.0D0 * ETA * (4.0D0 * XI - 1.0D0) * (1.0D0 - ETA)
        GRAD(1, 5, I) = ETA * (4.0D0 * XI - 1.0D0) * (2.0D0 * ETA - 1.0D0)
        GRAD(1, 6, I) = 4.0D0 * ETA * (1.0D0 - 2.0D0 * XI) * (2.0D0 * ETA &
            & - 1.0D0)
        GRAD(1, 7, I) = (4.0D0 * XI - 3.0D0) * (2.0D0 * ETA - 1.0D0) * ETA
        GRAD(1, 8, I) = 4.0D0 * ETA * (4.0D0 * XI - 3.0D0) * (1.0D0 - ETA)
        GRAD(1, 9, I) = 16.0D0 * (1.0D0 - 2.0D0 * XI) * (1.0D0 - ETA) * ETA
        !
        GRAD(2, 1, I) = (2.0D0 * XI - 1.0D0) * (4.0D0 * ETA - 3.0D0) * (XI &
            & - 1.0D0)
        GRAD(2, 2, I) = 4.0D0 * XI * (1.0D0 - XI) * (4.0D0 * ETA - 3.0D0)
        GRAD(2, 3, I) = XI * (2.0D0 * XI - 1.0D0) * (4.0D0 * ETA - 3.0D0)
        GRAD(2, 4, I) = 4.0D0 * XI * (2.0D0 * XI - 1.0D0) * (1.0D0 - &
            & 2.0D0 * ETA)
        GRAD(2, 5, I) = XI * (2.0D0 * XI - 1.0D0) * (4.0D0 * ETA - 1.0D0)
        GRAD(2, 6, I) = 4.0D0 * XI * (1.0D0 - XI) * (4.0D0 * ETA - 1.0D0)
        GRAD(2, 7, I) = (XI - 1.0D0) * (2.0D0 * XI - 1.0D0) * (4.0D0 * ETA &
            & - 1.0D0)
        GRAD(2, 8, I) = 4.0D0 * (2.0D0 * XI - 1.0D0) * (XI - 1.0D0) * &
            & (1.0D0 - 2.0D0 * ETA)
        GRAD(2, 9, I) = 16.0D0 * XI * (1.0D0 - XI) * (1.0D0 - 2.0D0 * ETA)
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SF2DG09
    !
END MODULE SHAPE_2D_MOD
