! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE SHAPE_3D_MOD
!
! Shape function library for 3D elements. Contains full descriptions of 4-node
!   (linear) and 10-node (quadratic) tetrahedral elements, as well as 8-node
!   cubic elements. Currently only 10-node tetrahedral elements are supported
!   elsewhere in the code.
!
! Contains subroutines:
! SFDER_HPAR: Compute quadrature quantities given a set of local coordinates.
!   For 10-node tetrahedral elements specfically:
! SF10T_EVAL_VEC:Evaluate shape fcns for vector quantity at array of pts.
! T10_SHAPE_HPAR: Shape functions.
! T10_DERIV_HPAR: Shape function derivatives.
!   For 4-node tetrahedral elements specifically:
! SF4T_EVAL_VEC: Evaluate shape fcns. for vector quantity at array of pts.
! T4_SHAPE_HPAR: Shape functions.
! T4_DERIV_HPAR: Shape function derivatives.
!   For 8-node cubic elements specifically:
! SF8B_EVAL_VEC:Evaluate shape fcns for vector quantity at array of pts.
! B8_SHAPE_HPAR: Shape functions.
! B8_DERIV_HPAR: Shape function derivatives.
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
!
! From libfepx:
!
USE READ_INPUT_MOD, ONLY: NNPE, EL_SUB1, EL_SUP1, KDIM1
!
IMPLICIT NONE
!
CONTAINS
    !
    SUBROUTINE SFDER_HPAR(LOC0, LOC1, LOC2, COORDS, DNDX, DNDY, DNDZ, DET, &
        & IJAC11, IJAC12, IJAC13, IJAC21, IJAC22, IJAC23, IJAC31, IJAC32, &
        & IJAC33)
    !
    ! Compute quadrature quantities given a set of local coordinates.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! LOC0,LOC1,LOC2 : Coordinates of point (in reference element)
    ! COORDS : Real coordinates of element's nodes
    ! DET : Determinant of jacobian
    ! DNDX,DNDY,DNDZ : Derivatives of mapped shape functions
    ! IJACxx : Components of inverse jacobian
    !
    REAL(RK), INTENT(IN) :: LOC0
    REAL(RK), INTENT(IN) :: LOC1
    REAL(RK), INTENT(IN) :: LOC2
    REAL(RK), INTENT(IN) :: COORDS(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: DNDX(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: DNDY(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: DNDZ(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: DET(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: IJAC11(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: IJAC12(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: IJAC13(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: IJAC21(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: IJAC22(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: IJAC23(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: IJAC31(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: IJAC32(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: IJAC33(EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    INTEGER :: I0
    INTEGER :: I1
    INTEGER :: I2
    !
    REAL(RK) :: DND1
    REAL(RK) :: DND2
    REAL(RK) :: DND3
    !
    REAL(RK) :: DNDA(0:NNPE)
    REAL(RK) :: DNDB(0:NNPE)
    REAL(RK) :: DNDC(0:NNPE)
    !
    REAL(RK) :: JAC11(EL_SUB1:EL_SUP1)
    REAL(RK) :: JAC12(EL_SUB1:EL_SUP1)
    REAL(RK) :: JAC13(EL_SUB1:EL_SUP1)
    REAL(RK) :: JAC21(EL_SUB1:EL_SUP1)
    REAL(RK) :: JAC22(EL_SUB1:EL_SUP1)
    REAL(RK) :: JAC23(EL_SUB1:EL_SUP1)
    REAL(RK) :: JAC31(EL_SUB1:EL_SUP1)
    REAL(RK) :: JAC32(EL_SUB1:EL_SUP1)
    REAL(RK) :: JAC33(EL_SUB1:EL_SUP1)
    !
    INTEGER :: IO
    !
    !---------------------------------------------------------------------------
    !
    ! Evaluate shape function derivatives.
    !
    ! CALL B8_DERIV_HPAR(LOC0, LOC1, LOC2, DNDA, DNDB, DNDC)
    ! CALL T4_DERIV_HPAR(LOC0, LOC1, LOC2, DNDA, DNDB, DNDC)
    CALL T10_DERIV_HPAR(LOC0, LOC1, LOC2, DNDA, DNDB, DNDC)
    !
    ! Initialize Jacobian matrix
    !
    JAC11 = 0.0D0
    JAC12 = 0.0D0
    JAC13 = 0.0D0
    JAC21 = 0.0D0
    JAC22 = 0.0D0
    JAC23 = 0.0D0
    JAC31 = 0.0D0
    JAC32 = 0.0D0
    JAC33 = 0.0D0
    !
    DO I = 0, NNPE
        !
        I0 = 3 * I
        I1 = I0 + 1
        I2 = I1 + 1
        !
        DND1 = DNDA(I)
        DND2 = DNDB(I)
        DND3 = DNDC(I)
        !
        ! DO J = EL_SUB1,EL_SUP1 (legacy?)
        IJAC11(:) = COORDS(I0, :)
        IJAC22(:) = COORDS(I1, :)
        IJAC33(:) = COORDS(I2, :)
        ! END DO
        !
        JAC11 = JAC11 + IJAC11 * DND1
        JAC21 = JAC21 + IJAC11 * DND2
        JAC31 = JAC31 + IJAC11 * DND3
        JAC12 = JAC12 + IJAC22 * DND1
        JAC22 = JAC22 + IJAC22 * DND2
        JAC32 = JAC32 + IJAC22 * DND3
        JAC13 = JAC13 + IJAC33 * DND1
        JAC23 = JAC23 + IJAC33 * DND2
        JAC33 = JAC33 + IJAC33 * DND3
        !
    END DO
    !
    ! Determinant of the Jacobian matrix
    !
    IJAC11 = JAC11 * JAC22 * JAC33
    IJAC12 = JAC12 * JAC23 * JAC31
    IJAC13 = JAC13 * JAC21 * JAC32
    IJAC21 = JAC11 * JAC23 * JAC32
    IJAC22 = JAC12 * JAC21 * JAC33
    IJAC23 = JAC13 * JAC22 * JAC31
    !
    DET = IJAC11 + IJAC12 + IJAC13
    DET = DET - (IJAC21 + IJAC22 + IJAC23)
    !
    ! Inverse of the Jacobian matrix
    !
    IJAC11 = JAC22 * JAC33 - JAC23 * JAC32
    IJAC21 = JAC23 * JAC31 - JAC21 * JAC33
    IJAC31 = JAC21 * JAC32 - JAC22 * JAC31
    IJAC12 = JAC13 * JAC32 - JAC12 * JAC33
    IJAC22 = JAC11 * JAC33 - JAC13 * JAC31
    IJAC32 = JAC31 * JAC12 - JAC11 * JAC32
    IJAC13 = JAC12 * JAC23 - JAC13 * JAC22
    IJAC23 = JAC13 * JAC21 - JAC11 * JAC23
    IJAC33 = JAC11 * JAC22 - JAC12 * JAC21
    !
    IJAC11 = IJAC11 / DET
    IJAC12 = IJAC12 / DET
    IJAC13 = IJAC13 / DET
    IJAC21 = IJAC21 / DET
    IJAC22 = IJAC22 / DET
    IJAC23 = IJAC23 / DET
    IJAC31 = IJAC31 / DET
    IJAC32 = IJAC32 / DET
    IJAC33 = IJAC33 / DET
    !
    ! Shape function derivatives
    !
    DO I = 0, NNPE
        !
        DND1 = DNDA(I)
        DND2 = DNDB(I)
        DND3 = DNDC(I)
        !
        JAC11 = IJAC11 * DND1 + IJAC12 * DND2 + IJAC13 * DND3
        JAC22 = IJAC21 * DND1 + IJAC22 * DND2 + IJAC23 * DND3
        JAC33 = IJAC31 * DND1 + IJAC32 * DND2 + IJAC33 * DND3
        !
        DNDX(I, :) = JAC11
        DNDY(I, :) = JAC22
        DNDZ(I, :) = JAC33
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SFDER_HPAR
    !
    !===========================================================================
    !
    SUBROUTINE SF10T_EVAL_VEC(VEC, NPT, PNT, VAL)
    !
    ! Evaluate sf functions for vector valued quantity at array of points.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! VEC: Values of vector at nodes.
    ! PNT: Array of points in reference element.
    ! NPT: Number of points
    ! VAL: Values of vector at points
    ! Note: we need explicit interface to use assumed shape arrays as arguments,
    !   i.e. PNT(:), VAL(:); unless this is put in a module or an interface, we
    !   can still use the assumed-size arrays.
    !
    REAL(RK), INTENT(IN) :: VEC(3, 10)
    INTEGER, INTENT(IN) :: NPT
    REAL(RK), INTENT(IN) :: PNT(*)
    REAL(RK), INTENT(OUT) :: VAL(*)
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: XDOF
    INTEGER :: YDOF
    INTEGER :: ZDOF
    REAL(RK) :: X1
    REAL(RK) :: X2
    REAL(RK) :: X3
    REAL(RK) :: SF(10)
    REAL(RK) :: TMPVEC(3)
    !
    !---------------------------------------------------------------------------
    !
    XDOF = -2
    DO I = 1, NPT
        !
        XDOF = XDOF + 3
        YDOF = XDOF + 1
        ZDOF = YDOF + 1
        !
        X1 = PNT(XDOF)
        X2 = PNT(YDOF)
        X3 = PNT(ZDOF)
        !
        SF( 1) = 2.0D0 * (X1 + X2 + X3 - 1.0D0) * (X1 + X2 + X3 - 0.5D0)
        SF( 2) = -4.0D0 * (X1 + X2 + X3 - 1.0D0) * X1
        SF( 3) = 2.0D0 * X1 * (X1 - 0.5D0)
        SF( 4) = 4.0D0 * X2 * X1
        SF( 5) = 2.0D0 * X2 * (X2 - 0.5D0)
        SF( 6) = -4.0D0 * (X1 + X2 + X3 - 1.0D0) * X2
        SF( 7) = -4.0D0 * (X1 + X2 + X3 - 1.0D0) * X3
        SF( 8) = 4.0D0 * X1 * X3
        SF( 9) = 4.0D0 * X2 * X3
        SF(10) = 2.0D0 * X3 * (X3 - 0.5D0)
        !
        TMPVEC = MATMUL(VEC, SF)
        !
        VAL(XDOF) = TMPVEC(1)
        VAL(YDOF) = TMPVEC(2)
        VAL(ZDOF) = TMPVEC(3)
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SF10T_EVAL_VEC
    !
    !===========================================================================
    !
    SUBROUTINE T10_SHAPE_HPAR(LOC0, LOC1, LOC2, SHAPE)
    !
    ! Shape functions for 10-node tetrahedron.
    !
    ! Node ordering: Looking down 3-axis (here's the coordinate system):
    !   2
    !   |
    !   3--1
    !
    ! Top:
    !   10
    !
    ! Middle:
    !   9
    !   7 8
    !
    ! Bottom:
    !   5
    !   6 4
    !   1 2 3
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! LOC0,LOC1,LOC2: Coordinates of point
    ! SHAPE: Shape function value at point, distributed over all elements
    !
    REAL(RK), INTENT(IN) :: LOC0
    REAL(RK), INTENT(IN) :: LOC1
    REAL(RK), INTENT(IN) :: LOC2
    REAL(RK), INTENT(OUT) :: SHAPE(0:NNPE)
    !
    !---------------------------------------------------------------------------
    !
    SHAPE(0) = 2.0D0 * (LOC0 + LOC1 + LOC2 - 1.0D0) * (LOC0 + LOC1 + LOC2 - &
        & 0.5D0)
    SHAPE(1) = -4.0D0 * (LOC0 + LOC1 + LOC2 - 1.0D0) * LOC0
    SHAPE(2) = 2.0D0 * LOC0 * (LOC0 - 0.5D0)
    SHAPE(3) = 4.0D0 * LOC1 * LOC0
    SHAPE(4) = 2.0D0 * LOC1 * (LOC1 - 0.5D0)
    SHAPE(5) = -4.0D0 * (LOC0 + LOC1 + LOC2 - 1.0D0) * LOC1
    SHAPE(6) = -4.0D0 * (LOC0 + LOC1 + LOC2 - 1.0D0) * LOC2
    SHAPE(7) = 4.0D0 * LOC0 * LOC2
    SHAPE(8) = 4.0D0 * LOC1 * LOC2
    SHAPE(9) = 2.0D0 * LOC2 * (LOC2 - 0.5D0)
    !
    RETURN
    !
    END SUBROUTINE T10_SHAPE_HPAR
    !
    !===========================================================================
    !
    SUBROUTINE T10_DERIV_HPAR(LOC0, LOC1, LOC2, DNDA, DNDB, DNDC)
    !
    ! Shape function derivatives for 10-node tetrahedron.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! LOC0,LOC1,LOC2: Coordinates of point
    ! DNDA,DNDB,DNDC: shape function derivatives at point, distributed over all
    !   elements
    !
    REAL(RK), INTENT(IN) :: LOC0
    REAL(RK), INTENT(IN) :: LOC1
    REAL(RK), INTENT(IN) :: LOC2
    REAL(RK), INTENT(OUT) :: DNDA(0:NNPE)
    REAL(RK), INTENT(OUT) :: DNDB(0:NNPE)
    REAL(RK), INTENT(OUT) :: DNDC(0:NNPE)
    !
    !---------------------------------------------------------------------------
    !
    DNDA(0) = 4.0D0 * (LOC0 + LOC1 + LOC2) - 3.0D0
    DNDA(1) = -4.0D0 * (2.0D0 * LOC0 + LOC1 + LOC2 - 1.0D0)
    DNDA(2) = 4.0D0 * LOC0 - 1.0D0
    DNDA(3) = 4.0D0 * LOC1
    DNDA(4) = 0.0D0
    DNDA(5) = -4.0D0 * LOC1
    DNDA(6) = -4.0D0 * LOC2
    DNDA(7) = 4.0D0 * LOC2
    DNDA(8) = 0.0D0
    DNDA(9) = 0.0D0
    !
    DNDB(0) = 4.0D0 * (LOC0 + LOC1 + LOC2) - 3.0D0
    DNDB(1) = -4.0D0 * LOC0
    DNDB(2) = 0.0D0
    DNDB(3) = 4.0D0 * LOC0
    DNDB(4) = 4.0D0 * LOC1 - 1.0D0
    DNDB(5) = -4.0D0 * (LOC0 + 2.0D0 * LOC1 + LOC2 - 1.0D0)
    DNDB(6) = -4.0D0 * LOC2
    DNDB(7) = 0.0D0
    DNDB(8) = 4.0D0 * LOC2
    DNDB(9) = 0.0D0
    !
    DNDC(0) = 4.0D0 * (LOC0 + LOC1 + LOC2) - 3.0D0
    DNDC(1) = -4.0D0 * LOC0
    DNDC(2) = 0.0D0
    DNDC(3) = 0.0D0
    DNDC(4) = 0.0D0
    DNDC(5) = -4.0D0 * LOC1
    DNDC(6) = -4.0D0 * (LOC0 + LOC1 + 2.0D0 * LOC2 - 1.0D0)
    DNDC(7) = 4.0D0 * LOC0
    DNDC(8) = 4.0D0 * LOC1
    DNDC(9) = 4.0D0 * LOC2 - 1.0D0
    !
    RETURN
    !
    END SUBROUTINE T10_DERIV_HPAR
    !
    !===========================================================================
    !
    SUBROUTINE SF4T_EVAL_VEC(VEC, NPT, PNT, VAL)
    !
    ! Evaluate sf functions for vector valued quantity at array of points.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! VEC: Values of vector at nodes.
    ! PNT: Array of points in reference element.
    ! NPT: Number of points
    ! VAL: Values of vector at points
    ! Note: we need explicit interface to use assumed shape arrays as arguments,
    !   i.e. PNT(:), VAL(:); unless this is put in a module or an interface, we
    !   can still use the assumed-size arrays.
    !
    REAL(RK), INTENT(IN) :: VEC(3, 4)
    INTEGER, INTENT(IN) :: NPT
    REAL(RK), INTENT(IN) :: PNT(*)
    REAL(RK), INTENT(OUT):: VAL(*)
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: XDOF
    INTEGER :: YDOF
    INTEGER :: ZDOF
    REAL(RK) :: X1
    REAL(RK) :: X2
    REAL(RK) :: X3
    REAL(RK) :: SF(4)
    REAL(RK) :: TMPVEC(3)
    !
    !---------------------------------------------------------------------------
    !
    XDOF = -2
    DO I = 1, NPT
        !
        XDOF = XDOF + 3
        YDOF = XDOF + 1
        ZDOF = YDOF + 1
        !
        X1 = PNT(XDOF)
        X2 = PNT(YDOF)
        X3 = PNT(ZDOF)
        !
        SF(1) = X1
        SF(2) = X2
        SF(3) = X3
        SF(4) = 1 - (X1 + X2 + X3)
        !
        TMPVEC = MATMUL(VEC, SF)
        !
        VAL(XDOF) = TMPVEC(1)
        VAL(YDOF) = TMPVEC(2)
        VAL(ZDOF) = TMPVEC(3)
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SF4T_EVAL_VEC
    !
    !===========================================================================
    !
    SUBROUTINE T4_SHAPE_HPAR(LOC0, LOC1, LOC2, SHAPE)
    !
    ! Shape functions for 4-node tetrahedron.
    !
    ! Node ordering: Looking down 3-axis (here's the coordinate system):
    !   2
    !   |
    !   3--1
    !
    ! Top:
    !   2
    !
    ! Bottom:
    !   3
    !   4 1
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! LOC0,LOC1,LOC2: Coordinates of point
    ! SHAPE: Shape function value at point, distributed over all elements
    !
    REAL(RK), INTENT(IN) :: LOC0
    REAL(RK), INTENT(IN) :: LOC1
    REAL(RK), INTENT(IN) :: LOC2
    REAL(RK), INTENT(OUT) :: SHAPE(0:NNPE)
    !
    !---------------------------------------------------------------------------
    !
    SHAPE(0) = LOC0
    SHAPE(1) = LOC1
    SHAPE(2) = LOC2
    SHAPE(3) = 1.0 - (LOC0 + LOC1 + LOC2)
    !
    RETURN
    !
    END SUBROUTINE T4_SHAPE_HPAR
    !
    !===========================================================================
    !
    SUBROUTINE T4_DERIV_HPAR(LOC0, LOC1, LOC2, DNDA, DNDB, DNDC)
    !
    ! Shape function derivatives for 4-node tetrahedron.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! LOC0,LOC1,LOC2: Coordinates of point
    ! DNDA,DNDB,DNDC: shape function derivatives at point, distributed over all
    !   elements
    !
    REAL(RK), INTENT(IN) :: LOC0
    REAL(RK), INTENT(IN) :: LOC1
    REAL(RK), INTENT(IN) :: LOC2
    REAL(RK), INTENT(OUT) :: DNDA(0:NNPE)
    REAL(RK), INTENT(OUT) :: DNDB(0:NNPE)
    REAL(RK), INTENT(OUT) :: DNDC(0:NNPE)
    !
    !---------------------------------------------------------------------------
    !
    DNDA(0) = 1.0D0
    DNDA(1) = 0.0D0
    DNDA(2) = 0.0D0
    DNDA(3) = -1.0D0
    !
    DNDB(0) = 0.0D0
    DNDB(1) = 0.0D0
    DNDB(2) = 1.0D0
    DNDB(3) = -1.0D0
    !
    DNDC(0) = 0.0D0
    DNDC(1) = 1.0D0
    DNDC(2) = 0.0D0
    DNDC(3) = -1.0D0
    !
    RETURN
    !
    END SUBROUTINE T4_DERIV_HPAR
    !
    !===========================================================================
    !
    SUBROUTINE SF8B_EVAL_VEC(VEC, NPT, PNT, VAL)
    !
    ! Evaluate sf functions for vector valued quantity at array of points.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments.
    ! VEC: Values of vector at nodes.
    ! PNT: Array of points in reference element.
    ! NPT: Number of points
    ! VAL: Values of vector at points
    ! Note: we need explicit interface to use assumed shape arrays as arguments,
    !   i.e. PNT(:), VAL(:); unless this is put in a module or an interface, we
    !   can still use the assumed-size arrays.
    !
    REAL(RK), INTENT(IN) :: VEC(3, 8)
    INTEGER, INTENT(IN) :: NPT
    REAL(RK), INTENT(IN) :: PNT(*)
    REAL(RK), INTENT(OUT):: VAL(*)
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: XDOF
    INTEGER :: YDOF
    INTEGER :: ZDOF
    REAL(RK) :: X1
    REAL(RK) :: X2
    REAL(RK) :: X3
    REAL(RK) :: SF(8)
    REAL(RK) :: TMPVEC(3)
    !
    !---------------------------------------------------------------------------
    !
    XDOF = -2
    DO I = 1, NPT
        !
        XDOF = XDOF + 3
        YDOF = XDOF + 1
        ZDOF = YDOF + 1
        !
        X1 = PNT(XDOF)
        X2 = PNT(YDOF)
        X3 = PNT(ZDOF)
        !
        SF(1) = -(X3 - 1.0) * (X2 - 1.0) * (X1 - 1.0)
        SF(2) =  (X3 - 1.0) * (X2 - 1.0) * X1
        SF(3) = -(X3 - 1.0) * X2 * X1
        SF(4) =  (X3 - 1.0) * X2 * (X1 - 1.0)
        SF(5) =  X3 * (X2 - 1.0) * (X1 - 1.0)
        SF(6) = -X3 * (X2 - 1.0) * X1
        SF(7) =  X3 * X2 * X1
        SF(8) = -X3 * X2 * (X1 - 1.0)
        !
        TMPVEC = MATMUL(VEC, SF)
        !
        VAL(XDOF) = TMPVEC(1)
        VAL(YDOF) = TMPVEC(2)
        VAL(ZDOF) = TMPVEC(3)
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE SF8B_EVAL_VEC
    !
    !===========================================================================
    !
    SUBROUTINE B8_SHAPE_HPAR(LOC0, LOC1, LOC2, SHAPE)
    !
    !     Shape functions for 8-node brick.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments.
    ! LOC0,LOC1,LOC2: Coordinates of point
    ! SHAPE: Shape function value at point, distributed over all elements
    !
    REAL(RK), INTENT(IN) :: LOC0
    REAL(RK), INTENT(IN) :: LOC1
    REAL(RK), INTENT(IN) :: LOC2
    REAL(RK), INTENT(OUT) :: SHAPE(0:NNPE)
    !
    !---------------------------------------------------------------------------
    !
    SHAPE(0) = -(LOC2 - 1.0) * (LOC1 - 1.0) * (LOC0 - 1.0)
    SHAPE(1) =  (LOC2 - 1.0) * (LOC1 - 1.0) * LOC0
    SHAPE(2) = -(LOC2 - 1.0) * LOC1 * LOC0
    SHAPE(3) =  (LOC2 - 1.0) * LOC1 * (LOC0 - 1.0)
    SHAPE(4) =  LOC2 * (LOC1 - 1.0) * (LOC0 - 1.0)
    SHAPE(5) = -LOC2 * (LOC1 - 1.0) * LOC0
    SHAPE(6) =  LOC2 * LOC1 * LOC0
    SHAPE(7) = -LOC2 * LOC1 * (LOC0 - 1.0)
    !
    RETURN
    !
    END SUBROUTINE B8_SHAPE_HPAR
    !
    !===========================================================================
    !
      SUBROUTINE B8_DERIV_HPAR(LOC0, LOC1, LOC2, DNDA, DNDB, DNDC)
    !
    ! Shape function derivatives for 8-node brick.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! LOC0,LOC1,LOC2: Coordinates of point
    ! DNDA,DNDB,DNDC: Shape function derivatives at point, distributed over all
    !   elements
    !
    REAL(RK), INTENT(IN) :: LOC0
    REAL(RK), INTENT(IN) :: LOC1
    REAL(RK), INTENT(IN) :: LOC2
    REAL(RK), INTENT(OUT) :: DNDA(0:NNPE)
    REAL(RK), INTENT(OUT) :: DNDB(0:NNPE)
    REAL(RK), INTENT(OUT) :: DNDC(0:NNPE)
    !
    !---------------------------------------------------------------------------
    !
    DNDA(0) = -(LOC2 - 1.0) * (LOC1 - 1.0)
    DNDA(1) =  (LOC2 - 1.0) * (LOC1 - 1.0)
    DNDA(2) = -(LOC2 - 1.0) * LOC1
    DNDA(3) =  (LOC2 - 1.0) * LOC1
    DNDA(4) =  LOC2 * (LOC1 - 1.0)
    DNDA(5) = -LOC2 * (LOC1 - 1.0)
    DNDA(6) =  LOC2 * LOC1
    DNDA(7) = -LOC2 * LOC1
    !
    DNDB(0) = -(LOC2 - 1.0) * (LOC0 - 1.0)
    DNDB(1) =  (LOC2 - 1.0) * LOC0
    DNDB(2) = -(LOC2 - 1.0) * LOC0
    DNDB(3) =  (LOC2 - 1.0) * (LOC0 - 1.0)
    DNDB(4) =  LOC2 * (LOC0 - 1.0)
    DNDB(5) = -LOC2 * LOC0
    DNDB(6) =  LOC2 * LOC0
    DNDB(7) = -LOC2 * (LOC0 - 1.0)
    !
    DNDC(0) = -(LOC1 - 1.0) * (LOC0 - 1.0)
    DNDC(1) =  (LOC1 - 1.0) * LOC0
    DNDC(2) = -LOC1 * LOC0
    DNDC(3) =  LOC1 * (LOC0 - 1.0)
    DNDC(4) =  (LOC1 - 1.0) * (LOC0 - 1.0)
    DNDC(5) = -(LOC1 - 1.0) * LOC0
    DNDC(6) =  LOC1 * LOC0
    DNDC(7) = -LOC1 * (LOC0 - 1.0)
    !
    RETURN
    !
    END SUBROUTINE B8_DERIV_HPAR
    !
END MODULE SHAPE_3D_MOD
