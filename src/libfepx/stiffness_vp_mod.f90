! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE STIFFNESS_VP_MOD
!
! Module for elemental stiffness matrix for the viscoplastic problem
!
! Contains subroutines:
! ELEMENT_STIF_VP: Form elemental stiffness matrix for viscoplastic problem
! ADD_TO_STIFFNESS: Add component to stiffness matrix
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
!
! From libfepx:
!
USE CONVERGENCE_MOD, ONLY: CV_OPTIONS
USE DIMENSIONS_MOD
USE MATERIAL_MATRIX_VP_MOD
USE MATRIX_OPERATIONS_MOD
USE QUADRATURE_MOD
USE READ_INPUT_MOD
USE SHAPE_3D_MOD
!
IMPLICIT NONE
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: ELEMENT_STIF_VP
PUBLIC :: ADD_TO_STIFFNESS
!
CONTAINS
    !
    SUBROUTINE ELEMENT_STIF_VP(ITYPE, GSTIFF, GCOORDS, GVEL, PSCALE, PCNST, &
        & EPSEFF)
    !
    !     Form elemental stiffness matrix for viscoplastic problem.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! ITYPE: 0 = isotropic, 1 = anisotropic
    ! GSTIFF:
    ! GCOORDS:
    ! GVEL:
    ! PSCALE:
    ! PCNST:
    ! EPSEFF:
    !
    INTEGER :: ITYPE
    REAL(RK) :: GSTIFF(0:KDIM1, 0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: GCOORDS(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: GVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: PSCALE(EL_SUB1:EL_SUP1)
    REAL(RK) :: PCNST(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: EPSEFF(EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    INTEGER :: IQPT
    INTEGER :: I
    INTEGER :: J
    INTEGER :: I1
    INTEGER :: I2
    INTEGER :: I3
    INTEGER :: J1
    INTEGER :: J2
    INTEGER :: J3
    INTEGER :: IER
    REAL(RK) :: WT
    REAL(RK) :: DNDX(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDY(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDZ(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: DET(EL_SUB1:EL_SUP1)
    REAL(RK) :: S11(EL_SUB1:EL_SUP1)
    REAL(RK) :: S12(EL_SUB1:EL_SUP1)
    REAL(RK) :: S13(EL_SUB1:EL_SUP1)
    REAL(RK) :: S21(EL_SUB1:EL_SUP1)
    REAL(RK) :: S22(EL_SUB1:EL_SUP1)
    REAL(RK) :: S23(EL_SUB1:EL_SUP1)
    REAL(RK) :: S31(EL_SUB1:EL_SUP1)
    REAL(RK) :: S32(EL_SUB1:EL_SUP1)
    REAL(RK) :: S33(EL_SUB1:EL_SUP1)
    REAL(RK) :: T11(EL_SUB1:EL_SUP1)
    REAL(RK) :: MMTX(EL_SUB1:EL_SUP1)
    REAL(RK) :: SCLFAC(EL_SUB1:EL_SUP1)
    REAL(RK) :: C(TVEC, TVEC, EL_SUB1:EL_SUP1)
    REAL(RK) :: XNI(3, 5, EL_SUB1:EL_SUP1)
    REAL(RK) :: XNJ(5, 3, EL_SUB1:EL_SUP1)
    REAL(RK) :: TEMP1(3, 5, EL_SUB1:EL_SUP1)
    REAL(RK) :: TEMP2(3, 3, EL_SUB1:EL_SUP1)
    REAL(RK) :: LOC0, LOC1, LOC2
    !
    !---------------------------------------------------------------------------
    !
    PCNST = 0.0D0
    MMTX = 0.0D0
    PSCALE = 0.0D0
    !
    LOC0 = 0.25D0
    LOC1 = 0.25D0
    LOC2 = 0.25D0
    !
    CALL SFDER_HPAR(LOC0, LOC1, LOC2, GCOORDS, DNDX, DNDY, DNDZ, DET, S11, &
        & S12, S13, S21, S22, S23, S31, S32, S33)
    !
    CALL MATERIAL_MATRIX_VP(ITYPE, C, DNDX, DNDY, DNDZ, GVEL, SCLFAC, EPSEFF)
    !
    T11 = 0.0D0
    !
    DO IQPT = 0, NQPT1
        !
        LOC0 = QPLOC(0, IQPT)
        LOC1 = QPLOC(1, IQPT)
        LOC2 = QPLOC(2, IQPT)
        !
        WT = WTQP(0, IQPT)
        !
        CALL SFDER_HPAR(LOC0, LOC1, LOC2, GCOORDS, DNDX, DNDY, DNDZ, DET, S11, &
            & S12, S13, S21, S22, S23, S31, S32, S33)
        !
        DET = DET * WT
        PSCALE = PSCALE + SCLFAC * WT
        MMTX = MMTX + DET
        !
        DO I = 0, NNPE
            !
            I1 = 3 * I
            I2 = I1 + 1
            I3 = I2 + 1
            !
            PCNST(I1, :) = PCNST(I1, :) - DNDX(I, :) * DET
            PCNST(I2, :) = PCNST(I2, :) - DNDY(I, :) * DET
            PCNST(I3, :) = PCNST(I3, :) - DNDZ(I, :) * DET
            !
            XNI(1, 1, :) = DNDX(I, :)
            XNI(2, 1, :) = -DNDY(I, :)
            XNI(3, 1, :) = 0.0D0
            XNI(1, 2, :) = -DNDX(I, :) / 3.0D0
            XNI(2, 2, :) = -DNDY(I, :) / 3.0D0
            XNI(3, 2, :) = 2.0D0 * DNDZ(I, :) / 3.0D0
            XNI(1, 3, :) = 0.5D0 * DNDY(I, :)
            XNI(2, 3, :) = 0.5D0 * DNDX(I, :)
            XNI(3, 3, :) = 0.0D0
            XNI(1, 4, :) = 0.5D0 * DNDZ(I, :)
            XNI(2, 4, :) = 0.0D0
            XNI(3, 4, :) = 0.5D0 * DNDX(I, :)
            XNI(1, 5, :) = 0.0D0
            XNI(2, 5, :) = 0.5D0 * DNDZ(I, :)
            XNI(3, 5, :) = 0.5D0 * DNDY(I, :)
            !
            TEMP1 = 0.0D0
            !
            CALL GEN_MATRIX_MULT(TEMP1, XNI, C, IER)
            !
            DO J = 0, I
                !
                J1 = 3 * J
                J2 = J1 + 1
                J3 = J2 + 1
                !
                XNJ(1, 1, :) = DNDX(J, :)
                XNJ(1, 2, :) = -DNDY(J, :)
                XNJ(1, 3, :) = 0.0D0
                XNJ(2, 1, :) = -DNDX(J, :) / 3.0D0
                XNJ(2, 2, :) = -DNDY(J, :) / 3.0D0
                XNJ(2, 3, :) = 2.0D0 * DNDZ(J, :) / 3.0D0
                XNJ(3, 1, :) = 0.5D0 * DNDY(J, :)
                XNJ(3, 2, :) = 0.5D0 * DNDX(J, :)
                XNJ(3, 3, :) = 0.0D0
                XNJ(4, 1, :) = 0.5D0 * DNDZ(J, :)
                XNJ(4, 2, :) = 0.0D0
                XNJ(4, 3, :) = 0.5D0 * DNDX(J, :)
                XNJ(5, 1, :) = 0.0D0
                XNJ(5, 2, :) = 0.5D0 * DNDZ(J, :)
                XNJ(5, 3, :) = 0.5D0 * DNDY(J, :)
                !
                TEMP2 = 0.0D0
                !
                CALL GEN_MATRIX_MULT(TEMP2, TEMP1, XNJ, IER)
                !
                S11 = TEMP2(1, 1, :) * DET
                S12 = TEMP2(1, 2, :) * DET
                S13 = TEMP2(1, 3, :) * DET
                S21 = TEMP2(2, 1, :) * DET
                S22 = TEMP2(2, 2, :) * DET
                S23 = TEMP2(2, 3, :) * DET
                S31 = TEMP2(3, 1, :) * DET
                S32 = TEMP2(3, 2, :) * DET
                S33 = TEMP2(3, 3, :) * DET
                !
                CALL ADD_TO_STIFFNESS(GSTIFF, S11, I1, J1)
                CALL ADD_TO_STIFFNESS(GSTIFF, S22, I2, J2)
                CALL ADD_TO_STIFFNESS(GSTIFF, S33, I3, J3)
                CALL ADD_TO_STIFFNESS(GSTIFF, S12, I1, J2)
                CALL ADD_TO_STIFFNESS(GSTIFF, S13, I1, J3)
                CALL ADD_TO_STIFFNESS(GSTIFF, S21, I2, J1)
                CALL ADD_TO_STIFFNESS(GSTIFF, S23, I2, J3)
                CALL ADD_TO_STIFFNESS(GSTIFF, S31, I3, J1)
                CALL ADD_TO_STIFFNESS(GSTIFF, S32, I3, J2)
                !
            END DO
            !
        END DO !NNPE
        !
    END DO !NQPT
    !
    PSCALE = PSCALE / MMTX
    T11 = CV_OPTIONS%PACC * PSCALE
    !
    DO I = 0, NNPE
        !
        I1 = 3 * I
        I2 = I1 + 1
        I3 = I2 + 1
        !
        DO J = 0, I
            !
            J1 = 3 * J
            J2 = J1 + 1
            J3 = J2 + 1
            !
            S11 = PCNST(I1, :) * PCNST(J1, :) * T11
            S12 = PCNST(I1, :) * PCNST(J2, :) * T11
            S13 = PCNST(I1, :) * PCNST(J3, :) * T11
            S21 = PCNST(I2, :) * PCNST(J1, :) * T11
            S22 = PCNST(I2, :) * PCNST(J2, :) * T11
            S23 = PCNST(I2, :) * PCNST(J3, :) * T11
            S31 = PCNST(I3, :) * PCNST(J1, :) * T11
            S32 = PCNST(I3, :) * PCNST(J2, :) * T11
            S33 = PCNST(I3, :) * PCNST(J3, :) * T11
            !
            CALL ADD_TO_STIFFNESS(GSTIFF, S11, I1, J1)
            CALL ADD_TO_STIFFNESS(GSTIFF, S22, I2, J2)
            CALL ADD_TO_STIFFNESS(GSTIFF, S33, I3, J3)
            CALL ADD_TO_STIFFNESS(GSTIFF, S12, I1, J2)
            CALL ADD_TO_STIFFNESS(GSTIFF, S13, I1, J3)
            CALL ADD_TO_STIFFNESS(GSTIFF, S21, I2, J1)
            CALL ADD_TO_STIFFNESS(GSTIFF, S23, I2, J3)
            CALL ADD_TO_STIFFNESS(GSTIFF, S31, I3, J1)
            CALL ADD_TO_STIFFNESS(GSTIFF, S32, I3, J2)
            !
        END DO
        !
    END DO
    !
    END SUBROUTINE ELEMENT_STIF_VP
    !
    !===========================================================================
    !
    SUBROUTINE ADD_TO_STIFFNESS(GSTIFF, SIJ, I,J)
    !
    ! Add component to stiffness matrix.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! GSTIFF:
    ! SIJ:
    ! I:
    ! J:
    !
    REAL(RK), INTENT(INOUT) :: GSTIFF(0:KDIM1, 0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: SIJ(EL_SUB1:EL_SUP1)
    INTEGER :: I
    INTEGER :: J
    !
    !---------------------------------------------------------------------------
    !
    IF (I .GE. J) THEN
        !
        IF (I .EQ. J) THEN
            !
            GSTIFF(I, J, :) = GSTIFF(I, J, :) + SIJ
            !
        ELSE
            !
            GSTIFF(I, J, :) = GSTIFF(I, J, :) + SIJ
            GSTIFF(J, I, :) = GSTIFF(J, I, :) + SIJ
            !
        END IF
        !
    END IF
    !
    END SUBROUTINE ADD_TO_STIFFNESS
    !
END MODULE STIFFNESS_VP_MOD
