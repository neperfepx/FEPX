! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE STIFFNESS_EVPS_MOD
!
! Module for the elemental stiffness matrix for the elastic viscoplastic problem
!
! Contains subroutines:
! ELEMENT_STIF_EVPS: Construct elemental stiffness matrix for EVPS problem
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
!
! From libfepx:
!
USE DIMENSIONS_MOD
USE MATERIAL_MATRIX_EVPS_MOD
USE MATRIX_OPERATIONS_MOD
USE MICROSTRUCTURE_MOD
USE QUADRATURE_MOD
USE READ_INPUT_MOD
USE SHAPE_3D_MOD
USE STIFFNESS_VP_MOD
USE UNITS_MOD
!
! From libparallel
!
USE PARALLEL_MOD
!
IMPLICIT NONE
!
CONTAINS
    !
    SUBROUTINE ELEMENT_STIF_EVPS(GSTIFF, GTANSTIFF, F_VEC, GCOORDS, GVEL, &
        & C0_ANGS, C_ANGS, SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, &
        & E_BAR_VEC, E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, JITER_STATE, KEINV, &
        & WTS, EPSEFF, DTIME, INCR, CONVERGED_SOLUTION, AUTO_TIME, NR)
    !
    ! Construct elemental stiffness matrix for EVPS problem
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! GSTIFF: (enters = 0)
    ! GTANSTIFF: (enters = 0)
    ! F_VEC:
    ! GCOORDS:
    ! GVEL:
    ! C0_ANGS:
    ! C_ANGS:
    ! SIG_VEC_N:
    ! SIG_VEC:
    ! CRSS_N:
    ! CRSS:
    ! RSTAR_N:
    ! RSTAR:
    ! E_BAR_VEC:
    ! E_ELAS_KK_BAR:
    ! E_ELAS_KK:
    ! SIG_KK:
    ! JITER_STATE:
    ! KEINV:
    ! WTS:
    ! EPSEFF:
    ! DTIME:
    ! INCR:
    ! CONVERGED_SOLUTION:
    ! AUTO_TIME:
    ! NR:
    !
    REAL(RK), INTENT(INOUT) :: GSTIFF(0:KDIM1, 0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: GTANSTIFF(0:KDIM1, 0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: F_VEC(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: GCOORDS(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: GVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(OUT) :: SIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    REAL(RK), INTENT(OUT) :: SIG_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    REAL(RK), INTENT(IN) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    REAL(RK), INTENT(IN) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(OUT) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK_BAR(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(OUT) :: SIG_KK(EL_SUB1:EL_SUP1, 0:NQPT1)
    INTEGER, INTENT(OUT)  :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: KEINV(0:TVEC1, 1:NUMPHASES)
    REAL(RK), INTENT(IN) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: EPSEFF(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: DTIME
    INTEGER :: INCR
    LOGICAL, INTENT(INOUT) :: CONVERGED_SOLUTION
    INTEGER :: AUTO_TIME
    LOGICAL,  INTENT(IN) :: NR
    !
    ! Locals:
    !
    INTEGER :: IQPT
    INTEGER :: M_EL
    INTEGER :: IPHASE
    INTEGER :: NUMIND
    INTEGER, POINTER :: INDICES(:)=>NULL()
    REAL(RK), POINTER :: E_BAR_VEC_TMP(:, :, :)=>NULL()
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    INTEGER :: L
    INTEGER :: I1
    INTEGER :: I2
    INTEGER :: I3
    INTEGER :: J1
    INTEGER :: J2
    INTEGER :: J3
    INTEGER :: IER
    REAL(RK) :: WT
    REAL(RK) :: DNDX(0:NNPE, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: DNDY(0:NNPE, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: DNDZ(0:NNPE, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: DET(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: DETV(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: LOC0
    REAL(RK) :: LOC1
    REAL(RK) :: LOC2
    REAL(RK) :: S11(EL_SUB1:EL_SUP1)
    REAL(RK) :: S12(EL_SUB1:EL_SUP1)
    REAL(RK) :: S13(EL_SUB1:EL_SUP1)
    REAL(RK) :: S21(EL_SUB1:EL_SUP1)
    REAL(RK) :: S22(EL_SUB1:EL_SUP1)
    REAL(RK) :: S23(EL_SUB1:EL_SUP1)
    REAL(RK) :: S31(EL_SUB1:EL_SUP1)
    REAL(RK) :: S32(EL_SUB1:EL_SUP1)
    REAL(RK) :: S33(EL_SUB1:EL_SUP1)
    REAL(RK) :: C(TVEC, TVEC, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: C_TAN(TVEC, TVEC, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: F(TVEC, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: XNI  (3, 5, EL_SUB1:EL_SUP1)
    REAL(RK) :: XNJ  (5, 3, EL_SUB1:EL_SUP1)
    REAL(RK) :: TEMP1(3, 5, EL_SUB1:EL_SUP1)
    REAL(RK) :: TEMP2(3, 3, EL_SUB1:EL_SUP1)
    REAL(RK) :: TEMP3(3, 5, EL_SUB1:EL_SUP1)
    REAL(RK) :: TEMP4(3, 3, EL_SUB1:EL_SUP1)
    REAL(RK) :: BULK_FAC1(EL_SUB1:EL_SUP1)
    REAL(RK) :: BULK_FAC2(EL_SUB1:EL_SUP1)
    REAL(RK) :: FTEMP(3, EL_SUB1:EL_SUP1)
    INTEGER :: MY_PHASE(EL_SUB1:EL_SUP1)
    !
    !---------------------------------------------------------------------------
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    !
    F_VEC = 0.0D0
    !
    ! PRD: Update internal variables at QPs.
    !
    SIG_VEC_N = GSIG_VEC_N
    E_ELAS_KK_BAR = GELA_KK_BAR(0, :, :)
    !
    DO IQPT = 0, NQPT1
        !
        ! Coordinates in the parent element
        !
        LOC0 = QPLOC(0, IQPT)
        LOC1 = QPLOC(1, IQPT)
        LOC2 = QPLOC(2, IQPT)
        !
        ! Weight
        !
        WT = WTQP(0, IQPT)
        !
        ! Compute elastic strain: SIG_VEC_N --> E_BAR_VEC
        !
        DO IPHASE = 1, NUMPHASES
            !
            CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
            !
            INDICES = INDICES + EL_SUB1
            !
            IF (ASSOCIATED(E_BAR_VEC_TMP)) THEN
                !
                DEALLOCATE(E_BAR_VEC_TMP)
                !
            END IF
            !
            ALLOCATE(E_BAR_VEC_TMP(0:TVEC1, 0:NGRAIN1, 0:(NUMIND - 1)))
            !
            CALL VEC_D_VEC5(KEINV(:, IPHASE), SIG_VEC_N(:, :, INDICES, IQPT), &
                & E_BAR_VEC_TMP, NGRAIN, NUMIND)
            !
            E_BAR_VEC(:, :, INDICES, IQPT)  =E_BAR_VEC_TMP
            !
            DEALLOCATE(E_BAR_VEC_TMP)
            DEALLOCATE(INDICES)
            !
        END DO !NUMPHASES
        !
        ! Compute quadrature quantities given a set of local coordinates
        CALL SFDER_HPAR(LOC0, LOC1, LOC2, GCOORDS, DNDX(:, :, IQPT), &
            & DNDY(:, :, IQPT), DNDZ(:, :, IQPT), DET(:, IQPT), S11, S12, S13, &
            & S21, S22, S23, S31, S32, S33)
        !
    END DO !NQPT
    !
    CALL MATERIAL_MATRIX_EVPS(C, C_TAN, F, DETV, DNDX, DNDY, DNDZ, GVEL, &
        & C0_ANGS, C_ANGS, SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, &
        & KEINV, E_BAR_VEC, E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, JITER_STATE, &
        & WTS, EPSEFF, DTIME, INCR, CONVERGED_SOLUTION, AUTO_TIME, NR)
    !
    IF (.NOT. CONVERGED_SOLUTION .AND. AUTO_TIME .EQ. 1) RETURN
    !
    DO IQPT = 0, NQPT1
        !
        DO IPHASE = 1, NUMPHASES
            !
            WHERE (MY_PHASE == IPHASE)
                !
                BULK_FAC1 = CRYSTAL_PARM(8,IPHASE) / DETV(:, IQPT) * DTIME
                BULK_FAC2 = CRYSTAL_PARM(8,IPHASE) / DETV(:, IQPT) * &
                    & E_ELAS_KK_BAR(:, IQPT)
                !
            END WHERE
            !
        END DO
        !
        IF (MINVAL(DET) .LT. 0.0D0 ) THEN
            !
            DO I = EL_SUB1, EL_SUP1
                !
                IF (DET(I, IQPT) .LT. 0.0D0) THEN
                    !
                    WRITE(DFLT_U, *) 'Error  :       . Element: ',i,', &
                        &determinant: ',DET(I, IQPT)
                    !
                END IF
                !
            END DO
            !
            CALL PAR_QUIT('Error  :       . ELEMENT_STIF_EVPS: Negative &
                &Jacobian(s)', ABORT = .TRUE.)
            !
        END IF
        !
        DET(:, IQPT) = DET(:, IQPT) * WTQP(0, IQPT)
        !
        DO I = 0, NNPE
            !
            I1 = 3 * I
            I2 = I1 + 1
            I3 = I2 + 1
            !
            XNI(1, 1, :) = DNDX(I, :,IQPT)
            XNI(2, 1, :) = -DNDY(I, :,IQPT)
            XNI(3, 1, :) = 0.0D0
            XNI(1, 2, :) = -DNDX(I, :,IQPT) / 3.0D0
            XNI(2, 2, :) = -DNDY(I, :,IQPT) / 3.0D0
            XNI(3, 2, :) = 2.0D0 * DNDZ(I, :,IQPT) / 3.0D0
            XNI(1, 3, :) = 0.5D0 * DNDY(I, :,IQPT)
            XNI(2, 3, :) = 0.5D0 * DNDX(I, :,IQPT)
            XNI(3, 3, :) = 0.0D0
            XNI(1, 4, :) = 0.5D0 * DNDZ(I, :,IQPT)
            XNI(2, 4, :) = 0.0D0
            XNI(3, 4, :) = 0.5D0 * DNDX(I, :,IQPT)
            XNI(1, 5, :) = 0.0D0
            XNI(2, 5, :) = 0.5D0 * DNDZ(I, :,IQPT)
            XNI(3, 5, :) = 0.5D0 * DNDY(I, :,IQPT)
            !
            FTEMP = 0.0D0
            !
            !  RC 6/24/2016: Reordered for better memory striding
            !
            DO L = 1, 5
                !
                DO K = 1, 3
                    !
                    FTEMP(K, :) = FTEMP(K, :) + XNI(K, L, :) * F(L, :, IQPT)
                    !
                END DO
                !
            END DO
            !
            FTEMP(1, :) = FTEMP(1, :) - DNDX(I, :, IQPT) * BULK_FAC2
            FTEMP(2, :) = FTEMP(2, :) - DNDY(I, :, IQPT) * BULK_FAC2
            FTEMP(3, :) = FTEMP(3, :) - DNDZ(I, :, IQPT) * BULK_FAC2
            !
            F_VEC(I1, :) = F_VEC(I1, :) + FTEMP(1, :) * DET(:,IQPT)
            F_VEC(I2, :) = F_VEC(I2, :) + FTEMP(2, :) * DET(:,IQPT)
            F_VEC(I3, :) = F_VEC(I3, :) + FTEMP(3, :) * DET(:,IQPT)
            !
            TEMP1 = 0.0D0
            TEMP3 = 0.0D0
            !
            CALL GEN_MATRIX_MULT(TEMP1, XNI, C(:, :, :, IQPT), IER)
            !
            IF (NR) CALL GEN_MATRIX_MULT(TEMP3, XNI, C_TAN(:,:,:,IQPT), IER)
            !
            DO J = 0, I
                !
                J1 = 3 * J
                J2 = J1 + 1
                J3 = J2 + 1
                !
                XNJ(1, 1, :) = DNDX(J, :, IQPT)
                XNJ(1, 2, :) = -DNDY(J, :, IQPT)
                XNJ(1, 3, :) = 0.0D0
                XNJ(2, 1, :) = -DNDX(J, :, IQPT) / 3.0D0
                XNJ(2, 2, :) = -DNDY(J, :, IQPT) / 3.0D0
                XNJ(2, 3, :) = 2.0D0 * DNDZ(J, :, IQPT)  /3.0D0
                XNJ(3, 1, :) = 0.5D0 * DNDY(J, :, IQPT)
                XNJ(3, 2, :) = 0.5D0 * DNDX(J, :, IQPT)
                XNJ(3, 3, :) = 0.0D0
                XNJ(4, 1, :) = 0.5D0 * DNDZ(J, :, IQPT)
                XNJ(4, 2, :) = 0.0D0
                XNJ(4, 3, :) = 0.5D0 * DNDX(J, :, IQPT)
                XNJ(5, 1, :) = 0.0D0
                XNJ(5, 2, :) = 0.5D0 * DNDZ(J, :, IQPT)
                XNJ(5, 3, :) = 0.5D0 * DNDY(J, :, IQPT)
                !
                TEMP2 = 0.0D0
                TEMP4 = 0.0D0
                !
                CALL GEN_MATRIX_MULT(TEMP2, TEMP1, XNJ, IER)
                !
                IF (NR) CALL GEN_MATRIX_MULT(TEMP4, TEMP3, XNJ, IER)
                !
                ! Assemble secant + volumentric stiffness
                !
                S11 = TEMP2(1, 1, :) * DET(:, IQPT) + DNDX(I, :, IQPT) * &
                    & DNDX(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                S12 = TEMP2(1, 2, :) * DET(:, IQPT) + DNDX(I, :, IQPT) * &
                    & DNDY(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                S13 = TEMP2(1, 3, :) * DET(:, IQPT) + DNDX(I, :, IQPT) * &
                    & DNDZ(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                S21 = TEMP2(2, 1, :) * DET(:, IQPT) + DNDY(I, :, IQPT) * &
                    & DNDX(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                S22 = TEMP2(2, 2, :) * DET(:, IQPT) + DNDY(I, :, IQPT) * &
                    & DNDY(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                S23 = TEMP2(2, 3, :) * DET(:, IQPT) + DNDY(I, :, IQPT) * &
                    & DNDZ(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                S31 = TEMP2(3, 1, :) * DET(:, IQPT) + DNDZ(I, :, IQPT) * &
                    & DNDX(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                S32 = TEMP2(3, 2, :) * DET(:, IQPT) + DNDZ(I, :, IQPT) * &
                    & DNDY(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                S33 = TEMP2(3, 3, :) * DET(:, IQPT) + DNDZ(I, :, IQPT) * &
                    & DNDZ(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
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
                ! Assemble tangent + volumentric stiffness
                !
                IF (NR) THEN
                    !
                    S11 = TEMP4(1, 1, :) * DET(:, IQPT) + DNDX(I, :, IQPT) * &
                        & DNDX(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                    S12 = TEMP4(1, 2, :) * DET(:, IQPT) + DNDX(I, :, IQPT) * &
                        & DNDY(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                    S13 = TEMP4(1, 3, :) * DET(:, IQPT) + DNDX(I, :, IQPT) * &
                        & DNDZ(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                    S21 = TEMP4(2, 1, :) * DET(:, IQPT) + DNDY(I, :, IQPT) * &
                        & DNDX(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                    S22 = TEMP4(2, 2, :) * DET(:, IQPT) + DNDY(I, :, IQPT) * &
                        & DNDY(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                    S23 = TEMP4(2, 3, :) * DET(:, IQPT) + DNDY(I, :, IQPT) * &
                        & DNDZ(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                    S31 = TEMP4(3, 1, :) * DET(:, IQPT) + DNDZ(I, :, IQPT) * &
                        & DNDX(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                    S32 = TEMP4(3, 2, :) * DET(:, IQPT) + DNDZ(I, :, IQPT) * &
                        & DNDY(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                    S33 = TEMP4(3, 3, :) * DET(:, IQPT) + DNDZ(I, :, IQPT) * &
                        & DNDZ(J, :, IQPT) * BULK_FAC1 * DET(:, IQPT)
                    !
                    CALL ADD_TO_STIFFNESS(GTANSTIFF, S11, I1, J1)
                    CALL ADD_TO_STIFFNESS(GTANSTIFF, S22, I2, J2)
                    CALL ADD_TO_STIFFNESS(GTANSTIFF, S33, I3, J3)
                    CALL ADD_TO_STIFFNESS(GTANSTIFF, S12, I1, J2)
                    CALL ADD_TO_STIFFNESS(GTANSTIFF, S13, I1, J3)
                    CALL ADD_TO_STIFFNESS(GTANSTIFF, S21, I2, J1)
                    CALL ADD_TO_STIFFNESS(GTANSTIFF, S23, I2, J3)
                    CALL ADD_TO_STIFFNESS(GTANSTIFF, S31, I3, J1)
                    CALL ADD_TO_STIFFNESS(GTANSTIFF, S32, I3, J2)
                    !
                END IF ! NR
                !
            END DO
            !
        END DO !NNPE
        !
    END DO !NQPT1
    !
    END SUBROUTINE ELEMENT_STIF_EVPS
    !
END MODULE STIFFNESS_EVPS_MOD
