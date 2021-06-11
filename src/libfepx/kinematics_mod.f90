! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE KINEMATICS_MOD
!
! Calculations for kinematic (and adjacent) variables
!
! Contains subroutines:
! CALC_PLASTIC_WORK: Computes elemental plastic work for a time step.
! CALC_TOTAL_WORK: Computes elemental total work for a time step.
! DEFRATE: Compute deformation rate tensor
! DP_WP_HAT: Compute DP_HAT and WP_HAT
! EFF_DEF: Compute effective deformation rate
! FIND_WP_HAT: Compute the plastic spin in the intermediate config, `WP_HAT'
! PLASTICVELGRADSYMSKW: Plastic vel. grad. and derivative variables
! VEL_GRADIENT: Compute velocity gradient
!
! From libf95:
!
USE LIBF95, RK=>REAL_KIND
!
! From libfepx:
!
USE DIMENSIONS_MOD
USE MATRIX_OPERATIONS_MOD
USE MICROSTRUCTURE_MOD
USE READ_INPUT_MOD
!
IMPLICIT NONE
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: CALC_PLASTIC_WORK
PUBLIC :: CALC_TOTAL_WORK
PUBLIC :: DEFRATE
PUBLIC :: DP_WP_HAT
PUBLIC :: EFF_DEF
PUBLIC :: FIND_WP_HAT
PUBLIC :: PLASTICVELGRADSYMSKW
PUBLIC :: VEL_GRADIENT
!
CONTAINS
    !
    SUBROUTINE CALC_PLASTIC_WORK(DTIME, DP_HAT, C_ANGS, S_AVG_3X3, EL_WORKP_N, &
        & EL_WORKP_RATE_N, EL_WORKP)
    !
    ! Calculates elemental plastic work for a given time step.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! DTIME: Time step for current increment
    ! DP_HAT: Elemental plastic deformation rate tensor (5 vector)
    ! C_ANGS: Elemental orientation rotation matrix
    ! S_AVG_3X3: Elemental Cauchy stress tensor (3x3 matrix)
    ! EL_WORKP_N: Previous step's elemental plastic work
    ! EL_WORKP_RATE_N: Previous step's elemental plastic work rate
    ! EL_WORKP: Elemental plastic work
    !
    REAL(RK), INTENT(IN) :: DTIME
    REAL(RK), INTENT(IN) :: DP_HAT(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: S_AVG_3X3(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: EL_WORKP_N(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: EL_WORKP_RATE_N(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT)   :: EL_WORKP(EL_SUB1:EL_SUP1)
    !
    ! Locals:
    ! EL_WORKP_RATE: Elemental plastic work rate
    ! EL_WORKP_STEP: Elemental plastic work for the current step
    ! I: Looping indices
    ! QR5X5: 5x5 rotation matrix
    ! SIGM: Mean stress
    ! SIGDEV: Deviatoric stress
    ! SIGDEV5: Deviatoric stress (5 vector)
    ! DP_HAT_SAM: Pl. def. rate tensor in sample basis
    !
    REAL(RK) :: EL_WORKP_RATE(EL_SUB1:EL_SUP1)
    REAL(RK) :: EL_WORKP_STEP(EL_SUB1:EL_SUP1)
    INTEGER  :: I
    REAL(RK) :: QR5X5(0:TVEC1, 0:TVEC1)
    REAL(RK) :: SIGM
    REAL(RK) :: SIGDEV(0:DIMS1, 0:DIMS1)
    REAL(RK) :: SIGDEV5(0:TVEC1)
    REAL(RK) :: DP_HAT_SAM(0:TVEC1)
    !
    !---------------------------------------------------------------------------
    !
    ! Initialize
    !
    EL_WORKP_RATE = 0.0D0
    EL_WORKP_STEP = 0.0D0
    EL_WORKP = 0.0D0
    !
    DO I = EL_SUB1, EL_SUP1
        !
        ! Initialize looped variables
        !
        SIGM = 0.0D0
        SIGDEV = 0.0D0
        SIGDEV5 = 0.0D0
        DP_HAT_SAM = 0.0D0
        !
        ! Use transpose of crys-to-sample transformation:
        ! Lattice deform (below) transposes input (usually intented to go
        ! sample-to-crys). Transpose will let lattice_deform go crys-to-sample.
        !
        CALL ROT_MAT_SYMM_SER(TRANSPOSE(C_ANGS(:, :, 0, I)), QR5X5)
        !
        ! Calculate plastic work rate (tensor inner product of deviatoric cauchy
        !   stress and plastic deformation rate tensor). Since DP tensor is in 5
        !   vector form, easiest to convert 3x3 stress to deviatoric 6 vector,
        !   then to deviatoric 5 vector
        !
        ! First, construct deviatoric stress tensor (3x3)
        !
        SIGM = (S_AVG_3X3(0, 0, I) + S_AVG_3X3(1, 1, I) + S_AVG_3X3(2, 2, I)) &
            & / 3.0D0
        SIGDEV = S_AVG_3X3(:, :, I)
        SIGDEV(0, 0) = SIGDEV(0, 0) - SIGM
        SIGDEV(1, 1) = SIGDEV(1, 1) - SIGM
        SIGDEV(2, 2) = SIGDEV(2, 2) - SIGM
        !
        ! Next, construct deviatoric 5 vector for stress (SIGDEV5)
        !   Ordering (11-22), 33, 12, 13, 23 (with proper scalings)
        !
        CALL MAT_VEC_SYMM_SER(SIGDEV, SIGDEV5)
        !
        ! Rotate DP_HAT to sample reference frame, find plastic work rate
        !
        CALL LATTICE_DEFORM_SER(QR5X5, DP_HAT(:, I), DP_HAT_SAM)
        EL_WORKP_RATE(I) = (DP_HAT_SAM(0) * SIGDEV5(0)) + &
            & (DP_HAT_SAM(1) * SIGDEV5(1)) + (DP_HAT_SAM(2) * SIGDEV5(2)) + &
            & (DP_HAT_SAM(3) * SIGDEV5(3)) + (DP_HAT_SAM(4) * SIGDEV5(4))
        !
        ! Calculate work over step (trapezoidal time integration)
        !   Use previous work rate EL_WORKP_RATE_N
        !
        EL_WORKP_STEP(I) = DTIME * 0.5D0 * &
            & (EL_WORKP_RATE(I) + (EL_WORKP_RATE_N(I)))
        !
        ! Calculate cumulative work at current step
        !
        EL_WORKP(I) = EL_WORKP_STEP(I) + EL_WORKP_N(I)
        !
    ENDDO
    !
    ! Update previous variables
    !
    EL_WORKP_N = 0.0D0
    EL_WORKP_RATE_N = 0.0D0
    !
    ! All values below VTINY are forced to zero
    !
    WHERE (EL_WORKP .LE. VTINY)
        !
        EL_WORKP(:) = 0.0D0
        !
    END WHERE
    !
    WHERE (EL_WORKP_RATE .LE. VTINY)
        !
        EL_WORKP_RATE(:) = 0.0D0
        !
    END WHERE
    !
    EL_WORKP_N = EL_WORKP
    EL_WORKP_RATE_N = EL_WORKP_RATE
    !
    END SUBROUTINE CALC_PLASTIC_WORK
    !
    !===========================================================================
    !
    SUBROUTINE CALC_TOTAL_WORK(DTIME, D, S_AVG_3X3, EL_WORK_N, &
        & EL_WORK_RATE_N, EL_WORK)
    !
    ! Calculates elemental total work for a given time step.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! DTIME: Time step for current increment
    ! D: Elemental total deformation rate tensor (3x3 matrix)
    ! S_AVG_3X3: Elemental Cauchy stress tensor (3x3 matrix)
    ! EL_WORK_N: Previous step's elemental total work
    ! EL_WORK_RATE_N: Previous step's elemental total work rate
    ! EL_WORK: Elemental total work
    !
    REAL(RK), INTENT(IN) :: DTIME
    REAL(RK), INTENT(IN) :: D(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: S_AVG_3X3(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: EL_WORK_N(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: EL_WORK_RATE_N(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT)   :: EL_WORK(EL_SUB1:EL_SUP1)
    !
    ! Locals:
    ! EL_WORK_RATE: Elemental total work rate
    ! EL_WORK_STEP: Elemental total work for the current step
    ! I: Generic looping index
    !
    REAL(RK) :: EL_WORK_RATE(EL_SUB1:EL_SUP1)
    REAL(RK) :: EL_WORK_STEP(EL_SUB1:EL_SUP1)
    INTEGER  :: I
    !
    !---------------------------------------------------------------------------
    !
    ! Initialize
    !
    EL_WORK_RATE = 0.0D0
    EL_WORK_STEP = 0.0D0
    EL_WORK = 0.0D0
    !
    DO I = EL_SUB1, EL_SUP1
        !
        ! Calculate total work (tensor inner product of cauchy stress and
        !   deformation rate tensor, here both 3x3 matrices)
        !
        EL_WORK_RATE(I) = (D(0, 0, I) * S_AVG_3X3(0, 0, I)) + &
            & (D(0, 1, I) * S_AVG_3X3(0, 1, I)) + &
            & (D(0, 2, I) * S_AVG_3X3(0, 2, I)) + &
            & (D(1, 0, I) * S_AVG_3X3(1, 0, I)) + &
            & (D(1, 1, I) * S_AVG_3X3(1, 1, I)) + &
            & (D(1, 2, I) * S_AVG_3X3(1, 2, I)) + &
            & (D(2, 0, I) * S_AVG_3X3(2, 0, I)) + &
            & (D(2, 1, I) * S_AVG_3X3(2, 1, I)) + &
            & (D(2, 2, I) * S_AVG_3X3(2, 2, I))
        !
        ! Calculate work over step (trapezoidal time integration)
        !   Use previous work rate EL_WORK_RATE_N
        !
        EL_WORK_STEP(I) = DTIME * 0.5D0 * &
            & (EL_WORK_RATE(I) + (EL_WORK_RATE_N(I)))
        !
        ! Calculate cumulative work at current step
        !
        EL_WORK(I) = EL_WORK_STEP(I) + EL_WORK_N(I)
        !
    ENDDO
    !
    ! Update previous variables
    !
    EL_WORK_N = 0.0D0
    EL_WORK_RATE_N = 0.0D0
    !
    ! All values below VTINY are forced to zero
    !
    WHERE (EL_WORK .LE. VTINY)
        !
        EL_WORK(:) = 0.0D0
        !
    END WHERE
    !
    WHERE (EL_WORK_RATE .LE. VTINY)
        !
        EL_WORK_RATE(:) = 0.0D0
        !
    END WHERE
    !
    EL_WORK_N = EL_WORK
    EL_WORK_RATE_N = EL_WORK_RATE
    !
    END SUBROUTINE CALC_TOTAL_WORK
    !
    !===========================================================================
    !
    SUBROUTINE DEFRATE(D, DNDX, DNDY, DNDZ, GVEL)
    !
    ! Calculate deformation rate tensor
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! D:
    ! DNDX:
    ! DNDY:
    ! DNDZ:
    ! GVEL:
    !
    REAL(RK) :: D(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDX(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDY(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDZ(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: GVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: I1
    INTEGER :: I2
    INTEGER :: I3
    REAL(RK) :: DIVV(EL_SUB1:EL_SUP1)
    !
    !---------------------------------------------------------------------------
    !
    D = 0.0D0
    !
    DO I = 0, NNPE
        !
        I1 = 3 * I
        I2 = I1 + 1
        I3 = I2 + 1
        !
        D(0, 0, :) = D(0, 0, :) + DNDX(I, :) * GVEL(I1, :)
        D(1, 1, :) = D(1, 1, :) + DNDY(I, :) * GVEL(I2, :)
        D(2, 2, :) = D(2, 2, :) + DNDZ(I, :) * GVEL(I3, :)
        D(1, 0, :) = D(1, 0, :) + DNDX(I, :) * GVEL(I2, :) &
            & + DNDY(I, :) * GVEL(I1, :)
        D(2, 0, :) = D(2, 0, :) + DNDX(I, :) * GVEL(I3, :)&
            & + DNDZ(I, :) * GVEL(I1, :)
        D(2, 1, :) = D(2, 1, :) + DNDY(I, :) * GVEL(I3, :)&
            & + DNDZ(I, :) * GVEL(I2, :)
        !
    END DO
    !
    D(1, 0, :) = 0.5D0 * D(1, 0, :)
    D(2, 0, :) = 0.5D0 * D(2, 0, :)
    D(2, 1, :) = 0.5D0 * D(2, 1, :)
    !
    D(0, 1, :) = D(1, 0, :)
    D(0, 2, :) = D(2, 0, :)
    D(1, 2, :) = D(2, 1, :)
    !
    DIVV = D(0, 0, :) + D(1, 1, :) + D(2, 2, :)
    DIVV = DIVV / 3.0D0
    !
    D(0, 0, :) = D(0, 0, :) - DIVV
    D(1, 1, :) = D(1, 1, :) - DIVV
    D(2, 2, :) = D(2, 2, :) - DIVV
    !
    END SUBROUTINE DEFRATE
    !
    !===========================================================================
    !
    SUBROUTINE DP_WP_HAT(P_HAT_VEC, DP_HAT, WP_HAT, E_ELAS, E_BAR, W_VEC_LAT, &
        & GDOT, N_SLIP, DT, N, M, NUMIND, INDICES)
    !
    ! Add descriptions here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER, INTENT(IN)   :: N_SLIP, N, M, NUMIND, INDICES(1:NUMIND)
    REAL(RK), INTENT(OUT) :: DP_HAT(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: WP_HAT(0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: P_HAT_VEC(0:TVEC1,0:MAXSLIP1)
    REAL(RK), INTENT(IN)  :: DT
    REAL(RK), INTENT(IN)  :: E_ELAS(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: E_BAR(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: W_VEC_LAT(0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: GDOT(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    INTEGER  :: I, ISLIP
    REAL(RK) :: DP_HAT_TMP(0:TVEC1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: WP_HAT_TMP(0:DIMS1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: E_ELAS_TMP(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: E_BAR_TMP(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: W_VEC_LAT_TMP(0:DIMS1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: GDOT_TMP(0:MAXSLIP1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: P_HAT(0:DIMS1, 0:DIMS1, 0:MAXSLIP1)
    REAL(RK) :: X (0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: EE(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(NUMIND - 1))
    !
    !---------------------------------------------------------------------------
    !
    E_ELAS_TMP=E_ELAS(:, :, :, INDICES)
    E_BAR_TMP=E_BAR(:, :, :, INDICES)
    W_VEC_LAT_TMP=W_VEC_LAT(:, :, INDICES)
    GDOT_TMP=GDOT(:, :, INDICES)
    !
    CALL VEC_MAT_SYMM(P_HAT_VEC, P_HAT, N_SLIP)
    DP_HAT_TMP = 0.0D0
    !
    CALL MAT_X_MAT3(E_ELAS_TMP, E_BAR_TMP, EE, N, NUMIND)
    !
    WP_HAT_TMP(0, :, :) = W_VEC_LAT_TMP(0, :, :) + &
        & 0.5 / DT * (EE(1, 0, :, :) - EE(0, 1, :, :))
    WP_HAT_TMP(1, :, :) = W_VEC_LAT_TMP(1, :, :) + &
        & 0.5 / DT * (EE(2, 0, :, :) - EE(0, 2, :, :))
    WP_HAT_TMP(2, :, :) = W_VEC_LAT_TMP(2, :, :) + &
        & 0.5 / DT * (EE(2, 1, :, :) - EE(1, 2, :, :))
    !
    DO ISLIP = 0, N_SLIP - 1
        !
        CALL MAT_X_MATS3(E_ELAS_TMP, P_HAT(0, 0, ISLIP), X, N, NUMIND)
        !
        WP_HAT_TMP(0, :, :) = WP_HAT_TMP(0, :, :) - &
            & GDOT_TMP(ISLIP, :, :) * (X(1, 0, :, :) - X(0, 1, :, :))
        WP_HAT_TMP(1, :, :) = WP_HAT_TMP(1, :, :) - &
            & GDOT_TMP(ISLIP, :, :) * (X(2, 0, :, :) - X(0, 2, :, :))
        WP_HAT_TMP(2, :, :) = WP_HAT_TMP(2, :, :) - &
            & GDOT_TMP(ISLIP, :, :) * (X(2, 1, :, :) - X(1, 2, :, :))
        !
        DO I = 0, TVEC1
            !
            DP_HAT_TMP(I, :, :) = DP_HAT_TMP(I, :, :) + &
                & GDOT_TMP(ISLIP, :, :) * P_HAT_VEC(I, ISLIP)
            !
        ENDDO
        !
    ENDDO
    !
    DP_HAT(:,:,INDICES) = DP_HAT_TMP
    WP_HAT(:, :, INDICES) = WP_HAT_TMP
    !
    RETURN
    !
    END SUBROUTINE DP_WP_HAT
    !
    !===========================================================================
    !
    SUBROUTINE EFF_DEF(DEFF, D, M)
    !
    ! Compute the effective deformation rate from the deformation rate tensor.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! DEFF:
    ! D:
    ! M:
    REAL(RK), INTENT(OUT) :: DEFF(0:(M - 1))
    REAL(RK), INTENT(IN) :: D(0:DIMS1, 0:DIMS1, 0:(M - 1))
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    !
    REAL(RK), PARAMETER :: TWOTHIRDS = 2.0D0 / 3.0D0
    INTEGER :: I
    !
    !---------------------------------------------------------------------------
    !
    ! Effective deformation rate
    DO I = 0, M - 1
        !
        DEFF(I) = DSQRT(TWOTHIRDS * (D(0, 0, I) * D(0, 0, I) + D(1, 1, I) * &
            & D(1, 1, I) + D(2, 2, I) * D(2, 2, I) + 2.0D0 * ( D(0, 1, I) * &
            & D(0, 1, I) + D(0, 2, I) * D(0, 2, I) + D(1, 2, I) * D(1, 2, I))))
        !
    END DO
    !
    END SUBROUTINE EFF_DEF
    !
    !===========================================================================
    !
    SUBROUTINE FIND_WP_HAT(WP_HAT, E_ELAS, E_BAR, W_VEC_GRN, GDOT, QR5X5, &
        & DT, N, M)
    !
    ! Compute the plastic spin in the intermediate config, `WP_HAT'.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    REAL(RK), INTENT(OUT) :: WP_HAT(0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: E_ELAS(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: E_BAR(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: W_VEC_GRN(0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: GDOT(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: QR5X5(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: DT
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    !
    INTEGER :: I, ISLIP, IPHASE, NUMIND, N_SLIP
    INTEGER, POINTER :: INDICES(:) => NULL()
    !
    REAL(RK), POINTER :: P_HAT_VEC(:,:) => NULL()
    REAL(RK) :: EE(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: DP_HAT(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: TEMP(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: DP_HAT_TENS(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: X(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    INTEGER  :: MY_PHASE(0:(M-1))
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    !
    CALL MAT_X_MAT3(E_ELAS, E_BAR, EE, N, M)
    !
    WP_HAT(0, :, :) = W_VEC_GRN(0, :, :) + 0.5 / &
        & DT * (EE(1, 0, :, :) - EE(0, 1, :, :))
    WP_HAT(1, :, :) = W_VEC_GRN(1, :, :) + 0.5 / &
        & DT * (EE(2, 0, :, :) - EE(0, 2, :, :))
    WP_HAT(2, :, :) = W_VEC_GRN(2, :, :) + 0.5 / &
        & DT * (EE(2, 1, :, :) - EE(1, 2, :, :))
    !
    DP_HAT = 0.0D0
    !
    DO IPHASE = 1, NUMPHASES
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE), DEV=P_HAT_VEC)
        N_SLIP = CTYPE(IPHASE)%NUMSLIP
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
        !
        DO ISLIP = 0, (N_SLIP - 1)
            !
            DO I = 0, TVEC1
                !
                DP_HAT(I, :, INDICES) = DP_HAT(I, :, INDICES) + &
                    & GDOT(ISLIP, :, INDICES) * &
                    & P_HAT_VEC(I + 1, ISLIP + 1)
                !
            ENDDO
            !
        ENDDO
        !
        DEALLOCATE(P_HAT_VEC)
        DEALLOCATE(INDICES)
        !
    ENDDO
    !
    CALL MAT_X_VEC5(QR5X5, DP_HAT, TEMP, N, M)
    !
    CALL VEC_MAT_SYMM_GRN(TEMP, DP_HAT_TENS, N, M)
    !
    CALL MAT_X_MAT3(E_ELAS, DP_HAT_TENS, X, N, M)
    !
    WP_HAT(0, :, :) = WP_HAT(0, :, :) - X(1, 0, :, :) + X(0, 1, :, :)
    WP_HAT(1, :, :) = WP_HAT(1, :, :) - X(2, 0, :, :) + X(0, 2, :, :)
    WP_HAT(2, :, :) = WP_HAT(2, :, :) - X(2, 1, :, :) + X(1, 2, :, :)
    !
    RETURN
    !
    END SUBROUTINE FIND_WP_HAT
    !
    !===========================================================================
    !
    SUBROUTINE PLASTICVELGRADSYMSKW(DP_HAT,WP_HAT,DPEFF,GDOT,M_EL)
    !
    ! Calculate plastic velocity gradient and derivative variables.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! DP_HAT: Plastic deformation rate tensor
    ! WP_HAT: Plastic spin rate tensor
    ! DPEFF: Effective (equivalent) plastic deformation rate
    ! GDOT: Slip system shear rates
    ! M_EL: Number of elements
    !
    REAL(RK), INTENT(OUT) ::  DP_HAT(0:TVEC1,0:(M_EL - 1))
    REAL(RK), INTENT(OUT) ::  WP_HAT(0:DIMS1,0:(M_EL - 1))
    REAL(RK), INTENT(OUT) ::  DPEFF(0:(M_EL - 1))
    REAL(RK), INTENT(IN)  ::  GDOT(0:MAXSLIP1, 0:(M_EL - 1))
    INTEGER,  INTENT(IN)  ::  M_EL
    !
    ! Locals:
    !
    INTEGER :: MY_PHASE(0:(M_EL - 1))
    INTEGER :: ISLIP, I, IPHASE, NUMIND, N_SLIP
    INTEGER, POINTER  :: INDICES(:)  => NULL()
    REAL(RK), POINTER :: P_HAT_VEC(:, :) => NULL()
    REAL(RK), POINTER :: Q_HAT_VEC(:, :) => NULL()
    REAL(RK), ALLOCATABLE :: DP_HAT_TMP(:, :, :), WP_HAT_TMP(:, :, :)
    REAL(RK) :: TMP(0:(M_EL - 1))
    !
    REAL(RK), PARAMETER :: TWOTHIRDS = 2.0D0 / 3.0D0
    !
    !---------------------------------------------------------------------------
    !
    ! Initialize
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    DPEFF = 0.0D0
    DP_HAT = 0.0D0
    WP_HAT = 0.0D0
    !
    DO IPHASE = 1, NUMPHASES
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE), DEV = P_HAT_VEC, SKW = Q_HAT_VEC)
        !
        N_SLIP = CTYPE(IPHASE)%NUMSLIP
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
        !
        ALLOCATE(DP_HAT_TMP(0:TVEC1,0:(NUMIND-1),0:(N_SLIP - 1)))
        ALLOCATE(WP_HAT_TMP(0:DIMS1,0:(NUMIND-1),0:(N_SLIP - 1)))
        DP_HAT_TMP = 0.0D0
        WP_HAT_TMP = 0.0D0
        !
        DO ISLIP = 0, (N_SLIP - 1)
            !
            DO I = 0, TVEC1
                !
                DP_HAT_TMP(I, :, ISLIP) = DP_HAT_TMP(I, :, ISLIP) + &
                    & GDOT(ISLIP, INDICES) * P_HAT_VEC(I + 1, ISLIP + 1)
                DP_HAT(I,INDICES) = DP_HAT(I, INDICES) + DP_HAT_TMP(I, :, ISLIP)
                !
            END DO
            !
            DO I = 0, DIMS1
                !
                WP_HAT_TMP(I, :, ISLIP) = WP_HAT_TMP(I, :, ISLIP) + &
                    & GDOT(ISLIP, INDICES) * Q_HAT_VEC(I + 1, ISLIP + 1)
                WP_HAT(I, INDICES) = WP_HAT(I, INDICES) + &
                    & WP_HAT_TMP(I, :, ISLIP)
                !
            END DO
            !
        END DO
        !
        DEALLOCATE(P_HAT_VEC, Q_HAT_VEC)
        DEALLOCATE(INDICES)
        DEALLOCATE(DP_HAT_TMP, WP_HAT_TMP)
        !
    END DO
    !
    TMP = TWOTHIRDS * SUM(DP_HAT * DP_HAT, 1)
    DPEFF = DSQRT(TMP)
    !
    WHERE (TMP .LE. VTINY)
        !
        DPEFF(:) = 0.0D0
        !
    ENDWHERE
    !
    RETURN
    !
    END SUBROUTINE PLASTICVELGRADSYMSKW
    !
    !===========================================================================
    !
    SUBROUTINE VEL_GRADIENT(VGRAD, DNDX, DNDY, DNDZ, GVEL)
    !
    ! Compute velocity gradient.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! VGRAD:
    ! DNDX:
    ! DNDY:
    ! DNDZ:
    ! GVEL:
    !
    REAL(RK), INTENT(OUT) :: VGRAD(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: DNDX(0:NNPE,  EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: DNDY(0:NNPE,  EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: DNDZ(0:NNPE,  EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: GVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: I1
    INTEGER :: I2
    INTEGER :: I3
    !
    !---------------------------------------------------------------------------
    !
    VGRAD = 0.0D0
    !
    DO I = 0, NNPE
        !
        I1 = 3 * I
        I2 = I1 + 1
        I3 = I2 + 1
        !
        VGRAD(0, 0, :) = VGRAD(0, 0, :) + DNDX(I, :) * GVEL(I1, :)
        VGRAD(0, 1, :) = VGRAD(0, 1, :) + DNDY(I, :) * GVEL(I1, :)
        VGRAD(0, 2, :) = VGRAD(0, 2, :) + DNDZ(I, :) * GVEL(I1, :)
        VGRAD(1, 0, :) = VGRAD(1, 0, :) + DNDX(I, :) * GVEL(I2, :)
        VGRAD(1, 1, :) = VGRAD(1, 1, :) + DNDY(I, :) * GVEL(I2, :)
        VGRAD(1, 2, :) = VGRAD(1, 2, :) + DNDZ(I, :) * GVEL(I2, :)
        VGRAD(2, 0, :) = VGRAD(2, 0, :) + DNDX(I, :) * GVEL(I3, :)
        VGRAD(2, 1, :) = VGRAD(2, 1, :) + DNDY(I, :) * GVEL(I3, :)
        VGRAD(2, 2, :) = VGRAD(2, 2, :) + DNDZ(I, :) * GVEL(I3, :)
        !
    END DO
    !
    END SUBROUTINE VEL_GRADIENT
    !
END MODULE KINEMATICS_MOD
