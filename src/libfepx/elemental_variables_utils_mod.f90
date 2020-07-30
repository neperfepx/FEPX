! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE ELEMENTAL_VARIABLES_UTILS_MOD
!
! Module for calculating various elemental values
!
! Contians subroutines
! PLASTICVELGRADSYMSKW: Plastic vel. grad. and derivative variables
! CALC_TOTAL_WORK: Computes elemental total work for a time step.
! CALC_PLASTIC_WORK: Computes elemental plastic work for a time step.
!
USE INTRINSICTYPESMODULE, RK=>REAL_KIND
USE UTILSCRYSTALMODULE
USE DIMSMODULE
USE MICROSTRUCTURE_MOD
USE READ_INPUT_MOD
USE STRESS_STRAIN_MOD
!
IMPLICIT NONE
!
! Public
!
PUBLIC
!
CONTAINS
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
    INTEGER :: ISLIP, I, J, K, IPHASE, NUMIND, N_SLIP
    INTEGER, POINTER  :: INDICES(:)  => NULL()
    REAL(RK), POINTER :: P_HAT_VEC(:, :) => NULL()
    REAL(RK), POINTER :: Q_HAT_VEC(:, :) => NULL()
    REAL(RK), ALLOCATABLE :: DP_HAT_TMP(:, :, :), WP_HAT_TMP(:, :, :)
    REAL(RK) :: TMP(0:(M_EL - 1))
    !
    REAL(RK), PARAMETER :: TWOTHIRDS = 2.0_RK / 3.0_RK
    !
    !---------------------------------------------------------------------------
    !
    ! Initialize
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    DPEFF = 0.0_RK
    DP_HAT = 0.0_RK
    WP_HAT = 0.0_RK
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
        DP_HAT_TMP = 0.0_RK
        WP_HAT_TMP = 0.0_RK
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
    DPEFF = SQRT(TMP)
    !
    WHERE (TMP .LE. VTINY)
        !
        DPEFF(:) = 0.0_RK
        !
    ENDWHERE
    !
    RETURN
    !
    END SUBROUTINE PLASTICVELGRADSYMSKW
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
    EL_WORK_RATE = 0.0_RK
    EL_WORK_STEP = 0.0_RK
    EL_WORK = 0.0_RK
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
        EL_WORK_STEP(I) = DTIME * 0.5_RK * &
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
    EL_WORK_N = 0.0_RK
    EL_WORK_RATE_N = 0.0_RK
    !
    ! All values below VTINY are forced to zero
    !
    WHERE (EL_WORK .LE. VTINY)
        !
        EL_WORK(:) = 0.0_RK
        !
    END WHERE
    !
    WHERE (EL_WORK_RATE .LE. VTINY)
        !
        EL_WORK_RATE(:) = 0.0_RK
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
    EL_WORKP_RATE = 0.0_RK
    EL_WORKP_STEP = 0.0_RK
    EL_WORKP = 0.0_RK
    !
    DO I = EL_SUB1, EL_SUP1
        !
        ! Initialize looped variables
        !
        SIGM = 0.0_RK
        SIGDEV = 0.0_RK
        SIGDEV5 = 0.0_RK
        DP_HAT_SAM = 0.0_RK
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
            & / 3.0_RK
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
        EL_WORKP_STEP(I) = DTIME * 0.5_RK * &
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
    EL_WORKP_N = 0.0_RK
    EL_WORKP_RATE_N = 0.0_RK
    !
    ! All values below VTINY are forced to zero
    !
    WHERE (EL_WORKP .LE. VTINY)
        !
        EL_WORKP(:) = 0.0_RK
        !
    END WHERE
    !
    WHERE (EL_WORKP_RATE .LE. VTINY)
        !
        EL_WORKP_RATE(:) = 0.0_RK
        !
    END WHERE
    !
    EL_WORKP_N = EL_WORKP
    EL_WORKP_RATE_N = EL_WORKP_RATE
    !
    END SUBROUTINE CALC_PLASTIC_WORK
    !
END MODULE ELEMENTAL_VARIABLES_UTILS_MOD
