! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE POLYCRYSTAL_RESPONSE_EVPS_MOD
!
! Module hangling elastic-viscoplastic response for polycrystals.
!
! Contains subroutines:
! POLYCRYSTAL_RESPONSE_EVPS: EVPS response for polycrystals
! SOLVE_STATE_DEV_EVPS: Solves for the deviatoric state
! SOLVE_STATE_VOL_EVPS: Solves for the volumetric state
! POLYCRYSTAL_RESPONSE_EVPS_QP: EVPS response for polycrystal at quad points
! SOLVE_STATE_DEV_EVPS_QP: Solves for the deviatoric state at quad points
! SOLVE_STATE_VOL_EVPS_QP: Solves for the volumetric state at quad points
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
!
! From libfepx:
!
USE CONVERGENCE_MOD, ONLY: CV_OPTIONS
USE DIMENSIONS_MOD
USE MATRIX_OPERATIONS_MOD
USE MICROSTRUCTURE_MOD
USE READ_INPUT_MOD
USE RSTARN_SOLVE_MOD
USE STRESS_SOLVE_EVPS_MOD, ONLY: STRESS_SOLVE_EVPS
USE STRESS_SOLVE_VP_MOD
USE UNITS_MOD
!
IMPLICIT NONE
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: POLYCRYSTAL_RESPONSE_EVPS
PUBLIC :: POLYCRYSTAL_RESPONSE_EVPS_QP
!
CONTAINS
    !
    SUBROUTINE POLYCRYSTAL_RESPONSE_EVPS(D_VEC, W_VEC, C0_ANGS, C_ANGS, &
        & SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, E_BAR_VEC, EPSEFF, &
        & D_KK, SIG_KK, E_ELAS_KK_BAR, E_ELAS_KK, JITER_STATE, KEINV, INCR, &
        & DTIME, CONVERGED_SOLUTION, AUTO_TIME)
    !
    ! EVPS response for polycrystals
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! D_VEC:
    ! W_VEC:
    ! C0_ANGS:
    ! C_ANGS:
    ! SIG_VEC_N:
    ! SIG_VEC:
    ! CRSS_N:
    ! CRSS:
    ! RSTAR_N:
    ! RSTAR:
    ! E_BAR_VEC:
    ! EPSEFF:
    ! D_KK:
    ! SIG_KK:
    ! E_ELAS_KK_BAR:
    ! E_ELAS_KK:
    ! JITER_STATE:
    ! KEINV:
    ! INCR:
    ! DTIME:
    ! CONVERGED_SOLUTION:
    ! AUTO_TIME:
    !
    REAL(RK), INTENT(IN) :: D_VEC(0:TVEC1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(IN) :: W_VEC(0:DIMS1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(IN) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: C_ANGS (0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(IN) :: SIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
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
    REAL(RK), INTENT(IN) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    REAL(RK), INTENT(IN) :: EPSEFF(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(IN) :: D_KK(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(OUT) :: SIG_KK(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(IN) :: E_ELAS_KK_BAR(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK(EL_SUB1:EL_SUP1, 0:NQPT1)
    INTEGER, INTENT(OUT) :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: KEINV(0:TVEC1, 1:NUMPHASES)
    INTEGER :: INCR
    REAL(RK) :: DTIME
    LOGICAL, INTENT(INOUT) :: CONVERGED_SOLUTION
    INTEGER :: AUTO_TIME
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: M_EL
    REAL(RK) :: D_VEC_GRN(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: W_VEC_GRN(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: EPSEFF_LAT(0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    !
    !---------------------------------------------------------------------------
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    ! This part should be cleaned when NGRAIN1 stuff is removed
    ! Spread {d}_sm & {w}_sm to all grains in aggregate (Taylor Assumption)
    !
    DO J = 0, NQPT1
        !
        DO I = 0, TVEC1
            !
            D_VEC_GRN(I, :, :, J) = SPREAD(D_VEC(I, :, J), DIM = 1, &
                & NCOPIES = NGRAIN)
            !
        END DO
        !
        DO I = 0, DIMS1
            !
            W_VEC_GRN(I, :, :, J) = SPREAD(W_VEC(I, :, J), DIM = 1, &
                & NCOPIES = NGRAIN)
            !
        END DO
        !
        ! Spread over grains: EPSEFF --> EPSEFF_LAT
        !
        EPSEFF_LAT(:, :, J) = SPREAD(EPSEFF(:, J), DIM = 1, NCOPIES = NGRAIN)
        !
    END DO
    !
    ! Solve for state
    !
    ! Deviatoric
    !
    CALL SOLVE_STATE_DEV_EVPS(D_VEC_GRN, W_VEC_GRN,  C0_ANGS, C_ANGS, &
        & SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, E_BAR_VEC, KEINV, &
        & EPSEFF_LAT, JITER_STATE, INCR, DTIME, CONVERGED_SOLUTION, AUTO_TIME)
    !
    IF (.NOT. CONVERGED_SOLUTION .AND. AUTO_TIME .EQ. 1) RETURN
    !
    ! Volumetric
    !
    DO J = 0, NQPT1
        !
        CALL SOLVE_STATE_VOL_EVPS(E_ELAS_KK_BAR(:, J), E_ELAS_KK(:, J), &
            & D_KK(:, J), SIG_KK(:, J), DTIME, M_EL)
        !
    END DO
    !
    END SUBROUTINE POLYCRYSTAL_RESPONSE_EVPS
    !
    !===========================================================================
    !
    SUBROUTINE SOLVE_STATE_DEV_EVPS(D_VEC, W_VEC, C0_ANGS, C_ANGS, SIG_LAT_N, &
        & SIG_LAT, CRSS_N, CRSS, RSTAR_N, RSTAR, E_BAR_VEC, KEINV, EPSEFF, &
        & JITER_STATE, INCR, DTIME, CONVERGED_SOLUTION, AUTO_TIME)
    !
    ! Solves for the deviatoric state
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! D_VEC:
    ! W_VEC:
    ! C0_ANGS:
    ! C_ANGS:
    ! SIG_LAT_N:
    ! SIG_LAT:
    ! CRSS_N:
    ! CRSS:
    ! RSTAR_N:
    ! RSTAR:
    ! E_BAR_VEC:
    ! KEINV:
    ! EPSEFF:
    ! JITER_STATE:
    ! INCR:
    ! DTIME:
    ! CONVERGED_SOLUTION:
    ! AUTO_TIME:
    !
    REAL(RK), INTENT(IN) :: D_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(IN) :: W_VEC(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(IN) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: C_ANGS (0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(IN) :: SIG_LAT_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    REAL(RK), INTENT(OUT) :: SIG_LAT(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    REAL(RK), INTENT(IN) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    REAL(RK), INTENT(IN) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(IN) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    REAL(RK), INTENT(IN) :: KEINV(0:TVEC1, 1:NUMPHASES)
    REAL(RK), INTENT(IN) :: EPSEFF(0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    INTEGER, INTENT(OUT) :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER :: INCR
    REAL(RK) :: DTIME
    LOGICAL, INTENT(INOUT) :: CONVERGED_SOLUTION
    INTEGER :: AUTO_TIME
    !
    ! Locals:
    !
    LOGICAL, PARAMETER :: VP_LOG = .FALSE.
    LOGICAL :: DONE(0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    LOGICAL :: CONVERGED_NEWTON(0:NQPT1)
    LOGICAL :: CONVERGED_STATE(0:NQPT1)
    INTEGER :: ITER_STATE
    INTEGER :: M_EL
    INTEGER :: ISLIP
    INTEGER :: N_SLIP
    REAL(RK) :: QR5X5(0:TVEC1, 0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: QR3X3(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: W_VEC_LAT(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: WP_HAT(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: D_VEC_LAT(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: E_BAR_VEC_R(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: CRSS_0(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: NORM_S_0(0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: NORM_S(0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: DIFF_NORM_S(0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: DIFF_CRSS  (0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: D_RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    INTEGER :: IPHASE
    INTEGER :: I
    INTEGER :: J
    INTEGER :: MY_PHASE(0:(EL_SUP1-EL_SUB1))
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    CONVERGED_NEWTON = .TRUE.
    CONVERGED_STATE  = .TRUE.
    !
    JITER_STATE  = 0
    !
    CRSS = SPREAD(CRSS_N, 4, NQPT)
    !
    DONE = .FALSE.
    !
    ! Estimate for the stresses :
    ! For INCR = 1: based on viscoplastic solution
    ! For INCR > 1: based on previous solution (unrotated)
    !
    DO I = 0, NQPT1
        !
        IF (INCR .EQ. 1) THEN
            !
            ! C_ANGS [3x3] --> QR5X5 [5x5]
            !
            CALL ROT_MAT_SYMM(C_ANGS(:, :, : ,: ,I), QR5X5, NGRAIN, M_EL)
            !
            ! D_VEC(sample coords) --> D_VEC_LAT(crystal coords)
            ! {D_VEC_LAT} = [QR5X5]'{D_VEC}
            !
            CALL LATTICE_DEFORM(QR5X5, D_VEC(:, :, :, I), &
                & D_VEC_LAT(:, :, :, I), NGRAIN, M_EL)
            !
            ! Estimate of SIG_LAT (visco-plastic solution)
            !
            CALL STRESS_SOLVE_VP(SIG_LAT(:, :, :, I), D_VEC_LAT(:, :, :, I), &
                & CRSS(:, :, :, I), EPSEFF(:, :, I), VP_LOG)
            !
            ! Use a fraction of the viscoplastic solution as initial guess
            !
            SIG_LAT(:, :, :, I) = 0.5 * SIG_LAT(:, :, :, I)
            !
        ELSE
            !
            ! Use value from previous increment
            !
            SIG_LAT(:, :, :, I) = SIG_LAT_N(:, :, :, I)
            !
        END IF
        !
        ! First (Forward Euler) estimate for CRSS and RSTAR (c).
        !
        ! C_ANGS [3x3] --> QR3X3 [3x3]
        !
        CALL ROT_MAT_SKEW(C_ANGS(:, :, :, :, I), QR3X3, NGRAIN, M_EL)
        !
        ! W_VEC(sample coords) --> W_VEC_LAT(crystal coords)
        !
        CALL LATTICE_SPIN(QR3X3, W_VEC(:, :, :, I), W_VEC_LAT(:, :, :, I), &
            & NGRAIN, M_EL)
        !
        ! Estimate CRSS (g), RSTAR (R*) and D_RSTAR (dR*)
        ! Calculate also C_ANGS: [c]=[c_0]*[R*]
        !
        ACCUMSHEAR = GACCUMSHEAR(:, :, :, I)
        !
        CALL RSTARN_SOLVE(CRSS_N, CRSS(:, :, :, I), RSTAR_N, &
            & RSTAR(:, :, :, :, I), C0_ANGS,  C_ANGS(:, :, :, :, I), &
            & SIG_LAT(:, :, :, I), W_VEC_LAT(:, :, :, I), DTIME, &
            & EPSEFF(:, :, I), DONE(:, :, I), D_RSTAR(:, :, :, :, I), &
            & UPD_EULER_FWD)
        !
        ! Compute 2-norm for array of 5-vectors
        !
        CALL NORM_VEC(NORM_S_0(:, :, I), SIG_LAT(:, :, :, I), NGRAIN, M_EL)
        !
        ! Update CRSS_0
        !
        CRSS_0(:, :, :, I) = CRSS(:, :, :, I)
        !
    END DO !NQPT
    !
    ! Iterate for the material state
    !
    DO I = 0, NQPT1
        !
        ITER_STATE = 1
        ACCUMSHEAR = GACCUMSHEAR(:, :, :, I)
        !
        DO WHILE ((ANY(.NOT. DONE(:, :, I))) .AND. (ITER_STATE .LE. &
            & CV_OPTIONS% SX_MAX_ITERS_STATE))
            !
            ! C_ANGS [3x3] --> QR5X5 [5x5]
            !
            CALL ROT_MAT_SYMM(C_ANGS(:,:, :, :, I), QR5X5, NGRAIN, M_EL)
            !
            ! D_VEC(sample coords) --> D_VEC_LAT(crystal coords)
            !
            CALL LATTICE_DEFORM(QR5X5, D_VEC(:, :, :, I), &
                & D_VEC_LAT(:, :, :, I), NGRAIN, M_EL)
            !
            ! D_RSTAR [3x3] --> QR5X5 [5x5]
            !
            CALL ROT_MAT_SYMM(D_RSTAR(:,:, :, :, I), QR5X5, NGRAIN, M_EL)
            !
            ! Apply dR* to E_BAR_VEC --> E_BAR_VEC_R
            !
            CALL LATTICE_DEFORM(QR5X5, E_BAR_VEC(:, :, :, I), E_BAR_VEC_R, &
                & NGRAIN, M_EL)
            !
            ! --> SIG_LAT
            !
            CALL STRESS_SOLVE_EVPS(SIG_LAT(:, :, :, I), D_VEC_LAT(:, :, :, I), &
                & W_VEC_LAT(:, :, :, I), E_BAR_VEC_R, CRSS(:, :, :, I), KEINV, &
                & DTIME, WP_HAT(:, :, :, I), ITER_STATE, DONE(:, :, I), &
                & CONVERGED_NEWTON(I))
            !
            IF (.NOT. CONVERGED_NEWTON(I) .AND. AUTO_TIME .EQ. 1) GO TO 100
            !
            ! C_ANGS [3x3] --> QR3X3 [3x3]
            !
            CALL ROT_MAT_SKEW(C_ANGS(:,:, :, :, I), QR3X3, NGRAIN, M_EL)
            !
            ! W_VEC(sample coords) --> W_VEC_LAT(crystal coords)
            !
            CALL LATTICE_SPIN(QR3X3, W_VEC(:, :, :, I), W_VEC_LAT(:, :, :, I), &
                & NGRAIN, M_EL)
            !
            ! Calculate CRSS (g), RSTAR (R*) and D_RSTAR (dR*)
            ! Calculate also C_ANGS: [c]=[c_0]*[R*]
            !
            CALL RSTARN_SOLVE(CRSS_N, CRSS(:, :, :, I), RSTAR_N, &
                & RSTAR(:,:, :, :, I), C0_ANGS, C_ANGS(:,:, :, :, I), &
                & SIG_LAT(:, :, :, I), WP_HAT(:, :, :, I), DTIME, &
                & EPSEFF(:, :, I), DONE(:, :, I), D_RSTAR(:,:, :, :, I), &
                & UPD_EULER_BWD)
            !
            CALL NORM_VEC(NORM_S(:, :, I), SIG_LAT(:, :, :, I), NGRAIN, M_EL)
            !
            ! DEB: This section was originally dones with nested `WHERE'
            !   constructs, but was changed because the AIX compiler rejected
            !   them, although the CM compiler had no problem.
            !
            DO ISLIP = 0,MAXSLIP1
                !
                WHERE (.NOT. DONE(:, :, I))
                    !
                    DIFF_CRSS(ISLIP,:, :, I) = DABS(CRSS(ISLIP,:, :, I) - &
                        & CRSS_0(ISLIP,:, :, I))
                    !
                END WHERE
                !
            END DO
            !
            WHERE (.NOT. DONE(:, :, I))
                !
                DIFF_NORM_S(:, :, I) = DABS(NORM_S(:, :, I) - NORM_S_0(:, :, I))
                !
            END WHERE
            !
            !
            DO J = 0, NGRAIN - 1
                !
                DO IPHASE = 1, NUMPHASES
                    !
                    CALL CRYSTALTYPEGET(CTYPE(IPHASE))
                    !
                    N_SLIP = CTYPE(IPHASE)%NUMSLIP
                    !
                    DO ISLIP = 0, N_SLIP - 1
                        !
                        ! Currently have an initial CRSS_0 scaling factor from
                        !   the simple latent hardening model but might end up
                        !   getting rid of it later versions
                        !
                        WHERE((.NOT. DONE(J, :, I))  .AND. &
                            & (MY_PHASE .EQ. IPHASE) .AND. &
                            & (DIFF_NORM_S(J, :, I) .LT. (TOLER_STATE * &
                                & CRYSTAL_PARM(3, IPHASE))) .AND. &
                            & (DIFF_CRSS(ISLIP, J, :, I) .LT. (TOLER_STATE * &
                                & CRYSTAL_PARM(3, IPHASE))))
                            !
                            DONE(J, :, I) = .TRUE.
                            JITER_STATE(J, :) = ITER_STATE
                            !
                        END WHERE
                        !
                    END DO !N_SLIP
                    !
                END DO !NUMPHASES
                !
            END DO !NGRAIN
            !
            DO ISLIP = 0, MAXSLIP1
                !
                WHERE (.NOT. DONE(:, :, I))
                    !
                    NORM_S_0(:, :, I) = NORM_S(:, :, I)
                    CRSS_0(ISLIP,:, :, I) = CRSS(ISLIP,:, :, I)
                    !
                END WHERE
                !
            END DO
            !
            ITER_STATE = ITER_STATE + 1
            !
        END DO ! DO WHILE
        !
    END DO !NQPT
    !
    DO I = 0,NQPT1
        !
        IF (ANY(.NOT. DONE(:, :, I))) THEN
            !
            CONVERGED_STATE(I) = .FALSE.
            WRITE(DFLT_U, '(A)') 'Warning:       . Not all crystals converged.'
            !WRITE(DFLT_U, *) 'Warning:       . Crystals = ', NGRAIN * M_EL, ', &
            !    &converged = ', count(DONE(:, :, I))
            !
        END IF
        !
    END DO
    !
    100 IF (ANY(.NOT. CONVERGED_NEWTON) .OR. ANY(.NOT. CONVERGED_STATE)) &
        & CONVERGED_SOLUTION = .FALSE.
    !
    END SUBROUTINE SOLVE_STATE_DEV_EVPS
    !
    !===========================================================================
    !
    SUBROUTINE SOLVE_STATE_VOL_EVPS(E_ELAS_KK_BAR, E_ELAS_KK, D_KK, SIG_KK, &
        & DTIME, M)
    !
    ! Solves for the volumetric state
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! E_ELAS_KK_BAR:
    ! E_ELAS_KK:
    ! D_KK:
    ! SIG_KK:
    ! DTIME:
    ! M:
    !
    REAL(RK), INTENT(IN) :: E_ELAS_KK_BAR(0:(M - 1))
    REAL(RK), INTENT(OUT) :: E_ELAS_KK(0:(M - 1))
    REAL(RK), INTENT(IN) :: D_KK(0:(M - 1))
    REAL(RK), INTENT(OUT) :: SIG_KK(0:(M - 1))
    REAL(RK), INTENT(IN) :: DTIME
    INTEGER :: M
    !
    ! Locals:
    !
    INTEGER :: MY_PHASE(0:M - 1)
    INTEGER :: IPHASE
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    !
    ! Compute volumetric response in sample reference frame
    !
    E_ELAS_KK = E_ELAS_KK_BAR + DTIME * D_KK
    !
    DO IPHASE = 1, NUMPHASES
        !
        WHERE (MY_PHASE .EQ. IPHASE)
            !
            SIG_KK = 3.0D0 * CRYSTAL_PARM(8,IPHASE) * E_ELAS_KK
            !
        END WHERE
        !
    END DO !NUMPHASES
    !
    END SUBROUTINE SOLVE_STATE_VOL_EVPS
    !
    !===========================================================================
    !
    SUBROUTINE POLYCRYSTAL_RESPONSE_EVPS_QP(D_VEC, W_VEC, C0_ANGS, C_ANGS, &
        & SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, E_BAR_VEC, &
        & EPSEFF, D_KK, SIG_KK, E_ELAS_KK_BAR, E_ELAS_KK, JITER_STATE, KEINV, &
        & INCR, DTIME, CONVERGED_SOLUTION, AUTO_TIME)
    !
    ! EVPS response for polycrystal at quad points
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! D_VEC:
    ! W_VEC:
    ! C0_ANGS:
    ! C_ANGS:
    ! SIG_VEC_N:
    ! SIG_VEC:
    ! CRSS_N:
    ! CRSS:
    ! RSTAR_N:
    ! RSTAR:
    ! E_BAR_VEC:
    ! EPSEFF:
    ! D_KK:
    ! SIG_KK:
    ! E_ELAS_KK_BAR:
    ! E_ELAS_KK:
    ! JITER_STATE:
    ! KEINV:
    ! INCR:
    ! DTIME:
    ! CONVERGED_SOLUTION:
    ! AUTO_TIME:
    !
    REAL(RK), INTENT(IN) :: D_VEC(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: W_VEC(0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: SIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: SIG_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: EPSEFF(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: D_KK(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: SIG_KK(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: E_ELAS_KK_BAR(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK(EL_SUB1:EL_SUP1)
    INTEGER, INTENT(OUT) :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: KEINV(0:TVEC1, 1:NUMPHASES)
    INTEGER :: INCR
    REAL(RK) :: DTIME
    LOGICAL, INTENT(INOUT) :: CONVERGED_SOLUTION
    INTEGER :: AUTO_TIME
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: M_EL
    REAL(RK) :: D_VEC_GRN(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: W_VEC_GRN(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: EPSEFF_LAT(0:NGRAIN1, EL_SUB1:EL_SUP1)
    !
    !---------------------------------------------------------------------------
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    ! This part should be cleaned when NGRAIN1 stuff is removed
    ! Spread {d}_sm & {w}_sm to all grains in aggregate (Taylor Assumption)
    !
    DO I = 0, TVEC1
        !
        D_VEC_GRN(I, :, :) = SPREAD(D_VEC(I, :), DIM = 1, NCOPIES = NGRAIN)
        !
    END DO
    !
    DO I = 0, DIMS1
        !
        W_VEC_GRN(I, :, :) = SPREAD(W_VEC(I, :), DIM = 1, NCOPIES = NGRAIN)
        !
    END DO
    !
    ! Spread over grains: EPSEFF --> EPSEFF_LAT
    !
    EPSEFF_LAT = SPREAD(EPSEFF, DIM = 1, NCOPIES = NGRAIN)
    !
    ! Solve for State.
    !
    ! Deviatoric
    !
    CALL SOLVE_STATE_DEV_EVPS_QP(D_VEC_GRN, W_VEC_GRN, C0_ANGS, C_ANGS, &
        & SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, E_BAR_VEC, KEINV, &
        & EPSEFF_LAT, JITER_STATE, INCR, DTIME, CONVERGED_SOLUTION, AUTO_TIME)
    !
    IF (.NOT. CONVERGED_SOLUTION .AND. AUTO_TIME .EQ. 1) RETURN
    !
    ! Volumetric
    !
    CALL SOLVE_STATE_VOL_EVPS_QP(E_ELAS_KK_BAR, E_ELAS_KK, D_KK, SIG_KK, &
        & DTIME, M_EL)
    !
    END SUBROUTINE POLYCRYSTAL_RESPONSE_EVPS_QP
    !
    !===========================================================================
    !
    SUBROUTINE SOLVE_STATE_DEV_EVPS_QP(D_VEC, W_VEC, C0_ANGS, C_ANGS, &
        & SIG_LAT_N, SIG_LAT, CRSS_N, CRSS, RSTAR_N, RSTAR, E_BAR_VEC, KEINV, &
        & EPSEFF, JITER_STATE, INCR, DTIME, CONVERGED_SOLUTION, AUTO_TIME)
    !
    ! Solves for the deviatoric state at quad points
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! D_VEC:
    ! W_VEC:
    ! C0_ANGS:
    ! C_ANGS:
    ! SIG_LAT_N:
    ! SIG_LAT:
    ! CRSS_N:
    ! CRSS:
    ! RSTAR_N:
    ! RSTAR:
    ! E_BAR_VEC:
    ! KEINV:
    ! EPSEFF:
    ! JITER_STATE:
    ! INCR:
    ! DTIME:
    ! CONVERGED_SOLUTION:
    ! AUTO_TIME:
    !
    REAL(RK), INTENT(IN) :: D_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: W_VEC(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: SIG_LAT_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: SIG_LAT(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: KEINV(0:TVEC1, 1:NUMPHASES)
    REAL(RK), INTENT(IN) :: EPSEFF(0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER, INTENT(OUT) :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER :: INCR
    REAL(RK) :: DTIME
    LOGICAL, INTENT(INOUT) :: CONVERGED_SOLUTION
    INTEGER :: AUTO_TIME
    !
    ! Locals:
    !
    LOGICAL, PARAMETER :: VP_LOG = .FALSE.
    LOGICAL :: DONE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    LOGICAL :: CONVERGED_NEWTON
    LOGICAL :: CONVERGED_STATE
    INTEGER :: ITER_STATE
    INTEGER :: M_EL
    INTEGER :: ISLIP
    INTEGER :: N_SLIP
    REAL(RK) :: QR5X5(0:TVEC1, 0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: QR3X3(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: W_VEC_LAT(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: WP_HAT(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: D_VEC_LAT(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: E_BAR_VEC_R(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: CRSS_0(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: NORM_S_0(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: NORM_S(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DIFF_NORM_S(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DIFF_CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: D_RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER :: IPHASE
    INTEGER :: K
    INTEGER :: MY_PHASE(0:(EL_SUP1-EL_SUB1))
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    CONVERGED_NEWTON = .TRUE.
    CONVERGED_STATE  = .TRUE.
    !
    JITER_STATE  = 0
    !
    CRSS = CRSS_N
    DONE = .FALSE.
    !
    ! Estimate for the stresses :
    ! For INCR = 1 --> based on viscoplastic solution
    ! For INCR > 1 --> based on previous solution (unrotated)
    !
    IF (INCR .EQ. 1) THEN
        !
        ! C_ANGS [3x3] --> QR5X5 [5x5]
        !
        CALL ROT_MAT_SYMM(C_ANGS, QR5X5, NGRAIN, M_EL)
        !
        ! D_VEC(sample coords) --> D_VEC_LAT(crystal coords)
        ! {D_VEC_LAT} = [QR5X5]'{D_VEC}
        !
        CALL LATTICE_DEFORM(QR5X5, D_VEC, D_VEC_LAT, NGRAIN, M_EL)
        !
        ! Estimate of SIG_LAT (visco-plastic solution)
        !
        CALL STRESS_SOLVE_VP(SIG_LAT, D_VEC_LAT, CRSS, EPSEFF, VP_LOG)
        !
        ! Use a fraction of the viscoplastic solution as initial guess
        !
        SIG_LAT = 0.5 * SIG_LAT
        !
    ELSE
        !
        ! Use value from previous increment
        !
        SIG_LAT = SIG_LAT_N
        !
    END IF
    !
    ! First (Forward Euler) estimate for CRSS and RSTAR (c)
    !
    ! C_ANGS [3x3] --> QR3X3 [3x3]
    !
    CALL ROT_MAT_SKEW(C_ANGS, QR3X3, NGRAIN, M_EL)
    !
    ! W_VEC(sample coords) --> W_VEC_LAT(crystal coords)
    !
    CALL LATTICE_SPIN(QR3X3, W_VEC, W_VEC_LAT, NGRAIN, M_EL)
    !
    ! Estimate CRSS (g), RSTAR (R*) and D_RSTAR (dR*)
    ! Calculate also C_ANGS: [c]=[c_0]*[R*]
    !
    CALL RSTARN_SOLVE(CRSS_N, CRSS, RSTAR_N, RSTAR, C0_ANGS, C_ANGS, &
        & SIG_LAT, W_VEC_LAT, DTIME, EPSEFF, done, D_RSTAR, UPD_EULER_FWD)
    !
    ! Compute 2-norm for array of 5-vectors
    !
    CALL NORM_VEC(NORM_S_0, SIG_LAT, NGRAIN, M_EL)
    !
    ! Update CRSS_0
    !
    CRSS_0 = CRSS
    !
    ! Iterate for the material state
    !
    ITER_STATE = 1
    !
    DO WHILE ((ANY(.NOT. DONE)) .AND. (ITER_STATE .LE. &
        & CV_OPTIONS%SX_MAX_ITERS_STATE))
        !
        ! C_ANGS [3x3] --> QR5X5 [5x5]
        !
        CALL ROT_MAT_SYMM(C_ANGS, QR5X5, NGRAIN, M_EL)
        !
        ! D_VEC(sample coords) --> D_VEC_LAT(crystal coords)
        !
        CALL LATTICE_DEFORM(QR5X5, D_VEC, D_VEC_LAT, NGRAIN, M_EL)
        !
        ! D_RSTAR [3x3] --> QR5X5 [5x5]
        !
        CALL ROT_MAT_SYMM(D_RSTAR, QR5X5, NGRAIN, M_EL)
        !
        ! Apply dR* to E_BAR_VEC --> E_BAR_VEC_R
        !
        CALL LATTICE_DEFORM(QR5X5, E_BAR_VEC, E_BAR_VEC_R, NGRAIN, M_EL)
        !
        ! --> SIG_LAT
        !
        CALL STRESS_SOLVE_EVPS(SIG_LAT, D_VEC_LAT, W_VEC_LAT, E_BAR_VEC_R, &
            & CRSS, KEINV, DTIME, WP_HAT, ITER_STATE, DONE, CONVERGED_NEWTON)
        !
        IF (.NOT. CONVERGED_NEWTON .AND. AUTO_TIME .EQ. 1) GO TO 100
        !
        ! C_ANGS [3x3] --> QR3X3 [3x3]
        !
        CALL ROT_MAT_SKEW(C_ANGS, QR3X3, NGRAIN, M_EL)
        !
        ! W_VEC(sample coords) --> W_VEC_LAT(crystal coords)
        !
        CALL LATTICE_SPIN(QR3X3, W_VEC, W_VEC_LAT, NGRAIN, M_EL)
        !
        ! Calculate CRSS (g), RSTAR (R*) and D_RSTAR (dR*)
        ! Calculate also C_ANGS: [c]=[c_0]*[R*]
        !
        CALL RSTARN_SOLVE(CRSS_N, CRSS, RSTAR_N, RSTAR, C0_ANGS, C_ANGS, &
            & SIG_LAT, WP_HAT, DTIME, EPSEFF, DONE, D_RSTAR, UPD_EULER_BWD)
        !
        CALL NORM_VEC(NORM_S, SIG_LAT, NGRAIN, M_EL)
        !
        ! DEB This section was originally dones with nested `WHERE' constructs,
        !   but was changed because the AIX compiler rejected them, although the
        !   CM compiler had no problem.
        !
        DO ISLIP = 0, MAXSLIP1
            !
            WHERE (.NOT. DONE)
                !
                DIFF_CRSS(ISLIP, :, :) = DABS(CRSS(ISLIP, : ,:) - &
                    & CRSS_0(ISLIP, :, :))
                !
            END WHERE
            !
        END DO
        !
        WHERE (.NOT. DONE)
            !
            DIFF_NORM_S = DABS(NORM_S - NORM_S_0)
            !
        END WHERE
        !
        !
        DO K = 0, NGRAIN - 1
            !
            DO IPHASE = 1, NUMPHASES
                !
                CALL CRYSTALTYPEGET(CTYPE(IPHASE))
                !
                N_SLIP = CTYPE(IPHASE)%NUMSLIP
                !
                DO ISLIP = 0, N_SLIP - 1
                    !
                    ! Currently have an initial CRSS_0 scaling factor from the
                    !   simple latent hardening model but might end up getting
                    !   rid of it later versions
                    WHERE ((.NOT. DONE(K, :)) .AND. (MY_PHASE .EQ. IPHASE) &
                        & .AND. (DIFF_NORM_S(K, :) .LT. (TOLER_STATE * &
                            & CRYSTAL_PARM(3, IPHASE))) &
                        & .AND. (DIFF_CRSS(ISLIP, K, :) .LT. (TOLER_STATE * &
                            & CRYSTAL_PARM(3, IPHASE))))
                        !
                        DONE(K, :) = .TRUE.
                        JITER_STATE(K, :) = ITER_STATE
                        !
                    END WHERE
                    !
                END DO !N_SLIP
                !
            END DO !NUMPHASES
            !
        END DO !NGRAIN
        !
        DO ISLIP = 0, MAXSLIP1
            !
            WHERE (.NOT. DONE)
                !
                NORM_S_0 = NORM_S
                CRSS_0(ISLIP, :, :) = CRSS(ISLIP, :, :)
                !
            END WHERE
            !
        END DO
        !
        ITER_STATE = ITER_STATE + 1
        !
    END DO ! DO WHILE
    !
    IF (ANY(.NOT. DONE)) THEN
        !
        CONVERGED_STATE = .FALSE.
        WRITE(DFLT_U, '(A)') 'Warning:       . Not all crystals converged.'
        !WRITE(DFLT_U, *) 'Warning:       . Crystals = ', NGRAIN * M_EL, ', &
        !    &converged = ', COUNT(DONE)
        !
    END IF
    !
    100 IF ((.NOT. CONVERGED_NEWTON) .OR. (.NOT. CONVERGED_STATE)) &
        & CONVERGED_SOLUTION = .FALSE.
    !
    END SUBROUTINE SOLVE_STATE_DEV_EVPS_QP
    !
    !===========================================================================
    !
    SUBROUTINE SOLVE_STATE_VOL_EVPS_QP(E_ELAS_KK_BAR, E_ELAS_KK, D_KK, SIG_KK, &
        & DTIME, M)
    !
    ! Solves for the volumetric state at quad points
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! E_ELAS_KK_BAR:
    ! E_ELAS_KK:
    ! D_KK:
    ! SIG_KK:
    ! DTIME:
    ! M:
    !
    REAL(RK), INTENT(IN) :: E_ELAS_KK_BAR(0:(M - 1))
    REAL(RK), INTENT(OUT) :: E_ELAS_KK(0:(M - 1))
    REAL(RK), INTENT(IN) :: D_KK(0:(M - 1))
    REAL(RK), INTENT(OUT) :: SIG_KK(0:(M - 1))
    REAL(RK), INTENT(IN) :: DTIME
    INTEGER :: M
    !
    ! Locals:
    !
    INTEGER :: MY_PHASE(0:M - 1)
    INTEGER :: IPHASE
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    !
    ! Compute volumetric response in sample reference frame
    !
    E_ELAS_KK = E_ELAS_KK_BAR + DTIME * D_KK
    !
    DO IPHASE = 1, NUMPHASES
        !
        WHERE (MY_PHASE .EQ. IPHASE)
            !
            SIG_KK = 3.0D0 * CRYSTAL_PARM(8,IPHASE) * E_ELAS_KK
            !
        END WHERE
        !
    END DO !NUMPHASES
    !
    END SUBROUTINE SOLVE_STATE_VOL_EVPS_QP
    !
END MODULE POLYCRYSTAL_RESPONSE_EVPS_MOD
