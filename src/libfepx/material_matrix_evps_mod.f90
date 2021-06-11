! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE MATERIAL_MATRIX_EVPS_MOD
!
! Matrial matrix for the EVPS solution
!
! Contains subroutines:
! MATERIAL_MATRIX_EVPS: Material matrix for the EVPS solution
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
!
! From libfepx:
!
USE ANISO_EVPS_MOD
USE DIMENSIONS_MOD
USE KINEMATICS_MOD
USE MATRIX_OPERATIONS_MOD
USE MICROSTRUCTURE_MOD
USE POLYCRYSTAL_RESPONSE_EVPS_MOD
USE READ_INPUT_MOD
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
PUBLIC :: MATERIAL_MATRIX_EVPS
! MPK: Should these be public from here? These subroutines are located elsewhere
PUBLIC :: VEL_GRADIENT
PUBLIC :: EFF_DEF
!
CONTAINS
    !
    SUBROUTINE MATERIAL_MATRIX_EVPS(STIF, TAN_STIF, FE, DETV, DNDX, DNDY, &
        & DNDZ, GVEL, C0_ANGS, C_ANGS, SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, &
        & RSTAR_N, RSTAR, KEINV, E_BAR_VEC, E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, &
        & JITER_STATE, WTS, EPSEFF, DTIME, INCR, CONVERGED_SOLUTION, &
        & AUTO_TIME, NR)
    !
    ! Matrial matrix for EVPS solution
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! STIF:
    ! TAN_STIF:
    ! FE:
    ! DETV:
    ! DNDX:
    ! DNDY:
    ! DNDZ:
    ! GVEL:
    ! C0_ANGS:
    ! C_ANGS:
    ! SIG_VEC_N:
    ! SIG_VEC:
    ! CRSS_N:
    ! CRSS:
    ! RSTAR_N:
    ! RSTAR:
    ! KEINV:
    ! E_BAR_VEC:
    ! E_ELAS_KK_BAR:
    ! E_ELAS_KK:
    ! SIG_KK
    ! JITER_STATE:
    ! WTS:
    ! EPSEFF:
    ! DTIME:
    ! INCR:
    ! CONVERGED_SOLUTION:
    ! AUTO_TIME:
    ! NR:
    !
    REAL(RK), INTENT(OUT) :: STIF(TVEC, TVEC, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(OUT) :: TAN_STIF(TVEC, TVEC, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(OUT) :: FE(TVEC, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(OUT) :: DETV(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(IN) :: DNDX(0:NNPE, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(IN) :: DNDY(0:NNPE, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(IN) :: DNDZ(0:NNPE, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(IN) :: GVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
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
    REAL(RK), INTENT(IN) :: KEINV(0:TVEC1, 1:NUMPHASES)
    REAL(RK), INTENT(IN) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    REAL(RK), INTENT(IN) :: E_ELAS_KK_BAR(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(OUT) :: SIG_KK(EL_SUB1:EL_SUP1, 0:NQPT1)
    INTEGER, INTENT(OUT) :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: EPSEFF(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: DTIME
    INTEGER :: INCR
    LOGICAL, INTENT(INOUT) :: CONVERGED_SOLUTION
    INTEGER :: AUTO_TIME
    LOGICAL, INTENT(IN) :: NR
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: M_EL
    REAL(RK) :: D(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: W(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: VGRAD(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: W_VEC(0:DIMS1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: D_VEC(0:TVEC1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: D_KK(EL_SUB1:EL_SUP1, 0:NQPT1)
    !
    !---------------------------------------------------------------------------
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    DO I = 0, NQPT1
        !
        ! Compute velocity gradient (VGRAD) and its symm (d) and skew (w) parts
        !
        CALL VEL_GRADIENT(VGRAD, DNDX(:, :, I), DNDY(:, :, I), DNDZ(:, :, I), &
            & GVEL)
        !
        ! Compute:
        ! D_KK: Mean/volumetric/spherical part of the symmetric part of the
        !   velocity gradient
        ! D: Deviatoric part of the symmetric part of the velocity gradient
        ! W: Skew part of the velocity gradient
        !
        CALL SYMM_VGR(D, D_KK(:, I), VGRAD, M_EL)
        CALL SKEW_VGR(W, VGRAD, M_EL)
        !
        ! Convert:
        ! D [3x3] to D_VEC {5}
        ! W [3x3] to W_VEC {3}
        !
        CALL MAT_VEC_SYMM(D, D_VEC(:, :, I), M_EL)
        CALL MAT_VEC_SKEW(W, W_VEC(:, :, I), M_EL)
        !
        ! Calculate
        ! D, DTIME to EPSEFF
        !
        CALL EFF_DEF(EPSEFF(:, I), D, M_EL)
        !
    END DO
    !
    ! Variables @(t+dt):
    ! RSTAR [3x3]
    ! C_ANGS [3x3]
    ! CRSS (1)
    ! SIG_VEC {5}
    ! SIG_KK (1)
    ! E_ELAS_KK (1)
    !
    CALL POLYCRYSTAL_RESPONSE_EVPS(D_VEC, W_VEC, C0_ANGS, C_ANGS, SIG_VEC_N, &
        & SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, E_BAR_VEC, EPSEFF, D_KK, &
        & SIG_KK, E_ELAS_KK_BAR, E_ELAS_KK, JITER_STATE, KEINV, INCR, DTIME, &
        & CONVERGED_SOLUTION, AUTO_TIME)
    !
    IF (.NOT. CONVERGED_SOLUTION .AND. AUTO_TIME .EQ. 1) RETURN
    !
    DO I = 0, NQPT1
        !
        CALL ANISO_EVPS(STIF(:, :, :, I), TAN_STIF(:, :, :, I), FE(:, :, I), &
            & DETV(:, I), C_ANGS(:, :, :, :, I), SIG_VEC(:, :, :, I), &
            & CRSS(:, :, :, I), RSTAR_N, RSTAR(:, :, :, :, I), &
            & E_BAR_VEC(:, :, :, I), WTS, W_VEC(:, :, I), E_ELAS_KK(:, I), &
            & KEINV, DTIME, NR)
        !
    END DO
    !
    END SUBROUTINE MATERIAL_MATRIX_EVPS
    
END MODULE MATERIAL_MATRIX_EVPS_MOD
