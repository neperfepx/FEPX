! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE HARDENING_MOD
!
! Module containing definition of hardening evolution equations
!
! Contains subroutines:
! HARD_LAW: Evaluate hardening rate or derivative of hardening rate
! ISOTROPIC_HARDENING: Isotropic hardening assumption
! CYCLIC_HARDENING: Cyclic hardening assumption
! ANISOTROPIC_HARDENING: Anisotropic hardening assumption
! CALCULATE_SHRATE: Calculate the shear rate
!
! From libf95:
!
USE LIBF95, ONLY: RK=>REAL_KIND
!
! From libfepx:
!
USE DIMENSIONS_MOD
USE MATRIX_OPERATIONS_MOD
USE MICROSTRUCTURE_MOD
USE READ_INPUT_MOD
!
! From libparallel:
!
USE PARALLEL_MOD, ONLY: QUIT=>PAR_QUIT
!
IMPLICIT NONE
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: HARD_LAW_VOCE
PUBLIC :: EVAL_RATE
PUBLIC :: EVAL_RATE_DER
PUBLIC :: HARD_LAW
!
INTEGER, PARAMETER :: HARD_LAW_VOCE = 2
INTEGER, PARAMETER :: EVAL_RATE = 1
INTEGER, PARAMETER :: EVAL_RATE_DER = 2
INTEGER :: N_SLIP
!
CONTAINS
    !
    SUBROUTINE HARD_LAW(FUNC, DFUNC, CRSS, CRSS_SAT, SHEAR, SHRATE, EPSEFF, &
        & ICODE, IHARD, N, M)
    !
    ! Evaluate hardening rate or derivative of hardening rate
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! FUNC: Array of computed hardening rates
    ! DFUNC: Array of computed hardening rate derivatives
    ! CRSS: Array of critical resolved shear stress (hardness)
    ! CRSS_SAT: Array of satuRATIOn CRSS
    ! SHEAR: Shear values
    ! SHRATE: Shear rates
    ! EPSEFF:
    ! ICODE: Hardening law
    ! IHARD: Flag determining calculation of rate or derivative of rate
    ! N, M: Array dimensions
    !
    REAL(RK), INTENT(OUT) :: FUNC(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: DFUNC(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: CRSS(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: CRSS_SAT(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: SHRATE(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: SHEAR(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: EPSEFF(0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: ICODE
    INTEGER, INTENT(IN) :: IHARD
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    ! MPK: Do any of these need to be defined here? They have been moved to the
    !   individual subroutines...
    ! 
    REAL(RK) :: MYSIGN(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: RATIO(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    INTEGER :: NUMIND
    INTEGER :: N_SLIP
    INTEGER, POINTER :: INDICES(:)=>NULL()
    INTEGER :: IPHASE
    INTEGER :: ISLIP
    INTEGER :: MY_PHASE(0:(M - 1))
    REAL(RK) :: C3
    REAL(RK) :: C2
    REAL(RK) :: C10
    !
    !---------------------------------------------------------------------------
    !
    IF (IHARD /= HARD_LAW_VOCE) THEN
        !
        CALL PAR_QUIT('Error  :     > Voce hardening must be used.')
        !
    END IF
    !
    ! RC: Initialize these to zero or else the terms that are not used in a dual
    !   phase material get assigned random values.
    !
    FUNC = 0.0D0
    DFUNC = 0.0D0
    !
    ! RC: In future updates this will need to be set up such that it makes a
    !   CALL for either a isotropic case or latent hardening case
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    !
    MYSIGN = 1.0D0
    !
    DO ISLIP=0,MAXSLIP1
        !
        WHERE((CRSS_SAT - CRSS(ISLIP,:,:)) .LE. 0.0D0)
            !
            MYSIGN(ISLIP,:,:) = 0.0D0
            !
        END WHERE
        !
    END DO
    !
    ! Options for which hardening model to use
    !
    IF (OPTIONS%HARD_TYPE .EQ. 'isotropic') THEN
        !
        CALL ISOTROPIC_HARDENING(FUNC, DFUNC, CRSS, CRSS_SAT, SHRATE, ICODE, &
            & IHARD, N, M, MYSIGN, MY_PHASE)
        !
    ELSEIF (OPTIONS%HARD_TYPE .EQ. 'cyclic_isotropic') THEN
        !
        CALL CYCLIC_HARDENING(FUNC, DFUNC, CRSS, CRSS_SAT, SHEAR, SHRATE, &
            & EPSEFF, ICODE, IHARD, N, M, MYSIGN, MY_PHASE)
        !
    ELSEIF (OPTIONS%HARD_TYPE .EQ. 'cyclic_anisotropic') THEN
        !
        CALL PAR_QUIT('Error  :     > Anisotropic cyclic hardening is not currently supported.')
        !
    ELSEIF (OPTIONS%HARD_TYPE .EQ. 'latent') THEN
        !
        CALL ANISOTROPIC_HARDENING(FUNC, DFUNC, CRSS, CRSS_SAT,ICODE, IHARD, &
            & N, M, MYSIGN, MY_PHASE)
        !
    ELSE
        !
        CALL PAR_QUIT('Error  :     > Invalid hardening model option input.')
        !
    END IF
    !
    END SUBROUTINE HARD_LAW
    !
    !===========================================================================
    !
    SUBROUTINE ISOTROPIC_HARDENING(FUNC, DFUNC, CRSS, CRSS_SAT, SHRATE, ICODE, &
        & IHARD, N, M, MYSIGN, MY_PHASE)
    !
    ! Isotropic hardening assumption
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! FUNC: Array of computed hardening rates
    ! DFUNC: Array of computed hardening rate derivatives
    ! CRSS: Array of critical resolved shear stress (hardness)
    ! CRSS_SAT: Array of satuRATIOn CRSS
    ! SHRATE: Shear rates
    ! ICODE: Hardening law
    ! IHARD: Flag determining calculation of rate or derivative of rate
    ! N, M: Array dimensions
    ! MYSIGN:
    ! MY_PHASE:
    !
    REAL(RK), INTENT(OUT) :: FUNC(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: DFUNC(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: CRSS(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: CRSS_SAT(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: SHRATE(0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: ICODE
    INTEGER, INTENT(IN) :: IHARD
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    REAL(RK), INTENT(IN) :: MYSIGN(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: MY_PHASE(0:(M - 1))
    !
    ! Locals:
    !
    REAL(RK) :: RATIO(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    INTEGER :: NUMIND
    INTEGER :: N_SLIP
    INTEGER, POINTER :: INDICES(:)=>NULL()
    INTEGER :: IPHASE
    INTEGER :: ISLIP
    REAL(RK) :: C3
    REAL(RK) :: C2
    REAL(RK) :: C10
    !
    !---------------------------------------------------------------------------
    !
    DO IPHASE = 1, NUMPHASES
        !
        C2 = CRYSTAL_PARM(2, IPHASE)
        C3 = CRYSTAL_PARM(3, IPHASE)
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE))
        !
        N_SLIP = CTYPE(IPHASE)%NUMSLIP
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
        !
        ! RC: All of the slip systems are just looped over instead of just the
        !   one. Since only the isotropic CASE is being examined here it should
        !   be possible to do only one computation of the slip system and THEN
        !   assign all the other slip systems that value. Future update should
        !   look into that implementation.
        !
        DO ISLIP = 0, N_SLIP-1
            !
            IF (ISLIP==0) THEN
                !
                RATIO(ISLIP, :, INDICES) = (CRSS_SAT(:, INDICES) - &
                    & CRSS(ISLIP, :, INDICES)) / (CRSS_SAT(:, INDICES) - C3)
                !
                SELECT CASE(ICODE)
                !
                CASE (EVAL_RATE)
                    !
                    FUNC(ISLIP,:, INDICES) = SHRATE(:, INDICES) * &
                        & MYSIGN(ISLIP, :, INDICES) * C2 * &
                        & RATIO(ISLIP, :, INDICES) **N_VOCE(IPHASE)
                    !
                CASE (EVAL_RATE_DER)
                    !
                    DFUNC(ISLIP, : , INDICES) = -C2 /(CRSS_SAT(:, INDICES) - &
                        & CRYSTAL_PARM(3, IPHASE)) * SHRATE(:, INDICES) * &
                        & MYSIGN(ISLIP, :, INDICES) * N_VOCE(IPHASE)
                    !
                CASE DEFAULT
                    !
                END SELECT
                !
            ELSE
                !
                SELECT CASE(ICODE)
                !
                CASE (EVAL_RATE)
                    !
                    FUNC(ISLIP, :, INDICES) = FUNC(0, :, INDICES)
                    !
                CASE (EVAL_RATE_DER)
                    !
                    DFUNC(ISLIP, :, INDICES) = DFUNC(0, :, INDICES)
                    !
                CASE DEFAULT
                    !
                END SELECT
                !
            END IF
            !
        END DO !N_SLIP
        !
        DEALLOCATE(INDICES)
        !
    END DO !NUMPHASES
    !
    END SUBROUTINE ISOTROPIC_HARDENING
    !
    !===========================================================================
    !
    SUBROUTINE CYCLIC_HARDENING(FUNC, DFUNC, CRSS, CRSS_SAT, SHEAR, SHRATE, &
        & EPSEFF, ICODE, IHARD, N, M, MYSIGN, MY_PHASE)
    !
    ! Cyclic hardening assumption
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! FUNC: Array of computed hardening rates
    ! DFUNC: Array of computed hardening rate derivatives
    ! CRSS: Array of critical resolved shear stress (hardness)
    ! CRSS_SAT: Array of satuRATIOn CRSS
    ! SHEAR: Shear values
    ! SHRATE: Shear rates
    ! EPSEFF:
    ! ICODE: Hardening law
    ! IHARD: Flag determining calculation of rate or derivative of rate
    ! N, M: Array dimensions
    ! MYSIGN:
    ! MY_PHASE:
    !
    REAL(RK), INTENT(OUT) :: FUNC(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: DFUNC(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: CRSS(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: CRSS_SAT(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: SHRATE(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: SHEAR(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: EPSEFF(0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: ICODE
    INTEGER, INTENT(IN) :: IHARD
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    REAL(RK), INTENT(IN) :: MYSIGN(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: MY_PHASE(0:(M - 1))
    !
    ! Locals:
    !
    REAL(RK) :: RATIO(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: RATIO_SH(0:(N - 1), 0:(M - 1))
    REAL(RK) :: ACTIVE_SHRATE(0:(N - 1), 0:(M - 1))
    REAL(RK) :: SHEAR_CRIT(0:(N - 1), 0:(M - 1))
    INTEGER :: NUMIND
    INTEGER :: N_SLIP
    INTEGER :: N_GT1
    INTEGER :: N_LT0
    INTEGER, POINTER :: INDICES(:)=>NULL()
    INTEGER :: IPHASE
    INTEGER :: ISLIP
    REAL(RK) :: C3
    REAL(RK) :: C2
    REAL(RK) :: C10
    REAL(RK) :: SHR_MAX(0:(N - 1), 0:(M - 1))
    REAL(RK) :: SHR_MIN(0:(N - 1), 0:(M - 1))
    REAL(RK), PARAMETER :: MACH_EPS = 2.22D-16
    !
    !---------------------------------------------------------------------------
    !
    ACTIVE_SHRATE(:,:) = 0.0D0
    SHEAR_CRIT = 0.0D0
    N_LT0 = 0
    !
    DO IPHASE = 1, NUMPHASES
        !
        C2 = CRYSTAL_PARM(2,IPHASE)
        C3 = CRYSTAL_PARM(3,IPHASE)
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE))
        N_SLIP = CTYPE(IPHASE)%NUMSLIP
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
        !
        ! RC: All of the slip systems are just looped over instead of just the
        !   one. Since only the isotropic case is being examined here it should
        !   be possible to DO only one computation of the slip system and then
        !   assign all the other slip systems that value. Future update should
        !   look into that implementation.
        !
        RATIO_SH(:, INDICES) = DABS(CRSS(0, :, INDICES) / CRSS_SAT(:, INDICES))
        !
        SHEAR_CRIT(:, INDICES) = CYCLIC_PARM(0, IPHASE) * &
            & (RATIO_SH(:, INDICES) ** CYCLIC_PARM(1, IPHASE))
        !
        WHERE(SHEAR_CRIT(:, INDICES) .LT. MACH_EPS)
            !
            SHEAR_CRIT(:, INDICES) = 0.0D0
            !N_LT0 = N_LT0 + 1
            !
        END WHERE
        !
        IF (ANY(DABS(SHEAR_CRIT) .GE. 1.0D0)) THEN
            !
            N_GT1 = COUNT((DABS(SHEAR_CRIT) .GE. 1.0D0))
            !print *, 'Number greater than 1: ', N_GT1
            CALL PAR_QUIT('Error  :     > `SHEAR_CRIT` is greater than 1.')
            !
        END IF
        !
        DO ISLIP = 0, N_SLIP-1
            !
            WHERE (ACCUMSHEAR(ISLIP, 0, INDICES + EL_SUB1) .GE. &
                & SHEAR_CRIT(0, INDICES))
                !
                ACTIVE_SHRATE(0, INDICES) = ACTIVE_SHRATE(0, INDICES) + &
                    & DABS(SHEAR(ISLIP, 0, INDICES))
                !
            END WHERE
            !
        END DO
        !
        SHR_MIN = 1.0D-6 * EPSEFF
        SHR_MAX = 1.0D1 * EPSEFF
        !
        WHERE (ACTIVE_SHRATE .LE. SHR_MIN)
            !
            ACTIVE_SHRATE = SHR_MIN
            !
        END WHERE
        !
        WHERE (ACTIVE_SHRATE .GE. SHR_MAX)
            !
            ACTIVE_SHRATE = SHR_MAX
            !
        END WHERE
        !
        DO ISLIP = 0, N_SLIP - 1
            !
            IF (ISLIP==0) THEN
                !
                RATIO(ISLIP,:, INDICES) = (CRSS_SAT(:, INDICES) - &
                    & CRSS(ISLIP,:, INDICES)) / (CRSS_SAT(:, INDICES) - C3)
                !
                SELECT CASE(ICODE)
                !
                CASE (EVAL_RATE)
                    !
                    FUNC(ISLIP, :, INDICES) = ACTIVE_SHRATE(:, INDICES) * &
                        & MYSIGN(ISLIP, :, INDICES) * C2 * &
                        & RATIO(ISLIP, :, INDICES) ** N_VOCE(IPHASE)
                    !
                CASE (EVAL_RATE_DER)
                    !
                    DFUNC(ISLIP, :, INDICES) = -C2 / (CRSS_SAT(:, INDICES) - &
                        & CRYSTAL_PARM(3, IPHASE)) * ACTIVE_SHRATE(:, INDICES) &
                        & * MYSIGN(ISLIP, :, INDICES) * N_VOCE(IPHASE)
                    !
                CASE DEFAULT
                    !
                END SELECT
                !
            ELSE
                !
                SELECT CASE(ICODE)
                !
                CASE (EVAL_RATE)
                    !
                    FUNC(ISLIP, :, INDICES) = FUNC(0, :, INDICES)
                    !
                CASE (EVAL_RATE_DER)
                    !
                    DFUNC(ISLIP, :, INDICES) = DFUNC(0, :, INDICES)
                    !
                CASE DEFAULT
                    !
                END SELECT
                !
            END IF
            !
        END DO !N_SLIP
        !
        DEALLOCATE(INDICES)
        !
    END DO !NUMPHASES
    !
    END SUBROUTINE CYCLIC_HARDENING
    !
    !===========================================================================
    !
    SUBROUTINE ANISOTROPIC_HARDENING(FUNC, DFUNC, CRSS, CRSS_SAT, ICODE, &
        & IHARD, N, M, MYSIGN, MY_PHASE)
    !
    ! Anisotropic hardening assumption
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! FUNC: Array of computed hardening rates
    ! DFUNC: Array of computed hardening rate derivatives
    ! CRSS: Array of critical resolved shear stress (hardness)
    ! CRSS_SAT: Array of satuRATIOn CRSS
    ! ICODE: Hardening law
    ! IHARD: Flag determining calculation of rate or derivative of rate
    ! N, M: Array dimensions
    ! MYSIGN:
    ! MY_PHASE:
    !
    REAL(RK), INTENT(OUT) :: FUNC(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: DFUNC(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: CRSS(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: CRSS_SAT(0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: ICODE
    INTEGER, INTENT(IN) :: IHARD
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    REAL(RK), INTENT(IN) :: MYSIGN(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: MY_PHASE(0:(M - 1))
    !
    ! Locals:
    !
    REAL(RK) :: RATIO(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: SHRATE(0:MAXSLIP1)
    INTEGER :: NUMIND
    INTEGER :: NSLIP
    INTEGER :: XTYPE
    INTEGER, ALLOCATABLE :: INDICES(:)
    INTEGER :: IPHASE
    INTEGER :: IELMS
    INTEGER :: DIM
    INTEGER :: INDEX
    REAL(RK) :: C3
    REAL(RK) :: C2
    REAL(RK) :: C10
    LOGICAL :: LOGVEC(0:(M - 1))
    !
    !---------------------------------------------------------------------------
    !
    SHRATE = 0.0D0
    !
    DO IPHASE = 1, NUMPHASES
        !
        C2 = CRYSTAL_PARM(2,IPHASE)
        C3 = CRYSTAL_PARM(3,IPHASE)
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE))
        !
        NSLIP = CTYPE(IPHASE)%NUMSLIP - 1
        !
        IF (OPTIONS%HARD_TYPE .EQ. "test") THEN
            !
            ! RC: Set it outside the available crystal types to trigger the
            !   default calculation of the isotropic case
            !
            XTYPE = 4
            !
        ELSE
            !
            XTYPE = CTYPE(IPHASE)%CLASS
            !
        END IF
        !
        LOGVEC = (MY_PHASE .EQ. IPHASE)
        ! CALL INDEXFinder(LOGVEC, INDICES)
        DIM = SIZE(LOGVEC)
        !
        DO IELMS = 0, DIM - 1
            !
            ! Find where the element is the same phase IPHASE
            !
            IF (LOGVEC(IELMS)) THEN
                !
                RATIO(0:NSLIP, 0, IELMS) = (CRSS_SAT(0, IELMS) - &
                    & CRSS(0:NSLIP, 0, IELMS)) / (CRSS_SAT(0, IELMS) - C3)
                !
                CALL CALCULATE_SHRATE(SHRATE(:), XTYPE, IELMS)
                !
                SELECT CASE(ICODE)
                !
                CASE (EVAL_RATE)
                    !
                    ! All of the slip system FUNC and DFUNC values for the
                    !   element are calculated in one go
                    !
                    FUNC(0:NSLIP, 0, IELMS) = SHRATE(0:NSLIP) * &
                        & MYSIGN(0:NSLIP, 0, IELMS) * C2 * &
                        & RATIO(0:NSLIP, 0, IELMS) ** N_VOCE(IPHASE)
                    !
                CASE (EVAL_RATE_DER)
                    !
                    DFUNC(0:NSLIP, 0, IELMS) = -C2 /(CRSS_SAT(0, IELMS) - &
                        & CRYSTAL_PARM(3, IPHASE)) * SHRATE(0:NSLIP) * &
                        & MYSIGN(0:NSLIP, 0, IELMS) * N_VOCE(IPHASE)
                    !
                CASE DEFAULT
                    !
                END SELECT
                !
            END IF
            !
        END DO !N_SLIP
        !
        ! DEALLOCATE(INDICES)
        !
    END DO !NUMPHASES
    !
    END SUBROUTINE ANISOTROPIC_HARDENING
    !
    !===========================================================================
    !
    SUBROUTINE CALCULATE_SHRATE(SHRATE, XTYPE, IND)
    !
    ! SHRATEcalculations are based upon the assumption that slip systems that
    !   share the same slip plane interact instead of the same slip direction.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! SHRATE: Shear rate
    ! XTYPE: Crystal type
    ! IND:
    REAL(RK), INTENT(OUT) :: SHRATE(0:MAXSLIP1)
    INTEGER, INTENT(IN) :: XTYPE
    INTEGER, INTENT(IN) :: IND
    !
    ! Locals:
    !
    REAL(RK) :: TMP(0:1)
    REAL(RK) :: TMPO(0:1)
    REAL(RK) :: ONES(0:MAXSLIP1,0:MAXSLIP1)
    !
    !---------------------------------------------------------------------------
    !
    ONES = 1.0D0
    !
    SELECT CASE(XTYPE) ! Values are hard wired
    !
    CASE (1) !FCC XTYPE
        !
        CALL MATRIX_VEC_MULT(FCC_H1, ABS(GAMMADOT(0:2, 0, IND + EL_SUB1)), &
            & SHRATE(0:2), 3)
        CALL MATRIX_VEC_MULT(FCC_H2, ABS(GAMMADOT(3:5, 0, IND + EL_SUB1)), &
            & SHRATE(3:5), 3)
        CALL MATRIX_VEC_MULT(FCC_H3, ABS(GAMMADOT(6:8, 0, IND + EL_SUB1)), &
            & SHRATE(6:8), 3)
        CALL MATRIX_VEC_MULT(FCC_H4, ABS(GAMMADOT(9:11, 0, IND + EL_SUB1)), &
            & SHRATE(9:11), 3)
        !
    CASE (3) !HCP XTYPE
        !
        CALL MATRIX_VEC_MULT(HCP_H1, ABS(GAMMADOT(0:2, 0, IND + EL_SUB1)), &
            & SHRATE(0:2), 3)
        !
        SHRATE(3:5) = HCP_VERT * ABS(GAMMADOT(3:5 ,0, IND + EL_SUB1))
        !
        CALL MATRIX_VEC_MULT(HCP_H2, ABS(GAMMADOT(6:7, 0, IND + EL_SUB1)), &
            & SHRATE(6:7), 2)
        CALL MATRIX_VEC_MULT(HCP_H3, ABS(GAMMADOT(8:9, 0, IND + EL_SUB1)), &
            & SHRATE(8:9), 2)
        CALL MATRIX_VEC_MULT(HCP_H4, ABS(GAMMADOT(10:11, 0, IND + EL_SUB1)), &
            & SHRATE(10:11), 2)
        CALL MATRIX_VEC_MULT(HCP_H5, ABS(GAMMADOT(12:13, 0, IND + EL_SUB1)), &
            & SHRATE(12:13), 2)
        CALL MATRIX_VEC_MULT(HCP_H6, ABS(GAMMADOT(14:15, 0, IND + EL_SUB1)), &
            & SHRATE(14:15), 2)
        CALL MATRIX_VEC_MULT(HCP_H7, ABS(GAMMADOT(16:17, 0, IND + EL_SUB1)), &
            & SHRATE(16:17), 2)
        !
    CASE (2) !BCC XTYPE
        !
        TMP = ABS(GAMMADOT((/ 0,9/), 0, IND+EL_SUB1))
        CALL MATRIX_VEC_MULT(BCC_H1, TMP, TMPO,2)
        SHRATE((/ 0,9/)) = TMPO
        !
        TMP = ABS(GAMMADOT((/ 1,7/), 0, IND+EL_SUB1))
        CALL MATRIX_VEC_MULT(BCC_H2, TMP, TMPO,2)
        SHRATE((/ 1,7/)) = TMPO
        !
        TMP = ABS(GAMMADOT((/ 2,5/), 0, IND+EL_SUB1))
        CALL MATRIX_VEC_MULT(BCC_H3, TMP, TMPO,2)
        SHRATE((/ 2,5/)) = TMPO
        !
        TMP = ABS(GAMMADOT((/ 3,6/), 0, IND+EL_SUB1))
        CALL MATRIX_VEC_MULT(BCC_H4, TMP, TMPO,2)
        SHRATE((/ 3,6/)) = TMPO
        !
        TMP = ABS(GAMMADOT((/ 4,10/), 0, IND+EL_SUB1))
        CALL MATRIX_VEC_MULT(BCC_H5, TMP, TMPO,2)
        SHRATE((/ 4,10/)) = TMPO
        !
        TMP = ABS(GAMMADOT((/ 8,11/), 0, IND+EL_SUB1))
        CALL MATRIX_VEC_MULT(BCC_H6, TMP, TMPO,2)
        SHRATE((/ 8,11/)) = TMPO
        !
    CASE DEFAULT
        !
        CALL MATRIX_VEC_MULT(ONES, ABS(GAMMADOT(:, 0, IND + EL_SUB1)), SHRATE, &
            & MAXSLIP1 + 1)
        !
    END SELECT
    !
    END SUBROUTINE CALCULATE_SHRATE
    !
END MODULE HARDENING_MOD
