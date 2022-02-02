! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE RSTARN_SOLVE_MOD
!
! Module to update the orientations (rstar).
!
! Contains subroutines:
! RSTARN_SOLVE: Add descriptions to these.
! SLIP_QNT:
! UPD_HARD:
! UPD_ROTN:
! GET_C:
!
! To do:
! - Handle masks (e.g. done/converged array) in helper subroutines, e.g. get_c
! - Make use of MATMUL in get_c
!
! From libf95:
!
USE LIBF95, ONLY: RK => REAL_KIND, STDERR => STANDARD_ERROR
!
USE DIMENSIONS_MOD
USE HARDENING_MOD
USE MATRIX_OPERATIONS_MOD
USE MICROSTRUCTURE_MOD
USE READ_INPUT_MOD
USE STRESS_SOLVE_EVPS_MOD
USE STRESS_SOLVE_VP_MOD
USE UNITS_MOD, ONLY: OUNITS, LOG_U, DFLT_U
!
! From libparallel:
!
USE PARALLEL_MOD, ONLY: PAR_MESSAGE
!
IMPLICIT NONE
!
PRIVATE
!
PUBLIC :: RSTARN_SOLVE, SLIP_QNT
PUBLIC :: UPD_EULER_FWD, UPD_EULER_BWD
!
! Flags and Identifiers
!
INTEGER, PARAMETER :: UPD_EULER_FWD = 1, UPD_EULER_BWD = 2
!
! Hardening model iteration limits
!
! Note: The equation should be quadratic for the usual hardening model,
!   and so it should only require one iteration.
!
! MAX_ITER_HARD is now an available user input option from READ_INPUT_MOD
! and is assigned from the options type in UPD_HARD.
!
REAL(RK), PARAMETER :: TOLER_HARD = 1.0D-4, LS_CUTOFF = 0.001D0
INTEGER :: MAX_ITER_HARD
!
CONTAINS
    !
    SUBROUTINE RSTARN_SOLVE(CRSS_N, CRSS, RSTAR_N, RSTAR, C_0, C, SIG_LAT, &
        & WP_HAT, DTIME, EPSEFF, DONE, D_RSTAR, IFLAG)
    !
    ! Add description here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER  :: IFLAG
    REAL(RK) :: DTIME
    !
    LOGICAL, INTENT(IN)     :: DONE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)    :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: RSTAR  (0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: D_RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: CRSS  (0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)    :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT)   :: C  (0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)    :: C_0(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)    :: SIG_LAT(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)    :: WP_HAT (0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)    :: EPSEFF(0:NGRAIN1, EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    REAL(RK) :: SHRATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: WP_SS(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SHEAR(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER  :: M_EL
    !
    !---------------------------------------------------------------------------
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    ! Compute gammadots (shrate) and wp_ss
    !
    CALL SLIP_QNT(WP_SS, SHEAR, SHRATE, SIG_LAT, CRSS, EPSEFF, NGRAIN, M_EL)
    !
    ! Update crss
    !
    CALL UPD_HARD(CRSS, CRSS_N, SHEAR, SHRATE, EPSEFF, DTIME, IFLAG, DONE, &
        & NGRAIN, M_EL)
    !
    ! Update rstar (R*) and d_rstar(dR*)
    !
    CALL UPD_ROTN(RSTAR, RSTAR_N, D_RSTAR, DTIME, WP_HAT, WP_SS, DONE, &
        & NGRAIN, M_EL)
    !
    ! Compute c = c_0 * R*
    !
    CALL GET_C(C, RSTAR, C_0, DONE, NGRAIN, M_EL)
    !
    RETURN
    !
    END SUBROUTINE RSTARN_SOLVE
    !
    !===========================================================================
    !
    SUBROUTINE SLIP_QNT(WP_SS, SHEAR, SHRATE, SIG_LAT, CRSS, EPSEFF, N, M)
    !
    ! Add description here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER :: N_SLIP, N, M
    REAL(RK), INTENT(OUT) :: WP_SS(0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: SHRATE(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: SHEAR(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: SIG_LAT(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: CRSS(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: EPSEFF(0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    INTEGER  :: ISLIP, I, IPHASE, IS, NUMIND
    !
    REAL(RK), POINTER :: PLOCAL(:,:) => NULL()
    REAL(RK), POINTER :: QLOCAL(:,:) => NULL()
    INTEGER,  POINTER :: INDICES(:) => NULL()
    REAL(RK) :: RSS(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: SHR_MIN(0:(N - 1), 0:(M - 1))
    REAL(RK) :: SHR_MAX(0:(N - 1), 0:(M - 1))
    REAL(RK) :: ANISO_M_TEMP(0:MAXSLIP1)
    INTEGER  :: MY_PHASE(0:(M - 1))
    !
    !---------------------------------------------------------------------------
    !
    ! Crystal_type of the elements (on this node)
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    !
    ! Initialize the outputs
    !
    SHRATE(:,:)  = 0.0D0
    WP_SS(:,:,:) = 0.0D0
    SHEAR(:,:,:) = 0.0D0
    !
    DO IPHASE = 1, NUMPHASES
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE), DEV = PLOCAL, SKW = QLOCAL)
        N_SLIP=CTYPE(IPHASE)%NUMSLIP
        !
        ! Compute shear on each slip system
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
        !
        DO ISLIP = 0, N_SLIP - 1
            !
            CALL SS_PROJECT(RSS(ISLIP,:,:), PLOCAL(:,ISLIP+1), SIG_LAT, &
                & N, M, NUMIND, INDICES)
            RSS(ISLIP,:,INDICES)=RSS(ISLIP,:,INDICES)/CRSS(ISLIP,:,INDICES)
            !
            WHERE (ABS(RSS(ISLIP,:,INDICES)) .LT. T_MIN(IPHASE))
                !
                RSS(ISLIP,:,INDICES) = 0.0D0
                !
            ENDWHERE
            !
            IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .FALSE.) THEN
                !
                CALL POWER_LAW(SHEAR(ISLIP, :, :), RSS(ISLIP, :, :),&
                    & CRYSTAL_PARM(0,IPHASE), CRYSTAL_PARM(1,IPHASE), &
                    & T_MIN(IPHASE), N, M, NUMIND, INDICES)
                !
            ELSE IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .TRUE.) THEN
                !
                IF (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 3) THEN
                    !
                    ANISO_M_TEMP(0:2)  = CRYS_OPTIONS%ANISO_M(IPHASE,1)
                    ANISO_M_TEMP(3:5)  = CRYS_OPTIONS%ANISO_M(IPHASE,2)
                    ANISO_M_TEMP(6:17) = CRYS_OPTIONS%ANISO_M(IPHASE,3)
                    !
                    CALL POWER_LAW(SHEAR(ISLIP, :, :), RSS(ISLIP, :, :), &
                        & ANISO_M_TEMP(ISLIP), CRYSTAL_PARM(1,IPHASE), &
                        & T_MIN(IPHASE), N, M, NUMIND, INDICES)
                    !
                ELSE IF (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 4) &
                    & THEN
                    !
                    ANISO_M_TEMP(0:1)  = CRYS_OPTIONS%ANISO_M(IPHASE, 1)
                    ANISO_M_TEMP(2:3)  = CRYS_OPTIONS%ANISO_M(IPHASE, 2)
                    ANISO_M_TEMP(4:5) = CRYS_OPTIONS%ANISO_M(IPHASE, 3)
                    ANISO_M_TEMP(6:9)  = CRYS_OPTIONS%ANISO_M(IPHASE, 4)
                    ANISO_M_TEMP(10:11)  = CRYS_OPTIONS%ANISO_M(IPHASE, 5)
                    ANISO_M_TEMP(12:15) = CRYS_OPTIONS%ANISO_M(IPHASE, 6)
                    ANISO_M_TEMP(16:17)  = CRYS_OPTIONS%ANISO_M(IPHASE, 7)
                    ANISO_M_TEMP(18:19)  = CRYS_OPTIONS%ANISO_M(IPHASE, 8)
                    ANISO_M_TEMP(20:23) = CRYS_OPTIONS%ANISO_M(IPHASE, 9)
                    ANISO_M_TEMP(24:31)  = CRYS_OPTIONS%ANISO_M(IPHASE, 10)
                    !
                    CALL POWER_LAW(SHEAR(ISLIP, :, :), RSS(ISLIP, :, :), &
                        & ANISO_M_TEMP(ISLIP), CRYSTAL_PARM(1,IPHASE), &
                        & T_MIN(IPHASE), N, M, NUMIND, INDICES)
                    !
                END IF
                !
            ENDIF
            !
        ENDDO !N_SLIP
        !
        GAMMADOT(0:N_SLIP-1,:,EL_SUB1+INDICES)=SHEAR(0:N_SLIP-1,:,INDICES)
        !
        WHERE(ABS(GAMMADOT(0:N_SLIP-1,:,EL_SUB1+INDICES)) .LE. 1D-30)
            !
            GAMMADOT(0:N_SLIP-1,:,EL_SUB1+INDICES) = 0.0D0
            !
        ENDWHERE
        !
        ! Compute net shearing rate on all slip systems
        !
        DO IS = 0, N_SLIP - 1
            !
            SHRATE(:, INDICES) = SHRATE(:, INDICES) + ABS(SHEAR(IS, :, INDICES))
            !
        ENDDO
        !
        ! Compute plastic spin
        !
        DO IS = 0, N_SLIP - 1
            !
            DO I = 0, DIMS1
                !
                WP_SS(I, :, INDICES) = WP_SS(I, :, INDICES) + &
                    & QLOCAL(I+1, IS+1) * SHEAR(IS, :, INDICES)
                !
            enddo
            !
        enddo !N_SLIP
        !
        DEALLOCATE(INDICES)
        DEALLOCATE(PLOCAL)
        DEALLOCATE(QLOCAL)
        !
    ENDDO !NUMPHASES
    !
    SHR_MIN = 1.0D-6 * EPSEFF
    SHR_MAX = 1.0D1 * EPSEFF
    !
    WHERE (SHRATE .LE. SHR_MIN) SHRATE = SHR_MIN
    WHERE (SHRATE .GE. SHR_MAX) SHRATE = SHR_MAX
    !
    RETURN
    !
    END SUBROUTINE SLIP_QNT
    !
    !===========================================================================
    !
    SUBROUTINE UPD_HARD(CRSS, CRSS_0, SHEAR, SHRATE, EPSEFF, DTIME, &
        & IFLAG, DONE, N, M )
    !
    ! Add description here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    ! CRSS: Update hardnesses
    ! CRSS_0: Hardnesses before timestep
    ! SHRATE: Shearing rates
    ! DTIME: Timestep
    ! IFLAG: Indicates whether to use forward or backward Euler
    ! DONE: Mask of elements not to process
    ! N, M: Number of grains and elements
    !
    INTEGER,  INTENT(IN) :: IFLAG, N, M
    REAL(RK), INTENT(IN) :: DTIME
    REAL(RK), INTENT(INOUT) :: CRSS(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)    :: CRSS_0(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)    :: SHRATE(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)    :: SHEAR(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)    :: EPSEFF(0:(N - 1), 0:(M - 1))
    LOGICAL,  INTENT(IN)    :: DONE(0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    INTEGER, POINTER :: INDICES(:) => NULL()
    INTEGER  :: ITER_HARD, NM, INEWTON, IHARD, IPHASE, ISLIP
    INTEGER  :: NUMIND, N_SLIP, I
    !
    LOGICAL  :: NEWTON_OK(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    LOGICAL  :: NEWTON_OK_ALL(0:(N - 1), 0:(M - 1))
    LOGICAL  :: CONVERGED(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    LOGICAL  :: CONVERGED_ALL(0:(N - 1), 0:(M - 1))
    !
    REAL(RK) :: HARD_RATE(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: DHARD_RATE(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: CRSS_SAT(0:(N - 1), 0:(M - 1))
    REAL(RK) :: CRSS_TMP(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: DEL_CRSS(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: RES_N(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: RES(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: RATIO_RES(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: FJ31(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: FJAC(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: XLAM(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    INTEGER  :: MY_PHASE(0:(N - 1), 0:(M - 1))
    !
    !---------------------------------------------------------------------------
    !
    ! Reference list of certain useful crystal parameters:
    ! h_0         = crystal_parm(2)
    ! g_0         = crystal_parm(3)
    ! g_s0        = crystal_parm(4)
    ! m_prime     = crystal_parm(5)
    ! gammadot_s0 = crystal_parm(6)
    !
    MY_PHASE(N-1,:) = PHASE(EL_SUB1:EL_SUP1)
    !
    MAX_ITER_HARD = OPTIONS%MAX_ITER_HARD_LIMIT
    !
    ! Hardwire the legacy flag passed to `hard_law' - is this still used?
    !
    IHARD = 2
    !
    DO IPHASE = 1, NUMPHASES
        !
        IF (OPTIONS%SAT_EVO .EQV. .TRUE.) THEN
            !
            WHERE (MY_PHASE .EQ. IPHASE)
                !
                WHERE(SHRATE .GT. 0.0D0)
                    !
                    CRSS_SAT = CRYSTAL_PARM(4, IPHASE) &
                        & * ((SHRATE / CRYSTAL_PARM(6, IPHASE)) &
                        & ** CRYSTAL_PARM(5, IPHASE))
                    !
                    ELSEWHERE
                    !
                    CRSS_SAT = CRYSTAL_PARM(4, IPHASE)
                    !
                END WHERE
                !
            END WHERE
            !
        ELSE IF (OPTIONS%SAT_EVO .EQV. .FALSE.) THEN
            !
            WHERE (MY_PHASE .EQ. IPHASE)
                !
                CRSS_SAT = CRYSTAL_PARM(4, IPHASE)
                !
            END WHERE
            !
        END IF
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE(N-1,:), INDICES)
        !
        IF (CRYSTAL_PARM(2,IPHASE) .EQ. 0.D0) THEN
            !
            CRSS(:,:,INDICES) = CRSS_0(:,:,INDICES)
            !
        ENDIF
        !
        DEALLOCATE(INDICES)
        !
    ENDDO !NUMPHASES
    !
    IF (IFLAG == UPD_EULER_FWD) THEN
        !
        ! Initial guess via Forward Euler
        !
        CALL HARD_LAW(HARD_RATE, DHARD_RATE, CRSS_0, CRSS_SAT, SHEAR, &
            & SHRATE, EPSEFF, 1, IHARD, N, M)
        !
        CRSS = CRSS_0 + DTIME * HARD_RATE
        !
        RETURN
        !
    END IF
    !
    ! Start Newton iteration - Backward Euler approx. to hardening law
    !
    ! The newton_ok and converged variables are for each individual slip system.
    ! The newton_ok_all and converged_all variables are for each individual
    ! element. It was done this way to reduce the number of changes needed to
    ! be made to the legacy code.
    !
    NEWTON_OK = .TRUE.
    NEWTON_OK_ALL= .TRUE.
    CONVERGED = .TRUE.
    CONVERGED_ALL= .TRUE.
    !
    ! Initializing all the converged indices slip systems for each individual
    ! phase to false. In a multi-phase system where the different phases have
    ! different number of slip systems. This allows the following case where
    ! line 1 is a system with 3 slip systems and line 2 only has one for the
    ! converged variable:
    !
    ! Line 1: FFF
    ! Line 2: FTT
    !
    ! So, in the convergence variable, only the indices that are within the
    ! range of that phases number of slip systems are false.
    !
    DO IPHASE = 1, NUMPHASES
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE))
        N_SLIP=CTYPE(IPHASE)%NUMSLIP
        !
        DO ISLIP=0,N_SLIP-1
            !
            WHERE (MY_PHASE .EQ. IPHASE) CONVERGED(ISLIP,:,:) = .FALSE.
            !
        ENDDO
        !
    ENDDO
    !
    ! Finding where done is true across the element setting converged
    ! true for that entire element
    !
    DO ISLIP = 0, MAXSLIP1
        !
        WHERE(DONE) CONVERGED(ISLIP,:,:) = .TRUE.
        !
    ENDDO
    !
    NM = N * M
    !
    DO ITER_HARD = 1, MAX_ITER_HARD
        !
        ! Should we add "where (.not. converged) clause here?
        ! Yes, but we need to pass mask argument to hard_law
        !
        CRSS_TMP = CRSS
        !
        CALL HARD_LAW(HARD_RATE, DHARD_RATE, CRSS_TMP, CRSS_SAT, SHEAR, &
            & SHRATE, EPSEFF, 1, IHARD, N, M)
        !
        RES_N = - (CRSS_TMP - CRSS_0 - DTIME * HARD_RATE)
        !
        ! Patch to avoid divide-by-zero exception
        !
        WHERE (ABS(RES_N) == 0.0D0)
            !
            CONVERGED = .TRUE.
            !
        ENDWHERE
        !
        CALL HARD_LAW(HARD_RATE, DHARD_RATE, CRSS_TMP, CRSS_SAT, SHEAR, &
            & SHRATE, EPSEFF, 2, IHARD, N, M)
        !
        FJAC = 1.0 - DTIME * DHARD_RATE
        DEL_CRSS = RES_N / FJAC
        !
        ! Find where the newton_ok criterion is not satisfied
        !
        WHERE ((FJAC .LT. VTINY) .AND. (.NOT. CONVERGED))  NEWTON_OK = .FALSE.
        !
        XLAM = CRSS_TMP + DEL_CRSS
        !
        CALL HARD_LAW(HARD_RATE, DHARD_RATE, XLAM, CRSS_SAT, SHEAR, &
            & SHRATE, EPSEFF, 1, IHARD, N, M)
        !
        RES = - (XLAM - CRSS_0 - DTIME * HARD_RATE)
        !
        WHERE ((.NOT. CONVERGED))
            !
            RATIO_RES = ABS(RES) / ABS(RES_N)
            !
        ELSEWHERE
            !
            RATIO_RES = 0.0D0
            !
        ENDWHERE
        !
        ! Line search
        !
        FJ31 = 1.0D0
        !
        DO WHILE (ANY(RATIO_RES .GT. 1.0D0 &
            & .AND. NEWTON_OK .AND. .NOT. CONVERGED))
            !
            WHERE (RATIO_RES .GT. 1.0D0 &
                & .AND. NEWTON_OK .AND. .NOT. CONVERGED) FJ31 = FJ31 * 0.5D0
            !
            ! Find where the newton_ok criterion is not satisfied
            !
            WHERE(FJ31 .LT. LS_CUTOFF) NEWTON_OK = .FALSE.
            !
            XLAM = CRSS_TMP + FJ31 * DEL_CRSS
            !
            CALL HARD_LAW(HARD_RATE, DHARD_RATE, XLAM, CRSS_SAT, SHEAR, &
                & SHRATE, EPSEFF, 1, IHARD, N, M)
            !
            WHERE (RATIO_RES .GT. 1.0D0 .AND. NEWTON_OK .AND. .NOT. CONVERGED)
                !
                RES = - (XLAM - CRSS_0 - DTIME * HARD_RATE)
                RATIO_RES = ABS(RES) / ABS(RES_N)
                !
            ENDWHERE
            !
        ENDDO
        !
        WHERE((NEWTON_OK) .AND. (.NOT.CONVERGED)) &
            & CRSS = CRSS_TMP + FJ31 * DEL_CRSS
        !
        ! Find slip systems that have reached convergence tolerance
        !
        DO IPHASE = 1, NUMPHASES
            !
            CALL CRYSTALTYPEGET(CTYPE(IPHASE))
            N_SLIP=CTYPE(IPHASE)%NUMSLIP
            !
            DO ISLIP = 0, N_SLIP-1
                !
                WHERE (MY_PHASE .EQ. IPHASE)
                    !
                    WHERE ( (DABS(RES(ISLIP,:,:)) .LT. &
                        & TOLER_HARD * CRYSTAL_PARM(3,IPHASE)) .AND. &
                        & (NEWTON_OK(ISLIP,:,:)) .AND. (.NOT. DONE) ) &
                        & CONVERGED(ISLIP,:,:) = .TRUE.
                    !
                ENDWHERE
                !
            ENDDO
            !
        ENDDO
        !
        ! Finds out which areas have not converged and which
        !   newton_ok are not okay
        !
        DO IPHASE = 1, NUMPHASES
            !
            CALL CRYSTALTYPEGET(CTYPE(IPHASE))
            N_SLIP=CTYPE(IPHASE)%NUMSLIP
            !
            DO ISLIP = 0, N_SLIP-1
                !
                ! Converged_all set to false only when one slip system is false
                !
                WHERE(.NOT. CONVERGED(ISLIP,:,:)) CONVERGED_ALL=.FALSE.
                !
                ! Newton_ok_all set to false only when one slip system is false
                !
                WHERE(.NOT. NEWTON_OK(ISLIP,:,:)) NEWTON_OK_ALL=.FALSE.
                !
            ENDDO
            !
            ! If all slip systems for a given element are converged, then
            ! update the per-element arrays as needed. We limit the first index
            ! to the slip systems available for the given phase as above.
            !
            DO I = 0, M - 1
                !
                IF (ALL(CONVERGED(0:N_SLIP-1,:,I))) CONVERGED_ALL(:,I) = .TRUE.
                IF (ALL(NEWTON_OK(0:N_SLIP-1,:,I))) NEWTON_OK_ALL(:,I) = .TRUE.
                !
            END DO
            !
        ENDDO
        !
        INEWTON = COUNT(.NOT. NEWTON_OK_ALL)
        IF ((COUNT(CONVERGED_ALL) + INEWTON) .EQ. NM) EXIT
        !
    ENDDO
    !
    IF (ITER_HARD > MAX_ITER_HARD) THEN
        !
        ! All processes with this condition will write this message. This could
        ! likely be cleaned up to avoid
        !
        WRITE(DFLT_U, '(A)') 'Warning:       . &
            &Update hardening iteration limit reached.'
        !
    END IF
    !
    ! Perhaps a count of "not converged" is a better approach?
    WHERE (.NOT. CONVERGED) CRSS = CRSS_0
    !
    RETURN
    !
    END SUBROUTINE UPD_HARD
    !
    !===========================================================================
    !
    SUBROUTINE UPD_ROTN(RSTAR, RSTAR_N, DR, TIME_STEP, WP_HAT, WP_SS, DONE, &
        & N, M)
    !
    ! Add description here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER :: N, M
    REAL(RK) :: TIME_STEP
    REAL(RK), INTENT(INOUT) :: RSTAR(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)    :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(INOUT) :: DR(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)    :: WP_HAT(0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)    :: WP_SS(0:DIMS1, 0:(N - 1), 0:(M - 1))
    LOGICAL, INTENT(IN)     :: DONE(0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    INTEGER  :: I, J, K
    REAL(RK) :: TH1(0:(N - 1), 0:(M - 1)), TH2(0:(N - 1), 0:(M - 1))
    REAL(RK) :: TH3(0:(N - 1), 0:(M - 1))
    REAL(RK) :: TAU(0:(N - 1), 0:(M - 1)), TAUA(0:(N - 1), 0:(M - 1))
    REAL(RK) :: R(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    !
    !---------------------------------------------------------------------------
    !
    ! Integrate the evolution equation for r*.
    !
    TH1 = (WP_HAT(0, :, :) - WP_SS(0, :, :)) * TIME_STEP
    TH2 = (WP_HAT(1, :, :) - WP_SS(1, :, :)) * TIME_STEP
    TH3 = (WP_HAT(2, :, :) - WP_SS(2, :, :)) * TIME_STEP
    !
    TAU = SQRT(TH1 * TH1 + TH2 * TH2 + TH3 * TH3)
    !
    WHERE (TAU .NE. 0.0)
        !
        TAUA = TAN(TAU / 2.0) / TAU
        TH1 = TAUA * TH1
        TH2 = TAUA * TH2
        TH3 = TAUA * TH3
        TAU = TAUA * TAU
        !
    ENDWHERE
    !
    TAU = 2.0 / (1.0 + TAU * TAU)
    !
    WHERE (.NOT. DONE)
        !
        DR(0, 0, :, :) = 1.0 - TAU * (TH1 * TH1 + TH2 * TH2)
        DR(0, 1, :, :) = - TAU * (TH1 + TH2 * TH3)
        DR(0, 2, :, :) = TAU * ( - TH2 + TH1 * TH3)
        DR(1, 0, :, :) = TAU * (TH1 - TH2 * TH3)
        DR(1, 1, :, :) = 1.0 - TAU * (TH1 * TH1 + TH3 * TH3)
        DR(1, 2, :, :) = - TAU * (TH3 + TH1 * TH2)
        DR(2, 0, :, :) = TAU * (TH2 + TH1 * TH3)
        DR(2, 1, :, :) = TAU * (TH3 - TH1 * TH2)
        DR(2, 2, :, :) = 1.0 - TAU * (TH2 * TH2 + TH3 * TH3)
        !
    ENDWHERE
    !
    R = 0.0D0
    !
    ! RC 3/24/2016: Reordered for better memory striding
    DO J = 0, DIMS1
        !
        DO K = 0, DIMS1
            !
            DO I = 0, DIMS1
                !
                R(I, J, :, :) = R(I, J, :, :) + &
                    & RSTAR_N(I, K, :, :) * DR(K, J, :, :)
                !
            ENDDO
            !
        ENDDO
        !
    ENDDO
    !
    WHERE (.NOT. DONE)
        !
        RSTAR(0, 0, :, :) = R(0, 0, :, :)
        RSTAR(0, 1, :, :) = R(0, 1, :, :)
        RSTAR(0, 2, :, :) = R(0, 2, :, :)
        RSTAR(1, 0, :, :) = R(1, 0, :, :)
        RSTAR(1, 1, :, :) = R(1, 1, :, :)
        RSTAR(1, 2, :, :) = R(1, 2, :, :)
        RSTAR(2, 0, :, :) = R(2, 0, :, :)
        RSTAR(2, 1, :, :) = R(2, 1, :, :)
        RSTAR(2, 2, :, :) = R(2, 2, :, :)
        !
    ENDWHERE
    !
    RETURN
    !
    END SUBROUTINE UPD_ROTN
    !
    !===========================================================================
    !
    SUBROUTINE GET_C(C, R, C_0, DONE, N, M)
    !
    ! Add description here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER :: N, M
    REAL(RK), INTENT(IN)  :: R(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M -1))
    REAL(RK), INTENT(IN)  :: C_0(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M -1))
    REAL(RK), INTENT(OUT) :: C(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M -1))
    LOGICAL, INTENT(IN)   :: DONE(0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    INTEGER  :: I, J, K
    REAL(RK) :: C_AUX(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M -1))
    !
    !---------------------------------------------------------------------------
    !
    ! Compute : c = c_0 * rstar
    !
    C_AUX = 0.0D0
    !
    ! RC 3/24/2016: Reordered for better memory striding
    DO J = 0, DIMS1
        !
        DO K = 0, DIMS1
            !
            DO I = 0, DIMS1
                !
                C_AUX(I, J, :, :) = C_AUX(I, J, :, :) + &
                    & C_0(I, K, :, :) * R(K, J, :, :)
                !
            ENDDO
            !
        ENDDO
        !
    ENDDO
    !
    WHERE (.NOT. DONE)
        !
        C(0, 0, :, :) = C_AUX(0, 0, :, :)
        C(0, 1, :, :) = C_AUX(0, 1, :, :)
        C(0, 2, :, :) = C_AUX(0, 2, :, :)
        C(1, 0, :, :) = C_AUX(1, 0, :, :)
        C(1, 1, :, :) = C_AUX(1, 1, :, :)
        C(1, 2, :, :) = C_AUX(1, 2, :, :)
        C(2, 0, :, :) = C_AUX(2, 0, :, :)
        C(2, 1, :, :) = C_AUX(2, 1, :, :)
        C(2, 2, :, :) = C_AUX(2, 2, :, :)
        !
    ENDWHERE
    !
    RETURN
    !
    END SUBROUTINE GET_C
    !
END MODULE RSTARN_SOLVE_MOD
