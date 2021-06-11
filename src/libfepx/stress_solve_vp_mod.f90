! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE STRESS_SOLVE_VP_MOD
!
! Routines for solving viscoplastic crystal stress equations
!
! Contains subroutines:
! STRESS_SOLVE_VP: Routine which takes care of scaling and initial guesses
! SCALE_DOWN_DEFR: Rescale the deformation rate to unit size
! COMPUTE_WORK: Compute virtual plastic work for array of deformation rates
!   applied to vertex stresses
! FIND_VERTEX: Select the vertex which maximizes the plastic work
! VERTEX_STRESS: Set initial guess (vertex stress) for nonlinear solver
! SCALE_STRESS: Rescale stress initial guess according to hardness
! SS_PROJECT: Compute inner product of array of tensors with a fixed tensor
! SOLVE_NEWTON_VP: Nonlinear solution of viscoplastic crystal stress equations
! GET_RES: Compute residual for nonlinear VP crystal stress equation
! FORM_CRYSTIF: Form single crystal stiffness matrix
! POWER_LAW: Power law for VP single crystal
! COMPLIANCE: Form crystal compliance matrices
! SOLVIT: Solve an array of symmetric positive definite 5X5 systems
! CHECK_DIAGONALS: Determine where diagonal elements are small
! SCALE_UP_SIGM: Rescale the stress after solution is found
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
USE UNITS_MOD
!
! From libparallel:
!
USE PARALLEL_MOD
!
IMPLICIT  NONE
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: STRESS_SOLVE_VP
PUBLIC :: SS_PROJECT
PUBLIC :: POWER_LAW
PUBLIC :: COMPLIANCE
PUBLIC :: SOLVIT
PUBLIC :: CHECK_DIAGONALS
!
CONTAINS
    !
    SUBROUTINE STRESS_SOLVE_VP(SIG, D_VEC_LAT, CRSS, EPSEFF, VP_LOG)
    !
    ! Driver routine which takes care of scaling and initial guesses
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! SIG:
    ! D_VEC_LAT:
    ! CRSS:
    ! EPSEFF:
    ! VP_LOG:
    !
    REAL(RK), INTENT(OUT) :: SIG(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: D_VEC_LAT(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: CRSS  (0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: EPSEFF(0:NGRAIN1, EL_SUB1:EL_SUP1)
    LOGICAL,  INTENT(IN) :: VP_LOG
    !
    ! Locals:
    !
    LOGICAL :: CONVERGED(0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER :: I_EDGE
    INTEGER :: IER
    INTEGER :: M_EL
    INTEGER :: NTRIALS
    INTEGER :: VERTEX(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DIRECTION(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: PLWORK(0:MAX_VERT1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG_T(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: CRSS_AVG (0:NGRAIN1, EL_SUB1:EL_SUP1)
    !
    !---------------------------------------------------------------------------
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    CRSS_AVG = 0.0D0
    !
    CALL SCALE_DOWN_DEFR(D_VEC_LAT, EPSEFF, NGRAIN, M_EL)
    CALL COMPUTE_WORK(PLWORK, D_VEC_LAT, NGRAIN, M_EL)
    CALL COMPUTE_AVG_CRSS(CRSS,CRSS_AVG,M_EL)
    !
    CONVERGED = .FALSE.
    !
    ! This loop used to be from 1 to N_EDGE, in order to try all vertices as
    !   initial guesses. However, this can be very time-consuming and doesn't
    !   usually help.
    !
    NTRIALS = 4
    !
    DO I_EDGE = 1, NTRIALS
        !
        CALL FIND_VERTEX(VERTEX, DIRECTION, PLWORK,  NGRAIN, M_EL)
        !
        CALL VERTEX_STRESS(SIG_T, VERTEX, DIRECTION, NGRAIN, M_EL)
        !
        CALL SCALE_STRESS(SIG_T, CRSS_AVG, NGRAIN, M_EL)
        !
        WHERE (.NOT. CONVERGED)
            !
            SIG(0, :, :) = SIG_T(0, :, :)
            SIG(1, :, :) = SIG_T(1, :, :)
            SIG(2, :, :) = SIG_T(2, :, :)
            SIG(3, :, :) = SIG_T(3, :, :)
            SIG(4, :, :) = SIG_T(4, :, :)
            !
        END WHERE
        !
        CALL SOLVE_NEWTON_VP(SIG, D_VEC_LAT, CRSS_AVG, IER, CV_OPTIONS%SX_TOL, &
            & CONVERGED, NGRAIN, M_EL, VP_LOG)
        !
        IF (IER .EQ. 0) GO TO 50
        !
    END DO
    !
    IF ((MYID .EQ. 0) .AND. VP_LOG) THEN
        !
        WRITE(DFLT_U, '(A)') 'Warning:       . STRESS_SOLVE_VP: '
        WRITE(DFLT_U, '(A)') 'Warning:       . ',COUNT(.NOT. converged),'&
            & elements did not converge'
        WRITE(DFLT_U, '(A)') 'Warning:       . after ', NTRIALS,' trials.'
        !
    END IF
    !
    WHERE (.NOT. CONVERGED)
        !
        SIG(0, :, :) = SIG_T(0, :, :)
        SIG(1, :, :) = SIG_T(1, :, :)
        SIG(2, :, :) = SIG_T(2, :, :)
        SIG(3, :, :) = SIG_T(3, :, :)
        SIG(4, :, :) = SIG_T(4, :, :)
        !
    END WHERE
    !
    50 CONTINUE
    !
    CALL SCALE_UP_SIGM(SIG, EPSEFF, NGRAIN, M_EL)
    !
    END SUBROUTINE STRESS_SOLVE_VP
    !
    !===========================================================================
    !
    SUBROUTINE SCALE_DOWN_DEFR(D_VEC, EPSEFF, N, M)
    !
    ! Rescale the deformation rate to unit size.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! D_VEC: Array of deformation rates as 5-vectors (input/output)
    ! EPSEFF: Effective deformation of these tensors (input)
    ! N: Number of grains
    ! M: Number of elements
    !
    REAL(RK), INTENT(INOUT) :: D_VEC(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: EPSEFF(0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    !
    INTEGER :: I
    !
    !---------------------------------------------------------------------------
    !
    !-tm from donald's verion:
    !
    DO I = 0, TVEC1
        !
        WHERE (EPSEFF > 0.0D0)
            !
            D_VEC(I, :, :) = D_VEC(I, :, :) / EPSEFF
            !
        END WHERE
        !
    END DO
    !
    END SUBROUTINE SCALE_DOWN_DEFR
    !
    !===========================================================================
    !
    SUBROUTINE COMPUTE_WORK(PLWORK, D_VEC, N, M)
    !
    ! Compute virtual plastic work for array of deformation rates applied to
    !   vertex stresses.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! PLWORK:
    ! D_VEC:
    ! N:
    ! M:
    !
    REAL(RK), INTENT(OUT) :: PLWORK(0:MAX_VERT1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: D_VEC(0:TVEC1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    !
    REAL(RK), POINTER :: SIG_FS(:,:)=>NULL()
    INTEGER :: MY_PHASE(0:(N - 1), 0:(M - 1))
    INTEGER :: N_EDGE
    INTEGER :: I
    INTEGER :: J
    INTEGER :: IPHASE
    !
    !---------------------------------------------------------------------------
    !
    ! tsh, 1/26/03
    MY_PHASE(N - 1, :) = PHASE(EL_SUB1:EL_SUP1)
    !
    PLWORK = 0.0D0
    !
    DO IPHASE = 1, NUMPHASES
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE), VERTICES = SIG_FS)
        !
        N_EDGE = CTYPE(IPHASE)%NUMVERTICES
        SIG_FS = SIG_FS * CRYSTAL_PARM(9,IPHASE)
        !
        DO I = 0, (N_EDGE - 1)
            !
            DO J = 0, TVEC1
                !
                WHERE (MY_PHASE .EQ. IPHASE)
                    !
                    PLWORK(I, :, :) = PLWORK(I, :, :) + SIG_FS(J + 1, I + 1) &
                        & * D_VEC(J, :, :)
                    !
                END WHERE
                !
            END DO
            !
        END DO
        !
        DEALLOCATE(SIG_FS)
        !
    END DO !NUMPHASES
    !
    END SUBROUTINE COMPUTE_WORK
    !
    !===========================================================================
    !
    SUBROUTINE FIND_VERTEX(VERTEX, DIRECTION, PLWORK, N, M)
    !
    ! Select the vertex which maximizes the plastic work.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! VERTEX: List of optimal vertex numbers for each grain (output)
    ! DIRECTION: Sign of the vertex (1 or -1) (output)
    ! PLWORK: Array of virtual plastic work for each grain and all vertices (in)
    ! N: Number of grains
    ! M: Number of elements
    !
    INTEGER, INTENT(OUT) :: VERTEX(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: DIRECTION(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(INOUT) :: PLWORK(0:MAX_VERT1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    !
    INTEGER, POINTER :: INDICES(:)=>NULL()
    REAL(RK) :: PA(0:(N - 1), 0:(M - 1))
    INTEGER :: I
    INTEGER :: J
    INTEGER :: IPHASE
    INTEGER :: NUMIND
    INTEGER :: N_EDGE
    INTEGER  :: MY_PHASE(0:(N - 1),0:(M - 1))
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(N - 1,:) = PHASE(EL_SUB1:EL_SUP1)
    !
    VERTEX = 0
    DIRECTION = -1.0D0
    DO IPHASE = 1, NUMPHASES
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE(N - 1, :), INDICES)
        !
        N_EDGE = CTYPE(IPHASE)%NUMVERTICES
        DO I = 0, (N_EDGE - 1)
            !
            PA(:,INDICES) = ABS(PLWORK(I, :, INDICES))
            !
            WHERE (PA .GT. DIRECTION)
                !
                DIRECTION = PA
                VERTEX = I
                !
            END WHERE
            !
        END DO
        !
        DEALLOCATE(INDICES)
        !
    END DO !NUMPHASES
    !
    ! RC 6/24/2016: Reordered loops for better memory striding
    !
    DO J = 0,(M - 1)
        !
        DO I = 0,(N - 1)
            !
            DIRECTION(I, J) = PLWORK(VERTEX(I, J), I, J) / DIRECTION(I, J)
            !
        END DO
        !
    END DO
    !
    ! To keep from the vertex being reused
    !
    !  RC 6/24/2016: Reordered loops for better memory striding
    DO J = 0,(M - 1)
        !
        DO I = 0,(N - 1)
            !
            PLWORK(VERTEX(I, J), I, J) = 0.0D0
            !
        END DO
        !
    END DO
    !
    END SUBROUTINE FIND_VERTEX
    !
    !===========================================================================
    !
    SUBROUTINE VERTEX_STRESS(SIG, VERTEX, DIRECTION, N, M)
    !
    !     Set initial guess (vertex stress) for nonlinear solver.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! SIG: Initial guesses for stresses (output)
    ! SIG_FS: Vertex stresses (input)
    ! VERTEX: Optimal vertex numbers
    ! DIRECTION: Sign to multiply vertex stress by
    ! N: Number of grains
    ! M: Number of elements
    !
    REAL(RK), INTENT(OUT) :: SIG(0:TVEC1, 0:(N - 1), 0:(M - 1))
    INTEGER,INTENT(IN) :: VERTEX(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: DIRECTION(0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    !
    INTEGER :: MY_PHASE(0:(N - 1),0:(M - 1))
    REAL(RK), POINTER :: SIG_FS(:,:) => NULL()
    INTEGER :: N_EDGE
    INTEGER :: I
    INTEGER :: J
    INTEGER :: IPHASE
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(N - 1,:) = PHASE(EL_SUB1:EL_SUP1)
    !
    DO IPHASE = 1, NUMPHASES
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE), VERTICES = SIG_FS)
        !
        N_EDGE = CTYPE(IPHASE)%NUMVERTICES
        ! RC 6/24/2016 reordered the WHERE loops as well so it only runs the
        !   parts needed
        !
        SIG_FS = SIG_FS * CRYSTAL_PARM(9, IPHASE)
        !
        ! RC 6/24/2016: Reordered loops for better memory striding
        DO J = 0, (N_EDGE - 1)
            !
            DO I = 0, TVEC1
                !
                WHERE (J .EQ. VERTEX)
                    !
                    WHERE (MY_PHASE .EQ. IPHASE)
                        !
                        SIG(I, :, :) = SIG_FS(I + 1, J + 1) * DIRECTION(:, :)
                        !
                    END WHERE
                    !
                END WHERE
                !
            END DO
            !
        END DO
        !
        DEALLOCATE(SIG_FS)
        !
    END DO !NUMPHASES
    !
    END SUBROUTINE VERTEX_STRESS
    !
    !===========================================================================
    !
    SUBROUTINE SCALE_STRESS(SIG, CRSS, N, M)
    !
    ! Rescale stress initial guess according to hardness.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! SIG: Initial guess for stress (input/output)
    ! CRSS: Crystal hardnesses
    ! N: Number of grains
    ! M: Number of elements
    !
    REAL(RK), INTENT(INOUT) :: SIG(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: CRSS(0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    !
    ! Locals:
    !
    INTEGER :: MY_PHASE(0:(N - 1), 0:(M - 1))
    INTEGER :: ISLIP
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    INTEGER :: IPHASE
    INTEGER :: NUMIND
    INTEGER :: N_SLIP
    INTEGER, POINTER :: INDICES(:)=>NULL()
    REAL(RK), POINTER  :: P_HAT_VEC(:, :)=>NULL()
    REAL(RK) :: TAUMAX(0:(N - 1), 0:(M - 1))
    REAL(RK) :: TAU(0:(N - 1), 0:(M - 1))
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(N - 1,:) = PHASE(EL_SUB1:EL_SUP1)
    !
    !-tm  from donald's verion:
    !
    TAUMAX = 1.0D-8  ! deb/ prevent div by 0
    !TAUMAX = 0.0d0
    !
    DO IPHASE = 1, NUMPHASES
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE), DEV = P_HAT_VEC)
        !
        N_SLIP = CTYPE(IPHASE)%NUMSLIP
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE(N - 1, :), INDICES)
        !
        DO ISLIP = 0, (N_SLIP - 1)
            !
            CALL SS_PROJECT(TAU, P_HAT_VEC(:, ISLIP + 1), SIG, N, M, NUMIND, &
                & INDICES)
            !
            TAU(:,INDICES) = DABS(TAU(:, INDICES))
            !
            ! maybe: WHERE ((TAU(:,INDICES) .GT. TAUMAX) TAUMAX=TAU
            WHERE ((MY_PHASE .EQ. IPHASE) .AND. (TAU .GT. TAUMAX))
                !
                TAUMAX = TAU
                !
            END WHERE
            !
        END DO !N_SLIP
        !
        DEALLOCATE(P_HAT_VEC)
        DEALLOCATE(INDICES)
        !
    END DO !NUMPHASES
    !
    TAU = CRSS / TAUMAX
    !
    DO K = 0, (M - 1)
        !
        DO J = 0, (N - 1)
            !
            DO I = 0, TVEC1
                !
                SIG(I, J, K) = SIG(I, J, K) * TAU(J, K)
                !
            END DO
            !
        END DO
        !
    END DO
    !
    END SUBROUTINE SCALE_STRESS
    !
    !===========================================================================
    !
    SUBROUTINE SS_PROJECT(PROJ, PLOCAL, TENSOR, N, M, NUMIND, INDICES)
    !
    ! Compute inner product of array of tensors with a fixed tensor.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! PROJ:
    ! PLOCAL:
    ! TENSOR:
    ! N:
    ! M:
    ! NUMIND:
    ! INDICES:
    !
    REAL(RK), INTENT(OUT) :: PROJ(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: PLOCAL(0:TVEC1)
    REAL(RK), INTENT(IN) :: TENSOR(0:TVEC1, 0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    INTEGER, INTENT(IN) :: NUMIND
    INTEGER, INTENT(IN) :: INDICES(1:NUMIND)
    !
    ! Locals:
    !
    INTEGER :: I
    REAL(RK) :: PROJ_TMP(0:(N - 1), 0:(NUMIND-1))
    !
    !---------------------------------------------------------------------------
    !
    PROJ_TMP = 0.0D0
    !
    DO I = 0, TVEC1
        !
        PROJ_TMP = PROJ_TMP + PLOCAL(I) * TENSOR(I, :, INDICES)
        !
    END DO
    !
    PROJ(:, INDICES) = PROJ_TMP
    !
    END SUBROUTINE SS_PROJECT
    !
    !===========================================================================
    !
    SUBROUTINE SOLVE_NEWTON_VP(SIG, D_VEC, CRSS, IRC, EPS, CONVERGED, N, M, &
        &VP_LOG)
    !
    ! Nonlinear solution of viscoplastic crystal stress equations.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! SIG: Initial guess for stresses and final solution (input/output)
    ! D_VEC: Deformation rate for which to solve
    ! CRSS: Crystal hardnesses
    ! IRC: Return flag
    ! EPS: Error tolerance for nonlinear solution
    ! CONVERGED: Array telling what crystals have already converged
    ! N: Number of grains
    ! M: Number of elements
    ! VP_LOG: Write viscoplastic convergence output to log files
    !
    REAL(RK), INTENT(INOUT) :: SIG(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: D_VEC(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: CRSS(0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(OUT) :: IRC
    REAL(RK), INTENT(IN) :: EPS
    LOGICAL, INTENT(INOUT) :: CONVERGED(0:(N - 1), 0:(M - 1))
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(IN) :: M
    LOGICAL, INTENT(IN) :: VP_LOG
    !
    ! Locals:
    !
    LOGICAL :: NEWTON_OK(0:(N - 1), 0:(M - 1))
    !
    INTEGER :: ITER
    INTEGER :: NM
    INTEGER :: INEWTON
    !
    REAL(RK) :: RES(0:(N - 1), 0:(M - 1))
    REAL(RK) :: RES_N(0:(N - 1), 0:(M - 1))
    REAL(RK) :: FACT(0:(N - 1), 0:(M - 1))
    REAL(RK) :: RATIO_RES(0:(N - 1), 0:(M - 1))
    REAL(RK) :: RES_AUX(0:(N - 1), 0:(M - 1))
    REAL(RK) :: SIG_0(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: RHS(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: RSS(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: SHEAR(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: STIF(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: DEL_S(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: XLAMBDA(0:TVEC1, 0:(N - 1), 0:(M - 1))
    !
    !---------------------------------------------------------------------------
    !
    NM = N * M
    IRC = 0
    SIG_0 = 0.0D0
    !
    NEWTON_OK = .TRUE.
    !
    DO ITER = 1, CV_OPTIONS%SX_MAX_ITERS_NEWTON
        !
        SIG_0 = SIG
        !
        CALL GET_RES(RES_N, RHS, RSS, SHEAR, SIG, D_VEC, CRSS, N, M)
        !
        CALL FORM_CRYSTIF(stif, RSS, SHEAR, CRSS, N, M)
        !
        DEL_S = RHS
        !
        CALL SOLVIT(STIF, DEL_S, N, M)
        !
        CALL CHECK_DIAGONALS(STIF, NEWTON_OK, N, M)
        !
        XLAMBDA = SIG_0 + DEL_S
        !
        CALL GET_RES(RES, RHS, RSS, SHEAR, XLAMBDA, D_VEC, CRSS, N, M)
        !
        FACT = 1.0D0
        RATIO_RES = RES / RES_N
        !
        DO WHILE (ANY(RATIO_RES .GT. 1.0D0 .AND. NEWTON_OK .AND. .NOT. &
            &CONVERGED))
            !
            WHERE (RATIO_RES .GT. 1.0D0 .AND. NEWTON_OK .AND. .NOT. CONVERGED)
                !
                FACT = FACT*0.5D0
            END WHERE
            !
            IF (ANY(FACT .LT. 0.001D0)) THEN
                !
                IF (VP_LOG .AND. (MYID .EQ. 0)) THEN
                    !
                    WRITE(DFLT_U, '(A)') 'Warning:       . SOLVE_NEWTON_VP: &
                        &Line search failure for ', COUNT(FACT .LT. 0.001D0), &
                        & ' grains.'
                    !
                END IF
                !
                WHERE (FACT .LT. 0.001D0)
                    !
                    NEWTON_OK = .FALSE.
                    !
                END WHERE
                !
            END IF
            !
            XLAMBDA(0, :, :) = SIG_0(0, :, :) + FACT * DEL_S(0, :, :)
            XLAMBDA(1, :, :) = SIG_0(1, :, :) + FACT * DEL_S(1, :, :)
            XLAMBDA(2, :, :) = SIG_0(2, :, :) + FACT * DEL_S(2, :, :)
            XLAMBDA(3, :, :) = SIG_0(3, :, :) + FACT * DEL_S(3, :, :)
            XLAMBDA(4, :, :) = SIG_0(4, :, :) + FACT * DEL_S(4, :, :)
            !
            CALL GET_RES(RES_AUX, RHS, RSS, SHEAR, XLAMBDA, D_VEC, CRSS, N, M)
            !
            WHERE(RATIO_RES .GT. 1.0D0 .AND. NEWTON_OK .AND. .NOT. CONVERGED)
                !
                RES = RES_AUX
                RATIO_RES = RES / RES_N
                !
            END WHERE
            !
        END DO
        !
        WHERE (NEWTON_OK .AND. .NOT. CONVERGED)
            !
            SIG(0, :, :) = SIG_0(0, :, :) + FACT * DEL_S(0, :, :)
            SIG(1, :, :) = SIG_0(1, :, :) + FACT * DEL_S(1, :, :)
            SIG(2, :, :) = SIG_0(2, :, :) + FACT * DEL_S(2, :, :)
            SIG(3, :, :) = SIG_0(3, :, :) + FACT * DEL_S(3, :, :)
            SIG(4, :, :) = SIG_0(4, :, :) + FACT * DEL_S(4, :, :)
            !
        END WHERE
        !
        WHERE (RES .LE. EPS .AND. NEWTON_OK)
            !
            CONVERGED = .TRUE.
            !
        END WHERE
        !
        ! Return if all grains have converged or solution is only an estimate.
        !
        INEWTON = COUNT(.NOT. NEWTON_OK)
        !
        IF (((COUNT(CONVERGED) + INEWTON) .EQ. NM) .AND. (MYID .EQ. 0) ) &
            & THEN
            !
            IF ((INEWTON .GT. 0) .AND. VP_LOG) THEN
                !
                WRITE(DFLT_U,'(A)') 'Info   :     > SOLVE_NEWTON_VP: &
                    &Converged = ', COUNT(CONVERGED), ' remaining = ', &
                    & INEWTON, ' minval res = ', &
                    & MINVAL(RES, MASK = CONVERGED), ' maxval res = ', &
                    & MAXVAL(RES, MASK = CONVERGED)
                !
            END IF
            !
            IRC = INEWTON
            !
            RETURN
            !
        END IF
        !
    END DO !NEWTON_ITERATIONS
    !
    IRC = -2
    !
    END SUBROUTINE SOLVE_NEWTON_VP
    !
    !===========================================================================
    !
    SUBROUTINE GET_RES(RES, RHS, RSS, SHEAR, SIG, D, CRSS, N, M)
    !
    ! Compute residual for nonlinear VP crystal stress equation.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! RES:
    ! RHS:
    ! RSS:
    ! SHEAR:
    ! SIG:
    ! D:
    ! CRSS:
    ! N:
    ! M:
    !
    REAL(RK), INTENT(OUT) :: RES(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: RHS(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: RSS(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: SHEAR(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: SIG(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: D(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: CRSS(0:(N - 1), 0:(M - 1))
    INTEGER :: N
    INTEGER :: M
    !
    !     Locals:
    !
    INTEGER :: MY_PHASE(0:(N - 1), 0:(M - 1))
    INTEGER :: N_SLIP
    INTEGER, POINTER :: INDICES(:)=>NULL()
    REAL(RK), POINTER :: P(:, :)=>NULL()
    INTEGER :: ISLIP
    INTEGER :: J
    INTEGER :: IPHASE
    INTEGER :: NUMIND
    REAL(RK) :: XM_FAKE
    !
    !---------------------------------------------------------------------------
    !
    ! tsh, 1/26/03
    MY_PHASE(N - 1,:) = PHASE(EL_SUB1:EL_SUP1)
    !
    RHS = -D
    RES = 0.0D0
    !
    DO IPHASE = 1, NUMPHASES
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE), DEV = P)
        !
        N_SLIP = CTYPE(IPHASE)%NUMSLIP
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE(N - 1, :), INDICES)
        !
        DO ISLIP = 0, (N_SLIP - 1)
            !
            CALL SS_PROJECT(RSS(ISLIP, :, :), P(:, ISLIP + 1), SIG, N, M, &
                & NUMIND, INDICES)
            !
            RSS(ISLIP, :, INDICES) = RSS(ISLIP, :, INDICES) / CRSS(:, INDICES)
            !
            WHERE (ABS(RSS(ISLIP, :, INDICES)) .LT. T_MIN(IPHASE))
                !
                RSS(ISLIP, :, INDICES) = 0.0D0
                !
            END WHERE
            !
            XM_FAKE = 0.4D0
            !XM_FAKE=0.02d0
            !
            CALL POWER_LAW(SHEAR(ISLIP, :, :), RSS(ISLIP, :, :), XM_FAKE, &
                & CRYSTAL_PARM(1, IPHASE), T_MIN(IPHASE), N, M, NUMIND, INDICES)
            !
            DO J = 0, TVEC1
                !
                RHS(J, :, INDICES) = RHS(J, :, INDICES) + P(J + 1, ISLIP + 1) &
                    & * SHEAR(ISLIP, :, INDICES)
                !
            END DO
            !
        END DO !N_SLIP
        !
        DEALLOCATE(P)
        DEALLOCATE(INDICES)
        !
    END DO !NUMPHASES
    !
    DO J = 0, DIMS1
        !
        RES = RES + RHS(J, :, :) ** 2.0D0
        !
    END DO
    !
    RES = SQRT(RES)
    !
    END SUBROUTINE GET_RES
    !
    !===========================================================================
    !
    SUBROUTINE FORM_CRYSTIF(STIF, RSS, SHEAR, CRSS, N, M)
    !
    ! Form single crystal stiffness matrix.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! STIF:
    ! RSS:
    ! SHEAR:
    ! CRSS:
    ! N:
    ! M:
    !
    REAL(RK), INTENT(OUT) :: STIF(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: RSS(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: SHEAR(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: CRSS(0:(N - 1), 0:(M - 1))
    INTEGER :: N
    INTEGER :: M
    !
    ! Locals:
    !
    INTEGER :: N_SLIP
    REAL(RK), POINTER :: P(:, :)=>NULL()
    INTEGER :: MY_PHASE(0:(N - 1),0:(M - 1))
    INTEGER :: ISLIP
    INTEGER :: J
    INTEGER :: K
    INTEGER :: IPHASE
    INTEGER :: NUMIND
    INTEGER, POINTER :: INDICES(:)=>NULL()
    REAL(RK) :: COMP(0:(N - 1), 0:(M - 1))
    REAL(RK) :: XM_FAKE
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(N - 1, :) = PHASE(EL_SUB1:EL_SUP1)
    !
    STIF = 0.0D0
    !
    DO IPHASE = 1, NUMPHASES
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE), DEV = P)
        !
        N_SLIP = CTYPE(IPHASE)%NUMSLIP
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE(N - 1, :), INDICES)
        !
        DO ISLIP = 0, (N_SLIP - 1)
            !
            XM_FAKE = 0.4D0
            !
            CALL COMPLIANCE(COMP, RSS(ISLIP, :, :), SHEAR(ISLIP, :, :), CRSS, &
                & XM_FAKE, T_MIN(IPHASE), N, M, NUMIND, INDICES)
            !
            DO J = 0, TVEC1
                !
                DO K = 0, TVEC1
                    !
                    STIF(J, K, :, INDICES) = STIF(J, K, :, INDICES) - &
                        & COMP(:,INDICES) * P(J + 1, ISLIP + 1) * &
                        & P(K + 1, ISLIP + 1)
                    !
                END DO
                !
            END DO
            !
        END DO !N_SLIP
        !
        DEALLOCATE(P)
        DEALLOCATE(INDICES)
        !
    END DO !NUMPHASES
    !
    END SUBROUTINE FORM_CRYSTIF
    !
    !===========================================================================
    !
    SUBROUTINE POWER_LAW(POWER, T, XM, A_0, T_MIN, N, M, NUMIND, INDICES)
    !
    ! Power law for VP single crystal.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! POWER:
    ! T:
    ! XM:
    ! A_0:
    ! T_MIN:
    ! N:
    ! M:
    ! NUMIND:
    ! INDICES:
    !
    REAL(RK), INTENT(OUT) :: POWER(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: T(0:(N - 1), 0:(M - 1))
    REAL(RK) :: XM
    REAL(RK) :: A_0
    REAL(RK) :: T_MIN
    INTEGER :: N
    INTEGER :: M
    INTEGER :: NUMIND
    INTEGER :: INDICES(1:NUMIND)
    !
    ! Locals:
    !
    REAL(RK) :: P
    REAL(RK) :: POWER_TMP(0:(N - 1),0:(NUMIND-1))
    REAL(RK) :: AT(0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: ALOG(0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: BLOG(0:(N - 1), 0:(NUMIND - 1))
    !
    !---------------------------------------------------------------------------
    !
    P = 1.0D0 / XM - 1.0D0
    AT = ABS(T(:,INDICES))
    !
    WHERE (AT .GT. T_MIN)
        !
        ALOG = DLOG(AT)
        BLOG = P * ALOG
        POWER_TMP = A_0 * T(:, INDICES) * DEXP(BLOG)
        !
    ELSE WHERE
        !
        POWER_TMP = 0.0D0
        !
    END WHERE
    !
    POWER(:, INDICES) = POWER_TMP
    !
    END SUBROUTINE POWER_LAW
    !
    !===========================================================================
    !
    SUBROUTINE COMPLIANCE(COMP, T, SHEAR, CRSS, XM, T_MIN, N, M, NUMIND, &
        & INDICES)
    !
    ! Form crystal compliance matrices
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! COMP:
    ! T:
    ! SHEAR:
    ! CRSS:
    ! XM:
    ! T_MIN:
    ! N:
    ! M:
    ! NUMIND:
    ! INDICES:
    !
    REAL(RK), INTENT(OUT) :: COMP(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: T(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: SHEAR(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: CRSS(0:(N - 1), 0:(M - 1))
    REAL(RK) :: T_MIN
    INTEGER :: N
    INTEGER :: M
    INTEGER :: NUMIND
    INTEGER :: INDICES(1:NUMIND)
    !
    ! Locals:
    !
    REAL(RK) :: XM
    REAL(RK) :: COMP_TMP(0:(N - 1), 0:(NUMIND-1))
    !
    !---------------------------------------------------------------------------
    !
    COMP_TMP = 0.0D0
    !
    WHERE (ABS(T(:,INDICES)) .GT. T_MIN)
        !
        COMP_TMP = SHEAR(:, INDICES) / (T(:, INDICES) * CRSS(:, INDICES) * XM)
        !
    END WHERE
    !
    COMP(:, INDICES) = COMP_TMP
    !
    END SUBROUTINE COMPLIANCE
    !
    !===========================================================================
    !
    SUBROUTINE SOLVIT(A, X, N, M)
    !
    ! Solve an array of symmetric positive definite 5X5 systems.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! A:
    ! X:
    ! N:
    ! M:
    !
    REAL(RK), INTENT(IN) :: A(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(INOUT) :: X(0:TVEC1, 0:(N - 1), 0:(M - 1))
    INTEGER :: N
    INTEGER :: M
    !
    ! Locals:
    !
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A11
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A21
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A22
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A31
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A32
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A33
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A41
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A42
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A43
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A44
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A51
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A52
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A53
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A54
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: A55
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: X1
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: X2
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: X3
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: X4
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: X5
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: V1
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: V2
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: V3
    REAL(RK), DIMENSION(0:(N - 1), 0:(M - 1)) :: V4
    !
    !---------------------------------------------------------------------------
    !
    A11 = A(0, 0, :, :)
    A21 = A(1, 0, :, :)
    A31 = A(2, 0, :, :)
    A41 = A(3, 0, :, :)
    A51 = A(4, 0, :, :)
    A22 = A(1, 1, :, :)
    A32 = A(2, 1, :, :)
    A42 = A(3, 1, :, :)
    A52 = A(4, 1, :, :)
    A33 = A(2, 2, :, :)
    A43 = A(3, 2, :, :)
    A53 = A(4, 2, :, :)
    A44 = A(3, 3, :, :)
    A54 = A(4, 3, :, :)
    A55 = A(4, 4, :, :)
    X1 = X(0, :, :)
    X2 = X(1, :, :)
    X3 = X(2, :, :)
    X4 = X(3, :, :)
    X5 = X(4, :, :)
    !
    ! **  A = LDL'.
    ! **  j = 1.
    !
    A21 = A21 / A11
    A31 = A31 / A11
    A41 = A41 / A11
    A51 = A51 / A11
    !
    ! **  j = 2.
    !
    V1 = A21 * A11
    A22 = A22 - A21 * V1
    A32 = (A32 - A31 * V1) / A22
    A42 = (A42 - A41 * V1) / A22
    A52 = (A52 - A51 * V1) / A22
    !
    ! **  j = 3.
    !
    V1 = A31 * A11
    V2 = A32 * A22
    A33 = A33 - A31 * V1 - A32 * V2
    A43 = (A43 - A41 * V1 - A42 * V2) / A33
    A53 = (A53 - A51 * V1 - A52 * V2) / A33
    !
    ! **  j = 4.
    !
    V1 = A41 * A11
    V2 = A42 * A22
    V3 = A43 * A33
    A44 = A44 - A41 * V1 - A42 * V2 - A43 * V3
    A54 = (A54 - A51 * V1 - A52 * V2 - A53 * V3) / A44
    !
    ! **  j = 5.
    !
    V1 = A51 * A11
    V2 = A52 * A22
    V3 = A53 * A33
    V4 = A54 * A44
    A55 = A55 - A51 * V1 - A52 * V2 - A53 * V3 - A54 * V4
    !
    ! **  Ly=b.
    !
    X2 = X2 - A21 * X1
    X3 = X3 - A31 * X1 - A32 * X2
    X4 = X4 - A41 * X1 - A42 * X2 - A43 * X3
    X5 = X5 - A51 * X1 - A52 * X2 - A53 * X3 - A54 * X4
    !
    ! **  Dz=y.
    !
    X1 = X1 / A11
    X2 = X2 / A22
    X3 = X3 / A33
    X4 = X4 / A44
    X5 = X5 / A55
    !
    ! **  L'x=z.
    !
    X4 = X4 - A54 * X5
    X3 = X3 - A43 * X4 - A53 * X5
    X2 = X2 - A32 * X3 - A42 * X4 - A52 * X5
    X1 = X1 - A21 * X2 - A31 * X3 - A41 * X4 - A51 * X5
    x(0, :, :) = X1
    x(1, :, :) = X2
    x(2, :, :) = X3
    x(3, :, :) = X4
    x(4, :, :) = X5
    !
    END SUBROUTINE SOLVIT
    !
    !===========================================================================
    !
    SUBROUTINE CHECK_DIAGONALS(STIF, NEWTON_OK, N, M)
    !
    ! Determine where diagonal elements are small.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! STIF:
    ! NEWTON_OK:
    ! N:
    ! M:
    !
    REAL(RK),  INTENT(IN) :: STIF(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    LOGICAL, INTENT(INOUT) :: NEWTON_OK(0:(N - 1), 0:(M - 1))
    INTEGER :: N
    INTEGER :: M
    !
    ! Locals:
    !
    INTEGER :: I
    !
    !---------------------------------------------------------------------------
    !
    DO I = 0, TVEC1
        !
        WHERE (ABS(STIF(I, I, :, :)) .LT. VTINY)
            !
            NEWTON_OK = .FALSE.
            !
        END WHERE
        !
    END DO
    !
    END SUBROUTINE CHECK_DIAGONALS
    !
    !===========================================================================
    !
    SUBROUTINE SCALE_UP_SIGM(SIG, EPSEFF, N, M)
    !
    ! Rescale the stress after solution is found
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! SIG: Stress (input/output)
    ! EPSEFF: Effective deformation rate
    ! N: Number of grains
    ! M: Number of elements
    !
    REAL(RK), INTENT(INOUT) :: SIG(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: EPSEFF(0:(N - 1), 0:(M - 1))
    INTEGER :: N
    INTEGER :: M
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: NUMIND
    INTEGER :: IPHASE
    INTEGER, POINTER :: INDICES(:)
    REAL(RK) :: SCALE(0:(N - 1), 0:(M - 1))
    INTEGER :: MY_PHASE(0:(M - 1))
    REAL(RK) :: XM_FAKE
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    !
    DO IPHASE = 1, NUMPHASES
        !
        !-tm
        XM_FAKE = 0.4D0
        !XM_FAKE=0.02d0
        !
        SCALE = EPSEFF ** XM_FAKE
        !scale = EPSEFF**CRYSTAL_PARM(0,IPHASE)
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
        !
        DO J = 0, (N - 1)
            !
            DO I = 0, TVEC1
                !
                SIG(I, J, INDICES) = SIG(I, J, INDICES) * SCALE(J, INDICES)
                !
            END DO
            !
        END DO
        !
        DEALLOCATE(INDICES)
        !
    END DO !NUMPHASES
    !
    END SUBROUTINE SCALE_UP_SIGM
    !
    !===========================================================================
    !
    SUBROUTINE COMPUTE_AVG_CRSS(CRSS, CRSS_AVG, M_EL)
    !
    ! Computes the average strength of the crystal slip systems
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! CRSS:
    ! CRSS_AVG:
    ! M_EL:
    !
    REAL(RK), INTENT(IN) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: CRSS_AVG(0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER, INTENT(IN) :: M_EL
    !
    ! Locals:
    !
    INTEGER, POINTER :: INDICES(:)=>NULL()
    INTEGER :: ISLIP
    INTEGER :: N_SLIP
    INTEGER :: IPHASE
    INTEGER :: NUMIND
    INTEGER :: MY_PHASE(0:M_EL-1)
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    !
    ! Goes through the total number of phases in the material
    !
    DO IPHASE = 1, NUMPHASES
        !
        ! Finds the numbers of slip systems
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE))
        !
        N_SLIP = CTYPE(IPHASE)%NUMSLIP
        !
        ! Finds the indices corresponding to the current phase the loop is on
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
        !
        ! Sums up the crystal slip strengths corresponding to the given indices
        !
        DO ISLIP = 0, N_SLIP - 1
            !
            CRSS_AVG(:, INDICES + EL_SUB1) = CRSS_AVG(:, INDICES + EL_SUB1) + &
                & CRSS(ISLIP, :, INDICES + EL_SUB1)
            !
        END DO
        !
        ! Calculates the average of these slip systems
        !
        CRSS_AVG(:, INDICES + EL_SUB1) = CRSS_AVG(:, INDICES + EL_SUB1) / N_SLIP
        !
        DEALLOCATE (INDICES)
        !
    END DO
    !
    END SUBROUTINE COMPUTE_AVG_CRSS
    !
END MODULE STRESS_SOLVE_VP_MOD
