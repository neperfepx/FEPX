! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE ITERATE_STRESS_EVPS_MOD
!
! Module handling the iteration for a single time increment for the EVPS soln.
!
! Contains subroutines:
! RECOVER_PRESSURE_EVPS: Compute pressure from velocity field
!
! Contains functions:
! ITMETHOD_EVPS: Driver for the EVPS iteration required for a single time
!   increment.
!
! From libfepx:
!
USE CONJUGATE_GRADIENT_MOD
USE CONVERGENCE_MOD, ONLY: CV_OPTIONS
USE DIMENSIONS_MOD
USE MATRIX_OPERATIONS_MOD
USE MICROSTRUCTURE_MOD
USE READ_INPUT_MOD
USE STIFFNESS_EVPS_MOD
USE UNITS_MOD
USE WRITE_OUTPUT_MOD
!
! From libparallel:
!
USE GATHER_SCATTER_MOD
USE PARALLEL_MOD
!
IMPLICIT NONE
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: ITMETHOD_EVPS
!
CONTAINS
    !
    INTEGER FUNCTION ITMETHOD_EVPS(BCS, PFORCE, VEL, ELPRESS, EVEL, DTRACE, &
        & C0_ANGS, C_ANGS, SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, &
        & KEINV, E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, JITER_STATE, WTS, EPSEFF, &
        & DTIME, INCR, E_BAR_VEC, CONVERGED_SOLUTION, AUTO_TIME, ITER)
    !
    ! Driver for the EVPS iteration required for a single time increment.
    !
    ! ITMETHOD_EVPS is set:
    !   >0 (1) IF the solution converged
    !   <0 (-1) IF the solution did not converged (either the maximum number of
    !       iterations was exceeded or the solution was badly diverging)
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! BCS:
    ! PFORCE:
    ! VEL:
    ! ELPRESS:
    ! EVEL:
    ! DTRACE: Gather/scatter trace for degrees of freedom
    ! C0_ANGS:
    ! C_ANGS:
    ! SIG_VEC_N: (enters=0 at the first increment for all the grains, exits =/0
    !   only for grain=0)
    ! SIG_VEC:
    ! CRSS_N:
    ! CRSS:
    ! RSTAR_N:
    ! RSTAR:
    ! KEINV:
    ! E_ELAS_KK_BAR:
    ! E_ELAS_KK:
    ! SIG_KK:
    ! JITER_STATE:
    ! WTS:
    ! EPSEFF:
    ! DTIME:
    ! INCR: Current load increment
    ! E_BAR_VEC:
    ! CONVERGED_SOLUTION:
    ! AUTO_TIME:
    ! ITER:
    !
    LOGICAL, INTENT(IN) :: BCS(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(IN) :: PFORCE(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(INOUT) :: VEL(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(OUT) :: ELPRESS(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: EVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    TYPE(TRACE) :: DTRACE
    REAL(RK), INTENT(IN) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: SIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    REAL(RK), INTENT(OUT) :: SIG_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    REAL(RK), INTENT(IN) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    REAL(RK), INTENT(IN) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: KEINV(0:TVEC1, 1:NUMPHASES)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK_BAR(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK), INTENT(OUT) :: SIG_KK(EL_SUB1:EL_SUP1, 0:NQPT1)
    INTEGER, INTENT(OUT) :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: EPSEFF(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: DTIME
    INTEGER, INTENT(IN) :: INCR
    REAL(RK), INTENT(OUT) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, &
        & 0:NQPT1)
    LOGICAL, INTENT(INOUT) :: CONVERGED_SOLUTION
    INTEGER :: AUTO_TIME
    INTEGER, INTENT(OUT) :: ITER
    !
    ! Locals:
    !
    INTEGER, PARAMETER :: NR_SLOPE_START = 2
    INTEGER :: ITMETHOD
    INTEGER :: IDIV
    INTEGER :: I
    INTEGER :: CG_ITER_OUT
    INTEGER :: IER
    INTEGER :: CG_MAX_ITERS
    REAL(RK) :: CG_TOL
    REAL(RK) :: NL_TOL_STRICT
    REAL(RK) :: NL_TOL_LOOSE
    REAL(RK) :: NL_TOL_MIN
    REAL(RK) :: PART_R_NORM
    REAL(RK) :: R_NORM
    REAL(RK) :: R_NORM_O
    REAL(RK) :: R_NORM_N
    REAL(RK) :: PART_RX_NORM
    REAL(RK) :: RX_NORM
    REAL(RK) :: PART_F_NORM
    REAL(RK) :: F_NORM
    REAL(RK) :: PART_DELU_NORM
    REAL(RK) :: DELU_NORM
    REAL(RK) :: DELU_NORM_O
    REAL(RK) :: PART_DELUX_NORM
    REAL(RK) :: DELUX_NORM
    REAL(RK) :: PART_U_NORM
    REAL(RK) :: U_NORM
    REAL(RK) :: D_NORM
    REAL(RK) :: VEL_O(DOF_SUB1:DOF_SUP1) ! Velocity for previous iteration
    REAL(RK) :: DELTA_VEL(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: VEL_SA(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: ESTIFF(0:KDIM1, 0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ETANSTIFF(0:KDIM1, 0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: F_VEC(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TCOORDS(DOF_SUB1:DOF_SUP1), ECOORDS(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: FORCE(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: EFORCE(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: RESID(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: ERESID(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TC_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: TRSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: GDIAG(DOF_SUB1:DOF_SUP1)
    INTEGER  :: M_EL
    !
    ! Local node number for mid-nodes correction
    ! E1, E2: end node
    ! M : mid-node
    !
    INTEGER :: E1
    INTEGER :: E2
    INTEGER :: M
    REAL(RK) :: E_ONES(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: G_ONES(DOF_SUB1:DOF_SUP1)
    LOGICAL :: NR, NR_ATTEMPT
    LOGICAL :: NR_SLOW
    LOGICAL :: ITER_CONVERGED
    INTEGER :: NR_ITER
    REAL(RK) :: NR_TOL_SWITCH_REF
    REAL(RK) :: NR_TOL_CONV
    REAL(RK) :: NR_TOL_SWITCH
    REAL(RK) :: NR_CONV
    REAL(RK) :: NR_CONV_O
    !
    !---------------------------------------------------------------------------
    !
    NR_TOL_SWITCH_REF = CV_OPTIONS%NR_TOL_SWITCH_REF
    NR_TOL_CONV = CV_OPTIONS%NR_TOL_CONV
    !
    CG_MAX_ITERS = CV_OPTIONS%CG_MAX_ITERS
    CG_TOL = CV_OPTIONS%CG_TOL
    !
    NL_TOL_STRICT = CV_OPTIONS%NL_TOL_STRICT
    NL_TOL_LOOSE = CV_OPTIONS%NL_TOL_LOOSE
    NL_TOL_MIN = CV_OPTIONS%NL_TOL_MIN
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    ! Initialize ones arrays
    !
    E_ONES = 1.0D0
    G_ONES = 0.0D0
    !
    ! Scatter ones arrays and store multiplicity
    !
    CALL PART_SCATTER(G_ONES, E_ONES, NODES, DTRACE)
    !
    ! Initialization
    R_NORM_O = 0.0D0
    DELU_NORM_O = 0.0D0
    VEL_O = VEL
    VEL_SA = VEL
    !
    ! Set ITMETHOD_EVPS=1 at the beginning
    !
    ITMETHOD_EVPS = 1
    CG_ITER_OUT = 0
    D_NORM = 0.0D0
    !
    ITER_CONVERGED = .FALSE.
    NR = .FALSE.
    NR_ATTEMPT = .FALSE.
    NR_SLOW = .FALSE.
    NR_ITER = 0
    NR_TOL_SWITCH = NR_TOL_SWITCH_REF
    NR_CONV = 1.0D0
    NR_CONV_O = 1.0D0
    IDIV = 0
    !
    ! Non linear iteration
    !
    NONLINEAR_ITERATION : DO ITER = 1, CV_OPTIONS%NL_MAX_ITERS
        !
        IF (MYID .EQ. 0) THEN
            !
            WRITE(DFLT_U,'(A,I0)') 'Info   :     > ITMETHOD_EVPS: Iteration ', &
                & ITER
            !
        END IF
        !
        ! TSH: Adjust mid-nodes: re-position the mid-node to the middle of the
        !   edge-nodes
        !
        CALL PART_GATHER(ECOORDS, COORDS, NODES, DTRACE)
        !
        M = 1
        E1 = 0
        E2 = 2
        ECOORDS(3 * M, EL_SUB1:EL_SUP1) = (ECOORDS(3 * E1, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2, EL_SUB1:EL_SUP1)) / 2.0D0
        ECOORDS(3 * M + 1, EL_SUB1:EL_SUP1) = &
            & (ECOORDS(3 * E1 + 1, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2 + 1, EL_SUB1:EL_SUP1)) / 2.0D0
        ECOORDS(3 * M + 2, EL_SUB1:EL_SUP1) = &
            & (ECOORDS(3 * E1 + 2, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2 + 2, EL_SUB1:EL_SUP1)) / 2.0D0
        !
        M = 3
        E1 = 2
        E2 = 4
        ECOORDS(3 * M, EL_SUB1:EL_SUP1) = (ECOORDS(3 * E1, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2, EL_SUB1:EL_SUP1)) / 2.0D0
        ECOORDS(3 * M + 1, EL_SUB1:EL_SUP1) = &
            & (ECOORDS(3 * E1 + 1, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2 + 1, EL_SUB1:EL_SUP1)) / 2.0D0
        ECOORDS(3 * M + 2, EL_SUB1:EL_SUP1) = &
            & (ECOORDS(3 * E1 + 2, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2 + 2, EL_SUB1:EL_SUP1)) / 2.0D0
        !
        M = 5
        E1 = 0
        E2 = 4
        ECOORDS(3 * M, EL_SUB1:EL_SUP1) = (ECOORDS(3 * E1, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2, EL_SUB1:EL_SUP1)) / 2.0D0
        ECOORDS(3 * M + 1, EL_SUB1:EL_SUP1) = &
            & (ECOORDS(3 * E1 + 1, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2 + 1, EL_SUB1:EL_SUP1)) / 2.0D0
        ECOORDS(3 * M + 2, EL_SUB1:EL_SUP1) = &
            & (ECOORDS(3 * E1 + 2, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2 + 2, EL_SUB1:EL_SUP1)) / 2.0D0
        !
        M = 6
        E1 = 0
        E2 = 9
        ECOORDS(3 * M, EL_SUB1:EL_SUP1) = (ECOORDS(3 * E1, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2, EL_SUB1:EL_SUP1)) / 2.0D0
        ECOORDS(3 * M + 1, EL_SUB1:EL_SUP1) = &
            & (ECOORDS(3 * E1 + 1, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2 + 1, EL_SUB1:EL_SUP1)) / 2.0D0
        ECOORDS(3 * M + 2, EL_SUB1:EL_SUP1) = &
            & (ECOORDS(3 * E1 + 2, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2 + 2, EL_SUB1:EL_SUP1)) / 2.0D0
        !
        M = 7
        E1 = 2
        E2 = 9
        ECOORDS(3*M, EL_SUB1:EL_SUP1) = (ECOORDS(3 * E1, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2, EL_SUB1:EL_SUP1)) / 2.0D0
        ECOORDS(3 * M + 1, EL_SUB1:EL_SUP1) = &
            & (ECOORDS(3 * E1 + 1, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2 + 1, EL_SUB1:EL_SUP1)) / 2.0D0
        ECOORDS(3 * M + 2, EL_SUB1:EL_SUP1) = &
            & (ECOORDS(3 * E1 + 2, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2 + 2, EL_SUB1:EL_SUP1)) / 2.0D0
        !
        M = 8
        E1 = 4
        E2 = 9
        ECOORDS(3 * M, EL_SUB1:EL_SUP1) = (ECOORDS(3 * E1, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2,EL_SUB1:EL_SUP1)) / 2.0D0
        ECOORDS(3 * M + 1, EL_SUB1:EL_SUP1) = &
            & (ECOORDS(3 * E1 + 1, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2 + 1, EL_SUB1:EL_SUP1)) / 2.0D0
        ECOORDS(3 * M + 2, EL_SUB1:EL_SUP1) = &
            & (ECOORDS(3 * E1 + 2, EL_SUB1:EL_SUP1) + &
            & ECOORDS(3 * E2 + 2, EL_SUB1:EL_SUP1)) / 2.0D0
        !
        ! Reset coordinates
        !
        COORDS = 0.0D0
        !
        ! Scatter e_coordinates
        !
        CALL PART_SCATTER(COORDS, ECOORDS, NODES, DTRACE)
        !
        ! Divide the coordinates by the multiplicity
        !
        COORDS = COORDS / G_ONES
        !
        ! Estimate the new coordinates based on the velocity VEL
        ! (VEL is the guess for the velocity field)
        !
        TCOORDS = COORDS + DTIME * VEL
        !
        CALL PART_GATHER(ECOORDS, TCOORDS, NODES, DTRACE)
        CALL PART_GATHER(EVEL, VEL, NODES, DTRACE)
        !
        ESTIFF = 0.0D0
        ETANSTIFF = 0.0D0
        EFORCE = 0.0D0
        FORCE  = PFORCE !PFORCE=0
        TC_ANGS = SPREAD(C_ANGS, 5, NQPT)
        TRSTAR = SPREAD(RSTAR, 5, NQPT)
        !
        ! Dnter the subroutine that, through iterations at the constitutive
        !   level , computes:
        ! The state of each element:
        ! e*_kk: E_ELAS_KK
        ! tau_kk: SIG_KK
        ! tau': SIG_VEC
        ! R*: RSTAR
        ! g: CRSS
        ! Matrices [S]: ESTIFF
        ! Vectors [S]{h}: F_VEC
        !
        ! E_BAR_VEC: OUT
        ! SIG_VEC_N: INOUT (enters=0 at the first increment, only grain=0 is
        !   updated)
        !
        CALL ELEMENT_STIF_EVPS( ESTIFF, ETANSTIFF, F_VEC, ECOORDS, EVEL, &
            & C0_ANGS, TC_ANGS, SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, &
            & TRSTAR, E_BAR_VEC, E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, &
            & JITER_STATE, KEINV, WTS, EPSEFF, DTIME, INCR, &
            & CONVERGED_SOLUTION, AUTO_TIME, NR)
        !
        ! Essentially renames F_VEC EFORCE, consider removing
        !
        DO I = 0, KDIM1 ! 29
            !
            EFORCE(I, :) = EFORCE(I, :) + F_VEC(I, :)
            !
        END DO
        !
        ! Compute elemental norms
        !
        ERESID = 0.0D0
        RESID = 0.0D0
        PART_R_NORM = 0.0D0
        R_NORM = 0.0D0
        PART_RX_NORM = 0.0D0
        RX_NORM = 0.0D0
        PART_F_NORM = 0.0D0
        F_NORM = 0.0D0
        !
        CALL GEN_MATRIX_VECTOR_MULT(ERESID, ESTIFF, EVEL, IER)
        !
        DO I = 0, KDIM1 ! 29
            !
            ERESID(I,:) = EFORCE(I,:) - ERESID(I,:)
            !
        END DO
        !
        ! Scatter elemental residuals to all nodes
        !
        CALL PART_SCATTER(RESID, ERESID, NODES, DTRACE)
        !
        ! Calculate residual and force magnitudes
        !
        PART_R_NORM = SUM(RESID * RESID)
        !
        CALL PAR_SUM(PART_R_NORM, R_NORM)
        !
        R_NORM = DSQRT(R_NORM)
        !
        IF (ITER .EQ. 1) THEN
            !
            R_NORM_N = R_NORM
            !
        END IF
        !
        PART_RX_NORM = MAXVAL(ABS(RESID))
        !
        CALL PAR_MAX(PART_RX_NORM, RX_NORM)
        !
        CALL PART_SCATTER(FORCE, EFORCE, NODES, DTRACE)
        !
        PART_F_NORM = SUM(FORCE * FORCE)
        !
        CALL PAR_SUM(PART_F_NORM, F_NORM)
        !
        F_NORM = DSQRT(F_NORM)
        !
        ! Solve for new velocity field
        ! Zero the residual where velocities are specified
        !
        WHERE (BCS)
            !
            RESID = 0.0D0
            !
        END WHERE
        !
        DELTA_VEL = 0.0D0
        !
        IF (NR) THEN ! Use Newton-Raphson
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A)', ADVANCE='NO') 'Info   :       . Solving NR &
                    &iteration... '
                !
            END IF
            !
            ITMETHOD = 1
            NR_ITER = NR_ITER+1
            !
            ! Form the diagonal part of the stiffness matrix
            !
            CALL ASSEMBLE_DIAGONALS(GDIAG, ETANSTIFF, KDIM, DOF_SUB1, &
                & DOF_SUP1, EL_SUB1, EL_SUP1, DTRACE, NODES)
            !
            ! Compute the velocity field (VEL) using the conjugate gradient
            !   method
            !
            CG_ITER_OUT = CG_SOLVER_EBE(DELTA_VEL, D_NORM, RESID, GDIAG, &
                & ETANSTIFF, BCS, KDIM, DOF_SUB1, DOF_SUP1, EL_SUB1, EL_SUP1, &
                & CG_MAX_ITERS, CG_TOL, DTRACE, NODES)
            !
        ELSE
            !
            ! Use Successive Approximations
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A)', ADVANCE='NO') 'Info   :       . Solving SA &
                    &iteration... '
                !
            END IF
            !
            ITMETHOD = 0
            !
            CALL ASSEMBLE_DIAGONALS(GDIAG, ESTIFF, KDIM, DOF_SUB1, DOF_SUP1, &
                & EL_SUB1, EL_SUP1, DTRACE, NODES)
            !
            ! Compute the velocity field (VEL) using the conjugate gradient
            !   method
            CG_ITER_OUT = CG_SOLVER_EBE(DELTA_VEL, D_NORM, RESID, GDIAG, &
                & ESTIFF, BCS, KDIM, DOF_SUB1, DOF_SUP1, EL_SUB1, EL_SUP1, &
                & CG_MAX_ITERS, CG_TOL, DTRACE, NODES)
            !
        END IF
        !
        VEL = VEL_O + DELTA_VEL
        !
        ! Calculate velocity norm
        !
        PART_DELU_NORM = 0.0D0
        DELU_NORM = 0.0D0
        PART_DELUX_NORM = 0.0D0
        DELUX_NORM = 0.0D0
        PART_U_NORM = 0.0D0
        U_NORM = 0.0D0
        !
        PART_DELU_NORM = SUM(DELTA_VEL * DELTA_VEL)
        CALL PAR_SUM(PART_DELU_NORM, DELU_NORM)
        DELU_NORM = DSQRT(DELU_NORM)
        !
        PART_DELUX_NORM = MAXVAL(ABS(DELTA_VEL))
        CALL PAR_MAX(PART_DELUX_NORM, DELUX_NORM)
        !
        PART_U_NORM = SUM(VEL_O * VEL_O)
        CALL PAR_SUM(PART_U_NORM, U_NORM)
        U_NORM = DSQRT(U_NORM)
        !
        ! Normalize velocity norms
        !
        DELU_NORM = DELU_NORM/U_NORM
        DELUX_NORM = DELUX_NORM/U_NORM
        !
        ! Newton-Raphson convergence parameters
        IF (NR_ITER .GE. NR_SLOPE_START) THEN
            !
            NR_CONV = LOG(DELU_NORM_O) - LOG(DELU_NORM)
            !
        END IF
        !
        IF (NR_ITER .EQ. NR_SLOPE_START) THEN
            !
            NR_CONV_O = NR_CONV
            !
        END IF
        !
        ! Successive Approximation convergence parameter
        !
        IF (.NOT. NR) THEN
            !
            IF (DELU_NORM .GT. DELU_NORM_O) THEN
                !
                IDIV = IDIV + 1
                !
            ELSE
                !
                IDIV = 0
                !
            END IF
            !
        END IF
        !
        ! Output convergence parameters
        !
        IF ((MYID .EQ. 0) .AND. PRINT_OPTIONS%PRINT_CONV) THEN
            !
            CALL WRITE_CONV_FILE_DATA(INCR, ITER, ITMETHOD, R_NORM, RX_NORM, &
                & F_NORM, DELU_NORM, DELUX_NORM, U_NORM, CG_ITER_OUT)
            !
        END IF
        !
        IF (MYID .EQ. 0) THEN
            !
            WRITE(DFLT_U,'(A,E10.4,A,I0,A)', ADVANCE='YES') 'R = ', DELU_NORM, &
                & ' (', CG_ITER_OUT, ' iters)'
            !
        END IF
        !
        ! Check convergence of the velocity solution
        ! Cases 1-4:  convergence, strict or loose
        ! Cases 5-6:  switch back to SA due to problems with NR
        ! Case 7:  failure
        !
        IF ((DELU_NORM < NL_TOL_STRICT) .AND. (ITER .GT. 1)) THEN
            !
            ! Case 1
            ! Solution converged
            !
            ITER_CONVERGED = .TRUE.
            !
            EXIT
            !
        ELSE IF ((DELU_NORM * U_NORM < NL_TOL_MIN * MAXDOF) .AND. &
            & (ITER .GT. 1)) THEN
            !
            ! Case 2
            ! Solution converged
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A)') 'Info   :       . Change in velocity is &
                    &below threshold value.'
                !
            END IF
            !
            ITER_CONVERGED = .TRUE.
            !
            EXIT
            !
        ELSE IF ((NR_ITER .GT. NR_SLOPE_START) .AND. (NR_CONV .LT. &
            & NR_TOL_CONV*NR_CONV_O) .AND. (DELU_NORM < NL_TOL_LOOSE)) THEN
            !
            ! Case 3
            ! Newton-Raphson slow to converge, but converged satisfactorily
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A)') 'Info   :       . Newton-Raphson is slow &
                    &to converge, but converged satisfactorily.'
                !
            END IF
            !
            ITER_CONVERGED = .TRUE.
            !
            EXIT
            !
        ELSE IF ((.NOT. NR) .AND. (NR_ATTEMPT) .AND. &
            & (DELU_NORM < NL_TOL_LOOSE)) THEN
            !
            ! Case 4
            ! SA slow to converge, but converged satisfactorily
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A)') 'Info   :       . Successive Approximations&
                    & is slow to converge, but converged satisfactorily.'
                !
            END IF
            !
            ITER_CONVERGED = .TRUE.
            !
            EXIT
            !
        ELSE IF (NR .AND. (R_NORM .GT. 1.1 * R_NORM_N)) THEN
            !
            ! Case 5
            ! Newton-Raphson solution diverging
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U, '(A)') 'Warning:     > Newton-Raphson solution &
                    &diverging.'
                !
            END IF
            !
            ! Revert to previous SA solution and use SA
            !
            NR = .FALSE.
            VEL = VEL_SA
            !
            ! Enforce stricter switch-over tolerance between methods
            !
            NR_TOL_SWITCH = NR_TOL_SWITCH / 10.0
            !
            ! Set flag that NR had been attempted
            !
            NR_ATTEMPT = .TRUE.
            !
            ! Reset number of NR iterations
            !
            NR_ITER = 0
            !
        ELSE IF ((NR_ITER .GT. NR_SLOPE_START) .AND. (NR_CONV .LT. &
            & NR_TOL_CONV*NR_CONV_O)) THEN
            !
            ! Case 6
            ! Newton-Raphson converging slowly
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U, '(A)') 'Warning:     > Newton-Raphson converging &
                    &slowly.'
                !
            END IF
            !
            ! Switch to SA
            !
            NR = .FALSE.
            !
            ! Set flag that NR had been attempted
            !
            NR_ATTEMPT = .TRUE.
            NR_SLOW = .TRUE.
            !
            ! Reset number of NR iterations
            !
            NR_ITER = 0
            !
        ELSE IF (IDIV .ge. 5) THEN
            !
            ! Case 7
            ! Successive approximations diverging
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U, '(A)') 'Warning:     > Iterations are diverging.'
                !
            END IF
            !
            ! Failure to converge
            !
            ITMETHOD_EVPS = -1
            !
            RETURN
            !
        END IF
        !
        ! Switch from SA to NR
        ! Removed requirement for iter .GT. 1
        IF ((DELU_NORM .LT. NR_TOL_SWITCH) .AND. (.NOT. NR_SLOW) .AND. &
            & (.NOT. NR)) THEN
            !
            NR = .TRUE.
            IDIV = 0
            !
            ! Save velocity field
            !
            VEL_SA = VEL
            !
        END IF
        !
        VEL_O = VEL
        R_NORM_O = R_NORM
        DELU_NORM_O = DELU_NORM
        !
    ENDDO NONLINEAR_ITERATION
    !
    IF (ITER_CONVERGED) THEN
        !
        IF (MYID .EQ. 0) THEN
            !
            WRITE(DFLT_U, '(A,I0,A)') 'Info   :     > Converged in ', iter,&
                & ' iterations'
            !
        END IF
        !
        CALL RECOVER_PRESSURE_EVPS(ELPRESS, E_ELAS_KK)
        !
    ELSE
        !
        ! Maximum number of iterations exceeded. Failure to converge.
        !
        ITMETHOD_EVPS = -1
        !
    END IF
    !
    END FUNCTION ITMETHOD_EVPS
    !
    !===========================================================================
    !
    SUBROUTINE RECOVER_PRESSURE_EVPS(ELPRESS, E_ELAS_KK)
    !
    ! Compute pressure from velocity field
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! ELPRESS:
    ! E_ELAS_KK:
    !
    REAL(RK), INTENT(OUT) :: ELPRESS(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: E_ELAS_KK(EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    INTEGER :: IPHASE
    INTEGER :: MY_PHASE(EL_SUB1:EL_SUP1)
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    !
    DO IPHASE = 1, NUMPHASES
        !
        WHERE (MY_PHASE .EQ. IPHASE)
            !
            ELPRESS = -CRYSTAL_PARM(8, IPHASE) * E_ELAS_KK
            !
        END WHERE
        !
    END DO
    !
    END SUBROUTINE RECOVER_PRESSURE_EVPS
    !
END MODULE ITERATE_STRESS_EVPS_MOD
