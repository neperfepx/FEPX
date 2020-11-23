! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE ITERATE_STRESS_VP_MOD
!
! Module handling the iteration for a single time increment for the viscoplastic
!   solution.
!
! Contains subroutines:
! RECOVER_PRESSURE_VP: Compute pressure from velocity field.
!
! Contains functions:
! ITMETHOD_VP: Driver for the viscoplastic iteration required for a single time
!   increment.
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
!
! From libfepx:
!
USE CONJUGATE_GRADIENT_MOD
USE CONVERGENCE_MOD, ONLY: CV_OPTIONS
USE DIMENSIONS_MOD
USE MICROSTRUCTURE_MOD
USE READ_INPUT_MOD
USE STIFFNESS_VP_MOD
USE SURFACE_MOD
USE UNITS_MOD
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
PUBLIC :: ITMETHOD_VP
!
CONTAINS
    !
    INTEGER FUNCTION ITMETHOD_VP(ITYPE, BCS, PFORCE, VEL, ELPRESS, EVEL, &
        & DOF_TRACE, NP_TRACE, QR5X5, WTS, EPSEFF, DTIME, INCR)
    !
    ! Driver for the viscoplastic iteration required for a single time
    !   increment
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! ITYPE:
    ! BCS:
    ! PFORCE:
    ! VEL:
    ! ELPRESS:
    ! EVEL:
    ! DOF_TRACE:
    ! NP_TRACE:
    ! QR5X5:
    ! WTS:
    ! EPSEFF:
    ! DTIME:
    ! INCR:
    !
    INTEGER :: ITYPE
    LOGICAL :: BCS(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: PFORCE(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: VEL(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: ELPRESS(EL_SUB1:EL_SUP1)
    REAL(RK) :: EVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    TYPE(TRACE) :: DOF_TRACE
    TYPE(TRACE) :: NP_TRACE
    REAL(RK) :: QR5X5(0:TVEC1, 0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: EPSEFF(EL_SUB1:EL_SUP1)
    REAL(RK) :: DTIME
    INTEGER  :: INCR
    !
    ! Locals:
    !
    INTEGER :: N_SLIP
    INTEGER :: N_EDGE
    REAL(RK) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER :: ITER
    INTEGER :: IDIV
    INTEGER :: I
    INTEGER :: J
    INTEGER :: CG_ITER_OUT
    INTEGER :: CG_MAX_ITERS
    REAL(RK) :: CG_TOL
    REAL(RK) :: EPS_O
    REAL(RK) :: D_NORM
    REAL(RK) :: PART_DELOMAX
    REAL(RK) :: DELOMAX
    REAL(RK) :: PART_U_NORM
    REAL(RK) :: U_NORM
    REAL(RK) :: PART_X_NORM
    REAL(RK) :: X_NORM
    REAL(RK) :: PART_P_NORM
    REAL(RK) :: P_NORM
    REAL(RK) :: ESTIFF(0:KDIM1, 0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ECOORDS(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: FORCE(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: EFORCE(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: VEL_O(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: DEL(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: GDIAG(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: PSCALE(EL_SUB1:EL_SUP1)
    REAL(RK) :: PCNST(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: EQPLAS_TR(EL_SUB1:EL_SUP1)
    !
    !---------------------------------------------------------------------------
    !
    ITMETHOD_VP = 1
    !
    CG_MAX_ITERS = CV_OPTIONS%CG_MAX_ITERS
    CG_TOL = CV_OPTIONS%CG_TOL
    !
    EPS_O = 1.0D30
    VEL_O = VEL
    !
    CALL PART_GATHER(EVEL, VEL, NODES, DOF_TRACE)
    CALL PART_GATHER(ECOORDS, COORDS, NODES, DOF_TRACE)
    !
    ! Non linear iteration loop
    !
    NONLINEAR_ITERATION : DO ITER = 1, CV_OPTIONS%NL_MAX_ITERS
        !
        IF (MYID .EQ. 0) THEN
            !
            WRITE(DFLT_U,'(A,I0)') 'Info   :     > ITMETHOD_VP: Iteration ', &
                & ITER
            WRITE(DFLT_U,'(A)',ADVANCE = 'NO') 'Info   :       . Solving NL &
                &iteration... '
            !
        END IF
        !
        ESTIFF = 0.0D0
        EFORCE = 0.0D0
        FORCE = PFORCE !PFORCE=0
        EQPLAS_TR = 0.0D0
        !
        CALL ELEMENT_STIF_VP(ITYPE, ESTIFF, ECOORDS, EVEL, PSCALE, PCNST, &
            & QR5X5, WTS, EQPLAS_TR, EPSEFF, DTIME, INCR)
        !
        DO I = 0, KDIM1
            !
            EFORCE(I, :) = EFORCE(I, :) - ELPRESS * PCNST(I, :)
            !
        END DO
        !
        ! EFORCE --> FORCE
        !
        CALL PART_SCATTER(FORCE, EFORCE, NODES, .FALSE., DOF_TRACE)
        !
        ! Zero the forces WHERE velocities are specified
        !
        WHERE (BCS)
            !
            FORCE = 0.0D0
            !
        END WHERE
        !
        ! Preconditioned conjugate gradient solver
        !
        ! Form the diagonal part of the stiffness matrix
        !
        CALL ASSEMBLE_DIAGONALS(GDIAG, ESTIFF, KDIM,DOF_SUB1, DOF_SUP1, &
            & EL_SUB1, EL_SUP1, DOF_TRACE, NODES)
        !
        ! Compute the velocity field (VEL) using the conjugate gradient method
        !
        CG_ITER_OUT = CG_SOLVER_EBE(VEL, D_NORM, FORCE, GDIAG, ESTIFF, BCS, &
            & KDIM, DOF_SUB1, DOF_SUP1, EL_SUB1, EL_SUP1, CG_MAX_ITERS, &
            & CG_TOL, DOF_TRACE, NODES)
        !
        ! VEL --> EVEL
        !
        CALL PART_GATHER(EVEL, VEL, NODES, DOF_TRACE)
        !
        CALL RECOVER_PRESSURE_VP(EVEL, PSCALE, ELPRESS, PCNST)
        !
        PART_U_NORM = SUM((VEL_O - VEL) * (VEL_O - VEL))
        !
        CALL PAR_SUM(PART_U_NORM, U_NORM)
        !
        U_NORM = DSQRT(U_NORM)
        !
        PART_X_NORM = MAXVAL(ABS(VEL_O - VEL))
        !
        CALL PAR_MAX(PART_X_NORM, X_NORM)
        !
        IF (MYID .EQ. 0) THEN
            !
            WRITE(DFLT_U,'(A,E10.4,A,I0,A)', ADVANCE='YES') 'R = ', U_NORM, &
                & ' (', CG_ITER_OUT, ' iters)'
            !
        END IF
        !
        VEL_O = VEL
        !
        IF (ITER .GT. 1) GO TO 20
        !
        IF (U_NORM .GT. EPS_O) THEN
            !
            IDIV = IDIV + 1
            !
        ELSE
            !
            IDIV = 0
            !
        END IF
        !
        IF (IDIV .GE. 5) THEN
            !
            ITMETHOD_VP = -1
            !
            RETURN
            !
        END IF
        !
        EPS_O = U_NORM
        !
    END DO NONLINEAR_ITERATION
    !
    ITMETHOD_VP = -1
    !
    RETURN
    !
    20 CONTINUE
    !
    ! DEB I don't know why one more iteration is required here, but it seems
    !   like it can't hurt.  Note that the stiffness is not recomputed.
    !
    FORCE = PFORCE
    !
    DO I = 0, KDIM1
        !
        EFORCE(I, :) = -ELPRESS * PCNST(I, :)
        !
    END DO
    !
    ! EFORCE --> FORCE
    !
    CALL PART_SCATTER(FORCE, EFORCE, NODES, .FALSE., DOF_TRACE)
    !
    ! Zero the forces where velocities are specified
    !
    WHERE (BCS) FORCE = 0.0D0
    !
    ! Compute the velocity field (VEL) using the conjugate gradient method
    !
    CG_ITER_OUT = CG_SOLVER_EBE(VEL, D_NORM, FORCE, GDIAG, ESTIFF, BCS, KDIM, &
        & DOF_SUB1, DOF_SUP1, EL_SUB1, EL_SUP1, CG_MAX_ITERS, CG_TOL, &
        & DOF_TRACE, NODES)
    !
    ! VEL --> EVEL
    !
    CALL PART_GATHER(EVEL, VEL, NODES, DOF_TRACE)
    !
    CALL RECOVER_PRESSURE_VP(EVEL, PSCALE, ELPRESS, PCNST)
    !
    PART_U_NORM = SUM((VEL_O - VEL) * (VEL_O - VEL))
    !
    CALL PAR_SUM(PART_U_NORM, U_NORM)
    !
    U_NORM = SQRT(U_NORM)
    DEL = ABS(VEL_O)
    PART_DELOMAX = MAXVAL(DEL)
    !
    CALL PAR_MAX(PART_DELOMAX, DELOMAX)
    !
    IF (DELOMAX .EQ. 0.0) DELOMAX = 1.0D0
    !
    DEL = ABS((VEL_O - VEL) / DELOMAX)
    PART_P_NORM = MAXVAL(DEL)
    !
    CALL PAR_MAX(PART_P_NORM, P_NORM)
    !
    PART_X_NORM = MAXVAL(ABS(VEL_O - VEL))
    !
    CALL PAR_MAX(PART_X_NORM, X_NORM)
    !
    ITER = ITER + 1
    !
    VEL_O = VEL
    !
    IF (MYID .EQ. 0)WRITE(DFLT_U,'(A,I0,A)') 'Info   :     > Converged in ',&
        & ITER-1, ' iterations'
    !
    IF (MYID .EQ. 0)WRITE(DFLT_U,'(A,A)') 'Info   :     > Ready for ',&
        &'anisotropic elasto-viscoplastic simulation.'
    !
    RETURN
    !
    END FUNCTION ITMETHOD_VP
    !
    !===========================================================================
    !
    SUBROUTINE RECOVER_PRESSURE_VP(EVEL, PSCALE, ELPRESS, PCNST)
    !
    ! Compute pressure from VELocity field.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! EVEL:
    ! PSCALE:
    ! ELPRESS:
    ! PCNST:
    !
    REAL(RK) :: EVEL (0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: PSCALE (EL_SUB1:EL_SUP1)
    REAL(RK) :: ELPRESS(EL_SUB1:EL_SUP1)
    REAL(RK) :: PCNST(0:KDIM1, EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    INTEGER :: I
    REAL(RK) :: SUM(EL_SUB1:EL_SUP1)
    !
    !---------------------------------------------------------------------------
    !
    SUM = 0.0D0
    !
    DO I = 0, KDIM1
        !
        SUM = SUM + PCNST(I, :) * EVEL(I, :)
        !
    END DO
    !
    ELPRESS = ELPRESS + CV_OPTIONS%PACC * PSCALE * SUM
    !
    RETURN
    !
    END SUBROUTINE RECOVER_PRESSURE_VP
    !
END MODULE ITERATE_STRESS_VP_MOD
