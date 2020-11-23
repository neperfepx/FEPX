! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE CONJUGATE_GRADIENT_MOD
!
! Routines for the preconditioned conjugate gradient solver
!
! Contains subroutines:
! ASSEMBLE_DIAGONALS: Form the diagonal part of the stiffness matrix for use in
!   preconditioning
!
! Contains functions:
! CG_SOLVER_EBE: Driver for element by element conjugate gradient solver
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
!
! From libfepx:
!
USE MATRIX_OPERATIONS_MOD
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
PUBLIC :: ASSEMBLE_DIAGONALS
PUBLIC :: CG_SOLVER_EBE
!
CONTAINS
    !
    SUBROUTINE ASSEMBLE_DIAGONALS(DIAGONALS, GSTIF, NNPE, NSUB, NSUP, ESUB, &
        & ESUP, DTRACE, NP)
    !
    ! Form the diagonal part of the stiffness matrix, for use in
    !   preconditioning.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! DIAGONALS: Diagonal matrix (result of this routine)
    ! GSTIF: Global stiffness matrix
    ! NNPE: Number of nodes per element
    ! NSUB, NSUP: Node range for this process
    ! ESUB, ESUP: Element range for this process
    ! DTRACE: Commmunication information for gather/scatter
    ! NP: Connectivity
    !
    REAL(RK) :: DIAGONALS(NSUB:NSUP)
    REAL(RK) :: GSTIF(0:(NNPE - 1), 0:(NNPE - 1), ESUB:ESUP)
    INTEGER :: NNPE
    INTEGER :: NSUB
    INTEGER :: NSUP
    INTEGER :: ESUB
    INTEGER :: ESUP
    TYPE(TRACE) :: DTRACE
    INTEGER :: NP(0:(NNPE-1), ESUB:ESUP)
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: J
    REAL(RK) :: EDIAGONALS(0:(NNPE - 1), ESUB:ESUP)
    REAL(RK) :: DIAGMIN
    REAL(RK) :: DIAGMAX
    !
    !---------------------------------------------------------------------------
    !
    ! RC 3/24/2016: Reordered for better memory striding
    !
    DO J = ESUB, ESUP
        !
        DO I = 0,(NNPE - 1)
            !
            EDIAGONALS(I, J) = GSTIF(I, I, J)
            !
        END DO
        !
    END DO
    !
    DIAGONALS = 0.0D0
    !
    CALL PART_SCATTER(DIAGONALS, EDIAGONALS, NP, .FALSE., DTRACE)
    !
    DIAGONALS = 1.0D0 / DIAGONALS
    !
    RETURN
    !
    END SUBROUTINE ASSEMBLE_DIAGONALS
    !
    !===========================================================================
    !
    INTEGER FUNCTION CG_SOLVER_EBE(SOL, RES_NORM, RHS, DIAGONALS, GSTIF, BCS, &
        & NNPE, NSUB, NSUP, ESUB, ESUP, MAX_ITER, TOL, DTRACE, NP)
    !
    ! Driver for element by element conjugate gradient solver
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! SOL:
    ! RES_NORM:
    ! RHS:
    ! DIAGONALS:
    ! GSTIF:
    ! BCS:
    ! NNPE:
    ! NSUB:
    ! NSUP:
    ! ESUB:
    ! ESUP:
    ! MAX_ITER:
    ! TOL:
    ! DTRACE:
    ! NP:
    !
    REAL(RK) :: SOL(NSUB:NSUP)
    REAL(RK) :: RES_NORM
    REAL(RK) :: RHS(NSUB:NSUP)
    REAL(RK) :: DIAGONALS(NSUB:NSUP)
    REAL(RK) :: GSTIF(0:(NNPE - 1), 0:(NNPE - 1), ESUB:ESUP)
    LOGICAL :: BCS(NSUB:NSUP)
    INTEGER :: NNPE
    INTEGER :: NSUB
    INTEGER :: NSUP
    INTEGER :: ESUB
    INTEGER :: ESUP
    INTEGER :: MAX_ITER
    REAL(RK) :: TOL
    TYPE(TRACE) :: DTRACE
    INTEGER :: NP(0:(NNPE-1), ESUB:ESUP)
    !
    ! Locals:
    !
    REAL(RK) :: PART_RES_NORM
    INTEGER :: ITER
    INTEGER :: N_ITER
    INTEGER :: IER
    INTEGER :: I
    INTEGER :: J
    REAL(RK) :: ELAPSED
    REAL(RK) :: ZU (NSUB:NSUP)
    REAL(RK) :: RU(NSUB:NSUP)
    REAL(RK) :: APU(NSUB:NSUP)
    REAL(RK) :: PU(NSUB:NSUP)
    REAL(RK) :: TEMP1(0:(NNPE - 1), ESUB:ESUP)
    REAL(RK) :: TEMP2(0:(NNPE - 1), ESUB:ESUP)
    REAL(RK) :: PART_ALPHA
    REAL(RK) :: ALPHA
    REAL(RK) :: BETA
    REAL(RK) :: PART_ERROR
    REAL(RK) :: ERROR
    REAL(RK) :: PART_XNUMER
    REAL(RK) :: XNUMER
    REAL(RK) :: XNUMER_O
    REAL(RK) :: PART_MAG
    REAL(RK) :: MAG
    !
    !---------------------------------------------------------------------------
    !
    ! CG initializations
    !
    CALL SPARSE_MATVEC_EBE(RU, SOL, TEMP1, TEMP2, GSTIF, BCS, NNPE, NSUB, &
        & NSUP, ESUB, ESUP, DTRACE, NP)
    !
    RU  = RHS - RU
    RHS = RHS * DIAGONALS
    ZU  = RU * DIAGONALS
    !
    PART_MAG = SUM(ZU * ZU)
    !
    CALL PAR_SUM(PART_MAG, MAG)
    !
    MAG = DSQRT(MAG)
    !
    ! CG iterations
    !
    XNUMER_O = 1.0D0
    PU = 0.0D0
    ERROR = 1.0D0
    N_ITER = 0
    !
    DO WHILE (ERROR .GT. TOL)
        !
        N_ITER = N_ITER + 1
        !
        IF (N_ITER .GT. MAX_ITER) THEN
            !
            CALL PAR_quit('ERROR  :       . CG_SOLVER_EBE: Convergence failure.')
            !
        END IF
        !
        PART_XNUMER = SUM(RU * DIAGONALS * RU)
        !
        CALL PAR_SUM(PART_XNUMER, XNUMER)
        !
        BETA = XNUMER / XNUMER_O
        !
        IF (N_ITER .EQ. 1) BETA = 0.0D0
        !
        XNUMER_O = XNUMER
        !
        PU = ZU + BETA * PU
        !
        CALL SPARSE_MATVEC_EBE(APU, PU, TEMP1, TEMP2, GSTIF, BCS, NNPE, NSUB, &
            & NSUP, ESUB, ESUP, DTRACE, NP)
        !
        PART_ALPHA = SUM(PU * APU)
        !
        CALL PAR_SUM(PART_ALPHA, ALPHA)
        !
        ALPHA = XNUMER / ALPHA
        !
        SOL = SOL + ALPHA * PU
        RU = RU - ALPHA * APU
        ZU = RU * DIAGONALS
        !
        PART_ERROR = SUM(ZU * ZU)
        !
        CALL PAR_SUM(PART_ERROR, ERROR)
        !
        ERROR = DSQRT(ERROR) / MAG
        !
    END DO
    !
    PART_RES_NORM = SUM(SOL * SOL)
    !
    CALL PAR_SUM(PART_RES_NORM, RES_NORM)
    !
    RES_NORM = (RES_NORM) ** 0.5D0
    !
    CG_SOLVER_EBE = N_ITER
    !
    RETURN
    !
    END FUNCTION CG_SOLVER_EBE
    !
END MODULE CONJUGATE_GRADIENT_MOD
