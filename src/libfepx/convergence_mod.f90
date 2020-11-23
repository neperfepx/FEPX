! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE CONVERGENCE_MOD
!
! Convergence tolerances: data and input routines
!
! Contains subroutines:
! Various subroutines related to convergence options
!
! Contains functions:
! CONVERGENCEKEYWORDINPUT: Keyword input for module
!
! From libf95:
!
USE LIBF95, RK=>REAL_KIND
!
IMPLICIT NONE
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: CONVOPTIONSTYPE
PUBLIC :: CONVERGENCEKEYWORDINPUT
PUBLIC :: CV_OPTIONS
!
TYPE CONVOPTIONSTYPE
    !
    ! Conjugate Gradient (CG)
    INTEGER :: CG_MAX_ITERS = 16000
    REAL(RK) :: CG_TOL = 1.0D-8
    !
    ! Global Nonlinear Iteration (NL)
    INTEGER :: NL_MAX_ITERS = 50
    REAL(RK) :: NL_TOL_STRICT = 5.0D-4
    REAL(RK) :: NL_TOL_LOOSE = 5.0D-4
    REAL(RK) :: NL_TOL_MIN = 1.0D-10
    !
    ! Newton-Raphson
    REAL(RK) :: NR_TOL_SWITCH_REF = 1.0D-2
    REAL(RK) :: NR_TOL_CONV = 0.2D0
    !
    ! PACC (pressure acceleration: obsolete?)
    REAL(RK) :: PACC = 1.0D0
    !
    ! Single Crystal/ State Calculations (SX)
    INTEGER :: SX_MAX_ITERS_STATE = 100
    INTEGER :: SX_MAX_ITERS_NEWTON = 100
    REAL(RK) :: SX_TOL = 1.0D-4
    !
END TYPE CONVOPTIONSTYPE
!
TYPE(CONVOPTIONSTYPE) :: CV_OPTIONS
!
CONTAINS
    !
    FUNCTION CONVERGENCEKEYWORDINPUT(CMDLINE, INUNIT, STATUS) RESULT(MATCH)
    !
    !  Keyword input for this module.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! CMDLINE: Input command line
    ! INUNIT: Input unit for further data if necessary
    ! STATUS: Command status if command line was matched
    ! MATCH: (Result) true if the command line was matched
    !
    CHARACTER(LEN=*), INTENT(IN) :: CMDLINE
    INTEGER, INTENT(IN) :: INUNIT
    INTEGER, INTENT(OUT) :: STATUS
    LOGICAL :: MATCH
    !
    ! Locals:
    !
    CHARACTER(LEN=*), PARAMETER :: MYNAME = 'CONVERGENCEKEYWORDINPUT'
    !
    !--------------*-------------------------------------------------------
    !
    STATUS = 0
    MATCH = .TRUE.
    !
    IF (EXECCOMMAND('cg_max_iters', CMDLINE, EXEC_CG_MAX_ITERS, STATUS)) RETURN
    IF (EXECCOMMAND('cg_tol', CMDLINE, EXEC_CG_TOL, STATUS)) RETURN
    IF (EXECCOMMAND('nl_max_iters', CMDLINE, EXEC_NL_MAX_ITERS, STATUS)) RETURN
    IF (EXECCOMMAND('nl_tol_strict', CMDLINE, &
        & EXEC_NL_TOL_STRICT, STATUS)) RETURN
    IF (EXECCOMMAND('nl_tol_loose', CMDLINE, EXEC_NL_TOL_LOOSE, STATUS)) RETURN
    IF (EXECCOMMAND('nl_tol_min', CMDLINE, EXEC_NL_TOL_MIN, STATUS)) RETURN
    IF (EXECCOMMAND('nr_tol_switch_ref', CMDLINE, &
        & EXEC_NR_TOL_SWITCH_REF, STATUS)) RETURN
    IF (EXECCOMMAND('nr_tol_conv', CMDLINE, EXEC_NR_TOL_CONV, STATUS)) RETURN
    IF (EXECCOMMAND('sx_max_iters_state', CMDLINE, &
        & EXEC_SX_MAX_ITERS_STATE, STATUS)) RETURN
    IF (EXECCOMMAND('sx_max_iters_newton', CMDLINE, &
        & EXEC_SX_MAX_ITERS_NEWTON, STATUS)) RETURN
    IF (EXECCOMMAND('sx_tol', CMDLINE, EXEC_SX_TOL, STATUS)) RETURN
    !
    ! No calls matched the keyword.
    !
    MATCH = .FALSE.
    !
    END FUNCTION CONVERGENCEKEYWORDINPUT
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_CG_MAX_ITERS(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A
    INTEGER, INTENT(OUT) :: S
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    READ(A, *, IOSTAT = S) CV_OPTIONS%CG_MAX_ITERS
    !
    END SUBROUTINE EXEC_CG_MAX_ITERS
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_CG_TOL(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A
    INTEGER, INTENT(OUT) :: S
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    READ(A, *, IOSTAT = S) CV_OPTIONS%CG_TOL
    !
    END SUBROUTINE EXEC_CG_TOL
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NL_MAX_ITERS(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A
    INTEGER,          INTENT(OUT) :: S
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    READ(A, *, IOSTAT = S) CV_OPTIONS%NL_MAX_ITERS
    !
    END SUBROUTINE EXEC_NL_MAX_ITERS
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NL_TOL_STRICT(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A
    INTEGER,          INTENT(OUT) :: S
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    READ(A, *, IOSTAT = S) CV_OPTIONS%NL_TOL_STRICT
    !
    END SUBROUTINE EXEC_NL_TOL_STRICT
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NL_TOL_LOOSE(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A
    INTEGER,          INTENT(OUT) :: S
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    READ(A, *, IOSTAT = S) CV_OPTIONS%NL_TOL_LOOSE
    !
    END SUBROUTINE EXEC_NL_TOL_LOOSE
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NL_TOL_MIN(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A
    INTEGER,          INTENT(OUT) :: S
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    READ(A, *, IOSTAT = S) CV_OPTIONS%NL_TOL_MIN
    !
    END SUBROUTINE EXEC_NL_TOL_MIN
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NR_TOL_SWITCH_REF(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A
    INTEGER,          INTENT(OUT) :: S
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    READ(A, *, IOSTAT = S) CV_OPTIONS%NR_TOL_SWITCH_REF
    !
    END SUBROUTINE EXEC_NR_TOL_SWITCH_REF
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NR_TOL_CONV(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A
    INTEGER,          INTENT(OUT) :: S
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    READ(A, *, IOSTAT = S) CV_OPTIONS%NR_TOL_CONV
    !
    END SUBROUTINE EXEC_NR_TOL_CONV
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_SX_MAX_ITERS_STATE(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A
    INTEGER,          INTENT(OUT) :: S
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    READ(A, *, IOSTAT = S) CV_OPTIONS%SX_MAX_ITERS_STATE
    !
    END SUBROUTINE EXEC_SX_MAX_ITERS_STATE
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_SX_MAX_ITERS_NEWTON(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A
    INTEGER,          INTENT(OUT) :: S
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    READ(A, *, IOSTAT = S) CV_OPTIONS%SX_MAX_ITERS_NEWTON
    !
    END SUBROUTINE EXEC_SX_MAX_ITERS_NEWTON
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_SX_TOL(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A
    INTEGER,          INTENT(OUT) :: S
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    READ(A, *, IOSTAT = S) CV_OPTIONS % SX_TOL
    !
    END SUBROUTINE EXEC_SX_TOL
    !
END MODULE CONVERGENCE_MOD
