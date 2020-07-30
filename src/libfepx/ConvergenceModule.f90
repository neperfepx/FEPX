! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE ConvergenceModule
  !
  !  Convergence tolerances: data and input routines
  !
  !  ==================== Other Modules (USE statements)
  !
  USE LibF95, RK=>REAL_KIND

  IMPLICIT NONE
  !
  !  ==================== Public Entities
  !
  !  variables, procedures, constants, derived types and namelist groups
  !
  PRIVATE   ! all objects are private unless declared otherwise

  PUBLIC :: ConvOptionsType
  PUBLIC :: ConvergenceKeywordInput
  PUBLIC :: cv_options
  !
  !  ==================== Derived Types
  !
  !  NOTE: defaults are set in the type definition
  !
  TYPE ConvOptionsType
    ! Conjugate Gradient (CG)
    INTEGER  :: cg_max_iters = 16000
    REAL(RK) :: cg_tol = 1e-8_RK
    ! Global Nonlinear Iteration (NL)
    INTEGER  :: nl_max_iters = 50
    REAL(RK) :: nl_tol_strict = 5e-4_RK
    REAL(RK) :: nl_tol_loose = 5e-4_RK
    REAL(RK) :: nl_tol_min = 1e-10_RK
    ! Newton-Raphson
    REAL(RK) :: nr_tol_switch_ref = 1.0e-2_RK
    REAL(RK) :: nr_tol_conv = 0.2_RK
    ! PACC (pressure acceleration: obsolete?)
    REAL(RK) :: pacc = 1.0_RK
    ! Single Crystal/ State Calculations (SX)
    INTEGER  :: sx_max_iters_state = 100
    INTEGER  :: sx_max_iters_newton = 100
    REAL(RK) :: sx_tol = 1.0e-4

  END TYPE ConvOptionsType
  !
  ! Module Data
  !
  TYPE(ConvOptionsType) :: cv_options
  !
CONTAINS ! ============================================= MODULE PROCEDURES
  !
  ! ==================================================== BEGIN:  ConvergenceKeywordInput
  !
  FUNCTION ConvergenceKeywordInput(cmdLine, inUnit, status) RESULT(match)
    !
    !  Keyword input for this module.
    !
    ! ========== Arguments
    !
    !  .     REQUIRED ARGS
    !  cmdLine -- input command line 
    !  inUnit  -- input unit for further data if necessary
    !  status  -- command status if command line was matched
    !
    !  .     OPTIONAL ARGS
    !
    !  .     FUNCTION RESULT
    !  match -- true if the command line was matched
    !
    CHARACTER(LEN=*), INTENT(IN)  :: cmdLine
    INTEGER,          INTENT(IN)  :: inUnit
    INTEGER,          INTENT(OUT) :: status
    !
    LOGICAL :: match
    !
    !
    ! ========== Locals
    !
    CHARACTER(LEN=*), PARAMETER :: myNAME = 'ConvergenceKeywordInput'
    !
    !--------------*-------------------------------------------------------
    !
    status = 0;  match  = .TRUE.
    !
    IF (ExecCommand('cg_max_iters', cmdLine, exec_cg_max_iters, status)) RETURN 
    IF (ExecCommand('cg_tol', cmdLine, exec_cg_tol, status)) RETURN 
    IF (ExecCommand('nl_max_iters', cmdLine, exec_nl_max_iters, status)) RETURN  
    IF (ExecCommand('nl_tol_strict', cmdLine, exec_nl_tol_strict, status)) RETURN 
    IF (ExecCommand('nl_tol_loose', cmdLine, exec_nl_tol_loose, status)) RETURN 
    IF (ExecCommand('nl_tol_min', cmdLine, exec_nl_tol_min, status)) RETURN 
    IF (ExecCommand('nr_tol_switch_ref', cmdLine, exec_nr_tol_switch_ref, status)) RETURN 
    IF (ExecCommand('nr_tol_conv', cmdLine, exec_nr_tol_conv, status)) RETURN 
    IF (ExecCommand('sx_max_iters_state', cmdLine, exec_sx_max_iters_state, status)) RETURN 
    IF (ExecCommand('sx_max_iters_newton', cmdLine, exec_sx_max_iters_newton, status)) RETURN 
    IF (ExecCommand('sx_tol', cmdLine, exec_sx_tol, status)) RETURN 
    !
    !  No calls matched the keyword.
    !
    match = .FALSE.
    !
  END FUNCTION ConvergenceKeywordInput
  !
  ! ================================
  !
  SUBROUTINE exec_cg_max_iters(a, s)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: a
    INTEGER,          INTENT(OUT) :: s
    !
    ! ============================== Executable Code
    !
    s = 0
    READ(a, *, IOSTAT=s) cv_options % cg_max_iters
    !
  END SUBROUTINE exec_cg_max_iters
  !
  ! ================================ cg_tol
  !
  SUBROUTINE exec_cg_tol(a, s)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: a
    INTEGER,          INTENT(OUT) :: s
    !
    ! ============================== Executable Code
    !
    s = 0
    READ(a, *, IOSTAT=s) cv_options % cg_tol
    !
  END SUBROUTINE exec_cg_tol
  !
  ! ================================ nl_max_iters
  !
  SUBROUTINE exec_nl_max_iters(a, s)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: a
    INTEGER,          INTENT(OUT) :: s
    !
    ! ============================== Executable Code
    !
    s = 0
    READ(a, *, IOSTAT=s) cv_options % nl_max_iters
    !
  END SUBROUTINE exec_nl_max_iters
  !
  ! ================================ nl_tol_strict
  !
  SUBROUTINE exec_nl_tol_strict(a, s)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: a
    INTEGER,          INTENT(OUT) :: s
    !
    ! ============================== Executable Code
    !
    s = 0
    READ(a, *, IOSTAT=s) cv_options % nl_tol_strict
    !
  END SUBROUTINE exec_nl_tol_strict
  !
  ! ================================ nl_tol_loose
  !
  SUBROUTINE exec_nl_tol_loose(a, s)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: a
    INTEGER,          INTENT(OUT) :: s
    !
    ! ============================== Executable Code
    !
    s = 0
    READ(a, *, IOSTAT=s) cv_options % nl_tol_loose
    !
  END SUBROUTINE exec_nl_tol_loose
  !
  ! ================================ nl_tol_min
  !
  SUBROUTINE exec_nl_tol_min(a, s)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: a
    INTEGER,          INTENT(OUT) :: s
    !
    ! ============================== Executable Code
    !
    s = 0
    READ(a, *, IOSTAT=s) cv_options % nl_tol_min
    !
  END SUBROUTINE exec_nl_tol_min
  !
  ! ================================ nr_tol_switch_ref
  !
  SUBROUTINE exec_nr_tol_switch_ref(a, s)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: a
    INTEGER,          INTENT(OUT) :: s
    !
    ! ============================== Executable Code
    !
    s = 0
    READ(a, *, IOSTAT=s) cv_options % nr_tol_switch_ref
    !
  END SUBROUTINE exec_nr_tol_switch_ref
  !
  ! ================================ nr_tol_conv
  !
  SUBROUTINE exec_nr_tol_conv(a, s)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: a
    INTEGER,          INTENT(OUT) :: s
    !
    ! ============================== Executable Code
    !
    s = 0
    READ(a, *, IOSTAT=s) cv_options % nr_tol_conv
    !
  END SUBROUTINE exec_nr_tol_conv
  !
  ! ================================ sx_max_iters_state
  !
  SUBROUTINE exec_sx_max_iters_state(a, s)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: a
    INTEGER,          INTENT(OUT) :: s
    !
    ! ============================== Executable Code
    !
    s = 0
    READ(a, *, IOSTAT=s) cv_options % sx_max_iters_state
    !
  END SUBROUTINE exec_sx_max_iters_state
  !
  !
  ! ================================ sx_max_iters_newton
  !
  SUBROUTINE exec_sx_max_iters_newton(a, s)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: a
    INTEGER,          INTENT(OUT) :: s
    !
    ! ============================== Executable Code
    !
    s = 0
    READ(a, *, IOSTAT=s) cv_options % sx_max_iters_newton
    !
  END SUBROUTINE exec_sx_max_iters_newton
  !
  !
  ! ================================ sx_tol
  !
  SUBROUTINE exec_sx_tol(a, s)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: a
    INTEGER,          INTENT(OUT) :: s
    !
    ! ============================== Executable Code
    !
    s = 0
    READ(a, *, IOSTAT=s) cv_options % sx_tol
    !
  END SUBROUTINE exec_sx_tol
  !
  ! ====================================================   END:  ConvergenceKeywordInput
  !
END MODULE ConvergenceModule
