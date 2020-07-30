! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE OptimModule
  !
  !  Optimization routines.
  !
  !  ==================== Other Modules (USE statements)
  !
  USE IntrinsicTypesModule, RK => REAL_KIND
  USE ConstantsModule, GCONJ => RK_GOLD_CONJ

  IMPLICIT NONE
  !
  ! ==================== Public Entities
  !
  ! variables, procedures, constants, derived types and namelist groups
  !
  PRIVATE   ! all objects are private unless declared otherwise
  !
  PUBLIC :: minSearchGolden, solveNewton1D
  PUBLIC :: NEWTON_RES, NEWTON_DER
  !
  ! ==================== Module Data
  !
  ! Flags for newton solver
  !
  INTEGER, PARAMETER :: NEWTON_RES = 0, NEWTON_DER = 1
  !
CONTAINS ! ============================================= MODULE PROCEDURES
  !
  ! ==================================================== BEGIN:  minSearchGolden
  !
  SUBROUTINE minSearchGolden(f, d, interval, fmin, xmin, tolx)
    !
    ! Use golden section search to find a minimum.
    !
    ! In the future, more sophisticated convergence criteria could be specified.
    ! The current implementation only accepts an absolute tolerance on x, but it
    ! would be nice to have relative tolerances both on x and on the value.  For 
    ! that reason, the tolerance arguemnt is made optional so that other tolerance
    ! settings can be easily incorporated.
    !
    ! * ARGS *
    !
    ! f is a function
    !      the function to be minimized
    ! d is a vector
    !      data (fixed parameters) to be passed to f
    ! interval is a 2-vector
    !          giving endpoints of interval containing min; this is updated
    !          in this routine
    ! fmin is a scalar, the minimum value found
    ! xmin is a scalar, the location of the minimizing point
    ! tolx is a scalar, the absolute tolerance on x
    !
    INTERFACE
      FUNCTION f(x, d) RESULT(res)
        USE IntrinsicTypesModule, RK => REAL_KIND
        REAL(RK), INTENT(IN)  :: x, d(:)
        REAL(RK) :: res
      END FUNCTION f
    END INTERFACE
    !
    REAL(RK), INTENT(IN)  :: d(:)
    REAL(RK), INTENT(IN OUT) :: interval(2)
    REAL(RK), INTENT(OUT) :: fmin, xmin
    REAL(RK), INTENT(IN) :: tolx
    !
    OPTIONAL :: tolx
    !
    !  ========== Locals
    !
    REAL(RK) :: l_k, r_k, ak, bk, dk, fa, fb
    !
    !  ================================================== Executable Code
    !
    l_k = interval(1)
    r_k = interval(2)
    dk = r_k - l_k
    ak = r_k - GCONJ * dk
    bk = l_k + GCONJ * dk
    fa = f(ak, d)
    fb = f(bk, d)
    fmin = MIN(fa, fb)
    !
    ! Loop until tolerance is reached.
    !
    ! * Interval is reduced by a constant factor each time.
    !
    DO WHILE (dk > tolx)
      IF (fa <= fb) THEN
        ! l_k = l_k
        r_k = bk
        dk = r_k - l_k
        ak = r_k - GCONJ * dk
        bk = ak
        fb = fa
        fa = f(ak, d)
      ELSE
        l_k = ak
        ! r_k = r_k
        dk = r_k - l_k
        ak = bk
        bk = l_k + GCONJ * dk
        fa = fb
        fb = f(bk, d)
      END IF
    END DO
    !
    ! Update output variables on return.
    !
    fmin = MIN(fa, fb)
    IF (fa <= fb) THEN
      interval(1) = l_k
      xmin = ak
      interval(2) = bk
    ELSE
      interval(1) = ak
      xmin = bk
      interval(2) = r_k
    END IF
    !
  END SUBROUTINE minSearchGolden
  !
  ! ====================================================   END:  minSearchGolden
  ! ==================================================== BEGIN:  solveNewton
  !
  SUBROUTINE solveNewton1D(f, d, x, tolx, status)
    !
    ! Newton solver for nonlinear equations
    !
    ! f is a function
    !      the residual to make zero; it also returns the
    !      the derivative of the residual on request
    ! d is an array
    !      fixed parameters to pass to f and fprime
    ! x is a scalar
    !      on input, the initial guess; on ouput, the solution, if converged
    ! tolx is a scalar
    !      the solution tolerance on x
    ! status 
    !      flag indicating success(==0)|failure(/=0)
    !
    INTERFACE
      FUNCTION f(x, d, der) RESULT(res)
        !
        ! x -- the point to evaluate
        ! d -- array of fixed data
        ! der -- flag indicating return value of residual or derivative
        !
        USE IntrinsicTypesModule, RK => REAL_KIND
        REAL(RK), INTENT(IN)  :: x, d(:)
        INTEGER :: der
        REAL(RK) :: res
      END FUNCTION f

    END INTERFACE

    REAL(RK), INTENT(IN) :: d(:), tolx
    REAL(RK), INTENT(IN OUT) :: x
    INTEGER, INTENT(OUT) :: status
    !
    !  ========== Locals
    !
    ! hardwire max iterations to prevent infinite loops
    !
    INTEGER :: ITERMAX = 1000, N_INCREASED = 3
    INTEGER :: i
    REAL(RK) :: res, fprime, dx, dxold, n_inc
    !
    ! =================================================== Executable Code
    !
    status = 1  ! so that failure is set iterations reaching itermax
    !
    n_inc = 0
    dxold = HUGE(dx)
    DO i=1, ITERMAX

      res = f(x, d, NEWTON_RES)
      fprime = f(x, d, NEWTON_DER)
      IF (fprime == RK_ZERO) THEN
        ! zero derivative
        status = 2
        EXIT
      END IF
      dx = res/fprime
      x = x - dx
      dx = abs(dx)
      !
      IF (dx < tolx) THEN
        ! converged
        status = 0 
        EXIT
      END IF
      !
      ! Check for divergence
      !
      IF (dx < dxold) THEN
        n_inc = 0
      ELSE
        n_inc = n_inc + 1
        IF (n_inc >= N_INCREASED) THEN
          ! Diverged
          status = 1
          EXIT
        END IF
      END IF
      !
      dxold = dx

    END DO
    !
  END SUBROUTINE solveNewton1D
  !
  ! ====================================================   END:  solveNewton

  !
END MODULE OptimModule
