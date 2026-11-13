! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Time integration of a state vector over one increment, as a general numerical method (step
!> 6.4): x is any vector of state variables whose rates xdot(x) and Jacobian d(xdot)/dx a caller
!> can evaluate (the slip-system strengths of the saturation hardening laws; for a later law, the
!> variable it evolves, e.g. a dislocation density, from which it then computes the strengths).
!> Nothing here knows what x is.
!>
!>   forward Euler   x = x_prev + dt xdot(x_prev)
!>   backward Euler  r(x) = x - x_prev - dt xdot(x) = 0, by Newton's method with the full
!>                   Jacobian dr/dx = I - dt d(xdot)/dx and a line search halving the step.
!>
!> The backward-Euler solver uses reverse communication (as ODEPACK and ARPACK do): the caller
!> evaluates xdot and d(xdot)/dx at the point the solver asks for (backward_euler_t%x_eval) and
!> hands them back, until the solver is done. The rate evaluation therefore stays in the
!> caller's pure kernel, with no procedure argument and no class dispatch per point
!> (dev/conventions.md: polymorphism only at model selection). The solver state holds its
!> arrays and the work arrays of a step at the exact size of the state vector, allocated once by
!> backward_euler_alloc outside the per-point loops (no allocation per point, ROADMAP D100).
!> Use:
!>
!>     call backward_euler_alloc(be, n)                  ! once, outside the point loop
!>     ...
!>     call backward_euler_start(be, x_prev, x_start, dt, toler, max_iters)
!>     do while (.not. be%done)
!>       call <rates>(be%x_eval, xdot, jac)              ! the caller's law
!>       call backward_euler_step(be, xdot, jac)
!>     end do
!>     call backward_euler_finish(be, x, converged, num_iters)
!>
!> Source : port of v2-fixed hardening_mod1.f90::update_hardening (the backward-Euler Newton
!>          iteration of the strengths: a step halved while the residual of a component grows,
!>          failure below ls_cutoff = 0.001; convergence when |r_i| < toler_i for every
!>          component; at most max_iters iterations; components that did not converge keep
!>          x_prev), written for any vector (addendum_hardening.md §5). Deviations (addendum §5):
!>          Newton on the vector with the full Jacobian (v2: a scalar Newton per component).
!>          When the Jacobian is diagonal (decoupled equations) the components keep v2's
!>          separate line-search factors; otherwise one factor halves the whole step while the
!>          norm of r_i/toler_i grows (a per-component test would wrongly shorten a component
!>          whose residual grows only because a component it depends on has not converged yet).
!>          A lower-triangular Newton system is solved by forward substitution, any other by
!>          the normal equations (linalg_mod). As v2's update_hardening, each iteration takes
!>          the Newton step and then tests the residual after it, so the start estimate is
!>          never accepted without a step.
!>
!> Public types:
!>   backward_euler_t    the backward-Euler solver state (components at the type)
!>
!> Public routines:
!>   forward_euler           x = x_prev + dt xdot
!>   backward_euler_alloc    allocate a solver state for a state vector of a given size
!>   backward_euler_start    start the backward-Euler solve
!>   backward_euler_step     advance it with xdot and d(xdot)/dx at x_eval
!>   backward_euler_finish   the result (non-converged components reset to x_prev)
module time_integration_mod

  use kinds_mod, only: rk
  use linalg_mod, only: normal_equations_solve

  implicit none

  private

  public :: backward_euler_t, forward_euler, backward_euler_alloc, backward_euler_start, &
    & backward_euler_step, backward_euler_finish

  ! v2's smallest line-search factor (ls_cutoff).
  real(rk), parameter :: ls_cutoff = 1.0e-3_rk

  ! Solver stages: rates wanted at the iterate (residual and Newton step), or at a trial point.
  integer, parameter :: stage_iterate = 1, stage_trial = 2

  !> Backward-Euler solver state, every array of size n (backward_euler_alloc); only n,
  !> x_eval, done are read by the caller.
  !>
  !> Components:
  !>   n          integer          size of the state vector
  !>   x_prev     real(rk)(n)      state at the start of the increment
  !>   x          real(rk)(n)      current iterate
  !>   x_eval     real(rk)(n)      point at which the caller evaluates xdot and d(xdot)/dx next
  !>   toler      real(rk)(n)      convergence tolerance per component, |r_i| < toler_i
  !>   residual   real(rk)(n)      residual at the iterate
  !>   delta      real(rk)(n)      Newton step
  !>   fact       real(rk)(n)      line-search factors (all equal unless the Jacobian is
  !>                               diagonal)
  !>   conv, ok   logical(n)       per component: converged; Newton usable
  !>   mat        real(rk)(n, n)   work: the Newton matrix I - dt jac
  !>   mat_tm     real(rk)(n, n)   work: mat^T mat of the normal equations, then its factors
  !>   rhs        real(rk)(n)      work: the Newton right-hand side -r
  !>   residual_try  real(rk)(n)   work: the residual at a trial point
  !>   grew       logical(n)       work: components whose line search continues
  !>   dt         real(rk)         time increment
  !>   max_iters  integer          Newton iteration limit
  !>   iter       integer          Newton iterations performed
  !>   stage      integer          stage_iterate or stage_trial
  !>   lower      logical          the current Jacobian is lower triangular
  !>   diagonal   logical          the current Jacobian is diagonal (per-component line search)
  !>   done       logical          the solve has ended (call backward_euler_finish)
  type, public :: backward_euler_t
    integer :: n = 0
    real(rk), allocatable :: x_prev(:)
    real(rk), allocatable :: x(:)
    real(rk), allocatable :: x_eval(:)
    real(rk), allocatable :: toler(:)
    real(rk), allocatable :: residual(:)
    real(rk), allocatable :: delta(:)
    real(rk), allocatable :: fact(:)
    logical, allocatable :: conv(:)
    logical, allocatable :: ok(:)
    real(rk), allocatable :: mat(:, :)
    real(rk), allocatable :: mat_tm(:, :)
    real(rk), allocatable :: rhs(:)
    real(rk), allocatable :: residual_try(:)
    logical, allocatable :: grew(:)
    real(rk) :: dt = 0.0_rk
    integer :: max_iters = 0
    integer :: iter = 0
    integer :: stage = stage_iterate
    logical :: lower = .true.
    logical :: diagonal = .true.
    logical :: done = .false.
  end type backward_euler_t

contains

  !> Forward Euler: x = x_prev + dt xdot (xdot evaluated by the caller at x_prev).
  !>
  !> Arguments:
  !>   x_prev  (in)  real(rk)(n)  state at the start of the increment
  !>   xdot    (in)  real(rk)(n)  its rate at x_prev
  !>   dt      (in)  real(rk)     time increment
  !>   x       (out) real(rk)(n)  state at the end of the increment
  pure subroutine forward_euler(x_prev, xdot, dt, x)

    real(rk), intent(in) :: x_prev(:)
    real(rk), intent(in) :: xdot(:)
    real(rk), intent(in) :: dt
    real(rk), intent(out) :: x(:)

    x = x_prev + dt*xdot

  end subroutine forward_euler

  !> Allocate a solver state for a state vector of size n (its arrays and the work arrays of a
  !> step), once, outside the loops that solve (ROADMAP D100); a state already of size n is
  !> kept as it is.
  !>
  !> Arguments:
  !>   be  (inout) type(backward_euler_t)  solver state
  !>   n   (in)    integer                 size of the state vector (>= 1)
  pure subroutine backward_euler_alloc(be, n)

    type(backward_euler_t), intent(inout) :: be
    integer, intent(in) :: n

    ! Nothing to do when the size is already n.
    if (be%n == n .and. allocated(be%x)) return

    ! Every array at the new size.
    if (allocated(be%x)) then
      deallocate (be%x_prev, be%x, be%x_eval, be%toler, be%residual, be%delta, be%fact, &
        & be%conv, be%ok, be%mat, be%mat_tm, be%rhs, be%residual_try, be%grew)
    end if
    be%n = n
    allocate (be%x_prev(n), be%x(n), be%x_eval(n), be%toler(n), be%residual(n), be%delta(n), &
      & be%fact(n), be%conv(n), be%ok(n), be%mat(n, n), be%mat_tm(n, n), be%rhs(n), &
      & be%residual_try(n), be%grew(n))

  end subroutine backward_euler_alloc

  !> Start a backward-Euler solve: the first point the caller evaluates is the start estimate.
  !> Every component of the state is set (as for a new state).
  !>
  !> Arguments:
  !>   be         (inout) type(backward_euler_t)  solver state (backward_euler_alloc with the
  !>                                              size of x_prev)
  !>   x_prev     (in)    real(rk)(n)  state at the start of the increment
  !>   x_start    (in)    real(rk)(n)  start estimate (e.g. forward Euler, or the last solution)
  !>   dt         (in)    real(rk)     time increment
  !>   toler      (in)    real(rk)(n)  convergence tolerance per component, |r_i| < toler_i
  !>   max_iters  (in)    integer      Newton iteration limit
  pure subroutine backward_euler_start(be, x_prev, x_start, dt, toler, max_iters)

    type(backward_euler_t), intent(inout) :: be
    real(rk), intent(in) :: x_prev(:)
    real(rk), intent(in) :: x_start(:)
    real(rk), intent(in) :: dt
    real(rk), intent(in) :: toler(:)
    integer, intent(in) :: max_iters

    ! Size check (the state is allocated once, outside the point loops).
    if (.not. allocated(be%x) .or. be%n /= size(x_prev)) error stop &
      & 'backward_euler_start: the solver state is not allocated for this size'

    ! Initial state: rates wanted at the start estimate; the rest as in a new state.
    be%x_prev(:) = x_prev
    be%x(:) = x_start
    be%x_eval(:) = x_start
    be%toler(:) = toler
    be%residual(:) = 0.0_rk
    be%delta(:) = 0.0_rk
    be%fact(:) = 1.0_rk
    be%conv(:) = .false.
    be%dt = dt
    be%max_iters = max_iters
    be%iter = 0
    be%stage = stage_iterate
    be%lower = .true.
    be%diagonal = .true.
    be%ok(:) = .true.
    be%done = .false.

  end subroutine backward_euler_start

  !> Advance the backward-Euler solve with xdot and d(xdot)/dx evaluated at be%x_eval.
  !>   stage_iterate (x_eval = iterate): residual and convergence; stop when every component has
  !>     converged, one cannot continue, or the iteration limit is reached; else a Newton step
  !>     (I - dt jac) delta = -r and the full step as the first trial point. At the start
  !>     estimate no component counts as converged: as v2, the convergence test applies to the
  !>     residual after a Newton step.
  !>   stage_trial (x_eval = trial point): line search. Diagonal Jacobian: a component whose
  !>     |r_i| grew halves its factor (failure below ls_cutoff) and the trial is repeated;
  !>     otherwise the whole step is halved while the norm of r_i/toler_i grew. Otherwise the
  !>     unconverged, usable components take their steps and the new iterate is evaluated.
  !>
  !> Arguments:
  !>   be    (inout) type(backward_euler_t)  solver state
  !>   xdot  (in)    real(rk)(n)     rates at be%x_eval
  !>   jac   (in)    real(rk)(n, n)  d(xdot)/dx at be%x_eval (used at stage_iterate)
  !> Locals:
  !>   n             integer          size
  !>   a, ios        integer          diagonal index, linear-solve status
  !>   (the Newton matrix, its right-hand side, the trial residual and the line-search mask are
  !>   be's work arrays mat, rhs, residual_try, grew)
  pure subroutine backward_euler_step(be, xdot, jac)

    type(backward_euler_t), intent(inout) :: be
    real(rk), intent(in) :: xdot(:)
    real(rk), intent(in) :: jac(:, :)

    integer :: n, a, ios

    n = be%n
    select case (be%stage)

    case (stage_iterate)
      ! Residual at the iterate and convergence of each component.
      be%residual(:) = be%x - be%x_prev - be%dt*xdot
      be%conv(:) = abs(be%residual) < be%toler
      ! v2 tests the residual after a Newton step: nothing has converged at the start
      ! estimate.
      if (be%iter == 0) be%conv(:) = .false.
      if (all(be%conv) .or. .not. all(be%ok) .or. be%iter >= be%max_iters) then
        be%done = .true.
        return
      end if
      be%iter = be%iter + 1

      ! Newton step: (I - dt jac) delta = -r.
      be%mat(:, :) = -be%dt*jac
      do a = 1, n
        be%mat(a, a) = be%mat(a, a) + 1.0_rk
      end do
      be%lower = is_lower_triangular(be%mat)
      be%diagonal = be%lower .and. is_lower_triangular(transpose(be%mat))
      be%rhs(:) = -be%residual
      call solve_newton(be%mat, be%rhs, be%lower, be%delta, be%mat_tm, ios)
      if (ios /= 0) then
        be%ok(:) = .false.
        be%done = .true.
        return
      end if

      ! First trial point: the full step.
      be%fact(:) = 1.0_rk
      be%x_eval(:) = be%x + be%delta
      be%stage = stage_trial

    case (stage_trial)
      ! Line search: which components (or the whole step) must be shortened.
      be%residual_try(:) = be%x_eval - be%x_prev - be%dt*xdot
      if (be%diagonal) then
        be%grew(:) = abs(be%residual_try) > abs(be%residual) .and. be%ok .and. .not. be%conv
      else
        be%grew(:) = .not. (norm2(be%residual_try/be%toler) <= norm2(be%residual/be%toler)) &
          & .and. all(be%ok)
      end if

      ! Shorten and try again (a component below ls_cutoff fails).
      if (any(be%grew)) then
        where (be%grew) be%fact(:) = 0.5_rk*be%fact
        where (be%fact < ls_cutoff) be%ok(:) = .false.
        be%x_eval(:) = be%x + be%fact*be%delta
        return
      end if

      ! Accept the steps of the usable, unconverged components (v2); evaluate the new iterate.
      where (be%ok .and. .not. be%conv) be%x(:) = be%x + be%fact*be%delta
      be%x_eval(:) = be%x
      be%stage = stage_iterate

    end select

  end subroutine backward_euler_step

  !> The backward-Euler result: components that did not converge keep x_prev (v2).
  !>
  !> Arguments:
  !>   be         (in)  type(backward_euler_t)  solver state (be%done)
  !>   x          (out) real(rk)(n)  state at the end of the increment
  !>   converged  (out) logical      every component converged
  !>   num_iters  (out, optional) integer  Newton iterations performed
  pure subroutine backward_euler_finish(be, x, converged, num_iters)

    type(backward_euler_t), intent(in) :: be
    real(rk), intent(out) :: x(:)
    logical, intent(out) :: converged
    integer, intent(out), optional :: num_iters

    converged = all(be%conv)
    x = merge(be%x, be%x_prev, be%conv)
    if (present(num_iters)) num_iters = be%iter

  end subroutine backward_euler_finish

  !> Whether a square matrix is lower triangular (no nonzero entry above the diagonal).
  !>
  !> Arguments:
  !>   mat    (in)  real(rk)(:, :)  the matrix
  !>   lower  (out, result) logical
  !> Locals:
  !>   b      integer  column
  pure function is_lower_triangular(mat) result(lower)

    real(rk), intent(in) :: mat(:, :)
    logical :: lower

    integer :: b

    lower = .true.
    do b = 2, size(mat, 2)
      if (any(abs(mat(1:b - 1, b)) > 0.0_rk)) lower = .false.
    end do

  end function is_lower_triangular

  !> Solve the Newton system mat x = rhs: forward substitution when mat is lower triangular,
  !> else the normal equations (linalg_mod).
  !>
  !> Arguments:
  !>   mat     (in)  real(rk)(n, n)  the matrix
  !>   rhs     (in)  real(rk)(n)     right-hand side
  !>   lower   (in)  logical         mat is lower triangular
  !>   x       (out) real(rk)(n)     solution
  !>   mat_tm  (out) real(rk)(n, n)  work array of the normal equations (mat^T mat, then its
  !>                                 factors)
  !>   status  (out) integer         0, or nonzero for a singular matrix
  !> Locals:
  !>   a       integer               row
  pure subroutine solve_newton(mat, rhs, lower, x, mat_tm, status)

    real(rk), intent(in) :: mat(:, :)
    real(rk), intent(in) :: rhs(:)
    logical, intent(in) :: lower
    real(rk), intent(out) :: x(:)
    real(rk), intent(out) :: mat_tm(:, :)
    integer, intent(out) :: status

    integer :: a

    status = 0
    if (lower) then
      ! Forward substitution; a zero pivot is singular.
      do a = 1, size(rhs)
        if (.not. (abs(mat(a, a)) > 0.0_rk)) then
          status = 1
          return
        end if
        x(a) = (rhs(a) - dot_product(mat(a, 1:a - 1), x(1:a - 1)))/mat(a, a)
      end do
    else
      call normal_equations_solve(mat, rhs, x, mat_tm, status)
    end if

  end subroutine solve_newton

end module time_integration_mod
