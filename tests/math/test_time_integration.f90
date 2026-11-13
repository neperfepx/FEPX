! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of time_integration_mod (step 6.4): the general forward and backward Euler on
!> state vectors that are not hardening variables.
!>   - forward Euler: x_prev + dt xdot;
!>   - backward Euler, decoupled (diagonal Jacobian): xdot_i = -k_i x_i^2, closed-form roots
!>     x = (-1 + sqrt(1 + 4 dt k x_prev))/(2 dt k), random k, x_prev, dt (large dt too: the line
!>     search);
!>   - backward Euler, lower triangular: xdot_1 = -x_1^3, xdot_i = -x_i + x_1^2 (i > 1); the root
!>     from the scalar root of component 1 and then closed forms;
!>   - backward Euler, coupled (full Jacobian, normal equations): xdot = A x, random A; the root
!>     (I - dt A)^-1 x_prev, checked by its residual;
!>   - failure: an iteration limit of 1 on a nonlinear problem with a strict tolerance: not
!>     converged, the unconverged components keep x_prev; num_iters reported;
!>   - the start estimate already within the tolerance: one Newton step is still taken, as in
!>     v2 (the residual is tested after the step): a zero rate keeps x = x_prev in one
!>     iteration; a linear equation started within the tolerance of its root ends at the root
!>     to round-off;
!>   - the solver state, allocated once (backward_euler_alloc, step 9.15, D100) and reused by
!>     every solve above: a coupled solve with the reused state equals, bit for bit, the same
!>     solve with a new state, also after the state was allocated for another size and back.
!>
!> Locals:
!>   be, be_new    type(backward_euler_t)  solver state; a new one
!>   n             integer                 size
!>   x_prev, x, k, xdot, x_ref  real(rk)(n)  states, coefficients, rates, references
!>   x_new         real(rk)(n)             the solution with the new state
!>   a_mat, jac    real(rk)(n, n)          linear-system matrix, Jacobian
!>   dt, r         real(rk)                time increment, scalar root
!>   conv          logical                 converged
!>   i, trial, its integer                 indices, iterations
program test_time_integration

  use kinds_mod, only: rk
  use time_integration_mod, only: backward_euler_t, forward_euler, backward_euler_alloc, &
    & backward_euler_start, backward_euler_step, backward_euler_finish
  use tensor_test_utils_mod, only: seed_random
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary

  implicit none

  integer, parameter :: n = 7
  type(backward_euler_t) :: be, be_new
  real(rk) :: x_new(n)
  real(rk) :: x_prev(n), x(n), k(n), xdot(n), x_ref(n), a_mat(n, n), jac(n, n), dt, r
  logical :: conv
  integer :: i, trial, its

  call seed_random(20261010)
  call backward_euler_alloc(be, n)

  ! Forward Euler.
  call random_number(x_prev)
  call random_number(xdot)
  call forward_euler(x_prev, xdot, 0.3_rk, x)
  call assert_close(x, x_prev + 0.3_rk*xdot, 'forward_euler', rel_tol=1.0e-15_rk)

  do trial = 1, 20
    ! Decoupled: xdot_i = -k_i x_i^2 (dt up to 100: large steps need the line search).
    call random_number(k)
    k = 0.1_rk + k
    call random_number(x_prev)
    x_prev = 1.0_rk + 9.0_rk*x_prev
    call random_number(dt)
    dt = 10.0_rk**(4.0_rk*dt - 2.0_rk)
    call backward_euler_start(be, x_prev, x_prev, dt, spread(1.0e-12_rk, 1, n), 100)
    do while (.not. be%done)
      xdot = -k*be%x_eval(1:n)**2
      jac = 0.0_rk
      do i = 1, n
        jac(i, i) = -2.0_rk*k(i)*be%x_eval(i)
      end do
      call backward_euler_step(be, xdot, jac)
    end do
    call backward_euler_finish(be, x, conv)
    call assert_true(conv, 'decoupled_converged')
    x_ref = (-1.0_rk + sqrt(1.0_rk + 4.0_rk*dt*k*x_prev))/(2.0_rk*dt*k)
    call assert_close(x, x_ref, 'decoupled_closed_form', rel_tol=1.0e-10_rk)

    ! Lower triangular: xdot_1 = -x_1^3, xdot_i = -x_i + x_1^2.
    call backward_euler_start(be, x_prev, x_prev, dt, spread(1.0e-12_rk, 1, n), 100)
    do while (.not. be%done)
      xdot(1) = -be%x_eval(1)**3
      xdot(2:n) = -be%x_eval(2:n) + be%x_eval(1)**2
      jac = 0.0_rk
      jac(1, 1) = -3.0_rk*be%x_eval(1)**2
      do i = 2, n
        jac(i, i) = -1.0_rk
        jac(i, 1) = 2.0_rk*be%x_eval(1)
      end do
      call backward_euler_step(be, xdot, jac)
    end do
    call backward_euler_finish(be, x, conv)
    call assert_true(conv, 'lower_converged')
    r = scalar_cubic_root(x_prev(1), dt)
    call assert_close(x(1), r, 'lower_component_1', rel_tol=1.0e-10_rk)
    call assert_close(x(2:n), (x_prev(2:n) + dt*r**2)/(1.0_rk + dt), 'lower_components', &
      & rel_tol=1.0e-10_rk)

    ! Coupled linear: xdot = A x with a full A (stable: A = -(B^T B + I)/2).
    call random_number(a_mat)
    a_mat = -0.5_rk*(matmul(transpose(a_mat), a_mat))
    do i = 1, n
      a_mat(i, i) = a_mat(i, i) - 0.5_rk
    end do
    call backward_euler_start(be, x_prev, x_prev, dt, spread(1.0e-11_rk, 1, n), 100)
    do while (.not. be%done)
      xdot = matmul(a_mat, be%x_eval(1:n))
      call backward_euler_step(be, xdot, a_mat)
    end do
    call backward_euler_finish(be, x, conv, its)
    call assert_true(conv, 'coupled_converged')
    call assert_true(maxval(abs(x - x_prev - dt*matmul(a_mat, x))) < 1.0e-11_rk, &
      & 'coupled_residual')
  end do

  ! Failure: one iteration, strict tolerance; the unconverged components keep x_prev.
  k = 1.0_rk
  x_prev = 5.0_rk
  call backward_euler_start(be, x_prev, x_prev, 1.0_rk, spread(1.0e-14_rk, 1, n), 1)
  do while (.not. be%done)
    jac = 0.0_rk
    do i = 1, n
      jac(i, i) = -2.0_rk*be%x_eval(i)
    end do
    call backward_euler_step(be, -be%x_eval(1:n)**2, jac)
  end do
  call backward_euler_finish(be, x, conv, its)
  call assert_true(.not. conv, 'failure_not_converged')
  call assert_close(x, x_prev, 'failure_reset', rel_tol=0.0_rk)
  call assert_equal(its, 1, 'failure_iterations')

  ! The start estimate already converged (zero rate): one Newton step (as v2), x unchanged.
  x_prev = 2.0_rk
  call backward_euler_start(be, x_prev, x_prev, 1.0_rk, spread(1.0e-12_rk, 1, n), 10)
  do while (.not. be%done)
    call backward_euler_step(be, spread(0.0_rk, 1, n), spread(spread(0.0_rk, 1, n), 1, n))
  end do
  call backward_euler_finish(be, x, conv, its)
  call assert_true(conv, 'start_converged')
  call assert_equal(its, 1, 'start_one_step')
  call assert_close(x, x_prev, 'start_unchanged', rel_tol=0.0_rk)

  ! Linear xdot_i = -k_i x_i started within the tolerance of the root x_prev/(1 + dt k): one
  ! Newton step ends at the root to round-off (as v2).
  call random_number(k)
  k = 0.5_rk + k
  x_prev = 3.0_rk
  dt = 0.7_rk
  x_ref = x_prev/(1.0_rk + dt*k)
  call backward_euler_start(be, x_prev, x_ref + 0.5e-3_rk, dt, spread(1.0e-2_rk, 1, n), 10)
  do while (.not. be%done)
    jac = 0.0_rk
    do i = 1, n
      jac(i, i) = -k(i)
    end do
    call backward_euler_step(be, -k*be%x_eval(1:n), jac)
  end do
  call backward_euler_finish(be, x, conv, its)
  call assert_true(conv, 'linear_converged')
  call assert_equal(its, 1, 'linear_one_step')
  call assert_close(x, x_ref, 'linear_exact', rel_tol=1.0e-14_rk)

  ! The reused state against a new one (coupled solve, the normal equations), after the state
  ! was allocated for another size and back.
  call backward_euler_alloc(be, 3)
  call backward_euler_alloc(be, n)
  call random_number(a_mat)
  a_mat = -2.0_rk*a_mat
  call random_number(x_prev)
  call coupled_solve(be, x)
  call backward_euler_alloc(be_new, n)
  call coupled_solve(be_new, x_new)
  call assert_close(x, x_new, 'reused_state_bit_for_bit', abs_tol=0.0_rk)

  call test_summary()

contains

  !> The coupled solve xdot = A x (a_mat, x_prev, dt 0.5 of the host) with a given solver state.
  !>
  !> Arguments:
  !>   state   (inout) type(backward_euler_t)  solver state (allocated for n)
  !>   x_out   (out)   real(rk)(n)             the solution
  !> Locals:
  !>   ok      logical                         converged
  subroutine coupled_solve(state, x_out)

    type(backward_euler_t), intent(inout) :: state
    real(rk), intent(out) :: x_out(n)

    logical :: ok

    call backward_euler_start(state, x_prev, x_prev, 0.5_rk, spread(1.0e-11_rk, 1, n), 100)
    do while (.not. state%done)
      call backward_euler_step(state, matmul(a_mat, state%x_eval), a_mat)
    end do
    call backward_euler_finish(state, x_out, ok)
    call assert_true(ok, 'coupled_solve converged')

  end subroutine coupled_solve

  !> Root of x - x_prev + dt x^3 = 0 (x_prev > 0) by bisection on [0, x_prev], for the
  !> lower-triangular check (independent of the solver).
  !>
  !> Arguments:
  !>   xp, dtt  (in)  real(rk)  x_prev, dt
  !> Locals:
  !>   lo, hi   real(rk)  bracket
  !>   it       integer   iteration
  function scalar_cubic_root(xp, dtt) result(x_root)

    real(rk), intent(in) :: xp, dtt
    real(rk) :: x_root

    real(rk) :: lo, hi
    integer :: it

    ! The residual is increasing, negative at 0 and positive at x_prev.
    lo = 0.0_rk
    hi = xp
    do it = 1, 200
      x_root = 0.5_rk*(lo + hi)
      if (x_root - xp + dtt*x_root**3 > 0.0_rk) then
        hi = x_root
      else
        lo = x_root
      end if
    end do
    x_root = 0.5_rk*(lo + hi)

  end function scalar_cubic_root

end program test_time_integration
