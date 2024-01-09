! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module exec_type_mod

! Module to define types

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use types_mod, only : exec_type

  implicit none

  public

  contains

  subroutine exec_set_default(exec)

    type(exec_type), intent(out) :: exec

    exec%restart = .false.
    ! Convergence options
    ! Conjugate gradient values
    exec%check_necking = .false.
    exec%max_iter_hard_limit = 10
    exec%cg_max_iters = 16000
    exec%cg_tol = 1.0d-8
    ! Global nonlinear iteration
    exec%nl_max_iters = 50
    exec%nl_tol_strict = 5.0d-4
    exec%nl_tol_loose = 5.0d-4
    exec%nl_tol_min = 1.0d-10
    ! Newton-Raphson
    exec%nr_tol_switch_ref = 1.0d-2
    exec%nr_tol_conv= 0.2d0
    ! Single crystal / state calculations
    exec%sx_max_iters_state = 100
    exec%sx_max_iters_newton = 100
    exec%sx_tol = 1.0d-4
    exec%min_pert_frac = 1.0d-3

  end subroutine exec_set_default

end module exec_type_mod
