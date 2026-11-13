! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The isotropic viscoplastic start of a crystal-plasticity run (initialization
!> isotropic_viscoplastic, ROADMAP D78): the velocity guess of increment 1 from a linear,
!> incompressible, isotropic viscous problem on the initial configuration, solved by v2's three
!> augmented-Lagrangian (Uzawa) iterations with an element-constant pressure.
!>
!> Theory : dev/theory/addendum_isotropic_viscoplastic.md §1 (approved):
!>            {sigma'} = 2 mu [P_dev] {d},  mu = 1 MPa s                                  (I1)
!>            {g_e} = - sum_q w_q J_q [B_q]^T {m}  ({g_e}.{V_e} = - int_e div v dV)          (I2)
!>            V_e = sum_q w_q J_q,  r_e = mu / (6 V_e)                                     (I3)
!>            [k_e] = sum_q w_q J_q [B_q]^T 2 mu [P_dev] [B_q] + r_e {g_e}{g_e}^T            (I4)
!>            [K]{V}^k = - sum_e p_e^(k-1) {g_e} on the free dofs                          (I5)
!>            p_e^k = p_e^(k-1) + r_e {g_e}.{V_e}^k,  p_e^0 = 0,  k = 1, 2, 3               (I6)
!>          the stationarity conditions of the augmented Lagrangian of the incompressible
!>          viscous problem and its multiplier update (Fortin and Glowinski 1983; Glowinski and
!>          Le Tallec 1989). The velocity does not depend on mu (r_e is proportional to it).
!> Source : port of v2-fixed solveit_isovp_mod1.f90::solveit_vp (three iterations, no
!>          convergence test, the norm of the change logged, a NaN norm stops the run),
!>          solveit_isovp_mod2.f90::elt_stif_vp (pcnst = {g_e}, pscale = r_e, the penalty
!>          matrix pscale pcnst pcnst^T) and ::recover_pressure_vp (I6),
!>          solveit_isovp_mod3.f90::material_matrix_vp and solveit_isovp_mod4.f90::isotropic (the
!>          viscosity 1: v2's Hart model after its return is dead code, not ported); v2 writes
!>          the deviatoric law in its 5-vector basis (d11 - d22, d'33, d12, d13, d23) with
!>          mu diag(1, 3, 4, 4, 4), the same energy 2 mu d':d' as (I1) (tested).
!>          Deviations (addendum §2, §3; user decisions Q-I1..Q-I3, 2026-10-21): the solve is
!>          constrained by the multi-point constraints (v2 solves unconstrained, then projects);
!>          it runs for every kind of boundary condition (v2: the uniaxial driver only);
!>          r_e = mu / (6 V_e) for every element type (v2: (sum_q w_q) mu / V_e, 1/6 for its
!>          tetrahedron; hex20 is FEPX's); the element matrices and {g_e} are formed once (v2
!>          forms them at every iteration, with the same values: the problem is linear and the
!>          geometry fixed); the solve is in the residual form {r} = {f} - [K]{V},
!>          {V} <- {V} + {dV} (v2's CG mode "v" solves for {V} with the prescribed values held:
!>          the same system).
!>
!> Policy (layout of increment_mod): velocities are local dof vectors (owned nodes first, then
!> ghosts), valid on the ghosts; each element belongs to one rank, so the pressures are local.
!>
!> Public routines:
!>   iso_material_matrix   2 mu [P_dev], Mandel (I1) (pure)
!>   iso_element_system    [k_e], {g_e} and r_e of one element (I2)-(I4) (pure)
!>   isotropic_start_solve the three iterations (I5)-(I6) (collective)
module isotropic_start_mod

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use kinds_mod, only: rk, ik
  use element_mod, only: element_t, element_geometry_qpt, mandel_b_matrix
  use element_matrix_mod, only: element_work_t, element_work_alloc, element_stiffness
  use mesh_mod, only: mesh_t
  use halo_mod, only: dist_norm
  use error_mod, only: fatal_error
  use log_mod, only: log_info, log_int_text, log_real_text, log_count_text
  use ebe_mod, only: ebe_matvec, ebe_diagonal, ebe_assemble
  use pcg_mod, only: pcg_solve
  use constraint_mod, only: constraint_t, constraint_diagonal, constraint_diagonal_lumped, &
    & constraint_zero_mean

  implicit none

  private

  public :: iso_material_matrix, iso_element_system, isotropic_start_solve

  ! The viscosity of v2's isotropic problem (solveit_isovp_mod4.f90::isotropic) [MPa s].
  real(rk), parameter, public :: iso_viscosity = 1.0_rk
  ! The number of augmented-Lagrangian iterations (v2 solveit_vp).
  integer, parameter, public :: iso_num_iters = 3
  ! The penalty factor: r_e = iso_penalty_factor mu / V_e (v2: the sum of its tetrahedron
  ! weights, 1/6; Q-I3: the same for every element type) [-].
  real(rk), parameter, public :: iso_penalty_factor = 1.0_rk/6.0_rk

contains

  !> The isotropic viscous matrix 2 mu [P_dev] in Mandel form (I1): {sigma'} = 2 mu [P_dev] {d}.
  !>
  !> Arguments:
  !>   mu     (in)  real(rk)        viscosity [MPa s]
  !>   mat_s  (out) real(rk)(6, 6)  the matrix, sample frame [MPa s]
  !> Locals:
  !>   i, j   integer               rows and columns
  pure subroutine iso_material_matrix(mu, mat_s)

    real(rk), intent(in) :: mu
    real(rk), intent(out) :: mat_s(6, 6)

    integer :: i, j

    ! The identity minus (1/3) {m}{m}^T on the normal components, times 2 mu.
    mat_s = 0.0_rk
    do i = 1, 6
      mat_s(i, i) = 2.0_rk*mu
    end do
    do j = 1, 3
      do i = 1, 3
        mat_s(i, j) = mat_s(i, j) - 2.0_rk*mu/3.0_rk
      end do
    end do

  end subroutine iso_material_matrix

  !> Element matrix (I4), pressure vector (I2) and penalty (I3) of one element.
  !>
  !> Arguments:
  !>   elt      (in)  type(element_t)  element descriptor
  !>   coo      (in)  real(rk)(:, :)   (3, num_nodes) node coordinates, sample frame [mm]
  !>   mu       (in)  real(rk)         viscosity [MPa s]
  !>   k_e      (out) real(rk)(:, :)   (num_dofs, num_dofs) element matrix [N s/mm]
  !>   g_e      (out) real(rk)(:)      (num_dofs) {g_e}: {g_e}.{V_e} = - int_e div v dV [mm^2]
  !>   r_e      (out) real(rk)         penalty mu / (6 V_e) [MPa s/mm^3]
  !>   work     (inout) type(element_work_t)  element work arrays (element_work_alloc): also
  !>                                   dshape_s, d N_a / d x_i [1/mm], and b, the Mandel B
  !>                                   matrix [1/mm], of the pressure vector
  !> Locals:
  !>   mat_s     real(rk)(6, 6, num_qpts)  2 mu [P_dev] at every point [MPa s]
  !>   i_qpt     integer(ik)               quadrature point
  !>   det_jac   real(rk)                  det J [mm^3]
  !>   vol       real(rk)                  element volume V_e [mm^3]
  !>   i, j      integer(ik)               element dofs
  pure subroutine iso_element_system(elt, coo, mu, k_e, g_e, r_e, work)

    type(element_t), intent(in) :: elt
    real(rk), intent(in) :: coo(:, :)
    real(rk), intent(in) :: mu
    real(rk), intent(out) :: k_e(:, :)
    real(rk), intent(out) :: g_e(:)
    real(rk), intent(out) :: r_e
    type(element_work_t), intent(inout) :: work

    real(rk) :: mat_s(6, 6, elt%num_qpts)
    integer(ik) :: i_qpt
    real(rk) :: det_jac
    real(rk) :: vol
    integer(ik) :: i, j

    ! Deviatoric part: the element matrix of the material matrix 2 mu [P_dev] (R14 form).
    do i_qpt = 1_ik, elt%num_qpts
      call iso_material_matrix(mu, mat_s(:, :, i_qpt))
    end do
    call element_stiffness(elt, coo, mat_s, k_e, work)

    ! The pressure vector {g_e} = - sum_q w_q J_q [B_q]^T {m} (the trace rows of B) and the
    ! volume.
    g_e = 0.0_rk
    vol = 0.0_rk
    do i_qpt = 1_ik, elt%num_qpts
      call element_geometry_qpt(elt, coo, i_qpt, work%dshape_s, det_jac)
      call mandel_b_matrix(work%dshape_s, work%b)
      g_e = g_e - (elt%qpt_weight(i_qpt)*det_jac)*(work%b(1, :) + work%b(2, :) + work%b(3, :))
      vol = vol + elt%qpt_weight(i_qpt)*det_jac
    end do

    ! The penalty and its matrix r_e {g_e}{g_e}^T.
    r_e = iso_penalty_factor*mu/vol
    do j = 1_ik, elt%num_dofs
      do i = 1_ik, elt%num_dofs
        k_e(i, j) = k_e(i, j) + r_e*g_e(i)*g_e(j)
      end do
    end do

  end subroutine iso_element_system

  !> The three augmented-Lagrangian iterations (I5)-(I6) from the velocity guess, which holds the
  !> prescribed values (and satisfies the constraints, when present). Collective.
  !>
  !> Arguments:
  !>   mesh          (in)    type(mesh_t)       rank-local mesh
  !>   elt           (in)    type(element_t)    element descriptor
  !>   coo           (in)    real(rk)(:, :)     (3, num_nodes_loc) node coordinates [mm]
  !>   bc_mask       (in)    logical(:)         (3 num_nodes_loc) prescribed dofs
  !>   cg_tol        (in)    real(rk)           CG relative tolerance (the run's cg_tol) [-]
  !>   cg_max_iters  (in)    integer            CG iteration limit (the run's cg_max_iters)
  !>   vel           (inout) real(rk)(:)        (3 num_nodes_loc) the guess on entry, the
  !>                                            velocity of the third iteration on exit [mm/s]
  !>   cons          (in, optional) type(constraint_t)  constraints (used when active)
  !>   cg_precond    (in, optional) integer  CG preconditioner (pcg_precond_*, default Jacobi;
  !>                                         step 9.2)
  !> Locals:
  !>   use_cons      logical             constraints present and active
  !>   n_owned       integer(ik)         owned dofs
  !>   k_elt         real(rk)(:, :, :)   element matrices (I4)
  !>   g_elt         real(rk)(:, :)      element pressure vectors (I2)
  !>   r_elt         real(rk)(:)         element penalties (I3)
  !>   p_elt         real(rk)(:)         element pressures p_e [MPa]
  !>   f_elt         real(rk)(:, :)      element forces - p_e {g_e} [N]
  !>   f, k_vel, resid, dvel, diag, vel_prev  real(rk)(:)  global vectors
  !>   no_mask       logical(:)          no prescribed dofs (the product over all dofs)
  !>   iter, cg_iters, cg_status  integer   iteration, CG iterations, CG status
  !>   i_elt         integer(ik)         element
  !>   dv_norm       real(rk)            norm of the change of the iteration [mm/s]
  !>   ework         type(element_work_t)  element work arrays (D100)
  subroutine isotropic_start_solve(mesh, elt, coo, bc_mask, cg_tol, cg_max_iters, vel, cons, &
    & cg_precond)

    type(mesh_t), intent(in) :: mesh
    type(element_t), intent(in) :: elt
    real(rk), intent(in) :: coo(:, :)
    logical, intent(in) :: bc_mask(:)
    real(rk), intent(in) :: cg_tol
    integer, intent(in) :: cg_max_iters
    real(rk), intent(inout) :: vel(:)
    type(constraint_t), intent(in), optional :: cons
    integer, intent(in), optional :: cg_precond

    logical :: use_cons
    integer(ik) :: n_owned
    real(rk), allocatable :: k_elt(:, :, :), g_elt(:, :), r_elt(:), p_elt(:), f_elt(:, :)
    real(rk), allocatable :: f(:), k_vel(:), resid(:), dvel(:), diag(:), vel_prev(:)
    logical, allocatable :: no_mask(:)
    integer :: iter, cg_iters, cg_status
    integer(ik) :: i_elt
    real(rk) :: dv_norm
    type(element_work_t) :: ework

    n_owned = 3_ik*mesh%num_nodes_owned
    use_cons = .false.
    if (present(cons)) use_cons = cons%active
    allocate (k_elt(elt%num_dofs, elt%num_dofs, mesh%num_elts), &
      & g_elt(elt%num_dofs, mesh%num_elts), r_elt(mesh%num_elts), p_elt(mesh%num_elts), &
      & f_elt(elt%num_dofs, mesh%num_elts))
    allocate (f(size(vel)), k_vel(size(vel)), resid(size(vel)), dvel(size(vel)), &
      & diag(size(vel)), vel_prev(size(vel)), no_mask(size(vel)))
    no_mask = .false.

    ! Element matrices, pressure vectors and penalties, once (linear problem, fixed geometry).
    call element_work_alloc(ework, elt)
    do i_elt = 1_ik, mesh%num_elts
      call iso_element_system(elt, coo(:, mesh%elt_nodes_loc(:, i_elt)), iso_viscosity, &
        & k_elt(:, :, i_elt), g_elt(:, i_elt), r_elt(i_elt), ework)
    end do

    ! The preconditioner: the assembled diagonal (constrained: the exact reduced one or the
    ! lumped periodic one, as increment_mod).
    if (use_cons) then
      if (cons%diag_exact) then
        call constraint_diagonal(cons, mesh, k_elt, diag)
      else
        call constraint_diagonal_lumped(cons, mesh, k_elt, diag)
      end if
    else
      call ebe_diagonal(mesh, k_elt, diag)
    end if

    ! The iterations (v2: three, no convergence test).
    p_elt = 0.0_rk
    do iter = 1, iso_num_iters
      ! Forces of the current pressures, - p_e {g_e} (v2 eforce = - elpress pcnst), and the
      ! residual {f} - [K]{V}, zero at the prescribed dofs.
      do i_elt = 1_ik, mesh%num_elts
        f_elt(:, i_elt) = -p_elt(i_elt)*g_elt(:, i_elt)
      end do
      call ebe_assemble(mesh, f_elt, f)
      call ebe_matvec(mesh, k_elt, no_mask, vel, k_vel)
      resid = f - k_vel
      where (bc_mask) resid = 0.0_rk

      ! The change: [K]{dV} = {r} on the free dofs (constrained: the reduced system).
      dvel = 0.0_rk
      if (use_cons) then
        call pcg_solve(mesh, k_elt, bc_mask, diag, resid, dvel, cg_tol, cg_max_iters, &
          & cg_iters, cg_status, cons)
      else
        call pcg_solve(mesh, k_elt, bc_mask, diag, resid, dvel, cg_tol, cg_max_iters, &
          & cg_iters, cg_status, precond=cg_precond)
      end if
      if (cg_status == 2) call fatal_error('isotropic_start_solve', 'initialization ' // &
        & 'isotropic_viscoplastic: a nodal block of the block-Jacobi preconditioner is not ' // &
        & 'positive definite')
      if (cg_status /= 0) call fatal_error('isotropic_start_solve', 'initialization ' // &
        & 'isotropic_viscoplastic: the conjugate gradient reached cg_max_iters')
      vel_prev = vel
      vel = vel + dvel
      ! Periodic conditions: the rigid translation fixed by the zero-mean fluctuation (Q-P3).
      if (use_cons) then
        if (cons%zero_mean) call constraint_zero_mean(cons, mesh, vel)
      end if
      dvel = vel - vel_prev

      ! The pressure update (I6) (v2 recover_pressure_vp).
      do i_elt = 1_ik, mesh%num_elts
        p_elt(i_elt) = p_elt(i_elt) + r_elt(i_elt)*dot_product(g_elt(:, i_elt), &
          & vel(mesh%elt_dofs_loc(:, i_elt)))
      end do

      ! The norm of the change (v2: absolute, all dofs; a NaN stops the run).
      dv_norm = dist_norm(dvel, n_owned)
      if (.not. ieee_is_finite(dv_norm)) call fatal_error('isotropic_start_solve', &
        & 'initialization isotropic_viscoplastic: the velocity change is not finite')
      call log_info('Iteration ' // log_int_text(iter) // ': |dV| = ' // &
        & log_real_text(dv_norm, 4) // ' mm/s (' // log_count_text(cg_iters, 'CG iteration') // &
        & ')', 2)
    end do

  end subroutine isotropic_start_solve

end module isotropic_start_mod
