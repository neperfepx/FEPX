! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> One increment of the recoupled formulation: the global iteration for the end-of-increment
!> velocities (successive approximation and Newton-Raphson, D43, step 6.3), then the state at
!> the converged velocities.
!>
!> Theory : recoupled_formulation.md §6-7: at each iteration the kinematics are evaluated in the
!>          end-of-increment configuration x = x0 + dt v (v2's choice), the material returns
!>          [s] and {h} (R10)-(R12), the element matrices are (R14), and the velocity change
!>          solves [K]{dV} = {F_h} - [K]{V} on the free dofs.
!> Source : port of v2 solveit_evp_mod1.f90::solveit_evp (successive approximation: step 4.4,
!>          D43; Newton-Raphson and the switching: step 6.3): the residual form ("dv" solve);
!>          the norms (velocity change relative to the previous velocities, Euclidean over all
!>          dofs); SA iterations solve [K]{dV} = {F_h} - [K]{V}, NR iterations the same secant
!>          residual with the tangent matrices [K_tan] (from [s_tan], §4); the switch from SA to
!>          NR when the relative change falls below the switch tolerance (global_nr_switch at
!>          the start of the increment), saving the SA velocities; the NR rate (log of the
!>          ratio of consecutive changes from NR iteration global_nr_slope_start = 2 on); and
!>          v2's cases, in v2's order:
!>            1  relative change < global_tol, iteration > 1: converged;
!>            2  absolute change < global_tol_abs 3 num_nodes, iteration > 1: converged;
!>            3  NR slow (rate below global_nr_slow_ratio times the reference rate) and
!>               relative change < global_tol_slow: converged;
!>            4  SA after an NR attempt and relative change < global_tol_slow: converged;
!>            5  NR diverging (residual norm above 1.1 times that of iteration 1, prescribed
!>               dofs included, as v2): back to the saved SA velocities and to SA, switch
!>               tolerance divided by 10;
!>            6  NR slow: back to SA for the rest of the increment;
!>            7  SA change growing 5 times in a row: failure;
!>            8  (FEPX, D78) a residual or velocity-change norm that is not finite (NaN or
!>               infinite, e.g. from failed point solves): failure (diverging) before the tests
!>               above, which a NaN would pass (its comparisons are false; v2 has no such
!>               test);
!>          the iteration limit global_max_iters; v2's defaults (global_max_iters 50,
!>          global_tol 5e-4, global_tol_abs 1e-10, global_tol_slow 5e-4, global_nr_switch
!>          1e-2, global_nr_slow_ratio 0.2: exec_type_mod.f90; global_nr_slope_start 2:
!>          types_mod.f90); then the state at the converged velocities (v2 update_state_evp).
!>          With the elastic model [s_tan] = [s], so NR iterations equal SA iterations.
!>          Mid-node linearization (D26) is applied by the caller to the start-of-increment
!>          coordinates (v2 resets the mid-nodes at every iteration, which leaves them unchanged
!>          after the first).
!>          Constraints v = T v_hat + g (step 8.3a, MPCs; v2 solveit_evp with "MPC"): the
!>          residual and every norm stay those of the full dof vector (v2), the velocity change
!>          solves the reduced system (pcg_mod with the constraints, the diagonal by
!>          constraint_mod::constraint_diagonal), so that dV = T dV_hat keeps the constraints that
!>          the initial velocities satisfy (dev/theory/addendum_periodic_mpc.md §3.1).
!>          Periodic conditions (step 8.3b, §3.3): the diagonal is the exact one or the lumped one
!>          (constraint_diagonal_lumped, by cons%diag_exact), and after each correction the velocity
!>          is shifted by the uniform translation that gives the zero-mean fluctuation
!>          (constraint_zero_mean, Q-P3); the norms use the shifted correction.
!>
!> Names (ROADMAP D50): the keywords and fields carry FEPX's names; v2's are nl_max_iters,
!> nl_tol_strict, nl_tol_min, nl_tol_loose, nr_tol_switch_ref, nr_tol_conv and nr_slope_start for
!> global_max_iters, global_tol, global_tol_abs, global_tol_slow, global_nr_switch,
!> global_nr_slow_ratio and global_nr_slope_start.
!>
!> Public types:
!>   phase_model_t     a material model and the local elements of its phase
!>   increment_opts_t  iteration parameters (components at the type)
!>   increment_info_t  iteration report (components at the type)
!> Public parameters:
!>   increment_converged, increment_max_iters, increment_diverging, increment_point_failed,
!>   increment_cg_failed  status codes (increment_point_failed is no longer set: material points
!>   whose solve did not converge are counted in increment_info_t%num_point_failed and the
!>   iteration continues, as v2 with its default auto_time 0; ROADMAP D49)
!> Public routines:
!>   increment_solve      solve one increment (collective)
!>   midnode_linearize    reset the edge nodes of 10-node tets to their edge midpoints (D26;
!>                        collective)
module increment_mod

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use kinds_mod, only: rk, ik
  use element_mod, only: element_t
  use element_matrix_mod, only: element_work_t, element_work_alloc, element_kinematics_qpt, &
    & element_stiffness, element_force_h
  use material_model_mod, only: material_model_t
  use material_state_mod, only: material_state_t, material_state_update_derived
  use mesh_mod, only: mesh_t
  use parallel_mod, only: par_sum, par_max, par_wtime
  use halo_mod, only: dist_norm, halo_forward, halo_reverse_add
  use tet10_mod, only: tet10_gmsh_type
  use error_mod, only: fatal_error
  use log_mod, only: log_info, log_warning, log_int_text, log_real_text, log_count_text
  use ebe_mod, only: ebe_matvec, ebe_diagonal, ebe_assemble
  use pcg_mod, only: pcg_solve, pcg_precond_jacobi
  use constraint_mod, only: constraint_t, constraint_diagonal, constraint_diagonal_lumped, &
    & constraint_zero_mean

  implicit none

  private

  public :: phase_model_t, increment_opts_t, increment_info_t, increment_solve, midnode_linearize
  public :: increment_converged, increment_max_iters, increment_diverging, &
    & increment_point_failed, increment_cg_failed

  ! Edges of the 10-node tet in Gmsh order (D35): edge node 4 + k joins corners
  ! tet10_edges(1:2, k).
  integer(ik), parameter :: tet10_edges(2, 6) = reshape([1_ik, 2_ik, 2_ik, 3_ik, 3_ik, 1_ik, &
    & 4_ik, 1_ik, 4_ik, 3_ik, 4_ik, 2_ik], [2, 6])

  ! Status codes.
  integer, parameter :: increment_converged = 0, increment_max_iters = 1, &
    & increment_diverging = 2, increment_point_failed = 3, increment_cg_failed = 4

  !> A material model and the local elements of its phase.
  !>
  !> Components:
  !>   model  class(material_model_t), allocatable  the model
  !>   elts   integer(ik)(:)                         local elements of the phase
  type :: phase_model_t
    class(material_model_t), allocatable :: model
    integer(ik), allocatable :: elts(:)
  end type phase_model_t

  !> Iteration parameters (v2 defaults).
  !>
  !> Components:
  !>   global_max_iters       integer   global iteration limit
  !>   global_tol             real(rk)  relative velocity-change tolerance [-]
  !>   global_tol_abs         real(rk)  absolute velocity-change tolerance, per dof [mm/s]
  !>   cg_tol                 real(rk)  CG relative tolerance [-]
  !>   cg_max_iters           integer   CG iteration limit
  !>   cg_precond             integer   CG preconditioner (pcg_mod pcg_precond_*; step 9.2)
  !>   global_tol_slow        real(rk)  relative tolerance accepted after a slow NR or SA (v2
  !>                                    cases 3, 4) [-]
  !>   global_nr_switch       real(rk)  relative change below which SA switches to NR [-]
  !>   global_nr_slow_ratio   real(rk)  NR is slow when its rate falls below this fraction of
  !>                                    the reference rate [-]
  !>   global_nr_slope_start  integer   NR iteration from which the rate is measured (v2: 2,
  !>                                    internal)
  type :: increment_opts_t
    integer :: global_max_iters = 50
    real(rk) :: global_tol = 5.0e-4_rk
    real(rk) :: global_tol_abs = 1.0e-10_rk
    real(rk) :: cg_tol = 1.0e-8_rk
    integer :: cg_max_iters = 16000
    integer :: cg_precond = pcg_precond_jacobi
    real(rk) :: global_tol_slow = 5.0e-4_rk
    real(rk) :: global_nr_switch = 1.0e-2_rk
    real(rk) :: global_nr_slow_ratio = 0.2_rk
    integer :: global_nr_slope_start = 2
  end type increment_opts_t

  !> Iteration report.
  !>
  !> Components:
  !>   status       integer   increment_converged or a failure code
  !>   num_iters    integer   global iterations performed
  !>   cg_iters     integer   CG iterations, summed over the global iterations
  !>   delu_norm    real(rk)  last relative velocity change [-]
  !>   time_total   real(rk)  wall time of the increment (rank 0's clock) [s]
  !>   time_elements real(rk) wall time of the kinematics, point solves and element matrices [s]
  !>   time_cg      real(rk)  wall time of the CG solves (incl. their diagonal) [s]
  !>   num_point_failed  integer(ik)  material point solves that did not converge, summed over
  !>                     the ranks, at the last evaluation (the converged velocities); the points
  !>                     keep their last iterate (D49)
  !>   max_point_failed  integer(ik)  the largest such count over the global iterations
  !>   num_nr_iters      integer      NR iterations performed (all attempts)
  !>   num_nr_switches   integer      switches from SA to NR
  !>   num_nr_fallbacks  integer      returns from NR to SA (v2 cases 5 and 6)
  !>   exit_case         integer      v2 case that ended the iteration (1-4 converged, 7
  !>                                  diverging, 8 not finite; 0 iteration limit or CG
  !>                                  failure)
  !>   iter_method  character(2)(:)   (num_iters) method of each global iteration, 'SA' or 'NR'
  !>                                  (the convergence record, v2 post.conv; step 6.5b)
  !>   iter_norms   real(rk)(:, :)    (6, num_iters) norms of each global iteration, v2's
  !>                                  columns: r_norm (residual 2-norm, prescribed dofs
  !>                                  included) [N], rx_norm (largest |residual|) [N], f_norm
  !>                                  (2-norm of {F_h}) [N], delu_norm (2-norm of the velocity
  !>                                  change over that of the previous velocities) [-],
  !>                                  delux_norm (largest |change| over the same norm) [-],
  !>                                  u_norm (2-norm of the previous velocities) [mm/s]
  !>   iter_cg      integer(:)        (num_iters) CG iterations of each global iteration
  type :: increment_info_t
    integer :: status = increment_converged
    integer(ik) :: num_point_failed = 0_ik
    integer(ik) :: max_point_failed = 0_ik
    integer :: num_nr_iters = 0
    integer :: num_nr_switches = 0
    integer :: num_nr_fallbacks = 0
    integer :: exit_case = 0
    integer :: num_iters = 0
    integer :: cg_iters = 0
    real(rk) :: delu_norm = 0.0_rk
    real(rk) :: time_total = 0.0_rk
    real(rk) :: time_elements = 0.0_rk
    real(rk) :: time_cg = 0.0_rk
    character(len=2), allocatable :: iter_method(:)
    real(rk), allocatable :: iter_norms(:, :)
    integer, allocatable :: iter_cg(:)
  end type increment_info_t

contains

  !> Solve one increment: the end-of-increment velocities by successive approximation, then
  !> the state (stress, elastic strain, ...) at those velocities. Collective.
  !>
  !> Arguments:
  !>   mesh        (in)    type(mesh_t)            rank-local mesh
  !>   elt         (in)    type(element_t)         element descriptor
  !>   phases      (in)    type(phase_model_t)(:)  models and their elements
  !>   coo0_s      (in)    real(rk)(:, :)          (3, num_nodes_loc) start-of-increment
  !>                                               coordinates, sample frame [mm]
  !>   bc_mask     (in)    logical(:)              (3 num_nodes_loc) prescribed dofs
  !>   dt          (in)    real(rk)                time increment [s]
  !>   opts        (in)    type(increment_opts_t)  iteration parameters
  !>   state_prev  (in)    type(material_state_t)  state at the start of the increment
  !>   state       (inout) type(material_state_t)  state at the end of the increment (out)
  !>   vel         (inout) real(rk)(:)             (3 num_nodes_loc) velocities: in, the initial
  !>                                               guess holding the prescribed values; out, the
  !>                                               converged velocities (owned and ghost) [mm/s]
  !>   info        (out)   type(increment_info_t)  iteration report
  !>   cons        (in, optional) type(constraint_t)  velocity constraints (used when active;
  !>                                               vel must satisfy them on entry)
  !> Locals:
  !>   use_cons     logical                  constraints present and active
  !>   n_owned      integer(ik)              owned dofs
  !>   k_elt, k_tan real(rk)(:, :, :)        secant and tangent element matrices
  !>   f_elt        real(rk)(:, :)           element vectors {f_h}
  !>   f_h, k_vel   real(rk)(:)              assembled {F_h} and [K]{V}
  !>   resid, dvel  real(rk)(:)              residual and velocity change
  !>   diag         real(rk)(:)              assembled diagonal of the matrix solved
  !>   vel_prev     real(rk)(:)              velocities of the previous iteration
  !>   vel_sa       real(rk)(:)              velocities saved at the switch to NR (v2 vel_sa)
  !>   no_mask      logical(:)               all false: an unmasked product for the residual
  !>                                         norm, prescribed dofs included (v2 r_norm)
  !>   delu_norm, delu_prev, u_norm, dv_norm  real(rk)  velocity-change norms
  !>   r_norm, r_norm_ini  real(rk)          residual norms, this and the first iteration [N]
  !>   rx_norm, f_norm  real(rk)             largest |residual| and 2-norm of {F_h} [N]
  !>   delux_norm   real(rk)                 largest |velocity change| over u_norm [-]
  !>   nr_tol_switch  real(rk)               current SA-to-NR switch tolerance
  !>   nr_conv, nr_conv_prev  real(rk)       NR rate and its reference (v2 nr_conv, nr_conv_prev)
  !>   use_nr       logical                  this iteration is NR
  !>   nr_attempted, nr_isslow  logical      NR was tried; NR was found slow (v2 flags)
  !>   num_nr_now   integer                  NR iterations of the current attempt (v2
  !>                                         num_nr_iters)
  !>   num_sa_growing  integer               consecutive SA iterations with a growing change
  !>   iter, cg_iters, cg_status  integer    iteration, CG iterations and status
  !>   num_failed   integer(ik)              failed point solves (all ranks)
  !>   t0, t1       real(rk)                 wall-clock times [s]
  subroutine increment_solve(mesh, elt, phases, coo0_s, bc_mask, dt, opts, state_prev, state, &
    & vel, info, cons)

    type(mesh_t), intent(in) :: mesh
    type(element_t), intent(in) :: elt
    type(phase_model_t), intent(in) :: phases(:)
    real(rk), intent(in) :: coo0_s(:, :)
    logical, intent(in) :: bc_mask(:)
    real(rk), intent(in) :: dt
    type(increment_opts_t), intent(in) :: opts
    type(material_state_t), intent(in) :: state_prev
    type(material_state_t), intent(inout) :: state
    real(rk), intent(inout) :: vel(:)
    type(increment_info_t), intent(out) :: info
    type(constraint_t), intent(in), optional :: cons

    logical :: use_cons
    integer(ik) :: n_owned
    real(rk), allocatable :: k_elt(:, :, :), k_tan(:, :, :), f_elt(:, :)
    real(rk), allocatable :: f_h(:), k_vel(:), resid(:), dvel(:), diag(:), vel_prev(:), vel_sa(:)
    logical, allocatable :: no_mask(:)
    real(rk) :: delu_norm, delu_prev, u_norm, dv_norm, r_norm, r_norm_ini, nr_tol_switch
    real(rk) :: rx_norm, f_norm, delux_norm
    real(rk) :: nr_conv, nr_conv_prev
    logical :: use_nr, nr_attempted, nr_isslow
    integer :: num_nr_now, num_sa_growing, iter, cg_iters, cg_status
    integer(ik) :: num_failed
    real(rk) :: t0, t1

    info%time_total = par_wtime()
    n_owned = 3_ik*mesh%num_nodes_owned
    use_cons = .false.
    if (present(cons)) use_cons = cons%active
    allocate (k_elt(elt%num_dofs, elt%num_dofs, mesh%num_elts), f_elt(elt%num_dofs, mesh%num_elts))
    allocate (k_tan(elt%num_dofs, elt%num_dofs, mesh%num_elts))
    allocate (f_h(size(vel)), k_vel(size(vel)), resid(size(vel)), dvel(size(vel)), &
      & diag(size(vel)), vel_prev(size(vel)), vel_sa(size(vel)), no_mask(size(vel)))
    no_mask = .false.
    allocate (info%iter_method(opts%global_max_iters), &
      & info%iter_norms(6, opts%global_max_iters), info%iter_cg(opts%global_max_iters))

    ! v2's initialization of the iteration (solveit_evp).
    delu_prev = 0.0_rk
    r_norm_ini = 0.0_rk
    vel_sa = vel
    use_nr = .false.
    nr_attempted = .false.
    nr_isslow = .false.
    num_nr_now = 0
    nr_tol_switch = opts%global_nr_switch
    nr_conv = 1.0_rk
    nr_conv_prev = 1.0_rk
    num_sa_growing = 0
    info%status = increment_max_iters

    do iter = 1, opts%global_max_iters
      ! Element matrices and vectors at the current velocities (tangent matrices for NR).
      t0 = par_wtime()
      if (use_nr) then
        call element_system(mesh, elt, phases, coo0_s, dt, state_prev, state, vel, k_elt, &
          & f_elt, num_failed, k_tan)
      else
        call element_system(mesh, elt, phases, coo0_s, dt, state_prev, state, vel, k_elt, &
          & f_elt, num_failed)
      end if
      t1 = par_wtime()
      info%time_elements = info%time_elements + (t1 - t0)
      ! Points that did not converge keep their last iterate and are counted (D49; v2 warns
      ! and continues).
      info%max_point_failed = max(info%max_point_failed, num_failed)

      ! Secant residual {F_h} - [K]{V}: its norm over all dofs (prescribed included, v2
      ! r_norm), then zero at the prescribed dofs for the solve.
      call ebe_assemble(mesh, f_elt, f_h)
      call ebe_matvec(mesh, k_elt, no_mask, vel, k_vel)
      resid = f_h - k_vel
      r_norm = dist_norm(resid, n_owned)
      if (iter == 1) r_norm_ini = r_norm
      ! Largest |residual| and the norm of {F_h} (v2 rx_norm, f_norm; the convergence record).
      call par_max(max(0.0_rk, maxval(abs(resid(1:n_owned)))), rx_norm)
      f_norm = dist_norm(f_h, n_owned)
      where (bc_mask) resid = 0.0_rk
      ! Case 8: a residual that is not finite is a failure (before the solve, whose CG would
      ! return at once with a NaN error and a zero change that passes case 1).
      if (.not. ieee_is_finite(r_norm)) then
        call log_warning('The residual of the global iteration is not finite.', 2)
        info%exit_case = 8
        info%status = increment_diverging
        exit
      end if

      ! Velocity change: [K]{dV} = residual (SA) or [K_tan]{dV} = residual (NR).
      t0 = par_wtime()
      dvel = 0.0_rk
      if (use_nr) then
        num_nr_now = num_nr_now + 1
        info%num_nr_iters = info%num_nr_iters + 1
        if (use_cons) then
          if (cons%diag_exact) then
            call constraint_diagonal(cons, mesh, k_tan, diag)
          else
            call constraint_diagonal_lumped(cons, mesh, k_tan, diag)
          end if
          call pcg_solve(mesh, k_tan, bc_mask, diag, resid, dvel, opts%cg_tol, &
            & opts%cg_max_iters, cg_iters, cg_status, cons)
        else
          call ebe_diagonal(mesh, k_tan, diag)
          call pcg_solve(mesh, k_tan, bc_mask, diag, resid, dvel, opts%cg_tol, &
            & opts%cg_max_iters, cg_iters, cg_status, precond=opts%cg_precond)
        end if
      else
        if (use_cons) then
          if (cons%diag_exact) then
            call constraint_diagonal(cons, mesh, k_elt, diag)
          else
            call constraint_diagonal_lumped(cons, mesh, k_elt, diag)
          end if
          call pcg_solve(mesh, k_elt, bc_mask, diag, resid, dvel, opts%cg_tol, &
            & opts%cg_max_iters, cg_iters, cg_status, cons)
        else
          call ebe_diagonal(mesh, k_elt, diag)
          call pcg_solve(mesh, k_elt, bc_mask, diag, resid, dvel, opts%cg_tol, &
            & opts%cg_max_iters, cg_iters, cg_status, precond=opts%cg_precond)
        end if
      end if
      t1 = par_wtime()
      info%time_cg = info%time_cg + (t1 - t0)
      info%cg_iters = info%cg_iters + cg_iters
      if (cg_status /= 0) then
        ! (Status 2: a nodal block of the block-Jacobi preconditioner not positive definite.)
        if (cg_status == 2) call log_warning('A nodal block of the block-Jacobi ' // &
          & 'preconditioner is not positive definite.', 2)
        info%status = increment_cg_failed
        exit
      end if
      vel_prev = vel
      vel = vel + dvel
      ! Periodic conditions: the rigid translation fixed by the zero-mean fluctuation (Q-P3);
      ! the correction is the projected one.
      if (use_cons) then
        if (cons%zero_mean) then
          call constraint_zero_mean(cons, mesh, vel)
          dvel = vel - vel_prev
        end if
      end if

      ! Norms of the change, relative to the previous velocities (v2: over all dofs).
      dv_norm = dist_norm(dvel, n_owned)
      u_norm = dist_norm(vel_prev, n_owned)
      delu_norm = dv_norm/max(u_norm, tiny(1.0_rk))
      info%num_iters = iter
      info%delu_norm = delu_norm

      ! Convergence record of the iteration (v2 write_conv_file_data, written after the norms).
      call par_max(max(0.0_rk, maxval(abs(dvel(1:n_owned)))), delux_norm)
      delux_norm = delux_norm/max(u_norm, tiny(1.0_rk))
      info%iter_method(iter) = merge('NR', 'SA', use_nr)
      info%iter_norms(:, iter) = [r_norm, rx_norm, f_norm, delu_norm, delux_norm, u_norm]
      info%iter_cg(iter) = cg_iters
      ! The iteration on the log (step 9.7, D94, Q60 (B)): method, R = delu_norm (the relative
      ! velocity change, v2's R) and the CG iterations of its solve.
      call log_info('Iteration ' // log_int_text(iter) // ' (' // info%iter_method(iter) // &
        & '): R = ' // log_real_text(delu_norm, 4) // ' (' // &
        & log_count_text(cg_iters, 'CG iteration') // ')', 2)

      ! NR rate (v2 nr_conv): its reference at NR iteration global_nr_slope_start, then the log of
      ! the ratio of consecutive changes.
      if (use_nr) then
        if (num_nr_now == opts%global_nr_slope_start) then
          nr_conv_prev = nr_conv
        else if (num_nr_now > opts%global_nr_slope_start) then
          nr_conv = log(delu_prev) - log(delu_norm)
        end if
      else
        ! SA: consecutive growing changes.
        if (delu_norm > delu_prev) then
          num_sa_growing = num_sa_growing + 1
        else
          num_sa_growing = 0
        end if
      end if

      ! Case 8: a velocity change that is not finite is a failure.
      if (.not. ieee_is_finite(dv_norm) .or. .not. ieee_is_finite(u_norm)) then
        call log_warning('The velocity change of the global iteration is not finite.', 2)
        info%exit_case = 8
        info%status = increment_diverging
        exit
      end if

      ! v2's cases, in v2's order.
      if (iter > 1 .and. delu_norm < opts%global_tol) then
        ! Case 1: converged.
        info%exit_case = 1
        info%status = increment_converged
        exit
      else if (iter > 1 .and. dv_norm < opts%global_tol_abs*3.0_rk*real(mesh%num_nodes_glob, rk)) &
        & then
        ! Case 2: converged (change below the absolute threshold).
        info%exit_case = 2
        info%status = increment_converged
        exit
      else if (num_nr_now > opts%global_nr_slope_start .and. &
        & nr_conv < opts%global_nr_slow_ratio*nr_conv_prev .and. &
        & delu_norm < opts%global_tol_slow) then
        ! Case 3: NR slow to converge, but converged satisfactorily.
        info%exit_case = 3
        info%status = increment_converged
        call log_info('Newton-Raphson is slow to converge, but converged satisfactorily', 2)
        exit
      else if (.not. use_nr .and. nr_attempted .and. delu_norm < opts%global_tol_slow) then
        ! Case 4: SA slow to converge (after an NR attempt), but converged satisfactorily.
        info%exit_case = 4
        info%status = increment_converged
        call log_info('Successive approximation is slow to converge, but converged ' // &
          & 'satisfactorily', 2)
        exit
      else if (use_nr .and. r_norm > 1.1_rk*r_norm_ini) then
        ! Case 5: NR diverging: back to the saved SA velocities and to SA, with a stricter
        ! switch tolerance.
        call log_warning('Newton-Raphson diverging: back to SA.', 2)
        vel = vel_sa
        use_nr = .false.
        nr_tol_switch = nr_tol_switch/10.0_rk
        nr_attempted = .true.
        num_nr_now = 0
        info%num_nr_fallbacks = info%num_nr_fallbacks + 1
      else if (num_nr_now > opts%global_nr_slope_start .and. &
        & nr_conv < opts%global_nr_slow_ratio*nr_conv_prev) then
        ! Case 6: NR converging slowly: SA for the rest of the increment.
        call log_warning('Newton-Raphson converging slowly: back to SA.', 2)
        use_nr = .false.
        nr_attempted = .true.
        nr_isslow = .true.
        num_nr_now = 0
        info%num_nr_fallbacks = info%num_nr_fallbacks + 1
      else if (num_sa_growing >= 5) then
        ! Case 7: SA diverging: failure.
        info%exit_case = 7
        info%status = increment_diverging
        exit
      end if

      ! Switch from SA to NR, saving the SA velocities (v2).
      if (delu_norm < nr_tol_switch .and. .not. nr_isslow .and. .not. use_nr) then
        use_nr = .true.
        num_sa_growing = 0
        vel_sa = vel
        info%num_nr_switches = info%num_nr_switches + 1
      end if
      delu_prev = delu_norm
    end do
    if (info%status /= increment_converged) then
      info%time_total = par_wtime() - info%time_total
      return
    end if

    ! State at the converged velocities (v2 update_state_evp).
    t0 = par_wtime()
    call element_system(mesh, elt, phases, coo0_s, dt, state_prev, state, vel, k_elt, f_elt, &
      & num_failed)
    info%num_point_failed = num_failed
    info%max_point_failed = max(info%max_point_failed, num_failed)
    info%time_elements = info%time_elements + (par_wtime() - t0)
    info%time_total = par_wtime() - info%time_total

  end subroutine increment_solve

  !> Element matrices and vectors at given velocities: kinematics at every point in the
  !> end-of-increment configuration x0 + dt v, point solves by each phase's model (which set
  !> the current-iterate state), the model-independent derived quantities (kinematics,
  !> strains, work, lattice rotation rate; material_state_update_derived), then (R14).
  !>
  !> Arguments:
  !>   mesh, elt, phases, coo0_s, dt, state_prev  (in)  as increment_solve
  !>   state       (inout) type(material_state_t)  current-iterate state (out)
  !>   vel         (in)    real(rk)(:)             (3 num_nodes_loc) velocities [mm/s]
  !>   k_elt       (out)   real(rk)(:, :, :)       element matrices [N s/mm]
  !>   f_elt       (out)   real(rk)(:, :)          element vectors {f_h} [N]
  !>   num_failed  (out)   integer(ik)             failed point solves, summed over the ranks
  !>   k_tan       (out, optional) real(rk)(:, :, :) tangent element matrices (R14 with
  !>                                               [s_tan], §4; Newton-Raphson) [N s/mm]
  !> Locals:
  !>   coo, vel_n  real(rk)(:, :)       (3, num_nodes_loc) end-of-increment coordinates [mm] and
  !>                                    velocities [mm/s]
  !>   defrate_s   real(rk)(:, :, :)    (6, num_qpts, num_elts) {d} [1/s]
  !>   spin_s      real(rk)(:, :, :)    (3, num_qpts, num_elts) w [1/s]
  !>   mat_s       real(rk)(:, :, :, :) (6, 6, num_qpts, num_elts) [s] [MPa s]
  !>   vec_h       real(rk)(:, :, :)    (6, num_qpts, num_elts) {h} [1/s]
  !>   mat_tan     real(rk)(:, :, :, :) (6, 6, num_qpts, num_elts) [s_tan] (with k_tan) [MPa s]
  !>   i_elt, i_qpt  integer(ik)        element and point
  !>   i_phase     integer              phase
  !>   failed, failed_all  integer(ik)  failed point solves (this rank's phase, all phases)
  !>   ecoo, evel  real(rk)(:, :)       (3, num_nodes) one element's coordinates [mm] and
  !>                                    velocities [mm/s], gathered once per element
  !>   a           integer(ik)          element node
  !>   ework       type(element_work_t) the element routines' work arrays (once per call, D100)
  subroutine element_system(mesh, elt, phases, coo0_s, dt, state_prev, state, vel, k_elt, f_elt, &
    & num_failed, k_tan)

    type(mesh_t), intent(in) :: mesh
    type(element_t), intent(in) :: elt
    type(phase_model_t), intent(in) :: phases(:)
    real(rk), intent(in) :: coo0_s(:, :)
    real(rk), intent(in) :: dt
    type(material_state_t), intent(in) :: state_prev
    type(material_state_t), intent(inout) :: state
    real(rk), intent(in) :: vel(:)
    real(rk), intent(out) :: k_elt(:, :, :)
    real(rk), intent(out) :: f_elt(:, :)
    integer(ik), intent(out) :: num_failed
    real(rk), intent(out), optional :: k_tan(:, :, :)

    real(rk), allocatable :: coo(:, :), vel_n(:, :)
    real(rk), allocatable :: defrate_s(:, :, :), spin_s(:, :, :), mat_s(:, :, :, :), vec_h(:, :, :)
    real(rk), allocatable :: mat_tan(:, :, :, :)
    real(rk), allocatable :: ecoo(:, :), evel(:, :)
    integer(ik) :: i_elt, i_qpt, failed, failed_all, a
    integer :: i_phase
    type(element_work_t) :: ework

    ! End-of-increment configuration; one element's arrays and the element work arrays.
    vel_n = reshape(vel, [3_ik, mesh%num_nodes_loc])
    coo = coo0_s + dt*vel_n
    allocate (ecoo(3, elt%num_nodes), evel(3, elt%num_nodes))
    call element_work_alloc(ework, elt)

    ! Kinematics at every point (the element's coordinates and velocities gathered once).
    allocate (defrate_s(6, elt%num_qpts, mesh%num_elts), spin_s(3, elt%num_qpts, mesh%num_elts))
    do i_elt = 1_ik, mesh%num_elts
      do a = 1_ik, elt%num_nodes
        ecoo(:, a) = coo(:, mesh%elt_nodes_loc(a, i_elt))
        evel(:, a) = vel_n(:, mesh%elt_nodes_loc(a, i_elt))
      end do
      do i_qpt = 1_ik, elt%num_qpts
        call element_kinematics_qpt(elt, ecoo, evel, i_qpt, defrate_s(:, i_qpt, i_elt), &
          & spin_s(:, i_qpt, i_elt), ework)
      end do
    end do

    ! Point solves, one call per phase.
    allocate (mat_s(6, 6, elt%num_qpts, mesh%num_elts), vec_h(6, elt%num_qpts, mesh%num_elts))
    failed_all = 0_ik
    if (present(k_tan)) allocate (mat_tan(6, 6, elt%num_qpts, mesh%num_elts))
    do i_phase = 1, size(phases)
      if (present(k_tan)) then
        call phases(i_phase)%model%solve_points(phases(i_phase)%elts, dt, defrate_s, spin_s, &
          & state_prev, state, mat_s, vec_h, failed, mat_tan)
      else
        call phases(i_phase)%model%solve_points(phases(i_phase)%elts, dt, defrate_s, spin_s, &
          & state_prev, state, mat_s, vec_h, failed)
      end if
      failed_all = failed_all + failed
    end do
    call par_sum(failed_all, num_failed)

    ! Model-independent derived quantities of every point (work, lattice rotation rate).
    call material_state_update_derived(state_prev, defrate_s, spin_s, dt, state)

    ! Element matrices and vectors (R14), the element's coordinates gathered once.
    do i_elt = 1_ik, mesh%num_elts
      do a = 1_ik, elt%num_nodes
        ecoo(:, a) = coo(:, mesh%elt_nodes_loc(a, i_elt))
      end do
      call element_stiffness(elt, ecoo, mat_s(:, :, :, i_elt), k_elt(:, :, i_elt), ework)
      call element_force_h(elt, ecoo, mat_s(:, :, :, i_elt), vec_h(:, :, i_elt), f_elt(:, i_elt), &
        & ework)
      if (present(k_tan)) call element_stiffness(elt, ecoo, mat_tan(:, :, :, i_elt), &
        & k_tan(:, :, i_elt), ework)
    end do

  end subroutine element_system

  !> Mid-node linearization (ROADMAP D26): every element places its edge nodes at the
  !> midpoints of their corner nodes; the values of all elements sharing a node (on every rank)
  !> are averaged, which leaves corner nodes unchanged. 10-node tets only. Collective.
  !>
  !> Source : port of v2 solveit_evp_mod1.f90::solveit_evp (start of each iteration: element
  !>          coordinates with mid-nodes reset, scattered with part_scatter and divided by the
  !>          node multiplicity g_ones).
  !>
  !> Arguments:
  !>   mesh   (in)    type(mesh_t)     rank-local mesh
  !>   elt    (in)    type(element_t)  element descriptor (must be tet10)
  !>   coo_s  (inout) real(rk)(:, :)   (3, num_nodes_loc) coordinates, sample frame [mm];
  !>                                   valid on the owned and ghost nodes on exit
  !> Locals:
  !>   sums    real(rk)(:)             (3 num_nodes_loc) sums of the element values
  !>   counts  real(rk)(:)             (num_nodes_loc) elements contributing to each node
  !>   ecoo    real(rk)(3, 10)         one element's coordinates, mid-nodes reset [mm]
  !>   i_elt, k, a  integer(ik)        element, edge, node
  subroutine midnode_linearize(mesh, elt, coo_s)

    type(mesh_t), intent(in) :: mesh
    type(element_t), intent(in) :: elt
    real(rk), intent(inout) :: coo_s(:, :)

    real(rk), allocatable :: sums(:), counts(:)
    real(rk) :: ecoo(3, 10)
    integer(ik) :: i_elt, k, a

    if (elt%gmsh_type /= tet10_gmsh_type) call fatal_error('midnode_linearize', &
      & 'midnode_linearization applies to 10-node tetrahedra only')

    ! Element values (edge nodes at the midpoints), summed per node with their counts.
    allocate (sums(3*mesh%num_nodes_loc), counts(mesh%num_nodes_loc))
    sums = 0.0_rk
    counts = 0.0_rk
    do i_elt = 1_ik, mesh%num_elts
      ecoo = coo_s(:, mesh%elt_nodes_loc(:, i_elt))
      do k = 1_ik, 6_ik
        ecoo(:, 4_ik + k) = 0.5_rk*(ecoo(:, tet10_edges(1, k)) + ecoo(:, tet10_edges(2, k)))
      end do
      do a = 1_ik, 10_ik
        associate (node => mesh%elt_nodes_loc(a, i_elt))
          sums(3_ik*node - 2_ik:3_ik*node) = sums(3_ik*node - 2_ik:3_ik*node) + ecoo(:, a)
          counts(node) = counts(node) + 1.0_rk
        end associate
      end do
    end do

    ! Totals on the owners, averages, then the owners' values on the ghosts.
    call halo_reverse_add(mesh, sums, 3_ik)
    call halo_reverse_add(mesh, counts, 1_ik)
    do a = 1_ik, mesh%num_nodes_owned
      if (counts(a) > 0.0_rk) coo_s(:, a) = sums(3_ik*a - 2_ik:3_ik*a)/counts(a)
    end do
    sums = reshape(coo_s, [3_ik*mesh%num_nodes_loc])
    call halo_forward(mesh, sums, 3_ik)
    coo_s = reshape(sums, [3_ik, mesh%num_nodes_loc])

  end subroutine midnode_linearize

end module increment_mod
