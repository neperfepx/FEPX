! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Preconditioned conjugate gradient on the element-by-element operator (ebe_mod), with
!> prescribed (Dirichlet) dofs held fixed; Jacobi preconditioner (v2's, the default), nodal
!> block-Jacobi, or none (step 9.2).
!>
!> Source : port of v2 conjugate_gradient_mod1.f90::cg_solver_ebe (case "general"): the same
!>          iteration, preconditioner (the inverse of the assembled diagonal), stopping test
!>          (preconditioned residual norm relative to its initial value below cg_tol) and
!>          iteration limit; the dot products are over the owned dofs and summed over the ranks
!>          (halo_mod::dist_dot, v2 par_sum). Defaults of the callers: cg_tol 1e-8,
!>          cg_max_iters 16000 (v2 exec_type_mod.f90). Deviations: the iteration-limit failure is
!>          returned as a status (v2 stops the run: the caller does, with its own message); a
!>          zero right-hand side (preconditioned residual exactly zero) returns the initial
!>          guess without iterating (v2 would divide 0 by 0 in its first step).
!>          With constraints v = T v_hat + g (step 8.3a; v2 cg_solver_ebe "MPC"/"PBC", R-P4
!>          Algorithm 1; dev/theory/addendum_periodic_mpc.md §3.1) the same iteration runs on the
!>          reduced system T^T K T dv_hat = T^T rhs: the right-hand side reduced by T^T, the
!>          products by constraint_mod::constraint_matvec, the residual zero at the prescribed
!>          dofs, the preconditioner the inverse of the reduced diagonal (zero at the constrained
!>          dofs, as v2's "where (gdiag /= 0)"), and the solution expanded, dv = T dv_hat (no
!>          offsets: a correction).
!>          Periodic conditions (step 8.3b, addendum §3.3, Q-P15): the rigid translations of the
!>          components without a velocity condition are a null space of T^T K T; the reduced
!>          right-hand side is projected onto the range (its component along each null vector
!>          n_k, the indicator of the independent dofs of component k, subtracted), so that the
!>          singular system is consistent to round-off and the PCG can meet its relative
!>          tolerance (Kaasschieter 1988, J. Comput. Appl. Math. 24:265-275; Bochev & Lehoucq
!>          2005, SIAM Review 47:50-66, §3). v2 has no projection (v2-fixed C74 adds it).
!>          Preconditioners (step 9.2, dev/theory/addendum_block_jacobi.md §2; Barrett et al.
!>          1994, Templates, §2.5 Fig. 2.5, §3.2, §3.2.1; Saad 2003 §9.2, §12.2): M = diag(K)
!>          (pcg_precond_jacobi, v2's), M = the nodal 3 x 3 diagonal blocks of K restricted to
!>          the free dofs (pcg_precond_block_jacobi; each block inverted once per solve by
!>          linalg_mod::ldlt_factor and ldlt_substitute, a port of v2
!>          solve_evp_vpstress_mod3.f90::solvit), or M = I (pcg_precond_none, plain CG, for
!>          debugging). Only the step z = M^-1 r depends on M; r.z is computed from z
!>          (bit-for-bit v2's r.(D^-1 r) for Jacobi). With constraints only Jacobi (Q-B2;
!>          config_mod rejects the others).
!>
!> Public routines:
!>   pcg_solve   solve [K] x = rhs on the free dofs
!> Public constants:
!>   pcg_precond_none, pcg_precond_jacobi, pcg_precond_block_jacobi   preconditioner codes
module pcg_mod

  use kinds_mod, only: rk, ik
  use mesh_mod, only: mesh_t
  use halo_mod, only: halo_forward, dist_dot
  use ebe_mod, only: ebe_matvec
  use parallel_mod, only: par_any
  use ebe_mod, only: ebe_block_diagonal
  use linalg_mod, only: ldlt_factor, ldlt_substitute
  use error_mod, only: internal_error
  use constraint_mod, only: constraint_t, constraint_matvec, constraint_reduce, &
    & constraint_expand

  implicit none

  private

  public :: pcg_solve
  public :: pcg_precond_none, pcg_precond_jacobi, pcg_precond_block_jacobi

  !> Preconditioner codes: none (M = I), Jacobi (M = diag(K)), nodal block-Jacobi.
  integer, parameter :: pcg_precond_none = 0, pcg_precond_jacobi = 1, &
    & pcg_precond_block_jacobi = 2

contains

  !> Solve [K] x = rhs for the free dofs (x unchanged at the prescribed dofs), by PCG.
  !>
  !> Arguments:
  !>   mesh          (in)    type(mesh_t)       rank-local mesh
  !>   k_elt         (in)    real(rk)(:, :, :)  element matrices (symmetric positive definite on
  !>                                            the free dofs once assembled)
  !>   bc_mask       (in)    logical(:)         (3 num_nodes_loc) prescribed dofs
  !>   diag          (in)    real(rk)(:)        (3 num_nodes_loc) assembled diagonal of [K]
  !>                                            (ebe_diagonal); read by Jacobi only
  !>   rhs           (in)    real(rk)(:)        (3 num_nodes_loc) right-hand side, owned dofs
  !>                                            (ignored at the prescribed dofs)
  !>   x             (inout) real(rk)(:)        (3 num_nodes_loc) initial guess (usually zero) and
  !>                                            solution, valid on the owned and ghost dofs on exit
  !>   cg_tol        (in)    real(rk)           relative tolerance (v2 default 1e-8) [-]
  !>   cg_max_iters  (in)    integer            iteration limit (v2 default 16000)
  !>   num_iters     (out)   integer            iterations performed
  !>   status        (out)   integer            0: converged; 1: iteration limit reached;
  !>                                            2: a nodal block not positive definite
  !>                                            (block-Jacobi; on every rank)
  !>   cons          (in, optional) type(constraint_t)  constraints (used when active): diag
  !>                                            is then the reduced diagonal
  !>                                            (constraint_diagonal), rhs a full residual
  !>                                            (reduced here), x on exit T x_hat
  !>   precond       (in, optional) integer     pcg_precond_* (default pcg_precond_jacobi);
  !>                                            with active constraints only Jacobi
  !> Locals:
  !>   use_cons  logical          constraints present and active
  !>   rhs_red   real(rk)(:)      right-hand side, reduced by T^T when constrained (periodic:
  !>                              projected onto the range)
  !>   null_k    real(rk)(:)      periodic: the null vector of component k (1 at its owned
  !>                              independent dofs)
  !>   k         integer          component
  !>   n_owned   integer(ik)      owned dofs (3 num_nodes_owned)
  !>   precond_kind  integer         preconditioner code (pcg_precond_*)
  !>   gdiag     real(rk)(:)      Jacobi: inverse of the diagonal; none: 1
  !>   binv      real(rk)(:, :, :) block-Jacobi: (3, 3, num_nodes_owned) inverse of each owned
  !>                              node's block restricted to its free components
  !>   ru        real(rk)(:)      residual rhs - [K] x
  !>   zu        real(rk)(:)      preconditioned residual
  !>   pu, apu   real(rk)(:)      search direction and [K] times it
  !>   mag       real(rk)         norm of the initial preconditioned residual
  !>   error     real(rk)         relative preconditioned residual norm
  !>   xnumer, xnumer_o  real(rk) r . (M^-1 r) of this and the previous iteration
  !>   alpha, beta       real(rk) CG step and direction coefficients
  subroutine pcg_solve(mesh, k_elt, bc_mask, diag, rhs, x, cg_tol, cg_max_iters, num_iters, &
    & status, cons, precond)

    type(mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: k_elt(:, :, :)
    logical, intent(in) :: bc_mask(:)
    real(rk), intent(in) :: diag(:)
    real(rk), intent(in) :: rhs(:)
    real(rk), intent(inout) :: x(:)
    real(rk), intent(in) :: cg_tol
    integer, intent(in) :: cg_max_iters
    integer, intent(out) :: num_iters
    integer, intent(out) :: status
    type(constraint_t), intent(in), optional :: cons
    integer, intent(in), optional :: precond

    logical :: use_cons
    integer :: precond_kind
    real(rk), allocatable :: binv(:, :, :)
    real(rk), allocatable :: rhs_red(:), null_k(:)
    integer :: k
    integer(ik) :: n_owned
    real(rk), allocatable :: gdiag(:), ru(:), zu(:), pu(:), apu(:)
    real(rk) :: mag, error, xnumer, xnumer_o, alpha, beta

    n_owned = 3_ik*mesh%num_nodes_owned
    allocate (gdiag(size(x)), ru(size(x)), zu(size(x)), pu(size(x)), apu(size(x)))
    num_iters = 0
    status = 0

    use_cons = .false.
    if (present(cons)) use_cons = cons%active
    precond_kind = pcg_precond_jacobi
    if (present(precond)) precond_kind = precond
    if (use_cons .and. precond_kind /= pcg_precond_jacobi) call internal_error('pcg_solve', &
      & 'only the Jacobi preconditioner is available with constraints')

    ! Preconditioner. Jacobi: the inverse assembled diagonal (constrained: zero at the
    ! constrained dofs, whose reduced diagonal is zero). None: 1 (z = r exactly). Block-Jacobi:
    ! the inverses of the nodal blocks (a block not positive definite on any rank: status 2).
    select case (precond_kind)
    case (pcg_precond_jacobi)
      if (use_cons) then
        gdiag = 0.0_rk
        where (.not. cons%is_constrained) gdiag = 1.0_rk/diag
      else
        gdiag = 1.0_rk/diag
      end if
    case (pcg_precond_none)
      gdiag = 1.0_rk
    case (pcg_precond_block_jacobi)
      call block_inverses(binv, status)
      if (status /= 0) return
    case default
      call internal_error('pcg_solve', 'unknown preconditioner code')
    end select

    ! Right-hand side (constrained: T^T rhs).
    rhs_red = rhs
    if (use_cons) call constraint_reduce(cons, n_owned, rhs_red)

    ! Periodic conditions: the reduced right-hand side projected onto the range (Q-P15).
    if (use_cons) then
      if (cons%zero_mean) then
        allocate (null_k(size(x)))
        do k = 1, 3
          if (cons%mean_fixed(k)) cycle
          null_k = 0.0_rk
          null_k(k::3) = 1.0_rk
          where (cons%is_constrained .or. bc_mask) null_k = 0.0_rk
          rhs_red = rhs_red - dist_dot(null_k, rhs_red, n_owned)/ &
            & dist_dot(null_k, null_k, n_owned)*null_k
        end do
      end if
    end if

    ! Initial residual ru = rhs - [K] x (zero at the prescribed dofs) and its preconditioned norm.
    call halo_forward(mesh, x, 3_ik)
    call apply_operator(x, apu)
    ru = rhs_red - apu
    where (bc_mask) ru = 0.0_rk
    call apply_preconditioner(ru, zu)
    mag = sqrt(dist_dot(zu, zu, n_owned))

    ! Zero right-hand side: the initial guess is the solution (constrained: expanded).
    if (.not. (mag > 0.0_rk)) then
      if (use_cons) call constraint_expand(cons, x, .false.)
      return
    end if

    ! CG iterations (v2 order of operations).
    xnumer_o = 1.0_rk
    pu = 0.0_rk
    error = 1.0_rk
    do while (error > cg_tol)
      num_iters = num_iters + 1
      if (num_iters > cg_max_iters) then
        status = 1
        num_iters = cg_max_iters
        exit
      end if

      ! Direction: pu = zu + beta pu, with beta = (r.M^-1 r) / (previous r.M^-1 r).
      xnumer = dist_dot(ru, zu, n_owned)
      beta = xnumer/xnumer_o
      if (num_iters == 1) beta = 0.0_rk
      xnumer_o = xnumer
      pu = zu + beta*pu

      ! Step: alpha = (r.M^-1 r) / (p.Kp).
      call halo_forward(mesh, pu, 3_ik)
      call apply_operator(pu, apu)
      alpha = xnumer/dist_dot(pu, apu, n_owned)

      ! Update of the solution and of the residual.
      x = x + alpha*pu
      ru = ru - alpha*apu
      call apply_preconditioner(ru, zu)
      error = sqrt(dist_dot(zu, zu, n_owned))/mag
    end do

    ! Solution valid on the ghost dofs too (constrained: expanded, dv = T dv_hat).
    call halo_forward(mesh, x, 3_ik)
    if (use_cons) call constraint_expand(cons, x, .false.)

  contains

    !> The operator: [K] v, or T^T [K] T v with constraints; zero at the prescribed dofs.
    !>
    !> Arguments:
    !>   v   (in)   real(rk)(:)  input, valid on the owned and ghost dofs
    !>   kv  (out)  real(rk)(:)  result, complete on the owned dofs
    subroutine apply_operator(v, kv)

      real(rk), intent(in) :: v(:)
      real(rk), intent(out) :: kv(:)

      if (use_cons) then
        call constraint_matvec(cons, mesh, k_elt, bc_mask, v, kv)
      else
        call ebe_matvec(mesh, k_elt, bc_mask, v, kv)
      end if

    end subroutine apply_operator

    !> The preconditioner step z = M^-1 r: Jacobi and none z = r gdiag (on all local dofs, as
    !> v2); block-Jacobi z_a = binv_a r_a on the owned nodes, zero on the ghosts (their values
    !> are not used: the direction is completed by halo_forward before a product).
    !>
    !> Arguments:
    !>   r  (in)   real(rk)(:)  residual, complete on the owned dofs, zero at the prescribed
    !>   z  (out)  real(rk)(:)  preconditioned residual
    !> Locals:
    !>   a     integer(ik)  owned node
    !>   base  integer(ik)  dof before node a's, 3 (a - 1)
    subroutine apply_preconditioner(r, z)

      real(rk), intent(in) :: r(:)
      real(rk), intent(out) :: z(:)

      integer(ik) :: a, base

      if (precond_kind == pcg_precond_block_jacobi) then
        z = 0.0_rk
        do a = 1_ik, mesh%num_nodes_owned
          base = 3_ik*(a - 1_ik)
          z(base + 1) = binv(1, 1, a)*r(base + 1) + binv(1, 2, a)*r(base + 2) + &
            & binv(1, 3, a)*r(base + 3)
          z(base + 2) = binv(2, 1, a)*r(base + 1) + binv(2, 2, a)*r(base + 2) + &
            & binv(2, 3, a)*r(base + 3)
          z(base + 3) = binv(3, 1, a)*r(base + 1) + binv(3, 2, a)*r(base + 2) + &
            & binv(3, 3, a)*r(base + 3)
        end do
      else
        z = r*gdiag
      end if

    end subroutine apply_preconditioner

    !> The inverses of the owned nodes' 3 x 3 diagonal blocks (ebe_block_diagonal), each
    !> restricted to the node's free components: the row and column of a prescribed component
    !> replaced by those of the identity, so that z is zero there (addendum_block_jacobi.md §2).
    !> Each inverse from one LDL^T factorization (ldlt_factor; the block is symmetric positive
    !> definite: a principal submatrix of [K] on the free dofs) and a substitution per unit
    !> vector (ldlt_substitute); zero when the factorization fails.
    !>
    !> Arguments:
    !>   inv     (out) real(rk)(:, :, :)  (3, 3, num_nodes_owned) inverses
    !>   status  (out) integer            0; 2: a block not positive definite on some rank
    !> Locals:
    !>   blocks  real(rk)(:)    (9 num_nodes_loc) assembled blocks
    !>   block   real(rk)(3,3)  one node's block, restricted, then its LDL^T factors
    !>   a       integer(ik)    owned node
    !>   base    integer(ik)    dof before node a's, 3 (a - 1)
    !>   k       integer        component
    !>   ldlt_status  integer   ldlt_factor status
    !>   bad, any_bad logical   a block not positive definite on this rank, on any rank
    subroutine block_inverses(inv, status)

      real(rk), allocatable, intent(out) :: inv(:, :, :)
      integer, intent(out) :: status

      real(rk), allocatable :: blocks(:)
      real(rk) :: block(3, 3)
      integer(ik) :: a, base
      integer :: k, ldlt_status
      logical :: bad, any_bad

      ! Assembled blocks.
      allocate (blocks(9_ik*mesh%num_nodes_loc), inv(3, 3, mesh%num_nodes_owned))
      call ebe_block_diagonal(mesh, k_elt, blocks)

      ! Per owned node: restrict to the free components, factorize, invert column by column.
      bad = .false.
      do a = 1_ik, mesh%num_nodes_owned
        base = 3_ik*(a - 1_ik)
        block = reshape(blocks(9_ik*(a - 1_ik) + 1:9_ik*a), [3, 3])
        do k = 1, 3
          if (bc_mask(base + k)) then
            block(k, :) = 0.0_rk
            block(:, k) = 0.0_rk
            block(k, k) = 1.0_rk
          end if
        end do
        call ldlt_factor(block, ldlt_status)
        if (ldlt_status /= 0) then
          bad = .true.
          inv(:, :, a) = 0.0_rk
          cycle
        end if
        do k = 1, 3
          inv(:, k, a) = 0.0_rk
          inv(k, k, a) = 1.0_rk
          call ldlt_substitute(block, inv(:, k, a))
        end do
      end do

      ! The same status on every rank (the iteration's sums are collective).
      call par_any(bad, any_bad)
      status = merge(2, 0, any_bad)

    end subroutine block_inverses

  end subroutine pcg_solve

end module pcg_mod
