! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Elastic-only material model (ROADMAP D13, D40): anisotropic linear elasticity with no slip,
!> i.e. the formulation of recoupled_formulation.md with [M*] = 0.
!>
!> Theory : recoupled_formulation.md, with d^p = 0 and [M*] = 0:
!>          (R1)  {d} = {e^e dot} + [W^]{e^e},  W^ built from w^p by (R8);
!>          (R2)  w^p = w + 1/2 (e^e dot e^e - e^e e^e dot)   (d^p = 0);
!>          (R9)  {e^e dot} = ({e^e} - {e^e_0})/dt, e^e_0 held in the sample frame;
!>          (R3)-(R4)  {tau} = [C_s]{e^e}, sigma = tau/beta, beta = det(I + e^e);
!>          (R11) [s] = (dt/beta)[C_s];  (R12) {h} = [W^]{e^e} - {e^e_0}/dt;
!>          tangent (§4, [M*_tan] = 0): [s_tan] = [s].
!>          Substituting (R9) in (R1): ((1/dt)[I] + [W^]){e^e} = {d} + {e^e_0}/dt, a linear 6 x 6
!>          system for a given w^p; then sigma = [s]({d} - {h}) exactly (R10).
!> Source : per-point algorithm: port of v2 solve_evp_mod2.f90::solve_evp_dev (§5): alternate
!>          the stress solve for the current w^p and the update of w^p (v2:
!>          solve_evp_rstar / find_wp_hat) until the stress norm changes by less than
!>          tol_state, at most max_iters_state times (v2 defaults 1e-5 and 100,
!>          types_mod.f90). Deviation (D40): v2 scales tol_state by the initial slip strength
!>          tau_0, which an elastic-only model does not have; here the threshold is relative,
!>          |norm(sigma_k) - norm(sigma_(k-1))| <= tol_state norm(sigma_k). The 6 x 6 system is
!>          solved as v2 solves its per-point systems: normal equations and LDL^T
!>          (linalg_mod). The orientation is updated with the lattice spin inside the loop, as
!>          v2 (step 5.2, replacing D40's start-of-increment stiffness): a forward-Euler estimate
!>          from the continuum spin before the loop, then after each stress solve the
!>          backward-Euler update from the new w^p (v2 solve_evp_rstar, 'euler_forward' /
!>          'euler_backward'; rotation_mod::rot_c2s_spin_update); each stress solve uses the
!>          stiffness of the current orientation. Deviation: the forward estimate expresses the
!>          spin in the start-of-increment crystal frame (v2: in that of its current results%ori,
!>          the previous global iterate), a difference of order (w dt)^2.
!>
!> Public types:
!>   elastic_model_t   the model (extends material_model_t; components at the type)
!>
!> Public routines:
!>   elastic_model_init   set up the model from a crystal-frame Mandel stiffness
!>   elastic_point_solve  the per-point kernel (pure)
module elastic_model_mod

  use kinds_mod, only: rk, ik
  use tensor_mod, only: matrix_from_mandel, mandel_norm, skew_from_matrix, spin_commutator
  use rotation_mod, only: transform_crystal_to_sample_mat6, rot_c2s_spin_update
  use linalg_mod, only: normal_equations_solve
  use material_state_mod, only: material_state_t
  use material_model_mod, only: material_model_t

  implicit none

  private

  public :: elastic_model_t, elastic_model_init, elastic_point_solve

  !> Elastic-only model.
  !>
  !> Components (in addition to material_model_t's num_slip = 0):
  !>   stiff_el_c       real(rk)(6, 6)  elastic stiffness [C_c], Mandel, crystal frame [MPa]
  !>   tol_state      real(rk)        relative stress-norm tolerance of the point loop [-]
  !>   max_iters_state  integer         iteration limit of the point loop
  type, extends(material_model_t) :: elastic_model_t
    real(rk) :: stiff_el_c(6, 6) = 0.0_rk
    real(rk) :: tol_state = 1.0e-5_rk
    integer :: max_iters_state = 100
  contains
    procedure :: init_state => elastic_init_state
    procedure :: solve_points => elastic_solve_points
  end type elastic_model_t

contains

  !> Set up the model.
  !>
  !> Arguments:
  !>   model            (out) type(elastic_model_t)  model
  !>   stiff_el_c       (in)  real(rk)(6, 6)         [C_c], Mandel, crystal frame [MPa]
  !>                                                 (crystal_mod stiffness_* builders)
  !>   tol_state      (in, optional) real(rk)      point-loop tolerance (default 1e-5, v2)
  !>   max_iters_state  (in, optional) integer       point-loop limit (default 100, v2)
  subroutine elastic_model_init(model, stiff_el_c, tol_state, max_iters_state)

    type(elastic_model_t), intent(out) :: model
    real(rk), intent(in) :: stiff_el_c(6, 6)
    real(rk), intent(in), optional :: tol_state
    integer, intent(in), optional :: max_iters_state

    ! No slip systems; the stiffness and the loop controls.
    model%num_slip = 0_ik
    model%stiff_el_c = stiff_el_c
    if (present(tol_state)) model%tol_state = tol_state
    if (present(max_iters_state)) model%max_iters_state = max_iters_state

  end subroutine elastic_model_init

  !> Initial state: nothing beyond the common fields (zero stress and elastic strain, set by
  !> material_state_alloc; orientations from the mesh, set by the caller).
  !>
  !> Arguments:
  !>   self   (in)    class(elastic_model_t)   model
  !>   elts   (in)    integer(ik)(:)           elements of the phase
  !>   state  (inout) type(material_state_t)   state
  subroutine elastic_init_state(self, elts, state)

    class(elastic_model_t), intent(in) :: self
    integer(ik), intent(in) :: elts(:)
    type(material_state_t), intent(inout) :: state

    ! Zero stress and elastic strain at the phase's points (crss has no rows).
    state%stress_s(:, :, elts) = 0.0_rk
    state%strain_el_s(:, :, elts) = 0.0_rk
    state%spin_lat_s(:, :, elts) = 0.0_rk
    if (self%num_slip /= 0_ik) state%crss(:, :, elts) = 0.0_rk

  end subroutine elastic_init_state

  !> Point solve at every quadrature point of the given elements (material_model_t contract).
  !>
  !> Arguments: see material_model_mod::solve_points_iface.
  !> Locals:
  !>   i, i_elt, i_qpt  integer(ik)  element list index, element, quadrature point
  !>   status           integer      point status (0: converged)
  subroutine elastic_solve_points(self, elts, dt, defrate_s, spin_s, state_prev, state, &
    & mat_secant_s, vec_h_s, num_failed, mat_tangent_s)

    class(elastic_model_t), intent(in) :: self
    integer(ik), intent(in) :: elts(:)
    real(rk), intent(in) :: dt
    real(rk), intent(in) :: defrate_s(:, :, :)
    real(rk), intent(in) :: spin_s(:, :, :)
    type(material_state_t), intent(in) :: state_prev
    type(material_state_t), intent(inout) :: state
    real(rk), intent(inout) :: mat_secant_s(:, :, :, :)
    real(rk), intent(inout) :: vec_h_s(:, :, :)
    integer(ik), intent(out) :: num_failed
    real(rk), intent(inout), optional :: mat_tangent_s(:, :, :, :)

    integer(ik) :: i, i_elt, i_qpt
    integer :: status

    num_failed = 0_ik
    do i = 1_ik, int(size(elts), ik)
      i_elt = elts(i)
      do i_qpt = 1_ik, state%num_qpts
        ! Per-point solve (stress, elastic strain, lattice spin and orientation at the end of
        ! the increment).
        call elastic_point_solve(self%stiff_el_c, state_prev%rot_c2s(:, :, i_qpt, i_elt), &
          & state_prev%strain_el_s(:, i_qpt, i_elt), defrate_s(:, i_qpt, i_elt), &
          & spin_s(:, i_qpt, i_elt), dt, self%tol_state, self%max_iters_state, &
          & state%stress_s(:, i_qpt, i_elt), state%strain_el_s(:, i_qpt, i_elt), &
          & state%spin_lat_s(:, i_qpt, i_elt), state%rot_c2s(:, :, i_qpt, i_elt), &
          & mat_secant_s(:, :, i_qpt, i_elt), vec_h_s(:, i_qpt, i_elt), status)
        if (status /= 0) num_failed = num_failed + 1_ik

        ! The spin of the orientation update (the last w^p of the loop); no slip spin.
        state%spin_hat_update_s(:, i_qpt, i_elt) = state%spin_lat_s(:, i_qpt, i_elt)
        state%spin_slip_update_s(:, i_qpt, i_elt) = 0.0_rk

        ! Tangent: equal to the secant ([M*_tan] = 0).
        if (present(mat_tangent_s)) then
          mat_tangent_s(:, :, i_qpt, i_elt) = mat_secant_s(:, :, i_qpt, i_elt)
        end if
      end do
    end do

  end subroutine elastic_solve_points

  !> Elastic-only point solve (module header): alternate the linear stress solve for the
  !> current lattice spin w^p and orientation, and the update of w^p by (R2) and of the
  !> orientation, until the stress norm settles.
  !>
  !> Theory : module header ((R1)-(R4), (R8)-(R12) with [M*] = 0).
  !> Source : port of v2 solve_evp_mod2.f90::solve_evp_dev loop structure (module header).
  !>
  !> Arguments:
  !>   stiff_el_c       (in)  real(rk)(6, 6)  [C_c], Mandel, crystal frame [MPa]
  !>   rot_c2s_prev     (in)  real(rk)(3, 3)  crystal-to-sample rotation, start of increment [-]
  !>   strain_el_s_prev (in)  real(rk)(6)     {e^e_0}, start of increment, sample frame [-]
  !>   defrate_s        (in)  real(rk)(6)     {d}, Mandel, sample frame [1/s]
  !>   spin_s           (in)  real(rk)(3)     w (w23 w13 w12), sample frame [1/s]
  !>   dt               (in)  real(rk)        time increment [s]
  !>   tol_state      (in)  real(rk)        relative stress-norm tolerance [-]
  !>   max_iters_state  (in)  integer         iteration limit
  !>   stress_s         (out) real(rk)(6)     Cauchy stress sigma, Mandel, sample frame [MPa]
  !>   strain_el_s      (out) real(rk)(6)     {e^e}, sample frame [-]
  !>   spin_lat_s       (out) real(rk)(3)     w^p (w23 w13 w12), sample frame [1/s]
  !>   rot_c2s          (out) real(rk)(3, 3)  crystal-to-sample rotation, end of increment [-]
  !>   mat_secant_s     (out) real(rk)(6, 6)  [s] = (dt/beta)[C_s], [C_s] of the final
  !>                                          orientation (v2) [MPa s]
  !>   vec_h_s          (out) real(rk)(6)     {h} = [W^]{e^e} - {e^e_0}/dt [1/s]
  !>   status           (out) integer         0: converged; 1: iteration limit; 2: singular
  !>                                          6 x 6 system
  !> Locals:
  !>   stiff_el_s       real(rk)(6, 6)        [C_s] of the current orientation (first: of the
  !>                                          start orientation) [MPa]
  !>   spin_mat         real(rk)(6, 6)        [W^] of the final w^p [1/s]
  !>   amat             real(rk)(6, 6)        (1/dt)[I] + [W^] [1/s]
  !>   amat_ta          real(rk)(6, 6)        A^T A of the normal equations, then its LDL^T
  !>                                          factors [1/s^2]
  !>   rhs              real(rk)(6)           {d} + {e^e_0}/dt [1/s]
  !>   strain_mat       real(rk)(3, 3)        e^e as a matrix [-]
  !>   strain_rate_mat  real(rk)(3, 3)        e^e dot as a matrix [1/s]
  !>   beta             real(rk)              det(I + e^e) [-]
  !>   norm_prev, norm  real(rk)              stress norms of consecutive iterations [MPa]
  !>   iter, i, solve_status  integer         iteration, index, linear-solve status
  pure subroutine elastic_point_solve(stiff_el_c, rot_c2s_prev, strain_el_s_prev, defrate_s, &
    & spin_s, dt, tol_state, max_iters_state, stress_s, strain_el_s, spin_lat_s, rot_c2s, &
    & mat_secant_s, vec_h_s, status)

    real(rk), intent(in) :: stiff_el_c(6, 6)
    real(rk), intent(in) :: rot_c2s_prev(3, 3)
    real(rk), intent(in) :: strain_el_s_prev(6)
    real(rk), intent(in) :: defrate_s(6)
    real(rk), intent(in) :: spin_s(3)
    real(rk), intent(in) :: dt
    real(rk), intent(in) :: tol_state
    integer, intent(in) :: max_iters_state
    real(rk), intent(out) :: stress_s(6)
    real(rk), intent(out) :: strain_el_s(6)
    real(rk), intent(out) :: spin_lat_s(3)
    real(rk), intent(out) :: rot_c2s(3, 3)
    real(rk), intent(out) :: mat_secant_s(6, 6)
    real(rk), intent(out) :: vec_h_s(6)
    integer, intent(out) :: status

    real(rk) :: stiff_el_s(6, 6), spin_mat(6, 6), amat(6, 6), amat_ta(6, 6), rhs(6)
    real(rk) :: strain_mat(3, 3), strain_rate_mat(3, 3), beta, norm_prev, norm
    integer :: iter, i, solve_status

    ! The right-hand side {d} + {e^e_0}/dt.
    rhs = defrate_s + strain_el_s_prev/dt

    ! Start: lattice spin = continuum spin, forward-Euler orientation estimate from it (v2
    ! 'euler_forward'), stress norm of the start-of-increment strain (start orientation).
    spin_lat_s = spin_s
    rot_c2s = rot_c2s_spin_update(rot_c2s_prev, rot_c2s_prev, spin_lat_s, dt)
    strain_el_s = strain_el_s_prev
    stiff_el_s = transform_crystal_to_sample_mat6(rot_c2s_prev, stiff_el_c)
    norm_prev = mandel_norm(matmul(stiff_el_s, strain_el_s_prev))
    stiff_el_s = transform_crystal_to_sample_mat6(rot_c2s, stiff_el_c)
    stress_s = 0.0_rk
    beta = 1.0_rk
    status = 1

    do iter = 1, max_iters_state
      ! Sample-frame stiffness of the current orientation (R4).
      stiff_el_s = transform_crystal_to_sample_mat6(rot_c2s, stiff_el_c)

      ! Stress solve for the current w^p: ((1/dt)[I] + [W^]){e^e} = rhs.
      amat = spin_commutator(spin_lat_s)
      do i = 1, 6
        amat(i, i) = amat(i, i) + 1.0_rk/dt
      end do
      call normal_equations_solve(amat, rhs, strain_el_s, amat_ta, solve_status)
      if (solve_status /= 0) then
        status = 2
        exit
      end if

      ! Kirchhoff and Cauchy stress (R3)-(R4).
      strain_mat = matrix_from_mandel(strain_el_s)
      beta = det_identity_plus(strain_mat)
      stress_s = matmul(stiff_el_s, strain_el_s)/beta

      ! Lattice spin update (R2) with d^p = 0: w^p = w + 1/2 (e dot e - e e dot).
      strain_rate_mat = matrix_from_mandel((strain_el_s - strain_el_s_prev)/dt)
      spin_lat_s = spin_s + skew_from_matrix(0.5_rk*(matmul(strain_rate_mat, strain_mat) - &
        & matmul(strain_mat, strain_rate_mat)))

      ! Orientation from the new lattice spin (v2 'euler_backward': the spin in the crystal
      ! frame of the current estimate, applied to the start orientation).
      rot_c2s = rot_c2s_spin_update(rot_c2s_prev, rot_c2s, spin_lat_s, dt)

      ! Convergence of the stress norm (v2's test, relative threshold).
      norm = mandel_norm(stress_s)
      if (abs(norm - norm_prev) <= tol_state*norm) then
        status = 0
        exit
      end if
      norm_prev = norm
    end do

    ! Secant matrix (R11) with the stiffness of the final orientation (v2 elt_stif_evp uses the
    ! updated results%ori; the stress is that of the last stress solve), and vector (R12) with
    ! the final w^p.
    stiff_el_s = transform_crystal_to_sample_mat6(rot_c2s, stiff_el_c)
    mat_secant_s = (dt/beta)*stiff_el_s
    spin_mat = spin_commutator(spin_lat_s)
    vec_h_s = matmul(spin_mat, strain_el_s) - strain_el_s_prev/dt

  end subroutine elastic_point_solve

  !> det(I + A) of a 3 x 3 matrix (rule of Sarrus).
  !>
  !> Arguments:
  !>   a    (in)  real(rk)(3, 3)  matrix
  !> Locals:
  !>   m    real(rk)(3, 3)        I + A
  pure function det_identity_plus(a) result(det)

    real(rk), intent(in) :: a(3, 3)
    real(rk) :: det

    real(rk) :: m(3, 3)

    m = a
    m(1, 1) = m(1, 1) + 1.0_rk
    m(2, 2) = m(2, 2) + 1.0_rk
    m(3, 3) = m(3, 3) + 1.0_rk
    det = m(1, 1)*(m(2, 2)*m(3, 3) - m(2, 3)*m(3, 2)) - m(1, 2)*(m(2, 1)*m(3, 3) - &
      & m(2, 3)*m(3, 1)) + m(1, 3)*(m(2, 1)*m(3, 2) - m(2, 2)*m(3, 1))

  end function det_identity_plus

end module elastic_model_mod
