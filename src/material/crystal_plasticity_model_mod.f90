! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Crystal plasticity material model (ROADMAP step 6.2): anisotropic elasticity and power-law
!> slip kinetics with the recoupled formulation of recoupled_formulation.md, solved per point
!> for the stress, the slip strengths (hardening_mod, step 6.4), the lattice spin and the
!> orientation; returns [s], {h} and [s_tan].
!>
!> Theory : recoupled_formulation.md:
!>          (R1)  {d} = {e^e dot} + {d^p} + [W^]{e^e},  W^ from w^p by (R8);
!>          (R2)  w^p = w - (e^e d^p - d^p e^e) + (e^e e^e_0 - e^e_0 e^e)/(2 dt)  (with (R9));
!>          (R3)-(R4)  {tau} = [C]{e^e}, sigma = tau/beta, beta = det(I + e^e);
!>          (R5)-(R7)  slip kinetics (slip_kinetics_mod);  (R9) e^e_0 held in the sample frame;
!>          (R10)-(R12) [s]^-1 = beta([C]^-1/dt + [M*]), {h} = [W^]{e^e} - {e^e_0}/dt;
!>          §4 tangent: [M*_tan] in place of [M*];
!>          §5 point equation R(tau) = (1/dt)[C]^-1{tau} + {d^p}(tau) + [W^][C]^-1{tau} - {d}
!>          - {e^e_0}/dt = 0 (6 equations, the volumetric part included).
!>          The lattice rotation rate is w^p minus the slip spin sum_a gammadot^a q^a
!>          (q = skw(s (x) m)), integrated as in v2 (rotation_mod::rot_c2s_spin_update_crystal).
!> Source : port of v2-fixed:
!>          - solve_evp_mod2.f90::solve_evp_dev (state loop): start from the previous
!>            increment's stress; forward-Euler estimates of the strengths and of the
!>            orientation (continuum spin and slip spin) from the slip rates of that stress
!>            ('euler_forward'); then alternate the stress solve and the backward-Euler updates
!>            of the strengths and the orientation ('euler_backward', solve_evp_rstar order)
!>            until the stress norm and every strength change by less than tol_state times the
!>            system-1 initial strength, at most point_max_iters_state times (v2 tol_state
!>            1e-5, 100; types_mod.f90); v2 accepts when one system's strength has settled, here
!>            all must (the same for the isotropic law, whose systems all change alike);
!>          - solve_evp_stress_mod2.f90::solve_stress_evp_newton (stress solve): Newton on the
!>            residual above in the crystal frame of the current orientation, with the Jacobian
!>            (1/dt)[S] + [W^][S] + sum_a (d gammadot^a/d tau^a) p p^T ([S] = [C]^-1, w^p frozen),
!>            solved by the normal equations (linalg_mod), a diagonal check of J^T J (1e-15,
!>            v2 vtiny), a line search halving the step while the residual norm grows (fail
!>            below a factor 0.001), convergence when the residual norm is at most point_tol
!>            (v2 default 1e-4, an absolute strain rate [1/s]; a line-search failure from a
!>            start that already meets it counts as converged, a deviation), at most
!>            point_max_iters_newton
!>            iterations (100); the slip rates of the Newton use the cutoff tiny(1.0)**m of v2's
!>            Newton (|tau/g| at or below it: no slip), w^p recomputed by (R2) at each residual
!>            (v2 dp_wp_hat_phase);
!>          - solve_evp_rstar_mod2.f90 (compute_sliprate with the cutoff t_min of set_tmin and
!>            |gammadot| <= 1e-30 set to 0; update_rstar; update_ori);
!>          - solveit_evp_mod4.f90::material_matrix_evp_aniso ([s], [s_tan], {h} from the final
!>            stress, orientation and w^p (kinematics_mod_bis.f90::find_wp_hat_phase), slip
!>            rates with the cutoff t_min).
!>          Deviations (ROADMAP D49, step 6.2): one 6-equation system (v2: 5 deviatoric + a
!>          separate volumetric update), so the state-loop norm is that of the full Kirchhoff
!>          6-vector (v2: deviatoric); the forward estimate uses the start-of-increment crystal
!>          frame (as the elastic model, elastic_model_mod); a point whose Newton or state loop
!>          does not converge keeps its last iterate and is reported (v2 replaced its stress by
!>          the average of the converged points of the rank, which depends on the number of
!>          ranks); m and gammadot_0 per slip system (D48); the threshold strength is the point's
!>          phase strength of system 1 (v2: that of the rank's first element, CHANGES.md C36).
!>
!>          State-loop damping (step 7.3, D62, dev/theory/addendum_state_damping.md):
!>          with irradiation softening only, the state x = (rho, H) is damped between state
!>          iterations, x_{k+1} = x_k + omega_k r_k with r_k = x~_k - x_k (x~_k: the
!>          backward-Euler update), omega_1 = 1 and the Aitken factor of Kuttler & Wall (2008),
!>          Comput. Mech. 43:61-72, Eqs. 33, 35 and 44 (recursion on omega kept, Remark 8),
!>          bounded to [0.001, 1] (an estimate <= 0, no overshoot, gives 1), reset to 1 at every point solve (their 'relaxation'); the convergence test uses the
!>          undamped strengths crss(x~_k). Every other phase: the loop above, unchanged.
!>
!> Names (ROADMAP D50): point_tol, point_max_iters_newton and point_max_iters_state are v2's
!> sx_tol, sx_max_iters_newton and sx_max_iters_state.
!>
!> Work arrays (ROADMAP D100): the per-point routines take a cp_work_t holding every array they
!> need at its exact size (slip-system and vertex arrays, the hardening work arrays with the
!> point's hardening values), allocated once per call of cp_solve_points (cp_work_alloc), so
!> that nothing is allocated per point.
!>
!> Public types:
!>   cp_model_t   the model (extends material_model_t; components at the type)
!>   cp_work_t    work arrays of the per-point routines (components at the type)
!>
!> Public routines:
!>   cp_model_init    set up the model of a phase
!>   cp_work_alloc    allocate the work arrays for a model
!>   cp_point_solve   the per-point kernel (pure)
!>   cp_vertex_estimate  v2's increment-1 viscoplastic stress estimate from a yield-surface
!>                    vertex (pure; initialization isotropic_viscoplastic, step 6.6)
!>   cp_residual      the point residual (R1) at a stress (pure; also for tests)
!>   cp_state_damping_factor  the bounded Aitken factor of the state-loop damping (pure;
!>                    also for tests; step 7.3, D62)
module crystal_plasticity_model_mod

  use kinds_mod, only: rk, ik
  use tensor_mod, only: matrix_from_mandel, mandel_norm, skew_from_matrix, spin_commutator, &
    & mandel_dev, equiv_strain, dev5_from_mandel, mandel_from_dev5, mat5_from_mat6
  use rotation_mod, only: transform_crystal_to_sample_mat6, transform_crystal_to_sample_mandel, &
    & transform_sample_to_crystal_mandel, transform_crystal_to_sample_skew, &
    & transform_sample_to_crystal_skew, rot_c2s_spin_update_crystal
  use linalg_mod, only: ldlt_factor, ldlt_substitute, normal_equations_solve
  use slip_kinetics_mod, only: slip_rss_cutoff, resolved_shear_stress, slip_power_law, &
    & plastic_compliance
  use hardening_mod, only: hardening_params_t, hardening_point_set, hardening_work_t, &
    & hardening_work_alloc, hardening_forward, &
    & hardening_backward, hardening_kocks_mecking, hardening_km_strengths, &
    & hardening_km_k2, hardening_km_rate, hardening_km_forward, hardening_km_backward, &
    & hardening_irr_t, hardening_irr_dmg_0, hardening_irr_strengths, hardening_irr_rate, &
    & hardening_irr_forward, hardening_irr_backward
  use material_state_mod, only: material_state_t
  use material_model_mod, only: material_model_t
  use error_mod, only: internal_error

  implicit none

  private

  public :: cp_model_t, cp_work_t, cp_model_init, cp_model_set_element_values, cp_work_alloc, &
    & cp_point_solve, cp_residual, cp_vertex_estimate, cp_state_damping_factor
  public :: cp_converged, cp_newton_failed, cp_state_not_converged, cp_singular, &
    & cp_hardening_failed

  ! Point-solve status codes (cp_point_solve).
  integer, parameter :: cp_converged = 0, cp_newton_failed = 1, cp_state_not_converged = 2, &
    & cp_singular = 3, cp_hardening_failed = 4

  ! v2's limits: the smallest diagonal of J^T J (general_mod.f90 vtiny), the smallest line-search
  ! factor, and the slip rates set to zero (compute_sliprate).
  real(rk), parameter :: diag_min = 1.0e-15_rk
  real(rk), parameter :: fact_min = 1.0e-3_rk
  real(rk), parameter :: sliprate_zero = 1.0e-30_rk

  ! v2's increment-1 viscoplastic estimate (solve_evp_vpstress): the rate sensitivity of its
  ! power law (xm_fake), the number of vertices tried (ntrials) and the fraction of the estimate
  ! that starts the point solve (solve_evp_dev).
  real(rk), parameter :: vertex_m = 0.4_rk
  integer, parameter :: vertex_num_trials = 4
  real(rk), parameter :: vertex_start_fraction = 0.5_rk

  ! State-loop damping (D62, addendum_state_damping.md §3-§4): the lower bound of the
  ! factor omega, and the relative size of the residual change below which the previous factor
  ! is kept (no information in the last two residuals).
  real(rk), parameter :: damp_omega_min = 0.001_rk
  real(rk), parameter :: damp_change_min = 1.0e-14_rk

  !> Crystal plasticity model of one phase.
  !>
  !> Components (in addition to material_model_t's num_slip):
  !>   stiff_el_c        real(rk)(6, 6)         elastic stiffness [C], Mandel, crystal frame [MPa]
  !>   compl_el_c        real(rk)(6, 6)         elastic compliance [S] = [C]^-1 [1/MPa]
  !>   schmid_sym_c      real(rk)(6, num_slip)  Schmid tensors p, Mandel, crystal frame
  !>   schmid_skew_c     real(rk)(3, num_slip)  skew parts q (w23 w13 w12), crystal frame
  !>   m                 real(rk)(num_slip)     rate sensitivities [-]
  !>   gammadot_0        real(rk)(num_slip)     reference slip rates [1/s]
  !>   tau_0             real(rk)(num_slip)     initial strengths [MPa]
  !>   rss_cutoff        real(rk)(num_slip)     cutoff t_min of |tau/g| of the slip rates and
  !>                                            matrices (slip_rss_cutoff) [-]
  !>   rss_cutoff_newton real(rk)(num_slip)     cutoff tiny(1.0)**m of the Newton's slip rates [-]
  !>   tol_state         real(rk)               state-loop tolerance, times tau_0(1)
  !>                                            (point_tol_state) [-]
  !>   max_iters_state   integer                state-loop limit (point_max_iters_state)
  !>   max_iters_newton  integer                Newton limit (point_max_iters_newton)
  !>   tol_newton        real(rk)               Newton residual tolerance (point_tol) [1/s]
  !>   hard              type(hardening_params_t)  hardening law and parameters (hardening_mod;
  !>                                            default: none, constant strengths)
  !>   vertices_c        real(rk)(6, :)         yield-surface vertices, one of each +- pair,
  !>                                            deviatoric Mandel, crystal frame (only with
  !>                                            initialization isotropic_viscoplastic, else
  !>                                            size 0): stresses computed from the phase's
  !>                                            systems and tau_0 (yield_vertices_compute) [MPa]
  !>   elt_tau_0         real(rk)(num_slip, :)  initial strengths of each local element from
  !>                                            simulation.opt (step 7.1, D59), or size 0 in
  !>                                            the second dimension: the phase's tau_0 [MPa]
  !>   elt_tau_s         real(rk)(num_slip, :)  the same for the saturation strength (or tau_s0)
  !>                                            [MPa]
  !>   elt_h_0           real(rk)(num_slip, :)  the same for the hardening rate coefficient
  !>                                            [MPa]
  !>   elt_rho_0         real(rk)(:)            initial network dislocation density of each
  !>                                            local element from simulation.opt (Kocks-Mecking,
  !>                                            step 7.2, D60), or size 0: the phase's km_rho_0
  !>                                            [1/mm^2]
  !>   irradiation       logical                irradiation softening on (step 7.3, D61;
  !>                                            with Kocks-Mecking only)
  !>   irr               type(hardening_irr_t)  its parameters (plane_c set by the caller)
  !>   elt_loop_density  real(rk)(:)            loop number density of each local element
  !>                                            from simulation.opt, or size 0: the phase's
  !>                                            [1/mm^3]
  !>   elt_loop_diameter real(rk)(:)            the same for the loop diameter [mm]
  !>   damp_state       logical                Aitken damping of the state loop with
  !>                                            irradiation (D62; default on; the tests turn
  !>                                            it off to compare with the plain loop)
  type, extends(material_model_t) :: cp_model_t
    real(rk) :: stiff_el_c(6, 6) = 0.0_rk
    real(rk) :: compl_el_c(6, 6) = 0.0_rk
    real(rk), allocatable :: schmid_sym_c(:, :)
    real(rk), allocatable :: schmid_skew_c(:, :)
    real(rk), allocatable :: m(:)
    real(rk), allocatable :: gammadot_0(:)
    real(rk), allocatable :: tau_0(:)
    real(rk), allocatable :: rss_cutoff(:)
    real(rk), allocatable :: rss_cutoff_newton(:)
    real(rk) :: tol_state = 1.0e-5_rk
    integer :: max_iters_state = 100
    integer :: max_iters_newton = 100
    real(rk) :: tol_newton = 1.0e-4_rk
    type(hardening_params_t) :: hard
    real(rk), allocatable :: vertices_c(:, :)
    real(rk), allocatable :: elt_tau_0(:, :)
    real(rk), allocatable :: elt_tau_s(:, :)
    real(rk), allocatable :: elt_h_0(:, :)
    real(rk), allocatable :: elt_rho_0(:)
    logical :: irradiation = .false.
    type(hardening_irr_t) :: irr
    real(rk), allocatable :: elt_loop_density(:)
    real(rk), allocatable :: elt_loop_diameter(:)
    logical :: damp_state = .true.
  contains
    procedure :: init_state => cp_init_state
    procedure :: solve_points => cp_solve_points
  end type cp_model_t

  !> Work arrays of the per-point routines of a phase, at the exact size (num_slip systems, the
  !> model's vertices), allocated once per call of cp_solve_points by cp_work_alloc (ROADMAP
  !> D100). Their values are not kept from one call to the next, except hardening%point, which
  !> the caller sets before each point solve (cp_work_alloc sets the phase's values).
  !>
  !> Components:
  !>   num_slip         integer                 number of slip systems N
  !>   hardening        type(hardening_work_t)  hardening work arrays and the point's hardening
  !>                                            values (hardening%point)
  !>   tau_0_pt         real(rk)(N)  cp_point_solve: initial strengths used [MPa]
  !>   rss, coef_sec, coef_tan  real(rk)(N)  cp_point_solve: resolved shear stresses [MPa],
  !>                                         secant and tangent coefficients [1/(MPa s)]
  !>   crss_iter        real(rk)(N)  cp_point_solve: strengths of the previous state
  !>                                 iteration [MPa]
  !>   newton_coef_tan  real(rk)(N)  newton_stress: tangent coefficients [1/(MPa s)]
  !>   residual_rss, residual_sliprate, residual_coef_sec  real(rk)(N)  newton_stress: the work
  !>                                 arrays of cp_residual [MPa], [1/s], [1/(MPa s)]
  !>   vertex_work      real(rk)(num_vertices)  cp_vertex_estimate: virtual plastic work of each
  !>                                 vertex not yet used [MPa/s] (size 0 without vertices)
  !>   vertex_m, vertex_crss  real(rk)(N)  cp_vertex_estimate: m = 0.4 and crss_1 for every
  !>                                 system [-], [MPa] (size 0 without vertices)
  !>   vertex_rss, vertex_rate, vertex_coef_s, vertex_coef_t  real(rk)(N)  its residual's
  !>                                 resolved shear stresses [MPa], slip rates [1/s], secant and
  !>                                 tangent coefficients [1/(MPa s)] (size 0 without vertices)
  type :: cp_work_t
    integer :: num_slip = 0
    type(hardening_work_t) :: hardening
    real(rk), allocatable :: tau_0_pt(:)
    real(rk), allocatable :: rss(:)
    real(rk), allocatable :: coef_sec(:)
    real(rk), allocatable :: coef_tan(:)
    real(rk), allocatable :: crss_iter(:)
    real(rk), allocatable :: newton_coef_tan(:)
    real(rk), allocatable :: residual_rss(:)
    real(rk), allocatable :: residual_sliprate(:)
    real(rk), allocatable :: residual_coef_sec(:)
    real(rk), allocatable :: vertex_work(:)
    real(rk), allocatable :: vertex_m(:)
    real(rk), allocatable :: vertex_crss(:)
    real(rk), allocatable :: vertex_rss(:)
    real(rk), allocatable :: vertex_rate(:)
    real(rk), allocatable :: vertex_coef_s(:)
    real(rk), allocatable :: vertex_coef_t(:)
  end type cp_work_t

contains

  !> Set up the model of a phase.
  !>
  !> Arguments:
  !>   model             (out) type(cp_model_t)  model
  !>   stiff_el_c        (in)  real(rk)(6, 6)    [C], Mandel, crystal frame [MPa]
  !>   schmid_sym_c      (in)  real(rk)(6, :)    Schmid tensors p of the active systems
  !>                                             (crystal_mod), crystal frame
  !>   schmid_skew_c     (in)  real(rk)(3, :)    their skew parts q, crystal frame
  !>   m, gammadot_0, tau_0  (in)  real(rk)(:)   per-system rate sensitivity [-], reference
  !>                                             slip rate [1/s], initial strength [MPa]
  !>   tol_newton        (in, optional) real(rk) point_tol (default 1e-4, v2) [1/s]
  !>   max_iters_newton  (in, optional) integer  point_max_iters_newton (default 100, v2)
  !>   max_iters_state   (in, optional) integer  point_max_iters_state (default 100, v2)
  !>   tol_state         (in, optional) real(rk) state-loop tolerance (default 1e-5, v2) [-]
  !>   hard              (in, optional) type(hardening_params_t)  hardening (default: none)
  !>   vertices_c        (in, optional) real(rk)(6, :)  yield-surface vertices (initialization
  !>                                             isotropic_viscoplastic; default: none) [MPa]
  !> Locals:
  !>   n       integer      number of slip systems
  !>   ios     integer      inversion status
  subroutine cp_model_init(model, stiff_el_c, schmid_sym_c, schmid_skew_c, m, gammadot_0, &
    & tau_0, tol_newton, max_iters_newton, max_iters_state, tol_state, hard, vertices_c, irr)

    type(cp_model_t), intent(out) :: model
    real(rk), intent(in) :: stiff_el_c(6, 6)
    real(rk), intent(in) :: schmid_sym_c(:, :)
    real(rk), intent(in) :: schmid_skew_c(:, :)
    real(rk), intent(in) :: m(:)
    real(rk), intent(in) :: gammadot_0(:)
    real(rk), intent(in) :: tau_0(:)
    real(rk), intent(in), optional :: tol_newton
    integer, intent(in), optional :: max_iters_newton
    integer, intent(in), optional :: max_iters_state
    real(rk), intent(in), optional :: tol_state
    type(hardening_params_t), intent(in), optional :: hard
    real(rk), intent(in), optional :: vertices_c(:, :)
    type(hardening_irr_t), intent(in), optional :: irr

    integer :: n, ios

    ! Slip systems and their parameters.
    n = size(m)
    model%num_slip = int(n, ik)
    model%stiff_el_c = stiff_el_c
    allocate (model%schmid_sym_c(6, n), model%schmid_skew_c(3, n))
    allocate (model%m(n), model%gammadot_0(n), model%tau_0(n))
    allocate (model%rss_cutoff(n), model%rss_cutoff_newton(n))
    model%schmid_sym_c = schmid_sym_c
    model%schmid_skew_c = schmid_skew_c
    model%m = m
    model%gammadot_0 = gammadot_0
    model%tau_0 = tau_0

    ! The two cutoffs of v2: t_min (set_tmin) and tiny**m (solve_stress_evp_newton).
    model%rss_cutoff = slip_rss_cutoff(m)
    model%rss_cutoff_newton = tiny(1.0_rk)**m

    ! Elastic compliance (the stiffness is positive definite, checked by the configuration).
    call invert_spd(stiff_el_c, model%compl_el_c, ios)
    if (ios /= 0) call internal_error('cp_model_init', 'the elastic stiffness is singular')

    ! Controls.
    if (present(tol_newton)) model%tol_newton = tol_newton
    if (present(max_iters_newton)) model%max_iters_newton = max_iters_newton
    if (present(max_iters_state)) model%max_iters_state = max_iters_state
    if (present(tol_state)) model%tol_state = tol_state
    if (present(hard)) model%hard = hard
    ! Irradiation softening (step 7.3): its parameters, with the slip planes N^a.
    if (present(irr)) then
      model%irradiation = .true.
      model%irr = irr
      if (size(model%irr%plane_c, 2) /= n) call internal_error('cp_model_init', &
        & 'the irradiation slip planes do not match the slip systems')
    end if

    ! Yield-surface vertices of the increment-1 estimate (initialization
    ! isotropic_viscoplastic).
    if (present(vertices_c)) then
      model%vertices_c = vertices_c
    else
      allocate (model%vertices_c(6, 0))
    end if

    ! No per-element values (cp_model_set_element_values).
    allocate (model%elt_tau_0(n, 0), model%elt_tau_s(n, 0), model%elt_h_0(n, 0))
    allocate (model%elt_rho_0(0), model%elt_loop_density(0), model%elt_loop_diameter(0))

  end subroutine cp_model_init

  !> Per-element values from simulation.opt (step 7.1, D59; addendum_hardening_modifiers.md §6):
  !> each given array replaces the phase's value for every local element (columns indexed by
  !> the local element number; elements of other phases unused).
  !>
  !> Arguments:
  !>   model  (inout) type(cp_model_t)        model
  !>   tau_0  (in, optional) real(rk)(num_slip, num_elts)  initial strengths [MPa]
  !>   tau_s  (in, optional) real(rk)(num_slip, num_elts)  saturation strengths (or tau_s0) [MPa]
  !>   h_0    (in, optional) real(rk)(num_slip, num_elts)  hardening rate coefficients [MPa]
  !>   rho_0  (in, optional) real(rk)(num_elts)  initial network dislocation densities
  !>                                             (Kocks-Mecking, step 7.2) [1/mm^2]
  !>   loop_density   (in, optional) real(rk)(num_elts)  loop number densities (irradiation,
  !>                                             step 7.3) [1/mm^3]
  !>   loop_diameter  (in, optional) real(rk)(num_elts)  loop diameters [mm]
  subroutine cp_model_set_element_values(model, tau_0, tau_s, h_0, rho_0, loop_density, &
    & loop_diameter)

    type(cp_model_t), intent(inout) :: model
    real(rk), intent(in), optional :: tau_0(:, :)
    real(rk), intent(in), optional :: tau_s(:, :)
    real(rk), intent(in), optional :: h_0(:, :)
    real(rk), intent(in), optional :: rho_0(:)
    real(rk), intent(in), optional :: loop_density(:)
    real(rk), intent(in), optional :: loop_diameter(:)

    if (present(tau_0)) model%elt_tau_0 = tau_0
    if (present(tau_s)) model%elt_tau_s = tau_s
    if (present(h_0)) model%elt_h_0 = h_0
    if (present(rho_0)) model%elt_rho_0 = rho_0
    if (present(loop_density)) model%elt_loop_density = loop_density
    if (present(loop_diameter)) model%elt_loop_diameter = loop_diameter

  end subroutine cp_model_set_element_values

  !> Allocate the work arrays of the per-point routines for a model, at the exact size (once
  !> per call of cp_solve_points, outside the point loops; ROADMAP D100); arrays already of the
  !> right size are kept. The point's hardening values are set to the phase's.
  !>
  !> Arguments:
  !>   work   (inout) type(cp_work_t)   work arrays
  !>   model  (in)    type(cp_model_t)  model of the phase
  !> Locals:
  !>   n, num_vertices  integer         slip systems, vertices (0 without a vertex table)
  pure subroutine cp_work_alloc(work, model)

    type(cp_work_t), intent(inout) :: work
    type(cp_model_t), intent(in) :: model

    integer :: n, num_vertices

    ! Sizes.
    n = int(model%num_slip)
    num_vertices = 0
    if (allocated(model%vertices_c)) num_vertices = size(model%vertices_c, 2)
    work%num_slip = n

    ! The hardening work arrays (the evolved variable's solver state per law; the damage tensor
    ! with irradiation, which needs Kocks-Mecking) and the point's values.
    call hardening_work_alloc(work%hardening, model%hard, n, &
      & model%irradiation .and. model%hard%law == hardening_kocks_mecking)

    ! The slip-system arrays of the point solve, Newton and the residual.
    call alloc_vector(work%tau_0_pt, n)
    call alloc_vector(work%rss, n)
    call alloc_vector(work%coef_sec, n)
    call alloc_vector(work%coef_tan, n)
    call alloc_vector(work%crss_iter, n)
    call alloc_vector(work%newton_coef_tan, n)
    call alloc_vector(work%residual_rss, n)
    call alloc_vector(work%residual_sliprate, n)
    call alloc_vector(work%residual_coef_sec, n)

    ! The vertex estimate's arrays (empty without vertices).
    call alloc_vector(work%vertex_work, num_vertices)
    call alloc_vector(work%vertex_m, merge(n, 0, num_vertices > 0))
    call alloc_vector(work%vertex_crss, merge(n, 0, num_vertices > 0))
    call alloc_vector(work%vertex_rss, merge(n, 0, num_vertices > 0))
    call alloc_vector(work%vertex_rate, merge(n, 0, num_vertices > 0))
    call alloc_vector(work%vertex_coef_s, merge(n, 0, num_vertices > 0))
    call alloc_vector(work%vertex_coef_t, merge(n, 0, num_vertices > 0))

  contains

    !> An allocatable vector of size m (kept if it already has it).
    !>
    !> Arguments:
    !>   vec  (inout) real(rk)(:), allocatable  the vector
    !>   m    (in)    integer                   its size
    pure subroutine alloc_vector(vec, m)

      real(rk), allocatable, intent(inout) :: vec(:)
      integer, intent(in) :: m

      if (allocated(vec)) then
        if (size(vec) == m) return
        deallocate (vec)
      end if
      allocate (vec(m))

    end subroutine alloc_vector

  end subroutine cp_work_alloc

  !> Initial state at the phase's points: zero stress, elastic strain and lattice spin; slip
  !> strengths tau_0, or the element's from simulation.opt (rows beyond num_slip zero); with
  !> Kocks-Mecking hardening (step 7.2), the density rho_0 (the phase's km_rho_0 or the
  !> element's from simulation.opt) and the strengths tau_0 + b mu_0 sqrt(h_n rho_0)
  !> (addendum_kocks_mecking.md §5); with irradiation (step 7.3), H_0 from the element's loop
  !> density and diameter (I4) and the strengths (I1) (addendum_irradiation.md §7).
  !> Orientations are set by the caller.
  !>
  !> Arguments:
  !>   self   (in)    class(cp_model_t)        model
  !>   elts   (in)    integer(ik)(:)           elements of the phase
  !>   state  (inout) type(material_state_t)   state (crss with at least num_slip rows)
  !> Locals:
  !>   i, i_qpt  integer(ik)  element list index, point
  !>   rho_0     real(rk)     initial density of the element [1/mm^2]
  !>   loop_density, loop_diameter  real(rk)  the element's loops [1/mm^3], [mm]
  subroutine cp_init_state(self, elts, state)

    class(cp_model_t), intent(in) :: self
    integer(ik), intent(in) :: elts(:)
    type(material_state_t), intent(inout) :: state

    integer(ik) :: i, i_qpt
    real(rk) :: rho_0, loop_density, loop_diameter

    ! Zero mechanical state.
    state%stress_s(:, :, elts) = 0.0_rk
    state%strain_el_s(:, :, elts) = 0.0_rk
    state%spin_lat_s(:, :, elts) = 0.0_rk

    ! Initial strengths.
    if (state%num_slip < self%num_slip) call internal_error('cp_init_state', &
      & 'the state has fewer slip-system rows than the model')
    do i = 1_ik, int(size(elts), ik)
      do i_qpt = 1_ik, state%num_qpts
        state%crss(:, i_qpt, elts(i)) = 0.0_rk
        if (size(self%elt_tau_0, 2) > 0) then
          state%crss(1:self%num_slip, i_qpt, elts(i)) = self%elt_tau_0(:, elts(i))
        else
          state%crss(1:self%num_slip, i_qpt, elts(i)) = self%tau_0
        end if
        ! Kocks-Mecking: the initial density and the strengths it gives.
        if (self%hard%law == hardening_kocks_mecking) then
          rho_0 = self%hard%km_rho_0
          if (size(self%elt_rho_0) > 0) rho_0 = self%elt_rho_0(elts(i))
          state%km_rho(i_qpt, elts(i)) = rho_0
          state%crss(1:self%num_slip, i_qpt, elts(i)) = hardening_km_strengths(self%hard, &
            & state%crss(1:self%num_slip, i_qpt, elts(i)), rho_0)
        end if
        ! Irradiation: H_0 from the element's loops, and the defect term of the strengths
        ! (the Kocks-Mecking strengths recomputed with it from tau_0).
        if (self%irradiation) then
          loop_density = self%irr%loop_density
          if (size(self%elt_loop_density) > 0) loop_density = self%elt_loop_density(elts(i))
          loop_diameter = self%irr%loop_diameter
          if (size(self%elt_loop_diameter) > 0) loop_diameter = &
            & self%elt_loop_diameter(elts(i))
          state%dmg_c(:, i_qpt, elts(i)) = hardening_irr_dmg_0(self%irr, loop_density, &
            & loop_diameter)
          if (size(self%elt_tau_0, 2) > 0) then
            state%crss(1:self%num_slip, i_qpt, elts(i)) = hardening_irr_strengths(self%hard, &
              & self%irr, self%elt_tau_0(:, elts(i)), rho_0, state%dmg_c(:, i_qpt, elts(i)))
          else
            state%crss(1:self%num_slip, i_qpt, elts(i)) = hardening_irr_strengths(self%hard, &
              & self%irr, self%tau_0, rho_0, state%dmg_c(:, i_qpt, elts(i)))
          end if
        end if
      end do
    end do

  end subroutine cp_init_state

  !> Point solve at every quadrature point of the given elements (material_model_t contract),
  !> strengths included (hardening_mod, step 6.4), and the slip output quantities of the state
  !> (step 6.5): slip rates and resolved shear stresses of the final stress, accumulated slip
  !> and plastic strain, plastic deformation rate and plastic spin.
  !>
  !> Theory : slip = slip_prev + gammadot dt (signed; backward Euler with the end rates);
  !>          strain_pl_c = strain_pl_c_prev + {d^p}_c dt with {d^p}_c = sum gammadot^a {p^a}
  !>          in the crystal frame (accumulated in the lattice frame and written in the sample
  !>          frame of the current orientation); {d^p}_s and the plastic spin
  !>          sum gammadot^a q^a in the sample frame of the end orientation.
  !> Source : port of v2-fixed finalize_res_mod1.f90::finalize_res (results%slip and
  !>          results%strain_pl += rate * dtime, with dp_hat from finalize_res_dphat_wphat in
  !>          the lattice frame) and write_res_mod1.f90 (strain_pl, defrate_pl and spinrate
  !>          transformed to the sample frame with the current orientation).
  !>
  !> Arguments: see material_model_mod::solve_points_iface.
  !> Locals:
  !>   i, i_elt, i_qpt  integer(ik)            element list index, element, point
  !>   status           integer                point status (cp_converged: converged)
  !>   n                integer                slip systems
  !>   sliprate         real(rk)(num_slip)     slip rates at the final stress [1/s]
  !>   rss              real(rk)(num_slip)     resolved shear stresses at the final stress [MPa]
  !>   tangent          real(rk)(6, 6)         [s_tan] of the point [MPa s]
  !>   defrate_pl_c     real(rk)(6)            {d^p}, Mandel, crystal frame [1/s]
  !>   tau_est_c        real(rk)(6)            increment-1 vertex estimate (initialization
  !>                                           isotropic_viscoplastic), crystal frame [MPa]
  !>   est_status       integer                its status (0: converged; else v2's fallback,
  !>                                           the vertex stress)
  !>   tau_0_pt         real(rk)(num_slip)     initial strengths of the point's element (the
  !>                                           phase's or from simulation.opt) [MPa]
  !>   spin_pl_c        real(rk)(3)            plastic spin sum gammadot q, crystal frame [1/s]
  !>   work             type(cp_work_t)        work arrays of the per-point routines (once per
  !>                                           call, D100), with the per-point hardening values
  !>                                           in work%hardening%point: the element's tau_s and
  !>                                           h_0, the accumulated slip and d_eq (step 7.1)
  subroutine cp_solve_points(self, elts, dt, defrate_s, spin_s, state_prev, state, &
    & mat_secant_s, vec_h_s, num_failed, mat_tangent_s)

    class(cp_model_t), intent(in) :: self
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
    integer :: status, n
    real(rk) :: sliprate(self%num_slip), rss(self%num_slip), tangent(6, 6), defrate_pl_c(6)
    real(rk) :: tau_est_c(6), tau_0_pt(self%num_slip)
    real(rk) :: spin_pl_c(3)
    integer :: est_status
    type(cp_work_t) :: work

    n = int(self%num_slip)
    num_failed = 0_ik
    call cp_work_alloc(work, self)
    do i = 1_ik, int(size(elts), ik)
      i_elt = elts(i)
      do i_qpt = 1_ik, state%num_qpts
        ! The density starts from its start-of-increment value (Kocks-Mecking; zero and unused
        ! otherwise).
        state%km_rho(i_qpt, i_elt) = state_prev%km_rho(i_qpt, i_elt)
        state%dmg_c(:, i_qpt, i_elt) = state_prev%dmg_c(:, i_qpt, i_elt)

        ! Per-point values (step 7.1): the element's tau_0, tau_s, h_0 (phase or simulation.opt),
        ! the start-of-increment accumulated slip and d_eq = sqrt(2/3 d':d') (cyclic hardening).
        tau_0_pt = self%tau_0
        if (size(self%elt_tau_0, 2) > 0) tau_0_pt = self%elt_tau_0(:, i_elt)
        call hardening_point_set(work%hardening%point, self%hard, &
          & slip_accum=state_prev%slip_accum(1:n, i_qpt, i_elt), &
          & defrate_eq=equiv_strain(defrate_s(:, i_qpt, i_elt)))
        if (size(self%elt_tau_s, 2) > 0) work%hardening%point%tau_s(:) = self%elt_tau_s(:, i_elt)
        if (size(self%elt_h_0, 2) > 0) work%hardening%point%h_0(:) = self%elt_h_0(:, i_elt)

        ! Per-point solve from the start-of-increment state (strengths included, step 6.4); in
        ! increment 1 with initialization isotropic_viscoplastic, from the fraction 0.5 of
        ! v2's viscoplastic vertex estimate (computed at every global iteration, as v2
        ! solve_evp_dev).
        if (state_prev%start_from_vertex) then
          call cp_vertex_estimate(self, state_prev%crss(1:n, i_qpt, i_elt), &
            & transform_sample_to_crystal_mandel(state_prev%rot_c2s(:, :, i_qpt, i_elt), &
            & defrate_s(:, i_qpt, i_elt)), tau_est_c, est_status, work)
          tau_est_c = vertex_start_fraction*tau_est_c
          call cp_point_solve(self, state_prev%crss(1:n, i_qpt, i_elt), &
            & state_prev%rot_c2s(:, :, i_qpt, i_elt), state_prev%stress_s(:, i_qpt, i_elt), &
            & state_prev%strain_el_s(:, i_qpt, i_elt), defrate_s(:, i_qpt, i_elt), &
            & spin_s(:, i_qpt, i_elt), dt, state%stress_s(:, i_qpt, i_elt), &
            & state%strain_el_s(:, i_qpt, i_elt), state%spin_lat_s(:, i_qpt, i_elt), &
            & state%rot_c2s(:, :, i_qpt, i_elt), state%crss(1:n, i_qpt, i_elt), &
            & mat_secant_s(:, :, i_qpt, i_elt), vec_h_s(:, i_qpt, i_elt), tangent, sliprate, &
            & status, rss_final=rss, &
            & spin_hat_update_s=state%spin_hat_update_s(:, i_qpt, i_elt), &
            & spin_slip_update_s=state%spin_slip_update_s(:, i_qpt, i_elt), &
            & tau_start_c=tau_est_c, tau_0_point=tau_0_pt, km_rho=state%km_rho(i_qpt, i_elt), &
            & dmg_c=state%dmg_c(:, i_qpt, i_elt), work=work)
        else
          call cp_point_solve(self, state_prev%crss(1:n, i_qpt, i_elt), &
            & state_prev%rot_c2s(:, :, i_qpt, i_elt), state_prev%stress_s(:, i_qpt, i_elt), &
            & state_prev%strain_el_s(:, i_qpt, i_elt), defrate_s(:, i_qpt, i_elt), &
            & spin_s(:, i_qpt, i_elt), dt, state%stress_s(:, i_qpt, i_elt), &
            & state%strain_el_s(:, i_qpt, i_elt), state%spin_lat_s(:, i_qpt, i_elt), &
            & state%rot_c2s(:, :, i_qpt, i_elt), state%crss(1:n, i_qpt, i_elt), &
            & mat_secant_s(:, :, i_qpt, i_elt), vec_h_s(:, i_qpt, i_elt), tangent, sliprate, &
            & status, rss_final=rss, &
            & spin_hat_update_s=state%spin_hat_update_s(:, i_qpt, i_elt), &
            & spin_slip_update_s=state%spin_slip_update_s(:, i_qpt, i_elt), &
            & tau_0_point=tau_0_pt, km_rho=state%km_rho(i_qpt, i_elt), &
            & dmg_c=state%dmg_c(:, i_qpt, i_elt), work=work)
        end if
        if (status /= cp_converged) num_failed = num_failed + 1_ik

        ! Kocks-Mecking (step 7.2): k_2 and the density rate at the slip rates of the final
        ! stress (zero at Gammadot = 0); km_rho was updated by the point solve.
        if (self%hard%law == hardening_kocks_mecking) then
          state%km_k2(i_qpt, i_elt) = hardening_km_k2(self%hard, sum(abs(sliprate)))
          state%km_rhorate(i_qpt, i_elt) = hardening_km_rate(self%hard, sum(abs(sliprate)), &
            & state%km_rho(i_qpt, i_elt))
        end if
        ! Irradiation (step 7.3): H dot at the slip rates of the final stress and the updated H.
        if (self%irradiation) call hardening_irr_rate(self%irr, sliprate, &
          & state%dmg_c(:, i_qpt, i_elt), state%dmgrate_c(:, i_qpt, i_elt))
        if (present(mat_tangent_s)) mat_tangent_s(:, :, i_qpt, i_elt) = tangent

        ! Slip quantities (step 6.5): rates and resolved shear stresses of the final stress,
        ! slip and plastic strain accumulated from the start of the increment (plastic strain in
        ! the crystal frame), plastic deformation rate and spin in the sample frame.
        state%sliprate(1:n, i_qpt, i_elt) = sliprate
        state%rss(1:n, i_qpt, i_elt) = rss
        state%slip(1:n, i_qpt, i_elt) = state_prev%slip(1:n, i_qpt, i_elt) + sliprate*dt
        state%slip_accum(1:n, i_qpt, i_elt) = state_prev%slip_accum(1:n, i_qpt, i_elt) + &
          & abs(sliprate)*dt
        defrate_pl_c = matmul(self%schmid_sym_c, sliprate)
        state%strain_pl_c(:, i_qpt, i_elt) = state_prev%strain_pl_c(:, i_qpt, i_elt) + &
          & defrate_pl_c*dt
        state%defrate_pl_s(:, i_qpt, i_elt) = transform_crystal_to_sample_mandel( &
          & state%rot_c2s(:, :, i_qpt, i_elt), defrate_pl_c)
        spin_pl_c = matmul(self%schmid_skew_c, sliprate)
        state%spin_pl_s(:, i_qpt, i_elt) = transform_crystal_to_sample_skew( &
          & state%rot_c2s(:, :, i_qpt, i_elt), spin_pl_c)
      end do
    end do

  end subroutine cp_solve_points

  !> Crystal plasticity point solve (module header): state loop of stress solves (Newton),
  !> strength and orientation updates, then [s], {h}, [s_tan] (R10)-(R12), §4.
  !>
  !> Theory : module header.
  !> Source : port of v2 solve_evp_dev, solve_stress_evp_newton, solve_evp_rstar and
  !>          material_matrix_evp_aniso (module header).
  !>
  !> Arguments:
  !>   model            (in)  type(cp_model_t)  model of the phase
  !>   crss_prev        (in)  real(rk)(num_slip) slip strengths at the start of the
  !>                                          increment [MPa]
  !>   rot_c2s_prev     (in)  real(rk)(3, 3)  orientation at the start of the increment [-]
  !>   stress_s_prev    (in)  real(rk)(6)     Cauchy stress at the start, sample frame [MPa]
  !>   strain_el_s_prev (in)  real(rk)(6)     {e^e_0}, start of the increment, sample frame [-]
  !>   defrate_s        (in)  real(rk)(6)     {d}, Mandel, sample frame [1/s]
  !>   spin_s           (in)  real(rk)(3)     w (w23 w13 w12), sample frame [1/s]
  !>   dt               (in)  real(rk)        time increment [s]
  !>   stress_s         (out) real(rk)(6)     Cauchy stress, sample frame [MPa]
  !>   strain_el_s      (out) real(rk)(6)     {e^e}, sample frame [-]
  !>   spin_lat_s       (out) real(rk)(3)     w^p (w23 w13 w12), sample frame [1/s]
  !>   rot_c2s          (out) real(rk)(3, 3)  orientation at the end of the increment [-]
  !>   crss             (out) real(rk)(num_slip) slip strengths at the end of the increment
  !>                                          (hardening_mod; unchanged with hardening none)
  !>                                          [MPa]
  !>   mat_secant_s     (out) real(rk)(6, 6)  [s], sample frame [MPa s]
  !>   vec_h_s          (out) real(rk)(6)     {h}, sample frame [1/s]
  !>   mat_tangent_s    (out) real(rk)(6, 6)  [s_tan], sample frame [MPa s]
  !>   sliprate         (out) real(rk)(num_slip) slip rates (cutoff t_min) [1/s]
  !>   num_state_iters  (out, optional) integer  state-loop iterations performed (tests)
  !>   rss_final        (out, optional) real(rk)(num_slip) resolved shear stresses tau^a of
  !>                                          the final Kirchhoff stress (those of sliprate)
  !>                                          [MPa]
  !>   spin_hat_update_s  (out, optional) real(rk)(3)  the spin w^p of the last orientation
  !>                                          update (w23 w13 w12), sample frame of the final
  !>                                          orientation [1/s] (v2 d_rstar_spin; step 6.5b)
  !>   tau_start_c      (in, optional) real(rk)(6)  deviatoric start of the Kirchhoff stress,
  !>                                          crystal frame of the start orientation [MPa]
  !>                                          (initialization isotropic_viscoplastic,
  !>                                          increment 1: the given fraction of
  !>                                          cp_vertex_estimate; it replaces the
  !>                                          deviatoric part of the previous stress, v2
  !>                                          solve_evp_dev)
  !>   spin_slip_update_s (out, optional) real(rk)(3)  the slip spin sum gammadot q of the
  !>                                          last orientation update, the same form [1/s]
  !>                                          (v2 d_rstar_slip); the update's lattice spin is
  !>                                          spin_hat_update_s - spin_slip_update_s
  !>   tau_0_point      (in, optional) real(rk)(num_slip)  initial strengths of the point's
  !>                                          element (default: the phase's tau_0) [MPa]
  !>   km_rho            (inout, optional) real(rk)  Kocks-Mecking network dislocation density:
  !>                                          start of the increment in, end out (default:
  !>                                          the phase's km_rho_0 at the start, step 7.2)
  !>                                          [1/mm^2]
  !>   dmg_c            (inout, optional) real(rk)(6)  irradiation damage tensor H (Mandel,
  !>                                          crystal frame; step 7.3): start of the increment
  !>                                          in, end out; used when the model has
  !>                                          irradiation [1/mm^2]
  !>   status           (out) integer         cp_converged, cp_newton_failed (a stress solve
  !>                                          did not converge), cp_state_not_converged,
  !>                                          cp_singular (a singular Newton system), or
  !>                                          cp_hardening_failed (a backward-Euler strength
  !>                                          update did not converge: those systems keep
  !>                                          their start strength, v2)
  !>   work             (inout) type(cp_work_t)  work arrays (cp_work_alloc), with the point's
  !>                                          hardening values in work%hardening%point (the
  !>                                          phase's, no accumulated slip, after cp_work_alloc);
  !>                                          its tau_0_pt, rss, coef_sec, coef_tan and
  !>                                          crss_iter are this routine's: the initial strengths
  !>                                          used, the resolved shear stresses and the secant
  !>                                          and tangent coefficients of the current stress, the
  !>                                          strengths of the previous state iteration
  !> Locals:
  !>   tau_c            real(rk)(6)           Kirchhoff stress, crystal frame of rot_c2s [MPa]
  !>   strain_mat       real(rk)(3, 3)        e^e as a matrix [-]
  !>   beta             real(rk)              det(I + e^e) [-]
  !>   strain_el_c, strain_el0_c  real(rk)(6) e^e and e^e_0, crystal frame [-]
  !>   defrate_c        real(rk)(6)           {d}, crystal frame [1/s]
  !>   spin_c           real(rk)(3)           w, crystal frame [1/s]
  !>   spin_hat_c       real(rk)(3)           w^p, crystal frame [1/s]
  !>   spin_lat_c       real(rk)(3)           lattice spin w^p - sum gammadot q, crystal frame
  !>   spin_hat_upd_c, spin_slip_upd_c  real(rk)(3)  w^p and sum gammadot q of the last
  !>                                          orientation update, crystal frame [1/s]
  !>   compl_s, mstar_s, mtan_s  real(rk)(6, 6)  [S], [M*], [M*_tan], sample frame
  !>   amat             real(rk)(6, 6)        a compliance to invert
  !>   norm_prev, norm  real(rk)              stress norms of consecutive state iterations [MPa]
  !>   threshold        real(rk)              state-loop threshold [MPa]
  !>   newton_status, iter, ios  integer      Newton status, state iteration, solve status
  !>   newton_failed    logical               a stress solve of the state loop failed
  !>   hard_ok, hard_failed  logical          the last strength update converged; one failed
  !>   dp_c             real(rk)(6)           {d^p} of the final slip rates, crystal frame [1/s]
  !>   spin_mat         real(rk)(6, 6)        [W^] of the final w^p, sample frame [1/s]
  !>   kocks_mecking    logical               the law is Kocks-Mecking (the density is
  !>                                          integrated and the strengths follow from it)
  !>   rho_prev, rho_pt real(rk)              density at the start and current [1/mm^2]
  !>   irradiation      logical               irradiation softening (with Kocks-Mecking)
  !>   dmg_prev, dmg_pt real(rk)(6)           damage tensor H at the start and current,
  !>                                          Mandel, crystal frame [1/mm^2]
  !>   dmg_ok           logical               the last H update converged
  !>   damp            logical               damp the state between state iterations
  !>                                          (irradiation and model%damp_state; D62)
  !>   rho_iter, dmg_iter  real(rk), real(rk)(6)  the state x_k whose strengths the current
  !>                                          stress solve uses (damped iterate) [1/mm^2]
  !>   resid, resid_prev  real(rk)(7)         state residuals r_k = x~_k - x_k of this and
  !>                                          the previous state iteration, (rho, H Mandel)
  !>                                          [1/mm^2]
  !>   omega            real(rk)              damping factor of the last update [-]
  pure subroutine cp_point_solve(model, crss_prev, rot_c2s_prev, stress_s_prev, &
    & strain_el_s_prev, defrate_s, spin_s, dt, stress_s, strain_el_s, spin_lat_s, rot_c2s, crss, &
    & mat_secant_s, vec_h_s, mat_tangent_s, sliprate, status, num_state_iters, rss_final, &
    & spin_hat_update_s, spin_slip_update_s, tau_start_c, tau_0_point, km_rho, dmg_c, work)

    type(cp_model_t), intent(in) :: model
    real(rk), intent(in) :: crss_prev(:)
    real(rk), intent(in) :: rot_c2s_prev(3, 3)
    real(rk), intent(in) :: stress_s_prev(6)
    real(rk), intent(in) :: strain_el_s_prev(6)
    real(rk), intent(in) :: defrate_s(6)
    real(rk), intent(in) :: spin_s(3)
    real(rk), intent(in) :: dt
    real(rk), intent(out) :: stress_s(6)
    real(rk), intent(out) :: strain_el_s(6)
    real(rk), intent(out) :: spin_lat_s(3)
    real(rk), intent(out) :: rot_c2s(3, 3)
    real(rk), intent(out) :: crss(:)
    real(rk), intent(out) :: mat_secant_s(6, 6)
    real(rk), intent(out) :: vec_h_s(6)
    real(rk), intent(out) :: mat_tangent_s(6, 6)
    real(rk), intent(out) :: sliprate(:)
    integer, intent(out) :: status
    integer, intent(out), optional :: num_state_iters
    real(rk), intent(out), optional :: rss_final(:)
    real(rk), intent(out), optional :: spin_hat_update_s(3)
    real(rk), intent(out), optional :: spin_slip_update_s(3)
    real(rk), intent(in), optional :: tau_start_c(6)
    real(rk), intent(in), optional :: tau_0_point(:)
    real(rk), intent(inout), optional :: km_rho
    real(rk), intent(inout), optional :: dmg_c(6)
    type(cp_work_t), intent(inout) :: work

    logical :: kocks_mecking, irradiation, dmg_ok, damp
    real(rk) :: rho_iter, dmg_iter(6), resid(7), resid_prev(7), omega
    real(rk) :: rho_prev, rho_pt, dmg_prev(6), dmg_pt(6)
    real(rk) :: spin_hat_upd_c(3), spin_slip_upd_c(3)
    real(rk) :: tau_c(6), strain_mat(3, 3), beta, strain_el_c(6), strain_el0_c(6)
    real(rk) :: defrate_c(6), spin_c(3), spin_hat_c(3), spin_lat_c(3)
    real(rk) :: dp_c(6), spin_mat(6, 6)
    real(rk) :: compl_s(6, 6), mstar_s(6, 6), mtan_s(6, 6), amat(6, 6)
    real(rk) :: norm_prev, norm, threshold
    integer :: newton_status, iter, ios
    logical :: newton_failed, hard_ok, hard_failed

    status = cp_converged
    newton_failed = .false.
    hard_failed = .false.

    ! Initial strengths of the point (step 7.1; its hardening values are in
    ! work%hardening%point).
    work%tau_0_pt(:) = model%tau_0
    if (present(tau_0_point)) work%tau_0_pt(:) = tau_0_point
    ! Kocks-Mecking: the density at the start of the increment (step 7.2).
    kocks_mecking = model%hard%law == hardening_kocks_mecking
    rho_prev = model%hard%km_rho_0
    if (present(km_rho)) rho_prev = km_rho
    rho_pt = rho_prev
    ! Irradiation: the damage tensor at the start of the increment (step 7.3; zero if not
    ! given).
    irradiation = model%irradiation .and. kocks_mecking
    dmg_prev = 0.0_rk
    if (present(dmg_c)) dmg_prev = dmg_c
    dmg_pt = dmg_prev

    ! Start: the previous increment's Kirchhoff stress, in the crystal frame of the start
    ! orientation (v2 reuses its lattice-frame stress as it is).
    beta = det_identity_plus(matrix_from_mandel(strain_el_s_prev))
    tau_c = transform_sample_to_crystal_mandel(rot_c2s_prev, beta*stress_s_prev)
    ! Increment 1 with initialization isotropic_viscoplastic: the given deviatoric start,
    ! with the spherical part of the previous stress (v2: sig_vec = 0.5 x the viscoplastic estimate; sig_kk from
    ! its own volumetric update).
    if (present(tau_start_c)) tau_c = tau_c - mandel_dev(tau_c) + mandel_dev(tau_start_c)

    ! Forward-Euler estimates (v2 'euler_forward'), from the slip rates of the start stress and
    ! strengths: the strengths (hardening_mod; Kocks-Mecking: the density, then the strengths
    ! from it), then the orientation from the continuum spin and the slip spin, in the start
    ! crystal frame.
    call slip_rates_tmin(model, tau_c, crss_prev, sliprate, work%rss, work%coef_sec, &
      & work%coef_tan)
    if (irradiation) then
      call hardening_km_forward(model%hard, rho_prev, sliprate, dt, rho_pt)
      call hardening_irr_forward(model%irr, dmg_prev, sliprate, dt, dmg_pt)
      crss = hardening_irr_strengths(model%hard, model%irr, work%tau_0_pt, rho_pt, dmg_pt)
    else if (kocks_mecking) then
      call hardening_km_forward(model%hard, rho_prev, sliprate, dt, rho_pt)
      crss = hardening_km_strengths(model%hard, work%tau_0_pt, rho_pt)
    else
      call hardening_forward(model%hard, crss_prev, work%tau_0_pt, sliprate, dt, crss, &
        & work%hardening)
    end if
    spin_hat_upd_c = transform_sample_to_crystal_skew(rot_c2s_prev, spin_s)
    spin_slip_upd_c = matmul(model%schmid_skew_c, sliprate)
    spin_lat_c = spin_hat_upd_c - spin_slip_upd_c
    rot_c2s = rot_c2s_spin_update_crystal(rot_c2s_prev, spin_lat_c, dt)

    ! State-loop damping (irradiation only, D62): the first iterate x_1 is the forward
    ! estimate, the factor starts at 1 for every point solve.
    damp = irradiation .and. model%damp_state
    rho_iter = rho_pt
    dmg_iter = dmg_pt
    resid_prev = 0.0_rk
    omega = 1.0_rk

    ! State loop: stress solve in the current frame, then the orientation from the new spin.
    threshold = model%tol_state*work%tau_0_pt(1)
    norm_prev = mandel_norm(tau_c)
    spin_hat_c = 0.0_rk
    status = cp_state_not_converged
    newton_status = cp_converged
    do iter = 1, model%max_iters_state
      ! Inputs in the crystal frame of the current orientation; e^e_0 held in the sample frame
      ! (recoupled_formulation.md §4).
      defrate_c = transform_sample_to_crystal_mandel(rot_c2s, defrate_s)
      spin_c = transform_sample_to_crystal_skew(rot_c2s, spin_s)
      strain_el0_c = transform_sample_to_crystal_mandel(rot_c2s, strain_el_s_prev)

      ! Stress solve (Newton); a failure keeps the last iterate (D49).
      call newton_stress(model, crss, defrate_c, spin_c, strain_el0_c, dt, tau_c, spin_hat_c, &
        & newton_status, work)
      if (newton_status == cp_singular) exit
      if (newton_status == cp_newton_failed) newton_failed = .true.

      ! Backward-Euler updates (v2 'euler_backward', solve_evp_rstar order), from the slip
      ! rates of the new stress with the cutoff t_min: the strengths (Kocks-Mecking: the
      ! density, then the strengths from it), then the orientation from w^p of the stress solve
      ! minus the slip spin.
      call slip_rates_tmin(model, tau_c, crss, sliprate, work%rss, work%coef_sec, &
        & work%coef_tan)
      work%crss_iter(:) = crss
      if (irradiation) then
        ! With damping the backward Euler starts from the damped iterate x_k (without, x_k
        ! is the previous update already).
        if (damp) then
          rho_pt = rho_iter
          dmg_pt = dmg_iter
        end if
        call hardening_km_backward(model%hard, rho_prev, sliprate, dt, rho_pt, hard_ok, &
          & work%hardening)
        call hardening_irr_backward(model%hard, model%irr, dmg_prev, sliprate, dt, dmg_pt, &
          & dmg_ok, work%hardening)
        hard_ok = hard_ok .and. dmg_ok
        crss = hardening_irr_strengths(model%hard, model%irr, work%tau_0_pt, rho_pt, dmg_pt)
      else if (kocks_mecking) then
        call hardening_km_backward(model%hard, rho_prev, sliprate, dt, rho_pt, hard_ok, &
          & work%hardening)
        crss = hardening_km_strengths(model%hard, work%tau_0_pt, rho_pt)
      else
        call hardening_backward(model%hard, crss_prev, work%tau_0_pt, sliprate, dt, crss, hard_ok, &
          & work%hardening)
      end if
      if (.not. hard_ok) hard_failed = .true.
      spin_hat_upd_c = spin_hat_c
      spin_slip_upd_c = matmul(model%schmid_skew_c, sliprate)
      spin_lat_c = spin_hat_upd_c - spin_slip_upd_c
      rot_c2s = rot_c2s_spin_update_crystal(rot_c2s_prev, spin_lat_c, dt)

      ! Convergence of the stress norm and of every strength (v2: strict "<", threshold
      ! tol_state tau_0(1); v2 accepts when one system's strength has settled, here all
      ! must, the same for the isotropic law, whose systems change alike).
      ! With damping crss are the undamped strengths crss(x~_k) and crss_iter crss(x_k): the
      ! test is on the undamped residual (D62 §5).
      norm = mandel_norm(tau_c)
      if (abs(norm - norm_prev) < threshold .and. all(abs(crss - work%crss_iter) < threshold)) then
        status = cp_converged
        exit
      end if
      norm_prev = norm

      ! Damped state for the next iteration (D62 §3; not after the last iteration, whose
      ! undamped update is kept): the residual r_k, the factor (1 at the first iteration,
      ! then Aitken), x_{k+1} = x_k + omega r_k, and its strengths.
      if (damp .and. iter < model%max_iters_state) then
        resid(1) = rho_pt - rho_iter
        resid(2:7) = dmg_pt - dmg_iter
        if (iter > 1) omega = cp_state_damping_factor(resid_prev, resid, omega)
        rho_iter = rho_iter + omega*resid(1)
        dmg_iter = dmg_iter + omega*resid(2:7)
        resid_prev = resid
        crss = hardening_irr_strengths(model%hard, model%irr, work%tau_0_pt, rho_iter, dmg_iter)
      end if
    end do
    if (present(num_state_iters)) num_state_iters = min(iter, model%max_iters_state)
    if (present(km_rho)) km_rho = rho_pt
    if (present(dmg_c) .and. irradiation) dmg_c = dmg_pt
    if (hard_failed) status = cp_hardening_failed
    if (newton_failed) status = cp_newton_failed
    if (newton_status == cp_singular) status = cp_singular

    ! Final state: elastic strain, stress and w^p in the sample frame of the final orientation
    ! (v2 material_matrix_evp_aniso, find_wp_hat_phase).
    strain_el_c = matmul(model%compl_el_c, tau_c)
    strain_mat = matrix_from_mandel(strain_el_c)
    beta = det_identity_plus(strain_mat)
    strain_el_s = transform_crystal_to_sample_mandel(rot_c2s, strain_el_c)
    stress_s = transform_crystal_to_sample_mandel(rot_c2s, tau_c)/beta
    call slip_rates_tmin(model, tau_c, crss, sliprate, work%rss, work%coef_sec, &
      & work%coef_tan)
    if (present(rss_final)) rss_final = work%rss
    strain_el0_c = transform_sample_to_crystal_mandel(rot_c2s, strain_el_s_prev)
    spin_c = transform_sample_to_crystal_skew(rot_c2s, spin_s)
    dp_c = matmul(model%schmid_sym_c, sliprate)
    spin_hat_c = plastic_spin_hat(spin_c, strain_el_c, strain_el0_c, dp_c, dt)
    spin_lat_s = transform_crystal_to_sample_skew(rot_c2s, spin_hat_c)

    ! The two parts of the last orientation update, sample frame (outputs rotrate_spin,
    ! rotrate_slip: the rotation vector of exp(w dt) in the crystal frame is w dt, and its
    ! sample-frame form is the same with the start or the final orientation).
    if (present(spin_hat_update_s)) spin_hat_update_s = &
      & transform_crystal_to_sample_skew(rot_c2s, spin_hat_upd_c)
    if (present(spin_slip_update_s)) spin_slip_update_s = &
      & transform_crystal_to_sample_skew(rot_c2s, spin_slip_upd_c)

    ! [M*], [M*_tan] and [S] in the sample frame; [s] = ((beta/dt)[S] + beta[M*])^-1 (R11) and
    ! its tangent form; {h} (R12).
    compl_s = transform_crystal_to_sample_mat6(rot_c2s, model%compl_el_c)
    mstar_s = transform_crystal_to_sample_mat6(rot_c2s, &
      & plastic_compliance(model%schmid_sym_c, work%coef_sec))
    mtan_s = transform_crystal_to_sample_mat6(rot_c2s, &
      & plastic_compliance(model%schmid_sym_c, work%coef_tan))
    amat = beta*(compl_s/dt + mstar_s)
    call invert_spd(amat, mat_secant_s, ios)
    if (ios /= 0) status = cp_singular
    amat = beta*(compl_s/dt + mtan_s)
    call invert_spd(amat, mat_tangent_s, ios)
    if (ios /= 0) status = cp_singular
    spin_mat = spin_commutator(spin_lat_s)
    vec_h_s = matmul(spin_mat, strain_el_s) - strain_el_s_prev/dt

  end subroutine cp_point_solve

  !> Newton solve of the point equation for the Kirchhoff stress at fixed orientation and
  !> strengths (v2 solve_stress_evp_newton; module header).
  !>
  !> Arguments:
  !>   model         (in)    type(cp_model_t)    model
  !>   crss          (in)    real(rk)(:)         strengths [MPa]
  !>   defrate_c     (in)    real(rk)(6)         {d}, crystal frame [1/s]
  !>   spin_c        (in)    real(rk)(3)         w, crystal frame [1/s]
  !>   strain_el0_c  (in)    real(rk)(6)         e^e_0, crystal frame [-]
  !>   dt            (in)    real(rk)            time increment [s]
  !>   tau_c         (inout) real(rk)(6)         Kirchhoff stress: start value in, solution (or
  !>                                             last iterate) out [MPa]
  !>   spin_hat_c    (out)   real(rk)(3)         w^p at the returned stress [1/s]
  !>   status        (out)   integer             cp_converged, cp_newton_failed or cp_singular
  !>   work          (inout) type(cp_work_t)     work arrays (cp_work_alloc): newton_coef_tan,
  !>                                             the tangent coefficients at tau_0 [1/(MPa s)];
  !>                                             residual_rss, residual_sliprate,
  !>                                             residual_coef_sec, cp_residual's
  !> Locals:
  !>   tau_0         real(rk)(6)      stress at the start of an iteration [MPa]
  !>   delta         real(rk)(6)      Newton step [MPa]
  !>   rhs, rhs_try  real(rk)(6)      residual vectors [1/s]
  !>   residual_norm_start  real(rk)    residual norm at tau_0 [1/s]
  !>   residual_norm        real(rk)    residual norm at the trial stress [1/s]
  !>   jac           real(rk)(6, 6)   Jacobian [1/(MPa s)]
  !>   jac_tj        real(rk)(6, 6)   J^T J of the normal equations, then its LDL^T factors
  !>   spin_mat      real(rk)(6, 6)   [W^] of w^p at tau_0 [1/s]
  !>   spin_try      real(rk)(3)      w^p at the trial stress [1/s]
  !>   fact          real(rk)         line-search factor
  !>   iter, i, ios  integer          iteration, index, solve status
  pure subroutine newton_stress(model, crss, defrate_c, spin_c, strain_el0_c, dt, tau_c, &
    & spin_hat_c, status, work)

    type(cp_model_t), intent(in) :: model
    real(rk), intent(in) :: crss(:)
    real(rk), intent(in) :: defrate_c(6)
    real(rk), intent(in) :: spin_c(3)
    real(rk), intent(in) :: strain_el0_c(6)
    real(rk), intent(in) :: dt
    real(rk), intent(inout) :: tau_c(6)
    real(rk), intent(out) :: spin_hat_c(3)
    integer, intent(out) :: status
    type(cp_work_t), intent(inout) :: work

    real(rk) :: tau_0(6), delta(6), rhs(6), rhs_try(6), residual_norm_start, residual_norm
    real(rk) :: jac(6, 6), jac_tj(6, 6), spin_mat(6, 6), spin_try(3)
    real(rk) :: fact
    integer :: iter, i, ios

    status = cp_newton_failed
    spin_hat_c = 0.0_rk
    do iter = 1, model%max_iters_newton
      ! Residual and w^p at the current stress; Jacobian with w^p frozen (v2 matrix_fjac).
      tau_0 = tau_c
      call cp_residual(model, tau_0, crss, defrate_c, spin_c, strain_el0_c, dt, rhs, &
        & residual_norm_start, spin_hat_c, work%newton_coef_tan, work%residual_rss, &
        & work%residual_sliprate, work%residual_coef_sec)
      spin_mat = spin_commutator(spin_hat_c)
      jac = model%compl_el_c/dt + matmul(spin_mat, model%compl_el_c) + &
        & plastic_compliance(model%schmid_sym_c, work%newton_coef_tan)

      ! Newton step from the normal equations, with v2's check of the diagonal of J^T J.
      do i = 1, 6
        if (dot_product(jac(:, i), jac(:, i)) < diag_min) then
          status = cp_singular
          return
        end if
      end do
      call normal_equations_solve(jac, rhs, delta, jac_tj, ios)
      if (ios /= 0) then
        status = cp_singular
        return
      end if

      ! Line search: halve the step while the residual norm grows; below a factor 0.001 the
      ! solve fails and the stress stays at the iteration's start value (v2), unless that start
      ! already meets the tolerance (deviation, D49: a stress that solves the equations to
      ! round-off, e.g. re-entered from the state loop, cannot be improved, and v2 would report
      ! a failure). A residual that is not a number (overflow of the power law) counts as
      ! growing (v2 would accept it).
      fact = 1.0_rk
      call cp_residual(model, tau_0 + delta, crss, defrate_c, spin_c, strain_el0_c, dt, &
        & rhs_try, residual_norm, spin_try, work%newton_coef_tan, work%residual_rss, &
        & work%residual_sliprate, work%residual_coef_sec)
      do while (.not. (residual_norm <= residual_norm_start))
        fact = 0.5_rk*fact
        if (fact < fact_min) then
          tau_c = tau_0
          if (.not. (residual_norm_start > model%tol_newton)) status = cp_converged
          return
        end if
        call cp_residual(model, tau_0 + fact*delta, crss, defrate_c, spin_c, strain_el0_c, &
          & dt, rhs_try, residual_norm, spin_try, work%newton_coef_tan, work%residual_rss, &
          & work%residual_sliprate, work%residual_coef_sec)
      end do

      ! Accept the step; converged when its residual norm is at most point_tol.
      tau_c = tau_0 + fact*delta
      spin_hat_c = spin_try
      if (.not. (residual_norm > model%tol_newton)) then
        status = cp_converged
        return
      end if
    end do

  end subroutine newton_stress

  !> Residual of the point equation (R1) at a Kirchhoff stress, with the slip rates of the
  !> Newton cutoff and w^p by (R2) (v2 residual, dp_wp_hat_phase):
  !> rhs = {d} - ({e^e} - {e^e_0})/dt - {d^p} - [W^]{e^e}, {e^e} = [S]{tau}.
  !>
  !> Theory : recoupled_formulation.md (R1), (R2), (R5)-(R9).
  !> Source : port of v2 solve_evp_stress_mod3.f90::residual and kinematics_mod.f90::
  !>          dp_wp_hat_phase.
  !>
  !> Arguments:
  !>   model         (in)  type(cp_model_t)  model
  !>   tau_c         (in)  real(rk)(6)       Kirchhoff stress, crystal frame [MPa]
  !>   crss          (in)  real(rk)(:)       strengths [MPa]
  !>   defrate_c     (in)  real(rk)(6)       {d}, crystal frame [1/s]
  !>   spin_c        (in)  real(rk)(3)       w, crystal frame [1/s]
  !>   strain_el0_c  (in)  real(rk)(6)       e^e_0, crystal frame [-]
  !>   dt            (in)  real(rk)          time increment [s]
  !>   rhs           (out) real(rk)(6)       residual vector [1/s]
  !>   residual_norm (out) real(rk)          its Euclidean norm [1/s]
  !>   spin_hat_c    (out) real(rk)(3)       w^p, crystal frame [1/s]
  !>   coef_tan      (out) real(rk)(:)       tangent coefficients d gammadot/d tau [1/(MPa s)]
  !>   rss, sliprate, coef_sec  (out) real(rk)(:)  the caller's work arrays (ROADMAP D100):
  !>                                         resolved shear stresses [MPa], slip rates (Newton
  !>                                         cutoff) [1/s], secant coefficients [1/(MPa s)]
  !> Locals:
  !>   strain_el_c   real(rk)(6)             e^e [-]
  !>   dp_c          real(rk)(6)             {d^p} [1/s]
  !>   spin_mat      real(rk)(6, 6)          [W^] of w^p [1/s]
  pure subroutine cp_residual(model, tau_c, crss, defrate_c, spin_c, strain_el0_c, dt, rhs, &
    & residual_norm, spin_hat_c, coef_tan, rss, sliprate, coef_sec)

    type(cp_model_t), intent(in) :: model
    real(rk), intent(in) :: tau_c(6)
    real(rk), intent(in) :: crss(:)
    real(rk), intent(in) :: defrate_c(6)
    real(rk), intent(in) :: spin_c(3)
    real(rk), intent(in) :: strain_el0_c(6)
    real(rk), intent(in) :: dt
    real(rk), intent(out) :: rhs(6)
    real(rk), intent(out) :: residual_norm
    real(rk), intent(out) :: spin_hat_c(3)
    real(rk), intent(out) :: coef_tan(:)
    real(rk), intent(out) :: rss(:)
    real(rk), intent(out) :: sliprate(:)
    real(rk), intent(out) :: coef_sec(:)

    real(rk) :: strain_el_c(6), dp_c(6), spin_mat(6, 6)

    ! Elastic strain and slip (Newton cutoff tiny**m, v2 solve_stress_evp_newton).
    strain_el_c = matmul(model%compl_el_c, tau_c)
    rss = resolved_shear_stress(model%schmid_sym_c, tau_c)
    call slip_power_law(rss, crss, model%m, model%gammadot_0, model%rss_cutoff_newton, &
      & sliprate, coef_sec, coef_tan)
    dp_c = matmul(model%schmid_sym_c, sliprate)

    ! w^p (R2) and the residual (R1).
    spin_hat_c = plastic_spin_hat(spin_c, strain_el_c, strain_el0_c, dp_c, dt)
    spin_mat = spin_commutator(spin_hat_c)
    rhs = defrate_c - (strain_el_c - strain_el0_c)/dt - dp_c - matmul(spin_mat, strain_el_c)
    residual_norm = norm2(rhs)

  end subroutine cp_residual

  !> w^p by (R2) with the backward-difference elastic strain rate:
  !> w^p = w - (e d^p - d^p e) + (e e_0 - e_0 e)/(2 dt).
  !>
  !> Theory : recoupled_formulation.md (R2) and its note (v2 find_wp_hat).
  !>
  !> Arguments:
  !>   spin        (in)  real(rk)(3)  w (w23 w13 w12) [1/s]
  !>   strain_el   (in)  real(rk)(6)  e^e [-]
  !>   strain_el0  (in)  real(rk)(6)  e^e_0 [-]
  !>   dp          (in)  real(rk)(6)  d^p [1/s]
  !>   dt          (in)  real(rk)     time increment [s]
  !>   spin_hat    (out, result) real(rk)(3)  w^p [1/s]
  !>   (all in one frame)
  !> Locals:
  !>   e, e0, p    real(rk)(3, 3)  the tensors as matrices
  pure function plastic_spin_hat(spin, strain_el, strain_el0, dp, dt) result(spin_hat)

    real(rk), intent(in) :: spin(3)
    real(rk), intent(in) :: strain_el(6)
    real(rk), intent(in) :: strain_el0(6)
    real(rk), intent(in) :: dp(6)
    real(rk), intent(in) :: dt
    real(rk) :: spin_hat(3)

    real(rk) :: e(3, 3), e0(3, 3), p(3, 3)

    e = matrix_from_mandel(strain_el)
    e0 = matrix_from_mandel(strain_el0)
    p = matrix_from_mandel(dp)
    spin_hat = spin - skew_from_matrix(matmul(e, p) - matmul(p, e)) + &
      & skew_from_matrix((matmul(e, e0) - matmul(e0, e))/(2.0_rk*dt))

  end function plastic_spin_hat

  !> v2's increment-1 stress estimate (initialization isotropic_viscoplastic): the stress of
  !> the rigid-viscoplastic crystal (power law with m = 0.4) for the deviatoric deformation rate scaled to
  !> unit equivalent rate, solved by Newton's method from the yield-surface vertex that does the
  !> most plastic work on that rate, then scaled back. The caller starts the point solve from
  !> 0.5 times the result (vertex_start_fraction).
  !>
  !> Theory : with d' the deviatoric deformation rate in the crystal frame and d_eq =
  !>          sqrt(2/3 d':d'), solve for a deviatoric tau the 5 equations (deviatoric basis,
  !>          addendum_5vectors.md)
  !>            r(tau) = sum_a gammadot^a(tau) p^a - d'/d_eq = 0,
  !>            gammadot^a = gammadot_0^a |t^a|^(1/0.4 - 1) t^a,  t^a = p^a:tau / crss_1
  !>          (zero for |t^a| <= t_min, the cutoff of the phase's own m; crss_1 the strength of
  !>          system 1) by Newton's method, J = sum_a (dgammadot^a/dtau^a) p^a (x) p^a (the
  !>          tangent compliance with m = 0.4, slip_kinetics_mod), with v2's line search
  !>          (halve the step while |r| grows; fail below a factor fact_min) and the stop
  !>          |r| <= point_tol after a step, at most point_max_iters_newton iterations. The start
  !>          is the vertex sigma_v sgn(w_v) with the largest |w_v|, w_v = sigma_v : d'/d_eq
  !>          (first on ties), the vertices being stresses computed from the phase's systems and
  !>          tau_0 (D57; v2's tables, in units of crss_1, are not used); a vertex is not used
  !>          twice; up to 4
  !>          vertices are tried; if none converges, the last vertex stress is kept. The result
  !>          is scaled by d_eq^0.4.
  !> Source : port of v2-fixed solve_evp_vpstress_mod1.f90::solve_evp_vpstress and
  !>          solve_evp_vpstress_mod2/3.f90 (scale_down_defr, compute_work, compute_avg_crss,
  !>          find_vertex, vertex_stress, solve_newton_vp, get_res, form_crystif, power_law,
  !>          compliance, solvit, check_diagonals, scale_up_sigm), with CHANGES.md C43 (the
  !>          residual norm of all five components) and C44 (vertex per phase); the vertices
  !>          computed by yield_vertices_compute (addendum_yield_vertices.md, D57) in place of v2's
  !>          tables (removed in step 9.4, D56). v2's "average" strength is that
  !>          of system 1 (compute_avg_crss), and its work uses the initial strength of system 1
  !>          (equal to the current one in increment 1); its power law uses the phase's
  !>          gammadot_0 (here per system, equal when uniform) and the cutoff t_min of the real m.
  !>          Deviations: a zero deformation rate gives a zero estimate (v2 divides 0 by 0 in
  !>          find_vertex); the Newton system is solved with J (positive semidefinite) and the
  !>          step negated (v2: -J and the step), and a nonpositive pivot of the LDL^T
  !>          factorization ends the trial (v2 detects only a zero diagonal, check_diagonals,
  !>          which is also tested here).
  !>
  !> Arguments:
  !>   model      (in)  type(cp_model_t)    model (vertices_c set)
  !>   crss       (in)  real(rk)(num_slip)  slip strengths [MPa]
  !>   defrate_c  (in)  real(rk)(6)         {d}, Mandel, crystal frame [1/s]
  !>   tau_dev_c  (out) real(rk)(6)         the estimate, deviatoric, Mandel, crystal frame [MPa]
  !>   status     (out) integer             0: a trial converged; 1: none did (the last vertex
  !>                                        stress, as v2); 2: no vertex table
  !>   work       (inout) type(cp_work_t)   work arrays (cp_work_alloc): vertex_work, the virtual
  !>                                        plastic work of each vertex not yet used; vertex_m
  !>                                        and vertex_crss, m = 0.4 and crss_1 for every system;
  !>                                        vertex_rss, vertex_rate, vertex_coef_s, vertex_coef_t,
  !>                                        vertex_residual's
  !> Locals:
  !>   deq          real(rk)                 equivalent deformation rate d_eq [1/s]
  !>   d_unit       real(rk)(6)              d'/d_eq [-]
  !>   crss_1       real(rk)                 strength of system 1 [MPa]
  !>   best         real(rk)                 largest |work|
  !>   sign_v       real(rk)                 sign of the chosen vertex's work
  !>   tau_vertex   real(rk)(6)              the trial's vertex stress [MPa]
  !>   tau, tau_iter, tau_try  real(rk)(6)   iterate, its value at the iteration start, a
  !>                                         line-search point [MPa]
  !>   step         real(rk)(6)              Newton step [MPa]
  !>   res_iter, res  real(rk)               residual norms at tau_iter and at the new point
  !>   fact         real(rk)                 line-search factor
  !>   r5           real(rk)(5)              residual, deviatoric 5-vector
  !>   jac5         real(rk)(5, 5)           Jacobian, deviatoric 5-vector basis (then its LDL^T
  !>                                         factors)
  !>   x5           real(rk)(5)              solution of jac5 x5 = r5
  !>   i_trial, iter, i, v, ios  integer     trial, Newton iteration, index, vertex, status
  !>   newton_ok, converged  logical         the trial's Newton is usable; it converged
  pure subroutine cp_vertex_estimate(model, crss, defrate_c, tau_dev_c, status, work)

    type(cp_model_t), intent(in) :: model
    real(rk), intent(in) :: crss(:)
    real(rk), intent(in) :: defrate_c(6)
    real(rk), intent(out) :: tau_dev_c(6)
    integer, intent(out) :: status
    type(cp_work_t), intent(inout) :: work

    real(rk) :: deq, d_unit(6), crss_1, best, sign_v
    real(rk) :: tau_vertex(6), tau(6), tau_iter(6), tau_try(6), step(6), res_iter, res, fact
    real(rk) :: r5(5), jac5(5, 5), x5(5)
    integer :: i_trial, iter, i, v, ios
    logical :: newton_ok, converged

    tau_dev_c = 0.0_rk
    status = 2
    if (size(model%vertices_c, 2) == 0) return

    ! Deviatoric rate scaled to unit equivalent rate (v2 scale_down_defr); a zero rate gives a
    ! zero estimate.
    status = 0
    deq = equiv_strain(defrate_c)
    if (.not. (deq > 0.0_rk)) return
    d_unit = mandel_dev(defrate_c)/deq

    ! Strength of system 1 for every system, m = 0.4, and the virtual plastic work of each
    ! vertex (v2 compute_avg_crss, compute_work); the vertices are stresses (D57).
    crss_1 = crss(1)
    work%vertex_crss(:) = crss_1
    work%vertex_m(:) = vertex_m
    work%vertex_work(:) = matmul(d_unit, model%vertices_c)

    ! Trials from the vertices of largest work (v2 solve_evp_vpstress).
    converged = .false.
    tau_vertex = 0.0_rk
    do i_trial = 1, vertex_num_trials
      ! The vertex of largest |work| (first on ties), its sign; not used again (find_vertex,
      ! vertex_stress).
      best = -1.0_rk
      v = 1
      do i = 1, size(work%vertex_work)
        if (abs(work%vertex_work(i)) > best) then
          best = abs(work%vertex_work(i))
          v = i
        end if
      end do
      sign_v = work%vertex_work(v)/best
      work%vertex_work(v) = 0.0_rk
      tau_vertex = sign_v*model%vertices_c(:, v)
      tau = tau_vertex

      ! Newton's method with v2's line search (solve_newton_vp).
      newton_ok = .true.
      do iter = 1, model%max_iters_newton
        ! Residual and Jacobian at the iteration start; the step from J x = r.
        tau_iter = tau
        call vertex_residual(tau_iter, res_iter, r5, jac5, work%vertex_rss, work%vertex_rate, &
          & work%vertex_coef_s, work%vertex_coef_t)
        do i = 1, 5
          if (abs(jac5(i, i)) < diag_min) newton_ok = .false.
        end do
        if (.not. newton_ok) exit
        call ldlt_factor(jac5, ios)
        if (ios /= 0) newton_ok = .false.
        if (.not. newton_ok) exit
        x5 = r5
        call ldlt_substitute(jac5, x5)
        step = -mandel_from_dev5(x5)

        ! Line search: halve the step while the residual norm grows.
        fact = 1.0_rk
        tau_try = tau_iter + step
        call vertex_residual(tau_try, res, r5, jac5, work%vertex_rss, work%vertex_rate, &
          & work%vertex_coef_s, work%vertex_coef_t)
        do while (res > res_iter)
          fact = 0.5_rk*fact
          if (fact < fact_min) then
            newton_ok = .false.
            exit
          end if
          tau_try = tau_iter + fact*step
          call vertex_residual(tau_try, res, r5, jac5, work%vertex_rss, work%vertex_rate, &
            & work%vertex_coef_s, work%vertex_coef_t)
        end do
        if (.not. newton_ok) exit
        tau = tau_try

        ! Converged when the residual after the step is at most point_tol.
        if (res <= model%tol_newton) then
          converged = .true.
          exit
        end if
      end do
      if (converged) exit
    end do

    ! No trial converged: the last vertex stress (v2); then scaled back by d_eq^0.4
    ! (scale_up_sigm).
    if (.not. converged) then
      tau = tau_vertex
      status = 1
    end if
    tau_dev_c = tau*deq**vertex_m

  contains

    !> The residual of the estimate at a stress, its norm and its Jacobian (v2 get_res,
    !> form_crystif), in the deviatoric 5-vector basis.
    !>
    !> Arguments:
    !>   tau_at  (in)  real(rk)(6)     stress [MPa]
    !>   norm    (out) real(rk)        |r| [1/s]
    !>   res5    (out) real(rk)(5)     r [1/s]
    !>   j5      (out) real(rk)(5, 5)  J [1/(MPa s)]
    !>   rss, rate   (out) real(rk)(num_slip)  work arrays (the host's work%vertex_rss,
    !>                                 vertex_rate): resolved shear stresses [MPa], slip rates
    !>                                 [1/s]
    !>   coef_s, coef_t  (out) real(rk)(num_slip)  work arrays (work%vertex_coef_s,
    !>                                 vertex_coef_t): secant and tangent coefficients
    !>                                 [1/(MPa s)]
    !> Locals:
    !>   dp          real(rk)(6)  sum gammadot p [1/s]
    pure subroutine vertex_residual(tau_at, norm, res5, j5, rss, rate, coef_s, coef_t)

      real(rk), intent(in) :: tau_at(6)
      real(rk), intent(out) :: norm
      real(rk), intent(out) :: res5(5)
      real(rk), intent(out) :: j5(5, 5)
      real(rk), intent(out) :: rss(:)
      real(rk), intent(out) :: rate(:)
      real(rk), intent(out) :: coef_s(:)
      real(rk), intent(out) :: coef_t(:)

      real(rk) :: dp(6)

      ! Slip rates of the m = 0.4 power law with the strength crss_1 and the phase's cutoff.
      rss = resolved_shear_stress(model%schmid_sym_c, tau_at)
      call slip_power_law(rss, work%vertex_crss, work%vertex_m, model%gammadot_0, &
        & model%rss_cutoff, rate, coef_s, coef_t)

      ! r = sum gammadot p - d'/d_eq (all five deviatoric components, C43); J = [M*_tan].
      dp = matmul(model%schmid_sym_c, rate)
      res5 = dev5_from_mandel(dp - d_unit)
      norm = norm2(res5)
      j5 = mat5_from_mat6(plastic_compliance(model%schmid_sym_c, coef_t))

    end subroutine vertex_residual

  end subroutine cp_vertex_estimate

  !> Slip rates with the cutoff t_min, set to zero at or below 1e-30 (v2 compute_sliprate), and
  !> the secant and tangent coefficients with the same cutoff.
  !>
  !> Arguments:
  !>   model     (in)  type(cp_model_t)  model
  !>   tau_c     (in)  real(rk)(6)       Kirchhoff stress, crystal frame [MPa]
  !>   crss      (in)  real(rk)(:)       strengths [MPa]
  !>   sliprate  (out) real(rk)(:)       slip rates [1/s]
  !>   rss       (out) real(rk)(:)       resolved shear stresses [MPa]
  !>   coef_sec, coef_tan  (out) real(rk)(:)  secant and tangent coefficients [1/(MPa s)]
  pure subroutine slip_rates_tmin(model, tau_c, crss, sliprate, rss, coef_sec, coef_tan)

    type(cp_model_t), intent(in) :: model
    real(rk), intent(in) :: tau_c(6)
    real(rk), intent(in) :: crss(:)
    real(rk), intent(out) :: sliprate(:)
    real(rk), intent(out) :: rss(:)
    real(rk), intent(out) :: coef_sec(:)
    real(rk), intent(out) :: coef_tan(:)

    rss = resolved_shear_stress(model%schmid_sym_c, tau_c)
    call slip_power_law(rss, crss, model%m, model%gammadot_0, model%rss_cutoff, sliprate, &
      & coef_sec, coef_tan)
    where (.not. (abs(sliprate) > sliprate_zero)) sliprate = 0.0_rk

  end subroutine slip_rates_tmin

  !> Inverse of a symmetric positive definite 6 x 6 matrix: one LDL^T factorization
  !> (ldlt_factor), then column by column (ldlt_substitute); zero if a is not positive definite.
  !>
  !> Arguments:
  !>   a       (in)  real(rk)(6, 6)  matrix
  !>   a_inv   (out) real(rk)(6, 6)  inverse
  !>   status  (out) integer         0, or nonzero if a is not positive definite
  !> Locals:
  !>   i       integer         column
  !>   fact    real(rk)(6, 6)  LDL^T factors of a
  pure subroutine invert_spd(a, a_inv, status)

    real(rk), intent(in) :: a(6, 6)
    real(rk), intent(out) :: a_inv(6, 6)
    integer, intent(out) :: status

    integer :: i
    real(rk) :: fact(6, 6)

    ! Factorization; a failure gives a zero inverse.
    fact = a
    call ldlt_factor(fact, status)
    if (status /= 0) then
      a_inv = 0.0_rk
      return
    end if

    ! The columns of the identity, solved in place.
    do i = 1, 6
      a_inv(:, i) = 0.0_rk
      a_inv(i, i) = 1.0_rk
      call ldlt_substitute(fact, a_inv(:, i))
    end do

  end subroutine invert_spd

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

  !> Bounded Aitken damping factor of the state loop (step 7.3, D62,
  !> dev/theory/addendum_state_damping.md §3-§4).
  !>
  !> Theory : addendum_state_damping.md §3 (revision 3): omega^A = -omega_prev r_prev.(r - r_prev) /
  !>          |r - r_prev|^2; omega = 1 if omega^A <= 0 (lambda >= 1: no overshoot, damping cannot
  !>          help), else min(1, max(0.001, omega^A)); the previous factor when
  !>          |r - r_prev| <= 1e-14 |r|. In the linear scalar case omega^A = 1/(1 - lambda),
  !>          lambda the slope of the undamped map.
  !> Source : Kuttler & Wall (2008), Comput. Mech. 43:61-72, Eq. 44 (vector Aitken factor
  !>          after Irons & Tuck 1969, with the recursion on omega, Remark 8).
  !>
  !> Arguments:
  !>   resid_prev  (in) real(rk)(:)  residual r_{k-1} of the previous state iteration
  !>   resid       (in) real(rk)(:)  residual r_k of this iteration (same size and units)
  !>   omega_prev  (in) real(rk)     factor applied to r_{k-1} [-]
  !> Result:
  !>   omega       real(rk)          factor to apply to r_k, in [0.001, 1] [-]
  !> Locals (explicit loops in the order of the dot products: no array temporary, D100):
  !>   i           integer    component
  !>   change      real(rk)   r_k,i - r_{k-1},i
  !>   change_sq   real(rk)   |r_k - r_{k-1}|^2
  !>   resid_sq    real(rk)   |r_k|^2
  !>   prev_change real(rk)   r_{k-1} . (r_k - r_{k-1})
  pure function cp_state_damping_factor(resid_prev, resid, omega_prev) result(omega)

    real(rk), intent(in) :: resid_prev(:)
    real(rk), intent(in) :: resid(:)
    real(rk), intent(in) :: omega_prev
    real(rk) :: omega

    integer :: i
    real(rk) :: change, change_sq, resid_sq, prev_change

    ! Change of the residual between the two iterations: its square, |r_k|^2 and
    ! r_{k-1} . (r_k - r_{k-1}).
    change_sq = 0.0_rk
    resid_sq = 0.0_rk
    prev_change = 0.0_rk
    do i = 1, size(resid)
      change = resid(i) - resid_prev(i)
      change_sq = change_sq + change*change
      resid_sq = resid_sq + resid(i)*resid(i)
      prev_change = prev_change + resid_prev(i)*change
    end do

    ! No information (the residual did not change): keep the previous factor.
    if (change_sq <= (damp_change_min**2)*resid_sq) then
      omega = omega_prev
      return
    end if

    ! Aitken factor (Eq. 44); a non-positive estimate (lambda >= 1, no overshoot) gives 1, any
    ! other is bounded to [omega_min, 1].
    omega = -omega_prev*prev_change/change_sq
    if (omega <= 0.0_rk) then
      omega = 1.0_rk
    else
      omega = min(1.0_rk, max(damp_omega_min, omega))
    end if

  end function cp_state_damping_factor

end module crystal_plasticity_model_mod
