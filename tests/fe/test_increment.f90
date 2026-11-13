! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> MPI test of increment_mod and face_force_mod (1, 2 and 4 ranks) on the unit cube sx2.msh:
!>   - one increment of an isotropic elastic cube (cubic constants with c44 = (c11 - c12)/2)
!>     under uniaxial_minimal along z at v2's elastic-case rate and increment (strain
!>     2.5e-4): the successive approximation converges, every element's stress is uniaxial
!>     (sigma_zz = E eps to O(eps), the other components small), the elastic strain is
!>     sigma_zz / E along z and -nu sigma_zz / E across;
!>   - face forces: F(z1) = sigma_zz A(z1) along z, F(z0) = -F(z1), forces on the lateral
!>     faces near zero, areas near 1;
!>   - mid-node linearization: leaves the straight-edged mesh unchanged (to the file's 12
!>     decimals); mid-nodes displaced by a smooth field are reset to the same coordinates as
!>     those of the unmodified mesh (the reset depends on the corners only), on the owned and
!>     the ghost nodes of every rank;
!>   - step 6.3, Newton-Raphson: with the elastic model the tangent equals the secant (D43), so
!>     an increment forced to switch to NR (global_nr_switch huge) gives bit-identical
!>     velocities and the same number of iterations as SA only (global_nr_switch 0); a plastic
!>     increment (fcc crystal plasticity, tau_0 10 MPa, so the cube yields during the
!>     increment) converged tightly (global_tol 1e-10) with NR uses NR iterations and at most 30
!>     iterations (19 measured; SA alone needs 449 and agrees to 3.4e-9: measured once, too slow
!>     for the suite); restarted from its solution with global_tol 1e-12 it moves less than 1e-8
!>     of the largest velocity; SA alone does not converge in 30 iterations.
!>
!> Source : FEPX (ROADMAP steps 4.4 and 6.3). The first argument is the mesh directory
!>          (tests/regression/cases/meshes). Tolerances: the end-of-increment configuration and
!>          the finite-strain factors make the solution differ from small-strain uniaxial stress
!>          by O(eps) = 2.5e-4 relative; the checks use 2e-3 (relative to sigma_zz).
!>
!> Locals:
!>   meshes        character(512)        mesh directory
!>   mesh          type(mesh_t)          rank-local mesh
!>   elt           type(element_t)       element descriptor
!>   model         type(elastic_model_t) the elastic model
!>   phases        type(phase_model_t)(1) the phase
!>   state_prev, state  type(material_state_t)  start and end states
!>   specs         type(bc_spec_t)(:)    preset conditions
!>   bc_mask       logical(:)            prescribed dofs
!>   bc_vel, vel   real(rk)(:)           prescribed and solved velocities [mm/s]
!>   opts          type(increment_opts_t) iteration parameters
!>   info          type(increment_info_t) iteration report
!>   coo, coo_mod  real(rk)(:, :)        coordinates, original and with displaced mid-nodes [mm]
!>   stress_c      real(rk)(:, :)        (6, num_elts) centroid stresses [MPa]
!>   loads, areas  real(rk)(:, :), (:)   face-set forces [N] and areas [mm^2]
!>   sig_zz, eps   real(rk)              expected stress [MPa] and strain [-]
!>   worst_*       real(rk)              largest deviations found (this rank)
!>   all_worst     real(rk)(3)           the same, over all ranks
!>   i_elt, i, j   integer(ik)           indices
!>   unit_log      integer               scratch log unit
!>   info_sa, info_nr  type(increment_info_t)  reports of the SA-only and NR runs (step 6.3)
!>   vel_sa, vel_nr    real(rk)(:)           their velocities [mm/s]
!>   cp_model      type(cp_model_t)      the crystal plasticity model of the plastic increment
!>   crys          type(crystal_t)       its slip geometry
!>   n_slip        integer               its number of slip systems
program test_increment

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_init, par_finalize, par_max
  use mesh_mod, only: mesh_t
  use mesh_blocks_mod, only: mesh_load_blocks
  use element_mod, only: element_t, element_init
  use crystal_mod, only: stiffness_cubic
  use material_state_mod, only: material_state_t, material_state_alloc
  use elastic_model_mod, only: elastic_model_t, elastic_model_init
  use crystal_mod, only: crystal_t, crystal_init, crystal_fcc
  use crystal_plasticity_model_mod, only: cp_model_t, cp_model_init
  use bc_mod, only: bc_spec_t, bc_preset_specs, bc_resolve, bc_preset_minimal
  use increment_mod, only: phase_model_t, increment_opts_t, increment_info_t, increment_solve, &
    & midnode_linearize, increment_converged
  use face_force_mod, only: face_forces
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary

  implicit none

  real(rk), parameter :: young = 200.0e3_rk, nu = 0.3_rk
  real(rk), parameter :: rate = 0.01_rk, dt = 0.025_rk

  character(len=512) :: meshes
  type(mesh_t) :: mesh
  type(element_t) :: elt
  type(elastic_model_t) :: model
  type(phase_model_t) :: phases(1)
  type(material_state_t) :: state_prev, state
  type(bc_spec_t), allocatable :: specs(:)
  logical, allocatable :: bc_mask(:)
  real(rk), allocatable :: bc_vel(:), vel(:)
  type(increment_opts_t) :: opts
  type(increment_info_t) :: info
  real(rk), allocatable :: coo(:, :), coo_mod(:, :), stress_c(:, :), loads(:, :), areas(:)
  real(rk) :: sig_zz, eps, lam, mu
  real(rk) :: worst_axial, worst_lateral, worst_strain, worst_mid, all_worst(3)
  integer(ik) :: i_elt, i, j
  integer :: unit_log
  type(increment_info_t) :: info_sa, info_nr
  real(rk), allocatable :: vel_sa(:), vel_nr(:)
  type(cp_model_t) :: cp_model
  type(crystal_t) :: crys
  integer :: n_slip

  call par_init()
  call get_command_argument(1, meshes)
  open (newunit=unit_log, status='scratch', action='readwrite', form='formatted')
  call log_set_unit(unit_log)
  call mesh_load_blocks(mesh, trim(meshes) // '/sx2.msh')
  call log_set_unit(output_unit)
  close (unit_log)
  call element_init(elt, mesh%elt_type)

  ! --- Mid-node linearization ---

  ! Unchanged on a straight-edged mesh, to the 12 decimals of Neper's coordinates (measured
  ! 5.0e-13: the stored mid-nodes are rounded midpoints).
  coo = mesh%coo_s
  call midnode_linearize(mesh, elt, coo)
  call par_max(maxval(abs(coo - mesh%coo_s)), worst_mid)
  call assert_true(worst_mid < 1.0e-12_rk, &
    & 'midnode: straight-edged mesh unchanged (file precision)')

  ! Mid-nodes (nodes 5-10 of the elements) displaced by a smooth field are restored.
  coo_mod = mesh%coo_s
  do i_elt = 1_ik, mesh%num_elts
    do i = 5_ik, 10_ik
      j = mesh%elt_nodes_loc(i, i_elt)
      coo_mod(:, j) = mesh%coo_s(:, j) + 1.0e-3_rk*[sin(7.0_rk*mesh%coo_s(2, j)), &
        & cos(5.0_rk*mesh%coo_s(3, j)), sin(3.0_rk*mesh%coo_s(1, j))]
    end do
  end do
  call midnode_linearize(mesh, elt, coo_mod)
  call par_max(maxval(abs(coo_mod - coo)), worst_mid)
  call assert_true(worst_mid < 1.0e-15_rk, 'midnode: displaced mid-nodes reset as the ' // &
    & 'unmodified mesh (owned and ghost)')

  ! --- One elastic increment ---

  ! Isotropic model (cubic with c44 = (c11 - c12)/2), state with the mesh orientations.
  lam = young*nu/((1.0_rk + nu)*(1.0_rk - 2.0_rk*nu))
  mu = young/(2.0_rk*(1.0_rk + nu))
  call elastic_model_init(model, stiffness_cubic(lam + 2.0_rk*mu, lam, mu))
  allocate (phases(1)%model, source=model)
  phases(1)%elts = [(i_elt, i_elt = 1_ik, mesh%num_elts)]
  call material_state_alloc(state_prev, elt%num_qpts, mesh%num_elts, 0_ik)
  do i_elt = 1_ik, mesh%num_elts
    do i = 1_ik, elt%num_qpts
      state_prev%rot_c2s(:, :, i, i_elt) = mesh%rot_c2s(:, :, i_elt)
    end do
  end do
  call phases(1)%model%init_state(phases(1)%elts, state_prev)
  state = state_prev

  ! uniaxial_minimal along z, the moving face at rate times the gage length (1).
  call bc_preset_specs(bc_preset_minimal, 3, rate*1.0_rk, specs)
  allocate (bc_mask(3*mesh%num_nodes_loc), bc_vel(3*mesh%num_nodes_loc))
  call bc_resolve(mesh, specs, bc_mask, bc_vel)
  vel = bc_vel

  ! Solve.
  opts%cg_tol = 1.0e-12_rk
  call increment_solve(mesh, elt, phases, mesh%coo_s, bc_mask, dt, opts, state_prev, state, &
    & vel, info)
  call assert_equal(info%status, increment_converged, 'increment converged')
  call assert_true(info%num_iters >= 2 .and. info%num_iters <= 5, 'converged in 2 to 5 iterations')

  ! Uniaxial stress in every element (centroid point), elastic strain.
  eps = rate*dt
  sig_zz = young*eps
  worst_axial = 0.0_rk
  worst_lateral = 0.0_rk
  worst_strain = 0.0_rk
  allocate (stress_c(6, mesh%num_elts))
  do i_elt = 1_ik, mesh%num_elts
    stress_c(:, i_elt) = state%stress_s(:, elt%qpt_centroid, i_elt)
    worst_axial = max(worst_axial, abs(stress_c(3, i_elt) - sig_zz)/sig_zz)
    worst_lateral = max(worst_lateral, maxval(abs(stress_c([1, 2, 4, 5, 6], i_elt)))/sig_zz)
    worst_strain = max(worst_strain, maxval(abs(state%strain_el_s(:, elt%qpt_centroid, i_elt) &
      & - [-nu*eps, -nu*eps, eps, 0.0_rk, 0.0_rk, 0.0_rk]))/eps)
  end do
  call par_max(worst_axial, all_worst(1))
  call par_max(worst_lateral, all_worst(2))
  call par_max(worst_strain, all_worst(3))
  worst_axial = all_worst(1)
  worst_lateral = all_worst(2)
  worst_strain = all_worst(3)
  call assert_true(worst_axial < 2.0e-3_rk, 'sigma_zz = E eps')
  call assert_true(worst_lateral < 2.0e-3_rk, 'other stress components small')
  call assert_true(worst_strain < 2.0e-3_rk, 'elastic strain (eps, -nu eps)')

  ! Face forces at the end-of-increment configuration.
  coo = mesh%coo_s + dt*reshape(vel, [3_ik, mesh%num_nodes_loc])
  allocate (loads(3, size(mesh%fasets)), areas(size(mesh%fasets)))
  call face_forces(mesh, elt, coo, stress_c, loads, areas)
  do i = 1_ik, size(mesh%fasets, kind=ik)
    call assert_close(areas(i), 1.0_rk, 'area ' // mesh%fasets(i)%label, abs_tol=2.0e-3_rk, &
      & rel_tol=0.0_rk)
    select case (mesh%fasets(i)%label)
    case ('z1')
      call assert_close(loads(3, i), sig_zz*areas(i), 'force z1 = sigma_zz A', &
        & abs_tol=2.0e-3_rk*sig_zz, rel_tol=0.0_rk)
      call assert_true(maxval(abs(loads(1:2, i))) < 2.0e-3_rk*sig_zz, 'z1: no lateral force')
    case ('z0')
      call assert_close(loads(3, i), -sig_zz, 'force z0 = -sigma_zz A', &
        & abs_tol=2.0e-3_rk*sig_zz, rel_tol=0.0_rk)
    case default
      call assert_true(maxval(abs(loads(:, i))) < 2.0e-3_rk*sig_zz, 'lateral face ' // &
        & mesh%fasets(i)%label // ': force near zero')
    end select
  end do

  ! --- Step 6.3: NR with the elastic model equals SA (tangent = secant, D43) ---
  vel_sa = bc_vel
  state = state_prev
  opts%global_nr_switch = 0.0_rk
  call increment_solve(mesh, elt, phases, mesh%coo_s, bc_mask, dt, opts, state_prev, state, &
    & vel_sa, info_sa)
  vel_nr = bc_vel
  state = state_prev
  opts%global_nr_switch = huge(1.0_rk)
  call increment_solve(mesh, elt, phases, mesh%coo_s, bc_mask, dt, opts, state_prev, state, &
    & vel_nr, info_nr)
  call assert_equal(info_sa%num_nr_iters, 0, 'elastic SA only: no NR iteration')
  call assert_true(info_nr%num_nr_iters > 0, 'elastic forced NR: NR iterations')
  call assert_equal(info_nr%num_iters, info_sa%num_iters, 'elastic NR = SA: iterations')
  call assert_true(.not. any(abs(vel_nr - vel_sa) > 0.0_rk), 'elastic NR = SA: velocities ' // &
    & '(bit-identical)')

  ! --- Step 6.3: a plastic increment, NR against SA only ---
  ! fcc crystal plasticity with tau_0 = 10 MPa (the elastic stress of the increment is about
  ! 50 MPa), m 0.02, gammadot_0 1.
  call crystal_init(crys, crystal_fcc)
  n_slip = crys%num_slip
  call cp_model_init(cp_model, stiffness_cubic(lam + 2.0_rk*mu, lam, mu), crys%schmid_sym_c, &
    & crys%schmid_skew_c, spread(0.02_rk, 1, n_slip), spread(1.0_rk, 1, n_slip), &
    & spread(10.0_rk, 1, n_slip))
  deallocate (phases(1)%model)
  allocate (phases(1)%model, source=cp_model)
  call material_state_alloc(state_prev, elt%num_qpts, mesh%num_elts, int(n_slip, ik))
  do i_elt = 1_ik, mesh%num_elts
    do i = 1_ik, elt%num_qpts
      state_prev%rot_c2s(:, :, i, i_elt) = mesh%rot_c2s(:, :, i_elt)
    end do
  end do
  call phases(1)%model%init_state(phases(1)%elts, state_prev)
  ! Converged by the relative test only (case 2's absolute threshold, 1.6e-7 mm/s here, would
  ! stop the runs early: the velocities are about 1e-2 mm/s).
  opts%global_tol = 1.0e-10_rk
  opts%global_tol_abs = 1.0e-30_rk

  ! NR (default switch): converged with NR iterations, in at most 30 iterations.
  vel_nr = bc_vel
  state = state_prev
  opts%global_max_iters = 50
  opts%global_nr_switch = 1.0e-2_rk
  call increment_solve(mesh, elt, phases, mesh%coo_s, bc_mask, dt, opts, state_prev, state, &
    & vel_nr, info_nr)
  call assert_equal(info_nr%status, increment_converged, 'plastic NR converged')
  call assert_true(info_nr%num_nr_iters > 0, 'plastic: NR iterations used')
  call assert_true(info_nr%num_iters <= 30, 'plastic: NR converged in at most 30 iterations')

  ! Restarted from its solution with a 100 times tighter tolerance, NR hardly moves: the
  ! solution is the fixed point of the secant residual.
  vel_sa = vel_nr
  state = state_prev
  opts%global_tol = 1.0e-12_rk
  call increment_solve(mesh, elt, phases, mesh%coo_s, bc_mask, dt, opts, state_prev, state, &
    & vel_sa, info_sa)
  call assert_equal(info_sa%status, increment_converged, 'plastic NR restart converged')
  call par_max(maxval(abs(vel_sa - vel_nr)), worst_mid)
  call par_max(maxval(abs(vel_nr)), all_worst(1))
  call assert_true(worst_mid < 1.0e-8_rk*all_worst(1), 'plastic: NR solution is the fixed point')

  ! SA alone, capped at 30 iterations, does not reach the tolerance (it needs about 450).
  vel_sa = bc_vel
  state = state_prev
  opts%global_tol = 1.0e-10_rk
  opts%global_max_iters = 30
  opts%global_nr_switch = 0.0_rk
  call increment_solve(mesh, elt, phases, mesh%coo_s, bc_mask, dt, opts, state_prev, state, &
    & vel_sa, info_sa)
  call assert_true(info_sa%status /= increment_converged .and. info_sa%num_nr_iters == 0, &
    & 'plastic: SA alone not converged in 30 iterations')
  if (mesh%num_nodes_owned >= 0_ik) write (output_unit, '(a,i0,a,i0,a,es10.3)') &
    & 'plastic increment: NR ', info_nr%num_iters, ' iterations (', info_nr%num_nr_iters, &
    & ' NR); restart change ', worst_mid/all_worst(1)

  call test_summary(verbose=(mesh%num_nodes_owned >= 0_ik))
  call par_finalize()

end program test_increment
