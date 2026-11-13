! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> MPI test of the time stepping of step 5.2 (1, 2 and 4 ranks) on the unit cube sx2.msh with
!> every node prescribed, an HCP elastic crystal (Ti constants, C33 independent) in the mesh
!> orientation: a homogeneous pre-strain increment, then a rigid rotation of the strained body
!> in num_rot increments, with the geometry updated after each increment (x = x + dt v) as the
!> fepx driver does. Checks after the rotation (Q = exp(W dt) the rotation of one increment,
!> Q^N that of the whole rotation):
!>   - geometry: the coordinates are Q^N times the pre-strained coordinates (1e-12 mm);
!>   - lattice rotation by the continuum spin: rot_c2s = Q^N rot_c2s_0 at every quadrature
!>     point (1e-4; the lattice also turns with the elastic stretch, (R2), about e^e theta, and
!>     the spin of each increment's velocities is sin(theta)/dt, not theta/dt);
!>   - elastic strain carried with the lattice: the centroid stress is Q^N sigma_0 Q^N^T, sigma_0
!>     the stress after the pre-strain (1e-2 of max |sigma_0|: the rate of a finite rotation
!>     increment, evaluated in the end configuration, has a spurious symmetric part theta^2/(2
!>     dt), which accumulates to num_rot theta^2/2 = 1.5e-6 strain, 0.15 % of the pre-strain;
!>     rotating the stress changes it by about 6 %, so a missing or wrong-signed rotation fails).
!>
!> Source : FEPX (ROADMAP step 5.2). The first argument is the mesh directory
!>          (tests/regression/cases/meshes).
!>
!> Locals:
!>   num_rot, theta, dt_rot  integer, real(rk)  rotation increments, angle per increment [rad],
!>                                    time increment [s]
!>   axis          real(rk)(3)        rotation axis (unit)
!>   defrate       real(rk)(3, 3)     pre-strain rate (symmetric) [1/s]
!>   meshes        character(512)     mesh directory
!>   mesh          type(mesh_t)       rank-local mesh
!>   elt           type(element_t)    element descriptor
!>   model         type(elastic_model_t)  the elastic model
!>   phases        type(phase_model_t)(1) the phase
!>   state_prev, state, state0  type(material_state_t)  start, end and pre-strained states
!>   bc_mask       logical(:)         prescribed dofs (all)
!>   vel           real(rk)(:)        nodal velocities [mm/s]
!>   opts          type(increment_opts_t)  iteration parameters
!>   info          type(increment_info_t)  iteration report
!>   coo, coo0     real(rk)(:, :)     current and pre-strained coordinates [mm]
!>   q_step, q_all real(rk)(3, 3)     rotation of one increment, of all of them
!>   q6            real(rk)(6, 6)     Mandel rotation of q_all
!>   worst         real(rk)(3)        largest deviations (coordinates, orientation, stress), this
!>                                    rank; all_worst over the ranks
!>   sig_scale     real(rk)           max |sigma_0| over the ranks [MPa]
!>   i_incr, n, k  integer            increment, node, dof
!>   i_elt, i_qpt  integer(ik)        element, point
!>   unit_log      integer            scratch log unit
program test_rigid_rotation

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_init, par_finalize, par_max
  use mesh_mod, only: mesh_t
  use mesh_blocks_mod, only: mesh_load_blocks
  use element_mod, only: element_t, element_init
  use crystal_mod, only: stiffness_hexagonal
  use tensor_mod, only: skew_from_matrix, mandel_rotation
  use rotation_mod, only: rotation_exp_skew
  use material_state_mod, only: material_state_t, material_state_alloc
  use elastic_model_mod, only: elastic_model_t, elastic_model_init
  use increment_mod, only: phase_model_t, increment_opts_t, increment_info_t, increment_solve, &
    & increment_converged
  use assert_mod, only: assert_true, assert_equal, test_summary

  implicit none

  integer, parameter :: num_rot = 300
  real(rk), parameter :: theta = 1.0e-4_rk, dt_rot = 1.0e-3_rk
  real(rk), parameter :: defrate(3, 3) = 1.0e-3_rk*reshape([1.0_rk, 0.2_rk, 0.0_rk, &
    & 0.2_rk, -0.3_rk, 0.1_rk, 0.0_rk, 0.1_rk, -0.4_rk], [3, 3])

  character(len=512) :: meshes
  type(mesh_t) :: mesh
  type(element_t) :: elt
  type(elastic_model_t) :: model
  type(phase_model_t) :: phases(1)
  type(material_state_t) :: state_prev, state, state0
  logical, allocatable :: bc_mask(:)
  real(rk), allocatable :: vel(:), coo(:, :), coo0(:, :)
  type(increment_opts_t) :: opts
  type(increment_info_t) :: info
  real(rk) :: axis(3), q_step(3, 3), q_all(3, 3), q6(6, 6), worst(3), all_worst(3), sig_scale
  integer :: i_incr, n, unit_log
  integer(ik) :: i_elt, i_qpt

  call par_init()
  call get_command_argument(1, meshes)
  open (newunit=unit_log, status='scratch', action='readwrite', form='formatted')
  call log_set_unit(unit_log)
  call mesh_load_blocks(mesh, trim(meshes) // '/sx2.msh')
  call log_set_unit(output_unit)
  close (unit_log)
  call element_init(elt, mesh%elt_type)

  ! HCP model (Ti constants, C33 independent), state with the mesh orientation.
  call elastic_model_init(model, stiffness_hexagonal(162400.0_rk, 92000.0_rk, 69000.0_rk, &
    & 180700.0_rk, 46700.0_rk))
  allocate (phases(1)%model, source=model)
  phases(1)%elts = [(i_elt, i_elt = 1_ik, mesh%num_elts)]
  call material_state_alloc(state, elt%num_qpts, mesh%num_elts, 0_ik)
  do i_elt = 1_ik, mesh%num_elts
    do i_qpt = 1_ik, elt%num_qpts
      state%rot_c2s(:, :, i_qpt, i_elt) = mesh%rot_c2s(:, :, i_elt)
    end do
  end do
  call phases(1)%model%init_state(phases(1)%elts, state)

  ! Every dof prescribed; tight solver tolerances.
  allocate (bc_mask(3*mesh%num_nodes_loc), vel(3*mesh%num_nodes_loc))
  bc_mask = .true.
  opts%cg_tol = 1.0e-12_rk
  coo = mesh%coo_s

  ! --- Homogeneous pre-strain: v = D x, one increment of 1 s ---

  do n = 1, int(mesh%num_nodes_loc)
    vel(3*n - 2:3*n) = matmul(defrate, coo(:, n))
  end do
  call advance(1.0_rk)
  state0 = state
  coo0 = coo

  ! --- Rigid rotation: v = (Q - I) x / dt, so that x = Q x after each increment ---

  axis = [1.0_rk, 2.0_rk, 3.0_rk]/sqrt(14.0_rk)
  q_step = rotation_exp_skew(skew_from_matrix(theta*cross_matrix(axis)))
  q_all = reshape([1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk], &
    & [3, 3])
  do i_incr = 1, num_rot
    do n = 1, int(mesh%num_nodes_loc)
      vel(3*n - 2:3*n) = (matmul(q_step, coo(:, n)) - coo(:, n))/dt_rot
    end do
    call advance(dt_rot)
    q_all = matmul(q_step, q_all)
  end do

  ! --- Checks ---

  ! Geometry.
  worst(1) = maxval(abs(coo - matmul(q_all, coo0)))

  ! Lattice orientation at every point; stress at the centroid points.
  q6 = mandel_rotation(q_all)
  worst(2) = 0.0_rk
  worst(3) = 0.0_rk
  sig_scale = 0.0_rk
  do i_elt = 1_ik, mesh%num_elts
    do i_qpt = 1_ik, elt%num_qpts
      worst(2) = max(worst(2), maxval(abs(state%rot_c2s(:, :, i_qpt, i_elt) - &
        & matmul(q_all, state0%rot_c2s(:, :, i_qpt, i_elt)))))
    end do
    worst(3) = max(worst(3), maxval(abs(state%stress_s(:, elt%qpt_centroid, i_elt) - &
      & matmul(q6, state0%stress_s(:, elt%qpt_centroid, i_elt)))))
    sig_scale = max(sig_scale, maxval(abs(state0%stress_s(:, elt%qpt_centroid, i_elt))))
  end do
  do n = 1, 3
    call par_max(worst(n), all_worst(n))
  end do
  call par_max(sig_scale, worst(1))
  sig_scale = worst(1)
  call assert_true(sig_scale > 50.0_rk, 'pre-strain gives stresses of order 100 MPa')
  call assert_true(all_worst(1) < 1.0e-12_rk, 'coordinates rotated rigidly')
  call assert_true(all_worst(2) < 1.0e-4_rk, 'lattice rotated with the body')
  call assert_true(all_worst(3) < 1.0e-2_rk*sig_scale, 'stress rotated with the lattice')

  call test_summary()
  call par_finalize()

contains

  !> One increment with the prescribed velocities vel, then the geometry update.
  !>
  !> Arguments:
  !>   dt  (in)  real(rk)  time increment [s]
  subroutine advance(dt)

    real(rk), intent(in) :: dt

    ! Solve from the current coordinates, then x = x + dt v.
    state_prev = state
    call increment_solve(mesh, elt, phases, coo, bc_mask, dt, opts, state_prev, state, vel, &
      & info)
    call assert_equal(info%status, increment_converged, 'increment converged')
    coo = coo + dt*reshape(vel, [3_ik, mesh%num_nodes_loc])

  end subroutine advance

  !> Skew matrix [a]x of a vector (a x v = [a]x v).
  !>
  !> Arguments:
  !>   a  (in)  real(rk)(3)  vector
  pure function cross_matrix(a) result(m)

    real(rk), intent(in) :: a(3)
    real(rk) :: m(3, 3)

    m = reshape([0.0_rk, a(3), -a(2), -a(3), 0.0_rk, a(1), a(2), -a(1), 0.0_rk], [3, 3])

  end function cross_matrix

end program test_rigid_rotation
