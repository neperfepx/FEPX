! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Lamé verification (dev/theory/addendum_lame.md, approved; ROADMAP step 4.4, D43): one elastic
!> increment of a quarter annulus in plane strain with a prescribed inner radial displacement and
!> a free outer surface, against the analytic field (L3)-(L8), for nu = 0.3, 0.49, 0.499, on
!> tet10 and hex20 meshes at two refinements (tests/verification/cases/meshes/annulus_*.msh).
!>
!> Checks, per mesh and nu:
!>   - the increment converges;
!>   - stress error: max over the elements of |sigma - sigma_analytic| at the centroid point
!>     (every Mandel component), relative to max |sigma_analytic| (sigma_theta at r = a);
!>   - displacement error: max over the nodes of |u - u_analytic|, relative to u_a;
!>   - each error is below its tolerance (per mesh and nu) and is smaller on the finer mesh.
!> The measured errors are printed (one line per mesh and nu); tolerances and measurements are
!> documented in tests/verification/README.md.
!>
!> Set-up (addendum §3): v_z = 0 on z0 and z1 (plane strain), v_x = 0 on x0, v_y = 0 on y0,
!> v = v_a (x, y)/r on the inner surface; isotropic elasticity through the cubic model with
!> c11 = lambda + 2 mu, c12 = lambda, c44 = mu; one increment dt = 1 s with u_a/a = 1e-6 (the
!> end-of-increment configuration then differs from the reference one by 1e-6).
!>
!> Source : FEPX (ROADMAP step 4.4). The first argument is the mesh directory
!>          (tests/verification/cases/meshes).
!>
!> Locals:
!>   meshes      character(512)   mesh directory
!>   names       character(16)(4) mesh names (coarse and fine, tet10 and hex20)
!>   nus         real(rk)(3)      Poisson's ratios
!>   err_s, err_u  real(rk)(4, 3) stress and displacement errors per mesh and nu
!>   i_mesh, i_nu  integer        indices
!>   converged   logical          whether a run converged
!>   label       character(96)    check label
!>
!> The module lame_run_mod below holds one run (module procedure).
!>
!> lame_run_mod public routines: lame_run.
module lame_run_mod

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use constants_mod, only: sqrt2
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_max
  use mesh_mod, only: mesh_t
  use mesh_blocks_mod, only: mesh_load_blocks
  use element_mod, only: element_t, element_init
  use crystal_mod, only: stiffness_cubic
  use material_state_mod, only: material_state_t, material_state_alloc
  use elastic_model_mod, only: elastic_model_t, elastic_model_init
  use increment_mod, only: phase_model_t, increment_opts_t, increment_info_t, increment_solve, &
    & increment_converged

  implicit none

  private

  public :: lame_run

contains

  !> Run one increment on a mesh for one nu and return the errors against (L3)-(L8).
  !>
  !> Arguments:
  !>   mesh_file   (in)  character  mesh file
  !>   nu          (in)  real(rk)   Poisson's ratio [-]
  !>   err_stress  (out) real(rk)   max centroid stress error / max |sigma_analytic| [-]
  !>   err_disp    (out) real(rk)   max nodal displacement error / u_a [-]
  !>   converged   (out) logical    whether the increment converged
  !> Locals:
  !>   young, a, b, dt, u_a  real(rk)  modulus [MPa], radii [mm], increment [s], inner
  !>                                   displacement [mm]
  !>   lam, mu, dd, cfac     real(rk)  Lamé constants [MPa], D = (1 - 2 nu) a^2 + b^2 [mm^2],
  !>                                   E u_a a / ((1 + nu) D) [MPa]
  !>   mesh, elt, model, phases, state_prev, state, opts, info   the problem
  !>   bc_mask, vel          logical(:), real(rk)(:)  prescribed dofs and velocities
  !>   i_set, i, n, i_elt, i_qpt  integer(ik)  indices
  !>   x                     real(rk)(3)     a position [mm]
  !>   r, c, s               real(rk)        radius [mm], cos and sin of the angle
  !>   srr, stt, szz         real(rk)        analytic stresses [MPa]
  !>   sig_ref               real(rk)(6)     analytic stress, Mandel [MPa]
  !>   u_ref                 real(rk)(3)     analytic displacement [mm]
  !>   err_s_loc, err_u_loc  real(rk)        this rank's errors
  !>   unit_log              integer         scratch log unit
  subroutine lame_run(mesh_file, nu, err_stress, err_disp, converged)

    character(len=*), intent(in) :: mesh_file
    real(rk), intent(in) :: nu
    real(rk), intent(out) :: err_stress, err_disp
    logical, intent(out) :: converged

    real(rk), parameter :: young = 200.0e3_rk, a = 1.0_rk, b = 2.0_rk, dt = 1.0_rk
    real(rk), parameter :: u_a = 1.0e-6_rk*a
    real(rk) :: lam, mu, dd, cfac
    type(mesh_t) :: mesh
    type(element_t) :: elt
    type(elastic_model_t) :: model
    type(phase_model_t) :: phases(1)
    type(material_state_t) :: state_prev, state
    type(increment_opts_t) :: opts
    type(increment_info_t) :: info
    logical, allocatable :: bc_mask(:)
    real(rk), allocatable :: vel(:)
    integer(ik) :: i_set, i, n, i_elt, i_qpt
    real(rk) :: x(3), r, c, s, srr, stt, szz, sig_ref(6), u_ref(3), err_s_loc, err_u_loc
    integer :: unit_log

    ! Mesh and element.
    open (newunit=unit_log, status='scratch', action='readwrite', form='formatted')
    call log_set_unit(unit_log)
    call mesh_load_blocks(mesh, mesh_file)
    call log_set_unit(output_unit)
    close (unit_log)
    call element_init(elt, mesh%elt_type)

    ! Isotropic model and state (identity orientation from the mesh).
    lam = young*nu/((1.0_rk + nu)*(1.0_rk - 2.0_rk*nu))
    mu = young/(2.0_rk*(1.0_rk + nu))
    call elastic_model_init(model, stiffness_cubic(lam + 2.0_rk*mu, lam, mu))
    allocate (phases(1)%model, source=model)
    phases(1)%elts = [(i_elt, i_elt = 1_ik, mesh%num_elts)]
    call material_state_alloc(state_prev, elt%num_qpts, mesh%num_elts, 0_ik)
    do i_elt = 1_ik, mesh%num_elts
      do i_qpt = 1_ik, elt%num_qpts
        state_prev%rot_c2s(:, :, i_qpt, i_elt) = mesh%rot_c2s(:, :, i_elt)
      end do
    end do
    call phases(1)%model%init_state(phases(1)%elts, state_prev)
    state = state_prev

    ! Conditions (addendum §3): plane strain, symmetry, radial inner velocity.
    allocate (bc_mask(3*mesh%num_nodes_loc), vel(3*mesh%num_nodes_loc))
    bc_mask = .false.
    vel = 0.0_rk
    do i_set = 1_ik, size(mesh%nsets, kind=ik)
      do i = 1_ik, size(mesh%nsets(i_set)%nodes, kind=ik)
        n = mesh%nsets(i_set)%nodes(i)
        select case (mesh%nsets(i_set)%label)
        case ('z0', 'z1')
          bc_mask(3*n) = .true.
        case ('x0')
          bc_mask(3*n - 2) = .true.
        case ('y0')
          bc_mask(3*n - 1) = .true.
        case ('inner')
          r = norm2(mesh%coo_s(1:2, n))
          bc_mask(3*n - 2:3*n - 1) = .true.
          vel(3*n - 2:3*n - 1) = (u_a/dt)*mesh%coo_s(1:2, n)/r
        end select
      end do
    end do
    ! Symmetry planes win where they meet the inner surface (their velocity is zero there up
    ! to round-off of the coordinates).
    do i_set = 1_ik, size(mesh%nsets, kind=ik)
      do i = 1_ik, size(mesh%nsets(i_set)%nodes, kind=ik)
        n = mesh%nsets(i_set)%nodes(i)
        if (mesh%nsets(i_set)%label == 'x0') vel(3*n - 2) = 0.0_rk
        if (mesh%nsets(i_set)%label == 'y0') vel(3*n - 1) = 0.0_rk
      end do
    end do

    ! One increment.
    opts%cg_tol = 1.0e-12_rk
    opts%global_tol = 1.0e-10_rk
    call increment_solve(mesh, elt, phases, mesh%coo_s, bc_mask, dt, opts, state_prev, state, &
      & vel, info)
    converged = info%status == increment_converged

    ! Stress error at the centroid points (reference configuration positions).
    dd = (1.0_rk - 2.0_rk*nu)*a**2 + b**2
    cfac = young*u_a*a/((1.0_rk + nu)*dd)
    err_s_loc = 0.0_rk
    do i_elt = 1_ik, mesh%num_elts
      x = matmul(mesh%coo_s(:, mesh%elt_nodes_loc(:, i_elt)), elt%shape_qpt(:, elt%qpt_centroid))
      r = norm2(x(1:2))
      c = x(1)/r
      s = x(2)/r
      srr = cfac*(1.0_rk - b**2/r**2)
      stt = cfac*(1.0_rk + b**2/r**2)
      szz = 2.0_rk*nu*cfac
      sig_ref = [srr*c**2 + stt*s**2, srr*s**2 + stt*c**2, szz, 0.0_rk, 0.0_rk, &
        & sqrt2*(srr - stt)*s*c]
      err_s_loc = max(err_s_loc, maxval(abs(state%stress_s(:, elt%qpt_centroid, i_elt) - &
        & sig_ref)))
    end do
    call par_max(err_s_loc, err_stress)
    err_stress = err_stress/(cfac*(1.0_rk + b**2/a**2))

    ! Displacement error at the owned nodes: u = dt v against u_r(r) (x, y)/r.
    err_u_loc = 0.0_rk
    do n = 1_ik, mesh%num_nodes_owned
      r = norm2(mesh%coo_s(1:2, n))
      u_ref = [((1.0_rk - 2.0_rk*nu)*u_a*a/dd*r + u_a*a*b**2/(dd*r))*mesh%coo_s(1:2, n)/r, &
        & 0.0_rk]
      err_u_loc = max(err_u_loc, maxval(abs(dt*vel(3*n - 2:3*n) - u_ref)))
    end do
    call par_max(err_u_loc, err_disp)
    err_disp = err_disp/u_a

  end subroutine lame_run

end module lame_run_mod

program test_lame

  use kinds_mod, only: rk
  use parallel_mod, only: par_init, par_finalize, my_rank, root_rank
  use assert_mod, only: assert_true, test_summary
  use lame_run_mod, only: lame_run

  implicit none

  ! Tolerances (relative), about 1.5 times the errors measured on 2026-10-09 (rows: the meshes
  ! of names; columns: the nu of nus; tests/verification/README.md). The stress errors at
  ! nu -> 1/2 are large (volumetric locking of the fully integrated elements, Q17): these
  ! tolerances record the present behaviour, so that any change is seen.
  real(rk), parameter :: tol_stress(4, 3) = reshape([ &
    & 1.7e-2_rk, 5.0e-3_rk, 2.1e-2_rk, 6.5e-3_rk, &
    & 1.8e-1_rk, 5.4e-2_rk, 3.0e-1_rk, 9.5e-2_rk, &
    & 1.7e+0_rk, 4.8e-1_rk, 3.0e+0_rk, 9.3e-1_rk], [4, 3])
  real(rk), parameter :: tol_disp(4, 3) = reshape([ &
    & 1.6e-3_rk, 2.7e-4_rk, 3.9e-5_rk, 3.0e-6_rk, &
    & 4.9e-3_rk, 7.2e-4_rk, 2.6e-5_rk, 2.4e-6_rk, &
    & 1.9e-2_rk, 3.2e-3_rk, 2.1e-5_rk, 2.0e-6_rk], [4, 3])

  character(len=512) :: meshes
  character(len=16), parameter :: names(4) = [character(len=16) :: 'annulus_tet10_n4', &
    & 'annulus_tet10_n8', 'annulus_hex20_n4', 'annulus_hex20_n8']
  real(rk), parameter :: nus(3) = [0.3_rk, 0.49_rk, 0.499_rk]
  real(rk) :: err_s(4, 3), err_u(4, 3)
  integer :: i_mesh, i_nu
  logical :: converged
  character(len=96) :: label

  call par_init()
  call get_command_argument(1, meshes)

  ! Every mesh and nu.
  do i_mesh = 1, 4
    do i_nu = 1, 3
      call lame_run(trim(meshes) // '/' // trim(names(i_mesh)) // '.msh', nus(i_nu), &
        & err_s(i_mesh, i_nu), err_u(i_mesh, i_nu), converged)
      write (label, '(a,a,f5.3)') trim(names(i_mesh)), ' nu ', nus(i_nu)
      if (my_rank == root_rank) write (*, '(a,a,es10.3,a,es10.3)') trim(label), &
        & ': stress error ', err_s(i_mesh, i_nu), ', displacement error ', err_u(i_mesh, i_nu)
      call assert_true(converged, trim(label) // ': converged')
      call assert_true(err_s(i_mesh, i_nu) < tol_stress(i_mesh, i_nu), trim(label) // &
        & ': stress error below tolerance')
      call assert_true(err_u(i_mesh, i_nu) < tol_disp(i_mesh, i_nu), trim(label) // &
        & ': displacement error below tolerance')
    end do
  end do

  ! Refinement lowers the errors.
  do i_nu = 1, 3
    do i_mesh = 1, 3, 2
      write (label, '(a,a,f5.3)') trim(names(i_mesh)), ' nu ', nus(i_nu)
      call assert_true(err_s(i_mesh + 1, i_nu) < err_s(i_mesh, i_nu), trim(label) // &
        & ': stress error decreases with refinement')
      call assert_true(err_u(i_mesh + 1, i_nu) < err_u(i_mesh, i_nu), trim(label) // &
        & ': displacement error decreases with refinement')
    end do
  end do

  call test_summary(verbose=(my_rank == root_rank))
  call par_finalize()

end program test_lame
