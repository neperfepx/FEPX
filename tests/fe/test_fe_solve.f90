! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> MPI tests of the element-by-element operations (ebe_mod), the Jacobi-PCG (pcg_mod) and the
!> Dirichlet conditions (bc_mod), on 1, 2 and 4 ranks.
!>
!> Source : FEPX (ROADMAP step 4.2).
!>
!> Usage: mpiexec -n N test_fe_solve <regression mesh dir> (tests/regression/cases/meshes).
!>
!> Mesh: poly_part1.msh (unit cube, 656 tet10), distributed with partitioning "even"; material:
!> isotropic [s] = [C] (E = 200 GPa, nu = 0.3; dt = 1 s, beta = 1) at every point, so that the
!> element matrices are the elastic stiffness. Checks:
!>   - matrix-vector product: for x = sin(0.1 id + k) (global node id, component k), [K]x on the
!>     owned dofs equals a serial assembly of the whole mesh read on every rank (relative 1e-12),
!>     with zero at the prescribed dofs when a mask is given;
!>   - presets (grip, symmetry, minimal along x, y, z; triaxial): the number of prescribed dofs
!>     over all ranks equals the union of the conditions' node sets in the whole mesh, and every
!>     prescribed value is the condition's; gage lengths are 1 mm; the triaxial strain-control
!>     conditions (step 5.3, D45) with a free face: faces, components and values;
!>   - solves (symmetry and minimal presets, axes x, y, z): [K]v = 0 with the prescribed
!>     velocities gives the exact homogeneous uniaxial-stress field v_a = e x_a, v_b = -nu e x_b,
!>     v_c = -nu e x_c (strain rate e = 1e-3/s) at every node: 1e-10 e with cg_tol 1e-12, and
!>     1e-6 e with v2's default cg_tol 1e-8 (converged within the iteration limit); the grip
!>     preset converges and keeps its prescribed values;
!>   - PCG edge cases: a zero right-hand side returns the initial guess with 0 iterations; an
!>     iteration limit of 1 returns status 1.
!>
!> Locals:
!>   ework  type(element_work_t)  element work arrays (element_work_alloc; D100)
!>   young, nu, rate  real(rk)        Young's modulus [MPa], Poisson's ratio, strain rate [1/s]
!>   meshes     character(512)        regression mesh directory
!>   mesh       type(mesh_t)          rank-local mesh
!>   ref        type(mesh_file_t)     whole mesh (serial references)
!>   elt        type(element_t)       element descriptor
!>   mat_iso    real(rk)(6,6)         isotropic [s] [MPa s]
!>   mat_s      real(rk)(:,:,:)       [s] at the quadrature points of an element
!>   k_elt      real(rk)(:,:,:)       local element matrices
!>   diag       real(rk)(:)           assembled diagonal
!>   unit_log   integer               scratch unit receiving the log
program test_fe_solve

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_init, par_finalize, par_sum
  use mesh_file_mod, only: mesh_file_t, mesh_file_load
  use mesh_mod, only: mesh_t
  use mesh_blocks_mod, only: mesh_load_blocks
  use element_mod, only: element_t, element_init
  use element_matrix_mod, only: element_stiffness
  use ebe_mod, only: ebe_matvec, ebe_diagonal
  use pcg_mod, only: pcg_solve
  use bc_mod, only: bc_spec_t, bc_preset_specs, bc_resolve, bc_gage_length, bc_preset_grip, &
    & bc_preset_minimal, bc_preset_triaxial, bc_triaxial_specs
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary
  use element_matrix_mod, only: element_work_t, element_work_alloc

  implicit none

  type(element_work_t) :: ework
  real(rk), parameter :: young = 200.0e3_rk
  real(rk), parameter :: nu = 0.3_rk
  real(rk), parameter :: rate = 1.0e-3_rk

  character(len=512) :: meshes
  type(mesh_t) :: mesh
  type(mesh_file_t) :: ref
  type(element_t) :: elt
  real(rk) :: mat_iso(6, 6)
  real(rk), allocatable :: mat_s(:, :, :), k_elt(:, :, :), diag(:)
  integer :: unit_log

  call par_init()
  call get_command_argument(1, meshes)

  ! Mesh (log to a scratch unit) and the whole mesh for the serial references.
  open (newunit=unit_log, status='scratch', action='readwrite', form='formatted')
  call log_set_unit(unit_log)
  call mesh_load_blocks(mesh, trim(meshes) // '/poly_part1.msh')
  call mesh_file_load(ref, trim(meshes) // '/poly_part1.msh')
  call log_set_unit(output_unit)
  close (unit_log)

  ! Isotropic [s] at every point and the element matrices.
  mat_iso = isotropic_mandel(young, nu)
  call element_init(elt, mesh%elt_type)
  allocate (mat_s(6, 6, elt%num_qpts))
  mat_s = spread(mat_iso, 3, int(elt%num_qpts))
  call build_matrices()

  call check_matvec()
  call check_presets()
  call check_solves()
  call check_pcg_edges()

  call test_summary()
  call par_finalize()

contains

  !> Element matrices of the local elements and the assembled diagonal.
  !>
  !> Locals:
  !>   i_elt  integer(ik)  local element
  subroutine build_matrices()

    integer(ik) :: i_elt

    allocate (k_elt(elt%num_dofs, elt%num_dofs, mesh%num_elts), diag(3*mesh%num_nodes_loc))
    do i_elt = 1_ik, mesh%num_elts
      call element_work_alloc(ework, elt)
      call element_stiffness(elt, mesh%coo_s(:, mesh%elt_nodes_loc(:, i_elt)), mat_s, &
        & k_elt(:, :, i_elt), ework)
    end do
    call ebe_diagonal(mesh, k_elt, diag)

  end subroutine build_matrices

  !> Matrix-vector product against a serial assembly of the whole mesh.
  !>
  !> Locals:
  !>   x, y          real(rk)(:)       local vectors
  !>   x_ref, y_ref  real(rk)(:,:)     (3, num_nodes) whole-mesh vectors
  !>   k_ref         real(rk)(:,:)     element matrix of the whole mesh's element
  !>   mask          logical(:)        prescribed dofs (every third dof)
  !>   i, e, a, k    integer           indices
  !>   ok            logical           all owned entries agree
  !>   scale         real(rk)          max |y_ref|
  subroutine check_matvec()

    real(rk), allocatable :: x(:), y(:), x_ref(:, :), y_ref(:, :), k_ref(:, :), xe(:)
    logical, allocatable :: mask(:)
    integer :: i, e, k
    integer(ik) :: g
    logical :: ok
    real(rk) :: scale

    ! Input vector from the global node ids.
    allocate (x(3*mesh%num_nodes_loc), y(3*mesh%num_nodes_loc), mask(3*mesh%num_nodes_loc))
    do i = 1, int(mesh%num_nodes_loc)
      do k = 1, 3
        x(3*(i - 1) + k) = sin(0.1_rk*mesh%node_id_glob(i) + k)
      end do
    end do

    ! Serial reference on the whole mesh.
    allocate (x_ref(3, ref%num_nodes), y_ref(3, ref%num_nodes), k_ref(elt%num_dofs, elt%num_dofs))
    do g = 1_ik, ref%num_nodes
      x_ref(:, g) = [(sin(0.1_rk*g + k), k = 1, 3)]
    end do
    y_ref = 0.0_rk
    do e = 1, int(ref%num_elts)
      call element_work_alloc(ework, elt)
      call element_stiffness(elt, ref%coo(:, ref%elt_nodes(:, e)), mat_s, k_ref, ework)
      xe = reshape(x_ref(:, ref%elt_nodes(:, e)), [elt%num_dofs])
      xe = matmul(k_ref, xe)
      y_ref(:, ref%elt_nodes(:, e)) = y_ref(:, ref%elt_nodes(:, e)) + reshape(xe, &
        & [3, int(elt%num_nodes)])
    end do
    scale = maxval(abs(y_ref))

    ! Distributed product, no prescribed dofs.
    mask = .false.
    call ebe_matvec(mesh, k_elt, mask, x, y)
    ok = .true.
    do i = 1, int(mesh%num_nodes_owned)
      if (any(abs(y(3*i - 2:3*i) - y_ref(:, mesh%node_id_glob(i))) > 1e-12_rk*scale)) ok = .false.
    end do
    call assert_true(ok, 'matvec: owned entries = serial assembly')

    ! With a mask: zero at the masked dofs, unchanged elsewhere.
    mask = .false.
    mask(1::3) = .true.
    call ebe_matvec(mesh, k_elt, mask, x, y)
    ok = .true.
    do i = 1, int(mesh%num_nodes_owned)
      if (abs(y(3*i - 2)) > 0.0_rk) ok = .false.
      if (any(abs(y(3*i - 1:3*i) - y_ref(2:3, mesh%node_id_glob(i))) > 1e-12_rk*scale)) &
        & ok = .false.
    end do
    call assert_true(ok, 'matvec: zero at the prescribed dofs')

  end subroutine check_matvec

  !> Prescribed dofs of the presets against the whole mesh's node sets; gage lengths.
  !>
  !> Locals:
  !>   preset, axis, i_spec, i_set, i  integer  loop indices
  !>   specs       type(bc_spec_t)(:)  conditions
  !>   bc_mask     logical(:)          local prescribed dofs
  !>   bc_vel      real(rk)(:)         local prescribed values
  !>   expected    logical(:,:)        (3, num_nodes) prescribed in the whole mesh
  !>   num_owned, num_all  integer     prescribed owned dofs on this rank, over all ranks
  !>   ok          logical             values agree
  subroutine check_presets()

    type(bc_spec_t), allocatable :: specs(:)
    logical, allocatable :: bc_mask(:), expected(:, :)
    real(rk), allocatable :: bc_vel(:), expected_vel(:, :)
    integer :: preset, axis, i_spec, i_set, i, num_owned, num_all
    logical :: ok

    allocate (bc_mask(3*mesh%num_nodes_loc), bc_vel(3*mesh%num_nodes_loc), &
      & expected(3, ref%num_nodes), expected_vel(3, ref%num_nodes))
    ! (Allocated empty first: gfortran 15 -O3 warns of an unset descriptor otherwise.)
    allocate (specs(0))
    do axis = 1, 3
      call assert_close(bc_gage_length(mesh, axis), 1.0_rk, 'gage length 1', abs_tol=1e-12_rk)
    end do

    do preset = bc_preset_grip, bc_preset_triaxial
      do axis = 1, 3
        if (preset == bc_preset_triaxial .and. axis > 1) cycle
        call bc_preset_specs(preset, axis, 7.0_rk, specs)
        call bc_resolve(mesh, specs, bc_mask, bc_vel)

        ! Expected prescribed dofs and values in the whole mesh.
        expected = .false.
        expected_vel = 0.0_rk
        do i_spec = 1, size(specs)
          do i_set = 1, size(ref%nsets)
            if (ref%nsets(i_set)%label == specs(i_spec)%nset) then
              expected(specs(i_spec)%comp, ref%nsets(i_set)%nodes) = .true.
              expected_vel(specs(i_spec)%comp, ref%nsets(i_set)%nodes) = specs(i_spec)%value
            end if
          end do
        end do

        ! Count over the ranks (owned dofs) and values.
        num_owned = count(bc_mask(1:3*mesh%num_nodes_owned))
        call par_sum(num_owned, num_all)
        call assert_equal(num_all, count(expected), 'preset prescribed dof count')
        ok = .true.
        do i = 1, int(mesh%num_nodes_loc)
          if (any(bc_mask(3*i - 2:3*i) .neqv. expected(:, mesh%node_id_glob(i)))) ok = .false.
          if (any(abs(bc_vel(3*i - 2:3*i) - expected_vel(:, mesh%node_id_glob(i))) > 0.0_rk)) &
            & ok = .false.
        end do
        call assert_true(ok, 'preset prescribed dofs and values')
      end do
    end do

    ! Triaxial strain control (D45): symmetry faces x0 y0 z0, then the prescribed faces x1 and
    ! z1 with their velocities; the free face y1 has no condition.
    call bc_triaxial_specs([1.0_rk, 2.0_rk, 3.0_rk], [.true., .false., .true.], specs)
    call assert_equal(size(specs), 5, 'triaxial: 3 symmetry and 2 moving faces')
    ok = specs(1)%nset == 'x0' .and. specs(1)%comp == 1 .and. specs(2)%nset == 'y0' .and. &
      & specs(2)%comp == 2 .and. specs(3)%nset == 'z0' .and. specs(3)%comp == 3 .and. &
      & all(abs([specs(1)%value, specs(2)%value, specs(3)%value]) < tiny(1.0_rk))
    call assert_true(ok, 'triaxial: symmetry faces with zero normal velocity')
    ok = specs(4)%nset == 'x1' .and. specs(4)%comp == 1 .and. specs(5)%nset == 'z1' .and. &
      & specs(5)%comp == 3
    call assert_true(ok, 'triaxial: moving faces x1 and z1, none on y1')
    call assert_close([specs(4)%value, specs(5)%value], [1.0_rk, 3.0_rk], &
      & 'triaxial: moving-face velocities')

  end subroutine check_presets

  !> Uniaxial-stress solves against the exact homogeneous field.
  !>
  !> Locals:
  !>   preset, axis, i, k, b, c, i_tol  integer  loop indices and lateral axes
  !>   specs      type(bc_spec_t)(:)  conditions
  !>   bc_mask    logical(:)          prescribed dofs
  !>   bc_vel, rhs, delta, vel, y  real(rk)(:)  vectors
  !>   exact      real(rk)(3)         exact velocity of a node
  !>   num_iters, status  integer     PCG results
  !>   err        real(rk)            largest velocity error over the owned nodes
  !>   tol        real(rk)            cg_tol of the solve
  subroutine check_solves()

    type(bc_spec_t), allocatable :: specs(:)
    logical, allocatable :: bc_mask(:)
    real(rk), allocatable :: bc_vel(:), rhs(:), delta(:), vel(:), y(:)
    real(rk) :: exact(3), err, tol
    integer :: preset, axis, i, i_tol, num_iters, status

    allocate (bc_mask(3*mesh%num_nodes_loc), bc_vel(3*mesh%num_nodes_loc), &
      & rhs(3*mesh%num_nodes_loc), delta(3*mesh%num_nodes_loc), vel(3*mesh%num_nodes_loc), &
      & y(3*mesh%num_nodes_loc))
    ! (Allocated empty first: gfortran 15 -O3 warns of an unset descriptor otherwise.)
    allocate (specs(0))

    do preset = bc_preset_grip, bc_preset_minimal
      do axis = 1, 3
        do i_tol = 1, 2
          tol = merge(1e-12_rk, 1e-8_rk, i_tol == 1)
          call bc_preset_specs(preset, axis, rate*bc_gage_length(mesh, axis), specs)
          call bc_resolve(mesh, specs, bc_mask, bc_vel)

          ! [K](v_bc + delta) = 0 on the free dofs: [K] delta = -[K] v_bc.
          call ebe_matvec(mesh, k_elt, bc_mask, bc_vel, y)
          rhs = -y
          delta = 0.0_rk
          call pcg_solve(mesh, k_elt, bc_mask, diag, rhs, delta, tol, 16000, num_iters, status)
          call assert_equal(status, 0, 'solve converged')
          vel = bc_vel + delta
          call assert_true(all(abs(vel - bc_vel) <= 0.0_rk .or. .not. bc_mask), &
            & 'solve keeps the prescribed values')
          if (preset == bc_preset_grip) cycle

          ! Exact homogeneous uniaxial-stress field.
          err = 0.0_rk
          do i = 1, int(mesh%num_nodes_owned)
            exact = -nu*rate*mesh%coo_s(:, i)
            exact(axis) = rate*mesh%coo_s(axis, i)
            err = max(err, maxval(abs(vel(3*i - 2:3*i) - exact)))
          end do
          call assert_close(err, 0.0_rk, 'solve: exact uniaxial field', &
            & abs_tol=merge(1e-10_rk, 1e-6_rk, i_tol == 1)*rate)
        end do
      end do
    end do

  end subroutine check_solves

  !> PCG with a zero right-hand side and with an iteration limit of 1.
  !>
  !> Locals:
  !>   bc_mask   logical(:)    no prescribed dofs except the corner conditions
  !>   rhs, x    real(rk)(:)   vectors
  !>   specs     type(bc_spec_t)(:)  minimal conditions
  !>   bc_vel    real(rk)(:)   prescribed values
  !>   num_iters, status  integer  PCG results
  subroutine check_pcg_edges()

    type(bc_spec_t), allocatable :: specs(:)
    logical, allocatable :: bc_mask(:)
    real(rk), allocatable :: rhs(:), x(:), bc_vel(:)
    integer :: num_iters, status

    allocate (bc_mask(3*mesh%num_nodes_loc), bc_vel(3*mesh%num_nodes_loc), &
      & rhs(3*mesh%num_nodes_loc), x(3*mesh%num_nodes_loc))
    ! (Allocated empty first: gfortran 15 -O3 warns of an unset descriptor otherwise.)
    allocate (specs(0))
    call bc_preset_specs(bc_preset_minimal, 3, rate, specs)
    call bc_resolve(mesh, specs, bc_mask, bc_vel)

    ! Zero right-hand side.
    rhs = 0.0_rk
    x = 0.0_rk
    call pcg_solve(mesh, k_elt, bc_mask, diag, rhs, x, 1e-8_rk, 16000, num_iters, status)
    call assert_equal(num_iters, 0, 'pcg zero rhs: no iterations')
    call assert_equal(status, 0, 'pcg zero rhs: status 0')
    call assert_true(all(abs(x) < tiny(1.0_rk)), 'pcg zero rhs: zero solution')

    ! Iteration limit.
    rhs = 1.0_rk
    x = 0.0_rk
    call pcg_solve(mesh, k_elt, bc_mask, diag, rhs, x, 1e-12_rk, 1, num_iters, status)
    call assert_equal(status, 1, 'pcg iteration limit: status 1')

  end subroutine check_pcg_edges

  !> Isotropic Mandel stiffness lambda {I}{I}^T + 2 mu [I].
  !>
  !> Arguments:
  !>   e_mod, poisson  (in)  real(rk)  Young's modulus [MPa], Poisson's ratio [-]
  !> Locals:
  !>   lam, mu   real(rk)  Lame constants [MPa]
  !>   i, j      integer   indices
  pure function isotropic_mandel(e_mod, poisson) result(stiff)

    real(rk), intent(in) :: e_mod, poisson
    real(rk) :: stiff(6, 6)

    real(rk) :: lam, mu
    integer :: i, j

    lam = e_mod*poisson/((1.0_rk + poisson)*(1.0_rk - 2.0_rk*poisson))
    mu = e_mod/(2.0_rk*(1.0_rk + poisson))
    stiff = 0.0_rk
    do i = 1, 3
      do j = 1, 3
        stiff(i, j) = lam
      end do
    end do
    do i = 1, 6
      stiff(i, i) = stiff(i, i) + 2.0_rk*mu
    end do

  end function isotropic_mandel

end program test_fe_solve
