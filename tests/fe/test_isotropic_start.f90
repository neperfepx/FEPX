! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> MPI tests of the isotropic viscoplastic start (isotropic_start_mod; ROADMAP D78,
!> dev/theory/addendum_isotropic_viscoplastic.md §1, §5), on 1, 2 and 4 ranks.
!>
!> Source : FEPX (ROADMAP step 8.3b, D78).
!>
!> Usage: mpiexec -n N test_isotropic_start <verification mesh dir>
!>        (tests/verification/cases/meshes).
!>
!> Checks:
!>   - (I1): d . (2 mu [P_dev]) d equals v2's energy mu [(d11 - d22)^2 + 3 d'33^2 +
!>     4 (d12^2 + d13^2 + d23^2)] (solveit_isovp_mod3.f90, the 5-vector basis of
!>     solveit_isovp_mod2.f90::elt_stif_vp) for random rates (relative 1e-12);
!>   - (I2)-(I4) on random curved tet10 and hex20 elements: the deviatoric part of [k_e] equals
!>     v2's sum_q w_q J_q [B5_q]^T mu diag(1, 3, 4, 4, 4) [B5_q] and {g_e} equals v2's pcnst
!>     (relative 1e-12); [k_e] symmetric; [k_e] annihilates rigid motions; its penalty part
!>     annihilates isochoric affine fields; {g_e}.{V_e} = - V_e tr(L) for an affine field
!>     (relative 1e-12); r_e = mu / (6 V_e);
!>   - the solve on the periodic 2 x 1 x 1 box (small strain, 11 imposed at a, the other
!>     components free): with a homogeneous isotropic material every iterate is affine, and the
!>     three augmented-Lagrangian iterations give exactly D22 = D33 = -a/18 and zero shear
!>     (the lateral rates a/2, a/6, -a/18 of the iterations k = 1, 2, 3 from
!>     4 (b_k - a) + sum_(j <= k) (a + 2 b_j) = 0; the addendum's §5): Fbar-dot from the drive
!>     pairs to 1e-8 of a, the fluctuation v - Fbar-dot x (mean removed) below 1e-8 of the
!>     largest velocity, the periodicity relations kept.
!>
!> Locals:
!>   ework  type(element_work_t)  element work arrays (element_work_alloc; D100)
!>   meshes      character(512)  verification mesh directory
!>   unit_log    integer         scratch unit receiving the log
program test_isotropic_start

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds_mod, only: rk, ik
  use log_mod, only: log_set_unit
  use parallel_mod, only: par_init, par_finalize, par_max, par_min
  use mesh_mod, only: mesh_t
  use mesh_blocks_mod, only: mesh_load_blocks
  use element_mod, only: element_t, element_init, element_geometry_qpt, mandel_b_matrix
  use tet10_mod, only: tet10_gmsh_type, tet10_node_coo_ref
  use hex20_mod, only: hex20_gmsh_type, hex20_node_coo_ref
  use constraint_mod, only: constraint_t, constraint_set_rates, constraint_project, &
    & constraint_zero_mean, constraint_macro_rate
  use periodic_mod, only: periodic_build
  use isotropic_start_mod, only: iso_material_matrix, iso_element_system, &
    & isotropic_start_solve
  use tensor_test_utils_mod, only: seed_random, random_matrix, determinant3
  use assert_mod, only: assert_true, assert_close, test_summary
  use element_matrix_mod, only: element_work_t, element_work_alloc

  implicit none

  type(element_work_t) :: ework
  character(len=512) :: meshes
  integer :: unit_log

  call par_init()
  call seed_random(20261021)
  call get_command_argument(1, meshes)
  open (newunit=unit_log, status='scratch', action='readwrite', form='formatted')
  call log_set_unit(unit_log)

  call check_material()
  call check_element(tet10_gmsh_type, 'tet10')
  call check_element(hex20_gmsh_type, 'hex20')
  call check_periodic_solve(trim(meshes) // '/periodic_box_2x1x1.msh')

  call log_set_unit(output_unit)
  close (unit_log)
  call test_summary()
  call par_finalize()

contains

  !> (I1) against v2's 5-vector energy on random rates.
  !>
  !> Locals:
  !>   mu       real(rk)        viscosity [MPa s]
  !>   mat      real(rk)(6, 6)  2 mu [P_dev]
  !>   amat     real(rk)(3, 3)  random matrix (its symmetric part is the rate)
  !>   d        real(rk)(3, 3)  the rate [1/s]
  !>   dm       real(rk)(6)     Mandel rate
  !>   q        real(rk)(5)     v2's 5-vector (d11 - d22, d'33, d12, d13, d23)
  !>   e_new, e_v2  real(rk)    the two energies
  !>   i_try    integer         trial
  subroutine check_material()

    real(rk), parameter :: mu = 1.7_rk
    real(rk) :: mat(6, 6), amat(3, 3), d(3, 3), dm(6), q(5), e_new, e_v2
    integer :: i_try
    logical :: ok

    call iso_material_matrix(mu, mat)
    call assert_close(maxval(abs(mat - transpose(mat))), 0.0_rk, 'I1 symmetric', abs_tol=0.0_rk)
    call assert_close(maxval(abs(matmul(mat, [1.0_rk, 1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, &
      & 0.0_rk]))), 0.0_rk, 'I1 annihilates the spherical part', abs_tol=1.0e-14_rk)
    ok = .true.
    do i_try = 1, 20
      amat = random_matrix()
      d = 0.5_rk*(amat + transpose(amat))
      dm = [d(1, 1), d(2, 2), d(3, 3), sqrt(2.0_rk)*d(2, 3), sqrt(2.0_rk)*d(1, 3), &
        & sqrt(2.0_rk)*d(1, 2)]
      q = [d(1, 1) - d(2, 2), d(3, 3) - (d(1, 1) + d(2, 2) + d(3, 3))/3.0_rk, d(1, 2), &
        & d(1, 3), d(2, 3)]
      e_new = dot_product(dm, matmul(mat, dm))
      e_v2 = mu*(q(1)**2 + 3.0_rk*q(2)**2 + 4.0_rk*(q(3)**2 + q(4)**2 + q(5)**2))
      if (abs(e_new - e_v2) > 1.0e-12_rk*abs(e_v2)) ok = .false.
    end do
    call assert_true(ok, 'I1 energy equals v2''s 5-vector energy (20 random rates)')

  end subroutine check_material

  !> (I2)-(I4) on a random curved element of one type, against v2's formulas.
  !>
  !> Arguments:
  !>   gmsh_type  (in)  integer(ik)  element type
  !>   name       (in)  character    its name (labels)
  !> Locals:
  !>   mu        real(rk)              viscosity [MPa s]
  !>   elt       type(element_t)       descriptor
  !>   coo       real(rk)(:, :)        element node coordinates [mm]
  !>   k_e, k_dev, k_v2  real(rk)(:, :)  element matrix, its deviatoric part, v2's
  !>   g_e, g_v2 real(rk)(:)           pressure vector, v2's pcnst
  !>   r_e       real(rk)              penalty
  !>   b5        real(rk)(:, :)        v2's 5 x num_dofs matrix at a point
  !>   c5        real(rk)(5, 5)        v2's mu diag(1, 3, 4, 4, 4)
  !>   dshape_s  real(rk)(:, :)        shape derivatives [1/mm]
  !>   det_jac, vol  real(rk)          det J, volume
  !>   vel       real(rk)(:)           element velocities
  !>   lmat, wmat  real(rk)(3, 3)      a velocity gradient, a skew matrix
  !>   i_qpt, a, i  integer(ik)        point, node, dof
  subroutine check_element(gmsh_type, name)

    integer(ik), intent(in) :: gmsh_type
    character(len=*), intent(in) :: name

    real(rk), parameter :: mu = 1.0_rk
    type(element_t) :: elt
    real(rk), allocatable :: coo(:, :), k_e(:, :), k_dev(:, :), k_v2(:, :), g_e(:), g_v2(:)
    real(rk), allocatable :: b5(:, :), dshape_s(:, :), vel(:)
    real(rk) :: r_e, c5(5, 5), det_jac, vol, lmat(3, 3), wmat(3, 3), scale
    integer(ik) :: i_qpt, a, i

    call element_init(elt, gmsh_type)
    coo = random_element(elt)
    allocate (k_e(elt%num_dofs, elt%num_dofs), g_e(elt%num_dofs), vel(elt%num_dofs))
    call element_work_alloc(ework, elt)
    call iso_element_system(elt, coo, mu, k_e, g_e, r_e, ework)

    ! v2's matrices: the 5-vector form and pcnst (solveit_isovp_mod2.f90::elt_stif_vp).
    allocate (k_v2(elt%num_dofs, elt%num_dofs), g_v2(elt%num_dofs))
    allocate (b5(5, elt%num_dofs), dshape_s(3, elt%num_nodes))
    c5 = 0.0_rk
    c5(1, 1) = mu
    c5(2, 2) = 3.0_rk*mu
    c5(3, 3) = 4.0_rk*mu
    c5(4, 4) = 4.0_rk*mu
    c5(5, 5) = 4.0_rk*mu
    k_v2 = 0.0_rk
    g_v2 = 0.0_rk
    vol = 0.0_rk
    do i_qpt = 1_ik, elt%num_qpts
      call element_geometry_qpt(elt, coo, i_qpt, dshape_s, det_jac)
      b5 = 0.0_rk
      do a = 1_ik, elt%num_nodes
        i = 3_ik*(a - 1_ik)
        b5(1, i + 1) = dshape_s(1, a)
        b5(1, i + 2) = -dshape_s(2, a)
        b5(2, i + 1) = -dshape_s(1, a)/3.0_rk
        b5(2, i + 2) = -dshape_s(2, a)/3.0_rk
        b5(2, i + 3) = 2.0_rk*dshape_s(3, a)/3.0_rk
        b5(3, i + 1) = 0.5_rk*dshape_s(2, a)
        b5(3, i + 2) = 0.5_rk*dshape_s(1, a)
        b5(4, i + 1) = 0.5_rk*dshape_s(3, a)
        b5(4, i + 3) = 0.5_rk*dshape_s(1, a)
        b5(5, i + 2) = 0.5_rk*dshape_s(3, a)
        b5(5, i + 3) = 0.5_rk*dshape_s(2, a)
        g_v2(i + 1:i + 3) = g_v2(i + 1:i + 3) - elt%qpt_weight(i_qpt)*det_jac*dshape_s(:, a)
      end do
      k_v2 = k_v2 + elt%qpt_weight(i_qpt)*det_jac*matmul(transpose(b5), matmul(c5, b5))
      vol = vol + elt%qpt_weight(i_qpt)*det_jac
    end do

    ! The deviatoric part and the pressure vector against v2's.
    k_dev = k_e
    do i = 1_ik, elt%num_dofs
      k_dev(:, i) = k_dev(:, i) - r_e*g_e*g_e(i)
    end do
    scale = maxval(abs(k_v2))
    call assert_close(maxval(abs(k_dev - k_v2))/scale, 0.0_rk, name // ': deviatoric [k_e] ' // &
      & 'equals v2''s 5-vector form', abs_tol=1.0e-12_rk)
    call assert_close(maxval(abs(g_e - g_v2))/maxval(abs(g_v2)), 0.0_rk, name // &
      & ': {g_e} equals v2''s pcnst', abs_tol=1.0e-12_rk)
    call assert_close(r_e, mu/(6.0_rk*vol), name // ': r_e = mu / (6 V_e)', rel_tol=1.0e-13_rk)
    call assert_close(maxval(abs(k_e - transpose(k_e)))/maxval(abs(k_e)), 0.0_rk, name // &
      & ': [k_e] symmetric', abs_tol=1.0e-14_rk)

    ! Rigid motions: a translation and a rotation give zero forces.
    do a = 1_ik, elt%num_nodes
      vel(3*a - 2:3*a) = [0.3_rk, -1.1_rk, 0.7_rk]
    end do
    call assert_close(maxval(abs(matmul(k_e, vel)))/(scale*maxval(abs(vel))), 0.0_rk, name // &
      & ': [k_e] annihilates a translation', abs_tol=1.0e-12_rk)
    wmat = random_matrix()
    wmat = 0.5_rk*(wmat - transpose(wmat))
    do a = 1_ik, elt%num_nodes
      vel(3*a - 2:3*a) = matmul(wmat, coo(:, a))
    end do
    call assert_close(maxval(abs(matmul(k_e, vel)))/(scale*maxval(abs(vel))), 0.0_rk, name // &
      & ': [k_e] annihilates a rotation', abs_tol=1.0e-12_rk)

    ! Affine fields: {g_e}.{V_e} = - V_e tr(L); an isochoric one is not penalized.
    lmat = random_matrix()
    do a = 1_ik, elt%num_nodes
      vel(3*a - 2:3*a) = matmul(lmat, coo(:, a))
    end do
    call assert_close(dot_product(g_e, vel), -vol*(lmat(1, 1) + lmat(2, 2) + lmat(3, 3)), &
      & name // ': {g_e}.{V_e} = - V_e tr(L)', rel_tol=1.0e-12_rk, abs_tol=1.0e-14_rk)
    lmat(3, 3) = -lmat(1, 1) - lmat(2, 2)
    do a = 1_ik, elt%num_nodes
      vel(3*a - 2:3*a) = matmul(lmat, coo(:, a))
    end do
    call assert_close(abs(dot_product(g_e, vel))/(maxval(abs(g_e))*maxval(abs(vel))), 0.0_rk, &
      & name // ': an isochoric field is not penalized', abs_tol=1.0e-13_rk)

  end subroutine check_element

  !> The three iterations on the periodic box against their exact affine solution.
  !>
  !> Arguments:
  !>   file  (in)  character  mesh file
  !> Locals:
  !>   rate     real(rk)            imposed D11 = a [1/s]
  !>   mesh     type(mesh_t)        rank-local mesh
  !>   elt      type(element_t)     descriptor
  !>   cons     type(constraint_t)  periodicity relations (save: gfortran 15 -O3 otherwise
  !>                                warns of uninitialized descriptors, ROADMAP section G)
  !>   imposed  logical(9)          imposed components (11)
  !>   rates    real(rk)(9)         imposed rates [1/s]
  !>   vel      real(rk)(:)         velocities [mm/s]
  !>   bc_mask  logical(:)          no velocity conditions
  !>   fdot     real(rk)(3, 3)      Fbar-dot from the drive pairs [1/s]
  !>   expect   real(rk)(3, 3)      the exact Fbar-dot [1/s]
  !>   fluct    real(rk)(3)         v - Fbar-dot x at a node [mm/s]
  !>   lo, hi   real(rk)(3)         smallest and largest fluctuation components, this rank
  !>   lo_all, hi_all  real(rk)(3)  the same over the ranks
  !>   vmax, vmax_all  real(rk)     largest |v|, this rank and over the ranks [mm/s]
  !>   i        integer(ik)         node
  !>   k        integer             component
  subroutine check_periodic_solve(file)

    character(len=*), intent(in) :: file

    real(rk), parameter :: rate = 1.0e-3_rk
    type(mesh_t) :: mesh
    type(element_t) :: elt
    type(constraint_t), save :: cons
    logical :: imposed(9)
    real(rk) :: rates(9), fdot(3, 3), expect(3, 3), fluct(3), lo(3), hi(3), lo_all(3), &
      & hi_all(3), vmax, vmax_all
    real(rk), allocatable :: vel(:)
    logical, allocatable :: bc_mask(:)
    integer(ik) :: i
    integer :: k

    call mesh_load_blocks(mesh, file)
    call element_init(elt, mesh%elt_type)
    imposed = .false.
    imposed(1) = .true.
    rates = 0.0_rk
    rates(1) = rate
    call periodic_build(cons, mesh, elt, .true., imposed, .true.)
    call constraint_set_rates(cons, rates)

    ! The guess of the driver: zero, projected onto the relations (the offsets), zero mean.
    allocate (vel(3*mesh%num_nodes_loc), bc_mask(3*mesh%num_nodes_loc))
    vel = 0.0_rk
    bc_mask = .false.
    call constraint_project(cons, vel)
    call constraint_zero_mean(cons, mesh, vel)
    call isotropic_start_solve(mesh, elt, mesh%coo_s, bc_mask, 1.0e-13_rk, 16000, vel, cons)

    ! Fbar-dot: D11 = a, D22 = D33 = -a/18, no shear.
    fdot = constraint_macro_rate(cons, vel)
    expect = 0.0_rk
    expect(1, 1) = rate
    expect(2, 2) = -rate/18.0_rk
    expect(3, 3) = -rate/18.0_rk
    call assert_close(maxval(abs(fdot - expect))/rate, 0.0_rk, 'periodic box: Fbar-dot = ' // &
      & 'diag(a, -a/18, -a/18) after the three iterations', abs_tol=1.0e-8_rk)

    ! Affine: the fluctuation v - Fbar-dot x is a uniform translation (its spread over the
    ! nodes of all ranks below 1e-8 of the largest velocity).
    lo = huge(1.0_rk)
    hi = -huge(1.0_rk)
    vmax = 0.0_rk
    do i = 1_ik, mesh%num_nodes_owned
      fluct = vel(3*i - 2:3*i) - matmul(fdot, mesh%coo_s(:, i))
      lo = min(lo, fluct)
      hi = max(hi, fluct)
      vmax = max(vmax, maxval(abs(vel(3*i - 2:3*i))))
    end do
    do k = 1, 3
      call par_min(lo(k), lo_all(k))
      call par_max(hi(k), hi_all(k))
    end do
    call par_max(vmax, vmax_all)
    call assert_close(maxval(hi_all - lo_all)/vmax_all, 0.0_rk, 'periodic box: the solution is ' // &
      & 'affine (homogeneous isotropic material)', abs_tol=1.0e-8_rk)

  end subroutine check_periodic_solve

  !> A random curved element: an affine image of the reference element close to the identity,
  !> with the edge nodes moved by up to 3 % of the reference size (as test_element_matrix).
  !>
  !> Arguments:
  !>   elt  (in)  type(element_t)  descriptor
  !> Locals:
  !>   amat   real(rk)(3,3)  affine map
  !>   move   real(rk)(3)    edge-node displacement
  !>   a, first  integer     node, first edge node
  !>   scale  real(rk)       reference size
  function random_element(elt) result(coo)

    type(element_t), intent(in) :: elt
    real(rk), allocatable :: coo(:, :)

    real(rk) :: amat(3, 3), move(3), scale
    integer :: a, first

    do
      amat = 0.3_rk*random_matrix()
      amat(1, 1) = amat(1, 1) + 1.0_rk
      amat(2, 2) = amat(2, 2) + 1.0_rk
      amat(3, 3) = amat(3, 3) + 1.0_rk
      if (determinant3(amat) > 0.5_rk) exit
    end do
    if (elt%gmsh_type == tet10_gmsh_type) then
      coo = matmul(amat, tet10_node_coo_ref)
      first = 5
      scale = 1.0_rk
    else
      coo = matmul(amat, hex20_node_coo_ref)
      first = 9
      scale = 2.0_rk
    end if
    do a = first, int(elt%num_nodes)
      call random_number(move)
      coo(:, a) = coo(:, a) + 0.03_rk*scale*(2.0_rk*move - 1.0_rk)
    end do

  end function random_element

end program test_isotropic_start
