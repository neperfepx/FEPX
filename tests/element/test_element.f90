! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of element_mod: element descriptor, Jacobian and its inverse, sample-frame
!> shape-function derivatives, Mandel B matrix, element volume, face normals; patch tests.
!>
!> Source : FEPX (ROADMAP step 3.4).
!>
!> Usage: test_element <element fixture dir> <regression mesh dir>
!>        (tests/element/fixtures, tests/regression/cases/meshes).
!>
!> Checks:
!>   - descriptor: sizes, quadrature rules, shape functions and reference derivatives at every
!>     volume and face quadrature point equal the element module's (tet10, hex20);
!>   - reference: det J and d N/dx, dy, dz at the 15 quadrature points of three tet10 elements
!>     (scaled reference, affine, curved) equal FEPX's own (fixture geometry_reference.txt,
!>     written by dump_geometry.f90 in step 9.13; relative 1e-13; until step 9.13 v2-fixed's
!>     shape_3d_mod::sfder_hpar, to the same tolerance);
!>   - Jacobian on random curved elements (both types, 50 each): J from element_jacobian equals
!>     the explicit sum, det J equals an independent determinant, J J^-1 = I (1e-12); an
!>     inverted element (mirror image x -> -x) has det J < 0 and a zero inverse at every point;
!>   - chain rule: for a random nodal field f, sum_a d N_a/dx_i f_a dotted with the
!>     finite-difference change of x(xi) equals the finite-difference change of f(xi) along each
!>     reference direction (1e-7 relative);
!>   - B matrix against explicit 3x3 algebra: [B]{V} equals the Mandel vector of sym(grad v)
!>     with grad v(i, j) = sum_a V_a(i) d N_a/dx_j, random velocities (1e-13);
!>   - displacement patch test: nodal values of a linear field v = v0 + G x give grad v = G and
!>     [B]{V} = {sym G} at every quadrature point of curved elements (1e-12); of a quadratic
!>     field, the exact gradient at every quadrature point of straight-edged (affine) elements;
!>   - volume: element volumes add up to 1 mm^3 on the Neper meshes sx2.msh, poly_part1.msh
!>     (tet10) and hex_cube.msh (hex20) (1e-12); the volume of a curved tet10 equals
!>     (1/3) sum over its four faces of the face integral of x . n dA with the face elements
!>     (exact for both rules here; 1e-13); affine hex20 volume = |det A| 8;
!>   - quadrature in the sample frame: on an affine tet10, the volume integral of
!>     d/dx (x^(p+1) y^q z^r)/(p+1) = x^p y^q z^r (p + q + r <= 4) equals the surface integral
!>     over the faces (divergence theorem; both exact) (1e-12 relative);
!>   - force patch test (constant stress, poly_part1.msh): the assembled nodal forces
!>     sum_e integral B^T {sigma} dV equal, at every node, the face integrals of N_a sigma n dA
!>     over the fasets, i.e. zero at interior nodes and the consistent traction at boundary nodes
!>     (1e-12 of the largest force).
!>
!> Locals:
!>   element_dir, mesh_dir  character(512)  command-line directories
!>   tet, hex               type(element_t) descriptors
!>   tet_faces              integer(6, 4)   the faces of a tet10 (Gmsh order, outward normals)
program test_element

  use kinds_mod, only: rk, ik
  use element_mod, only: element_t, element_init, element_jacobian, jacobian_det_inv, &
    & shape_deriv_sample, element_geometry_qpt, mandel_b_matrix, element_volume, &
    & face_area_normal
  use tet10_mod, only: tet10_gmsh_type, tet10_node_coo_ref, tet10_qpt_coo, &
    & tet10_qpt_weight, tet10_shape, tet10_shape_deriv, tet10_face_qpt_coo, &
    & tet10_face_qpt_weight, tet10_face_shape, tet10_face_shape_deriv
  use hex20_mod, only: hex20_gmsh_type, hex20_node_coo_ref, hex20_qpt_coo, &
    & hex20_qpt_weight, hex20_shape, hex20_shape_deriv, hex20_face_qpt_coo, &
    & hex20_face_shape, hex20_face_shape_deriv
  use tensor_mod, only: mandel_from_matrix
  use mesh_file_mod, only: mesh_file_t, mesh_file_load
  use tensor_test_utils_mod, only: seed_random, random_matrix, random_symmetric, determinant3
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary

  implicit none

  character(len=512) :: element_dir, mesh_dir
  type(element_t) :: tet, hex
  integer(ik), parameter :: tet_faces(6, 4) = reshape([1_ik, 3_ik, 2_ik, 7_ik, 6_ik, 5_ik, &
    & 1_ik, 2_ik, 4_ik, 5_ik, 10_ik, 8_ik, 2_ik, 3_ik, 4_ik, 6_ik, 9_ik, 10_ik, &
    & 1_ik, 4_ik, 3_ik, 8_ik, 9_ik, 7_ik], [6, 4])

  call get_command_argument(1, element_dir)
  call get_command_argument(2, mesh_dir)
  call seed_random(20261009)

  call element_init(tet, tet10_gmsh_type)
  call element_init(hex, hex20_gmsh_type)

  call check_descriptor()
  call check_reference(trim(element_dir) // '/geometry_reference.txt')
  call check_jacobian(tet, 'tet10')
  call check_jacobian(hex, 'hex20')
  call check_chain_rule(tet, 'tet10')
  call check_chain_rule(hex, 'hex20')
  call check_b_matrix(tet, 'tet10')
  call check_b_matrix(hex, 'hex20')
  call check_patch_displacement(tet, 'tet10')
  call check_patch_displacement(hex, 'hex20')
  call check_volumes()
  call check_divergence_quadrature()
  call check_mesh_volume(trim(mesh_dir) // '/sx2.msh', tet, 'sx2')
  call check_mesh_volume(trim(mesh_dir) // '/poly_part1.msh', tet, 'poly')
  call check_mesh_volume(trim(element_dir) // '/hex_cube.msh', hex, 'hex_cube')
  call check_patch_force(trim(mesh_dir) // '/poly_part1.msh')

  call test_summary()

contains

  !> The descriptors against the element modules' tables.
  !>
  !> Locals:
  !>   q  integer  quadrature point
  subroutine check_descriptor()

    integer :: q

    ! tet10.
    call assert_equal(tet%name, 'tet10', 'tet name')
    call assert_equal(tet%num_nodes, 10_ik, 'tet num_nodes')
    call assert_equal(tet%num_dofs, 30_ik, 'tet num_dofs')
    call assert_equal(tet%num_qpts, 15_ik, 'tet num_qpts')
    call assert_close(tet%qpt_coo(:, tet%qpt_centroid), [0.25_rk, 0.25_rk, 0.25_rk], &
      & 'tet output point at the centroid', abs_tol=1.0e-15_rk, rel_tol=0.0_rk)
    call assert_close(hex%qpt_coo(:, hex%qpt_centroid), [0.0_rk, 0.0_rk, 0.0_rk], &
      & 'hex output point at the centroid', abs_tol=1.0e-15_rk, rel_tol=0.0_rk)
    call assert_equal(tet%face_num_nodes, 6_ik, 'tet face_num_nodes')
    call assert_equal(tet%face_num_qpts, 7_ik, 'tet face_num_qpts')
    call assert_close(reshape(tet%qpt_coo, [45]), reshape(tet10_qpt_coo, [45]), 'tet qpt_coo')
    call assert_close(tet%qpt_weight, tet10_qpt_weight, 'tet qpt_weight')
    call assert_close(tet%face_qpt_weight, tet10_face_qpt_weight, 'tet face_qpt_weight')
    do q = 1, 15
      call assert_close(tet%shape_qpt(:, q), tet10_shape(tet10_qpt_coo(:, q)), 'tet shape_qpt')
      call assert_close(reshape(tet%dshape_qpt(:, :, q), [30]), &
        & reshape(tet10_shape_deriv(tet10_qpt_coo(:, q)), [30]), 'tet dshape_qpt')
    end do
    do q = 1, 7
      call assert_close(tet%face_shape_qpt(:, q), tet10_face_shape(tet10_face_qpt_coo(:, q)), &
        & 'tet face_shape_qpt')
      call assert_close(reshape(tet%face_dshape_qpt(:, :, q), [12]), &
        & reshape(tet10_face_shape_deriv(tet10_face_qpt_coo(:, q)), [12]), &
        & 'tet face_dshape_qpt')
    end do

    ! hex20.
    call assert_equal(hex%name, 'hex20', 'hex name')
    call assert_equal(hex%num_nodes, 20_ik, 'hex num_nodes')
    call assert_equal(hex%num_dofs, 60_ik, 'hex num_dofs')
    call assert_equal(hex%num_qpts, 27_ik, 'hex num_qpts')
    call assert_equal(hex%face_num_nodes, 8_ik, 'hex face_num_nodes')
    call assert_equal(hex%face_num_qpts, 9_ik, 'hex face_num_qpts')
    call assert_close(hex%qpt_weight, hex20_qpt_weight, 'hex qpt_weight')
    do q = 1, 27
      call assert_close(hex%shape_qpt(:, q), hex20_shape(hex20_qpt_coo(:, q)), 'hex shape_qpt')
      call assert_close(reshape(hex%dshape_qpt(:, :, q), [60]), &
        & reshape(hex20_shape_deriv(hex20_qpt_coo(:, q)), [60]), 'hex dshape_qpt')
    end do
    do q = 1, 9
      call assert_close(hex%face_shape_qpt(:, q), hex20_face_shape(hex20_face_qpt_coo(:, q)), &
        & 'hex face_shape_qpt')
      call assert_close(reshape(hex%face_dshape_qpt(:, :, q), [16]), &
        & reshape(hex20_face_shape_deriv(hex20_face_qpt_coo(:, q)), [16]), &
        & 'hex face_dshape_qpt')
    end do

  end subroutine check_descriptor

  !> det J and sample-frame derivatives against the reference fixture.
  !>
  !> Arguments:
  !>   path  (in)  character  fixture file
  !> Locals:
  !>   unit, ios, num, e, q, a  integer   file unit, status, counts and indices
  !>   word        character(16)          section name
  !>   coo         real(rk)(3,10,3)       element coordinates [mm]
  !>   line        real(rk)(31)           det J and the 30 derivatives of the reference
  !>   dshape_s    real(rk)(3,10)         FEPX derivatives [1/mm]
  !>   det_jac     real(rk)               FEPX det J [mm^3]
  subroutine check_reference(path)

    character(len=*), intent(in) :: path

    integer :: unit, ios, num, e, q, a
    character(len=16) :: word
    real(rk) :: coo(3, 10, 3), line(31), dshape_s(3, 10), det_jac

    open (newunit=unit, file=path, status='old', action='read', iostat=ios)
    call assert_equal(ios, 0, 'reference geometry fixture open')
    if (ios /= 0) return

    ! Element coordinates.
    read (unit, *) word, num
    call assert_equal(num, 3, 'reference geometry elements')
    do e = 1, 3
      do a = 1, 10
        read (unit, *) coo(:, a, e)
      end do
    end do

    ! det J and derivatives at each quadrature point.
    read (unit, *) word, num
    call assert_equal(num, 15, 'reference geometry qpts')
    do e = 1, 3
      do q = 1, 15
        read (unit, *) line
        call element_geometry_qpt(tet, coo(:, :, e), int(q, ik), dshape_s, det_jac)
        call assert_close(det_jac, line(1), 'reference det J', rel_tol=1e-13_rk)
        call assert_close(dshape_s(1, :), line(2:11), 'reference dN/dx', abs_tol=1e-13_rk*maxval( &
          & abs(line(2:31))))
        call assert_close(dshape_s(2, :), line(12:21), 'reference dN/dy', abs_tol=1e-13_rk*maxval( &
          & abs(line(2:31))))
        call assert_close(dshape_s(3, :), line(22:31), 'reference dN/dz', abs_tol=1e-13_rk*maxval( &
          & abs(line(2:31))))
      end do
    end do
    close (unit)

  end subroutine check_reference

  !> Jacobian, determinant and inverse on random curved elements; an inverted element.
  !>
  !> Arguments:
  !>   elt    (in)  type(element_t)  descriptor
  !>   label  (in)  character        check label prefix
  !> Locals:
  !>   i_elt, q, a     integer         loop indices
  !>   coo             real(rk)(:,:)   element coordinates [mm]
  !>   jac, jac_inv    real(rk)(3,3)   Jacobian and inverse
  !>   jac_sum         real(rk)(3,3)   explicit sum over the nodes
  !>   det_jac         real(rk)        det J
  !>   ident           real(rk)(3,3)   identity
  subroutine check_jacobian(elt, label)

    type(element_t), intent(in) :: elt
    character(len=*), intent(in) :: label

    integer :: i_elt, q, a
    real(rk), allocatable :: coo(:, :)
    real(rk) :: jac(3, 3), jac_inv(3, 3), jac_sum(3, 3), det_jac, ident(3, 3)

    ident = 0.0_rk
    ident(1, 1) = 1.0_rk
    ident(2, 2) = 1.0_rk
    ident(3, 3) = 1.0_rk
    do i_elt = 1, 50
      coo = random_element(elt)
      do q = 1, int(elt%num_qpts)
        jac = element_jacobian(coo, elt%dshape_qpt(:, :, q))
        ! Explicit sum J(i, j) = sum_a x_a(i) dN_a/dxi_j.
        jac_sum = 0.0_rk
        do a = 1, int(elt%num_nodes)
          jac_sum = jac_sum + spread(coo(:, a), 2, 3)*spread(elt%dshape_qpt(:, a, q), 1, 3)
        end do
        call assert_close(reshape(jac, [9]), reshape(jac_sum, [9]), label // ' J sum', &
          & abs_tol=1e-14_rk)
        call jacobian_det_inv(jac, det_jac, jac_inv)
        call assert_close(det_jac, determinant3(jac), label // ' det J', rel_tol=1e-13_rk)
        call assert_true(det_jac > 0.0_rk, label // ' random element valid')
        call assert_close(reshape(matmul(jac, jac_inv), [9]), reshape(ident, [9]), &
          & label // ' J J^-1 = I', abs_tol=1e-12_rk)
        call assert_close(reshape(matmul(jac_inv, jac), [9]), reshape(ident, [9]), &
          & label // ' J^-1 J = I', abs_tol=1e-12_rk)
      end do
    end do

    ! Inverted element: the mirror image x -> -x of a valid element (det J < 0 everywhere).
    coo = random_element(elt)
    coo(1, :) = -coo(1, :)
    do q = 1, int(elt%num_qpts)
      call jacobian_det_inv(element_jacobian(coo, elt%dshape_qpt(:, :, q)), det_jac, jac_inv)
      call assert_true(det_jac < 0.0_rk, label // ' inverted element det J < 0')
      call assert_true(all(abs(jac_inv) < tiny(1.0_rk)), label // &
        & ' inverted element zero inverse')
    end do

  end subroutine check_jacobian

  !> Chain rule against finite differences of the isoparametric map.
  !>
  !> Arguments:
  !>   elt    (in)  type(element_t)  descriptor
  !>   label  (in)  character        check label prefix
  !> Locals:
  !>   i_elt, i_pt, j   integer          loop indices
  !>   coo              real(rk)(:,:)    element coordinates [mm]
  !>   f                real(rk)(:)      nodal field values
  !>   xi, xi_p, xi_m   real(rk)(3)      point and shifted points
  !>   dshape_s         real(rk)(:,:)    sample-frame derivatives
  !>   grad             real(rk)(3)      grad f at xi
  !>   dx               real(rk)(3)      central difference of x(xi)
  !>   df               real(rk)         central difference of f(xi)
  !>   jac, jac_inv     real(rk)(3,3)    Jacobian and inverse
  !>   det_jac, h       real(rk)         det J, step
  subroutine check_chain_rule(elt, label)

    type(element_t), intent(in) :: elt
    character(len=*), intent(in) :: label

    integer :: i_elt, i_pt, j
    real(rk), allocatable :: coo(:, :), f(:), dshape_s(:, :)
    real(rk) :: xi(3), xi_p(3), xi_m(3), grad(3), dx(3), df, jac(3, 3), jac_inv(3, 3)
    real(rk) :: det_jac, h

    h = 1e-5_rk
    allocate (f(elt%num_nodes))
    do i_elt = 1, 20
      coo = random_element(elt)
      call random_number(f)
      do i_pt = 1, 10
        xi = random_ref_point(elt)
        jac = element_jacobian(coo, ref_deriv(elt, xi))
        call jacobian_det_inv(jac, det_jac, jac_inv)
        dshape_s = shape_deriv_sample(ref_deriv(elt, xi), jac_inv)
        grad = matmul(dshape_s, f)
        do j = 1, 3
          xi_p = xi
          xi_p(j) = xi(j) + h
          xi_m = xi
          xi_m(j) = xi(j) - h
          dx = matmul(coo, ref_shape(elt, xi_p) - ref_shape(elt, xi_m))
          df = dot_product(f, ref_shape(elt, xi_p) - ref_shape(elt, xi_m))
          call assert_close(dot_product(grad, dx), df, label // ' chain rule vs FD', &
            & abs_tol=1e-7_rk*max(abs(df), 1e-3_rk))
        end do
      end do
    end do

  end subroutine check_chain_rule

  !> Mandel B matrix against explicit 3x3 algebra.
  !>
  !> Arguments:
  !>   elt    (in)  type(element_t)  descriptor
  !>   label  (in)  character        check label prefix
  !> Locals:
  !>   i_elt, q, a     integer         loop indices
  !>   coo             real(rk)(:,:)   element coordinates [mm]
  !>   vel             real(rk)(:,:)   (3, num_nodes) nodal velocities
  !>   dshape_s        real(rk)(:,:)   sample-frame derivatives
  !>   b               real(rk)(:,:)   B matrix
  !>   velgrad         real(rk)(3,3)   grad v(i, j) = sum_a V_a(i) dN_a/dx_j
  !>   det_jac         real(rk)        det J
  subroutine check_b_matrix(elt, label)

    type(element_t), intent(in) :: elt
    character(len=*), intent(in) :: label

    integer :: i_elt, q
    real(rk), allocatable :: coo(:, :), vel(:, :), dshape_s(:, :), b(:, :)
    real(rk) :: velgrad(3, 3), det_jac

    allocate (vel(3, elt%num_nodes), dshape_s(3, elt%num_nodes))
    do i_elt = 1, 20
      coo = random_element(elt)
      call random_number(vel)
      vel = 2.0_rk*vel - 1.0_rk
      do q = 1, int(elt%num_qpts)
        call element_geometry_qpt(elt, coo, int(q, ik), dshape_s, det_jac)
        if (allocated(b)) deallocate (b)
        allocate (b(6, 3*size(dshape_s, 2)))
        call mandel_b_matrix(dshape_s, b)
        call assert_equal(size(b, 2), int(elt%num_dofs), label // ' B columns')
        velgrad = matmul(vel, transpose(dshape_s))
        call assert_close(matmul(b, reshape(vel, [elt%num_dofs])), mandel_from_matrix(velgrad), &
          & label // ' B V = {sym grad v}', abs_tol=1e-13_rk*maxval(abs(velgrad)))
      end do
    end do

  end subroutine check_b_matrix

  !> Displacement patch test: linear fields exact on curved elements, quadratic fields exact on
  !> straight-edged elements.
  !>
  !> Arguments:
  !>   elt    (in)  type(element_t)  descriptor
  !>   label  (in)  character        check label prefix
  !> Locals:
  !>   i_elt, q, a     integer         loop indices
  !>   coo             real(rk)(:,:)   element coordinates [mm]
  !>   vel             real(rk)(:,:)   nodal values of the field
  !>   dshape_s        real(rk)(:,:)   sample-frame derivatives
  !>   grad_g          real(rk)(3,3)   linear field gradient G
  !>   v0              real(rk)(3)     linear field offset
  !>   hess            real(rk)(3,3,3) quadratic field: v_i = sum H(i,j,k) x_j x_k
  !>   xq              real(rk)(3)     sample coordinates of a quadrature point
  !>   exact           real(rk)(3,3)   exact gradient of the quadratic field at xq
  !>   det_jac         real(rk)        det J
  subroutine check_patch_displacement(elt, label)

    type(element_t), intent(in) :: elt
    character(len=*), intent(in) :: label

    integer :: i_elt, q, a, i, j
    real(rk), allocatable :: coo(:, :), vel(:, :), dshape_s(:, :), bmat(:, :)
    real(rk) :: grad_g(3, 3), v0(3), hess(3, 3, 3), xq(3), exact(3, 3), det_jac

    allocate (vel(3, elt%num_nodes), dshape_s(3, elt%num_nodes))
    do i_elt = 1, 20
      ! Linear field on a curved element.
      coo = random_element(elt)
      grad_g = random_matrix()
      call random_number(v0)
      do a = 1, int(elt%num_nodes)
        vel(:, a) = v0 + matmul(grad_g, coo(:, a))
      end do
      do q = 1, int(elt%num_qpts)
        call element_geometry_qpt(elt, coo, int(q, ik), dshape_s, det_jac)
        call assert_close(reshape(matmul(vel, transpose(dshape_s)), [9]), &
          & reshape(grad_g, [9]), label // ' linear field: grad exact', abs_tol=1e-12_rk)
        if (allocated(bmat)) deallocate (bmat)
        allocate (bmat(6, 3*size(dshape_s, 2)))
        call mandel_b_matrix(dshape_s, bmat)
        call assert_close(matmul(bmat, reshape(vel, [elt%num_dofs])), &
          & mandel_from_matrix(grad_g), label // ' linear field: B V = {sym G}', &
          & abs_tol=1e-12_rk)
      end do

      ! Quadratic field on a straight-edged (affine) element.
      coo = affine_element(elt)
      call random_number(hess)
      hess = 2.0_rk*hess - 1.0_rk
      do a = 1, int(elt%num_nodes)
        do i = 1, 3
          vel(i, a) = dot_product(coo(:, a), matmul(hess(i, :, :), coo(:, a)))
        end do
      end do
      do q = 1, int(elt%num_qpts)
        call element_geometry_qpt(elt, coo, int(q, ik), dshape_s, det_jac)
        xq = matmul(coo, elt%shape_qpt(:, q))
        do i = 1, 3
          do j = 1, 3
            exact(i, j) = dot_product(hess(i, j, :) + hess(i, :, j), xq)
          end do
        end do
        call assert_close(reshape(matmul(vel, transpose(dshape_s)), [9]), reshape(exact, [9]), &
          & label // ' quadratic field (affine element): grad exact', abs_tol=1e-11_rk)
      end do
    end do

  end subroutine check_patch_displacement

  !> Element volumes: curved tet10 against the surface integral (1/3) x . n dA over its faces;
  !> affine hex20 against |det A| 8.
  !>
  !> Locals:
  !>   i_elt, f, q      integer          loop indices
  !>   coo              real(rk)(:,:)    element coordinates [mm]
  !>   face_coo         real(rk)(3,6)    face node coordinates [mm]
  !>   vol_surf         real(rk)         (1/3) surface integral of x . n [mm^3]
  !>   xq, normal       real(rk)(3)      face point and area-weighted normal
  !>   amat             real(rk)(3,3)    affine map
  subroutine check_volumes()

    integer :: i_elt, f, q
    real(rk), allocatable :: coo(:, :)
    real(rk) :: face_coo(3, 6), vol_surf, xq(3), normal(3), amat(3, 3)

    ! Curved tet10: divergence theorem with the face elements (both rules exact here).
    do i_elt = 1, 50
      coo = random_element(tet)
      vol_surf = 0.0_rk
      do f = 1, 4
        face_coo = coo(:, tet_faces(:, f))
        do q = 1, int(tet%face_num_qpts)
          xq = matmul(face_coo, tet%face_shape_qpt(:, q))
          normal = face_area_normal(face_coo, tet%face_dshape_qpt(:, :, q))
          vol_surf = vol_surf + tet%face_qpt_weight(q)*dot_product(xq, normal)/3.0_rk
        end do
      end do
      call assert_close(element_volume(tet, coo), vol_surf, 'tet10 curved volume = surface', &
        & rel_tol=1e-13_rk)
    end do

    ! Affine hex20: volume |det A| times the reference volume 8.
    do i_elt = 1, 20
      amat = random_affine()
      coo = matmul(amat, hex20_node_coo_ref)
      call assert_close(element_volume(hex, coo), 8.0_rk*abs(determinant3(amat)), &
        & 'hex20 affine volume', rel_tol=1e-13_rk)
    end do

  end subroutine check_volumes

  !> Volume quadrature in the sample frame against the divergence theorem on an affine tet10:
  !> integral of x^p y^q z^r dV = surface integral of (x^(p+1) y^q z^r / (p+1)) n_x dA.
  !>
  !> Locals:
  !>   p, q, r, f, iq   integer          exponents and loop indices
  !>   coo              real(rk)(:,:)    element coordinates [mm]
  !>   face_coo         real(rk)(3,6)    face node coordinates [mm]
  !>   vol_int, surf_int real(rk)        the two integrals
  !>   xq, normal       real(rk)(3)      point, area-weighted normal
  !>   det_jac          real(rk)         det J
  !>   jac_inv          real(rk)(3,3)    unused inverse
  subroutine check_divergence_quadrature()

    integer :: p, q, r, f, iq
    real(rk), allocatable :: coo(:, :)
    real(rk) :: face_coo(3, 6), vol_int, surf_int, xq(3), normal(3), det_jac, jac_inv(3, 3)

    coo = affine_element(tet)
    do p = 0, 4
      do q = 0, 4 - p
        do r = 0, 4 - p - q
          ! Volume integral with the 15-point rule.
          vol_int = 0.0_rk
          do iq = 1, int(tet%num_qpts)
            xq = matmul(coo, tet%shape_qpt(:, iq))
            call jacobian_det_inv(element_jacobian(coo, tet%dshape_qpt(:, :, iq)), det_jac, &
              & jac_inv)
            vol_int = vol_int + tet%qpt_weight(iq)*det_jac*xq(1)**p*xq(2)**q*xq(3)**r
          end do
          ! Surface integral with the 7-point face rule over the four faces.
          surf_int = 0.0_rk
          do f = 1, 4
            face_coo = coo(:, tet_faces(:, f))
            do iq = 1, int(tet%face_num_qpts)
              xq = matmul(face_coo, tet%face_shape_qpt(:, iq))
              normal = face_area_normal(face_coo, tet%face_dshape_qpt(:, :, iq))
              surf_int = surf_int + tet%face_qpt_weight(iq)*normal(1)*xq(1)**(p + 1)* &
                & xq(2)**q*xq(3)**r/real(p + 1, rk)
            end do
          end do
          call assert_close(vol_int, surf_int, 'tet10 volume quadrature = divergence theorem', &
            & abs_tol=1e-12_rk*max(abs(surf_int), 1.0_rk))
        end do
      end do
    end do

  end subroutine check_divergence_quadrature

  !> Element volumes of a Neper mesh of the unit cube add up to 1 mm^3; det J > 0 everywhere.
  !>
  !> Arguments:
  !>   path   (in)  character        mesh file
  !>   elt    (in)  type(element_t)  descriptor of the mesh's element type
  !>   label  (in)  character        check label prefix
  !> Locals:
  !>   mesh      type(mesh_file_t)   parsed mesh
  !>   i_elt, q  integer             loop indices
  !>   vol       real(rk)            total volume [mm^3]
  !>   min_det   real(rk)            smallest det J
  !>   dshape_s  real(rk)(:,:)       sample-frame derivatives
  !>   det_jac   real(rk)            det J
  subroutine check_mesh_volume(path, elt, label)

    character(len=*), intent(in) :: path
    type(element_t), intent(in) :: elt
    character(len=*), intent(in) :: label

    type(mesh_file_t) :: mesh
    integer :: i_elt, q
    real(rk) :: vol, min_det, det_jac
    real(rk), allocatable :: dshape_s(:, :)

    call mesh_file_load(mesh, path)
    call assert_equal(mesh%elt_type, elt%gmsh_type, label // ' element type')
    allocate (dshape_s(3, elt%num_nodes))
    vol = 0.0_rk
    min_det = huge(1.0_rk)
    do i_elt = 1, int(mesh%num_elts)
      vol = vol + element_volume(elt, mesh%coo(:, mesh%elt_nodes(:, i_elt)))
      do q = 1, int(elt%num_qpts)
        call element_geometry_qpt(elt, mesh%coo(:, mesh%elt_nodes(:, i_elt)), int(q, ik), &
          & dshape_s, det_jac)
        min_det = min(min_det, det_jac)
      end do
    end do
    call assert_close(vol, 1.0_rk, label // ' mesh volume 1', rel_tol=1e-12_rk)
    call assert_true(min_det > 0.0_rk, label // ' det J > 0 everywhere')

  end subroutine check_mesh_volume

  !> Force patch test: constant stress; assembled volume forces sum_e integral B^T {sigma} dV
  !> equal the face integrals of N_a sigma n dA over the fasets, at every node.
  !>
  !> Arguments:
  !>   path  (in)  character  mesh file (tet10, with $Fasets covering the whole boundary)
  !> Locals:
  !>   mesh             type(mesh_file_t)  parsed mesh
  !>   sigma            real(rk)(3,3)      constant stress [MPa]
  !>   sigma_m          real(rk)(6)        its Mandel vector [MPa]
  !>   force_vol        real(rk)(:,:)      (3, num_nodes) assembled volume forces [N]
  !>   force_surf       real(rk)(:,:)      (3, num_nodes) assembled face tractions [N]
  !>   i_elt, q, a, i_set, i_face  integer  loop indices
  !>   coo              real(rk)(:,:)      element coordinates [mm]
  !>   dshape_s         real(rk)(:,:)      sample-frame derivatives
  !>   fe               real(rk)(:)        element force vector
  !>   face_coo         real(rk)(3,6)      face coordinates [mm]
  !>   normal           real(rk)(3)        area-weighted normal [mm^2]
  !>   det_jac          real(rk)           det J
  subroutine check_patch_force(path)

    character(len=*), intent(in) :: path

    type(mesh_file_t) :: mesh
    real(rk) :: sigma(3, 3), sigma_m(6), face_coo(3, 6), normal(3), det_jac
    real(rk), allocatable :: force_vol(:, :), force_surf(:, :), coo(:, :), dshape_s(:, :), fe(:)
    real(rk), allocatable :: bmat(:, :)
    integer :: i_elt, q, a, i_set, i_face

    call mesh_file_load(mesh, path)
    sigma = 100.0_rk*random_symmetric()
    sigma_m = mandel_from_matrix(sigma)
    allocate (force_vol(3, mesh%num_nodes), force_surf(3, mesh%num_nodes))
    allocate (dshape_s(3, tet%num_nodes), fe(tet%num_dofs))
    force_vol = 0.0_rk
    force_surf = 0.0_rk

    ! Volume: element force vectors integral B^T {sigma} dV, assembled.
    do i_elt = 1, int(mesh%num_elts)
      coo = mesh%coo(:, mesh%elt_nodes(:, i_elt))
      fe = 0.0_rk
      do q = 1, int(tet%num_qpts)
        call element_geometry_qpt(tet, coo, int(q, ik), dshape_s, det_jac)
        if (.not. allocated(bmat)) allocate (bmat(6, 3*size(dshape_s, 2)))
        call mandel_b_matrix(dshape_s, bmat)
        fe = fe + tet%qpt_weight(q)*det_jac*matmul(transpose(bmat), sigma_m)
      end do
      do a = 1, int(tet%num_nodes)
        force_vol(:, mesh%elt_nodes(a, i_elt)) = force_vol(:, mesh%elt_nodes(a, i_elt)) + &
          & fe(3*a - 2:3*a)
      end do
    end do

    ! Surface: face integrals of N_a sigma n dA over every faset face, assembled.
    do i_set = 1, size(mesh%fasets)
      do i_face = 1, size(mesh%fasets(i_set)%elts)
        face_coo = mesh%coo(:, mesh%fasets(i_set)%nodes(:, i_face))
        do q = 1, int(tet%face_num_qpts)
          normal = face_area_normal(face_coo, tet%face_dshape_qpt(:, :, q))
          do a = 1, int(tet%face_num_nodes)
            force_surf(:, mesh%fasets(i_set)%nodes(a, i_face)) = &
              & force_surf(:, mesh%fasets(i_set)%nodes(a, i_face)) + &
              & tet%face_qpt_weight(q)*tet%face_shape_qpt(a, q)*matmul(sigma, normal)
          end do
        end do
      end do
    end do

    call assert_close(reshape(force_vol, [3*mesh%num_nodes]), &
      & reshape(force_surf, [3*mesh%num_nodes]), 'patch test: volume forces = face tractions', &
      & abs_tol=1e-12_rk*maxval(abs(force_surf)))
    call assert_close(sum(force_vol, dim=2), [0.0_rk, 0.0_rk, 0.0_rk], &
      & 'patch test: forces in equilibrium', abs_tol=1e-12_rk*maxval(abs(force_surf)))

  end subroutine check_patch_force

  !> A random affine map close to the identity (det > 0.5).
  !>
  !> Locals:
  !>   amat  real(rk)(3,3)  the map
  function random_affine() result(amat)

    real(rk) :: amat(3, 3)

    do
      amat = 0.3_rk*random_matrix()
      amat(1, 1) = amat(1, 1) + 1.0_rk
      amat(2, 2) = amat(2, 2) + 1.0_rk
      amat(3, 3) = amat(3, 3) + 1.0_rk
      if (determinant3(amat) > 0.5_rk) exit
    end do

  end function random_affine

  !> A random straight-edged element: affine image of the reference element.
  !>
  !> Arguments:
  !>   elt  (in)  type(element_t)  descriptor
  !> Locals:
  !>   amat   real(rk)(3,3)  affine map
  !>   shift  real(rk)(3)    translation
  !>   a      integer        node
  function affine_element(elt) result(coo)

    type(element_t), intent(in) :: elt
    real(rk), allocatable :: coo(:, :)

    real(rk) :: amat(3, 3), shift(3)
    integer :: a

    amat = random_affine()
    call random_number(shift)
    allocate (coo(3, elt%num_nodes))
    do a = 1, int(elt%num_nodes)
      coo(:, a) = matmul(amat, ref_node(elt, a)) + shift
    end do

  end function affine_element

  !> A random curved element: an affine element whose edge nodes are moved by up to 3 % of the
  !> reference size (det J stays positive; checked by the callers).
  !>
  !> Arguments:
  !>   elt  (in)  type(element_t)  descriptor
  !> Locals:
  !>   move   real(rk)(3)    displacement of an edge node
  !>   a      integer        node
  !>   first  integer        first edge node (5 or 9)
  !>   scale  real(rk)       reference size (1 or 2)
  function random_element(elt) result(coo)

    type(element_t), intent(in) :: elt
    real(rk), allocatable :: coo(:, :)

    real(rk) :: move(3), scale
    integer :: a, first

    coo = affine_element(elt)
    first = merge(5, 9, elt%gmsh_type == tet10_gmsh_type)
    scale = merge(1.0_rk, 2.0_rk, elt%gmsh_type == tet10_gmsh_type)
    do a = first, int(elt%num_nodes)
      call random_number(move)
      coo(:, a) = coo(:, a) + 0.03_rk*scale*(2.0_rk*move - 1.0_rk)
    end do

  end function random_element

  !> Reference coordinates of node a.
  !>
  !> Arguments:
  !>   elt  (in)  type(element_t)  descriptor
  !>   a    (in)  integer          node
  function ref_node(elt, a) result(xi)

    type(element_t), intent(in) :: elt
    integer, intent(in) :: a
    real(rk) :: xi(3)

    if (elt%gmsh_type == tet10_gmsh_type) then
      xi = tet10_node_coo_ref(:, a)
    else
      xi = hex20_node_coo_ref(:, a)
    end if

  end function ref_node

  !> A random point of the reference element (tet: by sorting; hex: uniform in [-1, 1]^3).
  !>
  !> Arguments:
  !>   elt  (in)  type(element_t)  descriptor
  !> Locals:
  !>   u  real(rk)(3)  uniform numbers
  function random_ref_point(elt) result(xi)

    type(element_t), intent(in) :: elt
    real(rk) :: xi(3)

    real(rk) :: u(3)

    call random_number(u)
    if (elt%gmsh_type == tet10_gmsh_type) then
      u = [minval(u), u(1) + u(2) + u(3) - minval(u) - maxval(u), maxval(u)]
      xi = [u(1), u(2) - u(1), u(3) - u(2)]
    else
      xi = 2.0_rk*u - 1.0_rk
    end if

  end function random_ref_point

  !> Shape functions of the element type at a reference point.
  !>
  !> Arguments:
  !>   elt  (in)  type(element_t)  descriptor
  !>   xi   (in)  real(rk)(3)      reference point
  function ref_shape(elt, xi) result(shape)

    type(element_t), intent(in) :: elt
    real(rk), intent(in) :: xi(3)
    real(rk) :: shape(elt%num_nodes)

    if (elt%gmsh_type == tet10_gmsh_type) then
      shape = tet10_shape(xi)
    else
      shape = hex20_shape(xi)
    end if

  end function ref_shape

  !> Reference derivatives of the element type at a reference point.
  !>
  !> Arguments:
  !>   elt  (in)  type(element_t)  descriptor
  !>   xi   (in)  real(rk)(3)      reference point
  function ref_deriv(elt, xi) result(dshape)

    type(element_t), intent(in) :: elt
    real(rk), intent(in) :: xi(3)
    real(rk) :: dshape(3, elt%num_nodes)

    if (elt%gmsh_type == tet10_gmsh_type) then
      dshape = tet10_shape_deriv(xi)
    else
      dshape = hex20_shape_deriv(xi)
    end if

  end function ref_deriv

end program test_element
