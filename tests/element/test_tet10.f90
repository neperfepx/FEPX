! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of tet10_mod (10-node tetrahedron and its 15-point quadrature rule).
!>
!> Source : FEPX (ROADMAP step 3.1).
!>
!> Usage: test_tet10 <element fixture dir> <io fixture dir> <regression mesh dir>
!>        (tests/element/fixtures, tests/io/fixtures, tests/regression/cases/meshes).
!>
!> Quadrature (15 points):
!>   - reference: points and weights equal FEPX's own (fixture tet10_reference.txt, written by
!>     dump_tet10.f90 in step 9.13; README.md) exactly (until step 9.13 the fixture was v2-fixed's
!>     quadrature_init values, also equal exactly);
!>   - weights positive, sum 1/6 (the reference volume); points in the closed reference element;
!>   - exactness: every monomial xi1^p xi2^q xi3^r with p + q + r <= 5 is integrated to
!>     p! q! r! / (p + q + r + 3)! (relative 1e-14); degree exactly 5: some degree-6 monomial is
!>     not (relative error > 1e-6).
!> Shape functions:
!>   - reference: values and derivatives equal the fixture's at the 15 quadrature points, the
!>     10 nodes and 5 other points (1e-15; until step 9.13 v2-fixed's t10_shape_hpar /
!>     t10_deriv_hpar after v2's node renumbering, also to 1e-15);
!>   - N_a(node b) = delta_ab at the reference node coordinates;
!>   - at 200 random points of the element: partition of unity (sum N = 1, sum dN = 0) and
!>     quadratic completeness: sum_a N_a q(xi_a) = q(xi) and sum_a dN_a q(xi_a) = grad q for the
!>     10 monomials of degree <= 2 (1e-13); derivatives against central finite differences
!>     (exact for quadratics up to rounding; 1e-9);
!>   - the corners are the nodes with integer reference coordinates (1-4) and every other node
!>     is the midpoint of two corners.
!> Node order = Gmsh order (ROADMAP D35):
!>   - corners 1-4 at the origin and the unit vectors; edge nodes 5-10 on the edges 1-2, 2-3,
!>     3-1, 4-1, 4-3, 4-2 (Gmsh manual "Node ordering");
!>   - tiny.msh element 1 is the reference element itself (corners at the origin and the unit
!>     vectors, edge nodes at the midpoints): after reading, its node coordinates equal
!>     tet10_node_coo_ref exactly;
!>   - Neper meshes sx2.msh and poly_part1.msh (unit cube; nodes in the file's order):
!>     every element edge node is at the midpoint of its two corners (Neper's elements are
!>     straight-edged; 1e-11 mm, as the coordinates are written with 12 decimals), det(J) > 0
!>     at every quadrature point, the volume integrated with the rule equals the corner
!>     tetrahedron's volume (relative 1e-8: the 5e-13 mm coordinate rounding gives relative
!>     differences up to about 1e-10 on the smallest elements), and the volumes add up to
!>     1 mm^3 (relative 1e-12).
!> Face element (6-node triangle) and face rule (7 points; step 3.3):
!>   - reference: the rule's points and weights equal the fixture's exactly, the shape
!>     functions and derivatives at 16 points to 1e-15 (until step 9.13 v2-fixed's qp2d, wt2d,
!>     sf2d06 and sf2dg06, also equal);
!>   - the rule equals its closed form (Felippa Ch. 24, Fig. 24.3) evaluated here with sqrt(15)
!>     (points 2e-16, weights relative 4e-16: 1-2 ulp), weights sum to 1/2, exactness for
!>     x^p y^q with p + q <= 5 (relative 1e-15), not 6;
!>   - Kronecker property, partition of unity, quadratic completeness (6 monomials), finite
!>     differences at 200 random points; node order: corners (0,0), (1,0), (0,1), edge nodes on
!>     1-2, 2-3, 3-1 (Gmsh manual, Triangle6);
!>   - the tet10 shape functions on the face xi3 = 0 equal the face shape functions at the face
!>     nodes 1, 2, 3, 5, 6, 7 and vanish at the other nodes (100 random points);
!>   - Neper fasets of sx2.msh and poly_part1.msh (each line reversed by the reader): every
!>     face's edge nodes 4, 5, 6 at the midpoints of 1-2, 2-3, 3-1 (1e-11 mm), the surface
!>     normal dx/dxi1 x dx/dxi2 at every face quadrature point along the normal of the cube face
!>     named by the faset (relative tangential part < 1e-10), and each faset's area integrated
!>     with the face rule equal to 1 mm^2 (relative 1e-12), and pointing out of the body
!>     (x0: -x, ..., z1: +z), the orientation the reader assumes (ROADMAP D37). (sx2.msh
!>     replaced the single-cell sx.msh, whose x0, y0, z0 faces Neper orients into the body:
!>     ROADMAP Q12.)
!>
!> Locals:
!>   element_dir, io_dir, mesh_dir  character(512)  command-line directories
!>   corner         integer(ik)(4)    corner nodes of the element order
!>   edge_corners   integer(ik)(2, 10)  the two corners of each edge node (0 for corners)
program test_tet10

  use kinds_mod, only: rk, ik
  use tet10_mod, only: tet10_num_nodes, tet10_node_coo_ref, tet10_num_qpts, tet10_qpt_coo, &
    & tet10_qpt_weight, tet10_shape, tet10_shape_deriv, tet10_face_num_nodes, &
    & tet10_face_node_coo_ref, tet10_face_num_qpts, tet10_face_qpt_coo, tet10_face_qpt_weight, &
    & tet10_face_shape, tet10_face_shape_deriv
  use mesh_file_mod, only: mesh_file_t, mesh_file_load
  use tensor_test_utils_mod, only: seed_random, determinant3
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary

  implicit none

  character(len=512) :: element_dir, io_dir, mesh_dir
  integer(ik), parameter :: corner(4) = [1_ik, 2_ik, 3_ik, 4_ik]
  integer(ik) :: edge_corners(2, tet10_num_nodes)

  call get_command_argument(1, element_dir)
  call get_command_argument(2, io_dir)
  call get_command_argument(3, mesh_dir)
  call seed_random(20261009)

  call check_reference(trim(element_dir) // '/tet10_reference.txt')
  call check_quadrature()
  call check_node_layout()
  call check_shape()
  call check_gmsh_order()
  call check_tiny(trim(io_dir) // '/tiny.msh')
  call check_neper_mesh(trim(mesh_dir) // '/sx2.msh', 'sx2')
  call check_neper_mesh(trim(mesh_dir) // '/poly_part1.msh', 'poly')
  call check_face_quadrature()
  call check_face_shape()
  call check_neper_fasets(trim(mesh_dir) // '/sx2.msh', 'sx2', .true.)
  call check_neper_fasets(trim(mesh_dir) // '/poly_part1.msh', 'poly', .true.)

  call test_summary()

contains

  !> Compare the quadrature rule and the shape functions with the reference fixture.
  !>
  !> Arguments:
  !>   path  (in)  character  fixture file
  !> Locals:
  !>   unit, ios    integer         file unit, I/O status
  !>   word         character(16)   section name
  !>   num          integer         section size
  !>   i_qpt, i_pt  integer         loop indices
  !>   qpt_line     real(rk)(4)     xi and weight of a quadrature point
  !>   pt_line      real(rk)(43)    xi, shape, d/dxi1, d/dxi2, d/dxi3 at a point
  !>   dshape       real(rk)(3,10)  FEPX derivatives
  !>   tri_qpt_line real(rk)(3)     xi and weight of a face quadrature point
  !>   tri_line     real(rk)(20)    xi, face shape, d/dxi1, d/dxi2 at a point
  !>   face_dshape  real(rk)(2,6)   FEPX face derivatives
  subroutine check_reference(path)

    character(len=*), intent(in) :: path

    integer :: unit, ios, num, i_qpt, i_pt
    character(len=16) :: word
    real(rk) :: qpt_line(4), pt_line(43), dshape(3, 10)
    real(rk) :: tri_qpt_line(3), tri_line(20), face_dshape(2, 6)

    open (newunit=unit, file=path, status='old', action='read', iostat=ios)
    call assert_equal(ios, 0, 'reference fixture open')
    if (ios /= 0) return

    ! Quadrature points and weights: exactly the reference's.
    read (unit, *) word, num
    call assert_equal(trim(word), 'qpts', 'reference qpts header')
    call assert_equal(num, int(tet10_num_qpts), 'reference num_qpts')
    do i_qpt = 1, num
      read (unit, *) qpt_line
      call assert_close(tet10_qpt_coo(:, i_qpt), qpt_line(1:3), 'reference qpt coo')
      call assert_close(tet10_qpt_weight(i_qpt), qpt_line(4), 'reference qpt weight')
    end do

    ! Shape functions and derivatives at the fixture points.
    read (unit, *) word, num
    call assert_equal(trim(word), 'points', 'reference points header')
    call assert_equal(num, 30, 'reference num points')
    do i_pt = 1, num
      read (unit, *) pt_line
      dshape = tet10_shape_deriv(pt_line(1:3))
      call assert_close(tet10_shape(pt_line(1:3)), pt_line(4:13), 'reference shape', &
        & abs_tol=1e-15_rk)
      call assert_close(dshape(1, :), pt_line(14:23), 'reference dshape/dxi1', abs_tol=1e-15_rk)
      call assert_close(dshape(2, :), pt_line(24:33), 'reference dshape/dxi2', abs_tol=1e-15_rk)
      call assert_close(dshape(3, :), pt_line(34:43), 'reference dshape/dxi3', abs_tol=1e-15_rk)
    end do

    ! Face rule: exactly the reference's.
    read (unit, *) word, num
    call assert_equal(trim(word), 'tri_qpts', 'reference tri_qpts header')
    call assert_equal(num, int(tet10_face_num_qpts), 'reference face num_qpts')
    do i_qpt = 1, num
      read (unit, *) tri_qpt_line
      call assert_close(tet10_face_qpt_coo(:, i_qpt), tri_qpt_line(1:2), 'reference face qpt coo')
      call assert_close(tet10_face_qpt_weight(i_qpt), tri_qpt_line(3), 'reference face qpt weight')
    end do

    ! Face shape functions and derivatives.
    read (unit, *) word, num
    call assert_equal(trim(word), 'tri_points', 'reference tri_points header')
    call assert_equal(num, 16, 'reference num tri points')
    do i_pt = 1, num
      read (unit, *) tri_line
      face_dshape = tet10_face_shape_deriv(tri_line(1:2))
      call assert_close(tet10_face_shape(tri_line(1:2)), tri_line(3:8), &
        & 'reference face shape', abs_tol=1e-15_rk)
      call assert_close(face_dshape(1, :), tri_line(9:14), 'reference face dshape/dxi1', &
        & abs_tol=1e-15_rk)
      call assert_close(face_dshape(2, :), tri_line(15:20), 'reference face dshape/dxi2', &
        & abs_tol=1e-15_rk)
    end do
    close (unit)

  end subroutine check_reference

  !> Weights, point locations and polynomial exactness of the 15-point rule.
  !>
  !> Locals:
  !>   p, q, r       integer   monomial exponents
  !>   i_qpt         integer   quadrature point
  !>   approx, exact real(rk)  rule value and exact integral of a monomial
  !>   err6          real(rk)  largest relative error over the degree-6 monomials
  !>   xi            real(rk)(3)  point coordinates
  subroutine check_quadrature()

    integer :: p, q, r, i_qpt
    real(rk) :: approx, exact, err6, xi(3)

    ! Weights and points.
    call assert_true(all(tet10_qpt_weight > 0.0_rk), 'qpt weights positive')
    call assert_close(sum(tet10_qpt_weight), 1.0_rk/6.0_rk, 'qpt weights sum 1/6', rel_tol=1e-15_rk)
    do i_qpt = 1, int(tet10_num_qpts)
      xi = tet10_qpt_coo(:, i_qpt)
      call assert_true(all(xi > -1e-15_rk) .and. sum(xi) < 1.0_rk + 1e-15_rk, 'qpt inside')
    end do

    ! Exactness up to degree 5; error at degree 6.
    err6 = 0.0_rk
    do p = 0, 6
      do q = 0, 6 - p
        do r = 0, 6 - p - q
          approx = 0.0_rk
          do i_qpt = 1, int(tet10_num_qpts)
            xi = tet10_qpt_coo(:, i_qpt)
            approx = approx + tet10_qpt_weight(i_qpt)*xi(1)**p*xi(2)**q*xi(3)**r
          end do
          exact = gamma(real(p + 1, rk))*gamma(real(q + 1, rk))*gamma(real(r + 1, rk)) &
            & /gamma(real(p + q + r + 4, rk))
          if (p + q + r <= 5) then
            call assert_close(approx, exact, 'qpt exact degree <= 5', rel_tol=1e-14_rk)
          else
            err6 = max(err6, abs(approx - exact)/exact)
          end if
        end do
      end do
    end do
    call assert_true(err6 > 1e-6_rk, 'qpt not exact at degree 6')

  end subroutine check_quadrature

  !> Corners and edge nodes of the element order; fills edge_corners.
  !>
  !> Locals:
  !>   a, i, j   integer(ik)  node indices
  !>   found     integer      number of corner pairs whose midpoint is node a
  subroutine check_node_layout()

    integer(ik) :: a, i, j
    integer :: found

    ! Corners: the nodes with integer reference coordinates.
    do a = 1_ik, tet10_num_nodes
      call assert_true(all(abs(tet10_node_coo_ref(:, a) - anint(tet10_node_coo_ref(:, a))) &
        & < 1e-15_rk) .eqv. any(corner == a), 'node layout corner')
    end do

    ! Every other node is the midpoint of exactly one pair of corners.
    edge_corners = 0_ik
    do a = 1_ik, tet10_num_nodes
      if (any(corner == a)) cycle
      found = 0
      do i = 1_ik, 3_ik
        do j = i + 1_ik, 4_ik
          if (all(abs(tet10_node_coo_ref(:, a) - 0.5_rk*(tet10_node_coo_ref(:, corner(i)) + &
            & tet10_node_coo_ref(:, corner(j)))) < 1e-15_rk)) then
            found = found + 1
            edge_corners(:, a) = [corner(i), corner(j)]
          end if
        end do
      end do
      call assert_equal(found, 1, 'node layout edge node')
    end do

  end subroutine check_node_layout

  !> Kronecker property, partition of unity, quadratic completeness, finite differences.
  !>
  !> Locals:
  !>   a, b, k, i_pt  integer       loop indices
  !>   xi, xi_h       real(rk)(3)   point, shifted point
  !>   shape          real(rk)(10)  shape functions
  !>   dshape, dfd    real(rk)(3,10)  derivatives, finite-difference derivatives
  !>   delta          real(rk)(10)  Kronecker column
  !>   h              real(rk)      finite-difference step
  !>   q_nodes        real(rk)(10)  monomial values at the nodes
  !>   q_val          real(rk)      monomial value at xi
  !>   q_grad         real(rk)(3)   monomial gradient at xi
  subroutine check_shape()

    integer :: a, b, k, i_pt
    real(rk) :: xi(3), xi_h(3), shape(10), dshape(3, 10), dfd(3, 10), delta(10), h
    real(rk) :: q_nodes(10), q_val, q_grad(3)

    ! N_a(node b) = delta_ab.
    do b = 1, int(tet10_num_nodes)
      delta = 0.0_rk
      delta(b) = 1.0_rk
      call assert_close(tet10_shape(tet10_node_coo_ref(:, b)), delta, 'shape kronecker', &
        & abs_tol=1e-15_rk)
    end do

    ! Random points of the element.
    h = 1e-5_rk
    do i_pt = 1, 200
      xi = random_point()
      shape = tet10_shape(xi)
      dshape = tet10_shape_deriv(xi)

      ! Partition of unity.
      call assert_close(sum(shape), 1.0_rk, 'shape partition of unity', abs_tol=1e-14_rk)
      call assert_close(sum(dshape, dim=2), [0.0_rk, 0.0_rk, 0.0_rk], 'dshape sum zero', &
        & abs_tol=1e-13_rk)

      ! Quadratic completeness, monomial by monomial.
      do k = 1, 10
        do a = 1, int(tet10_num_nodes)
          call monomial(k, tet10_node_coo_ref(:, a), q_nodes(a), q_grad)
        end do
        call monomial(k, xi, q_val, q_grad)
        call assert_close(dot_product(shape, q_nodes), q_val, 'shape complete', &
          & abs_tol=1e-13_rk)
        call assert_close(matmul(dshape, q_nodes), q_grad, 'dshape complete', abs_tol=1e-13_rk)
      end do

      ! Central finite differences.
      do k = 1, 3
        xi_h = xi
        xi_h(k) = xi(k) + h
        dfd(k, :) = tet10_shape(xi_h)
        xi_h(k) = xi(k) - h
        dfd(k, :) = (dfd(k, :) - tet10_shape(xi_h))/(2.0_rk*h)
      end do
      call assert_close(reshape(dshape, [30]), reshape(dfd, [30]), 'dshape finite difference', &
        & abs_tol=1e-9_rk)
    end do

  end subroutine check_shape

  !> The element node order against Gmsh's documented order (ROADMAP D35).
  !>
  !> Locals:
  !>   gmsh_edge  integer(ik)(2, 6)  corners of the edges of edge nodes 5-10 (Gmsh manual)
  !>   i_edge     integer(ik)        edge index
  !>   unit_vec   real(rk)(3,4)      reference coordinates of corners 1-4
  subroutine check_gmsh_order()

    integer(ik), parameter :: gmsh_edge(2, 6) = reshape([1_ik, 2_ik, 2_ik, 3_ik, 3_ik, 1_ik, &
      & 4_ik, 1_ik, 4_ik, 3_ik, 4_ik, 2_ik], [2, 6])
    real(rk), parameter :: unit_vec(3, 4) = reshape([0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, &
      & 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk], [3, 4])
    integer(ik) :: i_edge

    ! Corners 1-4 at the origin and the unit vectors (Gmsh's reference tetrahedron).
    call assert_close(reshape(tet10_node_coo_ref(:, 1:4), [12]), reshape(unit_vec, [12]), &
      & 'gmsh order corners')

    ! Edge node 4 + k on Gmsh edge k (1-2, 2-3, 3-1, 4-1, 4-3, 4-2).
    do i_edge = 1_ik, 6_ik
      call assert_true(all(edge_corners(:, 4_ik + i_edge) == gmsh_edge(:, i_edge)) .or. &
        & all(edge_corners(:, 4_ik + i_edge) == gmsh_edge(2:1:-1, i_edge)), &
        & 'gmsh order edge node')
    end do

  end subroutine check_gmsh_order

  !> tiny.msh element 1 is the reference element: its coordinates after reading.
  !>
  !> Arguments:
  !>   path  (in)  character  mesh file
  !> Locals:
  !>   mesh  type(mesh_file_t)  parsed mesh
  subroutine check_tiny(path)

    character(len=*), intent(in) :: path

    type(mesh_file_t) :: mesh

    call mesh_file_load(mesh, path)
    call assert_close(reshape(mesh%coo(:, mesh%elt_nodes(:, 1)), [30]), &
      & reshape(tet10_node_coo_ref, [30]), 'tiny element 1 is the reference element')

  end subroutine check_tiny

  !> Geometry of every element of a Neper mesh of the unit cube (straight-edged elements).
  !>
  !> Arguments:
  !>   path   (in)  character  mesh file
  !>   label  (in)  character  check label prefix
  !> Locals:
  !>   mesh          type(mesh_file_t)  parsed mesh
  !>   i_elt, i_qpt  integer(ik)        loop indices
  !>   a             integer(ik)        element node
  !>   coo           real(rk)(3,10)     element node coordinates [mm]
  !>   jac           real(rk)(3,3)      jac(i, j) = d x_i / d xi_j [mm]
  !>   det_jac       real(rk)           det(jac) [mm^3]
  !>   vol, vol_lin  real(rk)           element volume with the rule, of the corner tetrahedron
  !>   vol_total     real(rk)           mesh volume [mm^3]
  !>   min_det       real(rk)           smallest det(jac) over the mesh
  !>   max_mid_err   real(rk)           largest edge-node distance from its edge midpoint [mm]
  !>   max_vol_err   real(rk)           largest relative difference vol / vol_lin - 1
  subroutine check_neper_mesh(path, label)

    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: label

    type(mesh_file_t) :: mesh
    integer(ik) :: i_elt, i_qpt, a
    real(rk) :: coo(3, 10), jac(3, 3), det_jac, vol, vol_lin, vol_total
    real(rk) :: min_det, max_mid_err, max_vol_err

    call mesh_file_load(mesh, path)
    call assert_equal(mesh%num_elt_nodes, tet10_num_nodes, label // ' tet10 mesh')

    vol_total = 0.0_rk
    min_det = huge(1.0_rk)
    max_mid_err = 0.0_rk
    max_vol_err = 0.0_rk
    do i_elt = 1_ik, mesh%num_elts
      coo = mesh%coo(:, mesh%elt_nodes(:, i_elt))

      ! Edge nodes at the midpoints of their corners.
      do a = 1_ik, tet10_num_nodes
        if (any(corner == a)) cycle
        max_mid_err = max(max_mid_err, maxval(abs(coo(:, a) - 0.5_rk*(coo(:, &
          & edge_corners(1, a)) + coo(:, edge_corners(2, a))))))
      end do

      ! Jacobian determinant at the quadrature points and integrated volume.
      vol = 0.0_rk
      do i_qpt = 1_ik, tet10_num_qpts
        jac = matmul(coo, transpose(tet10_shape_deriv(tet10_qpt_coo(:, i_qpt))))
        det_jac = determinant3(jac)
        min_det = min(min_det, det_jac)
        vol = vol + tet10_qpt_weight(i_qpt)*det_jac
      end do

      ! Volume of the corner tetrahedron.
      vol_lin = determinant3(reshape([coo(:, 2) - coo(:, 1), coo(:, 3) - coo(:, 1), &
        & coo(:, 4) - coo(:, 1)], [3, 3]))/6.0_rk
      max_vol_err = max(max_vol_err, abs(vol/vol_lin - 1.0_rk))
      vol_total = vol_total + vol
    end do

    call assert_close(max_mid_err, 0.0_rk, label // ' edge nodes at midpoints', abs_tol=1e-11_rk)
    call assert_true(min_det > 0.0_rk, label // ' det(J) > 0 at every quadrature point')
    call assert_close(max_vol_err, 0.0_rk, label // ' element volume = corner tetrahedron', &
      & abs_tol=1e-8_rk)
    call assert_close(vol_total, 1.0_rk, label // ' mesh volume 1', rel_tol=1e-12_rk)

  end subroutine check_neper_mesh

  !> A random point of the reference element (uniform in the unit cube, folded into the
  !> tetrahedron by sorting: xi = (u1, u2 - u1, u3 - u2) with u sorted).
  !>
  !> Locals:
  !>   u  real(rk)(3)  sorted uniform numbers
  !>   t  real(rk)     swap
  function random_point() result(xi)

    real(rk) :: xi(3)

    real(rk) :: u(3), t

    call random_number(u)
    if (u(1) > u(2)) then
      t = u(1); u(1) = u(2); u(2) = t
    end if
    if (u(2) > u(3)) then
      t = u(2); u(2) = u(3); u(3) = t
    end if
    if (u(1) > u(2)) then
      t = u(1); u(1) = u(2); u(2) = t
    end if
    xi = [u(1), u(2) - u(1), u(3) - u(2)]

  end function random_point

  !> Value and gradient of the k-th monomial of degree <= 2:
  !> 1, x, y, z, x^2, y^2, z^2, xy, yz, xz (x, y, z = xi1, xi2, xi3).
  !>
  !> Arguments:
  !>   k     (in)   integer      monomial index 1..10
  !>   xi    (in)   real(rk)(3)  point
  !>   val   (out)  real(rk)     value
  !>   grad  (out)  real(rk)(3)  gradient
  subroutine monomial(k, xi, val, grad)

    integer, intent(in) :: k
    real(rk), intent(in) :: xi(3)
    real(rk), intent(out) :: val, grad(3)

    grad = 0.0_rk
    select case (k)
    case (1)
      val = 1.0_rk
    case (2:4)
      val = xi(k - 1)
      grad(k - 1) = 1.0_rk
    case (5:7)
      val = xi(k - 4)**2
      grad(k - 4) = 2.0_rk*xi(k - 4)
    case (8)
      val = xi(1)*xi(2)
      grad = [xi(2), xi(1), 0.0_rk]
    case (9)
      val = xi(2)*xi(3)
      grad = [0.0_rk, xi(3), xi(2)]
    case default
      val = xi(1)*xi(3)
      grad = [xi(3), 0.0_rk, xi(1)]
    end select

  end subroutine monomial

  !> The face rule: closed form, weights, exactness.
  !>
  !> Locals:
  !>   r15, g1, g2    real(rk)       sqrt(15), (6 - r15)/21, (6 + r15)/21
  !>   coo_ref        real(rk)(2,7)  closed-form points (v2's order)
  !>   w_ref          real(rk)(7)    closed-form weights
  !>   p, q, n        integer        exponents, point
  !>   approx, exact  real(rk)       rule value and exact integral
  !>   err6           real(rk)       largest relative error at degree 6
  subroutine check_face_quadrature()

    real(rk) :: r15, g1, g2, coo_ref(2, 7), w_ref(7), approx, exact, err6
    integer :: p, q, n

    ! Closed form (Felippa Ch. 24, Fig. 24.3; weights halved for the area 1/2).
    r15 = sqrt(15.0_rk)
    g1 = (6.0_rk - r15)/21.0_rk
    g2 = (6.0_rk + r15)/21.0_rk
    coo_ref = reshape([1.0_rk/3.0_rk, 1.0_rk/3.0_rk, 1.0_rk - 2.0_rk*g2, g2, g2, &
      & 1.0_rk - 2.0_rk*g2, g2, g2, 1.0_rk - 2.0_rk*g1, g1, g1, 1.0_rk - 2.0_rk*g1, g1, g1], [2, 7])
    w_ref = [9.0_rk/80.0_rk, (155.0_rk + r15)/2400.0_rk, (155.0_rk + r15)/2400.0_rk, &
      & (155.0_rk + r15)/2400.0_rk, (155.0_rk - r15)/2400.0_rk, (155.0_rk - r15)/2400.0_rk, &
      & (155.0_rk - r15)/2400.0_rk]
    call assert_equal(int(tet10_face_num_qpts), 7, 'face num_qpts')
    call assert_close(reshape(tet10_face_qpt_coo, [14]), reshape(coo_ref, [14]), &
      & 'face qpt closed form', abs_tol=2e-16_rk)
    ! The weights evaluated here differ from the correctly rounded literals by up to 1 ulp.
    call assert_close(tet10_face_qpt_weight, w_ref, 'face weight closed form', rel_tol=4e-16_rk)

    ! Weights.
    call assert_true(all(tet10_face_qpt_weight > 0.0_rk), 'face weights positive')
    call assert_close(sum(tet10_face_qpt_weight), 0.5_rk, 'face weights sum 1/2', &
      & rel_tol=1e-15_rk)

    ! Exactness up to degree 5, not at degree 6.
    err6 = 0.0_rk
    do p = 0, 6
      do q = 0, 6 - p
        approx = 0.0_rk
        do n = 1, int(tet10_face_num_qpts)
          approx = approx + tet10_face_qpt_weight(n)*tet10_face_qpt_coo(1, n)**p* &
            & tet10_face_qpt_coo(2, n)**q
        end do
        exact = gamma(real(p + 1, rk))*gamma(real(q + 1, rk))/gamma(real(p + q + 3, rk))
        if (p + q <= 5) then
          call assert_close(approx, exact, 'face qpt exact degree <= 5', rel_tol=1e-15_rk)
        else
          err6 = max(err6, abs(approx - exact)/exact)
        end if
      end do
    end do
    call assert_true(err6 > 1e-6_rk, 'face qpt not exact at degree 6')

  end subroutine check_face_quadrature

  !> Face shape functions: node order, Kronecker property, partition of unity, completeness,
  !> finite differences, and the restriction of the tet10 functions to the face xi3 = 0.
  !>
  !> Locals:
  !>   face_nodes     integer(ik)(6)  tet10 nodes on the face xi3 = 0, in face order
  !>   coo_ref        real(rk)(2,6)   expected face node coordinates (Gmsh Triangle6)
  !>   a, b, k, i_pt  integer         loop indices
  !>   xi, xi_h       real(rk)(2)     point, shifted point
  !>   shape, delta   real(rk)(6)     face shape functions, Kronecker column
  !>   dshape, dfd    real(rk)(2,6)   derivatives, finite-difference derivatives
  !>   q_nodes        real(rk)(6)     monomial values at the face nodes
  !>   q_val, q_grad  real(rk), real(rk)(2)  monomial value and gradient
  !>   q_grad3        real(rk)(3)     monomial gradient in 3D (monomial())
  !>   u              real(rk)(3)     uniform random numbers
  !>   vol_shape      real(rk)(10)    tet10 shape functions on the face
  !>   on_face        logical(10)     whether each tet10 node is on the face
  !>   h              real(rk)        finite-difference step
  subroutine check_face_shape()

    integer(ik), parameter :: face_nodes(6) = [1_ik, 2_ik, 3_ik, 5_ik, 6_ik, 7_ik]
    real(rk), parameter :: coo_ref(2, 6) = reshape([0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, &
      & 1.0_rk, 0.5_rk, 0.0_rk, 0.5_rk, 0.5_rk, 0.0_rk, 0.5_rk], [2, 6])
    integer :: a, b, k, i_pt
    real(rk) :: xi(2), xi_h(2), shape(6), delta(6), dshape(2, 6), dfd(2, 6), q_nodes(6)
    real(rk) :: q_val, q_grad(2), q_grad3(3), vol_shape(10), h, u(3)
    logical :: on_face(10)

    ! Node order (Gmsh Triangle6).
    call assert_equal(int(tet10_face_num_nodes), 6, 'face num_nodes')
    call assert_close(reshape(tet10_face_node_coo_ref, [12]), reshape(coo_ref, [12]), &
      & 'face node coordinates')

    ! N_a(node b) = delta_ab.
    do b = 1, 6
      delta = 0.0_rk
      delta(b) = 1.0_rk
      call assert_close(tet10_face_shape(tet10_face_node_coo_ref(:, b)), delta, &
        & 'face shape kronecker', abs_tol=1e-15_rk)
    end do

    on_face = .false.
    on_face(face_nodes) = .true.
    h = 1e-5_rk
    do i_pt = 1, 200
      ! A random point of the reference triangle (a random tet point with xi3 dropped would
      ! not be uniform; fold the unit square instead).
      call random_number(u)
      xi = u(1:2)
      if (sum(xi) > 1.0_rk) xi = 1.0_rk - xi
      shape = tet10_face_shape(xi)
      dshape = tet10_face_shape_deriv(xi)

      ! Partition of unity.
      call assert_close(sum(shape), 1.0_rk, 'face partition of unity', abs_tol=1e-14_rk)
      call assert_close(sum(dshape, dim=2), [0.0_rk, 0.0_rk], 'face dshape sum zero', &
        & abs_tol=1e-13_rk)

      ! Quadratic completeness (monomials 1, x, y, x^2, y^2, xy: those of monomial() without
      ! xi3, i.e. numbers 1, 2, 3, 5, 6, 8).
      do k = 1, 10
        if (any(k == [4, 7, 9, 10])) cycle
        do a = 1, 6
          call monomial(k, [tet10_face_node_coo_ref(:, a), 0.0_rk], q_nodes(a), q_grad3)
        end do
        call monomial(k, [xi, 0.0_rk], q_val, q_grad3)
        q_grad = q_grad3(1:2)
        call assert_close(dot_product(shape, q_nodes), q_val, 'face shape complete', &
          & abs_tol=1e-13_rk)
        call assert_close(matmul(dshape, q_nodes), q_grad, 'face dshape complete', &
          & abs_tol=1e-13_rk)
      end do

      ! Central finite differences.
      do k = 1, 2
        xi_h = xi
        xi_h(k) = xi(k) + h
        dfd(k, :) = tet10_face_shape(xi_h)
        xi_h(k) = xi(k) - h
        dfd(k, :) = (dfd(k, :) - tet10_face_shape(xi_h))/(2.0_rk*h)
      end do
      call assert_close(reshape(dshape, [12]), reshape(dfd, [12]), &
        & 'face dshape finite difference', abs_tol=1e-9_rk)

      ! Restriction of the tet10 functions to the face xi3 = 0.
      if (i_pt <= 100) then
        vol_shape = tet10_shape([xi, 0.0_rk])
        call assert_close(vol_shape(face_nodes), shape, 'face = tet10 on xi3 = 0', &
          & abs_tol=1e-15_rk)
        call assert_true(all(abs(pack(vol_shape, .not. on_face)) < 1e-15_rk), &
          & 'tet10 off-face nodes vanish on xi3 = 0')
      end if
    end do

  end subroutine check_face_shape

  !> Neper fasets: Gmsh order after the reader's reversal, normal direction, faset areas, and
  !> optionally outward orientation.
  !>
  !> Arguments:
  !>   path           (in)  character  mesh file
  !>   label          (in)  character  check label prefix
  !>   check_outward  (in)  logical    also check that the normals point out of the body
  !> Locals:
  !>   mesh             type(mesh_file_t)  parsed mesh
  !>   i_faset, i_face  integer            loop indices
  !>   i_qpt            integer(ik)        face quadrature point
  !>   coo              real(rk)(3,6)      face node coordinates [mm]
  !>   tang             real(rk)(3,2)      dx/dxi1, dx/dxi2 [mm]
  !>   normal           real(rk)(3)        dx/dxi1 x dx/dxi2 [mm^2]
  !>   outward          real(rk)(3)        outward normal of the cube face of the faset
  !>   area             real(rk)           faset area [mm^2]
  !>   max_mid_err      real(rk)           largest edge-node distance from its midpoint [mm]
  !>   max_tangential   real(rk)           largest |normal - (normal . outward) outward|/|normal|
  !>   min_outward      real(rk)           smallest (normal . outward)/|normal|
  subroutine check_neper_fasets(path, label, check_outward)

    character(len=*), intent(in) :: path
    character(len=*), intent(in) :: label
    logical, intent(in) :: check_outward

    type(mesh_file_t) :: mesh
    integer :: i_faset, i_face
    integer(ik) :: i_qpt
    real(rk) :: coo(3, 6), tang(3, 2), normal(3), outward(3), area
    real(rk) :: max_mid_err, max_tangential, min_outward

    call mesh_file_load(mesh, path)
    call assert_equal(size(mesh%fasets), 6, label // ' 6 fasets')
    max_mid_err = 0.0_rk
    max_tangential = 0.0_rk
    min_outward = huge(1.0_rk)
    do i_faset = 1, size(mesh%fasets)
      ! Outward normal of the cube face named by the label (x0, x1, y0, y1, z0, z1).
      outward = 0.0_rk
      outward(index('xyz', mesh%fasets(i_faset)%label(1:1))) = &
        & merge(-1.0_rk, 1.0_rk, mesh%fasets(i_faset)%label(2:2) == '0')
      call assert_equal(int(mesh%fasets(i_faset)%num_face_nodes), 6, label // ' tri6 fasets')
      area = 0.0_rk
      do i_face = 1, size(mesh%fasets(i_faset)%elts)
        coo = mesh%coo(:, mesh%fasets(i_faset)%nodes(:, i_face))
        ! Edge nodes 4, 5, 6 at the midpoints of 1-2, 2-3, 3-1.
        max_mid_err = max(max_mid_err, maxval(abs(coo(:, 4) - 0.5_rk*(coo(:, 1) + coo(:, 2)))), &
          & maxval(abs(coo(:, 5) - 0.5_rk*(coo(:, 2) + coo(:, 3)))), &
          & maxval(abs(coo(:, 6) - 0.5_rk*(coo(:, 3) + coo(:, 1)))))
        ! Normal at each face quadrature point and area.
        do i_qpt = 1_ik, tet10_face_num_qpts
          tang = matmul(coo, transpose(tet10_face_shape_deriv(tet10_face_qpt_coo(:, i_qpt))))
          normal = cross(tang(:, 1), tang(:, 2))
          min_outward = min(min_outward, dot_product(normal, outward)/norm2(normal))
          max_tangential = max(max_tangential, norm2(normal - dot_product(normal, outward)* &
            & outward)/norm2(normal))
          ! (normal . outward is +-|normal| on every face; its sign is checked below.)
          area = area + tet10_face_qpt_weight(i_qpt)*norm2(normal)
        end do
      end do
      call assert_close(area, 1.0_rk, label // ' faset area ' // mesh%fasets(i_faset)%label, &
        & rel_tol=1e-12_rk)
    end do
    call assert_close(max_mid_err, 0.0_rk, label // ' faset edge nodes at midpoints', &
      & abs_tol=1e-11_rk)
    if (check_outward) call assert_true(min_outward > 0.0_rk, label // &
      & ' faset normals point outward')
    call assert_close(max_tangential, 0.0_rk, label // ' faset normals along the face normal', &
      & abs_tol=1e-10_rk)

  end subroutine check_neper_fasets

  !> Cross product of two 3-vectors.
  !>
  !> Arguments:
  !>   u, v  (in)  real(rk)(3)  vectors
  pure function cross(u, v) result(w)

    real(rk), intent(in) :: u(3), v(3)
    real(rk) :: w(3)

    w = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1)]

  end function cross

end program test_tet10
