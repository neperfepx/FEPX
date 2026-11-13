! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of hex20_mod (20-node serendipity hexahedron and its 27-point quadrature rule).
!>
!> Source : FEPX (ROADMAP step 3.2); checks of dev/theory/addendum_hex20.md §1-§3.
!>
!> Usage: test_hex20 <element fixture dir> (tests/element/fixtures).
!>
!> Quadrature (27 points):
!>   - the one-dimensional points are the zeros of P3(x) = (5x^3 - 3x)/2 and the weights
!>     2 / ((1 - x^2) P3'(x)^2) (Abramowitz & Stegun Eq. 25.4.29; 1e-14, as the formula amplifies
!>     the rounding of x); the 27 weights equal the products of 5/9 and 8/9 (125/729, ...;
!>     4e-16); the points are the tensor product in the order xi1 fastest, then xi2, then xi3;
!>   - weights positive, sum 8 (the reference volume); points inside the element;
!>   - exactness: every xi1^p xi2^q xi3^r with p, q, r <= 5 is integrated exactly (1e-14);
!>     not for exponent 6.
!> Shape functions:
!>   - N_a(node b) = delta_ab at the reference node coordinates;
!>   - at 200 random points of [-1, 1]^3: partition of unity (sum N = 1, sum dN = 0);
!>     serendipity completeness: sum_a N_a q(xi_a) = q(xi) and sum_a dN_a q(xi_a) = grad q for
!>     the 20 monomials of superlinear degree <= 2 (Arnold & Awanou 2011; 1e-13); not for
!>     xi1^2 xi2^2 (outside the serendipity space); derivatives against central finite
!>     differences (1e-8).
!> Node order = Gmsh order (ROADMAP D35, addendum §1):
!>   - corners 1-8 at the sign pattern of the addendum's table ((-1,-1,-1), (1,-1,-1), ...);
!>     (x2 - x1) x (x4 - x1) . (x5 - x1) > 0;
!>   - every other node is the midpoint of exactly one edge (two corners differing in one
!>     coordinate), edge node 8 + k on the k-th edge of Gmsh's list 1-2, 1-4, 1-5, 2-3, 2-6,
!>     3-4, 3-7, 4-8, 5-6, 5-8, 6-7, 7-8 (Gmsh manual §10.2.1).
!>   - Neper mesh hex_cube.msh (unit cube, 125 elements of Gmsh type 17, read by
!>     mesh_file_load): every edge node at the midpoint of its corners (1e-11 mm), det(J) > 0
!>     at every quadrature point, element volume = (x2 - x1) x (x4 - x1) . (x5 - x1) (the
!>     elements are rectangular boxes; relative 1e-10), total volume 1 mm^3 (relative 1e-12).
!> Face element (8-node quadrangle) and face rule (9 points; step 3.3):
!>   - the rule: 3 x 3 tensor product of the same one-dimensional points, xi1 fastest, weights
!>     exactly products of 5/9 and 8/9 (25/81, 40/81, 64/81; 4e-16), sum 4, exactness for
!>     x^p y^q with p, q <= 5 (1e-14), not for exponent 6;
!>   - node order (Gmsh Quadrangle8): corners (-1,-1), (1,-1), (1,1), (-1,1), edge nodes on
!>     1-2, 2-3, 3-4, 4-1;
!>   - Kronecker property, partition of unity, completeness on the 8 two-dimensional
!>     serendipity monomials (1, x, y, x^2, xy, y^2, x^2 y, x y^2) and not on x^2 y^2, finite
!>     differences, at 200 random points;
!>   - the hex20 shape functions on the face xi3 = -1 equal the face shape functions at the
!>     face nodes 1, 2, 3, 4, 9, 12, 14, 10 and vanish at the other nodes (100 random points);
!>   - reference: values and derivatives equal FEPX's own at the 8 nodes and 4 other points
!>     (fixture quad8_reference.txt, written by dump_quad8.f90 in step 9.13; 1e-15; until step
!>     9.13 v2-fixed's sf2d08 / sf2dg08 under xi = 2 s - 1 and v2's node order, also 1e-15).
!>
!> Locals:
!>   element_dir    character(512)      fixture directory (command-line argument)
!>   edge_corners   integer(ik)(2, 20)  the two corners of each edge node (0 for corners)
program test_hex20

  use kinds_mod, only: rk, ik
  use hex20_mod, only: hex20_num_nodes, hex20_node_coo_ref, hex20_num_qpts, hex20_qpt_coo, &
    & hex20_qpt_weight, hex20_shape, hex20_shape_deriv, hex20_face_num_nodes, &
    & hex20_face_node_coo_ref, hex20_face_num_qpts, hex20_face_qpt_coo, hex20_face_qpt_weight, &
    & hex20_face_shape, hex20_face_shape_deriv
  use mesh_file_mod, only: mesh_file_t, mesh_file_load, elt_type_hex20
  use tensor_test_utils_mod, only: seed_random, determinant3
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary

  implicit none

  character(len=512) :: element_dir
  integer(ik) :: edge_corners(2, hex20_num_nodes)

  call get_command_argument(1, element_dir)
  call seed_random(20261009)

  call check_quadrature()
  call check_node_layout()
  call check_gmsh_order()
  call check_shape()
  call check_neper_mesh(trim(element_dir) // '/hex_cube.msh')
  call check_face_quadrature()
  call check_face_shape()
  call check_face_reference(trim(element_dir) // '/quad8_reference.txt')

  call test_summary()

contains

  !> The 27-point rule: one-dimensional Gauss-Legendre points and weights, tensor order,
  !> weights, exactness.
  !>
  !> Locals:
  !>   weight_num     integer(3)   one-dimensional weights times 9 (5, 8, 5)
  !>   x, w           real(rk)(3)  one-dimensional points and A&S weights at those points
  !>   n, i, j, k     integer      point and tensor indices
  !>   p, q, r        integer      monomial exponents
  !>   approx, exact  real(rk)     rule value and exact integral
  !>   err6           real(rk)     largest error with an exponent 6
  subroutine check_quadrature()

    integer, parameter :: weight_num(3) = [5, 8, 5]
    real(rk) :: x(3), w(3), approx, exact, err6
    integer :: n, i, j, k, p, q, r

    call assert_equal(int(hex20_num_qpts), 27, 'hex num_qpts')

    ! One-dimensional rule from the first three points (xi1 varies fastest).
    x = hex20_qpt_coo(1, 1:3)
    do i = 1, 3
      ! Zeros of P3 and the Gauss weights of A&S 25.4.29.
      call assert_close(0.5_rk*(5.0_rk*x(i)**3 - 3.0_rk*x(i)), 0.0_rk, 'hex gauss point P3 = 0', &
        & abs_tol=1e-15_rk)
      w(i) = 2.0_rk/((1.0_rk - x(i)**2)*(0.5_rk*(15.0_rk*x(i)**2 - 3.0_rk))**2)
    end do
    call assert_true(x(1) < x(2) .and. x(2) < x(3), 'hex gauss points increasing')
    call assert_close(x, [-sqrt(0.6_rk), 0.0_rk, sqrt(0.6_rk)], 'hex gauss points sqrt(3/5)', &
      & abs_tol=2e-16_rk)

    ! Tensor product: point i + 3(j - 1) + 9(k - 1).
    do k = 1, 3
      do j = 1, 3
        do i = 1, 3
          n = i + 3*(j - 1) + 9*(k - 1)
          call assert_close(hex20_qpt_coo(:, n), [x(i), x(j), x(k)], 'hex qpt tensor order')
          ! The A&S weights evaluated at the rounded points (the formula amplifies the
          ! rounding of x: 1e-14) and the exact products of 5/9 and 8/9 (2 ulp).
          call assert_close(hex20_qpt_weight(n), w(i)*w(j)*w(k), 'hex qpt weight product', &
            & rel_tol=1e-14_rk)
          call assert_close(hex20_qpt_weight(n), real(weight_num(i)*weight_num(j)*weight_num(k), &
            & rk)/729.0_rk, 'hex qpt weight exact', rel_tol=4e-16_rk)
        end do
      end do
    end do

    ! Weights and points.
    call assert_true(all(hex20_qpt_weight > 0.0_rk), 'hex weights positive')
    call assert_close(sum(hex20_qpt_weight), 8.0_rk, 'hex weights sum 8', rel_tol=1e-15_rk)
    call assert_true(all(abs(hex20_qpt_coo) < 1.0_rk), 'hex qpts inside')

    ! Exactness for exponents <= 5; an exponent 6 is not integrated exactly.
    err6 = 0.0_rk
    do p = 0, 6
      do q = 0, 6
        do r = 0, 6
          approx = 0.0_rk
          do n = 1, int(hex20_num_qpts)
            approx = approx + hex20_qpt_weight(n)*hex20_qpt_coo(1, n)**p*hex20_qpt_coo(2, n)**q* &
              & hex20_qpt_coo(3, n)**r
          end do
          exact = integral_1d(p)*integral_1d(q)*integral_1d(r)
          if (max(p, q, r) <= 5) then
            call assert_close(approx, exact, 'hex qpt exact', abs_tol=1e-14_rk)
          else
            err6 = max(err6, abs(approx - exact))
          end if
        end do
      end do
    end do
    call assert_true(err6 > 1e-3_rk, 'hex qpt not exact for exponent 6')

  end subroutine check_quadrature

  !> Corners and edge nodes of the element order; fills edge_corners.
  !>
  !> Locals:
  !>   a, b, c  integer(ik)  nodes
  !>   found    integer      corner pairs whose midpoint is node a
  subroutine check_node_layout()

    integer(ik) :: a, b, c
    integer :: found

    ! Corners (all coordinates +-1) are nodes 1-8; the others have one coordinate 0.
    do a = 1_ik, hex20_num_nodes
      call assert_true(all(abs(abs(hex20_node_coo_ref(:, a)) - 1.0_rk) < 1e-15_rk) .eqv. &
        & (a <= 8_ik), 'hex layout corner')
    end do

    ! Each edge node is the midpoint of exactly one pair of corners, which share an edge.
    edge_corners = 0_ik
    do a = 9_ik, hex20_num_nodes
      found = 0
      do b = 1_ik, 7_ik
        do c = b + 1_ik, 8_ik
          if (all(abs(hex20_node_coo_ref(:, a) - 0.5_rk*(hex20_node_coo_ref(:, b) + &
            & hex20_node_coo_ref(:, c))) < 1e-15_rk)) then
            found = found + 1
            edge_corners(:, a) = [b, c]
            call assert_equal(count(abs(hex20_node_coo_ref(:, b) - hex20_node_coo_ref(:, c)) &
              & > 1.0_rk), 1, 'hex layout edge of the element')
          end if
        end do
      end do
      call assert_equal(found, 1, 'hex layout edge node')
    end do

  end subroutine check_node_layout

  !> The element node order against the Gmsh manual (addendum §1).
  !>
  !> Locals:
  !>   corner_coo  real(rk)(3, 8)     corner coordinates of the addendum's table
  !>   gmsh_edge   integer(ik)(2, 12) the corners of the edges of nodes 9-20
  !>   k           integer(ik)        edge index
  !>   x           real(rk)(3, 8)     corner coordinates
  subroutine check_gmsh_order()

    real(rk), parameter :: corner_coo(3, 8) = reshape([ &
      & -1.0_rk, -1.0_rk, -1.0_rk, 1.0_rk, -1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, &
      & -1.0_rk, 1.0_rk, -1.0_rk, -1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, &
      & 1.0_rk, 1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, 1.0_rk], [3, 8])
    integer(ik), parameter :: gmsh_edge(2, 12) = reshape([1_ik, 2_ik, 1_ik, 4_ik, 1_ik, 5_ik, &
      & 2_ik, 3_ik, 2_ik, 6_ik, 3_ik, 4_ik, 3_ik, 7_ik, 4_ik, 8_ik, 5_ik, 6_ik, 5_ik, 8_ik, &
      & 6_ik, 7_ik, 7_ik, 8_ik], [2, 12])
    integer(ik) :: k
    real(rk) :: x(3, 8)

    ! Corners.
    x = hex20_node_coo_ref(:, 1:8)
    call assert_close(reshape(x, [24]), reshape(corner_coo, [24]), 'gmsh order hex corners')
    call assert_true(dot_product(cross(x(:, 2) - x(:, 1), x(:, 4) - x(:, 1)), &
      & x(:, 5) - x(:, 1)) > 0.0_rk, 'gmsh order hex right-handed')

    ! Edge node 8 + k on the k-th Gmsh edge.
    do k = 1_ik, 12_ik
      call assert_true(all(edge_corners(:, 8_ik + k) == gmsh_edge(:, k)), &
        & 'gmsh order hex edge node')
    end do

  end subroutine check_gmsh_order

  !> Kronecker property, partition of unity, serendipity completeness, finite differences.
  !>
  !> Locals:
  !>   a, b, m, k, i_pt  integer        loop indices
  !>   xi, xi_h          real(rk)(3)    point, shifted point
  !>   shape             real(rk)(20)   shape functions
  !>   dshape, dfd       real(rk)(3,20) derivatives, finite-difference derivatives
  !>   delta             real(rk)(20)   Kronecker column
  !>   h                 real(rk)       finite-difference step
  !>   q_nodes           real(rk)(20)   monomial values at the nodes
  !>   q_val, q_grad     real(rk), real(rk)(3)  monomial value and gradient at xi
  !>   err_q22           real(rk)       largest interpolation error of xi1^2 xi2^2
  !>   expo              integer(3, 20) exponents of the serendipity monomials
  subroutine check_shape()

    integer, parameter :: expo(3, 20) = reshape([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, &
      & 2, 0, 0, 0, 2, 0, 0, 0, 2, 1, 1, 0, 0, 1, 1, 1, 0, 1, &
      & 2, 1, 0, 2, 0, 1, 1, 2, 0, 0, 2, 1, 1, 0, 2, 0, 1, 2, &
      & 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2], [3, 20])
    integer :: a, b, m, k, i_pt
    real(rk) :: xi(3), xi_h(3), shape(20), dshape(3, 20), dfd(3, 20), delta(20), h
    real(rk) :: q_nodes(20), q_val, q_grad(3), err_q22

    ! N_a(node b) = delta_ab.
    do b = 1, int(hex20_num_nodes)
      delta = 0.0_rk
      delta(b) = 1.0_rk
      call assert_close(hex20_shape(hex20_node_coo_ref(:, b)), delta, 'hex shape kronecker', &
        & abs_tol=1e-15_rk)
    end do

    ! Random points of [-1, 1]^3.
    h = 1e-5_rk
    err_q22 = 0.0_rk
    do i_pt = 1, 200
      call random_number(xi)
      xi = 2.0_rk*xi - 1.0_rk
      shape = hex20_shape(xi)
      dshape = hex20_shape_deriv(xi)

      ! Partition of unity.
      call assert_close(sum(shape), 1.0_rk, 'hex partition of unity', abs_tol=1e-14_rk)
      call assert_close(sum(dshape, dim=2), [0.0_rk, 0.0_rk, 0.0_rk], 'hex dshape sum zero', &
        & abs_tol=1e-13_rk)

      ! Serendipity completeness, monomial by monomial.
      do m = 1, 20
        do a = 1, int(hex20_num_nodes)
          call monomial(expo(:, m), hex20_node_coo_ref(:, a), q_nodes(a), q_grad)
        end do
        call monomial(expo(:, m), xi, q_val, q_grad)
        call assert_close(dot_product(shape, q_nodes), q_val, 'hex shape complete', &
          & abs_tol=1e-13_rk)
        call assert_close(matmul(dshape, q_nodes), q_grad, 'hex dshape complete', &
          & abs_tol=1e-13_rk)
      end do

      ! xi1^2 xi2^2 is outside the serendipity space.
      do a = 1, int(hex20_num_nodes)
        call monomial([2, 2, 0], hex20_node_coo_ref(:, a), q_nodes(a), q_grad)
      end do
      call monomial([2, 2, 0], xi, q_val, q_grad)
      err_q22 = max(err_q22, abs(dot_product(shape, q_nodes) - q_val))

      ! Central finite differences (cubic terms: truncation error ~ h^2).
      do k = 1, 3
        xi_h = xi
        xi_h(k) = xi(k) + h
        dfd(k, :) = hex20_shape(xi_h)
        xi_h(k) = xi(k) - h
        dfd(k, :) = (dfd(k, :) - hex20_shape(xi_h))/(2.0_rk*h)
      end do
      call assert_close(reshape(dshape, [60]), reshape(dfd, [60]), 'hex dshape finite difference', &
        & abs_tol=1e-8_rk)
    end do
    call assert_true(err_q22 > 1e-2_rk, 'hex shape not complete for xi1^2 xi2^2')

  end subroutine check_shape

  !> Geometry of every element of the Neper hex mesh of the unit cube.
  !>
  !> Arguments:
  !>   path  (in)  character  mesh file
  !> Locals:
  !>   mesh          type(mesh_file_t)  parsed mesh
  !>   i_elt, i_qpt  integer(ik)        loop indices
  !>   a             integer(ik)        element node
  !>   coo           real(rk)(3,20)     element node coordinates [mm]
  !>   jac           real(rk)(3,3)      jac(i, j) = d x_i / d xi_j [mm]
  !>   det_jac       real(rk)           det(jac) [mm^3]
  !>   vol, vol_box  real(rk)           element volume with the rule, of the box [mm^3]
  !>   vol_total     real(rk)           mesh volume [mm^3]
  !>   min_det       real(rk)           smallest det(jac)
  !>   max_mid_err   real(rk)           largest edge-node distance from its edge midpoint [mm]
  !>   max_vol_err   real(rk)           largest relative difference vol / vol_box - 1
  subroutine check_neper_mesh(path)

    character(len=*), intent(in) :: path

    type(mesh_file_t) :: mesh
    integer(ik) :: i_elt, i_qpt, a
    real(rk) :: coo(3, 20), jac(3, 3), det_jac, vol, vol_box, vol_total
    real(rk) :: min_det, max_mid_err, max_vol_err

    call mesh_file_load(mesh, path)
    call assert_equal(mesh%elt_type, elt_type_hex20, 'hex mesh element type')
    call assert_equal(mesh%num_elts, 125_ik, 'hex mesh num_elts')

    vol_total = 0.0_rk
    min_det = huge(1.0_rk)
    max_mid_err = 0.0_rk
    max_vol_err = 0.0_rk
    do i_elt = 1_ik, mesh%num_elts
      coo = mesh%coo(:, mesh%elt_nodes(:, i_elt))

      ! Edge nodes at the midpoints of their corners.
      do a = 9_ik, hex20_num_nodes
        max_mid_err = max(max_mid_err, maxval(abs(coo(:, a) - 0.5_rk*(coo(:, &
          & edge_corners(1, a)) + coo(:, edge_corners(2, a))))))
      end do

      ! Jacobian determinant at the quadrature points and integrated volume.
      vol = 0.0_rk
      do i_qpt = 1_ik, hex20_num_qpts
        jac = matmul(coo, transpose(hex20_shape_deriv(hex20_qpt_coo(:, i_qpt))))
        det_jac = determinant3(jac)
        min_det = min(min_det, det_jac)
        vol = vol + hex20_qpt_weight(i_qpt)*det_jac
      end do

      ! Volume of the box spanned by the edges from corner 1.
      vol_box = dot_product(cross(coo(:, 2) - coo(:, 1), coo(:, 4) - coo(:, 1)), &
        & coo(:, 5) - coo(:, 1))
      max_vol_err = max(max_vol_err, abs(vol/vol_box - 1.0_rk))
      vol_total = vol_total + vol
    end do

    call assert_close(max_mid_err, 0.0_rk, 'hex mesh edge nodes at midpoints', abs_tol=1e-11_rk)
    call assert_true(min_det > 0.0_rk, 'hex mesh det(J) > 0 at every quadrature point')
    call assert_close(max_vol_err, 0.0_rk, 'hex mesh element volume = box', abs_tol=1e-10_rk)
    call assert_close(vol_total, 1.0_rk, 'hex mesh volume 1', rel_tol=1e-12_rk)

  end subroutine check_neper_mesh

  !> Exact integral of x^p over [-1, 1].
  !>
  !> Arguments:
  !>   p  (in)  integer  exponent
  pure function integral_1d(p) result(val)

    integer, intent(in) :: p
    real(rk) :: val

    val = merge(2.0_rk/real(p + 1, rk), 0.0_rk, mod(p, 2) == 0)

  end function integral_1d

  !> Value and gradient of the monomial xi1^e1 xi2^e2 xi3^e3.
  !>
  !> Arguments:
  !>   expo  (in)   integer(3)   exponents
  !>   xi    (in)   real(rk)(3)  point
  !>   val   (out)  real(rk)     value
  !>   grad  (out)  real(rk)(3)  gradient
  !> Locals:
  !>   i, j  integer  coordinate indices
  subroutine monomial(expo, xi, val, grad)

    integer, intent(in) :: expo(3)
    real(rk), intent(in) :: xi(3)
    real(rk), intent(out) :: val, grad(3)

    integer :: i, j

    val = product(xi**expo)
    do i = 1, 3
      grad(i) = 0.0_rk
      if (expo(i) == 0) cycle
      grad(i) = real(expo(i), rk)*xi(i)**(expo(i) - 1)
      do j = 1, 3
        if (j /= i) grad(i) = grad(i)*xi(j)**expo(j)
      end do
    end do

  end subroutine monomial

  !> Cross product of two 3-vectors.
  !>
  !> Arguments:
  !>   u, v  (in)  real(rk)(3)  vectors
  pure function cross(u, v) result(w)

    real(rk), intent(in) :: u(3), v(3)
    real(rk) :: w(3)

    w = [u(2)*v(3) - u(3)*v(2), u(3)*v(1) - u(1)*v(3), u(1)*v(2) - u(2)*v(1)]

  end function cross

  !> The 9-point face rule.
  !>
  !> Locals:
  !>   weight_num     integer(3)   one-dimensional weights times 9
  !>   n, i, j, p, q  integer      point, tensor indices, exponents
  !>   approx, exact  real(rk)     rule value and exact integral
  !>   err6           real(rk)     largest error with an exponent 6
  subroutine check_face_quadrature()

    integer, parameter :: weight_num(3) = [5, 8, 5]
    integer :: n, i, j, p, q
    real(rk) :: approx, exact, err6

    call assert_equal(int(hex20_face_num_qpts), 9, 'face num_qpts')
    do j = 1, 3
      do i = 1, 3
        n = i + 3*(j - 1)
        call assert_close(hex20_face_qpt_coo(:, n), [hex20_qpt_coo(1, i), hex20_qpt_coo(1, j)], &
          & 'face qpt tensor order')
        call assert_close(hex20_face_qpt_weight(n), real(weight_num(i)*weight_num(j), rk)/81.0_rk, &
          & 'face qpt weight exact', rel_tol=4e-16_rk)
      end do
    end do
    call assert_close(sum(hex20_face_qpt_weight), 4.0_rk, 'face weights sum 4', rel_tol=1e-15_rk)

    err6 = 0.0_rk
    do p = 0, 6
      do q = 0, 6
        approx = 0.0_rk
        do n = 1, int(hex20_face_num_qpts)
          approx = approx + hex20_face_qpt_weight(n)*hex20_face_qpt_coo(1, n)**p* &
            & hex20_face_qpt_coo(2, n)**q
        end do
        exact = integral_1d(p)*integral_1d(q)
        if (max(p, q) <= 5) then
          call assert_close(approx, exact, 'face qpt exact', abs_tol=1e-14_rk)
        else
          err6 = max(err6, abs(approx - exact))
        end if
      end do
    end do
    call assert_true(err6 > 1e-3_rk, 'face qpt not exact for exponent 6')

  end subroutine check_face_quadrature

  !> Face shape functions: node order, Kronecker property, partition of unity, completeness,
  !> finite differences, and the restriction of the hex20 functions to the face xi3 = -1.
  !>
  !> Locals:
  !>   coo_ref        real(rk)(2,8)   expected face node coordinates (Gmsh Quadrangle8)
  !>   face_nodes     integer(ik)(8)  hex20 nodes on the face xi3 = -1, in face order
  !>   expo           integer(3,8)    exponents of the 2D serendipity monomials (third 0)
  !>   a, b, m, k, i_pt  integer      loop indices
  !>   xi, xi_h       real(rk)(2)     point, shifted point
  !>   shape, delta   real(rk)(8)     face shape functions, Kronecker column
  !>   dshape, dfd    real(rk)(2,8)   derivatives, finite-difference derivatives
  !>   q_nodes        real(rk)(8)     monomial values at the face nodes
  !>   q_val          real(rk)        monomial value
  !>   q_grad         real(rk)(3)     monomial gradient (3D form)
  !>   vol_shape      real(rk)(20)    hex20 shape functions on the face
  !>   on_face        logical(20)     whether each hex20 node is on the face
  !>   err_q22        real(rk)        largest interpolation error of x^2 y^2
  !>   h              real(rk)        finite-difference step
  subroutine check_face_shape()

    real(rk), parameter :: coo_ref(2, 8) = reshape([-1.0_rk, -1.0_rk, 1.0_rk, -1.0_rk, 1.0_rk, &
      & 1.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, -1.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, -1.0_rk, &
      & 0.0_rk], [2, 8])
    integer(ik), parameter :: face_nodes(8) = [1_ik, 2_ik, 3_ik, 4_ik, 9_ik, 12_ik, 14_ik, 10_ik]
    integer, parameter :: expo(3, 8) = reshape([0, 0, 0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 1, 1, 0, &
      & 0, 2, 0, 2, 1, 0, 1, 2, 0], [3, 8])
    integer :: a, b, m, k, i_pt
    real(rk) :: xi(2), xi_h(2), shape(8), delta(8), dshape(2, 8), dfd(2, 8), q_nodes(8)
    real(rk) :: q_val, q_grad(3), vol_shape(20), err_q22, h
    logical :: on_face(20)

    ! Node order (Gmsh Quadrangle8).
    call assert_equal(int(hex20_face_num_nodes), 8, 'face num_nodes')
    call assert_close(reshape(hex20_face_node_coo_ref, [16]), reshape(coo_ref, [16]), &
      & 'face node coordinates')

    ! N_a(node b) = delta_ab.
    do b = 1, 8
      delta = 0.0_rk
      delta(b) = 1.0_rk
      call assert_close(hex20_face_shape(hex20_face_node_coo_ref(:, b)), delta, &
        & 'face shape kronecker', abs_tol=1e-15_rk)
    end do

    on_face = .false.
    on_face(face_nodes) = .true.
    h = 1e-5_rk
    err_q22 = 0.0_rk
    do i_pt = 1, 200
      call random_number(xi)
      xi = 2.0_rk*xi - 1.0_rk
      shape = hex20_face_shape(xi)
      dshape = hex20_face_shape_deriv(xi)

      ! Partition of unity.
      call assert_close(sum(shape), 1.0_rk, 'face partition of unity', abs_tol=1e-14_rk)
      call assert_close(sum(dshape, dim=2), [0.0_rk, 0.0_rk], 'face dshape sum zero', &
        & abs_tol=1e-13_rk)

      ! Serendipity completeness, and x^2 y^2 not reproduced.
      do m = 1, 8
        do a = 1, 8
          call monomial(expo(:, m), [hex20_face_node_coo_ref(:, a), 0.0_rk], q_nodes(a), q_grad)
        end do
        call monomial(expo(:, m), [xi, 0.0_rk], q_val, q_grad)
        call assert_close(dot_product(shape, q_nodes), q_val, 'face shape complete', &
          & abs_tol=1e-13_rk)
        call assert_close(matmul(dshape, q_nodes), q_grad(1:2), 'face dshape complete', &
          & abs_tol=1e-13_rk)
      end do
      do a = 1, 8
        call monomial([2, 2, 0], [hex20_face_node_coo_ref(:, a), 0.0_rk], q_nodes(a), q_grad)
      end do
      call monomial([2, 2, 0], [xi, 0.0_rk], q_val, q_grad)
      err_q22 = max(err_q22, abs(dot_product(shape, q_nodes) - q_val))

      ! Central finite differences.
      do k = 1, 2
        xi_h = xi
        xi_h(k) = xi(k) + h
        dfd(k, :) = hex20_face_shape(xi_h)
        xi_h(k) = xi(k) - h
        dfd(k, :) = (dfd(k, :) - hex20_face_shape(xi_h))/(2.0_rk*h)
      end do
      call assert_close(reshape(dshape, [16]), reshape(dfd, [16]), &
        & 'face dshape finite difference', abs_tol=1e-8_rk)

      ! Restriction of the hex20 functions to the face xi3 = -1.
      if (i_pt <= 100) then
        vol_shape = hex20_shape([xi, -1.0_rk])
        call assert_close(vol_shape(face_nodes), shape, 'face = hex20 on xi3 = -1', &
          & abs_tol=1e-15_rk)
        call assert_true(all(abs(pack(vol_shape, .not. on_face)) < 1e-15_rk), &
          & 'hex20 off-face nodes vanish on xi3 = -1')
      end if
    end do
    call assert_true(err_q22 > 1e-2_rk, 'face shape not complete for x^2 y^2')

  end subroutine check_face_shape

  !> Compare the face element with the reference fixture.
  !>
  !> Arguments:
  !>   path  (in)  character  fixture file
  !> Locals:
  !>   unit, ios, num, i_pt  integer  file unit, I/O status, point count, point index
  !>   word          character(16)   section name
  !>   line          real(rk)(26)    xi1 xi2, 8 values, 8 d/dxi1, 8 d/dxi2
  !>   xi            real(rk)(2)     reference coordinates of the point
  !>   dshape        real(rk)(2,8)   FEPX derivatives
  subroutine check_face_reference(path)

    character(len=*), intent(in) :: path

    integer :: unit, ios, num, i_pt
    character(len=16) :: word
    real(rk) :: line(26), xi(2), dshape(2, 8)

    open (newunit=unit, file=path, status='old', action='read', iostat=ios)
    call assert_equal(ios, 0, 'reference quad8 fixture open')
    if (ios /= 0) return
    read (unit, *) word, num
    call assert_equal(trim(word), 'quad_points', 'reference quad8 fixture header')
    call assert_equal(num, 12, 'reference quad8 fixture num points')
    do i_pt = 1, num
      read (unit, *) line
      xi = line(1:2)
      dshape = hex20_face_shape_deriv(xi)
      call assert_close(hex20_face_shape(xi), line(3:10), 'reference face shape', &
        & abs_tol=1e-15_rk)
      call assert_close(dshape(1, :), line(11:18), 'reference face dshape/dxi1', &
        & abs_tol=1e-15_rk)
      call assert_close(dshape(2, :), line(19:26), 'reference face dshape/dxi2', &
        & abs_tol=1e-15_rk)
    end do
    close (unit)

  end subroutine check_face_reference

end program test_hex20
