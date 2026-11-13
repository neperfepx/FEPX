! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> 20-node (quadratic serendipity) hexahedron: element node order (Gmsh's), reference node
!> coordinates, the shape functions and their derivatives with respect to the reference
!> coordinates, and the element's quadrature rule (each element module owns its rule: ROADMAP
!> D36).
!>
!> Theory : dev/theory/addendum_hex20.md (approved): §1 node order and reference element
!>          [-1, 1]^3, §2 serendipity shape functions (Zienkiewicz, Taylor & Zhu 2013; Bathe
!>          2014; their span is the serendipity space S_2 of Arnold & Awanou 2011). For node a
!>          with reference coordinates xi_a:
!>            corner (all |xi_ia| = 1):
!>              N_a = 1/8 (1 + xi1 xi1a)(1 + xi2 xi2a)(1 + xi3 xi3a)(xi1 xi1a + xi2 xi2a
!>                    + xi3 xi3a - 2)
!>            edge node with xi_ka = 0 (the other two coordinates i, j = +-1):
!>              N_a = 1/4 (1 - xi_k^2)(1 + xi_i xi_ia)(1 + xi_j xi_ja)
!>          §3 quadrature: tensor product of the 3-point Gauss-Legendre rule on [-1, 1]
!>          (Abramowitz & Stegun 1964, Eq. 25.4.29: points 0 and +-sqrt(3/5), weights 8/9 and
!>          5/9); exact for xi1^p xi2^q xi3^r with p, q, r <= 5; weights sum to 8 (the reference
!>          volume); point order xi1 fastest, then xi2, then xi3.
!> Source : new in FEPX (v2 has no 20-node hexahedron); node order: Gmsh reference manual
!>          (Gmsh 4.15.2) §10.2.1, Hexahedron20 diagram, checked by test_hex20 on a Neper mesh.
!>
!> Element node order = Gmsh order (ROADMAP D35; conventions §3a), 20-node hexahedron (Gmsh
!> type 17). Reference coordinates xi = (xi1, xi2, xi3) in [-1, 1]^3 along Gmsh's u (corner
!> 1 -> 2), v (1 -> 4) and w (1 -> 5); the element seen from above (looking down the xi3 axis),
!> one layer at a time:
!>
!>   Layer xi3 = -1 (bottom)        Layer xi3 = 0 (middle)      Layer xi3 = +1 (top)
!>
!>   xi2                            xi2                         xi2
!>    ^                              ^                           ^
!>    4-----14-----3                16-----------15              8-----20-----7
!>    |            |                 |            |              |            |
!>    |            |                 |            |              |            |
!>   10           12                 |            |             18           19
!>    |            |                 |            |              |            |
!>    |            |                 |            |              |            |
!>    1------9-----2 --> xi1        11-----------13 --> xi1      5-----17-----6 --> xi1
!>
!>   corners 1-4 COUNTERCLOCKWISE   edge nodes of the vertical   corners 5-8 above 1-4
!>   edge nodes 9, 10, 12, 14       edges: 11 (1-5), 13 (2-6),   edge nodes 17, 18, 19, 20
!>                                  15 (3-7), 16 (4-8)
!>
!>   node  xi              node  xi            on edge   node  xi            on edge
!>    1   (-1, -1, -1)       9   ( 0, -1, -1)    1-2      15   ( 1,  1,  0)    3-7
!>    2   ( 1, -1, -1)      10   (-1,  0, -1)    1-4      16   (-1,  1,  0)    4-8
!>    3   ( 1,  1, -1)      11   (-1, -1,  0)    1-5      17   ( 0, -1,  1)    5-6
!>    4   (-1,  1, -1)      12   ( 1,  0, -1)    2-3      18   (-1,  0,  1)    5-8
!>    5   (-1, -1,  1)      13   ( 1, -1,  0)    2-6      19   ( 1,  0,  1)    6-7
!>    6   ( 1, -1,  1)      14   ( 0,  1, -1)    3-4      20   ( 0,  1,  1)    7-8
!>    7   ( 1,  1,  1)
!>    8   (-1,  1,  1)
!>
!> Numbering rule:
!>   - corners 1-8 first: the bottom face 1, 2, 3, 4 counterclockwise seen from above
!>     (right-hand rule: (x2 - x1) x (x4 - x1) points towards x5), so det(J) > 0 for a valid
!>     element; then the top face 5-8, each above the corner numbered 4 less;
!>   - then one node per edge, in Gmsh's edge order "by lowest corner, then by the other
!>     corner": 1-2, 1-4, 1-5, 2-3, 2-6, 3-4, 3-7, 4-8, 5-6, 5-8, 6-7, 7-8. This is NOT the
!>     Abaqus/CalculiX C3D20 order (bottom face around, top face around, then the verticals).
!>
!> Face element (the boundary faces of a hex20 mesh; ROADMAP D36): 8-node serendipity
!> quadrangle, Gmsh order (Gmsh type 16; Gmsh reference manual §10.2.1, Quadrangle8), reference
!> coordinates (xi1, xi2) in [-1, 1]^2:
!>
!>   xi2
!>    ^
!>    4-----7-----3
!>    |           |
!>    |           |
!>    8           6
!>    |           |
!>    |           |
!>    1-----5-----2 --> xi1
!>
!>    node  xi                node  xi         on edge
!>     1   (-1, -1) corner      5   ( 0, -1)    1-2
!>     2   ( 1, -1) corner      6   ( 1,  0)    2-3
!>     3   ( 1,  1) corner      7   ( 0,  1)    3-4
!>     4   (-1,  1) corner      8   (-1,  0)    4-1
!>
!> Corners 1-4 counterclockwise in (xi1, xi2), then the edge nodes of edges 1-2, 2-3, 3-4, 4-1.
!> The faces of a hex20 are such quadrangles: e.g. the face xi3 = -1 is element nodes 1, 2, 3, 4,
!> 9, 12, 14, 10 (in this order its normal points into the element). Shape functions: the 2D
!> serendipity functions (addendum_hex20.md H2-H4), corner N_a = 1/4 (1 + xi1 xi1a)(1 + xi2 xi2a)
!> (xi1 xi1a + xi2 xi2a - 1), edge node with xi_ka = 0: N_a = 1/2 (1 - xi_k^2)(1 + xi_j xi_ja);
!> they equal v2's shape_2d_mod.f90::sf2d08 (unused in v2; on [0, 1]^2, xi = 2 s - 1, v2's node
!> order 1 2 3 4 5 6 7 8 = Gmsh 1 5 2 6 3 7 4 8), cross-checked by test_hex20. Neper writes no
!> $Fasets for hex meshes (2026-10-09), so these faces are not read from files yet; the reader's
!> reversal of face lines (mesh_file_mod) is assumed to give this order, with the normal
!> pointing out of the body, as for triangles (ROADMAP D37, Q12; to be checked when Neper
!> writes them). Face quadrature: 3 x 3 Gauss-Legendre (addendum_hex20.md §3),
!> xi1 fastest; weights sum to 4.
!>
!> Public parameters:
!>   hex20_gmsh_type     integer(ik)            Gmsh element type code (17)
!>   hex20_num_nodes     integer(ik)            number of nodes (20)
!>   hex20_node_coo_ref  real(rk)(3, 20)        reference coordinates of the nodes [-]
!>   hex20_num_qpts      integer(ik)            number of quadrature points (27)
!>   hex20_qpt_centroid  integer(ik)            the point at the element centre (0, 0, 0) (14):
!>                                              elemental output (v2 manual, output.rst;
!>                                              dev/cfg_keywords.md §6)
!>   hex20_qpt_coo       real(rk)(3, 27)        reference coordinates of the quadrature points [-]
!>   hex20_qpt_weight    real(rk)(27)           quadrature weights (sum 8) [-]
!>   hex20_face_num_nodes     integer(ik)       number of face element nodes (8)
!>   hex20_face_node_coo_ref  real(rk)(2, 8)    reference coordinates of the face nodes [-]
!>   hex20_face_num_qpts      integer(ik)       number of face quadrature points (9)
!>   hex20_face_qpt_coo       real(rk)(2, 9)    face reference coordinates of the points [-]
!>   hex20_face_qpt_weight    real(rk)(9)       face quadrature weights (sum 4) [-]
!>
!> Public routines:
!>   hex20_shape        shape functions at one point
!>   hex20_shape_deriv  their derivatives with respect to the reference coordinates at one point
!>   hex20_face_shape        face element (8-node quadrangle) shape functions at one point
!>   hex20_face_shape_deriv  their derivatives with respect to the face reference coordinates
module hex20_mod

  use kinds_mod, only: rk, ik

  implicit none

  private

  public :: hex20_gmsh_type, hex20_num_nodes, hex20_node_coo_ref
  public :: hex20_num_qpts, hex20_qpt_coo, hex20_qpt_weight, hex20_qpt_centroid
  public :: hex20_shape, hex20_shape_deriv
  public :: hex20_face_num_nodes, hex20_face_node_coo_ref
  public :: hex20_face_num_qpts, hex20_face_qpt_coo, hex20_face_qpt_weight
  public :: hex20_face_shape, hex20_face_shape_deriv

  ! Gmsh element type code.
  integer(ik), parameter :: hex20_gmsh_type = 17_ik

  ! Number of nodes, and number of corners (the first nodes).
  integer(ik), parameter :: hex20_num_nodes = 20_ik
  integer(ik), parameter :: num_corners = 8_ik

  ! Reference coordinates of the nodes (one column per node, Gmsh order).
  real(rk), parameter :: hex20_node_coo_ref(3, hex20_num_nodes) = reshape([ &
    ! corners 1-8
    & -1.0_rk, -1.0_rk, -1.0_rk, &
    & 1.0_rk, -1.0_rk, -1.0_rk, &
    & 1.0_rk, 1.0_rk, -1.0_rk, &
    & -1.0_rk, 1.0_rk, -1.0_rk, &
    & -1.0_rk, -1.0_rk, 1.0_rk, &
    & 1.0_rk, -1.0_rk, 1.0_rk, &
    & 1.0_rk, 1.0_rk, 1.0_rk, &
    & -1.0_rk, 1.0_rk, 1.0_rk, &
    ! edge nodes 9-20 (edges 1-2, 1-4, 1-5, 2-3, 2-6, 3-4, 3-7, 4-8, 5-6, 5-8, 6-7, 7-8)
    & 0.0_rk, -1.0_rk, -1.0_rk, &
    & -1.0_rk, 0.0_rk, -1.0_rk, &
    & -1.0_rk, -1.0_rk, 0.0_rk, &
    & 1.0_rk, 0.0_rk, -1.0_rk, &
    & 1.0_rk, -1.0_rk, 0.0_rk, &
    & 0.0_rk, 1.0_rk, -1.0_rk, &
    & 1.0_rk, 1.0_rk, 0.0_rk, &
    & -1.0_rk, 1.0_rk, 0.0_rk, &
    & 0.0_rk, -1.0_rk, 1.0_rk, &
    & -1.0_rk, 0.0_rk, 1.0_rk, &
    & 1.0_rk, 0.0_rk, 1.0_rk, &
    & 0.0_rk, 1.0_rk, 1.0_rk], [3, int(hex20_num_nodes)])

  ! For each edge node (9-20), the reference coordinate that is 0 on its edge (the edge
  ! direction): 1, 2 or 3.
  integer(ik), parameter :: edge_dir(num_corners + 1:hex20_num_nodes) = &
    & [1_ik, 2_ik, 3_ik, 2_ik, 3_ik, 1_ik, 3_ik, 3_ik, 1_ik, 2_ik, 2_ik, 1_ik]

  ! Quadrature rule (3 x 3 x 3 Gauss-Legendre, addendum_hex20.md §3).

  ! Number of quadrature points.
  integer(ik), parameter :: hex20_num_qpts = 27_ik

  ! The quadrature point at the element centre (0, 0, 0) (elemental output).
  integer(ik), parameter :: hex20_qpt_centroid = 14_ik

  ! 3-point Gauss-Legendre rule on [-1, 1]: points -sqrt(3/5), 0, +sqrt(3/5) and their weights
  ! 5/9, 8/9, 5/9 (sqrt(3/5) to 34 digits; the compiler rounds it to the nearest real).
  real(rk), parameter :: gauss3_point = 0.7745966692414833770358530799564799_rk
  real(rk), parameter :: gauss3_coo(3) = [-gauss3_point, 0.0_rk, gauss3_point]
  real(rk), parameter :: gauss3_weight(3) = [5.0_rk/9.0_rk, 8.0_rk/9.0_rk, 5.0_rk/9.0_rk]

  ! Tensor-product indices of each quadrature point: point n uses gauss3_coo(i) along xi1,
  ! gauss3_coo(j) along xi2 and gauss3_coo(k) along xi3, n = i + 3(j - 1) + 9(k - 1)
  ! (xi1 fastest).
  integer, parameter :: hex_i(27) = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, &
    & 1, 2, 3, 1, 2, 3, 1, 2, 3]
  integer, parameter :: hex_j(27) = [1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, &
    & 1, 1, 1, 2, 2, 2, 3, 3, 3]
  integer, parameter :: hex_k(27) = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
    & 3, 3, 3, 3, 3, 3, 3, 3, 3]

  ! Quadrature points (one column per point).
  real(rk), parameter :: hex20_qpt_coo(3, hex20_num_qpts) = transpose(reshape( &
    & [gauss3_coo(hex_i), gauss3_coo(hex_j), gauss3_coo(hex_k)], [int(hex20_num_qpts), 3]))

  ! Weights: products of the one-dimensional weights.
  real(rk), parameter :: hex20_qpt_weight(hex20_num_qpts) = &
    & gauss3_weight(hex_i)*gauss3_weight(hex_j)*gauss3_weight(hex_k)

  ! Face element (8-node quadrangle): number of nodes, of corners, reference coordinates (Gmsh
  ! order) and, for edge nodes 5-8, the coordinate that is 0 on the edge (1 or 2).
  integer(ik), parameter :: hex20_face_num_nodes = 8_ik
  integer(ik), parameter :: num_face_corners = 4_ik
  real(rk), parameter :: hex20_face_node_coo_ref(2, hex20_face_num_nodes) = reshape([ &
    & -1.0_rk, -1.0_rk, &
    & 1.0_rk, -1.0_rk, &
    & 1.0_rk, 1.0_rk, &
    & -1.0_rk, 1.0_rk, &
    & 0.0_rk, -1.0_rk, &
    & 1.0_rk, 0.0_rk, &
    & 0.0_rk, 1.0_rk, &
    & -1.0_rk, 0.0_rk], [2, int(hex20_face_num_nodes)])
  integer(ik), parameter :: face_edge_dir(num_face_corners + 1:hex20_face_num_nodes) = &
    & [1_ik, 2_ik, 1_ik, 2_ik]

  ! Face quadrature rule: 3 x 3 tensor product of the same one-dimensional rule, point
  ! i + 3(j - 1) (xi1 fastest).
  integer(ik), parameter :: hex20_face_num_qpts = 9_ik
  integer, parameter :: quad_i(9) = [1, 2, 3, 1, 2, 3, 1, 2, 3]
  integer, parameter :: quad_j(9) = [1, 1, 1, 2, 2, 2, 3, 3, 3]
  real(rk), parameter :: hex20_face_qpt_coo(2, hex20_face_num_qpts) = transpose(reshape( &
    & [gauss3_coo(quad_i), gauss3_coo(quad_j)], [int(hex20_face_num_qpts), 2]))
  real(rk), parameter :: hex20_face_qpt_weight(hex20_face_num_qpts) = &
    & gauss3_weight(quad_i)*gauss3_weight(quad_j)

contains

  !> Shape functions of the 20-node hexahedron at one point.
  !>
  !> Theory : addendum_hex20.md §2 (formulas in the module header).
  !> Source : new in FEPX (no v2 hexahedron).
  !>
  !> Arguments:
  !>   xi     (in)  real(rk)(3)   reference coordinates of the point [-]
  !>   shape  (out) real(rk)(20)  shape function of each element node (Gmsh order) at the point [-]
  !> Locals:
  !>   a        integer(ik)   node
  !>   k        integer(ik)   edge direction of an edge node
  !>   lin      real(rk)(3)   1 + xi_i xi_ia for each coordinate i [-]
  pure function hex20_shape(xi) result(shape)

    real(rk), intent(in) :: xi(3)
    real(rk) :: shape(hex20_num_nodes)

    integer(ik) :: a, k
    real(rk) :: lin(3)

    ! Corners: 1/8 (1 + xi1 xi1a)(1 + xi2 xi2a)(1 + xi3 xi3a)(xi . xi_a - 2).
    do a = 1_ik, num_corners
      lin = 1.0_rk + xi*hex20_node_coo_ref(:, a)
      shape(a) = 0.125_rk*lin(1)*lin(2)*lin(3)* &
        & (dot_product(xi, hex20_node_coo_ref(:, a)) - 2.0_rk)
    end do

    ! Edge nodes: 1/4 (1 - xi_k^2) times the two other linear factors.
    do a = num_corners + 1_ik, hex20_num_nodes
      k = edge_dir(a)
      lin = 1.0_rk + xi*hex20_node_coo_ref(:, a)
      lin(k) = 1.0_rk - xi(k)**2
      shape(a) = 0.25_rk*lin(1)*lin(2)*lin(3)
    end do

  end function hex20_shape

  !> Derivatives of the shape functions of the 20-node hexahedron with respect to the reference
  !> coordinates, at one point.
  !>
  !> Theory : addendum_hex20.md §2, differentiated:
  !>            corner: dN_a/dxi_i = 1/8 xi_ia (prod over j /= i of (1 + xi_j xi_ja))
  !>                                 (2 xi_i xi_ia + sum over j /= i of xi_j xi_ja - 1)
  !>            edge node (direction k): dN_a/dxi_k = -1/2 xi_k (prod over j /= k of
  !>                                 (1 + xi_j xi_ja)); for i /= k: dN_a/dxi_i =
  !>                                 1/4 (1 - xi_k^2) xi_ia (1 + xi_j xi_ja), j the third index.
  !> Source : new in FEPX (no v2 hexahedron).
  !>
  !> Arguments:
  !>   xi      (in)  real(rk)(3)      reference coordinates of the point [-]
  !>   dshape  (out) real(rk)(3, 20)  dshape(i, a) = d N_a / d xi_i at the point [-]
  !> Locals:
  !>   a, i, k  integer(ik)   node, coordinate, edge direction
  !>   coo_a    real(rk)(3)   reference coordinates of node a
  !>   lin      real(rk)(3)   1 + xi_i xi_ia for each coordinate i [-]
  !>   dlin     real(rk)(3)   d lin_i / d xi_i [-]
  !>   proj     real(rk)      xi . xi_a [-]
  pure function hex20_shape_deriv(xi) result(dshape)

    real(rk), intent(in) :: xi(3)
    real(rk) :: dshape(3, hex20_num_nodes)

    integer(ik) :: a, i, k
    real(rk) :: coo_a(3), lin(3), dlin(3), proj

    ! Corners: product rule on lin(1) lin(2) lin(3) (proj - 2).
    do a = 1_ik, num_corners
      coo_a = hex20_node_coo_ref(:, a)
      lin = 1.0_rk + xi*coo_a
      proj = dot_product(xi, coo_a)
      do i = 1_ik, 3_ik
        ! d/dxi_i: xi_ia (product of the other two factors) ((proj - 2) + lin_i) / 8.
        dshape(i, a) = 0.125_rk*coo_a(i)*product(lin, mask=[1, 2, 3] /= i)* &
          & (proj - 2.0_rk + lin(i))
      end do
    end do

    ! Edge nodes: the factor in the edge direction k is 1 - xi_k^2, the others are linear.
    do a = num_corners + 1_ik, hex20_num_nodes
      coo_a = hex20_node_coo_ref(:, a)
      k = edge_dir(a)
      lin = 1.0_rk + xi*coo_a
      dlin = coo_a
      lin(k) = 1.0_rk - xi(k)**2
      dlin(k) = -2.0_rk*xi(k)
      do i = 1_ik, 3_ik
        ! d/dxi_i: derivative of factor i times the other two, over 4.
        dshape(i, a) = 0.25_rk*dlin(i)*product(lin, mask=[1, 2, 3] /= i)
      end do
    end do

  end function hex20_shape_deriv

  !> Shape functions of the face element (8-node serendipity quadrangle) at one point.
  !>
  !> Theory : module header (face element; addendum_hex20.md H2-H4).
  !> Source : equal to v2 shape_2d_mod.f90::sf2d08 under xi = 2 s - 1 and Gmsh node numbering
  !>          (written here in the [-1, 1]^2 form, like hex20_shape).
  !>
  !> Arguments:
  !>   xi     (in)  real(rk)(2)  face reference coordinates of the point [-]
  !>   shape  (out) real(rk)(8)  shape function of each face node (Gmsh order) at the point [-]
  !> Locals:
  !>   a        integer(ik)   node
  !>   k        integer(ik)   edge direction of an edge node
  !>   lin      real(rk)(2)   1 + xi_i xi_ia for each coordinate i [-]
  pure function hex20_face_shape(xi) result(shape)

    real(rk), intent(in) :: xi(2)
    real(rk) :: shape(hex20_face_num_nodes)

    integer(ik) :: a, k
    real(rk) :: lin(2)

    ! Corners: 1/4 (1 + xi1 xi1a)(1 + xi2 xi2a)(xi . xi_a - 1).
    do a = 1_ik, num_face_corners
      lin = 1.0_rk + xi*hex20_face_node_coo_ref(:, a)
      shape(a) = 0.25_rk*lin(1)*lin(2)*(dot_product(xi, hex20_face_node_coo_ref(:, a)) - 1.0_rk)
    end do

    ! Edge nodes: 1/2 (1 - xi_k^2) times the other linear factor.
    do a = num_face_corners + 1_ik, hex20_face_num_nodes
      k = face_edge_dir(a)
      lin = 1.0_rk + xi*hex20_face_node_coo_ref(:, a)
      lin(k) = 1.0_rk - xi(k)**2
      shape(a) = 0.5_rk*lin(1)*lin(2)
    end do

  end function hex20_face_shape

  !> Derivatives of the face element (8-node quadrangle) shape functions with respect to the face
  !> reference coordinates, at one point.
  !>
  !> Theory : derivatives of the functions of hex20_face_shape:
  !>            corner: dN_a/dxi_i = 1/4 xi_ia (1 + xi_j xi_ja)(xi . xi_a - 1 + 1 + xi_i xi_ia),
  !>                    j the other index;
  !>            edge node (direction k): dN_a/dxi_k = -xi_k (1 + xi_j xi_ja);
  !>                    dN_a/dxi_j = 1/2 (1 - xi_k^2) xi_ja.
  !> Source : equal to v2 shape_2d_mod.f90::sf2dg08 under xi = 2 s - 1 (factor 1/2) and Gmsh
  !>          node numbering.
  !>
  !> Arguments:
  !>   xi      (in)  real(rk)(2)     face reference coordinates of the point [-]
  !>   dshape  (out) real(rk)(2, 8)  dshape(i, a) = d N_a / d xi_i at the point [-]
  !> Locals:
  !>   a, k     integer(ik)   node, edge direction
  !>   coo_a    real(rk)(2)   reference coordinates of node a
  !>   lin      real(rk)(2)   1 + xi_i xi_ia for each coordinate i [-]
  !>   dlin     real(rk)(2)   d lin_i / d xi_i [-]
  !>   proj     real(rk)      xi . xi_a [-]
  pure function hex20_face_shape_deriv(xi) result(dshape)

    real(rk), intent(in) :: xi(2)
    real(rk) :: dshape(2, hex20_face_num_nodes)

    integer(ik) :: a, k
    real(rk) :: coo_a(2), lin(2), dlin(2), proj

    ! Corners: product rule on lin(1) lin(2) (proj - 1).
    do a = 1_ik, num_face_corners
      coo_a = hex20_face_node_coo_ref(:, a)
      lin = 1.0_rk + xi*coo_a
      proj = dot_product(xi, coo_a)
      dshape(1, a) = 0.25_rk*coo_a(1)*lin(2)*(proj - 1.0_rk + lin(1))
      dshape(2, a) = 0.25_rk*coo_a(2)*lin(1)*(proj - 1.0_rk + lin(2))
    end do

    ! Edge nodes: the factor in the edge direction k is 1 - xi_k^2, the other is linear.
    do a = num_face_corners + 1_ik, hex20_face_num_nodes
      coo_a = hex20_face_node_coo_ref(:, a)
      k = face_edge_dir(a)
      lin = 1.0_rk + xi*coo_a
      dlin = coo_a
      lin(k) = 1.0_rk - xi(k)**2
      dlin(k) = -2.0_rk*xi(k)
      dshape(1, a) = 0.5_rk*dlin(1)*lin(2)
      dshape(2, a) = 0.5_rk*dlin(2)*lin(1)
    end do

  end function hex20_face_shape_deriv

end module hex20_mod
