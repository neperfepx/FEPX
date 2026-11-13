! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> 10-node (quadratic) tetrahedron: element node order (Gmsh's), reference node coordinates, the
!> shape functions and their derivatives with respect to the reference coordinates, and the
!> element's quadrature rule (each element module owns its rule: ROADMAP D36).
!>
!> Theory : standard isoparametric quadratic (Lagrange) tetrahedron. With the volume
!>          (barycentric) coordinates L = (1 - xi1 - xi2 - xi3, xi1, xi2, xi3), a corner node
!>          at L_i = 1 has N = L_i (2 L_i - 1) and an edge node between corners i and j has
!>          N = 4 L_i L_j. The functions are the unique quadratic polynomials with
!>          N_a(node b) = delta_ab (checked by test_tet10).
!>          Quadrature: 15-point rule of degree 5, Keast, P. (1986), "Moderate-degree tetrahedral
!>          quadrature formulas", Comput. Methods Appl. Mech. Engrg. 55:339-348 (as cited by
!>          v2); it integrates every polynomial of degree <= 5 exactly. Its points form four
!>          orbits of the tetrahedron's symmetry group (barycentric coordinates L):
!>            4 points  L = permutations of (1/3, 1/3, 1/3, 0)       weight w1
!>            1 point   L = (1/4, 1/4, 1/4, 1/4)                     weight w2
!>            4 points  L = permutations of (1/11, 1/11, 1/11, 8/11) weight w3
!>            6 points  L = permutations of (a, a, b, b), a + b = 1/2 weight w4
!>          The weights include the reference volume 1/6 (they sum to 1/6), so that the
!>          integral over the element of f = sum over i_qpt of tet10_qpt_weight(i_qpt) f det(J).
!> Source : port of v2 shape_3d_mod.f90::t10_shape_hpar (shape functions) and t10_deriv_hpar
!>          (derivatives): the same expressions, renumbered from v2's node order to Gmsh's
!>          (ROADMAP D35); quadrature: port of v2 quadrature_mod.f90::quadrature_init (3D part:
!>          qploc, wtqp), the same values typed as in v2, in v2's point order (quadrature-point
!>          data of FEPX and v2 are then point-for-point comparable). Node order: Gmsh
!>          reference manual, "Node ordering", 10-node tetrahedron (element type 11); checked by
!>          test_tet10 on Neper meshes.
!>
!> Element node order = Gmsh order (ROADMAP D35; conventions §3a), 10-node tetrahedron
!> (Gmsh type 11). FEPX uses the node order of the mesh file internally (v2 renumbered the
!> nodes when reading). Reference coordinates xi = (xi1, xi2, xi3); the element seen from above
!> (looking down the xi3 axis), one layer at a time:
!>
!>   Layer xi3 = 0 (base face)      Layer xi3 = 1/2            Layer xi3 = 1 (apex)
!>
!>   xi2                            xi2                        xi2
!>    ^                              ^                          ^
!>    3                              |                          |
!>    | \                            9                          |
!>    |   \                          | \                        |
!>    7     6                        |   \                      |
!>    |       \                      |     \                    |
!>    |         \                    |       \                  |
!>    1-----5-----2 --> xi1          8-------10 --> xi1         4 --> xi1
!>
!>   5, 6, 7 COUNTERCLOCKWISE        8, 9, 10 CLOCKWISE
!>   (midpoints of 1-2, 2-3, 3-1)    (above corners 1, 3, 2:
!>                                    midpoints of 4-1, 4-3, 4-2)
!>
!>   node  xi                 node  xi                on edge
!>    1   (0, 0, 0) corner      5   (1/2, 0,   0  )    1-2
!>    2   (1, 0, 0) corner      6   (1/2, 1/2, 0  )    2-3
!>    3   (0, 1, 0) corner      7   (0,   1/2, 0  )    3-1
!>    4   (0, 0, 1) corner      8   (0,   0,   1/2)    4-1
!>                              9   (0,   1/2, 1/2)    4-3
!>                             10   (1/2, 0,   1/2)    4-2
!>
!> Numbering rule:
!>   - corners 1-4 first; 1, 2, 3 are counterclockwise seen from 4 (right-hand rule:
!>     (x2 - x1) x (x3 - x1) points towards x4), so det(J) > 0 for a valid element;
!>   - then one node per edge, in Gmsh's edge order: the base edges going around, 1-2, 2-3, 3-1
!>     (counterclockwise), then the edges from the apex 4-1, 4-3, 4-2 (4-3 before 4-2: the
!>     middle layer runs clockwise seen from above).
!> v2 numbered the same nodes layer by layer (base 1-6, middle 7-9, apex 10); Gmsh nodes 1-10
!> are v2 nodes 1 3 5 10 2 4 6 7 9 8.
!>
!> Face element (the boundary faces of a tet10 mesh, read from $Fasets; ROADMAP D36, D37):
!> 6-node triangle, Gmsh order (Gmsh type 9; Gmsh reference manual §10.2.1, Triangle6), reference
!> coordinates (xi1, xi2), corners (0,0), (1,0), (0,1):
!>
!>   xi2
!>    ^
!>    3
!>    | \
!>    |   \
!>    6     5
!>    |       \
!>    |         \
!>    1-----4-----2 --> xi1
!>
!>    node  xi            node  xi            on edge
!>     1   (0, 0) corner    4   (1/2, 0  )    1-2
!>     2   (1, 0) corner    5   (1/2, 1/2)    2-3
!>     3   (0, 1) corner    6   (0,   1/2)    3-1
!>
!> Corners 1, 2, 3 counterclockwise in (xi1, xi2), then the edge nodes of edges 1-2, 2-3, 3-1
!> (counterclockwise). The faces of a tet10 are such triangles: e.g. the face xi3 = 0 is element
!> nodes 1, 2, 3, 5, 6, 7 (in this order its normal points into the element). Faces read from
!> $Fasets are assumed to have the normal (x2 - x1) x (x3 - x1) pointing out of the body, as
!> Neper writes them (not checked; mesh_file_mod, ROADMAP D37).
!> Shape functions: port of v2 shape_2d_mod.f90::sf2d06 and
!> sf2dg06 (same expressions; Gmsh node a is v2 node 5 1 3 6 2 4), with L = (1 - xi1 - xi2, xi1,
!> xi2): corner N = L (2L - 1), edge node N = 4 L_i L_j.
!> Face quadrature: the 7-point rule of degree 5 (Radon 1948) in its closed form, Felippa,
!> Introduction to Finite Element Methods, Ch. 24, Fig. 24.3 (refs/felippa_ifem_ch24.pdf):
!> points (g, g, 1 - 2g) in triangular coordinates for g = (6 -+ sqrt 15)/21, and the centroid;
!> weights (155 -+ sqrt 15)/1200 and 9/40 for unit area, halved here for the reference triangle
!> (area 1/2, so the weights sum to 1/2). v2's point order (quadrature_mod.f90::quadrature_init,
!> qp2d, wt2d; v2-fixed CHANGES.md C23 gave v2 the same full-precision values).
!>
!> Public parameters:
!>   tet10_gmsh_type     integer(ik)            Gmsh element type code (11)
!>   tet10_num_nodes     integer(ik)            number of nodes (10)
!>   tet10_node_coo_ref  real(rk)(3, 10)        reference coordinates of the nodes [-]
!>   tet10_num_qpts      integer(ik)            number of quadrature points (15)
!>   tet10_qpt_centroid  integer(ik)            the point at the centroid (orbit 2, point 5):
!>                                              elemental output (v2 manual, output.rst;
!>                                              dev/cfg_keywords.md §6)
!>   tet10_qpt_coo       real(rk)(3, 15)        reference coordinates of the quadrature points [-]
!>   tet10_qpt_weight    real(rk)(15)           quadrature weights (sum 1/6) [-]
!>
!>   tet10_face_num_nodes     integer(ik)       number of face element nodes (6)
!>   tet10_face_node_coo_ref  real(rk)(2, 6)    reference coordinates of the face nodes [-]
!>   tet10_face_num_qpts      integer(ik)       number of face quadrature points (7)
!>   tet10_face_qpt_coo       real(rk)(2, 7)    face reference coordinates of the points [-]
!>   tet10_face_qpt_weight    real(rk)(7)       face quadrature weights (sum 1/2) [-]
!>
!> Public routines:
!>   tet10_shape        shape functions at one point
!>   tet10_shape_deriv  their derivatives with respect to the reference coordinates at one point
!>   tet10_face_shape        face element (6-node triangle) shape functions at one point
!>   tet10_face_shape_deriv  their derivatives with respect to the face reference coordinates
module tet10_mod

  use kinds_mod, only: rk, ik

  implicit none

  private

  public :: tet10_gmsh_type, tet10_num_nodes, tet10_node_coo_ref
  public :: tet10_num_qpts, tet10_qpt_coo, tet10_qpt_weight, tet10_qpt_centroid
  public :: tet10_shape, tet10_shape_deriv
  public :: tet10_face_num_nodes, tet10_face_node_coo_ref
  public :: tet10_face_num_qpts, tet10_face_qpt_coo, tet10_face_qpt_weight
  public :: tet10_face_shape, tet10_face_shape_deriv

  ! Gmsh element type code.
  integer(ik), parameter :: tet10_gmsh_type = 11_ik

  ! Number of nodes.
  integer(ik), parameter :: tet10_num_nodes = 10_ik

  ! Reference coordinates of the nodes (one column per node, Gmsh order).
  real(rk), parameter :: tet10_node_coo_ref(3, tet10_num_nodes) = reshape([ &
    & 0.0_rk, 0.0_rk, 0.0_rk, &
    & 1.0_rk, 0.0_rk, 0.0_rk, &
    & 0.0_rk, 1.0_rk, 0.0_rk, &
    & 0.0_rk, 0.0_rk, 1.0_rk, &
    & 0.5_rk, 0.0_rk, 0.0_rk, &
    & 0.5_rk, 0.5_rk, 0.0_rk, &
    & 0.0_rk, 0.5_rk, 0.0_rk, &
    & 0.0_rk, 0.0_rk, 0.5_rk, &
    & 0.0_rk, 0.5_rk, 0.5_rk, &
    & 0.5_rk, 0.0_rk, 0.5_rk], [3, int(tet10_num_nodes)])

  ! Quadrature rule (Keast 1986, 15 points, degree 5; v2's values and point order).

  ! Number of quadrature points.
  integer(ik), parameter :: tet10_num_qpts = 15_ik

  ! The quadrature point at the centroid (orbit 2) (elemental output).
  integer(ik), parameter :: tet10_qpt_centroid = 5_ik

  ! Coordinate values of the rule, typed as in v2 (named for the orbit they belong to):
  !   tet_third   1/3                       (orbit of (1/3, 1/3, 1/3, 0))
  !   tet_quarter 1/4                       (centroid)
  !   tet_a3      1/11, tet_b3  8/11        (orbit of (1/11, 1/11, 1/11, 8/11))
  !   tet_a4      a,    tet_b4  b = 1/2 - a (orbit of (a, a, b, b))
  real(rk), parameter :: tet_third = 0.333333333333333333_rk
  real(rk), parameter :: tet_quarter = 0.25_rk
  real(rk), parameter :: tet_a3 = 0.909090909090909091e-1_rk
  real(rk), parameter :: tet_b3 = 0.727272727272727273_rk
  real(rk), parameter :: tet_a4 = 0.665501535736642813e-1_rk
  real(rk), parameter :: tet_b4 = 0.433449846426335728_rk

  ! Weights of the four orbits, typed as in v2.
  real(rk), parameter :: tet_w1 = 0.602678571428571597e-2_rk
  real(rk), parameter :: tet_w2 = 0.302836780970891856e-1_rk
  real(rk), parameter :: tet_w3 = 0.116452490860289742e-1_rk
  real(rk), parameter :: tet_w4 = 0.109491415613864534e-1_rk

  ! Quadrature points (one column per point, v2's order).
  real(rk), parameter :: tet10_qpt_coo(3, tet10_num_qpts) = reshape([ &
    ! orbit 1: (1/3, 1/3, 1/3, 0), points 1-4
    & tet_third, tet_third, tet_third, &
    & tet_third, tet_third, 0.0_rk, &
    & tet_third, 0.0_rk, tet_third, &
    & 0.0_rk, tet_third, tet_third, &
    ! orbit 2: the centroid, point 5
    & tet_quarter, tet_quarter, tet_quarter, &
    ! orbit 3: (1/11, 1/11, 1/11, 8/11), points 6-9
    & tet_a3, tet_a3, tet_a3, &
    & tet_a3, tet_a3, tet_b3, &
    & tet_a3, tet_b3, tet_a3, &
    & tet_b3, tet_a3, tet_a3, &
    ! orbit 4: (a, a, b, b), points 10-15
    & tet_a4, tet_a4, tet_b4, &
    & tet_a4, tet_b4, tet_a4, &
    & tet_a4, tet_b4, tet_b4, &
    & tet_b4, tet_a4, tet_a4, &
    & tet_b4, tet_a4, tet_b4, &
    & tet_b4, tet_b4, tet_a4], [3, int(tet10_num_qpts)])

  ! Quadrature weights (v2's order).
  real(rk), parameter :: tet10_qpt_weight(tet10_num_qpts) = [ &
    & tet_w1, tet_w1, tet_w1, tet_w1, &
    & tet_w2, &
    & tet_w3, tet_w3, tet_w3, tet_w3, &
    & tet_w4, tet_w4, tet_w4, tet_w4, tet_w4, tet_w4]

  ! Face element (6-node triangle): number of nodes and reference coordinates (Gmsh order).
  integer(ik), parameter :: tet10_face_num_nodes = 6_ik
  real(rk), parameter :: tet10_face_node_coo_ref(2, tet10_face_num_nodes) = reshape([ &
    & 0.0_rk, 0.0_rk, &
    & 1.0_rk, 0.0_rk, &
    & 0.0_rk, 1.0_rk, &
    & 0.5_rk, 0.0_rk, &
    & 0.5_rk, 0.5_rk, &
    & 0.0_rk, 0.5_rk], [2, int(tet10_face_num_nodes)])

  ! Face quadrature rule (7 points, degree 5; Felippa Ch. 24, Fig. 24.3; v2's point order).
  ! Coordinate values (34 digits; the compiler rounds them to the nearest real):
  !   tri_third   1/3                       (centroid)
  !   tri_a1      1 - 2 g2, tri_b1  g2      g2 = (6 + sqrt 15)/21 (points 2-4)
  !   tri_a2      1 - 2 g1, tri_b2  g1      g1 = (6 - sqrt 15)/21 (points 5-7)
  integer(ik), parameter :: tet10_face_num_qpts = 7_ik
  real(rk), parameter :: tri_third = 0.33333333333333333333333333333333333_rk
  real(rk), parameter :: tri_a1 = 0.059715871789769820459117580973104799_rk
  real(rk), parameter :: tri_b1 = 0.47014206410511508977044120951344760_rk
  real(rk), parameter :: tri_a2 = 0.79742698535308732239802527616975234_rk
  real(rk), parameter :: tri_b2 = 0.10128650732345633880098736191512383_rk

  ! Weights: 9/80 at the centroid, (155 + sqrt 15)/2400 at points 2-4, (155 - sqrt 15)/2400 at
  ! points 5-7.
  real(rk), parameter :: tri_w0 = 9.0_rk/80.0_rk
  real(rk), parameter :: tri_w1 = 0.066197076394253090368824693916576000_rk
  real(rk), parameter :: tri_w2 = 0.062969590272413576297841972750090667_rk

  ! Face quadrature points (one column per point).
  real(rk), parameter :: tet10_face_qpt_coo(2, tet10_face_num_qpts) = reshape([ &
    & tri_third, tri_third, &
    & tri_a1, tri_b1, &
    & tri_b1, tri_a1, &
    & tri_b1, tri_b1, &
    & tri_a2, tri_b2, &
    & tri_b2, tri_a2, &
    & tri_b2, tri_b2], [2, int(tet10_face_num_qpts)])

  ! Face quadrature weights.
  real(rk), parameter :: tet10_face_qpt_weight(tet10_face_num_qpts) = &
    & [tri_w0, tri_w1, tri_w1, tri_w1, tri_w2, tri_w2, tri_w2]

contains

  !> Shape functions of the 10-node tetrahedron at one point.
  !>
  !> Theory : module header (N = L_i (2 L_i - 1) at corners, 4 L_i L_j at edge nodes).
  !> Source : port of v2 shape_3d_mod.f90::t10_shape_hpar (same expressions, Gmsh node
  !>          numbering; v2's loc0, loc1, loc2 = xi(1), xi(2), xi(3)).
  !>
  !> Arguments:
  !>   xi     (in)  real(rk)(3)   reference coordinates of the point [-]
  !>   shape  (out) real(rk)(10)  shape function of each element node (Gmsh order) at the point [-]
  !> Locals:
  !>   sum_xi  real(rk)  xi1 + xi2 + xi3 (= 1 - L_1) [-]
  pure function tet10_shape(xi) result(shape)

    real(rk), intent(in) :: xi(3)
    real(rk) :: shape(tet10_num_nodes)

    real(rk) :: sum_xi

    ! Sum of the reference coordinates (L_1 = 1 - sum_xi).
    sum_xi = xi(1) + xi(2) + xi(3)

    ! Corner nodes 1-4: L (2L - 1) (v2 nodes 1, 3, 5, 10).
    shape(1) = 2.0_rk*(sum_xi - 1.0_rk)*(sum_xi - 0.5_rk)
    shape(2) = 2.0_rk*xi(1)*(xi(1) - 0.5_rk)
    shape(3) = 2.0_rk*xi(2)*(xi(2) - 0.5_rk)
    shape(4) = 2.0_rk*xi(3)*(xi(3) - 0.5_rk)

    ! Edge nodes 5-10: 4 L_i L_j (v2 nodes 2, 4, 6, 7, 9, 8).
    shape(5) = -4.0_rk*(sum_xi - 1.0_rk)*xi(1)
    shape(6) = 4.0_rk*xi(2)*xi(1)
    shape(7) = -4.0_rk*(sum_xi - 1.0_rk)*xi(2)
    shape(8) = -4.0_rk*(sum_xi - 1.0_rk)*xi(3)
    shape(9) = 4.0_rk*xi(2)*xi(3)
    shape(10) = 4.0_rk*xi(1)*xi(3)

  end function tet10_shape

  !> Derivatives of the shape functions of the 10-node tetrahedron with respect to the
  !> reference coordinates, at one point.
  !>
  !> Theory : derivatives of the functions of tet10_shape.
  !> Source : port of v2 shape_3d_mod.f90::t10_deriv_hpar (same expressions, Gmsh node
  !>          numbering; v2's dnda, dndb, dndc = dshape(1, :), dshape(2, :), dshape(3, :)).
  !>
  !> Arguments:
  !>   xi      (in)  real(rk)(3)      reference coordinates of the point [-]
  !>   dshape  (out) real(rk)(3, 10)  dshape(i, a) = d N_a / d xi_i at the point [-]
  !> Locals:
  !>   sum_xi  real(rk)  xi1 + xi2 + xi3 [-]
  pure function tet10_shape_deriv(xi) result(dshape)

    real(rk), intent(in) :: xi(3)
    real(rk) :: dshape(3, tet10_num_nodes)

    real(rk) :: sum_xi

    ! Sum of the reference coordinates.
    sum_xi = xi(1) + xi(2) + xi(3)

    ! Derivatives with respect to xi1 (v2 dnda, renumbered).
    dshape(1, 1) = 4.0_rk*sum_xi - 3.0_rk
    dshape(1, 2) = 4.0_rk*xi(1) - 1.0_rk
    dshape(1, 3) = 0.0_rk
    dshape(1, 4) = 0.0_rk
    dshape(1, 5) = -4.0_rk*(2.0_rk*xi(1) + xi(2) + xi(3) - 1.0_rk)
    dshape(1, 6) = 4.0_rk*xi(2)
    dshape(1, 7) = -4.0_rk*xi(2)
    dshape(1, 8) = -4.0_rk*xi(3)
    dshape(1, 9) = 0.0_rk
    dshape(1, 10) = 4.0_rk*xi(3)

    ! Derivatives with respect to xi2 (v2 dndb, renumbered).
    dshape(2, 1) = 4.0_rk*sum_xi - 3.0_rk
    dshape(2, 2) = 0.0_rk
    dshape(2, 3) = 4.0_rk*xi(2) - 1.0_rk
    dshape(2, 4) = 0.0_rk
    dshape(2, 5) = -4.0_rk*xi(1)
    dshape(2, 6) = 4.0_rk*xi(1)
    dshape(2, 7) = -4.0_rk*(xi(1) + 2.0_rk*xi(2) + xi(3) - 1.0_rk)
    dshape(2, 8) = -4.0_rk*xi(3)
    dshape(2, 9) = 4.0_rk*xi(3)
    dshape(2, 10) = 0.0_rk

    ! Derivatives with respect to xi3 (v2 dndc, renumbered).
    dshape(3, 1) = 4.0_rk*sum_xi - 3.0_rk
    dshape(3, 2) = 0.0_rk
    dshape(3, 3) = 0.0_rk
    dshape(3, 4) = 4.0_rk*xi(3) - 1.0_rk
    dshape(3, 5) = -4.0_rk*xi(1)
    dshape(3, 6) = 0.0_rk
    dshape(3, 7) = -4.0_rk*xi(2)
    dshape(3, 8) = -4.0_rk*(xi(1) + xi(2) + 2.0_rk*xi(3) - 1.0_rk)
    dshape(3, 9) = 4.0_rk*xi(2)
    dshape(3, 10) = 4.0_rk*xi(1)

  end function tet10_shape_deriv

  !> Shape functions of the face element (6-node triangle) at one point.
  !>
  !> Theory : module header (face element: N = L (2L - 1) at corners, 4 L_i L_j at edge nodes).
  !> Source : port of v2 shape_2d_mod.f90::sf2d06 (same expressions, Gmsh node numbering; v2's
  !>          xi, eta, zeta = xi(1), xi(2), 1 - xi(1) - xi(2)).
  !>
  !> Arguments:
  !>   xi     (in)  real(rk)(2)  face reference coordinates of the point [-]
  !>   shape  (out) real(rk)(6)  shape function of each face node (Gmsh order) at the point [-]
  !> Locals:
  !>   zeta   real(rk)  1 - xi1 - xi2 (the volume coordinate of corner 1) [-]
  pure function tet10_face_shape(xi) result(shape)

    real(rk), intent(in) :: xi(2)
    real(rk) :: shape(tet10_face_num_nodes)

    real(rk) :: zeta

    ! Third volume coordinate.
    zeta = 1.0_rk - xi(1) - xi(2)

    ! Corners 1-3 (v2 nodes 5, 1, 3).
    shape(1) = (2.0_rk*zeta - 1.0_rk)*zeta
    shape(2) = (2.0_rk*xi(1) - 1.0_rk)*xi(1)
    shape(3) = (2.0_rk*xi(2) - 1.0_rk)*xi(2)

    ! Edge nodes 4-6 (v2 nodes 6, 2, 4).
    shape(4) = 4.0_rk*xi(1)*zeta
    shape(5) = 4.0_rk*xi(2)*xi(1)
    shape(6) = 4.0_rk*xi(2)*zeta

  end function tet10_face_shape

  !> Derivatives of the face element (6-node triangle) shape functions with respect to the face
  !> reference coordinates, at one point.
  !>
  !> Theory : derivatives of the functions of tet10_face_shape.
  !> Source : port of v2 shape_2d_mod.f90::sf2dg06 (same expressions, Gmsh node numbering).
  !>
  !> Arguments:
  !>   xi      (in)  real(rk)(2)     face reference coordinates of the point [-]
  !>   dshape  (out) real(rk)(2, 6)  dshape(i, a) = d N_a / d xi_i at the point [-]
  !> Locals:
  !>   zeta    real(rk)  1 - xi1 - xi2 [-]
  pure function tet10_face_shape_deriv(xi) result(dshape)

    real(rk), intent(in) :: xi(2)
    real(rk) :: dshape(2, tet10_face_num_nodes)

    real(rk) :: zeta

    ! Third volume coordinate.
    zeta = 1.0_rk - xi(1) - xi(2)

    ! Derivatives with respect to xi1 (v2 grad(1, :), renumbered).
    dshape(1, 1) = -4.0_rk*zeta + 1.0_rk
    dshape(1, 2) = 4.0_rk*xi(1) - 1.0_rk
    dshape(1, 3) = 0.0_rk
    dshape(1, 4) = 4.0_rk*zeta - 4.0_rk*xi(1)
    dshape(1, 5) = 4.0_rk*xi(2)
    dshape(1, 6) = -4.0_rk*xi(2)

    ! Derivatives with respect to xi2 (v2 grad(2, :), renumbered).
    dshape(2, 1) = -4.0_rk*zeta + 1.0_rk
    dshape(2, 2) = 0.0_rk
    dshape(2, 3) = 4.0_rk*xi(2) - 1.0_rk
    dshape(2, 4) = -4.0_rk*xi(1)
    dshape(2, 5) = 4.0_rk*xi(1)
    dshape(2, 6) = 4.0_rk*zeta - 4.0_rk*xi(2)

  end function tet10_face_shape_deriv

end module tet10_mod
