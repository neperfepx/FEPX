! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Writes FEPX's element geometry (src/element/element_mod.f90::element_geometry_qpt: Jacobian
!> determinant and the shape-function derivatives with respect to x, y, z) at the 15 quadrature
!> points of three 10-node tetrahedra, as the fixture geometry_reference.txt of
!> tests/element/test_element.f90.
!> Not part of FEPX; built against the objects of an FEPX build (README.md in this directory).
!>
!> Source : FEPX (ROADMAP step 9.13, D91: replaces v2_dump_geometry.f90, which wrote v2-fixed's
!>          values for the same elements).
!>
!> Elements (node coordinates in Gmsh order, from the reference coordinates xi_a):
!>   1. x = 2 xi (scaled reference element);
!>   2. x = A xi + b (affine, A below);
!>   3. element 2 with its edge nodes 5-10 moved by 0.05 (sin 3a, cos 2a, sin(a + 1)) (curved).
!>
!> Output: "elements 3", then for each element 10 lines "x y z" (Gmsh order); "qpts 15", then for
!> each element and quadrature point one line: det J, then d N_a/dx (10), d N_a/dy (10),
!> d N_a/dz (10) in Gmsh order (31 values, es25.17).
!>
!> Locals:
!>   amat      real(rk)(3, 3)      matrix A of the affine map [mm]
!>   bvec      real(rk)(3)         offset b of the affine map [mm]
!>   tet       type(element_t)     10-node tetrahedron descriptor
!>   coo       real(rk)(3, 10, 3)  node coordinates of the three elements [mm]
!>   dshape_s  real(rk)(3, 10)     shape-function derivatives with respect to x [1/mm]
!>   det_jac   real(rk)            Jacobian determinant [mm^3]
!>   a, e, q   integer(ik)         node, element and quadrature point indices
program dump_geometry

  use kinds_mod, only: rk, ik
  use element_mod, only: element_t, element_init, element_geometry_qpt
  use tet10_mod, only: tet10_gmsh_type, tet10_num_nodes, tet10_node_coo_ref, tet10_num_qpts

  implicit none

  integer(ik), parameter :: num_elts = 3
  real(rk), parameter :: amat(3, 3) = reshape([1.0_rk, 0.2_rk, 0.1_rk, 0.3_rk, 1.5_rk, &
    & 0.2_rk, 0.1_rk, 0.4_rk, 0.8_rk], [3, 3])
  real(rk), parameter :: bvec(3) = [0.5_rk, -0.2_rk, 1.0_rk]
  type(element_t) :: tet
  real(rk) :: coo(3, tet10_num_nodes, num_elts), dshape_s(3, tet10_num_nodes), det_jac
  integer(ik) :: a, e, q

  call element_init(tet, tet10_gmsh_type)

  ! Element coordinates, Gmsh order.
  do a = 1, tet10_num_nodes
    coo(:, a, 1) = 2.0_rk*tet10_node_coo_ref(:, a)
    coo(:, a, 2) = matmul(amat, tet10_node_coo_ref(:, a)) + bvec
    coo(:, a, 3) = coo(:, a, 2)
    if (a >= 5) coo(:, a, 3) = coo(:, a, 3) + 0.05_rk*[sin(3.0_rk*a), cos(2.0_rk*a), &
      & sin(a + 1.0_rk)]
  end do
  write (*, '(a,1x,i0)') 'elements', num_elts
  do e = 1, num_elts
    do a = 1, tet10_num_nodes
      write (*, '(3es25.17)') coo(:, a, e)
    end do
  end do

  ! det J and the derivatives at each quadrature point.
  write (*, '(a,1x,i0)') 'qpts', tet10_num_qpts
  do e = 1, num_elts
    do q = 1, tet10_num_qpts
      call element_geometry_qpt(tet, coo(:, :, e), q, dshape_s, det_jac)
      write (*, '(31es25.17)') det_jac, dshape_s(1, :), dshape_s(2, :), dshape_s(3, :)
    end do
  end do

end program dump_geometry
