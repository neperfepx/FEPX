! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Writes FEPX's 8-node quadrilateral face shape functions and derivatives of the 20-node
!> hexahedron (src/element/hex20_mod.f90::hex20_face_shape, hex20_face_shape_deriv) to standard
!> output, as the fixture quad8_reference.txt of tests/element/test_hex20.f90.
!> Not part of FEPX; built against the objects of an FEPX build (README.md in this directory).
!>
!> Source : FEPX (ROADMAP step 9.13, D91: replaces v2_dump_quad8.f90, which wrote v2-fixed's
!>          values on v2's reference square [0, 1]^2 in v2's node order).
!>
!> Output: "quad_points 12", then one line per point: xi1 xi2 (reference coordinates in
!> [-1, 1]^2), the 8 shape functions, the 8 derivatives with respect to xi1 and the 8 with respect
!> to xi2 (26 values, es25.17, Gmsh face node order). Points: the 8 face nodes and 4 others.
!>
!> Locals:
!>   others  real(rk)(2, 4)  reference coordinates of the 4 other points [-]
!>   dshape  real(rk)(2, 8)  shape-function derivatives with respect to xi [-]
!>   i       integer(ik)     point index
program dump_quad8

  use kinds_mod, only: rk, ik
  use hex20_mod, only: hex20_face_num_nodes, hex20_face_node_coo_ref, hex20_face_shape, &
    & hex20_face_shape_deriv

  implicit none

  real(rk), parameter :: others(2, 4) = reshape([-0.8_rk, -0.6_rk, 0.4_rk, -0.2_rk, -0.3_rk, &
    & 0.8_rk, 0.0_rk, 0.0_rk], [2, 4])
  real(rk) :: dshape(2, hex20_face_num_nodes)
  integer(ik) :: i

  ! The face nodes, then the other points.
  write (*, '(a,1x,i0)') 'quad_points', hex20_face_num_nodes + 4
  do i = 1, hex20_face_num_nodes
    dshape = hex20_face_shape_deriv(hex20_face_node_coo_ref(:, i))
    write (*, '(26es25.17)') hex20_face_node_coo_ref(:, i), &
      & hex20_face_shape(hex20_face_node_coo_ref(:, i)), dshape(1, :), dshape(2, :)
  end do
  do i = 1, 4
    dshape = hex20_face_shape_deriv(others(:, i))
    write (*, '(26es25.17)') others(:, i), hex20_face_shape(others(:, i)), dshape(1, :), &
      & dshape(2, :)
  end do

end program dump_quad8
