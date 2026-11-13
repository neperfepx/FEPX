! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Writes FEPX's 15-point tetrahedron quadrature rule, its 10-node tetrahedron shape functions
!> and derivatives, and the 7-point rule and 6-node triangle of its faces (src/element/tet10_mod.f90)
!> to standard output, as the fixture tet10_reference.txt of tests/element/test_tet10.f90.
!> Not part of FEPX; built against the objects of an FEPX build (README.md in this directory).
!>
!> Source : FEPX (ROADMAP step 9.13, D91: replaces v2_dump_tet10.f90, which wrote v2-fixed's
!>          values in v2's node order).
!>
!> Output (es25.17, Gmsh node order throughout):
!>   "qpts 15", then one line per quadrature point: xi1 xi2 xi3 weight;
!>   "points 30", then one line per point: xi1 xi2 xi3, the 10 shape functions, then the 10
!>   derivatives with respect to xi1, the 10 with respect to xi2 and the 10 with respect to xi3
!>   (43 values). Points: the 15 quadrature points, the 10 element nodes and 5 other points inside
!>   the element.
!>   "tri_qpts 7", then one line per face quadrature point: xi1 xi2 weight;
!>   "tri_points 16", then one line per point: xi1 xi2, the 6 face shape functions, the 6
!>   derivatives with respect to xi1 and the 6 with respect to xi2 (20 values). Points: the 7
!>   face quadrature points, the 6 face nodes and 3 others.
!>
!> Locals:
!>   others       real(rk)(3, 5)  reference coordinates of the 5 other points [-]
!>   tri_others   real(rk)(2, 3)  reference coordinates of the 3 other face points [-]
!>   i            integer(ik)     point index
program dump_tet10

  use kinds_mod, only: rk, ik
  use tet10_mod, only: tet10_num_nodes, tet10_node_coo_ref, tet10_num_qpts, tet10_qpt_coo, &
    & tet10_qpt_weight, tet10_shape, tet10_shape_deriv, tet10_face_num_nodes, &
    & tet10_face_node_coo_ref, tet10_face_num_qpts, tet10_face_qpt_coo, tet10_face_qpt_weight, &
    & tet10_face_shape, tet10_face_shape_deriv

  implicit none

  real(rk), parameter :: others(3, 5) = reshape([ &
    & 0.1_rk, 0.2_rk, 0.3_rk, 0.25_rk, 0.25_rk, 0.25_rk, 0.6_rk, 0.1_rk, 0.05_rk, &
    & 0.05_rk, 0.7_rk, 0.2_rk, 0.15_rk, 0.05_rk, 0.75_rk], [3, 5])
  real(rk), parameter :: tri_others(2, 3) = reshape([0.1_rk, 0.2_rk, 0.6_rk, 0.3_rk, 0.15_rk, &
    & 0.7_rk], [2, 3])
  integer(ik) :: i

  ! Quadrature rule.
  write (*, '(a,1x,i0)') 'qpts', tet10_num_qpts
  do i = 1, tet10_num_qpts
    write (*, '(4es25.17)') tet10_qpt_coo(:, i), tet10_qpt_weight(i)
  end do

  ! Shape functions and derivatives at the quadrature points, the nodes and the other points.
  write (*, '(a,1x,i0)') 'points', tet10_num_qpts + tet10_num_nodes + 5
  do i = 1, tet10_num_qpts
    call dump(tet10_qpt_coo(:, i))
  end do
  do i = 1, tet10_num_nodes
    call dump(tet10_node_coo_ref(:, i))
  end do
  do i = 1, 5
    call dump(others(:, i))
  end do

  ! Face rule and 6-node triangle.
  write (*, '(a,1x,i0)') 'tri_qpts', tet10_face_num_qpts
  do i = 1, tet10_face_num_qpts
    write (*, '(3es25.17)') tet10_face_qpt_coo(:, i), tet10_face_qpt_weight(i)
  end do
  write (*, '(a,1x,i0)') 'tri_points', tet10_face_num_qpts + tet10_face_num_nodes + 3
  do i = 1, tet10_face_num_qpts
    call dump_tri(tet10_face_qpt_coo(:, i))
  end do
  do i = 1, tet10_face_num_nodes
    call dump_tri(tet10_face_node_coo_ref(:, i))
  end do
  do i = 1, 3
    call dump_tri(tri_others(:, i))
  end do

contains

  !> Writes one line of the tetrahedron at a point.
  !>
  !> Arguments:
  !>   xi      (in)  real(rk)(3)  reference coordinates [-]
  !> Locals:
  !>   dshape  real(rk)(3, 10)    shape-function derivatives with respect to xi [-]
  subroutine dump(xi)
    real(rk), intent(in) :: xi(3)
    real(rk) :: dshape(3, tet10_num_nodes)
    dshape = tet10_shape_deriv(xi)
    write (*, '(43es25.17)') xi, tet10_shape(xi), dshape(1, :), dshape(2, :), dshape(3, :)
  end subroutine dump

  !> Writes one line of the face triangle at a point.
  !>
  !> Arguments:
  !>   xi      (in)  real(rk)(2)  face reference coordinates [-]
  !> Locals:
  !>   dshape  real(rk)(2, 6)     face shape-function derivatives with respect to xi [-]
  subroutine dump_tri(xi)
    real(rk), intent(in) :: xi(2)
    real(rk) :: dshape(2, tet10_face_num_nodes)
    dshape = tet10_face_shape_deriv(xi)
    write (*, '(20es25.17)') xi, tet10_face_shape(xi), dshape(1, :), dshape(2, :)
  end subroutine dump_tri

end program dump_tet10
