! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module shape_3d_mod

! Shape function library for 3d elements. Contains full descriptions of 4-node
!   (linear) and 10-node (quadratic) tetrahedral elements, as well as 8-node
!   cubic elements. Currently only 10-node tetrahedral elements are supported
!   elsewhere in the code.

! Contains subroutines:
! sfder_hpar: Compute quadrature quantities given a set of local coordinates.
!   For 10-node tetrahedral elements specfically:
! sf10t_eval_vec:Evaluate shape fcns for vector quantity at array of pts.
! t10_shape_hpar: Shape functions.
! t10_deriv_hpar: Shape function derivatives.
!   For 4-node tetrahedral elements specifically:
! sf4t_eval_vec: Evaluate shape fcns. for vector quantity at array of pts.
! t4_shape_hpar: Shape functions.
! t4_deriv_hpar: Shape function derivatives.
!   For 8-node cubic elements specifically:
! sf8b_eval_vec:Evaluate shape fcns for vector quantity at array of pts.
! b8_shape_hpar: Shape functions.
! b8_deriv_hpar: Shape function derivatives.

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod

  implicit none

contains

  subroutine sfder_hpar(loc0, loc1, loc2, coo, dndx, dndy, dndz, det, &
      & ijac11, ijac12, ijac13, ijac21, ijac22, ijac23, ijac31, ijac32, &
      & ijac33)

    ! Compute quadrature quantities given a set of local coordinates.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! loc0,loc1,loc2 : Coordinates of point (in reference element)
    ! coo : Real coordinates of element's nodes
    ! det : Determinant of jacobian
    ! dndx,dndy,dndz : Derivatives of mapped shape functions
    ! IJACxx : Components of inverse jacobian

    real(rk), intent(in) :: loc0
    real(rk), intent(in) :: loc1
    real(rk), intent(in) :: loc2
    real(rk), intent(in) :: coo(kdim, elt_sub:elt_sup)
    real(rk), intent(out) :: dndx(ndim, elt_sub:elt_sup)
    real(rk), intent(out) :: dndy(ndim, elt_sub:elt_sup)
    real(rk), intent(out) :: dndz(ndim, elt_sub:elt_sup)
    real(rk), intent(out) :: det(elt_sub:elt_sup)
    real(rk), intent(out) :: ijac11(elt_sub:elt_sup)
    real(rk), intent(out) :: ijac12(elt_sub:elt_sup)
    real(rk), intent(out) :: ijac13(elt_sub:elt_sup)
    real(rk), intent(out) :: ijac21(elt_sub:elt_sup)
    real(rk), intent(out) :: ijac22(elt_sub:elt_sup)
    real(rk), intent(out) :: ijac23(elt_sub:elt_sup)
    real(rk), intent(out) :: ijac31(elt_sub:elt_sup)
    real(rk), intent(out) :: ijac32(elt_sub:elt_sup)
    real(rk), intent(out) :: ijac33(elt_sub:elt_sup)

    ! Locals:

    integer :: i
    integer :: i0
    integer :: i1
    integer :: i2

    real(rk) :: dnd1
    real(rk) :: dnd2
    real(rk) :: dnd3

    real(rk) :: dnda(ndim)
    real(rk) :: dndb(ndim)
    real(rk) :: dndc(ndim)

    real(rk) :: jac11(elt_sub:elt_sup)
    real(rk) :: jac12(elt_sub:elt_sup)
    real(rk) :: jac13(elt_sub:elt_sup)
    real(rk) :: jac21(elt_sub:elt_sup)
    real(rk) :: jac22(elt_sub:elt_sup)
    real(rk) :: jac23(elt_sub:elt_sup)
    real(rk) :: jac31(elt_sub:elt_sup)
    real(rk) :: jac32(elt_sub:elt_sup)
    real(rk) :: jac33(elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    ! Evaluate shape function derivatives.

    ! call b8_deriv_hpar(loc0, loc1, loc2, dnda, dndb, dndc)
    ! call t4_deriv_hpar(dnda, dndb, dndc)
    call t10_deriv_hpar(loc0, loc1, loc2, dnda, dndb, dndc)

    ! Initialize Jacobian matrix

    jac11 = 0.0d0
    jac12 = 0.0d0
    jac13 = 0.0d0
    jac21 = 0.0d0
    jac22 = 0.0d0
    jac23 = 0.0d0
    jac31 = 0.0d0
    jac32 = 0.0d0
    jac33 = 0.0d0

    do i = 1, ndim
      i0 = 3*(i - 1) + 1
      i1 = i0 + 1
      i2 = i1 + 1

      dnd1 = dnda(i)
      dnd2 = dndb(i)
      dnd3 = dndc(i)

      ! do j = elt_sub,elt_sup (legacy?)
      ijac11(:) = coo(i0, :)
      ijac22(:) = coo(i1, :)
      ijac33(:) = coo(i2, :)
      ! end do

      jac11 = jac11 + ijac11*dnd1
      jac21 = jac21 + ijac11*dnd2
      jac31 = jac31 + ijac11*dnd3
      jac12 = jac12 + ijac22*dnd1
      jac22 = jac22 + ijac22*dnd2
      jac32 = jac32 + ijac22*dnd3
      jac13 = jac13 + ijac33*dnd1
      jac23 = jac23 + ijac33*dnd2
      jac33 = jac33 + ijac33*dnd3
    end do

    ! Determinant of the Jacobian matrix

    ijac11 = jac11*jac22*jac33
    ijac12 = jac12*jac23*jac31
    ijac13 = jac13*jac21*jac32
    ijac21 = jac11*jac23*jac32
    ijac22 = jac12*jac21*jac33
    ijac23 = jac13*jac22*jac31

    det = ijac11 + ijac12 + ijac13
    det = det - (ijac21 + ijac22 + ijac23)

    ! Inverse of the Jacobian matrix

    ijac11 = jac22*jac33 - jac23*jac32
    ijac21 = jac23*jac31 - jac21*jac33
    ijac31 = jac21*jac32 - jac22*jac31
    ijac12 = jac13*jac32 - jac12*jac33
    ijac22 = jac11*jac33 - jac13*jac31
    ijac32 = jac31*jac12 - jac11*jac32
    ijac13 = jac12*jac23 - jac13*jac22
    ijac23 = jac13*jac21 - jac11*jac23
    ijac33 = jac11*jac22 - jac12*jac21

    ijac11 = ijac11/det
    ijac12 = ijac12/det
    ijac13 = ijac13/det
    ijac21 = ijac21/det
    ijac22 = ijac22/det
    ijac23 = ijac23/det
    ijac31 = ijac31/det
    ijac32 = ijac32/det
    ijac33 = ijac33/det

    ! Shape function derivatives

    do i = 1, ndim
      dnd1 = dnda(i)
      dnd2 = dndb(i)
      dnd3 = dndc(i)

      jac11 = ijac11*dnd1 + ijac12*dnd2 + ijac13*dnd3
      jac22 = ijac21*dnd1 + ijac22*dnd2 + ijac23*dnd3
      jac33 = ijac31*dnd1 + ijac32*dnd2 + ijac33*dnd3

      dndx(i, :) = jac11
      dndy(i, :) = jac22
      dndz(i, :) = jac33
    end do

    return

  end subroutine sfder_hpar

  !===========================================================================

  subroutine sf10t_eval_vec(vec, npt, pnt, val)

    ! Evaluate sf functions for vector valued quantity at array of points.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! vec: Values of vector at nodes.
    ! pnt: Array of points in reference element.
    ! npt: Number of points
    ! val: Values of vector at points
    ! Note: we need explicit interface to use assumed shape arrays as arguments,
    !   i.e. pnt(:), val(:); unless this is put in a module or an interface, we
    !   can still use the assumed-size arrays.

    real(rk), intent(in) :: vec(3, 10)
    integer, intent(in) :: npt
    real(rk), intent(in) :: pnt(*)
    real(rk), intent(out) :: val(*)

    ! Locals:

    integer :: i
    integer :: xdof
    integer :: ydof
    integer :: zdof
    real(rk) :: x1
    real(rk) :: x2
    real(rk) :: x3
    real(rk) :: sf(10)
    real(rk) :: tmpvec(3)

    !---------------------------------------------------------------------------

    xdof = -2
    do i = 1, npt
      xdof = xdof + 3
      ydof = xdof + 1
      zdof = ydof + 1

      x1 = pnt(xdof)
      x2 = pnt(ydof)
      x3 = pnt(zdof)

      sf(1) = 2.0d0*(x1 + x2 + x3 - 1.0d0)*(x1 + x2 + x3 - 0.5d0)
      sf(2) = -4.0d0*(x1 + x2 + x3 - 1.0d0)*x1
      sf(3) = 2.0d0*x1*(x1 - 0.5d0)
      sf(4) = 4.0d0*x2*x1
      sf(5) = 2.0d0*x2*(x2 - 0.5d0)
      sf(6) = -4.0d0*(x1 + x2 + x3 - 1.0d0)*x2
      sf(7) = -4.0d0*(x1 + x2 + x3 - 1.0d0)*x3
      sf(8) = 4.0d0*x1*x3
      sf(9) = 4.0d0*x2*x3
      sf(10) = 2.0d0*x3*(x3 - 0.5d0)

      tmpvec = matmul(vec, sf)

      val(xdof) = tmpvec(1)
      val(ydof) = tmpvec(2)
      val(zdof) = tmpvec(3)
    end do

    return

  end subroutine sf10t_eval_vec

  !===========================================================================

  subroutine t10_shape_hpar(loc0, loc1, loc2, shape)

    ! Shape functions for 10-node tetrahedron.

    ! Node ordering: Looking down 3-axis (here's the coordinate system):
    !   2
    !   |
    !   3--1

    ! Top:
    !   10

    ! Middle:
    !   9
    !   7 8

    ! Bottom:
    !   5
    !   6 4
    !   1 2 3

    !---------------------------------------------------------------------------

    ! Arguments:
    ! loc0,loc1,loc2: Coordinates of point
    ! shape: Shape function value at point, distributed over all elements

    real(rk), intent(in) :: loc0
    real(rk), intent(in) :: loc1
    real(rk), intent(in) :: loc2
    real(rk), intent(out) :: shape(ndim)

    !---------------------------------------------------------------------------

    shape(1) = 2.0d0*(loc0 + loc1 + loc2 - 1.0d0)*(loc0 + loc1 + loc2 - &
        & 0.5d0)
    shape(2) = -4.0d0*(loc0 + loc1 + loc2 - 1.0d0)*loc0
    shape(3) = 2.0d0*loc0*(loc0 - 0.5d0)
    shape(4) = 4.0d0*loc1*loc0
    shape(5) = 2.0d0*loc1*(loc1 - 0.5d0)
    shape(6) = -4.0d0*(loc0 + loc1 + loc2 - 1.0d0)*loc1
    shape(7) = -4.0d0*(loc0 + loc1 + loc2 - 1.0d0)*loc2
    shape(8) = 4.0d0*loc0*loc2
    shape(9) = 4.0d0*loc1*loc2
    shape(10) = 2.0d0*loc2*(loc2 - 0.5d0)

    return

  end subroutine t10_shape_hpar

  !===========================================================================

  subroutine t10_deriv_hpar(loc0, loc1, loc2, dnda, dndb, dndc)

    ! Shape function derivatives for 10-node tetrahedron.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! loc0,loc1,loc2: Coordinates of point
    ! dnda,dndb,dndc: shape function derivatives at point, distributed over all
    !   elements

    real(rk), intent(in) :: loc0
    real(rk), intent(in) :: loc1
    real(rk), intent(in) :: loc2
    real(rk), intent(out) :: dnda(ndim)
    real(rk), intent(out) :: dndb(ndim)
    real(rk), intent(out) :: dndc(ndim)

    !---------------------------------------------------------------------------

    dnda(1) = 4.0d0*(loc0 + loc1 + loc2) - 3.0d0
    dnda(2) = -4.0d0*(2.0d0*loc0 + loc1 + loc2 - 1.0d0)
    dnda(3) = 4.0d0*loc0 - 1.0d0
    dnda(4) = 4.0d0*loc1
    dnda(5) = 0.0d0
    dnda(6) = -4.0d0*loc1
    dnda(7) = -4.0d0*loc2
    dnda(8) = 4.0d0*loc2
    dnda(9) = 0.0d0
    dnda(10) = 0.0d0

    dndb(1) = 4.0d0*(loc0 + loc1 + loc2) - 3.0d0
    dndb(2) = -4.0d0*loc0
    dndb(3) = 0.0d0
    dndb(4) = 4.0d0*loc0
    dndb(5) = 4.0d0*loc1 - 1.0d0
    dndb(6) = -4.0d0*(loc0 + 2.0d0*loc1 + loc2 - 1.0d0)
    dndb(7) = -4.0d0*loc2
    dndb(8) = 0.0d0
    dndb(9) = 4.0d0*loc2
    dndb(10) = 0.0d0

    dndc(1) = 4.0d0*(loc0 + loc1 + loc2) - 3.0d0
    dndc(2) = -4.0d0*loc0
    dndc(3) = 0.0d0
    dndc(4) = 0.0d0
    dndc(5) = 0.0d0
    dndc(6) = -4.0d0*loc1
    dndc(7) = -4.0d0*(loc0 + loc1 + 2.0d0*loc2 - 1.0d0)
    dndc(8) = 4.0d0*loc0
    dndc(9) = 4.0d0*loc1
    dndc(10) = 4.0d0*loc2 - 1.0d0

    return

  end subroutine t10_deriv_hpar

  !===========================================================================

  subroutine sf4t_eval_vec(vec, npt, pnt, val)

    ! Evaluate sf functions for vector valued quantity at array of points.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! vec: Values of vector at nodes.
    ! pnt: Array of points in reference element.
    ! npt: Number of points
    ! val: Values of vector at points
    ! Note: we need explicit interface to use assumed shape arrays as arguments,
    !   i.e. pnt(:), val(:); unless this is put in a module or an interface, we
    !   can still use the assumed-size arrays.

    real(rk), intent(in) :: vec(3, 4)
    integer, intent(in) :: npt
    real(rk), intent(in) :: pnt(*)
    real(rk), intent(out):: val(*)

    ! Locals:

    integer :: i
    integer :: xdof
    integer :: ydof
    integer :: zdof
    real(rk) :: x1
    real(rk) :: x2
    real(rk) :: x3
    real(rk) :: sf(4)
    real(rk) :: tmpvec(3)

    !---------------------------------------------------------------------------

    xdof = -2
    do i = 1, npt
      xdof = xdof + 3
      ydof = xdof + 1
      zdof = ydof + 1

      x1 = pnt(xdof)
      x2 = pnt(ydof)
      x3 = pnt(zdof)

      sf(1) = x1
      sf(2) = x2
      sf(3) = x3
      sf(4) = 1 - (x1 + x2 + x3)

      tmpvec = matmul(vec, sf)

      val(xdof) = tmpvec(1)
      val(ydof) = tmpvec(2)
      val(zdof) = tmpvec(3)
    end do

    return

  end subroutine sf4t_eval_vec

  !===========================================================================

  subroutine t4_shape_hpar(loc0, loc1, loc2, shape)

    ! Shape functions for 4-node tetrahedron.

    ! Node ordering: Looking down 3-axis (here's the coordinate system):
    !   2
    !   |
    !   3--1

    ! Top:
    !   2

    ! Bottom:
    !   3
    !   4 1

    !---------------------------------------------------------------------------

    ! Arguments:
    ! loc0,loc1,loc2: Coordinates of point
    ! shape: Shape function value at point, distributed over all elements

    real(rk), intent(in) :: loc0
    real(rk), intent(in) :: loc1
    real(rk), intent(in) :: loc2
    real(rk), intent(out) :: shape(ndim)

    !---------------------------------------------------------------------------

    shape(1) = loc0
    shape(2) = loc1
    shape(3) = loc2
    shape(4) = 1.0 - (loc0 + loc1 + loc2)

    return

  end subroutine t4_shape_hpar

  !===========================================================================

  subroutine t4_deriv_hpar(dnda, dndb, dndc)

    ! Shape function derivatives for 4-node tetrahedron.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! dnda,dndb,dndc: shape function derivatives at point, distributed over all
    !   elements

    real(rk), intent(out) :: dnda(ndim)
    real(rk), intent(out) :: dndb(ndim)
    real(rk), intent(out) :: dndc(ndim)

    !---------------------------------------------------------------------------

    dnda(1) = 1.0d0
    dnda(2) = 0.0d0
    dnda(3) = 0.0d0
    dnda(4) = -1.0d0

    dndb(1) = 0.0d0
    dndb(2) = 0.0d0
    dndb(3) = 1.0d0
    dndb(4) = -1.0d0

    dndc(1) = 0.0d0
    dndc(2) = 1.0d0
    dndc(3) = 0.0d0
    dndc(4) = -1.0d0

    return

  end subroutine t4_deriv_hpar

  !===========================================================================

  subroutine sf8b_eval_vec(vec, npt, pnt, val)

    ! Evaluate sf functions for vector valued quantity at array of points.

    !---------------------------------------------------------------------------

    ! Arguments.
    ! vec: Values of vector at nodes.
    ! pnt: Array of points in reference element.
    ! npt: Number of points
    ! val: Values of vector at points
    ! Note: we need explicit interface to use assumed shape arrays as arguments,
    !   i.e. pnt(:), val(:); unless this is put in a module or an interface, we
    !   can still use the assumed-size arrays.

    real(rk), intent(in) :: vec(3, 8)
    integer, intent(in) :: npt
    real(rk), intent(in) :: pnt(*)
    real(rk), intent(out):: val(*)

    ! Locals:

    integer :: i
    integer :: xdof
    integer :: ydof
    integer :: zdof
    real(rk) :: x1
    real(rk) :: x2
    real(rk) :: x3
    real(rk) :: sf(8)
    real(rk) :: tmpvec(3)

    !---------------------------------------------------------------------------

    xdof = -2
    do i = 1, npt
      xdof = xdof + 3
      ydof = xdof + 1
      zdof = ydof + 1

      x1 = pnt(xdof)
      x2 = pnt(ydof)
      x3 = pnt(zdof)

      sf(1) = -(x3 - 1.0)*(x2 - 1.0)*(x1 - 1.0)
      sf(2) = (x3 - 1.0)*(x2 - 1.0)*x1
      sf(3) = -(x3 - 1.0)*x2*x1
      sf(4) = (x3 - 1.0)*x2*(x1 - 1.0)
      sf(5) = x3*(x2 - 1.0)*(x1 - 1.0)
      sf(6) = -x3*(x2 - 1.0)*x1
      sf(7) = x3*x2*x1
      sf(8) = -x3*x2*(x1 - 1.0)

      tmpvec = matmul(vec, sf)

      val(xdof) = tmpvec(1)
      val(ydof) = tmpvec(2)
      val(zdof) = tmpvec(3)
    end do

    return

  end subroutine sf8b_eval_vec

  !===========================================================================

  subroutine b8_shape_hpar(loc0, loc1, loc2, shape)

    !     Shape functions for 8-node brick.

    !---------------------------------------------------------------------------

    ! Arguments.
    ! loc0,loc1,loc2: Coordinates of point
    ! shape: Shape function value at point, distributed over all elements

    real(rk), intent(in) :: loc0
    real(rk), intent(in) :: loc1
    real(rk), intent(in) :: loc2
    real(rk), intent(out) :: shape(ndim)

    !---------------------------------------------------------------------------

    shape(1) = -(loc2 - 1.0)*(loc1 - 1.0)*(loc0 - 1.0)
    shape(2) = (loc2 - 1.0)*(loc1 - 1.0)*loc0
    shape(3) = -(loc2 - 1.0)*loc1*loc0
    shape(4) = (loc2 - 1.0)*loc1*(loc0 - 1.0)
    shape(5) = loc2*(loc1 - 1.0)*(loc0 - 1.0)
    shape(6) = -loc2*(loc1 - 1.0)*loc0
    shape(7) = loc2*loc1*loc0
    shape(8) = -loc2*loc1*(loc0 - 1.0)

    return

  end subroutine b8_shape_hpar

  !===========================================================================

  subroutine b8_deriv_hpar(loc0, loc1, loc2, dnda, dndb, dndc)

    ! Shape function derivatives for 8-node brick.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! loc0,loc1,loc2: Coordinates of point
    ! dnda,dndb,dndc: Shape function derivatives at point, distributed over all
    !   elements

    real(rk), intent(in) :: loc0
    real(rk), intent(in) :: loc1
    real(rk), intent(in) :: loc2
    real(rk), intent(out) :: dnda(ndim)
    real(rk), intent(out) :: dndb(ndim)
    real(rk), intent(out) :: dndc(ndim)

    !---------------------------------------------------------------------------

    dnda(1) = -(loc2 - 1.0)*(loc1 - 1.0)
    dnda(2) = (loc2 - 1.0)*(loc1 - 1.0)
    dnda(3) = -(loc2 - 1.0)*loc1
    dnda(4) = (loc2 - 1.0)*loc1
    dnda(5) = loc2*(loc1 - 1.0)
    dnda(6) = -loc2*(loc1 - 1.0)
    dnda(7) = loc2*loc1
    dnda(8) = -loc2*loc1

    dndb(1) = -(loc2 - 1.0)*(loc0 - 1.0)
    dndb(2) = (loc2 - 1.0)*loc0
    dndb(3) = -(loc2 - 1.0)*loc0
    dndb(4) = (loc2 - 1.0)*(loc0 - 1.0)
    dndb(5) = loc2*(loc0 - 1.0)
    dndb(6) = -loc2*loc0
    dndb(7) = loc2*loc0
    dndb(8) = -loc2*(loc0 - 1.0)

    dndc(1) = -(loc1 - 1.0)*(loc0 - 1.0)
    dndc(2) = (loc1 - 1.0)*loc0
    dndc(3) = -loc1*loc0
    dndc(4) = loc1*(loc0 - 1.0)
    dndc(5) = (loc1 - 1.0)*(loc0 - 1.0)
    dndc(6) = -(loc1 - 1.0)*loc0
    dndc(7) = loc1*loc0
    dndc(8) = -loc1*(loc0 - 1.0)

    return

  end subroutine b8_deriv_hpar

end module shape_3d_mod
