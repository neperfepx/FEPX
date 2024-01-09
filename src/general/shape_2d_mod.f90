! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module shape_2d_mod

! Shape function library for "2d" elements (faces of 3d elements). Contains full
!   deescriptions of 3-node (linear) and 6-node (quadratic) triangular
!   elements, and 4-node (linear), 8-node (quadratic) and 9-node quadrilateral
!   elements.

! Contains subroutines:
! sf2d: Evaluates 2d shape functions at a list of points
! sf2dg: Evaluates 2d shape function gradients at a list of points
! sf2d03: 3-node triangular shape function values
! sf2dg03: 3-node triangular shape function gradients
! sf2d04: 4-node quadrilateral shape function values
! sf2dg04: 4-node quadrilateral shape function gradients
! sf2d06: 6-node triangular shape function values
! sf2dg06: 6-node triangular shape function gradients
! sf2d08: 8-node quadrilateral shape function values
! sf2dg08: 8-node quadrilateral shape function gradients
! sf2d09: 9-node quadrilateral shape function values
! sf2dg09: 9-node quadrilateral shape function gradients

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod

  implicit none

contains

  subroutine sf2d(itype, npts, pts, vals, mv1, istat)

    ! Evaluate shape functions at a list of points.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! itype: Element type
    ! npts: Number of points
    ! pts: Points
    ! vals: Point values
    ! mv1:
    ! istat: Error flagging

    integer :: itype
    integer :: npts
    real(rk) :: pts(2, *)
    real(rk) :: vals(mv1, *)
    integer :: mv1
    integer :: istat

    !---------------------------------------------------------------------------

    istat = 0

    if (itype .eq. 3) then
      call sf2d03(npts, pts, vals, mv1)

    else if (itype .eq. 4) then
      call sf2d04(npts, pts, vals, mv1)

    else if (itype .eq. 6) then
      call sf2d06(npts, pts, vals, mv1)

    else if (itype .eq. 8) then
      call sf2d08(npts, pts, vals, mv1)

    else if (itype .eq. 9) then
      call sf2d09(npts, pts, vals, mv1)

    else
      ! Unrecognized element type.

      istat = 1
    end if

    return

  end subroutine sf2d

  !===========================================================================

  subroutine sf2dg(itype, npts, pts, grads, mg1, istat)

    ! Evaluate shape function gradients at a list of points.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! itype: Element type
    ! npts: Number of points
    ! pts: Points
    ! grads: Point gradients
    ! mg1:
    ! istat: Error flagging

    integer :: itype
    integer :: npts
    real(rk) :: pts(2, *)
    real(rk) :: grads(2, mg1, *)
    integer :: mg1
    integer :: istat

    !---------------------------------------------------------------------------

    istat = 0

    if (itype .eq. 3) then
      call sf2dg03(npts, pts, grads, mg1)

    else if (itype .eq. 4) then
      call sf2dg04(npts, pts, grads, mg1)

    else if (itype .eq. 6) then
      call sf2dg06(npts, pts, grads, mg1)

    else if (itype .eq. 8) then
      call sf2dg08(npts, pts, grads, mg1)

    else if (itype .eq. 9) then
      call sf2dg09(npts, pts, grads, mg1)

    else
      ! Unrecognized element type.

      istat = 1
    end if

    return

  end subroutine sf2dg

  !===========================================================================

  subroutine sf2d03(npts, pnts, value, mv1)

    ! 3-node triangular shape function values

    !---------------------------------------------------------------------------

    ! Arguments:
    ! npts: Number of points
    ! pnts: Points
    ! value: Point values
    ! mv1:

    integer :: npts
    real(rk) :: pnts(2, *)
    real(rk) :: value(mv1, *)
    integer :: mv1

    ! Locals:

    integer :: i
    real(rk) :: xi
    real(rk) :: eta
    real(rk) :: zeta

    !---------------------------------------------------------------------------

    do i = 1, npts
      xi = pnts(1, i)
      eta = pnts(2, i)
      zeta = 1.0d0 - xi - eta

      ! Nodal locations:

      ! 2
      ! 31

      value(1, i) = xi
      value(2, i) = eta
      value(3, i) = zeta
    end do

    return

  end subroutine sf2d03

  !===========================================================================

  subroutine sf2dg03(npts, pnts, grad, mg1)

    ! 3-node triangular shape function gradients

    !---------------------------------------------------------------------------

    ! Arguments:
    ! npts: Number of points
    ! pnts: Points
    ! grad: Point gradients
    ! mg1:

    integer :: npts
    real(rk) :: pnts(2, *)
    real(rk) :: grad(2, mg1, *)
    integer :: mg1

    ! Locals:

    integer :: i
    real(rk) :: xi
    real(rk) :: eta
    real(rk) :: zeta

    !---------------------------------------------------------------------------

    do i = 1, npts
      xi = pnts(1, i)
      eta = pnts(2, i)
      zeta = 1.0d0 - xi - eta

      grad(1, 1, i) = 1.0d0
      grad(1, 2, i) = 0.0d0
      grad(1, 3, i) = -1.0d0

      grad(2, 1, i) = 0.0d0
      grad(2, 2, i) = 1.0d0
      grad(2, 3, i) = -1.0d0
    end do

    return

  end subroutine sf2dg03

  !===========================================================================

  subroutine sf2d04(npts, pnts, value, mv1)

    ! 4-node quadrilateral shape function values

    !---------------------------------------------------------------------------

    ! Arguments:
    ! npts: Number of points
    ! pnts: Points
    ! value: Point values
    ! mv1:

    integer :: npts
    real(rk) :: pnts(2, *)
    real(rk) :: value(mv1, *)
    integer :: mv1

    ! Locals:

    integer :: i
    real(rk) :: xi
    real(rk) :: eta

    !---------------------------------------------------------------------------

    do i = 1, npts
      xi = pnts(1, i)
      eta = pnts(2, i)

      value(1, i) = (1.0d0 - xi)*(1.0d0 - eta)
      value(2, i) = xi*(1.0d0 - eta)
      value(3, i) = xi*eta
      value(4, i) = (1.0d0 - xi)*eta
    end do

    return

  end subroutine sf2d04

  !===========================================================================

  subroutine sf2dg04(npts, pnts, grad, mg1)

    ! 4-node quadrilateral shape function gradients

    !---------------------------------------------------------------------------

    ! Arguments:
    ! npts: Number of points
    ! pnts: Points
    ! grad: Point gradients
    ! mg1:

    integer :: npts
    real(rk) :: pnts(2, *)
    real(rk) :: grad(2, mg1, *)
    integer :: mg1

    ! Locals:

    integer :: i
    real(rk) :: xi
    real(rk) :: eta

    !---------------------------------------------------------------------------

    do i = 1, npts
      xi = pnts(1, i)
      eta = pnts(2, i)

      grad(1, 1, i) = eta - 1.0d0
      grad(1, 2, i) = 1.0d0 - eta
      grad(1, 3, i) = eta
      grad(1, 4, i) = -eta

      grad(2, 1, i) = xi - 1.0d0
      grad(2, 2, i) = -xi
      grad(2, 3, i) = xi
      grad(2, 4, i) = 1.0d0 - xi
    end do

    return

  end subroutine sf2dg04

  !===========================================================================

  subroutine sf2d06(npts, pnts, value, mv1)

    ! 6-node triangular shape function values

    !---------------------------------------------------------------------------

    ! Arguements:
    ! npts: Number of points
    ! pnts: Points
    ! value: Point values
    ! mv1:

    integer :: npts
    real(rk) :: pnts(2, *)
    real(rk) :: value(mv1, *)
    integer :: mv1

    ! Locals:

    integer :: i
    real(rk) :: xi
    real(rk) :: eta
    real(rk) :: zeta

    !---------------------------------------------------------------------------

    do i = 1, npts
      xi = pnts(1, i)
      eta = pnts(2, i)
      zeta = 1.0d0 - xi - eta

      ! Nodal locations:
      ! 3
      ! 42
      ! 561

      value(1, i) = (2.0d0*xi - 1.0d0)*xi
      value(2, i) = 4.0d0*eta*xi
      value(3, i) = (2.0d0*eta - 1.0d0)*eta
      value(4, i) = 4.0d0*eta*zeta
      value(5, i) = (2.0d0*zeta - 1.0d0)*zeta
      value(6, i) = 4.0d0*xi*zeta
    end do

    return

  end subroutine sf2d06

  !===========================================================================

  subroutine sf2dg06(npts, pnts, grad, mg1)

    ! 6-node triangular shape function gradients

    !---------------------------------------------------------------------------

    ! Arguments:
    ! npts: Number of points
    ! pnts: Points
    ! grad: Point gradients
    ! mg1:

    integer :: npts
    real(rk) :: pnts(2, *)
    real(rk) :: grad(2, mg1, *)
    integer :: mg1

    ! Locals:

    integer :: i
    real(rk) :: xi
    real(rk) :: eta
    real(rk) :: zeta

    !---------------------------------------------------------------------------

    do i = 1, npts
      xi = pnts(1, i)
      eta = pnts(2, i)
      zeta = 1.0d0 - xi - eta

      grad(1, 1, i) = 4.0d0*xi - 1.0d0
      grad(1, 2, i) = 4.0d0*eta
      grad(1, 3, i) = 0.0d0
      grad(1, 4, i) = -4.0d0*eta
      grad(1, 5, i) = -4.0d0*zeta + 1.0d0
      grad(1, 6, i) = 4.0d0*zeta - 4.0d0*xi

      grad(2, 1, i) = 0.0d0
      grad(2, 2, i) = 4.0d0*xi
      grad(2, 3, i) = 4.0d0*eta - 1.0d0
      grad(2, 4, i) = 4.0d0*zeta - 4.0d0*eta
      grad(2, 5, i) = -4.0d0*zeta + 1.0d0
      grad(2, 6, i) = -4.0d0*xi
    end do

    return

  end subroutine sf2dg06

  !===========================================================================

  subroutine sf2d08(npts, pnts, value, mv1)

    ! 8-node quadrilateral shape function values

    !---------------------------------------------------------------------------

    ! Arguments:
    ! npts: Number of points
    ! pnts: Points
    ! value: Point values
    ! mv1:

    integer :: npts
    real(rk) :: pnts(2, *)
    real(rk) :: value(mv1, *)
    integer :: mv1

    ! Locals:

    integer :: i
    real(rk) :: xi
    real(rk) :: eta

    !---------------------------------------------------------------------------

    do i = 1, npts
      xi = pnts(1, i)
      eta = pnts(2, i)

      value(1, i) = (1.0d0 - xi)*(1.0d0 - eta)*(-2.0d0*xi - 2.0d0 &
                                               & *eta + 1.0d0)
      value(2, i) = 4.0d0*xi*(1.0d0 - xi)*(1.0d0 - eta)
      value(3, i) = xi*(1.0d0 - eta)*(2.0d0*xi - 2.0d0*eta - &
          & 1.0d0)
      value(4, i) = 4.0d0*xi*eta*(1.0d0 - eta)
      value(5, i) = xi*eta*(2.0d0*xi + 2.0d0*eta - 3.0d0)
      value(6, i) = 4.0d0*xi*eta*(1.0d0 - xi)
      value(7, i) = eta*(1.0d0 - xi)*(-2.0d0*xi + 2.0d0*eta - &
          & 1.0d0)
      value(8, i) = 4.0d0*eta*(1.0d0 - xi)*(1.0d0 - eta)
    end do

    return

  end subroutine sf2d08

  !===========================================================================

  subroutine sf2dg08(npts, pnts, grad, mg1)

    ! 8-node quadrilateral shape function gradients

    !---------------------------------------------------------------------------

    ! Arguments:
    ! npts: Number of points
    ! pnts: Points
    ! grad: Point gradients
    ! mg1:

    integer :: npts
    real(rk) :: pnts(2, *)
    real(rk) :: grad(2, mg1, *)
    integer :: mg1

    ! Locals:

    integer :: i
    real(rk) :: xi
    real(rk) :: eta

    !---------------------------------------------------------------------------

    do i = 1, npts
      xi = pnts(1, i)
      eta = pnts(2, i)

      grad(1, 1, i) = -1.0d0*(1.0d0 - eta)*(3.0d0 - 4.0d0*xi - &
          & 2.0d0*eta)
      grad(1, 2, i) = 4.0d0*(1.0d0 - eta)*(1.0d0 - 2.0d0*xi)
      grad(1, 3, i) = (1.0d0 - eta)*(4.0d0*xi - 2.0d0*eta - 1.0d0)
      grad(1, 4, i) = 4.0d0*eta*(1.0d0 - eta)
      grad(1, 5, i) = eta*(4.0d0*xi + 2.0d0*eta - 3.0d0)
      grad(1, 6, i) = 4.0d0*eta*(1.0d0 - 2.0d0*xi)
      grad(1, 7, i) = -eta*(-4.0d0*xi + 2.0d0*eta + 1.0d0)
      grad(1, 8, i) = -4.0d0*eta*(1.0d0 - eta)

      grad(2, 1, i) = -1.0d0*(1.0d0 - xi)*(3.0d0 - 4.0d0*eta - &
          & 2.0d0*xi)
      grad(2, 2, i) = -4.0d0*xi*(1.0d0 - xi)
      grad(2, 3, i) = -xi*(2.0d0*xi - 4.0d0*eta + 1.0d0)
      grad(2, 4, i) = 4.0d0*xi*(1.0d0 - 2.0d0*eta)
      grad(2, 5, i) = xi*(2.0d0*xi + 4.0d0*eta - 3.0d0)
      grad(2, 6, i) = 4.0d0*xi*(1.0d0 - xi)
      grad(2, 7, i) = (1.0d0 - xi)*(-2.0d0*xi + 4.0d0*eta - 1.0d0)
      grad(2, 8, i) = 4.0d0*(1.0d0 - xi)*(1.0d0 - 2.0d0*eta)
    end do

    return

  end subroutine sf2dg08

  !===========================================================================

  subroutine sf2d09(npts, pnts, value, mv1)

    ! 9-node quadrilateral shape function values

    !---------------------------------------------------------------------------

    ! Arguments:
    ! npts: Number of points
    ! pnts: Points
    ! value: Point values
    ! mv1:

    integer :: npts
    real(rk) :: pnts(2, *)
    real(rk) :: value(mv1, *)
    integer :: mv1

    ! Locals:

    integer :: i
    real(rk) :: xi
    real(rk) :: eta

    !---------------------------------------------------------------------------

    do i = 1, npts
      xi = pnts(1, i)
      eta = pnts(2, i)

      value(1, i) = (2.0d0*xi - 1.0d0)*(2.0d0*eta - 1.0d0)*(xi - &
          & 1.0d0)*(eta - 1.0d0)
      value(2, i) = 4.0d0*xi*(1.0d0 - xi)*(eta - 1.0d0)*(2.0d0* &
          & eta - 1.0d0)
      value(3, i) = xi*(2.0d0*xi - 1.0d0)*(2.0d0*eta - 1.0d0)* &
          & (eta - 1.0d0)
      value(4, i) = 4.0d0*xi*eta*(2.0d0*xi - 1.0d0)*(1.0d0 - &
          & eta)
      value(5, i) = xi*eta*(2.0d0*xi - 1.0d0)*(2.0d0*eta - &
          & 1.0d0)
      value(6, i) = 4.0d0*xi*eta*(1.0d0 - xi)*(2.0d0*eta - &
          & 1.0d0)
      value(7, i) = (xi - 1.0d0)*(2.0d0*xi - 1.0d0)*(2.0d0*eta - &
          & 1.0d0)*eta
      value(8, i) = 4.0d0*(2.0d0*xi - 1.0d0)*(xi - 1.0d0)* &
          & (1.0d0 - eta)*eta
      value(9, i) = 16.0d0*xi*eta*(1.0d0 - xi)*(1.0d0 - eta)
    end do

    return

  end subroutine sf2d09

  !===========================================================================

  subroutine sf2dg09(npts, pnts, grad, mg1)

    ! 9-node quadrilateral shape function gradients

    !---------------------------------------------------------------------------

    ! Arguments:
    ! npts: Number of points
    ! pnts: Points
    ! grad: Point gradients
    ! mg1:

    integer :: npts
    real(rk) :: pnts(2, *)
    real(rk) :: grad(2, mg1, *)
    integer :: mg1

    ! Locals:

    integer :: i
    real(rk) :: xi
    real(rk) :: eta

    !---------------------------------------------------------------------------

    do i = 1, npts
      xi = pnts(1, i)
      eta = pnts(2, i)

      grad(1, 1, i) = (4.0d0*xi - 3.0d0)*(2.0d0*eta - 1.0d0)* &
          & (eta - 1.0d0)
      grad(1, 2, i) = 4.0d0*(1.0d0 - 2.0d0*xi)*(eta - 1.0d0)* &
          & (2.0d0*eta - 1.0d0)
      grad(1, 3, i) = (4.0d0*xi - 1.0d0)*(2.0d0*eta - 1.0d0)* &
          & (eta - 1.0d0)
      grad(1, 4, i) = 4.0d0*eta*(4.0d0*xi - 1.0d0)*(1.0d0 - eta)
      grad(1, 5, i) = eta*(4.0d0*xi - 1.0d0)*(2.0d0*eta - 1.0d0)
      grad(1, 6, i) = 4.0d0*eta*(1.0d0 - 2.0d0*xi)*(2.0d0*eta &
          & - 1.0d0)
      grad(1, 7, i) = (4.0d0*xi - 3.0d0)*(2.0d0*eta - 1.0d0)*eta
      grad(1, 8, i) = 4.0d0*eta*(4.0d0*xi - 3.0d0)*(1.0d0 - eta)
      grad(1, 9, i) = 16.0d0*(1.0d0 - 2.0d0*xi)*(1.0d0 - eta)*eta

      grad(2, 1, i) = (2.0d0*xi - 1.0d0)*(4.0d0*eta - 3.0d0)*(xi &
          & - 1.0d0)
      grad(2, 2, i) = 4.0d0*xi*(1.0d0 - xi)*(4.0d0*eta - 3.0d0)
      grad(2, 3, i) = xi*(2.0d0*xi - 1.0d0)*(4.0d0*eta - 3.0d0)
      grad(2, 4, i) = 4.0d0*xi*(2.0d0*xi - 1.0d0)*(1.0d0 - &
          & 2.0d0*eta)
      grad(2, 5, i) = xi*(2.0d0*xi - 1.0d0)*(4.0d0*eta - 1.0d0)
      grad(2, 6, i) = 4.0d0*xi*(1.0d0 - xi)*(4.0d0*eta - 1.0d0)
      grad(2, 7, i) = (xi - 1.0d0)*(2.0d0*xi - 1.0d0)*(4.0d0*eta &
          & - 1.0d0)
      grad(2, 8, i) = 4.0d0*(2.0d0*xi - 1.0d0)*(xi - 1.0d0)* &
          & (1.0d0 - 2.0d0*eta)
      grad(2, 9, i) = 16.0d0*xi*(1.0d0 - xi)*(1.0d0 - 2.0d0*eta)
    end do

    return

  end subroutine sf2dg09

end module shape_2d_mod
