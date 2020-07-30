! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE shape_2_mod
  !
  ! Shape function library - 2d.
  !
  USE IntrinsicTypesModule, RK=>REAL_KIND
  !
  IMPLICIT NONE
  !
CONTAINS
  !     
  !
  !***********************************************************************
  !
  SUBROUTINE sf2d(itype, npts, pts, vals, MV1, istat)
    !
    !     Evaluate shape functions at a list of points.
    !
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !     Arguments:
    !
    INTEGER :: itype, npts, MV1, istat
    REAL(RK) :: pts(2,*), vals(MV1, *)
    !
    !----------------------------------------------------------------------
    !
    istat = 0
    !
    IF (itype .eq. 3) THEN
       !
       call sf2d03(npts, pts, vals, MV1)
       !
    ELSEIF (itype .eq. 4) THEN
       !
       call sf2d04(npts, pts, vals, MV1)
       !
    ELSEIF (itype .eq. 6) THEN
       !
       call sf2d06(npts, pts, vals, MV1)
       !
    ELSEIF (itype .eq. 8) THEN
       !
       call sf2d08(npts, pts, vals, MV1)
       !
    ELSEIF (itype .eq. 9) THEN
       !
       call sf2d09(npts, pts, vals, MV1)
       !
    ELSE
       !
       !     Unrecognized element type.
       !
       istat = 1
       !
    ENDIF
    !
    RETURN
  END SUBROUTINE sf2d
  !
  !***********************************************************************
  !
  SUBROUTINE sf2dg(itype, npts, pts, grads, MG1, istat)
    !
    !     Evaluate shape function gradients at a list of points.
    !
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    !     Arguments:
    !
    INTEGER :: itype, npts, MG1, istat
    REAL(RK) :: pts(2,*), grads(2,MG1, *)
    !
    !----------------------------------------------------------------------
    !
    istat = 0
    !
    IF (itype .eq. 3) THEN
       !
       call sf2dg03(npts, pts, grads, MG1)
       !
    ELSEIF (itype .eq. 4) THEN
       !
       call sf2dg04(npts, pts, grads, MG1)
       !
    ELSEIF (itype .eq. 6) THEN
       !
       call sf2dg06(npts, pts, grads, MG1)
       !
    ELSEIF (itype .eq. 8) THEN
       !
       call sf2dg08(npts, pts, grads, MG1)
       !
    ELSEIF (itype .eq. 9) THEN
       !
       call sf2dg09(npts, pts, grads, MG1)
       !
    ELSE
       !
       !     Unrecognized element type.
       !
       istat = 1
       !
    ENDIF
    !
    RETURN
  END SUBROUTINE sf2dg
  !
  !***********************************************************************
  !
  !     ***** t3 ******
  !
  SUBROUTINE sf2d03(npts, pnts, value, MV1)
    !
    !     three node triangular element (values)
    !     
    IMPLICIT NONE
    !
    INTEGER  :: npts, MV1
    REAL(RK) :: value(MV1,*), pnts(2, *)
    !
    INTEGER  :: i
    REAL(RK) ::xi, eta, zeta
    !
    do i=1, npts
       xi  = pnts(1,i)
       eta = pnts(2,i)
       zeta = 1.d0 - xi - eta
       !
       ! nodal locations:
       !
       ! 2
       ! 31
       !
       value(1,i) = xi
       value(2,i) = eta
       value(3,i) = zeta
    enddo
    !
    RETURN
  END SUBROUTINE sf2d03
  !
  !***********************************************************************
  !
  SUBROUTINE sf2dg03(npts, pnts, grad, MG1)
    !
    !     three node tirangle element (gradients)
    !     
    IMPLICIT NONE
    !
    INTEGER  :: npts, MG1
    REAL(RK) :: grad(2,MG1,*), pnts(2,*)
    !
    INTEGER  :: i
    REAL(RK) :: xi, eta, zeta
    !
    do i=1, npts
       xi  = pnts(1,i)
       eta = pnts(2,i)
       zeta = 1.0d0 - xi - eta
       !
       grad(1,1,i) = 1.0d0
       grad(1,2,i) = 0.0d0
       grad(1,3,i) = -1.0d0
       !
       grad(2,1,i) = 0.d0
       grad(2,2,i) = 1.d0
       grad(2,3,i) = -1.d0
    enddo
    !
    RETURN
  END SUBROUTINE sf2dg03
  !
  !***********************************************************************
  !
  !     ***** q4 ******
  !
  SUBROUTINE sf2d04(npts, pnts, value, MV1)
    !
    !     four node quadrilateral element (values)
    !     
    IMPLICIT NONE
    !
    INTEGER  :: npts, MV1
    REAL(RK) :: value(MV1,*), pnts(2, *)
    !
    INTEGER  :: i
    REAL(RK) :: xi, eta
    !
    do i=1, npts
       xi  = pnts(1,i)
       eta = pnts(2,i)
       !
       value(1,i) = (1.0d0-xi)*(1.0d0-eta)
       value(2,i) = xi*(1.0d0-eta)
       value(3,i) = xi*eta
       value(4,i) = (1.0d0-xi)*eta
    enddo
    !
    RETURN
  END SUBROUTINE sf2d04
  !
  !***********************************************************************
  !
  SUBROUTINE sf2dg04(npts, pnts, grad, MG1)
    !
    !     four node quadrilateral element (gradients)
    !     
    IMPLICIT NONE
    !
    INTEGER  :: npts, MG1
    REAL(RK) :: grad(2,MG1,*), pnts(2,*)
    !
    INTEGER  :: i
    REAL(RK) :: xi, eta
    !
    do i=1, npts
       xi  = pnts(1,i)
       eta = pnts(2,i)
       !
       grad(1,1,i) = eta-1.0d0
       grad(1,2,i) = 1.0d0-eta
       grad(1,3,i) = eta
       grad(1,4,i) = -eta
       !
       grad(2,1,i) =xi-1.0d0
       grad(2,2,i) =-xi
       grad(2,3,i) =xi
       grad(2,4,i) =1.0d0-xi
    enddo
    !
    RETURN
  END SUBROUTINE sf2dg04
  !
  !***********************************************************************
  !
  !     ****** t6 ******
  !
  SUBROUTINE sf2d06(npts, pnts, value, MV1)
    !
    !     six node triangular element (values)
    !     
    IMPLICIT NONE
    !
    INTEGER  :: npts, MV1
    REAL(RK) :: value(MV1,*), pnts(2, *)
    !
    INTEGER  :: i
    REAL(RK) :: xi, eta, zeta
    !
    do i=1, npts
       xi  = pnts(1,i)
       eta = pnts(2,i)
       zeta = 1.0d0 - xi - eta
       ! nodal locations:
       !
       ! 3
       ! 42
       ! 561
       !
       !
       value(1,i) = (2.0d0*xi-1.0d0)*xi
       value(2,i) = 4.0d0*eta*xi
       value(3,i) = (2.0d0*eta-1.0d0)*eta
       value(4,i) = 4.0d0*eta*zeta
       value(5,i) = (2.0d0*zeta-1.0d0)*zeta
       value(6,i) = 4.0d0*xi*zeta
    enddo
    !
    RETURN
  END SUBROUTINE sf2d06
  !
  !***********************************************************************
  !
  SUBROUTINE sf2dg06(npts, pnts, grad, MG1)
    !
    !     six node triangular element (gradients)
    !     
    IMPLICIT NONE
    !
    INTEGER  :: npts, MG1
    REAL(RK) :: grad(2,MG1,*), pnts(2,*)
    !
    INTEGER  :: i
    REAL(RK) :: xi, eta, zeta
    !
    do i=1, npts
       xi  = pnts(1,i)
       eta = pnts(2,i)
       zeta = 1.0d0 - xi - eta
       !
       grad(1,1,i) = 4.0d0*xi-1.0d0
       grad(1,2,i) = 4.0d0*eta
       grad(1,3,i) = 0.0d0
       grad(1,4,i) = -4.0d0*eta
       grad(1,5,i) = -4.0d0*zeta+1.0d0
       grad(1,6,i) = 4.0d0*zeta-4.0d0*xi
       !
       grad(2,1,i) = 0.0d0
       grad(2,2,i) = 4.0d0*xi
       grad(2,3,i) = 4.0d0*eta-1.0d0
       grad(2,4,i) = 4.0d0*zeta-4.0d0*eta
       grad(2,5,i) = -4.0d0*zeta+1.0d0
       grad(2,6,i) = -4.0d0*xi
    enddo
    !
    RETURN
  END SUBROUTINE sf2dg06
  !
  !***********************************************************************
  !
  !     ***** q8 *****
  !
  SUBROUTINE sf2d08(npts, pnts, value, MV1)
    !
    !     eight-node quadilateral ( serendipity ) element (values)
    !
    IMPLICIT NONE
    !
    INTEGER  :: npts, MV1
    REAL(RK) :: value(MV1,*), pnts(2, *)
    !
    INTEGER  :: i
    REAL(RK) :: xi, eta
    !
    do i=1, npts
       xi  = pnts(1,i)
       eta = pnts(2,i)
       !
       value(1,i) = (1.0-xi)*(1.0-eta)*(-2.0*xi-2.0*eta+1.0)
       value(2,i) = 4.0*xi*(1.0-xi)*(1.0-eta)
       value(3,i) = xi*(1.0-eta)*(2.0*xi-2.0*eta-1.0)
       value(4,i) = 4.0*xi*eta*(1.0-eta)
       value(5,i) = xi*eta*(2.0*xi+2.0*eta-3.0)
       value(6,i) = 4.0*xi*eta*(1.0-xi)
       value(7,i) = eta*(1.0-xi)*(-2.0*xi+2.0*eta-1.0)
       value(8,i) = 4.0*eta*(1.0-xi)*(1.0-eta)
    enddo
    !
    RETURN
  END SUBROUTINE sf2d08
  !
  !***********************************************************************
  !
  SUBROUTINE sf2dg08(npts, pnts, grad, MG1)
    !
    !     eight-node quadilateral ( serendipity ) element (gradients)
    !
    IMPLICIT NONE
    !
    INTEGER  :: npts, MG1
    REAL(RK) :: grad(2,MG1,*), pnts(2,*)
    !
    INTEGER  :: i
    REAL(RK) :: xi, eta
    !
    do i=1, npts
       xi  = pnts(1,i)
       eta = pnts(2,i)
       !
       grad(1,1,i) = -1.0*(1.0-eta)*(3.0-4.0*xi-2.0*eta)
       grad(1,2,i) = 4.0*(1.0-eta)*(1.0-2.0*xi)
       grad(1,3,i) = (1.0-eta)*(4.0*xi-2.0*eta-1.0)
       grad(1,4,i) = 4.0*eta*(1.0-eta)
       grad(1,5,i) = eta*(4.0*xi+2.0*eta-3.0)
       grad(1,6,i) = 4.0*eta*(1.0-2.0*xi)
       grad(1,7,i) = -eta*(-4.0*xi+2.0*eta+1.0)
       grad(1,8,i) = -4.0*eta*(1.0-eta)
       !
       grad(2,1,i) = -1.0*(1.0-xi)*(3.0-4.0*eta-2.0*xi)
       grad(2,2,i) = -4.0*xi*(1.0-xi)
       grad(2,3,i) = -xi*(2.0*xi-4.0*eta+1.0)
       grad(2,4,i) = 4.0*xi*(1.0-2.0*eta)
       grad(2,5,i) = xi*(2.0*xi+4.0*eta-3.0)
       grad(2,6,i) = 4.0*xi*(1.0-xi)
       grad(2,7,i) = (1.0-xi)*(-2.0*xi+4.0*eta-1.0)
       grad(2,8,i) = 4.0*(1.0-xi)*(1.0-2.0*eta)
    enddo
    !
    RETURN
  END SUBROUTINE sf2dg08
  !
  !***********************************************************************
  !
  !     ***** q9 *****
  !
  SUBROUTINE sf2d09(npts, pnts, value, MV1)
    !
    !     nine-node quadilateral element ( lagrangian ) (values)
    ! 
    IMPLICIT NONE
    !
    INTEGER  :: npts, MV1
    REAL(RK) :: value(MV1,*), pnts(2, *)
    !
    INTEGER  :: i
    REAL(RK) :: xi, eta
    !
    do i=1, npts
       xi  = pnts(1,i)
       eta = pnts(2,i)
       !
       value(1,i) =  (2.0*xi-1.0)*(2.0*eta-1.0)*(xi-1.0)*(eta-1.0)
       value(2,i) =  4.0*xi*(1.0-xi)*(eta-1.0)*(2.0*eta-1.0)
       value(3,i) =  xi*(2.0*xi-1.0)*(2.0*eta-1.0)*(eta-1.0)
       value(4,i) =  4.0*xi*eta*(2.0*xi-1.0)*(1.0-eta)
       value(5,i) =  xi*eta*(2.0*xi-1.0)*(2.0*eta-1.0)
       value(6,i) =  4.0*xi*eta*(1.0-xi)*(2.0*eta-1.0)
       value(7,i) =  (xi-1.0)*(2.0*xi-1.0)*(2.0*eta-1.0)*eta
       value(8,i) =  4.0*(2.0*xi-1.0)*(xi-1.0)*(1.0-eta)*eta
       value(9,i) =  16.0*xi*eta*(1.0-xi)*(1.0-eta)
    enddo
    !
    RETURN
  END SUBROUTINE sf2d09
  !
  !***********************************************************************
  !
  SUBROUTINE sf2dg09(npts, pnts, grad, MG1)
    !
    !     nine-node quadilateral element ( lagrangian ) (gradients)
    ! 
    IMPLICIT NONE
    !
    INTEGER  :: npts, MG1
    REAL(RK) :: grad(2,MG1,*), pnts(2,*)
    !
    INTEGER  :: i
    REAL(RK) :: xi, eta
    !
    do i=1, npts
       xi  = pnts(1,i)
       eta = pnts(2,i)
       !
       grad(1,1,i) =  (4.0*xi-3.0)*(2.0*eta-1.0)*(eta-1.0)
       grad(1,2,i) =  4.0*(1.0-2.0*xi)*(eta-1.0)*(2.0*eta-1.0)
       grad(1,3,i) =  (4.0*xi-1.0)*(2.0*eta-1.0)*(eta-1.0)
       grad(1,4,i) =  4.0*eta*(4.0*xi-1.0)*(1.0-eta)
       grad(1,5,i) =  eta*(4.0*xi-1.0)*(2.0*eta-1.0)
       grad(1,6,i) =  4.0*eta*(1.0-2.0*xi)*(2.0*eta-1.0)
       grad(1,7,i) =  (4.0*xi-3.0)*(2.0*eta-1.0)*eta
       grad(1,8,i) =  4.0*eta*(4.0*xi-3.0)*(1.0-eta)
       grad(1,9,i) =  16.0*(1.0-2.0*xi)*(1.0-eta)*eta
       !
       grad(2,1,i) =  (2.0*xi-1.0)*(4.0*eta-3.0)*(xi-1.0)
       grad(2,2,i) =  4.0*xi*(1.0-xi)*(4.0*eta-3.0)
       grad(2,3,i) =  xi*(2.0*xi-1.0)*(4.0*eta-3.0)
       grad(2,4,i) =  4.0*xi*(2.0*xi-1.0)*(1.0-2.0*eta)
       grad(2,5,i) =  xi*(2.0*xi-1.0)*(4.0*eta-1.0)
       grad(2,6,i) =  4.0*xi*(1.0-xi)*(4.0*eta-1.0)
       grad(2,7,i) =  (xi-1.0)*(2.0*xi-1.0)*(4.0*eta-1.0)
       grad(2,8,i) =  4.0*(2.0*xi-1.0)*(xi-1.0)*(1.0-2.0*eta)
       grad(2,9,i) =  16.0*xi*(1.0-xi)*(1.0-2.0*eta)
    enddo
    !
    RETURN
  END SUBROUTINE sf2dg09
  !
  !
END MODULE shape_2_mod
