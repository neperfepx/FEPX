! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE shape_3d_mod
  !
  ! Shape function library - 3d.
  !
  USE IntrinsicTypesModule, RK=>REAL_KIND
  !
  IMPLICIT NONE
  !
CONTAINS
!
!     Shape function routines for 4-node, 10-node tetrahedrons and 8-node brick
!     
!     4-NODE TETRAHEDRAON
!     --------------------------
!     sf4t_eval_vec: evaluate shape functions for vector valued quantity
!     .            : at array of points
!     T4_shape_hpar: shape functions for 4-node tetrahedron
!     T4_deriv_hpar: shape function derivatives for 4-node brick
!
!     10-NODE TETRAHEDRAON
!     --------------------------
!     sf10t_eval_vec: evaluate shape functions for vector valued quantity
!                   : at array of points
!     T10_shape_hpar: shape functions for 10-node tetrahedron
!     T10_deriv_hpar: shape function derivatives for 10-node tetrahedron
!
!     8-NODE BRICK
!     --------------------------
!     sf8b_eval_vec: evaluate shape functions for vector valued quantity
!     .            : at array of points
!     B8_shape_hpar: shape functions for 8-node brick
!     B8_deriv_hpar: shape function derivatives for 8-node brick
!
!     sfder_hpar   : compute jacobian determinant and global derivatives
!
!
!***********************************************************************
!
      SUBROUTINE sf10t_eval_vec(vec, npt, pnt, val)
!
!     Evaluate sf functions for vector valued quantity at array of points.
!
!---------------------------------------------------------------------72
!
      USE IntrinsicTypesModule, RK=>REAL_KIND
!
      IMPLICIT NONE
!
!     Arguments.
!
!     vec  -- values of vector at nodes
!     pnt  -- array of points in reference element
!     npt  -- number of points
!     val  -- values of vector at points
!
!     NOTE:  need explicit interface to use assumed shape arrays
!     .      as arguments, i.e. pnt(:), val(:); unless this is
!     .      put in a module or an interface, we can still use
!     .      the assumed-size arrays;
!
      INTEGER, INTENT(IN)   :: npt
      REAL(RK), INTENT(IN)  :: vec(3, 10), pnt(*)
      REAL(RK), INTENT(OUT) :: val(*)
!
!     Locals.
!
      INTEGER  :: i, xdof, ydof, zdof
      REAL(RK) :: x1, x2, x3, sf(10), tmpvec(3)
!
!---------------------------------------------------------------------72
!

      xdof = -2
      do i=1, npt
!
        xdof = xdof + 3; ydof = xdof + 1; zdof = ydof + 1
!
        x1 = pnt(xdof); x2 = pnt(ydof); x3 = pnt(zdof)
!        
        sf( 1) = 2.0d0*(x1+x2+x3-1.0d0)*(x1+x2+x3-0.5d0)
        sf( 2) = -4.0d0*(x1+x2+x3-1.0d0)*x1
        sf( 3) = 2.0d0*x1*(x1-0.5d0)
        sf( 4) = 4.0d0*x2*x1
        sf( 5) = 2.0d0*x2*(x2-0.5d0)
        sf( 6) = -4.0d0*(x1+x2+x3-1.0d0)*x2
        sf( 7) = -4.0d0*(x1+x2+x3-1.0d0)*x3
        sf( 8) = 4.0d0*x1*x3
        sf( 9) = 4.0d0*x2*x3
        sf(10) = 2.0d0*x3*(x3-0.5d0)   

!
        tmpvec = MATMUL(vec, sf)
!
        val(xdof) = tmpvec(1)  
        val(ydof) = tmpvec(2)  
        val(zdof) = tmpvec(3)
!
      enddo
!
      RETURN
      END SUBROUTINE
!
!**********************************************************************
!
      SUBROUTINE T10_shape_hpar(loc0, loc1, loc2, shape)
!
!----------------------------------------------------------------------
!
!     Shape functions for 10-node tetrahedron.
!    
!     looking down 3-axis
!     2
!     |
!     3--1
!    
!     top:
!    
!    
!     10
!     middle:
!    
!     9
!     78
!     bottom:
!     5
!     64
!     123
!    
!----------------------------------------------------------------------
!
      USE READ_INPUT_MOD, ONLY: NNPE
      USE IntrinsicTypesModule, RK=>REAL_KIND
!
      IMPLICIT NONE
!
!     Arguments.
!
!     loc0,loc1,loc2: coordinates of point
!     shape: shape function value at point, distributed over all elements
!
      REAL(RK)  ::  loc0
      REAL(RK)  ::  loc1
      REAL(RK)  ::  loc2
!
      REAL(RK)  ::  shape(0:nnpe)
!
!----------------------------------------------------------------------
!

      shape(0) = 2.0d0*(loc0+loc1+loc2-1.0d0)*(loc0+loc1+loc2-0.5d0)
      shape(1) = -4.0d0*(loc0+loc1+loc2-1.0d0)*loc0
      shape(2) = 2.0d0*loc0*(loc0-0.5d0)
      shape(3) = 4.0d0*loc1*loc0
      shape(4) = 2.0d0*loc1*(loc1-0.5d0)
      shape(5) = -4.0d0*(loc0+loc1+loc2-1.0d0)*loc1
      shape(6) = -4.0d0*(loc0+loc1+loc2-1.0d0)*loc2
      shape(7) = 4.0d0*loc0*loc2
      shape(8) = 4.0d0*loc1*loc2
      shape(9) = 2.0d0*loc2*(loc2-0.5d0)

      RETURN
      END SUBROUTINE
!
!**********************************************************************
!
      SUBROUTINE T10_deriv_hpar(loc0, loc1, loc2, dnda, dndb, dndc)
!
!     Shape function derivatives for 10-node tetrahedron.
!
!----------------------------------------------------------------------
!
!     Parameters and common blocks.
!     
      USE READ_INPUT_MOD, ONLY : nnpe, el_sub1, el_sup1
!
      IMPLICIT NONE
!
!     Arguments.
!
!     loc0,loc1,loc2: coordinates of point
!     dnda,dndb,dndc: shape function derivatives at point, distributed 
!     *             : over all elements
!
      REAL(RK)  ::  loc0
      REAL(RK)  ::  loc1
      REAL(RK)  ::  loc2
!
      REAL(RK)  ::  dnda(0:nnpe)
      REAL(RK)  ::  dndb(0:nnpe)
      REAL(RK)  ::  dndc(0:nnpe)
!
!----------------------------------------------------------------------
! 

      dnda(0) = 4.0d0*(loc0+loc1+loc2)-3.0d0
      dnda(1) = -4.0d0*(2.0d0*loc0+loc1+loc2-1.0d0)
      dnda(2) = 4.0d0*loc0-1.0d0
      dnda(3) = 4.0d0*loc1
      dnda(4) = 0.0d0
      dnda(5) = -4.0d0*loc1
      dnda(6) = -4.0d0*loc2
      dnda(7) = 4.0d0*loc2
      dnda(8) = 0.0d0
      dnda(9) = 0.0d0

      dndb(0) = 4.0d0*(loc0+loc1+loc2)-3.0d0
      dndb(1) = -4.0d0*loc0
      dndb(2) = 0.0d0
      dndb(3) = 4.0d0*loc0
      dndb(4) = 4.0d0*loc1-1.0d0
      dndb(5) = -4.0d0*(loc0+2.0d0*loc1+loc2-1.0d0)
      dndb(6) = -4.0d0*loc2
      dndb(7) = 0.0d0
      dndb(8) = 4.0d0*loc2
      dndb(9) = 0.0d0

      dndc(0) = 4.0d0*(loc0+loc1+loc2)-3.0d0
      dndc(1) = -4.0d0*loc0
      dndc(2) = 0.0d0
      dndc(3) = 0.0d0
      dndc(4) = 0.0d0
      dndc(5) = -4.0d0*loc1
      dndc(6) = -4.0d0*(loc0+loc1+2.0d0*loc2-1.0d0)
      dndc(7) = 4.0d0*loc0
      dndc(8) = 4.0d0*loc1
      dndc(9) = 4.0d0*loc2-1.0d0

      RETURN
      END SUBROUTINE
!     
!
!***********************************************************************
!
      SUBROUTINE sf4t_eval_vec(vec, npt, pnt, val)
!
!     Evaluate sf functions for vector valued quantity at array of points.
!
!---------------------------------------------------------------------72
!
      USE IntrinsicTypesModule, RK=>REAL_KIND
!
      IMPLICIT NONE
!
!     Arguments.
!
!     vec  -- values of vector at nodes
!     pnt  -- array of points in reference element
!     npt  -- number of points
!     val  -- values of vector at points
!
!     NOTE:  need explicit interface to use assumed shape arrays
!     .      as arguments, i.e. pnt(:), val(:); unless this is
!     .      put in a module or an interface, we can still use
!     .      the assumed-size arrays;
!
      INTEGER, INTENT(IN):: npt
      REAL(RK), INTENT(IN):: vec(3, 4), pnt(*)
      REAL(RK), INTENT(OUT):: val(*)
!
!     Locals.
!
      INTEGER  :: i, xdof, ydof, zdof
      REAL(RK) :: x1, x2, x3, sf(4), tmpvec(3)
!
!---------------------------------------------------------------------72
!

      xdof = -2
      do i=1, npt
!
        xdof = xdof + 3; ydof = xdof + 1; zdof = ydof + 1
!
        x1 = pnt(xdof); x2 = pnt(ydof); x3 = pnt(zdof)
!        
        sf(1) = x1
        sf(2) = x2
        sf(3) = x3
        sf(4) = 1 - (x1 + x2 + x3)

!
        tmpvec = MATMUL(vec, sf)
!
        val(xdof) = tmpvec(1)  
        val(ydof) = tmpvec(2)  
        val(zdof) = tmpvec(3)
!
      enddo
!
      RETURN
      END SUBROUTINE
!
!**********************************************************************
!
      SUBROUTINE T4_shape_hpar(loc0, loc1, loc2, shape)
!
!     Shape functions for 4-node tetrahedron.
!    
!     looking down 3-axis
!     2
!     |
!     3--1
!    
!     top:
!    
!    
!     2
!     bottom:
!     3
!     | \
!     4--1
!    
!----------------------------------------------------------------------
!
      USE READ_INPUT_MOD, ONLY: NNPE
!
      IMPLICIT NONE
!
!     Arguments.
!
!     loc0,loc1,loc2: coordinates of point
!     shape: shape function value at point, distributed over all elements
!
      REAL(RK)  ::  loc0
      REAL(RK)  ::  loc1
      REAL(RK)  ::  loc2
!
      REAL(RK)  ::  shape(0:nnpe)
!
!----------------------------------------------------------------------
!
      shape(0) = loc0
      shape(1) = loc1
      shape(2) = loc2
      shape(3) = 1.0 - (loc0 + loc1 + loc2)

      RETURN
      END SUBROUTINE
!
!**********************************************************************
!
      SUBROUTINE T4_deriv_hpar(loc0, loc1, loc2, dnda, dndb, dndc)
!
!     Shape function derivatives for 4-node tetrahedron.
!
!----------------------------------------------------------------------
!
!     Parameters and common blocks.
!     
      USE READ_INPUT_MOD, ONLY : nnpe, el_sub1, el_sup1
!
      IMPLICIT NONE
!
!     Arguments.
!
!     loc0,loc1,loc2: coordinates of point
!     dnda,dndb,dndc: shape function derivatives at point, distributed 
!     *             : over all elements
!
      REAL(RK)  ::  loc0
      REAL(RK)  ::  loc1
      REAL(RK)  ::  loc2
!
      REAL(RK)  ::  dnda(0:nnpe)
      REAL(RK)  ::  dndb(0:nnpe)
      REAL(RK)  ::  dndc(0:nnpe)
!
!----------------------------------------------------------------------
! 
      dnda(0) = 1.0d0
      dnda(1) = 0.0d0
      dnda(2) = 0.0d0
      dnda(3) = -1.0d0

      dndb(0) = 0.0d0
      dndb(1) = 0.0d0
      dndb(2) = 1.0d0
      dndb(3) = -1.0d0

      dndc(0) = 0.0d0
      dndc(1) = 1.0d0
      dndc(2) = 0.0d0
      dndc(3) = -1.0d0

      RETURN
      END SUBROUTINE
!     
!***********************************************************************
!
      SUBROUTINE sf8b_eval_vec(vec, npt, pnt, val)
!
!     Evaluate sf functions for vector valued quantity at array of points.
!
!---------------------------------------------------------------------72
!
      USE IntrinsicTypesModule, RK=>REAL_KIND
!
      IMPLICIT NONE
!
!     Arguments.
!
!     vec  -- values of vector at nodes
!     pnt  -- array of points in reference element
!     npt  -- number of points
!     val  -- values of vector at points
!
!     NOTE:  need explicit interface to use assumed shape arrays
!     .      as arguments, i.e. pnt(:), val(:); unless this is
!     .      put in a module or an interface, we can still use
!     .      the assumed-size arrays;
!
      INTEGER, INTENT(IN):: npt
      REAL(RK), INTENT(IN):: vec(3, 8), pnt(*)
      REAL(RK), INTENT(OUT):: val(*)
!
!     Locals.
!
      INTEGER :: i, xdof, ydof, zdof
      REAL(RK) :: x1, x2, x3, sf(8), tmpvec(3)
!
!---------------------------------------------------------------------72
!


      xdof = -2
      do i=1, npt
!
        xdof = xdof + 3; ydof = xdof + 1; zdof = ydof + 1
!
        x1 = pnt(xdof); x2 = pnt(ydof); x3 = pnt(zdof)
!        
        sf(1) = -(x3 - 1.0) * (x2 - 1.0) * (x1 - 1.0)
        sf(2) =  (x3 - 1.0) * (x2 - 1.0) * x1
        sf(3) = -(x3 - 1.0) * x2 * x1
        sf(4) =  (x3 - 1.0) * x2 * (x1 - 1.0)
        sf(5) =  x3 * (x2 - 1.0) * (x1 - 1.0)
        sf(6) = -x3 * (x2 - 1.0) * x1
        sf(7) =  x3 * x2 * x1
        sf(8) = -x3 * x2 * (x1 - 1.0)
!
        tmpvec = MATMUL(vec, sf)
!
        val(xdof) = tmpvec(1)  
        val(ydof) = tmpvec(2)  
        val(zdof)= tmpvec(3)
!
      enddo
!
      RETURN
      END SUBROUTINE
!
!**********************************************************************
!
      SUBROUTINE B8_shape_hpar(loc0, loc1, loc2, shape)
!
!     Shape functions for 8-node brick.
!
!----------------------------------------------------------------------
!
      USE READ_INPUT_MOD, ONLY: NNPE
!
      IMPLICIT NONE
!
!     Arguments.
!
!     loc0,loc1,loc2: coordinates of point
!     shape: shape function value at point, distributed over all elements
!
      REAL(RK)  ::  loc0
      REAL(RK)  ::  loc1
      REAL(RK)  ::  loc2
!
      REAL(RK)  ::  shape(0:nnpe)
!
!----------------------------------------------------------------------
!
      shape(0) = -(loc2 - 1.0) * (loc1 - 1.0) * (loc0 - 1.0)
      shape(1) =  (loc2 - 1.0) * (loc1 - 1.0) * loc0
      shape(2) = -(loc2 - 1.0) * loc1 * loc0
      shape(3) =  (loc2 - 1.0) * loc1 * (loc0 - 1.0)
      shape(4) =  loc2 * (loc1 - 1.0) * (loc0 - 1.0)
      shape(5) = -loc2 * (loc1 - 1.0) * loc0
      shape(6) =  loc2 * loc1 * loc0
      shape(7) = -loc2 * loc1 * (loc0 - 1.0)

      RETURN
      END SUBROUTINE
!
!**********************************************************************
!
      SUBROUTINE B8_deriv_hpar(loc0, loc1, loc2, dnda, dndb, dndc)
!
!     Shape function derivatives for 8-node brick.
!
!----------------------------------------------------------------------
!
!     Parameters and common blocks.
!     
      USE READ_INPUT_MOD, ONLY : nnpe, el_sub1, el_sup1
!
      IMPLICIT NONE
!
!     Arguments.
!
!     loc0,loc1,loc2: coordinates of point
!     dnda,dndb,dndc: shape function derivatives at point, distributed 
!     *             : over all elements
!
      REAL(RK)  ::  loc0
      REAL(RK)  ::  loc1
      REAL(RK)  ::  loc2
!
      REAL(RK)  ::  dnda(0:nnpe)
      REAL(RK)  ::  dndb(0:nnpe)
      REAL(RK)  ::  dndc(0:nnpe)
!
!----------------------------------------------------------------------
! 
      dnda(0) = -(loc2 - 1.0) * (loc1 - 1.0)
      dnda(1) =  (loc2 - 1.0) * (loc1 - 1.0)
      dnda(2) = -(loc2 - 1.0) * loc1
      dnda(3) =  (loc2 - 1.0) * loc1
      dnda(4) =  loc2 * (loc1 - 1.0)
      dnda(5) = -loc2 * (loc1 - 1.0)
      dnda(6) =  loc2 * loc1
      dnda(7) = -loc2 * loc1
  
      dndb(0) = -(loc2 - 1.0) * (loc0 - 1.0)
      dndb(1) =  (loc2 - 1.0) * loc0
      dndb(2) = -(loc2 - 1.0) * loc0
      dndb(3) =  (loc2 - 1.0) * (loc0 - 1.0)
      dndb(4) =  loc2 * (loc0 - 1.0)
      dndb(5) = -loc2 * loc0
      dndb(6) =  loc2 * loc0
      dndb(7) = -loc2 * (loc0 - 1.0)
 
      dndc(0) = -(loc1 - 1.0) * (loc0 - 1.0)
      dndc(1) =  (loc1 - 1.0) * loc0
      dndc(2) = -loc1 * loc0
      dndc(3) =  loc1 * (loc0 - 1.0)
      dndc(4) =  (loc1 - 1.0) * (loc0 - 1.0)
      dndc(5) = -(loc1 - 1.0) * loc0
      dndc(6) =  loc1 * loc0
      dndc(7) = -loc1 * (loc0 - 1.0)

      RETURN
      END SUBROUTINE
!     
!**********************************************************************
!     
      SUBROUTINE sfder_hpar(&
  &   loc0, loc1, loc2, coords, dndx, dndy, dndz, det,&
  &   ijac11, ijac12, ijac13, ijac21, ijac22, ijac23,&
  &   ijac31, ijac32, ijac33&
  &   )
!     
!     Compute quadrature quantities given a set of local coordinates.
!     
!----------------------------------------------------------------------
!     
      USE READ_INPUT_MOD, ONLY : el_sub1, el_sup1, kdim1, nnpe
      use units_mod
!
      IMPLICIT NONE
!     
!     Arguments.
!     
!     loc0,loc1,loc2: coordinates of point (in reference element)
!     coords        : real coordinates of element's nodes
!     det           : determinant of jacobian 
!     dndx,dndy,dndz: derivatives of mapped shape functions 
!     ijac{i}{j}    : components of inverse jacobian
!     
      REAL(RK)  ::  loc0
      REAL(RK)  ::  loc1
      REAL(RK)  ::  loc2

      REAL(RK)  ::  coords(0:kdim1, el_sub1:el_sup1), det(el_sub1:el_sup1)
      
      REAL(RK)  ::  dndx(0:nnpe, el_sub1:el_sup1)
      REAL(RK)  ::  dndy(0:nnpe, el_sub1:el_sup1)
      REAL(RK)  ::  dndz(0:nnpe, el_sub1:el_sup1)
      
      REAL(RK)  ::  ijac11(el_sub1:el_sup1), ijac12(el_sub1:el_sup1)
      REAL(RK)  ::  ijac13(el_sub1:el_sup1), ijac21(el_sub1:el_sup1)
      REAL(RK)  ::  ijac22(el_sub1:el_sup1), ijac23(el_sub1:el_sup1)
      REAL(RK)  ::  ijac31(el_sub1:el_sup1), ijac32(el_sub1:el_sup1)
      REAL(RK)  ::  ijac33(el_sub1:el_sup1)
!     
!     Locals.
!     
      INTEGER i, j, k, i0, i1, i2
!     
      REAL(RK)  ::  dnd1
      REAL(RK)  ::  dnd2
      REAL(RK)  ::  dnd3
!     
      REAL(RK)  ::  dnda(0:nnpe)
      REAL(RK)  ::  dndb(0:nnpe)
      REAL(RK)  ::  dndc(0:nnpe)

      REAL(RK)  ::  jac11(el_sub1:el_sup1)
      REAL(RK)  ::  jac12(el_sub1:el_sup1)
      REAL(RK)  ::  jac13(el_sub1:el_sup1)
!     
      REAL(RK)  ::  jac21(el_sub1:el_sup1)
      REAL(RK)  ::  jac22(el_sub1:el_sup1)
      REAL(RK)  ::  jac23(el_sub1:el_sup1)
!     
      REAL(RK)  ::  jac31(el_sub1:el_sup1)
      REAL(RK)  ::  jac32(el_sub1:el_sup1)
      REAL(RK)  ::  jac33(el_sub1:el_sup1)
!
      integer io
!     
!----------------------------------------------------------------------
!     
!     Evaluate shape function derivatives
!     
!     hardwired
!      call B8_deriv_hpar(loc0, loc1, loc2, dnda, dndb, dndc)
!      call T4_deriv_hpar(loc0, loc1, loc2, dnda, dndb, dndc)
      call T10_deriv_hpar(loc0, loc1, loc2, dnda, dndb, dndc)
!     
!     Zero the jacobian matrix
!     
      jac11 = 0.0
      jac12 = 0.0
      jac13 = 0.0
      jac21 = 0.0
      jac22 = 0.0
      jac23 = 0.0
      jac31 = 0.0
      jac32 = 0.0
      jac33 = 0.0

      do i = 0, nnpe
        i0 = 3 * i
        i1 = i0 + 1
        i2 = i1 + 1

        dnd1 = dnda(i)
        dnd2 = dndb(i)
        dnd3 = dndc(i)

!        do j = el_sub1,el_sup1
        ijac11(:) = coords(i0, :)
        ijac22(:) = coords(i1, :)
        ijac33(:) = coords(i2, :)
!        enddo

        jac11 = jac11 + ijac11 * dnd1
        jac21 = jac21 + ijac11 * dnd2
        jac31 = jac31 + ijac11 * dnd3
        jac12 = jac12 + ijac22 * dnd1
        jac22 = jac22 + ijac22 * dnd2
        jac32 = jac32 + ijac22 * dnd3
        jac13 = jac13 + ijac33 * dnd1
        jac23 = jac23 + ijac33 * dnd2
        jac33 = jac33 + ijac33 * dnd3
      end do
!     
!     Determinant of the jacobian matrix
!     
      ijac11 = jac11 * jac22 * jac33
      ijac12 = jac12 * jac23 * jac31
      ijac13 = jac13 * jac21 * jac32
      ijac21 = jac11 * jac23 * jac32
      ijac22 = jac12 * jac21 * jac33
      ijac23 = jac13 * jac22 * jac31
      
      det = ijac11 + ijac12 + ijac13
      det = det - (ijac21 + ijac22 + ijac23)
!     
!     Inverse of the jacobian matrix
!     
      ijac11 = jac22 * jac33 - jac23 * jac32
      ijac21 = jac23 * jac31 - jac21 * jac33
      ijac31 = jac21 * jac32 - jac22 * jac31
      ijac12 = jac13 * jac32 - jac12 * jac33
      ijac22 = jac11 * jac33 - jac13 * jac31
      ijac32 = jac31 * jac12 - jac11 * jac32
      ijac13 = jac12 * jac23 - jac13 * jac22
      ijac23 = jac13 * jac21 - jac11 * jac23
      ijac33 = jac11 * jac22 - jac12 * jac21
      
      ijac11 = ijac11 / det
      ijac12 = ijac12 / det
      ijac13 = ijac13 / det
      ijac21 = ijac21 / det
      ijac22 = ijac22 / det
      ijac23 = ijac23 / det
      ijac31 = ijac31 / det
      ijac32 = ijac32 / det
      ijac33 = ijac33 / det
!     
!     Finally the shape function derivatives
!     
      do i = 0, nnpe
        dnd1 = dnda(i)
        dnd2 = dndb(i)
        dnd3 = dndc(i)
        
        jac11 = ijac11 * dnd1 + ijac12 * dnd2 + ijac13 * dnd3
        jac22 = ijac21 * dnd1 + ijac22 * dnd2 + ijac23 * dnd3
        jac33 = ijac31 * dnd1 + ijac32 * dnd2 + ijac33 * dnd3
        
        dndx(i, :) = jac11
        dndy(i, :) = jac22
        dndz(i, :) = jac33
      end do
      
      RETURN
      END SUBROUTINE
!
!
END MODULE shape_3d_mod
