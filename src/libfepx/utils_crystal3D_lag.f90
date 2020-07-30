! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE UtilsCrystalModule
  !
  !  *** Program Unit:  Module
  !  ***    Unit Name:  UtilsCrystalModule
  !  ***  Description:
  !
  !  This module ....
  !
  !  *** Use Statements:
  !
  USE IntrinsicTypesModule, RK=>REAL_KIND
  USE DimsModule
  USE READ_INPUT_MOD
  USE quadrature_mod
  USE shape_3d_mod
  !
  !  *** End:
  !
  IMPLICIT NONE
  !
CONTAINS 
!
!     Crystal utility routines:
!
!     symm_vgr: Compute symmetric part of an array of velocity gradients.
!     skew_vgr: Compute skew part of an array of velocity gradients.
!     mat_vec_symm_ser: Convert 3x3 symmetric matrix to 5-vector.
!     mat_vec_symm: Convert array of 3x3 symmetric matrices to
!     *           : array of 5-vectors.
!     vec_mat_symm: Convert 5-vector to symmetric matrix.
!     vec_mat_symm_grn: Convert 5-vector to symmetric matrix across
!     *               : grains and elements.
!     mat_vec_skew_ser: Convert skew 3x3 matrix to 3-vector.
!     mat_vec_skew: Convert array of skew matrices to array of 3-vectors.
!     vec_mat_skew: Convert array of 3-vectors to array of skew matrices.
!     vec_mat_skew_grn: Convert 3-vectors to skew matrices across
!     *               : grains and elements.
!     rot_mat_symm: Convert rotation matrix c to an operator on 5-vectors
!     *           : which is equivalent to a change of coordinates.
!     rot_mat_skew: Construct 3x3 matrix acting on skew 3-vectors which
!     *           : is equivalent to a change of coordinates.
!     lattice_deform: Change coordinates from sample to lattice (for 5-vectors).
!     lattice_spin: Convert a skew 3-vector from sample to lattice coordinates.
!     mat_x_mat3: Matrix multiplication for arrays of 3x3 matrices. (c = a*b)
!     mat_x_mat5: Matrix multiplication for arrays of 5x5 matrices.
!     mat_x_matt3: Matrix multiplication by transpose. (3x3)
!     mat_x_matt5: Matrix multiplication by transpose for arrays of 5x5 matrices.
!     mat_x_mats3: Multiply array of 3x3 matrices by a fixed 3x3 matrix, tranposed.
!     vec_d_vec5: Multiply diagonal matrix times array of vectors. (5 dim)
!     mat_x_vec5: Multiply array of matrices times array of vectors. (5 dim)
!     norm_vec: Compute 2-norm for array of 5-vectors.
!     determinant_grn: Compute determinant of array of 3x3 tensors.
!     calc_elvol: Compute element volume
!
!**********************************************************************
!
      SUBROUTINE symm_vgr_ser(d, dkk, vgrad)
!
!     Compute symmetric part of velocity gradients.  
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!     Arguments:
!
!     d    : deviatoric part of velocity gradients (output)
!     dkk  : trace of velocity gradients           (output)
!     vgrad: array of velocity gradients           (input)
!     m    : number of elements                    (input)
!
      INTEGER m
!
      REAL(RK) vgrad(0:DIMS1, 0:DIMS1)
      REAL(RK) d(0:DIMS1, 0:DIMS1)
      REAL(RK) dkk
!
!----------------------------------------------------------------------
!c
      d(0, 0) = vgrad(0, 0)
      d(1, 1) = vgrad(1, 1)
      d(2, 2) = vgrad(2, 2)
      d(1, 0) = 0.5 * (vgrad(1, 0) + vgrad(0, 1))
      d(2, 0) = 0.5 * (vgrad(2, 0) + vgrad(0, 2))
      d(2, 1) = 0.5 * (vgrad(2, 1) + vgrad(1, 2))
      d(0, 1) = d(1, 0)
      d(0, 2) = d(2, 0)
      d(1, 2) = d(2, 1)

      dkk = d(0, 0) + d(1, 1) + d(2, 2)

      d(0, 0) = d(0, 0) - dkk / 3.0
      d(1, 1) = d(1, 1) - dkk / 3.0
      d(2, 2) = d(2, 2) - dkk / 3.0

      RETURN
      
      END SUBROUTINE symm_vgr_ser
!
!**********************************************************************
!
      SUBROUTINE symm_vgr(d, dkk, vgrad, m)
!
!     Compute symmetric part of an array of velocity gradients.  
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     d    : deviatoric part of velocity gradients (output)
!     dkk  : trace of velocity gradients           (output)
!     vgrad: array of velocity gradients           (input)
!     m    : number of elements                    (input)
!
      INTEGER m
!
      REAL(RK) vgrad(0:DIMS1, 0:DIMS1, 0:(m -1))
      REAL(RK) d(0:DIMS1, 0:DIMS1, 0:(m - 1))
      REAL(RK) dkk(0:(m - 1))
!
!----------------------------------------------------------------------
!
      d(0, 0, :) = vgrad(0, 0, :)
      d(1, 1, :) = vgrad(1, 1, :)
      d(2, 2, :) = vgrad(2, 2, :)
      d(1, 0, :) = 0.5 * (vgrad(1, 0, :) + vgrad(0, 1, :))
      d(2, 0, :) = 0.5 * (vgrad(2, 0, :) + vgrad(0, 2, :))
      d(2, 1, :) = 0.5 * (vgrad(2, 1, :) + vgrad(1, 2, :))
      d(0, 1, :) = d(1, 0, :)
      d(0, 2, :) = d(2, 0, :)
      d(1, 2, :) = d(2, 1, :)

      dkk = d(0, 0, :) + d(1, 1, :) + d(2, 2, :)

      d(0, 0, :) = d(0, 0, :) - dkk / 3.0
      d(1, 1, :) = d(1, 1, :) - dkk / 3.0
      d(2, 2, :) = d(2, 2, :) - dkk / 3.0

      RETURN
      !
      END SUBROUTINE symm_vgr
!
!**********************************************************************
!
      SUBROUTINE skew_vgr(w, vgrad, m)
!
!     Compute skew part of an array of velocity gradients.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     w    : array of skew matrices (output)
!     vgrad: array of velocity gradients
!     m    : number of elements
!
      INTEGER m
!
      REAL(RK)&
  &   vgrad(0:DIMS1, 0:DIMS1, 0:(m - 1)),&
  &   w    (0:DIMS1, 0:DIMS1, 0:(m - 1))
!
!     Locals:
!
!      INTEGER   i, j
!
!----------------------------------------------------------------------
!
      w(:, :, :) = 0.0_RK
 
      w(0, 1, :) = 0.5 * (vgrad(0, 1, :) - vgrad(1, 0, :))
      w(0, 2, :) = 0.5 * (vgrad(0, 2, :) - vgrad(2, 0, :))
      w(1, 2, :) = 0.5 * (vgrad(1, 2, :) - vgrad(2, 1, :))
      w(1, 0, :) = - w(0, 1, :)
      w(2, 0, :) = - w(0, 2, :)
      w(2, 1, :) = - w(1, 2, :)

      RETURN
      
      END SUBROUTINE skew_vgr
!
!**********************************************************************
!
      SUBROUTINE mat_vec_symm_ser(mat, vec)
!
!     Convert 3x3 symmetric matrix to 5-vector.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     mat: matrix
!     vec: vector
!
      REAL(RK) mat(0:DIMS1, 0:DIMS1)
      REAL(RK) vec(0:TVEC1)
!
!     Locals:
!
      REAL(RK) sqr2, sqr32
!
!----------------------------------------------------------------------
!
      sqr2  = dsqrt(2.d0)
      sqr32 = dsqrt(3.d0 / 2.d0)
!
      vec(0) = (mat(0, 0) - mat(1, 1)) / sqr2
      vec(1) = mat(2, 2) * sqr32
      vec(2) = mat(1, 0) * sqr2
      vec(3) = mat(2, 0) * sqr2
      vec(4) = mat(2, 1) * sqr2

      RETURN
      
      END SUBROUTINE mat_vec_symm_ser
!
!**********************************************************************
!
      SUBROUTINE mat_vec_symm(mat, vec, m)
!
!     Convert array of 3x3 symmetric matrices to array of 5-vectors.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
      INTEGER   m
      REAL(RK)    mat(0:DIMS1, 0:DIMS1, 0:(m - 1))
      REAL(RK)    vec(0:TVEC1, 0:(m - 1))
!
!     Locals:
      REAL(RK)    sqr2, sqr32
!
!----------------------------------------------------------------------
!
      sqr2=dsqrt(2.0d0)
      sqr32=dsqrt(1.5d0)
!
      vec(0, :) = (mat(0, 0, :) - mat(1, 1, :)) / sqr2
      vec(1, :) = mat(2, 2, :) * sqr32
      vec(2, :) = mat(1, 0, :) * sqr2
      vec(3, :) = mat(2, 0, :) * sqr2
      vec(4, :) = mat(2, 1, :) * sqr2

      RETURN
      
      END SUBROUTINE mat_vec_symm
!
!**********************************************************************
!
      SUBROUTINE vec5_vec6(vec, mat)
!
!     Convert 5-vector to 6-vector of symmetric matrix. 
!
!     Note: Returns the upper triangle in the order of: 11 22 33 23 13 12
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     vec: 5-vector of deviatoric part of second-order tensor - inner product preserved
!     mat: 6-vecotr of upper triangle of symmetric tensor - no inner prodct preservation
!
      REAL(RK)  mat(0:TVEC)
      REAL(RK)  vec(0:TVEC1)
!
!     Locals:
!
      REAL(RK)  sqr2, sqr23
!
!----------------------------------------------------------------------
!
      sqr2=dsqrt(2.0d0)
      sqr23=dsqrt(2.0d0/3.0d0)
!
!     diagonal -> 11 22 33
      mat(0) = 0.5 * (sqr2 * vec(0) - sqr23 * vec(1))
      mat(1) =-0.5 * (sqr2 * vec(0) + sqr23 * vec(1))
      mat(2) = vec(1) * sqr23
!
!     off-diagonal -> 23 13 12
      mat(3) = vec(4) / sqr2
      mat(4) = vec(3) / sqr2
      mat(5) = vec(2) / sqr2

      RETURN
      
      END SUBROUTINE vec5_vec6
!
!**********************************************************************
!
      SUBROUTINE vec_mat_symm_ser(vec, mat)
!
!     Convert 5-vector to symmetric matrix.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     vec: array of vectors
!     mat: array of 3x3 symmetric matrices
!
      REAL(RK)  mat(0:DIMS1, 0:DIMS1)
      REAL(RK)  vec(0:TVEC1)
!
!     Locals:
!
      REAL(RK)    sqr2, sqr23
!
!----------------------------------------------------------------------
!
      sqr2=dsqrt(2.0d0)
      sqr23=dsqrt(2.0d0/3.0d0)
!
      mat(0, 0) = 0.5 * (sqr2 * vec(0) - sqr23 * vec(1))
      mat(1, 1) =-0.5 * (sqr2 * vec(0) + sqr23 * vec(1))
      mat(2, 2) = vec(1) * sqr23
      mat(1, 0) = vec(2) / sqr2
      mat(2, 0) = vec(3) / sqr2
      mat(2, 1) = vec(4) / sqr2
              
      mat(0, 1) = mat(1, 0)
      mat(0, 2) = mat(2, 0)
      mat(1, 2) = mat(2, 1)

      RETURN
      
      END SUBROUTINE vec_mat_symm_ser
!
!**********************************************************************
!
      SUBROUTINE vec_mat_symm(vec, mat, m)
!
!     Convert 5-vector to symmetric matrix.
!     {5} --> [3x3]sym
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     vec: array of vectors
!     mat: array of 3x3 symmetric matrices
!     m: number of elements
!
      INTEGER , INTENT(IN)  :: m
      REAL(RK), INTENT(IN)  :: vec(0:TVEC1, 0:(m - 1))      
      REAL(RK), INTENT(OUT) :: mat(0:DIMS1, 0:DIMS1, 0:(m - 1))
!
!     Locals:
!
      REAL(RK) :: sqr2, sqr23
!
!----------------------------------------------------------------------
!
      sqr2=dsqrt(2.0d0)
      sqr23=dsqrt(2.0d0/3.0d0)
!
      mat(0, 0, :) = 0.5 * (sqr2 * vec(0, :) - sqr23 * vec(1, :))
      mat(1, 1, :) =-0.5 * (sqr2 * vec(0, :) + sqr23 * vec(1, :))
      mat(2, 2, :) = vec(1, :) * sqr23
      mat(1, 0, :) = vec(2, :) / sqr2
      mat(2, 0, :) = vec(3, :) / sqr2
      mat(2, 1, :) = vec(4, :) / sqr2
 
      mat(0, 1, :) = mat(1, 0, :)
      mat(0, 2, :) = mat(2, 0, :)
      mat(1, 2, :) = mat(2, 1, :)

      RETURN
      
      END SUBROUTINE vec_mat_symm
!
!**********************************************************************
!
      SUBROUTINE vec_mat_symm_grn(vec, mat, n, m)
!
!     Convert 5-vector to symmetric matrix across grains and elements.
!     {5} --> [3x3]sym
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     vec: array of vectors
!     mat: array of matrices
!     n: number of grains
!     m: number of elements
!
      INTEGER , INTENT(IN)  :: n, m
      REAL(RK), INTENT(IN)  :: vec(0:TVEC1, 0:(n - 1), 0:(m - 1))      
      REAL(RK), INTENT(OUT) :: mat(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      REAL(RK) :: sqr2, sqr23
!
!----------------------------------------------------------------------
!
      sqr2=dsqrt(2.0d0)
      sqr23=dsqrt(2.0d0/3.0d0)
!
      mat(0, 0, :, :) = 0.5 *&
  &                     (sqr2 * vec(0, :, :) - sqr23 * vec(1, :, :))
      mat(1, 1, :, :) =-0.5 *&
  &                     (sqr2 * vec(0, :, :) + sqr23 * vec(1, :, :))
      mat(2, 2, :, :) = vec(1, :, :) * sqr23
      mat(1, 0, :, :) = vec(2, :, :) / sqr2
      mat(2, 0, :, :) = vec(3, :, :) / sqr2
      mat(2, 1, :, :) = vec(4, :, :) / sqr2
 
      mat(0, 1, :, :) = mat(1, 0, :, :)
      mat(0, 2, :, :) = mat(2, 0, :, :)
      mat(1, 2, :, :) = mat(2, 1, :, :)
 
      RETURN
      
      END SUBROUTINE vec_mat_symm_grn
!
!**********************************************************************
!
      SUBROUTINE mat_vec_skew_ser(mat, vec)
!
!     Convert skew 3x3 matrix to 3-vector.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     mat: matrix
!     vec: vector
!
      REAL(RK) mat(0:DIMS1, 0:DIMS1)
      REAL(RK) vec(0:DIMS1)
!
!----------------------------------------------------------------------
!
      vec(0) = mat(1, 0)
      vec(1) = mat(2, 0)
      vec(2) = mat(2, 1)
 
      RETURN
      
      END SUBROUTINE mat_vec_skew_ser
!
!**********************************************************************
!
      SUBROUTINE mat_vec_skew(mat, vec, m)
!
!     Convert array of skew matrices to array of 3-vectors.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     mat: skew matrices
!     vec: 3-vectors
!     m: number of elements
!
      INTEGER m
      REAL(RK)  mat(0:DIMS1, 0:DIMS1, 0:(m - 1))
      REAL(RK)  vec(0:DIMS1, 0:(m - 1))
!
!----------------------------------------------------------------------
!
      vec(0, :) = mat(1, 0, :)
      vec(1, :) = mat(2, 0, :)
      vec(2, :) = mat(2, 1, :)
 
      RETURN
      
      END SUBROUTINE mat_vec_skew
!
!**********************************************************************
!
      SUBROUTINE vec_mat_skew(vec, mat, m)
!
!     Convert array of 3-vectors to array of skew matrices.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
      INTEGER m
      REAL(RK)  mat(0:DIMS1, 0:DIMS1, 0:(m - 1))
      REAL(RK)  vec(0:DIMS1, 0:(m - 1))
!
!----------------------------------------------------------------------
!
      mat(0, 0, :) = 0.
      mat(1, 1, :) = 0.
      mat(2, 2, :) = 0.
 
      mat(1, 0, :) = vec(0, :)
      mat(2, 0, :) = vec(1, :)
      mat(2, 1, :) = vec(2, :)
 
      mat(0, 1, :) = - vec(0, :)
      mat(0, 2, :) = - vec(1, :)
      mat(1, 2, :) = - vec(2, :)

      RETURN
      
      END SUBROUTINE vec_mat_skew
!
!**********************************************************************
!
      SUBROUTINE vec_mat_skew_grn(vec, mat, n, m)
!
!     Convert 3-vectors to skew matrices across grains and elements.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     vec: array of 3-vectors (input)
!     mat: array of 3x3 skew matrices (output)
!     n: number of grains
!     m: number of elements
!
      INTEGER   n, m
      REAL(RK)    mat(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK)    vec(0:DIMS1, 0:(n - 1), 0:(m - 1))
!
!----------------------------------------------------------------------
!
      mat(0, 0, :, :) = 0.
      mat(1, 1, :, :) = 0.
      mat(2, 2, :, :) = 0.
 
      mat(1, 0, :, :) = vec(0, :, :)
      mat(2, 0, :, :) = vec(1, :, :)
      mat(2, 1, :, :) = vec(2, :, :)
 
      mat(0, 1, :, :) = - vec(0, :, :)
      mat(0, 2, :, :) = - vec(1, :, :)
      mat(1, 2, :, :) = - vec(2, :, :)
 
      RETURN
      
      END SUBROUTINE vec_mat_skew_grn
!
!**********************************************************************
!
      SUBROUTINE rot_mat_symm(c, qr5x5, n, m)
!
!     Convert rotation matrix c to an operator on 5-vectors which is
!     equivalent to a change of coordinates.
!
!     c     --> qr5x5
!     [3x3] --> [5x5] 
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
      INTEGER n, m
      REAL(RK), INTENT(IN)  :: c(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT) :: qr5x5(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      REAL(RK) sqr3
!
      REAL(RK) c11(0:(n - 1), 0:(m - 1)), c12(0:(n - 1), 0:(m - 1))
      REAL(RK) c13(0:(n - 1), 0:(m - 1)), c21(0:(n - 1), 0:(m - 1))
      REAL(RK) c22(0:(n - 1), 0:(m - 1)), c23(0:(n - 1), 0:(m - 1))
      REAL(RK) c31(0:(n - 1), 0:(m - 1)), c32(0:(n - 1), 0:(m - 1))
      REAL(RK) c33(0:(n - 1), 0:(m - 1))
!
!----------------------------------------------------------------------
!
      sqr3=dsqrt(3.0d0)
!
!     construct 5x5 crystal rotation matrix: 
!                                 T
!             [A_sm]=[C][A_lat][C] <=>  {A_sm} = [Q]{A_lat}
!     with: {A}={()/sqr2,sqr32*(),sqr2*(),sqr2*(),sqr2*()}
!
      c11 = c(0, 0, :, :)
      c21 = c(1, 0, :, :)
      c31 = c(2, 0, :, :)
      c12 = c(0, 1, :, :)
      c22 = c(1, 1, :, :)
      c32 = c(2, 1, :, :)
      c13 = c(0, 2, :, :)
      c23 = c(1, 2, :, :)
      c33 = c(2, 2, :, :)

      qr5x5(0, 0, :, :)  =  0.5d0 * (c11 * c11 - c12 * c12 -&
  &                                           c21 * c21 + c22 * c22)
      qr5x5(0, 1, :, :)  =  sqr3 / 2.d0 * (c13 * c13 - c23 * c23)
      qr5x5(0, 2, :, :)  =  c11 * c12 - c21 * c22
      qr5x5(0, 3, :, :)  =  c11 * c13 - c21 * c23
      qr5x5(0, 4, :, :)  =  c12 * c13 - c22 * c23
      qr5x5(1, 0, :, :)  =  sqr3 / 2.d0 * (c31 * c31 - c32 * c32)
      qr5x5(1, 1, :, :)  =  1.5d0 * c33 * c33 - 0.5d0
      qr5x5(1, 2, :, :)  =  sqr3 * c31 * c32
      qr5x5(1, 3, :, :)  =  sqr3 * c31 * c33
      qr5x5(1, 4, :, :)  =  sqr3 * c32 * c33
      qr5x5(2, 0, :, :)  =  c11 * c21 - c12 * c22
      qr5x5(2, 1, :, :)  =  sqr3 * c13 * c23
      qr5x5(2, 2, :, :)  =  c11 * c22 + c12 * c21
      qr5x5(2, 3, :, :)  =  c11 * c23 + c13 * c21
      qr5x5(2, 4, :, :)  =  c12 * c23 + c13 * c22
      qr5x5(3, 0, :, :)  =  c11 * c31 - c12 * c32
      qr5x5(3, 1, :, :)  =  sqr3 * c13 * c33
      qr5x5(3, 2, :, :)  =  c11 * c32 + c12 * c31
      qr5x5(3, 3, :, :)  =  c11 * c33 + c13 * c31
      qr5x5(3, 4, :, :)  =  c12 * c33 + c13 * c32
      qr5x5(4, 0, :, :)  =  c21 * c31 - c22 * c32
      qr5x5(4, 1, :, :)  =  sqr3 * c23 * c33
      qr5x5(4, 2, :, :)  =  c21 * c32 + c22 * c31
      qr5x5(4, 3, :, :)  =  c21 * c33 + c23 * c31
      qr5x5(4, 4, :, :)  =  c22 * c33 + c23 * c32

      RETURN
      
      END SUBROUTINE rot_mat_symm
!
!**********************************************************************
!
      SUBROUTINE rot_mat_symm_ser(c, qr5x5)
!
!     Convert rotation matrix c to an operator on 5-vectors which is
!     equivalent to a change of coordinates.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
      REAL(RK)  c(0:DIMS1, 0:DIMS1)
      REAL(RK)  qr5x5(0:TVEC1, 0:TVEC1)
!
!     Locals:
!
      REAL(RK) sqr3
!     PARAMETER (sqr3=dsqrt(3.0d0))
!
      REAL(RK) c11, c12
      REAL(RK) c13, c21
      REAL(RK) c22, c23
      REAL(RK) c31, c32
      REAL(RK) c33
!
!----------------------------------------------------------------------
!
      sqr3=dsqrt(3.0d0)
!
!     construct 5x5 crystal rotation matrix: 
!                                 T
!             [A_sm]=[C][A_lat][C] <=>  {A_sm} = [Q]{A_lat}
!     with: {A}={()/sqr2,sqr32*(),sqr2*(),sqr2*(),sqr2*()}
!
      c11 = c(0, 0)
      c21 = c(1, 0)
      c31 = c(2, 0)
      c12 = c(0, 1)
      c22 = c(1, 1)
      c32 = c(2, 1)
      c13 = c(0, 2)
      c23 = c(1, 2)
      c33 = c(2, 2)

      qr5x5(0, 0)  =  0.5d0 * (c11 * c11 - c12 * c12 -&
  &                                           c21 * c21 + c22 * c22)
      qr5x5(0, 1)  =  sqr3 / 2.d0 * (c13 * c13 - c23 * c23)
      qr5x5(0, 2)  =  c11 * c12 - c21 * c22
      qr5x5(0, 3)  =  c11 * c13 - c21 * c23
      qr5x5(0, 4)  =  c12 * c13 - c22 * c23
      qr5x5(1, 0)  =  sqr3 / 2.d0 * (c31 * c31 - c32 * c32)
      qr5x5(1, 1)  =  1.5d0 * c33 * c33 - 0.5d0
      qr5x5(1, 2)  =  sqr3 * c31 * c32
      qr5x5(1, 3)  =  sqr3 * c31 * c33
      qr5x5(1, 4)  =  sqr3 * c32 * c33
      qr5x5(2, 0)  =  c11 * c21 - c12 * c22
      qr5x5(2, 1)  =  sqr3 * c13 * c23
      qr5x5(2, 2)  =  c11 * c22 + c12 * c21
      qr5x5(2, 3)  =  c11 * c23 + c13 * c21
      qr5x5(2, 4)  =  c12 * c23 + c13 * c22
      qr5x5(3, 0)  =  c11 * c31 - c12 * c32
      qr5x5(3, 1)  =  sqr3 * c13 * c33
      qr5x5(3, 2)  =  c11 * c32 + c12 * c31
      qr5x5(3, 3)  =  c11 * c33 + c13 * c31
      qr5x5(3, 4)  =  c12 * c33 + c13 * c32
      qr5x5(4, 0)  =  c21 * c31 - c22 * c32
      qr5x5(4, 1)  =  sqr3 * c23 * c33
      qr5x5(4, 2)  =  c21 * c32 + c22 * c31
      qr5x5(4, 3)  =  c21 * c33 + c23 * c31
      qr5x5(4, 4)  =  c22 * c33 + c23 * c32

      RETURN
      
      END SUBROUTINE rot_mat_symm_ser
!
!**********************************************************************
!
      SUBROUTINE rot_mat_skew(c, qr3x3, n, m)
!
!     Construct 3x3 matrix acting on skew 3-vectors which is equivalent
!     to a change of coordinates.
!     c [3x3] --> qr3x3 [3x3]
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     c: orientation matrix
!     qr3x3: 3x3 matrix which acts on skew vectors
!     n: number of grains
!     m: number of elements
!
      INTEGER , INTENT(IN)  :: n, m
      REAL(RK), INTENT(IN)  :: c(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT) :: qr3x3(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      REAL(RK) :: c11(0:(n - 1), 0:(m - 1)), c12(0:(n - 1), 0:(m - 1))
      REAL(RK) :: c13(0:(n - 1), 0:(m - 1)), c21(0:(n - 1), 0:(m - 1))
      REAL(RK) :: c22(0:(n - 1), 0:(m - 1)), c23(0:(n - 1), 0:(m - 1))
      REAL(RK) :: c31(0:(n - 1), 0:(m - 1)), c32(0:(n - 1), 0:(m - 1))
      REAL(RK) :: c33(0:(n - 1), 0:(m - 1))
!
!----------------------------------------------------------------------
!
!     Construct 3X3 rotation matrix for skew 2nd order tensors
!     [W]_sm = [c] [W]_lat [c]'  <=>  {W}_sm = [qr3x3] {W}_lat

      c11 = c(0, 0, :, :)
      c21 = c(1, 0, :, :)
      c31 = c(2, 0, :, :)
      c12 = c(0, 1, :, :)
      c22 = c(1, 1, :, :)
      c32 = c(2, 1, :, :)
      c13 = c(0, 2, :, :)
      c23 = c(1, 2, :, :)
      c33 = c(2, 2, :, :)
 
      qr3x3(0, 0, :, :) = c22 * c11 - c21 * c12
      qr3x3(0, 1, :, :) = c23 * c11 - c21 * c13
      qr3x3(0, 2, :, :) = c23 * c12 - c22 * c13
      qr3x3(1, 0, :, :) = c32 * c11 - c31 * c12
      qr3x3(1, 1, :, :) = c33 * c11 - c31 * c13
      qr3x3(1, 2, :, :) = c33 * c12 - c32 * c13
      qr3x3(2, 0, :, :) = c32 * c21 - c31 * c22
      qr3x3(2, 1, :, :) = c33 * c21 - c31 * c23
      qr3x3(2, 2, :, :) = c33 * c22 - c32 * c23

      RETURN
      
      END SUBROUTINE rot_mat_skew
!
!**********************************************************************
!
      SUBROUTINE rot_mat_skew_ser(c, qr3x3)
!
!     Construct 3x3 matrix acting on skew 3-vectors which is equivalent
!     to a change of coordinates.
!     c [3x3] --> qr3x3 [3x3]
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     c: orientation matrix
!     qr3x3: 3x3 matrix which acts on skew vectors
!
      REAL(RK), INTENT(IN)  :: c(0:DIMS1, 0:DIMS1)
      REAL(RK), INTENT(OUT) :: qr3x3(0:DIMS1, 0:DIMS1)
!
!     Locals:
!
      REAL(RK) :: c11
      REAL(RK) :: c12
      REAL(RK) :: c13
      REAL(RK) :: c21
      REAL(RK) :: c22
      REAL(RK) :: c23
      REAL(RK) :: c31
      REAL(RK) :: c32
      REAL(RK) :: c33
!
!----------------------------------------------------------------------
!
!     Construct 3X3 rotation matrix for skew 2nd order tensors
!     [W]_sm = [c] [W]_lat [c]'  <=>  {W}_sm = [qr3x3] {W}_lat

      c11 = c(0, 0)
      c21 = c(1, 0)
      c31 = c(2, 0)
      c12 = c(0, 1)
      c22 = c(1, 1)
      c32 = c(2, 1)
      c13 = c(0, 2)
      c23 = c(1, 2)
      c33 = c(2, 2)
 
      qr3x3(0, 0) = c22 * c11 - c21 * c12
      qr3x3(0, 1) = c23 * c11 - c21 * c13
      qr3x3(0, 2) = c23 * c12 - c22 * c13
      qr3x3(1, 0) = c32 * c11 - c31 * c12
      qr3x3(1, 1) = c33 * c11 - c31 * c13
      qr3x3(1, 2) = c33 * c12 - c32 * c13
      qr3x3(2, 0) = c32 * c21 - c31 * c22
      qr3x3(2, 1) = c33 * c21 - c31 * c23
      qr3x3(2, 2) = c33 * c22 - c32 * c23

      RETURN
      
      END SUBROUTINE rot_mat_skew_ser
!
!**********************************************************************
!
      SUBROUTINE lattice_deform_ser(qr5x5, vec, vec_lat)
!
!     Change coordinates from sample to lattice (for 5-vectors).
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
      REAL(RK)  qr5x5(0:TVEC1, 0:TVEC1)
      REAL(RK)  vec(0:TVEC1)
      REAL(RK)  vec_lat(0:TVEC1)
!
!     Locals:
!
      INTEGER  i, j
!
!----------------------------------------------------------------------
!
!  Rotate {vec}_sm -> {vec}_lat : {v_lat} = [qr5x5]' {v_sm}

      vec_lat = 0.0

      do i = 0, TVEC1
         do j = 0, TVEC1
            vec_lat(i) = vec_lat(i) +&
  &                              qr5x5(j, i) * vec(j)
         end do
      end do

      RETURN
      
      END SUBROUTINE lattice_deform_ser
!
!**********************************************************************
!
      SUBROUTINE lattice_deform(qr5x5, vec, vec_lat, n, m)
!
!     Change coordinates from sample to lattice (for 5-vectors).
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
      INTEGER, INTENT(IN)   :: n, m
      REAL(RK), INTENT(IN)  :: qr5x5(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: vec(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT) :: vec_lat(0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER :: i, j
!
!----------------------------------------------------------------------
!
!  Rotate {vec}_sm -> {vec}_lat : {v_lat} = [qr5x5]' {v_sm}

      vec_lat = 0.0

      do i = 0, TVEC1
         do j = 0, TVEC1
            vec_lat(i, :, :) = vec_lat(i, :, :) +&
  &                              qr5x5(j, i, :, :) * vec(j, :, :)
         end do
      end do

      RETURN
      
      END SUBROUTINE lattice_deform
!
!**********************************************************************
!
      SUBROUTINE lattice_spin(qr3x3, w_vec, w_vec_lat, n, m)
!
!     Convert a skew 3-vector from sample to lattice coordinates.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
      INTEGER, INTENT(IN)   :: n, m
      REAL(RK), INTENT(IN)  :: qr3x3(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: w_vec(0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT) :: w_vec_lat(0:DIMS1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER :: i, j
!
!----------------------------------------------------------------------
!
!     Rotate {w_vec}_sm -> {w_vec}_lat : {w_lat} = [qr3x3]'{w_sm}

      w_vec_lat = 0.0

      do i = 0, DIMS1
         do j = 0, DIMS1
            w_vec_lat(i, :, :) = w_vec_lat(i, :, :) +&
  &                               qr3x3(j, i, :, :) * w_vec(j, :, :)
         enddo
      enddo
!
      RETURN
      
      END SUBROUTINE lattice_spin
!
!**********************************************************************
!
      SUBROUTINE lattice_spin_ser(qr3x3, w_vec, w_vec_lat)
!
!     Convert a skew 3-vector from sample to lattice coordinates.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
      REAL(RK), INTENT(IN)  :: qr3x3(0:DIMS1, 0:DIMS1)
      REAL(RK), INTENT(IN)  :: w_vec(0:DIMS1)
      REAL(RK), INTENT(OUT) :: w_vec_lat(0:DIMS1)
!
!     Locals:
!
      INTEGER :: i, j
!
!----------------------------------------------------------------------
!
!     Rotate {w_vec}_sm -> {w_vec}_lat : {w_lat} = [qr3x3]'{w_sm}

      w_vec_lat = 0.0

      do i = 0, DIMS1
         do j = 0, DIMS1
            w_vec_lat(i) = w_vec_lat(i) +&
  &                               qr3x3(j, i) * w_vec(j)
         enddo
      enddo
!
      RETURN
      
      END SUBROUTINE lattice_spin_ser
!
!**********************************************************************
!
      SUBROUTINE vec_d_vec5(d, b, c, n, m)
!
!     c = d*b where c and b are 5-vectors and d is a diagonal 5x5 matrix
!     represented in compact form
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     d   : a fixed diagonal matrix (in 5-vector form)
!     b   : input matrix            (in 5-vector form)
!     c   : output  matrix (d*b)    (in 5-vector form)
!     n   : number of grains
!     m   : number of elements
!
      INTEGER n, m
!
      REAL(RK) d(0:TVEC1)
      REAL(RK) b(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK) c(0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
      INTEGER k
!
!----------------------------------------------------------------------
!
      DO k = 0, TVEC1
        c(k, :, :) = d(k)*b(k, :, :)
      ENDDO

      RETURN
      
      END SUBROUTINE vec_d_vec5
!
!**********************************************************************
!
!
      SUBROUTINE mat_x_mat3(a, b, c, n, m)
!
!     Matrix multiplication for arrays of 3x3 matrices. (c = a*b)
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     a,b: input matrices
!     c: output (a*b)
!     n: number of grains
!     m: number of elements
!
      INTEGER n, m
      REAL(RK) a(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK) b(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK) c(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
      INTEGER i, j, k
!
!----------------------------------------------------------------------
!
      c = 0.0

      do j = 0, DIMS1
         do k = 0, DIMS1
            do i = 0, DIMS1
               c(i, j, :, :) = c(i, j, :, :) +&
  &                                   a(i, k, :, :) * b(k, j, :, :)
            enddo
         enddo
      enddo
   
      RETURN
      
      END SUBROUTINE mat_x_mat3
!
!**********************************************************************
!
      SUBROUTINE mat_x_mat5(a, b, c, n, m)
!
!     Matrix multiplication for arrays of 5x5 matrices.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     a,b: arrays of 5x5 matrices (input)
!     c: array of 5x5 matrices (c=a*b) (output)
!     n: number of grains
!     m: number of elements
!
      INTEGER n, m
      REAL(RK) a(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK) b(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK) c(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
      INTEGER i, j, k
!
!----------------------------------------------------------------------
!
      c = 0.0_RK

      do j = 0, TVEC1
         do k = 0, TVEC1
            do i = 0, TVEC1
               c(i, j, :, :) = c(i, j, :, :) +&
  &                                   a(i, k, :, :) * b(k, j, :, :)
            enddo
         enddo
      enddo
   
      RETURN
      
      END SUBROUTINE mat_x_mat5
!
!**********************************************************************
!
      SUBROUTINE mat_x_matt3(a, b, c, n, m)
!
!     Matrix multiplication by transpose. (3x3)
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     a,b: arrays of 3x3 matrices
!     c: array of 3x3 matrices, c=a*b^t
!     n: number of grains
!     m: number of elements
!
      INTEGER n, m
!
      REAL(RK) a(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK) b(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK) c(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
      INTEGER i, j, k
!
!----------------------------------------------------------------------
!
      c = 0.0_RK

      do j = 0, DIMS1
         do k = 0, DIMS1
            do i = 0, DIMS1
               c(i, j, :, :) = c(i, j, :, :) +&
  &                                   a(i, k, :, :) * b(j, k, :, :)
            enddo
         enddo
      enddo
   
      RETURN
      
      END SUBROUTINE mat_x_matt3
!
!**********************************************************************
!
      SUBROUTINE matt_x_mat3(a, b, c, n, m)
!
!     Matrix multiplication by transpose. (3x3)
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     a,b: arrays of 3x3 matrices
!     c: array of 3x3 matrices, c=a^t*b
!     n: number of grains
!     m: number of elements
!
      INTEGER n, m
!
      REAL(RK) a(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK) b(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK) c(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
      INTEGER i, j, k
!
!----------------------------------------------------------------------
!
      c = 0.0_RK

      do j = 0, DIMS1
         do i = 0, DIMS1
            do k = 0, DIMS1
               c(i, j, :, :) = c(i, j, :, :) +&
  &                                   a(k, i, :, :) * b(k, j, :, :)
            enddo
         enddo
      enddo
   
      RETURN
      
      END SUBROUTINE matt_x_mat3
!
!**********************************************************************
!
      SUBROUTINE mat_x_matt5(a, b, c, n, m)
!
!     Matrix multiplication by transpose for arrays of 5x5 matrices.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     a,b: arrays of 5x5 matrices (input)
!     c: array of 5x5 matrices; c=a*b^t (output)
!     n: number of grains
!     m: number of elements
!
      INTEGER n, m
      REAL(RK)  a(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  b(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  c(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
      INTEGER  i, j, k
!
!----------------------------------------------------------------------
!
      c = 0.0_RK

      do k = 0, TVEC1
         do j = 0, TVEC1
            do i = 0, TVEC1
               c(i, j, :, :) = c(i, j, :, :) +&
  &                                   a(i, k, :, :) * b(j, k, :, :)
            enddo
         enddo
      enddo
   
      RETURN
      
      END SUBROUTINE mat_x_matt5
!
!**********************************************************************
!
      SUBROUTINE mat_x_mats3(a, b, c, n, m)
!
!     Multiply array of 3x3 matrices by a fixed 3x3 matrix, tranposed.
!
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!
!     Arguments:
!
!     a: array of matrices to be multiplied
!     b: fixed 3x3 matrix
!     c: array of 3x3 matrices, c=a*b^t
!     n: number of grains
!     m: number of elements
!
      INTEGER n, m
      REAL(RK)  a(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  b(0:DIMS1, 0:DIMS1)
      REAL(RK)  c(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
      INTEGER i, j, k
!
!----------------------------------------------------------------------
!
      c = 0.0_RK

      do k = 0, DIMS1
         do j = 0, DIMS1
            do i = 0, DIMS1
               c(i, j, :, :) = c(i, j, :, :) +&
  &                                   a(i, k, :, :) * b(j, k)
            enddo
         enddo
      enddo
   
      RETURN
      
      END SUBROUTINE mat_x_mats3
!
!**********************************************************************
!
      SUBROUTINE mat_x_vec5_ser(matrix, vector, product)
!
!     Multiply array of matrices times array of vectors. (5 dim)
!     
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     matrix: array of matrices (5x5)
!     vector: array of 5-vectors
!     product: array of 5-vectors; product=matrix * vector
!
      REAL(RK)    matrix(0:TVEC1, 0:TVEC1)
      REAL(RK)    vector(0:TVEC1)
      REAL(RK)    product(0:TVEC1)
!
!     Locals:
      INTEGER   i, j
!
!----------------------------------------------------------------------
!
      product = 0.0_RK
 
      do j = 0, TVEC1
         do i = 0, TVEC1
            product(i) = product(i) +&
  &                            matrix(i, j) * vector(j)
         enddo
      enddo
 
      RETURN
      
      END SUBROUTINE mat_x_vec5_ser
!
!**********************************************************************
!
      SUBROUTINE mat_x_vec5(matrix, vector, product, n, m)
!
!     Multiply array of matrices times array of vectors. (5 dim)
!     
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     matrix: array of matrices (5x5)
!     vector: array of 5-vectors
!     product: array of 5-vectors; product=matrix * vector
!     n: number of grains
!     m: number of elements
!
      INTEGER, INTENT(IN)   :: n, m
      REAL(RK), INTENT(IN)  :: matrix(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: vector(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT) :: product(0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
      INTEGER  ::  i, j
!
!----------------------------------------------------------------------
!
      product = 0.0_RK
 
      do j = 0, TVEC1
         do i = 0, TVEC1
            product(i, :, :) = product(i, :, :) +&
  &                            matrix(i, j, :, :) * vector(j, :, :)
         enddo
      enddo
 
      RETURN
      
      END SUBROUTINE mat_x_vec5
!
!**********************************************************************
!
      SUBROUTINE norm_vec(norm, vector, n, m)
!
!     Compute 2-norm for array of 5-vectors.
!
!----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
      INTEGER   n, m
      REAL(RK)    norm(0:(n - 1), 0:(m - 1))
      REAL(RK)    vector(0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER   i
!
!----------------------------------------------------------------------
!
      norm = 0.0_RK
      do i = 0, TVEC1
         norm = norm + vector(i, :, :) * vector(i, :, :)
      enddo
      norm = dsqrt(norm)
 
      RETURN
      
      END SUBROUTINE norm_vec
!
!**********************************************************************
!
      SUBROUTINE determinant_grn(tensor, determ, n, m)
!
!     Compute determinant of array of 3x3 tensors.
!     
!-----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
!     tensor: array of 3x3 matrices
!     determ: array of scalars (determinant of tensor)
!     n: number of grains
!     m: number of elements
!
      INTEGER   n, m
      REAL(RK)    determ(0:(n - 1), 0:(m - 1))
      REAL(RK)    tensor(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
      INTEGER   i
      REAL(RK)    det11(0:(n - 1), 0:(m - 1)), det12(0:(n - 1), 0:(m - 1))
      REAL(RK)    det13(0:(n - 1), 0:(m - 1)), det21(0:(n - 1), 0:(m - 1))
      REAL(RK)    det22(0:(n - 1), 0:(m - 1)), det23(0:(n - 1), 0:(m - 1))
!
!-----------------------------------------------------------------------
!
!     Determinant of a second order tensor (3x3 matrix)
 
      det11 = tensor(0, 0, :, :)* tensor(1, 1, :, :)* tensor(2, 2, :, :)
      det12 = tensor(0, 1, :, :)* tensor(1, 2, :, :)* tensor(2, 0, :, :)
      det13 = tensor(0, 2, :, :)* tensor(1, 0, :, :)* tensor(2, 1, :, :)
      det21 = tensor(0, 0, :, :)* tensor(1, 2, :, :)* tensor(2, 1, :, :)
      det22 = tensor(0, 1, :, :)* tensor(1, 0, :, :)* tensor(2, 2, :, :)
      det23 = tensor(0, 2, :, :)* tensor(1, 1, :, :)* tensor(2, 0, :, :)
      
      determ = det11 + det12 + det13 - det21 - det22 - det23
 
      RETURN
      
      END SUBROUTINE determinant_grn
      
!
!**********************************************************************
!
      SUBROUTINE find_indices(numind, target, vector, indices)
!
!     Find the members of vector (integer) which are equal to target (integer), 
!     and return indices and numind
!
!----------------------------------------------------------------------
!
      USE DimsModule
!
      IMPLICIT NONE
!
!
!     Arguments:
!
!
      INTEGER :: target, vector(:)
      INTEGER, pointer :: indices(:)
      INTEGER :: numind
!
!     Locals:
!
      integer :: k, lb, ub, i

!
!----------------------------------------------------------------------
!
!-tm:
      NULLIFY(indices)
!
!      see if indices has been allocated. if so, deallocate
      if (associated(indices)) then
         deallocate(indices)
      endif
!      find how many indices there are
      numind=count(vector.EQ.target)
!      allocate indices to be numind long
      allocate(indices(1:numind))
!      find which indices of vector equal target, and save them in indices
      k=0
      lb=lbound(vector,1) ! always equals 1
      ub=ubound(vector,1)
      do i=lb,ub
         if (vector(i).EQ.target) then
            k=k+1
            indices(k)=i
         endif
      enddo
      ! 0-based
      indices=indices-1
      
      RETURN
      
      END SUBROUTINE find_indices
!      
!**********************************************************************
!
      SUBROUTINE calc_elvol(elvol, ecoords)
!
!     Compute elemental volumes
!     
!-----------------------------------------------------------------------
!
      IMPLICIT  NONE
!
!     Arguments:
!
!     elvol: array of elemental volumes
!     ecoords: coordinates of elemental nodal points
!
      REAL(RK), INTENT(OUT) :: elvol(el_sub1:el_sup1)
      REAL(RK), INTENT(IN)  :: ecoords(0:kdim1, el_sub1:el_sup1)
!
!     Locals:
!
      REAL(RK)  ::  loc0
      REAL(RK)  ::  loc1
      REAL(RK)  ::  loc2
      REAL(RK)  ::  wt(el_sub1:el_sup1)

      REAL(RK)  ::  det(el_sub1:el_sup1)
      
      REAL(RK)  ::  dndx(0:nnpe, el_sub1:el_sup1)
      REAL(RK)  ::  dndy(0:nnpe, el_sub1:el_sup1)
      REAL(RK)  ::  dndz(0:nnpe, el_sub1:el_sup1)
      
      REAL(RK)  ::  s11(el_sub1:el_sup1), s12(el_sub1:el_sup1)
      REAL(RK)  ::  s13(el_sub1:el_sup1), s21(el_sub1:el_sup1)
      REAL(RK)  ::  s22(el_sub1:el_sup1), s23(el_sub1:el_sup1)
      REAL(RK)  ::  s31(el_sub1:el_sup1), s32(el_sub1:el_sup1)
      REAL(RK)  ::  s33(el_sub1:el_sup1)

      INTEGER   ::  i
!
!-----------------------------------------------------------------------
!
      elvol = 0.0_RK

      do i = 0, nqpt1
         loc0 = qploc(0, i)
         loc1 = qploc(1, i)
         loc2 = qploc(2, i)
         wt = wtqp(0, i)

         call sfder_hpar(loc0, loc1, loc2, ecoords, dndx, dndy,&
              &      dndz, det, s11, s12, s13, s21, s22, s23, s31, s32, s33)
         
         elvol = elvol + det * wt

      end do

      RETURN

    END SUBROUTINE calc_elvol     
!
!**********************************************************************
!
      
END MODULE UtilsCrystalModule
