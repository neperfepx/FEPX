! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE Tensor3DModule 
!
!  *** Program Unit:  Module
!  ***    Unit Name:  Tensor3DModule
!  ***  Description:
!
!  This module handles operations on 3D tensors.
!
!  *** Use Statements:
!
USE IntrinsicTypesModule, &
     &  RK=>REAL_KIND, IK=>INTEGER_KIND, LK=>LOGICAL_KIND
USE ConstantsModule, ONLY:  &
     &  ZERO   => RK_ZERO,     ONE    => RK_ONE, TWO => RK_TWO, &
     &  HALF   => RK_ONE_HALF, THIRD  => RK_ONE_THIRD, &
     &  ROOT_2 => RK_ROOT_2,   ROOT_6 =>RK_ROOT_6
!
!  *** End:
!
IMPLICIT NONE
!
PRIVATE 
!
!  Decomposition conventions.
!
INTEGER, PARAMETER :: DECOMP_MPSIM=0, DECOMP_FEMEVPS=1
INTEGER :: DECOMP_DFLT = DECOMP_MPSIM
!!!INTEGER :: DECOMP_DFLT = DECOMP_FEMEVPS
!
!  Constants.
!
REAL(RK), PARAMETER :: SQ2_I = ONE/ROOT_2, SQ6_I = ONE/ROOT_6
REAL(RK), PARAMETER :: TWOSQ6_I = TWO*ROOT_6
!
!--------------Public Entities
!
!  *** Public Data:
!
!  Decomposition flags:
!
!  DECOMP_MPSIM, 
!  DECOMP_FEMEVPS -- conventions for five vector and skew part decompositions
!
PUBLIC :: DECOMP_MPSIM, DECOMP_FEMEVPS
!
!  *** Public Procedures:
!
PUBLIC :: Tensor3DDecompose, Tensor3Dcompose
!
!  *** End:
!
!--------------*-------------------------------------------------------
!
CONTAINS 
!
!  *** Program Unit:  subroutine
!  ***    Unit Name:  Tensor3DDecompose
!
!  *** Unit Declaration: 
!
SUBROUTINE Tensor3DDecompose(mat, DEV, SKW, SPH, &
     &   DECOMP)
  !
  !  ***  Description:  
  !
  !  Decompose matrix into deviatoric, skew, or spherical parts.
  !
  !  Note that the three components of the decomposition are
  !  orthogonal, but the underlying basis is not orthonormal
  !  due to scaling.
  !
  !  We are using the basis advocated by Tome' and Kocks, 1985.
  !
  !  Note:  dev(2) = sqrt(3/2)*mat(3,3), when mat is deviatoric
  !
  !  *** Arguments:
  !
  !  mat -- array of 3x3 matrices
  !
  REAL(RK), INTENT(IN) :: mat(:, :, :)
  !
  !  DEV -- array of 5-vectors representing symmetric, traceles parts
  !         (MPSIM convention)
  !
  REAL(RK), INTENT(OUT), OPTIONAL :: DEV(:, :)
  !
  !  SKW -- array of 3-vectors representing the axial vectors of 
  !         the skew parts (MPSIM convention)
  !
  REAL(RK), INTENT(OUT), OPTIONAL :: SKW(:, :)
  !
  !  SPH -- array of scalars representing one third of the trace
  !
  REAL(RK), INTENT(OUT), OPTIONAL :: SPH(:)
  !
  !
  !  DECOMP -- flag indicating decomposition convention
  !
  INTEGER, INTENT(IN), OPTIONAL :: DECOMP
  !
  !  *** Locals:
  !
  INTEGER ::  decomp_conv
  !
  !  *** End:
  !
  !--------------*-------------------------------------------------------
  !
  decomp_conv = DECOMP_DFLT
  if (PRESENT(DECOMP)) then
     decomp_conv = DECOMP
  end if
  !
  if (PRESENT(DEV)) then
     !
     SELECT CASE(decomp_conv)
        !
     CASE (DECOMP_MPSIM)
        !
        DEV(1,:) = ( mat(2,2,:) - mat(1,1,:) )*SQ2_I
        DEV(2,:) = ( mat(3,3,:) + mat(3,3,:) - mat(1,1,:) - mat(2,2,:) ) * SQ6_I
        !
        DEV(3,:) = SQ2_I *(mat(2,3,:)+mat(3,2,:))
        DEV(4,:) = SQ2_I *(mat(1,3,:)+mat(3,1,:))
        DEV(5,:) = SQ2_I *(mat(1,2,:)+mat(2,1,:))
        !
     CASE (DECOMP_FEMEVPS)
        !
        DEV(1,:) = ( mat(1,1,:) - mat(2,2,:) )*SQ2_I
        DEV(2,:) = ( mat(3,3,:) + mat(3,3,:) - mat(1,1,:) - mat(2,2,:) ) * SQ6_I
        !
        DEV(3,:) = SQ2_I *(mat(1,2,:)+mat(2,1,:))
        DEV(4,:) = SQ2_I *(mat(1,3,:)+mat(3,1,:))
        DEV(5,:) = SQ2_I *(mat(2,3,:)+mat(3,2,:))
        !
     CASE DEFAULT 
        !  Return error status
     END SELECT
     !
  end if
  !
  if (PRESENT(SKW)) then
     !
     SELECT CASE(decomp_conv)
        !
     CASE (DECOMP_MPSIM)
        !
        SKW(1,:) = HALF*(mat(3,2,:) - mat(2,3,:))
        SKW(2,:) = HALF*(mat(1,3,:) - mat(3,1,:))
        SKW(3,:) = HALF*(mat(2,1,:) - mat(1,2,:))
        !
     CASE (DECOMP_FEMEVPS)
        !
        SKW(1,:) =  HALF*(mat(2,1,:) - mat(1,2,:))
        SKW(2,:) = -HALF*(mat(1,3,:) - mat(3,1,:))
        SKW(3,:) =  HALF*(mat(3,2,:) - mat(2,3,:))
        !
     CASE DEFAULT 
        !  Return error status
     END SELECT
     !
  end if
  !
  if (PRESENT(SPH)) then
     SPH = ( mat(1,1,:) + mat(2,2,:) + mat(3,3,:) ) * THIRD
  end if
  !
END SUBROUTINE Tensor3DDecompose
!
!  *** Program Unit:  subroutine
!  ***    Unit Name:  Tensor3DCompose
!
!  *** Unit Declaration: 
!
SUBROUTINE Tensor3DCompose(mat, DEV, SKW, SPH)
  !
  !  ***  Description:  
  !
  !  Form  matrix from deviatoric, skew, or spherical parts.
  !  If no parts are passed, the matix is zeroed.
  !
  !  We are using the basis advocated by Tome' and Kocks, 1985.
  !
  !  *** Arguments:
  !
  !  mat -- the resulting array of 3x3 matrices 
  !
  REAL(RK), INTENT(OUT) :: mat(:, :, :)
  !
  !  DEV -- array of 5-vectors representing symmetric, traceles parts
  !
  REAL(RK), INTENT(IN), OPTIONAL :: DEV(:, :)
  !
  !  SKW -- array of 3-vectors representing the axial vectors of 
  !         the skew parts 
  !
  REAL(RK), INTENT(IN), OPTIONAL :: SKW(:, :)
  !
  !  SPH -- array of scalars representing one third of the trace
  !
  REAL(RK), INTENT(IN), OPTIONAL :: SPH(:)
  !
  !  *** Locals:
  !

  !  *** End:
  !
  !--------------*-------------------------------------------------------
  !
  mat = ZERO
  !
  if (PRESENT(DEV)) then
     !
     mat(1,1,:) = - DEV(1,:)*SQ2_I - DEV(2,:)*SQ6_I
     mat(2,2,:) = + DEV(1,:)*SQ2_I - DEV(2,:)*SQ6_I
     mat(3,3,:) =                + DEV(2,:)*TWOSQ6_I
     !     
     mat(2,3,:) = DEV(3,:) * SQ2_I
     mat(3,2,:) = mat(2,3,:)
     !     
     mat(1,3,:) = DEV(4,:) * SQ2_I
     mat(3,1,:) = mat(1,3,:)
     !     
     mat(1,2,:) = DEV(5,:) * SQ2_I
     mat(2,1,:) = mat(1,2,:)
     !
  end if
  !
  if (PRESENT(SKW)) then
     !
     mat(2,3,:) = mat(2,3,:) - SKW(1,:) 
     mat(3,2,:) = mat(3,2,:) + SKW(1,:) 
     !
     mat(1,3,:) = mat(1,3,:) + SKW(2,:)
     mat(3,1,:) = mat(3,1,:) - SKW(2,:)
     !
     mat(1,2,:) = mat(1,2,:) - SKW(3,:)
     mat(2,1,:) = mat(2,1,:) + SKW(3,:)
     !
  end if
  !
  if (PRESENT(SPH)) then
     !
     mat(1,1,:) = mat(1,1,:) + SPH
     mat(2,2,:) = mat(2,2,:) + SPH
     mat(3,3,:) = mat(3,3,:) + SPH
     !
  end if
  !
END SUBROUTINE Tensor3DCompose
END MODULE Tensor3DModule 
!
