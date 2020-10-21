! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE ORIENTATION_CONVERSION_MOD
!
! Module containing orientation conversions for input and printing.
!
! Note: The conversions here assume no convention ("active" or "passive"). The
! default assumption in FEPX is a "passive" convention. If an "active" 
! convention is chosen by the user, the internal rotation matrix is transposed 
! when converted from it's original orientation convention, and again before 
! output to files. FEPX assumes "active" to mean "sample-to-crystal" 
! transformation, and "passive" to mean "crystal-to-sample" transformation. 
!
! All conversions are summarized in a paper by Rowenhorst et. al.: 
!     doi:10.1088/0965-0393/23/8/083501, 2015.
!
! Note: The following subroutines rely on other conversions:
! RODRIGUES_TO_ROT_MATS: First converts to axis-angle, then to rotation matrix.
! ROT_MATS_TO_AXIS_ANGLE: First converts to quaternions, then to axis-angle.
! ROT_MATS_TO_RODRIGUES: First converts to quaternions, then to Rodrigues.
!
! Contains subroutines:
! AXIS_ANGLE_TO_ROT_MATS: Axis-angle (degrees) to rotation matrices.
! EULER_BUNGE_TO_ROT_MATS: Euler-Bunge angles (degrees) to rotation matrices.
! EULER_KOCKS_TO_ROT_MATS: Euler-Kocks angles (degrees) to rotation matrices.
! RODRIGUES_TO_ROT_MATS: Rodrigues vectors to rotation matrices.
! ROT_MATS_TO_AXIS_ANGLE: Rotation matrices to axis-angle (degrees)
! ROT_MATS_TO_EULER_BUNGE: Rotation matrices to Euler-Bunge angles (degrees).
! ROT_MATS_TO_EULER_KOCKS: Rotation matrices to Euler-Kocks angles (degrees).
! ROT_MATS_TO_RODRIGUES: Rotation matrices to rodrigues vectors.
! ROT_MATS_TO_QUATERNIONS: Rotation matrices to quaternions.
! QUATERNIONS_TO_ROT_MATS: Quaternions to rotation matrices.
!
USE INTRINSICTYPESMODULE, ONLY: RK=>REAL_KIND
!
USE CONSTANTSMODULE, ONLY: RK_PI, RK_PI_OVER_180
USE DIMSMODULE
!
IMPLICIT NONE
!
! Public
!
PUBLIC
!
CONTAINS
    !
    SUBROUTINE AXIS_ANGLE_TO_ROT_MATS(N, M, AXIS, ANGLE, R)
        !
        ! Convert axis-angle pairs (degrees) to rotation matrices.
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments:
        ! N: Number of orientations per element (legacy)
        ! M: Number of elements
        ! AXIS: Array of normalized axes
        ! ANGLE: Array of angles (degrees)
        ! R: Array of rotation matrices
        !
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: M
        REAL(RK), INTENT(IN) :: AXIS(0:DIMS1, 0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(IN) :: ANGLE(0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: R(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
        !
        ! Locals:
        ! COM: Cosine of angle
        ! SOM: Sine of angle
        ! I, J: Looping indices
        !
        REAL(RK) :: ANGLERAD
        REAL(RK) :: COM
        REAL(RK) :: SOM
        INTEGER :: I, J
        !
        !-----------------------------------------------------------------------
        !
        R = 0.0_RK
        !
        DO I = 0, (N - 1)
            !
            DO J = 0, (M - 1)
                !
                ANGLERAD = 0.0_RK
                COM = 0.0_RK
                SOM = 0.0_RK
                !
                ANGLERAD = ANGLE(I, J) * RK_PI_OVER_180
                COM = COS(ANGLERAD)
                SOM = SIN(ANGLERAD)
                !
                ! Check if returned cos/sin values are near machine epsilon 
                IF (ABS(COM) .LT. VTINY) COM = 0.0_RK
                IF (ABS(SOM) .LT. VTINY) SOM = 0.0_RK
                !
                R(0, 0, I, J) = COM + (1.0_RK - COM) * (AXIS(0, I, J) ** 2)
                R(0, 1, I, J) = (1.0_RK - COM) * AXIS(0, I, J) * &
                    & AXIS(1, I, J) + SOM * AXIS(2, I, J)
                R(0, 2, I, J) = (1.0_RK - COM) * AXIS(0, I, J) * &
                    & AXIS(2, I, J) - SOM * AXIS(1, I, J)
                R(1, 0, I, J) = (1.0_RK - COM) * AXIS(0, I, J) * &
                    & AXIS(1, I, J) - SOM * AXIS(2, I, J)
                R(1, 1, I, J) = COM + (1.0_RK - COM) * (AXIS(1, I, J) ** 2)
                R(1, 2, I, J) = (1.0_RK - COM) * AXIS(1, I, J) * &
                    & AXIS(2, I, J) + SOM * AXIS(0, I, J)
                R(2, 0, I, J) = (1.0_RK - COM) * AXIS(0, I, J) * &
                    & AXIS(2, I, J) + SOM * AXIS(1, I, J)
                R(2, 1, I, J) = (1.0_RK - COM) * AXIS(1, I, J) * &
                    & AXIS(2, I, J) - SOM * AXIS(0, I, J)
                R(2, 2, I, J) = COM + (1.0_RK - COM) * (AXIS(2, I, J) ** 2)
                !
            ENDDO
            !
        ENDDO
        !
    END SUBROUTINE AXIS_ANGLE_TO_ROT_MATS
    !
    !===========================================================================
    !
    SUBROUTINE EULER_BUNGE_TO_ROT_MATS(N, M, PSI1, PHI, PSI2, R)
        !
        ! Convert Euler-Bunge angles (degrees) to rotation matrices.
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments:
        ! N: Number of orientations per element (legacy)
        ! M: Number of elements
        ! PSI1, PHI, PSI2: Euler-Bunge angles
        ! R: Array of rotation matrices
        !
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: M
        REAL(RK), INTENT(IN) :: PSI1(0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(IN) :: PHI(0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(IN) :: PSI2(0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: R(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
        !
        ! Locals:
        ! PSI1RAD, PHIRAD, PSI2RAD: Angles in radians
        ! SPS1: Sine of psi1
        ! CPS1: Cosine of psi1
        ! SPH: Sine of phi
        ! CPH: Cosine of phi
        ! SPS2: Sine of psi2
        ! CPS2: Cosine of psi2
        ! I/J: Looping indices
        !
        REAL(RK) :: PSI1RAD(0:(N - 1), 0:(M - 1))
        REAL(RK) :: PHIRAD(0:(N - 1), 0:(M - 1))
        REAL(RK) :: PSI2RAD(0:(N - 1), 0:(M - 1))
        REAL(RK) :: SPS1(0:(N - 1), 0:(M - 1))
        REAL(RK) :: CPS1(0:(N - 1), 0:(M - 1))
        REAL(RK) :: SPH(0:(N - 1), 0:(M - 1))
        REAL(RK) :: CPH(0:(N - 1), 0:(M - 1))
        REAL(RK) :: SPS2(0:(N - 1), 0:(M - 1))
        REAL(RK) :: CPS2(0:(N - 1), 0:(M - 1))
        INTEGER  :: I, J
        !
        !-----------------------------------------------------------------------
        !
        R = 0.0_RK
        !
        PSI1RAD = 0.0_RK
        PHIRAD = 0.0_RK
        PSI2RAD = 0.0_RK
        SPS1 = 0.0_RK
        CPS1 = 0.0_RK
        SPH = 0.0_RK
        CPH = 0.0_RK
        SPS2 = 0.0_RK
        CPS2 = 0.0_RK
        !
        ! Convert from degrees to radians
        !
        PSI1RAD = PSI1 * RK_PI_OVER_180
        PHIRAD = PHI * RK_PI_OVER_180
        PSI2RAD = PSI2 * RK_PI_OVER_180
        !
        SPS1 = SIN(PSI1RAD)
        CPS1 = COS(PSI1RAD)
        SPH = SIN(PHIRAD)
        CPH = COS(PHIRAD)
        SPS2 = SIN(PSI2RAD)
        CPS2 = COS(PSI2RAD)
        !
        ! Check if returned cos/sin values are near machine epsilon 
        DO I = 0, (N - 1)
            !
            DO J = 0, (M - 1)
                !
                IF (ABS(SPS1(I,J)) .LT. VTINY) SPS1(I,J) = 0.0_RK
                IF (ABS(CPS1(I,J)) .LT. VTINY) CPS1(I,J) = 0.0_RK
                IF (ABS(SPH(I,J))  .LT. VTINY) SPH(I,J)  = 0.0_RK
                IF (ABS(CPH(I,J))  .LT. VTINY) CPH(I,J)  = 0.0_RK
                IF (ABS(SPS2(I,J)) .LT. VTINY) SPS2(I,J) = 0.0_RK
                IF (ABS(CPS2(I,J)) .LT. VTINY) CPS2(I,J) = 0.0_RK
                !
            END DO
            !
        END DO
        !
        R(0, 0, :, :) = CPS1 * CPS2 - SPS1 * CPH * SPS2
        R(0, 1, :, :) = SPS1 * CPS2 + CPS1 * CPH * SPS2
        R(0, 2, :, :) = SPH * SPS2
        R(1, 0, :, :) = -CPS1 * SPS2 - SPS1 * CPH * CPS2
        R(1, 1, :, :) = -SPS1 * SPS2 + CPS1 * CPH * CPS2
        R(1, 2, :, :) = SPH * CPS2
        R(2, 0, :, :) = SPS1 * SPH
        R(2, 1, :, :) = -CPS1 * SPH
        R(2, 2, :, :) = CPH
        !
        !
        RETURN
        !
    END SUBROUTINE EULER_BUNGE_TO_ROT_MATS
    !
    !===========================================================================
    !
    SUBROUTINE EULER_KOCKS_TO_ROT_MATS(N, M, PSI, THE, PHI, R)
        !
        ! Convert Euler-Kocks angles (degrees) to rotation matrices.
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments:
        ! N: Number of orientations per element (legacy)
        ! M: Number of elements
        ! PSI, THE, PHI: Euler-Kocks angles (degrees)
        ! R: Array of rotation matrices
        !
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: M
        REAL(RK), INTENT(IN) :: PSI(0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(IN) :: THE(0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(IN) :: PHI(0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: R(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
        !
        ! Locals:
        ! PSIRAD, THERAD, PHIRAD: Angles in radians
        ! SPS: Sine of psi
        ! CPS: Cosine of psi
        ! STH: Sine of theta
        ! CTH: Cosine of theta
        ! SPH: Sine of phi
        ! CPH: Cosine of phi
        ! I/J: Looping indices
        !
        REAL(RK) :: PSIRAD(0:(N - 1), 0:(M - 1))
        REAL(RK) :: THERAD(0:(N - 1), 0:(M - 1))
        REAL(RK) :: PHIRAD(0:(N - 1), 0:(M - 1))
        REAL(RK) :: SPS(0:(N - 1), 0:(M - 1))
        REAL(RK) :: CPS(0:(N - 1), 0:(M - 1))
        REAL(RK) :: STH(0:(N - 1), 0:(M - 1))
        REAL(RK) :: CTH(0:(N - 1), 0:(M - 1))
        REAL(RK) :: SPH(0:(N - 1), 0:(M - 1))
        REAL(RK) :: CPH(0:(N - 1), 0:(M - 1))
        INTEGER  :: I, J
        !
        !-----------------------------------------------------------------------
        !
        R = 0.0_RK
        !
        PSIRAD = 0.0_RK
        THERAD = 0.0_RK
        PHIRAD = 0.0_RK
        SPS = 0.0_RK
        CPS = 0.0_RK
        STH = 0.0_RK
        CTH = 0.0_RK
        SPH = 0.0_RK
        CPH = 0.0_RK
        !
        ! Convert from degrees to radians
        !
        PSIRAD = PSI * RK_PI_OVER_180
        THERAD = THE * RK_PI_OVER_180
        PHIRAD = PHI * RK_PI_OVER_180
        !
        SPS = SIN(PSIRAD)
        CPS = COS(PSIRAD)
        STH = SIN(THERAD)
        CTH = COS(THERAD)
        SPH = SIN(PHIRAD)
        CPH = COS(PHIRAD)
        !
        ! Check if returned cos/sin values are near machine epsilon 
        DO I = 0, (N - 1)
            !
            DO J = 0, (M - 1)
                !
                IF (ABS(SPS(I,J)) .LT. VTINY) SPS(I,J) = 0.0_RK
                IF (ABS(CPS(I,J)) .LT. VTINY) CPS(I,J) = 0.0_RK
                IF (ABS(STH(I,J)) .LT. VTINY) STH(I,J) = 0.0_RK
                IF (ABS(CTH(I,J)) .LT. VTINY) CTH(I,J) = 0.0_RK
                IF (ABS(SPH(I,J)) .LT. VTINY) SPH(I,J) = 0.0_RK
                IF (ABS(CPH(I,J)) .LT. VTINY) CPH(I,J) = 0.0_RK
                !
            END DO
            !
        END DO
        !
        R(0, 0, :, :) = -SPS * SPH - CPS * CPH * CTH
        R(0, 1, :, :) = CPS * SPH - SPS * CPH * CTH
        R(0, 2, :, :) = CPH * STH
        R(1, 0, :, :) = CPH * SPS - SPH * CPS * CTH
        R(1, 1, :, :) = -CPS * CPH - SPS * SPH * CTH
        R(1, 2, :, :) = SPH * STH
        R(2, 0, :, :) = CPS * STH
        R(2, 1, :, :) = SPS * STH
        R(2, 2, :, :) = CTH
        !
    END SUBROUTINE EULER_KOCKS_TO_ROT_MATS
    !
    !===========================================================================
    !
    SUBROUTINE RODRIGUES_TO_ROT_MATS(N, M, RODS, R)
        !
        ! Convert Rodrigues vectors to rotation matrices.
        ! Note: first converts Rodrigues vector to axis-angle parameterization,
        ! then to rotation matrix.
        !
        !-----------------------------------------------------------------------
        !
        ! Arguments:
        ! N: Number of orientations per element (legacy)
        ! M: Number of elements
        ! RODS: Array of Rodrigues vectors
        ! R: Array of rotation matrices
        !
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: M
        REAL(RK), INTENT(IN) :: RODS(0:DIMS1, 0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: R(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
        !
        ! Locals:
        ! I, J: Looping indices
        ! ANGLE: Angle of rotation about axis
        ! AXIS: Axis of rotation
        ! COM: Cosine of angle
        ! SOM: Sine of angle
        !
        INTEGER :: I, J
        REAL(RK) :: ANGLE
        REAL(RK) :: AXIS(0:DIMS1)
        REAL(RK) :: COM
        REAL(RK) :: SOM
        !
        !-----------------------------------------------------------------------
        !
        R = 0.0_RK
        !
        DO I = 0, (N - 1)
            !
            DO J = 0, (M - 1)
                !
                ANGLE = 0.0_RK
                AXIS = 0.0_RK
                COM = 0.0_RK
                SOM = 0.0_RK
                !
                ! Convert to axis-angle representation
                !
                ANGLE = 2.0_RK * ATAN(NORM2(RODS(:, I, J)))
		        !
                IF (ANGLE .GT. VTINY) THEN
                    !
                    AXIS = RODS(:, I, J) / NORM2(RODS(:, I, J))
                    !
                ELSE
                    !
                    AXIS(0) = 1.0_RK
                    !
                END IF
                !
                COM = COS(ANGLE)
                SOM = SIN(ANGLE)
                !
                ! Check if returned cos/sin values are near machine epsilon 
                IF (ABS(COM) .LT. VTINY) COM = 0.0_RK
                IF (ABS(SOM) .LT. VTINY) SOM = 0.0_RK
                !
                ! Then from axis-angle to rotation matrix
                !
                R(0, 0, I, J) = COM + (1.0_RK - COM) * (AXIS(0) ** 2)
                R(0, 1, I, J) = (1.0_RK - COM) * AXIS(0) * AXIS(1) + &
                    & SOM * AXIS(2)
                R(0, 2, I, J) = (1.0_RK - COM) * AXIS(0) * AXIS(2) - &
                    & SOM * AXIS(1)
                R(1, 0, I, J) = (1.0_RK - COM) * AXIS(0) * AXIS(1) - &
                    & SOM * AXIS(2)
                R(1, 1, I, J) = COM + (1.0_RK - COM) * (AXIS(1) ** 2)
                R(1, 2, I, J) = (1.0_RK - COM) * AXIS(1) * AXIS(2) + &
                    & SOM * AXIS(0)
                R(2, 0, I, J) = (1.0_RK - COM) * AXIS(0) * AXIS(2) + &
                    & SOM * AXIS(1)
                R(2, 1, I, J) = (1.0_RK - COM) * AXIS(1) * AXIS(2) - &
                    & SOM * AXIS(0)
                R(2, 2, I, J) = COM + (1.0_RK - COM) * (AXIS(2) ** 2)
                !
            ENDDO
            !
        ENDDO
        !
    END SUBROUTINE RODRIGUES_TO_ROT_MATS
    !
    !===========================================================================
    !
    SUBROUTINE ROT_MATS_TO_AXIS_ANGLE(N, M, R, AXIS, ANGLE)
        !
        ! Convert rotation matrices to axis-angle (degrees)
        !
        !-----------------------------------------------------------------------
        !
        ! Arugments:
        ! N: Number of orientations per element (legacy)
        ! M: Number of elements
        ! R: Array of rotation matrices
        ! AXIS: Array of axes
        ! ANGLE: Array of rotations (degrees)
        !
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: M
        REAL(RK), INTENT(IN) :: R(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: AXIS(0:DIMS1, 0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: ANGLE(0:(N - 1), 0:(M - 1))
        !
        ! Locals:
        ! Q0, Q1, Q2, Q3: Individual components of quaternion
        ! S: Normalization factor for axis
        ! I, J: Looping indices
	    !
        REAL(RK) :: Q0
        REAL(RK) :: Q1
        REAL(RK) :: Q2
        REAL(RK) :: Q3
        REAL(RK) :: S
        INTEGER :: I, J
        !
        !-----------------------------------------------------------------------
        !
        ANGLE = 0.0_RK
        AXIS = 0.0_RK
        !
        DO I = 0, (N - 1)
            !
            DO J = 0, (M - 1)
                !
                ! Convert first to quaternion
                !
                Q0 = 0.0_RK
                Q1 = 0.0_RK
                Q2 = 0.0_RK
                Q3 = 0.0_RK
                !
                Q0 = 0.5_RK * SQRT(1 + R(0, 0, I, J) + R(1, 1, I, J) + &
                    & R(2, 2, I, J))
                Q1 = (-0.5_RK) * SQRT(1 + R(0, 0, I, J) - R(1, 1, I, J) - &
                    & R(2, 2, I, J))
                Q2 = (-0.5_RK) * SQRT(1 - R(0, 0, I, J) + R(1, 1, I, J) - &
                    & R(2, 2, I, J))
                Q3 = (-0.5_RK) * SQRT(1 - R(0, 0, I, J) - R(1, 1, I, J) + &
                    & R(2, 2, I, J))
                !
                IF (R(2, 1, I, J) .LT. R(1, 2, I, J)) THEN
                    !
                    Q1 = -Q1
                    !
                ENDIF
                !
                IF (R(0, 2, I, J) .LT. R(2, 0, I, J)) THEN
                    !
                    Q2 = -Q2
                    !
                ENDIF
                !
                IF (R(1, 0, I, J) .LT. R(0, 1, I, J)) THEN
                    !
                    Q3 = -Q3
                    !
                ENDIF
                !
                ! Then from quaternion to axis-angle
                !
                ANGLE(I, J) = 2.0_RK * ACOS(Q0)
                !
                IF (ANGLE(I, J) .NE. 0.0_RK) THEN
                    !
                    IF (Q0 .NE. 0.0_RK) THEN
                        !
                        S = SQRT((Q1 ** 2) + (Q2 ** 2) + (Q3 ** 2))
                        AXIS(0, I, J) = Q1 / S
                        AXIS(1, I, J) = Q2 / S
                        AXIS(2, I, J) = Q3 / S
                        !
                        IF (Q0 .LT. 0.0_RK) THEN
                            !
                            AXIS = -AXIS
                            !
                        ENDIF
                        !
                    ELSEIF (Q0 .EQ. 0.0_RK) THEN
                        !
                        ANGLE(I, J) = RK_PI
                        AXIS(0, I, J) = Q1
                        AXIS(1, I, J) = Q2
                        AXIS(2, I, J) = Q3
                        !
                    ENDIF
                    !
                ELSEIF (ANGLE(I, J) .EQ. 0.0_RK) THEN
                    !
                    AXIS(0, I, J) = 0.0_RK
                    AXIS(1, I, J) = 0.0_RK
                    AXIS(2, I, J) = 1.0_RK
                    !
                ENDIF
                !
                ! Convert from radians to degrees
                !
                ANGLE(I, J) = ANGLE(I, J) / RK_PI_OVER_180
                !
            ENDDO
            !
        ENDDO
        !
    END SUBROUTINE ROT_MATS_TO_AXIS_ANGLE
    !
    !===========================================================================
    !
    SUBROUTINE ROT_MATS_TO_EULER_BUNGE(N, M, R, PSI1, PHI, PSI2)
        !
        ! Convert rotation matrices to Euler-Bunge angles (degrees).
        !
        !-----------------------------------------------------------------------
        !
        ! Arugments:
        ! N: Number of orientations per element (legacy)
        ! M: Number of elements
        ! R: Array of rotation matrices
        ! PSI1, PHI, PSI2: Arrays of Euler angles
        !
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: M
        REAL(RK), INTENT(IN) :: R(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: PSI1(0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: PHI(0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: PSI2(0:(N - 1), 0:(M - 1))
        !
        ! Locals:
        ! XI: Internal scaling factor
        !
        REAL(RK) :: XI(0:(N - 1), 0:(M - 1))
        !
        !-----------------------------------------------------------------------
        !
        PSI1 = 0.0_RK
        PHI = 0.0_RK
        PSI2 = 0.0_RK
        XI = 0.0_RK
        !
        XI = 1.0_RK / SQRT(1 - (R(2, 2, :, :) ** 2))
        !
        WHERE (ABS(R(2, 2, :, :)) .NE. 1.0_RK)
            !
            PSI1 = ATAN2(R(2, 0, :, :) * XI, -R(2, 1, :, :) * XI)
            PHI = ACOS(R(2, 2, :, :))
            PSI2 = ATAN2(R(0, 2, :, :) * XI, R(1, 2, :, :) * XI)
            !
        ELSEWHERE
            !
            PSI1 = ATAN2(R(0, 1, :, :), R(0, 0, :, :))
            PHI = (RK_PI/2.0_RK) * (1.0_RK - R(2, 2, :, :))
            PSI2 = 0.0_RK
            !
        ENDWHERE
        !
        ! Convert from radians to degrees
        !
        PSI1 = PSI1 / RK_PI_OVER_180
        PHI = PHI / RK_PI_OVER_180
        PSI2 = PSI2 / RK_PI_OVER_180
        !
    END SUBROUTINE ROT_MATS_TO_EULER_BUNGE
    !
    !===========================================================================
    !
    SUBROUTINE ROT_MATS_TO_EULER_KOCKS(N, M, R, PSI, THE, PHI)
        !
        ! Convert rotation matrices to Euler-Kocks angles (degrees).
        !
        !-----------------------------------------------------------------------
        !
        ! Arugments:
        ! N: Number of orientations per element (legacy)
        ! M: Number of elements
        ! R: Array of rotation matrices
        ! PSI, THE, PHI: Arrays of Euler-Kocks angles
        !
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: M
        REAL(RK), INTENT(IN) :: R(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: PSI(0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: THE(0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: PHI(0:(N - 1), 0:(M - 1))
        !
        ! Locals:
        ! STH: Sine of theta
        !
        REAL(RK) :: STH(0:(N - 1), 0:(M - 1))
        !
        !-----------------------------------------------------------------------
        !
        PSI = 0.0_RK
        THE = 0.0_RK
        PHI = 0.0_RK
        !
        STH = 0.0_RK
        !
        THE = ACOS(R(2, 2, :, :))
        !
        WHERE (ABS(R(2, 2, :, :)) .NE. 1.)
            !
            STH = SIN(THE)
            PSI = ATAN2(R(2, 1, :, :) / STH, R(2, 0, :, :) / STH)
            PHI = ATAN2(R(1, 2, :, :) / STH, R(0, 2, :, :) / STH)
            !
        ELSEWHERE
            !
            PSI = 0.0_RK
            PHI = ATAN2(-R(1, 0, :, :), -R(0, 0, :, :))
            !
        ENDWHERE
        !
        ! Convert from radians to degrees
        !
        PSI = PSI / RK_PI_OVER_180
        THE = THE / RK_PI_OVER_180
        PHI = PHI / RK_PI_OVER_180
        !
    END SUBROUTINE ROT_MATS_TO_EULER_KOCKS
    !
    !===========================================================================
    !
    SUBROUTINE ROT_MATS_TO_RODRIGUES(N, M, R, RODS)
        !
        ! Convert rotation matrices to Rodrigues vectors.
        ! Note: first converts rotation matrices to quaternions and then to 
        ! Rodrigues vectors.
        !
        !-----------------------------------------------------------------------
        !
        ! Arugments:
        ! N: Number of orientations per element (legacy)
        ! M: Number of elements
        ! R: Array of rotation matrices
        ! RODS: Array of Rodrigues vectors
        !
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: M
        REAL(RK), INTENT(IN) :: R(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: RODS(0:DIMS1, 0:(N - 1), 0:(M - 1))
        !
        ! Locals:
        ! Q0, Q1, Q2, Q3: Individual components of quaternion
        ! S: Normalization factor for axis
        ! I, J: Looping indices
        !
        REAL(RK) :: Q0
        REAL(RK) :: Q1
        REAL(RK) :: Q2
        REAL(RK) :: Q3
        REAL(RK) :: S
        REAL(RK) :: T
        INTEGER :: I, J
        !
        !-----------------------------------------------------------------------
        !
        RODS = 0.0_RK
        !
        DO I = 0, (N - 1)
            !
            DO J = 0, (M - 1)
                !
                ! Convert first to quaternions
                !
                Q0 = 0.0_RK
                Q1 = 0.0_RK
                Q2 = 0.0_RK
                Q3 = 0.0_RK
                !
                ! Take ABS() of quantity within SQRT() to avoid NaNs
                Q0 = 0.5_RK * SQRT(ABS(1 + R(0, 0, I, J) + R(1, 1, I, J) + &
                    & R(2, 2, I, J)))
                Q1 = (-0.5_RK) * SQRT(ABS(1 + R(0, 0, I, J) - R(1, 1, I, J) - &
                    & R(2, 2, I, J)))        
                Q2 = (-0.5_RK) * SQRT(ABS(1 - R(0, 0, I, J) + R(1, 1, I, J) - &
                    & R(2, 2, I, J)))
                Q3 = (-0.5_RK) * SQRT(ABS(1 - R(0, 0, I, J) - R(1, 1, I, J) + &
                    & R(2, 2, I, J)))
                !
                IF (R(2, 1, I, J) .LT. R(1, 2, I, J)) THEN
                    !
                    Q1 = -Q1
                    !
                ENDIF
                !
                IF (R(0, 2, I, J) .LT. R(2, 0, I, J)) THEN
                    !
                    Q2 = -Q2
                    !
                ENDIF
                !
                IF (R(1, 0, I, J) .LT. R(0, 1, I, J)) THEN
                    !
                    Q3 = -Q3
                    !
                ENDIF
                !
                ! Then from quaternions to Rodrigues
                !
                S = SQRT((Q1 ** 2) + (Q2 ** 2) + (Q3 ** 2))
                T = TAN(ACOS(Q0))
                !
                IF (S .LE. VTINY) THEN
                    !
                    RODS(0, I, J) = 0.0_RK
                    RODS(1, I, J) = 0.0_RK
                    RODS(2, I, J) = 0.0_RK
                    !
                ELSE                  
                    RODS(0, I, J) = (Q1 / S) * T
                    RODS(1, I, J) = (Q2 / S) * T
                    RODS(2, I, J) = (Q3 / S) * T
                    !
                END IF
                !
            ENDDO
            !
        ENDDO
        !
    END SUBROUTINE ROT_MATS_TO_RODRIGUES
    !
    !===========================================================================
    !
    SUBROUTINE ROT_MATS_TO_QUATERNIONS(N, M, R, QUAT)
        !
        ! Convert rotation matrices to quaternions.
        !
        !-----------------------------------------------------------------------
        !
        ! Arugments:
        ! N: Number of orientations per element (legacy)
        ! M: Number of elements
        ! R: Array of rotation matrices
        ! QUAT: Array of quaternions
        !
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: M
        REAL(RK), INTENT(IN) :: R(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: QUAT(0:DIMS, 0:(N - 1), 0:(M - 1))
        !
        ! Locals:
        ! Q0, Q1, Q2, Q3: Individual components of quaternion
        ! I, J: Looping indices
	    !
        REAL(RK) :: Q0
        REAL(RK) :: Q1
        REAL(RK) :: Q2
        REAL(RK) :: Q3
        INTEGER :: I, J
        !
        !-----------------------------------------------------------------------
        !
        QUAT = 0.0_RK
        !
        DO I = 0, (N - 1)
            !
            DO J = 0, (M - 1)
                !
                Q0 = 0.0_RK
                Q1 = 0.0_RK
                Q2 = 0.0_RK
                Q3 = 0.0_RK
                !
                Q0 = 0.5_RK * SQRT(1 + R(0, 0, I, J) + R(1, 1, I, J) + &
                    & R(2, 2, I, J))
                Q1 = (-0.5_RK) * SQRT(1 + R(0, 0, I, J) - R(1, 1, I, J) - &
                    & R(2, 2, I, J))
                Q2 = (-0.5_RK) * SQRT(1 - R(0, 0, I, J) + R(1, 1, I, J) - &
                    & R(2, 2, I, J))
                Q3 = (-0.5_RK) * SQRT(1 - R(0, 0, I, J) - R(1, 1, I, J) + &
                    & R(2, 2, I, J))
                !
                IF (R(2, 1, I, J).LT.R(1, 2, I, J)) THEN
                    !
                    Q1 = -Q1
                    !
                ENDIF
                !
                IF (R(0, 2, I, J).LT.R(2, 0, I, J)) THEN
                    !
                    Q2 = -Q2
                    !
                ENDIF
                !
                IF (R(1, 0, I, J).LT.R(0, 1, I, J)) THEN
                    !
                    Q3 = -Q3
                    !
                ENDIF
                !
                QUAT(0, I, J) = Q0
                QUAT(1, I, J) = Q1
                QUAT(2, I, J) = Q2
                QUAT(3, I, J) = Q3
                QUAT(:, I, J) = QUAT(:, I, J) / NORM2(QUAT(:, I, J))
                !
            ENDDO
            !
        ENDDO
        !
    END SUBROUTINE ROT_MATS_TO_QUATERNIONS
    !
    !===========================================================================
    !
    SUBROUTINE QUATERNIONS_TO_ROT_MATS(N, M, QUAT, R)
        !
        ! Convert quaternions to rotation matrices.
        !
        !-----------------------------------------------------------------------
        !
        ! Arugments:
        ! N: Number of orientations per element (legacy)
        ! M: Number of elements
        ! QUAT: Array of quaternions
        ! R: Array of rotation matrices
        !
        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(IN) :: M
        REAL(RK), INTENT(IN) :: QUAT(0:DIMS, 0:(N - 1), 0:(M - 1))
        REAL(RK), INTENT(INOUT) :: R(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
        !
        ! Locals:
        ! QBAR: Quaternion magnitude
        ! I, J: Looping indices
        !
        REAL(RK) :: QBAR
        INTEGER :: I, J
        !
        !-----------------------------------------------------------------------
        !
        R = 0.0_RK
        !
        DO I = 0, (N - 1)
            !
            DO J = 0, (M - 1)
                !
                QBAR = 0.0_RK
                QBAR = (QUAT(0, I, J) ** 2) - ((QUAT(1, I, J)  ** 2) + &
                    & (QUAT(2, I, J)  ** 2) + (QUAT(3, I, J)  ** 2))
                !
                R(0, 0, I, J) = QBAR + (2 * (QUAT(1, I, J)**2))
                R(0, 1, I, J) = 2 * ((QUAT(1, I, J) * QUAT(2, I, J)) + &
                    & (QUAT(0, I, J) * QUAT(3, I, J)))
                R(0, 2, I, J) = 2 * ((QUAT(1, I, J) * QUAT(3, I, J)) - &
                    & (QUAT(0, I, J) * QUAT(2, I, J)))
                R(1, 0, I, J) = 2 * ((QUAT(1, I, J) * QUAT(2, I, J)) - &
                    & (QUAT(0, I, J) * QUAT(3, I, J)))
                R(1, 1, I, J) = QBAR + (2 * (QUAT(2, I, J)**2))
                R(1, 2, I, J) = 2 * ((QUAT(2, I, J) * QUAT(3, I, J)) + &
                    & (QUAT(0, I, J) * QUAT(1, I, J)))
                R(2, 0, I, J) = 2 * ((QUAT(1, I, J) * QUAT(3, I, J)) + &
                    & (QUAT(0, I, J) * QUAT(2, I, J)))
                R(2, 1, I, J) = 2 * ((QUAT(2, I, J) * QUAT(3, I, J)) - &
                    & (QUAT(0, I, J) * QUAT(1, I, J)))
                R(2, 2, I, J) = QBAR + (2 * (QUAT(3, I, J)**2))
                !
            ENDDO
            !
        ENDDO
        !
    END SUBROUTINE QUATERNIONS_TO_ROT_MATS
    !
END MODULE ORIENTATION_CONVERSION_MOD
