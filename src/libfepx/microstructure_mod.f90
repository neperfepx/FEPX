! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE MICROSTRUCTURE_MOD
!
! Module to handle reading of simulation.grains file, simulations.oris file,
!   assignment of orientations to elements, initialization of hardening values,
!
! Contains subroutines:
! ASSIGN_ANGLES_PHASES: Assign orientations to the elements and compute rotation
!   matrices
! SET_TMIN: Set minimum RSS to prevent overflow
! ALLOCATE_GAMMADOT: Allocate memory for shear rate
! ALLOCATE_GACCUMSHEAR: Allocate memory for accumulated shear
! ALLOCATE_CRSS_N: Allocate memory for CRSS of previous step
! CRSS_N_INITIALIZE: Initialize the crystal slip system hardness
! SET_FCC_BLOCK_MATRICES: Define hardening interaction matrix for FCC crystal
! SET_HCP_BLOCK_MATRICES: Define hardening interaction matrix for HCP crystal
! SET_BCC_BLOCK_MATRICES: Define hardening interaction matrix for BCC crystal
! PROCESS_MATERIAL_PARAMETERS:  Routine to process all crystal parameter info
!   from data read in from the configuration file.
!
! From libf95:
!
USE LIBF95, ONLY: RK => REAL_KIND, STANDARD_OUTPUT
!
! From libfepx:
!
USE CRYSTAL_TYPE_MOD
USE DIMENSIONS_MOD
USE MATRIX_OPERATIONS_MOD
USE ORIENTATION_CONVERSION_MOD
USE PARALLEL_MOD
USE READ_INPUT_MOD
USE UNITS_MOD
!
IMPLICIT NONE
!
PRIVATE
!
! Public subroutines
!
PUBLIC :: ALLOCATE_CRSS_N
PUBLIC :: ALLOCATE_GACCUMSHEAR
PUBLIC :: ALLOCATE_GAMMADOT
PUBLIC :: ASSIGN_ANGLES_PHASES
PUBLIC :: CRSS_N_INITIALIZE
!PUBLIC :: READ_GRAIN_ANGLES
!PUBLIC :: READ_GRAINS_PHASES
PUBLIC :: PROCESS_MATERIAL_PARAMETERS
PUBLIC :: SET_TMIN
PUBLIC :: CRYSTALTYPEGET, CRYSTALTYPECREATE ! From CRYSTAL_TYPE_MOD
!
! Public variables
!
PUBLIC :: ACCUMSHEAR
PUBLIC :: ACCUMSHEAR_CEN
PUBLIC :: CRYSTAL_PARM
PUBLIC :: CTYPE
PUBLIC :: ELAS_COEFFS
PUBLIC :: GACCUMSHEAR
PUBLIC :: GAMMADOT
PUBLIC :: NUMPHASES
PUBLIC :: N_VOCE
PUBLIC :: PHASE
PUBLIC :: T_MIN
PUBLIC :: FCC_H1, FCC_H2, FCC_H3, FCC_H4
PUBLIC :: BCC_H1, BCC_H2, BCC_H3, BCC_H4, BCC_H5, BCC_H6
PUBLIC :: HCP_H1, HCP_H2, HCP_H3, HCP_H4, HCP_H5, HCP_H6, HCP_H7, HCP_VERT
PUBLIC :: CYCLIC_PARM
!
! Orientation options (public):
!
PUBLIC :: ORIENTATION_OPTIONS

!
!  Crystal types:
!
TYPE(CRYSTALTYPETYPE), ALLOCATABLE :: CTYPE(:)
!
!  Crystal parameters:
!
REAL(RK), ALLOCATABLE :: CRYSTAL_PARM(:, :), N_VOCE(:), T_MIN(:)
REAL(RK), ALLOCATABLE :: ELAS_COEFFS(:, :)
!
!  Number of phases:
!
INTEGER :: NUMPHASES
!
INTEGER :: NUMEL, NUMSUB !, NUM_GRAINS
!INTEGER, ALLOCATABLE :: PHASE(:)
REAL(RK), ALLOCATABLE :: GAMMADOT(:, :, :), GACCUMSHEAR(:, :, :, :)
REAL(RK), ALLOCATABLE :: ACCUMSHEAR(:, :, :), ACCUMSHEAR_CEN(:, :, :), &
    & CYCLIC_PARM(:, :)
!
REAL(RK) :: FCC_H1(0:2, 0:2), FCC_H2(0:2, 0:2), FCC_H3(0:2, 0:2), &
    & FCC_H4(0:2, 0:2)
REAL(RK) :: BCC_H1(0:1, 0:1), BCC_H2(0:1, 0:1), BCC_H3(0:1, 0:1), &
    & BCC_H4(0:1, 0:1), BCC_H5(0:1, 0:1), BCC_H6(0:1, 0:1)
REAL(RK) :: HCP_H1(0:2, 0:2), HCP_H2(0:1, 0:1), HCP_H3(0:1, 0:1), &
    & HCP_H4(0:1, 0:1), HCP_H5(0:1, 0:1), HCP_H6(0:1, 0:1), HCP_H7(0:1, 0:1), &
    & HCP_VERT(0:2)
!
! Locals
!
!INTEGER, ALLOCATABLE :: UESUB(:)
!REAL(RK), ALLOCATABLE :: GRAIN_ORIENTATION(:,:)
CHARACTER(LEN = 128) :: MESSAGE
!
CONTAINS
    !
    SUBROUTINE ASSIGN_ANGLES_PHASES(C0_ANGS, CRSS_N, RSTAR_N, WTS)
    !
    ! Assign orientations to the elements and compute rotation matrices
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Fortran unit number
    ! C0_ANGS: Array of initial rotation tensors
    ! CRSS_N: crystal slip system hardnesses
    ! RSTAR_N: Array of rotation matrices representing change in
    !   orientation, initialized to identity in this routine
    ! WTS: weights
    !
    !INTEGER :: IO
    REAL(RK) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    INTEGER  :: MY_PHASE(0:NGRAIN1, 0:(EL_SUP1-EL_SUB1))
    INTEGER  :: MY_UESUB(EL_SUB1:EL_SUP1)
    !
    REAL(RK) :: ANGLE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: AXIS(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: PSI1(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: PHI(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: PSI2(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: PSI(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: THE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: RODS(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: QUAT(0:DIMS, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    !
    INTEGER  :: I, J, NUM, I1, J1, IPHASE, NUMIND, IGRAIN, M
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(0,:) = PHASE(EL_SUB1:EL_SUP1)
    MY_UESUB = UESUB(EL_SUB1:EL_SUP1)
    !
    M = EL_SUP1 - EL_SUB1 + 1
    !
    ! Weights
    !
    WTS = 1.0D0 / NGRAIN
    !
    ! Check whether or not element orientations are to be used. Assign
    ! accordingly either per-grain or per-element.
    IF (ELEMENT_ORIS .EQV. .FALSE.) THEN
        !
        ! Assign the orientations to each element FROM GRAINS
        !
        DO I = 0, (NUM_GRAINS - 1)
            !
            DO IGRAIN = 0, NGRAIN1
                !
                IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
                    & .EQ. 'axis-angle') THEN
                    !
                    WHERE (I+1 .EQ. MY_UESUB) ! I+1 as UESUB is 1-indexed
                        !
                        AXIS(0, IGRAIN, :) = GRAIN_ORIENTATION(0, I)
                        AXIS(1, IGRAIN, :) = GRAIN_ORIENTATION(1, I)
                        AXIS(2, IGRAIN, :) = GRAIN_ORIENTATION(2, I)
                        ANGLE(IGRAIN, :) = GRAIN_ORIENTATION(3, I)
                        !
                    ENDWHERE
                    !
                ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
                    & .EQ. 'euler-bunge') THEN
                    !
                    WHERE (I+1 .EQ. MY_UESUB)
                        !
                        PSI1(IGRAIN, :) = GRAIN_ORIENTATION(0, I)
                        PHI(IGRAIN, :) = GRAIN_ORIENTATION(1, I)
                        PSI2(IGRAIN, :) = GRAIN_ORIENTATION(2, I)
                        !
                    ENDWHERE
                    !
                ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
                    & .EQ. 'euler-kocks') THEN
                    !
                    WHERE (I+1 .EQ. MY_UESUB)
                        !
                        PSI(IGRAIN, :) = GRAIN_ORIENTATION(0, I)
                        THE(IGRAIN, :) = GRAIN_ORIENTATION(1, I)
                        PHI(IGRAIN, :) = GRAIN_ORIENTATION(2, I)
                        !
                    ENDWHERE
                    !
                ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
                    & .EQ. 'rodrigues') THEN
                    !
                    WHERE (I+1 .EQ. MY_UESUB)
                        !
                        RODS(0, IGRAIN, :) = GRAIN_ORIENTATION(0, I)
                        RODS(1, IGRAIN, :) = GRAIN_ORIENTATION(1, I)
                        RODS(2, IGRAIN, :) = GRAIN_ORIENTATION(2, I)
                        !
                    ENDWHERE
                    !
                ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
                    & .EQ. 'quaternion') THEN
                    !
                    WHERE (I+1 .EQ. MY_UESUB)
                        !
                        QUAT(0, IGRAIN, :) = GRAIN_ORIENTATION(0, I)
                        QUAT(1, IGRAIN, :) = GRAIN_ORIENTATION(1, I)
                        QUAT(2, IGRAIN, :) = GRAIN_ORIENTATION(2, I)
                        QUAT(3, IGRAIN, :) = GRAIN_ORIENTATION(3, I)
                        !
                    ENDWHERE
                    !
                ENDIF
                !
            ENDDO
            !
        ENDDO
        !
    ELSE IF (ELEMENT_ORIS .EQV. .TRUE.) THEN
        !
        ! Assign the orientations to each element DIRECTLY
        !
        DO I = 0, (NUMELM - 1)
            !
            DO IGRAIN = 0, NGRAIN1
                !
                IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
                    & .EQ. 'axis-angle') THEN
                    !
                    IF ((I .GE. EL_SUB1) .AND. (I .LE. EL_SUP1)) THEN
                        !
                        AXIS(0, IGRAIN, I) = GRAIN_ORIENTATION(0, I)
                        AXIS(1, IGRAIN, I) = GRAIN_ORIENTATION(1, I)
                        AXIS(2, IGRAIN, I) = GRAIN_ORIENTATION(2, I)
                        ANGLE(IGRAIN, I) = GRAIN_ORIENTATION(3, I)
                        !
                    END IF
                    !
                ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
                    & .EQ. 'euler-bunge') THEN
                    !
                    IF ((I .GE. EL_SUB1) .AND. (I .LE. EL_SUP1)) THEN
                        !
                        PSI1(IGRAIN, I) = GRAIN_ORIENTATION(0, I)
                        PHI(IGRAIN, I) = GRAIN_ORIENTATION(1, I)
                        PSI2(IGRAIN, I) = GRAIN_ORIENTATION(2, I)
                        !
                    END IF
                    !
                ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
                    & .EQ. 'euler-kocks') THEN
                    !
                    IF ((I .GE. EL_SUB1) .AND. (I .LE. EL_SUP1)) THEN
                        !
                        PSI(IGRAIN, I) = GRAIN_ORIENTATION(0, I)
                        THE(IGRAIN, I) = GRAIN_ORIENTATION(1, I)
                        PHI(IGRAIN, I) = GRAIN_ORIENTATION(2, I)
                        !
                    END IF
                    !
                ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
                    & .EQ. 'rodrigues') THEN
                    !
                    IF ((I .GE. EL_SUB1) .AND. (I .LE. EL_SUP1)) THEN
                        !
                        RODS(0, IGRAIN, I) = GRAIN_ORIENTATION(0, I)
                        RODS(1, IGRAIN, I) = GRAIN_ORIENTATION(1, I)
                        RODS(2, IGRAIN, I) = GRAIN_ORIENTATION(2, I)
                        !
                    END IF
                    !
                ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
                    & .EQ. 'quaternion') THEN
                    !
                    IF ((I .GE. EL_SUB1) .AND. (I .LE. EL_SUP1)) THEN
                        !
                        QUAT(0, IGRAIN, I) = GRAIN_ORIENTATION(0, I)
                        QUAT(1, IGRAIN, I) = GRAIN_ORIENTATION(1, I)
                        QUAT(2, IGRAIN, I) = GRAIN_ORIENTATION(2, I)
                        QUAT(3, IGRAIN, I) = GRAIN_ORIENTATION(3, I)
                        !
                    END IF
                    !
                ENDIF
                !
            ENDDO
            !
        ENDDO
        !
    END IF
    !
    ! Find initial orientation matrix, c0_angs
    !
    ! Determine parameterization from orientation options
    !
    IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION .EQ. &
        & 'axis-angle') THEN
        !
        CALL AXIS_ANGLE_TO_ROT_MATS(NGRAIN, M, AXIS, ANGLE, c0_angs)
        !
    ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION .EQ. &
        & 'euler-bunge') THEN
        !
        CALL EULER_BUNGE_TO_ROT_MATS(NGRAIN, M, PSI1, PHI, PSI2, &
            & c0_angs)
        !
    ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION .EQ. &
        & 'euler-kocks') THEN
        !
        CALL EULER_KOCKS_TO_ROT_MATS(NGRAIN, M, PSI, THE, PHI, &
            & c0_angs)
        !
    ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION .EQ. &
        & 'rodrigues') THEN
        !
        CALL RODRIGUES_TO_ROT_MATS(NGRAIN, M, RODS, c0_angs)
        !
    ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION .EQ. &
        & 'quaternion') THEN
        !
        CALL QUATERNIONS_TO_ROT_MATS(NGRAIN, M, QUAT, c0_angs)
        !
    ENDIF
    !
    ! Determine passive (C2S) or active (S2C) from orientation options
    !
    IF (ORIENTATION_OPTIONS%ORIENTATION_CONVENTION .EQ. &
        & 'passive') THEN
        !
        ! Don't do anything - FEPX assumes passive convention!!!
        !
    ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_CONVENTION .EQ. &
        & 'active') THEN
        !
        DO J = EL_SUB1, EL_SUP1
            !
            DO IGRAIN = 0, NGRAIN1
                !
                C0_ANGS(:, :, IGRAIN, J) = TRANSPOSE(C0_ANGS(:, :, &
                    & IGRAIN, J))
                !
            ENDDO
            !
        ENDDO
        !
    ENDIF
    !
    ! Initialize hardnesses
    !
    CALL CRSS_N_INITIALIZE(CRSS_N, MY_PHASE)
    !
    ! Initialize RSTAR_N (identity matrix)
    !
    RSTAR_N = 0.0D0
    RSTAR_N(0, 0, :, :) = 1.0D0
    RSTAR_N(1, 1, :, :) = 1.0D0
    RSTAR_N(2, 2, :, :) = 1.0D0
    !
    RETURN
    !
    END SUBROUTINE ASSIGN_ANGLES_PHASES
    !
    !===========================================================================
    !
    SUBROUTINE SET_TMIN(T_MIN, XM)
    !
    ! Set minimum tau (resolved shear stress) from rate dependence, to prevent
    !   overflow.
    !
    !----------------------------------------------------------------------
    !
    ! Arguments:
    ! T_MIN: Minimum resolved shear stress
    ! XM: ??
    !
    REAL(RK), INTENT(OUT) :: T_MIN
    REAL(RK), INTENT(IN)  :: XM
    !
    ! Locals:
    ! XNN: ??
    !
    REAL(RK) :: XNN
    !
    !----------------------------------------------------------------------
    !
    XNN = 1.0 / XM - 1.0
    !
    IF (XNN .GT. 1.0) THEN
        !
        T_MIN = DEXP((-321.D0 / XNN + 2.D0) * DLOG(10.D0))
        !
    ELSE
        !
        T_MIN = DEXP(-50.D0)
        !
    ENDIF
    !
    END SUBROUTINE SET_TMIN
    !
    !===========================================================================
    !
    SUBROUTINE ALLOCATE_GAMMADOT(NGRAIN1)
    !
    ! Allocates memory for shear rate
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NGRAIN1: Number of grains per element (minus 1, obsolete, always 0)
    !
    INTEGER, INTENT(IN) :: NGRAIN1
    !
    !---------------------------------------------------------------------------
    !
    ALLOCATE(GAMMADOT(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1))
    !
    ! Initialize
    !
    GAMMADOT = 0.0D0
    !
    END SUBROUTINE ALLOCATE_GAMMADOT
    !
    !===========================================================================
    !
    SUBROUTINE ALLOCATE_GACCUMSHEAR(NGRAIN1, NQPT1)
    !
    ! Allocate memory for accumulated shear
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NGRAIN1: Number of grains per element (minus 1, obsolete, always 0)
    ! NQPT1: Number of quadrature points per element (minus 1)
    !
    INTEGER, INTENT(IN) :: NGRAIN1
    INTEGER, INTENT(IN) :: NQPT1
    !
    !---------------------------------------------------------------------------
    !
    ALLOCATE(GACCUMSHEAR(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1))
    ALLOCATE(ACCUMSHEAR(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1))
    ALLOCATE(ACCUMSHEAR_CEN(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1))
    !
    ! Initialize
    !
    GACCUMSHEAR = 0.0D0
    ACCUMSHEAR = 0.0D0
    ACCUMSHEAR_CEN = 0.0D0
    !
    END SUBROUTINE ALLOCATE_GACCUMSHEAR
    !
    !===========================================================================
    !
    SUBROUTINE ALLOCATE_CRSS_N(CRSS_N, NGRAIN1)
    !
    ! Allocate memory for CRSS of previous step
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments
    ! CRSS_N: CRSS of previous step
    ! NGRAIN1: Number of grains per element (minus 1, obsolete, always 0)
    !
    REAL(RK), INTENT(INOUT), ALLOCATABLE :: CRSS_N(:,:,:)
    INTEGER, INTENT(IN) :: NGRAIN1
    !
    ALLOCATE(CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1))
    !
    CRSS_N = 0.0D0
    !
    END SUBROUTINE ALLOCATE_CRSS_N
    !
    !===========================================================================
    !
    SUBROUTINE CRSS_N_INITIALIZE(CRSS_N,MY_PHASE)
    !
    ! Initialize the crystal slip system hardness
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! CRSS_N: Crystal slip system hardnesses
    ! MY_PHASE: Crystal phase currently being initialized on this processor
    !
    REAL(RK), INTENT(OUT) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER,  INTENT(IN)  :: MY_PHASE(0:(EL_SUP1-EL_SUB1))
    !
    ! Locals:
    ! INDICES: (??)
    ! ISLIP: Loop index over number of slip systems
    ! N_SLIP: Number of slip systems for the crystal phase and type
    ! IPHASE: Loop index over number of crystal phases
    ! NUMIND: Number of indices returned by FIND_INDICES (??)
    !
    INTEGER, POINTER :: INDICES(:) => NULL()
    INTEGER :: ISLIP, N_SLIP, IPHASE, NUMIND
    !
    !---------------------------------------------------------------------------
    !
    ! Initialize
    CRSS_N= 0.0D0
    !
    ! Assign the crss based on the crystal type
    DO IPHASE = 1, NUMPHASES
        !
        ! Find the numbers of slip systems
        CALL CRYSTALTYPEGET(CTYPE(IPHASE))
        N_SLIP=CTYPE(IPHASE)%NUMSLIP
        !
        ! Finds the indices corresponding to the current phase the loop is on
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
        !
        ! Initialize the crystal slip strengths corresponding to the slip number
        DO ISLIP = 0, N_SLIP-1
            !
            CRSS_N(ISLIP,:,EL_SUB1 + INDICES)=CRYSTAL_PARM(9,IPHASE)
            !        
        ENDDO
        !
    ENDDO
    !
    END SUBROUTINE CRSS_N_INITIALIZE
    !
    !===========================================================================
    !
    SUBROUTINE SET_FCC_BLOCK_MATRICES(H1,H2,H3,H4,DIAG)
    !
    ! Define hardening interaction matrix for FCC crystal type
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! H1, H2, H3, H4: Hardening interaction values
    ! DIAG: Matrix diagonal value
    !
    REAL(RK), INTENT(IN) :: H1,H2,H3,H4
    REAL(RK), INTENT(IN) :: DIAG
    !
    !---------------------------------------------------------------------------
    !
    FCC_H1 = RESHAPE((/ DIAG, H1, H1, H1, DIAG, H1, H1, H1, DIAG /), (/3, 3/))
    FCC_H2 = RESHAPE((/ DIAG, H2, H2, H2, DIAG, H2, H2, H2, DIAG /), (/3, 3/))
    FCC_H3 = RESHAPE((/ DIAG, H3, H3, H3, DIAG, H3, H3, H3, DIAG /), (/3, 3/))
    FCC_H4 = RESHAPE((/ DIAG, H4, H4, H4, DIAG, H4, H4, H4, DIAG /), (/3, 3/))
    !
    END SUBROUTINE SET_FCC_BLOCK_MATRICES
    !
    !===========================================================================
    !
    SUBROUTINE SET_HCP_BLOCK_MATRICES(H1,H2,H3,H4,H5,H6,H7,DIAG)
    !
    ! Define hardening interaction matrix for HCP crystal type
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! H1, H2, H3, H4, H5, H6, H7: Hardening interaction values
    ! DIAG: Matrix diagonal value
    !
    REAL(RK), INTENT(IN) :: H1,H2,H3,H4,H5,H6,H7
    REAL(RK), INTENT(IN) :: DIAG
    !
    !---------------------------------------------------------------------------
    !
    HCP_H1 = RESHAPE((/ DIAG, H1, H1, H1, DIAG, H1, H1, H1, DIAG /), (/3, 3/))
    HCP_VERT = DIAG
    HCP_H2 = RESHAPE((/ DIAG, H2, H2, DIAG /), (/2,2/))
    HCP_H3 = RESHAPE((/ DIAG, H3, H3, DIAG /), (/2,2/))
    HCP_H4 = RESHAPE((/ DIAG, H4, H4, DIAG /), (/2,2/))
    HCP_H5 = RESHAPE((/ DIAG, H5, H5, DIAG /), (/2,2/))
    HCP_H6 = RESHAPE((/ DIAG, H6, H6, DIAG /), (/2,2/))
    HCP_H7 = RESHAPE((/ DIAG, H7, H7, DIAG /), (/2,2/))
    !
    END SUBROUTINE SET_HCP_BLOCK_MATRICES
    !
    !===========================================================================
    !
    SUBROUTINE SET_BCC_BLOCK_MATRICES (H1,H2,H3,H4,H5,H6,DIAG)
    !
    ! Define hardening interaction matrix for BCC crystal type
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! H1, H2, H3, H4, H5, H6: Hardening interaction values
    ! DIAG: Matrix diagonal value
    !
    REAL(RK), INTENT(IN) :: H1,H2,H3,H4,H5,H6
    REAL(RK), INTENT(IN) :: DIAG
    !
    !---------------------------------------------------------------------------
    !
    BCC_H1 = RESHAPE((/ DIAG, H1, H1, DIAG /), (/2,2/))
    BCC_H2 = RESHAPE((/ DIAG, H2, H2, DIAG /), (/2,2/))
    BCC_H3 = RESHAPE((/ DIAG, H3, H3, DIAG /), (/2,2/))
    BCC_H4 = RESHAPE((/ DIAG, H4, H4, DIAG /), (/2,2/))
    BCC_H5 = RESHAPE((/ DIAG, H5, H5, DIAG /), (/2,2/))
    BCC_H6 = RESHAPE((/ DIAG, H6, H6, DIAG /), (/2,2/))
    !
    END SUBROUTINE SET_BCC_BLOCK_MATRICES
    !
    !===========================================================================
    !
    SUBROUTINE PROCESS_MATERIAL_PARAMETERS(KELAS, KEINV)
    !
    ! Routine to process all crystal parameter info from data read in from the
    !   configuration file.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! KELAS: Anisotropic elastic stiffness tensor - crystal type dependent.
    ! KEINV: Inverse of anisotropic elastic stiffness tensor.
    !
    REAL(RK), POINTER, INTENT(INOUT) :: KELAS(:, :), KEINV(:, :)
    !
    ! Locals:
    ! I/IPHASE: Generic loop index values
    ! NUMSLIP: Number of restricted slip systems for crystal phase.
    ! MAXVERT: Maximum number of single crystal yield surface vertices.
    ! HCP_RATIOS: Temporary array for extracting HCP specific data.
    ! C11,C12,C13,C33,C44: Crystal elastic constants for phase.
    ! DIAG,H1-H7: Slip family scaling terms for hardening interaction matrix.
    !
    INTEGER           :: I, IPHASE, NUMSLIP, MAXVERT
    REAL(RK)          :: HCP_RATIOS(0:2)
    REAL(RK)          :: C11, C12, C13, C33, C44
    REAL(RK)          :: DIAG, H1, H2, H3, H4, H5, H6, H7
    !
    ! Notes:
    ! The elastic constants for each phase (Cij) are written in the Strength of
    ! Materials (SOM) convention and C44 is multiplied by 2 after being read in
    ! from the *.config file. The scaling terms for the hardening interaction
    ! matrix (DIAG, H1-H7) follow the implementation of anisotropic (latent)
    ! hardening into FEpX by Carson et. al. in:
    !
    !   Characterizing heterogeneous intragranular deformations in poly-
    !   -crystalline solids using diffraction-based and mechanics-based metrics
    !   https://doi.org/10.1088/1361-651X/aa6dc5
    !
    ! Similarly, the cyclic hardening parameters from CYCLIC_PARM follow the 
    ! implementation of Turkmen's model into FEpX by Turkmen et. al in:
    !
    !   A formulation for the pseudo‐saturation behavior observed during 
    !   variable amplitude multiaxial cyclic plasticity
    !   https://doi.org/10.1063/1.1766796
    !
    !---------------------------------------------------------------------------
    !
    ! Initalize variables.
    NUMSLIP    = 0
    MAXVERT    = 0
    NUMPHASES  = CRYS_OPTIONS%NUMBER_OF_PHASES
    !
    ! Print initial message to console.
    IF (MYID .EQ. 0) THEN
        !
        WRITE(DFLT_U, '(A)') 'Info   :   - Material parameters:'
        WRITE(DFLT_U, '(A,I0)') 'Info   :     > number_of_phases ', NUMPHASES
        !
    END IF
    !
    ! Allocate crystal parameter arrays.
    ALLOCATE(CRYSTAL_PARM(0:N_PARM-1, 1:NUMPHASES))
    ALLOCATE(CTYPE(1:NUMPHASES))
    ALLOCATE(N_VOCE(1:NUMPHASES))
    ALLOCATE(KELAS(0:TVEC1, 1:NUMPHASES))
    ALLOCATE(KEINV(0:TVEC1, 1:NUMPHASES))
    ALLOCATE(ELAS_COEFFS(0:3, 1:NUMPHASES))
    ALLOCATE(T_MIN(1:NUMPHASES))
    ALLOCATE(CYCLIC_PARM(0:1, 1:NUMPHASES))
    !
    !---------------------------------------------------------------------------
    !
    ! Loop over each phase and push the material data into the proper arrays.
    !
    DO IPHASE = 1, NUMPHASES
        !
        ! Assign standard crystal parameters for each phase.
        !
        CRYSTAL_PARM(0, IPHASE) = CRYS_OPTIONS%M(IPHASE)
        CRYSTAL_PARM(1, IPHASE) = CRYS_OPTIONS%GAMMADOT_0(IPHASE)
        CRYSTAL_PARM(2, IPHASE) = CRYS_OPTIONS%H_0(IPHASE)
        CRYSTAL_PARM(3, IPHASE) = CRYS_OPTIONS%G_0(IPHASE)
        CRYSTAL_PARM(4, IPHASE) = CRYS_OPTIONS%G_S0(IPHASE)
        CRYSTAL_PARM(5, IPHASE) = CRYS_OPTIONS%M_PRIME(IPHASE)
        CRYSTAL_PARM(6, IPHASE) = CRYS_OPTIONS%GAMMADOT_S0(IPHASE)
        CRYSTAL_PARM(7, IPHASE) = HUGE(0.0D0) ! Legacy fix to remove shear mod.
        N_VOCE(IPHASE)          = CRYS_OPTIONS%N(IPHASE)
        !
        ! Compute bulk modulus for the phase from elastic stiffness tensor.
        !
        SELECT CASE (CRYS_OPTIONS%CRYSTAL_TYPE(IPHASE))
            !
            CASE (3) ! Only need to handle HCP differently due to decoupling.
                !
                ! Assign stiffness tensor values for the phase.
                !
                C11 = CRYS_OPTIONS%C11(IPHASE)
                C12 = CRYS_OPTIONS%C12(IPHASE)
                C13 = CRYS_OPTIONS%C13(IPHASE)
                C44 = CRYS_OPTIONS%C44(IPHASE)
                !
                ! Compute C33 from decoupled hexagonal values.
                !
                C33 = C11 + C12 - C13
                !
                ELAS_COEFFS(:, IPHASE) = (/ C11, C12, C13, C44 /)
                KELAS(:, IPHASE) = (/ (C11 - C12), &
                             & 1/3.*(C11 + C12 - 4.*C13 + 2.*C33), &
                             & 2.*C44, 2.*C44, (C11 - C12) /)
                !
                CRYSTAL_PARM(8, IPHASE) = (2.*C11 + 2.*C12 + 4.*C13 + C33) / 9.
                !
                NUMSLIP = 18
                !
            CASE DEFAULT ! Handle the FCC and BCC construction the same way.
                ! 
                ! Assign stiffness tensor values for the phase.
                !
                C11 = CRYS_OPTIONS%C11(IPHASE)
                C12 = CRYS_OPTIONS%C12(IPHASE)
                C13 = CRYS_OPTIONS%C12(IPHASE)
                C44 = CRYS_OPTIONS%C44(IPHASE)
                !
                ELAS_COEFFS(:, IPHASE) = (/ C11, C12, C12, C44 /)
                KELAS(:, IPHASE) = (/ C11-C12, C11-C12, 2*C44, 2*C44, 2*C44 /)
                !
                CRYSTAL_PARM(8, IPHASE) = (C11 + 2 * C12) / 3.
                !
                IF (NUMSLIP .LE. 1) THEN ! Legacy code addition. Needed?
                    !
                    NUMSLIP = 12
                    !                
                END IF
                !
        END SELECT
        !
        ! Construct the inverse of KELAS.
        !
        KEINV(:,IPHASE) = 1.0D0 / KELAS(:,IPHASE)
        !
        ! Finish assigning standard crystal parameters for each phase.
        !
        CRYSTAL_PARM(9, IPHASE)  = CRYSTAL_PARM(3, IPHASE)
        CRYSTAL_PARM(10, IPHASE) = 0.0D0
        !
        ! Check if any values in this phase are invalid and if so quit.
        ! DEBUG: This check needs to be reconciled better - JC
        !
        IF ((ANY(CRYSTAL_PARM(:, IPHASE) .LT. 0.0D0)) .AND. &
            & (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .FALSE.)) THEN
            !
            CALL PAR_QUIT('Error  :     > Number of phases&
                    & does not match input data.')
        END IF
        !
        ! If current phase is HCP, handle specific phase data.
        !
        IF (CRYS_OPTIONS%CRYSTAL_TYPE(IPHASE) .EQ. 3) THEN
            HCP_RATIOS = 1.0D0
            !
            HCP_RATIOS(0) = CRYS_OPTIONS%C_OVER_A(IPHASE)
            HCP_RATIOS(1) = CRYS_OPTIONS%PYRAMIDAL_TO_BASAL(IPHASE)
            HCP_RATIOS(2) = CRYS_OPTIONS%PRISMATIC_TO_BASAL(IPHASE)
            !
            CRYSTAL_PARM(11, IPHASE) = HCP_RATIOS(0)
        END IF
        !
        ! Make the CrystalType object (CTYPE).
        !
        SELECT CASE (CRYS_OPTIONS%CRYSTAL_TYPE(IPHASE))
            !
            CASE (1) ! FCC
                !
                CTYPE(IPHASE) = CRYSTALTYPECREATE(CLASS_FCC,&
                             ! & VERTEX_FILE='vertices-fcc',&
                              & DECOMP=DECOMP_FEMEVPS)
                !
            CASE (2) ! BCC
                !
                CTYPE(IPHASE) = CRYSTALTYPECREATE(CLASS_BCC,&
                             ! & VERTEX_FILE='vertices-fcc',&
                              & DECOMP=DECOMP_FEMEVPS)
                !
            CASE (3) ! HCP
                !
                CTYPE(IPHASE) = CRYSTALTYPECREATE(CLASS_HCP,&
                              & C_OVER_A=HCP_RATIOS(0),&
                              & HRATIO_HCP=HCP_RATIOS(1),&
                              & HRATIO_HCP_PRISM=HCP_RATIOS(2),&
                             ! & VERTEX_FILE='vertices-hcp',&
                              & DECOMP=DECOMP_FEMEVPS)
                !
        END SELECT
        !
        ! Get maximum number of SCYS vertices.
        !
        MAXVERT = MAX(MAXVERT,CTYPE(IPHASE)%NUMVERTICES)
        !
        ! Set minimum value for tau to prevent overflow.
        !
        CALL SET_TMIN(T_MIN(IPHASE), CRYSTAL_PARM(0, IPHASE))
        !
        ! Handle latent hardening parameters - if available.
        !
        IF (OPTIONS%HARD_TYPE .EQ. 'latent') THEN
            !
            ! Assign block matrix values for the phase.
            DIAG = CRYS_OPTIONS%LATENT_PARAMETERS(IPHASE,1)
            H1   = CRYS_OPTIONS%LATENT_PARAMETERS(IPHASE,2)
            H2   = CRYS_OPTIONS%LATENT_PARAMETERS(IPHASE,3)
            H3   = CRYS_OPTIONS%LATENT_PARAMETERS(IPHASE,4)
            H4   = CRYS_OPTIONS%LATENT_PARAMETERS(IPHASE,5)
            H5   = CRYS_OPTIONS%LATENT_PARAMETERS(IPHASE,6)
            H6   = CRYS_OPTIONS%LATENT_PARAMETERS(IPHASE,7)
            H7   = CRYS_OPTIONS%LATENT_PARAMETERS(IPHASE,8)
            !
            SELECT CASE (CRYS_OPTIONS%CRYSTAL_TYPE(IPHASE))  
                !
                CASE (1) ! FCC          
                    ! 
                    CALL SET_FCC_BLOCK_MATRICES(H1, H2, H3, H4, DIAG)
                    !
                CASE (2) ! BCC
                    !
                    CALL SET_BCC_BLOCK_MATRICES(H1, H2, H3, H4, H5, H6, DIAG)
                    !
                CASE (3) ! HCP
                    !
                    CALL SET_HCP_BLOCK_MATRICES(H1, H2, H3, H4,&
                                              & H5, H6, H7, DIAG)
                    !
                CASE DEFAULT
                    !
                    CALL PAR_QUIT('Error  :     > Invalid crystal type&
                                 & in READ_SLIP_DATA.')
                    !
            END SELECT
            !
        END IF
        !
        ! Check if any latent parameters are invalid and if so quit.
        !
        IF (OPTIONS%HARD_TYPE .EQ. 'latent') THEN
            !
            SELECT CASE (CRYS_OPTIONS%CRYSTAL_TYPE(IPHASE))
                !
                CASE (1) ! FCC
                    !
                    IF (ANY(CRYS_OPTIONS%LATENT_PARAMETERS(IPHASE,1:5) .LT. 0.0D0)) THEN
                        !            
                        CALL PAR_QUIT('Error  :     > Invalid latent hardening&
                               & parameters provided for phase.')
                        !        
                    END IF
                    !
                CASE (2) ! BCC
                    !
                    IF (ANY(CRYS_OPTIONS%LATENT_PARAMETERS(IPHASE,1:7) .LT. 0.0D0)) THEN
                        !            
                        CALL PAR_QUIT('Error  :     > Invalid latent hardening&
                               & parameters provided for phase.')
                        !        
                    END IF
                    !
                CASE (3) ! HCP
                    !
                    IF (ANY(CRYS_OPTIONS%LATENT_PARAMETERS(IPHASE,1:8) .LT. 0.0D0)) THEN
                        !            
                        CALL PAR_QUIT('Error  :     > Invalid latent hardening&
                               & parameters provided for phase.')
                        !        
                    END IF
                    !
            END SELECT
            !
        END IF
        !
        ! Handle cyclic hardening parameters - if available.
        !
        IF (OPTIONS%HARD_TYPE .EQ. 'cyclic') THEN
            !
            ! Assign values to CYCLIC_PARM array.
            CYCLIC_PARM(0, IPHASE) = CRYS_OPTIONS%CYCLIC_PARAMETER_A(IPHASE)
            CYCLIC_PARM(1, IPHASE) = CRYS_OPTIONS%CYCLIC_PARAMETER_C(IPHASE)
            !
        END IF
        !
        ! Check if any cyclic parameters are invalid and if so quit.
        !
        IF (OPTIONS%HARD_TYPE .EQ. 'cyclic') THEN
            !
            IF (ANY(CYCLIC_PARM(:, IPHASE) .LT. 0.0D0)) THEN
                !            
                CALL PAR_QUIT('Error  :     > Invalid cyclic hardening parameters&
                    & provided for phase.')
                !        
            END IF
            !
        END IF
        !
        ! Print out crystal parameters to the console for user confirmation.
        !
        IF (MYID .EQ. 0) THEN
            !
            ! Need logic for crystal type - should be SELECT CASE
            IF (CRYS_OPTIONS%CRYSTAL_TYPE(IPHASE) .EQ. 1) THEN
                !
                WRITE(DFLT_U, '(A,I0,A)') 'Info   :     > phase ', &
                    & IPHASE, ' - crystal type: FCC'
                !
            ELSE IF (CRYS_OPTIONS%CRYSTAL_TYPE(IPHASE) .EQ. 2) THEN
                !
                WRITE(DFLT_U, '(A,I0,A)') 'Info   :     > phase ', &
                    & IPHASE, ' - crystal type: BCC'
                !
            ELSE
                !
                WRITE(DFLT_U, '(A,I0,A)') 'Info   :     > phase ', &
                    & IPHASE, ' - crystal type: HCP'
                !
            END IF
            !
            IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .FALSE.) THEN
                !
                WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > m           ', &
                    & CRYSTAL_PARM(0, IPHASE)
                !
            ELSE IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .TRUE.) THEN
                !
                WRITE(DFLT_U, '(A,3(E14.6))') 'Info   :     > m           ', &
                    & CRYS_OPTIONS%ANISO_M(IPHASE,1:3)
                !
            END IF
            !
            WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > gammadot_0  ', CRYSTAL_PARM(1, IPHASE)
            WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > h_0         ', CRYSTAL_PARM(2, IPHASE)
            !
            IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .FALSE.) THEN
                !
                WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > g_0         ', &
                    & CRYSTAL_PARM(3, IPHASE)
                !
            ELSE IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .TRUE.) THEN
                !
                WRITE(DFLT_U, '(A,3(E14.6))') 'Info   :     > g_0         ', &
                    & CRYSTAL_PARM(3, IPHASE), &
                    & HCP_RATIOS(2) * CRYSTAL_PARM(3, IPHASE), &
                    & HCP_RATIOS(1) * CRYSTAL_PARM(3, IPHASE)
                !
            END IF
            !
            WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > g_s0        ', CRYSTAL_PARM(4, IPHASE)
            WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > m_prime     ', CRYSTAL_PARM(5, IPHASE)
            WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > gammadot_s0 ', CRYSTAL_PARM(6, IPHASE)
            WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > n           ', N_VOCE(IPHASE)
            WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > c11         ', C11
            WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > c12         ', C12
            !
            IF (CRYS_OPTIONS%CRYSTAL_TYPE(IPHASE) .EQ. 3) THEN            
                !
                WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > c13         ', C13
                !
            ENDIF
            !            
            WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > c44         ', C44
            !
            ! Print out cyclic hardening parameters - if available.
            !
            IF (OPTIONS%HARD_TYPE .EQ. 'cyclic') THEN
                !
                WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > cyclic_param_a ',&
                    & CYCLIC_PARM(0, IPHASE)
                WRITE(DFLT_U, '(A,E14.6)') 'Info   :     > cyclic_param_c ',&
                    &CYCLIC_PARM(1, IPHASE)
                !
            END IF
            !
            ! Print out latent hardening parameters - if available.
            !
            IF (OPTIONS%HARD_TYPE .EQ. 'latent') THEN
                !
                SELECT CASE (CRYS_OPTIONS%CRYSTAL_TYPE(IPHASE))  
                    !
                    CASE (1) ! FCC 
                        !
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > diag        ', DIAG
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h1          ', H1
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h2          ', H2
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h3          ', H3
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h4          ', H4
                        !
                    CASE (2) ! BCC
                        !
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > diag        ', DIAG
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h1          ', H1
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h2          ', H2
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h3          ', H3
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h4          ', H4
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h5          ', H5 
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h6          ', H6
                        !
                    CASE (3) ! HCP
                        !
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > diag        ', DIAG
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h1          ', H1
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h2          ', H2
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h3          ', H3
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h4          ', H4
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h5          ', H5 
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h6          ', H6 
                        WRITE(DFLT_U,'(A,E14.6,A)') 'Info   :     > h7          ', H7
                        !
                END SELECT
                !
            END IF
            !
        END IF
        ! 
    END DO
    !
    ! Set the maximum number of slip systems.
    !
    CALL SET_MAXSLIP(NUMSLIP, MAXVERT)
    !
    END SUBROUTINE PROCESS_MATERIAL_PARAMETERS
    !
END MODULE MICROSTRUCTURE_MOD
