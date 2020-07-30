! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE BOUNDARY_CONDITIONS_MOD
!
! Module to calculate or read (from file) essential boundary conditions.
!
! Contains subroutines:
! CALC_BCS: Calculate grip, symmetry, or initial triaxial boundary conditions.
! READ_BCS: Read boundary conditions from file.
!
USE INTRINSICTYPESMODULE, ONLY: RK=>REAL_KIND
!
USE READ_INPUT_MOD, ONLY: COORDS, DOF_SUB1, DOF_SUP1, NP_SUB1, NP_SUP1
USE PARALLEL_MOD, ONLY: PAR_MAX, PAR_MIN, PAR_QUIT, PAR_MESSAGE, MYID
USE SIMULATION_CONFIGURATION_MOD, ONLY: BCS_OPTIONS, OPTIONS
USE UNITS_MOD
!
IMPLICIT NONE
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: CALC_BCS, READ_BCS
!
CONTAINS
    !
    SUBROUTINE CALC_BCS(GLOBAL_BCS, GLOBAL_VEL, GLOBAL_COORDS, GLOBAL_FORCE)
    !
    ! Calculate automatic grip or symmetry boundary conditions.
    ! Subroutine also initializes triaxial loading boundary conditions.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! GLOBAL_BCS: Global D.O.F. array indicating applied velocity bc's.
    ! GLOBAL_VEL: Global D.O.F. array storing velocity components.
    ! GLOBAL_COORDS: Global coordinates for the spatial mesh.
    !
    LOGICAL, INTENT(INOUT)  :: GLOBAL_BCS(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(INOUT) :: GLOBAL_VEL(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(IN)    :: GLOBAL_COORDS(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(INOUT) :: GLOBAL_FORCE(DOF_SUB1:DOF_SUP1)
    !
    ! Locals:
    ! CONTROL_FACE: Surface plane defined spatially opposite to LOADING_FACE.
    ! SYMMETRY_FACE_1: Surface plane orthogonal to CONTROL_FACE.
    ! SYMMETRY_FACE_2: Surface plane orthogonal to CONTROL_FACE.
    ! COORD_OFFSET: Offset that denotes which [X Y Z] components we loop over.
    ! SYMM_OFFSET_1: Offset that denotes which [X Y Z] components we loop over.
    ! SYMM_OFFSET_2: Offset that denotes which [X Y Z] components we loop over.
    ! I/J/II: Generic loop index values.
    ! MIN_SPATIAL_DIM: Global minimum [X Y Z] spatial component (CONTROL FACE).
    ! MAX_SPATIAL_DIM: Global maximum [X Y Z] spatial component (LOADING FACE).
    ! MIN_SYMM_DIM_1: Global minimum [X Y Z] spatial component (SYMM_FACE_1).
    ! MIN_SYMM_DIM_2: Global minimum [X Y Z] spatial component (SYMM_FACE_2).
    ! MIN_(X/Y/Z)_DIM: Global minimum [X Y Z] spatial component (ALL FACES).
    ! MAX_(X/Y/Z)_DIM: Global maximum [X Y Z] spatial component (ALL FACES).
    ! SWAP: Temporary variable for exchanging of MIN/MAX spatial dimensions.
    ! TEMP: Temporary array that holds [X Y Z] component for evaluation.
    ! MAX/MIN_FACE_DIM: Global max/min [X Y Z] spatial component (MINIMAL BCS).
    ! COORD_OFFSET_1/2: Offsets that denote the face normal directions of above.
    ! ROTATION_DIR: Direction that rigid body rotation is constrained in.
    !
    INTEGER                 :: CONTROL_FACE, COORD_OFFSET
    INTEGER                 :: SYMMETRY_FACE_1, SYMMETRY_FACE_2
    INTEGER                 :: SYMM_OFFSET_1, SYMM_OFFSET_2
    INTEGER                 :: I, J, II
    REAL(RK)                :: MIN_SPATIAL_DIM, MAX_SPATIAL_DIM
    REAL(RK)                :: MIN_SYMM_DIM_1, MIN_SYMM_DIM_2
    REAL(RK)                :: MIN_X_DIM, MAX_X_DIM, MIN_Y_DIM, MAX_Y_DIM
    REAL(RK)                :: MIN_Z_DIM, MAX_Z_DIM
    REAL(RK)                :: SWAP
    REAL(RK), ALLOCATABLE   :: TEMP(:)
    REAL(RK)                :: MAX_FACE_DIM, MIN_FACE_DIM
    INTEGER                 :: COORD_OFFSET_1, COORD_OFFSET_2, ROTATION_DIR
    !
    ! Parsed Values:
    ! LOADING_FACE: User-defined face which nodal velocities are applied.
    ! LOAD_DIRECTION: User-defined direction in which LOADING_FACE is displaced.
    ! BOUNDARY_CONDITIONS: User-defined boundary condition type.
    ! CONTROL_TYPE: User-defined deformation control type.
    ! STRAIN_RATE: User-defined engineering strain rate (1/s)
    ! NODAL_VELOCITY: Strain-rate computed velocity applied to LOADING_FACE.
    ! GAGE_LENGTH: Spatial distance between LOADING_FACE and CONTROL_FACE.
    !
    INTEGER                 :: LOADING_FACE, LOAD_DIRECTION, BOUNDARY_CONDITIONS
    INTEGER                 :: CONTROL_TYPE
    REAL(RK)                :: NODAL_VELOCITY, STRAIN_RATE, GAGE_LENGTH
    !
    ! Notes: 
    ! The mesh face set is always defined from [1 2 3 4 5 6] as
    ! [x_min x_max y_min y_max z_min z_max]. This should be consistent
    ! with the face set provided during mesh generation. The loading 
    ! directions are defined along the sample axes [X Y Z] as [0 1 2].
    !
    ! This subroutine by default assumes that displacement occurs in a 
    ! positive axis direction (that is, a domain loading in tension or 
    ! sheared along a positive coordinate axis). In order to perform 
    ! compression or negative axis direction shearing tests, the user 
    ! must sign the velocity in their *.loads or *.disp files.
    !
    !---------------------------------------------------------------------------
    !
    ! Initialize variables.
    GLOBAL_BCS = .FALSE.
    GLOBAL_VEL = 0.0_RK
    GLOBAL_FORCE = 0.0_RK
    !
    ! Read in boundary condition type from the options module.
    BOUNDARY_CONDITIONS = BCS_OPTIONS%BOUNDARY_CONDITIONS
    CONTROL_TYPE        = OPTIONS%DEF_CONTROL_BY
    !
    ! Error handling to ensure CONTROL_TYPE AND BOUNDARY_CONDITIONS agree.
    IF ((CONTROL_TYPE .EQ. 1) .AND. (BOUNDARY_CONDITIONS .EQ. 3)) THEN
        ! Uniaxial load target and triaxial bcs set
        CALL PAR_QUIT('Error  :     > Parameters "DEF_CONTROL_BY" and&
            & "BOUNDARY_CONDITIONS" do not agree.')
        !
    ELSE IF ((CONTROL_TYPE .EQ. 2) .AND. (BOUNDARY_CONDITIONS .EQ. 3)) THEN
        ! Uniaxial strain target and triaxial bcs set
        CALL PAR_QUIT('Error  :     > Parameters "DEF_CONTROL_BY" and&
            & "BOUNDARY_CONDITIONS" do not agree.')
        !
    ELSE IF ((CONTROL_TYPE .EQ. 3) .AND. (BOUNDARY_CONDITIONS .NE. 3)) THEN
        ! Triaxial CSR and uniaxial bcs set
        CALL PAR_QUIT('Error  :     > Parameters "DEF_CONTROL_BY" and&
            & "BOUNDARY_CONDITIONS" do not agree.')
        !
    ELSE IF ((CONTROL_TYPE .EQ. 4) .AND. (BOUNDARY_CONDITIONS .NE. 3)) THEN
        ! Triaxial CLR and uniaxial bcs set
        CALL PAR_QUIT('Error  :     > Parameters "DEF_CONTROL_BY" and&
            & "BOUNDARY_CONDITIONS" do not agree.')
        !
    ELSE
        ! Print notification to processor 0.
        IF (MYID .EQ. 0) THEN
            !
            WRITE(DFLT_U, '(A)') 'Info   :   - Initializing boundary conditions...'
            !
        END IF
        !         
    END IF
    !
    SELECT CASE (BOUNDARY_CONDITIONS)
    !
    !---------------------------------------------------------------------------
    !
    CASE (1) ! Begin Case (UNIAXIAL_GRIP)
        !
        LOADING_FACE    = BCS_OPTIONS%LOADING_FACE
        LOAD_DIRECTION  = BCS_OPTIONS%LOADING_DIRECTION
        STRAIN_RATE     = BCS_OPTIONS%STRAIN_RATE
        !
        ! Determine the control face based on the defined loading face.
        SELECT CASE (LOADING_FACE)
            !
            CASE (1)
                !
                CONTROL_FACE = 2
                COORD_OFFSET = 0
                !
            CASE (2)
                !
                CONTROL_FACE = 1
                COORD_OFFSET = 0
                !
            CASE (3)
                !
                CONTROL_FACE = 4
                COORD_OFFSET = 1
                !
            CASE (4)
                !
                CONTROL_FACE = 3
                COORD_OFFSET = 1
                !
            CASE (5)
                !
                CONTROL_FACE = 6
                COORD_OFFSET = 2
                !
            CASE (6)
                CONTROL_FACE = 5
                COORD_OFFSET = 2
                !
            CASE DEFAULT
                !
                CALL PAR_QUIT('Error  :     > Control face determination failed.')
                !
        END SELECT
        !
        ! Locate the extrema for our coordinate of interest in spatial mesh.
        ALLOCATE(TEMP(0:SIZE(GLOBAL_COORDS)/3))
        TEMP = 0.0_RK
        J = 0
        !
        DO I = DOF_SUB1, DOF_SUP1, 3
            !
            TEMP(J) = GLOBAL_COORDS(I + COORD_OFFSET)
            J = J + 1
            !
        END DO
        !
        CALL PAR_MAX(MAXVAL(TEMP), MAX_SPATIAL_DIM)
        CALL PAR_MIN(MINVAL(TEMP), MIN_SPATIAL_DIM)
        !
        DEALLOCATE(TEMP)
        !
        ! Convert the strain rate to a velocity.
        GAGE_LENGTH = MAX_SPATIAL_DIM - MIN_SPATIAL_DIM
        NODAL_VELOCITY = STRAIN_RATE * GAGE_LENGTH
        !
        ! If LOADING_FACE is a minimum swap MIN/MAX_SPATIAL_DIMS.
        ! Also, sign the NODAL_VELOCITY to always have positive displacement.
        IF (MODULO(LOADING_FACE,2) .EQ. 1) THEN
            !
            SWAP = MAX_SPATIAL_DIM
            MAX_SPATIAL_DIM = MIN_SPATIAL_DIM
            MIN_SPATIAL_DIM = SWAP
            NODAL_VELOCITY = -1 * NODAL_VELOCITY
            !
        END IF
        !
        ! Loop over global coords as triplets and determine if a node is on
        ! a BC surface. If it is apply the appropriate values.
        DO I = DOF_SUB1, DOF_SUP1, 3
            ! If the coordinate component is on the CONTROL_FACE.
            IF (GLOBAL_COORDS(I+COORD_OFFSET) .EQ. MIN_SPATIAL_DIM) THEN
                !
                GLOBAL_BCS(I)   = .TRUE.
                GLOBAL_BCS(I+1) = .TRUE.
                GLOBAL_BCS(I+2) = .TRUE.
                !
                GLOBAL_VEL(I)   = 0.0_RK
                GLOBAL_VEL(I+1) = 0.0_RK
                GLOBAL_VEL(I+2) = 0.0_RK
                !
            ! If the coordinate component is on the LOADING_FACE.
            ELSE IF ((GLOBAL_COORDS(I+COORD_OFFSET)) .EQ. MAX_SPATIAL_DIM) THEN
                !
                GLOBAL_BCS(I)   = .TRUE.
                GLOBAL_BCS(I+1) = .TRUE.
                GLOBAL_BCS(I+2) = .TRUE.
                !
                GLOBAL_VEL(I + LOAD_DIRECTION) = NODAL_VELOCITY
                !
            END IF
            !
        END DO
        !
    ! End Case (UNIAXIAL_GRIP)
    !
    !---------------------------------------------------------------------------
    !
    CASE (2) ! Begin Case (UNIAXIAL_SYMMETRY)
        !
        LOAD_DIRECTION  = BCS_OPTIONS%LOADING_DIRECTION
        STRAIN_RATE     = BCS_OPTIONS%STRAIN_RATE
        !
        ! Determine face information based on user-defined loading direction.
        SELECT CASE (LOAD_DIRECTION)
            !
            CASE (0) ! X_MAX
                !
                CONTROL_FACE = 1
                LOADING_FACE = 2
                SYMMETRY_FACE_1 = 3
                SYMMETRY_FACE_2 = 5
                SYMM_OFFSET_1 = 1
                SYMM_OFFSET_2 = 2
                !
            CASE (1) ! Y_MAX
                ! 
                CONTROL_FACE = 3
                LOADING_FACE = 4
                SYMMETRY_FACE_1 = 1
                SYMMETRY_FACE_2 = 5
                SYMM_OFFSET_1 = 0
                SYMM_OFFSET_2 = 2
                !
            CASE (2) ! Z_MAX
                !
                CONTROL_FACE = 5
                LOADING_FACE = 6
                SYMMETRY_FACE_1 = 1
                SYMMETRY_FACE_2 = 3
                SYMM_OFFSET_1 = 0
                SYMM_OFFSET_2 = 1
                !
            CASE DEFAULT
                !
                CALL PAR_QUIT('Error  :     > Control face determination failed.')
                !
        END SELECT
        !
        ! Override LOADING_FACE from options to ensure correct behavior in
        ! loading control routines.
        BCS_OPTIONS%LOADING_FACE = LOADING_FACE
        !
        ! Locate the extrema for our coordinate of interest in spatial mesh.
        ALLOCATE(TEMP(0:SIZE(GLOBAL_COORDS)/3))
        TEMP = 0.0_RK
        ! First, consider the loading and control faces.
        J = 0
        DO I = DOF_SUB1, DOF_SUP1, 3
            !
            TEMP(J) = GLOBAL_COORDS(I + LOAD_DIRECTION)
            J = J + 1
            !
        END DO
        !
        CALL PAR_MAX(MAXVAL(TEMP), MAX_SPATIAL_DIM)
        CALL PAR_MIN(MINVAL(TEMP), MIN_SPATIAL_DIM)
        !
        ! Convert the strain rate to a velocity.
        GAGE_LENGTH = MAX_SPATIAL_DIM - MIN_SPATIAL_DIM
        NODAL_VELOCITY = STRAIN_RATE * GAGE_LENGTH
        ! Next, consider the first symmetry face.
        J = 0
        DO I = DOF_SUB1, DOF_SUP1, 3
            !
            TEMP(J) = GLOBAL_COORDS(I + SYMM_OFFSET_1)
            J = J + 1
            !
        END DO
        !
        CALL PAR_MIN(MINVAL(TEMP), MIN_SYMM_DIM_1)
        !
        ! Then, consider the second symmetry face.
        J = 0
        DO I = DOF_SUB1, DOF_SUP1, 3
            !
            TEMP(J) = GLOBAL_COORDS(I + SYMM_OFFSET_2)
            J = J + 1
            !
        END DO
        !
        CALL PAR_MIN(MINVAL(TEMP), MIN_SYMM_DIM_2)
        !
        DEALLOCATE(TEMP)
        !
        ! If LOADING_FACE is a minimum swap MIN/MAX_SPATIAL_DIMS.
        ! Also, sign the NODAL_VELOCITY to always have positive displacement.
        IF (MODULO(LOADING_FACE,2) .EQ. 1) THEN
            !
            SWAP = MAX_SPATIAL_DIM
            MAX_SPATIAL_DIM = MIN_SPATIAL_DIM
            MIN_SPATIAL_DIM = SWAP
            NODAL_VELOCITY = -1 * NODAL_VELOCITY
            !
        END IF
        !
        ! Loop over global coords as triplets and determine if a node is on
        ! a BC surface. If it is apply the appropriate values.
        ! Note that multiple DO loops are utilized to handle edge cases.
        !
        DO I = DOF_SUB1, DOF_SUP1, 3
            ! If the coordinate component is on the SYMMETRY_FACE_1.
            IF ((GLOBAL_COORDS(I+SYMM_OFFSET_1)) .EQ. MIN_SYMM_DIM_1) THEN
                !
                GLOBAL_BCS(I+SYMM_OFFSET_1) = .TRUE.
                !
                GLOBAL_VEL(I+SYMM_OFFSET_1) = 0.0_RK
                !
            END IF
        END DO
        !
        DO I = DOF_SUB1, DOF_SUP1, 3
            ! If the coordinate component is on the SYMMETRY_FACE_2.
            IF ((GLOBAL_COORDS(I+SYMM_OFFSET_2)) .EQ. MIN_SYMM_DIM_2) THEN
                !
                GLOBAL_BCS(I+SYMM_OFFSET_2) = .TRUE.
                !
                GLOBAL_VEL(I+SYMM_OFFSET_2) = 0.0_RK
                !
            END IF
        END DO
        !
        DO I = DOF_SUB1, DOF_SUP1, 3
            ! If the coordinate component is on the CONTROL_FACE.
            IF (GLOBAL_COORDS(I+LOAD_DIRECTION) .EQ. MIN_SPATIAL_DIM) THEN
                !
                GLOBAL_BCS(I + LOAD_DIRECTION) = .TRUE.
                !
                GLOBAL_VEL(I + LOAD_DIRECTION) = 0.0_RK
                !
            ! If the coordinate component is on the LOADING_FACE.
            ELSE IF ((GLOBAL_COORDS(I+LOAD_DIRECTION)) .EQ. MAX_SPATIAL_DIM) THEN
                !
                GLOBAL_BCS(I)   = .TRUE.
                GLOBAL_BCS(I+1) = .TRUE.
                GLOBAL_BCS(I+2) = .TRUE.
                !
                GLOBAL_VEL(I + LOAD_DIRECTION) = NODAL_VELOCITY
                !
            END IF
            !
        END DO
        !
    ! End Case (UNIAXIAL_SYMMETRY)
    !
    !---------------------------------------------------------------------------
    !
    CASE (3) ! Begin Case (TRIAXIAL)
        !
        ! Note:
        ! Actual boundary conditions (i.e. velocities) are calculated in
        ! the corresponding triaxial driver modules. Initial conditions
        ! are expected to be such that all surface normal DOFs are prescribed
        ! and set to zero.
        !
        ! Initialize dimension variables.
        MIN_X_DIM = 0.0_RK
        MAX_X_DIM = 0.0_RK
        MIN_Y_DIM = 0.0_RK
        MAX_Y_DIM = 0.0_RK
        MIN_Z_DIM = 0.0_RK
        MAX_Z_DIM = 0.0_RK
        !
        ! Locate the extrema for our coordinate of interest in spatial mesh.
        ALLOCATE(TEMP(0:SIZE(GLOBAL_COORDS)/3))
        TEMP = 0.0_RK
        !
        ! Loop over all load directions to find surface nodes.
        DO LOAD_DIRECTION = 0, 2
            J = 0
            DO I = DOF_SUB1, DOF_SUP1, 3
                !
                TEMP(J) = GLOBAL_COORDS(I + LOAD_DIRECTION)
                J = J + 1
                !
            END DO
            !
            IF (LOAD_DIRECTION .EQ. 0) THEN
                !
                CALL PAR_MIN(MINVAL(TEMP), MIN_X_DIM)
                CALL PAR_MAX(MAXVAL(TEMP), MAX_X_DIM)
                !
            ELSE IF (LOAD_DIRECTION .EQ. 1) THEN
                !
                CALL PAR_MIN(MINVAL(TEMP), MIN_Y_DIM)
                CALL PAR_MAX(MAXVAL(TEMP), MAX_Y_DIM)
                !
            ELSE IF (LOAD_DIRECTION .EQ. 2) THEN
                !
                CALL PAR_MIN(MINVAL(TEMP), MIN_Z_DIM)
                CALL PAR_MAX(MAXVAL(TEMP), MAX_Z_DIM)
                !
            END IF
            !
        END DO
        !
        DEALLOCATE(TEMP)
        !
        ! Loop over global coords as triplets and determine if a node is on
        ! a BC surface. If it is apply the appropriate values.
        ! Note that multiple DO loops are utilized to handle edge cases.
        !
        ! First, check the X faces.
        DO I = DOF_SUB1, DOF_SUP1, 3
            !
            IF (GLOBAL_COORDS(I + 0) .EQ. MIN_X_DIM) THEN
                !
                GLOBAL_BCS(I + 0) = .TRUE.
                !
                GLOBAL_VEL(I + 0) = 0.0_RK
                !
            ELSE IF ((GLOBAL_COORDS(I + 0)) .EQ. MAX_X_DIM) THEN
                !
                GLOBAL_BCS(I + 0) = .TRUE.
                !
                GLOBAL_VEL(I + 0) = 0.0_RK
                !
            END IF
            !
        END DO
        ! Next, check the Y faces.
        DO I = DOF_SUB1, DOF_SUP1, 3
            !
            IF (GLOBAL_COORDS(I + 1) .EQ. MIN_Y_DIM) THEN
                !
                GLOBAL_BCS(I + 1) = .TRUE.
                !
                GLOBAL_VEL(I + 1) = 0.0_RK
                !
            ELSE IF ((GLOBAL_COORDS(I + 1)) .EQ. MAX_Y_DIM) THEN
                !
                GLOBAL_BCS(I + 1) = .TRUE.
                !
                GLOBAL_VEL(I + 1) = 0.0_RK
                !
            END IF
            !
        END DO
        ! Finally, check the Z faces.
        DO I = DOF_SUB1, DOF_SUP1, 3
            !
            IF (GLOBAL_COORDS(I + 2) .EQ. MIN_Z_DIM) THEN
                !
                GLOBAL_BCS(I + 2) = .TRUE.
                !
                GLOBAL_VEL(I + 2) = 0.0_RK
                !
            ELSE IF ((GLOBAL_COORDS(I + 2)) .EQ. MAX_Z_DIM) THEN
                !
                GLOBAL_BCS(I + 2) = .TRUE.
                !
                GLOBAL_VEL(I + 2) = 0.0_RK
                !
            END IF
            !
        END DO        
        !
    ! End Case (TRIAXIAL)
    !
    !---------------------------------------------------------------------------
    !
    CASE (4) ! Begin Case (MINIMAL)
        !
        ! Note: This boundary condition is another implementation of uniaxial
        ! grip however it has `minimal` nodal constraints. The loading and 
        ! control face are constrained in the normal directions and two nodes
        ! are constrained to prevent rigid body motion. All other DOFs are free.
        !
        LOAD_DIRECTION  = BCS_OPTIONS%LOADING_DIRECTION
        STRAIN_RATE     = BCS_OPTIONS%STRAIN_RATE
        ! 
        ! Determine face information based on user-defined loading direction.
        SELECT CASE (LOAD_DIRECTION)
            !
            CASE (0) ! X_MAX
                !
                CONTROL_FACE = 1 ! x0
                LOADING_FACE = 2 ! x1
                !
            CASE (1) ! Y_MAX
                !
                CONTROL_FACE = 3 ! y0
                LOADING_FACE = 4 ! y1
                !
            CASE (2) ! Z_MAX
                !
                CONTROL_FACE = 5 ! z0
                LOADING_FACE = 6 ! z1
                !
            CASE DEFAULT
                !
                CALL PAR_QUIT('Error  :     > Control face determination failed.')
                !
        END SELECT
        !
        ! Override LOADING_FACE from options to ensure correct behavior in
        ! loading control routines.
        BCS_OPTIONS%LOADING_FACE = LOADING_FACE
        !
        ! Initialize dimension variables.
        MIN_X_DIM = 0.0_RK
        MAX_X_DIM = 0.0_RK
        MIN_Y_DIM = 0.0_RK
        MAX_Y_DIM = 0.0_RK
        MIN_Z_DIM = 0.0_RK
        MAX_Z_DIM = 0.0_RK
        !
        ! Locate the extrema for our coordinate of interest in spatial mesh.
        ALLOCATE(TEMP(0:SIZE(GLOBAL_COORDS)/3))
        TEMP = 0.0_RK
        !
        ! Loop over all load directions to find surface nodes.
        DO I = 0, 2
            !
            J = 0
            DO II = DOF_SUB1, DOF_SUP1, 3
                !
                TEMP(J) = GLOBAL_COORDS(II + LOAD_DIRECTION)
                J = J + 1
                !
            END DO
            !
            IF (I .EQ. 0) THEN
                !
                CALL PAR_MIN(MINVAL(TEMP), MIN_X_DIM)
                CALL PAR_MAX(MAXVAL(TEMP), MAX_X_DIM)
                !
            ELSE IF (I .EQ. 1) THEN
                !
                CALL PAR_MIN(MINVAL(TEMP), MIN_Y_DIM)
                CALL PAR_MAX(MAXVAL(TEMP), MAX_Y_DIM)
                !
            ELSE IF (I .EQ. 2) THEN
                !
                CALL PAR_MIN(MINVAL(TEMP), MIN_Z_DIM)
                CALL PAR_MAX(MAXVAL(TEMP), MAX_Z_DIM)
                !
            END IF
            !
        END DO
        !
        DEALLOCATE(TEMP)
        !
        ! Set the min/max spatial dim variables here to avoid issues
        ! When I set them in the DO loop above it had unexpected behavior - JC
        IF (LOAD_DIRECTION .EQ. 0) THEN
            !
            MIN_SPATIAL_DIM = MIN_X_DIM
            MAX_SPATIAL_DIM = MAX_X_DIM
            ! Set the nodal velocity while we are here
            GAGE_LENGTH = MAX_X_DIM - MIN_X_DIM
            NODAL_VELOCITY = STRAIN_RATE * GAGE_LENGTH  
            !
        ELSE IF (LOAD_DIRECTION .EQ. 1) THEN
            !
            MIN_SPATIAL_DIM = MIN_Y_DIM
            MAX_SPATIAL_DIM = MAX_Y_DIM
            !
            GAGE_LENGTH = MAX_Y_DIM - MIN_Y_DIM
            NODAL_VELOCITY = STRAIN_RATE * GAGE_LENGTH 
            !
        ELSE IF (LOAD_DIRECTION .EQ. 2) THEN
            !
            MIN_SPATIAL_DIM = MIN_Z_DIM
            MAX_SPATIAL_DIM = MAX_Z_DIM
            !
            GAGE_LENGTH = MAX_Z_DIM - MIN_Z_DIM
            NODAL_VELOCITY = STRAIN_RATE * GAGE_LENGTH
            !
        END IF
        !
        ! Loop over global coords as triplets and determine if a node is on
        ! a BC surface. If it is apply the appropriate values.
        ! First constrain the normal directions on the control and loading face.
        DO I = DOF_SUB1, DOF_SUP1, 3
            ! If the coordinate component is on the CONTROL_FACE.
            IF (GLOBAL_COORDS(I+LOAD_DIRECTION) .EQ. MIN_SPATIAL_DIM) THEN
                !
                GLOBAL_BCS(I + LOAD_DIRECTION) = .TRUE.
                !
                GLOBAL_VEL(I + LOAD_DIRECTION) = 0.0_RK
                !
            ! If the coordinate component is on the LOADING_FACE.
            ELSE IF ((GLOBAL_COORDS(I+LOAD_DIRECTION)) .EQ. MAX_SPATIAL_DIM) THEN
                !
                GLOBAL_BCS(I + LOAD_DIRECTION) = .TRUE.
                !
                GLOBAL_VEL(I + LOAD_DIRECTION) = NODAL_VELOCITY
                !
            END IF
            !
        END DO
        !
        ! Next, locate the appropriate two nodes to constrain to prevent RBM.
        ! Node 1 is ALWAYS fixed where the minimum X, Y, and Z faces converge.
        DO I = DOF_SUB1, DOF_SUP1, 3
            !
            IF (GLOBAL_COORDS(I+0) .EQ. MIN_X_DIM) THEN
                !
                IF (GLOBAL_COORDS(I+1) .EQ. MIN_Y_DIM) THEN
                    !
                    IF (GLOBAL_COORDS(I+2) .EQ. MIN_Z_DIM) THEN
                        !
                        ! If we reach here then fixed this nodes values.
                        GLOBAL_BCS(I+0) = .TRUE.
                        GLOBAL_BCS(I+1) = .TRUE.
                        GLOBAL_BCS(I+2) = .TRUE.
                        !
                        GLOBAL_VEL(I+0) = 0.0_RK
                        GLOBAL_VEL(I+1) = 0.0_RK
                        GLOBAL_VEL(I+2) = 0.0_RK
                        !
                    END IF
                    !
                END IF
                !
            END IF
            !
        END DO
        !
        ! Node 2 is additionally constrained in the direction that would cause rigid
        ! body rotation about the loading direction axis. It is constrained
        ! on the maximum face on the same edge as Node 1.
        !
        ! Abstract the RHS variables we need to check to be load direction dependent.
        IF (LOAD_DIRECTION .EQ. 0) THEN ! X-direction
            !
            MAX_FACE_DIM    = MAX_Z_DIM ! z1
            MIN_FACE_DIM    = MIN_Y_DIM ! y0
            !
            COORD_OFFSET_1 = 2 ! z-face normal
            COORD_OFFSET_2 = 1 ! y-face normal
            ROTATION_DIR   = 1 ! Fixed in Y direction
            !
        ELSE IF (LOAD_DIRECTION .EQ. 1) THEN ! Y-direction
            !
            MAX_FACE_DIM    = MAX_Z_DIM ! z1
            MIN_FACE_DIM    = MIN_X_DIM ! x0
            !
            COORD_OFFSET_1 = 2 ! z-face normal
            COORD_OFFSET_2 = 0 ! x-face normal
            ROTATION_DIR   = 0 ! Fixed in X direction
            !
        ELSE IF (LOAD_DIRECTION .EQ. 2) THEN ! Z-direction
            !
            MAX_FACE_DIM    = MAX_Y_DIM ! y1
            MIN_FACE_DIM    = MIN_X_DIM ! x0
            !
            COORD_OFFSET_1 = 1 ! y-face normal
            COORD_OFFSET_2 = 0 ! x-face normal
            ROTATION_DIR   = 0 ! Fixed in X direction
            !
        END IF
        !
        DO I = DOF_SUB1, DOF_SUP1, 3
            ! Check if we are on the control face
            IF (GLOBAL_COORDS(I+LOAD_DIRECTION) .EQ. MIN_SPATIAL_DIM) THEN
                ! Check if we are on MIN_FACE
                IF (GLOBAL_COORDS(I+COORD_OFFSET_2) .EQ. MIN_FACE_DIM) THEN
                    ! Check if we are on MAX_FACE
                    IF (GLOBAL_COORDS(I+COORD_OFFSET_1) .EQ. MAX_FACE_DIM) THEN
                        !
                        ! If we reach here then fixed this nodes values.
                        GLOBAL_BCS(I+ROTATION_DIR) = .TRUE.
                        !
                        GLOBAL_VEL(I+ROTATION_DIR) = 0.0_RK
                        !
                    END IF
                    !
                END IF
                !
            END IF
            !
        END DO
        !
    ! End Case (MINIMAL)
    !
    !---------------------------------------------------------------------------
    !
    CASE DEFAULT
        ! Exit the program if input is invalid.
        CALL PAR_QUIT('Error  :     > User-defined BOUNDARY_CONDITIONS parameter&
            & is invalid.')
        !
    ! End Case (DEFAULT)
    !
    !---------------------------------------------------------------------------
    !
    END SELECT
    !
    END SUBROUTINE CALC_BCS
    !
    !===========================================================================
    !
    SUBROUTINE READ_BCS(BCS, VELOCITY, FORCE)
    !
    ! Read boundary conditions from file.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! BCS: Array indicating essential boundary conditions (T/F).
    ! VELOCITY: Array containing nodal velocities.
    ! FORCE: Forces for force boundary conditions (not yet implemented).
    !
    LOGICAL, INTENT(INOUT)  ::  BCS(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(INOUT) ::  VELOCITY(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(INOUT) ::  FORCE(DOF_SUB1:DOF_SUP1)
    !
    ! Locals:
    ! NUMUBC, I: Indices for number of boundary conditions (from file) and loop.
    ! II: Dummy index for substring parsing.
    ! J: Index to loop over parsed string constraint directions.
    ! INODE: Index of node with boundary condition.
    ! K: Index in global array for INODE.
    ! NSPACE, NCOMMA: Parsed number of space and comma delimiter for a line.
    ! NVALS: Number of unique substrings to be handled per line.
    ! IOSTATUS, IERR: Error handling for FORTRAN read to check for EOF.
    ! VAL: Applied velocity in line specified directions (from file).
    ! IOFILE: Default boundary conditions file name that is searched for.
    ! MESSAGE: Statement passed to processor 0 to notify user of read-in.
    ! DIR: Direction to constrain node on current line.
    ! IARRAY: Temporary substring array for line parsing and var assignment.
    !
    INTEGER  :: NUMUBC, I, II, J, INODE, K, NSPACE, NCOMMA, NVALS
    INTEGER  :: IOSTATUS, IERR
    REAL(RK) :: VAL
    CHARACTER(LEN=256) :: IOFILE, MESSAGE
    CHARACTER(LEN=256) :: LINE
    CHARACTER(LEN=5) :: DIR
    CHARACTER(LEN=8) :: IARRAY(5)
    !
    !---------------------------------------------------------------------------
    !
    ! Initialize variables.
    BCS = .FALSE.
    VELOCITY = 0.0_RK
    FORCE = 0.0_RK
    !
    ! If user does not provide a specific file name to read in check default.
    IF (BCS_OPTIONS%BCS_FILE .EQ. '') THEN
        !
        IF (MYID .EQ. 0) THEN
            !
            WRITE(DFLT_U, '(A)') 'Info   :   [i] Parsing bcs file `simulation.bcs`...'
            IOFILE = 'simulation.bcs'
            !
        END IF
        !
    ELSE IF (LEN_TRIM(BCS_OPTIONS%BCS_FILE) .GT. 0) THEN
        !
        IOFILE = BCS_OPTIONS%BCS_FILE
        !
    ELSE
        !
        CALL PAR_QUIT('Error  :     > Failure to locate boundary conditions (*.bcs) file.')
        !    
    END IF
    !
    ! Look into usage of STANDARD_INPUT from LIBF95 OR UNITS MOD
    OPEN(IUNITS(BCS_U), FILE = IOFILE, STATUS='OLD', ACTION='READ', IOSTAT=IOSTATUS)
    !
    IF (IOSTATUS .NE. 0) THEN
        !
        CALL PAR_QUIT('Error  :     > Failure to open boundary conditions (*.bcs) file.')
        !
    ELSE
        !
        IF (MYID .EQ. 0) THEN
            !
            WRITE(DFLT_U, '(A,A,A)') 'Info   :   [i] Parsing bcs file `', TRIM(ADJUSTL(IOFILE)),'`...'
            !
        END IF
        !
    END IF
    !
    ! Read the number of lines in the file to set NUMUBC
    ! This was previously read in from the *.bcs file explicitly
    !
    NUMUBC = 0
    !
    DO ! Count total number of lines in file
        !
        READ(IUNITS(BCS_U), *, IOSTAT = IERR)
        !
        IF (IERR .NE. 0) EXIT
        NUMUBC = NUMUBC + 1
        !
    END DO
    !
    ! Rewind the file to read in the orientations to an allocated array
    REWIND (IUNITS(BCS_U))
    !
    DO I = 1, NUMUBC
        !
        READ(IUNITS(BCS_U), '(A)') LINE
        !
        ! Trim the full string into NVALS number of substrings to parse
        NSPACE = COUNT( (/ (LINE(II:II), II=1, LEN_TRIM(LINE)) /) == " ")
        NCOMMA = COUNT( (/ (LINE(II:II), II=1, LEN_TRIM(LINE)) /) == ",")
        NVALS = (NSPACE + NCOMMA) + 1
        !
        ! Internal read of line to store substrings
        READ(LINE, *) IARRAY(1:NVALS)
        !
        ! Loop over the number of directions on this node line and assign BCs
        DO J = 2, (NVALS - 1)
            !
            ! Reread the substring array and assign values
            READ(IARRAY(1), *) INODE ! INODE is fixed over outer index I
            INODE = INODE - 1 ! Need to shift the value to read in 1-indexed BCs
            READ(IARRAY(J), *) DIR
            READ(IARRAY(NVALS), *) VAL ! VAL is fixed over outer index I
            !
            ! Assign node direction/value pair into corresponding DOF array location
            IF ( (INODE .GE. NP_SUB1) .AND. (INODE .LE. NP_SUP1) ) THEN
                !
                IF ((TRIM(DIR) .EQ. 'X') .OR. (TRIM(DIR) .EQ. 'x')) THEN
                    !
                    K = 3 * INODE
                    !
                ELSE IF ((TRIM(DIR) .EQ. 'Y') .OR. (TRIM(DIR) .EQ. 'y')) THEN
                    !
                    K = 3 * INODE + 1
                    !
                ELSE IF ((TRIM(DIR) .EQ. 'Z') .OR. (TRIM(DIR) .EQ. 'z')) THEN
                    !
                    K = 3 * INODE + 2
                    !
                ELSE
                    !
                    CALL PAR_QUIT('Fatal error: Unknown direction character&
                        & provided in user-defined boundary conditions.')
                    !
                END IF
                !
                ! Check that the node has not already been assigned a value
                IF ( BCS(K) .NEQV. .FALSE. ) THEN
                    !
                    WRITE(MESSAGE, '(a,i0,a)') 'Fatal error: Node number ', INODE, &
                        & 'has been assigned multiple boundary conditions.'
                    !
                    CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                    !
                ELSE ! If there are no issues, assign the BC
                    !
                    BCS(K) = .TRUE.
                    VELOCITY(K) = VAL
                    !
                END IF
                !
            ENDIF
            !
        END DO
        !
    END DO
    !
    !
    IF (MYID .EQ. 0) THEN
        !
        WRITE(DFLT_U, '(A,A,A)') 'Info   :   [i] Parsed file `', TRIM(ADJUSTL(IOFILE)),'`.'
        !
    END IF
    !
    CLOSE(UNIT = IUNITS(BCS_U))
    !
    END SUBROUTINE READ_BCS
    !
END MODULE BOUNDARY_CONDITIONS_MOD
