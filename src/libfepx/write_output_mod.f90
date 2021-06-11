! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE WRITE_OUTPUT_MOD
!
! Module for printing variables to file at the end of a step.
!
! Contains subroutines:
!
! General printing and handling of field variable data:
! PRINT_STEP: Print variables on a given load step
!
! Handling the post.report, post.forceX, and post.conv files
! WRITE_REPORT_FILE_HEADER: Writes preamble information about the simulation
! WRITE_REPORT_FILE_OUTPUT_FILES: Writes the files requested to be printed
! WRITE_REPORT_FILE_COMPLETE_STEPS: Writes the LAST completed step number
! WRITE_FORCE_FILE_HEADERS: Writes the file headers for table formatting
! WRITE_FORCE_FILE_DATA: Writes the surface forces [X Y Z] for all surfaces
! WRITE_CONV_FILE_HEADERS: Writes the file headers for the table formatting
! WRITE_CONV_FILE_DATA: Writes the various convergence statistics
!
! Writing restart files from the various drivers:
! WRITE_RESTART_FIELD: Writes field data for restarting a simulation.
! WRITE_UNIAXIAL_RESTART: Writes required information for the uniaxial restart
! WRITE_TRIAXCSR_RESTART: Writes required information for the triaxcsr restart
! WRITE_TRIAXCLR_RESTART: Writes required information for the triaxclr restart
!
! To do:
! - is PROB_TYPE necessary for PRINT_STEP?
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
!
! From libfepx:
!
USE DIMENSIONS_MOD
USE MATRIX_OPERATIONS_MOD
USE MICROSTRUCTURE_MOD
USE ORIENTATION_CONVERSION_MOD
USE READ_INPUT_MOD
USE SURFACE_MOD, ONLY: NSURFACES
USE UNITS_MOD
!
! From libparallel
!
USE PARALLEL_MOD
!
IMPLICIT NONE
!
! Public
!
PUBLIC
!
CONTAINS
    !
    SUBROUTINE PRINT_STEP(STEP, PROB_TYPE, CRD, VEL, ORIENT, HARD, ELVOL,&
        & ELAS, S3X3, DEFF, EQSTRAIN, DPEFF, EQPLSTRAIN, VGRAD, DP_HAT, &
        & WP_HAT, GAMMA, GAMMADOT, WORK, WORKP, D_TOT, WORKRATE, WORKRATEP, &
        & EQELSTRAIN, PLSTRAIN, TOTSTRAIN)
    !
    ! Print variables on a given load step
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! STEP: Step number
    ! PROB_TYPE: Type of problem (iso(0), aniso-vp(1), aniso-evps(2)) Necessary?
    ! CRD: Nodal coordinates
    ! VEL: Nodal velocities
    ! ORIENT: Elemental orientation tensors (to be converted)
    ! HARD: Elemental CRSS values
    ! ELVOL: Elemental volume
    ! ELAS: Elemental elastic strain tensors (upper triangle, i.e. 6-vector)
    ! S3X3: Elemental stress tensors (3x3 matrix)
    ! DEFF: Elemental effective deformation rate (scalar)
    ! EQSTRAIN: Elemental equivalent strain (scalar)
    ! DPEFF: Elemental effective plastic deformation rate (scalar)
    ! EQPLSTRAIN: Elemental equivalent plastic strain (scalar)
    ! VGRAD: Elemental velocity gradient (3x3 matrix)
    ! DP_HAT: Elemental plastic deformation rate (5-vector)
    ! WP_HAT: Elemental plastic spin rate (3-vector)
    ! GAMMA: Elemental accumulated shears
    ! WORK: Elemental total work (scalar)
    ! WORKP: Elemental plastic work (scalar)
    ! D_TOT: Elemental total deformation rate tensor (3x3 matrix)
    ! WORKRATE: Elemental total work rate (scalar)
    ! WORKRATEP: Elemental plastic work rate (scalar)
    ! EQELSTRAIN: Elemental equivalent elastic strain (scalar)
    ! PLSTRAIN: Elemental plastic strain tensor (5-vector)
    ! TOTSTRAIN: Elemental total strain tensor (3x3 matrix)
    !
    INTEGER, INTENT(IN) :: STEP
    INTEGER, INTENT(IN) :: PROB_TYPE
    REAL(RK), INTENT(IN) :: CRD(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(IN) :: VEL(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(IN) :: ORIENT(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: HARD(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: ELVOL(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: ELAS(0:TVEC, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: S3X3(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: EQSTRAIN(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: DEFF(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: DPEFF(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: EQPLSTRAIN(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: VGRAD (0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: DP_HAT(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: WP_HAT(0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: GAMMA(0:MAXSLIP1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: GAMMADOT(0:MAXSLIP1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: WORK(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: WORKP(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: D_TOT(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: WORKRATE(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: WORKRATEP(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: EQELSTRAIN(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: PLSTRAIN(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN), OPTIONAL :: TOTSTRAIN(0:DIMS1, 0:DIMS1, &
        & EL_SUB1:EL_SUP1)
    !
    ! Locals:
    ! IO: File control integer
    ! I, J, IGRAIN: Looping indices
    ! M_EL: Number of elements
    ! P: Passive (1) or active (-1) orientation convention
    ! ANGLE, AXIS: Angle-axis arrays
    ! PSI1, PHI PSI2: Euler-Bunge arrays
    ! PSI, THE, (PHI): Euler-Kocks arrays
    ! RODS: Rodrigues' vectors
    ! QUAT: Quaternions
    ! ELAS_SAM: Elastic strain in the sample basis
    ! EQSTRESS: Equivalent (von Mises) stress
    !
    INTEGER :: IO
    INTEGER :: I, J, IGRAIN
    INTEGER :: M_EL
    LOGICAL :: PFLAG_NODE = .TRUE.
    LOGICAL :: PFLAG_ELEM = .TRUE.
    REAL(RK) :: ORIENT_INT(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ANGLE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: AXIS(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: PSI1(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: PHI(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: PSI2(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: PSI(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: THE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: RODS(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: QUAT(0:DIMS, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ELAS_SAM(0:5)
    REAL(RK) :: PLSTRAINMAT(0:DIMS1, 0:DIMS1)
    REAL(RK) :: PLSTRAIN6(0:5)
    REAL(RK) :: PLSTRAIN6_SAM(0:5)
    REAL(RK) :: EQSTRESS(EL_SUB1:EL_SUP1)
    REAL(RK) :: DP_HAT_MAT(0:5)
    REAL(RK) :: DP_HAT_SAM(0:5)
    REAL(RK) :: WP_HAT_SAM(0:DIMS1)
    REAL(RK) :: QR5X5(0:TVEC1, 0:TVEC1)
    REAL(RK) :: QR3X3(0:DIMS1, 0:DIMS1)
    ! CHARACTER(LEN=16), ALLOCATABLE :: MESH_RESULTS(:)
    CHARACTER(LEN=16), ALLOCATABLE :: NODE_RESULTS(:)
    CHARACTER(LEN=16), ALLOCATABLE :: ELEM_RESULTS(:)
    !
    !----------------------------------------------------------------------
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    ! Allocate results string arrays with "headers"
    !
    ! ALLOCATE(MESH_RESULTS(1))
    ALLOCATE(NODE_RESULTS(1))
    ALLOCATE(ELEM_RESULTS(1))
    !
    ! MESH_RESULTS(1) = 'results_mesh'
    NODE_RESULTS(1) = 'results_nodes'
    ELEM_RESULTS(1) = 'results_elements'
    !
    ! Nodal coordinates
    !
    IF (PRINT_OPTIONS%PRINT_COORDINATES) THEN
        !
        IO = OUNITS(COORDS_U)
        WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, DOF_SUB1 + 1, DOF_SUP1 + 1
        WRITE(IO, '(3(e13.7,1x))') (CRD(J), J = DOF_SUB1, DOF_SUP1)
        !
        CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'coo', NODE_RESULTS, PFLAG_NODE)
        !
    END IF
    !
    ! Nodal displacements
    !
    IF (PRINT_OPTIONS%PRINT_DISPLACEMENTS) THEN
        !
        IO = OUNITS(DISP_U)
        WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, DOF_SUB1 + 1, DOF_SUP1 + 1
        WRITE(IO, '(3(e13.7,1x))') ((CRD(J) - COORDS_ORIG(J)), J = DOF_SUB1, &
            & DOF_SUP1)
        !
        CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'disp', NODE_RESULTS, PFLAG_NODE)
        !
    END IF
    !
    ! Nodal velocities
    !
    IF (PRINT_OPTIONS%PRINT_VELOCITIES) THEN
        !
        IO = OUNITS(VEL_U)
        WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, DOF_SUB1 + 1, DOF_SUP1 + 1
        WRITE(IO, '(3(e13.7,1x))') (VEL(J), J = DOF_SUB1, DOF_SUP1)
        !
        CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'vel', NODE_RESULTS, PFLAG_NODE)
        !
    END IF
    !
    ! Orientations
    !
    IF (PRINT_OPTIONS%PRINT_ORIENTATIONS) THEN
        !
        IO = OUNITS(ANGLES_U)
        !
        CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'ori', ELEM_RESULTS, PFLAG_ELEM)
        !
        ! Write header
        !
        WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
        !
        ! Determine passive (C2S) or active (S2C) from orientation options
        !
        ORIENT_INT = 0.0D0
        IF (ORIENTATION_OPTIONS%ORIENTATION_CONVENTION .EQ. &
            & 'passive') THEN
            !
            ! Don't do anything - FEPX assumes passive convention!!!
            !
            ORIENT_INT = ORIENT
            !
        ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_CONVENTION .EQ. &
            & 'active') THEN
            !
            DO J = EL_SUB1, EL_SUP1
                !
                DO IGRAIN = 0, NGRAIN1
                    !
                    ORIENT_INT(:, :, IGRAIN, J) = TRANSPOSE(ORIENT(:, :, &
                        & IGRAIN, J))
                    !
                END DO
                !
            END DO
            !
        END IF
        !
        ! Determine parameterization from orientation options
        !
        IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION .EQ. &
            & 'axis-angle') THEN
            !
            CALL ROT_MATS_TO_AXIS_ANGLE(NGRAIN, M_EL, ORIENT_INT, AXIS, ANGLE)
            !
            DO J = EL_SUB1, EL_SUP1
                !
                DO IGRAIN = 0, NGRAIN1
                    !
                    WRITE(IO, '(4f13.7)') AXIS(0, IGRAIN, J), &
                        & AXIS(1, IGRAIN, J), AXIS(2, IGRAIN, J), &
                        & ANGLE(IGRAIN, J)
                    !
                END DO
                !
            END DO
            !
        ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION .EQ. &
            & 'euler-bunge') THEN
            !
            CALL ROT_MATS_TO_EULER_BUNGE(NGRAIN, M_EL, ORIENT_INT, PSI1, PHI, &
                & PSI2)
            !
            DO J = EL_SUB1, EL_SUP1
                !
                DO IGRAIN = 0, NGRAIN1
                    !
                    WRITE(IO, '(4(e13.7,1x))') PSI1(IGRAIN, J), &
                        & PHI(IGRAIN, J), PSI2(IGRAIN, J)
                    !
                END DO
                !
            END DO
            !
        ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION .EQ. &
            & 'euler-kocks') THEN
            !
            CALL ROT_MATS_TO_EULER_KOCKS(NGRAIN, M_EL, ORIENT_INT, PSI, THE, &
                & PHI)
            !
            DO J = EL_SUB1, EL_SUP1
                !
                DO IGRAIN = 0, NGRAIN1
                    !
                    WRITE(IO, '(4(e13.7,1x))') PSI(IGRAIN, J), THE(IGRAIN, J), &
                        & PHI(IGRAIN, J)
                    !
                END DO
                !
            END DO
            !
        ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION .EQ. &
            & 'rodrigues') THEN
            !
            CALL ROT_MATS_TO_RODRIGUES(NGRAIN, M_EL, ORIENT_INT, RODS)
            !
            DO J = EL_SUB1, EL_SUP1
                !
                DO IGRAIN = 0, NGRAIN1
                    !
                    WRITE(IO, '(4(e13.7,1x))') (RODS(I, IGRAIN, J), I = 0, &
                        & DIMS1)
                    !
                END DO
                !
            END DO
            !
        ELSEIF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION .EQ. &
            & 'quaternion') THEN
            !
            CALL ROT_MATS_TO_QUATERNIONS(NGRAIN, M_EL, ORIENT_INT, QUAT)
            !
            DO J = EL_SUB1, EL_SUP1
                !
                DO IGRAIN = 0, NGRAIN1
                    !
                    WRITE(IO, '(4(e13.7,1x))') (QUAT(I, IGRAIN, J), I = 0, DIMS)
                    !
                END DO
                !
            END DO
            !
        END IF
        !
    END IF
    !
    ! CRSS Values
    !
    IF (PRINT_OPTIONS%PRINT_CRSS) THEN
        !
        IO=OUNITS(CRSS_U)
        !
        CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'crss', ELEM_RESULTS, PFLAG_ELEM)
        !
        WRITE(IO, '(A2,(3(I0,1X)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
        !
        DO J = EL_SUB1, EL_SUP1
            !
            DO IGRAIN = 0, NGRAIN1
                !
                ! First, check if the element is FCC/BCC.
                IF ((CRYS_OPTIONS%CRYSTAL_TYPE(PHASE(J)) .EQ. 1) .OR. &
                    & (CRYS_OPTIONS%CRYSTAL_TYPE(PHASE(J)) .EQ. 2)) THEN
                    !
                    ! Next, check if the hardening is isotropic or anisotropic.
                    IF ((OPTIONS%HARD_TYPE .EQ. 'isotropic') .OR. &
                        & (OPTIONS%HARD_TYPE .EQ. 'cyclic_isotropic')) THEN
                        !
                        WRITE(IO, '(1(E13.7,1X))') (HARD(0, IGRAIN, J))
                        !
                    ELSE IF (OPTIONS%HARD_TYPE .EQ. 'latent') THEN
                        !
                        WRITE(IO, '(12(E13.7,1X))') HARD(0:11, IGRAIN, J)
                        !
                    ELSE
                        !
                        CALL PAR_QUIT('Error  :     > Invalid hardening type &
                            &provided.')
                        !
                    END IF
                    !
                ! Else, check if the element is HCP.
                ELSE IF (CRYS_OPTIONS%CRYSTAL_TYPE(PHASE(J)) .EQ. 3) THEN
                    !
                    ! Next, check if the hardening is isotropic or anisotropic.
                    IF ((OPTIONS%HARD_TYPE .EQ. 'isotropic') .OR. &
                        & (OPTIONS%HARD_TYPE .EQ. 'cyclic_isotropic')) THEN
                        !
                        WRITE(IO, '(3(E13.7,1X))') (HARD(0, IGRAIN, J)), &
                            & CRYS_OPTIONS%PRISMATIC_TO_BASAL(PHASE(J)) &
                            & * (HARD(0, IGRAIN, J)), &
                            & CRYS_OPTIONS%PYRAMIDAL_TO_BASAL(PHASE(J)) &
                            & * (HARD(0, IGRAIN, J))
                        !
                    ELSE IF (OPTIONS%HARD_TYPE .EQ. 'latent') THEN
                        !
                        WRITE(IO, '(18(E13.7,1X))') &
                            & HARD(0:2, IGRAIN, J), &
                            & CRYS_OPTIONS%PRISMATIC_TO_BASAL(PHASE(J)) &
                            & * HARD(3:5, IGRAIN, J), &
                            & CRYS_OPTIONS%PYRAMIDAL_TO_BASAL(PHASE(J)) &
                            & * HARD(6:17, IGRAIN, J)
                        !
                    ELSE
                        !
                        CALL PAR_QUIT('Error  :     > Invalid hardening type &
                            &provided.')
                        !
                    END IF
                    !
                END IF
                !
            END DO
            !
        END DO
        !
    END IF
    !
    ! Element volume (scalar)
    !
    IF (PRINT_OPTIONS%PRINT_ELVOL) THEN
        !
        IO = OUNITS(ELVOL_U)
        !
        CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'elt-vol', ELEM_RESULTS, &
            & PFLAG_ELEM)
        !
        WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
        !
        DO I = EL_SUB1, EL_SUP1
            !
            WRITE(IO, '(e13.7)') ELVOL(I)
            !
        END DO
        !
    END IF
    !
    !
    IF (PROB_TYPE == ANISOTROPIC_EVPS) THEN
        !
        ! Elemental elastic strain tensors (6-vector)
        !
        IF (PRINT_OPTIONS%PRINT_STRAIN_EL) THEN
            !
            IO = OUNITS(ELSTRAIN_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'strain-el', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                DO IGRAIN = 0, NGRAIN1
                    !
                    ELAS_SAM = 0.0D0
                    CALL VEC6_CRYS_TO_VEC6_SAM(ELAS(:, IGRAIN, I), &
                        & ORIENT(:, :, IGRAIN, I), ELAS_SAM)
                    !
                    ! Written as 11, 22, 33, 23, 13, 12 in sample basis
                    !
                    WRITE(IO, '(6(e13.7,1x))') ELAS_SAM(0), ELAS_SAM(3), &
                        & ELAS_SAM(5), ELAS_SAM(4), &
                        & ELAS_SAM(2), ELAS_SAM(1)
                    !
                END DO
                !
            END DO
            !
        END IF
        !
        ! Elemental plastic strain tensors
        !
        IF (PRINT_OPTIONS%PRINT_STRAIN_PL) THEN
            !
            IO = OUNITS(PLSTRAIN_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'strain-pl', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                ! Convert from 5-vec to 3x3 matrix, then to 6-vector
                !
                PLSTRAINMAT = 0.0D0
                PLSTRAIN6 = 0.0D0
                CALL VEC_MAT_SYMM_SER(PLSTRAIN(:, I), PLSTRAINMAT)
                PLSTRAIN6 = (/ PLSTRAINMAT(0, 0), PLSTRAINMAT(0, 1), &
                    & PLSTRAINMAT(0, 2), PLSTRAINMAT(1, 1), PLSTRAINMAT(1, 2), &
                    & PLSTRAINMAT(2, 2) /)
                !
                ! Transform from crystal to sample basis
                !
                CALL VEC6_CRYS_TO_VEC6_SAM(PLSTRAIN6, &
                    & ORIENT(:, :, 0, I), PLSTRAIN6_SAM)
                !
                ! Written as 11, 22, 33, 23, 13, 12 in sample basis
                !
                WRITE(IO, '(6(e13.7,1x))') PLSTRAIN6_SAM(0), PLSTRAIN6_SAM(3), &
                    & PLSTRAIN6_SAM(5), PLSTRAIN6_SAM(4), &
                    & PLSTRAIN6_SAM(2), PLSTRAIN6_SAM(1)
                !
            END DO
            !
        END IF
        !
        ! Elemental total strain tensors
        !
        IF (PRINT_OPTIONS%PRINT_STRAIN_TOT) THEN
            !
            IO = OUNITS(TOTSTRAIN_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'strain', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                !
                ! Written as 11, 22, 33, 23, 13, 12 in sample basis
                !
                WRITE(IO, '(6(e13.7,1x))') TOTSTRAIN(0, 0, I), &
                    & TOTSTRAIN(1, 1 ,I), TOTSTRAIN(2, 2, I), &
                    & TOTSTRAIN(1, 2, I), TOTSTRAIN(0, 2, I), &
                    & TOTSTRAIN(0, 1, I)
                !
            END DO
            !
        END IF
        !
        ! Elemental stress tensors (3x3 matrix)
        !
        IF (PRINT_OPTIONS%PRINT_STRESS) THEN
            !
            IO = OUNITS(STRESS_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'stress', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                ! Written as 11, 22, 33, 23, 13, 12 in sample basis
                !
                WRITE(IO, '(6(e13.7,1x))') S3X3(0,0,I), S3X3(1,1,I), &
                    & S3X3(2,2,I), S3X3(1,2,I), S3X3(0,2,I), S3X3(0,1,I)
                !
            END DO
            !
        END IF
        !
        ! Elemental shear rates (gammadot)
        !
        IF (PRINT_OPTIONS%PRINT_GAMMADOT) THEN
            !
            IO = OUNITS(GAMMADOT_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'sliprate', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                DO IGRAIN = 0, NGRAIN1
                    !
                    WRITE(IO, '(18(e13.7,1x))') (GAMMADOT(J, IGRAIN, I), &
                        & J = 0, MAXSLIP1)
                    !
                END DO
                !
            END DO
            !
        END IF
        !
        ! Elemental equivalent stress (scalar)
        !
        IF (PRINT_OPTIONS%PRINT_EQSTRESS) THEN
            !
            CALL STRESS_EQUIV_3X3(S3X3, EQSTRESS)
            !
            IO = OUNITS(EQSTRESS_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'stress-eq', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                WRITE(IO, '(e13.7)') EQSTRESS(I)
                !
            END DO
            !
        END IF
        !
        !
        ! Elemental effective deformation rate (scalar)
        !
        IF (PRINT_OPTIONS%PRINT_DEFF) THEN
            !
            IO = OUNITS(DEFF_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'defrate-eq', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                WRITE(IO, '(e13.7)') DEFF(I)
                !
            END DO
            !
        END IF
        !
        ! Elemental equivalent strain (scalar)
        !
        IF (PRINT_OPTIONS%PRINT_EQSTRAIN) THEN
            !
            IO = OUNITS(EQSTRAIN_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'strain-eq', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                WRITE(IO, '(e13.7)') EQSTRAIN(I)
                !
            END DO
            !
        END IF
        !
        ! Elemental effective plastic deformation rate (scalar)
        !
        IF (PRINT_OPTIONS%PRINT_DPEFF) THEN
            !
            IO = OUNITS(DPEFF_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'defrate-pl-eq', &
                & ELEM_RESULTS, PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                WRITE(IO, '(e13.7)') DPEFF(I)
                !
            END DO
            !
        END IF
        !
        ! Elemental equivalent plastic strain (scalar)
        !
        IF (PRINT_OPTIONS%PRINT_EQPLSTRAIN) THEN
            !
            IO = OUNITS(EQPLSTRAIN_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'strain-pl-eq', &
                & ELEM_RESULTS, PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                WRITE(IO, '(e13.7)') EQPLSTRAIN(I)
                !
            END DO
            !
        END IF
        !
        ! Elemental equivalent plastic strain (scalar)
        !
        IF (PRINT_OPTIONS%PRINT_EQELSTRAIN) THEN
            !
            IO = OUNITS(EQELSTRAIN_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'strain-el-eq', &
                & ELEM_RESULTS, PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                WRITE(IO, '(e13.7)') EQELSTRAIN(I)
                !
            END DO
            !
        END IF
        !
        ! Elemental velocity gradient (3x3 matrix)
        !
        IF (PRINT_OPTIONS%PRINT_VGRAD) THEN
            !
            IO = OUNITS(VGRAD_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'velgrad', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                WRITE(IO, '(9(e13.7,1x))') VGRAD(0,0,I), VGRAD(0,1,I), &
                    & VGRAD(0,2,I), VGRAD(1,0,I), VGRAD(1,1,I), VGRAD(1,2,I), &
                    & VGRAD(2,0,I), VGRAD(2,1,I), VGRAD(2,2,I)
                !
            END DO
            !
        END IF
        !
        ! Elemental plastic deformation rate (6-vector)
        !
        IF (PRINT_OPTIONS%PRINT_DPHAT) THEN
            !
            IO = OUNITS(DPHAT_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'defrate-pl', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                DP_HAT_MAT = 0.0D0
                DP_HAT_SAM = 0.0D0
                !
                ! Convert the transpose of the elemental rotation matrix into
                ! an operator on 5-vectors. The transpose outputs QR5X5 in
                ! a crystal-to-sample form.
                !
                CALL ROT_MAT_SYMM_SER(TRANSPOSE(ORIENT(:, :, 0, I)), QR5X5)
                !
                ! LATTICE_DEFORM transforms DP_HAT from crystal to sample frame.
                ! This is because ORIENT was transposed (above).
                ! 5-vector order is maintained here as (11-22), 33, 12, 13, 23
                ! with the proper scalings.
                !
                CALL LATTICE_DEFORM_SER(QR5X5, DP_HAT(:, I), DP_HAT_SAM)
                !
                ! Convert 5-vector with scaling into proper 6-vector form.
                !
                CALL VEC5_VEC6(DP_HAT_SAM, DP_HAT_MAT)
                !
                ! Written as 11, 22, 33, 23, 13, 12 in sample basis.
                !
                WRITE(IO, '(6(e13.7,1x))') DP_HAT_MAT(0), DP_HAT_MAT(1), &
                    & DP_HAT_MAT(2), DP_HAT_MAT(3), &
                    & DP_HAT_MAT(4), DP_HAT_MAT(5)
                !
            END DO
            !
        END IF
        !
        ! Elemental plastic spin rate (3-vector)
        !
        IF (PRINT_OPTIONS%PRINT_WPHAT) THEN
            !
            IO = OUNITS(WPHAT_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'spinrate', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                ! Convert the transpose of the elemental rotation matrix into
                ! an operator on skew 3-vectors. The transpose outputs QR3X3 in
                ! a crystal-to-sample form.
                !
                CALL ROT_MAT_SKEW_SER(TRANSPOSE(ORIENT(:, :, 0, I)), QR3X3)
                !
                ! LATTICE_SPIN transforms WP_HAT from crystal to sample frame.
                ! This is because ORIENT was transposed (above).
                ! Skew 3-vector order is maintained here as 21 31 32.
                !
                CALL LATTICE_SPIN_SER(QR3X3, WP_HAT(:, I), WP_HAT_SAM)
                !
                ! Negative values are output here in order to obtain desired
                ! order while maintaining skew-symmetry constraints.
                !
                ! Written as 12, 13, 23 in sample basis.
                !
                WRITE(IO, '(3(e13.7,1x))') - WP_HAT_SAM(0), - WP_HAT_SAM(1), &
                    & - WP_HAT_SAM(2)
                !
            END DO
            !
        END IF
        !
        ! Elemental accumulated shears
        !
        IF (PRINT_OPTIONS%PRINT_GAMMA) THEN
            !
            IO = OUNITS(GAMMA_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'slip', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                DO IGRAIN = 0, NGRAIN1
                    !
                    WRITE(IO, '(18(e13.7,1x))') (GAMMA(J, IGRAIN, I), &
                        & J = 0, MAXSLIP1)
                    !
                END DO
                !
            END DO
            !
        END IF
        !
        ! Elemental total work (scalar)
        !
        IF (PRINT_OPTIONS%PRINT_WORK) THEN
            !
            IO = OUNITS(WORK_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'work', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                WRITE(IO, '(e13.7)') WORK(I)
                !
            END DO
            !
        END IF
        !
        ! Elemental plastic work (scalar)
        !
        IF (PRINT_OPTIONS%PRINT_WORKP) THEN
            !
            IO = OUNITS(WORKP_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'work-pl', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                WRITE(IO, '(e13.7)') WORKP(I)
                !
            END DO
            !
        END IF
        !
        ! Elemental total deformation rate tensor (3x3 matrix)
        !
        IF (PRINT_OPTIONS%PRINT_DEFRATE) THEN
            !
            IO = OUNITS(DEFRATE_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'defrate', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                ! Written as 11, 22, 33, 23, 13, 12 in sample basis
                !
                WRITE(IO, '(9(e13.7,1x))') D_TOT(0,0,I), D_TOT(1,1,I), &
                    & D_TOT(2,2,I), D_TOT(1,2,I), D_TOT(0,2,I), D_TOT(0,1,I)
                !
            END DO
            !
        END IF
        !
        ! Elemental total work rate (scalar)
        !
        IF (PRINT_OPTIONS%PRINT_WORKRATE) THEN
            !
            IO = OUNITS(WORKRATE_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'workrate', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                WRITE(IO, '(e13.7)') WORKRATE(I)
                !
            END DO
            !
        END IF
        !
        ! Elemental plastic work rate (scalar)
        !
        IF (PRINT_OPTIONS%PRINT_WORKRATEP) THEN
            !
            IO = OUNITS(WORKRATEP_U)
            !
            CALL ADD_TO_REPORT_OUTPUT_FILES(STEP, 'workrate-pl', ELEM_RESULTS, &
                & PFLAG_ELEM)
            !
            WRITE(IO, '(a2,(3(i0,1x)))') '% ', STEP, EL_SUB1 + 1, EL_SUP1 + 1
            !
            DO I = EL_SUB1, EL_SUP1
                !
                WRITE(IO, '(e13.7)') WORKRATEP(I)
                !
            END DO
            !
        END IF
        !
        ! Write output files to the post.result files
        !
        CALL WRITE_REPORT_FILE_OUTPUT_FILES(STEP, NODE_RESULTS, PFLAG_NODE)
        CALL WRITE_REPORT_FILE_OUTPUT_FILES(STEP, ELEM_RESULTS, PFLAG_ELEM)
        PFLAG_NODE = .FALSE.
        PFLAG_ELEM = .FALSE.
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE PRINT_STEP
    !
    !===========================================================================
    !
    SUBROUTINE WRITE_REPORT_FILE_HEADER(PART_ARRAY)
    !
    ! Write a post.report file containing required information to facilitate
    ! automated post-processing via Neper.
    !
    ! This prints the number of nodes, elements, partitions, elements by part,
    ! nodes by part, number of slip systems, and orientation definition.
    !
    ! Arugments:
    ! PART_ARRAY: Gathered array from fepx.f90 that contains partition info.
    !
    INTEGER, DIMENSION(0:1,0:NUMPROCS-1), INTENT(IN) :: PART_ARRAY
    !
    ! Locals:
    ! I: Looping index
    !
    INTEGER :: I
    !---------------------------------------------------------------------------
    !
    ! The file is opened in fepx.f90 so just begin printing from master proc.
    ! All printed data should be public from the top-level so no arguments.
    IF (MYID .EQ. 0) THEN
        !
        ! Print the number of nodes, elements, and partitions
        !
        WRITE(OUNITS(REPORT_U), '(A,I0)') 'number_of_nodes ', NUMNP
        WRITE(OUNITS(REPORT_U), '(A,I0)') 'number_of_elements ', NUMELM
        WRITE(OUNITS(REPORT_U), '(A,I0)') 'number_of_partitions ', NUMPROCS
        !
        ! Print number of elements on each partition
        !
        WRITE(OUNITS(REPORT_U), '(A)', ADVANCE='NO') &
            & 'number_of_elements_bypartition '
        !
        DO I = 0, NUMPROCS-1
            !
            WRITE(OUNITS(REPORT_U), '(I0, A)', ADVANCE='NO') PART_ARRAY(0,I), &
                & ' '
            !
        END DO
        !
        ! Print number of nodes on each partition
        !
        WRITE(OUNITS(REPORT_U), '(/,A)', ADVANCE='NO') &
            & 'number_of_nodes_bypartition '
        !
        DO I = 0, NUMPROCS-1
            !
            WRITE(OUNITS(REPORT_U), '(I0, A)', ADVANCE='NO') PART_ARRAY(1,I), &
                & ' '
            !
        END DO
        !
        ! Print number of phases
        !
        WRITE(OUNITS(REPORT_U), '(/,A,I0)', ADVANCE='NO') 'number_of_phases ', &
            & NUMPHASES
        !
        ! Print number of slip systems and the orientation definition
        !
        WRITE(OUNITS(REPORT_U), '(/,A)', ADVANCE='NO') 'number_of_slip_systems '
        !
        DO I = 1, NUMPHASES
            !
            WRITE(OUNITS(REPORT_U), '(I0, A)', ADVANCE='NO') CTYPE(I)%NUMSLIP, &
                & ' '
            !
        END DO
        !
        WRITE(OUNITS(REPORT_U), '(/,A,A,A,A)') 'orientation_definition ', &
            & TRIM(ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION), ':', &
            & TRIM(ORIENTATION_OPTIONS%ORIENTATION_CONVENTION)
        !
    END IF
    !
    END SUBROUTINE WRITE_REPORT_FILE_HEADER
    !
    !===========================================================================
    !
    SUBROUTINE WRITE_REPORT_FILE_OUTPUT_FILES(CURRENT_STEP, INPUT_STRING, PFLAG)
    !
    ! This writes the file names of the user-defined output files from the
    ! simulation.config.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! INPUT_STRING: Character array denoting file name to print.
    ! CURRENT_STEP: Current step at which we are printing information for.
    ! PFLAG: Determine if printing should happen or not
    !
    CHARACTER(LEN=*), INTENT(IN) :: INPUT_STRING(:)
    INTEGER, INTENT(IN) :: CURRENT_STEP
    LOGICAL, INTENT(IN) :: PFLAG
    !
    ! Locals:
    ! I: Looping variable for printing of INPUT_STRING
    !
    INTEGER :: I
    !
    !---------------------------------------------------------------------------
    !
    ! Only print to the post.report file from the root processor
    !
    IF (MYID .EQ. 0) THEN
        !
        ! If the first step is being printed then print to the report file.
        ! If the "first" step after restarting a simulation, also print.
        IF ((PFLAG .EQV. .TRUE.) .OR. ((OPTIONS%RESTART .EQV. .TRUE.) .AND. &
            & (CURRENT_STEP .EQ. OPTIONS%RESTART_INITIAL_STEP))) THEN
            !
            WRITE(OUNITS(REPORT_U), '(25(A,1X))') &
                & (TRIM(INPUT_STRING(I)),I = 1,SIZE(INPUT_STRING))
            !
        END IF
        !
    END IF
    !
    END SUBROUTINE WRITE_REPORT_FILE_OUTPUT_FILES
    !
    !===========================================================================
    !
    SUBROUTINE ADD_TO_REPORT_OUTPUT_FILES(CURRENT_STEP, INPUT_STRING, &
        &STRING_ARRAY, PFLAG)
    !
    ! This writes the file names of the user-defined output files from the
    ! simulation.config.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! INPUT_STRING: Character array denoting file name to print.
    ! CURRENT_STEP: Current step at which we are printing information for.
    !
    CHARACTER(LEN=*), INTENT(IN) :: INPUT_STRING
    CHARACTER(LEN=16), INTENT(INOUT), ALLOCATABLE :: STRING_ARRAY(:)
    INTEGER, INTENT(IN) :: CURRENT_STEP
    LOGICAL, INTENT(IN) :: PFLAG
    !
    CHARACTER(LEN=16), ALLOCATABLE :: TEMP_ARRAY(:)
    INTEGER :: ARRAYSIZE, I
    !
    !---------------------------------------------------------------------------
    !
    ! Only print to the post.report file from the root processor
    !
    IF (MYID .EQ. 0) THEN
        !
        ! If the first step is being printed then print to the report file.
        ! If the "first" step after restarting a simulation, also print.
        IF ((PFLAG .EQV. .TRUE.) .OR. ((OPTIONS%RESTART .EQV. .TRUE.) .AND. &
            & (CURRENT_STEP .EQ. OPTIONS%RESTART_INITIAL_STEP))) THEN
            !
            ! Append the INPUT_STRING to the present array and reallocate.
            IF(ALLOCATED(STRING_ARRAY)) THEN
                !
                ARRAYSIZE = SIZE(STRING_ARRAY)
                ALLOCATE(TEMP_ARRAY(ARRAYSIZE+1))
                !
                DO I = 1, ARRAYSIZE
                    !
                    TEMP_ARRAY(I) = STRING_ARRAY(I)
                    !
                END DO
                !
                TEMP_ARRAY(ARRAYSIZE+1) = INPUT_STRING
                DEALLOCATE(STRING_ARRAY)
                CALL MOVE_ALLOC(TEMP_ARRAY, STRING_ARRAY)
                !
            END IF
            !
        END IF
        !
    END IF
    !
    END SUBROUTINE ADD_TO_REPORT_OUTPUT_FILES
    !
    !===========================================================================
    !
    SUBROUTINE WRITE_REPORT_FILE_COMPLETE_STEPS(FINAL_STEP, PRINT_ARRAY)
    !
    ! This writes the completed number of steps from the driver as the current
    ! simulation is either finishing successfully or terminated early.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! FINAL_STEP: Final completed step in this simulation
    ! PRINT_ARRAY: Array of integers denoting which steps are to be printed
    !
    INTEGER, INTENT(IN) :: FINAL_STEP
    INTEGER, INTENT(IN) :: PRINT_ARRAY(:)
    !
    ! Locals:
    ! I: Generic looping index
    !
    INTEGER :: I
    !
    !---------------------------------------------------------------------------
    !
    ! If the first step is being printed then print to the report file.
    IF (MYID .EQ. 0) THEN
        !
        ! Print the completed number of steps
        !
        WRITE(OUNITS(REPORT_U), '(A,I0)') 'number_of_steps ', &
            & FINAL_STEP
        !
        ! Print the IDs of the steps that were requested to be printed
        !
        WRITE(OUNITS(REPORT_U), '(A)', ADVANCE='NO') 'printed_steps '
        !
        DO I = 1, SIZE(PRINT_ARRAY,1)
            !
            IF (PRINT_ARRAY(I) .EQ. 0) WRITE(OUNITS(REPORT_U), &
                & '(I0, A)', ADVANCE='NO') I, ' '
            !
        END DO
        !
        ! Close the report file before ending the process
        !
        CLOSE(OUNITS(REPORT_U))
        !
    END IF
    !
    END SUBROUTINE WRITE_REPORT_FILE_COMPLETE_STEPS
    !
    !===========================================================================
    !
    SUBROUTINE WRITE_FORCE_FILE_HEADERS(DRIVER_TYPE)
    !
    ! This writes the file headers for surface force files. Consistent format
    ! across all drivers is maintained here. Assumes that the CALL is wrapped
    ! in `IF (PRINT_OPTIONS%PRINT_FORCES) THEN` logic from where it is called.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! DRIVER_TYPE: Flag that denotes if the calling driver is uni- or triaxial
    !   (1) = uniaxial, (2) = triaxial
    !
    INTEGER :: DRIVER_TYPE
    !
    ! Locals:
    ! FORCE_HEADERU1/2: Headers for the uniaxial drivers
    ! FORCE_HEADERT1/2: Headers for the triaxial drivers
    !
    CHARACTER(LEN=*), PARAMETER :: FORCE_HEADERU1 = &
        &   '% step     INCR     Fx             Fy            &
        & Fz             area A         TIME'
    CHARACTER(LEN=*), PARAMETER :: FORCE_HEADERU2 = &
        &   '% -------- -------- -------------- --------------&
        & -------------- -------------- --------------'
    !
    CHARACTER(LEN=*), PARAMETER :: FORCE_HEADERT1 = &
        &   '% step     incr     Fx             Fy            &
        & Fz             area A         time           length'
    CHARACTER(LEN=*), PARAMETER :: FORCE_HEADERT2 = &
        &   '% -------- -------- -------------- --------------&
        & -------------- -------------- -------------- -----------------'
    !
    !---------------------------------------------------------------------------
    !
    ! If restartinbg with appending then we don't need the headers
    !
    IF (OPTIONS%RESTART_FILE_HANDLING .EQ. 0) THEN
        !
        RETURN
        !
    END IF
    !
    ! Check the input DRIVER_TYPE and print the correct headers
    !
    IF (DRIVER_TYPE .EQ. 1) THEN ! uniaxial
        !
        ! Print the headers to the files
        WRITE(OUNITS(FORCE_U1),'(a)') FORCE_HEADERU1
        WRITE(OUNITS(FORCE_U1),'(a)') FORCE_HEADERU2
        WRITE(OUNITS(FORCE_U2),'(a)') FORCE_HEADERU1
        WRITE(OUNITS(FORCE_U2),'(a)') FORCE_HEADERU2
        WRITE(OUNITS(FORCE_U3),'(a)') FORCE_HEADERU1
        WRITE(OUNITS(FORCE_U3),'(a)') FORCE_HEADERU2
        WRITE(OUNITS(FORCE_U4),'(a)') FORCE_HEADERU1
        WRITE(OUNITS(FORCE_U4),'(a)') FORCE_HEADERU2
        WRITE(OUNITS(FORCE_U5),'(a)') FORCE_HEADERU1
        WRITE(OUNITS(FORCE_U5),'(a)') FORCE_HEADERU2
        WRITE(OUNITS(FORCE_U6),'(a)') FORCE_HEADERU1
        WRITE(OUNITS(FORCE_U6),'(a)') FORCE_HEADERU2
        !
    ELSE IF (DRIVER_TYPE .EQ. 2) THEN ! triaxial
        !
        ! Print the headers to the files
        WRITE(OUNITS(FORCE_U1),'(a)') FORCE_HEADERT1
        WRITE(OUNITS(FORCE_U1),'(a)') FORCE_HEADERT2
        WRITE(OUNITS(FORCE_U2),'(a)') FORCE_HEADERT1
        WRITE(OUNITS(FORCE_U2),'(a)') FORCE_HEADERT2
        WRITE(OUNITS(FORCE_U3),'(a)') FORCE_HEADERT1
        WRITE(OUNITS(FORCE_U3),'(a)') FORCE_HEADERT2
        WRITE(OUNITS(FORCE_U4),'(a)') FORCE_HEADERT1
        WRITE(OUNITS(FORCE_U4),'(a)') FORCE_HEADERT2
        WRITE(OUNITS(FORCE_U5),'(a)') FORCE_HEADERT1
        WRITE(OUNITS(FORCE_U5),'(a)') FORCE_HEADERT2
        WRITE(OUNITS(FORCE_U6),'(a)') FORCE_HEADERT1
        WRITE(OUNITS(FORCE_U6),'(a)') FORCE_HEADERT2
        !
    END IF
    !
    END SUBROUTINE WRITE_FORCE_FILE_HEADERS
    !
    !===========================================================================
    !
    SUBROUTINE WRITE_FORCE_FILE_DATA(DRIVER_TYPE, ISTEP, INCR, LOAD_ARRAY, &
        & AREA, TIME, LENGTH)
    !
    ! This writes the actual data for surface force files. Consistent format
    ! across all drivers is maintained here. Assumes that the CALL is wrapped
    ! in `IF (PRINT_OPTIONS%PRINT_FORCES) THEN` logic from where it is called.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! DRIVER_TYPE: Flag that denotes if the calling driver is uni- or triaxial
    !   (1) = uniaxial, (2) = triaxial
    ! ISTEP: Current timestep being printed
    ! INCR: Current total increment being printed
    ! LOAD_ARRAY: Contains the [X Y Z] loads on all surfaces
    ! AREA: Current surface face areas bring printed
    ! TIME: Current time value being printed
    ! LENGTH: Current length of the mesh edges (triaxial only)
    !
    INTEGER :: DRIVER_TYPE
    INTEGER :: ISTEP
    INTEGER :: INCR
    REAL(RK) :: LOAD_ARRAY(:,:)
    REAL(RK) :: AREA(:)
    REAL(RK) :: TIME
    REAL(RK) :: LENGTH(:)
    !
    !---------------------------------------------------------------------------
    !
    ! Check the input DRIVER_TYPE and print the correct headers
    IF (DRIVER_TYPE .EQ. 1) THEN ! uniaxial
        !
        ! Print the data to the files
        WRITE(OUNITS(FORCE_U1), '(2(I8),5(E15.5))') ISTEP, INCR, &
            & LOAD_ARRAY(1,:), AREA(1), TIME
        WRITE(OUNITS(FORCE_U2), '(2(I8),5(E15.5))') ISTEP, INCR, &
            & LOAD_ARRAY(2,:), AREA(2), TIME
        WRITE(OUNITS(FORCE_U3), '(2(I8),5(E15.5))') ISTEP, INCR, &
            & LOAD_ARRAY(3,:), AREA(3), TIME
        WRITE(OUNITS(FORCE_U4), '(2(I8),5(E15.5))') ISTEP, INCR, &
            & LOAD_ARRAY(4,:), AREA(4), TIME
        WRITE(OUNITS(FORCE_U5), '(2(I8),5(E15.5))') ISTEP, INCR, &
            & LOAD_ARRAY(5,:), AREA(5), TIME
        WRITE(OUNITS(FORCE_U6), '(2(I8),5(E15.5))') ISTEP, INCR, &
            & LOAD_ARRAY(6,:), AREA(6), TIME
        !
    ELSE IF (DRIVER_TYPE .EQ. 2) THEN ! triaxial
        !
        ! Print the data to the files
        WRITE(OUNITS(FORCE_U1), '(2(I8),5(E15.5),E18.8)') ISTEP, INCR, &
            & LOAD_ARRAY(1,:), AREA(1), TIME, LENGTH(1)
        WRITE(OUNITS(FORCE_U2), '(2(I8),5(E15.5),E18.8)') ISTEP, INCR, &
            & LOAD_ARRAY(2,:), AREA(2), TIME, LENGTH(1)
        WRITE(OUNITS(FORCE_U3), '(2(I8),5(E15.5),E18.8)') ISTEP, INCR, &
            & LOAD_ARRAY(3,:), AREA(3), TIME, LENGTH(2)
        WRITE(OUNITS(FORCE_U4), '(2(I8),5(E15.5),E18.8)') ISTEP, INCR, &
            & LOAD_ARRAY(4,:), AREA(4), TIME, LENGTH(2)
        WRITE(OUNITS(FORCE_U5), '(2(I8),5(E15.5),E18.8)') ISTEP, INCR, &
            & LOAD_ARRAY(5,:), AREA(5), TIME, LENGTH(3)
        WRITE(OUNITS(FORCE_U6), '(2(I8),5(E15.5),E18.8)') ISTEP, INCR, &
            & LOAD_ARRAY(6,:), AREA(6), TIME, LENGTH(3)
        !
    END IF
    !
    END SUBROUTINE WRITE_FORCE_FILE_DATA
    !
    !===========================================================================
    !
    SUBROUTINE WRITE_CONV_FILE_HEADERS()
    !
    ! This writes the file headers for convergence report. Consistent format
    ! across all drivers is maintained here. Assumes that the CALL is wrapped
    ! in `IF (PRINT_OPTIONS%PRINT_CONV) THEN` logic from where it is called.
    !
    !---------------------------------------------------------------------------
    !
    ! Print the headers
    WRITE(OUNITS(CONV_U),'(a)') '%   INCR     iter       NR     r_norm&
        &        rx_norm       f_norm        delu_norm     delux_norm&
        &    u_norm      cg_iter'
    !
    END SUBROUTINE WRITE_CONV_FILE_HEADERS
    !
    !===========================================================================
    !
    SUBROUTINE WRITE_CONV_FILE_DATA(INCR, ITER, ITMETHOD, R_NORM, RX_NORM, &
        & F_NORM, DELU_NORM, DELUX_NORM, U_NORM, CG_ITER_OUT)
    !
    ! This writes the file data for convergence report. This is only called
    ! from the ITMETHOD_EVPS subroutine. Assumes that the CALL is wrapped in
    ! `IF (PRINT_OPTIONS%PRINT_CONV) THEN` logic from where it is called.
    !
    ! Arguments:
    ! INCR: Current total increment value
    ! ITER: Current iteration within step
    ! ITMETHOD: Flag denoting successive approx. (0) or newton-raphson (1)
    ! R_NORM: Residual norm -> sqrt(sum(resid * resid))
    ! RX_NORM: Maximal absolute residual -> maxval(abs(resid))
    ! F_NORM: Force norm -> sqrt(sum(force * force))
    ! DELU_NORM: Change in velocity norm -> sqrt(sum(delta_vel * delta_vel))
    ! DELUX_NORM: Maximal absolute change in velocity -> maxval(abs(delta_vel))
    ! U_NORM: Velocity norm -> sqrt(sum(vel_o * vel_o))
    ! CG_ITER_OUT: Number of conjugate-gradient (CG) iterations for this INCR
    !
    INTEGER  :: INCR, ITER, ITMETHOD, CG_ITER_OUT
    REAL(RK) :: R_NORM, RX_NORM, F_NORM, DELU_NORM, DELUX_NORM, U_NORM
    !---------------------------------------------------------------------------
    !
    ! Print the data to the files
    WRITE(OUNITS(CONV_U),'(I8,1X,I8,1X,I8,1X,6D14.4,1X,I8)') INCR, ITER, &
        & ITMETHOD, R_NORM, RX_NORM, F_NORM, DELU_NORM, DELUX_NORM, U_NORM, &
        & CG_ITER_OUT
    !
    END SUBROUTINE WRITE_CONV_FILE_DATA
    !
    !===========================================================================
    !
    SUBROUTINE WRITE_RESTART_FIELD(VELOCITY, C0_ANGS, C_ANGS, &
        & RSTAR, RSTAR_N, WTS, CRSS, CRSS_N, &
        & E_ELAS_KK_BAR, SIG_VEC_N, EQSTRAIN, EQPLSTRAIN, GAMMA, &
        & EL_WORK_N, EL_WORKP_N, EL_WORK_RATE_N, EL_WORKP_RATE_N, PLSTRAIN, &
        & TOTSTRAIN, ISTEP)
    !
    !  Writes field data for restarting simulation
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    REAL(RK), INTENT(IN) :: VELOCITY(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(IN) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: E_ELAS_KK_BAR(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: SIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: EQSTRAIN(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: EQPLSTRAIN(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: GAMMA(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: EL_WORK_N(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: EL_WORKP_N(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: EL_WORK_RATE_N(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: EL_WORKP_RATE_N(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: PLSTRAIN(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: TOTSTRAIN(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    INTEGER, INTENT(IN) :: ISTEP
    !
    ! Locals:
    !
    INTEGER :: MYUNIT
    CHARACTER(LEN=8) :: CHARID ! assumes less than 10,000 processes
    INTEGER :: RST_NUM
    LOGICAL :: FILE_EXISTS
    CHARACTER(LEN=8) :: RST_NUM_STR
    CHARACTER(LEN=64) :: FILENAME
    !
    INTRINSIC :: TRIM
    !
    !---------------------------------------------------------------------------
    !
    MYUNIT = NEWUNITNUMBER()
    WRITE(CHARID, '(I0)') MYID + 1
    !
    RST_NUM = 1000
    FILE_EXISTS = .FALSE.
    !
    DO WHILE (.NOT. FILE_EXISTS)
        !
        WRITE(RST_NUM_STR, '(I0)') RST_NUM
        !
        FILENAME = 'rst'//TRIM(RST_NUM_STR)//'.control'
        !
        INQUIRE(FILE = FILENAME, EXIST = FILE_EXISTS)
        RST_NUM = RST_NUM - 1
        !
        ! If we are on the first simulation, non-restart
        !
        IF (RST_NUM .EQ. -2) THEN
            !
            FILE_EXISTS = .TRUE.
            RST_NUM = -2
            !
        END IF
        !
    END DO
    !
    IF (ISTEP .EQ. 1) THEN
        !
        RST_NUM = RST_NUM + 2
        WRITE(RST_NUM_STR, '(I0)') RST_NUM
        !
    END IF
    !
    FILENAME = 'rst'//TRIM(RST_NUM_STR)//'.field.core'//TRIM(CHARID)
    OPEN(UNIT = MYUNIT, FILE = FILENAME, FORM = 'UNFORMATTED', ACTION = 'WRITE')
    !
    ! Velocity and coordinates.
    !
    WRITE(MYUNIT) COORDS
    WRITE(MYUNIT) VELOCITY
    !
    ! Orientations, weights and hardnesses.
    !
    WRITE(MYUNIT) C0_ANGS
    WRITE(MYUNIT) C_ANGS
    WRITE(MYUNIT) RSTAR
    WRITE(MYUNIT) RSTAR_N
    WRITE(MYUNIT) WTS
    WRITE(MYUNIT) CRSS
    WRITE(MYUNIT) CRSS_N
    !
    ! Elastic Strains.
    !
    WRITE(MYUNIT) GELA_KK_BAR
    WRITE(MYUNIT) GSIG_VEC_N
    WRITE(MYUNIT) PELA_KK_BAR
    WRITE(MYUNIT) PSIG_VEC_N
    WRITE(MYUNIT) E_ELAS_KK_BAR
    WRITE(MYUNIT) SIG_VEC_N
    !
    ! Equivalent Strains
    !
    WRITE(MYUNIT) EQSTRAIN
    WRITE(MYUNIT) EQPLSTRAIN
    WRITE(MYUNIT) GAMMA
    !
    ! Total work, plastic work, and rates.
    !
    WRITE(MYUNIT) EL_WORK_N
    WRITE(MYUNIT) EL_WORKP_N
    WRITE(MYUNIT) EL_WORK_RATE_N
    WRITE(MYUNIT) EL_WORKP_RATE_N
    !
    ! Other integrated quantities
    !
    WRITE(MYUNIT) PLSTRAIN
    WRITE(MYUNIT) TOTSTRAIN
    !
    CLOSE(MYUNIT)
    !
    END SUBROUTINE WRITE_RESTART_FIELD
    !
    !===========================================================================
    !
    SUBROUTINE WRITE_UNIAXIAL_RESTART(INCR, TIME, LOAD, AREA, &
        & AREA0, CURRENT_STEP, PREVIOUS_LOAD, STEP_COMPLETE, DTIME_OLD, &
        & PREV_STRAIN, CURR_STRAIN, ISTEP)
    !
    ! WRITE uniaxial control restart information.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! INCR: Current step increment (additive).
    ! TIME: Current step total time.
    ! LOAD: Current load on all surfaces in the mesh.
    ! PREV_LOAD: Previous step load value.
    ! AREA/AREA0: Current and initial surface areas of the mesh at current step.
    !
    INTEGER,  INTENT(IN) :: INCR
    REAL(RK), INTENT(IN) :: TIME
    REAL(RK), INTENT(IN) :: LOAD(NSURFACES,3)
    REAL(RK), INTENT(IN) :: AREA(NSURFACES)
    REAL(RK), INTENT(IN) :: AREA0(NSURFACES)
    INTEGER,  INTENT(IN) :: CURRENT_STEP
    REAL(RK), INTENT(IN) :: PREVIOUS_LOAD(:)
    LOGICAL,  INTENT(IN) :: STEP_COMPLETE
    REAL(RK), INTENT(IN) :: DTIME_OLD
    REAL(RK), INTENT(IN) :: PREV_STRAIN
    REAL(RK), INTENT(IN) :: CURR_STRAIN
    INTEGER, INTENT(IN) :: ISTEP
    !
    ! Locals:
    ! MYUNIT: Current unit number to open restart file.
    ! ISURF: Generic loop index to loop over mesh surfaces.
    !
    INTEGER :: MYUNIT
    INTEGER :: RST_NUM
    LOGICAL :: FILE_EXISTS
    CHARACTER(LEN=8) :: RST_NUM_STR
    CHARACTER(LEN=64) :: FILENAME
    INTEGER :: ISURF
    !
    !---------------------------------------------------------------------------
    !
    IF (MYID .EQ. 0) THEN
        !
        MYUNIT = NEWUNITNUMBER()
        !
        RST_NUM = 1000
        FILE_EXISTS = .FALSE.
        !
        DO WHILE (.NOT. FILE_EXISTS)
            !
            WRITE(RST_NUM_STR, '(I0)') RST_NUM
            !
            FILENAME = 'rst'//TRIM(RST_NUM_STR)//'.control'
            !
            INQUIRE(FILE = FILENAME, EXIST = FILE_EXISTS)
            RST_NUM = RST_NUM - 1
            !
            ! If we are on the first simulation, non-restart
            !
            IF (RST_NUM .EQ. -2) THEN
                !
                FILE_EXISTS = .TRUE.
                RST_NUM = -2
                !
            END IF
            !
        END DO
        !
        IF (ISTEP .EQ. 1) THEN
            !
            RST_NUM = RST_NUM + 2
            WRITE(RST_NUM_STR, '(I0)') RST_NUM
            !
        END IF
        !
        FILENAME = 'rst'//TRIM(RST_NUM_STR)//'.control'
        OPEN(UNIT = MYUNIT, FILE = FILENAME, FORM = 'UNFORMATTED', &
            & ACTION = 'WRITE')
        !
        !
        WRITE(MYUNIT) CURRENT_STEP
        WRITE(MYUNIT) PREVIOUS_LOAD
        WRITE(MYUNIT) STEP_COMPLETE
        WRITE(MYUNIT) DTIME_OLD
        WRITE(MYUNIT) INCR
        WRITE(MYUNIT) TIME
        !
        DO ISURF = 1, NSURFACES
            !
            WRITE(MYUNIT) LOAD(ISURF,:)
            !
        END DO
        !
        WRITE(MYUNIT) AREA
        WRITE(MYUNIT) AREA0
        WRITE(MYUNIT) PREV_STRAIN
        WRITE(MYUNIT) CURR_STRAIN
        !
        CLOSE(MYUNIT)
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE WRITE_UNIAXIAL_RESTART
    !
    !===========================================================================
    !
    SUBROUTINE WRITE_TRIAXCSR_RESTART(ISTEP,CURR_LOAD,PREV_LOAD,STEP_COMPLETE, &
        & DTIME,INCR,TIME,SURF_LOAD_ARRAY,AREA,AREA0,LENGTH,LENGTH0,CURR_VEL, &
        & S_PERT_MAG, T_PERT_MAG)
    !
    ! Write uniaxial control restart information.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! Needs to be defined - JC
    !
    LOGICAL, INTENT(IN)  :: STEP_COMPLETE
    INTEGER, INTENT(IN)  :: ISTEP
    INTEGER, INTENT(IN)  :: INCR
    REAL(RK), INTENT(IN) :: CURR_LOAD(3), PREV_LOAD(3)
    REAL(RK), INTENT(IN) :: DTIME, TIME
    REAL(RK), INTENT(IN) :: SURF_LOAD_ARRAY(NSURFACES,3)
    REAL(RK), INTENT(IN) :: AREA(NSURFACES), AREA0(NSURFACES)
    REAL(RK), INTENT(IN) :: LENGTH(3), LENGTH0(3)
    REAL(RK), INTENT(IN) :: CURR_VEL(3)
    REAL(RK), INTENT(IN) :: S_PERT_MAG, T_PERT_MAG
    !
    ! Locals:
    ! MYUNIT: Current unit number to open restart file.
    ! ISURF: Generic loop index to loop over mesh surfaces.
    !
    INTEGER :: MYUNIT
    INTEGER :: RST_NUM
    LOGICAL :: FILE_EXISTS
    CHARACTER(LEN=8) :: RST_NUM_STR
    CHARACTER(LEN=64) :: FILENAME
    INTEGER :: ISURF
    !
    !---------------------------------------------------------------------------
    !
    IF (MYID .EQ. 0) THEN
        !
        MYUNIT = NEWUNITNUMBER()
        !
        RST_NUM = 1000
        FILE_EXISTS = .FALSE.
        !
        DO WHILE (.NOT. FILE_EXISTS)
            !
            WRITE(RST_NUM_STR, '(I0)') RST_NUM
            !
            FILENAME = 'rst'//TRIM(RST_NUM_STR)//'.control'
            !
            INQUIRE(FILE = FILENAME, EXIST = FILE_EXISTS)
            RST_NUM = RST_NUM - 1
            !
            ! If we are on the first simulation, non-restart
            !
            IF (RST_NUM .EQ. -2) THEN
                !
                FILE_EXISTS = .TRUE.
                RST_NUM = -2
                !
            END IF
            !
        END DO
        !
        IF (ISTEP - 1 .EQ. 1) THEN
            !
            RST_NUM = RST_NUM + 2
            WRITE(RST_NUM_STR, '(I0)') RST_NUM
            !
        END IF
        !
        FILENAME = 'rst'//TRIM(RST_NUM_STR)//'.control'
        OPEN(UNIT = MYUNIT, FILE = FILENAME, FORM = 'UNFORMATTED', &
            & ACTION = 'WRITE')
        !
        WRITE(MYUNIT) ISTEP
        WRITE(MYUNIT) CURR_LOAD
        WRITE(MYUNIT) PREV_LOAD
        WRITE(MYUNIT) STEP_COMPLETE
        WRITE(MYUNIT) DTIME
        WRITE(MYUNIT) INCR
        WRITE(MYUNIT) TIME
        !
        DO ISURF = 1, NSURFACES
            !
            WRITE(MYUNIT) SURF_LOAD_ARRAY(ISURF,:)
            !
        END DO
        !
        WRITE(MYUNIT) AREA
        WRITE(MYUNIT) AREA0
        WRITE(MYUNIT) LENGTH
        WRITE(MYUNIT) LENGTH0
        WRITE(MYUNIT) CURR_VEL
        WRITE(MYUNIT) S_PERT_MAG
        WRITE(MYUNIT) T_PERT_MAG
        !
        CLOSE(MYUNIT)
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE WRITE_TRIAXCSR_RESTART
    !
    !===========================================================================
    !
    SUBROUTINE WRITE_TRIAXCLR_RESTART(ISTEP, CURR_LOAD, PREV_LOAD, &
        & FIRST_INCR_IN_STEP, INCR, TIME, SURF_LOAD_ARRAY, AREA, AREA0, &
        & LENGTH, LENGTH0, CURR_VEL, PREV_ACTION, CURR_ACTION, &
        & INITIAL_LOAD_DWELL_VEL, INITIAL_UNLOAD_DWELL_VEL)
    !
    !  Write TriaxCLR restart information.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    LOGICAL, INTENT(IN) :: FIRST_INCR_IN_STEP
    INTEGER, INTENT(IN) :: ISTEP
    INTEGER, INTENT(IN) :: INCR
    INTEGER, INTENT(IN) :: PREV_ACTION, CURR_ACTION
    REAL(RK), INTENT(IN) :: CURR_LOAD(3)
    REAL(RK), INTENT(IN) :: PREV_LOAD(3)
    REAL(RK), INTENT(IN) :: TIME
    REAL(RK), INTENT(IN) :: SURF_LOAD_ARRAY(NSURFACES,3)
    REAL(RK), INTENT(IN) :: AREA(NSURFACES)
    REAL(RK), INTENT(IN) :: AREA0(NSURFACES)
    REAL(RK), INTENT(IN) :: LENGTH(3), LENGTH0(3)
    REAL(RK), INTENT(IN) :: CURR_VEL(3)
    REAL(RK), INTENT(IN) :: INITIAL_LOAD_DWELL_VEL(3)
    REAL(RK), INTENT(IN) :: INITIAL_UNLOAD_DWELL_VEL(3)
    !
    ! Locals:
    !
    INTEGER :: MYUNIT
    INTEGER :: RST_NUM
    LOGICAL :: FILE_EXISTS
    CHARACTER(LEN=8) :: RST_NUM_STR
    CHARACTER(LEN=64) :: FILENAME
    INTEGER :: ISURF
    !
    !---------------------------------------------------------------------------
    !
    IF (MYID .EQ. 0) THEN
        !
        MYUNIT = NEWUNITNUMBER()
        !
        RST_NUM = 1000
        FILE_EXISTS = .FALSE.
        !
        DO WHILE (.NOT. FILE_EXISTS)
            !
            WRITE(RST_NUM_STR, '(I0)') RST_NUM
            !
            FILENAME = 'rst'//TRIM(RST_NUM_STR)//'.control'
            !
            INQUIRE(FILE = FILENAME, EXIST = FILE_EXISTS)
            RST_NUM = RST_NUM - 1
            !
            ! If we are on the first simulation, non-restart
            !
            IF (RST_NUM .EQ. -2) THEN
                !
                FILE_EXISTS = .TRUE.
                RST_NUM = -2
                !
            END IF
            !
        END DO
        !
        IF (ISTEP - 1 .EQ. 1) THEN
            !
            RST_NUM = RST_NUM + 2
            WRITE(RST_NUM_STR, '(I0)') RST_NUM
            !
        END IF
        !
        FILENAME = 'rst'//TRIM(RST_NUM_STR)//'.control'
        OPEN(UNIT = MYUNIT, FILE = FILENAME, FORM = 'UNFORMATTED', &
            & ACTION = 'WRITE')
        !
        WRITE(MYUNIT) ISTEP
        WRITE(MYUNIT) CURR_LOAD
        WRITE(MYUNIT) PREV_LOAD
        WRITE(MYUNIT) FIRST_INCR_IN_STEP
        WRITE(MYUNIT) INCR
        WRITE(MYUNIT) TIME
        !
        DO ISURF = 1,NSURFACES
            !
            WRITE(MYUNIT) SURF_LOAD_ARRAY(ISURF, :)
            !
        END DO
        WRITE(MYUNIT) AREA
        WRITE(MYUNIT) AREA0
        WRITE(MYUNIT) LENGTH
        WRITE(MYUNIT) LENGTH0
        WRITE(MYUNIT) CURR_VEL
        WRITE(MYUNIT) PREV_ACTION
        WRITE(MYUNIT) CURR_ACTION
        WRITE(MYUNIT) INITIAL_LOAD_DWELL_VEL
        WRITE(MYUNIT) INITIAL_UNLOAD_DWELL_VEL
        !
        CLOSE(MYUNIT)
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE WRITE_TRIAXCLR_RESTART
    !
END MODULE WRITE_OUTPUT_MOD
