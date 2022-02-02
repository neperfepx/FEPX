! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE READ_INPUT_MOD
!
! Module to handle all generic file read-in functionality. Specific file
! handling (restart, BCs, and powder-diffraction) are maintained in their
! own respective modules for ease of modification. Simulation configuration
! parameters (and their associated scope) is now stored here.
!
! Contains subroutines:
!
! For reading in the mesh file:
! ALLOCATE_MSH: Allocates mesh arrays based on the problem size.
! GET_MSH_SIZE: Scrapes the msh file and retrieve the number of nodes and elems.
! READ_SPATIAL_MSH: Parent subroutine to handle msh file parsing process.
! READ_RESTART_FIELD: Reads field data for restarting a simulation.
!
! For reading in simulation configuration parameters:
! INITIALIZE_OPTIONS - Initializes options and input.
! PROCESS_INPUT - Processes options and input.
!
! Contains helper subroutines:
!
! For parsing the mesh (and related external) file:
! GET_NODE_INFO: Gets the number of nodes in the mesh.
! GET_ELEM_INFO: Gets the number of elements in the mesh and partitions.
! READ_MESH_FORMAT: Reads in gmsh file format - 2.2 0 8 only.
! READ_MESH_VERSION: Reads in gmsh file version.
! READ_NODES: Read in node ID and [X Y Z] spatial coordinates.
! READ_ELEMENTS: Read in ALL elements and only store 3D elements.
! READ_NSETS: Read in node sets - currently not stored for read-in BCs.
! READ_FASETS: Read in surface face node sets (2D elements).
! READ_NODEPARTITIONS: Read in per-node partition distribution. -OPTIONAL
! READ_PHYSICALNAMES: Read in physical names and do not store (gmsh data only).
! READ_ELSETORIENTATIONS: Read in grain orientations.
! READ_ELEMENTORIENTATIONS: Read in element orientations.
! READ_GROUPS: Read in grain phases (assumes single-phase by default). -OPTIONAL
!
! For parsing the simulation configuration file:
! Standard options:
!  EXEC_MAX_INCR
!  EXEC_MAX_TOTAL_TIME
!  EXEC_DEF_CONTROL_BY
!  EXEC_CHECK_NECKING
!  EXEC_LOAD_TOL
!  EXEC_DTIME_FACTOR
!  EXEC_RESTART
!  EXEC_HARD_TYPE
!  EXEC_READ_ORI_FROM_FILE
!  EXEC_READ_PHASE_FROM_FILE
!  EXEC_MAX_ITER_HARD_LIMIT
! Printing options:
!  EXEC_PRINT
!  EXEC_SUPPRESS
! Boundary conditions options:
!  EXEC_READ_BCS_FROM_FILE
!  EXEC_BOUNDARY_CONDITIONS
!  EXEC_LOADING_FACE
!  EXEC_LOADING_DIRECTION
!  EXEC_STRAIN_RATE
!  EXEC_LOAD_RATE
! Crystal parameter options:
!  EXEC_NUMBER_OF_PHASES
!  EXEC_CRYSTAL_TYPE
!  EXEC_PHASEEXEC_PHASE
!  EXEC_M
!  EXEC_GAMMADOT_0
!  EXEC_H_0
!  EXEC_G_0
!  EXEC_G_S0
!  EXEC_M_PRIME
!  EXEC_GAMMADOT_S0
!  EXEC_N
!  EXEC_C11
!  EXEC_C12
!  EXEC_C13
!  EXEC_C44
!  EXEC_C66
!  EXEC_C_OVER_A
!  EXEC_CYCLIC_A
!  EXEC_CYCLIC_C
!  EXEC_LATENT_PARAMETERS
!  EXEC_A_P
!  EXEC_F_P
!  EXEC_B_P
!  EXEC_R_P
! Uniaxial control options:
!  EXEC_NUMBER_OF_STRAIN_STEPS
!  EXEC_NUMBER_OF_STRAIN_RATE_JUMPS
!  EXEC_TARGET_STRAIN
!  EXEC_STRAIN_RATE_JUMP
!  EXEC_NUMBER_OF_LOAD_STEPS
!  EXEC_TARGET_LOAD
! Triaxial CSR options:
!  EXEC_MAX_BC_ITER
!  EXEC_MIN_PERT_FRAC
!  EXEC_LOAD_TOL_ABS
!  EXEC_LOAD_TOL_REL
!  EXEC_NUMBER_OF_CSR_LOAD_STEPS
!  EXEC_TARGET_CSR_LOAD
! Triaxial CLR options:
!  EXEC_MAX_STRAIN_INCR
!  EXEC_MAX_STRAIN
!  EXEC_MAX_EQSTRAIN
!  EXEC_NUMBER_OF_CLR_LOAD_STEPS
!  EXEC_TARGET_CLR_LOAD
!  EXEC_NUMBER_OF_LOAD_RATE_JUMPS
!  EXEC_LOAD_RATE_JUMP
!  EXEC_NUMBER_OF_DWELL_EPISODES
!  EXEC_DWELL_EPISODE
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
USE LIBF95, ONLY: EXECCOMMAND, NEWUNITNUMBER
!
! From libfepx:
!
USE CONVERGENCE_MOD, ONLY: CONVERGENCEKEYWORDINPUT
USE DIMENSIONS_MOD
USE SURFACE_MOD
!
! From libparallel:
!
USE PARALLEL_MOD
!
IMPLICIT NONE
!
! Private
!
INTEGER, PARAMETER, PRIVATE :: MAX_PATH_LEN=1024
INTEGER, PARAMETER, PRIVATE :: MAX_FILE_LEN=256
INTEGER, PRIVATE :: IOERR
INTEGER, PRIVATE :: PRIV_INPUNIT
!
! Public
!
PUBLIC
!
! Mesh format
!
CHARACTER(LEN=5)   :: MESH_VERSION = "2.2"
!
! Element/DOF array bounds
!
INTEGER, PARAMETER :: NDIM  = 10 ! 10-noded tetrahedra
INTEGER, PARAMETER :: NNPE  = NDIM - 1
INTEGER, PARAMETER :: KDIM  = 3*NDIM
INTEGER, PARAMETER :: KDIM1 = KDIM - 1
!
! Array sizes and ranges.
!
INTEGER :: NUMELM, NUMNP, MAXEL, MAXEL1, MAXNP, MAXNP1, MAXDOF, MAXDOF1
INTEGER :: EL_SUB1, EL_SUP1, NP_SUB1, NP_SUP1, DOF_SUB1, DOF_SUP1
INTEGER :: EL_STARTID = 0
!
! Connectivities and coordinates.
!
INTEGER, ALLOCATABLE :: NP(:,:), NODES(:,:)
INTEGER, ALLOCATABLE :: UESUB(:), PHASE(:)
REAL(RK), ALLOCATABLE :: COORDS(:), ELEMENT_CRDS(:,:), &
    & COORDS_ORIG(:)
!
! Orientation options
!
LOGICAL :: ELEMENT_ORIS = .FALSE.
TYPE ORIENTATION_TYPE
    !
    CHARACTER(LEN=50) :: ORIENTATION_PARAMETERIZATION
    CHARACTER(LEN=50) :: ORIENTATION_CONVENTION
    !
END TYPE ORIENTATION_TYPE
TYPE(ORIENTATION_TYPE) :: ORIENTATION_OPTIONS
!
! Microstructure options
!
INTEGER :: NUM_GRAINS, GRAIN_ID_MAX
INTEGER, ALLOCATABLE :: GRAIN_IDS(:), GRAIN_IDS_INV(:)
REAL(RK), ALLOCATABLE :: GRAIN_ORIENTATION(:,:)
!
! Partition flags
!
LOGICAL :: READ_ELEM_PART = .FALSE.
LOGICAL :: READ_NODE_PART = .FALSE.
!
! Elastic strain/stress and post_update_n internal variables
!
REAL(RK), ALLOCATABLE :: PELA_KK_BAR(:)
REAL(RK), ALLOCATABLE :: PSIG_VEC_N(:,:,:)
REAL(RK), ALLOCATABLE :: PCRSS_N(:,:,:)
REAL(RK), ALLOCATABLE :: PRSTAR_N(:,:,:,:)
!
REAL(RK), ALLOCATABLE :: GELA_KK_BAR(:,:,:)
REAL(RK), ALLOCATABLE :: GSIG_VEC_N(:,:,:,:)
!
! Internal variables needed for work related calculations
!
REAL(RK), ALLOCATABLE :: EL_WORK_N(:)
REAL(RK), ALLOCATABLE :: EL_WORKP_N(:)
REAL(RK), ALLOCATABLE :: EL_WORK_RATE_N(:)
REAL(RK), ALLOCATABLE :: EL_WORKP_RATE_N(:)
REAL(RK), ALLOCATABLE :: EL_WORK(:)
REAL(RK), ALLOCATABLE :: EL_WORKP(:)
REAL(RK), ALLOCATABLE :: D_TOT(:, :, :)
!
! Configuration types
!
TYPE OPTIONS_TYPE
    !
    INTEGER :: MAX_INCR
    REAL(RK) :: MAX_TOTAL_TIME
    INTEGER :: DEF_CONTROL_BY
    LOGICAL :: CHECK_NECKING
    REAL(RK) :: LOAD_TOL
    REAL(RK) :: DTIME_FACTOR
    LOGICAL :: RESTART
    INTEGER :: RESTART_FILE_HANDLING
    INTEGER :: RESTART_INITIAL_STEP
    LOGICAL :: SAT_EVO
    LOGICAL :: PRECIP_HARD
    CHARACTER(LEN=18) :: HARD_TYPE
    LOGICAL :: READ_ORI_FROM_FILE
    CHARACTER(LEN=MAX_FILE_LEN) :: ORI_FILE
    LOGICAL :: READ_PHASE_FROM_FILE
    CHARACTER(LEN=MAX_FILE_LEN) :: PHASE_FILE
    INTEGER :: MAX_ITER_HARD_LIMIT
    !
END TYPE OPTIONS_TYPE
!
TYPE PRINT_OPTIONS_TYPE
    !
    LOGICAL :: PRINT_COORDINATES
    LOGICAL :: PRINT_DISPLACEMENTS
    LOGICAL :: PRINT_VELOCITIES
    LOGICAL :: PRINT_ORIENTATIONS
    LOGICAL :: PRINT_CRSS
    LOGICAL :: PRINT_EQSTRESS
    LOGICAL :: PRINT_STRESS
    LOGICAL :: PRINT_STRAIN_TOT
    LOGICAL :: PRINT_STRAIN_EL
    LOGICAL :: PRINT_STRAIN_PL
    LOGICAL :: PRINT_GAMMADOT
    LOGICAL :: PRINT_DEFF
    LOGICAL :: PRINT_DPEFF
    LOGICAL :: PRINT_ELVOL
    LOGICAL :: PRINT_EQSTRAIN
    LOGICAL :: PRINT_EQELSTRAIN
    LOGICAL :: PRINT_EQPLSTRAIN
    LOGICAL :: PRINT_VGRAD
    LOGICAL :: PRINT_DPHAT
    LOGICAL :: PRINT_WPHAT
    LOGICAL :: PRINT_GAMMA
    LOGICAL :: PRINT_RESTART
    LOGICAL :: PRINT_FORCES
    LOGICAL :: PRINT_CONV
    LOGICAL :: PRINT_WORK
    LOGICAL :: PRINT_WORKP
    LOGICAL :: PRINT_DEFRATE
    LOGICAL :: PRINT_WORKRATE
    LOGICAL :: PRINT_WORKRATEP
    !
END TYPE PRINT_OPTIONS_TYPE
!
TYPE BOUNDARY_CONDITIONS_TYPE
    !
    LOGICAL  :: READ_BCS_FROM_FILE
    CHARACTER(LEN=MAX_FILE_LEN) :: BCS_FILE
    INTEGER  :: BOUNDARY_CONDITIONS
    INTEGER  :: LOADING_FACE
    INTEGER  :: LOADING_DIRECTION
    REAL(RK) :: STRAIN_RATE
    REAL(RK) :: LOAD_RATE
    !
END TYPE BOUNDARY_CONDITIONS_TYPE
!
TYPE CRYS_OPTIONS_TYPE
    !
    INTEGER :: NUMBER_OF_PHASES
    INTEGER :: PHASE
    ! Required crystal parameters for all phases.
    INTEGER,  ALLOCATABLE :: CRYSTAL_TYPE(:)
    REAL(RK), ALLOCATABLE :: M(:)
    REAL(RK), ALLOCATABLE :: GAMMADOT_0(:)
    REAL(RK), ALLOCATABLE :: H_0(:)
    REAL(RK), ALLOCATABLE :: G_0(:)
    REAL(RK), ALLOCATABLE :: G_S0(:)
    REAL(RK), ALLOCATABLE :: M_PRIME(:)
    REAL(RK), ALLOCATABLE :: GAMMADOT_S0(:)
    REAL(RK), ALLOCATABLE :: N(:)
    REAL(RK), ALLOCATABLE :: C11(:)
    REAL(RK), ALLOCATABLE :: C12(:)
    REAL(RK), ALLOCATABLE :: C13(:) ! HCP and BCT only
    REAL(RK), ALLOCATABLE :: C44(:)
    REAL(RK), ALLOCATABLE :: C66(:) ! BCT only
    REAL(RK), ALLOCATABLE :: C_OVER_A(:) ! HCP and BCT only
    REAL(RK), ALLOCATABLE :: PRISMATIC_TO_BASAL(:) ! HCP only
    REAL(RK), ALLOCATABLE :: PYRAMIDAL_TO_BASAL(:) ! HCP only
    REAL(RK), ALLOCATABLE :: HRATIO_BCT_A(:) ! BCT only
    REAL(RK), ALLOCATABLE :: HRATIO_BCT_B(:) ! BCT only
    REAL(RK), ALLOCATABLE :: HRATIO_BCT_C(:) ! BCT only
    REAL(RK), ALLOCATABLE :: HRATIO_BCT_D(:) ! BCT only
    REAL(RK), ALLOCATABLE :: HRATIO_BCT_E(:) ! BCT only
    REAL(RK), ALLOCATABLE :: HRATIO_BCT_F(:) ! BCT only
    REAL(RK), ALLOCATABLE :: HRATIO_BCT_G(:) ! BCT only
    REAL(RK), ALLOCATABLE :: HRATIO_BCT_H(:) ! BCT only
    REAL(RK), ALLOCATABLE :: HRATIO_BCT_I(:) ! BCT only
    ! Hardening-related parameters.
    REAL(RK), ALLOCATABLE :: CYCLIC_A(:)
    REAL(RK), ALLOCATABLE :: CYCLIC_C(:)
    REAL(RK), ALLOCATABLE :: LATENT_PARAMETERS(:,:)
    ! Precipitation hardening parameters
    REAL(RK), ALLOCATABLE :: A_P(:)
    REAL(RK), ALLOCATABLE :: F_P(:)
    REAL(RK), ALLOCATABLE :: B_P(:)
    REAL(RK), ALLOCATABLE :: R_P(:)
    !
    ! Rate dependence options for HCP materials
    LOGICAL,  ALLOCATABLE :: USE_ANISO_M(:)
    REAL(RK), ALLOCATABLE :: ANISO_M(:, :)
    !
END TYPE CRYS_OPTIONS_TYPE
!
TYPE UNIAXIAL_CONTROL_TYPE
    !
    ! Uniaxial strain target parameters.
    INTEGER :: NUMBER_OF_STRAIN_STEPS
    INTEGER :: NUMBER_OF_STRAIN_RATE_JUMPS
    INTEGER :: CURRENT_STRAIN_TARGET
    INTEGER :: CURRENT_STRAIN_JUMP
    REAL(RK), ALLOCATABLE :: TARGET_STRAIN(:,:)
    REAL(RK), ALLOCATABLE :: STRAIN_RATE_JUMP(:,:)
    CHARACTER(LEN=15) :: TEMP
    ! Uniaxial load target parameters.
    INTEGER :: NUMBER_OF_LOAD_STEPS
    INTEGER :: CURRENT_LOAD_TARGET
    REAL(RK), ALLOCATABLE :: TARGET_LOAD(:,:)
    !
END TYPE UNIAXIAL_CONTROL_TYPE
!
TYPE TRIAXCSR_OPTIONS_TYPE
    !
    INTEGER  :: MAX_BC_ITER
    REAL(RK) :: MAX_STRAIN
    REAL(RK) :: MAX_EQSTRAIN
    REAL(RK) :: MIN_PERT_FRAC
    REAL(RK) :: LOAD_TOL_ABS
    REAL(RK) :: LOAD_TOL_REL
    CHARACTER(LEN=15) :: TEMP_CSR
    INTEGER :: NUMBER_OF_CSR_LOAD_STEPS
    INTEGER :: CURRENT_CSR_LOAD_TARGET
    REAL(RK), ALLOCATABLE :: TARGET_CSR_LOAD(:,:)
    !
END TYPE TRIAXCSR_OPTIONS_TYPE
!
TYPE TRIAXCLR_OPTIONS_TYPE
    !
    REAL(RK) :: LOAD_TOL_ABS
    INTEGER  :: MAX_BC_ITER
    REAL(RK) :: MAX_STRAIN_INCR
    REAL(RK) :: MAX_STRAIN
    REAL(RK) :: MAX_EQSTRAIN
    CHARACTER(LEN=15) :: TEMP_CLR
    INTEGER :: NUMBER_OF_CLR_LOAD_STEPS
    INTEGER :: CURRENT_CLR_LOAD_TARGET
    REAL(RK), ALLOCATABLE :: TARGET_CLR_LOAD(:,:)
    INTEGER :: NUMBER_OF_LOAD_RATE_JUMPS
    INTEGER :: NUMBER_OF_DWELL_EPISODES
    INTEGER :: CURRENT_LOAD_RATE_JUMP
    INTEGER :: CURRENT_DWELL_EPISODE
    REAL(RK), ALLOCATABLE :: LOAD_RATE_JUMP(:,:)
    REAL(RK), ALLOCATABLE :: DWELL_EPISODE(:,:)
    !
END TYPE TRIAXCLR_OPTIONS_TYPE
!
TYPE(OPTIONS_TYPE) :: OPTIONS
TYPE(PRINT_OPTIONS_TYPE) :: PRINT_OPTIONS
TYPE(BOUNDARY_CONDITIONS_TYPE) :: BCS_OPTIONS
TYPE(CRYS_OPTIONS_TYPE) :: CRYS_OPTIONS
TYPE(UNIAXIAL_CONTROL_TYPE) :: UNIAXIAL_OPTIONS
TYPE(TRIAXCSR_OPTIONS_TYPE) :: TRIAXCSR_OPTIONS
TYPE(TRIAXCLR_OPTIONS_TYPE) :: TRIAXCLR_OPTIONS
!
!  Deformation control :
!
INTEGER, PARAMETER :: UNIAXIAL_LOAD_TARGET = 1
INTEGER, PARAMETER :: UNIAXIAL_STRAIN_TARGET = 2
INTEGER, PARAMETER :: TRIAXIAL_CONSTANT_STRAIN_RATE = 3
INTEGER, PARAMETER :: TRIAXIAL_CONSTANT_LOAD_RATE = 4
!
CONTAINS
    !
    SUBROUTINE ALLOCATE_MSH(STATUS)
    !
    ! Allocate mesh arrays according to problem size.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! STATUS: Returns status of allocation to main program.
    !
    INTEGER, INTENT(INOUT) :: STATUS
    !
    !---------------------------------------------------------------------------
    !
    STATUS = 0
    !
    ALLOCATE (NP(0:NNPE, EL_SUB1:EL_SUP1), &
        & NODES(0:KDIM1, EL_SUB1:EL_SUP1), &
        & COORDS(DOF_SUB1:DOF_SUP1), &
        & COORDS_ORIG(DOF_SUB1:DOF_SUP1), &
        & ELEMENT_CRDS(0:KDIM1, EL_SUB1:EL_SUP1), &
        & UESUB(0:NUMELM-1), &
        & PHASE(0:NUMELM-1), &
        & STAT=STATUS)
    !
    ! Initialize these arrays here before the spatial mesh is parsed.
    !
    ! Note: UESUB stores the 1-indexed elset for a element of ID:(0:NUMELM-1)
    ! and we set phase to 1 to assume a single-phase simulation by default.
    !
    UESUB = 0
    PHASE = 1
    !
    RETURN
    !
    END SUBROUTINE ALLOCATE_MSH
    !
    !===========================================================================
    !
    SUBROUTINE GET_MSH_SIZE(IO)
    !
    ! Scrapes $Nodes and $Elements fields for the number of nodes and 3D elems
    ! to determine problem size. Also, attempt to retrieve embedded mesh
    ! partitioning info from $Elements and $NodePartitions (if available).
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! EOFSTAT: Value that confirms if FORTRAN is at the end of a file.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: EOFSTAT
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Initialize EOFSTAT and loop over the entire file until EOF is reached.
    EOFSTAT = 0
    !
    DO WHILE (EOFSTAT .EQ. 0)
        !
        ! Read in line and determine which section is to be parsed.
        ! The below line is a temporary fix to avoid program hangs - JC
        LINE = ""
        !
        READ(IO, '(A)', IOSTAT = EOFSTAT) LINE
        !
        SELECT CASE (LINE)
            !
            CASE('$Nodes')
                !
                BACKSPACE(IO)
                CALL GET_NODE_INFO(IO)
                !
            CASE('$Elements')
                !
                BACKSPACE(IO)
                CALL GET_ELEM_INFO(IO)
                !
            CASE('$NodePartitions')
                !
                BACKSPACE(IO)
                CALL READ_NODEPARTITIONS(IO)
                !
            ! END CASE
            !
        END SELECT
        !
    END DO
    !
    REWIND (IO)
    !
    RETURN
    !
    END SUBROUTINE GET_MSH_SIZE
    !
    !===========================================================================
    !
    SUBROUTINE GET_NODE_INFO(IO)
    !
    ! Read the $Nodes field and retrieves the total number of nodes in mesh.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I: Generic looping index.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, I
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$Nodes')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &nodes.')
        !
    END IF
    !
    ! Read in the number of nodes in the section (and set global NUMNP).
    READ(IO, *) NUMNP
    !
    DO I = 1, NUMNP
        !
        READ(IO, *)
        !
    END DO
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndNodes')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &nodes.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE GET_NODE_INFO
    !
    !===========================================================================
    !
    SUBROUTINE GET_ELEM_INFO(IO)
    !
    ! Read the $Elements field and retrieves the total number of nodes in mesh.
    ! Also, retrieves elemental partition information if present in the mesh.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I: Generic looping index.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! NSPACE/NVALS: Storage values used for dividing a line into substrings.
    ! ELMTYPE: Gmsh defined value that determines the element type.
    ! NELM: Locally determined number of 3D elements in mesh.
    ! MINBND/MAXBND: Array bounds to assist in trimming element partition array.
    ! NPARTS: Number of partitions read in from the mesh.
    ! IPART: Partition loop index to check which partition is being checked.
    ! LINE: Input line on current record to be parsed.
    ! IARRAY: Substring array for internal read parsing.
    ! TEMP: Temporary array storage for element paritions including non-3D.
    ! ELEM_PARTS: Array storing the trimmed element paritions.
    ! TEMPVAL: Array of EL_SUB1/EL_SUP1 bounds for each partition read in.
    ! CHG_INDEX: Index in array immediately after a change in value occurs.
    ! MASK: Masking array to trim TEMP into ELEM_PARTS.
    !
    INTEGER :: IERR, I
    INTEGER :: NLINES, NSPACE, NVALS, ELMTYPE, NELM
    INTEGER :: MINBND, MAXBND, NPARTS, IPART
    CHARACTER(LEN=256) :: LINE
    CHARACTER(LEN=12)  :: IARRAY(16)
    INTEGER, ALLOCATABLE :: TEMP(:), ELEM_PARTS(:), TEMPVAL(:,:), CHG_INDEX(:)
    LOGICAL, ALLOCATABLE :: MASK(:)
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$Elements')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &elements.')
        !
    END IF
    !
    ! Read in the number of elements in the section (and set NLINES).
    READ(IO, *) NLINES
    !
    ! Initialize the 3D elem counter and partition array.
    NELM = 0
    !
    ALLOCATE(TEMP(1:NLINES))
    ALLOCATE(MASK(1:NLINES))
    TEMP = -1
    !
    DO I = 1, NLINES
        !
        READ(IO, '(A)') LINE
        !
        ! Trim the full string into NVALS number of substrings to parse.
        NSPACE = COUNT(TRANSFER(LINE, 'A', LEN_TRIM(LINE)) == " ")
        NVALS = NSPACE + 1
        IARRAY = ""
        !
        ! Internal read of line to store substrings.
        READ(LINE, *) IARRAY(1:NVALS)
        READ(IARRAY(2), *) ELMTYPE
        !
        IF (ELMTYPE .EQ. 11) THEN
            !
            if (EL_STARTID .EQ. 0) THEN
            !
                READ(IARRAY(1), *) EL_STARTID
            !
            END IF

            ! Assign data from IARRAY(6) the TAG3 value which contain partition.
            READ(IARRAY(6), *) TEMP(I)
            !
            ! Increment NELM by 1 to ensure NUMELM is assigned correctly.
            NELM = NELM + 1
            !
        END IF
        !
    END DO
    !
    ! Check to see if any ELMTYPE=11 have been parsed (if not, quit)
    IF (NELM .EQ. 0) THEN
        !
        CALL PAR_QUIT('Error  :     > Mesh requires second-order elements.')
        !
    END IF
    !
    ! Set the global NUMELM value.
    NUMELM = NELM
    !
    ! Reallocate TEMP array to properly sized array (1:NUMELM).
    MASK = TEMP .EQ. -1
    MINBND = COUNT(MASK) + 1
    MAXBND = SIZE(TEMP)
    !
    ELEM_PARTS = TEMP(MINBND:MAXBND)
    DEALLOCATE(MASK)
    DEALLOCATE(TEMP)
    !
    ! Set a flag if we read in paritions and strip the ELEM_PARTS array down.
    IF (MAXVAL(ELEM_PARTS) .GT. 0) THEN
        !
        READ_ELEM_PART = .TRUE.
        NPARTS = MAXVAL(ELEM_PARTS)
        ALLOCATE(TEMPVAL(1:NPARTS, 2))
        ALLOCATE(CHG_INDEX(NPARTS-1))
        IPART = 1
        !
        ! We can split ELEM_PARTS in 0-indexed EL_SUB1 and EL_SUP1 here.
        DO I = 1, SIZE(ELEM_PARTS)
            !
            ! A change has been detected so do something.
            IF (ELEM_PARTS(I) .NE. IPART) THEN
                !
                CHG_INDEX(IPART) = I
                IF (IPART .NE. NPARTS) IPART = IPART + 1
                !
            END IF
            !
        END DO
        !
        ! Build bound array for NPARTS -> (1,:) is EL_SUB1, (2,:) is EL_SUP1.
        TEMPVAL = 0
        DO I = 1, NPARTS
            !
            ! Need to handle the edge cases differently.
            IF (I .EQ. 1) THEN
                !
                TEMPVAL(1,1) = 1
                TEMPVAL(1,2) = CHG_INDEX(I) - 1
                !
            ELSE IF (I .EQ. NPARTS) THEN
                !
                TEMPVAL(I,1) = CHG_INDEX(I-1)
                TEMPVAL(I,2) = SIZE(ELEM_PARTS)
                !
            ELSE
                !
                TEMPVAL(I,1) = CHG_INDEX(I-1)
                TEMPVAL(I,2) = CHG_INDEX(I) - 1
                !
            END IF
            !
        END DO
        !
    END IF
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndElements')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &elements.')
        !
    END IF
    !
    ! Do not need to maintain element partitions for now so deallocate.
    IF (ALLOCATED(ELEM_PARTS)) DEALLOCATE(ELEM_PARTS)
    IF (ALLOCATED(CHG_INDEX))  DEALLOCATE(CHG_INDEX)
    IF (ALLOCATED(TEMPVAL))    DEALLOCATE(TEMPVAL)
    !
    RETURN
    !
    END SUBROUTINE GET_ELEM_INFO
    !
    !===========================================================================
    !
    SUBROUTINE READ_SPATIAL_MSH(IO, CONSOLE_UNIT)
    !
    ! Parse the .msh file and have each processor store only what it needs.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    INTEGER, INTENT(IN) :: CONSOLE_UNIT
    !
    ! Locals:
    ! EOFSTAT: Value that confirms if FORTRAN is at the end of a file.
    ! I: Generic looping index.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: EOFSTAT, I
    CHARACTER(LEN=256) :: LINE
    !
    ! Notes:
    ! The .msh file follows the Gmsh format 2.2 and is parsed in chunks
    ! and the field headers ($FieldName) are read to determine which
    ! helper subroutines to call.
    !
    !---------------------------------------------------------------------------
    !
    ! Initialize EOFSTAT and loop over the entire file until EOF is reached.
    EOFSTAT = 0
    !
    DO WHILE (EOFSTAT .EQ. 0)
        !
        ! Read in line and determine which section is to be parsed.
        ! The below line is a temporary fix to avoid program hangs - JC
        LINE = ""
        !
        READ(IO, '(A)', IOSTAT = EOFSTAT) LINE
        !
        SELECT CASE (LINE)
            !
            CASE('$MeshFormat')
                !
                BACKSPACE(IO)
                CALL READ_MESH_FORMAT(IO)
                !
            CASE('$MeshVersion')
                !
                BACKSPACE(IO)
                CALL READ_MESH_VERSION(IO)
                !
            CASE('$Nodes')
                !
                BACKSPACE(IO)
                CALL READ_NODES(IO)
                !
            CASE('$Elements')
                !
                BACKSPACE(IO)
                CALL READ_ELEMENTS(IO)
                !
            CASE('$NSets') ! Nodal values for BCs assignments (OPTIONAL).
                !
                BACKSPACE(IO)
                CALL READ_NSETS(IO)
                !
            CASE('$Fasets')
                !
                BACKSPACE(IO)
                CALL READ_FASETS(IO)
                !
                DO I = 1, NSURFACES
                    !
                    CALL PART_GATHER(SURFACES(I)%CRDS, COORDS, &
                        &SURFACES(I)%CONN3D, SURFACES(I)%TR)
                    !
                ENDDO
                !
            CASE('$PhysicalNames')
                !
                BACKSPACE(IO)
                CALL READ_PHYSICALNAMES(IO)
                !
            CASE('$ElsetOrientations') ! Per-grain orientations.
                !
                BACKSPACE(IO)
                CALL READ_ELSETORIENTATIONS(IO)
                !
            CASE('$ElementOrientations') ! Per-element orientations (OPTIONAL).
                !
                BACKSPACE(IO)
                CALL READ_ELEMENTORIENTATIONS(IO)
                !
            CASE('$Groups') ! Grain/phase assignment for multiphase (OPTIONAL).
                !
                BACKSPACE(IO)
                CALL READ_GROUPS(IO, CONSOLE_UNIT)
                !
            ! END CASE
            !
        END SELECT
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE READ_SPATIAL_MSH
    !
    !===========================================================================
    !
    SUBROUTINE READ_RESTART_FIELD(VELOCITY, C0_ANGS, C_ANGS, &
        & RSTAR, RSTAR_N, WTS, CRSS, CRSS_N, &
        & E_ELAS_KK_BAR, SIG_VEC_N, EQSTRAIN, EQPLSTRAIN, GAMMA, &
        & EL_WORK_N, EL_WORKP_N, EL_WORK_RATE_N, EL_WORKP_RATE_N, PLSTRAIN, &
        & TOTSTRAIN)
    !
    !  Reads field data for restarting simulation
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    REAL(RK), INTENT(OUT) :: VELOCITY(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(OUT) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK_BAR(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: SIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: EQSTRAIN(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: EQPLSTRAIN(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: GAMMA(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: EL_WORK_N(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: EL_WORKP_N(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: EL_WORK_RATE_N(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: EL_WORKP_RATE_N(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: PLSTRAIN(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: TOTSTRAIN(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
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
    !----------------------------------------------------------------------
    !
    MYUNIT = NEWUNITNUMBER()
    WRITE(CHARID, '(I0)') MYID + 1
    !
    RST_NUM = 1000
    FILE_EXISTS = .FALSE.
    !
    ! Find max value N of rstN.control
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
        IF (RST_NUM .EQ. -2) THEN
            !
            CALL PAR_QUIT('Error  :     > Restart control file not found.')
            !
        END IF
        !
    END DO
    !
    RST_NUM = RST_NUM + 1
    WRITE(RST_NUM_STR, '(I0)') RST_NUM
    !
    FILENAME = 'rst'//TRIM(RST_NUM_STR)//'.field.core'//TRIM(CHARID)
    !
    OPEN(UNIT = MYUNIT, FILE = FILENAME, FORM = 'UNFORMATTED', ACTION = 'READ')
    !
    ! Velocity and coordinates.
    !
    READ(MYUNIT) COORDS
    READ(MYUNIT) VELOCITY
    !
    ! Orientations, weights and hardnesses.
    !
    READ(MYUNIT) C0_ANGS
    READ(MYUNIT) C_ANGS
    READ(MYUNIT) RSTAR
    READ(MYUNIT) RSTAR_N
    READ(MYUNIT) WTS
    READ(MYUNIT) CRSS
    READ(MYUNIT) CRSS_N
    !
    ! Elastic Strains.
    !
    READ(MYUNIT) GELA_KK_BAR
    READ(MYUNIT) GSIG_VEC_N
    READ(MYUNIT) PELA_KK_BAR
    READ(MYUNIT) PSIG_VEC_N
    READ(MYUNIT) E_ELAS_KK_BAR
    READ(MYUNIT) SIG_VEC_N
    !
    ! Equivalent Strains.
    !
    READ(MYUNIT) EQSTRAIN
    READ(MYUNIT) EQPLSTRAIN
    READ(MYUNIT) GAMMA
    !
    ! Total work, plastic work, and rates.
    !
    READ(MYUNIT) EL_WORK_N
    READ(MYUNIT) EL_WORKP_N
    READ(MYUNIT) EL_WORK_RATE_N
    READ(MYUNIT) EL_WORKP_RATE_N
    !
    ! Other integrated quantities.
    !
    READ(MYUNIT) PLSTRAIN
    READ(MYUNIT) TOTSTRAIN
    !
    CLOSE(MYUNIT)
    !
    END SUBROUTINE READ_RESTART_FIELD
    !
    !===========================================================================
    !
    SUBROUTINE READ_MESH_FORMAT(IO)
    !
    ! Read the mesh file format and confirm it is of correct type.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! NSPACE/NVALS: Storage values used for dividing a line into substrings.
    ! FILE_TYPE: Gmsh defined file type of mesh - must be '0'.
    ! DATA_SIZE: Gmsh defined data size of mesh - must be '8'.
    ! MESH_FORMAT: Gmsh defined mesh file format - must be '2.2'.
    ! IARRAY: Substring array for internal read parsing.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, NSPACE, NVALS
    INTEGER :: FILE_TYPE, DATA_SIZE
    CHARACTER(LEN=3)   :: MESH_FORMAT
    CHARACTER(LEN=12)  :: IARRAY(16)
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$MeshFormat')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &mesh format.')
        !
    END IF
    !
    ! Read in the mesh format string.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    ! Trim the full string into NVALS number of substrings to parse.
    NSPACE = COUNT(TRANSFER(LINE, 'A', LEN_TRIM(LINE)) == " ")
    NVALS = NSPACE + 1
    IARRAY = ""
    !
    ! Internal read of line to store substrings.
    READ(LINE, *) IARRAY(1:NVALS)
    READ(IARRAY(1), *) MESH_FORMAT
    READ(IARRAY(2), *) FILE_TYPE
    READ(IARRAY(3), *) DATA_SIZE
    !
    IF (MESH_FORMAT .NE. '2.2') &
        &CALL PAR_QUIT('Error  :     > Incorrect mesh format version provided.')
    IF (FILE_TYPE .NE. 0) &
        &CALL PAR_QUIT('Error  :     > Incorrect mesh file type provided.')
    IF (DATA_SIZE .NE. 8) &
        &CALL PAR_QUIT('Error  :     > Incorrect mesh data size provided.')
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndMeshFormat')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &mesh format.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_MESH_FORMAT
    !
    !===========================================================================
    !
    SUBROUTINE READ_MESH_VERSION(IO)
    !
    ! Read the mesh file version
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$MeshVersion')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &mesh version.')
        !
    END IF
    !
    ! Read in the mesh version string.
    READ(IO, '(A)', IOSTAT = IERR) MESH_VERSION
    !
    IF (MESH_VERSION(1:3) .NE. '2.2') &
        &CALL PAR_QUIT('Error  :     > Incorrect mesh version provided.')
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndMeshVersion')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &mesh version.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_MESH_VERSION
    !
    !===========================================================================
    !
    SUBROUTINE READ_NODES(IO)
    !
    ! Read in node ID and coordinates and store per processor.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! INODE: Node 1-indexed ID value.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! I: Generic loop index.
    ! K1/K2/K3: DOF mapping values for the node coordinates.
    ! X/Y/Z: Cartesian coordinate value of the node.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER  :: IERR, INODE, NLINES
    INTEGER  :: I, K1, K2, K3
    REAL(RK) :: X, Y, Z
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$Nodes')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &nodes.')
        !
    END IF
    !
    ! Read in the number of nodes in the section (and set NLINES).
    READ(IO, *) NLINES
    !
    DO I = 1, NLINES
        !
        READ(IO, *) INODE, X, Y, Z
        !
        ! Need to shift by 1 to store zero-indexed node ID for internal use.
        INODE = INODE - 1
        !
        ! Store the nodal coords only if within local processor range.
        IF ((INODE .GE. NP_SUB1) .AND. (INODE .LE. NP_SUP1)) THEN
            !
            K1 = 3 * INODE
            K2 = K1 + 1
            K3 = K2 + 1
            !
            COORDS(K1) = X
            COORDS(K2) = Y
            COORDS(K3) = Z
            !
        END IF
        !
    END DO
    !
    ! Make copy of original coordinates.
    !
    COORDS_ORIG = COORDS
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndNodes')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &nodes.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_NODES
    !
    !===========================================================================
    !
    SUBROUTINE READ_ELEMENTS(IO)
    !
    ! Read in all elements [0D -> 3D] and only store 3D elements within
    ! local processor range.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I/J: Generic loop index.
    ! J1/J2/K3/K1/K2/K3: DOF mapping values for the element nodes.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! NSPACE/NVALS: Storage values used for dividing a line into substrings.
    ! ELMTYPE: Gmsh defined value that determines the element type.
    ! NELM: Locally determined number of 3D elements in mesh.
    ! NTAGS: Number of tags in the 3D element line - must be 3.
    ! GRAIN_ID: Grain (or elset) ID - 1-indexed.
    ! NODES_FE: Local element nodal connectivity array.
    ! IARRAY: Substring array for internal read parsing.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, I, J, J1, J2, J3, K1, K2, K3
    INTEGER :: NLINES, NSPACE, NVALS, ELMTYPE, NELM, NTAGS, GRAIN_ID, GRAIN_POS
    INTEGER :: NODES_FE(0:NNPE)
    CHARACTER(LEN=12) :: IARRAY(16)
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$Elements')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &elements.')
        !
    END IF
    !
    ! Read in the number of elements in the section (and set NLINES).
    READ(IO, *) NLINES
    ALLOCATE(GRAIN_IDS(0:NLINES-1))
    !
    NELM = 0
    NUM_GRAINS = 0
    !
    DO I = 1, NLINES
        !
        READ(IO, '(A)') LINE
        !
        ! Trim the full string into NVALS number of substrings to parse.
        NSPACE = COUNT(TRANSFER(LINE, 'A', LEN_TRIM(LINE)) == " ")
        NVALS = NSPACE + 1
        IARRAY = ""
        !
        ! Internal read of line to store substrings.
        READ(LINE, *) IARRAY(1:NVALS)
        READ(IARRAY(2), *) ELMTYPE
        !
        IF (ELMTYPE .EQ. 11) THEN
            !
            ! If 10-node tetrahedral element, extract information from the line.
            ! The line has format of:
            ! IELEM ELMTYPE NTAGS GRAIN_ID TAG2 TAG3 NODE0 ... NODE9
            READ(IARRAY(3), *) NTAGS
            READ(IARRAY(4), *) GRAIN_ID ! elset
            !
            ! Do we know this grain?
            GRAIN_POS = -1
            DO J = 0, NUM_GRAINS - 1
                !
                IF (GRAIN_IDS(J) .EQ. GRAIN_ID) THEN
                    !
                    GRAIN_POS = J
                    ! escape here
                END IF
                !
            END DO
            !
            IF (GRAIN_POS .EQ. -1) THEN
                !
                NUM_GRAINS = NUM_GRAINS + 1
                GRAIN_POS = NUM_GRAINS - 1
                GRAIN_IDS (GRAIN_POS) = GRAIN_ID
                !
            END IF
            !
            ! Set the UESUB value for this element (assigns elset/grain).
            UESUB(NELM) = GRAIN_POS
            !
            IF (NTAGS .EQ. 3) THEN
                !
                ! The .msh format modifies the local node order of an element.
                ! The remapping follows from Gmsh to FEPX:
                ! GMSH ORDER: 1 2 3 4 5 6 7 8 9 10
                ! FEPX ORDER: 1 3 5 10 2 4 6 7 9 8
                READ(IARRAY(7),*) NODES_FE(0)
                READ(IARRAY(8),*) NODES_FE(2)
                READ(IARRAY(9),*) NODES_FE(4)
                READ(IARRAY(10),*) NODES_FE(9)
                READ(IARRAY(11),*) NODES_FE(1)
                READ(IARRAY(12),*) NODES_FE(3)
                READ(IARRAY(13),*) NODES_FE(5)
                READ(IARRAY(14),*) NODES_FE(6)
                READ(IARRAY(15),*) NODES_FE(8)
                READ(IARRAY(16),*) NODES_FE(7)
                !
                ! Shift the array in order to maintain internal zero-indexing.
                NODES_FE(0:9) = NODES_FE(0:9) - 1
                !
                ! Store elements only on local processor range.
                IF ((NELM .GE. EL_SUB1) .AND. (NELM .LE. EL_SUP1)) THEN
                    !
                    NP(:, NELM) = NODES_FE
                    !
                    ! Map the local DOF to global DOF.
                    DO J = 0, NNPE
                        !
                        J1 = 3 * J
                        J2 = J1 + 1
                        J3 = J2 + 1
                        !
                        K1 = 3 * NODES_FE(J)
                        K2 = K1 + 1
                        K3 = K2 + 1
                        !
                        NODES(J1, NELM) = K1
                        NODES(J2, NELM) = K2
                        NODES(J3, NELM) = K3
                        !
                    END DO
                    !
                END IF
                !
            END IF
            !
            ! Increment NELM by 1 to ensure NP is assigned correctly
            NELM = NELM + 1
            !
        END IF
        !
    END DO
    !
    GRAIN_ID_MAX = 0
    !
    DO I = 0, NUM_GRAINS - 1
        !
        GRAIN_ID_MAX = MAX (GRAIN_IDS (I), GRAIN_ID_MAX)
        !
    END DO
    !
    ALLOCATE(GRAIN_IDS_INV(0:GRAIN_ID_MAX))
    GRAIN_IDS_INV = 0
    !
    DO I = 0, NUM_GRAINS - 1
        !
        GRAIN_IDS_INV(GRAIN_IDS(I)) = I
        !
    END DO
    !
    ! Confirm that NELM (local) and NUMELM (global) values match
    IF (NELM .NE. NUMELM) CALL PAR_QUIT('Error  :     > &
        &Number of elements in $Elements does not match problem size.')
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndElements')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &elements.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_ELEMENTS
    !
    !===========================================================================
    !
    SUBROUTINE READ_NSETS(IO)
    !
    ! Read in node sets and do not store - currently unused, but will be
    ! integrated into custom boundary conditions in the future.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I/J: Generic loop index.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! NSETS: Number of NSets within the field.
    ! NSETLABEL: NSet label.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, I, J, NLINES, NSETS
    CHARACTER(LEN=12) :: NSETLABEL
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$NSets')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &nsets.')
        !
    END IF
    !
    ! Read in the number of nsets in the section (and set NSETS).
    READ(IO, *) NSETS
    !
    ! Loop over the number of NSETS and extract information per set.
    DO I = 1, NSETS
        !
        ! Read in the NSet label.
        READ(IO, *) NSETLABEL
        !
        ! Read in the number of nodes in this set.
        READ(IO, *) NLINES
        !
        ! Loop over single set and do NOT store.
        DO J = 1, NLINES
            !
            READ(IO, *)
            !
        END DO
        !
    END DO
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndNSets')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &nsets.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_NSETS
    !
    !===========================================================================
    !
    SUBROUTINE READ_FASETS(IO)
    !
    ! Read in facesets and store the surface 2D elements. This subroutine also
    ! prepares the 3D connectivity for the surface nodes which is later used
    ! to integrate the loads on each surface.
    !
    ! Notes:
    ! This information should be readily extracted from the READ_ELEMENTS
    ! subroutine by searching for ELM_TYPE (9) and extracting, however, this
    ! would not provide faset labels or IDs.
    !
    ! Legacy note: type = 6 for triangle - tsh
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I/J/K: Generic loop index.
    ! TYPE/TYPE1/TYPE3: Parameters to define the internal element type/bounds.
    ! ELEM_ID: Surface element ID - 1-indexed.
    ! ELEM_NODES/ELEM_NODES_TMP: Surface element connectivity.
    ! NSEL: Number of elements on a given surface.
    ! SEMIN/SEMAX: Partitioned array bounds to distribute surface sections.
    ! STATUS: Confirms if allocation of surface section arrays was successful.
    ! J3/N3: 3D connectivity DOF mapping values.
    ! IERR: Value that confirms if a READ() fails.
    ! EL_DOF_MIN/EL_DOF_MAX: DOF mapping to scatter surface elements.
    ! FASETS: Number of fasets in the field to be parsed.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER, PARAMETER :: TYPE = 6, TYPE1 = TYPE - 1, TYPE3 = 3*TYPE -1
    INTEGER :: ELEM_ID, ELEM_NODES(0:TYPE1), ELEM_NODES_TMP(0:TYPE1)
    INTEGER :: I, J, K, NSEL, SEMIN, SEMAX, STATUS, J3, N3, IERR
    INTEGER :: EL_DOF_MIN, EL_DOF_MAX, FASETS
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$Fasets')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &fasets.')
        !
    END IF
    !
    ! Read in the number of fasets in the section (and set FASETS).
    READ(IO, *) FASETS
    !
    ! This needs to be set here, but is in the scope of SURFACE_MOD
    !   currently.
    NSURFACES = FASETS
    !
    ! Check if 6 surfaces are to be read in. This is hardwired for cubic geom.
    IF (FASETS .NE. 6) CALL PAR_QUIT('Error  :     > &
        &6 surfaces necessary in $Fasets.')
    !
    DO I = 1, FASETS
        !
        ! Read the faset label - unused currently.
        READ(IO,*) FASET(I)
        ! number of elements on the surface.
        READ(IO,*) NSEL
        !
        ! Should SEMIN/SEMAX use the the surface elements from the 3D elements
        ! that are stored on the same processor range instead of something else?
        CALL PAR_PARTITION(NSEL, NUMPROCS, MYID, SEMIN, SEMAX)
        STATUS = ALLOCATE_SURFACE_SECTION(TYPE, SEMIN, SEMAX, SURFACES(I))
        !
        IF (STATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Surface allocation error.')
            !
        ENDIF
        !
        DO J = 0, NSEL-1
            !
            ! For each element, read in the ID and 3D connectivity.
            READ(IO, *) ELEM_ID, ELEM_NODES_TMP
            !
            ! Shift the ID for msh file version > 2.2
            IF (MESH_VERSION .NE. '2.2') THEN
              !
              ELEM_ID = ELEM_ID - EL_STARTID + 1
              !
            END IF
            !
            ! Shift the ID to maintain internal zero-indexing.
            ELEM_ID = ELEM_ID - 1
            !
            IF (ELEM_ID .LT. 0 .OR. ELEM_ID .GE. NUMELM) CALL PAR_QUIT('Error  :     > &
                &Element index out of bounds')
            !
            ! The node values read also need to be shifted and reordered.
            ! GMSH ORDER: 1 2 3 4 5 6
            ! FEPX ORDER: 6 4 2 5 3 1
            ELEM_NODES(0) = ELEM_NODES_TMP(5) - 1
            ELEM_NODES(1) = ELEM_NODES_TMP(2) - 1
            ELEM_NODES(2) = ELEM_NODES_TMP(4) - 1
            ELEM_NODES(3) = ELEM_NODES_TMP(1) - 1
            ELEM_NODES(4) = ELEM_NODES_TMP(3) - 1
            ELEM_NODES(5) = ELEM_NODES_TMP(0) - 1
            !
            IF ((J <= SEMAX) .AND. (J>=SEMIN)) THEN
                !
                ! Create the 2D connectivity first.
                SURFACES(I)%ECONN(0,J) = 6 * ELEM_ID
                SURFACES(I)%ECONN(1,J) = 6 * ELEM_ID + 1
                SURFACES(I)%ECONN(2,J) = 6 * ELEM_ID + 2
                SURFACES(I)%ECONN(3,J) = 6 * ELEM_ID + 3
                SURFACES(I)%ECONN(4,J) = 6 * ELEM_ID + 4
                SURFACES(I)%ECONN(5,J) = 6 * ELEM_ID + 5
                !
                ! Store into nodes and ID into global connectivity?
                SURFACES(I)%CONN(:,J)=ELEM_NODES(:)
                SURFACES(I)%ELEM(J)=ELEM_ID
                !
            END IF
            !
        END DO
        !
        ! Create a 3d connectivity for the coordinate gather.
        !
        DO K = SEMIN, SEMAX
            !
            DO J = 0, TYPE1
                !
                J3 = 3 * J
                N3 = 3 * SURFACES(I)%CONN(J,K)
                SURFACES(I)%CONN3D(J3,K)   = N3
                SURFACES(I)%CONN3D(J3 + 1,K) = N3 + 1
                SURFACES(I)%CONN3D(J3 + 2,K) = N3 + 2
                !
            ENDDO
            !
        ENDDO
        !
        CALL PART_SCATTER_SETUP(0,  TYPE3, DOF_SUB1, DOF_SUP1,&
            & SEMIN, SEMAX, SURFACES(I)%CONN3D, SURFACES(I)%TR)
        EL_DOF_MIN = 6 * EL_SUB1
        EL_DOF_MAX = 6 * EL_SUP1 + 5
        CALL PART_SCATTER_SETUP(0,  5, EL_DOF_MIN, EL_DOF_MAX,&
            & SEMIN, SEMAX, SURFACES(I)%ECONN, SURFACES(I)%ETR)
        !
    ENDDO
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndFasets')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &fasets.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_FASETS
    !
    !===========================================================================
    !
    SUBROUTINE READ_NODEPARTITIONS(IO)
    !
    ! Read in node paritions and do not store - currently unused, but will be
    ! integrated in the future for utilizing imprted SCOTCH partitions.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I: Generic loop index.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! NPARTS: Number of partitions read in from the mesh.
    ! IPART: Partition loop index to check which partition is being checked.
    ! TEMP: Temporary storage for current node ID - not stored.
    ! NODE_PARTS: Array storing the nodal paritions.
    ! TEMPVAL: Array of NP_SUB1/NP_SUP1 bounds for each partition read in.
    ! CHG_INDEX: Index in array immediately after a change in value occurs.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, I, NLINES, NPARTS, IPART, TEMP
    INTEGER, ALLOCATABLE :: NODE_PARTS(:), TEMPVAL(:,:), CHG_INDEX(:)
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$NodePartitions')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &node partitions.')
        !
    END IF
    !
    ! Set the flag if this section is being read in
    READ_NODE_PART = .TRUE.
    !
    ! Read in the number of node partitions in the section (and set NLINES).
    READ(IO, *) NLINES
    !
    ! Allocate the nodal partition array -> N(1,:) is 1-index ID, N(2,:) is part.
    ALLOCATE(NODE_PARTS(1:NLINES))
    !
    ! Loop over the number of NLINES and do NOT store information.
    DO I = 1, NLINES
        !
        READ(IO, *) TEMP, NODE_PARTS(I)
        !
    END DO
    !
    NPARTS = MAXVAL(NODE_PARTS)
    ALLOCATE(TEMPVAL(1:NPARTS, 2))
    ALLOCATE(CHG_INDEX(NPARTS-1))
    IPART = 1
    !
    ! We can split ELEM_PARTS in 0-indexed EL_SUB1 and EL_SUP1 here.
    DO I = 1, SIZE(NODE_PARTS)
        !
        ! A change has been detected so do something
        IF (NODE_PARTS(I) .NE. IPART) THEN
            !
            CHG_INDEX(IPART) = I
            IF (IPART .NE. NPARTS) IPART = IPART + 1
            !
        END IF
        !
    END DO
    !
    ! Build the bound array for NPARTS -> (1,:) is NP_SUB1, (2,:) is NP_SUP1.
    TEMPVAL = 0
    DO I = 1, NPARTS
        !
        ! Handle the edge cases differently
        IF (I .EQ. 1) THEN
            !
            TEMPVAL(1,1) = 1
            TEMPVAL(1,2) = CHG_INDEX(I) - 1
            !
        ELSE IF (I .EQ. NPARTS) THEN
            !
            TEMPVAL(I,1) = CHG_INDEX(I-1)
            TEMPVAL(I,2) = SIZE(NODE_PARTS)
            !
        ELSE
            !
            TEMPVAL(I,1) = CHG_INDEX(I-1)
            TEMPVAL(I,2) = CHG_INDEX(I) - 1
            !
        END IF
        !
    END DO
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndNodePartitions')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &node partitions.')
        !
    END IF
    !
    ! Do not need to maintain node partitions for now so deallocate.
    DEALLOCATE(NODE_PARTS)
    DEALLOCATE(CHG_INDEX)
    DEALLOCATE(TEMPVAL)
    !
    RETURN
    !
    END SUBROUTINE READ_NODEPARTITIONS
    !
    !===========================================================================
    !
    SUBROUTINE READ_PHYSICALNAMES(IO)
    !
    ! Read in physical names and do not store - will never be used as these
    ! are for gmsh internal use only and must be stored to successfully open
    ! the msh file in any version of gmsh.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I: Generic loop index.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, I, NLINES
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$PhysicalNames')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &physical names.')
        !
    END IF
    !
    ! Read in the number of physical names in the section (and set NLINES).
    READ(IO, *) NLINES
    !
    ! Loop over the number of NLINES and do NOT store information.
    DO I = 1, NLINES
        !
        READ(IO, *)
        !
    END DO
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndPhysicalNames')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &physical names.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_PHYSICALNAMES
    !
    !===========================================================================
    !
    SUBROUTINE READ_ELSETORIENTATIONS(IO)
    !
    ! Read elset orientation information from either the msh file or an external
    ! simulation.ori file. If we are reading in orientations from an external
    ! file we must deallocate the current storage array and start over.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I/J: Generic loop index.
    ! S: Status value that read in string is a valid option.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! NSPACE/NVALS: Storage values used for dividing a line into substrings.
    ! ELSET_ID: Grain ID - 1-indexed.
    ! DELIM_POS: Line position where the ':' delimiter is found.
    ! IDEAL_DIM/IDEAL_DIM1: Array widths for storing orientations (3 or 4).
    ! PARM_STRING: String describing input orientation parameterization.
    ! CONV_STRING: String describing input orientation convention.
    ! ORI_STRING: Temp string used to split the 'descriptor:convention' line.
    ! IARRAY: Substring array for internal read parsing.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, I, J, S
    INTEGER :: NLINES, NSPACE, NVALS, ELSET_ID, DELIM_POS
    INTEGER, PARAMETER :: IDEAL_DIM  = 4
    INTEGER, PARAMETER :: IDEAL_DIM1 = IDEAL_DIM - 1
    REAL(RK) :: ORI(0:3)
    CHARACTER(LEN=50)  :: PARM_STRING, CONV_STRING, ORI_STRING
    CHARACTER(LEN=32)  :: IARRAY(16)
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! If per-element orientations have been read-in previously then skip this.
    IF (ELEMENT_ORIS .EQV. .TRUE.) THEN
        !
        ! Read the first line and confirm it is the correct record.
        READ(IO, '(A)', IOSTAT = IERR) LINE
        !
        IF ((IERR .LT. 0) .OR. (LINE .NE. '$ElsetOrientations')) THEN
            !
            CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
                &elset orientations.')
            !
        END IF
        !
        ! Extract the number of lines to skip
        READ(IO, '(A)', IOSTAT = IERR) LINE
        NSPACE = COUNT(TRANSFER(LINE, 'A', LEN_TRIM(LINE)) == " ")
        NVALS = NSPACE + 1
        IARRAY = ""
        !
        ! Internal read of the line into the primary substring arrays.
        READ(LINE, *) IARRAY(1:NVALS)
        READ(IARRAY(1), *) NLINES
        !
        ! Loop over NLINES and do not store read in values.
        DO I = 0, (NLINES - 1)
            !
            READ(IO, *) LINE
            !
        END DO
        !
        ! Read the end of section footer.
        READ(IO, '(A)', IOSTAT = IERR) LINE
        !
        IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndElsetOrientations')) THEN
            !
            CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
                &elset orientations.')
            !
        END IF
        !
        RETURN
        !
    END IF
    !
    ! Check if the array is already allocated from reading the mesh file.
    IF (ALLOCATED(GRAIN_ORIENTATION) .EQV. .TRUE.) THEN
        !
        DEALLOCATE(GRAIN_ORIENTATION)
        !
    END IF
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$ElsetOrientations')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &elset orientations.')
        !
    END IF
    !
    ! Read in the section information string and prepare to parse.
    !
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    ! Trim the full string into two substring arrays split by the space delim.
    NSPACE = COUNT(TRANSFER(LINE, 'A', LEN_TRIM(LINE)) == " ")
    NVALS = NSPACE + 1
    IARRAY = ""
    !
    ! Internal read of the line into the primary substring arrays.
    READ(LINE, *) IARRAY(1:NVALS)
    !
    ! Internal read of the primary substring arrays to get the number of
    ! lines to set to parse and prepare to parse the orientation description.
    READ(IARRAY(1), *) NLINES
    READ(IARRAY(2), '(A)') ORI_STRING
    !
    ! Find the position of the ":" character in ORI_STRING as FORTRAN
    ! refuses to allow the READ() function to handle this automatically.
    DELIM_POS = INDEX(ORI_STRING, ":")
    !
    ! Internal read to store the orientation parameterization and convention.
    READ(ORI_STRING(1:DELIM_POS-1), *) PARM_STRING
    READ(ORI_STRING(DELIM_POS+1:LEN_TRIM(ORI_STRING)), *) CONV_STRING
    !
    ! Set a per-element logical to false.
    ELEMENT_ORIS = .FALSE.
    !
    ! Parse the parameterization and confirm it is a valid option.
    !
    S = 0
    !
    READ(PARM_STRING, '(A)') &
        & ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION
    !
    ! Confirm that the rest of the string is a valid input.
    IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'axis-angle') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'euler-bunge') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'euler-kocks') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'rodrigues') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'quaternion') THEN
        !
        S = 0
        !
    ELSE ! Orientation parameterization defined is incorrect.
        !
        S = 1
        !
    END IF
    !
    IF (S .EQ. 1) THEN
        !
        CALL PAR_QUIT("Error  :     > ORIENTATION_PARAMETERIZATION &
            & contains an error or unexpected input type.")
        !
    END IF
    !
    ! Parse the convention and confirm it is a valid option.
    !
    S = 0
    !
    READ(CONV_STRING, '(A)') &
        & ORIENTATION_OPTIONS%ORIENTATION_CONVENTION
    !
    ! Confirm that the rest of the string is a valid input.
    IF (ORIENTATION_OPTIONS%ORIENTATION_CONVENTION &
        & .EQ. 'active') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_CONVENTION &
        & .EQ. 'passive') THEN
        !
        S = 0
        !
    ELSE ! Orientation convention defined is incorrect.
        !
        S = 1
        !
    END IF
    !
    IF (S .EQ. 1) THEN
        !
        CALL PAR_QUIT("Error  :     > ORIENTATION_CONVENTION &
            & contains an error or unexpected input type.")
        !
    END IF
    !
    ! Check orientation parameterization to set column number.
    IF ( (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'axis-angle') .OR. &
        & (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'quaternion') ) THEN
        !
        ALLOCATE(GRAIN_ORIENTATION(0:IDEAL_DIM1, 0:NLINES - 1))
        !
    ELSE ! All other parameterization have 3 values per orientation.
        !
        ALLOCATE(GRAIN_ORIENTATION(0:IDEAL_DIM1-1, 0:NLINES - 1))
        !
    END IF
    !
    IF (SIZE(GRAIN_ORIENTATION, 1) .EQ. 4) THEN ! Only quaternion or axis-angle.
        !
        DO I = 0, (NLINES - 1)
            !
            READ(IO, *) ELSET_ID, ORI(0:3)
            DO J = 0, 3
                !
                GRAIN_ORIENTATION(J,GRAIN_IDS_INV(ELSET_ID)) = ORI(J)
                !
            END DO
            !
        END DO
        !
    ELSE IF (SIZE(GRAIN_ORIENTATION,1) .EQ. 3) THEN
        !
        DO I = 0, (NLINES - 1)
            !
            READ(IO, *) ELSET_ID, ORI(0:2)
            DO J = 0, 2
                !
                GRAIN_ORIENTATION(J,GRAIN_IDS_INV(ELSET_ID)) = ORI(J)
                !
            END DO
            !
        END DO
        !
    ELSE
        !
        CALL PAR_QUIT("Error  :     > GRAIN_ORIENTATION array not initialized&
            & correctly.")
        !
    END IF
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndElsetOrientations')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &elset orientations.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_ELSETORIENTATIONS
    !
    !===========================================================================
    !
    SUBROUTINE READ_ELEMENTORIENTATIONS(IO)
    !
    ! Read elem orientation information from either the msh file or an external
    ! simulation.ori file. If we are reading in orientations from an external
    ! file we must deallocate the current storage array and start over.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I/J: Generic loop index.
    ! S: Status value that read in string is a valid option.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! NSPACE/NVALS: Storage values used for dividing a line into substrings.
    ! ELT: Element ID - 1-indexed.
    ! DELIM_POS: Line position where the ':' delimiter is found.
    ! IDEAL_DIM/IDEAL_DIM1: Array widths for storing orientations (3 or 4).
    ! PARM_STRING: String describing input orientation parameterization.
    ! CONV_STRING: String describing input orientation convention.
    ! ORI_STRING: Temp string used to split the 'descriptor:convention' line.
    ! IARRAY: Substring array for internal read parsing.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, I, J, S
    INTEGER :: NLINES, NSPACE, NVALS, ELT, DELIM_POS
    INTEGER, PARAMETER :: IDEAL_DIM  = 4
    INTEGER, PARAMETER :: IDEAL_DIM1 = IDEAL_DIM - 1
    CHARACTER(LEN=50)  :: PARM_STRING, CONV_STRING, ORI_STRING
    CHARACTER(LEN=32)  :: IARRAY(16)
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Check if the array is already allocated from reading the mesh file.
    IF (ALLOCATED(GRAIN_ORIENTATION) .EQV. .TRUE.) THEN
        !
        DEALLOCATE(GRAIN_ORIENTATION)
        !
    END IF

    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$ElementOrientations')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &element orientations.')
        !
    END IF
    !
    ! Read in the section information string and prepare to parse.
    !
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    ! Trim the full string into two substring arrays split by the space delim.
    NSPACE = COUNT(TRANSFER(LINE, 'A', LEN_TRIM(LINE)) == " ")
    NVALS = NSPACE + 1
    IARRAY = ""
    !
    ! Internal read of the line into the primary substring arrays.
    READ(LINE, *) IARRAY(1:NVALS)
    !
    ! Internal read of the primary substring arrays to get the number of
    ! lines to set to parse and prepare to parse the orientation description.
    READ(IARRAY(1), *) NLINES
    READ(IARRAY(2), '(A)') ORI_STRING
    !
    ! Find the position of the ":" character in ORI_STRING as FORTRAN
    ! refuses to allow the READ() function to handle this automatically.
    DELIM_POS = INDEX(ORI_STRING, ":")
    !
    ! Internal read to store the orientation parameterization and convention.
    READ(ORI_STRING(1:DELIM_POS-1), *) PARM_STRING
    READ(ORI_STRING(DELIM_POS+1:LEN_TRIM(ORI_STRING)), *) CONV_STRING
    !
    ! Set a per-element logical to true.
    ELEMENT_ORIS = .TRUE.
    !
    ! Parse the parameterization and confirm it is a valid option.
    !
    S = 0
    !
    READ(PARM_STRING, '(A)') &
        & ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION
    !
    ! Confirm that the rest of the string is a valid input.
    IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'axis-angle') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'euler-bunge') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'euler-kocks') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'rodrigues') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'quaternion') THEN
        !
        S = 0
        !
    ELSE ! Orientation parameterization defined is incorrect.
        !
        S = 1
        !
    END IF
    !
    IF (S .EQ. 1) THEN
        !
        CALL PAR_QUIT("Error  :     > ORIENTATION_PARAMETERIZATION &
            & contains an error or unexpected input type.")
        !
    END IF
    !
    ! Parse the convention and confirm it is a valid option.
    !
    S = 0
    !
    READ(CONV_STRING, '(A)') &
        & ORIENTATION_OPTIONS%ORIENTATION_CONVENTION
    !
    ! Confirm that the rest of the string is a valid input.
    IF (ORIENTATION_OPTIONS%ORIENTATION_CONVENTION &
        & .EQ. 'active') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_CONVENTION &
        & .EQ. 'passive') THEN
        !
        S = 0
        !
    ELSE ! Orientation convention defined is incorrect.
        !
        S = 1
        !
    END IF
    !
    IF (S .EQ. 1) THEN
        !
        CALL PAR_QUIT("Error  :     > ORIENTATION_CONVENTION &
            & contains an error or unexpected input type.")
        !
    END IF
    !
    ! Check orientation parameterization to set column number.
    IF ( (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'axis-angle') .OR. &
        & (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'quaternion') ) THEN
        !
        ALLOCATE(GRAIN_ORIENTATION(0:IDEAL_DIM1, 0:NLINES - 1))
        !
    ELSE ! All other parameterization have 3 values per orientation.
        !
        ALLOCATE(GRAIN_ORIENTATION(0:IDEAL_DIM1-1, 0:NLINES - 1))
        !
    END IF
    !
    IF (SIZE(GRAIN_ORIENTATION, 1) .EQ. 4) THEN ! Only quaternion or axis-angle.
        !
        DO I = 0, (NLINES - 1)
            !
            READ(IO, *) ELT, (GRAIN_ORIENTATION(J,I), J = 0, IDEAL_DIM1)
            !
        END DO
        !
    ELSE IF (SIZE(GRAIN_ORIENTATION,1) .EQ. 3) THEN
        !
        DO I = 0, (NLINES - 1)
            !
            READ(IO, *) ELT, (GRAIN_ORIENTATION(J,I), J = 0, IDEAL_DIM1-1)
            !
        END DO
        !
    ELSE
        !
        CALL PAR_QUIT("Error  :     > GRAIN_ORIENTATION array not initialized&
            & correctly.")
        !
    END IF
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndElementOrientations')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &orientations.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_ELEMENTORIENTATIONS
    !
    !===========================================================================
    !
    SUBROUTINE READ_GROUPS(IO, CONSOLE_UNIT)
    !
    ! Read grain and phase information from either the .msh file or an external
    ! simulation.phase file. If we are reading in phases from an external
    ! file we must deallocate the current storage array and start over.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    INTEGER, INTENT(IN) :: CONSOLE_UNIT
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I: Generic loop index.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! ELSET_ID: Grain ID - 1-indexed.
    ! ELSET_PHASE: Storage array to hold which phase is assigned to which grain.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: I, IERR, NLINES, ELSET_ID, PHASE_TMP
    INTEGER, ALLOCATABLE :: ELSET_PHASE(:)
    CHARACTER(LEN=256)   :: LINE
    !
    ! Notes:
    ! If the simulation is single-phase then we don't input any information
    ! on the grain phases anymore. Therefore we need to attempt to read in
    ! the next record on the mesh unit and see if it fails by EOF. This only
    ! works as $Groups will always be the last section in the mesh file iff
    ! it is present and generated via Neper.
    !
    ! Todo:
    ! It seems to be possible to define a number of grain/phase assignments
    ! that is out of bounds of the previously allocated array since we don't
    ! explicity enforce this.
    !----------------------------------------------------------------------
    !
    ! Attempt to read the next line (may be EOF, but unknown)
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF (IERR .LT. 0) THEN
        !
        IF (MYID .EQ. 0) THEN
            !
            WRITE(CONSOLE_UNIT, '(A)') 'Info   :   [i] No phase assignments &
                &found in `simulation.msh`.'
            !
        END IF
        !
        RETURN
        !
    ELSE IF (IERR .EQ. 0) THEN
        !
        ! Check that the line read in previously was $Groups.
        IF (LINE .NE. '$Groups') THEN
            !
            CALL PAR_QUIT('Error  :     > Phase assignment section header &
                &should be `$Groups`.')
            !
        END IF
        !
        ! Phase assignments are per elset so we need to do some mapping with
        ! temporary arrays.
        READ(IO, *) ! This skips the `elset` line
        READ(IO, *) NLINES
        !
        ! Check if the array is already allocated from the mesh read in
        IF (ALLOCATED(ELSET_PHASE) .EQV. .TRUE.) THEN
            !
            DEALLOCATE(ELSET_PHASE)
            !
        END IF
        !
        ! Allocate and fill an array with elset-based `grain/phase` pairs
        ! (1,:) is phase assign to grain ID (row index)
        ALLOCATE(ELSET_PHASE(0:NLINES - 1))
        !
        DO I = 0, NLINES - 1
            !
            READ(IO, *) ELSET_ID, PHASE_TMP
            !
            ELSET_PHASE(GRAIN_IDS_INV(ELSET_ID)) = PHASE_TMP
            !
        END DO
        !
        ! Testing ELSET_PHASE
        IF (MINVAL (ELSET_PHASE) .LE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > ELSET_PHASE is out of bounds.')
            !
        END IF
        !
        ! Now, loop over UESUB and reassign phase accordingly
        DO I = 0, NUMELM-1
            !
            ! Assign phase per element based on elset value of UESUB per element
            ! UESUB(I) returns an elset value that is 0-indexed that will grab
            ! the corresponding elset row in ELSET_PHASE and, in turn, the phase
            PHASE(I) = ELSET_PHASE(UESUB(I))

            !
        END DO
        !
        ! Testing PHASE
        IF (MINVAL (PHASE) .LE. 0)  THEN
            !
            CALL PAR_QUIT('Error  :     > Phase index is out of bounds.')
            !
        END IF
        !
        ! Now that PHASE is built there is a chance that the user supplied
        ! an incorrect `number_of_phases` parameter for what the .msh has
        IF (MAXVAL(ELSET_PHASE) .NE. CRYS_OPTIONS%NUMBER_OF_PHASES) THEN
            !
            CALL PAR_QUIT('Error  :     > Number of phases in mesh does not&
                & match the config file.')
            !
        END IF
        !
    ELSE
        ! If something goes wrong terminate the simulation.
        CALL PAR_QUIT('Error  :     > Phase assignment parsing failed.')
        !
    END IF
    !
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndGroups')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &groups.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_GROUPS
    !
    !===========================================================================
    !
    SUBROUTINE INITIALIZE_OPTIONS()
    !
    ! Initialize inputs and assign default values.
    !
    !---------------------------------------------------------------------------
    !
    ! Standard options
    OPTIONS%MAX_INCR = 50000
    OPTIONS%MAX_TOTAL_TIME = 12000.0D0
    OPTIONS%DEF_CONTROL_BY = UNIAXIAL_STRAIN_TARGET
    OPTIONS%CHECK_NECKING = .FALSE.
    OPTIONS%LOAD_TOL = 0.0
    OPTIONS%DTIME_FACTOR = 1.001D0
    OPTIONS%HARD_TYPE = 'isotropic' ! set what type of hardening module used
    OPTIONS%RESTART = .FALSE.
    OPTIONS%RESTART_FILE_HANDLING = -1
    OPTIONS%RESTART_INITIAL_STEP = -1
    OPTIONS%READ_ORI_FROM_FILE = .FALSE.
    OPTIONS%SAT_EVO = .FALSE.
    OPTIONS%PRECIP_HARD = .FALSE.
    OPTIONS%ORI_FILE = ''
    OPTIONS%READ_PHASE_FROM_FILE = .FALSE.
    OPTIONS%PHASE_FILE = ''
    OPTIONS%MAX_ITER_HARD_LIMIT = 10
    !
    ! Printing options
    PRINT_OPTIONS%PRINT_COORDINATES = .FALSE.
    PRINT_OPTIONS%PRINT_DISPLACEMENTS = .FALSE.
    PRINT_OPTIONS%PRINT_VELOCITIES = .FALSE.
    PRINT_OPTIONS%PRINT_ORIENTATIONS = .FALSE.
    PRINT_OPTIONS%PRINT_CRSS = .FALSE.
    PRINT_OPTIONS%PRINT_STRAIN_TOT = .FALSE.
    PRINT_OPTIONS%PRINT_STRAIN_EL = .FALSE.
    PRINT_OPTIONS%PRINT_STRAIN_PL = .FALSE.
    PRINT_OPTIONS%PRINT_STRESS = .FALSE.
    PRINT_OPTIONS%PRINT_GAMMADOT = .FALSE.
    PRINT_OPTIONS%PRINT_DEFF = .FALSE.
    PRINT_OPTIONS%PRINT_DPEFF = .FALSE.
    PRINT_OPTIONS%PRINT_ELVOL = .FALSE.
    PRINT_OPTIONS%PRINT_EQSTRESS = .FALSE.
    PRINT_OPTIONS%PRINT_EQSTRAIN = .FALSE.
    PRINT_OPTIONS%PRINT_EQELSTRAIN = .FALSE.
    PRINT_OPTIONS%PRINT_EQPLSTRAIN = .FALSE.
    PRINT_OPTIONS%PRINT_VGRAD = .FALSE.
    PRINT_OPTIONS%PRINT_DPHAT = .FALSE.
    PRINT_OPTIONS%PRINT_WPHAT = .FALSE.
    PRINT_OPTIONS%PRINT_GAMMA = .FALSE.
    PRINT_OPTIONS%PRINT_RESTART = .FALSE.
    PRINT_OPTIONS%PRINT_FORCES = .FALSE.
    PRINT_OPTIONS%PRINT_CONV = .FALSE.
    PRINT_OPTIONS%PRINT_WORK = .FALSE.
    PRINT_OPTIONS%PRINT_WORKP = .FALSE.
    PRINT_OPTIONS%PRINT_DEFRATE = .FALSE.
    PRINT_OPTIONS%PRINT_WORKRATE = .FALSE.
    PRINT_OPTIONS%PRINT_WORKRATEP = .FALSE.
    !
    ! Boundary condition options
    BCS_OPTIONS%READ_BCS_FROM_FILE = .FALSE.
    BCS_OPTIONS%BCS_FILE = ''
    BCS_OPTIONS%BOUNDARY_CONDITIONS = 0
    BCS_OPTIONS%LOADING_FACE = 0
    BCS_OPTIONS%LOADING_DIRECTION = 0
    !
    ! Crystal parameters
    CRYS_OPTIONS%NUMBER_OF_PHASES = 0
    CRYS_OPTIONS%PHASE = 0
    !
    ! Uniaxial control options
    UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_STEPS = 0
    UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS = 0
    UNIAXIAL_OPTIONS%NUMBER_OF_LOAD_STEPS = 0
    !
    ! Triaxial CSR options
    TRIAXCSR_OPTIONS%NUMBER_OF_CSR_LOAD_STEPS = 0
    TRIAXCSR_OPTIONS%MAX_BC_ITER = 10
    TRIAXCSR_OPTIONS%MIN_PERT_FRAC = 1.0D-3
    TRIAXCSR_OPTIONS%LOAD_TOL_ABS = 0.1D0
    TRIAXCSR_OPTIONS%LOAD_TOL_REL = 0.001D0
    TRIAXCSR_OPTIONS%MAX_STRAIN = 0.2D0
    TRIAXCSR_OPTIONS%MAX_EQSTRAIN = 0.2D0
    !
    ! Triaxial CLR options
    TRIAXCLR_OPTIONS%NUMBER_OF_CLR_LOAD_STEPS = 0
    TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS = 0
    TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES = 0
    TRIAXCLR_OPTIONS%MAX_BC_ITER = 10
    TRIAXCLR_OPTIONS%LOAD_TOL_ABS = 0.1D0
    TRIAXCLR_OPTIONS%MAX_STRAIN_INCR = 0.001D0
    TRIAXCLR_OPTIONS%MAX_STRAIN = 0.2D0
    TRIAXCLR_OPTIONS%MAX_EQSTRAIN = 0.2D0
    !
    END SUBROUTINE INITIALIZE_OPTIONS
    !
    !===========================================================================
    !
    LOGICAL FUNCTION PROCESS_INPUT(CMDLINE, INUNIT, STATUS)
    !
    !  Process options and input.
    !
    !---------------------------------------------------------------------------
    !
    CHARACTER(LEN=*), INTENT(IN) :: CMDLINE
    INTEGER, INTENT(IN)  :: INUNIT
    INTEGER, INTENT(OUT) :: STATUS
    !
    !---------------------------------------------------------------------------
    !
    STATUS = 0; PROCESS_INPUT = .TRUE.
    !
    PRIV_INPUNIT = INUNIT ! make available to callbacks
    !
    ! Standard options
    IF (EXECCOMMAND('max_incr', CMDLINE, &
        & EXEC_MAX_INCR, STATUS)) RETURN
    IF (EXECCOMMAND('max_total_time', CMDLINE, &
        & EXEC_MAX_TOTAL_TIME, STATUS)) RETURN
    IF (EXECCOMMAND('def_control_by', CMDLINE, &
        & EXEC_DEF_CONTROL_BY, STATUS)) RETURN
    IF (EXECCOMMAND('check_necking', CMDLINE, &
        & EXEC_CHECK_NECKING, STATUS)) RETURN
    IF (EXECCOMMAND('load_tol', CMDLINE, &
        & EXEC_LOAD_TOL, STATUS)) RETURN
    IF (EXECCOMMAND('dtime_factor', CMDLINE, &
        & EXEC_DTIME_FACTOR, STATUS)) RETURN
    IF (EXECCOMMAND('restart', CMDLINE, &
        & EXEC_RESTART, STATUS)) RETURN
    IF (EXECCOMMAND('hard_type', CMDLINE, &
        & EXEC_HARD_TYPE, STATUS)) RETURN
    IF (EXECCOMMAND('read_ori_from_file', CMDLINE, &
        & EXEC_READ_ORI_FROM_FILE, STATUS)) RETURN
    IF (EXECCOMMAND('read_phase_from_file', CMDLINE, &
        & EXEC_READ_PHASE_FROM_FILE, STATUS)) RETURN
    IF (EXECCOMMAND('max_iter_hard_limit', CMDLINE, &
        & EXEC_MAX_ITER_HARD_LIMIT, STATUS)) RETURN
    !
    ! Printing options
    IF (EXECCOMMAND('print', CMDLINE, &
        & EXEC_PRINT, STATUS)) RETURN
    IF (EXECCOMMAND('suppress', CMDLINE, &
        & EXEC_SUPPRESS, STATUS)) RETURN
    !
    ! Boundary conditions options
    IF (EXECCOMMAND('read_bcs_from_file', CMDLINE, &
        & EXEC_READ_BCS_FROM_FILE, STATUS)) RETURN
    IF (EXECCOMMAND('boundary_conditions', CMDLINE, &
        & EXEC_BOUNDARY_CONDITIONS, STATUS)) RETURN
    IF (EXECCOMMAND('loading_face', CMDLINE, &
        & EXEC_LOADING_FACE, STATUS)) RETURN
    IF (EXECCOMMAND('loading_direction', CMDLINE, &
        & EXEC_LOADING_DIRECTION, STATUS)) RETURN
    IF (EXECCOMMAND('strain_rate', CMDLINE, &
        & EXEC_STRAIN_RATE, STATUS)) RETURN
    IF (EXECCOMMAND('load_rate', CMDLINE, &
        & EXEC_LOAD_RATE, STATUS)) RETURN
    !
    ! Crystal parameter options
    IF (EXECCOMMAND('number_of_phases', CMDLINE, &
        & EXEC_NUMBER_OF_PHASES, STATUS)) RETURN
    IF (EXECCOMMAND('crystal_type', CMDLINE, &
        & EXEC_CRYSTAL_TYPE, STATUS)) RETURN
    IF (EXECCOMMAND('phase', CMDLINE, &
        & EXEC_PHASE, STATUS)) RETURN
    IF (EXECCOMMAND('m', CMDLINE, &
        & EXEC_M, STATUS)) RETURN
    IF (EXECCOMMAND('gammadot_0', CMDLINE, &
        & EXEC_GAMMADOT_0, STATUS)) RETURN
    IF (EXECCOMMAND('h_0', CMDLINE, &
        & EXEC_H_0, STATUS)) RETURN
    IF (EXECCOMMAND('g_0', CMDLINE, &
        & EXEC_G_0, STATUS)) RETURN
    IF (EXECCOMMAND('g_s0', CMDLINE, &
        & EXEC_G_S0, STATUS)) RETURN
    IF (EXECCOMMAND('m_prime', CMDLINE, &
        & EXEC_M_PRIME, STATUS)) RETURN
    IF (EXECCOMMAND('gammadot_s0', CMDLINE, &
        & EXEC_GAMMADOT_S0, STATUS)) RETURN
    IF (EXECCOMMAND('n', CMDLINE, &
        & EXEC_N, STATUS)) RETURN
    IF (EXECCOMMAND('c11', CMDLINE, &
        & EXEC_C11, STATUS)) RETURN
    IF (EXECCOMMAND('c12', CMDLINE, &
        & EXEC_C12, STATUS)) RETURN
    IF (EXECCOMMAND('c13', CMDLINE, &
        & EXEC_C13, STATUS)) RETURN
    IF (EXECCOMMAND('c44', CMDLINE, &
        & EXEC_C44, STATUS)) RETURN
    IF (EXECCOMMAND('c66', CMDLINE, &
        & EXEC_C66, STATUS)) RETURN
    IF (EXECCOMMAND('c_over_a', CMDLINE, &
        & EXEC_C_OVER_A, STATUS)) RETURN
    IF (EXECCOMMAND('cyclic_a', CMDLINE, &
        & EXEC_CYCLIC_A, STATUS)) RETURN
    IF (EXECCOMMAND('cyclic_c', CMDLINE, &
        & EXEC_CYCLIC_C, STATUS)) RETURN
    IF (EXECCOMMAND('latent_parameters', CMDLINE, &
        & EXEC_LATENT_PARAMETERS, STATUS)) RETURN
    IF (EXECCOMMAND('a_p', CMDLINE, &
        & EXEC_A_P, STATUS)) RETURN
    IF (EXECCOMMAND('f_p', CMDLINE, &
        & EXEC_F_P, STATUS)) RETURN
    IF (EXECCOMMAND('b_p', CMDLINE, &
        & EXEC_B_P, STATUS)) RETURN
    IF (EXECCOMMAND('r_p', CMDLINE, &
        & EXEC_R_P, STATUS)) RETURN
    !
    ! Uniaxial control options
    IF (EXECCOMMAND('number_of_strain_steps', CMDLINE, &
        & EXEC_NUMBER_OF_STRAIN_STEPS, STATUS)) RETURN
    IF (EXECCOMMAND('number_of_strain_rate_jumps', CMDLINE, &
        & EXEC_NUMBER_OF_STRAIN_RATE_JUMPS, STATUS)) RETURN
    IF (EXECCOMMAND('target_strain', CMDLINE, &
        & EXEC_TARGET_STRAIN, STATUS)) RETURN
    IF (EXECCOMMAND('strain_rate_jump', CMDLINE, &
        & EXEC_STRAIN_RATE_JUMP, STATUS)) RETURN
    IF (EXECCOMMAND('number_of_load_steps', CMDLINE, &
        & EXEC_NUMBER_OF_LOAD_STEPS, STATUS)) RETURN
    IF (EXECCOMMAND('target_load', CMDLINE, &
        & EXEC_TARGET_LOAD, STATUS)) RETURN
    !
    ! Triaxial CSR options
    IF (EXECCOMMAND('max_bc_iter', CMDLINE, &
        & EXEC_MAX_BC_ITER, STATUS)) RETURN
    IF (EXECCOMMAND('min_pert_frac', CMDLINE, &
        & EXEC_MIN_PERT_FRAC, STATUS)) RETURN
    IF (EXECCOMMAND('load_tol_abs', CMDLINE, &
        & EXEC_LOAD_TOL_ABS, STATUS)) RETURN
    IF (EXECCOMMAND('load_tol_rel', CMDLINE, &
        & EXEC_LOAD_TOL_REL, STATUS)) RETURN
    IF (EXECCOMMAND('number_of_csr_load_steps', CMDLINE, &
        & EXEC_NUMBER_OF_CSR_LOAD_STEPS, STATUS)) RETURN
    IF (EXECCOMMAND('target_csr_load', CMDLINE, &
        & EXEC_TARGET_CSR_LOAD, STATUS)) RETURN
    !
    ! Triaxial CLR options
    IF (EXECCOMMAND('max_strain_incr', CMDLINE, &
        & EXEC_MAX_STRAIN_INCR, STATUS)) RETURN
    IF (EXECCOMMAND('max_strain', CMDLINE, &
        & EXEC_MAX_STRAIN, STATUS)) RETURN
    IF (EXECCOMMAND('max_eqstrain', CMDLINE, &
        & EXEC_MAX_EQSTRAIN, STATUS)) RETURN
    IF (EXECCOMMAND('number_of_clr_load_steps', CMDLINE, &
        & EXEC_NUMBER_OF_CLR_LOAD_STEPS, STATUS)) RETURN
    IF (EXECCOMMAND('target_clr_load', CMDLINE, &
        & EXEC_TARGET_CLR_LOAD, STATUS)) RETURN
    IF (EXECCOMMAND('number_of_load_rate_jumps', CMDLINE, &
        & EXEC_NUMBER_OF_LOAD_RATE_JUMPS, STATUS)) RETURN
    IF (EXECCOMMAND('load_rate_jump', CMDLINE, &
        & EXEC_LOAD_RATE_JUMP, STATUS)) RETURN
    IF (EXECCOMMAND('number_of_dwell_episodes', CMDLINE, &
        & EXEC_NUMBER_OF_DWELL_EPISODES, STATUS)) RETURN
    IF (EXECCOMMAND('dwell_episode', CMDLINE, &
        & EXEC_DWELL_EPISODE, STATUS)) RETURN
    !
    IF (CONVERGENCEKEYWORDINPUT(CMDLINE, STATUS)) RETURN
    !
    ! No calls matched the keyword
    !
    PROCESS_INPUT = .FALSE.
    !
    END FUNCTION PROCESS_INPUT
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_MAX_INCR(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) OPTIONS%MAX_INCR
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_MAX_INCR
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_MAX_TOTAL_TIME(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) OPTIONS%MAX_TOTAL_TIME
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_MAX_TOTAL_TIME
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_DEF_CONTROL_BY(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    !
    SELECT CASE (TRIM(ADJUSTL(A)))
        !
        CASE ('uniaxial_load_target')
            !
            OPTIONS%DEF_CONTROL_BY = UNIAXIAL_LOAD_TARGET
            !
        CASE ('uniaxial_strain_target')
            !
            OPTIONS%DEF_CONTROL_BY = UNIAXIAL_STRAIN_TARGET
            !
        CASE ('triaxial_constant_strain_rate')
            !
            OPTIONS%DEF_CONTROL_BY = TRIAXIAL_CONSTANT_STRAIN_RATE
            !
        CASE ('triaxial_constant_load_rate')
            !
            OPTIONS%DEF_CONTROL_BY = TRIAXIAL_CONSTANT_LOAD_RATE
            !
        CASE DEFAULT
            !
            S = 1
            !
        ! END CASE
        !
    END SELECT
    !
    END SUBROUTINE EXEC_DEF_CONTROL_BY
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_CHECK_NECKING(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    !
    SELECT CASE (TRIM(ADJUSTL(A)))
        !
        CASE ('on','On','ON','true','True','TRUE')
            !
            OPTIONS%CHECK_NECKING = .TRUE.
            !
        CASE ('off','Off','OFF','false','False','FALSE')
            !
            OPTIONS%CHECK_NECKING = .FALSE.
            !
        CASE DEFAULT
            !
            S = 1
            !
        ! END CASE
        !
    END SELECT
    !
    END SUBROUTINE EXEC_CHECK_NECKING
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_LOAD_TOL(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) OPTIONS%LOAD_TOL
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE IF ((IOERR .EQ. 0) .AND. (OPTIONS%LOAD_TOL .GE. 0.0)) THEN
        !
        S = 0
        !
    ELSE
        !
        S = 1
        !
    END IF
    !
    END SUBROUTINE EXEC_LOAD_TOL
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_DTIME_FACTOR(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) OPTIONS%DTIME_FACTOR
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE IF ((IOERR .EQ. 0) .AND. (OPTIONS%DTIME_FACTOR .GE. 1.0)) THEN
        !
        S = 0
        !
    ELSE
        !
        S = 1
        !
    END IF
    !
    END SUBROUTINE EXEC_DTIME_FACTOR
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_RESTART(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    CHARACTER(LEN=128) :: FIELD_NAME
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *) FIELD_NAME
    !
    S = 0
    !
    ! Note that (0) denotes file appending while (1) denotes new files will
    ! be written. Appending will also create new files if they are not present.
    !
    SELECT CASE (TRIM(ADJUSTL(FIELD_NAME)))
        !
        CASE('append', 'Append', 'APPEND')
            !
            CALL PAR_QUIT('Error  :     > Unsupported restart option. Use &
                &"restart new_file".')
            !OPTIONS%RESTART = .TRUE.
            !OPTIONS%RESTART_FILE_HANDLING = 0
            !
        CASE('on', 'On', 'ON', 'new_file', 'New_file', 'New_File', 'NEW_FILE')
            !
            OPTIONS%RESTART = .TRUE.
            OPTIONS%RESTART_FILE_HANDLING = 1
            !
        CASE DEFAULT
            !
            S = 1
            !
        ! END CASE
        !
    END SELECT
    !
    END SUBROUTINE EXEC_RESTART
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_HARD_TYPE(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    CHARACTER(LEN=128) :: FIELD_NAME
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *) FIELD_NAME
    !
    S = 0
    !
    SELECT CASE (TRIM(ADJUSTL(FIELD_NAME)))
        !
        CASE('ISOTROPIC','Isotropic','isotropic')
            !
            OPTIONS%HARD_TYPE = 'isotropic'
            !
        CASE('LATENT','Latent','latent','ANISOTROPIC', &
            & 'Anisotropic','anisotropic')
            !
            OPTIONS%HARD_TYPE = 'latent'
            !
        CASE('test','TEST','Test')
            !
            OPTIONS%HARD_TYPE = 'test'
            !
        CASE('cyclic_isotropic','Cyclic_Isotropic','CYCLIC_ISOTROPIC')
            !
            OPTIONS%HARD_TYPE = 'cyclic_isotropic'
            !
        CASE('cyclic_anisotropic','Cyclic_Anisotropic','CYCLIC_ANISOTROPIC')
            !
            OPTIONS%HARD_TYPE = 'cyclic_anisotropic'
            !
        CASE DEFAULT
            !
            S = 1
            !
        ! END CASE
        !
    END SELECT
    !
    END SUBROUTINE EXEC_HARD_TYPE
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_READ_ORI_FROM_FILE(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    ! The external file must be named 'simulation.ori'.
    OPTIONS%READ_ORI_FROM_FILE = .TRUE.
    !
    S = 0
    !
    END SUBROUTINE EXEC_READ_ORI_FROM_FILE
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_READ_PHASE_FROM_FILE(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    ! The external file must be named 'simulation.phase'.
    OPTIONS%READ_PHASE_FROM_FILE = .TRUE.
    !
    S = 0
    !
    END SUBROUTINE EXEC_READ_PHASE_FROM_FILE
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_PRINT(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    CHARACTER(LEN=128) :: FIELD_NAME
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    !
    READ(A, *) FIELD_NAME
    !
    SELECT CASE (TRIM(ADJUSTL(FIELD_NAME)))
        !
        CASE ('coo')
            !
            PRINT_OPTIONS%PRINT_COORDINATES = .TRUE.
            !
        CASE ('disp')
            !
            PRINT_OPTIONS%PRINT_DISPLACEMENTS = .TRUE.
            !
        CASE ('vel')
            !
            PRINT_OPTIONS%PRINT_VELOCITIES = .TRUE.
            !
        CASE ('ori')
            !
            PRINT_OPTIONS%PRINT_ORIENTATIONS = .TRUE.
            !
        CASE ('crss')
            !
            PRINT_OPTIONS%PRINT_CRSS = .TRUE.
            !
        CASE ('strain')
            !
            PRINT_OPTIONS%PRINT_STRAIN_TOT = .TRUE.
            !
        CASE ('strain-el')
            !
            PRINT_OPTIONS%PRINT_STRAIN_EL = .TRUE.
            !
        CASE ('strain-pl')
            !
            PRINT_OPTIONS%PRINT_STRAIN_PL = .TRUE.
            !
        CASE ('stress')
            !
            PRINT_OPTIONS%PRINT_STRESS = .TRUE.
            !
        CASE ('sliprate')
            !
            PRINT_OPTIONS%PRINT_GAMMADOT = .TRUE.
            !
        CASE ('defrate-eq')
            !
            PRINT_OPTIONS%PRINT_DEFF = .TRUE.
            !
        CASE ('defrate-pl-eq')
            !
            PRINT_OPTIONS%PRINT_DPEFF = .TRUE.
            !
        CASE ('elt-vol')
            !
            PRINT_OPTIONS%PRINT_ELVOL = .TRUE.
            !
        CASE ('stress-eq')
            !
            PRINT_OPTIONS%PRINT_EQSTRESS = .TRUE.
            !
        CASE ('strain-eq')
            !
            PRINT_OPTIONS%PRINT_EQSTRAIN = .TRUE.
            !
        CASE ('strain-el-eq')
            !
            PRINT_OPTIONS%PRINT_EQELSTRAIN = .TRUE.
            !
        CASE ('strain-pl-eq')
            !
            PRINT_OPTIONS%PRINT_EQPLSTRAIN = .TRUE.
            !
        CASE ('velgrad')
            !
            PRINT_OPTIONS%PRINT_VGRAD = .TRUE.
            !
        CASE ('defrate-pl')
            !
            PRINT_OPTIONS%PRINT_DPHAT = .TRUE.
            !
        CASE ('spinrate')
            !
            PRINT_OPTIONS%PRINT_WPHAT = .TRUE.
            !
        CASE ('slip')
            !
            PRINT_OPTIONS%PRINT_GAMMA = .TRUE.
            !
        CASE ('restart')
            !
            PRINT_OPTIONS%PRINT_RESTART = .TRUE.
            !
        CASE ('forces')
            !
            PRINT_OPTIONS%PRINT_FORCES = .TRUE.
            !
        CASE ('convergence')
            !
            PRINT_OPTIONS%PRINT_CONV = .TRUE.
            !
        CASE ('work')
            !
            PRINT_OPTIONS%PRINT_WORK = .TRUE.
            !
        CASE ('work-pl')
            !
            PRINT_OPTIONS%PRINT_WORKP = .TRUE.
            !
        CASE ('defrate')
            !
            PRINT_OPTIONS%PRINT_DEFRATE = .TRUE.
            !
        CASE ('workrate')
            !
            PRINT_OPTIONS%PRINT_WORKRATE = .TRUE.
            !
        CASE ('workrate-pl')
            !
            PRINT_OPTIONS%PRINT_WORKRATEP = .TRUE.
            !
        CASE DEFAULT
            !
            S = 1
            !
        ! END CASE
        !
    END SELECT
    !
    END SUBROUTINE EXEC_PRINT
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_SUPPRESS(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    CHARACTER(LEN=128) :: FIELD_NAME
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    !
    READ(A, *) FIELD_NAME
    !
    SELECT CASE (TRIM(ADJUSTL(FIELD_NAME)))
        !
        CASE ('coo')
            !
            PRINT_OPTIONS%PRINT_COORDINATES = .FALSE.
            !
        CASE ('disp')
            !
            PRINT_OPTIONS%PRINT_DISPLACEMENTS = .FALSE.
            !
        CASE ('vel')
            !
            PRINT_OPTIONS%PRINT_VELOCITIES = .FALSE.
            !
        CASE ('ori')
            !
            PRINT_OPTIONS%PRINT_ORIENTATIONS = .FALSE.
            !
        CASE ('crss')
            !
            PRINT_OPTIONS%PRINT_CRSS = .FALSE.
            !
        CASE ('strain')
            !
            PRINT_OPTIONS%PRINT_STRAIN_TOT = .FALSE.
            !
        CASE ('strain-el')
            !
            PRINT_OPTIONS%PRINT_STRAIN_EL = .FALSE.
            !
        CASE ('strain-pl')
            !
            PRINT_OPTIONS%PRINT_STRAIN_PL = .FALSE.
            !
        CASE ('stress')
            !
            PRINT_OPTIONS%PRINT_STRESS = .FALSE.
            !
        CASE ('sliprate')
            !
            PRINT_OPTIONS%PRINT_GAMMADOT = .FALSE.
            !
        CASE ('defrate-eq')
            !
            PRINT_OPTIONS%PRINT_DEFF = .FALSE.
            !
        CASE ('defrate-pl-eq')
            !
            PRINT_OPTIONS%PRINT_DPEFF = .FALSE.
            !
        CASE ('stress-eq')
            !
            PRINT_OPTIONS%PRINT_EQSTRESS = .FALSE.
            !
        CASE ('strain-eq')
            !
            PRINT_OPTIONS%PRINT_EQSTRAIN = .FALSE.
            !
        CASE ('strain-el-eq')
            !
            PRINT_OPTIONS%PRINT_EQELSTRAIN = .FALSE.
            !
        CASE ('strain-pl-eq')
            !
            PRINT_OPTIONS%PRINT_EQPLSTRAIN = .FALSE.
            !
        CASE ('velgrad')
            !
            PRINT_OPTIONS%PRINT_VGRAD = .FALSE.
            !
        CASE ('defrate-pl')
            !
            PRINT_OPTIONS%PRINT_DPHAT = .FALSE.
            !
        CASE ('elt-vol')
            !
            PRINT_OPTIONS%PRINT_ELVOL = .FALSE.
            !
        CASE ('spinrate')
            !
            PRINT_OPTIONS%PRINT_WPHAT = .FALSE.
            !
        CASE ('slip')
            !
            PRINT_OPTIONS%PRINT_GAMMA = .FALSE.
            !
        CASE ('restart')
            !
            PRINT_OPTIONS%PRINT_RESTART = .FALSE.
            !
        CASE ('forces')
            !
            PRINT_OPTIONS%PRINT_FORCES = .FALSE.
            !
        CASE ('convergence')
            !
            PRINT_OPTIONS%PRINT_CONV = .FALSE.
            !
        CASE ('work')
            !
            PRINT_OPTIONS%PRINT_WORK = .FALSE.
            !
        CASE ('work-pl')
            !
            PRINT_OPTIONS%PRINT_WORKP = .FALSE.
            !
        CASE ('defrate')
            !
            PRINT_OPTIONS%PRINT_DEFRATE = .FALSE.
            !
        CASE ('workrate')
            !
            PRINT_OPTIONS%PRINT_WORKRATE = .FALSE.
            !
        CASE ('workrate-pl')
            !
            PRINT_OPTIONS%PRINT_WORKRATEP = .FALSE.
            !
        CASE DEFAULT
            !
            S = 1
            !
        ! END CASE
        !
    END SELECT
    !
    END SUBROUTINE EXEC_SUPPRESS
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_READ_BCS_FROM_FILE(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    ! The external file must be named 'simulation.bcs'.
    BCS_OPTIONS%READ_BCS_FROM_FILE = .TRUE.
    !
    S = 0
    !
    END SUBROUTINE EXEC_READ_BCS_FROM_FILE
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_BOUNDARY_CONDITIONS(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    !
    SELECT CASE (TRIM(ADJUSTL(A)))
        !
        CASE ('uniaxial_grip','Uniaxial_Grip','UNIAXIAL_GRIP')
            !
            BCS_OPTIONS%BOUNDARY_CONDITIONS = 1
            !
        CASE ('uniaxial_symmetry','Uniaxial_Symmetry','UNIAXIAL_SYMMETRY')
            !
            BCS_OPTIONS%BOUNDARY_CONDITIONS = 2
            !
        CASE ('triaxial','Triaxial','TRIAXIAL')
            !
            BCS_OPTIONS%BOUNDARY_CONDITIONS = 3
            !
        CASE ('uniaxial_minimal','Uniaxial_Minimal','UNIAXIAL_MINIMAL')
            ! Uniaxial grip, but with minimal constraints. Only constrains
            ! normal directions on control faces and two nodes to prevent RBM.
            !
            BCS_OPTIONS%BOUNDARY_CONDITIONS = 4
            !
        CASE DEFAULT
            !
            S = 1
            !
        ! END CASE
        !
    END SELECT
    !
    END SUBROUTINE EXEC_BOUNDARY_CONDITIONS
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_LOADING_FACE(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    !
    SELECT CASE (TRIM(ADJUSTL(A)))
        !
        CASE ('x_min','X_min','X_Min','X_MIN','x0','X0','1')
            !
            BCS_OPTIONS%LOADING_FACE = 1
            !
        CASE ('x_max','X_max','X_Max','X_MAX','x1','X1','2')
            !
            BCS_OPTIONS%LOADING_FACE = 2
            !
        CASE ('y_min','Y_min','Y_Min','Y_MIN','y0','Y0','3')
            !
            BCS_OPTIONS%LOADING_FACE = 3
            !
        CASE ('y_max','Y_max','Y_Max','Y_MAX','y1','Y1','4')
            !
            BCS_OPTIONS%LOADING_FACE = 4
            !
        CASE ('z_min','Z_min','Z_Min','Z_MIN','z0','Z0','5')
            !
            BCS_OPTIONS%LOADING_FACE = 5
            !
        CASE ('z_max','Z_max','Z_Max','Z_MAX','z1','Z1','6')
            !
            BCS_OPTIONS%LOADING_FACE = 6
            !
        CASE DEFAULT
            !
            S = 1
            !
        ! END CASE
        !
    END SELECT
    !
    END SUBROUTINE EXEC_LOADING_FACE
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_LOADING_DIRECTION(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    !
    SELECT CASE (TRIM(ADJUSTL(A)))
        !
        CASE ('X','x','+x','+X','0')
            !
            BCS_OPTIONS%LOADING_DIRECTION = 0
            !
        CASE ('Y','y','+y','+Y','1')
            !
            BCS_OPTIONS%LOADING_DIRECTION = 1
            !
        CASE ('Z','z','+z','+Z','2')
            !
            BCS_OPTIONS%LOADING_DIRECTION = 2
            !
        CASE DEFAULT
            !
            S = 1
            !
        ! END CASE
        !
    END SELECT
    !
    END SUBROUTINE EXEC_LOADING_DIRECTION
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_STRAIN_RATE(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) BCS_OPTIONS%STRAIN_RATE
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_STRAIN_RATE
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_LOAD_RATE(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) BCS_OPTIONS%LOAD_RATE
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_LOAD_RATE
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NUMBER_OF_PHASES(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%NUMBER_OF_PHASES
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        RETURN
        !
    END IF
    !
    ! Allocate all material data arrays
    ALLOCATE(CRYS_OPTIONS%CRYSTAL_TYPE(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%M(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%GAMMADOT_0(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%H_0(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%G_0(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%G_S0(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%M_PRIME(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%GAMMADOT_S0(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%N(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%C11(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%C12(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%C13(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%C44(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%C66(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ! HCP specific arrays
    ALLOCATE(CRYS_OPTIONS%C_OVER_A(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%PRISMATIC_TO_BASAL(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%PYRAMIDAL_TO_BASAL(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ! BCT specific arrays
    ALLOCATE(CRYS_OPTIONS%HRATIO_BCT_A(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%HRATIO_BCT_B(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%HRATIO_BCT_C(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%HRATIO_BCT_D(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%HRATIO_BCT_E(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%HRATIO_BCT_F(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%HRATIO_BCT_G(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%HRATIO_BCT_H(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%HRATIO_BCT_I(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ! Hardening specific arrays
    ALLOCATE(CRYS_OPTIONS%CYCLIC_A(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%CYCLIC_C(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%LATENT_PARAMETERS(1:CRYS_OPTIONS%NUMBER_OF_PHASES, &
                & 1:11))
    ! Precipitation based hardening arrays
    ALLOCATE(CRYS_OPTIONS%A_P(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%F_P(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%B_P(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%R_P(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    !
    ! Rate dependence options for HCP/BCT materials
    ALLOCATE(CRYS_OPTIONS%USE_ANISO_M(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%ANISO_M(1:CRYS_OPTIONS%NUMBER_OF_PHASES, 1:10))
    !
    CRYS_OPTIONS%CRYSTAL_TYPE = -1
    CRYS_OPTIONS%M = -1.0D0
    CRYS_OPTIONS%GAMMADOT_0 = -1.0D0
    CRYS_OPTIONS%H_0 = -1.0D0
    CRYS_OPTIONS%G_0 = -1.0D0
    CRYS_OPTIONS%G_S0 = -1.0D0
    CRYS_OPTIONS%M_PRIME = -1.0D0
    CRYS_OPTIONS%GAMMADOT_S0 = -1.0D0
    CRYS_OPTIONS%N = -1.0D0
    CRYS_OPTIONS%C11 = -1.0D0
    CRYS_OPTIONS%C12 = -1.0D0
    CRYS_OPTIONS%C13 = -1.0D0
    CRYS_OPTIONS%C44 = -1.0D0
    CRYS_OPTIONS%C66 = -1.0D0
    !
    CRYS_OPTIONS%C_OVER_A = -1.0D0
    CRYS_OPTIONS%PRISMATIC_TO_BASAL = -1.0D0
    CRYS_OPTIONS%PYRAMIDAL_TO_BASAL = -1.0D0
    !
    CRYS_OPTIONS%HRATIO_BCT_A = -1.0D0
    CRYS_OPTIONS%HRATIO_BCT_B = -1.0D0
    CRYS_OPTIONS%HRATIO_BCT_C = -1.0D0
    CRYS_OPTIONS%HRATIO_BCT_D = -1.0D0
    CRYS_OPTIONS%HRATIO_BCT_E = -1.0D0
    CRYS_OPTIONS%HRATIO_BCT_F = -1.0D0
    CRYS_OPTIONS%HRATIO_BCT_G = -1.0D0
    CRYS_OPTIONS%HRATIO_BCT_H = -1.0D0
    CRYS_OPTIONS%HRATIO_BCT_I = -1.0D0
    !
    CRYS_OPTIONS%CYCLIC_A = -1.0D0
    CRYS_OPTIONS%CYCLIC_C = -1.0D0
    CRYS_OPTIONS%LATENT_PARAMETERS = -1.0D0
    !
    CRYS_OPTIONS%A_P = -1.0D0
    CRYS_OPTIONS%F_P = -1.0D0
    CRYS_OPTIONS%B_P = -1.0D0
    CRYS_OPTIONS%R_P = -1.0D0
    !
    CRYS_OPTIONS%ANISO_M = -1.0D0
    CRYS_OPTIONS%USE_ANISO_M = .FALSE.
    !
    S = 0
    !
    END SUBROUTINE EXEC_NUMBER_OF_PHASES
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_CRYSTAL_TYPE(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    S = 0
    !
    SELECT CASE (TRIM(ADJUSTL(A)))
        !
        CASE ('FCC', 'Fcc', 'fcc')
            !
            CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) = 1
            !
        CASE ('BCC', 'Bcc', 'bcc')
            !
            CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) = 2
            !
        CASE ('HCP', 'Hcp', 'hcp')
            !
            CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) = 3
            !
        CASE ('BCT', 'Bct', 'bct')
            !
            CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) = 4
            !
        CASE DEFAULT
            !
            S = 1
            !
        ! END CASE
        !
    END SELECT
    !
    END SUBROUTINE EXEC_CRYSTAL_TYPE
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_PHASE(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%PHASE
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_PHASE
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_M(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    INTEGER :: NSPACE, NVALS, II
    REAL(RK) :: TEMP_M, TEMP_M1, TEMP_M2, TEMP_M3, TEMP_M4, TEMP_M5, TEMP_M6, &
        & TEMP_M7, TEMP_M8, TEMP_M9, TEMP_M10
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! If current phase is FCC/BCC read in singular `m' value
    IF ((CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 1) .OR. &
        & (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 2)) THEN
        !
        ! Restructured input to be less ambiguous edge cases on input
        !
        ! Attempt to read the line as-is to determine the number of inputs
        READ(A, '(A)') LINE
        !
        ! Trim the full string into `NVALS' number of inputs
        II = 1
        NSPACE = COUNT( (/ (LINE(II:II), II=1, LEN_TRIM(LINE)) /) == " ")
        NVALS = NSPACE + 1
        !
        SELECT CASE(NVALS)
            !
            CASE(1)
                !
                ! If one value, read in the isotropic `m'
                !
                READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%M(CRYS_OPTIONS%PHASE)
                !
                IF (IOERR .NE. 0) THEN
                    !
                    S = 1
                    RETURN
                    !
                ELSE
                    !
                    S = 0
                    !
                END IF
                !
            CASE DEFAULT
                !
                ! The input is not 1 value so exit!
                !
                CALL PAR_QUIT&
                    &('Error  :     > One value expected for&
                    & "m" in FCC or BCC phases.')
                !
        END SELECT
        !
    ! Else if current phase is HCP, attempt a few options for `m'
    ELSE IF (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 3) THEN
        !
        ! Restructured input to be less ambiguous edge cases on input
        !
        ! Attempt to read the line as-is to determine the number of inputs
        READ(A, '(A)') LINE
        !
        ! Trim the full string into `NVALS' number of inputs
        II = 1
        NSPACE = COUNT( (/ (LINE(II:II), II=1, LEN_TRIM(LINE)) /) == " ")
        NVALS = NSPACE + 1
        !
        SELECT CASE(NVALS)
            !
            CASE(1)
                !
                ! If one value, read in the isotropic `m'
                !
                READ(A, *, IOSTAT=IOERR) TEMP_M
                !
                IF (IOERR .NE. 0) THEN
                    !
                    S = 1
                    RETURN
                    !
                ELSE
                    !
                    S = 0
                    !
                    ! Assign singular value to appropriate array
                    CRYS_OPTIONS%M(CRYS_OPTIONS%PHASE) = TEMP_M
                    !
                END IF
                !
            CASE(3)
                !
                ! If three values, read in the anisotropic `m'
                !
                READ(A, *, IOSTAT=IOERR) TEMP_M1, TEMP_M2, TEMP_M3
                !
                IF (IOERR .NE. 0) THEN
                    !
                    S = 1
                    RETURN
                    !
                ELSE
                    !
                    S = 0
                    !
                    ! Enable logical and assign temp values into array for phase
                    CRYS_OPTIONS%USE_ANISO_M(CRYS_OPTIONS%PHASE) = .TRUE.
                    CRYS_OPTIONS%ANISO_M(CRYS_OPTIONS%PHASE, 1) = TEMP_M1
                    CRYS_OPTIONS%ANISO_M(CRYS_OPTIONS%PHASE, 2) = TEMP_M2
                    CRYS_OPTIONS%ANISO_M(CRYS_OPTIONS%PHASE, 3) = TEMP_M3
                    !
                END IF
                !
            CASE DEFAULT
                !
                ! The input is not 1 or 3 values so exit!
                !
                CALL PAR_QUIT&
                    &('Error  :     > One or three values expected for&
                    & "m" in HCP phases.')
                !
        END SELECT
        !
    ! Else if current phase is BCT, attempt a few options for `m'
    ELSE IF (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 4) THEN
        !
        ! Restructured input to be less ambiguous edge cases on input
        !
        ! Attempt to read the line as-is to determine the number of inputs
        READ(A, '(A)') LINE
        !
        ! Trim the full string into `NVALS' number of inputs
        II = 1
        NSPACE = COUNT( (/ (LINE(II:II), II=1, LEN_TRIM(LINE)) /) == " ")
        NVALS = NSPACE + 1
        !
        SELECT CASE(NVALS)
            !
            CASE(1)
                !
                ! If one value, read in the isotropic `m'
                !
                READ(A, *, IOSTAT=IOERR) TEMP_M
                !
                IF (IOERR .NE. 0) THEN
                    !
                    S = 1
                    RETURN
                    !
                ELSE
                    !
                    S = 0
                    !
                    ! Assign singular value to appropriate array
                    CRYS_OPTIONS%M(CRYS_OPTIONS%PHASE) = TEMP_M
                    !
                END IF
                !
            CASE(10)
                !
                ! If ten values, read in the anisotropic `m'
                !
                READ(A, *, IOSTAT=IOERR) TEMP_M1, TEMP_M2, TEMP_M3, TEMP_M4, &
                    & TEMP_M5, TEMP_M6, TEMP_M7, TEMP_M8, TEMP_M9, TEMP_M10
                !
                IF (IOERR .NE. 0) THEN
                    !
                    S = 1
                    RETURN
                    !
                ELSE
                    !
                    S = 0
                    !
                    ! Enable logical and assign temp values into array for phase
                    CRYS_OPTIONS%USE_ANISO_M(CRYS_OPTIONS%PHASE) = .TRUE.
                    CRYS_OPTIONS%ANISO_M(CRYS_OPTIONS%PHASE, 1) = TEMP_M1
                    CRYS_OPTIONS%ANISO_M(CRYS_OPTIONS%PHASE, 2) = TEMP_M2
                    CRYS_OPTIONS%ANISO_M(CRYS_OPTIONS%PHASE, 3) = TEMP_M3
                    CRYS_OPTIONS%ANISO_M(CRYS_OPTIONS%PHASE, 4) = TEMP_M4
                    CRYS_OPTIONS%ANISO_M(CRYS_OPTIONS%PHASE, 5) = TEMP_M5
                    CRYS_OPTIONS%ANISO_M(CRYS_OPTIONS%PHASE, 6) = TEMP_M6
                    CRYS_OPTIONS%ANISO_M(CRYS_OPTIONS%PHASE, 7) = TEMP_M7
                    CRYS_OPTIONS%ANISO_M(CRYS_OPTIONS%PHASE, 8) = TEMP_M7
                    CRYS_OPTIONS%ANISO_M(CRYS_OPTIONS%PHASE, 9) = TEMP_M8
                    CRYS_OPTIONS%ANISO_M(CRYS_OPTIONS%PHASE, 10) = TEMP_M10
                    !
                END IF
                !
            CASE DEFAULT
                !
                ! The input is not 1 or 10 values so exit!
                !
                CALL PAR_QUIT&
                    &('Error  :     > One or ten values expected for&
                    & "m" in BCT phases.')
                !
        END SELECT
        !
    ! Exception error handling for future proofing
    ELSE
        !
        CALL PAR_QUIT&
            &('Error  :     > Invalid "m" value(s) input for material phase.')
        !
    END IF
    !
    END SUBROUTINE EXEC_M
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_GAMMADOT_0(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER,          INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%GAMMADOT_0(CRYS_OPTIONS%PHASE)
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_GAMMADOT_0
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_H_0(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%H_0(CRYS_OPTIONS%PHASE)
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_H_0
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_G_0(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    INTEGER :: NSPACE, NVALS, II
    REAL(RK) :: PRISM_TEMP, PYRAM_TEMP
    REAL(RK) :: G0BCT2, G0BCT3, G0BCT4, G0BCT5, G0BCT6, G0BCT7, G0BCT8, &
        & G0BCT9, G0BCT10
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! If current phase is FCC/BCC read in singular `g_0' value
    IF ((CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 1) .OR. &
        & (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 2)) THEN
        !
        ! Restructured input to be less ambiguous edge cases on input
        !
        ! Attempt to read the line as-is to determine the number of inputs
        READ(A, '(A)') LINE
        !
        ! Trim the full string into `NVALS' number of inputs
        II = 1
        NSPACE = COUNT( (/ (LINE(II:II), II=1, LEN_TRIM(LINE)) /) == " ")
        NVALS = NSPACE + 1
        !
        SELECT CASE(NVALS)
            !
            CASE(1)
                !
                ! If one value, read in the isotropic `g0'
                !
                READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE)
                !
                IF (IOERR .NE. 0) THEN
                    !
                    S = 1
                    RETURN
                    !
                ELSE
                    !
                    S = 0
                    !
                END IF
                !
            CASE DEFAULT
                !
                ! The input is not 1 value so exit!
                !
                CALL PAR_QUIT&
                    &('Error  :     > One value expected for&
                    & "g_0" in FCC or BCC phases.')
                !
        END SELECT
        !
    ! Else if current phase is HCP, attempt to read in three values for `g_0'
    ELSE IF (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 3) THEN
        !
        ! Restructured input to be less ambiguous edge cases on input
        !
        ! Attempt to read the line as-is to determine the number of inputs
        READ(A, '(A)') LINE
        !
        ! Trim the full string into `NVALS' number of inputs
        II = 1
        NSPACE = COUNT( (/ (LINE(II:II), II=1, LEN_TRIM(LINE)) /) == " ")
        NVALS = NSPACE + 1
        !
        SELECT CASE(NVALS)
            !
            CASE(3)
                !
                ! If three values, read in the anisotropic `g_0'
                !
                READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE), &
                    & PRISM_TEMP, PYRAM_TEMP
                !
                IF (IOERR .NE. 0) THEN
                    !
                    S = 1
                    ! Explicitly handle this error by the user since it acts on
                    ! a newly deprecated feature.
                    CALL PAR_QUIT('Error  :     > Three values for `g_0` &
                        &expected for HCP materials.')
                    !
                ELSE
                    !
                    S = 0
                    !
                    ! Scale the temporary variables by read-in `g_0'
                    ! (basal strength) in order to retrieve the internal
                    ! slip family strength ratios
                    !
                    CRYS_OPTIONS%PRISMATIC_TO_BASAL(CRYS_OPTIONS%PHASE) = &
                        & PRISM_TEMP / CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE)
                    CRYS_OPTIONS%PYRAMIDAL_TO_BASAL(CRYS_OPTIONS%PHASE) = &
                        & PYRAM_TEMP / CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE)
                    !
                END IF
                !
            CASE DEFAULT
                !
                ! The input is not 3 values so exit!
                !
                CALL PAR_QUIT&
                    &('Error  :     > Three values expected for&
                    & "g_0" in HCP phases.')
                !
        END SELECT
    !
    ! Else if current phase is BCT, attempt to read in ten values for `g_0'
    ELSE IF (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 4) THEN
        !
        ! Restructured input to be less ambiguous edge cases on input
        !
        ! Attempt to read the line as-is to determine the number of inputs
        READ(A, '(A)') LINE
        !
        ! Trim the full string into `NVALS' number of inputs
        II = 1
        NSPACE = COUNT( (/ (LINE(II:II), II=1, LEN_TRIM(LINE)) /) == " ")
        NVALS = NSPACE + 1
        !
        SELECT CASE(NVALS)
            !
            CASE(10)
                !
                ! If three values, read in the anisotropic `g_0'
                !
                READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE), &
                    & G0BCT2, G0BCT3, G0BCT4, G0BCT5, G0BCT6, G0BCT7, G0BCT8, &
                    & G0BCT9, G0BCT10
                !
                IF (IOERR .NE. 0) THEN
                    !
                    S = 1
                    ! Explicitly handle this error by the user since it acts on
                    ! a newly deprecated feature.
                    CALL PAR_QUIT('Error  :     > Ten values for `g_0` &
                        &expected for BCT materials.')
                    !
                ELSE
                    !
                    S = 0
                    !
                    ! Scale the temporary variables by read-in `g_0'
                    ! (basal strength) in order to retrieve the internal
                    ! slip family strength ratios
                    !
                    CRYS_OPTIONS%HRATIO_BCT_A(CRYS_OPTIONS%PHASE) = &
                        & G0BCT2 / CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE)
                    CRYS_OPTIONS%HRATIO_BCT_B(CRYS_OPTIONS%PHASE) = &
                        & G0BCT3 / CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE)
                    CRYS_OPTIONS%HRATIO_BCT_C(CRYS_OPTIONS%PHASE) = &
                        & G0BCT4 / CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE)
                    CRYS_OPTIONS%HRATIO_BCT_D(CRYS_OPTIONS%PHASE) = &
                        & G0BCT5 / CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE)
                    CRYS_OPTIONS%HRATIO_BCT_E(CRYS_OPTIONS%PHASE) = &
                        & G0BCT6 / CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE)
                    CRYS_OPTIONS%HRATIO_BCT_F(CRYS_OPTIONS%PHASE) = &
                        & G0BCT7 / CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE)
                    CRYS_OPTIONS%HRATIO_BCT_G(CRYS_OPTIONS%PHASE) = &
                        & G0BCT8 / CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE)
                    CRYS_OPTIONS%HRATIO_BCT_H(CRYS_OPTIONS%PHASE) = &
                        & G0BCT9 / CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE)
                    CRYS_OPTIONS%HRATIO_BCT_I(CRYS_OPTIONS%PHASE) = &
                        & G0BCT10 / CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE)
                    !
                END IF
                !
            CASE DEFAULT
                !
                ! The input is not 3 values so exit!
                !
                CALL PAR_QUIT&
                    &('Error  :     > Ten values expected for&
                    & "g_0" in BCT phases.')
                !
        END SELECT
        !
    ! Exception error handling for future proofing
    ELSE
        !
        CALL PAR_QUIT&
            &('Error  :     > Invalid "g_0" value(s) input for material phase.')
        !
    END IF
    !
    END SUBROUTINE EXEC_G_0
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_G_S0(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%G_S0(CRYS_OPTIONS%PHASE)
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_G_S0
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_M_PRIME(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%M_PRIME(CRYS_OPTIONS%PHASE)
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_M_PRIME
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_GAMMADOT_S0(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%GAMMADOT_S0(CRYS_OPTIONS%PHASE)
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_GAMMADOT_S0
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_N(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%N(CRYS_OPTIONS%PHASE)
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_N
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_C11(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%C11(CRYS_OPTIONS%PHASE)
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_C11
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_C12(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%C12(CRYS_OPTIONS%PHASE)
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        RETURN
        !
    END IF
    !
    ! Handle dupulicate assignments here since C12 MUST be read in for FCC/BCC.
    IF (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 1) THEN
        !
        CRYS_OPTIONS%C13(CRYS_OPTIONS%PHASE) = &
                & CRYS_OPTIONS%C12(CRYS_OPTIONS%PHASE)
        !
    ELSE IF (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 2) THEN
        !
        CRYS_OPTIONS%C13(CRYS_OPTIONS%PHASE) = &
                & CRYS_OPTIONS%C12(CRYS_OPTIONS%PHASE)
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_C12
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_C13(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    ! Read in data and store to allocated array in proper phase location.
    ! Do not allow for non-HCP OR BCT phases to read c13 in!
    IF ((CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 3) .OR. &
        & (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 4)) THEN
        !
        READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%C13(CRYS_OPTIONS%PHASE)
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        !
    ELSE
        !
        CALL PAR_QUIT&
            &('Error  :     > C13 input for FCC or BCC crystals invalid.')
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_C13
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_C44(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    ! Read in data and store to allocated array in proper phase location
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%C44(CRYS_OPTIONS%PHASE)
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_C44
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_C66(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    ! Read in data and store to allocated array in proper phase location.
    ! Do not allow for non-BCT phase to read c66 in!
    IF (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 4) THEN
        !
        READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%C66(CRYS_OPTIONS%PHASE)
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        !
    ELSE
        !
        CALL PAR_QUIT&
            &('Error  :     > C66 input for FCC, BCC, or HCP crystals invalid.')
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_C66
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_C_OVER_A(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    ! Read in data and store to allocated array in proper phase location
    IF ((CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 3) .OR. &
        & (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 4)) THEN
        !
        READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%C_OVER_A(CRYS_OPTIONS%PHASE)
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        !
    ELSE
        !
        CALL PAR_QUIT&
            &('Error  :     > c_over_a input for FCC or BCC crystals invalid.')
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_C_OVER_A
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_CYCLIC_A(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    IF (OPTIONS%HARD_TYPE .EQ. 'cyclic_isotropic') THEN
        !
        READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%CYCLIC_A(CRYS_OPTIONS%PHASE)
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_CYCLIC_A
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_CYCLIC_C(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    IF (OPTIONS%HARD_TYPE .EQ. 'cyclic_isotropic') THEN
        !
        READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%CYCLIC_C(CRYS_OPTIONS%PHASE)
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_CYCLIC_C
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_LATENT_PARAMETERS(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    IF (OPTIONS%HARD_TYPE .EQ. 'latent') THEN
        !
        SELECT CASE (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE))
            !
            CASE (1) ! FCC
                !
                READ(A, *, IOSTAT=IOERR) &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 1), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 2), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 3), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 4), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 5)
                !
                IF (IOERR .NE. 0) THEN
                    !
                    S = 1
                    RETURN
                    !
                END IF
                !
            CASE (2) ! BCC
                !
                READ(A, *, IOSTAT=IOERR) &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 1), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 2), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 3), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 4), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 5), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 6), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 7)
                !
                IF (IOERR .NE. 0) THEN
                    !
                    S = 1
                    RETURN
                    !
                END IF
                !
            CASE (3) ! HCP
                !
                READ(A, *, IOSTAT=IOERR) &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 1), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 2), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 3), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 4), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 5), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 6), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 7), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 8)
                !
                IF (IOERR .NE. 0) THEN
                    !
                    S = 1
                    RETURN
                    !
                END IF
                !
            CASE (4) ! BCT
                !
                READ(A, *, IOSTAT=IOERR) &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 1), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 2), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 3), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 4), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 5), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 6), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 7), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 8), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 9), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 10), &
                    & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE, 11)
                !
                IF (IOERR .NE. 0) THEN
                    !
                    S = 1
                    RETURN
                    !
                END IF
                !
            CASE DEFAULT
                !
                CALL PAR_QUIT('Error  :     > Invalid latent hardening &
                    &parameters provided.')
                !
        END SELECT
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_LATENT_PARAMETERS
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_A_P(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%A_P(CRYS_OPTIONS%PHASE)
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_A_P
    !
    !===========================================================================
    !    !
    SUBROUTINE EXEC_F_P(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%F_P(CRYS_OPTIONS%PHASE)
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_F_P
    !
    !===========================================================================
    !    !
    SUBROUTINE EXEC_B_P(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%B_P(CRYS_OPTIONS%PHASE)
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_B_P
    !
    !===========================================================================
    !    !
    SUBROUTINE EXEC_R_P(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%R_P(CRYS_OPTIONS%PHASE)
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_R_P
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NUMBER_OF_STRAIN_STEPS(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_STEPS
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        RETURN
        !
    END IF
    !
    IF (UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_STEPS .LT. 1) THEN
        !
        CALL PAR_QUIT('Error  :     > Invalid number of strain steps defined.')
        !
    END IF
    !
    ! Allocate all target arrays
    ALLOCATE(UNIAXIAL_OPTIONS%TARGET_STRAIN(&
        & 1:UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_STEPS, 1:3))
    !
    UNIAXIAL_OPTIONS%TARGET_STRAIN = -1.0D0
    UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET = 1
    !
    S = 0
    !
    END SUBROUTINE EXEC_NUMBER_OF_STRAIN_STEPS
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NUMBER_OF_STRAIN_RATE_JUMPS(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        RETURN
        !
    END IF
    !
    IF (UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS .LT. 0) THEN
        !
        CALL PAR_QUIT('Error  :     > Invalid number of strain rate jumps &
            &defined.')
        S = 1
        !
    ELSE IF (UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS .EQ. 0) THEN
        !
        S = 0
        RETURN
        !
    ELSE
        !
        ! Allocate all target arrays
        ALLOCATE(UNIAXIAL_OPTIONS%STRAIN_RATE_JUMP&
            & (1:UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS, 1:2))
        !
        UNIAXIAL_OPTIONS%STRAIN_RATE_JUMP = -1.0D0
        UNIAXIAL_OPTIONS%CURRENT_STRAIN_JUMP = 1
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_NUMBER_OF_STRAIN_RATE_JUMPS
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_TARGET_STRAIN(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    IF ((UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_STEPS .EQ. 0) .OR. &
        & (UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET .GT. &
        & UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_STEPS)) THEN
        !
        CALL PAR_QUIT('Error  :     > Input "number_of_strain_steps" is &
            &invalid.')
        !
    ELSE
        !
        READ(A, *, IOSTAT=IOERR) &
            & UNIAXIAL_OPTIONS%TARGET_STRAIN(&
                & UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET, 1), &
            & UNIAXIAL_OPTIONS%TARGET_STRAIN(&
                & UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET, 2), &
            & UNIAXIAL_OPTIONS%TEMP
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        !
        IF (UNIAXIAL_OPTIONS%TEMP .EQ. 'print_data') THEN
            !
            UNIAXIAL_OPTIONS%TARGET_STRAIN(&
                &UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET, 3) = 0
            !
        ELSE IF (UNIAXIAL_OPTIONS%TEMP .EQ. 'suppress_data') THEN
            !
            UNIAXIAL_OPTIONS%TARGET_STRAIN(&
                & UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET, 3) = 1
            !
        ELSE
            !
            CALL PAR_QUIT('Fatal error: String for print control is invalid &
                &in *.config file.')
            !
        END IF
        !
        UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET = &
            & UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET + 1
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_TARGET_STRAIN
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_STRAIN_RATE_JUMP(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    IF ((UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS .GT. 0) .AND. &
        & (UNIAXIAL_OPTIONS%CURRENT_STRAIN_JUMP .LE. &
        & UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS)) THEN
        !
        READ(A, *, IOSTAT=IOERR) &
            & UNIAXIAL_OPTIONS%STRAIN_RATE_JUMP(&
                & UNIAXIAL_OPTIONS%CURRENT_STRAIN_JUMP, 1), &
            & UNIAXIAL_OPTIONS%STRAIN_RATE_JUMP(&
                & UNIAXIAL_OPTIONS%CURRENT_STRAIN_JUMP, 2)
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        !
        UNIAXIAL_OPTIONS%CURRENT_STRAIN_JUMP = &
            & UNIAXIAL_OPTIONS%CURRENT_STRAIN_JUMP + 1
        !
    ELSE IF (UNIAXIAL_OPTIONS%CURRENT_STRAIN_JUMP .GT. &
        & UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS) THEN
        !
        CALL PAR_QUIT('Error  :     > Invalid number of strain jumps defined.')
        S = 1
        !
    ELSE IF (UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS .EQ. 0) THEN
        !
        RETURN
        S = 0
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_STRAIN_RATE_JUMP
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NUMBER_OF_LOAD_STEPS(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) UNIAXIAL_OPTIONS%NUMBER_OF_LOAD_STEPS
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        RETURN
        !
    END IF
    !
    IF (UNIAXIAL_OPTIONS%NUMBER_OF_LOAD_STEPS .LT. 1) THEN
        !
        CALL PAR_QUIT('Error  :     > Invalid number of load steps defined.')
        !
    END IF
    !
    ! Allocate all target arrays
    ALLOCATE(UNIAXIAL_OPTIONS%TARGET_LOAD(&
        & 1:UNIAXIAL_OPTIONS%NUMBER_OF_LOAD_STEPS, 1:4))
    !
    UNIAXIAL_OPTIONS%TARGET_LOAD = -1.0D0
    UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET = 1
    !
    S = 0
    !
    END SUBROUTINE EXEC_NUMBER_OF_LOAD_STEPS
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_TARGET_LOAD(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    IF ((UNIAXIAL_OPTIONS%NUMBER_OF_LOAD_STEPS .EQ. 0) .OR. &
        & (UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET .GT. &
        & UNIAXIAL_OPTIONS%NUMBER_OF_LOAD_STEPS)) THEN
        !
        CALL PAR_QUIT('Error  :     > Input "number_of_load_steps" is invalid.')
        !
    ELSE
        !
        READ(A, *, IOSTAT=IOERR) &
            & UNIAXIAL_OPTIONS%TARGET_LOAD(&
                & UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET, 1), &
            & UNIAXIAL_OPTIONS%TARGET_LOAD(&
                & UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET, 2), &
            & UNIAXIAL_OPTIONS%TARGET_LOAD(&
                & UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET, 3), &
            & UNIAXIAL_OPTIONS%TEMP
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        !
        IF (UNIAXIAL_OPTIONS%TEMP .EQ. 'print_data') THEN
            !
            UNIAXIAL_OPTIONS%TARGET_LOAD(&
                & UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET, 4) = 0
            !
        ELSE IF (UNIAXIAL_OPTIONS%TEMP .EQ. 'suppress_data') THEN
            !
            UNIAXIAL_OPTIONS%TARGET_LOAD(&
                & UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET, 4) = 1
            !
        ELSE
            !
            CALL PAR_QUIT('Error  :     > Invalid printing options in &
                & simulation.config.')
            !
        END IF
        !
        UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET = &
            & UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET + 1
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_TARGET_LOAD
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_MAX_BC_ITER(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    INTEGER :: TEMP
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) TEMP
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        RETURN
        !
    END IF
    !
    TRIAXCSR_OPTIONS%MAX_BC_ITER = TEMP
    TRIAXCLR_OPTIONS%MAX_BC_ITER = TEMP
    !
    s = 0
    !
    END SUBROUTINE EXEC_MAX_BC_ITER
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_MIN_PERT_FRAC(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) TRIAXCSR_OPTIONS%MIN_PERT_FRAC
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_MIN_PERT_FRAC
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_LOAD_TOL_ABS(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    REAL(RK) :: TEMP
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) TEMP
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        RETURN
        !
    END IF
    !
    TRIAXCSR_OPTIONS%LOAD_TOL_ABS = TEMP
    TRIAXCLR_OPTIONS%LOAD_TOL_ABS = TEMP
    !
    S = 0
    !
    END SUBROUTINE EXEC_LOAD_TOL_ABS
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_LOAD_TOL_REL(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) TRIAXCSR_OPTIONS%LOAD_TOL_REL
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_LOAD_TOL_REL
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NUMBER_OF_CSR_LOAD_STEPS(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) TRIAXCSR_OPTIONS%NUMBER_OF_CSR_LOAD_STEPS
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        RETURN
        !
    END IF
    !
    IF (TRIAXCSR_OPTIONS%NUMBER_OF_CSR_LOAD_STEPS .LT. 1) THEN
        !
        CALL PAR_QUIT('Error  :     > Invalid number of load steps defined.')
        !
    END IF
    !
    ! Allocate all target arrays
    ALLOCATE(TRIAXCSR_OPTIONS%TARGET_CSR_LOAD&
        & (1:TRIAXCSR_OPTIONS%NUMBER_OF_CSR_LOAD_STEPS, 1:6))
    !
    TRIAXCSR_OPTIONS%TARGET_CSR_LOAD = -1.0D0
    TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET = 1
    !
    S = 0
    !
    END SUBROUTINE EXEC_NUMBER_OF_CSR_LOAD_STEPS
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_TARGET_CSR_LOAD(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    IF ((TRIAXCSR_OPTIONS%NUMBER_OF_CSR_LOAD_STEPS .EQ. 0) .OR. &
        & (TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET .GT. &
        & TRIAXCSR_OPTIONS%NUMBER_OF_CSR_LOAD_STEPS)) THEN
        !
        CALL PAR_QUIT&
            & ('Error  :     > Input "number_of_csr_load_steps" is invalid.')
        !
        !
    ELSE
        !
        READ(A, *, IOSTAT=IOERR) &
            & TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(&
                & TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET, 1), &
            & TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(&
                & TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET, 2), &
            & TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(&
                & TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET, 3), &
            & TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(&
                & TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET, 4), &
            & TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(&
                & TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET, 5), &
            & TRIAXCSR_OPTIONS%TEMP_CSR
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        !
        IF (TRIAXCSR_OPTIONS%TEMP_CSR .EQ. 'print_data') THEN
            !
            TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(&
                & TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET, 6) = 0
            !
        ELSE IF (TRIAXCSR_OPTIONS%TEMP_CSR .EQ. 'suppress_data') THEN
            !
            TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(&
                & TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET, 6) = 1
            !
        ELSE
            !
            CALL PAR_QUIT('Error  :     > Invalid printing options in &
                &simulation.config.')
            !
        END IF
        !
        TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET = &
            & TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET + 1
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_TARGET_CSR_LOAD
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_MAX_STRAIN_INCR(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) TRIAXCLR_OPTIONS%MAX_STRAIN_INCR
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_MAX_STRAIN_INCR
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_MAX_STRAIN(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    REAL(RK) :: TEMP
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) TEMP
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        RETURN
        !
    END IF
    !
    TRIAXCSR_OPTIONS%MAX_STRAIN = TEMP
    TRIAXCLR_OPTIONS%MAX_STRAIN = TEMP
    !
    S = 0
    !
    END SUBROUTINE EXEC_MAX_STRAIN
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_MAX_EQSTRAIN(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    REAL(RK) :: TEMP
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) TEMP
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        RETURN
        !
    END IF
    !
    TRIAXCSR_OPTIONS%MAX_EQSTRAIN = TEMP
    TRIAXCLR_OPTIONS%MAX_EQSTRAIN = TEMP
    !
    S = 0
    !
    END SUBROUTINE EXEC_MAX_EQSTRAIN
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NUMBER_OF_CLR_LOAD_STEPS(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) TRIAXCLR_OPTIONS%NUMBER_OF_CLR_LOAD_STEPS
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        RETURN
        !
    END IF
    !
    IF (TRIAXCLR_OPTIONS%NUMBER_OF_CLR_LOAD_STEPS .LT. 1) THEN
        !
        CALL PAR_QUIT('Error  :     > Invalid number of load steps defined.')
        !
    END IF
    !
    ! Allocate all target arrays
    ALLOCATE(TRIAXCLR_OPTIONS%TARGET_CLR_LOAD&
        & (1:TRIAXCLR_OPTIONS%NUMBER_OF_CLR_LOAD_STEPS, 1:5))
    ! Note that allocation for step that are dwell episodes is handled in
    ! driver_triaxclr_mod to avoid overhead here.
    !
    TRIAXCLR_OPTIONS%TARGET_CLR_LOAD = -1.0D0
    TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET = 1
    !
    S = 0
    !
    END SUBROUTINE EXEC_NUMBER_OF_CLR_LOAD_STEPS
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_TARGET_CLR_LOAD(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    IF ((TRIAXCLR_OPTIONS%NUMBER_OF_CLR_LOAD_STEPS .EQ. 0) .OR. &
        & (TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET .GT. &
        & TRIAXCLR_OPTIONS%NUMBER_OF_CLR_LOAD_STEPS)) THEN
        !
        CALL PAR_QUIT&
            & ('Error  :     > Input "number_of_clr_load_steps" is invalid.')
        !
        !
    ELSE
        !
        READ(A, *, IOSTAT=IOERR) &
            & TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(&
                & TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET, 1), &
            & TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(&
                & TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET, 2), &
            & TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(&
                & TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET, 3), &
            & TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(&
                & TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET, 4), &
            & TRIAXCLR_OPTIONS%TEMP_CLR
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        !
        IF (TRIAXCLR_OPTIONS%TEMP_CLR .EQ. 'print_data') THEN
            !
            TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(&
                & TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET, 5) = 0
            !
        ELSE IF (TRIAXCLR_OPTIONS%TEMP_CLR .EQ. 'suppress_data') THEN
            !
            TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(&
                & TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET, 5) &
            & = 1
            !
        ELSE
            !
            CALL PAR_QUIT('Error  :     > Invalid printing options in &
                &simulation.config')
            !
        END IF
        !
        TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET = &
            & TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET + 1
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_TARGET_CLR_LOAD
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NUMBER_OF_LOAD_RATE_JUMPS(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        RETURN
        !
    END IF
    !
    IF (TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS .LT. 0) THEN
        !
        CALL PAR_QUIT('Error  :     > Invalid number of load rate jumps &
            &defined.')
        S = 1
        !
    ELSE IF (TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS .EQ. 0) THEN
        !
        S = 0
        RETURN
        !
    ELSE
        !
        ! Allocate all target arrays
        ALLOCATE(TRIAXCLR_OPTIONS%LOAD_RATE_JUMP&
            & (1:TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS, 1:2))
        !
        TRIAXCLR_OPTIONS%LOAD_RATE_JUMP = -1.0D0
        TRIAXCLR_OPTIONS%CURRENT_LOAD_RATE_JUMP = 1
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_NUMBER_OF_LOAD_RATE_JUMPS
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_LOAD_RATE_JUMP(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    IF ((TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS .GT. 0) .AND. &
        & (TRIAXCLR_OPTIONS%CURRENT_LOAD_RATE_JUMP .LE. &
        & TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS)) THEN
        !
        READ(A, *, IOSTAT=IOERR) &
            & TRIAXCLR_OPTIONS%LOAD_RATE_JUMP(&
                & TRIAXCLR_OPTIONS%CURRENT_LOAD_RATE_JUMP, 1), &
            & TRIAXCLR_OPTIONS%LOAD_RATE_JUMP(&
                & TRIAXCLR_OPTIONS%CURRENT_LOAD_RATE_JUMP, 2)
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        !
        TRIAXCLR_OPTIONS%CURRENT_LOAD_RATE_JUMP = &
            & TRIAXCLR_OPTIONS%CURRENT_LOAD_RATE_JUMP + 1
        !
    ELSE IF (TRIAXCLR_OPTIONS%CURRENT_LOAD_RATE_JUMP .GT. &
        & TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS) THEN
        !
        CALL PAR_QUIT('Error  :     > Invalid number of load rate jumps &
            &defined.')
        S = 1
        !
    ELSE IF (TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS .EQ. 0) THEN
        !
        RETURN
        S = 0
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_LOAD_RATE_JUMP
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_NUMBER_OF_DWELL_EPISODES(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        RETURN
        !
    END IF
    !
    IF (TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES .LT. 0) THEN
        !
        CALL PAR_QUIT('Error  :     > Invalid number of dwell episodes &
            &defined.')
        S = 1
        !
    ELSE IF (TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES .EQ. 0) THEN
        !
        S = 0
        RETURN
        !
    ELSE
        !
        ! Allocate all target arrays
        ALLOCATE(TRIAXCLR_OPTIONS%DWELL_EPISODE&
            & (1:TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES, 1:4))
        !
        TRIAXCLR_OPTIONS%DWELL_EPISODE = -1.0D0
        TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE = 1
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_NUMBER_OF_DWELL_EPISODES
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_DWELL_EPISODE(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    IF ((TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES .GT. 0) .AND. &
        & (TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE .LE. &
        & TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES)) THEN
        !
        READ(A, *, IOSTAT=IOERR) &
            & TRIAXCLR_OPTIONS%DWELL_EPISODE(&
                & TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE, 1), &
            & TRIAXCLR_OPTIONS%DWELL_EPISODE(&
                & TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE, 2), &
            & TRIAXCLR_OPTIONS%DWELL_EPISODE(&
                & TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE, 3), &
            & TRIAXCLR_OPTIONS%TEMP_CLR
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        !
        IF (TRIAXCLR_OPTIONS%TEMP_CLR .EQ. 'print_data') THEN
            !
            TRIAXCLR_OPTIONS%DWELL_EPISODE(&
                & TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE, 4) = 0
            !
        ELSE IF (TRIAXCLR_OPTIONS%TEMP_CLR .EQ. 'suppress_data') THEN
            !
            TRIAXCLR_OPTIONS%DWELL_EPISODE(&
                & TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE, 4) = 1
            !
        ELSE
            !
            CALL PAR_QUIT('Error  :     > Invalid printing options in &
                &simulation.config.')
            !
        END IF
        !
        TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE = &
            & TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE + 1
        !
    ELSE IF (TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE .GT. &
        & TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES) THEN
        !
        CALL PAR_QUIT('Error  :     > Invalid number of dwell episodes &
            &defined.')
        S = 1
        !
    ELSE IF (TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES .EQ. 0) THEN
        !
        RETURN
        S = 0
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_DWELL_EPISODE
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_MAX_ITER_HARD_LIMIT(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) OPTIONS%MAX_ITER_HARD_LIMIT
    !
    IF (IOERR .NE. 0) THEN
        !
        S = 1
        !
    ELSE
        !
        S = 0
        !
    END IF
    !
    END SUBROUTINE EXEC_MAX_ITER_HARD_LIMIT
    !
END MODULE READ_INPUT_MOD
