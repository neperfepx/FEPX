! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
PROGRAM FEPX
!
! Full-field 3D polycrystal FEM analysis with anisotropic elasticity and
! viscoplasticity. Primary arrays are automatically allocated here to be
! used in the subsequent driver routines which perform the simulation.
!
!-------------------------------------------------------------------------------
!
! LibF95 modules
USE INTRINSIC_TYPES_MOD, RK=>REAL_KIND
USE FILES_MOD
USE TIMER_MOD
!
! LibParallel modules
USE PARALLEL_MOD
USE GATHER_SCATTER_MOD
!
! Driver modules
USE DRIVER_ISO_VP_MOD
USE DRIVER_TRIAXCSR_MOD
USE DRIVER_TRIAXCLR_MOD
USE DRIVER_UNIAXIAL_CONTROL_MOD
!
! Other FEPX modules
USE READ_INPUT_MOD
USE BOUNDARY_CONDITIONS_MOD
USE MICROSTRUCTURE_MOD
USE FIBER_AVERAGE_MOD
USE WRITE_OUTPUT_MOD
USE DIMENSIONS_MOD
USE UNITS_MOD
USE SURFACE_MOD
USE QUADRATURE_MOD
!
IMPLICIT NONE
!
! Local variables:
!
CHARACTER(LEN = 80) :: IOFILE
INTEGER :: IOSTATUS
INTEGER :: AVG_NP_PER_PROC
INTEGER, DIMENSION(8) :: TIMEVALUES
!
! Parallel variables:
!
TYPE(TRACE) NP_TRACE, DOF_TRACE
!
! Spatial calculation:
!
LOGICAL,  ALLOCATABLE :: BCS(:)
REAL(RK), ALLOCATABLE :: VELOCITY(:)
REAL(RK), ALLOCATABLE :: PFORCE(:)
!
! Crystal arrays:
!
REAL(RK), ALLOCATABLE :: C0_ANGS(:, :, :, :)
REAL(RK), ALLOCATABLE :: RSTAR_N(:, :, :, :)
REAL(RK), ALLOCATABLE :: WTS(:, :)
REAL(RK), POINTER ::  KELAS(:, :), KEINV(:, :)
REAL(RK), ALLOCATABLE :: CRSS_N(:, :, :)
!
! Miscellaneous:
!
REAL(RK) :: CLOCK_START
INTEGER :: DEF_CONTROL_BY
INTEGER :: AUTO_TIME
INTEGER :: M_EL
INTEGER :: ALSTATUS, NUM_ELM_PART, NUM_NODE_PART
INTEGER, ALLOCATABLE :: PART_INFO(:,:)
INTEGER, ALLOCATABLE :: GLOBAL_INFO(:,:)
!
!-------------------------------------------------------------------------------
!
! Initialize parallel libraries.
!
CALL PAR_INIT()
!
! Initialize clock value
!
CALL CPU_TIME(CLOCK_START)
!
! Print FEPX header to console.
!
IF (MYID .EQ. 0) THEN
    !
    CALL DATE_AND_TIME(VALUES = TIMEVALUES)
    WRITE(DFLT_U,'(A)')'==========================    &
        & F   E   P   X   =========================='
    WRITE(DFLT_U,'(A)')'Info   : A finite element software package for &
        & polycrystal plasticity.'
    WRITE(DFLT_U,'(A)')'Info   : Version 1.2.1-3'
    WRITE(DFLT_U,'(A,I0,A)')'Info   : Running on ', NUMPROCS, ' cores.'
    WRITE(DFLT_U,'(A)')'Info   : <https://fepx.info>'
    WRITE(DFLT_U,'(A)')'Info   : Copyright (C) 1996-2021, DPLab, ACME Lab.'
    WRITE(DFLT_U,'(A)')'Info   : &
        &---------------------------------------------------------------'
    WRITE(DFLT_U,'(A,I0,A,I0,A,I0,A,I0,A,I2.2)')'Info   : Start time: ',&
        & TIMEVALUES(1), '-', TIMEVALUES(2), '-', TIMEVALUES(3), ' at ' ,&
        & TIMEVALUES(5), ':', TIMEVALUES(6)
    WRITE(DFLT_U,'(A)')'Info   : Loading simulation...'
    !
ENDIF
!
! Set default option values.
!
CALL INITIALIZE_OPTIONS()
!
IF (MYID .EQ. 0) THEN
    !
    WRITE(DFLT_U, '(A)') "Info   :   [i] Parsing file `simulation.config'..."
    !
END IF
!
IOFILE = 'simulation.config'
!
! Open options file.
!
IF (IOFILE .NE. '0') THEN
    !
    OPEN(UNIT = IUNITS(TMPI1_U), FILE = IOFILE, STATUS = 'old', &
        & ACTION = 'READ', IOSTAT = IOSTATUS)
    !
    IF (IOSTATUS .NE. 0) THEN
        !
        CALL PAR_QUIT("Error  :     > Failure to open `simulation.config'.")
        !
    END IF
    !
    ! Read optional input
    !
    CALL COMMANDLOOP(IUNITS(TMPI1_U), PROCESS_INPUT, IOSTATUS)
    !
    SELECT CASE(IOSTATUS)
        !
        CASE (LOOPSTAT_EOF)
            !
            CALL PAR_QUIT("Error  :     > Failure to open `simulation.config'.")
            !
        CASE (LOOPSTAT_NZ)
            !
            CALL PAR_QUIT('Error  :     > Non-zero status returned while &
                &reading options file.')
            !
        CASE (LOOPSTAT_NOMATCH)
            !
            CALL PAR_QUIT('Error  :     > No match on keyword while reading &
                &options file')
            !
    END SELECT
    !
ENDIF
!
! Read mesh size parameters, keeping file open to read rest of mesh later.
!
IOFILE = 'simulation.msh'
!
OPEN(UNIT = IUNITS(MSH_U), FILE = IOFILE, STATUS = 'old', ACTION = 'READ', &
    & IOSTAT = IOSTATUS)
!
IF (IOSTATUS .NE. 0) THEN
    !
    CALL PAR_QUIT("Error  :     > Failure to open `simulation.msh'.")
    !
END IF
!
! Extract the number of elements and nodes from the mesh and rewind.
!
CALL GET_MSH_SIZE(IUNITS(MSH_U))
!
! Set sizes for automatic allocation.
!
MAXEL = NUMELM
MAXEL1 = MAXEL - 1
MAXNP = NUMNP
MAXNP1 = MAXNP - 1
MAXDOF = 3 * MAXNP
MAXDOF1 = MAXDOF - 1
!
! Partition arrays among processes.
!
! Note: the dof array was partitioned to keep the degrees of freedom on the same
!   processor as the corresponding NODES. This shouldn't be necessary, but there
!   was some kind of bug when I did it the other way. (deb)
!
CALL PAR_PARTITION(MAXEL, NUMPROCS, MYID, EL_SUB1, EL_SUP1)
CALL PAR_PARTITION(MAXNP, NUMPROCS, MYID, NP_SUB1, NP_SUP1)
DOF_SUB1 = 3 * NP_SUB1
DOF_SUP1 = 3 * NP_SUP1 + 2
!
! Check if the number of nodes per processor is small.
!
AVG_NP_PER_PROC = MAXNP / NUMPROCS
!
IF ((MYID .EQ. 0) .AND. (AVG_NP_PER_PROC .LE. 50)) THEN
    !
    WRITE(DFLT_U, '(A)') 'Warning:     > The average number of nodes per &
        &processor is small. MPI '
    WRITE(DFLT_U, '(A)') '               issues may cause this simulation &
        &to hang.'
    !
END IF
!
! Allocate DOF arrays and various crystal arrays.
!
ALLOCATE(BCS(DOF_SUB1:DOF_SUP1))
ALLOCATE(VELOCITY(DOF_SUB1:DOF_SUP1))
ALLOCATE(PFORCE(DOF_SUB1:DOF_SUP1))
ALLOCATE(C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1))
ALLOCATE(RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1))
ALLOCATE(WTS(0:NGRAIN1, EL_SUB1:EL_SUP1))
!
! Hardwire some parameters and retrieve the selected simulation type.
!
M_EL = EL_SUP1 - EL_SUB1 + 1
AUTO_TIME = 0   ! (0)=No, (1):Yes
DEF_CONTROL_BY = OPTIONS%DEF_CONTROL_BY
!
! Process crystal parameters from configuration file.
!
CALL PROCESS_MATERIAL_PARAMETERS(KELAS, KEINV)
!
! Write boundary condition information to console
!
IF (MYID .EQ. 0) THEN
    !
    WRITE(DFLT_U, '(A)') "Info   :   - Boundary conditions:"
    !
    ! Write boundary conditions type
    !
    IF (BCS_OPTIONS%BOUNDARY_CONDITIONS .EQ. 1) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :     > Uniaxial grip"
        !
    ELSE IF (BCS_OPTIONS%BOUNDARY_CONDITIONS .EQ. 2) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :     > Uniaxial symmetry"
        !
    ELSE IF (BCS_OPTIONS%BOUNDARY_CONDITIONS .EQ. 3) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :     > Triaxial"
        !
    ELSE IF (BCS_OPTIONS%BOUNDARY_CONDITIONS .EQ. 4) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :     > Uniaxial minimal"
        !
    END IF
    !
    ! Write control type
    !
    IF (OPTIONS%DEF_CONTROL_BY .EQ. UNIAXIAL_LOAD_TARGET) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :     > Load targeting, constant &
            &strain rate"
        !
    ELSE IF (OPTIONS%DEF_CONTROL_BY .EQ. UNIAXIAL_STRAIN_TARGET) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :     > Strain targeting, constant &
            &strain rate"
        !
    ELSE IF (OPTIONS%DEF_CONTROL_BY .EQ. TRIAXIAL_CONSTANT_STRAIN_RATE) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :     > Load targeting, constant &
            &strain rate"
        !
    ELSE IF (OPTIONS%DEF_CONTROL_BY .EQ. TRIAXIAL_CONSTANT_LOAD_RATE) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :     > Load targeting, constant load &
            &rate"
        !
    END IF
    !
    ! Write loading rates rates
    !
    IF ((OPTIONS%DEF_CONTROL_BY .EQ. UNIAXIAL_LOAD_TARGET) .OR. &
        & (OPTIONS%DEF_CONTROL_BY .EQ. UNIAXIAL_STRAIN_TARGET) .OR. &
        & (OPTIONS%DEF_CONTROL_BY .EQ. TRIAXIAL_CONSTANT_STRAIN_RATE)) THEN
        !
        WRITE(DFLT_U, '(A,E14.6)') "Info   :     > Strain rate: ", &
            & BCS_OPTIONS%STRAIN_RATE
        !
    ELSE IF (OPTIONS%DEF_CONTROL_BY .EQ. TRIAXIAL_CONSTANT_LOAD_RATE) THEN
        !
        WRITE(DFLT_U, '(A,E14.6)') "Info   :     > Load rate: ", &
            & BCS_OPTIONS%LOAD_RATE
        !
    END IF
    !
    ! Write loading directions
    !
    IF (BCS_OPTIONS%LOADING_DIRECTION .EQ. 0) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :     > Loading direction: x"
        !
    ELSE IF (BCS_OPTIONS%LOADING_DIRECTION .EQ. 1) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :     > Loading direction: y"
        !
    ELSE IF (BCS_OPTIONS%LOADING_DIRECTION .EQ. 2) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :     > Loading direction: z"
        !
    END IF
    !
    ! Write loading face (if uniaxial)
    !
    IF ((BCS_OPTIONS%BOUNDARY_CONDITIONS .EQ. 1) .OR. &
        & (BCS_OPTIONS%BOUNDARY_CONDITIONS .EQ. 2) .OR. &
        & (BCS_OPTIONS%BOUNDARY_CONDITIONS .EQ. 4)) THEN
        !
        IF (BCS_OPTIONS%LOADING_FACE .EQ. 1) THEN
            !
            WRITE(DFLT_U, '(A)') "Info   :     > Loading face: x0"
            !
        ELSE IF (BCS_OPTIONS%LOADING_FACE .EQ. 2) THEN
            !
            WRITE(DFLT_U, '(A)') "Info   :     > Loading face: x1"
            !
        ELSE IF (BCS_OPTIONS%LOADING_FACE .EQ. 3) THEN
            !
            WRITE(DFLT_U, '(A)') "Info   :     > Loading face: y0"
            !
        ELSE IF (BCS_OPTIONS%LOADING_FACE .EQ. 4) THEN
            !
            WRITE(DFLT_U, '(A)') "Info   :     > Loading face: y1"
            !
        ELSE IF (BCS_OPTIONS%LOADING_FACE .EQ. 5) THEN
            !
            WRITE(DFLT_U, '(A)') "Info   :     > Loading face: z0"
            !
        ELSE IF (BCS_OPTIONS%LOADING_FACE .EQ. 6) THEN
            !
            WRITE(DFLT_U, '(A)') "Info   :     > Loading face: z1"
            !
        END IF
        !
    END IF
    !
END IF
!
! Write deformation history stats to console
!
IF (MYID .EQ. 0) THEN
    !
    WRITE(DFLT_U, '(A)') "Info   :   - Deformation history:"
    !
    IF (OPTIONS%DEF_CONTROL_BY .EQ. UNIAXIAL_LOAD_TARGET) THEN
        !
        WRITE(DFLT_U, '(A, I0)') "Info   :     > Number of load steps: ", &
            & UNIAXIAL_OPTIONS%NUMBER_OF_LOAD_STEPS
        !
    ELSE IF (OPTIONS%DEF_CONTROL_BY .EQ. UNIAXIAL_STRAIN_TARGET) THEN
        !
        WRITE(DFLT_U, '(A, I0)') "Info   :     > Number of strain steps: ", &
            & UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_STEPS
        !
        IF (UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS .GT. 0) THEN
            !
            WRITE(DFLT_U, '(A, I0)') "Info   :     > Number of strain rate &
                &jumps: ", UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS
            !
        END IF
        !
    ELSE IF (OPTIONS%DEF_CONTROL_BY .EQ. TRIAXIAL_CONSTANT_STRAIN_RATE) THEN
        !
        WRITE(DFLT_U, '(A, I0)') "Info   :     > Number of load steps: ", &
            & TRIAXCSR_OPTIONS%NUMBER_OF_CSR_LOAD_STEPS
        !
    ELSE IF (OPTIONS%DEF_CONTROL_BY .EQ. TRIAXIAL_CONSTANT_LOAD_RATE) THEN
        !
        WRITE(DFLT_U, '(A, I0)') "Info   :     > Number of load steps: ", &
            & TRIAXCLR_OPTIONS%NUMBER_OF_CLR_LOAD_STEPS
        !
        IF (TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS .GT. 0) THEN
            !
            WRITE(DFLT_U, '(A, I0)') "Info   :     > Number of load rate &
                &jumps: ", TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS
            !
        END IF
        !
        IF (TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES .GT. 0) THEN
            !
            WRITE(DFLT_U, '(A, I0)') "Info   :     > Number of dwell &
                &episodes: ", TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES
            !
        END IF
        !
    END IF
    !
END IF
!
IF (MYID .EQ. 0) THEN
    !
    WRITE(DFLT_U, '(A)') "Info   :   [i] Parsed file `simulation.config'."
    !
END IF
!
! Allocate memory for spatial mesh and prepare to parse.
!
CALL ALLOCATE_MSH(ALSTATUS)
!
IF (ALSTATUS .NE. 0) THEN
    !
    CALL PAR_QUIT('Error  :     > Failed to allocate mesh.')
    !
ENDIF
!
IF (MYID .EQ. 0) THEN
    !
    WRITE(DFLT_U, '(A)') "Info   :   [i] Parsing file `simulation.msh'..."
    WRITE(DFLT_U, '(A)') 'Info   :   - Mesh parameters:'
    WRITE(DFLT_U, '(A,I0)') 'Info   :     > Node number: ', NUMNP
    WRITE(DFLT_U, '(A,I0)') 'Info   :     > Elt  number: ', NUMELM
    !
END IF
!
! Read in the spatial mesh via helper routines in READ_INPUT_MOD
!
CALL READ_SPATIAL_MSH(IUNITS(MSH_U), DFLT_U)
!
IF (MYID .EQ. 0) THEN
    !
    WRITE(DFLT_U, '(A)') "Info   :   [i] Parsed file `simulation.msh'."
    !
END IF
!
! Close the mesh unit - shows up as a comment in editors due to F77 formatting.
!
CLOSE(IUNITS(MSH_U))
!
! Open files for writing
!
CALL OPEN_OUTPUT_FILES(MYID)
!
! Check for external orientation or phase files
!
IF (OPTIONS%READ_ORI_FROM_FILE .EQV. .TRUE.) THEN
    !
    IOFILE = 'simulation.ori'
    !
    OPEN(UNIT = IUNITS(ORI_U), FILE = IOFILE, STATUS = 'OLD', ACTION = 'READ', &
        & IOSTAT = IOSTATUS)
    !
    IF (IOSTATUS .NE. 0) THEN
        !
        CALL PAR_QUIT("Error  :     > Failure to open `simulation.ori' file.")
        !
    END IF
    !
    IF (MYID .EQ. 0) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :   [i] Parsing file `simulation.ori'..."
        !
    END IF
    !
    ! Use the mesh parsing parent subroutine to parse the external .ori file.
    ELEMENT_ORIS = .FALSE.
    CALL READ_SPATIAL_MSH(IUNITS(ORI_U), DFLT_U)
    !
    IF (MYID .EQ. 0) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :   [i] Parsed file `simulation.ori'."
        !
    END IF
    !
    CLOSE(IUNITS(ORI_U))
    !
END IF
!
IF (OPTIONS%READ_PHASE_FROM_FILE .EQV. .TRUE.) THEN
    !
    IOFILE = 'simulation.phase'
    !
    OPEN(UNIT = IUNITS(PHASE_U), FILE = IOFILE, STATUS = 'old', &
        & ACTION = 'READ', IOSTAT = IOSTATUS)
    !
    IF (IOSTATUS .NE. 0) THEN
        !
        CALL PAR_QUIT("Error  :     > Failure to open `simulation.phase' file.")
        !
    END IF
    !
    IF (MYID .EQ. 0) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :   [i] Parsing file `simulation.phase'..."
        !
    END IF
    !
    CALL READ_GROUPS(IUNITS(PHASE_U), DFLT_U)
    !
    IF (MYID .EQ. 0) THEN
        !
        WRITE(DFLT_U, '(A)') "Info   :   [i] Parsed file `simulation.phase'."
        !
    END IF
    !
    CLOSE(IUNITS(PHASE_U))
    !
END IF
!
! Read-in or automatically calculated boundary conditions.
!
IF (BCS_OPTIONS%READ_BCS_FROM_FILE .EQV. .TRUE.) THEN
    !
    CALL READ_BCS(BCS, VELOCITY, PFORCE)
    !
ELSE IF (BCS_OPTIONS%READ_BCS_FROM_FILE .EQV. .FALSE.) THEN
    !
    CALL CALC_BCS(BCS, VELOCITY, COORDS, PFORCE)
    !
ELSE
    !
    CALL PAR_QUIT('Error  :     > Failure to initialize boundary conditions.')
    !
END IF
!
! Allocate slip arrays and assign element orientations as rotation matrices.
!
CALL ALLOCATE_CRSS_N(CRSS_N, NGRAIN1)
CALL ASSIGN_ANGLES_PHASES(C0_ANGS, CRSS_N, RSTAR_N, WTS)
CALL ALLOCATE_GAMMADOT(NGRAIN1)
CALL ALLOCATE_GACCUMSHEAR(NGRAIN1, NQPT1)
!
! Allocate memory to work related variables if requested or restart is enabled
!
IF ((PRINT_OPTIONS%PRINT_WORK) .OR. (PRINT_OPTIONS%PRINT_DEFRATE) &
    & .OR. (PRINT_OPTIONS%PRINT_RESTART)) THEN
    !
    ALLOCATE(D_TOT(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1))
    D_TOT = 0.0_RK
    !
END IF
!
IF ((PRINT_OPTIONS%PRINT_WORK) .OR. (PRINT_OPTIONS%PRINT_RESTART)) THEN
    !
    ALLOCATE(EL_WORK_N(EL_SUB1:EL_SUP1))
    ALLOCATE(EL_WORK_RATE_N(EL_SUB1:EL_SUP1))
    ALLOCATE(EL_WORK(EL_SUB1:EL_SUP1))
    !
    EL_WORK_N = 0.0_RK
    EL_WORK_RATE_N = 0.0_RK
    EL_WORK = 0.0_RK
    !
END IF
!
IF ((PRINT_OPTIONS%PRINT_WORKP) .OR. (PRINT_OPTIONS%PRINT_RESTART)) THEN
    !
    ALLOCATE(EL_WORKP_N(EL_SUB1:EL_SUP1))
    ALLOCATE(EL_WORKP_RATE_N(EL_SUB1:EL_SUP1))
    ALLOCATE(EL_WORKP(EL_SUB1:EL_SUP1))
    !
    EL_WORKP_N = 0.0_RK
    EL_WORKP_RATE_N = 0.0_RK
    EL_WORKP = 0.0_RK
    !
END IF
!
! Gather partition information for the post.report file.
!
ALLOCATE(PART_INFO(0:1, 0:NUMPROCS-1))
ALLOCATE(GLOBAL_INFO(0:1, 0:NUMPROCS-1))
PART_INFO = 0
GLOBAL_INFO = 0
!
NUM_ELM_PART  = (EL_SUP1 - EL_SUB1) + 1
NUM_NODE_PART = (NP_SUP1 - NP_SUB1) + 1
PART_INFO(0, MYID) = NUM_ELM_PART
PART_INFO(1, MYID) = NUM_NODE_PART
!
CALL PAR_GATHER(PART_INFO(0:1,MYID), GLOBAL_INFO, 2)
!
IF (MYID .EQ. 0) THEN
    !
    CALL WRITE_REPORT_FILE_HEADER(GLOBAL_INFO)
    !
    WRITE(DFLT_U, '(A)') 'Info   : Initializing simulation...'
    !
END IF
!
DEALLOCATE(PART_INFO)
DEALLOCATE(GLOBAL_INFO)
!
! Process load history from read-in options
!
SELECT CASE (DEF_CONTROL_BY)
    !
    CASE(UNIAXIAL_LOAD_TARGET, UNIAXIAL_STRAIN_TARGET)
        !
        CALL PROCESS_CTRL_DATA(VELOCITY)
        !
    CASE(TRIAXIAL_CONSTANT_STRAIN_RATE)
        !
        CALL PROCESS_CTRL_DATA_CSR(VELOCITY)
        !
    CASE(TRIAXIAL_CONSTANT_LOAD_RATE)
        !
        CALL PROCESS_CTRL_DATA_CLR
    !
END SELECT
!
! Allocate memory for elastic stress/strain state variables.
!
!CALL ALLOCATE_STRESS_STRAIN(ALSTATUS)
!
ALLOCATE(GELA_KK_BAR(0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1), &
    & GSIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1), &
    & STAT=ALSTATUS)
!
GELA_KK_BAR = 0.0_RK
GSIG_VEC_N  = 0.0_RK
!
IF (ALSTATUS .NE. 0) THEN
    !
    CALL PAR_QUIT('Error  :     > Failure to allocate stress/strain variables.')
    !
ENDIF
!
! Allocate memory for POST_UPDATE_N state variables.
!
!CALL ALLOCATE_POST_UPDATE_N(ALSTATUS)
!
ALLOCATE(PELA_KK_BAR(EL_SUB1:EL_SUP1), &
     & PSIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1), &
     & STAT=ALSTATUS)
!
PELA_KK_BAR = 0.0_RK
PSIG_VEC_N  = 0.0_RK
!
IF (ALSTATUS .NE. 0) THEN
    !
    CALL PAR_QUIT('Error  :     > Failure to allocate post_update_n variables.')
    !
ENDIF
!
! Initialize gauss quadrature points and weights.
!
CALL INITIALIZE()
CALL INIT_SURF()
!
! Prepare to use Mika's MPI gather/scatter routines.
!
IF (MYID .EQ. 0) THEN
    !
    WRITE(DFLT_U, '(A)') 'Info   :   - Initializing parallel gather/scatter &
        &routines...'
    !
END IF
!
CALL PART_SCATTER_SETUP(0, KDIM1, DOF_SUB1, DOF_SUP1, EL_SUB1, EL_SUP1, NODES, &
    & DOF_TRACE)
CALL PART_SCATTER_SETUP(0, NNPE, NP_SUB1, NP_SUP1, EL_SUB1, EL_SUP1, NP, &
    & NP_TRACE)
!
! Initialize fiber average routines.
!
IF (FIBER_AVERAGE_OPTIONS%RUN_FIBER_AVERAGE) THEN
    !
    IF (MYID .EQ. 0) THEN
        !
        WRITE(DFLT_U, '(A)') 'Info   :   - Initializing fiber averaging &
            &processing...'
        !
    END IF
    !
    CALL INITIALIZE_FIBER_AVERAGE(IUNITS(TMPI1_U))
    !
ENDIF
!
! Begin the deformation simulation proper.
!
! Note: the uniaxial deformation drivers obtain an initial velocity field
! guess from a viscoplastic isotropic solution. This is NOT utilized by
! the triaxial deformation drivers.
!
IF (((DEF_CONTROL_BY .EQ. UNIAXIAL_LOAD_TARGET) .OR. (DEF_CONTROL_BY .EQ. &
       & UNIAXIAL_STRAIN_TARGET)) .AND. (.NOT. OPTIONS%RESTART)) THEN
    !
    CALL DRIVER_VP_SOLVE(0, BCS, VELOCITY, PFORCE, DOF_TRACE, CRSS_N, C0_ANGS)
    !
ENDIF
!
! Solve the anisotropic elasto-viscoplastic problem with the selected driver.
!
SELECT CASE (DEF_CONTROL_BY)
    !
    CASE (UNIAXIAL_LOAD_TARGET, UNIAXIAL_STRAIN_TARGET)
        !
        CALL DRIVER_UNIAXIAL_CONTROL(BCS, VELOCITY, PFORCE, DOF_TRACE, &
            & C0_ANGS, CRSS_N, RSTAR_N, KEINV, WTS, AUTO_TIME, GAMMADOT, &
            & CLOCK_START)
        !
    CASE (TRIAXIAL_CONSTANT_STRAIN_RATE)
        !
        CALL DRIVER_TRIAX_CSR(BCS, VELOCITY, PFORCE, DOF_TRACE, C0_ANGS, &
            & CRSS_N, RSTAR_N, KEINV, WTS, AUTO_TIME, GAMMADOT, CLOCK_START)
        !
    CASE (TRIAXIAL_CONSTANT_LOAD_RATE)
        !
        CALL DRIVER_TRIAX_CLR( BCS, VELOCITY, PFORCE, DOF_TRACE, C0_ANGS, &
            & CRSS_N, RSTAR_N, KEINV, WTS, AUTO_TIME, GAMMADOT, CLOCK_START)
        !
    CASE DEFAULT
        !
        CALL PAR_QUIT('Error  :     > Invalid deformation control option.')
    !
END SELECT
!
CALL PAR_QUIT('Info   : Completed simulation successfully')
!
END PROGRAM FEPX
