! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
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
USE INTRINSICTYPESMODULE, RK=>REAL_KIND
USE FILESMODULE
!
! LibParallel modules
USE PARALLEL_MOD
USE GATHER_SCATTER
!
! Driver modules
USE DRIVERVPMODULE
USE DRIVER_TRIAXCSR_MOD
USE DRIVER_TRIAXCLR_MOD
USE DRIVER_UNIAXIAL_CONTROL_MOD
!
! Other FEPX modules
USE SIMULATION_CONFIGURATION_MOD
USE RESTARTMODULE
USE BOUNDARY_CONDITIONS_MOD
USE MICROSTRUCTURE_MOD
USE STRESS_STRAIN_MOD
USE POST_UPDATE_N_MOD
USE POWDER_DIFFRACTION_MOD
USE WRITE_OUTPUT_MOD
USE DIMSMODULE
USE UNITS_MOD
USE SURFACE_MOD
USE SURF_INFO_MOD
USE QUADRATURE_MOD
USE READ_INPUT_MOD
!
IMPLICIT NONE
!
! Local variables:
!
CHARACTER(LEN = 80) :: INFILE
INTEGER :: NARGS, IOSTATUS
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
INTEGER :: DEF_CONTROL_BY
INTEGER :: AUTO_TIME
INTEGER :: M_EL, I
INTEGER :: ALSTATUS, NUM_ELM_PART, NUM_NODE_PART
INTEGER, ALLOCATABLE :: PART_INFO(:,:)
!
!-------------------------------------------------------------------------------
!
! Initialize parallel libraries.
!
CALL PAR_INIT()
!
! Print FEPX header to console.
!
IF (MYID .EQ. 0) THEN
    !
    CALL DATE_AND_TIME(VALUES = TIMEVALUES)
    WRITE(DFLT_U,'(A)')'==========================    F   E   P   X   =========================='
    WRITE(DFLT_U,'(A)')'Info   : A finite element software package for polycrystal plasticity.'
    WRITE(DFLT_U,'(A)')'Info   : Version 1.0.0'
    WRITE(DFLT_U,'(A,I0,A)')'Info   : Running on ', NUMPROCS, ' processors.'
    WRITE(DFLT_U,'(A)')'Info   : <http://fepx.info>'
    WRITE(DFLT_U,'(A)')'Info   : Copyright (C) 1996-2020, DPLab, ACME Lab.'
    WRITE(DFLT_U,'(A)')'Info   : ---------------------------------------------------------------'
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
INFILE = 'simulation.config'
!
! Open options file.
!
IF (INFILE .NE. '0') THEN
    !
    OPEN(UNIT = IUNITS(TMPI1_U), FILE = INFILE, STATUS = 'old', &
        & ACTION = 'READ', IOSTAT = IOSTATUS)
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
INFILE = 'simulation.msh'
!
OPEN(UNIT = IUNITS(MSH_U), FILE = INFILE, STATUS = 'old', ACTION = 'READ', &
    & IOSTAT = IOSTATUS)
!
IF (IOSTATUS .NE. 0) THEN
    !
    CALL PAR_QUIT("Error  :     > Failure to open `simulation.msh' file.")
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
CALL READ_SPATIAL_MSH(IUNITS(MSH_U))
!
IF (MYID .EQ. 0) THEN
    !
    WRITE(DFLT_U, '(A)') "Info   :   [i] Parsed file `simulation.msh'."
    !
END IF
!
! Extract crystal parameters from configuration file.
!
CALL READ_MATERIAL_PARAMETERS(KELAS, KEINV)
!
IF (MYID .EQ. 0) THEN
    !
    WRITE(DFLT_U, '(A)') "Info   :   [i] Parsed file `simulation.config'."
    !
END IF
!
! Close the mesh unit - this shows up as a comment in text editors due to F77 formatting.
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
    INFILE = 'simulation.ori'
    !
    OPEN(UNIT = IUNITS(ORI_U), FILE = INFILE, STATUS = 'OLD', ACTION = 'READ', &
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
    CALL READ_SPATIAL_MSH(IUNITS(ORI_U))
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
    INFILE = 'simulation.phase'
    !
    OPEN(UNIT = IUNITS(PHASE_U), FILE = INFILE, STATUS = 'old', ACTION = 'READ', &
        & IOSTAT = IOSTATUS)
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
    CALL READ_GROUPS(IUNITS(PHASE_U))
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
! Allocate slip arrays and assign element orientations as rotation matrices.
!
CALL ALLOCATE_CRSS_N(CRSS_N, NGRAIN1)
CALL ASSIGN_ANGLES_PHASES(C0_ANGS, CRSS_N, RSTAR_N, WTS)
CALL ALLOCATE_GAMMADOT(NGRAIN1)
CALL ALLOCATE_GACCUMSHEAR(NGRAIN1, NQPT1)
!
! Gather partition information for the post.report file.
!
ALLOCATE(PART_INFO(0:1, 0:NUMPROCS-1))
PART_INFO = 0
!
NUM_ELM_PART  = (EL_SUP1 - EL_SUB1) + 1
NUM_NODE_PART = (NP_SUP1 - NP_SUB1) + 1
PART_INFO(0, MYID) = NUM_ELM_PART
PART_INFO(1, MYID) = NUM_NODE_PART
!
CALL PAR_GATHER(PART_INFO(0:1,MYID), PART_INFO, 2)
!
IF (MYID .EQ. 0) THEN
    !
    CALL WRITE_REPORT_FILE_HEADER(PART_INFO)
    !
    WRITE(DFLT_U, '(A)') 'Info   : Initializing simulation...'
    !
END IF
!
DEALLOCATE(PART_INFO)
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
! Read load history.
!
SELECT CASE (DEF_CONTROL_BY)
    !
    CASE(UNIAXIAL_LOAD_TARGET, UNIAXIAL_STRAIN_TARGET)
        !
        CALL READ_CTRL_DATA(VELOCITY)
        !
    CASE(TRIAXIAL_CONSTANT_STRAIN_RATE)
        !
        CALL READ_CTRL_DATA_CSR(VELOCITY)
        !
    CASE(TRIAXIAL_CONSTANT_LOAD_RATE)
        !
        CALL READ_CTRL_DATA_CLR
    !
END SELECT
!
! Allocate memory for elastic stress/strain state variables.
!
CALL ALLOCATE_STRESS_STRAIN(ALSTATUS)
!
IF (ALSTATUS .NE. 0) THEN
    !
    CALL PAR_QUIT('Error  :     > Failure to allocate stress/strain variables.')
    !
ENDIF
!
! Allocate memory for POST_UPDATE_N state variables.
!
CALL ALLOCATE_POST_UPDATE_N(ALSTATUS)
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
    WRITE(DFLT_U, '(A)') 'Info   :   - Initializing parallel gather/scatter routines...'
    !
END IF
!
CALL PART_SCATTER_SETUP(0, KDIM1, DOF_SUB1, DOF_SUP1, EL_SUB1, EL_SUP1, NODES, DOF_TRACE)
CALL PART_SCATTER_SETUP(0, NNPE, NP_SUB1, NP_SUP1, EL_SUB1, EL_SUP1, NP, NP_TRACE)
!
! Initialize powder diffraction routines.
!
IF (POWDER_DIFFRACTION_OPTIONS%RUN_POWDER_DIFFRACTION) THEN
    !
    IF (MYID .EQ. 0) THEN
        !
        WRITE(DFLT_U, '(A)') 'Info   :   - Initializing powder diffraction routines...'
        !
    END IF
    !
    CALL INITIALIZE_POWDER_DIFFRACTION(IUNITS(TMPI1_U), DOF_TRACE)
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
    CALL DRIVER_VP_SOLVE(0, BCS, VELOCITY, PFORCE, DOF_TRACE, NP_TRACE, &
        & CRSS_N, C0_ANGS, RSTAR_N, WTS)
    !
ENDIF
!
! Solve the anisotropic elastoviscoplastic problem with the selected driver.
!
SELECT CASE (DEF_CONTROL_BY)
    !
    CASE (UNIAXIAL_LOAD_TARGET, UNIAXIAL_STRAIN_TARGET)
        !
        CALL DRIVER_UNIAXIAL_CONTROL(BCS, VELOCITY, PFORCE, DOF_TRACE, &
            & NP_TRACE, C0_ANGS, CRSS_N, RSTAR_N, KEINV, WTS, AUTO_TIME, &
            & GAMMADOT)
        !
    CASE (TRIAXIAL_CONSTANT_STRAIN_RATE)
        !
        CALL DRIVER_TRIAX_CSR(BCS, VELOCITY, PFORCE, DOF_TRACE, NP_TRACE, &
            & C0_ANGS, CRSS_N, RSTAR_N, KEINV, WTS, AUTO_TIME, GAMMADOT)
        !
    CASE (TRIAXIAL_CONSTANT_LOAD_RATE)
        !
        CALL DRIVER_TRIAX_CLR( BCS, VELOCITY, PFORCE, DOF_TRACE, NP_TRACE, &
            & C0_ANGS, CRSS_N, RSTAR_N, KEINV, WTS, AUTO_TIME, GAMMADOT)
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
