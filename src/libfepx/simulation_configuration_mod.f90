! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE SIMULATION_CONFIGURATION_MOD
!
! Module for reading simulation configuration parameters
!
! Contains subroutines:
! INITIALIZE_OPTIONS - Initializes options and input.
! PROCESS_INPUT - Processes options and input.
! Various subroutines to execute specific input.
!
USE LIBF95, ONLY: EXECCOMMAND, NEWUNITNUMBER
USE INTRINSICTYPESMODULE, ONLY: RK=>REAL_KIND
!
USE PARALLEL_MOD, ONLY: PAR_MESSAGE, PAR_QUIT
USE CONVERGENCEMODULE, ONLY: CONVERGENCEKEYWORDINPUT
!
IMPLICIT NONE
!
! Private
!
INTEGER, PARAMETER, PRIVATE :: MAX_PATH_LEN=1024
INTEGER, PARAMETER, PRIVATE :: MAX_FILE_LEN=256
INTEGER, PRIVATE :: IOERR
CHARACTER(LEN=MAX_PATH_LEN), PRIVATE :: PRIV_FNAME
INTEGER, PRIVATE :: PRIV_INPUNIT
!
! Public
!
TYPE OPTIONS_TYPE
    !
    INTEGER  :: MAX_INCR
    REAL(RK) :: MAX_TOTAL_TIME
    INTEGER  :: DEF_CONTROL_BY
    LOGICAL  :: CHECK_NECKING
    REAL(RK) :: LOAD_TOL
    REAL(RK) :: DTIME_FACTOR
    LOGICAL  :: RESTART
    CHARACTER(LEN=MAX_FILE_LEN) :: RSFIELD_BASE_IN
    CHARACTER(LEN=MAX_FILE_LEN) :: RSFIELD_BASE_OUT
    CHARACTER(LEN=MAX_FILE_LEN) :: RSCTRL_IN
    CHARACTER(LEN=MAX_FILE_LEN) :: RSCTRL_OUT
    CHARACTER(LEN=12) :: HARD_TYPE
    CHARACTER(LEN=12) :: CYCLIC
    LOGICAL  :: READ_ORI_FROM_FILE
    CHARACTER(LEN=MAX_FILE_LEN) :: ORI_FILE
    LOGICAL  :: READ_PHASE_FROM_FILE
    CHARACTER(LEN=MAX_FILE_LEN) :: PHASE_FILE
    !
END TYPE OPTIONS_TYPE
!
TYPE PRINT_OPTIONS_TYPE
    !
    LOGICAL :: PRINT_COORDINATES
    LOGICAL :: PRINT_VELOCITIES
    LOGICAL :: PRINT_ORIENTATIONS
    LOGICAL :: PRINT_CRSS
    LOGICAL :: PRINT_EQSTRESS
    LOGICAL :: PRINT_STRESS
    LOGICAL :: PRINT_STRAIN
    LOGICAL :: PRINT_GAMMADOT
    LOGICAL :: PRINT_DEFF
    LOGICAL :: PRINT_DPEFF
    LOGICAL :: PRINT_EQSTRAIN
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
    REAL(RK), ALLOCATABLE :: C13(:)
    REAL(RK), ALLOCATABLE :: C44(:)
    ! HCP specific parameter data.
    REAL(RK), ALLOCATABLE :: C_OVER_A(:)
    REAL(RK), ALLOCATABLE :: PRISMATIC_TO_BASAL(:)
    REAL(RK), ALLOCATABLE :: PYRAMIDAL_TO_BASAL(:)
    ! Hardening related parameters.
    REAL(RK), ALLOCATABLE :: CYCLIC_PARAMETER_A(:)
    REAL(RK), ALLOCATABLE :: CYCLIC_PARAMETER_C(:)
    REAL(RK), ALLOCATABLE :: LATENT_PARAMETERS(:,:)
    !
END TYPE CRYS_OPTIONS_TYPE
!
TYPE UNIAXIAL_CONTROL_TYPE
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
TYPE POWDER_DIFFRACTION_OPTIONS_TYPE
    !
    LOGICAL :: RUN_POWDER_DIFFRACTION
    LOGICAL :: READ_VOLUME
    LOGICAL :: READ_INTERIOR
    CHARACTER(LEN=MAX_FILE_LEN) :: POWDER_DIFFRACTION_FILE
    CHARACTER(LEN=MAX_FILE_LEN) :: VOLUME_FILE
    CHARACTER(LEN=MAX_FILE_LEN) :: INTERIOR_FILE
    !
END TYPE POWDER_DIFFRACTION_OPTIONS_TYPE
!
TYPE(OPTIONS_TYPE) :: OPTIONS
TYPE(PRINT_OPTIONS_TYPE) :: PRINT_OPTIONS
TYPE(BOUNDARY_CONDITIONS_TYPE) :: BCS_OPTIONS
TYPE(CRYS_OPTIONS_TYPE) :: CRYS_OPTIONS
TYPE(UNIAXIAL_CONTROL_TYPE) :: UNIAXIAL_OPTIONS
TYPE(TRIAXCSR_OPTIONS_TYPE) :: TRIAXCSR_OPTIONS
TYPE(TRIAXCLR_OPTIONS_TYPE) :: TRIAXCLR_OPTIONS
TYPE(POWDER_DIFFRACTION_OPTIONS_TYPE) :: POWDER_DIFFRACTION_OPTIONS
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
    SUBROUTINE INITIALIZE_OPTIONS()
    !
    ! Initialize inputs and assign default values.
    !
    !---------------------------------------------------------------------------
    !
    ! Standard options
    OPTIONS%MAX_INCR = 50000
    OPTIONS%MAX_TOTAL_TIME = 12000.0
    OPTIONS%DEF_CONTROL_BY = UNIAXIAL_STRAIN_TARGET
    OPTIONS%CHECK_NECKING = .FALSE.
    OPTIONS%LOAD_TOL = 0.0
    OPTIONS%DTIME_FACTOR = 1.001
    OPTIONS%HARD_TYPE = 'isotropic' ! set what type of hardening module used
    OPTIONS%CYCLIC = 'isotropic' !Turkmen's model initial set
    OPTIONS%RESTART = .FALSE.
    OPTIONS%RSFIELD_BASE_IN = 'pre.restart.field'
    OPTIONS%RSFIELD_BASE_OUT = 'post.restart.field'
    OPTIONS%RSCTRL_IN = 'pre.restart.control'
    OPTIONS%RSCTRL_OUT = 'post.restart.control'
    OPTIONS%READ_ORI_FROM_FILE = .FALSE.
    OPTIONS%ORI_FILE = ''
    OPTIONS%READ_PHASE_FROM_FILE = .FALSE.
    OPTIONS%PHASE_FILE = ''
    !
    ! Printing options
    PRINT_OPTIONS%PRINT_COORDINATES = .FALSE.
    PRINT_OPTIONS%PRINT_VELOCITIES = .FALSE.
    PRINT_OPTIONS%PRINT_ORIENTATIONS = .FALSE.
    PRINT_OPTIONS%PRINT_CRSS = .FALSE.
    PRINT_OPTIONS%PRINT_STRAIN = .FALSE.
    PRINT_OPTIONS%PRINT_STRESS = .FALSE.
    PRINT_OPTIONS%PRINT_GAMMADOT = .FALSE.
    PRINT_OPTIONS%PRINT_DEFF = .FALSE.
    PRINT_OPTIONS%PRINT_DPEFF = .FALSE.
    PRINT_OPTIONS%PRINT_EQSTRESS = .FALSE.
    PRINT_OPTIONS%PRINT_EQSTRAIN = .FALSE.
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
    TRIAXCSR_OPTIONS%MIN_PERT_FRAC = 1e-3
    TRIAXCSR_OPTIONS%LOAD_TOL_ABS = 0.1
    TRIAXCSR_OPTIONS%LOAD_TOL_REL = 0.001
    TRIAXCSR_OPTIONS%MAX_STRAIN = 0.2
    TRIAXCSR_OPTIONS%MAX_EQSTRAIN = 0.2
    !
    ! Triaxial CLR options
    TRIAXCLR_OPTIONS%NUMBER_OF_CLR_LOAD_STEPS = 0
    TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS = 0
    TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES = 0
    TRIAXCLR_OPTIONS%MAX_BC_ITER = 10
    TRIAXCSR_OPTIONS%LOAD_TOL_ABS = 0.1
    TRIAXCLR_OPTIONS%MAX_STRAIN_INCR = 0.001
    TRIAXCLR_OPTIONS%MAX_STRAIN = 0.2
    TRIAXCLR_OPTIONS%MAX_EQSTRAIN = 0.2
    !
    ! Powder diffraction options
    POWDER_DIFFRACTION_OPTIONS%RUN_POWDER_DIFFRACTION = .FALSE.
    POWDER_DIFFRACTION_OPTIONS%READ_VOLUME = .FALSE.
    POWDER_DIFFRACTION_OPTIONS%READ_INTERIOR = .FALSE.
    POWDER_DIFFRACTION_OPTIONS%POWDER_DIFFRACTION_FILE = ''
    POWDER_DIFFRACTION_OPTIONS%VOLUME_FILE = ''
    POWDER_DIFFRACTION_OPTIONS%INTERIOR_FILE = ''
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
    IF (EXECCOMMAND('rsfield_base_in', CMDLINE, &
                & EXEC_RSFIELD_BASE_IN, STATUS)) RETURN
    IF (EXECCOMMAND('rsfield_base_out', CMDLINE, &
                & EXEC_RSFIELD_BASE_OUT, STATUS)) RETURN
    IF (EXECCOMMAND('rsctrl_in', CMDLINE, &
                & EXEC_RSCTRL_IN, STATUS)) RETURN
    IF (EXECCOMMAND('rsctrl_out', CMDLINE, &
                & EXEC_RSCTRL_OUT, STATUS)) RETURN
    IF (EXECCOMMAND('hard_type', CMDLINE, &
                & EXEC_HARD_TYPE, STATUS)) RETURN
    IF (EXECCOMMAND('read_ori_from_file', CMDLINE, &
                & EXEC_READ_ORI_FROM_FILE, STATUS)) RETURN
    IF (EXECCOMMAND('read_phase_from_file', CMDLINE, &
                & EXEC_READ_PHASE_FROM_FILE, STATUS)) RETURN
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
    IF (EXECCOMMAND('c_over_a', CMDLINE, &
                & EXEC_C_OVER_A, STATUS)) RETURN
    IF (EXECCOMMAND('prismatic_to_basal', CMDLINE, &
                & EXEC_PRISMATIC_TO_BASAL, STATUS)) RETURN
    IF (EXECCOMMAND('pyramidal_to_basal', CMDLINE, &
                & EXEC_PYRAMIDAL_TO_BASAL, STATUS)) RETURN
    IF (EXECCOMMAND('cyclic_parameter_a', CMDLINE, &
                & EXEC_CYCLIC_PARAMETER_A, STATUS)) RETURN
    IF (EXECCOMMAND('cyclic_parameter_c', CMDLINE, &
                & EXEC_CYCLIC_PARAMETER_C, STATUS)) RETURN
    IF (EXECCOMMAND('latent_parameters', CMDLINE, &
                & EXEC_LATENT_PARAMETERS, STATUS)) RETURN
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
    ! Powder diffraction options
    IF (EXECCOMMAND('run_powder_diffraction', CMDLINE, &
                & EXEC_RUN_POWDER_DIFFRACTION, STATUS)) RETURN
    IF (EXECCOMMAND('read_volume', CMDLINE, &
                & EXEC_READ_VOLUME, STATUS)) RETURN
    IF (EXECCOMMAND('read_interior', CMDLINE, &
                & EXEC_READ_INTERIOR, STATUS)) RETURN
    !
    IF (ConvergenceKeywordInput(CMDLINE, INUNIT, STATUS)) RETURN
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
    !
    !---------------------------------------------------------------------------
    !
    OPTIONS%RESTART  = .TRUE.
    S = 0
    !
    END SUBROUTINE exec_restart
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_RSFIELD_BASE_IN(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *) OPTIONS%RSFIELD_BASE_IN
    OPTIONS%RSFIELD_BASE_IN = ADJUSTL(OPTIONS%RSFIELD_BASE_IN)
    !
    S = 0
    !
    END SUBROUTINE EXEC_RSFIELD_BASE_IN
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_RSFIELD_BASE_OUT(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *) OPTIONS%RSFIELD_BASE_OUT
    OPTIONS%RSFIELD_BASE_OUT = ADJUSTL(OPTIONS%RSFIELD_BASE_OUT)
    !
    S = 0
    !
    END SUBROUTINE EXEC_RSFIELD_BASE_OUT
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_RSCTRL_IN(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *) OPTIONS%RSCTRL_IN
    OPTIONS%RSCTRL_IN = ADJUSTL(OPTIONS%RSCTRL_IN)
    !
    S = 0
    !
    END SUBROUTINE EXEC_RSCTRL_IN
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_RSCTRL_OUT(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *) OPTIONS%RSCTRL_OUT
    OPTIONS%RSCTRL_OUT = ADJUSTL(OPTIONS%RSCTRL_OUT)
    !
    S = 0
    !
    END SUBROUTINE EXEC_RSCTRL_OUT
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
        CASE('cyclic','Cyclic','CYCLIC')
            !
            OPTIONS%HARD_TYPE = 'cyclic'
            OPTIONS%CYCLIC = 'isotropic'
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
    INTEGER :: STATUS
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
        CASE ('vel')
            !
            PRINT_OPTIONS%PRINT_VELOCITIES = .TRUE.
            !
        CASE ('ori')
            !
            PRINT_OPTIONS%PRINT_ORIENTATIONS =.TRUE.
            !
        CASE ('crss')
            !
            PRINT_OPTIONS%PRINT_CRSS = .TRUE.
            !
        CASE ('strain-el')
            !
            PRINT_OPTIONS%PRINT_STRAIN = .TRUE.
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
        CASE ('stress-eq')
            !
            PRINT_OPTIONS%PRINT_EQSTRESS = .TRUE.
            !
        CASE ('strain-eq')
            !
            PRINT_OPTIONS%PRINT_EQSTRAIN = .TRUE.
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
    INTEGER :: STATUS
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
        CASE ('vel')
            !
            PRINT_OPTIONS%PRINT_VELOCITIES = .FALSE.
            !
        CASE ('ori')
            !
            PRINT_OPTIONS%PRINT_ORIENTATIONS =.FALSE.
            !
        CASE ('crss')
            !
            PRINT_OPTIONS%PRINT_CRSS = .FALSE.
            !
        CASE ('strain-el')
            !
            PRINT_OPTIONS%PRINT_STRAIN = .FALSE.
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
    BCS_OPTIONS%READ_BCS_FROM_FILE = .TRUE.
    READ(A, *, IOSTAT=IOERR) BCS_OPTIONS%BCS_FILE
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
        CASE ('minimal','Minimal','MINIMAL')
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
    ! HCP specific arrays
    ALLOCATE(CRYS_OPTIONS%C_OVER_A(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%PRISMATIC_TO_BASAL(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%PYRAMIDAL_TO_BASAL(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ! Hardening specific arrays
    ALLOCATE(CRYS_OPTIONS%CYCLIC_PARAMETER_A(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%CYCLIC_PARAMETER_C(1:CRYS_OPTIONS%NUMBER_OF_PHASES))
    ALLOCATE(CRYS_OPTIONS%LATENT_PARAMETERS(1:CRYS_OPTIONS%NUMBER_OF_PHASES, &
                & 1:8))
    !
    CRYS_OPTIONS%CRYSTAL_TYPE = -1
    CRYS_OPTIONS%M = -1.0_RK
    CRYS_OPTIONS%GAMMADOT_0 = -1.0_RK
    CRYS_OPTIONS%H_0 = -1.0_RK
    CRYS_OPTIONS%G_0 = -1.0_RK
    CRYS_OPTIONS%G_S0 = -1.0_RK
    CRYS_OPTIONS%M_PRIME = -1.0_RK
    CRYS_OPTIONS%GAMMADOT_S0 = -1.0_RK
    CRYS_OPTIONS%N = -1.0_RK
    CRYS_OPTIONS%C11 = -1.0_RK
    CRYS_OPTIONS%C12 = -1.0_RK
    CRYS_OPTIONS%C13 = -1.0_RK
    CRYS_OPTIONS%C44 = -1.0_RK
    !
    CRYS_OPTIONS%C_OVER_A = -1.0_RK
    CRYS_OPTIONS%PRISMATIC_TO_BASAL = -1.0_RK
    CRYS_OPTIONS%PYRAMIDAL_TO_BASAL = -1.0_RK
    !
    CRYS_OPTIONS%CYCLIC_PARAMETER_A = -1.0_RK
    CRYS_OPTIONS%CYCLIC_PARAMETER_C = -1.0_RK
    CRYS_OPTIONS%LATENT_PARAMETERS = -1.0_RK
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
        CASE ('FCC','Fcc','fcc')
            !            
            CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) = 1
            !        
        CASE ('BCC','Bcc','bcc')
            !            
            CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) = 2
            !
        CASE ('HCP','Hcp','hcp')
            !            
            CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) = 3
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
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%M(CRYS_OPTIONS%PHASE)
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
    !
    !---------------------------------------------------------------------------
    !
    READ(A, *, IOSTAT=IOERR) CRYS_OPTIONS%G_0(CRYS_OPTIONS%PHASE)
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
    ! Do not allow for non-HCP phase to read c13 in!
    IF (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 3) THEN
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
    SUBROUTINE EXEC_C_OVER_A(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    ! Read in data and store to allocated array in proper phase location
    IF (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 3) THEN
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
            &('Error  :     > c/a input for FCC or BCC crystals invalid.')
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_C_OVER_A
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_PRISMATIC_TO_BASAL(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    ! Read in data and store to allocated array in proper phase location
    IF (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 3) THEN
        !
        READ(A, *, IOSTAT=IOERR) &
            &CRYS_OPTIONS%PRISMATIC_TO_BASAL(CRYS_OPTIONS%PHASE)
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
        CALL PAR_QUIT('Error  :     > Hardening ratio input for FCC or BCC &
            &crystals invalid.')
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_PRISMATIC_TO_BASAL
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_PYRAMIDAL_TO_BASAL(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    ! Read in data and store to allocated array in proper phase location
    IF (CRYS_OPTIONS%CRYSTAL_TYPE(CRYS_OPTIONS%PHASE) .EQ. 3) THEN
        !
        READ(A, *, IOSTAT=IOERR) &
            &CRYS_OPTIONS%PYRAMIDAL_TO_BASAL(CRYS_OPTIONS%PHASE)
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
        CALL PAR_QUIT('Error  :     > Hardening ratio input for FCC or BCC &
            &crystals invalid.')
        !
    END IF
    !
    S = 0
    !
    END SUBROUTINE EXEC_PYRAMIDAL_TO_BASAL
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_CYCLIC_PARAMETER_A(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    IF (OPTIONS%HARD_TYPE .EQ. 'cyclic') THEN
        !
        READ(A, *, IOSTAT=IOERR) &
            &CRYS_OPTIONS%CYCLIC_PARAMETER_A(CRYS_OPTIONS%PHASE)
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
    END SUBROUTINE EXEC_CYCLIC_PARAMETER_A
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_CYCLIC_PARAMETER_C(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    IF (OPTIONS%HARD_TYPE .EQ. 'cyclic') THEN
        !
        READ(A, *, IOSTAT=IOERR) &
            &CRYS_OPTIONS%CYCLIC_PARAMETER_C(CRYS_OPTIONS%PHASE)
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
    END SUBROUTINE EXEC_CYCLIC_PARAMETER_C
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
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,1), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,2), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,3), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,4), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,5)
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
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,1), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,2), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,3), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,4), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,5), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,6), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,7)
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
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,1), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,2), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,3), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,4), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,5), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,6), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,7), &
                & CRYS_OPTIONS%LATENT_PARAMETERS(CRYS_OPTIONS%PHASE,8)
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
            CALL PAR_QUIT('Error  :     > Invalid latent hardening parameters &
                &provided.')
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
    ALLOCATE(UNIAXIAL_OPTIONS%TARGET_STRAIN(1:UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_STEPS, &
                & 1:3))
    !
    UNIAXIAL_OPTIONS%TARGET_STRAIN = -1.0_RK
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
        UNIAXIAL_OPTIONS%STRAIN_RATE_JUMP = -1.0_RK
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
        CALL PAR_QUIT('Error  :     > Input "number_of_strain_steps" is invalid.')
        !   
    ELSE
        !
        READ(A, *, IOSTAT=IOERR) &
            & UNIAXIAL_OPTIONS%TARGET_STRAIN(UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET,1), &
            & UNIAXIAL_OPTIONS%TARGET_STRAIN(UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET,2), &
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
            UNIAXIAL_OPTIONS%TARGET_STRAIN(UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET,3) &
            & = 0
            !
        ELSE IF (UNIAXIAL_OPTIONS%TEMP .EQ. 'suppress_data') THEN
            !            
            UNIAXIAL_OPTIONS%TARGET_STRAIN(UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET,3) &
            & = 1
            !
        ELSE
            !
            CALL PAR_QUIT('Fatal error: String for print control is invalid in *.config file.') 
            !
        END IF
        !
        UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET = UNIAXIAL_OPTIONS%CURRENT_STRAIN_TARGET + 1
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
            &UNIAXIAL_OPTIONS%STRAIN_RATE_JUMP(UNIAXIAL_OPTIONS%CURRENT_STRAIN_JUMP,1), &
            & UNIAXIAL_OPTIONS%STRAIN_RATE_JUMP(UNIAXIAL_OPTIONS%CURRENT_STRAIN_JUMP,2)
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        !
        UNIAXIAL_OPTIONS%CURRENT_STRAIN_JUMP = UNIAXIAL_OPTIONS%CURRENT_STRAIN_JUMP + 1
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
    ALLOCATE(UNIAXIAL_OPTIONS%TARGET_LOAD(1:UNIAXIAL_OPTIONS%NUMBER_OF_LOAD_STEPS, &
                & 1:4))
    !
    UNIAXIAL_OPTIONS%TARGET_LOAD = -1.0_RK
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
            & UNIAXIAL_OPTIONS%TARGET_LOAD(UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET,1), &
            & UNIAXIAL_OPTIONS%TARGET_LOAD(UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET,2), &
            & UNIAXIAL_OPTIONS%TARGET_LOAD(UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET,3), &
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
            UNIAXIAL_OPTIONS%TARGET_LOAD(UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET,4) &
            & = 0
            !
        ELSE IF (UNIAXIAL_OPTIONS%TEMP .EQ. 'suppress_data') THEN
            !            
            UNIAXIAL_OPTIONS%TARGET_LOAD(UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET,4) &
            & = 1
            !
        ELSE
            !
            CALL PAR_QUIT('Error  :     > Invalid printing options in simulation.cfg.')
            !
        END IF
        !
        UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET = UNIAXIAL_OPTIONS%CURRENT_LOAD_TARGET + 1
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
    TRIAXCSR_OPTIONS%TARGET_CSR_LOAD = -1.0_RK
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
            & TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET,1), &
            & TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET,2), &
            & TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET,3), &
            & TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET,4), &
            & TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET,5), &
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
            TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET,6) &
            & = 0
            !
        ELSE IF (TRIAXCSR_OPTIONS%TEMP_CSR .EQ. 'suppress_data') THEN
            !            
            TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET,6) &
            & = 1
            !
        ELSE
            !
            CALL PAR_QUIT('Error  :     > Invalid printing options in simulation.cfg.')
            !
        END IF
        !
        TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET = TRIAXCSR_OPTIONS%CURRENT_CSR_LOAD_TARGET + 1
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
    TRIAXCLR_OPTIONS%TARGET_CLR_LOAD = -1.0_RK
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
            & TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET,1), &
            & TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET,2), &
            & TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET,3), &
            & TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET,4), &
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
            TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET,5) &
            & = 0
            !
        ELSE IF (TRIAXCLR_OPTIONS%TEMP_CLR .EQ. 'suppress_data') THEN
            !            
            TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET,5) &
            & = 1
            !
        ELSE
            !
            CALL PAR_QUIT('Error  :     > Invalid printing options in simulation.cfg')
            !
        END IF
        !
        TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET = TRIAXCLR_OPTIONS%CURRENT_CLR_LOAD_TARGET + 1
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
        CALL PAR_QUIT('Error  :     > Invalid number of load rate jumps defined.')
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
        TRIAXCLR_OPTIONS%LOAD_RATE_JUMP = -1.0_RK
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
            & TRIAXCLR_OPTIONS%LOAD_RATE_JUMP(TRIAXCLR_OPTIONS%CURRENT_LOAD_RATE_JUMP,1), &
            & TRIAXCLR_OPTIONS%LOAD_RATE_JUMP(TRIAXCLR_OPTIONS%CURRENT_LOAD_RATE_JUMP,2)
        !
        IF (IOERR .NE. 0) THEN
            !
            S = 1
            RETURN
            !
        END IF
        ! 
        TRIAXCLR_OPTIONS%CURRENT_LOAD_RATE_JUMP = TRIAXCLR_OPTIONS%CURRENT_LOAD_RATE_JUMP + 1
        !
    ELSE IF (TRIAXCLR_OPTIONS%CURRENT_LOAD_RATE_JUMP .GT. &
        & TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS) THEN
        !
        CALL PAR_QUIT('Error  :     > Invalid number of load rate jumps defined.')
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
        CALL PAR_QUIT('Error  :     > Invalid number of dwell episodes defined.')
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
        TRIAXCLR_OPTIONS%DWELL_EPISODE = -1.0_RK
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
            & TRIAXCLR_OPTIONS%DWELL_EPISODE(TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE,1), &
            & TRIAXCLR_OPTIONS%DWELL_EPISODE(TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE,2), &
            & TRIAXCLR_OPTIONS%DWELL_EPISODE(TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE,3), &
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
            TRIAXCLR_OPTIONS%DWELL_EPISODE(TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE,4) &
            & = 0
            !
        ELSE IF (TRIAXCLR_OPTIONS%TEMP_CLR .EQ. 'suppress_data') THEN
            !            
            TRIAXCLR_OPTIONS%DWELL_EPISODE(TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE,4) &
            & = 1
            !
        ELSE
            !
            CALL PAR_QUIT('Error  :     > Invalid printing options in simulation.cfg.')
            !
        END IF
        !
        TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE = TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE + 1
        !
    ELSE IF (TRIAXCLR_OPTIONS%CURRENT_DWELL_EPISODE .GT. &
        & TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES) THEN
        !
        CALL PAR_QUIT('Error  :     > Invalid number of dwell episodes defined.')
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
    SUBROUTINE EXEC_RUN_POWDER_DIFFRACTION(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    !
    POWDER_DIFFRACTION_OPTIONS%RUN_POWDER_DIFFRACTION = .TRUE.
    READ(A, *, IOSTAT=IOERR) POWDER_DIFFRACTION_OPTIONS%POWDER_DIFFRACTION_FILE
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
    END SUBROUTINE EXEC_RUN_POWDER_DIFFRACTION
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_READ_VOLUME(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN) :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    POWDER_DIFFRACTION_OPTIONS%READ_VOLUME = .TRUE.
    READ(A, *, IOSTAT=IOERR) POWDER_DIFFRACTION_OPTIONS%VOLUME_FILE
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
    END SUBROUTINE EXEC_READ_VOLUME
    !
    !===========================================================================
    !
    SUBROUTINE EXEC_READ_INTERIOR(A, S)
    !
    CHARACTER(LEN=*), INTENT(IN)  :: A ! command argument string
    INTEGER, INTENT(OUT) :: S ! STATUS
    !
    !---------------------------------------------------------------------------
    !
    POWDER_DIFFRACTION_OPTIONS%READ_INTERIOR = .TRUE.
    READ(A, *, IOSTAT=IOERR) POWDER_DIFFRACTION_OPTIONS%INTERIOR_FILE
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
    END SUBROUTINE EXEC_READ_INTERIOR
    !
END MODULE SIMULATION_CONFIGURATION_MOD
