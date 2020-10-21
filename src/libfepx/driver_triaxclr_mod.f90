! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE DRIVER_TRIAXCLR_MOD
!
! Driver for EVPS simulation using triaxial load control at constant load rate
!
! Contains subroutines:
! DRIVER_TRIAX_CLR: Driver for triaxial CLR simulation
! READ_CTRL_DATA_CLR: Read input data for load control at constant load rate
! READ_TRIAXCLR_RESTART: Read restart files for triaxial CLR simulations
! PRINT_HEADERS: Print headers to output files
!
USE INTRINSICTYPESMODULE, RK=>REAL_KIND
USE GATHER_SCATTER
USE PARALLEL_MOD
!
USE DIMSMODULE
USE FIBER_AVERAGE_MOD, ONLY: RUN_FIBER_AVERAGE
USE DRIVER_UTILITIES_MOD
USE TIMERMODULE
USE READ_INPUT_MOD
USE SURF_INFO_MOD
USE UNITS_MOD
USE WRITE_OUTPUT_MOD
USE MICROSTRUCTURE_MOD
USE ELEMENTAL_VARIABLES_UTILS_MOD
USE MATRIX_OPERATIONS_MOD, ONLY: CALC_ELVOL
!
IMPLICIT NONE
!
! Private data
!
INTEGER,  PRIVATE :: NLOADS
REAL(RK), ALLOCATABLE, PRIVATE :: END_LOAD(:, :)
INTEGER,  ALLOCATABLE, PRIVATE :: CONTROL_DIR(:)
REAL(RK), ALLOCATABLE, PRIVATE :: RAMP_RATE(:), DWELL_TIME(:)
REAL(RK), ALLOCATABLE, PRIVATE :: TARGET_TIME_INCR(:)
REAL(RK), ALLOCATABLE, PRIVATE :: LOAD_TOL(:)
INTEGER,  ALLOCATABLE, PRIVATE :: PRINT_FLAG(:)
!
CONTAINS
    !
    SUBROUTINE DRIVER_TRIAX_CLR(BCS, VELOCITY, PFORCE, DTRACE, NTRACE, &
        & C0_ANGS, CRSS_N, RSTAR_N, KEINV, WTS, AUTO_TIME, GAMMADOT)
    !
    !---------------------------------------------------------------------------
    !
    ! Driver for triaxial CLR simulations
    !
    ! Arguments:
    ! BCS:
    ! VELOCITY:
    ! PFORCE:
    ! DTRACE: gather/scatter trace for degrees of freedom
    ! NTRACE: gather/scatter trace for nodal points
    ! C0_ANGS:
    ! CRSS_N:
    ! RSTAR_N:
    ! KEINV:
    ! WTS:
    ! AUTO_TIME:
    ! GAMMADOT:
    !
    LOGICAL :: BCS(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: VELOCITY(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: PFORCE(DOF_SUB1:DOF_SUP1)
    TYPE(TRACE) :: DTRACE
    TYPE(TRACE) :: NTRACE
    REAL(RK) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: KEINV(0:TVEC1, 1:NUMPHASES)
    REAL(RK) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER :: AUTO_TIME
    REAL(RK) :: GAMMADOT(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    !
    ! Locals:
    INTEGER, PARAMETER :: LOADING = 1
    INTEGER, PARAMETER :: UNLOADING = 2
    INTEGER, PARAMETER :: DWELLING = 3
    !
    ! fraction of TARGET_TIME_INCR to use on first dwell increement
    REAL, PARAMETER :: INITIAL_TIME_FRAC = 0.01
    ! tolerance on dwell TIME expressed as a fraction of TARGET_TIME_INCR
    REAL, PARAMETER :: TIME_TOL = 0.001
    CHARACTER(LEN = 128) :: MESSAGE
    LOGICAL :: CONVERGED_SOLUTION
    INTEGER :: ITMETHOD_EVPS
    INTEGER :: INCR
    INTEGER :: ITERNL
    INTEGER :: IER
    INTEGER :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER :: IOSTATUS
    REAL(RK) :: DTIME, DTIME_STEP, TIME
    REAL(RK) :: D(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: D_VEC(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TRIAL_VELOCITY(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: VGRAD (0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: EVEL(0:kdim1, EL_SUB1:EL_SUP1)
    REAL(RK) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TRSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ELAS_TOT6(0:TVEC , 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: S_AVG_3X3(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ELPRESS(EL_SUB1:EL_SUP1)
    REAL(RK) :: E_ELAS_KK_BAR(EL_SUB1:EL_SUP1)
    REAL(RK) :: E_ELAS_KK(EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG_KK(EL_SUB1:EL_SUP1)
    REAL(RK) :: STIF(TVEC, TVEC, EL_SUB1:EL_SUP1)
    REAL(RK) :: FE(TVEC, EL_SUB1:EL_SUP1)
    INTEGER :: M_EL, I, J, K
    INTEGER :: ISTEP, IPRINT, IPHASE, IDIR
    INTEGER :: INDX(1:(DOF_SUP1-DOF_SUB1+1)/3)
    INTEGER :: INDY(1:(DOF_SUP1-DOF_SUB1+1)/3)
    INTEGER :: INDZ(1:(DOF_SUP1-DOF_SUB1+1)/3)
    REAL(RK) :: SURF_LOAD_ARRAY(NSURFACES,3)
    REAL(RK) :: AREA0(NSURFACES), AREA(NSURFACES)
    REAL(RK) :: LENGTH(3), LENGTH0(3)
    REAL(RK) :: MACRO_ENG_STRAIN(3)
    INTEGER  :: MAX_BC_ITER
    REAL(RK) :: MAX_STRAIN_INCR
    REAL(RK) :: MAX_STRAIN
    LOGICAL  :: FIRST_INCR_IN_STEP
    LOGICAL  :: START_RELOAD, START_UNLOAD, START_DWELL
    INTEGER  :: PREV_ACTION, CURR_ACTION
    INTEGER  :: NINCR_STEP, INCR_COUNT
    INTEGER  :: PERT_DIR
    INTEGER  :: BC_ITER_1, BC_ITER_2
    REAL(RK) :: STEP_FRAC
    REAL(RK) :: INITIAL_LOAD(3), PREV_LOAD(3), CURR_LOAD(3)
    REAL(RK) :: TARGET_LOAD(3), TRIAL_LOAD(3)
    REAL(RK) :: I_END_LOAD, I1_END_LOAD
    REAL(RK) :: DWELL_TIME_REMAINING
    REAL(RK) :: INITIAL_VEL(3), PREV_VEL(3), CURR_VEL(3)
    REAL(RK) :: INITIAL_LOAD_DWELL_VEL(3), INITIAL_UNLOAD_DWELL_VEL(3)
    REAL(RK) :: DELTA_VEL(3)
    REAL(RK) :: VEL_SCALE_FACTOR
    REAL(RK) :: PERT_MAG, PERT_SIGN
    REAL(RK) :: A(3,3)
    REAL(RK) :: B(3)
    REAL(RK) :: GAMMA(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DP_HAT(0:TVEC1,EL_SUB1:EL_SUP1)
    REAL(RK) :: WP_HAT(0:DIMS1,EL_SUB1:EL_SUP1)
    REAL(RK) :: DEFF(EL_SUB1:EL_SUP1)
    REAL(RK) :: DPEFF(EL_SUB1:EL_SUP1)
    REAL(RK) :: EQSTRAIN(EL_SUB1:EL_SUP1)
    REAL(RK) :: EQPLSTRAIN(EL_SUB1:EL_SUP1)
    REAL(RK) :: EAVG, NUAVG
    REAL(RK) :: MAX_EQSTRAIN
    REAL(RK) :: CURR_EQSTRAIN
    CHARACTER(LEN=16) :: FIELD, TIME_STRING, DTIME_STRING
    REAL(RK) :: DKK_STEP
    !
    TYPE(TIMERTYPE) :: DRIVER_TIMER, LIGHTUP_TIMER
    !
    REAL(RK), ALLOCATABLE :: ECOORDS(:, :)
    REAL(RK), ALLOCATABLE :: ELVOL(:)
    REAL(RK), ALLOCATABLE :: ELVOL_0(:)
    !
    !---------------------------------------------------------------------------
    !
    NLOADS = SIZE(END_LOAD, 1)
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    ! Extract optional parameters
    !
    MAX_BC_ITER = TRIAXCLR_OPTIONS%MAX_BC_ITER
    MAX_STRAIN_INCR = TRIAXCLR_OPTIONS%MAX_STRAIN_INCR
    MAX_STRAIN = TRIAXCLR_OPTIONS%MAX_STRAIN
    MAX_EQSTRAIN = TRIAXCLR_OPTIONS%MAX_EQSTRAIN
    !
    ! Indices for x, y, and z degrees of freedom
    !
    DO I = 1, ((DOF_SUP1-DOF_SUB1+1)/3)
        !
        INDX(I) = DOF_SUB1+3*(I-1)
        INDY(I) = DOF_SUB1+3*(I-1)+1
        INDZ(I) = DOF_SUB1+3*(I-1)+2
        !
    ENDDO
    !
    ! Restarting options
    !
    IF (OPTIONS%RESTART) THEN
        !
        CALL READ_RESTART_FIELD(VELOCITY, C0_ANGS, C_ANGS, &
            & RSTAR, RSTAR_N, WTS, CRSS, CRSS_N, &
            & E_ELAS_KK_BAR, SIG_VEC_N, EQSTRAIN, EQPLSTRAIN, GAMMA, &
            & EL_WORK_N, EL_WORKP_N, EL_WORK_RATE_N, EL_WORKP_RATE_N)
        !
        CALL READ_TRIAXCLR_RESTART(ISTEP, CURR_LOAD, PREV_LOAD, &
            & FIRST_INCR_IN_STEP, INCR, TIME, SURF_LOAD_ARRAY, AREA, AREA0, &
            & LENGTH, LENGTH0, CURR_VEL, PREV_ACTION, CURR_ACTION, &
            & INITIAL_LOAD_DWELL_VEL, INITIAL_UNLOAD_DWELL_VEL)
        !
    ELSE
        !
        ! Initialize areas and load arrays
        !
        SURF_LOAD_ARRAY = 0.0_RK
        CURR_LOAD = 0.0_RK
        AREA = 0.0_RK
        AREA0 = 0.0_RK
        !
        ! Compute initial area (AREA0)
        !
        CALL COMPUTE_AREA(COORDS, AREA0)
        !
        !IF (MYID .EQ. 0) WRITE(DFLT_U,'(a25,6f12.6)') 'initial areas (AREA0):   ', AREA0
        !
        ! Compute initial mesh dimensions
        !
        CALL CALC_MESH_DIM(LENGTH, INDX, INDY, INDZ)
        !
        LENGTH0 = LENGTH
        !
        ! Initialize state
        !
        E_ELAS_KK_BAR = 0.0_RK
        SIG_VEC_N = 0.0_RK
        C_ANGS = C0_ANGS
        !
        ! Initialize elvol arrays (and associated) iff it needs to be printed
        !
        IF ((PRINT_OPTIONS%PRINT_ELVOL) .OR. &
            & (FIBER_AVERAGE_OPTIONS%RUN_FIBER_AVERAGE)) THEN
            !
            ALLOCATE(ELVOL(EL_SUB1:EL_SUP1))
            ALLOCATE(ELVOL_0(EL_SUB1:EL_SUP1))
            ALLOCATE(ECOORDS(0:KDIM1, EL_SUB1:EL_SUP1))
            !
            ELVOL = 0.0_RK
            ELVOL_0 = 0.0_RK
            CALL PART_GATHER(ECOORDS, COORDS, NODES, DTRACE)
            CALL CALC_ELVOL(ELVOL_0, ECOORDS)
            !
        END IF
        !
        ! Initialize integrated quantities
        !
        EQPLSTRAIN = 0.0_RK
        EQSTRAIN = 0.0_RK
        GAMMA = 0.0_RK
        !
        ! Initialize deformation control
        !
        INCR = 1
        ISTEP = 1
        TIME = 0.0_RK
        FIRST_INCR_IN_STEP = .TRUE.
        PREV_ACTION = DWELLING
        CURR_ACTION = DWELLING
        INITIAL_LOAD_DWELL_VEL = 0.0_RK
        INITIAL_UNLOAD_DWELL_VEL = 0.0_RK
        CURR_EQSTRAIN = 0.0_RK
        !
        ! Print initial values
        !
        IF (PRINT_OPTIONS%PRINT_ELVOL) THEN
            !
            CALL PRINT_STEP(0, 0, COORDS, VELOCITY, C0_ANGS, WTS, CRSS_N, &
                & ELVOL_0)
            !
            DEALLOCATE(ELVOL_0)
            !
        ELSE
            !
            CALL PRINT_STEP(0, 0, COORDS, VELOCITY, C0_ANGS, WTS, CRSS_N)
            !
        END IF
        !
    ENDIF
    !
    ! Initialize load control
    !
    START_RELOAD = .FALSE.
    START_UNLOAD = .FALSE.
    START_DWELL = .FALSE.
    !
    ! Estimate bulk elastic moduli
    ! Assume uniform texture and equal element volumes
    !
    CALL EST_AVG_MOD(EAVG, NUAVG)
    ! IF (MYID .EQ. 0) WRITE(DFLT_U,'(2e16.8)') EAVG, NUAVG
    !
    ! Print headers to debug files - commented out to avoid standard printing
    ! CALL PRINT_HEADERS()
    !
    ! Debug printing
    ! WRITE(*,*) '1', END_LOAD(1,1), END_LOAD(1,2), END_LOAD(1,3), CONTROL_DIR(1), &
    ! & RAMP_RATE(1), DWELL_TIME(1), TARGET_TIME_INCR(1), LOAD_TOL(1), PRINT_FLAG(1)
    !
    ! Print headers to output files
    IF (MYID .EQ. 0) THEN
        !
        ! Write the header to file post.conv
        IF (PRINT_OPTIONS%PRINT_CONV) THEN
            !
            CALL WRITE_CONV_FILE_HEADERS()
            !
        END IF
        !
        ! Write the header to file post.force1-6
        IF (PRINT_OPTIONS%PRINT_FORCES) THEN
            !
            CALL WRITE_FORCE_FILE_HEADERS(2)
            !
        END IF
        !
    END IF
    !
    ! Start timers
    !
    ! CALL TIMERINIT(DRIVER_TIMER, 'DRIVER_TIMER', .TRUE.)
    ! CALL TIMERINIT(LIGHTUP_TIMER, 'LIGHTUP_TIMER', .TRUE.)
    ! CALL TIMERSTART(DRIVER_TIMER)
    !
    ! IF (MYID .EQ. 0) WRITE(DFLT_U,'(A)') 'Info   : Running step 1...'
    !
    ! Time stepping loop
    !
    TIME_STEPPING : DO WHILE (.TRUE.)
        !
        ! Iterate to find new configuration and material state.
        !
        CONVERGED_SOLUTION = .TRUE.
        !
        ! Initialize gammadots
        !
        GAMMADOT = 0.0_RK
        !
        ! If first increment in step, calculate time increment, number of
        !   increments in the step, and initialize increment counter. Check for
        !   start of unloading, reloading, and dwell episodes.
        !
        IF (FIRST_INCR_IN_STEP) THEN
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A,I0,A)') 'Info   : Running step ', ISTEP, '...'
                !
            ENDIF
            !
            ! Get end loads for steps i and i-1
            !
            IDIR = ABS(CONTROL_DIR(ISTEP))
            I_END_LOAD = END_LOAD(ISTEP,IDIR)
            !
            IF (ISTEP .EQ. 1) THEN
                !
                I1_END_LOAD = 0.0_RK
                !
            ELSE
                !
                I1_END_LOAD = END_LOAD(ISTEP-1, IDIR)
                !
            ENDIF
            !
            ! Get previous and current actions
            !
            PREV_ACTION = CURR_ACTION
            !
            IF (CONTROL_DIR(ISTEP) .LT. 0) THEN
                !
                CURR_ACTION = DWELLING
                !
            ELSEIF (I_END_LOAD .GT. I1_END_LOAD) THEN
                !
                CURR_ACTION = LOADING
                !
            ELSEIF (I_END_LOAD .LT. I1_END_LOAD) THEN
                !
                CURR_ACTION = UNLOADING
                !
            ELSE
                !
                CALL PAR_QUIT('Error  :     > Invalid load sequence')
                !
            ENDIF
            !
            ! Check for start of loading episode
            !
            IF ((CURR_ACTION .EQ. LOADING) .AND. (PREV_ACTION .NE. LOADING) &
                & .AND. (ISTEP .NE. 1)) THEN
                !
                START_RELOAD = .TRUE.
                !
                IF (MYID .EQ. 0) THEN
                    !
                    WRITE(DFLT_U,'(A)') 'Info   :   - Starting reload episode'
                    !
                ENDIF
                !
            ENDIF
            !
            ! Check for start of unloading episode
            !
            IF ((CURR_ACTION .EQ. UNLOADING) .AND. &
                & (PREV_ACTION .NE. UNLOADING)) THEN
                !
                START_UNLOAD = .TRUE.
                !
                IF (MYID .EQ. 0) THEN
                    !
                    WRITE(DFLT_U,'(A)') 'Info   :   - Starting unload episode'
                    !
                ENDIF
                !
            ENDIF
            !
            ! Check for start of dwell episode
            !
            IF ((CURR_ACTION .EQ. DWELLING) .AND. (PREV_ACTION .NE. DWELLING)) &
                & THEN
                !
                START_DWELL = .TRUE.
                !
                IF (MYID .EQ. 0) THEN
                    !
                    WRITE(DFLT_U,'(A)') 'Info   :   - Starting dwell episode'
                    !
                ENDIF
                !
            ENDIF
            !
            ! Calculate time step and time increment
            !
            IF (CURR_ACTION .EQ. DWELLING) THEN
                !
                DTIME_STEP = DWELL_TIME(ISTEP)
                DWELL_TIME_REMAINING = DTIME_STEP
                !
                ! NINCR_STEP, DTIME, and INCR_COUNT are updated on each
                !   increment of the dwell episode
                !
                NINCR_STEP = 1
                !
                IF (START_DWELL) THEN
                    !
                    DTIME = INITIAL_TIME_FRAC * TARGET_TIME_INCR(ISTEP)
                    !
                ELSE
                    !
                    DTIME = MIN(MAX_STRAIN_INCR * LENGTH(IDIR) / &
                        & CURR_VEL(IDIR), TARGET_TIME_INCR(ISTEP), &
                        & DWELL_TIME_REMAINING)
                    !
                ENDIF
                !
                INCR_COUNT = 1
                !
            ELSE
                !
                DTIME_STEP = ABS((I_END_LOAD - I1_END_LOAD) / RAMP_RATE(ISTEP))
                !
                ! NINCR_STEP and DTIME are static for each load/unload step
                !
                NINCR_STEP = MAX(NINT(DTIME_STEP / TARGET_TIME_INCR(ISTEP)), 1)
                DTIME = DTIME_STEP / NINCR_STEP
                INCR_COUNT = 1
                DWELL_TIME_REMAINING = 0.0_RK
                !
            ENDIF
            !
        ENDIF ! First increment in step
        !
        ! Update TIME
        !
        TIME = TIME + DTIME
        !
        ! Calculate target loads
        !
        STEP_FRAC = REAL(INCR_COUNT) / REAL(NINCR_STEP)
        !
        IF (ISTEP .EQ. 1) THEN
            !
            TARGET_LOAD = END_LOAD(1, :) * STEP_FRAC
            !
        ELSE
            !
            TARGET_LOAD = END_LOAD(ISTEP, :) * STEP_FRAC + &
                & END_LOAD(ISTEP-1, :) * (1.0-STEP_FRAC)
            !
        ENDIF
        !
        IF (MYID .EQ. 0) THEN
            !
            WRITE(FIELD, '(F0.6)') TIME
            IF (FIELD(1:1) == '.') FIELD = '0' // FIELD
            TIME_STRING = FIELD
            WRITE(FIELD, '(F0.6)') DTIME
            IF (FIELD(1:1) == '.') FIELD = '0' // FIELD
            DTIME_STRING = FIELD
            !               
            WRITE(DFLT_U,'(A,I0,A,A,A,A,A)') 'Info   :   - &
                &Increment ', INCR, ': t = ', TRIM(TIME_STRING), '&
                & secs, dt = ', TRIM(DTIME_STRING), ' secs'
            !
        ENDIF
        !
        ! Iterate on surface velocities (single velocity mode - fixed from previous
        !   time increment)
        !
        ! Store initial load and update previous load
        !
        INITIAL_LOAD = CURR_LOAD
        PREV_LOAD = CURR_LOAD
        !
        ! Initialize previous velocity
        !
        PREV_VEL = 0
        !
        ! For starting reload apply zero velocity
        !
        IF (START_RELOAD) THEN
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A)') 'Info   :   - Zero-velocity iteration'
                !
            ENDIF
            !
            VELOCITY = 0.0
            !
            CALL VELOCITY_ITERATION(BCS, PFORCE, VELOCITY, ELPRESS, EVEL, &
                & CURR_LOAD, DTRACE, NTRACE, C0_ANGS, C_ANGS, SIG_VEC_N, &
                & SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, TRSTAR, KEINV, &
                & E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, JITER_STATE, WTS, DEFF, &
                & DTIME, INCR, E_BAR_VEC, CONVERGED_SOLUTION, AUTO_TIME, ITERNL)
            !
        ENDIF
        !
        ! Initialize vel for first time increment or start of reload episode
        !
        IF ((INCR .EQ. 1) .OR. START_RELOAD) THEN
            !
            CURR_VEL(1) = ((TARGET_LOAD(1) - CURR_LOAD(1)) / AREA0(4) &
                & - NUAVG * (TARGET_LOAD(2) - CURR_LOAD(2)) / AREA0(6) &
                & - NUAVG * (TARGET_LOAD(3) - CURR_LOAD(3)) / AREA0(2)) / &
                & EAVG * LENGTH(1) / DTIME
            CURR_VEL(2) = (-NUAVG * (TARGET_LOAD(1) - CURR_LOAD(1)) / AREA0(4) &
                & + (TARGET_LOAD(2) - CURR_LOAD(2)) / AREA0(6) &
                & - NUAVG * (TARGET_LOAD(3) - CURR_LOAD(3)) / AREA0(2)) / &
                & EAVG * LENGTH(2) / DTIME
            CURR_VEL(3) = (-NUAVG * (TARGET_LOAD(1) - CURR_LOAD(1)) / AREA0(4) &
                & - NUAVG * (TARGET_LOAD(2) - CURR_LOAD(2)) / AREA0(6) &
                & + (TARGET_LOAD(3) - CURR_LOAD(3)) / AREA0(2)) / EAVG * &
                & LENGTH(3) / DTIME
            !
            VELOCITY(INDX) = CURR_VEL(1) * COORDS(INDX) / LENGTH(1)
            VELOCITY(INDY) = CURR_VEL(2) * COORDS(INDY) / LENGTH(2)
            VELOCITY(INDZ) = CURR_VEL(3) * COORDS(INDZ) / LENGTH(3)
            !
        ENDIF
        !
        ! Initialize velocity for start of unload - back off yield surface
        !
        IF (START_UNLOAD) THEN
            !
            CURR_VEL = CURR_VEL / 10.0_RK
            VELOCITY = VELOCITY / 10.0_RK
            !
        ENDIF
        !
        ! Initialize velocity for start of dwell episode
        !
        IF (START_DWELL) THEN
            !
            IF (PREV_ACTION .EQ. LOADING) THEN
                !
                CURR_VEL = INITIAL_LOAD_DWELL_VEL
                !
            ELSEIF (PREV_ACTION .EQ. UNLOADING) THEN
                !
                CURR_VEL = INITIAL_UNLOAD_DWELL_VEL
                !
            ENDIF
            !
            VELOCITY(INDX) = CURR_VEL(1) * COORDS(INDX) / LENGTH(1)
            VELOCITY(INDY) = CURR_VEL(2) * COORDS(INDY) / LENGTH(2)
            VELOCITY(INDZ) = CURR_VEL(3) * COORDS(INDZ) / LENGTH(3)
            !
        ENDIF
        !
        ! Store the velocity at start of increment
        !
        INITIAL_VEL = CURR_VEL
        !
        ! Initialize iteration counter
        !
        BC_ITER_1 = 0
        !
        BC_ITERATION_1 : DO WHILE (.TRUE.)
            !
            BC_ITER_1 = BC_ITER_1 + 1
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U, '(A,I0)') 'Info   :     > TRIAL_BC1: Iteration ', BC_ITER_1
                !
            ENDIF
            !
            CALL VELOCITY_ITERATION(BCS, PFORCE, VELOCITY, ELPRESS, EVEL, &
                & CURR_LOAD, DTRACE, NTRACE, C0_ANGS, C_ANGS, SIG_VEC_N, &
                & SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, TRSTAR, KEINV, &
                & E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, JITER_STATE, WTS, DEFF, &
                & DTIME, INCR, E_BAR_VEC, CONVERGED_SOLUTION, AUTO_TIME, ITERNL)
            !
            ! Print BCS_ITER data
            !
            IF (MYID .EQ. 0) THEN
                !
                ! For debug purpose only:
                ! WRITE(OUNITS(BCS_ITER_1_U),'(3(i8),10(e15.5))') ISTEP, INCR, &
                ! & BC_ITER_1, CURR_VEL, CURR_LOAD, TARGET_LOAD, DTIME
                !
                IF (MINVAL(ABS(CURR_VEL)) .LE. 0.0010) THEN
                    !
                    WRITE(DFLT_U,'(A,3(E12.2))') 'Info   :       . Velocity:    ', CURR_VEL
                    !
                ELSE
                    !
                    WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       . Velocity:    ', CURR_VEL
                    !
                END IF
                !
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       . Load:        ', CURR_LOAD
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       . Target load: ', TARGET_LOAD
                !
            ENDIF
            !
            ! Check if force is within range
            !
            IF ((ABS(CURR_LOAD(IDIR)-TARGET_LOAD(IDIR)) .LT. &
                & MAX(LOAD_TOL(ISTEP), 0.1_RK * ABS(TARGET_LOAD(IDIR) - &
                & INITIAL_LOAD(IDIR)))) .OR. START_DWELL .OR. START_RELOAD &
                & .OR. (CURR_ACTION .EQ. UNLOADING)) THEN
                !
                ! Within tolerance - exit loop
                !
                EXIT
                !
            ELSE
                !
                ! Not within tolerance - update guess for surface velocities
                !
                IF (BC_ITER_1 .EQ. 1) THEN
                    !
                    IF (CURR_ACTION .EQ. DWELLING) THEN
                        !
                        VEL_SCALE_FACTOR = 0.9_RK
                        !
                    ELSE
                        !
                        VEL_SCALE_FACTOR = 1.1_RK
                        !
                    ENDIF
                    !
                ELSE
                    !
                    VEL_SCALE_FACTOR = (1 - PREV_VEL(IDIR) / CURR_VEL(IDIR)) * &
                       & (TARGET_LOAD(IDIR) - PREV_LOAD(IDIR)) / &
                       & (CURR_LOAD(IDIR) - PREV_LOAD(IDIR)) + &
                       & PREV_VEL(IDIR) / CURR_VEL(IDIR)
                    !
                ENDIF
                !
                PREV_LOAD = CURR_LOAD
                PREV_VEL  = CURR_VEL
                CURR_VEL = CURR_VEL*VEL_SCALE_FACTOR
                VELOCITY = VELOCITY*VEL_SCALE_FACTOR
                !
            ENDIF
            !
            ! Check if maximum number of iteration are exceeded
            !
            IF (BC_ITER_1 .EQ. MAX_BC_ITER) THEN
                !
                CALL PAR_QUIT ('Error  :     > Maximum number of boundary&
                    & condition iterations exceeded.')
                !
            ENDIF
            !
        ENDDO BC_ITERATION_1
        !
        ! Iterate on surface velocities (all three velocity modes)
        !
        ! Initialize iterataion counter
        !
        BC_ITER_2 = 0
        !
        ! Calculate perturbation magnitude
        !
        IF (START_DWELL) THEN
            !
            PERT_MAG = 1e-6
            !
        ELSE
            !
            PERT_MAG = MAX(0.1_RK * ABS(CURR_VEL(IDIR) - INITIAL_VEL(IDIR)), &
                & 0.01_RK * ABS(CURR_VEL(IDIR)))
            !
        ENDIF
        !
        BC_ITERATION_2 : DO WHILE (.TRUE.)
            !
            BC_ITER_2 = BC_ITER_2 + 1
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A,I0)') 'Info   :     > TRIAL_BC2: Iteration ', BC_ITER_2
                !
            ENDIF
            !
            ! Perturb velocity
            !
            DO PERT_DIR = 1, DIMS
                !
                IF (MYID .EQ. 0) THEN
                    !
                    SELECT CASE (PERT_DIR)
                        !
                        CASE(1)
                            !
                            WRITE(DFLT_U,'(A)') 'Info   :       . Perturbing x-velocity'
                            !
                        CASE(2)
                            !
                            WRITE(DFLT_U,'(A)') 'Info   :       . Perturbing y-velocity'
                            !
                        CASE(3)
                            !
                            WRITE(DFLT_U,'(A)') 'Info   :       . Perturbing z-velocity'
                        !
                    END SELECT
                    !
                ENDIF
                !
                ! Apply perturbation
                !
                TRIAL_VELOCITY = VELOCITY
                !
                IF (CURR_LOAD(PERT_DIR) .LT. TARGET_LOAD(PERT_DIR)) THEN
                    !
                    PERT_SIGN = 1.0_RK
                    !
                ELSE
                    !
                    PERT_SIGN = -1.0_RK
                    !
                ENDIF
                !
                SELECT CASE (PERT_DIR)
                    !
                    CASE (1)
                        !
                        TRIAL_VELOCITY(INDX) = VELOCITY(INDX) + PERT_SIGN * &
                            & PERT_MAG * COORDS(INDX) / LENGTH(1)
                        !
                    CASE (2)
                        !
                        TRIAL_VELOCITY(INDY) = VELOCITY(INDY) + PERT_SIGN * &
                            & PERT_MAG * COORDS(INDY) / LENGTH(2)
                        !
                    CASE (3)
                        !
                        TRIAL_VELOCITY(INDZ) = VELOCITY(INDZ) + PERT_SIGN * &
                            & PERT_MAG * COORDS(INDZ) / LENGTH(3)
                    !
                END SELECT
                !
                CALL VELOCITY_ITERATION(BCS, PFORCE, TRIAL_VELOCITY, ELPRESS, &
                    & EVEL, TRIAL_LOAD, DTRACE, NTRACE, C0_ANGS, C_ANGS, &
                    & SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, &
                    & TRSTAR, KEINV, E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, &
                    & JITER_STATE, WTS, DEFF, DTIME, INCR, E_BAR_VEC, &
                    & CONVERGED_SOLUTION, AUTO_TIME, ITERNL)
                !
                IF (MYID .EQ. 0) THEN
                    !
                    WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :     > Trial load: ', TRIAL_LOAD
                    !
                ENDIF
                !
                ! Compute matrix coefficients
                !
                A(1, PERT_DIR) = (TRIAL_LOAD(1) - CURR_LOAD(1)) / &
                    & (PERT_SIGN * PERT_MAG)
                A(2, PERT_DIR) = (TRIAL_LOAD(2) - CURR_LOAD(2)) / &
                    & (PERT_SIGN * PERT_MAG)
                A(3, PERT_DIR) = (TRIAL_LOAD(3) - CURR_LOAD(3)) / &
                    & (PERT_SIGN * PERT_MAG)
                !
            ENDDO ! Perturbate velocity
            !
            ! Update velocity field
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A)') 'Info   :     > Updating velocity field'
                !
            ENDIF
            !
            ! Solve for surface velocity increment
            !
            ! From right-hand side array
            !
            B = TARGET_LOAD - CURR_LOAD
            !
            ! Solve system of equations
            !
            CALL SOLVE_LIN_SYS_3(A, B, DELTA_VEL)
            !
            ! Apply new velocity boundary conditions
            !
            VELOCITY(INDX) = VELOCITY(INDX) + DELTA_VEL(1) * COORDS(INDX) / &
                & LENGTH(1)
            VELOCITY(INDY) = VELOCITY(INDY) + DELTA_VEL(2) * COORDS(INDY) / &
                & LENGTH(2)
            VELOCITY(INDZ) = VELOCITY(INDZ) + DELTA_VEL(3) * COORDS(INDZ) / &
                & LENGTH(3)
            !
            CURR_VEL = CURR_VEL + DELTA_VEL
            !
            CALL VELOCITY_ITERATION(BCS, PFORCE, VELOCITY, ELPRESS, EVEL, &
                & CURR_LOAD, DTRACE, NTRACE, C0_ANGS, C_ANGS, SIG_VEC_N, &
                & SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, TRSTAR, KEINV, &
                & E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, JITER_STATE, WTS, DEFF, &
                & DTIME, INCR, E_BAR_VEC, CONVERGED_SOLUTION, AUTO_TIME, ITERNL)
            !
            ! Print BCS_ITER data
            !
            IF (MYID .EQ. 0) THEN
                !
                ! For debug purposes only
                !WRITE(OUNITS(BCS_ITER_2_U),'(3(i8),11(e15.5))') ISTEP, INCR, &
                !    & BC_ITER_2, CURR_VEL, CURR_LOAD, TARGET_LOAD, DTIME, &
                !    & PERT_MAG
                !
                IF (MINVAL(ABS(CURR_VEL)) .LE. 0.0010) THEN
                    !
                    WRITE(DFLT_U,'(A,3(E12.2))') 'Info   :       . Velocity:    ', CURR_VEL
                    !
                ELSE
                    !
                    WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       . Velocity:    ', CURR_VEL
                    !
                END IF
                !
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       . Load:        ', CURR_LOAD
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       . Target load: ', TARGET_LOAD
                !
            ENDIF
            !
            ! Check if force is within range
            !
            IF (MAXVAL(ABS(CURR_LOAD-TARGET_LOAD)) .LT. LOAD_TOL(ISTEP)) THEN
                !
                ! Within tolerance - exit loop
                !
                EXIT
                !
            ENDIF
            !
            ! Check if maximum number of iteration are exceeded
            !
            IF (BC_ITER_2 .EQ. MAX_BC_ITER) THEN
                !
                CALL PAR_QUIT('Error  :     > Maximum number of boundary&
                    & condition iterations exceeded')
                !
            ENDIF
            !
            ! Update perturbation magnitude
            !
            PERT_MAG = MAXVAL(ABS(DELTA_VEL))
            !
        ENDDO BC_ITERATION_2
        !
        ! Update RSTAR and SIG_VEC_N (?)
        !
        RSTAR = TRSTAR
        !
        ! Write boundary condition iteration statistics - debug only
        !IF (MYID .EQ. 0) THEN
        !    !
        !    WRITE(OUNITS(BCS_ITER_LOG_U),'(4(i12))') ISTEP, INCR, BC_ITER_1, &
        !        & BC_ITER_2
        !    !
        !ENDIF
        !
        ! Update state (at center of each element), geometry and tool bc's
        !
        ! Update:
        ! - COORDS @(t+dt), using DT and VELOCITY @(t+dt)
        !
        ! element centroid:
        ! - SIG_VEC_N: 5-vec deviatoric kirchhoff stress @(t)
        ! - SIG_VEC: 5-vec deviatoric kirchhoff stress @(t+dt)
        ! - SIG_KK: volumetric kirchhoff stress @(t+dt)
        ! - E_ELAS_KK_BAR: volumetric lattice strain @(t)
        ! - E_ELAS_KK: volumetric lattice strain @(t+dt)
        !
        ! These are outputs of the subroutine/function but are useless and can
        !   be removed from the argument list
        ! E_BAR_VEC: OUT
        ! E_ELAS_KK_BAR: OUT
        ! SIG_VEC_N: OUT
        !
        CALL UPDATE_STATE_EVPS(VELOCITY, DTRACE, C0_ANGS, C_ANGS, SIG_VEC_N, &
            & SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, E_BAR_VEC, E_ELAS_KK_BAR, &
            & E_ELAS_KK, SIG_KK, KEINV, WTS, DEFF, DTIME, STIF, FE, D, D_VEC, &
            & VGRAD, CONVERGED_SOLUTION, AUTO_TIME)
        !
        ! Calculate stress, strain, load, and area
        !
        CALL CALC_STRESS_STRAIN(S_AVG_3X3, ELAS_TOT6, SURF_LOAD_ARRAY, AREA, &
            & SIG_VEC, SIG_KK, E_ELAS_KK, C_ANGS, KEINV, WTS)
        !
        CURR_LOAD(1) = SURF_LOAD_ARRAY(2,1) ! LEGACY FORMAT: 4,1
        CURR_LOAD(2) = SURF_LOAD_ARRAY(4,2) !                6,2
        CURR_LOAD(3) = SURF_LOAD_ARRAY(6,3) !                2,3
        !
        ! Calculate mesh dimensions
        !
        CALL CALC_MESH_DIM(LENGTH, INDX, INDY, INDZ)
        !
        ! Calculate macroscopic engineering strain
        !
        DO I = 1, 3
            !
            MACRO_ENG_STRAIN(I) = LENGTH(I) / LENGTH0(I) - 1.0_RK
            !
        ENDDO
        !
        ! Calculate the current increment macroscopic eqstrain
        !
        CURR_EQSTRAIN = (2. / 3. ) * SQRT( (3 * ((MACRO_ENG_STRAIN(1) ** 2) &
            & + (MACRO_ENG_STRAIN(2) ** 2) + (MACRO_ENG_STRAIN(3) ** 2))/2))
        !
        ! Print the macroscopic strain values to console for monitoring
        !
        IF (MYID .EQ. 0) THEN
            !
            WRITE(DFLT_U, '(A,F12.6)') 'Info   :       . Maximum eng. strain: ',&
                &maxval(abs(MACRO_ENG_STRAIN(:)))
            WRITE(DFLT_U, '(A,F12.6)') 'Info   :       . Current eqv. strain: ',&
                &CURR_EQSTRAIN
            !
        END IF
        !
        ! Calculate DWELL_TIME_REMAINING
        !
        IF (CURR_ACTION .EQ. DWELLING) THEN
            !
            DWELL_TIME_REMAINING = DWELL_TIME_REMAINING - DTIME
            !
        ENDIF
        !
        ! Write loads to post.force# files
        !
        IF (MYID .EQ. 0) THEN
            !
            IF (PRINT_OPTIONS%PRINT_FORCES) THEN
                !
                CALL WRITE_FORCE_FILE_DATA(2, ISTEP, INCR, &
                    & SURF_LOAD_ARRAY, AREA, TIME, LENGTH)
                !
            END IF
            !
        ENDIF
        !
        ! Calculate DPEFF
        !
        CALL PLASTICVELGRADSYMSKW(DP_HAT, WP_HAT, DPEFF, GAMMADOT(:, 0, :), &
            & M_EL)
        !
        ! Update EQPLSTRAIN
        !
        EQPLSTRAIN = EQPLSTRAIN + DPEFF * DTIME
        !
        ! Update EQSTRAIN
        !
        EQSTRAIN = EQSTRAIN + DEFF * DTIME
        !
        ! Reconstruct total deformation rate tensor if requested
        !
        IF ((PRINT_OPTIONS%PRINT_WORK) .OR. (PRINT_OPTIONS%PRINT_DEFRATE) .OR. &
            & (PRINT_OPTIONS%PRINT_RESTART)) THEN
            !
            ! Copy current deviatoric portion to new array
            !
            D_TOT = D
            ! 
            ! Reconstruct the total deformation rate tensor
            !
            DO I = EL_SUB1, EL_SUP1
                !
                DKK_STEP = VGRAD(0, 0, I) + VGRAD(1, 1, I) + VGRAD(2, 2, I)
                D_TOT(0, 0, I) = D_TOT(0, 0, I) + DKK_STEP / 3.0
                D_TOT(1, 1, I) = D_TOT(1, 1, I) + DKK_STEP / 3.0
                D_TOT(2, 2, I) = D_TOT(2, 2, I) + DKK_STEP / 3.0
                !
            END DO
            !
        END IF
        !
        ! Calculate work, plastic work
        !
        IF ((PRINT_OPTIONS%PRINT_WORK) .OR. (PRINT_OPTIONS%PRINT_RESTART)) THEN
            !
            CALL CALC_TOTAL_WORK(DTIME, D_TOT, S_AVG_3X3, EL_WORK_N, &
                & EL_WORK_RATE_N, EL_WORK)
            !
        END IF
        !
        IF ((PRINT_OPTIONS%PRINT_WORKP) .OR. (PRINT_OPTIONS%PRINT_RESTART)) THEN
            !
            CALL CALC_PLASTIC_WORK(DTIME, DP_HAT, C_ANGS, S_AVG_3X3, &
                & EL_WORKP_N, EL_WORKP_RATE_N, EL_WORKP)
            !
        END IF
        !
        IF ((PRINT_OPTIONS%PRINT_ELVOL) .OR. &
            & (FIBER_AVERAGE_OPTIONS%RUN_FIBER_AVERAGE)) THEN
            !
            CALL PART_GATHER(ECOORDS, COORDS, NODES, DTRACE)
            CALL CALC_ELVOL(ELVOL, ECOORDS)
            !
        END IF
        !
        ! Update slip system shear
        !
        GAMMA = GAMMA + GAMMADOT * DTIME
        !
        ! Store velocity if at beginning of dwell episode
        !
        IF (START_DWELL) THEN
            !
            IF (PREV_ACTION .EQ. LOADING) THEN
                !
                INITIAL_LOAD_DWELL_VEL = CURR_VEL
                !
            ELSEIF (PREV_ACTION .EQ. UNLOADING) THEN
                !
                INITIAL_UNLOAD_DWELL_VEL = CURR_VEL
                !
            ENDIF
            !
        ENDIF
        !
        ! Output computed quantities at end of step
        !
        IF (((CURR_ACTION .NE. DWELLING) .AND. (INCR_COUNT .EQ. NINCR_STEP)) &
            & .OR. ((CURR_ACTION .EQ. DWELLING) .AND. (DWELL_TIME_REMAINING &
            & .LT. TIME_TOL * TARGET_TIME_INCR(ISTEP))))  THEN
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A,I0)') 'Info   :     > Converged on increment ', INCR
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :     > Loads on X0: ', SURF_LOAD_ARRAY(1,:)
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       Loads on X1: ', SURF_LOAD_ARRAY(2,:)
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       Loads on Y0: ', SURF_LOAD_ARRAY(3,:)
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       Loads on Y1: ', SURF_LOAD_ARRAY(4,:)
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       Loads on Z0: ', SURF_LOAD_ARRAY(5,:)
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       Loads on Z1: ', SURF_LOAD_ARRAY(6,:)
                !
            ENDIF
            !
            IF (PRINT_FLAG(ISTEP) .EQ. 0) THEN
                !
                CALL PRINT_STEP(ISTEP, ANISOTROPIC_EVPS, COORDS, VELOCITY, &
                    & C_ANGS, WTS, CRSS, ELVOL, ELAS_TOT6, S_AVG_3X3, DEFF, &
                    & EQSTRAIN, DPEFF, EQPLSTRAIN, VGRAD, DP_HAT, WP_HAT, &
                    & GAMMA, GAMMADOT, EL_WORK, EL_WORKP, D_TOT)
                !
                IF (PRINT_OPTIONS%PRINT_RESTART) THEN
                    !
                    CALL WRITE_RESTART_FIELD(VELOCITY, C0_ANGS, C_ANGS, RSTAR, &
                        & RSTAR_N, WTS, CRSS, CRSS_N, E_ELAS_KK_BAR, &
                        & SIG_VEC_N, EQSTRAIN, EQPLSTRAIN, GAMMA, EL_WORK_N, &
                        & EL_WORKP_N, EL_WORK_RATE_N, EL_WORKP_RATE_N)
                    !
                    CALL WRITE_TRIAXCLR_RESTART(ISTEP + 1, CURR_LOAD, &
                        & PREV_LOAD, .TRUE., INCR + 1, TIME, SURF_LOAD_ARRAY, &
                        & AREA, AREA0, LENGTH, LENGTH0, CURR_VEL, PREV_ACTION, &
                        & CURR_ACTION, INITIAL_LOAD_DWELL_VEL, &
                        & INITIAL_UNLOAD_DWELL_VEL)
                    !
                ENDIF
                !
                IF (FIBER_AVERAGE_OPTIONS%RUN_FIBER_AVERAGE) THEN
                    !
                    ! CALL TIMERSTART(LIGHTUP_TIMER)
                    CALL RUN_FIBER_AVERAGE(ISTEP, C_ANGS, ELAS_TOT6, DPEFF, &
                        & CRSS, ELVOL)
                    ! CALL TIMERSTOP(LIGHTUP_TIMER)
                    !
                ENDIF
                !
            ENDIF
            !
            IF (ISTEP .EQ. NLOADS) THEN
                !
                ! Job finished successfully
                !
                ! CALL TIMERSTOP(DRIVER_TIMER)
                !
                !IF (MYID .EQ. 0) THEN
                    !
                    ! CALL TIMERWRITE(DRIVER_TIMER, DFLT_U)
                    ! CALL TIMERWRITE(LIGHTUP_TIMER, DFLT_U)
                    !
                !ENDIF
                !
                CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP, PRINT_FLAG)
                !
                CALL PAR_QUIT('Info   : Final step terminated. Simulation&
                    & completed successfully.')
                !
            ENDIF
            !
            ISTEP = ISTEP + 1
            FIRST_INCR_IN_STEP = .TRUE.
            INCR_COUNT = 1
            !
        ELSE
            !
            FIRST_INCR_IN_STEP = .FALSE.
            INCR_COUNT = INCR_COUNT + 1
            !
            ! Compute new time increment for dwell
            !
            IF (CURR_ACTION .EQ. DWELLING) THEN
                !
                DTIME = MIN(MAX_STRAIN_INCR*LENGTH(IDIR)/ABS(CURR_VEL(IDIR)), &
                    &TARGET_TIME_INCR(ISTEP), DWELL_TIME_REMAINING)
                NINCR_STEP = NINCR_STEP + 1
                !
            ENDIF
            !
        ENDIF
        !
        ! Check that maximum strain has not been exceeded
        !
        IF (MAXVAL(ABS(MACRO_ENG_STRAIN)) .GE. MAX_STRAIN) THEN
            !
            CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP, PRINT_FLAG)
            !
            CALL PAR_QUIT('Error  :     > Maximum eng. strain exceeded.')
            !
        ENDIF
        !
        ! Check that maximum eqstrain has not been exceeded
        !
        IF (CURR_EQSTRAIN .GE. MAX_EQSTRAIN) THEN
            !
            CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP, PRINT_FLAG)
            !
            CALL PAR_QUIT('Error  :     > Maximum eqv. strain exceeded.')
            !
        ENDIF
        !
        INCR = INCR + 1
        START_RELOAD = .FALSE.
        START_UNLOAD = .FALSE.
        START_DWELL  = .FALSE.
        !
    ENDDO TIME_STEPPING
    !
    RETURN
    !
    END SUBROUTINE DRIVER_TRIAX_CLR
    !
    !===========================================================================
    !
    SUBROUTINE READ_CTRL_DATA_CLR
    !
    ! Process input data for load control at constant load rate.
    !
    !---------------------------------------------------------------------------
    !
    ! Locals:
    ! I/II: Generic loop index values.
    ! ILOAD: Index value used to access TARGET_CLR_LOAD sequentially by row.
    ! IDWELL: Index value used to access DWELL_EPISODE sequentially by row.
    ! STEP_RAMP_RATE: Array used to store ramp rate jump values.
    !
    INTEGER :: NLOADS, I, II, ILOAD, IDWELL
    REAL(RK), ALLOCATABLE :: STEP_RAMP_RATE(:)
    !
    ! Notes:
    ! User-defined input (ramp rates) should ALWAYS be positive. Compression 
    ! and similar signed deformation control are handled internally within the
    ! primary driver subroutine.
    !
    !---------------------------------------------------------------------------
    !
    ! Set NLOADS to account for possible dwell episodes to adjust the array
    !
    NLOADS = TRIAXCLR_OPTIONS%NUMBER_OF_CLR_LOAD_STEPS + &
        & TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES
    !
    ! Allocate the necessary arrays
    ALLOCATE(END_LOAD(NLOADS,DIMS))
    ALLOCATE(CONTROL_DIR(NLOADS))
    ALLOCATE(RAMP_RATE(NLOADS))
    ALLOCATE(DWELL_TIME(NLOADS))
    ALLOCATE(TARGET_TIME_INCR(NLOADS))
    ALLOCATE(LOAD_TOL(NLOADS))
    ALLOCATE(PRINT_FLAG(NLOADS))
    !
    ALLOCATE(STEP_RAMP_RATE(TRIAXCLR_OPTIONS%NUMBER_OF_CLR_LOAD_STEPS))
    !
    ! Initialize arrays.
    !
    CONTROL_DIR = BCS_OPTIONS%LOADING_DIRECTION + 1
    RAMP_RATE = BCS_OPTIONS%LOAD_RATE
    DWELL_TIME = 0.0_RK
    TARGET_TIME_INCR = 0.0_RK
    LOAD_TOL = TRIAXCLR_OPTIONS%LOAD_TOL_ABS
    !
    STEP_RAMP_RATE = BCS_OPTIONS%LOAD_RATE
    !    
    ! Build load rate array, specifically handle load rate jumps.
    !
    IF (TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS .GT. 0) THEN
        !
        II = 1
        !
        DO I = 1, TRIAXCLR_OPTIONS%NUMBER_OF_CLR_LOAD_STEPS
            !
            IF (TRIAXCLR_OPTIONS%LOAD_RATE_JUMP(II,1) .EQ. I) THEN
                !
                STEP_RAMP_RATE(I) = &
                    & TRIAXCLR_OPTIONS%LOAD_RATE_JUMP(II,2)
                !
                IF (II.EQ.TRIAXCLR_OPTIONS%NUMBER_OF_LOAD_RATE_JUMPS) THEN
                    !
                    II = II
                    !           
                ELSE
                    !
                    II = II + 1
                    !
                ENDIF
                !
            ENDIF
            !
        END DO
        !
    ENDIF
    !
    ! Check that a dwell doesn't occur before the first loading step
    !
    IF (TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES .GT. 0) THEN
        !
        IF (TRIAXCLR_OPTIONS%DWELL_EPISODE(1,1) .EQ. 1) THEN
            !
            CALL PAR_QUIT('Error  :     > A dwell episode can not be applied to &
                &an unloaded specimen on the first step.')
            !
        ENDIF
        !
    ENDIF
    !
    ! Increment the first column (step index) of dwell episodes so it can be 
    ! accessed in a logical order to build the expected arrays properly
    !
    IF (TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES .GT. 0) THEN
        !
        DO I = 1, TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES
            !
            TRIAXCLR_OPTIONS%DWELL_EPISODE(I,1) = &
                & TRIAXCLR_OPTIONS%DWELL_EPISODE(I,1) + 1
            !
        END DO
        !
    ENDIF
    !
    ! Initialize additional index values (used to access load and dwell arrays)
    !
    ILOAD = 1
    IDWELL = 1
    !
    ! Build the input arrays for all loads (including possible dwells)
    !
    DO I = 1, NLOADS
        !
        ! Initial error handling
        !
        IF (TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(ILOAD,5) .EQ. -1) THEN
            !
            CALL PAR_QUIT('Error  :     > Number of load steps does not match&
                & defined number of targets in *.config file.')
            !
        ENDIF
        !
        ! Check if the current index value is a dwell episode
        !
        IF (TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES .GT. 0) THEN
            !
            IF (I .EQ. TRIAXCLR_OPTIONS%DWELL_EPISODE(IDWELL,1)) THEN
                !
                ! Assign the load values from the previous step
                !
                END_LOAD(I, 1) = END_LOAD(I-1, 1)
                END_LOAD(I, 2) = END_LOAD(I-1, 2)
                END_LOAD(I, 3) = END_LOAD(I-1, 3)
                !
                ! Set the ramp rate to zero for dwell
                !
                RAMP_RATE(I)  = 0.0_RK
                CONTROL_DIR(I) = -1 * CONTROL_DIR(I)
                !
                DWELL_TIME(I) = TRIAXCLR_OPTIONS%DWELL_EPISODE(IDWELL,2)
                TARGET_TIME_INCR(I) = TRIAXCLR_OPTIONS%DWELL_EPISODE(IDWELL,3)
                PRINT_FLAG(I) = TRIAXCLR_OPTIONS%DWELL_EPISODE(IDWELL,4)
                !
                ! Increment the index IDWELL to access next available row
                ! Check IF max is reached to avoid out of bounds array access
                !
                IF (IDWELL .EQ. TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES) THEN
                    !
                    IDWELL = IDWELL
                    !           
                ELSE
                    !
                    IDWELL = IDWELL + 1
                    !
                ENDIF
                !
            ELSE ! If not a dwell episode than process a typical load step
                !
                ! Assign the load values for this step
                !
                END_LOAD(I, 1) = TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(ILOAD, 1)
                END_LOAD(I, 2) = TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(ILOAD, 2)
                END_LOAD(I, 3) = TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(ILOAD, 3)
                !
                RAMP_RATE(I)  = STEP_RAMP_RATE(ILOAD)
                TARGET_TIME_INCR(I) = TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(ILOAD, 4)
                PRINT_FLAG(I) = TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(ILOAD, 5)
                !
                ! Increment the index ILOAD to access next available row
                ! Check if max is reached to avoid out of bounds array access
                !
                IF (ILOAD .EQ. TRIAXCLR_OPTIONS%NUMBER_OF_CLR_LOAD_STEPS) THEN
                    !
                    ILOAD = ILOAD
                    !           
                ELSE
                    !
                    ILOAD = ILOAD + 1
                    !
                ENDIF
                !
            ENDIF
            !
        ENDIF
        !
        IF (TRIAXCLR_OPTIONS%NUMBER_OF_DWELL_EPISODES .EQ. 0) THEN
            !
            END_LOAD(I,1) = TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(ILOAD,1)
            END_LOAD(I,2) = TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(ILOAD,2)
            END_LOAD(I,3) = TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(ILOAD,3)
            !
            RAMP_RATE(I)  = STEP_RAMP_RATE(ILOAD)
            TARGET_TIME_INCR(I) = TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(ILOAD,4)
            PRINT_FLAG(I) = TRIAXCLR_OPTIONS%TARGET_CLR_LOAD(ILOAD,5)
            !
            ! Increment the index ILOAD to access next available row
            ! Check if max is reached to avoid out of bounds array access
            !
            IF (ILOAD .EQ. TRIAXCLR_OPTIONS%NUMBER_OF_CLR_LOAD_STEPS) THEN
                !
                ILOAD = ILOAD
                !           
            ELSE
                !
                ILOAD = ILOAD + 1
                !
            ENDIF
            !
        ENDIF
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE READ_CTRL_DATA_CLR
    !
    !===========================================================================
    !
    SUBROUTINE READ_TRIAXCLR_RESTART(ISTEP, CURR_LOAD, PREV_LOAD, &
        & FIRST_INCR_IN_STEP, INCR, TIME, SURF_LOAD_ARRAY, AREA, AREA0, &
        & LENGTH, LENGTH0, CURR_VEL, PREV_ACTION, CURR_ACTION, &
        & INITIAL_LOAD_DWELL_VEL, INITIAL_UNLOAD_DWELL_VEL)
    !
    ! Read TriaxCLR restart information.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    LOGICAL, INTENT(OUT)  :: FIRST_INCR_IN_STEP
    INTEGER, INTENT(OUT)  :: ISTEP
    INTEGER, INTENT(OUT)  :: INCR
    INTEGER, INTENT(OUT)  :: PREV_ACTION, CURR_ACTION
    REAL(RK), INTENT(OUT) :: CURR_LOAD(3)
    REAL(RK), INTENT(OUT) :: PREV_LOAD(3)
    REAL(RK), INTENT(OUT) :: TIME
    REAL(RK), INTENT(OUT) :: SURF_LOAD_ARRAY(NSURFACES,3)
    REAL(RK), INTENT(OUT) :: AREA(NSURFACES)
    REAL(RK), INTENT(OUT) :: AREA0(NSURFACES)
    REAL(RK), INTENT(OUT) :: LENGTH(3), LENGTH0(3)
    REAL(RK), INTENT(OUT) :: CURR_VEL(3)
    REAL(RK), INTENT(OUT) :: INITIAL_LOAD_DWELL_VEL(3)
    REAL(RK), INTENT(OUT) :: INITIAL_UNLOAD_DWELL_VEL(3)
    !
    ! Locals:
    !
    INTEGER :: MYUNIT
    INTEGER :: ISURF
    !
    !---------------------------------------------------------------------------
    !
    MYUNIT = NEWUNITNUMBER()
    OPEN(UNIT = MYUNIT, FILE = TRIM(OPTIONS%RSCTRL_IN), &
        & FORM = 'UNFORMATTED', ACTION = 'READ')
    !
    READ(MYUNIT) ISTEP
    READ(MYUNIT) CURR_LOAD
    READ(MYUNIT) PREV_LOAD
    READ(MYUNIT) FIRST_INCR_IN_STEP
    READ(MYUNIT) INCR
    READ(MYUNIT) TIME
    !
    DO ISURF = 1,NSURFACES
        !
        READ(MYUNIT) SURF_LOAD_ARRAY(ISURF, :)
        !
    ENDDO
    !
    READ(MYUNIT) AREA
    READ(MYUNIT) AREA0
    READ(MYUNIT) LENGTH
    READ(MYUNIT) LENGTH0
    READ(MYUNIT) CURR_VEL
    READ(MYUNIT) PREV_ACTION
    READ(MYUNIT) CURR_ACTION
    READ(MYUNIT) INITIAL_LOAD_DWELL_VEL
    READ(MYUNIT) INITIAL_UNLOAD_DWELL_VEL
    !
    IF (MYID .EQ. 0) THEN
        !
        WRITE(DFLT_U,'(A)') 'Info   : Reading restart control information...'
        WRITE(DFLT_U,'(A)') 'Info   :   - Restart parameters:'
        WRITE(DFLT_U,'(A, I0)')       'Info   :     > Increment:     ', INCR
        WRITE(DFLT_U,'(A, I0)')       'Info   :     > Current Step:  ', ISTEP
        WRITE(DFLT_U,'(A, E14.6)')    'Info   :     > Current Time:  ', TIME
        WRITE(DFLT_U,'(A, 3(E14.6))') 'Info   :     > Current Load:  ', CURR_LOAD
        WRITE(DFLT_U,'(A, 3(E14.6))') 'Info   :     > Previous Load: ', PREV_LOAD
        !
    ENDIF
    !
    CLOSE(MYUNIT)
    !
    RETURN
    !
    END SUBROUTINE READ_TRIAXCLR_RESTART
    !
    !===========================================================================
    !
    SUBROUTINE PRINT_HEADERS()
    !
    ! Print output file headers. For debug purposes only - commented out above
    !
    !---------------------------------------------------------------------------
    !
    ! Locals:
    ! IOSTATUS: Returns value about file command success
    !
    INTEGER :: IOSTATUS
    !
    !---------------------------------------------------------------------------
    !
    IF (MYID .EQ. 0) THEN
        !
        ! Open post.bcs_iter_1 file
        !
        OPEN(OUNITS(BCS_ITER_1_U), FILE = 'post.bcs_iter_1', IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > IO Failure to open &
                &post.bcs_iter_1 file.')
            !
        ENDIF
        !
        ! Write the header to the post.bcs_iter_1 file
        !
        WRITE(OUNITS(BCS_ITER_1_U), '(a)') ' %   step    incr bc_iter_1 &
            &vel_x          vel_y          vel_z          &
            &load_x         load_y         load_z         &
            &target_load_x  target_load_y  target_load_z  &
            &dtime'
        !
        ! Open post.bcs_iter_2 file
        !
        OPEN(OUNITS(BCS_ITER_2_U), FILE = 'post.bcs_iter_2', IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > IO Failure to open &
                &post.BCS_iter_2 file.')
            !
        ENDIF
        !
        ! Write the header to the post.BCS_iter_2 file
        !
        WRITE(OUNITS(BCS_ITER_2_U),'(a)') ' %   step    incr bc_iter_2 &
            &vel_x          vel_y          vel_z          &
            &load_x         load_y         load_z         &
            &target_load_x  target_load_y  target_load_z  &
            &dtime         pert_mag'
        !
        ! Open post.bcs_iter_log file
        !
        OPEN(OUNITS(BCS_ITER_LOG_U), FILE = 'post.bcs_iter_log', &
            & IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > IO Failure to open &
                &post.bcs_iter_log file.')
            !
        ENDIF
        !
        ! Write the header to the post.bcs_iter_log file
        !
        WRITE(OUNITS(BCS_ITER_LOG_U),'(a)') '%       step        incr   bc_iter_1   bc_iter_2'
        !
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE PRINT_HEADERS
    !
END MODULE DRIVER_TRIAXCLR_MOD
