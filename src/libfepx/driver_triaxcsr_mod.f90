! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE DRIVER_TRIAXCSR_MOD
!
! Driver for triaxial load control with constant strain rate.
!
! Contains subroutines:
! DRIVER_TRIAX_CSR: Primary driver for triaxial CSR control simulations.
! PROCESS_CTRL_DATA_CSR: Read data for control and modify init. velocity field.
! DRIVER_TRIAXCSR_RESTART: Read triaxial CSR restart information.
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, RK=>REAL_KIND
USE LIBF95, ONLY: NEWUNITNUMBER
!
! From libfepx:
!
USE DIMENSIONS_MOD, ONLY: DIMS1, TVEC, TVEC1, MAXSLIP1, NGRAIN1, &
    & ANISOTROPIC_EVPS, VTINY
USE DRIVER_UTILITIES_MOD, ONLY: VELOCITY_ITERATION, CALC_MESH_DIM, EST_AVG_MOD,&
    & UPDATE_STATE_EVPS, CALC_STRESS_STRAIN
USE KINEMATICS_MOD, ONLY: PLASTICVELGRADSYMSKW, CALC_TOTAL_WORK,&
    & CALC_PLASTIC_WORK
USE MATRIX_OPERATIONS_MOD, ONLY: CALC_ELVOL, SOLVE_LIN_SYS_3, &
    & STRAIN_EQUIV_3X3, VEC6_MAT_SYMM
USE MICROSTRUCTURE_MOD, ONLY: NUMPHASES, GAMMADOT
USE READ_INPUT_MOD, ONLY: KDIM1, EL_SUB1, EL_SUP1, DOF_SUB1, DOF_SUP1, COORDS,&
    & OPTIONS, BCS_OPTIONS, UNIAXIAL_OPTIONS, TRIAXCSR_OPTIONS, PRINT_OPTIONS, &
    & NODES, READ_RESTART_FIELD, D_TOT, EL_WORK_N, EL_WORKP_N, EL_WORK_RATE_N, &
    & EL_WORKP_RATE_N, EL_WORK, EL_WORKP
USE SURFACE_MOD, ONLY: NSURFACES, COMPUTE_AREA
USE UNITS_MOD, ONLY: DFLT_U, FORCE_U1, FORCE_U2, FORCE_U3, FORCE_U4, FORCE_U5,&
    & FORCE_U6, CONV_U, OUNITS, BCS_ITER_1_U, REPORT_U
USE WRITE_OUTPUT_MOD, ONLY: PRINT_STEP, WRITE_REPORT_FILE_COMPLETE_STEPS, &
    & WRITE_FORCE_FILE_HEADERS, WRITE_FORCE_FILE_DATA, WRITE_CONV_FILE_HEADERS,&
    & WRITE_TRIAXCSR_RESTART, WRITE_RESTART_FIELD
!
! From libparallel:
!
USE GATHER_SCATTER_MOD, ONLY: TRACE, PART_GATHER
USE PARALLEL_MOD, ONLY: MYID, PAR_MESSAGE, PAR_QUIT
!
IMPLICIT NONE
!
! Private
!
INTEGER, PRIVATE :: NSTEPS
REAL(RK), ALLOCATABLE, PRIVATE :: TARGET_LOAD(:,:), TARGET_SIGN(:), &
    & VEL_FACTOR(:), TIME_STEP(:), TIME_STEP_MIN(:)
LOGICAL, ALLOCATABLE, PRIVATE :: IS_SIGNED(:)
INTEGER, ALLOCATABLE, PRIVATE :: PRINT_FLAG(:)
!
CONTAINS
    !
    SUBROUTINE DRIVER_TRIAX_CSR(BCS, VELOCITY, PFORCE, DTRACE, C0_ANGS, &
        & CRSS_N, RSTAR_N, KEINV, WTS, AUTO_TIME, GAMMADOT, CLOCK_START)
    !
    ! Primary driver for triaxial load control with constant strain rate.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! BCS: Global D.O.F. array indicating applied velocity boundary conditions.
    ! VELOCITY: Global D.O.F. array storing nodal velocities.
    ! PFORCE: Global D.O.F. array with tractions (not implemented?).
    ! DTRACE: Gather/scatter trace for degrees of freedom.
    ! C0_ANGS: Initial orientations, rotation matrix.
    ! CRSS_N: Resolved shear stresses on each slip system.
    ! RSTAR_N: Current orientations.
    ! KEINV: Inverse of the single crystal elasticity tensor.
    ! WTS: Elemental weights (??).
    ! AUTO_TIME: (??).
    ! CLOCK_START: Timer start time
    !
    LOGICAL :: BCS(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: VELOCITY(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: PFORCE  (DOF_SUB1:DOF_SUP1)
    TYPE(TRACE) :: DTRACE
    REAL(RK) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: KEINV(0:TVEC1,1:NUMPHASES)
    REAL(RK) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER  :: AUTO_TIME
    REAL(RK) :: GAMMADOT(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: CLOCK_START
    !
    ! Locals:
    ! Need to be defined - JC
    !
    LOGICAL  :: CONVERGED_SOLUTION
    LOGICAL  :: STEP_COMPLETE
    LOGICAL  :: IS_NECKING
    LOGICAL  :: IS_LIMIT_TRIPPED
    INTEGER  :: INCR
    INTEGER  :: ITERNL
    INTEGER  :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DTIME, DTIME_OLD, TIME
    REAL(RK) :: D(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: D_VEC(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TRIAL_VELOCITY(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: VGRAD(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: EVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: RSTAR (0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TRSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ELAS_TOT6(0:TVEC , 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ELAS_TOT_3X3(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: D_TOT(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TOTSTRAIN(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: S_AVG_3X3(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: E_ELAS_KK_BAR(EL_SUB1:EL_SUP1)
    REAL(RK) :: E_ELAS_KK(EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG_KK(EL_SUB1:EL_SUP1)
    REAL(RK) :: STIF(TVEC, TVEC, EL_SUB1:EL_SUP1)
    REAL(RK) :: FE(TVEC, EL_SUB1:EL_SUP1)
    INTEGER :: M_EL, I, II
    INTEGER :: ISTEP
    REAL(RK) :: SURF_LOAD_ARRAY(NSURFACES,3)
    REAL(RK) :: AREA0(NSURFACES), AREA(NSURFACES)
    REAL(RK) :: LENGTH(3), LENGTH0(3)
    REAL(RK) :: GAMMA(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DP_HAT(0:TVEC1,EL_SUB1:EL_SUP1)
    REAL(RK) :: PLSTRAIN(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK) :: WP_HAT(0:DIMS1,EL_SUB1:EL_SUP1)
    REAL(RK) :: DEFF(EL_SUB1:EL_SUP1)
    REAL(RK) :: DPEFF(EL_SUB1:EL_SUP1)
    REAL(RK) :: EQSTRAIN(EL_SUB1:EL_SUP1)
    REAL(RK) :: EQELSTRAIN(EL_SUB1:EL_SUP1)
    REAL(RK) :: EQPLSTRAIN (EL_SUB1:EL_SUP1)
    !
    ! OPTIONAL PARAMETERS.
    INTEGER  :: MAX_BC_ITER
    INTEGER  :: CONTROL_DIR
    REAL(RK) :: INITIAL_VEL
    REAL(RK) :: MIN_PERT_FRAC
    REAL(RK) :: LOAD_TOL_ABS
    REAL(RK) :: LOAD_TOL_REL
    REAL(RK) :: MAX_STRAIN
    !
    ! CONTROL VARIABLES.
    INTEGER  :: PDIR, SDIR, TDIR
    INTEGER  :: IND(1:(DOF_SUP1-DOF_SUB1+1)/3, 1:3)
    INTEGER  :: BC_ITER
    REAL(RK) :: EAVG, NUAVG
    REAL(RK) :: CURR_VEL(3)
    REAL(RK) :: CURR_LOAD(3), PREV_LOAD(3)
    REAL(RK) :: IDEAL_LOAD(3), TRIAL_LOAD(3)
    REAL(RK) :: LOAD_STEP
    REAL(RK) :: STEP_FRAC
    REAL(RK) :: S_PERT_MAG, T_PERT_MAG
    REAL(RK) :: S_PERT_SIGN, T_PERT_SIGN
    REAL(RK) :: S_DELTA_VEL, T_DELTA_VEL
    REAL(RK) :: S_LOAD_ERR, T_LOAD_ERR, MAX_LOAD_ERR
    REAL(RK) :: A(3,3)
    REAL(RK) :: B(3)
    REAL(RK) :: COEFFS(3)
    REAL(RK) :: GAGE_LENGTH
    REAL(RK) :: STRAIN_RATE
    REAL(RK) :: MACRO_ENG_STRAIN(3)
    REAL(RK) :: MAX_EQSTRAIN
    REAL(RK) :: CURR_EQSTRAIN
    CHARACTER(LEN=16) :: FIELD, TIME_STRING, DTIME_STRING
    REAL(RK) :: DKK_STEP
    !
    REAL(RK), ALLOCATABLE :: ECOORDS(:, :)
    REAL(RK), ALLOCATABLE :: ELVOL(:)
    REAL(RK), ALLOCATABLE :: ELVOL_0(:)
    !
    REAL(RK) :: CLOCK_END
    !
    !---------------------------------------------------------------------------
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    ! Extract input/optional parameters
    !
    MAX_BC_ITER   = TRIAXCSR_OPTIONS%MAX_BC_ITER
    MAX_STRAIN = TRIAXCSR_OPTIONS%MAX_STRAIN
    MAX_EQSTRAIN = TRIAXCSR_OPTIONS%MAX_EQSTRAIN
    STRAIN_RATE = BCS_OPTIONS%STRAIN_RATE
    CONTROL_DIR = BCS_OPTIONS%LOADING_DIRECTION + 1
    MIN_PERT_FRAC = TRIAXCSR_OPTIONS%MIN_PERT_FRAC
    LOAD_TOL_ABS  = TRIAXCSR_OPTIONS%LOAD_TOL_ABS
    LOAD_TOL_REL  = TRIAXCSR_OPTIONS%LOAD_TOL_REL
    !
    ! Designate primary (control), secondary, and tertiary directions
    !
    SELECT CASE (CONTROL_DIR)
        !
        CASE (1)
            !
            PDIR = 1
            SDIR = 2
            TDIR = 3
            !
        CASE (2)
            !
            PDIR = 2
            SDIR = 3
            TDIR = 1
            !
        CASE (3)
            !
            PDIR = 3
            SDIR = 1
            TDIR = 2
            !
        CASE DEFAULT
            !
            CALL PAR_QUIT('Error  :     > Invalid control direction provided.')
            !
    END SELECT
    !
    ! Store indices for x, y, and z degrees of freedom
    !
    DO I=1, ((DOF_SUP1-DOF_SUB1+1)/3)
        !
        IND(I,1)= DOF_SUB1+3*(I-1)
        IND(I,2)= DOF_SUB1+3*(I-1)+1
        IND(I,3)= DOF_SUB1+3*(I-1)+2
        !
    ENDDO
    !
    ! Convert input strain rate to initial velocity
    !
    CALL CALC_MESH_DIM(LENGTH, IND(:,1), IND(:,2), IND(:,3))
    !
    GAGE_LENGTH = LENGTH(CONTROL_DIR)
    INITIAL_VEL = STRAIN_RATE * GAGE_LENGTH
    !
    ! Initialization
    !
    IS_NECKING = .FALSE.
    IS_LIMIT_TRIPPED = .FALSE.
    !
    IF (OPTIONS%RESTART) THEN
        !
        CALL READ_RESTART_FIELD(VELOCITY, C0_ANGS, C_ANGS, &
            & RSTAR, RSTAR_N, WTS, CRSS, CRSS_N, &
            & E_ELAS_KK_BAR, SIG_VEC_N, EQSTRAIN, EQPLSTRAIN, GAMMA, &
            & EL_WORK_N, EL_WORKP_N, EL_WORK_RATE_N, EL_WORKP_RATE_N, &
            & PLSTRAIN, TOTSTRAIN)
        !
        CALL READ_TRIAXCSR_RESTART(ISTEP, CURR_LOAD, PREV_LOAD, STEP_COMPLETE, &
            & DTIME, INCR, TIME, SURF_LOAD_ARRAY, &
            & AREA, AREA0, LENGTH, LENGTH0, CURR_VEL, S_PERT_MAG, T_PERT_MAG, &
            & VELOCITY)
        !
    ELSE
        !
        ! Initialize areas and load arrays
        !
        SURF_LOAD_ARRAY = 0.0D0
        CURR_LOAD = 0.0D0
        AREA = 0.0D0
        AREA0 = 0.0D0
        !
        ! Compute initial area (area0)
        !
        CALL COMPUTE_AREA(COORDS, AREA0)
        !
        !IF (MYID .EQ. 0) THEN
            !
        !    WRITE(DFLT_U,'(A25,6F12.6)') 'initial areas (area0):   ', AREA0
            !
        !END IF
        !
        ! Compute initial mesh dimensions
        !
        CALL CALC_MESH_DIM(LENGTH, IND(:,1), IND(:,2), IND(:,3))
        LENGTH0 = LENGTH
        !
        ! Initialize state
        !
        E_ELAS_KK_BAR = 0.0D0
        SIG_VEC_N = 0.0D0
        C_ANGS = C0_ANGS
        !
        ! Initialize elvol arrays (and associated) iff it needs to be printed
        !
        IF (PRINT_OPTIONS%PRINT_ELVOL) THEN
            !
            ALLOCATE(ELVOL(EL_SUB1:EL_SUP1))
            ALLOCATE(ELVOL_0(EL_SUB1:EL_SUP1))
            ALLOCATE(ECOORDS(0:KDIM1, EL_SUB1:EL_SUP1))
            !
            ELVOL = 0.0D0
            ELVOL_0 = 0.0D0
            CALL PART_GATHER(ECOORDS, COORDS, NODES, DTRACE)
            CALL CALC_ELVOL(ELVOL_0, ECOORDS)
            !
        END IF
        !
        ! Initialize integrated quantities
        !
        GAMMA = 0.0D0
        EQPLSTRAIN = 0.0D0
        PLSTRAIN = 0.0D0
        D_TOT = 0.0D0
        TOTSTRAIN = 0.0D0
        EQSTRAIN = 0.0D0
        !
        ! Initialize deformation control
        !
        INCR = 0
        TIME = 0.0D0
        ISTEP = 1
        INCR = 0
        STEP_COMPLETE = .FALSE.
        DTIME = TIME_STEP(1)
        S_PERT_MAG = MIN_PERT_FRAC*ABS(INITIAL_VEL)
        T_PERT_MAG = MIN_PERT_FRAC*ABS(INITIAL_VEL)
        CURR_EQSTRAIN = 0.0D0
        !
        ! Estimate bulk elastic moduli
        ! Assume uniform texture and equal element volumes
        !
        CALL EST_AVG_MOD(EAVG, NUAVG)
        !
        ! Initialize velocity field using elasticity
        !
        CURR_VEL(1) = (TARGET_LOAD(1,1)/AREA0(4)-NUAVG*TARGET_LOAD(1,2)&
            & /AREA0(6)-NUAVG*TARGET_LOAD(1,3)/AREA0(2))*LENGTH0(1)/EAVG
        !
        CURR_VEL(2) = (-NUAVG*TARGET_LOAD(1,1)/AREA0(4)+TARGET_LOAD(1,2)&
            & /AREA0(6)-NUAVG*TARGET_LOAD(1,3)/AREA0(2))*LENGTH0(2)/EAVG
        !
        CURR_VEL(3) = (-NUAVG*TARGET_LOAD(1,1)/AREA0(4)-NUAVG*TARGET_LOAD(1,2)&
            & /AREA0(6)+TARGET_LOAD(1,3)/AREA0(2))*LENGTH0(3)/EAVG
        !
        CURR_VEL = CURR_VEL*ABS(INITIAL_VEL/CURR_VEL(PDIR))
        !
        VELOCITY(IND(:,1)) = CURR_VEL(1)*COORDS(IND(:,1))
        VELOCITY(IND(:,2)) = CURR_VEL(2)*COORDS(IND(:,2))
        VELOCITY(IND(:,3)) = CURR_VEL(3)*COORDS(IND(:,3))
        !
        ! Print initial values
        !
        IF (PRINT_OPTIONS%PRINT_ELVOL) THEN
            !
            CALL PRINT_STEP(0, 0, COORDS, VELOCITY, C0_ANGS, CRSS_N, ELVOL_0)
            !
            DEALLOCATE(ELVOL_0)
            !
        ELSE
            !
            CALL PRINT_STEP(0, 0, COORDS, VELOCITY, C0_ANGS, CRSS_N)
            !
        END IF
        !
    ENDIF
    !
    IF (MYID .EQ. 0) THEN
        !
        ! Write the header to file post.force1-6
        !
        IF (PRINT_OPTIONS%PRINT_FORCES) THEN
            !
            CALL WRITE_FORCE_FILE_HEADERS(2)
            !
            ! If virgin sample, write 0th step
            !
            IF (OPTIONS%RESTART .EQV. .FALSE.) THEN
                !
                CALL WRITE_FORCE_FILE_DATA(2, 0, 0, SURF_LOAD_ARRAY, AREA0, &
                    & 0.0D0, LENGTH0)
                !
            END IF
            !
        END IF
        !
        IF (PRINT_OPTIONS%PRINT_CONV) THEN
            !
            CALL WRITE_CONV_FILE_HEADERS()
            !
        END IF
        !
        ! Open post.bcs_iter file - for debug only
        ! OPEN (OUNITS(BCS_ITER_1_U), FILE='post.bcs_iter_1',IOSTAT=IOSTATUS)
        !
        ! IF (IOSTATUS /= 0) THEN
            !
            ! CALL PAR_QUIT('Error  :     > IO Failure to open &
            ! &post.bcs_iter file.')
            !
        ! END IF
        !
        ! Write the header to file post.bcs_iter
        !
        ! WRITE(OUNITS(BCS_ITER_1_U),'(a)')  '%   step    incr bc_iter    &
        ! &vel_x          vel_y          vel_z          &
        ! &load_x         load_y         load_z         &
        ! &target_load_x  target_load_y  target_load_z  &
        ! &dtime          s_pert_mag     t_pert_mag'
        !
    END IF
    !
    ! IF (MYID .EQ. 0) WRITE(DFLT_U, '(A)') 'Info   :   - Starting T-stepping &
    ! &loop for anisotropic viscoplastic solution'
    ! IF (MYID .EQ. 0) WRITE(DFLT_U,'(A)') 'Info   : Running step 1...'
    !
    ! Time stepping loop
    !
    TIME_STEPPING : DO
        !
        ! Iterate to find new configuration and material state.
        !
        INCR = INCR + 1
        !
        IF (STEP_COMPLETE .AND.(MYID .EQ. 0))  THEN
            !
            WRITE(DFLT_U,'(A,I0,A)') 'Info   : Running step ', ISTEP, '...'
            !
        END IF
        !
        ! Initialize shear rates
        !
        GAMMADOT = 0.0D0
        !
        CONVERGED_SOLUTION = .TRUE.
        !
        ! Calculate new time increment
        ! (based on time_increment from load_disp_control_mod)
        !
        IF (INCR .GT. 1) THEN
            !
            DTIME_OLD = DTIME
            DTIME = TIME_STEP(ISTEP)
            !
            ! Check if we are close to target and need to adjust the time step
            !
            TRIAL_LOAD = CURR_LOAD + (CURR_LOAD-PREV_LOAD)*DTIME/DTIME_OLD
            !
            IF (TARGET_SIGN(ISTEP)*(TRIAL_LOAD(PDIR)-TARGET_LOAD(ISTEP,PDIR)) &
                & .GE. 0.0) THEN
                !
                IF (ABS(TRIAL_LOAD(PDIR) - CURR_LOAD(PDIR)) .GT. VTINY) THEN
                    !
                    DTIME = (TARGET_LOAD(ISTEP,PDIR) - CURR_LOAD(PDIR)) / &
                        & (TRIAL_LOAD(PDIR) - CURR_LOAD(PDIR)) * DTIME * &
                        & OPTIONS%DTIME_FACTOR
                    !
                    IF (DTIME .GT. TIME_STEP(ISTEP)) THEN
                        !
                        DTIME = TIME_STEP(ISTEP)
                        !
                    END IF
                    !
                    IF (DTIME .LT. TIME_STEP_MIN(ISTEP)) THEN
                        !
                        DTIME = TIME_STEP_MIN(ISTEP)
                        !
                    END IF
                    !
                END IF
                !
            END IF
            !
        END IF
        !
        ! Update the total time
        !
        TIME = TIME + DTIME
        !
        IF (MYID .EQ. 0) THEN
            !
            WRITE(FIELD, '(F0.4)') TIME
            IF (FIELD(1:1) == '.') FIELD = '0' // FIELD
            TIME_STRING = FIELD
            WRITE(FIELD, '(F0.4)') DTIME
            IF (FIELD(1:1) == '.') FIELD = '0' // FIELD
            DTIME_STRING = FIELD
            !
            WRITE(DFLT_U,'(A,I0,A,A,A,A,A)') 'Info   :   - &
                &Increment ', INCR, ': t = ', TRIM(TIME_STRING), '&
                & secs, dt = ', TRIM(DTIME_STRING), ' secs'
            !
        END IF
        !
        ! Store previous load
        !
        PREV_LOAD = CURR_LOAD
        !
        ! Update initial guess of velocity field for load reversal
        !
        IF (STEP_COMPLETE) THEN
            !
            VELOCITY = VEL_FACTOR(ISTEP)*VELOCITY
            CURR_VEL = VEL_FACTOR(ISTEP)*CURR_VEL
            !
        END IF
        !
        ! Surface velocity iteration loop
        !
        BC_ITER = 0
        !
        BC_ITERATION : DO
            !
            ! Check number of boundary condition iterations
            !
            IF (BC_ITER .GE. MAX_BC_ITER) THEN
                !
                CALL PAR_QUIT('Error  :     > Maximum number of &
                    &boundary conditions reached.')
                !
            END IF
            !
            BC_ITER = BC_ITER + 1
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U, '(A,I0)') 'Info   :     > TRIAL_BC: Iteration ', &
                    & BC_ITER
                !
            END IF
            !
            ! Begin velocity iteration
            !
            CALL VELOCITY_ITERATION(BCS, PFORCE, VELOCITY, EVEL, CURR_LOAD, &
                & DTRACE, C0_ANGS, C_ANGS, SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, &
                & RSTAR_N, RSTAR, TRSTAR, KEINV, E_ELAS_KK_BAR, E_ELAS_KK, &
                & SIG_KK, JITER_STATE, WTS, DTIME, INCR, E_BAR_VEC,&
                & CONVERGED_SOLUTION, AUTO_TIME, ITERNL)
            !
            ! Check load tolerances
            !
            IF (ISTEP .EQ. 1) THEN
                !
                STEP_FRAC = CURR_LOAD(PDIR)/TARGET_LOAD(1,PDIR)
                IDEAL_LOAD = STEP_FRAC*(TARGET_LOAD(1,:))
                !
            ELSE
                !
                STEP_FRAC = (CURR_LOAD(PDIR)-TARGET_LOAD(ISTEP-1,PDIR))/&
                    & (TARGET_LOAD(ISTEP,PDIR)-TARGET_LOAD(ISTEP-1,PDIR))
                IDEAL_LOAD = STEP_FRAC*TARGET_LOAD(ISTEP,:) &
                    & + (1-STEP_FRAC)*TARGET_LOAD(ISTEP-1,:)
                !
            END IF
            !
            ! Write output to console and file
            !
            IF (MYID .EQ. 0) THEN
                !
                ! For debug only
                ! WRITE(OUNITS(BCS_ITER_1_U),'(3(i8),12(e15.5))') ISTEP, INCR, &
                ! & BC_ITER, CURR_VEL, CURR_LOAD, IDEAL_LOAD, DTIME, &
                ! & S_PERT_MAG, T_PERT_MAG
                !
                IF (MINVAL(ABS(CURR_VEL)) .LE. 0.0010) THEN
                    !
                    WRITE(DFLT_U,'(A,3(E12.2))') 'Info   :       . &
                        &Velocity:    ', CURR_VEL
                    !
                ELSE
                    !
                    WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       . &
                        &Velocity:    ', CURR_VEL
                    !
                END IF
                !
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       . &
                    &Load:        ', CURR_LOAD
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       . &
                    &Target load: ', IDEAL_LOAD
                !
            END IF
            !
            ! Check if loads are within tolerance
            !
            S_LOAD_ERR = ABS(CURR_LOAD(SDIR)-IDEAL_LOAD(SDIR))
            T_LOAD_ERR = ABS(CURR_LOAD(TDIR)-IDEAL_LOAD(TDIR))
            MAX_LOAD_ERR = MAX(LOAD_TOL_ABS,LOAD_TOL_REL*ABS(CURR_LOAD(PDIR)))
            !
            IF ((S_LOAD_ERR .LE. MAX_LOAD_ERR) .AND. &
                & (T_LOAD_ERR .LE. MAX_LOAD_ERR)) THEN
                !
                ! Load within tolerance
                !
                EXIT
                !
            END IF
            !
            ! Perturb secondary direction velocity
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A)') 'Info   :       . Perturbing surface &
                    &velocity in secondary direction'
                !
            END IF
            !
            IF (CURR_LOAD(SDIR) .LE. IDEAL_LOAD(SDIR)) THEN
                !
                S_PERT_SIGN = 1.0
            ELSE
                !
                S_PERT_SIGN = -1.0
                !
            END IF
            !
            TRIAL_VELOCITY = VELOCITY
            TRIAL_VELOCITY(IND(:,SDIR)) = VELOCITY(IND(:,SDIR)) &
                & + S_PERT_SIGN*S_PERT_MAG*COORDS(IND(:,SDIR))/LENGTH(SDIR)
            !
            CALL VELOCITY_ITERATION(BCS, PFORCE, TRIAL_VELOCITY, EVEL, &
                & TRIAL_LOAD, DTRACE, C0_ANGS, C_ANGS, SIG_VEC_N, SIG_VEC, &
                & CRSS_N, CRSS, RSTAR_N, RSTAR, TRSTAR, KEINV, E_ELAS_KK_BAR, &
                & E_ELAS_KK, SIG_KK, JITER_STATE, WTS, DTIME, INCR, E_BAR_VEC, &
                & CONVERGED_SOLUTION, AUTO_TIME, ITERNL)
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :     > Trial load: ', &
                    & TRIAL_LOAD
                !
            END IF
            !
            ! Compute matrix coefficients
            !
            A(1,1) = TRIAL_LOAD(1)-CURR_LOAD(1)
            A(2,1) = TRIAL_LOAD(2)-CURR_LOAD(2)
            A(3,1) = TRIAL_LOAD(3)-CURR_LOAD(3)
            !
            ! Perturb tertiary direction velocity
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A)') 'Info   :       . Perturbing surface &
                    &velocity in tertiary direction'
                !
            END IF
            !
            IF (CURR_LOAD(TDIR) .LE. IDEAL_LOAD(TDIR)) THEN
                !
                T_PERT_SIGN = 1.0
                !
            ELSE
                !
                T_PERT_SIGN = -1.0
                !
            END IF
            !
            TRIAL_VELOCITY = VELOCITY
            TRIAL_VELOCITY(IND(:,TDIR)) = VELOCITY(IND(:,TDIR)) &
                & + T_PERT_SIGN*T_PERT_MAG*COORDS(IND(:,TDIR))/LENGTH(TDIR)
            !
            CALL VELOCITY_ITERATION(BCS, PFORCE, TRIAL_VELOCITY, EVEL, &
                & TRIAL_LOAD, DTRACE, C0_ANGS, C_ANGS, SIG_VEC_N, SIG_VEC, &
                & CRSS_N, CRSS, RSTAR_N, RSTAR, TRSTAR, KEINV, E_ELAS_KK_BAR, &
                & E_ELAS_KK, SIG_KK, JITER_STATE, WTS, DTIME, INCR, E_BAR_VEC, &
                & CONVERGED_SOLUTION, AUTO_TIME, ITERNL)
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :     > Trial load: ', &
                    & TRIAL_LOAD
                !
            END IF
            !
            ! Compute matrix coefficients
            !
            A(1,2) = TRIAL_LOAD(1)-CURR_LOAD(1)
            A(2,2) = TRIAL_LOAD(2)-CURR_LOAD(2)
            A(3,2) = TRIAL_LOAD(3)-CURR_LOAD(3)
            !
            ! Compute change to surface velocities
            !
            IF (ISTEP .EQ. 1) THEN
                !
                LOAD_STEP = TARGET_LOAD(ISTEP,PDIR)
                !
                A(1,3) = -TARGET_LOAD(ISTEP,1)/LOAD_STEP
                A(2,3) = -TARGET_LOAD(ISTEP,2)/LOAD_STEP
                A(3,3) = -TARGET_LOAD(ISTEP,3)/LOAD_STEP
                !
                B(1) = -CURR_LOAD(1)
                B(2) = -CURR_LOAD(2)
                B(3) = -CURR_LOAD(3)
                !
            ELSE
                !
                LOAD_STEP = TARGET_LOAD(ISTEP,PDIR)-TARGET_LOAD(ISTEP-1,PDIR)
                !
                A(1,3) = (TARGET_LOAD(ISTEP-1,1)-TARGET_LOAD(ISTEP,1))/LOAD_STEP
                A(2,3) = (TARGET_LOAD(ISTEP-1,2)-TARGET_LOAD(ISTEP,2))/LOAD_STEP
                A(3,3) = (TARGET_LOAD(ISTEP-1,3)-TARGET_LOAD(ISTEP,3))/LOAD_STEP
                !
                B(1) = TARGET_LOAD(ISTEP-1,1)-CURR_LOAD(1)
                B(2) = TARGET_LOAD(ISTEP-1,2)-CURR_LOAD(2)
                B(3) = TARGET_LOAD(ISTEP-1,3)-CURR_LOAD(3)
                !
            END IF
            !
            ! Solve linear system of equations
            !
            CALL SOLVE_LIN_SYS_3(A, B, COEFFS)
            !
            ! Extract changes to velocity field
            !
            S_DELTA_VEL = COEFFS(1)*S_PERT_SIGN*S_PERT_MAG
            T_DELTA_VEL = COEFFS(2)*T_PERT_SIGN*T_PERT_MAG
            !
            ! Apply new surface velocities
            !
            VELOCITY(IND(:,SDIR)) = VELOCITY(IND(:,SDIR)) &
                & + S_DELTA_VEL*COORDS(IND(:,SDIR))/LENGTH(SDIR)
            VELOCITY(IND(:,TDIR)) = VELOCITY(IND(:,TDIR)) &
                & + T_DELTA_VEL*COORDS(IND(:,TDIR))/LENGTH(TDIR)
            !
            CURR_VEL(SDIR) = CURR_VEL(SDIR) + S_DELTA_VEL
            CURR_VEL(TDIR) = CURR_VEL(TDIR) + T_DELTA_VEL
            !
            ! Update perturbation magnitudes
            S_PERT_MAG = MAX(ABS(S_DELTA_VEL), MIN_PERT_FRAC*CURR_VEL(PDIR))
            T_PERT_MAG = MAX(ABS(T_DELTA_VEL), MIN_PERT_FRAC*CURR_VEL(PDIR))
            !
        END DO BC_ITERATION
        !
        ! Update RSTAR and SIG_VEC_N
        !
        RSTAR = TRSTAR
        !
        ! Update state (at center of each element), geometry
        !
        ! Update:
        ! COORDS @(t+dt), using DT and VELOCITY @(t+dt)
        !
        ! Element centroid:
        ! SIG_VEC_N: 5-vec deviatoric Kirchhoff stress @(t)
        ! SIG_VEC: 5-vec deviatoric Kirchhoff stress @(t+dt)
        ! SIG_KK: volumetric Kirchhoff stress @(t+dt)
        ! E_ELAS_KK_BAR: volumetric lattice strain @(t)
        ! E_ELAS_KK: volumetric lattice strain @(t+dt)
        !
        ! These are outputs of the subroutine/function but are useless and can
        ! be removed from the argument list
        ! E_BAR_VEC: OUT
        ! E_ELAS_KK_BAR: OUT
        ! SIG_VEC_N: OUT
        !
        CALL UPDATE_STATE_EVPS(VELOCITY, DTRACE, C0_ANGS, C_ANGS, SIG_VEC_N, &
            & SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, E_BAR_VEC, E_ELAS_KK_BAR, &
            & E_ELAS_KK, SIG_KK, KEINV, DEFF, DTIME, STIF, FE, D, D_VEC, &
            & VGRAD, CONVERGED_SOLUTION, AUTO_TIME)
        !
        ! Calculate stress, strain, load, and area
        !
        CALL CALC_STRESS_STRAIN(S_AVG_3X3, ELAS_TOT6, SURF_LOAD_ARRAY, AREA, &
            & SIG_VEC, SIG_KK, E_ELAS_KK, C_ANGS, KEINV, WTS)
        !
        CURR_LOAD(1) = SURF_LOAD_ARRAY(2,1) ! 4,1 LEGACY FACESET VALUES
        CURR_LOAD(2) = SURF_LOAD_ARRAY(4,2) ! 6,2
        CURR_LOAD(3) = SURF_LOAD_ARRAY(6,3) ! 2,3
        !
        ! Compute mesh dimensions
        !
        CALL CALC_MESH_DIM(LENGTH, IND(:,1), IND(:,2), IND(:,3))
        !
        ! Calculate macroscopic engineering strain (ported from TriaxCLR)
        !
        DO II = 1, 3
            !
            MACRO_ENG_STRAIN(II) = LENGTH(II) / LENGTH0(II) - 1.0D0
            !
        END DO
        !
        ! Calculate the current increment macroscopic eqstrain
        !
        CURR_EQSTRAIN = (2.0D0 / 3.0D0) * &
            & DSQRT((3.0D0 / 2.0D0) * ((MACRO_ENG_STRAIN(1) ** 2.0D0) &
            & + (MACRO_ENG_STRAIN(2) ** 2.0D0) &
            & + (MACRO_ENG_STRAIN(3) ** 2.0D0)))
        !
        ! Print the macroscopic strain values to console for monitoring
        !
        IF (MYID .EQ. 0) THEN
            !
            WRITE(DFLT_U, '(A,F12.6)') 'Info   :       . &
                &Maximum eng. strain: ', MAXVAL(ABS(MACRO_ENG_STRAIN(:)))
            WRITE(DFLT_U, '(A,F12.6)') 'Info   :       . &
                &Current eqv. strain: ', CURR_EQSTRAIN
            !
        END IF
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
        END IF
        !
        ! Calculate DPEFF
        !
        CALL PLASTICVELGRADSYMSKW(DP_HAT,WP_HAT,DPEFF,GAMMADOT(:,0,:),M_EL)
        !
        ! Update EQPLSTRAIN
        !
        EQPLSTRAIN=EQPLSTRAIN+DPEFF*DTIME
        !
        ! Update EQELSTRAIN
        !
        CALL VEC6_MAT_SYMM(ELAS_TOT6, ELAS_TOT_3X3, M_EL, NGRAIN1)
        CALL STRAIN_EQUIV_3X3(ELAS_TOT_3X3, EQELSTRAIN)
        !
        ! Update PLSTRAIN
        !
        PLSTRAIN = PLSTRAIN + DP_HAT * DTIME
        !
        ! Update slip system shear
        !
        GAMMA=GAMMA+GAMMADOT*DTIME
        !
        ! Reconstruct total deformation rate tensor if requested
        !
        IF ((PRINT_OPTIONS%PRINT_WORK) .OR. (PRINT_OPTIONS%PRINT_DEFRATE) .OR. &
            & (PRINT_OPTIONS%PRINT_WORKRATE) .OR. &
            & (PRINT_OPTIONS%PRINT_STRAIN_TOT) .OR. &
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
        ! Update TOTSTRAIN
        !
        TOTSTRAIN = TOTSTRAIN + D_TOT * DTIME
        !
        ! Update EQSTRAIN
        !
        CALL STRAIN_EQUIV_3X3(TOTSTRAIN, EQSTRAIN)
        !
        ! Calculate work, plastic work
        !
        IF ((PRINT_OPTIONS%PRINT_WORK) .OR. (PRINT_OPTIONS%PRINT_WORKRATE) &
            & .OR. (PRINT_OPTIONS%PRINT_RESTART)) THEN
            !
            CALL CALC_TOTAL_WORK(DTIME, D_TOT, S_AVG_3X3, EL_WORK_N, &
                & EL_WORK_RATE_N, EL_WORK)
            !
        END IF
        !
        IF ((PRINT_OPTIONS%PRINT_WORKP) .OR. (PRINT_OPTIONS%PRINT_WORKRATEP) &
            & .OR. (PRINT_OPTIONS%PRINT_RESTART)) THEN
            !
            CALL CALC_PLASTIC_WORK(DTIME, DP_HAT, C_ANGS, S_AVG_3X3, &
                & EL_WORKP_N, EL_WORKP_RATE_N, EL_WORKP)
            !
        END IF
        !
        IF (PRINT_OPTIONS%PRINT_ELVOL) THEN
            !
            CALL PART_GATHER(ECOORDS, COORDS, NODES, DTRACE)
            CALL CALC_ELVOL(ELVOL, ECOORDS)
            !
        END IF
        !
        ! Check if the target load is reached
        !
        STEP_COMPLETE = .FALSE.
        !
        !IF (TARGET_SIGN(ISTEP)*(CURR_LOAD(PDIR)-TARGET_LOAD(ISTEP,PDIR))&
        !    & + OPTIONS%LOAD_TOL .GE. 0.0) THEN
        IF (ABS(CURR_LOAD(PDIR) - TARGET_LOAD(ISTEP, PDIR)) .LE. &
            & TRIAXCSR_OPTIONS%LOAD_TOL_ABS) THEN
            !
            STEP_COMPLETE = .TRUE.
            !
        END IF
        !
        ! Check for necking
        !
        IF (OPTIONS%CHECK_NECKING .AND. ((CURR_LOAD(PDIR) - PREV_LOAD(PDIR))&
            & /CURR_VEL(PDIR) .LT. 0.0)) THEN
            !
            IS_NECKING = .TRUE.
            !
        END IF
        !
        ! Check time and increment count
        !
        IF ((TIME .GE. OPTIONS%MAX_TOTAL_TIME) .OR. &
            & (INCR .GE. OPTIONS%MAX_INCR)) THEN
            !
            IS_LIMIT_TRIPPED = .TRUE.
            !
        ENDIF
        !
        ! Output computed quantities.
        !
        IF (STEP_COMPLETE .OR. IS_NECKING .OR. IS_LIMIT_TRIPPED) THEN
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A,I0)') 'Info   :     > Converged on &
                    &increment ', INCR
                WRITE(DFLT_U,'(A,3(E14.6))') 'Info   :     > Loads on X0: ', &
                    &SURF_LOAD_ARRAY(1,:)
                WRITE(DFLT_U,'(A,3(E14.6))') 'Info   :       Loads on X1: ', &
                    & SURF_LOAD_ARRAY(2,:)
                WRITE(DFLT_U,'(A,3(E14.6))') 'Info   :       Loads on Y0: ', &
                    & SURF_LOAD_ARRAY(3,:)
                WRITE(DFLT_U,'(A,3(E14.6))') 'Info   :       Loads on Y1: ', &
                    & SURF_LOAD_ARRAY(4,:)
                WRITE(DFLT_U,'(A,3(E14.6))') 'Info   :       Loads on Z0: ', &
                    & SURF_LOAD_ARRAY(5,:)
                WRITE(DFLT_U,'(A,3(E14.6))') 'Info   :       Loads on Z1: ', &
                    & SURF_LOAD_ARRAY(6,:)
                !
            ENDIF
            !
            IF (PRINT_FLAG(ISTEP) .EQ. 0) THEN
                !
                CALL PRINT_STEP(ISTEP, ANISOTROPIC_EVPS, COORDS, VELOCITY, &
                    & C_ANGS, CRSS, ELVOL, ELAS_TOT6, S_AVG_3X3, DEFF, &
                    & EQSTRAIN, DPEFF, EQPLSTRAIN, VGRAD, DP_HAT, WP_HAT, &
                    & GAMMA, GAMMADOT, EL_WORK, EL_WORKP, D_TOT, &
                    & EL_WORK_RATE_N, EL_WORKP_RATE_N, EQELSTRAIN, PLSTRAIN, &
                    & TOTSTRAIN)
                !
                IF (PRINT_OPTIONS%PRINT_RESTART) THEN
                    !
                    CALL WRITE_RESTART_FIELD(VELOCITY, C0_ANGS, C_ANGS, &
                        & RSTAR, RSTAR_N, WTS, CRSS, CRSS_N, &
                        & E_ELAS_KK_BAR, SIG_VEC_N, EQSTRAIN, EQPLSTRAIN, &
                        & GAMMA, EL_WORK_N, EL_WORKP_N, EL_WORK_RATE_N, &
                        & EL_WORKP_RATE_N, PLSTRAIN, TOTSTRAIN, ISTEP)
                    !
                    CALL WRITE_TRIAXCSR_RESTART(ISTEP+1, CURR_LOAD, PREV_LOAD, &
                        & STEP_COMPLETE, DTIME, INCR, TIME, SURF_LOAD_ARRAY, &
                        & AREA, AREA0, LENGTH, LENGTH0, CURR_VEL, S_PERT_MAG, &
                        & T_PERT_MAG)
                    !
                END IF
                !
            ENDIF
            !
            ! Check that maximum strain has not been exceeded
            !
            IF (MAXVAL(ABS(MACRO_ENG_STRAIN)) .GE. MAX_STRAIN) THEN
                !
                ! MPK - 10/2021: Don't write that we have completed the final
                !   step. We are currently in the middle of a step.
                CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP - 1, PRINT_FLAG)
                !
                CALL PAR_QUIT('Info   :     > Maximum eng. strain exceeded.')
                !
            ENDIF
            !
            ! Check that maximum eqstrain has not been exceeded
            !
            IF (CURR_EQSTRAIN .GE. MAX_EQSTRAIN) THEN
                !
                ! MPK - 10/2021: Don't write that we have completed the final
                !   step. We are currently in the middle of a step.
                CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP - 1, PRINT_FLAG)
                !
                CALL PAR_QUIT('Info   :     > Maximum eqv. strain exceeded.')
                !
            ENDIF
            !
            IF (IS_NECKING) THEN
                !
                CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP, PRINT_FLAG)
                !
                CALL PAR_QUIT('Error  :     > Specimen is necking.')
                !
            ENDIF
            !
            IF (IS_LIMIT_TRIPPED) THEN
                !
                CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP, PRINT_FLAG)
                !
                CALL PAR_QUIT('Error  :     > Maximum time or maximum &
                    &increments exceeded.')
                !
            ENDIF
            !
            IF (ISTEP .EQ. NSTEPS) THEN
                !
                CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP, PRINT_FLAG)
                !
                ! Finalize clock values and print to console
                !
                IF (MYID .EQ. 0) THEN
                    !
                    CALL CPU_TIME(CLOCK_END)
                    WRITE(DFLT_U, '(A, F10.3, A)') 'Info   : Elapsed time:', &
                        & (CLOCK_END - CLOCK_START), ' secs.'
                    !
                END IF
                !
                CALL PAR_QUIT('Info   : Final step terminated. Simulation&
                    & completed successfully.')
                !
            END IF
            !
            ! Advance step count
            !
            ISTEP = ISTEP + 1
            !
        ENDIF
        !
    END DO TIME_STEPPING
    !
    RETURN
    !
    END SUBROUTINE DRIVER_TRIAX_CSR
    !
    !===========================================================================
    !
    SUBROUTINE PROCESS_CTRL_DATA_CSR(VELOCITY)
    !
    ! Read data for triaxial load control at constant strain rate.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! VELOCITY: Global D.O.F. array storing velocity components.
    !
    REAL(RK), INTENT(INOUT) :: VELOCITY(DOF_SUB1:DOF_SUP1)
    !
    ! Locals:
    ! I/II: Generic loop index values.
    ! PRIMARY_DIR: Loading direction value used to index load arrays.
    ! DIFF_STRESS: Load target difFErential calculated between steps.
    ! USER_STRAIN_RATE: User-defined initial strain rate from BCs input.
    ! STEP_STRAIN_RATE: Array used to store strain rate jump values.
    ! IS_SIGNED: Array to enable load reversal checks to sign VEL_FACTOR.
    !
    INTEGER :: I, II, PRIMARY_DIR
    REAL(RK) :: DIFF_STRESS, USER_STRAIN_RATE
    REAL(RK), ALLOCATABLE :: STEP_STRAIN_RATE(:)
    !
    ! Notes:
    ! User-defined input (strain rates) should ALWAYS be positive. Compression
    ! and similar signed deformation control are handled internally by way of
    ! the sign in the targeted load.
    !
    !---------------------------------------------------------------------------
    !
    ! Initialize definition variables.
    !
    NSTEPS = TRIAXCSR_OPTIONS%NUMBER_OF_CSR_LOAD_STEPS
    PRIMARY_DIR = BCS_OPTIONS%LOADING_DIRECTION + 1
    USER_STRAIN_RATE = BCS_OPTIONS%STRAIN_RATE
    DIFF_STRESS = 0.0D0
    !
    ! Allocate the arrays for the load target sequence.
    !
    ALLOCATE(TARGET_LOAD(NSTEPS,3))
    ALLOCATE(TARGET_SIGN(NSTEPS))
    ALLOCATE(VEL_FACTOR(NSTEPS))
    ALLOCATE(TIME_STEP(NSTEPS))
    ALLOCATE(TIME_STEP_MIN(NSTEPS))
    ALLOCATE(PRINT_FLAG(NSTEPS))
    !
    ALLOCATE(IS_SIGNED(NSTEPS))
    ALLOCATE(STEP_STRAIN_RATE(NSTEPS))
    !
    ! Initialize arrays.
    !
    IS_SIGNED = .FALSE.
    STEP_STRAIN_RATE = BCS_OPTIONS%STRAIN_RATE
    !
    ! Build step strain rate array, specifically handle strain rate jumps.
    !
    IF (UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS .GT. 0) THEN
        !
        II = 1
        !
        DO I = 1, NSTEPS
            !
            IF (UNIAXIAL_OPTIONS%STRAIN_RATE_JUMP(II,1) .EQ. I) THEN
                !
                STEP_STRAIN_RATE(I:NSTEPS) = &
                    & UNIAXIAL_OPTIONS%STRAIN_RATE_JUMP(II,2)
                !
                IF (II.EQ.UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS) THEN
                    !
                    II = II
                    !
                ELSE
                    !
                    II = II + 1
                    !
                END IF
                !
            END IF
            !
        END DO
        !
    END IF
    !
    ! Check if the number of strain jumps defined is above the number of
    ! jumps provided by the user.
    !
    IF (UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_RATE_JUMPS .GT. 0) THEN
        !
        IF (UNIAXIAL_OPTIONS%STRAIN_RATE_JUMP(II,1) .EQ. -1) THEN
            !
            CALL PAR_QUIT('Error  :     > Input `number_of_strain_rate_jumps`&
                & does not match number of defined `strain_rate_jumps`.')
            !
        END IF
        !
    END IF
    !
    ! Process all load steps and build control arrays.
    !
    DO I = 1, NSTEPS
        !
        ! Initial error handling.
        !
        IF (TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,6) .EQ. -1) THEN
            !
            CALL PAR_QUIT('Error  :     > Number of load steps does not match&
                & defined number of targets in *.config file.')
            !
        END IF
        !
        ! Handle the first step alone to avoid out of bounds on arrays.
        !
        IF (I .EQ. 1) THEN
            !
            ! First, assign the fixed values
            !
            TARGET_LOAD(I,1) = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,1)
            TARGET_LOAD(I,2) = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,2)
            TARGET_LOAD(I,3) = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,3)
            TIME_STEP(I) = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,4)
            TIME_STEP_MIN(I) = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,5)
            PRINT_FLAG(I) = INT(TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,6))
            !
            ! Next, determine the velocity factor for the first step.
            !
            VEL_FACTOR(I) = ABS(STEP_STRAIN_RATE(I) / USER_STRAIN_RATE)
            !
            ! Check if the vel factor needs to be signed for compression.
            ! Calculate the differential stress from zero in pdir.
            !
            DIFF_STRESS = (TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,PRIMARY_DIR))
            !
            IF (DIFF_STRESS .GT. 0.0D0) THEN
                !
                IS_SIGNED(I)   = .FALSE.
                TARGET_SIGN(I) = 1
                !VEL_FACTOR(I)  = 1
                !
            ELSE IF (DIFF_STRESS .LT. 0.0D0) THEN
                !
                IS_SIGNED(I) = .TRUE.
                TARGET_SIGN(I) = -1
                !VEL_FACTOR(I)  = 1
                VELOCITY = VELOCITY * ((-1)*(VEL_FACTOR(I)))
                !
            ELSE ! The differential change is zero, thus dt is zero.
                !
                CALL PAR_QUIT('Error  :     > Load differential is zero&
                    & between steps.')
                !
            END IF
            !
        ELSE IF (I .GT. 1) THEN
            !
            ! First, assign the fixed values for I > 1
            !
            TARGET_LOAD(I,1) = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,1)
            TARGET_LOAD(I,2) = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,2)
            TARGET_LOAD(I,3) = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,3)
            TIME_STEP(I) = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,4)
            TIME_STEP_MIN(I) = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,5)
            PRINT_FLAG(I) = INT(TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,6))
            !
            ! Next, determine the velocity factor for the step.
            !
            VEL_FACTOR(I) = STEP_STRAIN_RATE(I) / STEP_STRAIN_RATE(I-1)
            !
            ! Check if the factor needs to be signed for compression.
            ! Calculate the differential stress from previous.
            !
            DIFF_STRESS = (TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,PRIMARY_DIR)) &
                & - (TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I-1,PRIMARY_DIR))
            !
            IF (DIFF_STRESS .GT. 0.0D0) THEN
                !
                IS_SIGNED(I)   = .FALSE.
                TARGET_SIGN(I) = 1
                !
            ELSE IF (DIFF_STRESS .LT. 0.0D0) THEN
                !
                IS_SIGNED(I) = .TRUE.
                TARGET_SIGN(I) = -1
                !
            ELSE ! The differential change is zero, thus dt is zero.
                !
                CALL PAR_QUIT('Error  :     > Load differential is zero&
                    & between steps.')
                !
            END IF
            !
            ! Check if a load reversal occured since the last step.
            !
            IF (IS_SIGNED(I) .NEQV. IS_SIGNED(I-1)) THEN
                !
                VEL_FACTOR(I) = VEL_FACTOR(I) * (-1)
                !
            ELSE IF (IS_SIGNED(I) .EQV. IS_SIGNED(I-1)) THEN
                !
                VEL_FACTOR(I) = VEL_FACTOR(I) * (1)
                !
            ELSE
                !
                CALL PAR_QUIT('Error  :     > Load reversal check failed.')
                !
            END IF
            !
        END IF
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE PROCESS_CTRL_DATA_CSR
    !
    !===========================================================================
    !
    SUBROUTINE READ_TRIAXCSR_RESTART(ISTEP,CURR_LOAD,PREV_LOAD,STEP_COMPLETE, &
        & DTIME,INCR,TIME,SURF_LOAD_ARRAY,AREA,AREA0,LENGTH,LENGTH0,CURR_VEL, &
        & S_PERT_MAG, T_PERT_MAG, VELOCITY)
    !
    ! Read triaxial CSR restart information.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! Needs to be defined - JC
    !
    LOGICAL, INTENT(OUT) :: STEP_COMPLETE
    INTEGER, INTENT(OUT) :: ISTEP
    INTEGER, INTENT(OUT) :: INCR
    REAL(RK), INTENT(OUT) :: CURR_LOAD(3)
    REAL(RK), INTENT(OUT) :: PREV_LOAD(3)
    REAL(RK), INTENT(OUT) :: DTIME
    REAL(RK), INTENT(OUT) :: TIME
    REAL(RK), INTENT(OUT) :: SURF_LOAD_ARRAY(NSURFACES,3)
    REAL(RK), INTENT(OUT) :: AREA(NSURFACES)
    REAL(RK), INTENT(OUT) :: AREA0(NSURFACES)
    REAL(RK), INTENT(OUT) :: LENGTH(3)
    REAL(RK), INTENT(OUT) :: LENGTH0(3)
    REAL(RK), INTENT(OUT) :: CURR_VEL(3)
    REAL(RK), INTENT(OUT) :: S_PERT_MAG
    REAL(RK), INTENT(OUT) :: T_PERT_MAG
    REAL(RK), INTENT(INOUT) :: VELOCITY(DOF_SUB1:DOF_SUP1)
    !
    ! Locals:
    ! MYUNIT: Current unit number to open restart file.
    ! ISURF: Generic loop index to loop over mesh surfaces.
    !
    INTEGER :: MYUNIT
    INTEGER :: ISURF
    INTEGER :: RST_NUM
    LOGICAL :: FILE_EXISTS
    CHARACTER(LEN=8) :: RST_NUM_STR
    CHARACTER(LEN=64) :: FILENAME
    REAL(RK) :: DIFF_STRESS
    INTEGER :: I
    !
    !---------------------------------------------------------------------------
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
    FILENAME = 'rst'//TRIM(RST_NUM_STR)//'.control'
    !
    OPEN(UNIT = MYUNIT, FILE = FILENAME, FORM = 'UNFORMATTED', ACTION = 'READ')
    !
    READ(MYUNIT) ISTEP
    READ(MYUNIT) CURR_LOAD
    READ(MYUNIT) PREV_LOAD
    READ(MYUNIT) STEP_COMPLETE
    READ(MYUNIT) DTIME
    READ(MYUNIT) INCR
    READ(MYUNIT) TIME
    !
    DO ISURF = 1, NSURFACES
        !
        READ(MYUNIT) SURF_LOAD_ARRAY(ISURF,:)
        !
    ENDDO
    !
    READ(MYUNIT) AREA
    READ(MYUNIT) AREA0
    READ(MYUNIT) LENGTH
    READ(MYUNIT) LENGTH0
    READ(MYUNIT) CURR_VEL
    READ(MYUNIT) S_PERT_MAG
    READ(MYUNIT) T_PERT_MAG
    !
    IF (MYID .EQ. 0) THEN
        !
        WRITE(DFLT_U,'(A)') 'Info   : Reading restart control information...'
        WRITE(DFLT_U,'(A)') 'Info   :   - Previous simulation ended with final:'
        WRITE(DFLT_U,'(A, I0)') 'Info   :     > Increments: ', INCR
        WRITE(DFLT_U,'(A, I0)') 'Info   :     > Steps:      ', ISTEP - 1
        WRITE(DFLT_U,'(A, E14.4)') 'Info   :     > Time:  ', TIME
        WRITE(DFLT_U,'(A, 3(E14.4))') 'Info   :     > Normal loads:  ', &
            & CURR_LOAD
        !
    ENDIF
    !
    CLOSE(MYUNIT)
    !
    ! Reinitialize step and time if new files
    ! In the future, we will have to handle the APPEND case, which will continue
    ! sequentially, rather than reinitializing
    !
    IF (OPTIONS%RESTART_FILE_HANDLING .EQ. 1) THEN ! New files
        !
        ISTEP = 1
        STEP_COMPLETE = .TRUE.
        INCR = 0
        TIME = 0.0D0
        !
    !ELSE IF (OPTIONS%RESTART_FILE_HANDLING .EQ. 0) THEN ! Append
    END IF
    !
    ! Set a temporary variable to assist in printing to post.report file
    !
    OPTIONS%RESTART_INITIAL_STEP = ISTEP
    !
    ! Check for a reversal in loading direction upon restart, correct loading
    ! history if so. The loading history as assigned in PROCESS_CTRL_DATA_CSR
    ! assumes that we are starting from 0 load. Below rectifies any directional
    ! issues with loading due to this assumption.
    !
    IF ((CURR_LOAD(BCS_OPTIONS%LOADING_DIRECTION + 1) - &
        & PREV_LOAD(BCS_OPTIONS%LOADING_DIRECTION + 1)) * &
        & (TARGET_LOAD(1, BCS_OPTIONS%LOADING_DIRECTION + 1) - &
        & CURR_LOAD(BCS_OPTIONS%LOADING_DIRECTION + 1)) .LT. 0.0D0) THEN
        !
        ! Change velocity direction
        !
        VELOCITY = -1.0D0 * VELOCITY
        !
        ! And check for any load reversals on subsequent steps
        !
        IF (NSTEPS .GE. 2) THEN
            !
            ! Calculate TARGET_SIGN and IS_SIGNED for the first step
            !
            DIFF_STRESS = (TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(1, &
                & BCS_OPTIONS%LOADING_DIRECTION + 1) - &
                & CURR_LOAD(BCS_OPTIONS%LOADING_DIRECTION + 1))
            !
            IF (DIFF_STRESS .GT. 0.0D0) THEN
                !
                IS_SIGNED(1)   = .FALSE.
                TARGET_SIGN(1) = 1
                !
            ELSE IF (DIFF_STRESS .LT. 0.0D0) THEN
                !
                IS_SIGNED(1) = .TRUE.
                TARGET_SIGN(1) = -1

                !
            ELSE ! The differential change is zero, thus dt is zero.
                !
                CALL PAR_QUIT('Error  :     > Load differential is zero&
                    & between steps.')
                !
            END IF
            !
            ! Reset VEL_FACTOR for reassignment in the below loops
            !
            VEL_FACTOR = ABS(VEL_FACTOR)
            !
            ! Calculate VEL_FACTOR for each step after the first
            !
            DO I = 2, NSTEPS
                !
                DIFF_STRESS = &
                    & (TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I, &
                    & BCS_OPTIONS%LOADING_DIRECTION + 1)) &
                    & - (TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I - 1, &
                    & BCS_OPTIONS%LOADING_DIRECTION + 1))
                !
                IF (DIFF_STRESS .GT. 0.0D0) THEN
                    !
                    IS_SIGNED(I)   = .FALSE.
                    TARGET_SIGN(I) = 1
                    !
                ELSE IF (DIFF_STRESS .LT. 0.0D0) THEN
                    !
                    IS_SIGNED(I) = .TRUE.
                    TARGET_SIGN(I) = -1
                    !
                ELSE ! The differential change is zero, thus dt is zero.
                    !
                    CALL PAR_QUIT('Error  :     > Load differential is zero&
                        & between steps.')
                    !
                END IF
                !
                ! Check if a load reversal occured since the last step.
                !
                IF (IS_SIGNED(I) .NEQV. IS_SIGNED(I-1)) THEN
                    !
                    VEL_FACTOR(I) = VEL_FACTOR(I) * (-1)
                    !
                ELSE IF (IS_SIGNED(I) .EQV. IS_SIGNED(I-1)) THEN
                    !
                    VEL_FACTOR(I) = VEL_FACTOR(I) * (1)
                    !
                ELSE
                    !
                    CALL PAR_QUIT('Error  :     > Load reversal check failed.')
                    !
                END IF
                !
            END DO
            !
        END IF
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_TRIAXCSR_RESTART
    !
END MODULE DRIVER_TRIAXCSR_MOD
