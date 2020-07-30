! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE DRIVER_TRIAXCSR_MOD
!
! Driver for triaxial load control with constant strain rate.
!
! Contains subroutines:
! DRIVER_TRIAX_CSR: Primary driver for triaxial CSR control simulations.
! READ_CTRL_DATA_CSR: Read data for control and modify initial velocity field.
! READ_TRIAXCSR_RESTART: Read triaxial CSR control restart information.
! WRITE_TRIAXCSR_RESTART: Write triaxial CSR control restart information.
!
USE LIBF95, ONLY: NEWUNITNUMBER
USE INTRINSICTYPESMODULE, RK=>REAL_KIND
!
USE DIMSMODULE, ONLY: DIMS1, TVEC, TVEC1, MAXSLIP1, NGRAIN1, ANISOTROPIC_EVPS, &
    & VTINY
USE DRIVER_UTILITIES_MOD, ONLY: VELOCITY_ITERATION, CALC_MESH_DIM, EST_AVG_MOD,&
    & SOLVE_LIN_SYS_3, UPDATE_STATE_EVPS, CALC_STRESS_STRAIN
USE GATHER_SCATTER, ONLY: TRACE
USE ITMETHODEVPSMODULE, ONLY: ITMETHOD_EVPS
USE POWDER_DIFFRACTION_MOD, ONLY: RUN_POWDER_DIFFRACTION
USE MICROSTRUCTURE_MOD, ONLY: NUMPHASES, GAMMADOT
USE READ_INPUT_MOD, ONLY: KDIM1, EL_SUB1, EL_SUP1, DOF_SUB1, DOF_SUP1, COORDS
USE PARALLEL_MOD, ONLY: MYID, PAR_MESSAGE, PAR_QUIT
USE RESTARTMODULE, ONLY: RESTARTWRITEFIELD, RESTARTREADFIELD
USE SIMULATION_CONFIGURATION_MOD, ONLY: OPTIONS, BCS_OPTIONS, UNIAXIAL_OPTIONS,&
    & POWDER_DIFFRACTION_OPTIONS, TRIAXCSR_OPTIONS, PRINT_OPTIONS
USE SURF_INFO_MOD, ONLY: NSURFACES, COMPUTE_AREA
USE UNITS_MOD, ONLY: DFLT_U, FORCE_U1, FORCE_U2, FORCE_U3, FORCE_U4, FORCE_U5,&
    & FORCE_U6, CONV_U, OUNITS, BCS_ITER_1_U, REPORT_U
USE WRITE_OUTPUT_MOD, ONLY: PRINT_STEP, WRITE_REPORT_FILE_COMPLETE_STEPS, &
    & WRITE_FORCE_FILE_HEADERS, WRITE_FORCE_FILE_DATA, WRITE_CONV_FILE_HEADERS,&
    & WRITE_TRIAXCSR_RESTART
USE ELEMENTAL_VARIABLES_UTILS_MOD, ONLY: PLASTICVELGRADSYMSKW, CALC_TOTAL_WORK,&
    & CALC_PLASTIC_WORK
!
IMPLICIT NONE
!
! Private
!
INTEGER, PRIVATE :: NSTEPS
REAL(RK), ALLOCATABLE, PRIVATE :: TARGET_LOAD(:,:), TARGET_SIGN(:), &
    & VEL_FACTOR(:), TIME_STEP(:), TIME_STEP_MIN(:)
INTEGER, ALLOCATABLE, PRIVATE :: PRINT_FLAG(:)
!
CONTAINS
    !
    SUBROUTINE driver_triax_csr(bcs, velocity, pforce, dtrace, ntrace, &
        & c0_angs, crss_n, rstar_n, keinv, wts, auto_time, gammadot)
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
    ! NTRACE: Gather/scatter trace for nodal points.
    ! C0_ANGS: Initial orientations, rotation matrix.
    ! CRSS_N: Resolved shear stresses on each slip system.
    ! RSTAR_N: Current orientations.
    ! KEINV: Inverse of the single crystal elasticity tensor.
    ! WTS: Elemental weights (??).
    ! AUTO_TIME: (??).
    !
    LOGICAL  :: BCS(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: VELOCITY(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: PFORCE  (DOF_SUB1:DOF_SUP1)
    TYPE(TRACE) DTRACE
    TYPE(TRACE) NTRACE
    REAL(RK) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: KEINV(0:TVEC1,1:NUMPHASES)
    REAL(RK) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER  :: AUTO_TIME
    REAL(RK) :: GAMMADOT(0:MAXSLIP1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    !
    ! Locals:
    ! Need to be defined - JC
    !
    LOGICAL  :: CONVERGED_SOLUTION
    LOGICAL  :: REJECT_SOLUTION
    LOGICAL  :: STEP_COMPLETE
    LOGICAL  :: IS_NECKING
    LOGICAL  :: IS_LIMIT_TRIPPED
    INTEGER  :: ITMETHOD_EVPS
    INTEGER  :: INCR
    INTEGER  :: ITERNL
    INTEGER  :: IER
    INTEGER  :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER  :: IOSTATUS
    REAL(RK) :: DTIME, DTIME_OLD, TIME
    REAL(RK) :: D(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: D_VEC     (0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TRIAL_VELOCITY(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: VGRAD (0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: EVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ECOORDS(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: RSTAR (0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TRSTAR (0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG_VEC  (0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ELAS_TOT6(0:TVEC , 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: S_AVG_3X3(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ELPRESS        (EL_SUB1:EL_SUP1)
    REAL(RK) :: E_ELAS_KK_BAR  (EL_SUB1:EL_SUP1)
    REAL(RK) :: E_ELAS_KK      (EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG_KK         (EL_SUB1:EL_SUP1)
    REAL(RK) :: STIF(TVEC, TVEC, EL_SUB1:EL_SUP1)
    REAL(RK) :: FE        (TVEC, EL_SUB1:EL_SUP1)
    INTEGER  :: M_EL, I, J, K, II
    INTEGER  :: ISTEP, IPHASE    
    REAL(RK) :: SURF_LOAD_ARRAY(NSURFACES,3)
    REAL(RK) :: AREA0(NSURFACES), AREA(NSURFACES)
    REAL(RK) :: LENGTH(3), LENGTH0(3)
    REAL(RK) :: GAMMA(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DP_HAT(0:TVEC1,EL_SUB1:EL_SUP1)
    REAL(RK) :: WP_HAT(0:DIMS1,EL_SUB1:EL_SUP1)
    REAL(RK) :: DEFF       (EL_SUB1:EL_SUP1)
    REAL(RK) :: DPEFF      (EL_SUB1:EL_SUP1)
    REAL(RK) :: EQSTRAIN   (EL_SUB1:EL_SUP1)
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
    REAL(RK), ALLOCATABLE :: EL_WORK_N(:)
    REAL(RK), ALLOCATABLE :: EL_WORKP_N(:)
    REAL(RK), ALLOCATABLE :: EL_WORK_RATE_N(:)
    REAL(RK), ALLOCATABLE :: EL_WORKP_RATE_N(:)
    REAL(RK), ALLOCATABLE :: EL_WORK(:)
    REAL(RK), ALLOCATABLE :: EL_WORKP(:)
    REAL(RK), ALLOCATABLE :: D_TOT(:, :, :)
    !
    !----------------------------------------------------------------------
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
        CALL RESTARTREADFIELD(VELOCITY, C0_ANGS, C_ANGS, &
            & RSTAR, RSTAR_N, WTS, CRSS, CRSS_N, &
            & E_ELAS_KK_BAR, SIG_VEC_N, EQSTRAIN, EQPLSTRAIN, GAMMA)
        !       
        CALL READ_TRIAXCSR_RESTART(ISTEP, CURR_LOAD, PREV_LOAD, STEP_COMPLETE, &
            & DTIME, INCR, TIME, SURF_LOAD_ARRAY, &
            & AREA, AREA0, LENGTH, LENGTH0, CURR_VEL, S_PERT_MAG, T_PERT_MAG)
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
        E_ELAS_KK_BAR = 0.0_RK
        SIG_VEC_N = 0.0_RK
        C_ANGS = C0_ANGS
        !
        ! Initialize total deformation rate tensor (optional)
        !
        IF ((PRINT_OPTIONS%PRINT_WORK) .OR. (PRINT_OPTIONS%PRINT_DEFRATE)) THEN
            !
            ALLOCATE(D_TOT(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1))
            D_TOT = 0.0_RK
            !
        END IF
        !
        ! Initialize work arrays iff it needs to be printed
        !
        IF (PRINT_OPTIONS%PRINT_WORK) THEN
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
        ! Initialize work-pl arrays iff it needs to be printed
        !
        IF (PRINT_OPTIONS%PRINT_WORKP) THEN
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
        ! Initialize integrated quantities
        !
        GAMMA = 0.0_RK
        EQPLSTRAIN = 0.0_RK
        EQSTRAIN   = 0.0_RK
        !
        ! Initialize deformation control
        !
        INCR = 0
        TIME = 0.0_RK
        ISTEP = 1
        INCR = 0
        STEP_COMPLETE = .FALSE.
        DTIME = TIME_STEP(1)   
        S_PERT_MAG = MIN_PERT_FRAC*ABS(INITIAL_VEL)
        T_PERT_MAG = MIN_PERT_FRAC*ABS(INITIAL_VEL)
        CURR_EQSTRAIN = 0.0_RK
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
        CALL PRINT_STEP(0, 0, COORDS, VELOCITY, C0_ANGS, WTS, CRSS_N)
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
    ! IF (MYID .EQ. 0) WRITE(DFLT_U, '(A)') 'Info   :   - Starting T-stepping loop&
    ! & for anisotropic viscoplastic solution'
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
        GAMMADOT = 0.0_RK
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
                WRITE(DFLT_U, '(A,I0)') 'Info   :     > TRIAL_BC: Iteration ', BC_ITER   
                !         
            END IF
            !
            ! Begin velocity iteration
            !
            CALL VELOCITY_ITERATION(BCS, PFORCE, VELOCITY,& 
                & ELPRESS, EVEL, CURR_LOAD, DTRACE, NTRACE,&
                & C0_ANGS, C_ANGS, SIG_VEC_N, SIG_VEC, CRSS_N,&
                & CRSS, RSTAR_N, RSTAR, TRSTAR, KEINV,&
                & E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, JITER_STATE,&
                & WTS, DEFF, DTIME, INCR, E_BAR_VEC,&
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
                    WRITE(DFLT_U,'(A,3(E12.2))') 'Info   :       . Velocity:    ', CURR_VEL
                    !
                ELSE
                    !
                    WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       . Velocity:    ', CURR_VEL
                    !
                END IF
                !
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       . Load:        ', CURR_LOAD
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       . Target load: ', IDEAL_LOAD
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
                WRITE(DFLT_U,'(A)') 'Info   :       . Perturbing surface velocity&
                    & in secondary direction'
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
            CALL VELOCITY_ITERATION(BCS, PFORCE, TRIAL_VELOCITY,& 
                & ELPRESS, EVEL, TRIAL_LOAD, DTRACE, NTRACE,&
                & C0_ANGS, C_ANGS, SIG_VEC_N, SIG_VEC, CRSS_N,&
                & CRSS, RSTAR_N, RSTAR, TRSTAR, KEINV,&
                & E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, JITER_STATE,&
                & WTS, DEFF, DTIME, INCR, E_BAR_VEC,&
                & CONVERGED_SOLUTION, AUTO_TIME, ITERNL)
            !
            IF (MYID .EQ. 0) THEN 
                !            
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :     > Trial load: ', TRIAL_LOAD
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
                WRITE(DFLT_U,'(A)') 'Info   :       . Perturbing surface velocity&
                    & in tertiary direction'
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
            CALL VELOCITY_ITERATION(BCS, PFORCE, TRIAL_VELOCITY,& 
                & ELPRESS, EVEL, TRIAL_LOAD, DTRACE, NTRACE,&
                & C0_ANGS, C_ANGS, SIG_VEC_N, SIG_VEC, CRSS_N,&
                & CRSS, RSTAR_N, RSTAR, TRSTAR, KEINV,&
                & E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, JITER_STATE,&
                & WTS, DEFF, DTIME, INCR, E_BAR_VEC,&
                & CONVERGED_SOLUTION, AUTO_TIME, ITERNL)
            !
            IF (MYID .EQ. 0) THEN 
                !            
                WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :     > Trial load: ', TRIAL_LOAD 
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
            & E_ELAS_KK, SIG_KK, KEINV, WTS, DEFF, DTIME, STIF, FE, D, D_VEC, &
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
            MACRO_ENG_STRAIN(II) = LENGTH(II) / LENGTH0(II) - 1.0_RK
            !
        END DO
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
        ! Update EQSTRAIN
        !
        EQSTRAIN=EQSTRAIN+DEFF*DTIME
        !
        ! Update slip system shear
        !
        GAMMA=GAMMA+GAMMADOT*DTIME
        !
        ! Reconstruct total deformation rate tensor if requested
        !
        IF ((PRINT_OPTIONS%PRINT_WORK) .OR. (PRINT_OPTIONS%PRINT_DEFRATE)) THEN
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
        IF (PRINT_OPTIONS%PRINT_WORK) THEN
            !
            CALL CALC_TOTAL_WORK(DTIME, D_TOT, S_AVG_3X3, EL_WORK_N, &
                & EL_WORK_RATE_N, EL_WORK)
            !
        END IF
        !
        IF (PRINT_OPTIONS%PRINT_WORKP) THEN
            !
            CALL CALC_PLASTIC_WORK(DTIME, DP_HAT, C_ANGS, S_AVG_3X3, &
                & EL_WORKP_N, EL_WORKP_RATE_N, EL_WORKP)
            !
        END IF
        !
        ! Check if the target load is reached
        !
        STEP_COMPLETE = .FALSE.
        !
        IF (TARGET_SIGN(ISTEP)*(CURR_LOAD(PDIR)-TARGET_LOAD(ISTEP,PDIR))&
            & + OPTIONS%LOAD_TOL .GE. 0.0) THEN
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
                    & C_ANGS, WTS, CRSS, ELAS_TOT6, S_AVG_3X3, DEFF, EQSTRAIN, &
                    & DPEFF, EQPLSTRAIN, VGRAD, DP_HAT, WP_HAT, GAMMA, &
                    & GAMMADOT, EL_WORK, EL_WORKP, D_TOT)
                !
                IF (PRINT_OPTIONS%PRINT_RESTART) THEN
                    !
                    CALL RESTARTWRITEFIELD(VELOCITY, C0_ANGS, C_ANGS, &
                        & RSTAR, RSTAR_N, WTS, CRSS, CRSS_N, &
                        & E_ELAS_KK_BAR, SIG_VEC_N, EQSTRAIN, EQPLSTRAIN, GAMMA)
                    !
                    CALL WRITE_TRIAXCSR_RESTART(ISTEP+1, CURR_LOAD, PREV_LOAD, &
                        & STEP_COMPLETE, DTIME, INCR, TIME, SURF_LOAD_ARRAY, &
                        & AREA, AREA0, LENGTH, LENGTH0, CURR_VEL, S_PERT_MAG, &
                        & T_PERT_MAG)
                    !
                END IF
                !
                IF (POWDER_DIFFRACTION_OPTIONS%RUN_POWDER_DIFFRACTION) THEN
                    !   
                    CALL RUN_POWDER_DIFFRACTION(ISTEP, C_ANGS, ELAS_TOT6, DPEFF, CRSS)
                    !
                ENDIF
                !
            ENDIF
            !
            ! Check that maximum strain has not been exceeded
            !
            IF (MAXVAL(ABS(MACRO_ENG_STRAIN)) .GE. MAX_STRAIN) THEN
                !
                CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP)
                !
                CALL PAR_QUIT('Error  :     > Maximum eng. strain exceeded.')
                !
            ENDIF
            !
            ! Check that maximum eqstrain has not been exceeded
            !
            IF (CURR_EQSTRAIN .GE. MAX_EQSTRAIN) THEN
                !
                CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP)
                !
                CALL PAR_QUIT('Error  :     > Maximum eqv. strain exceeded.')
                !
            ENDIF
            !
            IF (IS_NECKING) THEN
                !
                CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP)
                !
                CALL PAR_QUIT('Error  :     > Specimen is necking.')
                !
            ENDIF
            !
            IF (IS_LIMIT_TRIPPED) THEN
                !
                CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP)
                !
                CALL PAR_QUIT('Error  :     > Maximum time or maximum &
                    &increments exceeded.')
                !     
            ENDIF
            !
            IF (ISTEP .EQ. NSTEPS) THEN
                !
                CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP)
                !
                CALL PAR_QUIT('Info   : Final step terminated. Simulation&
                    & completed successfully.')
                !       
            ENDIF
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
    SUBROUTINE READ_CTRL_DATA_CSR(VELOCITY)
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
    ! MESSAGE: Temporary string for passing to PAR_MESSAGE.
    !
    INTEGER  :: I, II, PRIMARY_DIR
    REAL(RK) :: DIFF_STRESS, USER_STRAIN_RATE
    REAL(RK), ALLOCATABLE :: STEP_STRAIN_RATE(:)    
    LOGICAL, ALLOCATABLE  :: IS_SIGNED(:)
    CHARACTER(LEN=256) :: MESSAGE
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
    DIFF_STRESS = 0.0_RK
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
            CALL PAR_QUIT('Error  :     > Inputted `number_of_strain_rate_jumps`&
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
            TIME_STEP(I)     = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,4)
            TIME_STEP_MIN(I) = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,5)
            PRINT_FLAG(I)    = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,6)
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
            IF (DIFF_STRESS .GT. 0.0_RK) THEN
                !
                IS_SIGNED(I)   = .FALSE.
                TARGET_SIGN(I) = 1
                !VEL_FACTOR(I)  = 1
                !
            ELSE IF (DIFF_STRESS .LT. 0.0_RK) THEN
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
            TIME_STEP(I)     = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,4)
            TIME_STEP_MIN(I) = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,5)
            PRINT_FLAG(I)    = TRIAXCSR_OPTIONS%TARGET_CSR_LOAD(I,6)
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
            IF (DIFF_STRESS .GT. 0.0_RK) THEN
                !
                IS_SIGNED(I)   = .FALSE.
                TARGET_SIGN(I) = 1
                !
            ELSE IF (DIFF_STRESS .LT. 0.0_RK) THEN
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
    END SUBROUTINE READ_CTRL_DATA_CSR
    !
    !===========================================================================
    !
    SUBROUTINE READ_TRIAXCSR_RESTART(ISTEP,CURR_LOAD,PREV_LOAD,STEP_COMPLETE, &
        & DTIME,INCR,TIME,SURF_LOAD_ARRAY,AREA,AREA0,LENGTH,LENGTH0,CURR_VEL, &
        & S_PERT_MAG, T_PERT_MAG)
    !
    ! Read triaxial CSR restart information.
    !
    !--------------------------------------------------------------------------- 
    !
    ! Arguments:
    ! Needs to be defined - JC
    !
    LOGICAL, INTENT(OUT)  :: STEP_COMPLETE
    INTEGER, INTENT(OUT)  :: ISTEP
    INTEGER, INTENT(OUT)  :: INCR
    REAL(RK), INTENT(OUT) :: CURR_LOAD(3), PREV_LOAD(3)  
    REAL(RK), INTENT(OUT) :: DTIME, TIME    
    REAL(RK), INTENT(OUT) :: SURF_LOAD_ARRAY(NSURFACES,3)
    REAL(RK), INTENT(OUT) :: AREA(NSURFACES), AREA0(NSURFACES)
    REAL(RK), INTENT(OUT) :: LENGTH(3), LENGTH0(3)
    REAL(RK), INTENT(OUT) :: CURR_VEL(3)
    REAL(RK), INTENT(OUT) :: S_PERT_MAG, T_PERT_MAG
    !
    ! Locals:
    ! MYUNIT: Current unit number to open restart file.
    ! ISURF: Generic loop index to loop over mesh surfaces.
    !
    INTEGER :: MYUNIT
    INTEGER :: ISURF
    !
    !---------------------------------------------------------------------------
    !
    MYUNIT = NEWUNITNUMBER()
    !
    OPEN(UNIT=MYUNIT, FILE=TRIM(OPTIONS%RSCTRL_IN), &
         &   FORM='UNFORMATTED', ACTION='READ')
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
        WRITE(DFLT_U,'(A)') 'Info   :   - Restart parameters:'
        WRITE(DFLT_U,'(A,I0)') 'Info   :     > Current  step: ', ISTEP
        WRITE(DFLT_U,'(A,I0)') 'Info   :     > Previous incr: ', INCR
        WRITE(DFLT_U,'(A,F12.6)') 'Info   :     > Previous time: ', TIME
        WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :     > Previous load: ', PREV_LOAD
        WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :     > Loads on X0:   ', SURF_LOAD_ARRAY(1,:)
        WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       Loads on X1:   ', SURF_LOAD_ARRAY(2,:)
        WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       Loads on Y0:   ', SURF_LOAD_ARRAY(3,:)
        WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       Loads on Y1:   ', SURF_LOAD_ARRAY(4,:)
        WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       Loads on Z0:   ', SURF_LOAD_ARRAY(5,:)
        WRITE(DFLT_U,'(A,3(F12.6))') 'Info   :       Loads on Z1:   ', SURF_LOAD_ARRAY(6,:)
!
!        WRITE(DFLT_U,'(a36)') 'Reading restart control information'
!        WRITE(DFLT_U,'(a16, i8)')       'current_step  : ', ISTEP
!        WRITE(DFLT_U,'(a16, 3(e14.6))') 'current_load  : ', CURR_LOAD
!        WRITE(DFLT_U,'(a16, 3(e14.6))') 'previous_load : ', PREV_LOAD
!        WRITE(DFLT_U,'(a16, e14.6)')    'dtime         : ', DTIME
!        WRITE(DFLT_U,'(a8, i8)')        'incr     : ', INCR
!        WRITE(DFLT_U,'(a8, e14.6)')     'time     : ', TIME
!        WRITE(DFLT_U,'(a8, 3(e14.6))')  'load     : ', SURF_LOAD_ARRAY(1,:)
        !
!        DO ISURF = 2, NSURFACES
            !        
!            WRITE(DFLT_U,'(a8, 3(e14.6))')  '        ', SURF_LOAD_ARRAY(ISURF,:)
            !       
!        ENDDO
        !
!        WRITE(DFLT_U,'(a8, 3(e14.6))')  'area     : ', AREA
!        WRITE(DFLT_U,'(a8, 3(e14.6))')  'area0    : ', AREA0
!        WRITE(DFLT_U,'(a8, 3(e14.6))')  'length   : ', LENGTH
!        WRITE(DFLT_U,'(a8, 3(e14.6))')  'length0  : ', LENGTH0
!        WRITE(DFLT_U,'(a8, 3(e14.6))')  'curr_vel : ', CURR_VEL
!        WRITE(DFLT_U,'(a16, e14.6)')    's_pert_mag   : ', S_PERT_MAG
!        WRITE(DFLT_U,'(a16, e14.6)')    't_pert_mag   : ', T_PERT_MAG
!        WRITE(DFLT_U,*) ' '
        !
    ENDIF
    !
    CLOSE(MYUNIT)
    !
    RETURN
    !
    END SUBROUTINE READ_TRIAXCSR_RESTART
    !
END MODULE DRIVER_TRIAXCSR_MOD
