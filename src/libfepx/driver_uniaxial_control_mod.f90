! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE DRIVER_UNIAXIAL_CONTROL_MOD
!
! Driver for uniaxial strain control with either load or strain targets.
!
! Contains subroutines:
! DRIVER_UNIAXIAL_CONTROL: Primary driver for uniaxial control simulations.
! PROCESS_CTRL_DATA: Read data for load control and modify init. velocity field.
! READ_UNIAXIAL_RESTART: Read uniaxial control restart information.
!
! Contains functions:
! CHECK_LIMIT: Check if specimen time or increment limit is tripped.
! CHECK_NECKING: Check if specimen is necking.
! GET_ALL_STEPS_COMPLETE: Check if all steps complete.
! GET_CURRENT_STEP: Returns current target index.
! GET_INITIAL_DTIME: Returns first time step.
! GET_PRINT_FLAG: Returns current print flag.
! GET_STEP_COMPLETE: Set the current step as complete.
! GET_TARGET_INDEX: Returns current target index.
! GET_TARGET_SURFACE: Returns current target surface.
! GET_VEL_FACTOR: Check if load reverses direction.
! IN_RANGE: Check if load is in specified range.
! TIME_INCREMENT: Return TIME increment according to loading step.
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
USE LIBF95, ONLY: NEWUNITNUMBER
!
! From libfepx:
!
USE DIMENSIONS_MOD, ONLY: DIMS1, TVEC, TVEC1, MAXSLIP1, NGRAIN1, &
    & ANISOTROPIC_EVPS
USE DRIVER_UTILITIES_MOD, ONLY: UPDATE_STATE_EVPS, CALC_STRESS_STRAIN
USE FIBER_AVERAGE_MOD, ONLY: RUN_FIBER_AVERAGE
USE ITERATE_STRESS_EVPS_MOD, ONLY: ITMETHOD_EVPS
USE KINEMATICS_MOD, ONLY: PLASTICVELGRADSYMSKW, CALC_TOTAL_WORK,&
    & CALC_PLASTIC_WORK
USE MATRIX_OPERATIONS_MOD, ONLY: CALC_ELVOL, STRAIN_EQUIV_3X3, VEC6_MAT_SYMM
USE MICROSTRUCTURE_MOD, ONLY: NUMPHASES, GAMMADOT, GACCUMSHEAR, ACCUMSHEAR_CEN
USE QUADRATURE_MOD, ONLY: NQPT1
USE READ_INPUT_MOD, ONLY: KDIM1, EL_SUB1, EL_SUP1, DOF_SUB1, DOF_SUP1, COORDS,&
    & OPTIONS, BCS_OPTIONS, UNIAXIAL_OPTIONS, FIBER_AVERAGE_OPTIONS, &
    & UNIAXIAL_LOAD_TARGET, UNIAXIAL_STRAIN_TARGET, PRINT_OPTIONS, NODES, &
    & READ_RESTART_FIELD, D_TOT, EL_WORK_N, EL_WORKP_N, EL_WORK_RATE_N, &
    & EL_WORKP_RATE_N, EL_WORK, EL_WORKP
USE SURFACE_MOD, ONLY: NSURFACES, COMPUTE_AREA
USE UNITS_MOD, ONLY: DFLT_U, FORCE_U1, FORCE_U2, FORCE_U3, FORCE_U4, FORCE_U5,&
    & FORCE_U6, CONV_U, OUNITS, REPORT_U
USE WRITE_OUTPUT_MOD, ONLY: PRINT_STEP, WRITE_REPORT_FILE_COMPLETE_STEPS, &
    & WRITE_FORCE_FILE_HEADERS, WRITE_FORCE_FILE_DATA, WRITE_CONV_FILE_HEADERS,&
    & WRITE_UNIAXIAL_RESTART, WRITE_RESTART_FIELD
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
INTEGER, PRIVATE :: NSTEPS, CURRENT_STEP, TARGET_SURFACE
REAL(RK), ALLOCATABLE, PRIVATE :: TARGET_LOAD(:), TARGET_SIGN(:), &
    & VEL_FACTOR(:), TIME_STEP(:), TIME_STEP_MIN(:)
REAL(RK), PRIVATE :: PREV_STRAIN, CURR_STRAIN
REAL(RK), ALLOCATABLE, PRIVATE :: STEP_STRAIN_RATE(:)
LOGICAL, ALLOCATABLE, PRIVATE :: IS_SIGNED(:)
INTEGER, ALLOCATABLE, PRIVATE :: TARGET_INCR(:), TARGET_INDEX(:), PRINT_FLAG(:)
REAL(RK), DIMENSION(3), PRIVATE :: &
    & CURRENT_LOAD  =(/0.0d0,0.0d0,0.0d0/),&
    & PREVIOUS_LOAD =(/0.0d0,0.0d0,0.0d0/)
LOGICAL, PRIVATE :: STEP_COMPLETE, ALL_STEPS_COMPLETE
REAL(RK), PRIVATE :: DTIME_OLD
!
! Public
!
PUBLIC
!
CONTAINS
    !
    SUBROUTINE DRIVER_UNIAXIAL_CONTROL(BCS, VELOCITY, PFORCE, DTRACE, C0_ANGS, &
        & CRSS_N, RSTAR_N, KEINV, WTS, AUTO_TIME, GAMMADOT, CLOCK_START)
    !
    ! Primary driver for uniaxial strain control with load or strain targeting.
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
    ! GAMMADOT: Elemental slip system shearing rates.
    ! COUNT_INIT: Timer start time
    !
    LOGICAL  :: BCS(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: VELOCITY(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: PFORCE(DOF_SUB1:DOF_SUP1)
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
    LOGICAL  :: REJECT_SOLUTION
    LOGICAL  :: IS_NECKING
    LOGICAL  :: IS_LIMIT_TRIPPED
    INTEGER  :: INCR
    INTEGER  :: ITERNL
    INTEGER  :: IER
    INTEGER  :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DTIME, DTIME_N, TIME
    REAL(RK) :: D(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: D_VEC(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK) :: VGRAD(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: EVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ELAS_TOT6(0:TVEC, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ELAS_TOT_3X3(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: D_TOT(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TOTSTRAIN(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: S_AVG_3X3(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: E_ELAS_KK_BAR(EL_SUB1:EL_SUP1)
    REAL(RK) :: E_ELAS_KK(EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG_KK(EL_SUB1:EL_SUP1)
    REAL(RK) :: CRSS_Q(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: ELPRESS_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: E_ELAS_KK_BAR_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: E_ELAS_KK_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: SIG_KK_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: SIG_VEC_Q(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: SIG_VEC_N_Q(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: E_BAR_VEC_Q(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: DEFF_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: STIF(TVEC, TVEC, EL_SUB1:EL_SUP1)
    REAL(RK) :: FE(TVEC, EL_SUB1:EL_SUP1)
    INTEGER  :: M_EL, I
    INTEGER  :: ISTEP
    REAL(RK) :: LOAD(NSURFACES,3)
    REAL(RK) :: AREA0(NSURFACES), AREA(NSURFACES)
    REAL(RK) :: GAMMA(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DP_HAT(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK) :: PLSTRAIN(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK) :: WP_HAT(0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DEFF(EL_SUB1:EL_SUP1)
    REAL(RK) :: DPEFF(EL_SUB1:EL_SUP1)
    REAL(RK) :: EQSTRAIN(EL_SUB1:EL_SUP1)
    REAL(RK) :: EQELSTRAIN(EL_SUB1:EL_SUP1)
    REAL(RK) :: EQPLSTRAIN(EL_SUB1:EL_SUP1)
    INTEGER  :: LOAD_INDEX, LOAD_SURFACE, NTSTEPS
    CHARACTER(LEN=16) :: FIELD, TIME_STRING, DTIME_STRING
    REAL(RK) :: IDUMMY(3)
    INTEGER  :: RINCR
    REAL(RK) :: DKK_STEP
    !
    REAL(RK), ALLOCATABLE :: ECOORDS(:, :)
    REAL(RK), ALLOCATABLE :: ELVOL(:)
    !
    REAL(RK) :: CLOCK_END
    !
    !---------------------------------------------------------------------------
    !
    ! Initialization
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    LOAD_SURFACE = GET_TARGET_SURFACE()
    IS_NECKING = .FALSE.
    NTSTEPS = 0
    !
    IF (OPTIONS%RESTART) THEN
        !
        CALL READ_RESTART_FIELD(VELOCITY, C0_ANGS, C_ANGS, &
            & RSTAR, RSTAR_N, WTS, CRSS, CRSS_N, &
            & E_ELAS_KK_BAR, SIG_VEC_N, EQSTRAIN, EQPLSTRAIN, GAMMA, &
            & EL_WORK_N, EL_WORKP_N, EL_WORK_RATE_N, EL_WORKP_RATE_N, &
            & PLSTRAIN, TOTSTRAIN)
        !
        CALL READ_UNIAXIAL_RESTART(INCR, TIME, LOAD, PREVIOUS_LOAD, AREA, &
            & AREA0, VELOCITY, TIME_STEP)
        !
        RINCR = INCR
        !
        ! Update initial guess of velocity field for load reversal
        ! STEP_COMPLETE resets after call to TIME_INCREMENT
        !
        IF (GET_STEP_COMPLETE()) THEN
            !
            VELOCITY = GET_VEL_FACTOR()*VELOCITY
            !
        ENDIF
        !
        DTIME_N = TIME_INCREMENT(LOAD(LOAD_SURFACE,:), INCR)
        !
    ELSE
        !
        ! Initialize areas and load arrays
        !
        LOAD=0.0D0
        AREA=0.0D0
        AREA0=0.0D0
        !
        ! Compute initial area (AREA0)
        !
        CALL COMPUTE_AREA(COORDS, AREA0)
        !
        ! Initialize state
        !
        E_ELAS_KK_BAR_Q = 0.0D0
        SIG_VEC_N_Q = 0.0D0
        E_ELAS_KK_BAR = 0.0D0
        SIG_VEC_N = 0.0D0
        C_ANGS = C0_ANGS
        !
        ! Initialize elvol array (and associated) iff it needs to be printed
        !
        IF ((PRINT_OPTIONS%PRINT_ELVOL) .OR. &
            & (FIBER_AVERAGE_OPTIONS%RUN_FIBER_AVERAGE)) THEN
            !
            ALLOCATE(ELVOL(EL_SUB1:EL_SUP1))
            ALLOCATE(ECOORDS(0:KDIM1, EL_SUB1:EL_SUP1))
            !
            ELVOL = 0.0D0
            !
        END IF
        !
        ! Initialize integrated quantities
        !
        EQSTRAIN = 0.0D0
        EQPLSTRAIN = 0.0D0
        PLSTRAIN = 0.0D0
        D_TOT = 0.0D0
        TOTSTRAIN = 0.0D0
        GAMMA = 0.0D0
        !
        ! Initialize deformation control
        !
        INCR = 0
        TIME = 0.0D0
        DTIME_N = GET_INITIAL_DTIME()
        !
    END IF
    !
    IF (MYID .EQ. 0) THEN
        !
        ! Write the header to file post.force#
        !
        IF (PRINT_OPTIONS%PRINT_FORCES) THEN
            !
            CALL WRITE_FORCE_FILE_HEADERS(1)
            !
            IF (OPTIONS%RESTART .EQV. .FALSE.) THEN
                !
                ! If virgin sample, write 0th step
                !
                IDUMMY = 0.0D0
                CALL WRITE_FORCE_FILE_DATA(1, 0, 0, LOAD, AREA0, 0.0D0, &
                    & IDUMMY)
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
    END IF
    !
    ISTEP = GET_CURRENT_STEP()
    !
    IF ((ISTEP .EQ. 1) .AND. (MYID .EQ. 0)) THEN
        !
        WRITE(DFLT_U,'(A,I0,A)') 'Info   : Running step ', ISTEP, '...'
        !
    END IF
    !
    ! Time stepping loop
    !
    TIME_STEPPING : DO
        !
        ! Iterate to find new configuration and material state
        !
        INCR = INCR + 1
        ISTEP = GET_CURRENT_STEP()
        !
        IF ((GET_STEP_COMPLETE()) .AND. (MYID .EQ. 0)) THEN
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
        ! Use the DTIME provided by TIME_INCREMENT in the previous increment
        !
        DTIME = DTIME_N
        !
        ! Update the total time
        !
        TIME = TIME + DTIME
        !
        IF (MYID .EQ. 0) THEN
            !
            IF (GET_STEP_COMPLETE() .OR. (INCR .EQ. 1)) THEN
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
            ELSEIF ((OPTIONS%RESTART .EQV. .TRUE.) .AND. &
                    & RINCR .EQ. (INCR - 1)) THEN
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
            ELSE
                !
                WRITE(FIELD, '(F0.4)') TIME
                IF (FIELD(1:1) == '.') FIELD = '0' // FIELD
                TIME_STRING = FIELD
                WRITE(FIELD, '(F0.4)') DTIME
                IF (FIELD(1:1) == '.') FIELD = '0' // FIELD
                DTIME_STRING = FIELD
                !
                WRITE(DFLT_U,'(/,A,I0,A,A,A,A,A)') 'Info   :   - &
                    &Increment ', INCR, ': t = ', TRIM(TIME_STRING), '&
                    & secs, dt = ', TRIM(DTIME_STRING), ' secs'
                !
            END IF
            !
        END IF
        !
        ! Update initial guess of velocity field for load reversal
        !
        IF (GET_STEP_COMPLETE()) THEN
            !
            VELOCITY = GET_VEL_FACTOR()*VELOCITY
            NTSTEPS = 0
            !
        END IF
        !
        IF((GET_VEL_FACTOR() .EQ. -1) .AND. (NTSTEPS .EQ. 0)) THEN
            !
            GACCUMSHEAR = 0.0D0
            ACCUMSHEAR_CEN = 0.0D0
            !
        END IF
        !
        NTSTEPS = NTSTEPS + 1
        !
        ! Compute:
        ! VELOCITY @(t+dt) using PCG solver
        !
        ! These are outputs of the subroutine/function but are useless and can
        ! be removed from the argument list
        ! E_BAR_VEC: OUT
        ! E_ELAS_KK_BAR: OUT
        ! SIG_VEC_N: INOUT (enters=0 at the first increment), could be OUT?
        !
        IER = ITMETHOD_EVPS(BCS, PFORCE, VELOCITY, ELPRESS_Q, EVEL, DTRACE, &
            & C0_ANGS, C_ANGS, SIG_VEC_N_Q, SIG_VEC_Q, CRSS_N, CRSS_Q, &
            & RSTAR_N, RSTAR, KEINV, E_ELAS_KK_BAR_Q, E_ELAS_KK_Q, SIG_KK_Q, &
            & JITER_STATE, WTS, DEFF_Q, DTIME, INCR, E_BAR_VEC_Q, &
            & CONVERGED_SOLUTION, AUTO_TIME, ITERNL)
        !
        IF (IER .LT. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to converge.')
            !
        END IF
        !
        REJECT_SOLUTION=.FALSE.
        !
        IF (.NOT. CONVERGED_SOLUTION) THEN
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(A)') 'Warning:     > AUTO_TIME .NE. 1 with&
                & CONVERGED_SOLUTION = FALSE'
                !
            END IF
            !
        END IF
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
        !
        ! Calculate stress, strain, load, and area
        !
        CALL CALC_STRESS_STRAIN(S_AVG_3X3, ELAS_TOT6, LOAD, AREA, SIG_VEC, &
            & SIG_KK, E_ELAS_KK, C_ANGS, KEINV, WTS)
        !
        ! Write loads to post.force# files
        !
        IF (MYID .EQ. 0) THEN
            !
            IF (PRINT_OPTIONS%PRINT_FORCES) THEN
                !
                IDUMMY = 0.0D0
                CALL WRITE_FORCE_FILE_DATA(1, ISTEP, INCR, &
                    & LOAD, AREA, TIME, IDUMMY)
                !
            END IF
            !
            WRITE(DFLT_U, '(A)', ADVANCE='NO') 'Info   :       . Loads on &
                &control surface:  '
            !
            DO I = 1, 3
                !
                WRITE(FIELD, '(F0.4)') LOAD(LOAD_SURFACE,I)
                IF (FIELD(1:1) == '.') FIELD = '0' // FIELD
                IF (FIELD(1:2) == '-.') FIELD = '-0.' // FIELD(3:)
                !
                WRITE(DFLT_U, '(A,A)', ADVANCE='NO') TRIM(FIELD), '  '
                !
            END DO
            !
        END IF
        !
        LOAD_INDEX = GET_TARGET_INDEX()
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
        GAMMA = GAMMA + GAMMADOT * DTIME
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
        IF ((PRINT_OPTIONS%PRINT_ELVOL) .OR. &
            & (FIBER_AVERAGE_OPTIONS%RUN_FIBER_AVERAGE)) THEN
            !
            CALL PART_GATHER(ECOORDS, COORDS, NODES, DTRACE)
            CALL CALC_ELVOL(ELVOL, ECOORDS)
            !
        END IF
        !
        ! Check if the target load is reached and find new DTIME
        !
        DTIME_N = TIME_INCREMENT(LOAD(LOAD_SURFACE, :), INCR)
        !
        ! Check for necking
        !
        IF (OPTIONS%CHECK_NECKING) THEN
            !
            IS_NECKING = CHECK_NECKING()
            !
        END IF
        !
        ! Check time and increment count
        !
        IS_LIMIT_TRIPPED = CHECK_LIMIT(TIME, INCR)
        !
        ! Output computed quantities
        !
        IF (GET_STEP_COMPLETE())  THEN
            !
            IF (MYID .EQ. 0) THEN
                !
                WRITE(DFLT_U,'(/,A,I0)') 'Info   :     > Converged on &
                    &increment ', INCR
                WRITE(DFLT_U,'(A,3(E14.6))') 'Info   :     > Loads on X0: ', &
                    & LOAD(1,:)
                WRITE(DFLT_U,'(A,3(E14.6))') 'Info   :       Loads on X1: ', &
                    & LOAD(2,:)
                WRITE(DFLT_U,'(A,3(E14.6))') 'Info   :       Loads on Y0: ', &
                    & LOAD(3,:)
                WRITE(DFLT_U,'(A,3(E14.6))') 'Info   :       Loads on Y1: ', &
                    & LOAD(4,:)
                WRITE(DFLT_U,'(A,3(E14.6))') 'Info   :       Loads on Z0: ', &
                    & LOAD(5,:)
                WRITE(DFLT_U,'(A,3(E14.6))') 'Info   :       Loads on Z1: ', &
                    & LOAD(6,:)
                !
            END IF
            !
            IF (GET_PRINT_FLAG() .EQ. 0) THEN
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
                    CALL WRITE_RESTART_FIELD(VELOCITY, C0_ANGS, C_ANGS, RSTAR, &
                        & RSTAR_N, WTS, CRSS, CRSS_N, E_ELAS_KK_BAR, SIG_VEC_N,&
                        & EQSTRAIN, EQPLSTRAIN, GAMMA, EL_WORK_N, EL_WORKP_N, &
                        & EL_WORK_RATE_N, EL_WORKP_RATE_N, PLSTRAIN, &
                        & TOTSTRAIN, ISTEP)
                    !
                    IF (OPTIONS%DEF_CONTROL_BY .EQ. UNIAXIAL_STRAIN_TARGET) THEN
                        !
                        IF (ISTEP .EQ. 1) THEN
                            !
                            PREV_STRAIN = 0.0D0
                            !
                        ELSE
                            !
                            PREV_STRAIN = &
                                & UNIAXIAL_OPTIONS%TARGET_STRAIN(ISTEP - 1,1)
                            !
                        END IF
                        !
                        CURR_STRAIN = UNIAXIAL_OPTIONS%TARGET_STRAIN(ISTEP, 1)
                        !
                    ELSE IF (OPTIONS%DEF_CONTROL_BY .EQ. &
                        & UNIAXIAL_LOAD_TARGET) THEN
                        !
                        PREV_STRAIN = 0.0D0
                        CURR_STRAIN = 0.0D0
                        !
                    END IF
                    !
                    CALL WRITE_UNIAXIAL_RESTART(INCR, TIME, LOAD, AREA, AREA0, &
                        & CURRENT_STEP, PREVIOUS_LOAD, STEP_COMPLETE, &
                        & DTIME_OLD, PREV_STRAIN, CURR_STRAIN, ISTEP)
                    !
                END IF
                !
                IF (FIBER_AVERAGE_OPTIONS%RUN_FIBER_AVERAGE) THEN
                    !
                    CALL RUN_FIBER_AVERAGE(ISTEP, C_ANGS, ELAS_TOT6, DPEFF, &
                        & CRSS, ELVOL)
                    !
                END IF
                !
            END IF
            !
            IF (IS_NECKING) THEN
                !
                CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP, PRINT_FLAG)
                !
                CALL PAR_QUIT('Error  :     > Specimen is necking.')
                !
            END IF
            !
            IF (IS_LIMIT_TRIPPED) THEN
                !
                CALL WRITE_REPORT_FILE_COMPLETE_STEPS(ISTEP, PRINT_FLAG)
                !
                CALL PAR_QUIT('Error  :     > Maximum time or maximum &
                    &increments exceeded.')
                !
            END IF
            !
            IF (GET_ALL_STEPS_COMPLETE()) THEN
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
        END IF
        !
    END DO TIME_STEPPING
    !
    RETURN
    !
    END SUBROUTINE DRIVER_UNIAXIAL_CONTROL
    !
    !===========================================================================
    !
    SUBROUTINE PROCESS_CTRL_DATA(VELOCITY)
    !
    ! Read data for load control and modify initial velocity field.
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
    ! DIFF_STRAIN: Strain target difFErential calculated between steps.
    ! DIFF_STRESS: Load target difFErential calculated between steps.
    ! USER_STRAIN_RATE: User-defined initial strain rate from BCs input.
    ! STEP_STRAIN_RATE: Array used to store strain rate jump values.
    ! IS_SIGNED: Array to enable load reversal checks to sign VEL_FACTOR.
    !
    INTEGER :: I, II
    REAL(RK) :: DIFF_STRAIN, DIFF_STRESS, USER_STRAIN_RATE
    !
    ! Notes:
    ! User-defined input (strain rates) should ALWAYS be positive. Compression
    ! and similar signed deformation control are handled internally by way of
    ! the sign in the targeted load or strain.
    !
    !---------------------------------------------------------------------------
    !
    SELECT CASE (OPTIONS%DEF_CONTROL_BY)
    !
    CASE (UNIAXIAL_LOAD_TARGET) ! Begin case (UNIAXIAL_LOAD_TARGET)
        !
        ! Initialize definition variables.
        !
        NSTEPS = UNIAXIAL_OPTIONS%NUMBER_OF_LOAD_STEPS
        TARGET_SURFACE = BCS_OPTIONS%LOADING_FACE
        USER_STRAIN_RATE = BCS_OPTIONS%STRAIN_RATE
        DIFF_STRESS = 0.0D0
        !
        ! Allocate the arrays for the load target sequence.
        !
        ALLOCATE(TARGET_LOAD(NSTEPS))
        ALLOCATE(TARGET_INDEX(NSTEPS))
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
                CALL PAR_QUIT('Error  :     > Input &
                    &`number_of_strain_rate_ jumps` does not match number of &
                    &defined `strain_rate_jumps`.')
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
            IF (UNIAXIAL_OPTIONS%TARGET_LOAD(I,4) .EQ. -1) THEN
                !
                CALL PAR_QUIT('Error  :     > Number of load steps does not &
                    &match defined number of targets in *.config file.')
                !
            END IF
            !
            ! Handle the first step alone to avoid out of bounds on arrays.
            !
            IF (I .EQ. 1) THEN
                !
                ! First, assign the fixed values
                !
                TARGET_LOAD(I)   = UNIAXIAL_OPTIONS%TARGET_LOAD(I,1)
                TARGET_INDEX(I)  = BCS_OPTIONS%LOADING_DIRECTION + 1
                TIME_STEP(I) = UNIAXIAL_OPTIONS%TARGET_LOAD(I,2)
                TIME_STEP_MIN(I) = UNIAXIAL_OPTIONS%TARGET_LOAD(I,3)
                PRINT_FLAG(I) = INT(UNIAXIAL_OPTIONS%TARGET_LOAD(I,4))
                !
                ! Next, determine the velocity factor for the first step.
                !
                VEL_FACTOR(I) = ABS(STEP_STRAIN_RATE(I) / USER_STRAIN_RATE)
                !
                ! Check if the vel factor needs to be signed for compression.
                ! Calculate the differential stress from zero.
                !
                DIFF_STRESS = (UNIAXIAL_OPTIONS%TARGET_LOAD(I,1))
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
                TARGET_LOAD(I)  = UNIAXIAL_OPTIONS%TARGET_LOAD(I,1)
                TARGET_INDEX(I) = BCS_OPTIONS%LOADING_DIRECTION + 1
                TIME_STEP(I) = UNIAXIAL_OPTIONS%TARGET_LOAD(I,2)
                TIME_STEP_MIN(I) = UNIAXIAL_OPTIONS%TARGET_LOAD(I,3)
                PRINT_FLAG(I) = INT(UNIAXIAL_OPTIONS%TARGET_LOAD(I,4))
                !
                ! Next, determine the velocity factor for the first step.
                !
                VEL_FACTOR(I) = STEP_STRAIN_RATE(I) / STEP_STRAIN_RATE(I-1)
                !
                ! Check if the factor needs to be signed for compression.
                ! Calculate the differential stress from previous.
                !
                DIFF_STRESS = (UNIAXIAL_OPTIONS%TARGET_LOAD(I,1)) &
                    & - (UNIAXIAL_OPTIONS%TARGET_LOAD(I-1,1))
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
    ! End case (UNIAXIAL_LOAD_TARGET)
    !
    CASE (UNIAXIAL_STRAIN_TARGET) ! Begin case (UNIAXIAL_STRAIN_TARGET)
        !
        ! Initialize definition variables.
        !
        NSTEPS = UNIAXIAL_OPTIONS%NUMBER_OF_STRAIN_STEPS
        TARGET_SURFACE = BCS_OPTIONS%LOADING_FACE
        USER_STRAIN_RATE = BCS_OPTIONS%STRAIN_RATE
        DIFF_STRAIN = 0.0D0
        !
        ! Allocate the arrays for the strain target sequence.
        !
        ALLOCATE(TARGET_INCR(NSTEPS))
        ALLOCATE(TARGET_INDEX(NSTEPS))
        ALLOCATE(TARGET_SIGN(NSTEPS))
        ALLOCATE(VEL_FACTOR(NSTEPS))
        ALLOCATE(TIME_STEP(NSTEPS))
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
                CALL PAR_QUIT('Error  :     > Input &
                    &`number_of_strain_rate_jumps` does not match number of &
                    &defined `strain_rate_jumps`.')
                !
            END IF
            !
        END IF
        !
        ! Process all strain steps and build control arrays.
        !
        DO I = 1, NSTEPS
            !
            ! Initial error handling.
            !
            IF (UNIAXIAL_OPTIONS%TARGET_STRAIN(I,3) .EQ. -1) THEN
                !
                CALL PAR_QUIT('Error  :     > Number of strain steps does not&
                    & match defined number of targets in *.config file.')
                !
            END IF
            !
            ! Handle the first step alone to avoid out of bounds on arrays.
            !
            IF (I .EQ. 1) THEN
                !
                ! First, assign the fixed values.
                !
                TARGET_INCR(I) = INT(UNIAXIAL_OPTIONS%TARGET_STRAIN(I,2))
                TARGET_INDEX(I) = BCS_OPTIONS%LOADING_DIRECTION + 1
                TARGET_SIGN(I) = 1
                PRINT_FLAG(I) = INT(UNIAXIAL_OPTIONS%TARGET_STRAIN(I,3))
                !
                ! Next, determine the velocity factor for the first step.
                !
                VEL_FACTOR(I) = ABS(STEP_STRAIN_RATE(I) / USER_STRAIN_RATE)
                !
                ! Check if the vel factor needs to be signed for compression.
                ! Calculate the differential strain from zero.
                !
                DIFF_STRAIN = (UNIAXIAL_OPTIONS%TARGET_STRAIN(I,1))
                !
                IF (DIFF_STRAIN .GT. 0.0D0) THEN
                    !
                    IS_SIGNED(I) = .FALSE.
                    !
                ELSE IF (DIFF_STRAIN .LT. 0.0D0) THEN
                    !
                    IS_SIGNED(I) = .TRUE.
                    VELOCITY = VELOCITY * (-1) ! Sign the entire velocity field.
                    !
                ELSE ! The differential change is zero, thus dt is zero.
                    !
                    CALL PAR_QUIT('Error  :     > Strain differential&
                        & is zero between steps.')
                    !
                END IF
                !
                ! Then, calculate the time differential for the first step.
                !
                TIME_STEP(I) = (ABS(DIFF_STRAIN)) / &
                    & (STEP_STRAIN_RATE(I) * TARGET_INCR(I))
                !
            ELSE IF (I .GT. 1) THEN
                !
                ! First, assign the fixed values for I > 1.
                !
                TARGET_INCR(I) = INT(UNIAXIAL_OPTIONS%TARGET_STRAIN(I,2)) + &
                    & TARGET_INCR(I-1)
                !
                TARGET_INDEX(I) = BCS_OPTIONS%LOADING_DIRECTION + 1
                TARGET_SIGN(I) = 1
                PRINT_FLAG(I) = INT(UNIAXIAL_OPTIONS%TARGET_STRAIN(I,3))
                !
                ! Next, determine the velocity factor for the first step.
                !
                VEL_FACTOR(I) = STEP_STRAIN_RATE(I) / STEP_STRAIN_RATE(I-1)
                !
                ! Check if the factor needs to be signed for compression.
                ! Calculate the differential strain from zero.
                !
                DIFF_STRAIN = (UNIAXIAL_OPTIONS%TARGET_STRAIN(I,1)) &
                    & - (UNIAXIAL_OPTIONS%TARGET_STRAIN(I-1,1))
                !
                IF (DIFF_STRAIN .GT. 0.0D0) THEN
                    !
                    IS_SIGNED(I) = .FALSE.
                    !
                ELSE IF (DIFF_STRAIN .LT. 0.0D0) THEN
                    !
                    IS_SIGNED(I) = .TRUE.
                    !
                ELSE ! The differential change is zero, thus dt is zero.
                    !
                    CALL PAR_QUIT('Error  :     > Strain differential&
                        & is zero between steps.')
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
                ! Then, calculate the time differential for the first step
                !
                TIME_STEP(I) = (ABS(DIFF_STRAIN)) / &
                    & (STEP_STRAIN_RATE(I)*UNIAXIAL_OPTIONS%TARGET_STRAIN(I,2))
                !
            END IF
            !
        END DO
        !
        ! Check that increments are monotonically increasing.
        !
        DO I = 2, NSTEPS
            !
            IF (TARGET_INCR(I) .LE. TARGET_INCR(I-1)) THEN
                !
                CALL PAR_QUIT('Error  :     > Increments are not&
                    & monotonically increasing.')
                !
            END IF
            !
        END DO
        !
    ! End case (UNIAXIAL_STRAIN_TARGET)
    !
    CASE DEFAULT
        !
        CALL PAR_QUIT('Error  :     > Invalid control option defined.')
        !
    END SELECT
    !
    ! Define initial values for driver.
    !
    CURRENT_STEP = 1
    STEP_COMPLETE = .FALSE.
    ALL_STEPS_COMPLETE = .FALSE.
    DTIME_OLD = TIME_STEP(1)
    !
    RETURN
    !
    END SUBROUTINE PROCESS_CTRL_DATA
    !
    !===========================================================================
    !
    SUBROUTINE READ_UNIAXIAL_RESTART(INCR, TIME, LOAD, PREVIOUS_LOAD, AREA, &
        & AREA0, VELOCITY, TIME_STEP)
    !
    ! Read uniaxial control restart information.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! INCR: Current step increment (additive).
    ! TIME: Current step total time.
    ! LOAD: Current load on all surfaces in the mesh.
    ! PREVIOUS_LOAD: Previous step load value.
    ! AREA/AREA0: Current and initial surface areas of the mesh at current step.
    !
    INTEGER,  INTENT(OUT) :: INCR
    REAL(RK), INTENT(OUT) :: TIME
    REAL(RK), INTENT(OUT) :: LOAD(NSURFACES,3)
    REAL(RK), INTENT(OUT) :: PREVIOUS_LOAD(3)
    REAL(RK), INTENT(OUT) :: AREA(NSURFACES)
    REAL(RK), INTENT(OUT) :: AREA0(NSURFACES)
    REAL(RK), INTENT(INOUT) :: VELOCITY(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(INOUT) :: TIME_STEP(NSTEPS)
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
    REAL(RK) :: DIFF_STRAIN
    REAL(RK) :: DIFF_STRESS
    INTEGER :: I
    !
    ! Notes:
    ! Current_load becomes the previous_load inside subroutine time_increment.
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
    READ(MYUNIT) CURRENT_STEP
    READ(MYUNIT) PREVIOUS_LOAD
    READ(MYUNIT) STEP_COMPLETE
    READ(MYUNIT) DTIME_OLD
    READ(MYUNIT) INCR
    READ(MYUNIT) TIME
    !
    DO ISURF = 1, NSURFACES
        !
        READ(MYUNIT) LOAD(ISURF,:)
        !
    ENDDO
    !
    READ(MYUNIT) AREA
    READ(MYUNIT) AREA0
    READ(MYUNIT) PREV_STRAIN
    READ(MYUNIT) CURR_STRAIN
    !
    ! Re-assign variables for new restart
    !
    STEP_COMPLETE = .TRUE.
    CURRENT_LOAD(1) = LOAD(2,1)
    CURRENT_LOAD(2) = LOAD(4,2)
    CURRENT_LOAD(3) = LOAD(6,3)
    !
    ! Print statistics to console
    !
    IF (MYID .EQ. 0) THEN
        !
        WRITE(DFLT_U,'(A)') 'Info   : Reading restart control information...'
        WRITE(DFLT_U,'(A)') 'Info   :   - Restart parameters:'
        WRITE(DFLT_U,'(A, I0)')       'Info   :     > Prior Increments: ', &
            & INCR
        WRITE(DFLT_U,'(A, I0)')       'Info   :     > Prior Steps:      ', &
            & CURRENT_STEP - 1
        WRITE(DFLT_U,'(A, E14.4)')    'Info   :     > Prior Time:   ', TIME
        WRITE(DFLT_U,'(A, 3(E14.4))') 'Info   :     > Current Load: ', &
            & CURRENT_LOAD
        WRITE(DFLT_U,'(A, 3(E14.4))') 'Info   :     > Previous Load:', &
            & PREVIOUS_LOAD
        !
    ENDIF
    !
    ! Reinitialize step and time if new files
    ! In the future, we will have to handle the APPEND case, which will continue
    ! sequentially, rather than reinitializing
    !
    IF (OPTIONS%RESTART_FILE_HANDLING .EQ. 1) THEN ! New files
        !
        CURRENT_STEP = 1
        STEP_COMPLETE = .FALSE.
        INCR = 0
        TIME = 0.0D0
        !
    !ELSE IF (OPTIONS%RESTART_FILE_HANDLING .EQ. 0) THEN ! Append
    END IF
    !
    ! Reevaluate time step for first step if strain targeting
    !
    IF (OPTIONS%DEF_CONTROL_BY .EQ. UNIAXIAL_STRAIN_TARGET) THEN
        !
        DIFF_STRAIN = (UNIAXIAL_OPTIONS%TARGET_STRAIN(1,1)) - (CURR_STRAIN)
        TIME_STEP(1) = (ABS(DIFF_STRAIN)) / (STEP_STRAIN_RATE(1) * &
            & TARGET_INCR(1))
        !
    END IF
    !
    ! Check for a reversal in loading direction upon restart, correct loading
    ! history if so. The loading history as assigned in PROCESS_CTRL_DATA
    ! assumes that we are starting from 0 load. Below rectifies any directional
    ! issues with loading due to this assumption.
    !
    IF ((OPTIONS%DEF_CONTROL_BY .EQ. UNIAXIAL_LOAD_TARGET) .AND. &
        & (((CURRENT_LOAD(BCS_OPTIONS%LOADING_DIRECTION + 1) - &
        & PREVIOUS_LOAD(BCS_OPTIONS%LOADING_DIRECTION + 1)) * &
        & (TARGET_LOAD(1) -  &
        & CURRENT_LOAD(BCS_OPTIONS%LOADING_DIRECTION + 1))) .LT. 0.0D0)) THEN
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
            DIFF_STRESS = UNIAXIAL_OPTIONS%TARGET_LOAD(1, 1) - &
                & CURRENT_LOAD(BCS_OPTIONS%LOADING_DIRECTION + 1)
            !
            IF (DIFF_STRESS .GT. 0.0D0) THEN
                !
                IS_SIGNED(1) = .FALSE.
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
                DIFF_STRESS = (UNIAXIAL_OPTIONS%TARGET_LOAD(I, 1)) &
                    & - (UNIAXIAL_OPTIONS%TARGET_LOAD(I - 1, 1))
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
    ELSE IF ((OPTIONS%DEF_CONTROL_BY .EQ. UNIAXIAL_STRAIN_TARGET) .AND. &
        & (((CURR_STRAIN - PREV_STRAIN) * &
        & (UNIAXIAL_OPTIONS%TARGET_STRAIN(1, 1) - CURR_STRAIN)) &
        & .LT. 0.0D0)) THEN
        !
        ! Change velocity direction
        !
        VELOCITY = -1.0D0 * VELOCITY
        !
        ! And check for any load reversals on subsequent steps
        !
        IF (NSTEPS .GE. 2) THEN
            !
            DIFF_STRAIN = UNIAXIAL_OPTIONS%TARGET_STRAIN(1, 1) - CURR_STRAIN
            !
            IF (DIFF_STRAIN .GT. 0.0D0) THEN
                !
                IS_SIGNED(1) = .FALSE.
                !
            ELSE IF (DIFF_STRAIN .LT. 0.0D0) THEN
                !
                IS_SIGNED(1) = .TRUE.
                !
            ELSE ! The differential change is zero, thus dt is zero.
                !
                CALL PAR_QUIT('Error  :     > Strain differential&
                    & is zero between steps.')
                !
            END IF
            !
            ! Reset VEL_FACTOR for reassignment in the below loops
            !
            VEL_FACTOR = ABS(VEL_FACTOR)
            !
            DO I = 2, NSTEPS
                !
                DIFF_STRAIN = (UNIAXIAL_OPTIONS%TARGET_STRAIN(I, 1)) &
                    & - (UNIAXIAL_OPTIONS%TARGET_STRAIN(I - 1, 1))
                !
                IF (DIFF_STRAIN .GT. 0.0D0) THEN
                    !
                    IS_SIGNED(I) = .FALSE.
                    !
                ELSE IF (DIFF_STRAIN .LT. 0.0D0) THEN
                    !
                    IS_SIGNED(I) = .TRUE.
                    !
                ELSE ! The differential change is zero, thus dt is zero.
                    !
                    CALL PAR_QUIT('Error  :     > Load differential is zero&
                        & between steps.')
                    !
                END IF
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
    CLOSE(MYUNIT)
    !
    RETURN
    !
    END SUBROUTINE READ_UNIAXIAL_RESTART
    !
    !===========================================================================
    !
    FUNCTION CHECK_LIMIT(TIME, INCR) RESULT(IS_LIMIT_TRIPPED)
    !
    ! Check if specimen time or increment limit is tripped.
    !
    !---------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! Arguments:
    ! TIME: Current step total time.
    ! INCR: Current step increment (additive).
    ! IS_LIMIT_TRIPPED: Flag for checking if simulation controls are tripped.
    !
    REAL(RK) :: TIME
    INTEGER  :: INCR
    LOGICAL  :: IS_LIMIT_TRIPPED
    !
    !---------------------------------------------------------------------------
    !
    IS_LIMIT_TRIPPED = .FALSE.
    !
    IF ((TIME .GE. OPTIONS%MAX_TOTAL_TIME) .OR. &
        & (INCR .GE. OPTIONS%MAX_INCR)) THEN
        !
        IS_LIMIT_TRIPPED = .TRUE.
        !
        IF (.NOT. STEP_COMPLETE) THEN
            !
            ! Update step because GET_PRINT_FLAG checks flag from previous step.
            !
            CURRENT_STEP = CURRENT_STEP + 1
            !
        ENDIF
        !
        STEP_COMPLETE = .TRUE.
        ALL_STEPS_COMPLETE = .TRUE.
        !
    ENDIF
    !
    RETURN
    !
    !---------------------------------------------------------------------------
    !
    END FUNCTION CHECK_LIMIT
    !
    !===========================================================================
    !
    FUNCTION CHECK_NECKING() RESULT(IS_NECKING)
    !
    ! Check if specimen is necking (that is, the load on the target surface is
    ! decreasing from the previous step).
    !
    !---------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! Arguments:
    ! IS_NECKING: Flag for checking if the specimen has necked.
    !
    LOGICAL :: IS_NECKING
    !
    ! Locals:
    ! INDEX: Current step loading direction.
    !
    INTEGER :: INDEX
    !
    !---------------------------------------------------------------------------
    !
    IS_NECKING = .FALSE.
    !
    INDEX = TARGET_INDEX(CURRENT_STEP)
    !
    IF ((.NOT. STEP_COMPLETE) .AND. &
        & ((CURRENT_LOAD(INDEX)-PREVIOUS_LOAD(INDEX))/TARGET_SIGN(CURRENT_STEP)&
        & .LT. 0.0)) THEN
        !
        IS_NECKING = .TRUE.
        STEP_COMPLETE = .TRUE.
        ALL_STEPS_COMPLETE = .TRUE.
        !
        ! Update step because GET_PRINT_FLAG checks flag from previous step.
        !
        CURRENT_STEP = CURRENT_STEP + 1
        !
    ENDIF
    !
    RETURN
    !
    !---------------------------------------------------------------------------
    !
    END FUNCTION CHECK_NECKING
    !
    !===========================================================================
    !
    LOGICAL FUNCTION GET_ALL_STEPS_COMPLETE()
    !
    ! Check if all steps complete.
    !
    !---------------------------------------------------------------------------
    !
    GET_ALL_STEPS_COMPLETE = ALL_STEPS_COMPLETE
    !
    END FUNCTION GET_ALL_STEPS_COMPLETE
    !
    !===========================================================================
    !
    INTEGER FUNCTION GET_CURRENT_STEP()
    !
    ! Returns current target index.
    !
    !---------------------------------------------------------------------------
    !
    GET_CURRENT_STEP = CURRENT_STEP
    !
    END FUNCTION GET_CURRENT_STEP
    !
    !===========================================================================
    !
    REAL(RK) FUNCTION GET_INITIAL_DTIME()
    !
    ! First time step (avoid CALL to TIME_INCREMENT).
    !
    !--------------------------------------------------------------------------
    !
    GET_INITIAL_DTIME = TIME_STEP(1)
    !
    END FUNCTION GET_INITIAL_DTIME
    !
    !===========================================================================
    !
    INTEGER FUNCTION GET_PRINT_FLAG()
    !
    ! Returns current print flag.
    !
    ! Notes:
    ! The current_step gets incremented in TIME_INCREMENT when step is complete.
    !
    !---------------------------------------------------------------------------
    !
    GET_PRINT_FLAG = PRINT_FLAG(CURRENT_STEP - 1)
    !
    END FUNCTION GET_PRINT_FLAG
    !
    !===========================================================================
    !
    LOGICAL FUNCTION GET_STEP_COMPLETE()
    !
    ! Set the current step as complete.
    !
    !---------------------------------------------------------------------------
    !
    GET_STEP_COMPLETE = STEP_COMPLETE
    !
    END FUNCTION GET_STEP_COMPLETE
    !
    !===========================================================================
    !
    INTEGER FUNCTION GET_TARGET_INDEX()
    !
    ! Returns current target index.
    !
    !---------------------------------------------------------------------------
    !
    GET_TARGET_INDEX = TARGET_INDEX(CURRENT_STEP)
    !
    END FUNCTION GET_TARGET_INDEX
    !
    !===========================================================================
    !
    INTEGER FUNCTION GET_TARGET_SURFACE()
    !
    ! Returns current target surface.
    !
    !---------------------------------------------------------------------------
    !
    GET_TARGET_SURFACE = TARGET_SURFACE
    !
    END FUNCTION GET_TARGET_SURFACE
    !
    !===========================================================================
    !
    REAL(RK) FUNCTION GET_VEL_FACTOR()
    !
    ! Check if load reverses direction.
    !
    !---------------------------------------------------------------------------
    !
    GET_VEL_FACTOR = VEL_FACTOR(CURRENT_STEP)
    !
    END FUNCTION GET_VEL_FACTOR
    !
    !===========================================================================
    !
    FUNCTION IN_RANGE(TRIAL_LOAD) RESULT(STATUS)
    !
    ! Check if load is in specified range.
    !
    !---------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! Arguments:
    ! TRIAL_LOAD: Estimated load value to determine if dt needs to be changed.
    ! STATUS: Flag for checking if the trial load is within the range.
    !
    REAL(RK), DIMENSION(3) :: TRIAL_LOAD
    LOGICAL :: STATUS
    !
    ! Locals:
    ! INDEX: Current step loading direction.
    ! DIFFERENCE: DifFErence between the trial load and the target load on step.
    ! SIGN: +/- value for handling load INCReases/decreases.
    !
    INTEGER :: INDEX
    REAL(RK) :: DIFFERENCE, SIGN
    !
    !---------------------------------------------------------------------------
    !
    STATUS = .FALSE.
    !
    INDEX = TARGET_INDEX(CURRENT_STEP)
    SIGN  = TARGET_SIGN(CURRENT_STEP)
    !
    DIFFERENCE = TRIAL_LOAD(INDEX) - TARGET_LOAD(CURRENT_STEP)
    DIFFERENCE = DIFFERENCE*SIGN + OPTIONS%LOAD_TOL
    !
    IF (DIFFERENCE > 0.0D0) THEN
        !
        STATUS = .TRUE.
        !
    ENDIF
    !
    RETURN
    !
    !---------------------------------------------------------------------------
    !
    END FUNCTION IN_RANGE
    !
    !===========================================================================
    !
    FUNCTION TIME_INCREMENT(LOAD, INCR) RESULT(DELTA_T)
    !
    ! Return time increment according to loading step.
    !
    !---------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! Arguments:
    ! INCR: Current step INCRement (additive).
    ! DELTA_T: TIME INCRement for load step.
    ! LOAD: Current step load value.
    !
    INTEGER  :: INCR
    REAL(RK) :: DELTA_T
    REAL(RK), DIMENSION(3) :: LOAD
    !
    ! Locals:
    ! INDEX: Generic loop index value.
    ! DLOAD: DifFErential load between current and previous TIME steps.
    ! TRIAL_LOAD: Estimated load value to determine if dt needs to be changed.
    ! TRIAL_DT: Estimated new TIME step value for the trial load.
    ! MIN_TIME_STEP: Current step minimum defined TIME INCRement.
    ! DELTA_LOAD: DifFErence between the trial load and current load at step.
    !
    INTEGER :: INDEX
    REAL(RK) :: DLOAD(3), TRIAL_LOAD(3), TRIAL_DT, MIN_TIME_STEP, DELTA_LOAD
    !
    ! Notes:
    ! Used for load control only!
    ! The suggested TIME step is given in the input data.
    ! Near the target load, we use a linear approximation to match the load,
    ! with a minimum value to insure we get there.
    !
    !---------------------------------------------------------------------------
    !
    DELTA_T = TIME_STEP(CURRENT_STEP)
    !
    PREVIOUS_LOAD = CURRENT_LOAD
    CURRENT_LOAD = LOAD
    !
    ! First, check if we are at target load or INCRement.
    !
    SELECT CASE (OPTIONS%DEF_CONTROL_BY)
        !
        CASE (UNIAXIAL_LOAD_TARGET)
            !
            STEP_COMPLETE = IN_RANGE(CURRENT_LOAD)
            !
        CASE (UNIAXIAL_STRAIN_TARGET)
            !
            STEP_COMPLETE = .FALSE.
            !
            IF (INCR .EQ. TARGET_INCR(CURRENT_STEP)) THEN
                !
                STEP_COMPLETE = .TRUE.
                !
            END IF
            !
        CASE DEFAULT
            !
            CALL PAR_QUIT('Error  :     > Invalid control option.')
            !
    END SELECT
    !
    IF (STEP_COMPLETE) THEN
        !
        CURRENT_STEP = CURRENT_STEP + 1
        !
        IF (CURRENT_STEP > NSTEPS) THEN
            !
            ALL_STEPS_COMPLETE = .TRUE.
            !
            RETURN
            !
        ENDIF
        !
        DELTA_T = TIME_STEP(CURRENT_STEP)
        !
    ENDIF
    !
    ! Now, check to see if we are close to target, and need to adjust the time
    ! step.
    !
    IF (OPTIONS%DEF_CONTROL_BY .EQ. UNIAXIAL_LOAD_TARGET) THEN
        !
        DLOAD = (CURRENT_LOAD - PREVIOUS_LOAD)/DTIME_OLD
        TRIAL_LOAD = CURRENT_LOAD + DLOAD*DELTA_T
        !
        IF (IN_RANGE(TRIAL_LOAD)) THEN
            !
            INDEX = TARGET_INDEX(CURRENT_STEP)
            !
            !LEGACY NOTE: -TM. FROM DONALD'S VERSION:
            DELTA_LOAD = TRIAL_LOAD (INDEX) - CURRENT_LOAD(INDEX)
            !
            IF (ABS(DELTA_LOAD) < 1.0D-16) THEN
                !
                DTIME_OLD = DELTA_T;
                !
                RETURN
                !
            ENDIF
            !
            TRIAL_DT = (TARGET_LOAD(CURRENT_STEP) - CURRENT_LOAD(INDEX))/&
                & (TRIAL_LOAD (INDEX) - CURRENT_LOAD(INDEX))
            !
            TRIAL_DT = TRIAL_DT * DELTA_T
            DELTA_T  = TRIAL_DT * OPTIONS%DTIME_FACTOR
            !
            IF (DELTA_T > TIME_STEP(CURRENT_STEP)) THEN
                !
                DELTA_T = TIME_STEP(CURRENT_STEP)
                !
            ENDIF
            !
            MIN_TIME_STEP = TIME_STEP_MIN(CURRENT_STEP)
            !
            IF ((DELTA_T) < MIN_TIME_STEP) THEN
                !
                DELTA_T = MIN_TIME_STEP
                !
            ENDIF
            !
        ENDIF
        !
    ENDIF
    !
    DTIME_OLD = DELTA_T
    !
    RETURN
    !
    END FUNCTION TIME_INCREMENT
    !
END MODULE DRIVER_UNIAXIAL_CONTROL_MOD
