! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE DRIVER_ISO_VP_MOD
!
! Driver for the isotropic viscoplastic solution.
!
! Contains subroutines:
! DRIVER_VP_SOLVE: Solves isotropic viscoplastic solution.
!
! From libf95:
!
USE LIBF95, ONLY: RK=>REAL_KIND
!
! From libfepx:
!
USE DIMENSIONS_MOD
USE ITERATE_STRESS_VP_MOD
USE MATRIX_OPERATIONS_MOD, ONLY: CALC_ELVOL
USE READ_INPUT_MOD
USE UNITS_MOD
USE WRITE_OUTPUT_MOD, ONLY: PRINT_STEP
!
! From libparallel:
!
USE GATHER_SCATTER_MOD
USE PARALLEL_MOD
!
IMPLICIT NONE
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: DRIVER_VP_SOLVE
!
CONTAINS
    !
    SUBROUTINE DRIVER_VP_SOLVE(ITYPE, BCS, VELOCITY, PFORCE, DOF_TRACE, CRSS, &
        & C0_ANGS)
    !
    ! Driver for viscoplastic solution.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! ITYPE: 0=isotropic viscoplastic,1=anisotropic viscoplastic
    ! BCS: Global boundary condition vector
    ! VELOCITY: Global velocity vector
    ! PFORCE: Global force vector
    ! DOF_TRACE: Gather/scatter trace for degrees of freedom
    ! CRSS: Elemental critical resolved shear stress
    ! C0_ANGS: Initial orientation matrices
    !
    INTEGER :: ITYPE
    LOGICAL :: BCS(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: VELOCITY(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: PFORCE(DOF_SUB1:DOF_SUP1)
    TYPE(TRACE) :: DOF_TRACE
    REAL(RK) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    INTEGER :: INCR
    INTEGER :: IER
    INTEGER :: M_EL
    REAL(RK) :: DTIME
    REAL(RK) :: ELPRESS(EL_SUB1:EL_SUP1)
    REAL(RK) :: EVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: EPSEFF(EL_SUB1:EL_SUP1)
    ! tm268_M13: REAL(RK) :: eqplas(EL_SUB1:EL_SUP1)
    ! SIG_DUMMY is currently just for passing to print_incr, but in the future
    !   one could store and print the grain stresses.
    !
    ! mpk: Are the following 4 variables necessary??
    !
    !REAL(RK) :: SIG_AVG(0:TVEC1, EL_SUB1:EL_SUP1)
    !REAL(RK) :: SIG_AVG_KK(EL_SUB1:EL_SUP1)
    !REAL(RK) :: SIG_DUMMY(0:TVEC, EL_SUB1:EL_SUP1)
    !REAL(RK) :: SDUMMY3X3(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    !
    REAL(RK), ALLOCATABLE :: ELVOL_0(:)
    !
    !----------------------------------------------------------------------
    ! 
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    IF (ITYPE .EQ. 1) CALL PAR_QUIT(&
        &'Error  :     > ANISOTROPIC_VP is no longer implemented.')
    !
    ! Isotropic viscoplastic solution
    !
    IF (MYID .EQ. 0)  WRITE(DFLT_U,'(A)') 'Info   :   - Initializing fields &
        &from isotropic viscoplastic solution'
    !
    ELPRESS = 0.0D0
    EVEL = 0.0D0
    DTIME = 0.0D0
    INCR = 0
    !
    IER = ITMETHOD_VP(ITYPE, BCS, PFORCE, VELOCITY, ELPRESS, EVEL, DOF_TRACE, &
        & EPSEFF)
    !      
    IF (IER .LT. 0) THEN
        !
        CALL PAR_QUIT('Error  :     > ITMETHOD_VP failed to converge.')
        !
    END IF
    !
    ! Output computed quantities
    !
    CALL PART_GATHER(ELEMENT_CRDS, COORDS, NODES, DOF_TRACE)
    !
    IF (PRINT_OPTIONS%PRINT_ELVOL) THEN
        !
        ALLOCATE(ELVOL_0(EL_SUB1:EL_SUP1))
        !
        CALL CALC_ELVOL(ELVOL_0, ELEMENT_CRDS)
        !
        CALL PRINT_STEP(INCR, ITYPE, COORDS, VELOCITY, C0_ANGS, CRSS, &
            & ELVOL_0)
        !
        DEALLOCATE(ELVOL_0)
        !
    ELSE
        !
        CALL PRINT_STEP(INCR, ITYPE, COORDS, VELOCITY, C0_ANGS, CRSS)
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE DRIVER_VP_SOLVE
    !
END MODULE DRIVER_ISO_VP_MOD
