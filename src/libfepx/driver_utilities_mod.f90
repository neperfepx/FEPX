! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE DRIVER_UTILITIES_MOD
!
! This module defines common subroutines called by the various driver modules.
!
! Contains subroutines:
! CALC_MESH_DIM: Calculate mesh dimensions.
! CALC_STRESS_STRAIN: Calculate elemental stress/strain, mesh surface areas, 
!   and macroscopic loads.
! EST_AVG_MOD: Estimate bulk elastic moduli. Assumes a uniform texture and 
!   equal valued element volumes.
! SOLVE_LIN_SYS_3: Solve a linear system of three equations.
! TEMP_UPDATE_STATE_EVPS: Temporarily updates state at center of element only.
! UPDATE_STATE_EVPS: Update crystal states for entire mesh.
! VELOCITY_ITERATION: Perform iteration on velocity field.
!
USE IntrinsicTypesModule, RK=>REAL_KIND
!
USE DIMSMODULE
USE GATHER_SCATTER
USE ITMETHODEVPSMODULE
USE KINEMATICSMODULE
USE READ_INPUT_MOD
USE MICROSTRUCTURE_MOD
USE PARALLEL_MOD
USE POLYCRYSTALRESPONSEEVPSMODULE
USE POST_UPDATE_N_MOD
USE QUADRATURE_MOD
USE RSTARN_SOLVE_LAG_MOD
USE SHAPE_3D_MOD
USE STRESS_STRAIN_MOD
USE SURF_INFO_MOD
USE UNITS_MOD
USE UTILSCRYSTALMODULE
!
IMPLICIT NONE
!
PRIVATE
!
PUBLIC :: CALC_MESH_DIM, CALC_STRESS_STRAIN, EST_AVG_MOD, SOLVE_LIN_SYS_3, &
    & TEMP_UPDATE_STATE_EVPS, UPDATE_STATE_EVPS, VELOCITY_ITERATION
!
CONTAINS
    !
    SUBROUTINE CALC_MESH_DIM(LENGTH, INDX, INDY, INDZ)
    !
    ! Calculate mesh dimensions for triaxial loading.
    !      
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! LENGTH:
    ! INDX:
    ! INDY:
    ! INDZ:
    !
    REAL(RK), INTENT(OUT) :: LENGTH(3)
    INTEGER, INTENT(IN)   :: INDX(1:(DOF_SUP1-DOF_SUB1+1)/3)
    INTEGER, INTENT(IN)   :: INDY(1:(DOF_SUP1-DOF_SUB1+1)/3)
    INTEGER, INTENT(IN)   :: INDZ(1:(DOF_SUP1-DOF_SUB1+1)/3)
    !
    ! Locals:
    ! PART_XMAX,YMAX,ZMAX:
    ! LENGTH_X,Y,Z:
    !
    REAL(RK) :: PART_XMAX, PART_YMAX, PART_ZMAX
    REAL(RK) :: LENGTH_X, LENGTH_Y, LENGTH_Z
    !      
    !---------------------------------------------------------------------------
    !     
    ! Locate maximum coordinate values of the domain.
    !
    PART_XMAX = MAXVAL(COORDS(INDX))
    PART_YMAX = MAXVAL(COORDS(INDY))
    PART_ZMAX = MAXVAL(COORDS(INDZ))
    !
    CALL PAR_MAX(PART_XMAX, LENGTH_X)
    CALL PAR_MAX(PART_YMAX, LENGTH_Y)
    CALL PAR_MAX(PART_ZMAX, LENGTH_Z)
    !
    ! Store returned maximums into output array.
    !
    LENGTH(1) = LENGTH_X
    LENGTH(2) = LENGTH_Y
    LENGTH(3) = LENGTH_Z
    !
    RETURN
    !
    END SUBROUTINE CALC_MESH_DIM
    !
    !===========================================================================
    !
    SUBROUTINE CALC_STRESS_STRAIN(S_AVG_3X3, ELAS_TOT6, LOAD, AREA, SIG_VEC, &
        & SIG_KK, E_ELAS_KK, C_ANGS, KEINV, WTS)
    !
    ! Calculate elemental stress/strain, mesh surface areas, macroscopic loads.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! Needs to be defined - JC
    !
    REAL(RK), INTENT(OUT) :: S_AVG_3X3(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: ELAS_TOT6(0:TVEC , 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: LOAD(NSURFACES,3)
    REAL(RK), INTENT(OUT) :: AREA(NSURFACES)
    REAL(RK), INTENT(INOUT) :: SIG_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: SIG_KK(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: E_ELAS_KK(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: KEINV(0:TVEC1,1:NUMPHASES)
    REAL(RK), INTENT(IN) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    !
    ! Locals:
    ! Needs to be defined - JC
    !
    REAL(RK), POINTER :: E_ELAS_VEC_DEV_TMP(:,:,:) => NULL()
    REAL(RK)  ::  E_ELAS_VEC_DEV(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK)  ::  E_ELAS_KK_GRN(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK)  ::  S_KK_GRN(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK)  ::  ELAS_T3X3(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK)  ::  V_TENSOR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK)  ::  DETERM_V(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK)  ::  QR5X5(0:TVEC1, 0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK)  ::  SIG_SM(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK)  ::  S_AVG(0:TVEC1,  EL_SUB1:EL_SUP1)
    REAL(RK)  ::  S_AVG_KK(EL_SUB1:EL_SUP1)
    REAL(RK), ALLOCATABLE :: SIG_AVG_ALL(:)
    INTEGER :: MY_PHASE(0:(EL_SUP1-EL_SUB1))
    INTEGER :: M_EL, EL_DOF_MIN, EL_DOF_MAX, IDOF, I, J, K
    INTEGER :: IPHASE, NUMIND
    INTEGER, POINTER :: INDICES(:) => NULL()
    ! 
    ! Notes:
    ! 'MY_PHASE' here is 0-based because it is used only once for the local
    ! variable(POINTER): e_elas_vec_dev_TMP(0:TVEC1,0:NGRAIN1,0:(NUMIND-1))
    ! which is zero-based.
    ! 
    ! SIG_VEC, SIG_KK input as Kirchhoff stresses and output as Cauchy stresses.
    !      
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    ! Deviatoric elastic lattice strain @(t+dt)
    ! SIG_VEC --> e_elas_vec_dev
    !
    DO IPHASE = 1, NUMPHASES
        !
        ! 'MY_PHASE' is 0-based and so are the INDICES
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
        !
        ! Add EL_SUB1 to INDICES
        !
        INDICES=INDICES+EL_SUB1
        !
        IF (ASSOCIATED(E_ELAS_VEC_DEV_TMP)) THEN
            !
            DEALLOCATE(E_ELAS_VEC_DEV_TMP)
            !        
        ENDIF
        !
        ALLOCATE(E_ELAS_VEC_DEV_TMP(0:TVEC1, 0:NGRAIN1, 0:(NUMIND - 1)))
        !
        CALL VEC_D_VEC5(KEINV(:,IPHASE), SIG_VEC(:,:,INDICES), &
            & E_ELAS_VEC_DEV_TMP, NGRAIN, NUMIND)
        !
        E_ELAS_VEC_DEV(:,:,INDICES)=E_ELAS_VEC_DEV_TMP
        !
        DEALLOCATE(E_ELAS_VEC_DEV_TMP)
        DEALLOCATE(INDICES)
        !
    END DO
    !
    ! Deviatoric elastic strain:
    ! 5-vector --> symmetric matrix (3x3)
    ! e_elas_vec_dev --> elas_t3x3
    !
    CALL VEC_MAT_SYMM_GRN(E_ELAS_VEC_DEV, ELAS_T3X3, NGRAIN, M_EL)
    !
    ! Volumetric elastic strain: SPREAD over grains
    ! E_ELAS_KK --> E_ELAS_KK_grn
    !
    E_ELAS_KK_GRN = SPREAD(E_ELAS_KK, DIM = 1, NCOPIES = NGRAIN)
    !
    ! Volumetric kirchhoff stress: SPREAD over grains
    ! SIG_KK --> s_kk_grn
    !
    S_KK_GRN = SPREAD(SIG_KK, DIM = 1, NCOPIES = NGRAIN)
    !
    ! Elastic strain tensor: deviatoric + volumetric
    !
    ELAS_T3X3(0,0,:,:) = ELAS_T3X3(0,0,:,:) + E_ELAS_KK_GRN/3.0_RK
    ELAS_T3X3(1,1,:,:) = ELAS_T3X3(1,1,:,:) + E_ELAS_KK_GRN/3.0_RK
    ELAS_T3X3(2,2,:,:) = ELAS_T3X3(2,2,:,:) + E_ELAS_KK_GRN/3.0_RK
    !
    ! Elastic strain tensor in 6-vec format
    !
    ELAS_TOT6(0,:,:) = ELAS_T3X3(0,0,:,:)
    ELAS_TOT6(1,:,:) = ELAS_T3X3(0,1,:,:)
    ELAS_TOT6(2,:,:) = ELAS_T3X3(0,2,:,:)
    ELAS_TOT6(3,:,:) = ELAS_T3X3(1,1,:,:)
    ELAS_TOT6(4,:,:) = ELAS_T3X3(1,2,:,:)
    ELAS_TOT6(5,:,:) = ELAS_T3X3(2,2,:,:)
    !
    ! Determinant of tensor V* in lattice axes.
    ! V*=I+e*
    !
    V_TENSOR=ELAS_T3X3
    !
    DO I = 0, DIMS1
        !
        V_TENSOR(I, I, :, :) = ELAS_T3X3(I, I, :, :) + 1.0_RK
        !
    END DO
    !
    ! DET(V*)
    !
    CALL DETERMINANT_GRN(V_TENSOR, DETERM_V, NGRAIN, M_EL)
    !
    ! Cauchy Stress at current configuration (crystal COORDS).
    ! Calculate deviatoric part: divide by DET(V*)
    ! RC 6/24/16: Reordered loop ordering for better memory striding
    !
    DO K = EL_SUB1,EL_SUP1
        !
        DO J = 0,NGRAIN1
            !
            DO I = 0,TVEC1
                !
                SIG_VEC(I, J, K) = SIG_VEC(I, J, K) / DETERM_V(J, K)
                !            
            END DO
            !
        END DO
        !
    END DO
    !
    ! Calculate volumetric part: divide by DET(V*)
    !
    S_KK_GRN = S_KK_GRN / DETERM_V
    !
    ! Cauchy stress in sample COORDS
    ! C_ANGS [3x3] --> qr5x5 [5x5]
    !
    CALL ROT_MAT_SYMM(C_ANGS, QR5X5, NGRAIN, M_EL)
    !
    ! deviatoric cauchy stress: SIG_VEC (crystal) --> sig_sm (sample)
    !
    CALL MAT_X_VEC5(QR5X5, SIG_VEC, SIG_SM, NGRAIN, M_EL)
    !
    ! Average Cauchy stress (sample COORDS)
    !
    S_AVG = 0.0_RK
    S_AVG_KK = 0.0_RK
    !
    ! Get weighted average
    !
    DO I = 0, TVEC1
        !    
        S_AVG(I, :) = SUM(SIG_SM(I, :, :) * WTS, DIM = 1)
        !      
    END DO
    !      
    S_AVG_KK = SUM(S_KK_GRN * WTS, DIM = 1)
    !
    ! Deviatoric stress:
    ! 5-vector --> symmetric matrix (3x3)
    ! s_avg    --> S_AVG_3X3
    !
    CALL VEC_MAT_SYMM(S_AVG,S_AVG_3X3,M_EL)
    !
    ! Total stress tensor: deviatoric + volumetric
    !
    S_AVG_3X3(0,0,:) = S_AVG_3X3(0,0,:) + S_AVG_KK/3.0_RK
    S_AVG_3X3(1,1,:) = S_AVG_3X3(1,1,:) + S_AVG_KK/3.0_RK
    S_AVG_3X3(2,2,:) = S_AVG_3X3(2,2,:) + S_AVG_KK/3.0_RK
    !
    ! Update surface information.
    !
    EL_DOF_MIN = 6*EL_SUB1
    EL_DOF_MAX = 6*EL_SUP1 + 5
    !    
    ALLOCATE(SIG_AVG_ALL(EL_DOF_MIN:EL_DOF_MAX))
    !
    DO I=EL_SUB1, EL_SUP1
        !    
        IDOF = 6*I
        SIG_AVG_ALL(IDOF)   = S_AVG_3X3(0,0,I)
        SIG_AVG_ALL(IDOF+1) = S_AVG_3X3(1,0,I)
        SIG_AVG_ALL(IDOF+2) = S_AVG_3X3(2,0,I)
        SIG_AVG_ALL(IDOF+3) = S_AVG_3X3(1,1,I)
        SIG_AVG_ALL(IDOF+4) = S_AVG_3X3(2,1,I)
        SIG_AVG_ALL(IDOF+5) = S_AVG_3X3(2,2,I)
        !      
    END DO
    !
    ! Update areas and loads
    !
    CALL UPD_SURF(OUNITS(DEBUG_U), COORDS, SIG_AVG_ALL, LOAD, AREA)
    !
    DEALLOCATE(SIG_AVG_ALL)
    !
    RETURN
    !    
    END SUBROUTINE CALC_STRESS_STRAIN
    !
    !===========================================================================
    !
    SUBROUTINE EST_AVG_MOD(E_AVG, NU_AVG)
    !
    ! Estimate Voigt-averaged bulk elastic moduli (Hosford p22)
    ! Assumes uniform texture and equal element volumes.
    !      
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! E_AVG:
    ! NU_AVG:
    !
    REAL(RK), INTENT(OUT) :: E_avg, nu_avg
    !
    ! Locals:
    ! MY_PHASE:
    ! I/IPHASE:
    ! PART_NUMEL_PHASE:
    ! NUMEL_PHASE:
    ! PHASE_FRAC:
    ! K_PHASE:
    ! E_PHASE:
    ! NU_PHASE:
    ! C11, C12, C13, C33, C44, C66:
    ! F, G, H:
    !
    INTEGER :: MY_PHASE(EL_SUB1:EL_SUP1)
    INTEGER :: i, IPHASE
    REAL(RK) :: part_numel_phase, numel_phase
    REAL(RK) :: phase_frac
    REAL(RK) :: K_phase, E_phase, nu_phase
    REAL(RK) :: c11, c12, c13, c33, c44, c66
    REAL(RK) :: F, G, H
    !      
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    !
    E_AVG = 0.0
    NU_AVG = 0.0
    !
    DO IPHASE = 1, NUMPHASES
        !
        PART_NUMEL_PHASE = 0.0_RK
        !        
        DO I = EL_SUB1, EL_SUP1
            !
            IF (MY_PHASE(I) .EQ. IPHASE) THEN
                !
                PART_NUMEL_PHASE = PART_NUMEL_PHASE + 1.0
                !          
            END IF
            !
        END DO
        !
        CALL PAR_SUM(PART_NUMEL_PHASE, NUMEL_PHASE)
        PHASE_FRAC = NUMEL_PHASE/REAL(NUMELM, RK) 
        !
        IF ((IPHASE .EQ. 1) .OR. (IPHASE .EQ. 2)) THEN
            !
            ! Cubic (FCC or BCC)
            !
            C11 = ELAS_COEFFS(0,IPHASE)
            C12 = ELAS_COEFFS(1,IPHASE)
            C44 = ELAS_COEFFS(3,IPHASE)
            E_PHASE = (C11-C12+3*C44)*(C11+2*C12)/(2*C11+3*C12+C44)
            !
        ELSE IF (IPHASE .EQ. 3) THEN
            !
            ! Hexagonal (HCP)
            !
            C11 = ELAS_COEFFS(0,IPHASE)
            C12 = ELAS_COEFFS(1,IPHASE)
            C13 = ELAS_COEFFS(2,IPHASE)
            C44 = ELAS_COEFFS(3,IPHASE)
            C33 = C11+C12-C13
            C66 = (C11-C12)/2
            !
            F = (2*C11+C33)/3
            G = (C12+2*C13)/3
            H = (2*C44+C66)/3
            E_PHASE = (F-G+3*H)*(F+2*G)/(2*F+3*G+H)
            !
        ELSE
            !
          CALL PAR_QUIT('Error  :     > Invalid crystal type.')
            !
        END IF
        !
        K_PHASE = CRYSTAL_PARM(8, IPHASE)
        NU_PHASE = (3*K_PHASE-E_PHASE)/(6*K_PHASE)
        !
        E_AVG = E_AVG + PHASE_FRAC*E_PHASE
        NU_AVG = NU_AVG + PHASE_FRAC*NU_PHASE
        !
    END DO
    !
    RETURN
    !  
    END SUBROUTINE EST_AVG_MOD
    !
    !===========================================================================
    !
    SUBROUTINE SOLVE_LIN_SYS_3(MAT, VEC, SOL)
    !
    ! Solve a linear system of three equations. Used for triaxial loading.
    !      
    !---------------------------------------------------------------------------
    !      
    ! Arguments:
    !
    REAL(RK), INTENT(IN) :: MAT(3,3)
    REAL(RK), INTENT(IN) :: VEC(3)
    REAL(RK), INTENT(OUT) :: SOL(3)
    !
    ! Locals:
    !
    REAL(RK) :: INV(3,3)
    REAL(RK) :: A, B, C, D, E, F, G, H, K
    REAL(RK) :: DET, MATNORM, INVNORM, COND
    !      
    !---------------------------------------------------------------------------
    !      
    A = MAT(1, 1)
    B = MAT(1, 2)
    C = MAT(1, 3)
    D = MAT(2, 1)
    E = MAT(2, 2)
    F = MAT(2, 3)
    G = MAT(3, 1)
    H = MAT(3, 2)
    K = MAT(3, 3)
    !
    DET = A * (E * K - F * H) - B * (D * K - F * G) + C * (D * H - E * G)
    !
    INV(1, 1) = (E * K - F * H) / DET
    INV(1, 2) = (C * H - B * K) / DET
    INV(1, 3) = (B * F - C * E) / DET
    INV(2, 1) = (F * G - D * K) / DET
    INV(2, 2) = (A * K - C * G) / DET
    INV(2, 3) = (C * D - A * F) / DET
    INV(3, 1) = (D * H - E * G) / DET
    INV(3, 2) = (B * G - A * H) / DET
    INV(3, 3) = (A * E - B * D) / DET
    !
    SOL = MATMUL(INV, VEC)
    !
    ! Find conditioning number
    !
    MATNORM = MAX(MAT(1, 1) + MAT(1, 2) + MAT(1, 3), &
         & MAT(2, 1) + MAT(2, 2) + MAT(2, 3), &
         & MAT(3, 1) + MAT(3, 2) + MAT(3, 3))
    INVNORM = MAX(INV(1, 1) + INV(1, 2) + INV(1, 3), &
         & INV(2, 1) + INV(2, 2) + INV(2, 3), &
         & INV(3, 1) + INV(3, 2) + INV(3, 3))
    COND = MATNORM * INVNORM
    !
    IF (COND .GT. 1e3) THEN
        !
        CALL PAR_QUIT('Error  :     > Matrix is poorly conditioned.')
        !
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE solve_lin_sys_3
    !
    !===========================================================================
    !
    SUBROUTINE TEMP_UPDATE_STATE_EVPS(VELOCITY, DTRACE, C0_ANGS, C_ANGS, &
        & SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, E_BAR_VEC, &
        & E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, KEINV, WTS, DTIME, &
        & CONVERGED_SOLUTION, AUTO_TIME)
    !
    ! Temporarily updates state at center of element only.
    ! This subroutine is employed when iterating on velocity boundary
    !   conditions.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    REAL(RK), INTENT(INOUT) :: VELOCITY(DOF_SUB1:DOF_SUP1)
    TYPE(TRACE) DTRACE
    REAL(RK), INTENT(IN) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: SIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: SIG_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK_BAR(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: SIG_KK(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: KEINV(0:TVEC1,1:NUMPHASES)
    REAL(RK), INTENT(IN) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: DTIME
    LOGICAL, INTENT(INOUT) :: CONVERGED_SOLUTION
    INTEGER, INTENT(IN) :: AUTO_TIME
    !
    ! Locals:
    !
    REAL(RK), POINTER :: E_BAR_VEC_TMP(:,:,:) => NULL()
    INTEGER :: I, J, M_EL, IQPT
    INTEGER :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: S11(EL_SUB1:EL_SUP1), S12(EL_SUB1:EL_SUP1), S13(EL_SUB1:EL_SUP1)
    REAL(RK) :: S21(EL_SUB1:EL_SUP1), S22(EL_SUB1:EL_SUP1), S23(EL_SUB1:EL_SUP1)
    REAL(RK) :: S31(EL_SUB1:EL_SUP1), S32(EL_SUB1:EL_SUP1), S33(EL_SUB1:EL_SUP1)
    REAL(RK) :: VGRAD(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: D(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: D_VEC(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DEFF(EL_SUB1:EL_SUP1)
    REAL(RK) :: W(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: W_VEC(0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDX(0:nnpe, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDY(0:nnpe, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDZ(0:nnpe, EL_SUB1:EL_SUP1)
    REAL(RK) :: LOC0
    REAL(RK) :: LOC1
    REAL(RK) :: LOC2
    REAL(RK) :: DET(EL_SUB1:EL_SUP1)
    REAL(RK) :: D_KK(EL_SUB1:EL_SUP1)
    REAL(RK) :: ECOORDS(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: EVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SQR2, SQR32
    REAL(RK) :: ALPHA(TVEC), ELAPSED
    INTEGER :: E1, E2, M, IPHASE, NUMIND
    INTEGER, POINTER :: INDICES(:) => NULL()
    INTEGER :: MY_PHASE(0:(EL_SUP1-EL_SUB1))
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    SQR2  = SQRT(2.0_RK)
    SQR32 = SQRT(1.5_RK)
    !
    ! COORDS --> ECOORDS [30 x m]
    ! velocity --> EVEL  [30 x m]
    !
    CALL PART_GATHER(ECOORDS, COORDS, NODES, DTRACE)
    CALL PART_GATHER(EVEL, VELOCITY, NODES, DTRACE)
    !
    ! Coordinates in the parent element
    !
    LOC0 = 0.25_RK
    LOC1 = 0.25_RK
    LOC2 = 0.25_RK
    !
    ! Initialize state variables at centroid
    !
    SIG_VEC_N = PSIG_VEC_N
    E_ELAS_KK_BAR = PELA_KK_BAR
    !
    ! Compute elastic strain: SIG_VEC_N --> E_BAR_VEC
    !
    DO IPHASE = 1, NUMPHASES
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
        INDICES = INDICES + EL_SUB1
        !
        IF (ASSOCIATED(E_BAR_VEC_TMP)) THEN
            !
            DEALLOCATE(E_BAR_VEC_TMP)
            !
        ENDIF
        !
        ALLOCATE(E_BAR_VEC_TMP(0:TVEC1, 0:NGRAIN1, 0:(NUMIND - 1)))
        !
        CALL VEC_D_VEC5(KEINV(:, IPHASE), SIG_VEC_N(:, :, INDICES), &
            & E_BAR_VEC_TMP, NGRAIN, NUMIND)
        !
        E_BAR_VEC(:, :, INDICES) = E_BAR_VEC_TMP
        !
        DEALLOCATE(E_BAR_VEC_TMP)
        DEALLOCATE(INDICES)
        !
    ENDDO
    !
    ! Compute quadrature quantities given a set of local coordinates
    !
    CALL SFDER_HPAR(LOC0, LOC1, LOC2, ECOORDS, DNDX, DNDY, DNDZ, DET, S11, &
        & S12, S13, S21, S22, S23, S31, S32, S33)
    !
    ! Compute velocity gradient (VGRAD) and his sym. (d) and skew (w) parts.
    !
    CALL VEL_GRADIENT(VGRAD, DNDX, DNDY, DNDZ, EVEL)
    !
    ! D_KK: mean/volumetric part of the symmetric part of the velocity gradient
    ! D: deviatoric part of the symmetric part of the velocity gradient
    ! W: skew part of the velocity gradient
    !
    CALL SYMM_VGR(D, D_KK, VGRAD, M_EL)
    CALL SKEW_VGR(W, VGRAD, M_EL)
    !
    ! D [3x3] --> D_VEC {5}
    ! W [3x3] --> W_VEC {3}
    !
    CALL MAT_VEC_SYMM(D, D_VEC, M_EL)
    CALL MAT_VEC_SKEW(W, W_VEC, M_EL)
    CALL EFF_DEF(DEFF, D, DTIME, M_EL)
    !
    ! Variables @(t+dt):
    ! RSTAR [3x3]
    ! C_ANGS [3x3]
    ! CRSS (1)
    ! SIG_VEC {5}
    ! SIG_KK (1)
    ! E_ELAS_KK (1)
    !
    CALL POLYCRYSTAL_RESPONSE_EVPS_QP(D_VEC, W_VEC, C0_ANGS, C_ANGS, &
        & SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, E_BAR_VEC, WTS, &
        & DEFF, D_KK, SIG_KK, E_ELAS_KK_BAR, E_ELAS_KK, JITER_STATE, KEINV, &
        & 9999, DTIME, CONVERGED_SOLUTION, AUTO_TIME)
    !
    IF (.NOT. CONVERGED_SOLUTION .AND. AUTO_TIME .EQ. 1) THEN
        !
        CALL PAR_QUIT('Error  :     > No converged solution found.')
        !
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE TEMP_UPDATE_STATE_EVPS
    !
    !===========================================================================
    !
    SUBROUTINE UPDATE_STATE_EVPS(VELOCITY, DTRACE, C0_ANGS, C_ANGS, &
        & SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, E_BAR_VEC, &
        & E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, KEINV, WTS, DEFF, DTIME, STIF, FE, &
        & D, D_VEC, VGRAD, CONVERGED_SOLUTION, AUTO_TIME)
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    REAL(RK), INTENT(INOUT) :: VELOCITY(DOF_SUB1:DOF_SUP1)
    TYPE(TRACE) DTRACE
    REAL(RK), INTENT(IN) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: SIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: SIG_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(INOUT) :: RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK_BAR(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: SIG_KK(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: KEINV(0:TVEC1, 1:NUMPHASES)
    REAL(RK), INTENT(IN) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: DEFF(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: DTIME
    REAL(RK), INTENT(OUT) :: STIF(TVEC, TVEC, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: FE(TVEC, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: D(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: D_VEC(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: VGRAD(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    LOGICAL, INTENT(INOUT) :: CONVERGED_SOLUTION
    INTEGER, INTENT(IN) :: AUTO_TIME
    !
    ! Locals:
    !
    REAL(RK), POINTER ::  E_BAR_VEC_TMP(:,:,:) => NULL()
    INTEGER  :: I, J, M_EL, IQPT
    INTEGER  :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: S11(EL_SUB1:EL_SUP1), S12(EL_SUB1:EL_SUP1), S13(EL_SUB1:EL_SUP1)
    REAL(RK) :: S21(EL_SUB1:EL_SUP1), S22(EL_SUB1:EL_SUP1), S23(EL_SUB1:EL_SUP1)
    REAL(RK) :: S31(EL_SUB1:EL_SUP1), S32(EL_SUB1:EL_SUP1), S33(EL_SUB1:EL_SUP1)
    REAL(RK) :: W(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: W_VEC(0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDX(0:nnpe, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDY(0:nnpe, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDZ(0:nnpe, EL_SUB1:EL_SUP1)
    REAL(RK) :: LOC0
    REAL(RK) :: LOC1
    REAL(RK) :: LOC2
    REAL(RK) :: DET(EL_SUB1:EL_SUP1)
    REAL(RK) :: D_KK(EL_SUB1:EL_SUP1)
    REAL(RK) :: ECOORDS(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: EVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SQR2, SQR32
    REAL(RK) :: ALPHA(TVEC), ELAPSED
    REAL(RK) ::  D_VEC_Q(0:TVEC1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) ::  W_VEC_Q(0:DIMS1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) ::  C_ANGS_Q(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) ::  RSTAR_Q(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) ::  SIG_VEC_N_Q(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) ::  SIG_VEC_Q(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) ::  E_BAR_VEC_Q(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) ::  CRSS_Q(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) ::  E_ELAS_KK_BAR_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) ::  E_ELAS_KK_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) ::  SIG_KK_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) ::  DEFF_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) ::  DEFF_G(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) ::  D_KK_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) ::  SHEAR(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) ::  SHRATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) ::  WP_SS(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    INTEGER :: E1, E2, M, IPHASE, NUMIND
    INTEGER, POINTER :: INDICES(:) => NULL()
    INTEGER :: MY_PHASE(0:(EL_SUP1-EL_SUB1))
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    SQR2  = SQRT(2.0_RK)
    SQR32 = SQRT(1.5_RK)
    !
    ! Compute coordinates @(t+dt) using the velocity @(t+dt)
    !
    COORDS = COORDS + VELOCITY * DTIME
    !
    ! COORDS --> ECOORDS [30 x m]
    ! VELOCITY --> EVEL  [30 x m]
    !
    CALL PART_GATHER(ECOORDS, COORDS, NODES, DTRACE)
    CALL PART_GATHER(EVEL, VELOCITY, NODES, DTRACE)
    !
    ! Looping over gauss quadrature points
    !
    SIG_VEC_N_Q(:, :, :, :) = GSIG_VEC_N(:, :, :, :)
    E_ELAS_KK_BAR_Q(:, :) = GELA_KK_BAR(0, :, :)
    C_ANGS_Q = SPREAD(C_ANGS, 5, NQPT)
    RSTAR_Q = SPREAD(RSTAR, 5, NQPT)
    !
    DO IQPT = 0, NQPT1
        !
        ! Coordinates in the parent element
        !
        LOC0 = QPLOC(0, IQPT)
        LOC1 = QPLOC(1, IQPT)
        LOC2 = QPLOC(2, IQPT)
        !
        ! Initialize internal variables at quadrature points
        ! Compute elastic strain: SIG_VEC_N --> E_BAR_VEC
        !
        DO IPHASE = 1, NUMPHASES
            !
            CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
            INDICES = INDICES + EL_SUB1
            !
            IF (ASSOCIATED(E_BAR_VEC_TMP)) THEN
                !
                DEALLOCATE(E_BAR_VEC_TMP)
                !
            ENDIF
            !
            ALLOCATE(E_BAR_VEC_TMP(0:TVEC1 ,0:NGRAIN1, 0:(NUMIND - 1)))
            !
            CALL VEC_D_VEC5(KEINV(:, IPHASE), &
                & SIG_VEC_N_Q(:, :, INDICES, IQPT), E_BAR_VEC_TMP, NGRAIN, &
                & NUMIND)
           E_BAR_VEC_Q(:, :, INDICES, IQPT) = E_BAR_VEC_TMP
           !
           DEALLOCATE(E_BAR_VEC_TMP)
           DEALLOCATE(INDICES)
           !
        ENDDO !NUMPHASES
        !
        ! Compute quadrature quantities given a set of local coordinates.
        !
        CALL SFDER_HPAR(LOC0, LOC1, LOC2, ECOORDS, DNDX, DNDY, DNDZ, DET, S11, &
            & S12, S13, S21, S22, S23, S31, S32, S33)
        !
        ! Compute velocity gradient (VGRAD) and its sym. (D) and skew (W) parts.
        !
        CALL VEL_GRADIENT(VGRAD, DNDX, DNDY, DNDZ, EVEL)
        !
        ! D_KK: mean/volumetric/ part of the symm part of the velocity gradient
        ! D: deviatoric part of the symmetric part of the velocity gradient
        ! W: skew part of the velocity gradient
        !
        CALL SYMM_VGR(D, D_KK_Q(:, IQPT), VGRAD, M_EL)
        CALL SKEW_VGR(W, VGRAD, M_EL)
        !
        ! D [3x3] --> D_VEC {5}
        ! W [3x3] --> W_VEC {3}
        !
        CALL MAT_VEC_SYMM(d, D_VEC_Q(:,:,IQPT), M_EL)
        CALL MAT_VEC_SKEW(w, W_VEC_Q(:,:,IQPT), M_EL)
        !
        CALL EFF_DEF(DEFF_Q(:,IQPT), D, DTIME, M_EL)
        !
        ! Variables @(t+dt):
        ! RSTAR [3x3]
        ! C_ANGS [3x3]
        ! CRSS (1)
        ! SIG_VEC {5}
        ! SIG_KK (1)
        ! E_ELAS_KK (1)
        !
    ENDDO ! loop over quad points
    !
    CALL POLYCRYSTAL_RESPONSE_EVPS(D_VEC_Q, W_VEC_Q, C0_ANGS, C_ANGS_Q, &
        & SIG_VEC_N_Q, SIG_VEC_Q, CRSS_N, CRSS_Q, RSTAR_N, RSTAR_Q, &
        & E_BAR_VEC_Q, WTS, DEFF_Q, D_KK_Q, SIG_KK_Q, E_ELAS_KK_BAR_Q, &
        & E_ELAS_KK_Q, JITER_STATE, KEINV, 9999, DTIME, CONVERGED_SOLUTION, &
        & AUTO_TIME)
    !
    IF (.NOT. CONVERGED_SOLUTION .AND. AUTO_TIME .EQ. 1) THEN
        !
        CALL PAR_QUIT('Error  :     > No converged solution found.')
        !
    ENDIF
    !
    DO I = 0, NQPT1
        !
        DEFF_G(:,:) = SPREAD(DEFF_Q(:, I), DIM = 1, NCOPIES = NGRAIN)
        !
        CALL SLIP_QNT(WP_SS, SHEAR, SHRATE, SIG_VEC_Q(:, :, :, I), &
            & CRSS_Q(:, :, :, I), DEFF_G, NGRAIN, M_EL)
        !
        GACCUMSHEAR(:, 0, :, I) = GACCUMSHEAR(:, 0, :, I) + &
            & DABS(SHEAR(: ,0, :)) * DTIME
        !
    ENDDO
    !
    ! Update _n variables
    !
    GSIG_VEC_N(:, :, :, :) = SIG_VEC_Q(:, :, :, :)
    GELA_KK_BAR(0, :, :) = E_ELAS_KK_Q(:, :)
    !
    ! Note: CRSS_N and RSTAR_N at centroid will be stored.
    ! Repeat the update at the center of the element for post-processing
    !
    ! Coordinates in the parent element
    !
    LOC0 = 0.25_RK
    LOC1 = 0.25_RK
    LOC2 = 0.25_RK
    !
    ! Initialize state variables at centroid
    !
    SIG_VEC_N = PSIG_VEC_N
    E_ELAS_KK_BAR = PELA_KK_BAR
    !
    ! Compute elastic strain: SIG_VEC_N --> E_BAR_VEC
    !
    DO IPHASE = 1, NUMPHASES
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
        INDICES=INDICES+EL_SUB1
        !
        IF (ASSOCIATED(E_BAR_VEC_TMP)) THEN
            !
            DEALLOCATE(E_BAR_VEC_TMP)
            !
        ENDIF
        !
        ALLOCATE(E_BAR_VEC_TMP(0:TVEC1, 0:NGRAIN1, 0:(NUMIND - 1)))
        CALL VEC_D_VEC5(KEINV(:, IPHASE), SIG_VEC_N(:, :, INDICES), &
            & E_BAR_VEC_TMP, NGRAIN, NUMIND)
        E_BAR_VEC(:, :, INDICES) = E_BAR_VEC_TMP
        !
        DEALLOCATE(E_BAR_VEC_TMP)
        DEALLOCATE(INDICES)
        !
    ENDDO
    !
    ! Compute quadrature quantities given a set of local coordinates.
    !
    CALL SFDER_HPAR(LOC0, LOC1, LOC2, ECOORDS, DNDX, DNDY, DNDZ, DET, S11, &
        & S12, S13, S21, S22, S23, S31, S32, S33)
    !
    ! Compute velocity gradient (VGRAD) and his sym. (d) and skew (w) parts.
    !
    CALL vel_gradient(VGRAD, DNDX, DNDY, DNDZ, EVEL)
    !
    ! D_KK: mean/volumetric part of the symmetric part of the velocity gradient
    ! D: deviatoric part of the symmetric part of the velocity gradient
    ! W: skew part of the velocity gradient
    !
    CALL SYMM_VGR(D, D_KK, VGRAD, M_EL)
    CALL SKEW_VGR(W, VGRAD, M_EL)
    !
    ! D [3x3] --> D_VEC {5}
    ! W [3x3] --> W_VEC {3}
    !
    CALL MAT_VEC_SYMM(D, D_VEC, M_EL)
    CALL MAT_VEC_SKEW(W, W_VEC, M_EL)
    !
    CALL EFF_DEF(DEFF, d, DTIME, M_EL)
    !
    ! Variables @(t+dt):
    ! RSTAR [3x3]
    ! C_ANGS [3x3]
    ! CRSS (1)
    ! SIG_VEC {5}
    ! SIG_KK (1)
    ! E_ELAS_KK (1)
    !
    ACCUMSHEAR = ACCUMSHEAR_CEN
    !
    CALL POLYCRYSTAL_RESPONSE_EVPS_QP(D_VEC, W_VEC, C0_ANGS, C_ANGS, &
        & SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, RSTAR, E_BAR_VEC, WTS, &
        & DEFF, D_KK, SIG_KK, E_ELAS_KK_BAR, E_ELAS_KK, JITER_STATE, KEINV, &
        & 9999, DTIME, CONVERGED_SOLUTION, AUTO_TIME)
    !
    IF (.NOT. CONVERGED_SOLUTION .AND. AUTO_TIME .EQ. 1) THEN
        !
        CALL PAR_QUIT('Error  :     > No converged solution found.')
        !
    ENDIF
    !
    DEFF_G(:,:) = SPREAD(DEFF, DIM = 1, NCOPIES = NGRAIN)
    !
    CALL SLIP_QNT(WP_SS, SHEAR, SHRATE, SIG_VEC, CRSS, DEFF_G, NGRAIN, M_EL )
    !
    ACCUMSHEAR_CEN(:, 0, :) = ACCUMSHEAR_CEN(:, 0, :) + &
        & DABS(SHEAR(:, 0, :)) * DTIME
    !
    ! Store internal _n variables for post processing
    !
    PSIG_VEC_N  = SIG_VEC
    PELA_KK_BAR = E_ELAS_KK
    !
    ! Update global internal variables for critial resovolved shear stress
    !   and R^*_n at the center of each element.
    CRSS_N  = CRSS
    RSTAR_N = RSTAR
    !
    ! Fix STIF & FE to be used for stress post-processing.
    !
    ALPHA(1) = SQR2
    ALPHA(2) = 1.0_RK / SQR32
    ALPHA(3) = 1.0_RK / SQR2
    ALPHA(4) = 1.0_RK / SQR2
    ALPHA(5) = 1.0_RK / SQR2
    !
    DO J = 1, TVEC
        !
        FE(J, :) = ALPHA(J) * FE(J, :)
        !
        DO I = 1, TVEC
            !
            STIF(I, J, :) = ALPHA(I) * STIF(I, J, :) * ALPHA(J)
            !
        ENDDO
        !
    ENDDO
    !
    RETURN
    !
    END SUBROUTINE UPDATE_STATE_EVPS
    !
    !===========================================================================
    !
    SUBROUTINE VELOCITY_ITERATION(BCS, PFORCE, VELOCITY, ELPRESS, EVEL, LOAD, &
        & DTRACE, NTRACE, C0_ANGS, C_ANGS, SIG_VEC_N, SIG_VEC, CRSS_N, CRSS, &
        & RSTAR_N, RSTAR, TRSTAR, KEINV, E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, &
        & JITER_STATE, WTS, EPSEFF, DTIME, INCR, E_BAR_VEC, &
        & CONVERGED_SOLUTION, AUTO_TIME, ITERNL)
    !
    ! Perform velocity iteration
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    

    
    LOGICAL, INTENT(IN) :: BCS(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(IN) :: PFORCE(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(INOUT) :: VELOCITY(DOF_SUB1:DOF_SUP1)
    REAL(RK), INTENT(OUT) :: ELPRESS(EL_SUB1:EL_SUP1)
    REAL(RK) :: EVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: LOAD(3)
    TYPE(TRACE) :: DTRACE
    TYPE(TRACE) :: NTRACE
    REAL(RK), INTENT(IN) :: C0_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: SIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: SIG_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: CRSS_N(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: CRSS (0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: TRSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, &
        & EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: KEINV(0:TVEC1, 1:NUMPHASES)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK_BAR(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: E_ELAS_KK(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: SIG_KK(EL_SUB1:EL_SUP1)
    INTEGER, INTENT(OUT) :: JITER_STATE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: EPSEFF(EL_SUB1:EL_SUP1)
    REAL(RK) :: DTIME
    INTEGER, INTENT(IN) :: INCR   ! current load INCRement
    REAL(RK), INTENT(OUT) :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    LOGICAL, INTENT(OUT) :: CONVERGED_SOLUTION
    INTEGER :: AUTO_TIME
    INTEGER, INTENT(OUT) :: ITERNL
    !
    ! Locals:
    !
    INTEGER  :: IER
    REAL(RK) :: TSIG_VEC_N(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TC_ANGS (0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TCOORDS(DOF_SUB1:DOF_SUP1)
    REAL(RK) :: S_AVG_3X3(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: ELAS_TOT6(0:TVEC , 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SURF_LOAD_ARRAY(NSURFACES, 3)
    REAL(RK) :: AREA(NSURFACES)
    REAL(RK) :: CRSS_Q(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: EPSEFF_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: ELPRESS_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: SIG_VEC_Q(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: SIG_VEC_N_Q(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: E_ELAS_KK_BAR_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: E_ELAS_KK_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: SIG_KK_Q(EL_SUB1:EL_SUP1, 0:NQPT1)
    REAL(RK) :: E_BAR_VEC_Q(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1, 0:NQPT1)
    !
    !---------------------------------------------------------------------------
    !
    ! Initialization
    !
    CONVERGED_SOLUTION = .TRUE.
    !
    ! Assign temporary variables for itmethod_evps and temp_update_state_evps
    !
    TRSTAR = RSTAR
    TSIG_VEC_N = SIG_VEC_N
    SIG_VEC_N_Q = SPREAD(SIG_VEC_N, 4, NQPT)
    E_ELAS_KK_BAR_Q = SPREAD(E_ELAS_KK_BAR, 2, NQPT)
    TC_ANGS = C_ANGS
    !
    ! Iterate on velocity field
    !
    IER = ITMETHOD_EVPS(BCS, PFORCE, VELOCITY, ELPRESS_Q, EVEL, DTRACE, &
        & NTRACE, C0_ANGS, C_ANGS, SIG_VEC_N_Q, SIG_VEC_Q, CRSS_N, CRSS_Q, &
        & RSTAR_N, TRSTAR, KEINV, E_ELAS_KK_BAR_Q, E_ELAS_KK_Q, SIG_KK_Q, &
        & JITER_STATE, WTS, EPSEFF_Q, DTIME, INCR, E_BAR_VEC_Q, &
        & CONVERGED_SOLUTION, AUTO_TIME, ITERNL)
    !
    ! Evaluate convergence of velocity field and material state
    !
    IF (IER .LT. 0) THEN
        !
        CALL PAR_QUIT('Error  :     > Failure to converge.')
        !
    ENDIF
    !
    IF (.NOT. CONVERGED_SOLUTION) WRITE(DFLT_U,'(A)') 'Warning:     > AUTO_TIME&
        & .ne. 1 with CONVERGED_SOLUTION = false'
    !
    ! Temporatily update coordinates
    !
    TCOORDS = COORDS
    COORDS = COORDS + VELOCITY * DTIME
    !
    CALL TEMP_UPDATE_STATE_EVPS(VELOCITY, DTRACE, C0_ANGS, TC_ANGS, &
        & TSIG_VEC_N, SIG_VEC, CRSS_N, CRSS, RSTAR_N, TRSTAR, E_BAR_VEC, &
        & E_ELAS_KK_BAR, E_ELAS_KK, SIG_KK, KEINV, WTS, DTIME, &
        & CONVERGED_SOLUTION, AUTO_TIME)
    !
    ! Compute load
    !
    CALL CALC_STRESS_STRAIN(S_AVG_3X3, ELAS_TOT6, SURF_LOAD_ARRAY, AREA, &
        & SIG_VEC, SIG_KK, E_ELAS_KK, TC_ANGS, KEINV, WTS)
    !
    LOAD(1) = SURF_LOAD_ARRAY(2,1) ! 4,1
    LOAD(2) = SURF_LOAD_ARRAY(4,2) ! 6,2
    LOAD(3) = SURF_LOAD_ARRAY(6,3) ! 2,3
    !
    ! Reset coordinates
    !
    COORDS = TCOORDS
    !
    RETURN
    ! 
    END SUBROUTINE VELOCITY_ITERATION
    !
END MODULE DRIVER_UTILITIES_MOD
