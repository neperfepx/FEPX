! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE STRESS_SOLVE_EVPS_MOD
!
! Add module description here.
!
! Contains subroutines:
! STRESS_SOLVE_EVPS: Add descriptions here.
! SOLVE_NEWTON_EVPS:
! MATRIX_FJAC:
! JACOB5X5:
! DP_WP_HAT:
! WP_HAT_MAT5X5:
! WP_HAT_MAT5X5_ALL:
! RESIDUAL:
! CHECK_DIAGONALS_EVPS:
! SYMMETRIZE_JAC:
!
USE INTRINSICTYPESMODULE, RK=>REAL_KIND
USE PARALLEL_MOD
!
USE CONVERGENCEMODULE, ONLY: CV_OPTIONS
USE DIMSMODULE
USE MICROSTRUCTURE_MOD
USE READ_INPUT_MOD
USE STRESSSOLVEVPMODULE
USE UNITS_MOD
USE MATRIX_OPERATIONS_MOD
!
IMPLICIT NONE
!
PRIVATE
!
PUBLIC :: STRESS_SOLVE_EVPS, WP_HAT_MAT5X5_ALL
!
CONTAINS
    !
    SUBROUTINE STRESS_SOLVE_EVPS(SIG, D_VEC_LAT, W_VEC_LAT, E_BAR_VEC, &
        & CRSS, KEINV, DTIME, WP_HAT, ITER_STATE, DONE, CONVERGED_NEWTON)
    !
    ! Add description here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER, INTENT(IN)     :: ITER_STATE      
    LOGICAL, INTENT(INOUT)  :: CONVERGED_NEWTON
    LOGICAL, INTENT(IN)     :: DONE(0:NGRAIN1, EL_SUB1:EL_SUP1)
    !
    REAL(RK), INTENT(IN)    :: DTIME
    REAL(RK), INTENT(IN)    :: KEINV(0:TVEC1,1:NUMPHASES)
    REAL(RK), INTENT(INOUT) :: SIG(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)    :: D_VEC_LAT(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)    :: W_VEC_LAT(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)    :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)    :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT)   :: WP_HAT(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    INTEGER  :: IER , NUM, M_EL
    INTEGER  :: JITER(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SIG0_AVG, SIG1_AVG, SIG2_AVG, SIG3_AVG, SIG4_AVG
    LOGICAL  :: CONVERGED(0:NGRAIN1, EL_SUB1:EL_SUP1)
    !
    !---------------------------------------------------------------------------
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    CONVERGED = .FALSE.
    WHERE(DONE) CONVERGED = .TRUE.
    !
    ! Solve NLE for stresses
    !
    CALL SOLVE_NEWTON_EVPS(SIG, D_VEC_LAT, E_BAR_VEC, W_VEC_LAT, &
        & CRSS, CV_OPTIONS % SX_MAX_ITERS_NEWTON, IER, &
        & KEINV, DTIME, WP_HAT, JITER, &
        & CONVERGED, DONE, NGRAIN, M_EL)
    !
    IF (IER .EQ. 0) RETURN
    !
    WRITE(DFLT_U, '(A,I0)') 'Warning:       . STRESS_SOLVE_EVPS: ', &
        & count(.not. converged), ' elements did not converge'
    WRITE(DFLT_U, '(A,I0)') 'Warning:       .                    in iter_state &
        &= ', iter_state
    WRITE(DFLT_U, '(A)') 'Warning:       . These elements will receive &
        &average stress values'
    !
    NUM = COUNT(CONVERGED)
    !
    ! Adjust values to avoid divide-by-zero exceptions - deb
    IF (NUM .EQ. 0) NUM = 1
    !
    SIG0_AVG = SUM(SIG(0, :, :), MASK = CONVERGED) / NUM
    SIG1_AVG = SUM(SIG(1, :, :), MASK = CONVERGED) / NUM
    SIG2_AVG = SUM(SIG(2, :, :), MASK = CONVERGED) / NUM
    SIG3_AVG = SUM(SIG(3, :, :), MASK = CONVERGED) / NUM
    SIG4_AVG = SUM(SIG(4, :, :), MASK = CONVERGED) / NUM
    !
    WHERE (.NOT. CONVERGED)
        !
        SIG(0, :, :) = SIG0_AVG
        SIG(1, :, :) = SIG1_AVG
        SIG(2, :, :) = SIG2_AVG
        SIG(3, :, :) = SIG3_AVG
        SIG(4, :, :) = SIG4_AVG
        !
    ENDWHERE
    !
    CONVERGED_NEWTON = .FALSE.
    !
    RETURN
    !    
    END SUBROUTINE STRESS_SOLVE_EVPS
    !
    !===========================================================================
    !
    SUBROUTINE SOLVE_NEWTON_EVPS(SIG, D_VEC_LAT, E_BAR_VEC, W_VEC_LAT, &
        & CRSS, MAX_ITER_NEWTON, IRC, KEINV, DT, WP_HAT_VEC, JITER, CONVERGED, &
        & DONE, N, M)
    !
    ! Add description here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER, INTENT(IN)  :: MAX_ITER_NEWTON, N, M
    INTEGER, INTENT(OUT) :: JITER(0:(N - 1), 0:(M - 1))
    !
    INTEGER, INTENT(OUT)   :: IRC
    LOGICAL, INTENT(INOUT) :: CONVERGED(0:(N - 1), 0:(M - 1))
    LOGICAL, INTENT(IN)    :: DONE(0:(N - 1), 0:(M - 1))
    !
    REAL(RK), INTENT(IN)    :: DT
    REAL(RK), INTENT(IN)    :: KEINV(0:TVEC1, 1:NUMPHASES)
    REAL(RK), INTENT(INOUT) :: SIG(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)    :: D_VEC_LAT(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)    :: E_BAR_VEC(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)    :: W_VEC_LAT(0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)    :: CRSS(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT)   :: WP_HAT_VEC(0:DIMS1, 0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    LOGICAL  :: NEWTON_OK(0:(N - 1), 0:(M - 1))
    !
    INTEGER  :: ITER_NEWTON, ISLIP, I, NM, INEWTON, N_SLIP
    !
    REAL(RK) :: C1, XN(NUMPHASES), XNN(NUMPHASES)
    REAL(RK) :: SIG_0(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: E_BAR(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: E_VEC(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: E_ELAS(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    !    
    REAL(RK) :: GDOT(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: DGDOT(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: DP_HAT_VEC(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: WP_X_E(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: DEL_S(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: FJ(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: RES(0:(N - 1), 0:(M - 1)), RES_N(0:(N - 1), 0:(M - 1))
    REAL(RK) :: FACT(0:(N - 1), 0:(M - 1))
    REAL(RK) :: XLAMBDA(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: FJAC(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    !      
    REAL(RK) :: WP_HAT_MATX(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: TAU(0:(N - 1), 0:(M - 1)), TAUA(0:(N - 1), 0:(M - 1))
    REAL(RK) :: FJ21(0:(N - 1), 0:(M - 1))
    REAL(RK) :: RATIO_RES(0:(N - 1), 0:(M - 1))
    REAL(RK) :: RES_AUX(0:(N - 1), 0:(M - 1))
    !
    REAL(RK), POINTER :: P_HAT_VEC(:,:) => NULL()
    REAL(RK), POINTER :: PPT(:,:,:) => NULL()
    REAL(RK), POINTER :: E_VEC_TMP(:,:,:) => NULL()
    REAL(RK), POINTER :: E_ELAS_TMP(:,:,:,:) => NULL()
    REAL(RK), POINTER :: WP_X_E_TMP(:,:,:) => NULL()
    REAL(RK), POINTER :: FJAC_TMP(:,:,:,:) => NULL()
    !   
    REAL(RK), PARAMETER :: ONE_DP = 1.0_RK
    REAL(RK) :: UFLOW = TINY(ONE_DP)
    REAL(RK) :: TOOSMALL(NUMPHASES)
    !
    INTEGER  :: MY_PHASE(0:(M - 1)), IPHASE, NUMIND, IN
    INTEGER, POINTER :: INDICES(:) => NULL()
    !
    REAL(RK) :: ANISO_M_TEMP(0:17)
    REAL(RK) :: AXNN(0:17,NUMPHASES)
    REAL(RK) :: AXN(0:17,NUMPHASES)
    REAL(RK) :: ATOOSMALL(0:17,NUMPHASES)
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    ! 
    IRC   = 0
    JITER = 0
    SIG_0 = 0.0_RK
    NEWTON_OK = .TRUE.
    !
    NM  = N * M
    C1  = 1.0_RK / DT
    !
    ! e_bar_vec {5} --> e_bar [3x3]sym
    !
    CALL VEC_MAT_SYMM_GRN(E_BAR_VEC, E_BAR, N, M)
    !
    GDOT = 0.0_RK
    DGDOT = 0.0_RK
    !      
    ! Begin iterations
    !   
    DO ITER_NEWTON = 1, MAX_ITER_NEWTON
        !
        SIG_0 = SIG
        !
        ! Elastic strains
        !
        DO IPHASE=1,NUMPHASES
            !
            CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
            CALL CRYSTALTYPEGET(CTYPE(IPHASE), DEV=P_HAT_VEC,PPTRANS=PPT)
            N_SLIP=CTYPE(IPHASE)%NUMSLIP
            !
            IF(ASSOCIATED(E_VEC_TMP)) THEN
                !
                DEALLOCATE(E_VEC_TMP)
                !
            ENDIF
            !
            ALLOCATE(E_VEC_TMP(0:TVEC1,0:(N-1),0:(NUMIND-1)))
            CALL VEC_D_VEC5(KEINV(:,IPHASE), SIG(:,:,INDICES), E_VEC_TMP, &
                & N, NUMIND)
            !
            E_VEC(:,:,INDICES)=E_VEC_TMP
            !
            IF(ASSOCIATED(E_ELAS_TMP)) THEN
                !
                DEALLOCATE(E_ELAS_TMP)
                !            
            ENDIF
            !
            ALLOCATE(E_ELAS_TMP(0:DIMS1,0:DIMS1,0:(N-1),0:(NUMIND-1)))
            CALL VEC_MAT_SYMM_GRN(E_VEC(:,:,INDICES), E_ELAS_TMP, N, NUMIND)
            !            
            E_ELAS(:,:,:,INDICES) = E_ELAS_TMP
            !
            ! Power law viscoplastic model
            !
            ! Added parameter check for anisotropic rate sensitivity
            IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .FALSE.) THEN
                !
                XNN(IPHASE) = 1.0_RK / CRYSTAL_PARM(0,IPHASE)
                XN(IPHASE)  = XNN(IPHASE) - 1.0_RK
                TOOSMALL(IPHASE) = UFLOW**CRYSTAL_PARM(0,IPHASE)
                !
            ELSE IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .TRUE.) THEN
                !
                ANISO_M_TEMP(0:2)  = CRYS_OPTIONS%ANISO_M(IPHASE,1)
                ANISO_M_TEMP(3:5)  = CRYS_OPTIONS%ANISO_M(IPHASE,2)
                ANISO_M_TEMP(6:17) = CRYS_OPTIONS%ANISO_M(IPHASE,3)
                !
                AXNN(:,IPHASE) = 1.0_RK / ANISO_M_TEMP(:)
                AXN(:,IPHASE)  = AXNN(:,IPHASE) - 1.0_RK
                ATOOSMALL(:,IPHASE) = UFLOW**ANISO_M_TEMP(:)
                !
            END IF
            !
            DO ISLIP = 0, N_SLIP - 1
                !         
                TAU(:,INDICES) = 0.0_RK
                !                
                DO I = 0, TVEC1
                    !
                    TAU(:,INDICES) = TAU(:,INDICES) + P_HAT_VEC(I+1,ISLIP+1) * &
                        & SIG(I,:,INDICES)/CRSS(ISLIP,:,INDICES)
                    !
                ENDDO 
                !
                TAUA(:,INDICES) = DABS(TAU(:,INDICES))
                !
                IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .FALSE.) THEN
                    !
                    WHERE (TAUA(:,INDICES) .LE. TOOSMALL(IPHASE)) &
                        & TAUA(:,INDICES) = 0.0_RK
                    !
                    FJ21(:,INDICES) = &
                        &CRYSTAL_PARM(1,IPHASE) * TAUA(:,INDICES)**XN(IPHASE)  
                    DGDOT(ISLIP, :,INDICES) = &
                        &FJ21(:,INDICES) * XNN(IPHASE) / CRSS(ISLIP,:,INDICES)
                    GDOT(ISLIP, :, INDICES) = FJ21(:,INDICES) * TAU(:,INDICES) 
                    !
                ELSE IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .TRUE.) THEN
                    !
                    WHERE (TAUA(:,INDICES) .LE. ATOOSMALL(ISLIP,IPHASE)) &
                        & TAUA(:,INDICES) = 0.0_RK
                    !
                    FJ21(:,INDICES) = CRYSTAL_PARM(1,IPHASE) * &
                        &TAUA(:,INDICES)**AXN(ISLIP,IPHASE)  
                    DGDOT(ISLIP, :,INDICES) = FJ21(:,INDICES) * &
                        &AXNN(ISLIP,IPHASE) / CRSS(ISLIP,:,INDICES)
                    GDOT(ISLIP, :, INDICES) = FJ21(:,INDICES) * TAU(:,INDICES)
                    !
                END IF
                !
            ENDDO !N_SLIP
            !
            ! Set up the Jacobian
            !        
            CALL DP_WP_HAT(P_HAT_VEC, DP_HAT_VEC, &
                & WP_HAT_VEC, E_ELAS, E_BAR, W_VEC_LAT, GDOT, &
                & N_SLIP, DT, N, M, NUMIND, INDICES)
            !
            CALL WP_HAT_MAT5X5(WP_HAT_VEC,WP_HAT_MATX, N, M, NUMIND, INDICES)
            !
            IF(ASSOCIATED(WP_X_E_TMP)) THEN
                !
                DEALLOCATE(WP_X_E_TMP)
                !
            ENDIF
            !
            ALLOCATE(WP_X_E_TMP(0:TVEC1, 0:(N - 1), 0:(NUMIND - 1)))
            CALL MAT_X_VEC5(WP_HAT_MATX(:, :, :, INDICES), &
                & E_VEC(:, :, INDICES), WP_X_E_TMP, N, NUMIND)
            !            
            WP_X_E(:, :, INDICES)=WP_X_E_TMP
            !
            IF(ASSOCIATED(FJAC_TMP)) THEN
                !               
                DEALLOCATE(FJAC_TMP)
                !
            ENDIF
            !
            ALLOCATE(FJAC_TMP(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(NUMIND-1)))
            CALL MATRIX_FJAC(FJAC_TMP,&
                & KEINV(:,IPHASE), WP_HAT_MATX(:, :, :, INDICES), C1, &
                & DGDOT(:, :, INDICES), PPT, N_SLIP,&
                & N, NUMIND)
            !
            FJAC(:, :, :, INDICES)=FJAC_TMP
            !
            ! Set up the system function (RHS)      
            !
            CALL RESIDUAL(RES_N, DEL_S, &
                & D_VEC_LAT, E_VEC, E_BAR_VEC, &
                & DP_HAT_VEC, WP_X_E, C1, N, M, NUMIND, INDICES)
            !
            DEALLOCATE(INDICES)
            DEALLOCATE(FJAC_TMP)
            DEALLOCATE(PPT)
            DEALLOCATE(WP_X_E_TMP)
            DEALLOCATE(E_VEC_TMP)
            DEALLOCATE(E_ELAS_TMP)
            DEALLOCATE(P_HAT_VEC)
            !
        ENDDO !NUMPHASES
        !
        ! Compute new iteration
        !
        CALL SYMMETRIZE_JAC(FJAC, DEL_S, N, M)
        CALL SOLVIT(FJAC, DEL_S, N, M)
        CALL CHECK_DIAGONALS_EVPS(FJAC, NEWTON_OK, DONE, N, M)
        !
        ! Trial residual
        !
        XLAMBDA = SIG_0 + DEL_S
        !
        DO IPHASE = 1, NUMPHASES
            !
            CALL CRYSTALTYPEGET(CTYPE(IPHASE), DEV=P_HAT_VEC)
            CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
            N_SLIP=CTYPE(IPHASE)%NUMSLIP
            !            
            IF(ASSOCIATED(E_VEC_TMP)) THEN
                !
                DEALLOCATE(E_VEC_TMP)
                !            
            ENDIF
            !
            ALLOCATE(E_VEC_TMP(0:TVEC1,0:(N-1),0:(NUMIND-1)))
            CALL VEC_D_VEC5(KEINV(:,IPHASE), XLAMBDA(:,:,INDICES), E_VEC_TMP, &
                & N, NUMIND)
            !
            E_VEC(:,:,INDICES) = E_VEC_TMP
            !
            IF(ASSOCIATED(E_ELAS_TMP)) THEN
                !
                DEALLOCATE(E_ELAS_TMP)
                !
            ENDIF
            !
            ALLOCATE(E_ELAS_TMP(0:DIMS1,0:DIMS1,0:(N-1),0:(NUMIND-1)))
            CALL VEC_MAT_SYMM_GRN(E_VEC(:, :, INDICES),E_ELAS_TMP, N, NUMIND)
            !            
            E_ELAS(:,:,:,INDICES) = E_ELAS_TMP
            !
            ! Added parameter check for anisotropic rate sensitivity
            DO ISLIP = 0, N_SLIP - 1
                !
                TAU(:,INDICES) = 0.0_RK
                !               
                DO I = 0, TVEC1
                    !
                    TAU(:,INDICES) = TAU(:,INDICES) + P_HAT_VEC(I+1,ISLIP+1) * &
                        & XLAMBDA(I,:,INDICES)/CRSS(ISLIP,:,INDICES)
                    !
                ENDDO
                !                
                TAUA(:,INDICES) = DABS(TAU(:,INDICES))
                !                
                DO IN = 0, (N-1)
                    !
                    IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .FALSE.) THEN
                        !
                        WHERE (TAUA(IN,INDICES) .LE. TOOSMALL(IPHASE)) &
                            & TAUA(IN,INDICES) = 0.0_RK
                        GDOT(ISLIP, IN, :) = CRYSTAL_PARM(1,IPHASE)*&
                            &TAU(IN,:)*TAUA(IN,:)**XN(IPHASE) 
                        !
                    ELSE IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .TRUE.) THEN
                        !
                        WHERE (TAUA(IN,INDICES).LE.ATOOSMALL(ISLIP,IPHASE)) &
                            & TAUA(IN,INDICES) = 0.0_RK
                        GDOT(ISLIP, IN, :) = CRYSTAL_PARM(1,IPHASE)*&
                            &TAU(IN,:)*TAUA(IN,:)**AXN(ISLIP,IPHASE) 
                        !
                    ENDIF
                    !
                ENDDO
                ! 
            ENDDO !N_SLIP
            !    
            CALL DP_WP_HAT(P_HAT_VEC, DP_HAT_VEC, &
                & WP_HAT_VEC, E_ELAS, E_BAR, W_VEC_LAT, GDOT, &
                & N_SLIP, DT, N, M, NUMIND, INDICES)
            !
            CALL WP_HAT_MAT5X5(WP_HAT_VEC,WP_HAT_MATX, N, M, NUMIND, INDICES)
            !
            IF(ASSOCIATED(WP_X_E_TMP)) THEN
                !
                DEALLOCATE(WP_X_E_TMP)
                !
            ENDIF
            !
            ALLOCATE(WP_X_E_TMP(0:TVEC1, 0:(N - 1), 0:(NUMIND - 1)))
            CALL MAT_X_VEC5(WP_HAT_MATX(:, :, :, INDICES), &
                & E_VEC(:, :, INDICES), WP_X_E_TMP, N, NUMIND)
            !            
            WP_X_E(:, :, INDICES)=WP_X_E_TMP
            !
            CALL RESIDUAL(RES, FJ,&
                & D_VEC_LAT, E_VEC, E_BAR_VEC,&
                & DP_HAT_VEC, WP_X_E, C1, N, M, NUMIND, INDICES)
            !
            DEALLOCATE(INDICES)
            DEALLOCATE(WP_X_E_TMP)
            DEALLOCATE(E_VEC_TMP)
            DEALLOCATE(E_ELAS_TMP)
            DEALLOCATE(P_HAT_VEC)
            !         
        ENDDO !NUMPHASES
        !
        ! Line Search.
        !
        FACT = 1.0_RK
        RATIO_RES = RES / RES_N
        !
        DO WHILE(ANY(RATIO_RES .GT. 1.0 .AND. NEWTON_OK .AND. .NOT. CONVERGED))
            !
            WHERE(RATIO_RES .GT. 1.0 .AND. NEWTON_OK .AND. .NOT. CONVERGED) &
                & FACT = FACT*0.5_RK
            !
            IF(ANY(FACT .LT. 0.001)) THEN
                !
                ! DEBUG: Update this print - JC
                print *, ' error in line search',count(fact .lt. 0.1), iter_newton
                WHERE(FACT .LT. 0.001) NEWTON_OK = .FALSE.
                !
            ENDIF
            !
            XLAMBDA(0, :, :) = SIG_0(0, :, :) + FACT * DEL_S(0, :, :)
            XLAMBDA(1, :, :) = SIG_0(1, :, :) + FACT * DEL_S(1, :, :)
            XLAMBDA(2, :, :) = SIG_0(2, :, :) + FACT * DEL_S(2, :, :)
            XLAMBDA(3, :, :) = SIG_0(3, :, :) + FACT * DEL_S(3, :, :)
            XLAMBDA(4, :, :) = SIG_0(4, :, :) + FACT * DEL_S(4, :, :)      
            !
            DO IPHASE = 1, NUMPHASES
                !
                CALL CRYSTALTYPEGET(CTYPE(IPHASE), DEV=P_HAT_VEC)
                CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
                N_SLIP=CTYPE(IPHASE)%NUMSLIP
                !
                IF(ASSOCIATED(E_VEC_TMP)) THEN
                    !                  
                    DEALLOCATE(E_VEC_TMP)
                    !                
                ENDIF
                !
                ALLOCATE(E_VEC_TMP(0:TVEC1,0:(N-1),0:(NUMIND-1)))
                CALL VEC_D_VEC5(KEINV(:,IPHASE), XLAMBDA(:,:,INDICES), E_VEC_TMP, N, NUMIND)
                !                
                E_VEC(:,:,INDICES) = E_VEC_TMP
                !
                IF(ASSOCIATED(E_ELAS_TMP)) THEN
                    !                  
                    DEALLOCATE(E_ELAS_TMP)
                    !                
                ENDIF
                !
                ALLOCATE(E_ELAS_TMP(0:DIMS1,0:DIMS1,0:(N-1),0:(NUMIND-1)))
                CALL VEC_MAT_SYMM_GRN(E_VEC(:,:,INDICES),E_ELAS_TMP, N,NUMIND)
                !                
                E_ELAS(:,:,:,INDICES)=E_ELAS_TMP
                !
                DO ISLIP = 0, N_SLIP - 1
                    !
                    TAU(:,INDICES) = 0.0_RK
                    !                  
                    DO I = 0, TVEC1
                        !
                        TAU(:,INDICES)=TAU(:,INDICES)+P_HAT_VEC(I+1,ISLIP+1)*&
                            & XLAMBDA(I, :, INDICES) / CRSS(ISLIP,:,INDICES)
                        !
                    ENDDO
                    ! 
                    TAUA(:,INDICES) = DABS(TAU(:,INDICES))
                    !
                    DO IN = 0, N-1
                        !
                        ! Added parameter check for anisotropic rate sensitivity
                        IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .FALSE.) THEN
                            !
                            WHERE(TAUA(IN,INDICES).LE.TOOSMALL(IPHASE)) &
                                & TAUA(IN,INDICES) = 0.0
                            GDOT(ISLIP, IN, :) = CRYSTAL_PARM(1,IPHASE)*&
                                &TAU(IN,:)*TAUA(IN,:)**XN(IPHASE)
                            !
                        ELSE IF &
                        & (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .TRUE.) THEN
                            !
                            WHERE(TAUA(IN,INDICES).LE.ATOOSMALL(ISLIP,IPHASE)) &
                                & TAUA(IN,INDICES) = 0.0
                            GDOT(ISLIP, IN, :) = CRYSTAL_PARM(1,IPHASE)*&
                                &TAU(IN,:)*TAUA(IN,:)**AXN(ISLIP,IPHASE)
                            !
                        END IF
                        !
                    ENDDO
                    !
                ENDDO !N_SLIP
                !
                CALL DP_WP_HAT(P_HAT_VEC, DP_HAT_VEC, &
                    & WP_HAT_VEC, E_ELAS, E_BAR, W_VEC_LAT, GDOT, &
                    & N_SLIP, DT, N, M, NUMIND, INDICES)
                !
                CALL WP_HAT_MAT5X5(WP_HAT_VEC,WP_HAT_MATX, N, M, NUMIND, &
                    & INDICES)
                !
                IF(ASSOCIATED(WP_X_E_TMP)) THEN
                    !
                    DEALLOCATE(WP_X_E_TMP)
                    !
                ENDIF
                !
                ALLOCATE(WP_X_E_TMP(0:TVEC1, 0:(N - 1), 0:(NUMIND - 1)))
                CALL MAT_X_VEC5(WP_HAT_MATX(:, :, :, INDICES), &
                    & E_VEC(:, :, INDICES), WP_X_E_TMP, N, NUMIND)
                !
                WP_X_E(:, :, INDICES) = WP_X_E_TMP
                !
                CALL RESIDUAL(RES_AUX, FJ,&
                    & D_VEC_LAT, E_VEC, E_BAR_VEC, DP_HAT_VEC, WP_X_E,&
                    & C1, N, M, NUMIND, INDICES)
                !
                DEALLOCATE(INDICES)
                DEALLOCATE(WP_X_E_TMP)
                DEALLOCATE(E_VEC_TMP)
                DEALLOCATE(E_ELAS_TMP)
                DEALLOCATE(P_HAT_VEC)
                !
            ENDDO !NUMPHASES
            !
            WHERE((RATIO_RES .GT. 1.0) .AND. (NEWTON_OK) .AND. (.NOT. CONVERGED))
                !
                RES = RES_AUX
                RATIO_RES = RES / RES_N
                !
            ENDWHERE
            !
        ENDDO !DO WHILE
        !
        ! Update stresses and check convergence.
        !
        WHERE((NEWTON_OK) .AND. (.NOT. CONVERGED))
            !
            SIG(0, :, :) = SIG_0(0, :, :) + FACT * DEL_S(0, :, :)
            SIG(1, :, :) = SIG_0(1, :, :) + FACT * DEL_S(1, :, :)
            SIG(2, :, :) = SIG_0(2, :, :) + FACT * DEL_S(2, :, :)
            SIG(3, :, :) = SIG_0(3, :, :) + FACT * DEL_S(3, :, :)
            SIG(4, :, :) = SIG_0(4, :, :) + FACT * DEL_S(4, :, :)
            !
        ENDWHERE
        !
        WHERE((RES .LE. CV_OPTIONS % SX_TOL).AND.(NEWTON_OK).AND.(.NOT. DONE)&
            & .AND.(.NOT. CONVERGED)) JITER = ITER_NEWTON
        WHERE((RES .LE. CV_OPTIONS % SX_TOL).AND.(NEWTON_OK).AND.(.NOT. DONE))&
            & CONVERGED = .TRUE.
        !
        ! Return if all grains have converged.
        !
        INEWTON = COUNT(.NOT. NEWTON_OK)
        !
        IF ((COUNT(CONVERGED) + INEWTON) .EQ. NM) THEN
            !
            if (INEWTON .GT. 0) THEN
                !
                WRITE(DFLT_U, '(A,I0,A,I0)') 'Warning:       . SOLVE_NEWTON_EVPS:&
                    & Converged = ', COUNT(CONVERGED), ' Remaining = ', INEWTON
                WRITE(DFLT_U, '(A,G12.5,2X,G12.5)') '                 Residual for&
                    & converged grains = ', MINVAL(RES, MASK=CONVERGED), &
                    & MAXVAL(RES, MASK=CONVERGED)
                !
            END IF
            !
            IRC = INEWTON
            !
            RETURN
            !
        ENDIF
        ! 
    ENDDO
    !
    IRC = -2
    !
    RETURN
    !
    END SUBROUTINE SOLVE_NEWTON_EVPS
    !
    !===========================================================================
    !
    SUBROUTINE MATRIX_FJAC(FJAC, KEINV, W_MATX, DTI, DGDOT, PPT, N_SLIP, N, M)
    !
    ! Add descriptions here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER, INTENT(IN)   :: N_SLIP, N, M
    !
    REAL(RK), INTENT(OUT) :: FJAC(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: PPT(0:TVEC1, 0:TVEC1, 0:MAXSLIP1)
    REAL(RK), INTENT(IN)  :: DTI
    REAL(RK), INTENT(IN)  :: KEINV(0:TVEC1)
    REAL(RK), INTENT(IN)  :: W_MATX(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: DGDOT(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    INTEGER :: ISLIP, I, J, IN, IM
    !
    !---------------------------------------------------------------------------
    !
    DO IM = 0, M-1
        !    
        DO IN = 0, N-1
            !
            DO J = 0, TVEC1
                !               
                DO I = 0, TVEC1          
                    !                
                    FJAC(I, J, IN, IM) = W_MATX(I, J, IN, IM) * KEINV(J)
                    !               
                ENDDO
                !
                FJAC(J, J, IN, IM) = FJAC(J, J, IN, IM) + KEINV(J) * DTI
                !
            ENDDO
            !
            DO ISLIP = 0, N_SLIP-1
                !            
                FJAC(:,:,IN,IM) = FJAC(:,:,IN,IM) + &
                    & DGDOT(ISLIP,IN,IM) * PPT(:,:,ISLIP)
                !            
            ENDDO
            !
        ENDDO
        !      
    ENDDO
    !
    RETURN
    !
    END SUBROUTINE MATRIX_FJAC
    !
    !===========================================================================
    !
    SUBROUTINE JACOB5X5(JAC, A, N, M)
    !
    ! Computes the 5x5 Jacobian for matrix `a'. Doesn't appear to be used?
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER, INTENT(IN)   :: N, M
    !
    REAL(RK), INTENT(IN)  :: A(0:DIMS1, 0:DIMS1, 0:DIMS1, &
        & 0:DIMS1, 0:(N-1), 0:(M-1))
    REAL(RK), INTENT(OUT) :: JAC(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    REAL(RK) :: SQR3, SQR3B2
    !
    !---------------------------------------------------------------------------
    !
    SQR3   = DSQRT(3.D0)
    SQR3B2 = SQR3 / 2.D0 
    !
    JAC(0, 0, :, :) = (A(0, 0, 0, 0, :, :) + A(1, 1, 1, 1, :, :) - &
        & A(0, 0, 1, 1, :, :) - A(1, 1, 0, 0, :, :)) / 2.
    JAC(0, 1, :, :) = SQR3B2 * (A(0, 0, 2, 2, :, :) - &
        & A(1, 1, 2, 2, :, :))
    JAC(0, 2, :, :) = A(0, 0, 0, 1, :, :) - A(1, 1, 0, 1, :, :)
    JAC(0, 3, :, :) = A(0, 0, 0, 2, :, :) - A(1, 1, 0, 2, :, :)
    JAC(0, 4, :, :) = A(0, 0, 1, 2, :, :) - A(1, 1, 1, 2, :, :)
    !
    JAC(1, 0, :, :) = SQR3B2 * (A(2, 2, 0, 0, :, :) -&
        & A(2, 2, 1, 1, :, :))
    JAC(1, 1, :, :) = 1.5 * A(2, 2, 2, 2, :, :) - &
        & 0.5*(A(2, 2, 0, 0, :, :) + A(2, 2, 1, 1, :, :) + A(2, 2, 2, 2, :, :))
    JAC(1, 2, :, :) = SQR3 * A(2, 2, 0, 1, :, :)
    JAC(1, 3, :, :) = SQR3 * A(2, 2, 0, 2, :, :)
    JAC(1, 4, :, :) = SQR3 * A(2, 2, 1, 2, :, :)
    !
    JAC(2, 0, :, :) = A(0, 1, 0, 0, :, :) - A(0, 1, 1, 1, :, :)
    JAC(2, 1, :, :) = SQR3 * A(0, 1, 2, 2, :, :)
    JAC(2, 2, :, :) = 2. * A(0, 1, 0, 1, :, :)
    JAC(2, 3, :, :) = 2. * A(0, 1, 0, 2, :, :)
    JAC(2, 4, :, :) = 2. * A(0, 1, 1, 2, :, :)
    !
    JAC(3, 0, :, :) = A(0, 2, 0, 0, :, :) - A(0, 2, 1, 1, :, :)
    JAC(3, 1, :, :) = SQR3 * A(0, 2, 2, 2, :, :)
    JAC(3, 2, :, :) = 2. * A(0, 2, 0, 1, :, :)
    JAC(3, 3, :, :) = 2. * A(0, 2, 0, 2, :, :)
    JAC(3, 4, :, :) = 2. * A(0, 2, 1, 2, :, :)
    !
    JAC(4, 0, :, :) = A(1, 2, 0, 0, :, :) - A(1, 2, 1, 1, :, :)
    JAC(4, 1, :, :) = SQR3 * A(1, 2, 2, 2, :, :)
    JAC(4, 2, :, :) = 2. * A(1, 2, 0, 1, :, :)
    JAC(4, 3, :, :) = 2. * A(1, 2, 0, 2, :, :)
    JAC(4, 4, :, :) = 2. * A(1, 2, 1, 2, :, :)
    !
    RETURN
    !
    END SUBROUTINE JACOB5X5
    !
    !===========================================================================
    !
    SUBROUTINE DP_WP_HAT(P_HAT_VEC, DP_HAT, WP_HAT, E_ELAS, E_BAR, W_VEC_LAT, &
        & GDOT, N_SLIP, DT, N, M, NUMIND, INDICES)
    !
    ! Add descriptions here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER, INTENT(IN)   :: N_SLIP, N, M, NUMIND, INDICES(1:NUMIND)
    REAL(RK), INTENT(OUT) :: DP_HAT(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: WP_HAT(0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: P_HAT_VEC(0:TVEC1,0:MAXSLIP1)
    REAL(RK), INTENT(IN)  :: DT
    REAL(RK), INTENT(IN)  :: E_ELAS(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: E_BAR(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: W_VEC_LAT(0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: GDOT(0:MAXSLIP1, 0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    INTEGER  :: I, ISLIP
    REAL(RK) :: DP_HAT_TMP(0:TVEC1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: WP_HAT_TMP(0:DIMS1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: E_ELAS_TMP(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: E_BAR_TMP(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: W_VEC_LAT_TMP(0:DIMS1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: GDOT_TMP(0:MAXSLIP1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: P_HAT(0:DIMS1, 0:DIMS1, 0:MAXSLIP1)
    REAL(RK) :: X (0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(NUMIND - 1))
    REAL(RK) :: EE(0:DIMS1, 0:DIMS1, 0:(N - 1), 0:(NUMIND - 1))
    !
    !---------------------------------------------------------------------------
    !    
    E_ELAS_TMP=E_ELAS(:, :, :, INDICES)
    E_BAR_TMP=E_BAR(:, :, :, INDICES)
    W_VEC_LAT_TMP=W_VEC_LAT(:, :, INDICES)
    GDOT_TMP=GDOT(:, :, INDICES)
    !    
    CALL VEC_MAT_SYMM(P_HAT_VEC, P_HAT, N_SLIP)
    DP_HAT_TMP = 0.0_RK
    !
    CALL MAT_X_MAT3(E_ELAS_TMP, E_BAR_TMP, EE, N, NUMIND)
    !
    WP_HAT_TMP(0, :, :) = W_VEC_LAT_TMP(0, :, :) + &
        & 0.5 / DT * (EE(1, 0, :, :) - EE(0, 1, :, :))
    WP_HAT_TMP(1, :, :) = W_VEC_LAT_TMP(1, :, :) + &
        & 0.5 / DT * (EE(2, 0, :, :) - EE(0, 2, :, :))
    WP_HAT_TMP(2, :, :) = W_VEC_LAT_TMP(2, :, :) + &
        & 0.5 / DT * (EE(2, 1, :, :) - EE(1, 2, :, :))
    !
    DO ISLIP = 0, N_SLIP - 1
        !
        CALL MAT_X_MATS3(E_ELAS_TMP, P_HAT(0, 0, ISLIP), X, N, NUMIND)
        !
        WP_HAT_TMP(0, :, :) = WP_HAT_TMP(0, :, :) - &
            & GDOT_TMP(ISLIP, :, :) * (X(1, 0, :, :) - X(0, 1, :, :))
        WP_HAT_TMP(1, :, :) = WP_HAT_TMP(1, :, :) - &
            & GDOT_TMP(ISLIP, :, :) * (X(2, 0, :, :) - X(0, 2, :, :))
        WP_HAT_TMP(2, :, :) = WP_HAT_TMP(2, :, :) - &
            & GDOT_TMP(ISLIP, :, :) * (X(2, 1, :, :) - X(1, 2, :, :))
        !
        DO I = 0, TVEC1
            !        
            DP_HAT_TMP(I, :, :) = DP_HAT_TMP(I, :, :) + &
                & GDOT_TMP(ISLIP, :, :) * P_HAT_VEC(I, ISLIP)
            !
        ENDDO
        !
    ENDDO
    !
    DP_HAT(:,:,INDICES) = DP_HAT_TMP
    WP_HAT(:, :, INDICES) = WP_HAT_TMP
    !
    RETURN
    !
    END SUBROUTINE DP_WP_HAT
    !
    !===========================================================================
    !
    SUBROUTINE WP_HAT_MAT5X5(WP_HAT, WP_HAT_MATX, N, M, NUMIND, INDICES)
    !
    ! Add descriptions here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER, INTENT(IN)   :: N, M, NUMIND, INDICES(1:NUMIND)
    REAL(RK), INTENT(IN)  :: WP_HAT(0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: WP_HAT_MATX(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    INTEGER  :: I, J
    REAL(RK) :: SQR3
    !
    !---------------------------------------------------------------------------
    !      
    SQR3 = DSQRT(3.D0)
    !
    WP_HAT_MATX(:,:,:,INDICES) = 0.0
    !
    WP_HAT_MATX(0, 2, :, INDICES) =   WP_HAT(0, :, INDICES) * 2.0D0
    WP_HAT_MATX(0, 3, :, INDICES) =   WP_HAT(1, :, INDICES)
    WP_HAT_MATX(0, 4, :, INDICES) = - WP_HAT(2, :, INDICES)
    !                                          
    WP_HAT_MATX(1, 3, :, INDICES) = - WP_HAT(1, :, INDICES) * SQR3
    WP_HAT_MATX(1, 4, :, INDICES) = - WP_HAT(2, :, INDICES) * SQR3
    !                                          
    WP_HAT_MATX(2, 0, :, INDICES) = - WP_HAT(0, :, INDICES) * 2.0D0
    WP_HAT_MATX(2, 3, :, INDICES) =   WP_HAT(2, :, INDICES)
    WP_HAT_MATX(2, 4, :, INDICES) =   WP_HAT(1, :, INDICES)
    !                                          
    WP_HAT_MATX(3, 0, :, INDICES) = - WP_HAT(1, :, INDICES)
    WP_HAT_MATX(3, 1, :, INDICES) =   WP_HAT(1, :, INDICES) * SQR3
    WP_HAT_MATX(3, 2, :, INDICES) = - WP_HAT(2, :, INDICES)
    WP_HAT_MATX(3, 4, :, INDICES) =   WP_HAT(0, :, INDICES)
    !                                          
    WP_HAT_MATX(4, 0, :, INDICES) =   WP_HAT(2, :, INDICES)
    WP_HAT_MATX(4, 1, :, INDICES) =   WP_HAT(2, :, INDICES) * SQR3
    WP_HAT_MATX(4, 2, :, INDICES) = - WP_HAT(1, :, INDICES)
    WP_HAT_MATX(4, 3, :, INDICES) = - WP_HAT(0, :, INDICES)
    !
    RETURN
    !
    END SUBROUTINE WP_HAT_MAT5X5
    !
    !===========================================================================
    !
    SUBROUTINE WP_HAT_MAT5X5_ALL(WP_HAT, WP_HAT_MATX, N, M)
    !
    ! Add descriptions here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER, INTENT(IN)   :: N, M
    REAL(RK), INTENT(IN)  :: WP_HAT(0:DIMS1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: WP_HAT_MATX(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    INTEGER  :: I, J
    REAL(RK) :: SQR3
    !
    !---------------------------------------------------------------------------
    !
    SQR3 = DSQRT(3.D0)
    !
    WP_HAT_MATX(:,:,:,:) = 0.0
    !
    WP_HAT_MATX(0, 2, :, :) =   WP_HAT(0, :, :) * 2.0D0
    WP_HAT_MATX(0, 3, :, :) =   WP_HAT(1, :, :)
    WP_HAT_MATX(0, 4, :, :) = - WP_HAT(2, :, :)
    !                                      
    WP_HAT_MATX(1, 3, :, :) = - WP_HAT(1, :, :) * SQR3
    WP_HAT_MATX(1, 4, :, :) = - WP_HAT(2, :, :) * SQR3
    !                                      
    WP_HAT_MATX(2, 0, :, :) = - WP_HAT(0, :, :) * 2.0D0
    WP_HAT_MATX(2, 3, :, :) =   WP_HAT(2, :, :)
    WP_HAT_MATX(2, 4, :, :) =   WP_HAT(1, :, :)
    !                                      
    WP_HAT_MATX(3, 0, :, :) = - WP_HAT(1, :, :)
    WP_HAT_MATX(3, 1, :, :) =   WP_HAT(1, :, :) * SQR3
    WP_HAT_MATX(3, 2, :, :) = - WP_HAT(2, :, :)
    WP_HAT_MATX(3, 4, :, :) =   WP_HAT(0, :, :)
    !                                      
    WP_HAT_MATX(4, 0, :, :) =   WP_HAT(2, :, :)
    WP_HAT_MATX(4, 1, :, :) =   WP_HAT(2, :, :) * SQR3
    WP_HAT_MATX(4, 2, :, :) = - WP_HAT(1, :, :)
    WP_HAT_MATX(4, 3, :, :) = - WP_HAT(0, :, :)   
    !
    RETURN
    !
    END SUBROUTINE WP_HAT_MAT5X5_ALL
    !
    !===========================================================================
    !
    SUBROUTINE RESIDUAL(RES, RHS, D_VEC_LAT, E_VEC, E_BAR_VEC, DP_HAT_VEC, &
        & WP_X_E, C1, N, M, NUMIND, INDICES)
    !
    ! Add descriptions here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER, INTENT(IN)   :: N, M, NUMIND, INDICES(1:NUMIND)
    !
    REAL(RK), INTENT(OUT) :: RHS(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(OUT) :: RES(0:(N - 1), 0:(NUMIND - 1))
    REAL(RK), INTENT(IN)  :: C1
    REAL(RK), INTENT(IN)  :: D_VEC_LAT(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: E_VEC(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: E_BAR_VEC(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: DP_HAT_VEC(0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN)  :: WP_X_E(0:TVEC1, 0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    INTEGER :: I
    !
    !---------------------------------------------------------------------------
    !
    RHS = D_VEC_LAT - C1 * (E_VEC - E_BAR_VEC) - DP_HAT_VEC - WP_X_E
    RES(:,INDICES) = 0.0_RK
    !
    DO I = 0, TVEC1
        !         
        RES(:,INDICES) = RES(:,INDICES) + RHS(I, :, INDICES)* RHS(I, :, INDICES)
        !      
    ENDDO
    !
    RES(:,INDICES) = DSQRT(RES(:,INDICES))
    !
    RETURN
    !
    END SUBROUTINE RESIDUAL
    !
    !===========================================================================
    !
    SUBROUTINE CHECK_DIAGONALS_EVPS(STIF, NEWTON_OK, DONE, N, M)
    !
    ! Add descriptions here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER, INTENT(IN)  :: N, M
    !
    LOGICAL, INTENT(OUT) :: NEWTON_OK(0:(N - 1), 0:(M - 1))
    LOGICAL, INTENT(IN)  :: DONE(0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(IN) :: STIF(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    INTEGER :: I
    !
    !---------------------------------------------------------------------------
    !
    DO I = 0, TVEC1
        !
        WHERE ( (ABS(STIF(I, I, :, :)) .LT. VTINY) .AND. (.NOT. DONE) ) &
            & NEWTON_OK = .FALSE.
        !    
    ENDDO
    !
    RETURN
    !
    END SUBROUTINE CHECK_DIAGONALS_EVPS
    !
    !===========================================================================
    !
    SUBROUTINE SYMMETRIZE_JAC(FJAC, DEL_S, N, M)
    !
    ! Add descriptions here.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER, INTENT(IN) :: N, M
    !
    REAL(RK), INTENT(INOUT) :: FJAC(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK), INTENT(INOUT) :: DEL_S(0:TVEC1, 0:(N - 1), 0:(M - 1))
    !
    ! Locals:
    !
    INTEGER  :: I, J, K
    !
    REAL(RK) :: TEMPJ(0:TVEC1, 0:TVEC1, 0:(N - 1), 0:(M - 1))
    REAL(RK) :: TEMPS(0:TVEC1, 0:(N - 1), 0:(M - 1))
    !
    !---------------------------------------------------------------------------
    !
    TEMPJ = 0.0_RK
    TEMPS = 0.0_RK
    !
    ! RC 6/24/2016: Reordered loops for better memory striding      
    DO J = 0, TVEC1
        !
        DO I = 0, TVEC1
            !
            DO K = 0, TVEC1
                !
                TEMPJ(I, J, :, :) = TEMPJ(I, J, :, :) + &
                    & FJAC(K, I, :, :) * FJAC(K, J, :, :)
                !
            ENDDO
            !
        ENDDO
        !
    ENDDO
    !
    DO I = 0, TVEC1
        !
        DO J = 0, TVEC1
            !
            TEMPS(I, :, :) = TEMPS(I, :, :) + FJAC(J, I, :, :) * DEL_S(J, :, :)
            !
        ENDDO
        !
    ENDDO
    !
    FJAC = TEMPJ
    DEL_S = TEMPS
    !
    RETURN
    !
    END SUBROUTINE SYMMETRIZE_JAC
    !
END MODULE STRESS_SOLVE_EVPS_MOD
