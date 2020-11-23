! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE ANISO_EVPS_MOD
!
! Module to compute the anisotropic elasto-viscoplastic solution.
!
! Contains subroutines:
! ANISO_EVPS: Compute the anisotropic elasto-viscoplastic solution.
!
! From libfepx:
!
USE DIMENSIONS_MOD
USE KINEMATICS_MOD, ONLY: FIND_WP_HAT
USE MATERIAL_MATRIX_VP_MOD
USE MATRIX_OPERATIONS_MOD
USE MICROSTRUCTURE_MOD
USE READ_INPUT_MOD
USE STRESS_SOLVE_EVPS_MOD
USE STRESS_SOLVE_VP_MOD
USE UNITS_MOD
!
IMPLICIT NONE
!
PRIVATE
!
PUBLIC :: ANISO_EVPS
!
CONTAINS
    !
    SUBROUTINE ANISO_EVPS(C, C_TAN, F, DETV, C_ANGS, SIG_VEC, CRSS, RSTAR_N, &
        & RSTAR, E_BAR_VEC, WTS, W_VEC, E_ELAS_KK, KEINV, DTIME, NR)
    !
    ! Compute the anisotropic elasto-viscoplastic solution.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    REAL(RK), INTENT(OUT) :: C(TVEC, TVEC, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: C_TAN(TVEC, TVEC, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: F(TVEC, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(OUT) :: DETV(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)  :: DTIME
    REAL(RK), INTENT(IN)  :: KEINV(0:TVEC1,1:NUMPHASES)
    REAL(RK), INTENT(IN)  :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)  :: SIG_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)  :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)  :: RSTAR_N(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)  :: RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)  :: E_BAR_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)  :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)  :: W_VEC(0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN)  :: E_ELAS_KK(EL_SUB1:EL_SUP1)
    LOGICAL,  INTENT(IN)  :: NR
    !
    ! Locals:
    !
    INTEGER  :: ISLIP, I, J, M_EL, K, N_SLIP
    INTEGER, POINTER :: INDICES(:) => NULL()
    !
    REAL(RK) :: DTIMEI, SQR2, SQR32
    REAL(RK) :: ALPHA(TVEC)
    REAL(RK) :: E_VEC(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: D_RSTAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: QR5X5(0:TVEC1, 0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: E_BAR_VEC_R(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: E_BAR_SM(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: E_VEC_SM(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: E_BAR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: E_ELAS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: W_VEC_GRN(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    !           
    REAL(RK) :: STIF_VP(0:TVEC1, 0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: STIF_EVP(0:TVEC1, 0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TAN_STIF_VP(0:TVEC1, 0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TAN_STIF_EVP(0:TVEC1, 0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: KEINV_ALL(0:TVEC1, 0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    !            
    REAL(RK) :: RSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: GDOT(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: COMP(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: TEMP(0:TVEC1, 0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: V_TENSOR(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: E_ELAS_KK_GRN(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: DETERM_V(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: WP_HAT_VEC(0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: WP_HAT_MAT(0:TVEC1, 0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: WP_X_E(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: F_ELAS(0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    !
    REAL(RK), POINTER :: E_VEC_TMP(:,:,:) => NULL()
    REAL(RK), POINTER :: P_HAT_VEC(:,:) => NULL()
    !
    INTEGER  :: IPHASE, NUMIND
    INTEGER  :: MY_PHASE(0:(EL_SUP1-EL_SUB1))
    REAL(RK) :: ANISO_M_TEMP(0:17)
    REAL(RK) :: ANISO_M_MIN
    !
    !---------------------------------------------------------------------------
    !
    MY_PHASE(:) = PHASE(EL_SUB1:EL_SUP1)
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    TAN_STIF_EVP = 0.0D0
    C_TAN = 0.0D0
    !
    ! Initialize factors to correct [c] & {f}
    !
    SQR2  = DSQRT(2.0D0)
    SQR32 = DSQRT(1.5D0)
    !
    ALPHA(1) = 1.0D0 / SQR2
    ALPHA(2) = SQR32
    ALPHA(3) = SQR2
    ALPHA(4) = SQR2
    ALPHA(5) = SQR2
    !
    ! Compute elastic strain: sig_vec --> e_vec
    !
    DO IPHASE=1,NUMPHASES
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
        INDICES=INDICES+EL_SUB1
        !
        IF (ASSOCIATED(E_VEC_TMP)) THEN
            !
            DEALLOCATE(E_VEC_TMP)
            !        
        ENDIF
        !
        ALLOCATE(E_VEC_TMP(0:TVEC1,0:NGRAIN1,0:NUMIND-1))
        !        
        CALL VEC_D_VEC5(KEINV(:,IPHASE), SIG_VEC(:,:,INDICES), E_VEC_TMP, &
            & NGRAIN, NUMIND)
        !
        E_VEC(:,:,INDICES)=E_VEC_TMP
        !        
        DEALLOCATE(INDICES)
        DEALLOCATE(E_VEC_TMP)
        !    
    ENDDO !NUMPHASES
    !
    ! Rotate e_bar_vec to current configuration (lattice axes).
    !
    CALL MATT_X_MAT3(RSTAR_N, RSTAR, D_RSTAR, NGRAIN, M_EL)
    !
    CALL ROT_MAT_SYMM(D_RSTAR, QR5X5, NGRAIN, M_EL)
    !
    CALL LATTICE_DEFORM(QR5X5, E_BAR_VEC, E_BAR_VEC_R, NGRAIN, M_EL)
    !
    ! Transform e_bar_vec_r & e_vec to sample axes: {.} = [Q] {.}.
    !
    CALL ROT_MAT_SYMM(C_ANGS, QR5X5, NGRAIN, M_EL)
    !
    CALL MAT_X_VEC5(QR5X5, E_BAR_VEC_R, E_BAR_SM, NGRAIN, M_EL)
    !
    CALL MAT_X_VEC5(QR5X5, E_VEC, E_VEC_SM, NGRAIN, M_EL)
    !
    ! Tensor form of e_bar_sm, e_vec_sm: {.} -> [ ].
    !
    CALL VEC_MAT_SYMM_GRN(E_BAR_SM, E_BAR, NGRAIN, M_EL)
    !
    CALL VEC_MAT_SYMM_GRN(E_VEC_SM, E_ELAS, NGRAIN, M_EL)
    !
    ! {w}_sm to all grains in aggregate (Taylor Assumption).
    !
    DO I = 0, DIMS1
        !
        W_VEC_GRN(I, :, :) = SPREAD(W_VEC(I, :), DIM = 1,  NCOPIES = NGRAIN)
        !    
    ENDDO
    !
    ! Compute elasto-visco-plastic crystal stiffness and 
    !   visco-plastic crystal compliance.
    !
    TAN_STIF_VP = 0.0D0
    !
    DTIMEI = 1. / DTIME
    !
    DO IPHASE=1,NUMPHASES
        !
        CALL CRYSTALTYPEGET(CTYPE(IPHASE), DEV=P_HAT_VEC)
        N_SLIP=CTYPE(IPHASE)%NUMSLIP
        !
        CALL FIND_INDICES(NUMIND, IPHASE, MY_PHASE, INDICES)
        INDICES=INDICES+EL_SUB1
        !
        DO ISLIP = 0, N_SLIP - 1
            !
            ! Compute the RSS
            CALL SS_PROJECT(RSS(ISLIP,:,:),P_HAT_VEC(:,ISLIP+1), SIG_VEC, &
                & NGRAIN, M_EL, NUMIND, INDICES-EL_SUB1)
            !
            RSS(ISLIP,:,INDICES)=RSS(ISLIP,:,INDICES) /CRSS(ISLIP,:,INDICES)
            !
            WHERE (ABS(RSS(ISLIP,:,INDICES)) .LT. T_MIN(IPHASE))
                !
                RSS(ISLIP,:,INDICES) = 0.0D0
                !
            ENDWHERE
            !
            ! Introduced `anisotropic' strain rate sensitivity option here
            ! for HCP material phases - JC
            !
            IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .FALSE.) THEN
                !
                CALL POWER_LAW(GDOT(ISLIP, :, :), RSS(ISLIP, :, :),&
                    & CRYSTAL_PARM(0,IPHASE), CRYSTAL_PARM(1,IPHASE),&
                    & T_MIN(IPHASE), NGRAIN, M_EL, NUMIND, INDICES-EL_SUB1)
                !
                CALL COMPLIANCE(COMP, RSS(ISLIP, :, :),&
                    & GDOT(ISLIP, :, :), CRSS(ISLIP,:,:), &
                    & CRYSTAL_PARM(0,IPHASE), T_MIN(IPHASE), NGRAIN, M_EL, &
                    & NUMIND, INDICES-EL_SUB1)
                !
            ELSE IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .TRUE.) THEN
                !
                ANISO_M_TEMP(0:2)  = CRYS_OPTIONS%ANISO_M(IPHASE,1)
                ANISO_M_TEMP(3:5)  = CRYS_OPTIONS%ANISO_M(IPHASE,2)
                ANISO_M_TEMP(6:17) = CRYS_OPTIONS%ANISO_M(IPHASE,3)
                !
                CALL POWER_LAW(GDOT(ISLIP, :, :), RSS(ISLIP, :, :), &
                    & ANISO_M_TEMP(ISLIP), CRYSTAL_PARM(1,IPHASE), &
                    & T_MIN(IPHASE), NGRAIN, M_EL, NUMIND, INDICES-EL_SUB1)
                !
                CALL COMPLIANCE(COMP, RSS(ISLIP, :, :),&
                    & GDOT(ISLIP, :, :), CRSS(ISLIP,:,:), ANISO_M_TEMP(ISLIP), &
                    & T_MIN(IPHASE), NGRAIN, M_EL, NUMIND, INDICES-EL_SUB1)
                !
            END IF
            !
            ! 
            ! In order to remove one dimension from my_phase (the dimension over
            ! `ngrain'), we add a loop over that dimension here. - hritz 9/15/05
            !
            ! Reordered loop ordering for better memory striding. - RC 6/24/16
            DO K = 0, NGRAIN1
                !
                DO J = 0, TVEC1
                    !
                    DO I = 0, TVEC1
                        !
                        WHERE (MY_PHASE .EQ. IPHASE)
                            !
                            TAN_STIF_VP(I, J, K, :) = TAN_STIF_VP(I, J, K, :) +&
                                & COMP(K,:) * &
                                & P_HAT_VEC(I+1, ISLIP+1) * &
                                & P_HAT_VEC(J+1, ISLIP+1)
                            !
                        ENDWHERE
                        !
                    ENDDO
                    !
                ENDDO
                !
            ENDDO !NGRAIN
            !
        ENDDO !N_SLIP
        !
        DEALLOCATE(P_HAT_VEC)
        DEALLOCATE(INDICES)
        !
    ENDDO !NUMPHASES
    !
    ! Transform tan_stif_vp to sample axes: [S]_sm = [Q] [S]_lat [Q]'.
    !
    CALL MAT_X_MAT5(QR5X5, TAN_STIF_VP, TEMP, NGRAIN, M_EL)
    !
    CALL MAT_X_MATT5(TEMP, QR5X5, TAN_STIF_VP, NGRAIN, M_EL)
    !
    ! Calculate secant moduli
    !
    DO IPHASE = 1, NUMPHASES
        !
        IF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .FALSE.) THEN
            !
            DO K = 0, NGRAIN1
                !
                DO J = 0, TVEC1
                    !
                    DO I = 0, TVEC1
                        !
                        WHERE (MY_PHASE .EQ. IPHASE)
                            !
                            STIF_VP(I, J, K, :) = &
                                & CRYSTAL_PARM(0,IPHASE)*TAN_STIF_VP(I, J, K, :)
                            !
                        ENDWHERE
                        !
                    ENDDO
                    !
                ENDDO
                !
            ENDDO
            !
        ELSEIF (CRYS_OPTIONS%USE_ANISO_M(IPHASE) .EQV. .TRUE.) THEN
            !
            ANISO_M_MIN = MINVAL(CRYS_OPTIONS%ANISO_M(IPHASE,:))
            !
            DO K = 0, NGRAIN1
                !
                DO J = 0, TVEC1
                    !
                    DO I = 0, TVEC1
                        !
                        WHERE (MY_PHASE .EQ. IPHASE)
                            !
                            STIF_VP(I, J, K, :) = &
                                & ANISO_M_MIN*TAN_STIF_VP(I, J, K, :)
                            !
                        ENDWHERE
                        !
                    ENDDO
                    !
                ENDDO
                !
            ENDDO
            !
        ENDIF
        !
    ENDDO !NUMPHASES 
    !
    ! Spread keinv to all crystals and transform it to sample axes.  
    ! deb - 6/11/2000
    !
    keinv_all = 0.0D0
    !
    ! In order to remove one dimension from my_phase (the dimension over
    ! `ngrain'), we add a loop over that dimension here. - hritz 9/15/05
    !
    DO IPHASE = 1, NUMPHASES
        !
        DO K = 0, NGRAIN1
            !
            DO I = 0, TVEC1
                !
                WHERE (MY_PHASE .EQ. IPHASE)
                    !
                    KEINV_ALL(I, I, K, :) = KEINV(I,IPHASE)
                    !
                ENDWHERE
                !
            ENDDO
            !
        ENDDO
        !
    ENDDO !NUMPHASES
    !
    CALL MAT_X_MAT5(QR5X5, KEINV_ALL, TEMP, NGRAIN, M_EL)
    !
    CALL MAT_X_MATT5(TEMP, QR5X5, KEINV_ALL, NGRAIN, M_EL)
    !
    ! Determinant of tensor V* .
    !
    V_TENSOR = E_ELAS
    E_ELAS_KK_GRN = SPREAD(E_ELAS_KK, DIM = 1, NCOPIES = NGRAIN)
    !
    DO I = 0, DIMS1
        !
        V_TENSOR(I, I, :, :) = V_TENSOR(I, I, :, :) +  E_ELAS_KK_GRN / 3. + 1.
        !
    ENDDO
    !
    CALL DETERMINANT_GRN(V_TENSOR, DETERM_V, NGRAIN, M_EL)
    !
    ! Compute elasto-visco-plastic crystal compliance.
    !
    DO J = 0, TVEC1
        !
        DO I = 0, TVEC1
            !
            STIF_EVP(I, J, :, :) = DETERM_V * STIF_VP(I, J, :, :) + &
                & DETERM_V * KEINV_ALL(I, J, :, :) * DTIMEI
            !
        ENDDO
        !
    ENDDO
    !
    IF (NR) THEN
        !
        DO I = 0, TVEC1
            !
            DO J = 0, TVEC1
                !
                TAN_STIF_EVP(I, J, :, :) = DETERM_V * TAN_STIF_VP(I, J, :, :) &
                    & + DETERM_V * KEINV_ALL(I, J, :, :) * DTIMEI
                !
            ENDDO
            !
        ENDDO
        !
    END IF
    !
    ! Compute elasto-visco-plastic crystal stiffness.
    !
    CALL INVERT5X5(STIF_EVP, NGRAIN, M_EL)
    !
    IF (NR) CALL INVERT5X5(TAN_STIF_EVP, NGRAIN, M_EL)
    !
    ! Force vector due to elastic terms
    !
    CALL FIND_WP_HAT(WP_HAT_VEC, E_ELAS, E_BAR, W_VEC_GRN, GDOT, QR5X5, &
        & DTIME, NGRAIN, M_EL)
    !
    CALL WP_HAT_MAT5X5_ALL(WP_HAT_VEC, WP_HAT_MAT, NGRAIN, M_EL)  
    !
    CALL MAT_X_VEC5(WP_HAT_MAT, E_VEC_SM, WP_X_E, NGRAIN, M_EL)
    !
    WP_X_E = WP_X_E - 1. / DTIME * E_BAR_SM
    !
    CALL MAT_X_VEC5(STIF_EVP, WP_X_E, F_ELAS, NGRAIN, M_EL)
    ! 
    ! Averaged values of stif_evp, f_elas
    !
    DO J = 0, TVEC1
        !
        F(J + 1, :) = SUM(F_ELAS(J, :, :) * WTS, DIM = 1)
        !
        DO I= 0, TVEC1
            !
            C(I + 1, J + 1, :) = SUM(STIF_EVP(I, J, :, :) * WTS,  DIM = 1)
            !
        ENDDO
        !
    ENDDO
    !
    IF (NR) THEN
        !
        DO J = 0, TVEC1
            !
            DO I = 0, TVEC1
                !
                C_TAN(I + 1, J + 1, :) = SUM(TAN_STIF_EVP(I, J, :, :) * WTS, &
                    & DIM = 1)
                !
            ENDDO
            !
        ENDDO
        !
    END IF
    !
    DETV = SUM(DETERM_V * WTS, DIM = 1)
    !
    ! Fix [c] & {f} to be consistent with FE equations.
    !
    DO J = 1, TVEC
        !
        F(J, :) = ALPHA(J) * F(J, :)
        !
        DO I = 1, TVEC
            !
            C(I, J, :) = ALPHA(I) * C(I, J, :) * ALPHA(J)
            !
        ENDDO
        !
    ENDDO
    !
    IF (NR) THEN
        !
        DO J = 1, TVEC
            !
            DO I = 1, TVEC
                !
                C_TAN(I, J, :) = ALPHA(I) * C_TAN(I, J, :) * ALPHA(J)
                !
            ENDDO
            !
        ENDDO
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE ANISO_EVPS
    !
END MODULE ANISO_EVPS_MOD
