! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE MATERIAL_MATRIX_VP_MOD
!
! Material matrix for the viscoplastic solution
!
! Contains subroutines:
! MATERIAL_MATRIX_VP: Material matrix for the viscoplastic solution
! ISOTROPIC: Return viscosity for the isotropic viscoplastic solution
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY: RK=>REAL_KIND
!
! From libfepx:
!
USE DIMENSIONS_MOD
USE KINEMATICS_MOD
USE READ_INPUT_MOD
USE STRESS_SOLVE_VP_MOD
!
! From libparallel:
!
USE PARALLEL_MOD
!
IMPLICIT  NONE
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: MATERIAL_MATRIX_VP
PUBLIC :: DEFRATE
!
CONTAINS
    !
    SUBROUTINE MATERIAL_MATRIX_VP(TYPE, STIF, DNDX, DNDY, DNDZ, GVEL, SCALE, &
        & QR5X5, WTS, EPSEFF, DTIME, INCR)
    !
    ! Material matrix for the viscoplastic solution
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! TYPE:
    ! STIF:
    ! DNDX:
    ! DNDY:
    ! DNDZ:
    ! GVEL:
    ! SCALE:
    ! QR5X5:
    ! WTS:
    ! EPSEFF:
    ! DTIME:
    ! INCR:
    !
    INTEGER :: TYPE
    REAL(RK) :: STIF(TVEC, TVEC, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDX(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDY(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: DNDZ(0:NNPE, EL_SUB1:EL_SUP1)
    REAL(RK) :: GVEL(0:KDIM1, EL_SUB1:EL_SUP1)
    REAL(RK) :: SCALE(EL_SUB1:EL_SUP1)
    REAL(RK) :: QR5X5(0:TVEC1, 0:TVEC1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: WTS(0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK) :: EPSEFF(EL_SUB1:EL_SUP1)
    REAL(RK) :: DTIME
    INTEGER :: INCR
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: M_EL
    REAL(RK) :: D(0:DIMS1, 0:DIMS1, EL_SUB1:EL_SUP1)
    REAL(RK) :: D_VEC(0:TVEC1, EL_SUB1:EL_SUP1)
    REAL(RK) :: CMU(EL_SUB1:EL_SUP1)
    REAL(RK) :: TEMP_K(EL_SUB1:EL_SUP1)
    REAL(RK) :: DSDT_ISO(EL_SUB1:EL_SUP1)
    REAL(RK) :: STATE_ISO(EL_SUB1:EL_SUP1)
    REAL(RK) :: CCONST(TVEC)
    !
    ! Data:
    !
    DATA CCONST /1.0d0, 3.0d0, 4.0d0, 4.0d0, 4.0d0/
    !
    !---------------------------------------------------------------------------
    !
    M_EL = EL_SUP1 - EL_SUB1 + 1
    !
    CALL DEFRATE(D, DNDX, DNDY, DNDZ, GVEL)
    !
    ! tm268_M13:
    !CALL eff_def(EPSEFF, eqplas, d, DTIME, M_EL)
    CALL EFF_DEF(EPSEFF, D, DTIME, M_EL)
    !
    ! hr-tm
    ! we aren't going through this portion to update for two-phase compatibility
    IF (TYPE .EQ. ANISOTROPIC_VP) THEN
        !
        CALL par_quit('Error  :     > ANISOTROPIC_VP is no longer implemented.')
        !
    ELSE IF (TYPE .EQ. ISOTROPIC_VP) THEN
        !
        ! dbg: Note the hardwired temperature and state variable.
        !
        TEMP_K = 273.0D0 + 400.0D0
        STATE_ISO = 20.0D6
        !
        CALL ISOTROPIC(EPSEFF, TEMP_K, STATE_ISO, CMU, DSDT_ISO)
        !
        STIF = 0.0D0
        STIF(1, 1, :) = CCONST(1) * CMU
        STIF(2, 2, :) = CCONST(2) * CMU
        STIF(3, 3, :) = CCONST(3) * CMU
        STIF(4, 4, :) = CCONST(4) * CMU
        STIF(5, 5, :) = CCONST(5) * CMU
        !
    END IF
    !
    SCALE = 0.0D0
    !
    DO I = 1, TVEC
        !
        SCALE = SCALE + STIF(I, I, :) / CCONST(I)
        !
    END DO
    !
    SCALE = SCALE / 5.0D0
    !
    END SUBROUTINE MATERIAL_MATRIX_VP
    !
    !===========================================================================
    !
    SUBROUTINE ISOTROPIC(DII, T, S_PA, VISC, DSDT)
    !
    ! Return viscosity for isotropic problem. Apparently, this routine was
    !   disabled and was replaced by a linear viscous problem with unit
    !   viscosity. MPK: Is this still true?
    !
    ! Computes the viscosity for an element based on the effective strain rate.
    !   Uses the fit of Dawson to Hart's model for pure Aluminum.
    !
    ! Modifications:
    ! The effective viscosity computation has been changed to:
    !   visc = (2nd invariant of sigma) / 3*(2nd invariant of strain rate)
    !   so as to be compatible with definition of strain rate invariant
    !   by ajb (9/25/87). Units are MPa, kJoule/mole-K, etc.
    !
    !
    !----------------------------------------------------------------------
    !
    ! Arguments:
    ! DII: The effective strain rate for the element
    ! T: Element temperature (F)
    ! S_PA: The value of the state variable (Pascals)
    ! VISC: The resulting viscosity
    ! DSDT: The derivative of the state variable
    !
    REAL(RK) :: DII(EL_SUB1:EL_SUP1)
    REAL(RK) :: T(EL_SUB1:EL_SUP1)
    REAL(RK) :: S_PA(EL_SUB1:EL_SUP1)
    REAL(RK) :: VISC(EL_SUB1:EL_SUP1)
    REAL(RK) :: DSDT(EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    ! Constants:
    !
    !REAL(RK) :: a0
    !PARAMETER ( a0     = 9.64d52     )
    REAL(RK), PARAMETER :: LOG_A0 = 122.0D0
    REAL(RK), PARAMETER :: C0 = 6.19D-9
    REAL(RK), PARAMETER :: QPR = 1.45D4
    REAL(RK), PARAMETER :: GE = 24.20D3
    REAL(RK), PARAMETER :: M = 7.80D0
    REAL(RK), PARAMETER :: F0 = 2.12D19
    REAL(RK), PARAMETER :: QR = QPR
    REAL(RK), PARAMETER :: SMALLM = 5.0D0
    REAL(RK), PARAMETER :: LAMBDA = 0.14D0
    REAL(RK), PARAMETER :: MP = 3.5D0
    REAL(RK), PARAMETER :: N = 6.0D0
    !
    REAL(RK) :: A(EL_SUB1:EL_SUP1)
    REAL(RK) :: S(EL_SUB1:EL_SUP1)
    REAL(RK) :: DSTAR(EL_SUB1:EL_SUP1)
    REAL(RK) :: SIGMA(EL_SUB1:EL_SUP1)
    REAL(RK) :: SIGMAP(EL_SUB1:EL_SUP1)
    REAL(RK) :: SIGMAV(EL_SUB1:EL_SUP1)
    !
    INTEGER :: I
    !
    !---------------------------------------------------------------------------
    !
    ! Place the state variable in Pascals
    !
    VISC = 1.0D0
    !
    RETURN
    !
    S = S_PA / 1.0D6
    !
    ! Compute flow stress in viscous (friction) element
    !   a = a0 * exp(-QPR/t)
    A = LOG_A0 + (-QPR / T)
    A = EXP(A)
    !
    SIGMAV = GE * (DII / A) ** (1.0 / M)
    !
    ! Compute flow stress in the plastic element
    !
    DSTAR  = F0 * ((S / GE) ** SMALLM) * EXP(-QR / T)
    SIGMAP = S * EXP(-(DSTAR/DII) ** LAMBDA )
    !
    ! Total flow stress
    !
    SIGMA = SIGMAV + SIGMAP
    !
    ! Compute the derivative of the state vaiable
    !
    DSDT = C0 * S * DII * ((GE / S) ** MP) * ((SIGMAP/S) ** N)
    !
    ! Convert to Pascals
    !
    !SIGMA = s * DII**0.10
    SIGMA = SIGMA * 1.0D6
    !
    ! Compute viscosity
    !
    VISC = SIGMA / (3.0D0 * DII)
    !
    print *,'SIGMA', MINVAL(SIGMA), MAXVAL(SIGMA)
    print *,'visc', MINVAL(VISC), MAXVAL(VISC)
    !
    RETURN
    !
    END SUBROUTINE ISOTROPIC
    !
END MODULE MATERIAL_MATRIX_VP_MOD
