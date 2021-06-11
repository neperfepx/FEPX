! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE FIBER_AVERAGE_MOD
!
! Performs "powder diffraction" analysis on the data at the end of each step in
! the load history. Essentially calculates fiber-averaged values of certain
! metrics. A "valid fiber element" is an element that belongs to a defined
! crystallographic fiber and meets the diffraction condition.
!
! Contains subroutines:
! RUN_FIBER_AVERAGE: Run fiber average routine
! INITIALIZE_FIBER_AVERAGE: Initialize fiber average routine
! PRINT_ELEMS: Print valid fiber elements for each processor
! PRINT_STATS: Print fiber averaging statistics
! INITIALIZE_CRYSTAL_SYMS: Initialize crystal symmetry operations
!
! From libf95:
!
USE LIBF95
!
! From libfepx:
!
USE DIMENSIONS_MOD
USE MATRIX_OPERATIONS_MOD
USE MICROSTRUCTURE_MOD
USE READ_INPUT_MOD
USE UNITS_MOD
!
! From libparallel:
!
USE GATHER_SCATTER_MOD
USE PARALLEL_MOD
!
IMPLICIT NONE
!
! Private data
!
! HKLFIB: 3 x NFIB array of crystal directions
! UVWFIB: 3 x NFIB array of sample directions
! PHASEFIB: NFIB array of phase for each fiber
! TOLFIB: NFIB array of angular tolerance for each fiber
! ELVOL: Element volumes
! ISINTERIOR: Logical array of whether each element is within mesh interior
!
INTEGER,  PRIVATE  :: NFIB, NFIB1
INTEGER,  ALLOCATABLE, PRIVATE  :: HKLFIB(:, :)
REAL(RK), ALLOCATABLE, PRIVATE  :: UVWFIB(:, :)
INTEGER,  ALLOCATABLE, PRIVATE  :: PHASEFIB(:)
REAL(RK), ALLOCATABLE, PRIVATE  :: TOLFIB(:)
!REAL(RK), ALLOCATABLE, PRIVATE  :: ELVOL(:)
LOGICAL,  ALLOCATABLE, PRIVATE  :: ISINTERIOR(:)
!
INTEGER,  PARAMETER, PRIVATE :: NSYMCUB = 24, NSYMCUB1 = 23
INTEGER,  PARAMETER, PRIVATE :: NSYMHEX = 12, NSYMHEX1 = 11
REAL(RK), PRIVATE :: RMATSYMCUB(0:DIMS1, 0:DIMS1, 0:NSYMCUB1)
REAL(RK), PRIVATE :: RMATSYMHEX(0:DIMS1, 0:DIMS1, 0:NSYMHEX1)
!
REAL(RK), PARAMETER, PRIVATE :: PI_OVER_180 = (4.0D0 * DATAN(1.0D0)) / 180.0D0
!
CONTAINS
    !
    SUBROUTINE RUN_FIBER_AVERAGE(ISTEP, C_ANGS, ELAS_TOT6, DPEFF, CRSS, ELVOL)
    !
    ! Run fiber average processing
    !
    !---------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! Arguments:
    ! ISTEP: Index of current step number
    ! C_ANGS: Current elemental orientations
    ! ELAS_TOT6: Current elemental strains (6-vector)
    ! DPEFF: Current effective plastic deformation rate
    ! CRSS: Current slip systen strength
    ! ELVOL: Element volumes
    !
    INTEGER,  INTENT(IN) :: ISTEP
    REAL(RK), INTENT(IN) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: ELAS_TOT6(0:TVEC , 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: DPEFF(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: ELVOL(EL_SUB1:EL_SUP1)
    !
    ! Locals:
    ! NSYMCUB/NSYMHEX: Number of symmetries for cubic and hexagonal crystals
    ! NSYMCUB1/NSYMHEX1: 0-indexed values of the above local variables
    ! NUM_VFE: NFIB1 array with number of "valid fiber elements" per-fiber
    ! VOL_VFE: NFIB1 array with total "valid fiber element" volumes per-fiber
    ! LS_AVG/LS_STD: Avg. and std. deviation of elastic lattice strain
    ! DPEFF_AVG/DPEFF_STD: Avg. and std. deviation of eff. plastic def. rate.
    ! SGD_AVG/SGD_STD: Avg. and std. deviation of gammadot summation
    ! CRSS_AVG/CRSS_STD: Avg. and std. deviation of CRSS per slip system
    ! P_NUM_VFE_INT: Number of "valid fiber elements" on local partition
    ! P_NUM_VFE: Local partition array of "valid fiber elements" 
    ! NUM_VFE_IFIB: Summation of "valid fiber elements" per-fiber
    ! P_VOL_VFE: Local partition array of "valid fiber elements'" volumes 
    ! VOL_VFD_FIB: Summation of "valid fiber elements'" volume per-fiber
    ! P_LS_SUM/LS_SUM: Local/global volume-weighted lattice strain summations
    ! P_DPEFF_SUM/DPEFF_SUM: Local/global volume-weighted eff. plastic def. 
    !   rate. summations
    ! P_SGD_SUM/SGD_SUM: Local/global volume-weighted gammadot summations
    ! P_CRSS_SUM/CRSS_SUM: Local/global volume-weighted CRSS per-slip system
    ! LS: Elemental lattice strains
    ! SGD: Elemental sum of gammadots
    ! I/J/K/KELEM/ISS/IFB/IPLN: Generic looping indices
    ! NSYM1: Temporary storage variable for NSYMCUB or NSYMHEX for looping
    ! NUNQPLN: Counter of unique plane normals
    ! ISUNIQUEPLN: Logical denoting a unique plane normal
    ! ISVFE: Logical mask of "valid fiber elements" on local core
    ! ALLPLN: All equivalent plane normals from symmetry
    ! RMATSYM: Storage array for symmetry rotation matrices RMATSYMCUB/HEX
    ! HKL: DIMS1 array of extracted HKL values for IFIB
    ! UVW: DIMS1 array of extracted UVW values for IFIB
    ! PLN: Lattice plane unit normal vectors in Cartesian crystal basis
    ! DOTPROD: Dot product result storage variable
    ! SVEC: Plane normal in the sample basis
    ! CVEC: Sample direction in the crystal basis
    ! EPS: Elemental elastic strains in sample basis
    !
    INTEGER, PARAMETER  :: NSYMCUB = 24, NSYMCUB1 = 23
    INTEGER, PARAMETER  :: NSYMHEX = 12, NSYMHEX1 = 11
    REAL(RK), PARAMETER :: TOL = 1.0D-4
    INTEGER  :: NUM_VFE(0:NFIB1) ! VFE = "Valid Fiber Element"
    REAL(RK) :: VOL_VFE(0:NFIB1)
    REAL(RK) :: LS_AVG(0:NFIB1), LS_STD(0:NFIB1)
    REAL(RK) :: DPEFF_AVG(0:NFIB1), DPEFF_STD(0:NFIB1)
    REAL(RK) :: SGD_AVG(0:NFIB1), SGD_STD(0:NFIB1)
    REAL(RK) :: CRSS_AVG(0:NFIB1,0:MAXSLIP1), CRSS_STD(0:NFIB1,0:MAXSLIP1)
    INTEGER  :: P_NUM_VFE_INT
    REAL(RK) :: P_NUM_VFE, NUM_VFE_IFIB ! NUM_VFE is real for use with par_sum
    REAL(RK) :: P_VOL_VFE, VOL_VFE_IFIB
    REAL(RK) :: P_LS_SUM, LS_SUM
    REAL(RK) :: P_DPEFF_SUM, DPEFF_SUM
    REAL(RK) :: P_SGD_SUM, SGD_SUM
    REAL(RK) :: P_CRSS_SUM(0:MAXSLIP1), CRSS_SUM(0:MAXSLIP1)
    REAL(RK) :: LS(EL_SUB1:EL_SUP1)
    REAL(RK) :: SGD(EL_SUB1:EL_SUP1)
    INTEGER :: I, J, K, KELEM, ISS, IFIB, IPLN
    INTEGER :: NSYM1
    INTEGER :: NUNQPLN
    LOGICAL :: ISUNIQUEPLN
    LOGICAL :: ISVFE(EL_SUB1:EL_SUP1)
    REAL(RK), ALLOCATABLE :: ALLPLN(:, :)
    REAL(RK), ALLOCATABLE :: RMATSYM(:, :, :)
    REAL(RK) :: HKL(0:DIMS1), UVW(0:DIMS1), PLN(0:DIMS1)
    REAL(RK) :: DOTPROD
    REAL(RK) :: SVEC(0:DIMS1), CVEC(0:DIMS1)
    REAL(RK) :: EPS(0:DIMS1, 0:DIMS1)
    !
    !---------------------------------------------------------------------------
    !
    ! Initialization
    !
    LS_AVG = 0.0D0
    DPEFF_AVG = 0.0D0
    SGD_AVG = 0.0D0
    CRSS_AVG = 0.0D0
    LS_STD = 0.0D0
    DPEFF_STD = 0.0D0
    SGD_STD = 0.0D0
    CRSS_STD = 0.0D0
    !
    ! Calculate sum of the gammadots
    !
    SGD = 0.0D0
    DO KELEM = EL_SUB1, EL_SUP1
        !
        DO ISS = 0, MAXSLIP1
            !
            SGD(KELEM) = SGD(KELEM) + ABS(GAMMADOT(ISS, 0, KELEM))
            !
        ENDDO
        !
    ENDDO
    !
    ! Loop over all fibers
    !
    DO IFIB = 0, NFIB1
        !
        ! Initialization
        !
        LS = 0.0D0
        P_NUM_VFE = 0.0D0
        P_VOL_VFE = 0.0D0
        P_LS_SUM = 0.0D0
        P_DPEFF_SUM = 0.0D0
        P_SGD_SUM = 0.0D0
        P_CRSS_SUM = 0.0D0
        NUM_VFE_IFIB = 0.0D0
        VOL_VFE_IFIB = 0.0D0
        LS_SUM = 0.0D0
        DPEFF_SUM = 0.0D0
        SGD_SUM = 0.0D0
        CRSS_SUM = 0.0D0
        !
        ! Extract HKL
        !
        HKL(0) = REAL(HKLFIB(0, IFIB), RK)
        HKL(1) = REAL(HKLFIB(1, IFIB), RK)
        HKL(2) = REAL(HKLFIB(2, IFIB), RK)
        !
        ! Extract sample direction
        !
        UVW(0) = REAL(UVWFIB(0, IFIB), RK)
        UVW(1) = REAL(UVWFIB(1, IFIB), RK)
        UVW(2) = REAL(UVWFIB(2, IFIB), RK)
        UVW = UVW / DSQRT(UVW(0)*UVW(0) + UVW(1)*UVW(1) + UVW(2)*UVW(2))
        !
        ! Find multiplicity
        !
        SELECT CASE (CTYPE(PHASEFIB(IFIB))%CLASS)
            !
            CASE (1, 2)
                !
                ! Cubic:
                ! construct lattice plane unit normal vectors (PLN) in Cartesian
                ! crystal basis
                !
                PLN = HKL/DSQRT(HKL(0)*HKL(0) + HKL(1)*HKL(1) + HKL(2)*HKL(2))
                NSYM1 = NSYMCUB1
                ALLOCATE(RMATSYM(0:DIMS1, 0:DIMS1, 0:NSYM1))
                RMATSYM = RMATSYMCUB
                !
            CASE (3)
                !
                ! Hexagonal:
                ! construct lattice plane unit normal vectors (PLN) in Cartesian
                ! crystal basis
                !
                PLN(0) = HKL(0)
                PLN(1) = (2*HKL(1)+HKL(0)) / DSQRT(3.0D0)
                PLN(2) = HKL(2)/CRYSTAL_PARM(11,PHASEFIB(IFIB))
                PLN = PLN / DSQRT(PLN(0)*PLN(0) + PLN(1)*PLN(1) + PLN(2)*PLN(2))
                NSYM1 = NSYMHEX1
                ALLOCATE(RMATSYM(0:DIMS1, 0:DIMS1, 0:NSYM1))
                RMATSYM = RMATSYMHEX
                !
            CASE DEFAULT
                !
                CALL PAR_QUIT('Error  :     > RUN_FIBER_AVERAGE: &
                    &Invalid crystal type')
            !
        END SELECT
        !
        ALLOCATE(ALLPLN(0:DIMS1, 0:NSYM1))
        ALLPLN = 0.0
        !
        ! Find equivalent plane unit normals
        !
        DO K = 0, NSYM1
            !
            DO I = 0, DIMS1
                !
                DO J = 0, DIMS1
                    !
                    ALLPLN(I, K) = ALLPLN(I, K) + RMATSYM(I, J, K) * PLN(J)
                    !
                ENDDO
                !
            ENDDO
            !
        ENDDO
        !
        ! Find unique plane unit normals
        !
        NUNQPLN = 1
        !
        DO I = 1, NSYM1 ! Start from 1, since 0-th entry is unique
            !
            ISUNIQUEPLN = .TRUE.
            !
            DO J = 0, (NUNQPLN - 1)
                !
                DOTPROD = ALLPLN(0, I) * ALLPLN(0, J) + &
                    & ALLPLN(1, I) * ALLPLN(1, J) + &
                    & ALLPLN(2, I) * ALLPLN(2, J)
                !
                IF (ABS(1.0D0 - ABS(DOTPROD)) .LT. TOL) ISUNIQUEPLN = .FALSE.
                !
            ENDDO
            !
            IF (ISUNIQUEPLN) THEN
                !
                NUNQPLN = NUNQPLN + 1
                ALLPLN(:, NUNQPLN - 1) = ALLPLN(:, I)
                !
            ENDIF
            !
        ENDDO
        !
        ! Calculate misorientation
        !
        ISVFE = .FALSE.
        DO KELEM = EL_SUB1, EL_SUP1 ! Loop over all elements
            !
            IF (ISINTERIOR(KELEM) .AND. (PHASE(KELEM) .EQ. PHASEFIB(IFIB))) THEN
                !
                DOTPROD = 0.0D0
                !
                DO IPLN = 0, (NUNQPLN - 1)
                    !
                    SVEC = 0.0D0
                    !
                    DO I = 0, DIMS1
                        !
                        DO J = 0, DIMS1
                            !
                            SVEC(I) = SVEC(I) + C_ANGS(I, J, 0, KELEM) * &
                                & ALLPLN(J, IPLN)
                            !
                        ENDDO
                        !
                    ENDDO
                    !
                    DOTPROD = MAX(DOTPROD, ABS(SVEC(0) * UVW(0) + &
                        & SVEC(1) * UVW(1) + SVEC(2) * UVW(2)))
                    !
                ENDDO
                !
                IF (DOTPROD .GT. TOLFIB(IFIB)) THEN
                    !
                    ! Crystal belongs to fiber
                    !
                    ISVFE(KELEM) = .TRUE.
                    !
                    ! Calculate lattice strain
                    !
                    ! Scattering vector in crystal coordinates
                    ! {n}_c = [R]^T {n}_s
                    !
                    CVEC = 0.0D0
                    !
                    DO I = 0, DIMS1
                        !
                        DO J = 0, DIMS1
                            !
                            CVEC(I) = CVEC(I) + C_ANGS(J, I, 0, KELEM) * UVW(J)
                            !
                        ENDDO
                        !
                    ENDDO
                    !
                    ! Strain in crystal coordinates
                    !
                    EPS(0, 0) = ELAS_TOT6(0, 0, KELEM)
                    EPS(0, 1) = ELAS_TOT6(1, 0, KELEM)
                    EPS(0, 2) = ELAS_TOT6(2, 0, KELEM)
                    EPS(1, 0) = ELAS_TOT6(1, 0, KELEM)
                    EPS(1, 1) = ELAS_TOT6(3, 0, KELEM)
                    EPS(1, 2) = ELAS_TOT6(4, 0, KELEM)
                    EPS(2, 0) = ELAS_TOT6(2, 0, KELEM)
                    EPS(2, 1) = ELAS_TOT6(4, 0, KELEM)
                    EPS(2, 2) = ELAS_TOT6(5, 0, KELEM)
                    !
                    DO I = 0, DIMS1
                        !
                        DO J = 0, DIMS1
                            !
                            LS(KELEM) = LS(KELEM) + CVEC(I)*EPS(I, J) * CVEC(J)
                            !
                        ENDDO
                        !
                    ENDDO
                    !
                    ! Sum some quantities
                    !
                    P_NUM_VFE = P_NUM_VFE + 1.0D0
                    P_VOL_VFE = P_VOL_VFE + ELVOL(KELEM)
                    !
                    ! The following sums are weighted by element volume
                    !
                    P_LS_SUM = P_LS_SUM + LS(KELEM)*ELVOL(KELEM)
                    P_DPEFF_SUM = P_DPEFF_SUM + DPEFF(KELEM)*ELVOL(KELEM)
                    P_SGD_SUM = P_SGD_SUM + SGD(KELEM)*ELVOL(KELEM)
                    !
                    ! Compute the elt-volume-weighted CRSS for each slip system
                    !
                    ! Compute as normal if FCC or BCC
                    IF ((CTYPE(PHASEFIB(IFIB))%CLASS .EQ. 1) .OR. &
                        & (CTYPE(PHASEFIB(IFIB))%CLASS .EQ. 2)) THEN
                        !
                        DO ISS = 0, MAXSLIP1
                            !
                            P_CRSS_SUM(ISS) = P_CRSS_SUM(ISS) &
                                & + CRSS(ISS, 0, KELEM) * ELVOL(KELEM)
                            !
                        ENDDO
                        !
                    ! Else, scale CRSS values by strength ratios for HCP
                    ELSE IF (CTYPE(PHASEFIB(IFIB))%CLASS .EQ. 3) THEN
                        !
                        DO ISS = 0, 2
                            !
                            P_CRSS_SUM(ISS) = P_CRSS_SUM(ISS) &
                                & + CRSS(ISS, 0, KELEM) * ELVOL(KELEM)
                            !
                        ENDDO
                        !
                        DO ISS = 3, 5
                            !
                            P_CRSS_SUM(ISS) = P_CRSS_SUM(ISS) &
                                & + CRYS_OPTIONS%PRISMATIC_TO_BASAL&
                                &(PHASEFIB(IFIB)) * CRSS(ISS, 0,  KELEM) &
                                & * ELVOL(KELEM)
                            !
                        ENDDO
                        !
                        DO ISS = 6, 17
                            !
                            P_CRSS_SUM(ISS) = P_CRSS_SUM(ISS) &
                                & + CRYS_OPTIONS%PYRAMIDAL_TO_BASAL&
                                &(PHASEFIB(IFIB)) * CRSS(ISS, 0,  KELEM) &
                                & * ELVOL(KELEM)
                            !
                        ENDDO
                        !
                    ENDIF
                    !
                ENDIF
                !
            ENDIF
            !
        ENDDO ! Loop over all elements
        !
        ! Calculate fiber-averaged statistics
        !
        CALL PAR_SUM(P_NUM_VFE, NUM_VFE_IFIB)
        CALL PAR_SUM(P_VOL_VFE, VOL_VFE_IFIB)
        CALL PAR_SUM(P_LS_SUM, LS_SUM)
        CALL PAR_SUM(P_DPEFF_SUM, DPEFF_SUM)
        CALL PAR_SUM(P_SGD_SUM, SGD_SUM)
        !
        ! Loop over the slip systems to recover individual CRSS values
        !
        DO ISS = 0, MAXSLIP1
            !
            CALL PAR_SUM(P_CRSS_SUM(ISS), CRSS_SUM(ISS))
            !
        ENDDO
        !
        NUM_VFE(IFIB) = NINT(NUM_VFE_IFIB)
        VOL_VFE(IFIB) = VOL_VFE_IFIB
        !
        IF (VOL_VFE_IFIB .GT. 0.0) THEN
            !
            LS_AVG(IFIB) = LS_SUM / VOL_VFE_IFIB
            DPEFF_AVG(IFIB) = DPEFF_SUM / VOL_VFE_IFIB
            SGD_AVG(IFIB) = SGD_SUM / VOL_VFE_IFIB
            !
            DO ISS = 0, MAXSLIP1
                !
                CRSS_AVG(IFIB,ISS) = CRSS_SUM(ISS) / VOL_VFE_IFIB
                !
            ENDDO
            !
        ENDIF
        !
        ! Calculate standard deviations
        !
        P_LS_SUM = 0.0D0
        P_DPEFF_SUM = 0.0D0
        P_SGD_SUM = 0.0D0
        P_CRSS_SUM = 0.0D0
        LS_SUM = 0.0D0
        DPEFF_SUM = 0.0D0
        SGD_SUM = 0.0D0
        CRSS_SUM = 0.0D0
        !
        IF (VOL_VFE_IFIB .GT. 0.0D0) THEN
            !
            DO KELEM = EL_SUB1, EL_SUP1
                !
                IF (ISVFE(KELEM)) THEN
                    !
                    P_LS_SUM = P_LS_SUM + (LS(KELEM) - LS_AVG(IFIB))**2 * &
                        & ELVOL(KELEM)
                    P_DPEFF_SUM = P_DPEFF_SUM + (DPEFF(KELEM) - &
                        & DPEFF_AVG(IFIB))**2 * ELVOL(KELEM)
                    P_SGD_SUM = P_SGD_SUM + (SGD(KELEM) - SGD_AVG(IFIB))**2 * &
                        & ELVOL(KELEM)
                    !
                    ! Compute as normal if FCC or BCC
                    IF ((CTYPE(PHASEFIB(IFIB))%CLASS .EQ. 1) .OR. &
                        & (CTYPE(PHASEFIB(IFIB))%CLASS .EQ. 2)) THEN
                        !
                        DO ISS = 0, MAXSLIP1
                            !
                            P_CRSS_SUM(ISS) = P_CRSS_SUM(ISS) &
                                & + (CRSS(ISS, 0,  KELEM) &
                                & - CRSS_AVG(IFIB,ISS))**2 * ELVOL(KELEM)
                            !
                        ENDDO
                        !
                    ! Else, scale CRSS values by strength ratios for HCP
                    ELSE IF (CTYPE(PHASEFIB(IFIB))%CLASS .EQ. 3) THEN
                        !
                        DO ISS = 0, 2
                            !
                            P_CRSS_SUM(ISS) = P_CRSS_SUM(ISS) &
                                & + (CRSS(ISS, 0,  KELEM) &
                                & - CRSS_AVG(IFIB,ISS))**2 * ELVOL(KELEM)
                            !
                        ENDDO
                        !
                        DO ISS = 3, 5
                            !
                            P_CRSS_SUM(ISS) = P_CRSS_SUM(ISS) &
                                & + CRYS_OPTIONS%PRISMATIC_TO_BASAL&
                                &(PHASEFIB(IFIB)) * (CRSS(ISS, 0,  KELEM) &
                                & - CRSS_AVG(IFIB,ISS))**2 * ELVOL(KELEM)
                            !
                        ENDDO
                        !
                        DO ISS = 6, 17
                            !
                            P_CRSS_SUM(ISS) = P_CRSS_SUM(ISS) &
                                & + CRYS_OPTIONS%PYRAMIDAL_TO_BASAL&
                                &(PHASEFIB(IFIB)) * (CRSS(ISS, 0,  KELEM) &
                                & - CRSS_AVG(IFIB,ISS))**2 * ELVOL(KELEM)
                            !
                        ENDDO
                        !
                    ENDIF
                    !
                ENDIF
                !
            ENDDO
            !
            CALL PAR_SUM(P_LS_SUM, LS_SUM)
            CALL PAR_SUM(P_DPEFF_SUM, DPEFF_SUM)
            CALL PAR_SUM(P_SGD_SUM, SGD_SUM)
            !
            DO ISS = 0, MAXSLIP1
                !
                CALL PAR_SUM(P_CRSS_SUM(ISS), CRSS_SUM(ISS))
                !
            ENDDO
            !
            LS_STD(IFIB) = DSQRT(LS_SUM/VOL_VFE_IFIB)
            DPEFF_STD(IFIB) = DSQRT(DPEFF_SUM/VOL_VFE_IFIB)
            SGD_STD(IFIB) = DSQRT(SGD_SUM/VOL_VFE_IFIB)
            !
            DO ISS = 0, MAXSLIP1
                !
                CRSS_STD(IFIB,ISS) = DSQRT(CRSS_SUM(ISS)/VOL_VFE_IFIB)
                !
            ENDDO
            !
        ENDIF
        !
        DEALLOCATE(ALLPLN)
        DEALLOCATE(RMATSYM)
        !
        ! Print list of valid fiber elements
        !
        P_NUM_VFE_INT = NINT(P_NUM_VFE)
        !
        CALL PRINT_ELEMS(ISTEP, IFIB, P_NUM_VFE_INT, ISVFE)
        !
    ENDDO ! Loop over all fibers
    !
    ! Print statsitics
    !
    CALL PRINT_STATS(LS_AVG, LS_STD, DPEFF_AVG, DPEFF_STD, SGD_AVG, SGD_STD, &
        & CRSS_AVG, CRSS_STD, NUM_VFE, VOL_VFE, ISTEP)
    !
    RETURN
    !
    END SUBROUTINE RUN_FIBER_AVERAGE
    !
    !===========================================================================
    !
    SUBROUTINE INITIALIZE_FIBER_AVERAGE(INPUT_UNIT)
    !
    ! Read and process input data in fiber average file, open files for
    ! output, initialize symmetry operations
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! INPUT_UNIT: Unit for input file
    !
    INTEGER :: INPUT_UNIT
    !
    ! Locals:
    ! IFIB/IELEM: Generic looping indices
    ! IVOL: Element volume read in from file `simulation.vol'
    ! IINT: Interior element logical read in from file `simulation.int'
    ! FILENAME: Output file name string for per-core data
    ! IOFILE: Input file names for data read in
    ! CHARID: Local processor ID (1-indexed) appended to FILENAME
    ! IOSTATUS: Returns value about file command success
    !
    INTEGER :: IFIB
    INTEGER :: IELEM
    INTEGER :: IERR
    INTEGER :: IINT
    CHARACTER*128 FILENAME, IOFILE
    CHARACTER*4 CHARID      ! Assumes less than 10,000 processes
    INTEGER  :: IOSTATUS
    !
    !---------------------------------------------------------------------------
    !
    ! Initialize crystal symmetries
    !
    CALL INITIALIZE_CRYSTAL_SYMS()
    !
    ! Read input data
    !
    IOFILE = 'simulation.fib'
    !
    OPEN(INPUT_UNIT, FILE = IOFILE, STATUS = 'OLD', ACTION = 'READ', &
        & IOSTAT = IOSTATUS)
    !
    IF (IOSTATUS .NE. 0) THEN
        !
        CALL PAR_QUIT("Error  :     > Failure to open `simulation.fib' file.")
        !
    END IF
    !
    ! Determine NFIB from the number of lines in the file
    !
    NFIB = 0
    !
    DO ! Count total number of lines in file
        !
        READ(INPUT_UNIT, *, IOSTAT = IERR)
        !
        IF (IERR .NE. 0) EXIT
        NFIB = NFIB + 1
        !
    END DO
    !
    ! Rewind the file to begin read-in of the fibers
    !
    REWIND (INPUT_UNIT)
    !
    NFIB1 = NFIB-1
    !
    ALLOCATE(HKLFIB(0:DIMS1, 0:NFIB1))
    ALLOCATE(UVWFIB(0:DIMS1, 0:NFIB1))
    ALLOCATE(PHASEFIB(0:NFIB1))
    ALLOCATE(TOLFIB(0:NFIB1))
    !ALLOCATE(ELVOL(EL_SUB1:EL_SUP1))
    ALLOCATE(ISINTERIOR(EL_SUB1:EL_SUP1))
    !
    ! Read-in the fibers per-line
    !
    DO IFIB = 0, NFIB1
        !
        ! Read in the entire line for the fiber
        !
        READ(INPUT_UNIT, *) HKLFIB(0,IFIB), HKLFIB(1,IFIB), HKLFIB(2,IFIB), &
            & UVWFIB(0,IFIB), UVWFIB(1,IFIB), UVWFIB(2,IFIB), PHASEFIB(IFIB), &
            & TOLFIB(IFIB)
        !
    END DO
    !
    TOLFIB = COS(TOLFIB * PI_OVER_180)
    !
    CLOSE(UNIT = INPUT_UNIT)
    !
    ! Read interior elements
    !
    IF (FIBER_AVERAGE_OPTIONS%READ_INTERIOR) THEN
        !
        ! Read interior file
        !
        IOFILE = 'simulation.int'
        !
        OPEN(INPUT_UNIT, FILE = IOFILE, STATUS = 'OLD', ACTION = 'READ')
        !
        ISINTERIOR = .FALSE.
        !
        DO IELEM = 0, MAXEL1
            !
            READ(INPUT_UNIT, *) IINT
            !
            IF ((IELEM .GE. EL_SUB1) .AND. (IELEM .LE. EL_SUP1)) THEN
                !
                IF (IINT .EQ. 1) ISINTERIOR(IELEM) = .TRUE.
                !
            ENDIF
            !
        ENDDO
        !
        CLOSE(UNIT = INPUT_UNIT)
        !
    ELSE
        !
        ! Use all elements for fiber averaging
        !
        ISINTERIOR = .TRUE.
        !
    ENDIF
    !
    ! Open output files
    !
    IF (MYID .EQ. 0) THEN
        !
        ! Open post.fib.strain-el-lat file
        !
        OPEN (OUNITS(LS_STAT_U), FILE = 'post.fib.strain-el-lat', &
            & IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.fib.'&
                &'strain-el-lat file.')
            !
        ENDIF
        !
        ! Open post.fib.defrate-pl-eq file
        !
        OPEN (OUNITS(DPEFF_STAT_U), FILE = 'post.fib.defrate-pl-eq', &
            & IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.fib.'&
                &'defrate-pl-eq file.')
            !
        ENDIF
        !
        ! Open post.fib.sliprate-sum file
        !
        OPEN (OUNITS(SGD_STAT_U), FILE = 'post.fib.sliprate-sum', &
            & IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.fib.'&
                &'sliprate-sum file.')
            !
        ENDIF
        !
        ! Open post.fib.crss file
        !
        OPEN (OUNITS(CRSS_STAT_U), FILE = 'post.fib.crss', &
            & IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.fib.'&
                &'crss file.')
            !
        ENDIF
        !
        ! Open post.fib.elt-stats file
        !
        OPEN (OUNITS(ELT_STAT_U), FILE = 'post.fib.elt-stats', &
            & IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.fib.'&
                &'elt-stats file.')
            !
        ENDIF
        !
    ENDIF
    !
    ! Open post.fib.elt.core# files
    ! Need to incremement MYID by 1 so the `core#` string is 1-indexed here
    !
    WRITE(CHARID,'(I0)') (MYID + 1)
    FILENAME = 'post.fib.elt.core'//CHARID
    !
    OPEN (OUNITS(VFE_U), FILE = FILENAME)
    !
    RETURN
    !
    END SUBROUTINE INITIALIZE_FIBER_AVERAGE
    !
    !===========================================================================
    !
    SUBROUTINE PRINT_ELEMS(ISTEP, IFIB, P_NUM_VFE_INT, ISVFE)
    !
    ! Print list of elements that meet diffraction condition or "valid fiber
    ! elements".
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! ISTEP: Index of current step number
    ! IFIB: Index of current fiber being printed
    ! P_NUM_VFE_INT: Number of "valid fiber elements" on local partition
    ! ISVFE: Logical mask of which elements in local partition are valid
    !
    INTEGER  :: ISTEP, IFIB, P_NUM_VFE_INT
    LOGICAL  :: ISVFE(EL_SUB1:EL_SUP1)
    !
    ! Locals:
    ! IO: Output file unit for printing
    ! KELEM: Generic looping index for printing all valid local elements
    !
    INTEGER  :: IO
    INTEGER  :: KELEM
    !
    !---------------------------------------------------------------------------
    !
    IO = OUNITS(VFE_U)
    !
    ! Need to 1-index the printed fiber IDs to be consistent with input
    !
    WRITE(IO,'(A2,5(I12))') '% ', ISTEP, (IFIB + 1), P_NUM_VFE_INT, &
        & EL_SUB1 + 1, EL_SUP1 + 1
    !
    ! Find and print valid fiber elements
    !
    DO KELEM = EL_SUB1, EL_SUP1
        !
        IF (ISVFE(KELEM)) THEN
            !
            WRITE(IO, '(I12)') KELEM
            !
        ENDIF
        !
    ENDDO
    !
    RETURN
    !
    END SUBROUTINE PRINT_ELEMS
    !
    !===========================================================================
    !
    SUBROUTINE PRINT_STATS(LS_AVG, LS_STD, DPEFF_AVG, DPEFF_STD, SGD_AVG, &
        & SGD_STD, CRSS_AVG, CRSS_STD, NUM_VFE, VOL_VFE, ISTEP)
    !
    ! Print statistics for elements that meet diffraction condition or "valid
    ! fiber elements"
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NUM_VFE: NFIB1 array with number of "valid fiber elements" per-fiber
    ! VOL_VFE: NFIB1 array with total "valid fiber element" volumes per-fiber
    ! LS_AVG/LS_STD: Avg. and std. deviation of elastic lattice strain
    ! DPEFF_AVG/DPEFF_STD: Avg. and std. deviation of eff. plastic def. rate.
    ! SGD_AVG/SGD_STD: Avg. and std. deviation of gammadot summation
    ! CRSS_AVG/CRSS_STD: Avg. and std. deviation of CRSS per slip system
    ! ISTEP: Index of current step number
    !
    INTEGER  :: NUM_VFE(0:NFIB1)
    REAL(RK) :: VOL_VFE(0:NFIB1)
    REAL(RK) :: LS_AVG(0:NFIB1), LS_STD(0:NFIB1)
    REAL(RK) :: DPEFF_AVG(0:NFIB1), DPEFF_STD(0:NFIB1)
    REAL(RK) :: SGD_AVG(0:NFIB1), SGD_STD(0:NFIB1)
    REAL(RK) :: CRSS_AVG(0:NFIB1,0:MAXSLIP1), CRSS_STD(0:NFIB1,0:MAXSLIP1)
    INTEGER  :: ISTEP
    !
    ! Locals:
    ! I: Generic looping index
    !
    INTEGER :: I
    !
    !------------------------------------------------------------------
    !
    IF (MYID .EQ. 0) THEN
        !
        ! Write output file headers per-step
        !
        WRITE(OUNITS(ELT_STAT_U),'(A2,2(I12))') '% ', ISTEP, NFIB
        WRITE(OUNITS(LS_STAT_U),'(A2,2(I12))') '% ', ISTEP, NFIB
        WRITE(OUNITS(DPEFF_STAT_U),'(A2,2(I12))') '% ', ISTEP, NFIB
        WRITE(OUNITS(SGD_STAT_U),'(A2,2(I12))') '% ', ISTEP, NFIB
        WRITE(OUNITS(CRSS_STAT_U),'(A2,2(I12))') '% ', ISTEP, NFIB
        !
        DO I = 0, NFIB1
            !
            WRITE(OUNITS(ELT_STAT_U), '(I8,E14.6)') NUM_VFE(I), VOL_VFE(I)
            WRITE(OUNITS(LS_STAT_U), '(2(E14.6))') LS_AVG(I), LS_STD(I)
            WRITE(OUNITS(DPEFF_STAT_U), '(2(E14.6))') DPEFF_AVG(I), DPEFF_STD(I)
            WRITE(OUNITS(SGD_STAT_U), '(2(E14.6))') SGD_AVG(I), SGD_STD(I)
            !
            SELECT CASE (CTYPE(PHASEFIB(I))%CLASS)
                !
                CASE (1, 2) ! FCC/BCC
                    ! 
                    WRITE(OUNITS(CRSS_STAT_U), '(24(E14.6))') &
                        & CRSS_AVG(I,0:11), CRSS_STD(I,0:11)
                    !
                CASE (3) ! HCP
                    !
                    WRITE(OUNITS(CRSS_STAT_U), '(36(E14.6))') &
                        & CRSS_AVG(I,0:17), CRSS_STD(I,0:17)
                    !
                CASE DEFAULT
                    !
                    CALL PAR_QUIT('Error  :     > PRINT_STATS:&
                        & Invalid crystal type')
                    !
                ! END CASE
                !
            END SELECT
            !
        END DO
        !
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE PRINT_STATS
    !
    !===========================================================================
    !
    SUBROUTINE INITIALIZE_CRYSTAL_SYMS()
    !
    ! Initialize symmetry operations
    !
    !---------------------------------------------------------------------------
    !
    ! Locals:
    ! Z,P1,M1,P32,M32,P12,M12: Constant scalar values for RMATs
    ! RMATSYMCUB_DAT: 1D array of RMAT input data for cubic symmetry
    ! RMATSYMHEX_DAT: 1D array of RMAT input data for hexagonal symmetry
    ! RMATSYMCUB: 3 x 3 x NSYMCUB rotation matrices for cubic symmetry
    ! RMATSYMHEX: 3 x 3 x NSYMHEX rotation matrices for hexagonal symmetry
    !
    REAL(RK), PARAMETER :: Z = 0.0D0
    REAL(RK), PARAMETER :: P1 = 1.0D0
    REAL(RK), PARAMETER :: M1 = -1.0D0
    REAL(RK), PARAMETER :: P32 = DSQRT(3.0D0) / 2.0D0
    REAL(RK), PARAMETER :: M32 = -P32
    REAL(RK), PARAMETER :: P12 = 0.5D0
    REAL(RK), PARAMETER :: M12 = -0.5D0
    !
    REAL(RK), PARAMETER, DIMENSION(216) :: RMATSYMCUB_DAT = (/&
        & P1,  Z,  Z,  Z, P1,  Z,  Z,  Z, P1, &
         & P1,  Z,  Z,  Z,  Z, M1,  Z, P1,  Z, &
         & P1,  Z,  Z,  Z, M1,  Z,  Z,  Z, M1, &
         & P1,  Z,  Z,  Z,  Z, P1,  Z, M1,  Z, &
         &  Z,  Z, P1,  Z, P1,  Z, M1,  Z,  Z, &
         & M1,  Z,  Z,  Z, P1,  Z,  Z,  Z, M1, &
         &  Z,  Z, M1,  Z, P1,  Z, P1,  Z,  Z, &
         &  Z, M1,  Z, P1,  Z,  Z,  Z,  Z, P1, &
         & M1,  Z,  Z,  Z, M1,  Z,  Z,  Z, P1, &
         &  Z, P1,  Z, M1,  Z,  Z,  Z,  Z, P1, &
         &  Z,  Z, P1, P1,  Z,  Z,  Z, P1,  Z, &
         &  Z, P1,  Z,  Z,  Z, P1, P1,  Z,  Z, &
         &  Z, M1,  Z,  Z,  Z, P1, M1,  Z,  Z, &
         &  Z,  Z, M1, M1,  Z,  Z,  Z, P1,  Z, &
         &  Z, M1,  Z,  Z,  Z, M1, P1,  Z,  Z, &
         &  Z,  Z, P1, M1,  Z,  Z,  Z, M1,  Z, &
         &  Z,  Z, M1, P1,  Z,  Z,  Z, M1,  Z, &
         &  Z, P1,  Z,  Z,  Z, M1, M1,  Z,  Z, &
         &  Z, P1,  Z, P1,  Z,  Z,  Z,  Z, M1, &
         &  Z, M1,  Z, M1,  Z,  Z,  Z,  Z, M1, &
         &  Z,  Z, P1,  Z, M1,  Z, P1,  Z,  Z, &
         &  Z,  Z, M1,  Z, M1,  Z, M1,  Z,  Z, &
         & M1,  Z,  Z,  Z,  Z, P1,  Z, P1,  Z, &
         & M1,  Z,  Z,  Z,  Z, M1,  Z, M1,  Z/)
    !
    REAL(RK), PARAMETER, DIMENSION(108) :: RMATSYMHEX_DAT = (/&
         &   P1,   Z,   Z,   Z,  P1,   Z,   Z,   Z,  P1, &
         &  P12, M32,   Z, P32, P12,   Z,   Z,   Z,  P1, &
         &  M12, M32,   Z, P32, M12,   Z,   Z,   Z,  P1, &
         &   M1,   Z,   Z,   Z,  M1,   Z,   Z,   Z,  P1, &
         &  M12, P32,   Z, M32, M12,   Z,   Z,   Z,  P1, &
         &  P12, P32,   Z, M32, P12,   Z,   Z,   Z,  P1, &
         &   P1,   Z,   Z,   Z,  M1,   Z,   Z,   Z,  M1, &
         &  P12, P32,   Z, P32, M12,   Z,   Z,   Z,  M1, &
         &  M12, P32,   Z, P32, P12,   Z,   Z,   Z,  M1, &
         &   M1,   Z,   Z,   Z,  P1,   Z,   Z,   Z,  M1, &
         &  M12, M32,   Z, M32, P12,   Z,   Z,   Z,  M1, &
         &  P12, M32,   Z, M32, M12,   Z,   Z,   Z,  M1/)
    !
    !----------------------------------------------------------------------
    !
    RMATSYMCUB = RESHAPE(SOURCE = RMATSYMCUB_DAT, SHAPE = (/3, 3, NSYMCUB/), &
         & ORDER = (/2, 1, 3/))

    RMATSYMHEX = RESHAPE(SOURCE = RMATSYMHEX_DAT, SHAPE = (/3, 3, NSYMHEX/), &
         & ORDER = (/2, 1, 3/))
    !
    RETURN
    !
    END SUBROUTINE INITIALIZE_CRYSTAL_SYMS
    !
END MODULE FIBER_AVERAGE_MOD
