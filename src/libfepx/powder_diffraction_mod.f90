! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE POWDER_DIFFRACTION_MOD
!
! Performs "powder diffraction" analysis on the data at the end of each step in
!   the load history. Essentially calculates fiber-averaged values of certain
!   metrics.
!
! Contains subroutines:
! RUN_POWDER_DIFFRACTION: Run powder diffraction routine
! INITIALIZE_POWDER_DIFFRACTION: Initialize powder diffraction routine
! PRINT_ELEMS: Print litup elements for each processor
! PRINT_STATS: Print lightup statistics
! INITIALIZE_CRYSTAL_SYMS: Initialize crystal symmetry operations
!
USE GATHER_SCATTER
USE LIBF95
USE PARALLEL_MOD
!
USE DIMSMODULE
USE READ_INPUT_MOD
USE MICROSTRUCTURE_MOD
USE SIMULATION_CONFIGURATION_MOD
USE UNITS_MOD
USE UTILSCRYSTALMODULE
!
IMPLICIT NONE
!
! Private data
!
! HKLFIB: 3 x NFIB array of crystal directions
! UVWFIB: 3 x NFIB array of sample directions
! PHASEFIB: NFIB array of phase for each fiber
! TOLFIB: NFIB array of angular TOLerance for each fiber
! ELVOL: element volumes
! ISINTERIOR: logial array of whether each element is in mesh interior
!
INTEGER,  PRIVATE  :: NFIB, NFIB1
INTEGER,  ALLOCATABLE, PRIVATE  :: HKLFIB(:, :)
REAL(RK), ALLOCATABLE, PRIVATE  :: UVWFIB(:, :)
INTEGER,  ALLOCATABLE, PRIVATE  :: PHASEFIB(:)
REAL(RK), ALLOCATABLE, PRIVATE  :: TOLFIB(:)
REAL(RK), ALLOCATABLE, PRIVATE  :: ELVOL(:)
LOGICAL,  ALLOCATABLE, PRIVATE  :: ISINTERIOR(:)
!
INTEGER,  PARAMETER, PRIVATE :: NSYMCUB = 24, NSYMCUB1 = 23
INTEGER,  PARAMETER, PRIVATE :: NSYMHEX = 12, NSYMHEX1 = 11
REAL(RK), PRIVATE :: RMATSYMCUB(0:DIMS1, 0:DIMS1, 0:NSYMCUB1)
REAL(RK), PRIVATE :: RMATSYMHEX(0:DIMS1, 0:DIMS1, 0:NSYMHEX1)
!
CONTAINS
    !
    SUBROUTINE RUN_POWDER_DIFFRACTION(ISTEP, C_ANGS, ELAS_TOT6, DPEFF, CRSS)
    !
    ! Run powder diffraction routine
    !
    !---------------------------------------------------------------------------
    !
    IMPLICIT  NONE
    !
    ! Arguments:
    ! ISTEP: Index of current step number
    ! C_ANGS: Current elemental orientations
    ! ELAS_TOT6: Current elemental strains (6-vector)
    ! DPEFF: Current effective plastic deformation rate
    ! CRSS: Current slip systen strength
    !
    INTEGER,  INTENT(IN) :: ISTEP
    REAL(RK), INTENT(IN) :: C_ANGS(0:DIMS1, 0:DIMS1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: ELAS_TOT6(0:TVEC , 0:NGRAIN1, EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: DPEFF(EL_SUB1:EL_SUP1)
    REAL(RK), INTENT(IN) :: CRSS(0:MAXSLIP1, 0:NGRAIN1, EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    INTEGER, PARAMETER  :: NSYMCUB = 24, NSYMCUB1 = 23
    INTEGER, PARAMETER  :: NSYMHEX = 12, NSYMHEX1 = 11
    REAL(RK), PARAMETER :: TOL = 1e-4
    INTEGER  :: NUM_LIT(0:NFIB1)
    REAL(RK) :: VOL_LIT(0:NFIB1)
    REAL(RK) :: LS_AVG(0:NFIB1), LS_STD(0:NFIB1)
    REAL(RK) :: DPEFF_AVG(0:NFIB1), DPEFF_STD(0:NFIB1)
    REAL(RK) :: SGD_AVG(0:NFIB1), SGD_STD(0:NFIB1)
    REAL(RK) :: CRSS_AVG(0:NFIB1), CRSS_STD(0:NFIB1)
    INTEGER  :: P_NUM_LIT_INT
    REAL(RK) :: P_NUM_LIT, NUM_LIT_IFIB ! NUM_LIT is real for use with par_sum
    REAL(RK) :: P_VOL_LIT, VOL_LIT_IFIB
    REAL(RK) :: P_LS_SUM, LS_SUM
    REAL(RK) :: P_DPEFF_SUM, DPEFF_SUM
    REAL(RK) :: P_SGD_SUM, P_CRSS_SUM
    REAL(RK) :: SGD_SUM, CRSS_SUM
    REAL(RK) :: LS(EL_SUB1:EL_SUP1)
    REAL(RK) :: SGD(EL_SUB1:EL_SUP1)
    INTEGER :: I, J, K, KELEM, ISS, IFIB, IPLN
    INTEGER :: NSYM1
    INTEGER :: NUNQPLN
    LOGICAL :: ISUNIQUEPLN
    LOGICAL :: ISLIT(EL_SUB1:EL_SUP1)
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
    LS_AVG = 0.0_RK
    DPEFF_AVG = 0.0_RK
    SGD_AVG = 0.0_RK
    CRSS_AVG = 0.0_RK
    LS_STD = 0.0_RK
    DPEFF_STD = 0.0_RK
    SGD_STD = 0.0_RK
    CRSS_STD = 0.0_RK
    !
    ! Calculate sum of the gammadots
    !
    SGD = 0.0_RK
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
        LS = 0.0_RK
        P_NUM_LIT = 0.0_RK
        P_VOL_LIT = 0.0_RK
        P_LS_SUM = 0.0_RK
        P_DPEFF_SUM = 0.0_RK
        P_SGD_SUM = 0.0_RK
        P_CRSS_SUM = 0.0_RK
        NUM_LIT_IFIB = 0.0_RK
        VOL_LIT_IFIB = 0.0_RK
        LS_SUM = 0.0_RK
        DPEFF_SUM = 0.0_RK
        SGD_SUM = 0.0_RK
        CRSS_SUM = 0.0_RK
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
                !   crystal basis
                !
                PLN = HKL/SQRT(HKL(0)*HKL(0) + HKL(1)*HKL(1) + HKL(2)*HKL(2))
                NSYM1 = NSYMCUB1
                ALLOCATE(RMATSYM(0:DIMS1, 0:DIMS1, 0:NSYM1))
                RMATSYM = RMATSYMCUB
                !
            CASE (3)
                !
                ! Hexagonal:
                ! construct lattice plane unit normal vectors (PLN) in Cartesian
                !   crystal basis
                !
                PLN(0) = HKL(0)
                PLN(1) = (2*HKL(1)+HKL(0))/RK_ROOT_3
                PLN(2) = HKL(2)/crystal_parm(11,PHASEFIB(IFIB))
                PLN = PLN/SQRT(PLN(0)*PLN(0) + PLN(1)*PLN(1) + PLN(2)*PLN(2))
                NSYM1 = NSYMHEX1
                ALLOCATE(RMATSYM(0:DIMS1, 0:DIMS1, 0:NSYM1))
                RMATSYM = RMATSYMHEX
                !
            CASE DEFAULT
                !
                CALL PAR_QUIT('Error  :     > RUN_POWDER_DIFFRACTION: &
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
                IF (ABS(1.0_RK - ABS(DOTPROD)) .LT. TOL) ISUNIQUEPLN = .FALSE.
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
        ISLIT = .FALSE.
        DO KELEM = EL_SUB1, EL_SUP1 ! Loop over all elements
            !
            IF (ISINTERIOR(KELEM) .AND. (PHASE(KELEM) .EQ. PHASEFIB(IFIB))) THEN
                !
                DOTPROD = 0.0_RK
                !
                DO IPLN = 0, (NUNQPLN - 1)
                    !
                    SVEC = 0.0_RK
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
                    ISLIT(KELEM) = .TRUE.
                    !
                    ! Calculate lattice strain
                    !
                    ! Scattering vector in crystal coordinates
                    ! {n}_c = [R]^T {n}_s
                    !
                    CVEC = 0.0_RK
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
                    P_NUM_LIT = P_NUM_LIT + 1.0_RK
                    P_VOL_LIT = P_VOL_LIT + ELVOL(KELEM)
                    !
                    ! The following sums are weighted by element volume
                    !
                    P_LS_SUM = P_LS_SUM + LS(KELEM)*ELVOL(KELEM)
                    P_DPEFF_SUM = P_DPEFF_SUM + DPEFF(KELEM)*ELVOL(KELEM)
                    P_SGD_SUM = P_SGD_SUM + SGD(KELEM)*ELVOL(KELEM)
                    !
                    DO ISS =0,MAXSLIP1
                        !
                        P_CRSS_SUM = P_CRSS_SUM + CRSS(ISS, 0, KELEM) * &
                            & ELVOL(KELEM)
                        !
                    ENDDO
                    !
                ENDIF
                !
            ENDIF
            !
        ENDDO ! Loop over all elements
        !
        ! Calculate fiber-averaged statistics
        !
        CALL PAR_SUM(P_NUM_LIT, NUM_LIT_IFIB)
        CALL PAR_SUM(P_VOL_LIT, VOL_LIT_IFIB)
        CALL PAR_SUM(P_LS_SUM, LS_SUM)
        CALL PAR_SUM(P_DPEFF_SUM, DPEFF_SUM)
        CALL PAR_SUM(P_SGD_SUM, SGD_SUM)
        CALL PAR_SUM(P_CRSS_SUM, CRSS_SUM)
        !
        NUM_LIT(IFIB) = NINT(NUM_LIT_IFIB)
        VOL_LIT(IFIB) = VOL_LIT_IFIB
        !
        IF (VOL_LIT_IFIB .GT. 0.0) THEN
            !
            LS_AVG(IFIB) = LS_SUM / VOL_LIT_IFIB
            DPEFF_AVG(IFIB) = DPEFF_SUM / VOL_LIT_IFIB
            SGD_AVG(IFIB) = SGD_SUM / VOL_LIT_IFIB
            CRSS_AVG(IFIB) = CRSS_SUM / VOL_LIT_IFIB / MAXSLIP1
            !
        ENDIF
        !
        ! Calculate standard deviations
        !
        P_LS_SUM = 0.0_RK
        P_DPEFF_SUM = 0.0_RK
        P_SGD_SUM = 0.0_RK
        P_CRSS_SUM = 0.0_RK
        LS_SUM = 0.0_RK
        DPEFF_SUM = 0.0_RK
        SGD_SUM = 0.0_RK
        CRSS_SUM = 0.0_RK
        !
        IF (VOL_LIT_IFIB .GT. 0.0) THEN
            !
            DO KELEM = EL_SUB1, EL_SUP1
                !
                IF (ISLIT(KELEM)) THEN
                    !
                    P_LS_SUM = P_LS_SUM + (LS(KELEM) - LS_AVG(IFIB))**2 * &
                        & ELVOL(KELEM)
                    P_DPEFF_SUM = P_DPEFF_SUM + (DPEFF(KELEM) - &
                        & DPEFF_AVG(IFIB))**2 * ELVOL(KELEM)
                    P_SGD_SUM = P_SGD_SUM + (SGD(KELEM) - SGD_AVG(IFIB))**2 * &
                        & ELVOL(KELEM)
                    !
                    DO ISS = 0, MAXSLIP1
                        !
                        P_CRSS_SUM = P_CRSS_SUM + (CRSS(ISS, 0,  KELEM) - &
                            & CRSS_AVG(IFIB))**2 * ELVOL(KELEM)
                        !
                    ENDDO
                    !
                ENDIF
                !
            ENDDO
            !
            CALL PAR_SUM(P_LS_SUM, LS_SUM)
            CALL PAR_SUM(P_DPEFF_SUM, DPEFF_SUM)
            CALL PAR_SUM(P_SGD_SUM, SGD_SUM)
            CALL PAR_SUM(P_CRSS_SUM, CRSS_SUM)
            !
            LS_STD(IFIB) = SQRT(LS_SUM/VOL_LIT_IFIB)
            DPEFF_STD(IFIB) = SQRT(DPEFF_SUM/VOL_LIT_IFIB)
            SGD_STD(IFIB) = SQRT(SGD_SUM/VOL_LIT_IFIB)
            CRSS_STD(IFIB) = SQRT(CRSS_SUM/VOL_LIT_IFIB/MAXSLIP1)
            !
        ENDIF
        !
        DEALLOCATE(ALLPLN)
        DEALLOCATE(RMATSYM)
        !
        ! Print list of lit elements
        !
        P_NUM_LIT_INT = NINT(P_NUM_LIT)
        !
        CALL PRINT_ELEMS(ISTEP, IFIB, P_NUM_LIT_INT, ISLIT)
        !
    ENDDO ! Loop over all fibers
    !
    ! Print statsitics
    !
    CALL PRINT_STATS(LS_AVG, LS_STD, DPEFF_AVG, DPEFF_STD, SGD_AVG, SGD_STD, &
        & CRSS_AVG, CRSS_STD, NUM_LIT, VOL_LIT)
    !
    RETURN
    !
    END SUBROUTINE RUN_POWDER_DIFFRACTION
    !
    !===========================================================================
    !
    SUBROUTINE INITIALIZE_POWDER_DIFFRACTION(INPUT_UNIT, DTRACE)
    !
    ! Read and process input data in powder diffraction file, open files for
    !   output, initialize symmetry operations
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! INPUT_UNIT: Unit for input file
    ! DTRACE: Gather/scatter traces
    !
    INTEGER  :: INPUT_UNIT
    TYPE(trace) :: DTRACE
    !
    ! Locals:
    !
    INTEGER  :: I, J, IELEM
    INTEGER  :: IVOL, IINT
    REAL(RK) :: ECOORDS(0:KDIM1, EL_SUB1:EL_SUP1)
    CHARACTER*128 FILENAME
    CHARACTER*4   CHARID      ! assumes less than 10,000 processes
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
    WRITE(DFLT_U, *) 'Info   :   - Powder diffraction options'
    OPEN(INPUT_UNIT, FILE = POWDER_DIFFRACTION_OPTIONS%POWDER_DIFFRACTION_FILE,&
                & STATUS = 'old', ACTION = 'read')
    !
    READ(INPUT_UNIT,*) NFIB
    !
    NFIB1 = NFIB-1
    !
    ALLOCATE(HKLFIB(0:DIMS1, 0:NFIB1))
    ALLOCATE(UVWFIB (0:DIMS1, 0:NFIB1))
    ALLOCATE(PHASEFIB (0:NFIB1))
    ALLOCATE(TOLFIB (0:NFIB1))
    ALLOCATE(ELVOL(EL_SUB1:EL_SUP1))
    ALLOCATE(ISINTERIOR(EL_SUB1:EL_SUP1))
    !
    DO I = 0, DIMS1
        !
        READ(INPUT_UNIT, *) (HKLFIB(I, J), J = 0, NFIB1)
        !
    ENDDO
    !
    DO I = 0, DIMS1
        !
        READ(INPUT_UNIT, *) (UVWFIB(I, J), J = 0, NFIB1)
        !
    ENDDO
    !
    READ(INPUT_UNIT,*) (PHASEFIB(J), J = 0, NFIB1)
    !
    READ(INPUT_UNIT,*) (TOLFIB(J), J = 0, NFIB1)
    TOLFIB = COS(TOLFIB * RK_PI_OVER_180)
    !
    CLOSE(unit = INPUT_UNIT)
    !
    ! Read element volumes
    !
    IF (POWDER_DIFFRACTION_OPTIONS%READ_VOLUME) THEN
        !
        ! Read volume file
        !
        OPEN(INPUT_UNIT, FILE = POWDER_DIFFRACTION_OPTIONS%VOLUME_FILE, &
            & STATUS = 'old', ACTION = 'read')
        !
        DO IELEM = 0, MAXEL1
            !
            READ(INPUT_UNIT, *) IVOL
            !
            IF ((IELEM .GE. EL_SUB1) .AND. (IELEM .LE. EL_SUP1)) THEN
                !
                ELVOL(IELEM) = IVOL
                !
            ENDIF
            !
        ENDDO
        !
        CLOSE(UNIT = INPUT_UNIT)
        !
    ELSE
        !
        ! Calculate element volumes
        !
        CALL PART_GATHER(ECOORDS, COORDS, NODES, DTRACE)
        CALL CALC_ELVOL(ELVOL, ECOORDS)
        !
    ENDIF
    !
    ! Read interior elements
    !
    IF (POWDER_DIFFRACTION_OPTIONS%READ_INTERIOR) THEN
        !
        ! Read interior file
        !
        OPEN(INPUT_UNIT, FILE = POWDER_DIFFRACTION_OPTIONS%INTERIOR_FILE, &
            & STATUS = 'old', ACTION = 'read')
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
        ! Use all elements for lightup script
        !
        ISINTERIOR = .TRUE.
        !
    ENDIF
    !
    ! Open output files
    !
    IF (MYID .EQ. 0) THEN
        !
        ! Open post.powderdiffraction.ls_avg file
        !
        OPEN (OUNITS(LS_AVG_U), FILE = 'post.powderdiffraction.ls_avg', &
            & IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.powder'&
                &'diffraction.ls_avg file.')
            !
        ENDIF
        !
        ! Open post.powderdiffraction.ls_std file
        !
        OPEN (OUNITS(LS_STD_U), FILE = 'post.powderdiffraction.ls_std', &
            & IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.powder'&
                &'diffraction.ls_std file.')
            !
        ENDIF
        !
        ! Open post.powderdiffraction.dpeff_avg file
        !
        OPEN (OUNITS(DPEFF_AVG_U), FILE = 'post.powderdiffraction.dpeff_avg', &
            & IOSTAT = IOSTATUS)
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.powder'&
                &'diffraction.dpeff_avg file.')
            !
        ENDIF
        !
        ! Open post.powderdiffraction.dpeff_std file
        !
        OPEN (OUNITS(DPEFF_STD_U), FILE = 'post.powderdiffraction.dpeff_std', &
            & IOSTAT = IOSTATUS)
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.powder'&
                &'diffraction.dpeff_std file.')
            !
        ENDIF
        !
        ! Open post.powderdiffraction.sgd_avg file
        !
        OPEN (OUNITS(SGD_AVG_U), FILE = 'post.powderdiffraction.sgd_avg', &
            & IOSTAT = IOSTATUS)
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.powder'&
                &'diffraction.sgd_avg file.')
            !
        ENDIF
        !
        ! Open post.powderdiffraction.sgd_std file
        !
        OPEN (OUNITS(SGD_STD_U), FILE = 'post.powderdiffraction.sgd_std', &
            & IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.powder'&
                &'diffraction.sgd_std file.')
            !
        ENDIF
        !
        ! Open post.powderdiffraction.crss_avg file
        !
        OPEN (OUNITS(CRSS_AVG_U), FILE = 'post.powderdiffraction.crss_avg', &
            & IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.powder'&
                &'diffraction.crss_avg file.')
            !
        ENDIF
        !
        ! Open post.powderdiffraction.crss_std file
        !
        OPEN (OUNITS(CRSS_STD_U), FILE = 'post.powderdiffraction.crss_std', &
            & IOSTAT = IOSTATUS)
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.powder'&
                &'diffraction.crss_std file.')
            !
        ENDIF
        !
        ! Open post.powderdiffraction.num_lit file
        !
        OPEN (OUNITS(NUM_LIT_U), FILE = 'post.powderdiffraction.num_lit', &
            & IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.powder'&
                &'diffraction.num_lit file.')
            !
        ENDIF
        !
        ! Open post.powderdiffraction.vol_lit file
        !
        OPEN (OUNITS(VOL_LIT_U), FILE = 'post.powderdiffraction.vol_lit', &
            & IOSTAT = IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > Failure to open post.powder'&
                &'diffraction.vol_lit file.')
            !
        ENDIF
        !
    ENDIF
    !
    ! Open post.powderdiffraction.lit_elems.# files
    !
    IF (MYID .LE. 9) THEN
        !
        WRITE(CHARID,'(i1)') MYID
        !
    ELSEIF (MYID .LE. 99) THEN
        !
        WRITE(CHARID,'(i2)') MYID
        !
    ELSEIF (MYID .LE. 999) THEN
        !
        WRITE(CHARID,'(i3)') MYID
        !
    ELSEIF (MYID .LE. 9999) THEN
        !
        WRITE(CHARID,'(i4)') MYID
        !
    ENDIF
    !
    FILENAME = 'post.powderdiffraction.lit_elems.'//CHARID
    !
    OPEN (OUNITS(LIT_ELEMS_U), FILE = FILENAME)
    !
    RETURN
    !
    END SUBROUTINE INITIALIZE_POWDER_DIFFRACTION
    !
    !===========================================================================
    !
    SUBROUTINE PRINT_ELEMS(ISTEP, IFIB, P_NUM_LIT, ISLIT)
    !
    ! Print list of elements that meet diffraction condition ("lit up" elements)
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER  :: ISTEP, IFIB, P_NUM_LIT
    LOGICAL  :: ISLIT(EL_SUB1:EL_SUP1)
    !
    ! Locals:
    !
    INTEGER  :: IO
    INTEGER  :: KELEM
    !
    !---------------------------------------------------------------------------
    !
    IO = OUNITS(LIT_ELEMS_U)
    !
    WRITE(IO, '(a2,4(i12))') '% ', ISTEP, IFIB, EL_SUB1, EL_SUP1
    WRITE(IO, '(i12)') P_NUM_LIT
    !
    ! Find and print lit elements
    !
    DO KELEM = EL_SUB1, EL_SUP1
        !
        IF (ISLIT(KELEM)) THEN
            !
            WRITE(IO, '(i12)') KELEM
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
        & SGD_STD, CRSS_AVG, CRSS_STD, NUM_LIT, VOL_LIT)
    !
    ! Print statistics for elements that meet diffraction condition
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    !
    INTEGER  :: NUM_LIT(0:NFIB1)
    REAL(RK) :: VOL_LIT(0:NFIB1)
    REAL(RK) :: LS_AVG(0:NFIB1), LS_STD(0:NFIB1)
    REAL(RK) :: DPEFF_AVG(0:NFIB1), DPEFF_STD(0:NFIB1)
    REAL(RK) :: SGD_AVG(0:NFIB1), SGD_STD(0:NFIB1)
    REAL(RK) :: CRSS_AVG(0:NFIB1), CRSS_STD(0:NFIB1)
    !
    ! Locals:
    !
    INTEGER :: I
    CHARACTER(128) :: REAL_STR, INT_STR
    !
    !------------------------------------------------------------------
    !
    IF (MYID .EQ. 0) THEN
        !
        WRITE(REAL_STR, '(a1,i8,a6)') '(', NFIB, 'e14.6)'
        WRITE(INT_STR,  '(a1,i8,a3)') '(', NFIB, 'i8)'
        !
        WRITE(OUNITS(LS_AVG_U), REAL_STR) LS_AVG
        WRITE(OUNITS(LS_STD_U), REAL_STR) LS_STD
        WRITE(OUNITS(DPEFF_AVG_U), REAL_STR) DPEFF_AVG
        WRITE(OUNITS(DPEFF_STD_U), REAL_STR) DPEFF_STD
        WRITE(OUNITS(SGD_AVG_U), REAL_STR) SGD_AVG
        WRITE(OUNITS(SGD_STD_U), REAL_STR) SGD_STD
        WRITE(OUNITS(CRSS_AVG_U), REAL_STR) CRSS_AVG
        WRITE(OUNITS(CRSS_STD_U), REAL_STR) CRSS_STD
        WRITE(OUNITS(NUM_LIT_U), INT_STR) NUM_LIT
        WRITE(OUNITS(VOL_LIT_U), REAL_STR) VOL_LIT
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
    !
    REAL(RK), PARAMETER :: Z = 0.0_RK
    REAL(RK), PARAMETER :: P1 = 1.0_RK
    REAL(RK), PARAMETER :: M1 = -1.0_RK
    REAL(RK), PARAMETER :: P32 = RK_ROOT_3/RK_TWO
    REAL(RK), PARAMETER :: M32 = -P32
    REAL(RK), PARAMETER :: P12 = 0.5_RK
    REAL(RK), PARAMETER :: M12 = -0.5_RK
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
    INTEGER :: I, K
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
END MODULE POWDER_DIFFRACTION_MOD
