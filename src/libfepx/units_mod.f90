! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE UNITS_MOD
!
! Module to handle fortran unit numbers for printing.
!
! Contains subroutines:
! OPEN_LOG_FILES: Opens log files.
! OPEN_OUTPUT_FILES: Opens output files.
!
! From libfepx:
!
USE READ_INPUT_MOD
USE SURFACE_MOD, ONLY: FASET
!
! From libparallel:
!
USE PARALLEL_MOD
!
IMPLICIT NONE
!
!  Output units.
!
INTEGER, PRIVATE, PARAMETER :: N_OUNITS = 60
INTEGER, PRIVATE, PARAMETER :: O_UNIT_1 = 10
INTEGER, PRIVATE, PARAMETER :: O_UNIT_MAX = O_UNIT_1 + N_OUNITS - 1
INTEGER, PRIVATE :: OU
!
INTEGER, PUBLIC, DIMENSION(N_OUNITS ):: OUNITS = &
    & (/(OU, OU = O_UNIT_1, O_UNIT_MAX) /)
!
!  Unit association parameters.
!
INTEGER, PUBLIC, PARAMETER :: &
    & ANGLES_U = 1, &
    & COORDS_U = 2, &
    & VEL_U = 3, &
    & STRESS_U = 4, &
    & ELSTRAIN_U = 5, &
    & DFLT_U = 6, &
    & EQSTRAIN_U = 7,&
    & DEFF_U = 8, &
    & GAMMADOT_U = 9, &
    & C0ANGS_U = 10, &
    & FORCE_U1 = 11, &
    & FORCE_U2 = 12, &
    & FORCE_U3 = 13, &
    & FORCE_U4 = 14, &
    & FORCE_U5 = 15, &
    & FORCE_U6 = 16, &
    & STATS_U = 17, &
    & LOG_U = 18, &
    & TIMER_U = 19, &
    & DEBUG_U = 20, &    !  debugging file
    & DPEFF_U = 21, &
    & EQPLSTRAIN_U = 22, &
    & CONV_U = 23, &
    & BCS_ITER_1_U = 24, &
    & BCS_ITER_2_U = 25, &
    & BCS_ITER_LOG_U = 26, &
    & LS_STAT_U = 27, &
    & DPEFF_STAT_U = 28, &
    & SGD_STAT_U = 29, &
    & CRSS_STAT_U = 30, &
    & ELT_STAT_U = 31, &
    & VFE_U = 32, &
    & VGRAD_U = 33, &
    & DPHAT_U = 34, &
    & WPHAT_U = 35, &
    & GAMMA_U = 36, &
    & CRSS_U = 37, &
    & EQSTRESS_U = 38, &
    & REPORT_U = 39, &
    & WORK_U = 40, &
    & WORKP_U = 41, &
    & DEFRATE_U = 42, &
    & ELVOL_U = 43, &
    & WORKRATE_U = 44, &
    & WORKRATEP_U = 45, &
    & EQELSTRAIN_U = 46, &
    & PLSTRAIN_U = 47, &
    & TOTSTRAIN_U = 48, &
    & DISP_U = 49
!
!  Input units.
!
INTEGER, PRIVATE, PARAMETER :: N_IUNITS = 10
INTEGER, PRIVATE, PARAMETER :: I_UNIT_1   = O_UNIT_MAX + 1
INTEGER, PRIVATE, PARAMETER :: I_UNIT_MAX = I_UNIT_1 + N_IUNITS - 1
INTEGER, PRIVATE :: IU
!
INTEGER, PUBLIC, DIMENSION(N_IUNITS) :: IUNITS = &
    & (/(IU, IU = I_UNIT_1, I_UNIT_MAX) /)
!
!  Unit association parameters.
!
INTEGER, PUBLIC, PARAMETER :: TMPI1_U = 1
INTEGER, PUBLIC, PARAMETER :: TMPI2_U = 2
INTEGER, PUBLIC, PARAMETER :: DRV_U = 3
INTEGER, PUBLIC, PARAMETER :: MSH_U = 4
INTEGER, PUBLIC, PARAMETER :: BCS_U = 5
INTEGER, PUBLIC, PARAMETER :: ORI_U = 6
INTEGER, PUBLIC, PARAMETER :: PHASE_U = 7
!
CONTAINS
    !
    SUBROUTINE OPEN_LOG_FILES(MYID)
    !
    ! Open log files for output.  Append process identifier string to file name.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! MYID: Processor number
    !
    INTEGER, INTENT(IN) :: MYID
    !
    ! Locals:
    !
    CHARACTER(LEN=64) :: FILENAME
    CHARACTER(LEN=8)  :: CHARID ! assumes less than 10,000 processes
    LOGICAL :: FILE_EXISTS
    INTEGER :: RST_NUM
    CHARACTER(LEN=8) :: RST_NUM_STR
    !
    !---------------------------------------------------------------------------
    !
    !  Set up process identifier string to append to filenames.
    !  A limit of 10,000 processors is imposed.
    !
    WRITE(CHARID, '(i0)') MYID
    !
    SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
        !
        CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
            !
            FILENAME = 'post.log.core'//CHARID
            OPEN(OUNITS(LOG_U), FILE = FILENAME)
            !
        CASE (0) ! APPEND
            !
            FILENAME = 'post.log.core'//CHARID
            OPEN(OUNITS(LOG_U), FILE = FILENAME, &
                & STATUS='UNKNOWN', ACCESS='APPEND')
            !
        CASE (1) ! NEW_FILE
            !
            ! Find max value N of rstN.control
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
                    CALL PAR_QUIT('Error  :     > Restart control file not &
                        & found.')
                    !
                END IF
                !
            END DO
            !
            ! Increase N to write new files for next simulation cycle
            !
            RST_NUM = RST_NUM + 2
            WRITE(RST_NUM_STR, '(I0)') RST_NUM
            FILENAME = 'post.log.rst'//TRIM(RST_NUM_STR)//'.core'//CHARID
            !
            OPEN(OUNITS(LOG_U), FILE = FILENAME)
            !
        CASE DEFAULT ! UNSUPPORTED VALUE
            !
            CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                & option provided.')
            !
    END SELECT
    !
    ! Debug files are suppressed for general users. Developer usage only.
    ! FILENAME = 'post.debug.core'//CHARID
    ! OPEN(OUNITS(DEBUG_U), FILE = FILENAME)
    !
    END SUBROUTINE OPEN_LOG_FILES
    !
    !===========================================================================
    !
    SUBROUTINE OPEN_OUTPUT_FILES(MYID)
    !
    ! Open files for output.  Append process identifier string to file name.
    !
    !---------------------------------------------------------------------------
    !
    ! Arugments
    ! MYID: Processor number
    !
    INTEGER, INTENT(IN) :: MYID
    !
    ! Locals:
    !
    INTEGER :: IOSTATUS
    CHARACTER(LEN=64):: FILENAME
    CHARACTER(LEN=8) :: CHARID ! assumes less than 10,000 processes
    CHARACTER(LEN=256) :: MESSAGE
    LOGICAL :: FILE_EXISTS
    INTEGER :: RST_NUM, I
    CHARACTER(LEN=8) :: RST_NUM_STR
    !
    !---------------------------------------------------------------------------
    !
    !  Set up process identifier string to append to FILENAMEs.
    !  A limit of 10,000 processors is imposed.
    !
    ! Print 1-indexed process numbers
    !
    WRITE(CHARID, '(i0)') MYID + 1
    !
    ! If restarting simulation, get new restart cycle number
    ! Find max value N of rstN.control
    !
    IF (OPTIONS%RESTART) THEN
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
                CALL PAR_QUIT('Error  :     > Restart control file not &
                    & found.')
                !
            END IF
            !
        END DO
        !
        ! Increase N to write new files for next simulation cycle
        !
        RST_NUM = RST_NUM + 2
        WRITE(RST_NUM_STR, '(I0)') RST_NUM
        !
    END IF
    !
    ! For all processes
    !
    IF (PRINT_OPTIONS%PRINT_COORDINATES) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.coo.core'//CHARID
                OPEN(OUNITS(COORDS_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.coo.core'//CHARID
                OPEN(OUNITS(COORDS_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.coo.rst'//TRIM(RST_NUM_STR)//'.core'//CHARID
                OPEN(OUNITS(COORDS_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_DISPLACEMENTS) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.disp.core'//CHARID
                OPEN(OUNITS(DISP_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.disp.core'//CHARID
                OPEN(OUNITS(DISP_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.disp.rst'//TRIM(RST_NUM_STR)//'.core'//CHARID
                OPEN(OUNITS(DISP_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_CRSS) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.crss.core'//CHARID
                OPEN(OUNITS(CRSS_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.crss.core'//CHARID
                OPEN(OUNITS(CRSS_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.crss.rst'//TRIM(RST_NUM_STR)//'.core'//CHARID
                OPEN(OUNITS(CRSS_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_DEFF) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.defrate-eq.core'//CHARID
                OPEN(OUNITS(DEFF_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.defrate-eq.core'//CHARID
                OPEN(OUNITS(DEFF_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.defrate-eq.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(DEFF_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_DPEFF) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.defrate-pl-eq.core'//CHARID
                OPEN(OUNITS(DPEFF_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.defrate-pl-eq.core'//CHARID
                OPEN(OUNITS(DPEFF_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.defrate-pl-eq.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(DPEFF_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_DPHAT) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.defrate-pl.core'//CHARID
                OPEN(OUNITS(DPHAT_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.defrate-pl.core'//CHARID
                OPEN(OUNITS(DPHAT_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.defrate-pl.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(DPHAT_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_ELVOL) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.elt-vol.core'//CHARID
                OPEN(OUNITS(ELVOL_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.elt-vol.core'//CHARID
                OPEN(OUNITS(ELVOL_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.elt-vol.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(ELVOL_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_EQPLSTRAIN) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.strain-pl-eq.core'//CHARID
                OPEN(OUNITS(EQPLSTRAIN_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.strain-pl-eq.core'//CHARID
                OPEN(OUNITS(EQPLSTRAIN_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.strain-pl-eq.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(EQPLSTRAIN_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_EQELSTRAIN) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.strain-el-eq.core'//CHARID
                OPEN(OUNITS(EQELSTRAIN_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.strain-el-eq.core'//CHARID
                OPEN(OUNITS(EQELSTRAIN_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.strain-el-eq.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(EQELSTRAIN_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_EQSTRAIN) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.strain-eq.core'//CHARID
                OPEN(OUNITS(EQSTRAIN_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.strain-eq.core'//CHARID
                OPEN(OUNITS(EQSTRAIN_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.strain-eq.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(EQSTRAIN_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_EQSTRESS) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.stress-eq.core'//CHARID
                OPEN(OUNITS(EQSTRESS_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.stress-eq.core'//CHARID
                OPEN(OUNITS(EQSTRESS_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.stress-eq.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(EQSTRESS_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_GAMMA) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.slip.core'//CHARID
                OPEN(OUNITS(GAMMA_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.slip.core'//CHARID
                OPEN(OUNITS(GAMMA_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.slip.rst'//TRIM(RST_NUM_STR)//'.core'//CHARID
                OPEN(OUNITS(GAMMA_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_GAMMADOT) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.sliprate.core'//CHARID
                OPEN(OUNITS(GAMMADOT_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.sliprate.core'//CHARID
                OPEN(OUNITS(GAMMADOT_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.sliprate.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(GAMMADOT_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_ORIENTATIONS) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.ori.core'//CHARID
                OPEN(OUNITS(ANGLES_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.ori.core'//CHARID
                OPEN(OUNITS(ANGLES_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.ori.rst'//TRIM(RST_NUM_STR)//'.core'//CHARID
                OPEN(OUNITS(ANGLES_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_STRAIN_TOT) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.strain.core'//CHARID
                OPEN(OUNITS(TOTSTRAIN_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.strain.core'//CHARID
                OPEN(OUNITS(TOTSTRAIN_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.strain.rst'//TRIM(RST_NUM_STR)//'.core'//CHARID
                OPEN(OUNITS(TOTSTRAIN_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_STRAIN_EL) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.strain-el.core'//CHARID
                OPEN(OUNITS(ELSTRAIN_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.strain-el.core'//CHARID
                OPEN(OUNITS(ELSTRAIN_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.strain-el.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(ELSTRAIN_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_STRAIN_PL) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.strain-pl.core'//CHARID
                OPEN(OUNITS(PLSTRAIN_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.strain-pl.core'//CHARID
                OPEN(OUNITS(PLSTRAIN_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.strain-pl.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(PLSTRAIN_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_STRESS) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.stress.core'//CHARID
                OPEN(OUNITS(STRESS_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.stress.core'//CHARID
                OPEN(OUNITS(STRESS_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.stress.rst'//TRIM(RST_NUM_STR)//'.core'//CHARID
                OPEN(OUNITS(STRESS_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_VELOCITIES) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.vel.core'//CHARID
                OPEN(OUNITS(VEL_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.vel.core'//CHARID
                OPEN(OUNITS(VEL_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.vel.rst'//TRIM(RST_NUM_STR)//'.core'//CHARID
                OPEN(OUNITS(VEL_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_VGRAD) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.velgrad.core'//CHARID
                OPEN(OUNITS(VGRAD_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.velgrad.core'//CHARID
                OPEN(OUNITS(VGRAD_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.velgrad.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(VGRAD_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_WPHAT) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.spinrate.core'//CHARID
                OPEN(OUNITS(WPHAT_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.spinrate.core'//CHARID
                OPEN(OUNITS(WPHAT_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.spinrate.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(WPHAT_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_WORK) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.work.core'//CHARID
                OPEN(OUNITS(WORK_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.work.core'//CHARID
                OPEN(OUNITS(WORK_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.work.rst'//TRIM(RST_NUM_STR)//'.core'//CHARID
                OPEN(OUNITS(WORK_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_WORKP) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.work-pl.core'//CHARID
                OPEN(OUNITS(WORKP_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.work-pl.core'//CHARID
                OPEN(OUNITS(WORKP_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.work-pl.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(WORKP_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_DEFRATE) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.defrate.core'//CHARID
                OPEN(OUNITS(DEFRATE_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.defrate.core'//CHARID
                OPEN(OUNITS(DEFRATE_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.defrate.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(DEFRATE_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_WORKRATE) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.workrate.core'//CHARID
                OPEN(OUNITS(WORKRATE_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.workrate.core'//CHARID
                OPEN(OUNITS(WORKRATE_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.workrate.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(WORKRATE_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    IF (PRINT_OPTIONS%PRINT_WORKRATEP) THEN
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.workrate-pl.core'//CHARID
                OPEN(OUNITS(WORKRATEP_U), FILE = FILENAME)
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.workrate-pl.core'//CHARID
                OPEN(OUNITS(WORKRATEP_U), FILE = FILENAME, &
                    & STATUS='UNKNOWN', ACCESS='APPEND')
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.workrate-pl.rst'//TRIM(RST_NUM_STR)&
                    &//'.core'//CHARID
                OPEN(OUNITS(WORKRATEP_U), FILE = FILENAME)
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    ! Only for main process
    !
    IF (MYID .EQ. 0) THEN
        !
        ! post.force.#
        !
        IF (PRINT_OPTIONS%PRINT_FORCES) THEN
            !
            SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
                !
                CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                    !
                    DO I = 1, 6
                        !
                        ! Loop over the 6 force units (offset needed for int)
                        !
                        OPEN(OUNITS((FORCE_U1-1) + I), &
                            & FILE=('post.force.' // FASET(I)), IOSTAT=IOSTATUS)
                        !
                        IF (IOSTATUS .NE. 0) THEN
                            !
                            WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO &
                                &Failure to open post.force.', &
                                & TRIM(FASET(I)) ,' file'
                            CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                            !
                        END IF
                        !
                    END DO
                    !
                CASE (0) ! APPEND
                    !
                    DO I = 1, 6
                        !
                        ! Loop over the 6 force units (offset needed for int)
                        !
                        OPEN(OUNITS((FORCE_U1-1) + I), &
                            & FILE=('post.force.' // FASET(I)), &
                            & IOSTAT=IOSTATUS, STATUS='UNKNOWN', &
                            & ACCESS='APPEND')
                        !
                        IF (IOSTATUS .NE. 0) THEN
                            !
                            WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO &
                                &Failure to open post.force.', &
                                & TRIM(FASET(I)) ,' file'
                            CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                            !
                        END IF
                        !
                    END DO
                    !
                CASE (1) ! NEW_FILE
                    !
                    DO I = 1, 6
                        !
                        FILENAME = 'post.force.rst'//TRIM(RST_NUM_STR)&
                            &//'.'//FASET(I)
                        OPEN(OUNITS((FORCE_U1-1) + I), &
                            & FILE=FILENAME, IOSTAT=IOSTATUS)
                        !
                        IF (IOSTATUS .NE. 0) THEN
                            !
                            WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO &
                                &Failure to open ', TRIM(FILENAME) ,' file.'
                            CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                            !
                        END IF
                        !
                    END DO
                    !
                CASE DEFAULT ! UNSUPPORTED VALUE
                    !
                    CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                        & option provided.')
                    !
            END SELECT
            !
        END IF
        !
        ! post.stats
        !
        ! OPEN(OUNITS(STATS_U), FILE='post.stats',IOSTAT=IOSTATUS)
        !
        !  IF (IOSTATUS .NE. 0) THEN
            !
            ! CALL PAR_QUIT('Error  :     > IO Failure to open post.stats file.')
            !
        ! END IF
        !
        ! post.conv
        !
        IF (PRINT_OPTIONS%PRINT_CONV) THEN
            !
            SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
                !
                CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                    !
                    FILENAME = 'post.conv'
                    OPEN(OUNITS(CONV_U), FILE = FILENAME, IOSTAT = IOSTATUS)
                    !
                    IF (IOSTATUS .NE. 0) THEN
                        !
                        WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO &
                            &Failure to open ', TRIM(FILENAME) ,' file.'
                        CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                        !
                    END IF
                    !
                CASE (0) ! APPEND
                    !
                    FILENAME = 'post.conv'
                    OPEN(OUNITS(CONV_U), FILE = FILENAME, IOSTAT = IOSTATUS)
                    !
                    IF (IOSTATUS .NE. 0) THEN
                        !
                        WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO &
                            &Failure to open ', TRIM(FILENAME) ,' file.'
                        CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                        !
                    END IF
                    !
                CASE (1) ! NEW_FILE
                    !
                    FILENAME = 'post.conv.rst'//TRIM(RST_NUM_STR)
                    OPEN(OUNITS(CONV_U), FILE = FILENAME, IOSTAT = IOSTATUS)
                    !
                    IF (IOSTATUS .NE. 0) THEN
                        !
                        WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO &
                            &Failure to open ', TRIM(FILENAME) ,' file.'
                        CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                        !
                    END IF
                    !
                CASE DEFAULT ! UNSUPPORTED VALUE
                    !
                    CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                        & option provided.')
                    !
            END SELECT
            !
        END IF
        !
        !
        ! post.report
        !
        SELECT CASE (OPTIONS%RESTART_FILE_HANDLING)
            !
            CASE (-1) ! NO RESTART - LEGACY BEHAVIOR
                !
                FILENAME = 'post.report'
                OPEN(OUNITS(REPORT_U), FILE = FILENAME, IOSTAT = IOSTATUS)
                !
                IF (IOSTATUS .NE. 0) THEN
                    !
                    WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO &
                        &Failure to open ', TRIM(FILENAME) ,' file.'
                    CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                    !
                END IF
                !
            CASE (0) ! APPEND
                !
                FILENAME = 'post.report'
                OPEN(OUNITS(REPORT_U), FILE = FILENAME, IOSTAT = IOSTATUS)
                !
                IF (IOSTATUS .NE. 0) THEN
                    !
                    WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO &
                        &Failure to open ', TRIM(FILENAME) ,' file.'
                    CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                    !
                END IF
                !
            CASE (1) ! NEW_FILE
                !
                FILENAME = 'post.report.rst'//TRIM(RST_NUM_STR)
                OPEN(OUNITS(REPORT_U), FILE = FILENAME, IOSTAT = IOSTATUS)
                !
                IF (IOSTATUS .NE. 0) THEN
                    !
                    WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO &
                        &Failure to open ', TRIM(FILENAME) ,' file.'
                    CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                    !
                END IF
                !
            CASE DEFAULT ! UNSUPPORTED VALUE
                !
                CALL PAR_QUIT('Error  :     > Invalid restart file handling&
                    & option provided.')
                !
        END SELECT
        !
    END IF
    !
    END SUBROUTINE OPEN_OUTPUT_FILES
    !
END MODULE UNITS_MOD
