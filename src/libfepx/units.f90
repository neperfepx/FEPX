! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
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
USE PARALLEL_MOD
USE SIMULATION_CONFIGURATION_MOD
USE SURF_INFO_MOD, ONLY: FASET
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
INTEGER, PUBLIC, PARAMETER ::&
    & ANGLES_U = 1,&
    & COORDS_U = 2,&
    & VEL_U = 3,&
    & STRESS_U = 4,&
    & STRAIN_U = 5,&
    & DFLT_U = 6,&
    & EQSTRAIN_U = 7,&
    & DEFF_U = 8,&
    & GAMMADOT_U = 9,&
    & C0ANGS_U = 10,&
    & FORCE_U1 = 11,&
    & FORCE_U2 = 12,&
    & FORCE_U3 = 13,&
    & FORCE_U4 = 14,&
    & FORCE_U5 = 15,&
    & FORCE_U6 = 16,&
    & STATS_U = 17,&
    & LOG_U = 18,&
    & TIMER_U = 19,&
    & DEBUG_U = 20,&    !  debugging file
    & DPEFF_U = 21,&
    & EQPLSTRAIN_U = 22,&
    & CONV_U = 23,&
    & BCS_ITER_1_U = 24,&
    & BCS_ITER_2_U = 25,&
    & BCS_ITER_LOG_U = 26, &
    & LS_AVG_U = 27,&
    & LS_STD_U = 28,&
    & DPEFF_AVG_U = 29,&
    & DPEFF_STD_U = 30,&
    & SGD_AVG_U = 31,&
    & SGD_STD_U = 32,&
    & CRSS_AVG_U = 33,&
    & CRSS_STD_U = 34,&
    & NUM_LIT_U = 35,&
    & VOL_LIT_U = 36,&
    & LIT_ELEMS_U = 37,&
    & VGRAD_U = 38,&
    & DPHAT_U = 39,&
    & WPHAT_U = 40,&
    & GAMMA_U = 41,&
    & CRSS_U = 42,&
    & EQSTRESS_U = 43,&
    & REPORT_U = 44,&
    & WORK_U = 45,&
    & WORKP_U = 46,&
    & DEFRATE_U = 47
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
    INTEGER :: IOSTATUS
    CHARACTER(LEN=64) :: FILENAME
    CHARACTER(LEN=8)  :: CHARID ! assumes less than 10,000 processes
    !
    !---------------------------------------------------------------------------
    !
    !  Set up process identifier string to append to filenames.
    !  A limit of 10,000 processors is imposed.
    !
    WRITE(CHARID, '(i0)') MYID
    !
    FILENAME = 'post.log.core'//CHARID
    OPEN(OUNITS(LOG_U), FILE = FILENAME)
    !
    ! Debug files are suppressed for general users. Developer usage only.
!    FILENAME = 'post.debug.core'//CHARID
!    OPEN(OUNITS(DEBUG_U), FILE = FILENAME)
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
    !
    !---------------------------------------------------------------------------
    !
    !  Set up process identifier string to append to FILENAMEs.
    !  A limit of 10,000 processors is imposed.
    !
    ! Print 1-indexed process numbers
    WRITE(CHARID, '(i0)') MYID + 1
    !
    ! For all processes
    !
    IF (PRINT_OPTIONS%PRINT_COORDINATES) THEN
        !
        FILENAME = 'post.coo.core'//CHARID
        OPEN(OUNITS(COORDS_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_CRSS) THEN
        !
        FILENAME = 'post.crss.core'//CHARID
        OPEN(OUNITS(CRSS_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_DEFF) THEN
        !
        FILENAME = 'post.defrate-eq.core'//CHARID
        OPEN(OUNITS(DEFF_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_DPEFF) THEN
        !
        FILENAME = 'post.defrate-pl-eq.core'//CHARID
        OPEN(OUNITS(DPEFF_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_DPHAT) THEN
        !
        FILENAME = 'post.defrate-pl.core'//CHARID
        OPEN(OUNITS(DPHAT_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_EQPLSTRAIN) THEN
        !
        FILENAME = 'post.strain-pl-eq.core'//CHARID
        OPEN(OUNITS(EQPLSTRAIN_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_EQSTRAIN) THEN
        !
        FILENAME = 'post.strain-eq.core'//CHARID
        OPEN(OUNITS(EQSTRAIN_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_EQSTRESS) THEN
        !
        FILENAME = 'post.stress-eq.core'//CHARID
        OPEN(OUNITS(EQSTRESS_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_GAMMA) THEN
        !
        FILENAME = 'post.slip.core'//CHARID
        OPEN(OUNITS(GAMMA_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_GAMMADOT) THEN
        !
        FILENAME = 'post.sliprate.core'//CHARID
        OPEN(OUNITS(GAMMADOT_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_ORIENTATIONS) THEN
        !
        FILENAME = 'post.ori.core'//CHARID
        OPEN(OUNITS(ANGLES_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_STRAIN) THEN
        !
        FILENAME = 'post.strain-el.core'//CHARID
        OPEN(OUNITS(STRAIN_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_STRESS) THEN
        !
        FILENAME = 'post.stress.core'//CHARID
        OPEN(OUNITS(STRESS_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_VELOCITIES) THEN
        !
        FILENAME = 'post.vel.core'//CHARID
        OPEN(OUNITS(VEL_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_VGRAD) THEN
        !
        FILENAME = 'post.velgrad.core'//CHARID
        OPEN(OUNITS(VGRAD_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_WPHAT) THEN
        !
        FILENAME = 'post.spinrate.core'//CHARID
        OPEN(OUNITS(WPHAT_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_WORK) THEN
        !
        FILENAME = 'post.work.core'//CHARID
        OPEN(OUNITS(WORK_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_WORKP) THEN
        !
        FILENAME = 'post.work-pl.core'//CHARID
        OPEN(OUNITS(WORKP_U), FILE = FILENAME)
        !
    ENDIF
    !
    IF (PRINT_OPTIONS%PRINT_DEFRATE) THEN
        !
        FILENAME = 'post.defrate.core'//CHARID
        OPEN(OUNITS(DEFRATE_U), FILE = FILENAME)
        !
    ENDIF
    !
    ! Only for main process
    !
    IF (MYID .EQ. 0) THEN
        !
        ! post.force#
        !
        IF (PRINT_OPTIONS%PRINT_FORCES) THEN
            !
            OPEN(OUNITS(FORCE_U1), FILE=('post.force.' // FASET(1)), &
                &IOSTAT=IOSTATUS)
            !
            IF (IOSTATUS .NE. 0) THEN
                !
                WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO Failure to open &
                    &post.force.', TRIM(FASET(1)) ,' file'
                CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                !
            ENDIF
            !
            OPEN(OUNITS(FORCE_U2), FILE=('post.force.' // FASET(2)), &
                &IOSTAT=IOSTATUS)
            !
            IF (IOSTATUS .NE. 0) THEN
                !
                WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO Failure to open &
                    &post.force.', TRIM(FASET(2)) ,' file'
                CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                !
            ENDIF
            !
            OPEN(OUNITS(FORCE_U3), FILE=('post.force.' // FASET(3)), &
                &IOSTAT=IOSTATUS)
            !
            IF (IOSTATUS .NE. 0) THEN
                !
                WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO Failure to open &
                    &post.force.', TRIM(FASET(3)) ,' file'
                CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                !
            ENDIF
            !
            OPEN(OUNITS(FORCE_U4), FILE=('post.force.' // FASET(4)), &
                &IOSTAT=IOSTATUS)
            !
            IF (IOSTATUS .NE. 0) THEN
                !
                WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO Failure to open &
                    &post.force.', TRIM(FASET(4)) ,' file'
                CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                !
            ENDIF
            !
            OPEN(OUNITS(FORCE_U5), FILE=('post.force.' // FASET(5)), &
                &IOSTAT=IOSTATUS)
            !
            IF (IOSTATUS .NE. 0) THEN
                !
                WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO Failure to open &
                    &post.force.', TRIM(FASET(5)) ,' file'
                CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                !
            ENDIF
            !
            OPEN(OUNITS(FORCE_U6), FILE=('post.force.' // FASET(6)), &
                &IOSTAT=IOSTATUS)
            !
            IF (IOSTATUS .NE. 0) THEN
                !
                WRITE(MESSAGE, '(a,a,a)') 'Error  :     > IO Failure to open &
                    &post.force.', TRIM(FASET(6)) ,' file'
                CALL PAR_QUIT(TRIM(ADJUSTL(MESSAGE)))
                !
            ENDIF
            !
        END IF
        !
        ! post.stats
        !
!        OPEN(OUNITS(STATS_U), FILE='post.stats',IOSTAT=IOSTATUS)
        !
!        IF (IOSTATUS .NE. 0) THEN
            !
!            CALL PAR_QUIT('Error  :     > IO Failure to open post.stats file.')
            !
!        ENDIF
        !
        IF (PRINT_OPTIONS%PRINT_CONV) THEN
            !
            OPEN(OUNITS(CONV_U), FILE='post.conv',IOSTAT=IOSTATUS)
            !
            IF (IOSTATUS .NE. 0) THEN
                !
                CALL PAR_QUIT('Error  :     > IO Failure to open post.conv &
                    &file.')
                !
            ENDIF
            !
        END IF
        !
        OPEN(OUNITS(REPORT_U), FILE='post.report',IOSTAT=IOSTATUS)
        !
        IF (IOSTATUS .NE. 0) THEN
            !
            CALL PAR_QUIT('Error  :     > IO Failure to open post.report file.')
            !
        ENDIF
        !
    ENDIF
    !
    END SUBROUTINE OPEN_OUTPUT_FILES
    !
END MODULE UNITS_MOD
