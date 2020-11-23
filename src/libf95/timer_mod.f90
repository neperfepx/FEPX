! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE TIMER_MOD
!
! This module handles the timer object. Uses CPU_TIME instead of SYSTEM_CLOCK.
!
! Contains subroutines:
! TIMERINIT: Initialization for the timer object
! TIMERACTIVATE: Activate timer
! TIMERSTART: Start the timer
! TIMERSTOP: Stop timer and compute elapsed time
! TIMERWRITE: Print statistics for each timer
!
IMPLICIT NONE
!
TYPE TIMERTYPE
    !
    CHARACTER(LEN=32) :: NAME = ''
    LOGICAL :: ACTIVE = .TRUE.
    LOGICAL :: FAILED = .FALSE.
    INTEGER :: STARTS = 0
    INTEGER :: STOPS = 0
    REAL :: TIME_BEG = 0.0
    REAL :: TIME_END = 0.0
    REAL :: TIME_TOTAL = 0.0
    !
END TYPE TIMERTYPE
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: TIMERTYPE
PUBLIC :: TIMERINIT
PUBLIC :: TIMERACTIVATE
PUBLIC :: TIMERSTART
PUBLIC :: TIMERSTOP
PUBLIC :: TIMERWRITE
!
CONTAINS
    !
    SUBROUTINE TIMERINIT(SELF, NAME, ACTIVATE)
    !
    !  Initialization for the timer object.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! SELF:
    ! NAME:
    ! ACTIVATE: Optional argument to activate the timer immediately
    !
    TYPE(TIMERTYPE), INTENT(IN OUT) :: SELF
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    LOGICAL, INTENT(IN), OPTIONAL :: ACTIVATE
    !
    !---------------------------------------------------------------------------
    !
    SELF%NAME = NAME
    !
    IF (PRESENT(ACTIVATE)) THEN
        !
        SELF%ACTIVE = ACTIVATE
        !
    END IF
    !
    END SUBROUTINE TIMERINIT
    !
    !===========================================================================
    !
    SUBROUTINE TIMERACTIVATE(TIMER, ACTIVATE)
    !
    !  Activate timer.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! TIMER:
    ! ACTIVATE:
    !
    TYPE(TIMERTYPE), INTENT(IN OUT) :: TIMER
    LOGICAL, INTENT(IN), OPTIONAL :: ACTIVATE
    !
    !---------------------------------------------------------------------------
    !
    IF (PRESENT(ACTIVATE)) THEN
        !
        TIMER%ACTIVE = ACTIVATE
        !
    ELSE
        !
        TIMER%ACTIVE = .TRUE.
        !
    END IF
    !
    END SUBROUTINE TIMERACTIVATE
    !
    !===========================================================================
    !
    SUBROUTINE TIMERSTART(TIMER)
    !
    ! Start the timer
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! TIMER:
    !
    TYPE(TIMERTYPE), INTENT(IN OUT) :: TIMER
    !
    !---------------------------------------------------------------------------
    !
    IF (TIMER%ACTIVE) THEN
        !
        CALL CPU_TIME(TIMER%TIME_BEG)
        !
        IF (TIMER%TIME_BEG < 0.0) THEN
            !
            TIMER%ACTIVE = .FALSE.
            TIMER%FAILED = .TRUE.
            !
        ELSE
            !
            TIMER%STARTS = TIMER%STARTS + 1
            !
        END IF
        !
    END IF
    !
    END SUBROUTINE TIMERSTART
    !
    !===========================================================================
    !
    SUBROUTINE TIMERSTOP(TIMER)
    !
    !  Stop timer and compute elapsed time.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! TIMER:
    !
    TYPE(TIMERTYPE), INTENT(IN OUT) :: TIMER
    !
    !---------------------------------------------------------------------------
    !
    IF (TIMER%ACTIVE) THEN
        !
        CALL CPU_TIME(TIMER%TIME_END)
        !
        TIMER%TIME_TOTAL = TIMER%TIME_TOTAL + TIMER%TIME_END - TIMER%TIME_BEG
        TIMER%STOPS = TIMER%STOPS + 1
        !
    END IF
    !
    END SUBROUTINE TIMERSTOP
    !
    !===========================================================================
    !
    SUBROUTINE TIMERWRITE(TIMER, UNIT, LEVEL)
    !
    !  Print statistics for each timer.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! TIMER:
    ! UNIT: Output unit
    ! LEVEL: (Optional) indentation level
    !
    TYPE(TIMERTYPE), INTENT(IN) :: TIMER
    INTEGER, INTENT(IN) :: UNIT
    INTEGER, INTENT(IN), OPTIONAL :: LEVEL
    !
    ! Locals:
    !
    CHARACTER(LEN=4) :: INDENT='    '
    INTEGER :: MYLEVEL
    CHARACTER(LEN=32) :: MYFMT
    !
    !---------------------------------------------------------------------------
    !
    WRITE(UNIT, '(a)') 'Statistics for timer:  '//TIMER%NAME
    !
    IF (TIMER%ACTIVE) THEN
        !
        WRITE(UNIT, '(a30,a)')  'timer is currently:  ', 'ACTIVE'
        !
    ELSE
        !
        WRITE(UNIT, '(a30,a)')  'timer is currently:  ', 'INACTIVE'
        !
    END IF
    !
    IF (TIMER%FAILED) THEN
        !
        WRITE(UNIT, '(a)')  INDENT // '*** timer FAILED'
        !
    END IF
    !
    WRITE(UNIT, '(a30,i0)') 'timer STARTS:  ', TIMER%STARTS
    WRITE(UNIT, '(a30,i0)') 'timer STOPS:  ', TIMER%STOPS
    WRITE(UNIT, '(a30,f0.4)') 'total time in seconds:  ', TIMER%TIME_TOTAL
    !
    IF (TIMER%STOPS > 0) THEN
        !
        WRITE(UNIT, '(a30,f0.4)') 'average time per call:  ', &
            & TIMER%TIME_TOTAL/REAL(TIMER%STOPS)
        !
    ELSE
        !
        WRITE(UNIT, '(a30,a)') 'average time per call:  ', 'NA'
        !
        RETURN
        !
    END IF
    !
    ! Summary line.
    !
    IF (PRESENT(LEVEL)) THEN
        !
        MYLEVEL = 2 * LEVEL
        !
    ELSE
        !
        MYLEVEL = 0
        !
    END IF
    !
    IF (MYLEVEL > 0) THEN
        !
        WRITE(MYFMT, '(a,i0,a)') '(a30,a,', MYLEVEL, 'x,f0.3,a,i0,a,f0.3)'
        !
    ELSE
        !
        MYFMT = '(a30,a,f0.3,a,i0,a,f0.3)'
        !
    END IF
    !
    WRITE(UNIT, MYFMT) 'SUMMARY:  ', TRIM(ADJUSTL(TIMER%NAME)) // '---', &
        & TIMER%TIME_TOTAL, '/', TIMER%STOPS, &
        & '=',TIMER%TIME_TOTAL/REAL(TIMER%STOPS)
    !
    END SUBROUTINE TIMERWRITE
    !
END MODULE TIMER_MOD
