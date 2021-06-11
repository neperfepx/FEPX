! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE FILES_MOD
!
! Short module to handle fortran unit numbers. It provides parameters
!   STANDARD_INPUT, STANDARD_OUTPUT and STANDARD_ERROR, and one function
!   NEWUNITNUMBER() which returns a free FORTRAN unit number.
!
! Contains functions:
! NEWUNITNUMBER: Return available unit number
! GETUNITNUMBER: Return unit number of file, if opened, or open a new file with
!   appropriate options.
!
USE INTRINSIC_TYPES_MOD, ONLY: IK=>INTEGER_KIND, LK=>LOGICAL_KIND
!
IMPLICIT NONE
!
! Private
!
PRIVATE
!
! Standard unit numbers:
!
INTEGER(IK), PARAMETER :: STANDARD_INPUT = 5
INTEGER(IK), PARAMETER :: STANDARD_OUTPUT = 6
INTEGER(IK), PARAMETER :: STANDARD_ERROR = 0
!
INTEGER(IK) :: MAX_CHECK=1000
INTEGER(IK) :: FIRST_UNIT=10
INTEGER(IK) :: SET_NUMBER_TO_CHECK = 1
INTEGER(IK) :: SET_FIRST_UNIT = 2  ! to be used for Set command
!
! Public:
!
PUBLIC :: STANDARD_INPUT
PUBLIC :: STANDARD_OUTPUT
PUBLIC :: STANDARD_ERROR
PUBLIC :: SET_NUMBER_TO_CHECK
PUBLIC :: SET_FIRST_UNIT
!
PUBLIC :: NEWUNITNUMBER
PUBLIC :: GETUNITNUMBER
!
CONTAINS 
    !
    FUNCTION NEWUNITNUMBER() RESULT(UNIT)
    !
    ! Return available unit number
    !
    !---------------------------------------------------------------------------
    !
    ! ARGUMENTS:
    ! UNIT: -1 (failed to find new unit) >0 (available Fortran unit number)
    !  
    INTEGER(IK) :: UNIT
    !
    ! Locals:
    !
    INTEGER(IK), SAVE :: UNIT_NUMBER = -1
    INTEGER(IK) :: I
    INTEGER(IK) :: IOS
    LOGICAL(LK) :: IS_OPEN
    !
    !---------------------------------------------------------------------------
    !
    IF (UNIT_NUMBER < 0) THEN
        !
        UNIT_NUMBER = FIRST_UNIT
        !
    END IF
    !
    UNIT = -1
    !
    DO I = 1, MAX_CHECK
        !
        INQUIRE(UNIT = UNIT_NUMBER, OPENED = IS_OPEN, IOSTAT = IOS)
        !
        IF (IS_OPEN) THEN
            !
            UNIT_NUMBER = UNIT_NUMBER + 1
            !
        ELSE
            !
            UNIT = UNIT_NUMBER
            UNIT_NUMBER = UNIT_NUMBER + 1
            !
            RETURN
            !
        END IF
        !
    END DO
    !
    RETURN
    !
    END FUNCTION NEWUNITNUMBER
    !
    !===========================================================================
    !
    FUNCTION GETUNITNUMBER(FILE, ACTION, FORM)
    !
    ! Return unit number of file, if opened, or open a new file with
    !   appropriate options.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! FILE: Name of file (possibly open)
    ! ACTION: (optional) defaults to 'READWRITE'
    ! FORM:(optional) defaults to 'FORMATTED'
    ! GETUNITNUMBER: Result
    !
    CHARACTER(LEN=*), INTENT(IN) :: FILE
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: ACTION
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: FORM
    INTEGER(IK) :: GETUNITNUMBER
    !
    ! Locals:
    !
    !LOGICAL :: isopen
    INTEGER(IK) :: UNIT, STAT
    CHARACTER(LEN=16) :: MYFORM
    CHARACTER(LEN=16) :: MYACTION
    CHARACTER(LEN=16) :: ANSWER
    !
    !---------------------------------------------------------------------------
    !
    GETUNITNUMBER = -1
    !
    INQUIRE(FILE = FILE, NUMBER = UNIT)
    !
    IF (UNIT == -1) THEN
        !
        ! Need to open file
        !
        IF (PRESENT(ACTION)) THEN
            !
            IF ( (ACTION == 'READ') .OR. (ACTION == 'WRITE') .OR. (ACTION == &
                & 'READWRITE') ) THEN
                !
                MYACTION = ACTION
                !
            ELSE
                !
                RETURN ! failure
                !
            END IF
            !
        ELSE
            !
            MYACTION = 'READWRITE'
            !
        END IF
        !
        IF (PRESENT(FORM)) THEN
            !
            IF (FORM == 'FORMATTED') THEN
                !
                MYFORM = 'FORMATTED'
                !
            ELSEIF (FORM == 'UNFORMATTED') THEN
                !
                MYFORM = 'UNFORMATTED'
                !
            ELSE
                !
                RETURN ! failure
                !
            END IF
        ELSE
            !
            MYFORM = 'FORMATTED'
            !
        END IF
        !
        UNIT = NEWUNITNUMBER()
        OPEN(UNIT = UNIT, FILE = FILE, ACTION = MYACTION, FORM = MYFORM, &
            & IOSTAT = STAT)
        !
        IF (STAT == 0) THEN
            !
            GETUNITNUMBER = UNIT
            !
        END IF
       !
    ELSE
        !
        ! File already open, make sure options are okay
        !
        IF (PRESENT(ACTION) ) THEN
            !
            IF (ACTION == 'READ') THEN
                !
                INQUIRE(UNIT = UNIT,  READ = ANSWER)
                !
            ELSEIF (ACTION == 'WRITE') THEN
                !
                INQUIRE(UNIT = UNIT,  WRITE = ANSWER)
                !
            ELSEIF (ACTION == 'READWRITE') THEN
                !
                INQUIRE(UNIT = UNIT,  READWRITE = ANSWER)
                !
            ELSE
                !
                RETURN
                !
            END IF
            !
        ELSE
            !
            INQUIRE(UNIT = UNIT,  READWRITE = ANSWER) ! check this
            !
        END IF
        !
        IF (ANSWER /= 'YES') THEN
            !
            RETURN
        END IF
        !
        IF (PRESENT(FORM)) THEN
            !
            IF (FORM == 'FORMATTED') THEN
                !
                INQUIRE(UNIT = UNIT,  FORMATTED = ANSWER)
                !
            ELSEIF (FORM == 'UNFORMATTED') THEN
                !
                INQUIRE(UNIT = UNIT,  UNFORMATTED = ANSWER)
                !
            ELSE
                !
                ANSWER = 'NO'
                !
            END IF
            !
        ELSE
            !
            INQUIRE(UNIT = UNIT,  FORMATTED = ANSWER)
            !
        END IF
        !
        IF (ANSWER /= 'YES') THEN
            !
            RETURN
            !
        END IF
        !
        GETUNITNUMBER = UNIT
        !
    END IF
    !
    END FUNCTION GETUNITNUMBER
    !
END MODULE FILES_MOD
