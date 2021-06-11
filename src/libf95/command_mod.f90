! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE COMMAND_MOD
!
! Tools for processing command strings.
!
! Contains subroutines:
! COMMANDLOOP: Process command strings
!
! Contains functions
! EXECCOMMAND: Call a callback function if a command line matches a command
!
USE FILES_MOD, ONLY:  STANDARD_OUTPUT, STANDARD_ERROR
USE STRINGS_MOD, ONLY: GETWORD, GETLINE, STANDARD_LINE_LENGTH
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: LOOPSTAT_EOF = -1
INTEGER, PARAMETER :: LOOPSTAT_END = 0
INTEGER, PARAMETER :: LOOPSTAT_NZ = 1
INTEGER, PARAMETER :: LOOPSTAT_NOMATCH = -2
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: EXECCOMMAND
PUBLIC :: COMMANDLOOP
PUBLIC :: LOOPSTAT_EOF
PUBLIC :: LOOPSTAT_END
PUBLIC :: LOOPSTAT_NZ
PUBLIC :: LOOPSTAT_NOMATCH
!
CONTAINS
    !
    FUNCTION EXECCOMMAND(COMMAND, CMDLINE, CALLBACK, STATUS) RESULT(MATCH)
    !
    ! Call a callback function if a command line matches a command.
    !
    ! If :arg:`CMDLINE` matches the command string :arg:`command`, then
    ! call the subroutine :arg:`callback`, which sets the value of :arg:`status`
    ! Otherwise, set :arg:`match` to false.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! COMMAND: The command being checked
    ! CMDLINE: String with command and arguments
    ! CALLBACK: Subroutine to call if command matches
    ! STATUS: Return status of callback, if executed
    ! MATCH: (Output) true if CMDLINE matches command
    !
    CHARACTER(LEN=*), INTENT(IN) :: COMMAND
    CHARACTER(LEN=*), INTENT(IN) :: CMDLINE
    !
    INTERFACE
        !
        SUBROUTINE CALLBACK(A, S)
        !
        CHARACTER(LEN=*), INTENT(IN) :: A
        INTEGER, INTENT(OUT) :: S
        !
        END SUBROUTINE CALLBACK
        !
    END INTERFACE
    !
    INTEGER, INTENT(OUT) :: STATUS
    LOGICAL :: MATCH
    !
    ! Locals:
    !
    CHARACTER(LEN=LEN(CMDLINE)) :: WORD
    CHARACTER(LEN=LEN(CMDLINE)) :: ARGS
    INTEGER :: ISTAT
    !
    !---------------------------------------------------------------------------
    !
    STATUS = 0
    ARGS = CMDLINE
    !
    CALL GETWORD(ARGS, WORD, ISTAT)
    !
    IF (WORD == COMMAND) THEN
        !
        MATCH = .TRUE.
        CALL CALLBACK(ARGS, STATUS)
        !
    ELSE
        !
        MATCH = .FALSE.
        !
    END IF
    !
    END FUNCTION EXECCOMMAND
    !
    !===========================================================================
    !
    SUBROUTINE COMMANDLOOP(INUNIT, KWPROCESS, MYSTATUS)
    !
    ! Process command strings.
    !
    ! This subroutine repeatedly reads a line of input and calls the keyword
    !   processing subroutine. The subroutine returns if any of these conditions
    !   are met:
    !   - the input units reaches EOF
    !   - a nonzero status is returned from the keyword processing command
    !   - any command line is not matched
    !
    !---------------------------------------------------------------------------
    !
    INTEGER, INTENT(IN) :: INUNIT
    !
    INTERFACE
        !
        LOGICAL FUNCTION KWPROCESS(CMDLINE, INUNIT, STATUS)
        !
        CHARACTER(LEN=*), INTENT(IN) :: CMDLINE
        INTEGER, INTENT(IN) :: INUNIT
        INTEGER, INTENT(OUT) :: STATUS
        !
        END FUNCTION KWPROCESS
        !
    END INTERFACE
    !
    INTEGER, INTENT(OUT) :: MYSTATUS
    !
    ! Locals:
    !
    CHARACTER(LEN=*), PARAMETER :: MYNAME = 'COMMANDLOOP'
    CHARACTER(LEN=STANDARD_LINE_LENGTH) :: CMDLINE
    INTEGER :: STATUS
    LOGICAL :: CMDMATCHED
    !
    !---------------------------------------------------------------------------
    !
    MYSTATUS = 0
    !
    DO
        !
        CALL GETLINE(INUNIT, CMDLINE, STATUS, COMMENT='#')
        !
        IF (STATUS /= 0) THEN
            !
            ! Used to set MYSTATUS to LOOPSTAT_EOF, but since we opt to remove
            !   the EOF string the only place status can be non-zero is at EOF
            !
            MYSTATUS = LOOPSTAT_END; EXIT
            !
        END IF
        !
        ! Now call subroutine that checks for matches.
        !
        CMDMATCHED = KWPROCESS(CMDLINE, INUNIT, STATUS)
        !
        IF (CMDMATCHED) THEN
            !
            ! Matched: check status
            !
            IF (STATUS /= 0) THEN
                !
                CALL WRITEMESSAGE(STANDARD_ERROR, 'Error  :     > No match on &
                    &command: '// CMDLINE)
                !
                MYSTATUS = LOOPSTAT_NZ; EXIT
            ELSE
                !
                !CALL WRITEMESSAGE(STANDARD_OUTPUT, 'processed command:  ' //&
                !   &CMDLINE)
                ! Above CALL suppressed in favor of displaying commands
                !   elsewhere - JC
            CYCLE
                !
            END IF
            !
        ELSE
            !
            ! No match: exit
            !
            CALL WRITEMESSAGE(STANDARD_ERROR, 'Error  :     > No match on &
                &command: '// CMDLINE)
            MYSTATUS = LOOPSTAT_NOMATCH; EXIT
            !
        END IF
        !
    END DO
    !
    CONTAINS
    !
    SUBROUTINE WRITEMESSAGE(U, S)
        !
        ! Write message S to unit U
        !
        INTEGER, INTENT(IN) :: U
        CHARACTER(LEN=*), INTENT(IN) :: S
        !
        WRITE(U, '(a)') TRIM(S)
        !
    END SUBROUTINE WRITEMESSAGE
    !
    END SUBROUTINE COMMANDLOOP
    !
END MODULE COMMAND_MOD
