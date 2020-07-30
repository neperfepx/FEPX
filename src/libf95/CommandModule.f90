! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE CommandModule
  !
  !  Tools for processing command strings.
  !
  USE StringsModule, ONLY: getWord, getLine, STANDARD_LINE_LENGTH
  USE FilesModule, ONLY:  STANDARD_OUTPUT, STANDARD_ERROR

  IMPLICIT NONE
  !
  ! ========== Data
  !
  ! Status values for CommandLoop
  !
  INTEGER, PARAMETER :: LOOPSTAT_EOF = -1, LOOPSTAT_END = 0, LOOPSTAT_NZ = 1, &
       &   LOOPSTAT_NOMATCH = -2
  !
  ! ========== Public Entities
  !
  PRIVATE   ! all objects are private unless declared otherwise
  !
  PUBLIC :: ExecCommand, CommandLoop
  PUBLIC :: LOOPSTAT_EOF, LOOPSTAT_END, LOOPSTAT_NZ, LOOPSTAT_NOMATCH
  !
CONTAINS
  !
  ! ==================================================== BEGIN:  ExecCommand
  !
  FUNCTION ExecCommand(command, cmdLine, callback, status) RESULT(Match)
    !
    ! Call a callback function if a command line matches a command.
    !
    ! If :arg:`cmdLine` matches the command string :arg:`command`, then
    ! call the subroutine :arg:`callback`, which sets the value of :arg:`status`.  
    ! Otherwise, set :arg:`match` to false.
    ! 
    CHARACTER(LEN=*), INTENT(IN) :: command, cmdLine
    !
    INTERFACE
      SUBROUTINE callback(a, s)
        CHARACTER(LEN=*), INTENT(IN)  :: a
        INTEGER,          INTENT(OUT) :: s
      END SUBROUTINE callback
    END INTERFACE
    !
    INTEGER, INTENT(OUT) :: status
    !  
    !  Result.
    !  
    ! command  -- the command being checked
    ! cmdLine  -- string with command and arguments
    ! callback -- subroutine to call if command matches
    ! status   -- return status of callback, if executed
    !
    LOGICAL :: Match
    !
    ! match -- true if cmdline matches command
    !
    ! ========== Locals
    !
    CHARACTER(LEN=LEN(cmdLine)) :: word, args
    INTEGER :: iStat
    !
    ! ================================================== Executable Code
    !
    status = 0
    args = cmdLine
    !
    CALL getWord(args, word, iStat)
    IF (word == command) THEN
      Match = .TRUE.
      CALL callback(args, status)
    ELSE
      Match = .FALSE.
    END IF
    !
  END FUNCTION ExecCommand
  !
  ! ====================================================   END:  ExecCommand
  ! ==================================================== BEGIN:  CommandLoop
  !
  SUBROUTINE CommandLoop(inUnit, kwProcess, myStatus)
    !
    !  Process command strings.
    !
    !  This subroutine repeatedly reads a line of input and calls the
    !  the keyword processing subroutine.  The subroutine returns if
    !  any of these conditions are met:
    !  * the input units reaches EOF
    !  * a nonzero status is returned from the keyword processing command
    !  * any command line is not matched
    !
    INTEGER, INTENT(IN) :: inUnit
    !
    INTERFACE
      LOGICAL FUNCTION kwProcess(cmdLine, inUnit, status)
        CHARACTER(LEN=*), INTENT(IN)  :: cmdLine
        INTEGER,          INTENT(IN)  :: inUnit
        INTEGER,          INTENT(OUT) :: status
      END FUNCTION kwProcess
    END INTERFACE
    !
    INTEGER, INTENT(OUT) :: myStatus
    !
    !  inUnit -- input unit to read lines from
    !
    ! ========== Locals
    !
    CHARACTER(LEN=*), PARAMETER :: myNAME = 'CommandLoop'
    !
    CHARACTER(LEN=STANDARD_LINE_LENGTH) :: cmdLine
    !
    INTEGER :: status
    LOGICAL :: hitEnd, cmdMatched
    !
    ! ================================================== Executable Code
    !
    myStatus = 0
    DO 
      CALL getLine(inUnit, cmdLine, status, COMMENT='#')
      IF (status /= 0) THEN
        ! Used to set mystatus to LOOPSTAT_EOF, but since we opt to remove the
        ! EOF string the only place status can be non-zero is at EOF - JC
        mystatus = LOOPSTAT_END; EXIT
      END IF
      !
      ! Removed in favor of using EOF to end input - JC
!      IF ((cmdLine == 'end') .OR. (cmdLine == 'end-input')) THEN
!        mystatus = LOOPSTAT_END; EXIT
!      END IF
      !
      !  Now call subroutine that checks for matches.
      !
      cmdMatched = kwProcess(cmdLine, inUnit, status)
      !
      IF (cmdMatched) THEN
        !
        !  Matched:  check status
        !
        IF (status /= 0) THEN
          CALL writeMessage(STANDARD_ERROR, '*** Error executing command:' // cmdLine)
          mystatus = LOOPSTAT_NZ; EXIT
        ELSE
          !CALL writeMessage(STANDARD_OUTPUT, 'processed command:  ' // cmdLine)
          ! Above CALL suppressed in favor of displaying commands elsewhere - JC
          CYCLE
        END IF
        !
      ELSE
        !
        !  No match:  exit
        !
        CALL writeMessage(STANDARD_ERROR, '*** No match on command:' // cmdLine)
        mystatus = LOOPSTAT_NOMATCH; EXIT
        !
      END IF
      !
    END DO
    !
  CONTAINS
    !
    SUBROUTINE writeMessage(u, s)
      ! write message s to unit u
      INTEGER, INTENT(IN) :: u
      CHARACTER(LEN=*), INTENT(IN) :: s
      !
      WRITE(u, '(a)') TRIM(s)
      !q
    END SUBROUTINE writeMessage
    !
  END SUBROUTINE CommandLoop
  !
  ! ====================================================   END:  CommandLoop
  !
END MODULE CommandModule
