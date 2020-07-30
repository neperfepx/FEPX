! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE StringsModule
  !
  !  String parameters and operations.
  !
  !  PUBLIC ENTITIES:
  !  
  !  STANDARD_LINE_LENGTH
  !  getLine(unit, line, status, ECHO, COMMENT)
  !  getWord(string, word, status)
  !
  !  PRIVATE ENTITIES:
  !  change_char()  ! may want to make public if use arises
  !
  USE IntrinsicTypesModule, IK => INTEGER_KIND
  USE ConstantsModule, ONLY:  ASCII_TAB, ASCII_SPACE

  IMPLICIT NONE
  !
  ! ========== Public Entities
  !
  PRIVATE   ! all objects are private unless declared otherwise

  PUBLIC :: STANDARD_LINE_LENGTH
  PUBLIC :: getLine, getWord
  !
  ! ========== Data
  !
  !  Standard character size for reading a line from a text file.
  !
  INTEGER, PARAMETER :: STANDARD_LINE_LENGTH = 256
  !
CONTAINS
  !
  SUBROUTINE getLine(unit, line, status, ECHO, COMMENT)
    !
    ! Read line of input, ignoring comment cards, and echoing output
    ! if desired.  Remove any leading whitespace.
    !
    !  * Arguments:
    ! 
    !  unit    : input unit
    !  line    : a line of input
    !  status  : return value;
    !          : <  0 --> EOF
    !          : == 0 --> success
    !          : >  0 --> other error reading file
    !  ECHO    : (optional) echo unit
    !  COMMENT : (optional) comment character
    !
    CHARACTER(LEN=*), INTENT(OUT) :: line
    CHARACTER(LEN=1), INTENT(IN)  :: COMMENT
    !
    INTEGER(IK), INTENT(IN)       :: unit, ECHO
    INTEGER(IK), INTENT(OUT)      :: status
    !
    OPTIONAL :: ECHO, COMMENT
    !
    !  * Locals.
    !
    INTEGER :: iost, nchar
    LOGICAL :: echo_on, is_open, commented
    !
    CHARACTER(LEN=8) :: writable
    !
    !----------------------------------------------------------------------
    !
    !  Make sure the echo unit is associated with a writable file.
    !  Otherwise, turn it off.
    !
    if (PRESENT(ECHO)) then
      if (ECHO >= 0) then    !  define IK_ZERO in constants_mod
        !
        INQUIRE(UNIT=ECHO, OPENED=is_open, WRITE=writable)
        if (is_open) then
          if (writable .eq. 'YES') then
            echo_on = .TRUE.
          else
            echo_on = .FALSE.
          endif
        else
          echo_on = .FALSE.
        endif
        !
      endif
    else
      echo_on = .FALSE.
    endif
    !
    !  Determine whether to check for comment characters.
    !
    if (PRESENT(COMMENT)) then
      commented = .TRUE.
    else
      commented = .FALSE.
    endif
    !
    !  Read lines, ignoring blank lines and comments until we find
    !  an input line or end-of-file.
    !
    READ_LOOP: do
      !
      !  Read line.
      !
      read(unit, '(a)', iostat=status) line
      if (status .ne. 0) then
        EXIT READ_LOOP
      endif
      !
      !  Echo the line as requested.
      !
      if (echo_on) then
        write(ECHO, '(a)', IOSTAT=iost) TRIM(line)
        if (iost /= 0) then
          echo_on = .FALSE.
        endif
      endif
      !
      !  Check for blank lines and comment lines. 
      !
      call change_char(line, ASCII_TAB, ASCII_SPACE)
      nchar = LEN_TRIM(line)
      if (nchar == 0) then
        CYCLE READ_LOOP
      endif
      !
      line = ADJUSTL(line)
      if (commented) then
        if (line(1:1) == COMMENT) then
          CYCLE READ_LOOP
        endif
      endif
      !
      !  If we are here, then we have a line of input.
      !
      status = 0
      EXIT READ_LOOP
      !
    enddo READ_LOOP
    !
  END SUBROUTINE getLine
  !
  !
  SUBROUTINE getWord(string, word, status)
    !
    !  Extract first word from a string.  Words are separated by
    !  the space character.
    !
    !  * Arguments.
    !
    !  string    : input string; on output, string with first word removed
    !  word      : first word of string
    !  status    : return value;
    !            :  0    --> success
    !            :  /= 0 --> failure/no more words
    !
    CHARACTER(LEN=*), INTENT(IN OUT) ::  string
    CHARACTER(LEN=*), INTENT   (OUT) ::  word
    !
    INTEGER, INTENT(OUT) :: status
    !
    ! * Locals.
    !
    INTEGER lstr, loc
    !
    !----------------------------------------------------------------------
    !
    !  Check for blank string---no more words.
    !  
    lstr = LEN_TRIM(string)
    if (lstr == 0) then
      status = 1
      word   = ASCII_SPACE
      RETURN
    endif
    !
    !  Locate the first blank character, put the word in `word' and
    !  take it out of `string'.
    !
    string = ADJUSTL(string)
    loc = SCAN(string, ASCII_SPACE) - 1
    word = string(1:loc)
    string(1:loc) = ASCII_SPACE
    string = ADJUSTL(string)
    status = 0

  END SUBROUTINE getWord
  !
  !
  SUBROUTINE change_char(string, from_char, to_char)
    !
    !  Change all occurences of a particular character to another.
    !
    IMPLICIT NONE
    !
    !  Arguments.
    !
    !  string    : character string to change
    !  from_char : character to be changed
    !  to_char   : new character
    !
    CHARACTER(LEN=*), INTENT(INOUT) :: string
    CHARACTER, INTENT(IN)    :: from_char, to_char
    !
    !  Locals.
    !
    INTEGER(IK) :: len_string, i
    !
    !----------------------------------------------------------------------
    !
    len_string = LEN(string)
    !
    do i=1, len_string
      if (string(i:i) == from_char) then
        string(i:i) = to_char
      endif
    enddo
    !
    RETURN
    !
  END SUBROUTINE change_char
  !
END MODULE StringsModule
