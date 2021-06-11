! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE STRINGS_MOD
!
! String parameters and operations.
!
! Contains subroutines:
! GETLINE: Read line of input, ignoring comment cards, and echoing output if
!   desired. Remove any leading whitespace.
! GETWORD: Extract first word from a string.  Words are separated by the space
!   character.
! CHANGE_CHAR: Change all occurences of a particular character to another.
!
USE INTRINSIC_TYPES_MOD, ONLY: IK=>INTEGER_KIND
!
IMPLICIT NONE
!
! Private
!
PRIVATE
!
CHARACTER, PARAMETER :: ASCII_TAB = ACHAR(9)
CHARACTER, PARAMETER :: ASCII_SPACE = ACHAR(32)
!
! Public
!
PUBLIC :: STANDARD_LINE_LENGTH
PUBLIC :: GETLINE, GETWORD
!
INTEGER, PARAMETER :: STANDARD_LINE_LENGTH = 256
!
CONTAINS
    !
    SUBROUTINE GETLINE(UNIT, LINE, STATUS, ECHO, COMMENT)
    !
    ! Read line of input, ignoring comment cards, and echoing output if desired.
    !   Remove any leading whitespace.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! UNIT: Input unit
    ! LINE: A line of input
    ! STATUS: Return value;
    !       : <  0 --> EOF
    !       : == 0 --> success
    !       : >  0 --> other error reading file
    ! ECHO: (Optional) echo unit
    ! COMMENT: (Optional) comment character
    !
    INTEGER(IK), INTENT(IN) :: UNIT
    CHARACTER(LEN=*), INTENT(OUT) :: LINE
    INTEGER(IK), INTENT(OUT) :: STATUS
    INTEGER(IK), INTENT(IN), OPTIONAL :: ECHO
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: COMMENT
    !
    ! Locals:
    !
    INTEGER :: IOST
    INTEGER :: NCHAR
    LOGICAL :: ECHO_ON
    LOGICAL :: IS_OPEN
    LOGICAL :: COMMENTED
    CHARACTER(LEN=8) :: WRITABLE
    !
    !---------------------------------------------------------------------------
    !
    ! Make sure the echo unit is associated with a writable file. Otherwise,
    !   turn it off.
    !
    IF (PRESENT(ECHO)) THEN
        !
        IF (ECHO >= 0) THEN ! define IK_ZERO in constants_mod
            !
            INQUIRE(UNIT = ECHO, OPENED = IS_OPEN, WRITE = WRITABLE)
            !
            IF (IS_OPEN) THEN
                !
                IF (WRITABLE .EQ. 'YES') THEN
                    !
                    ECHO_ON = .TRUE.
                    !
                ELSE
                    !
                    ECHO_ON = .FALSE.
                    !
                END IF
                !
            ELSE
                !
                ECHO_ON = .FALSE.
                !
            END IF
            !
        END IF
        !
    ELSE
        !
        ECHO_ON = .FALSE.
        !
    END IF
    !
    ! Determine whether to check for comment characters.
    !
    IF (PRESENT(COMMENT)) THEN
        !
        COMMENTED = .TRUE.
        !
    ELSE
        !
        COMMENTED = .FALSE.
        !
    END IF
    !
    ! Read lines, ignoring blank lines and comments until we find an input line
    !   or end-of-file.
    !
    READ_LOOP: DO
        !
        !  Read line.
        !
        READ(UNIT, '(a)', IOSTAT = STATUS) LINE
        !
        IF (STATUS .NE. 0) THEN
            !
            EXIT READ_LOOP
            !
        END IF
        !
        !  Echo the line as requested.
        !
        IF (ECHO_ON) THEN
            !
            WRITE(ECHO, '(a)', IOSTAT = IOST) TRIM(LINE)
            !
            IF (IOST /= 0) THEN
                !
                ECHO_ON = .FALSE.
                !
            END IF
            !
        END IF
        !
        !  Check for blank lines and comment lines.
        !
        CALL CHANGE_CHAR(LINE, ASCII_TAB, ASCII_SPACE)
        !
        NCHAR = LEN_TRIM(LINE)
        !
        IF (NCHAR == 0) THEN
            !
            CYCLE READ_LOOP
            !
        END IF
        !
        LINE = ADJUSTL(LINE)
        !
        IF (COMMENTED) THEN
            !
            IF (LINE(1:1) == COMMENT) THEN
                !
                CYCLE READ_LOOP
                !
            END IF
            !
        END IF
        !
        !  If we are here, then we have a line of input.
        !
        STATUS = 0
        EXIT READ_LOOP
        !
    END DO READ_LOOP
    !
    END SUBROUTINE GETLINE
    !
    !===========================================================================
    !
    SUBROUTINE GETWORD(STRING, WORD, STATUS)
    !
    ! Extract first word from a string.  Words are separated by the space
    !   character.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! STRING: Input string; on output, string with first word removed
    ! WORD: First word of string
    ! STATUS: return value;
    !       :  0    --> success
    !       :  /= 0 --> failure/no more words
    !
    CHARACTER(LEN=*), INTENT(IN OUT) ::  STRING
    CHARACTER(LEN=*), INTENT(OUT) ::  WORD
    INTEGER, INTENT(OUT) :: STATUS
    !
    ! Locals:
    !
    INTEGER :: LSTR
    INTEGER :: LOC
    !
    !---------------------------------------------------------------------------
    !
    ! Check for blank string---no more words.
    !  
    LSTR = LEN_TRIM(STRING)
    !
    IF (LSTR == 0) THEN
        !
        STATUS = 1
        WORD = ASCII_SPACE
        !
        RETURN
        !
    END IF
    !
    ! Locate the first blank character, put the word in `word' and take it out
    !   of `string'.
    !
    STRING = ADJUSTL(STRING)
    LOC = SCAN(STRING, ASCII_SPACE) - 1
    WORD = STRING(1:LOC)
    STRING(1:LOC) = ASCII_SPACE
    STRING = ADJUSTL(STRING)
    STATUS = 0
    !
    END SUBROUTINE GETWORD
    !
    !===========================================================================
    !
    SUBROUTINE CHANGE_CHAR(STRING, FROM_CHAR, TO_CHAR)
    !
    ! Change all occurences of a particular character to another.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! STRING: Character string to change
    ! FROM_CHAR: Character to be changed
    ! TO_CHAR: New character
    !
    CHARACTER(LEN=*), INTENT(INOUT) :: STRING
    CHARACTER, INTENT(IN) :: FROM_CHAR
    CHARACTER, INTENT(IN) :: TO_CHAR
    !
    ! Locals:
    !
    INTEGER(IK) :: LEN_STRING
    INTEGER(IK) :: I
    !
    !---------------------------------------------------------------------------
    !
    LEN_STRING = LEN(string)
    !
    DO I = 1, LEN_STRING
        !
        IF (STRING(I:I) == FROM_CHAR) THEN
            !
            STRING(I:I) = TO_CHAR
            !
        END IF
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE CHANGE_CHAR
    !
END MODULE STRINGS_MOD
