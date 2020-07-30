! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE FilesModule 
  !
  !  Short module to handle fortran unit numbers.  It provides
  !  parameters STANDARD_INPUT, STANDARD_OUTPUT and STANDARD_ERROR,
  !  and one function NewUnitNumber() which returns a free FORTRAN
  !  unit number.
  !  
  !  *** Use Statements:
  !
  USE IntrinsicTypesModule, &
       &  IK=>INTEGER_KIND, LK=>LOGICAL_KIND
  !
  !  *** End:
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  !  Standard unit numbers.
  !
  INTEGER(IK), PARAMETER ::&
       STANDARD_INPUT  = 5,    &
       STANDARD_OUTPUT = 6,    &
       STANDARD_ERROR  = 0
  !
  INTEGER(IK) :: MAX_CHECK=1000, FIRST_UNIT=10
  INTEGER(IK) :: SET_NUMBER_TO_CHECK=1, SET_FIRST_UNIT=2  ! to be used for Set command
  !
  !--------------* Public Data:
  !
  PUBLIC :: STANDARD_INPUT, STANDARD_OUTPUT, STANDARD_ERROR, &
       &    SET_NUMBER_TO_CHECK, SET_FIRST_UNIT
  !
  PUBLIC :: NewUnitNumber, GetUnitNumber
  !
CONTAINS 
  !
  FUNCTION NewUnitNumber() RESULT(unit)
    !
    ! Return available unit number.
    !
    !  Result:
    !  
    INTEGER(IK) :: unit
    !
    ! unit : -1 --> failed to find new unit
    !      : >0 --> an available fortran unit number
    !
    !  Locals:
    !
    INTEGER(IK), SAVE :: unit_number=-1
    INTEGER(IK)       :: i, ios
    !
    LOGICAL(LK) :: is_open
    !
    !  *** End:
    !
    !----------------------------------------------------------------------
    !
    if (unit_number < 0) then
       unit_number = FIRST_UNIT
    end if
    !
    unit = -1
    !
    do i=1, MAX_CHECK
       !
       INQUIRE(UNIT=unit_number, OPENED=is_open, IOSTAT=ios)
       !
       if (is_open) then
          unit_number = unit_number + 1
       else
          unit        = unit_number
          unit_number = unit_number + 1
          RETURN
       endif
       !
    enddo
    !
    RETURN
    !
  END FUNCTION NewUnitNumber

  !
  FUNCTION  GetUnitNumber(file, ACTION, FORM)
    !
    !  Return unit number of file, if opened, or open a new file
    !  with appropriate options.
    !
    !--------------*-------------------------------------------------------
    !
    !  Arguments.
    !
    !  file   -- name of file (possibly open)
    !  ACTION -- (optional) defaults to 'READWRITE'
    !  FORM   -- (optional) defaults to 'FORMATTED'
    !
    INTEGER(IK) :: GetUnitNumber  ! result
    !
    CHARACTER(LEN=*) :: file, ACTION, FORM
    !
    INTENT(IN   ) :: file, ACTION, FORM
    !
    OPTIONAL :: ACTION, FORM
    !
    !  Locals.
    !
    !LOGICAL :: isopen
    INTEGER(IK) :: unit, stat
    CHARACTER(LEN=16) :: myform, myaction, answer
    !
    !--------------*-------------------------------------------------------
    !
    GetUnitNumber = -1
    !
    INQUIRE(FILE=file, NUMBER=unit)
    !
    if (unit == -1) then
       !  
       !  Need to open file.
       !
       if (PRESENT(ACTION) ) then
         IF ( (ACTION == 'READ') .OR. (ACTION == 'WRITE') .OR. (ACTION == 'READWRITE') ) THEN
           myaction = ACTION
         ELSE
           RETURN ! failure
         END IF
       else
          myaction = 'READWRITE'
       end if
       !
       if (PRESENT(FORM)) then
          if (FORM == 'FORMATTED') then
             myform = 'FORMATTED'
          else if (FORM == 'UNFORMATTED') then
             myform = 'UNFORMATTED'
          else
             RETURN ! failure
          end if
       else 
          myform = 'FORMATTED'
       end if
       !
       unit = NewUnitNumber()
       OPEN(UNIT=unit, FILE=file, ACTION=myaction, FORM=myform, IOSTAT=stat)
       if (stat == 0) then
          GetUnitNumber = unit
       endif
       !
    else
       !
       !  File already open, make sure options are okay.
       !
       if (PRESENT(ACTION) ) then
          if (ACTION == 'READ') then
             INQUIRE(UNIT=unit,  READ=answer)
          else if (ACTION == 'WRITE') then
             INQUIRE(UNIT=unit,  WRITE=answer)
          else if (ACTION == 'READWRITE') then
             INQUIRE(UNIT=unit,  READWRITE=answer)
          else
             RETURN 
          end if
       else
          INQUIRE(UNIT=unit,  READWRITE=answer) ! check this
       end if
       !
       if (answer /= 'YES') then
          RETURN 
       end if
       !
       if (PRESENT(FORM)) then
          if (FORM == 'FORMATTED') then
             INQUIRE(UNIT=unit,  FORMATTED=answer)           
          else if (FORM == 'UNFORMATTED') then
             INQUIRE(UNIT=unit,  UNFORMATTED=answer)           
          else
             answer = 'NO'
          end if
       else 
          INQUIRE(UNIT=unit,  FORMATTED=answer)
       end if
       !
       if (answer /= 'YES') then
          RETURN 
       end if
       !
       GetUnitNumber = unit
       !
    end if
    !
  END FUNCTION GetUnitNumber
END MODULE FilesModule
