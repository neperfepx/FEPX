! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE F90UtilitiesModule 
!
!  *** Program Unit:  Module
!  ***    Unit Name:  F90UtilitiesModule
!  ***  Description:
!
!  This module contains additions to the basic
!  f90 intrinsic functions.
!
!  *** Use Statements:
!
USE IntrinsicTypesModule, &
     &  RK=>REAL_KIND, IK=>INTEGER_KIND, LK=>LOGICAL_KIND
!
!  *** End:
!
IMPLICIT NONE
!
PRIVATE
!
!  *** Public Procedures:
!
PUBLIC CArrayToString, StringToCArray
!
!  *** End:
!
!--------------*-------------------------------------------------------
!
CONTAINS 
!
!  *** Program Unit:  function
!  ***    Unit Name:  CArrayToString
!
!  *** Unit Declaration: 
!
FUNCTION CArrayToString(carray) RESULT(Result)
  !
  !  ***  Description:  
  !
  !  Return string from Character array.
  !
  !  *** Argument Declarations:
  !
  !  carray -- an array of characters
  !
  CHARACTER(LEN=1), INTENT(IN) :: carray(:)
  !
  !  *** Result:
  !
  !  Result -- a string with the same characters
  !  
  CHARACTER(LEN=SIZE(carray, DIM=1)), POINTER:: Result
  !
  !  *** Locals:
  !
  INTEGER :: i
  !
  !  *** End:

  !--------------*-------------------------------------------------------
  !
  ALLOCATE(Result)
  do i=1, SIZE(carray, DIM=1)
     Result(i:i) = carray(i)
  end do
  !
  RETURN
  !
END FUNCTION CArrayToString
!
!  *** Program Unit:  function
!  ***    Unit Name:  StringToCArray
!
!  *** Unit Declaration: 
!
FUNCTION StringToCArray(string) RESULT(Result)
  !
  !  ***  Description:  
  !
  !  Convert a string to an array of characters.
  !
  !  *** Argument Declarations:
  !
  CHARACTER(LEN=*), INTENT(IN) :: string
  !
  !  *** Result:
  !  
  CHARACTER(LEN=1), POINTER :: Result(:)
  !
  !  *** Locals:
  !
  INTEGER :: i
  !
  !  *** Intrinsics:
  !
  INTRINSIC :: LEN
  !
  !  *** End:
  !
  !--------------*-------------------------------------------------------
  !
  ALLOCATE(Result(LEN(string)))
  do i=1, LEN(string)
     Result(i) = string(i:i)
  end do
  !
END FUNCTION 
END MODULE F90UtilitiesModule 
!
