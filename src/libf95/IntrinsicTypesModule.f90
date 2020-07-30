! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE IntrinsicTypesModule
  !
  !  ***  Description:
  !
  !  KIND specifications for intrinsic types.
  !
  !  INTEGERS:
  !  --------
  !
  !  Typically, there are 1-byte, 2-byte, 4-byte and 8-byte
  !  integers with 4-byte being the default and having 9 digits
  !  of precision.  Here, we only use short and long (4-byte
  !  and 8-byte).
  !
  !  Public Parameters
  !  -----------------
  !  
  !  INTEGER_KIND_S : short integer kind (at least 9 digits)
  !  INTEGER_KIND_L : long integer kind (at least 15 digits)
  !  
  !  INTEGER_KIND   : default integer kind for this module (short)
  !
  !  *** End:
  !
  INTEGER, PUBLIC, PARAMETER :: &
       &  INTEGER_KIND_S = SELECTED_INT_KIND(9),&
       &  INTEGER_KIND_L = SELECTED_INT_KIND(15)
  !
  INTEGER, PUBLIC, PARAMETER :: &
       &  INTEGER_KIND   = INTEGER_KIND_S
  !
  !  ***  Description:
  !
  !  REALS:
  !  -----
  !  
  !  The real kinds are specified using the PRECISION and RANGE
  !  arguments to SELECTED_REAL_KIND; on seahag (alpha/linux),
  !  the preicison is 6 for single and 15 for double, and the
  !  ranges are 37 for single and 307 for double.  Quadruple
  !  precision is also available on many compilers.
  !
  !  Public Parameters
  !  -----------------
  !
  !  REAL_KIND_S : single precision 
  !  REAL_KIND_D : double precision 
  !  REAL_KIND_Q : quadruple precision
  !  
  !  REAL_KIND   : default precision for this module (double)
  !
  !  *** End:
  !
  INTEGER, PUBLIC, PARAMETER :: &
       &  REAL_KIND_S = SELECTED_REAL_KIND( 5,  30),&
       &  REAL_KIND_D = SELECTED_REAL_KIND(14, 200),&
       &  REAL_KIND_Q = SELECTED_REAL_KIND(25, 200)
  !
  INTEGER, PUBLIC, PARAMETER :: REAL_KIND = REAL_KIND_D
  !
  !  ***  Description:
  !
  !  LOGICALS:
  !  --------
  !
  !  Usually, many logical kinds are available, 1-byte, 2-byte,
  !  4-byte and 8-byte; default is usually 4-byte.  There is no
  !  fortran function SELECTED_LOGICAL_KIND, so here we just
  !  stick with the system default logical kind.
  !
  !  Public Parameters:
  !  -----------------
  !
  !  LOGICAL_KIND_D : system default logical kind (KIND(.TRUE.))
  !
  !  LOGICAL_KIND   : default logical kind for this module
  !
  !  *** End:
  !
  INTEGER, PUBLIC, PARAMETER :: LOGICAL_KIND_D = KIND(.TRUE.)
  !
  INTEGER, PUBLIC, PARAMETER :: LOGICAL_KIND = LOGICAL_KIND_D
  !
  !  Other functions available regarding real numeric types are:
  !
  !  EPSILON, HUGE, MAXEXPONENT, MINEXPONENT, PRECISION,
  !  RADIX, RANGE and TINY
  !
  !
END MODULE IntrinsicTypesModule
