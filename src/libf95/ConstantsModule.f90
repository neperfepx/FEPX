! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE ConstantsModule
  !
  !  This module provides commonly used constants of all types.
  !
  ! ==================== Other Modules
  !
  USE IntrinsicTypesModule, ONLY : &
       &  RK => REAL_KIND, IK => INTEGER_KIND, LK => LOGICAL_KIND
  !
  IMPLICIT NONE
  !
  PRIVATE RK, IK, LK
  !
  ! ==================== Public Entities
  !
  PUBLIC ! all public
  !
  ! ========== Reals
  !
  !  Values for real constants are prefixed with
  !  "RK_" to indicate the default "REAL KIND".
  !  Below, the values are of the right type because of the
  !  rules for intrinsic assignment.
  !
  !  --- Real form of integers
  !
  REAL(RK), PUBLIC, PARAMETER ::&
       &  RK_ZERO  = 0.0_RK, RK_ONE  = 1.0_RK, RK_TWO = 2.0_RK, &
       &  RK_THREE = 3.0_RK, RK_FOUR = 4.0_RK, RK_SIX = 6.0_RK, &
       &  RK_180   = 180.0_RK
  !
  !  --- Common fractions.
  !
  REAL(RK), PUBLIC, PARAMETER ::&
       &  RK_ONE_HALF = 0.5_RK, RK_ONE_THIRD = RK_ONE/RK_THREE
  !
  !  --- Things related to pi.
  !
  !  The values for pi and various other numbers were  computed to 50 
  !  decimal places using `bc'.  They could not be initialized here since 
  !  they involve intrinsic functions and would have to be set at run-time.
  !
  REAL(RK), PUBLIC, PARAMETER ::&
       &    RK_PI = 3.14159265358979323846264338327950288419716939937508_RK
  !
  REAL(RK), PUBLIC, PARAMETER ::&
       &       RK_TWO_PI = RK_PI*RK_TWO,  &
       &    RK_PI_OVER_2 = RK_PI/RK_TWO,  &
       &    RK_PI_OVER_3 = RK_PI/RK_THREE,&
       &    RK_PI_OVER_4 = RK_PI/RK_FOUR, &
       &    RK_PI_OVER_6 = RK_PI/RK_SIX,  &
       &  RK_PI_OVER_180 = RK_PI/RK_180,  &
       &  RK_180_OVER_PI = RK_180/RK_PI

  REAL(RK), PUBLIC, PARAMETER ::&             !  For convenience ...
       &  DEGREES_PER_RADIAN = RK_180_OVER_PI,&
       &  RADIANS_PER_DEGREE = RK_PI_OVER_180
  !
  !  --- Architecture-Related Values
  !
  REAL(RK), PUBLIC, PARAMETER :: RK_HUGE = HUGE(RK_ONE)
  !
  !  --- Algebraic Numbers
  !
  !  Various roots, golden section.
  !
  !
  REAL(RK), PUBLIC, PARAMETER ::&
       &   RK_ROOT_2 =&
       &      1.41421356237309504880168872420969807856967187537694_RK,&
       &   RK_ROOT_3 =&
       &      1.73205080756887729352744634150587236694280525381038_RK,&
       &   RK_ROOT_5 =&
       &      2.23606797749978969640917366873127623544061835961152_RK,&
       &   RK_ROOT_7 =&
       &      2.64575131106459059050161575363926042571025918308245_RK
  !
  REAL(RK), PUBLIC, PARAMETER ::&
       &   RK_ROOT_6        = RK_ROOT_2 * RK_ROOT_3,&
       &   RK_ROOT_14       = RK_ROOT_2 * RK_ROOT_7,&
       &   RK_ROOT_3_HALVES = RK_ROOT_3/RK_ROOT_2,&
       &   RK_GOLD          = RK_ONE_HALF * (RK_ROOT_5 + RK_ONE),&
       &   RK_GOLD_CONJ     = RK_ONE_HALF * (RK_ROOT_5 - RK_ONE)
  !
  ! ========== Integers
  !
  INTEGER(IK), PUBLIC, PARAMETER :: &
       &  IK_ZERO = 0_IK, IK_ONE = 1_IK, IK_TWO = 2_IK, IK_THREE = 3_IK,&
       &  IK_FOUR = 4_IK
  !
  ! ========== Logicals
  !
  LOGICAL(LK), PUBLIC, PARAMETER :: &
       &  TRUE = .TRUE., FALSE = .FALSE.
  !
  ! ========== Characters
  !
  !  ASCII numbers of various characters.
  !  
  CHARACTER, PUBLIC, PARAMETER ::&
       &  ASCII_TAB  = ACHAR(9),  ASCII_NULL  = ACHAR(0),&
       &  ASCII_BELL = ACHAR(7),  ASCII_BACK  = ACHAR(8),&
       &  ASCII_LF   = ACHAR(10), ASCII_VTAB  = ACHAR(11),&
       &  ASCII_CR   = ACHAR(13), ASCII_SPACE = ACHAR(32)
  !
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: WHITESPACE = ASCII_SPACE//ASCII_TAB
  !
END MODULE ConstantsModule
