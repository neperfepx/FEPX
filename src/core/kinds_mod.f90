! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Kind parameters used throughout FEPX.
!>
!> Source : dev/conventions.md §6 ("Kinds"). Contains no numerical method.
!>
!> Every FEPX source file takes its kinds from this module (`use kinds_mod, only: rk, ...`), so
!> that the kind choice is made in exactly one place. The names and values are those of the
!> conventions: rk => real64, ik => int32, i8 => int64.
!>
!> Public parameters:
!>   rk   integer  kind of all real values (IEEE binary64, "double precision")
!>   ik   integer  kind of default integers: indices, counts local to one rank (32 bit)
!>   i8   integer  kind of integers that may exceed 2^31 - 1, e.g. global counts (64 bit)
module kinds_mod

  use, intrinsic :: iso_fortran_env, only: real64, int32, int64

  implicit none

  private

  public :: rk, ik, i8

  integer, parameter :: rk = real64
  integer, parameter :: ik = int32
  integer, parameter :: i8 = int64

end module kinds_mod
