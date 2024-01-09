! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module general_mod

  use, intrinsic:: iso_fortran_env, only: rk => real64

  public

  ! mesh partitionning
  integer :: num_procs
  integer :: myid
  integer :: elt_sub, elt_sup, node_sub, node_sup, dof_sub, dof_sup
  integer, allocatable :: list_dof_sub(:), list_dof_sup(:)

  real(rk) :: vtiny = 1.0d-16
  ! FEM parameters and other stuff
  integer, parameter :: ndim = 10 ! 10-noded tetrahedra
  integer, parameter :: kdim = 3*ndim
  integer, parameter :: nqpt = 15
  real(rk) :: qploc(3, nqpt)
  real(rk) :: wtqp(3, nqpt)
  real(rk) :: qp2d(2, 7)
  real(rk) :: wt2d(7)
  integer, parameter :: sftype = 6
  real(rk) :: sfqp2d(sftype, 7)
  real(rk) :: sfgqp2d(2, sftype, 7)

  !> Central quad point inside an element (used for post-processing)
  integer :: cqpt = 5

end module general_mod
