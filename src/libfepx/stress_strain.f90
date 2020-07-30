! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
module stress_strain_mod
  ! 
  ! elastic strain and stress data
  use DimsModule
  USE IntrinsicTypesModule, RK=>REAL_KIND
  USE READ_INPUT_MOD
  use quadrature_mod
  !
  implicit none
  !
  !---------------------------------------------------------------------
  !
  public   ! all objects are public unless declared otherwise 
  !
  ! elastic strain internal variable
  REAL(RK), allocatable :: gela_kk_bar(:,:,:)
  ! stress data
  REAL(RK), allocatable  :: gsig_vec_n(:,:,:,:)
  
  !
contains
  !
  ! allocates memory for elastic strain and stress 
  !**********************************************************************
  !
  subroutine allocate_stress_strain(status)
    !
    implicit none
    ! 
    integer status 
    !
    !---------------------------------------------------------------------
    status = 0
    !
    ! allocate memory for variables 
    !
    allocate(gela_kk_bar(0:ngrain1, el_sub1:el_sup1, 0:nqpt1), &
             gsig_vec_n(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1), &
             stat=status)
    !
    ! initialize variables
    gela_kk_bar  = 0.0_RK
    gsig_vec_n   = 0.0_RK
    !
  end subroutine allocate_stress_strain
  !
end module stress_strain_mod
