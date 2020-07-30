! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
module post_update_n_mod
  ! 
  ! elastic strain and stress data
  use DimsModule
  USE IntrinsicTypesModule, RK=>REAL_KIND
  USE READ_INPUT_MOD
  !
  implicit none
  !
  !---------------------------------------------------------------------
  !
  public   ! all objects are public unless declared otherwise 
  !
  ! elastic strain internal variable
  REAL(RK), allocatable :: pela_kk_bar(:)
  ! stress data
  REAL(RK), allocatable  :: psig_vec_n(:,:,:)
  ! critical resolved stress
  REAL(RK), allocatable  :: pcrss_n(:,:,:)
  ! R^*_n 
  REAL(RK), allocatable  :: prstar_n(:,:,:,:)
  !
contains
  !
  ! allocates memory for elastic strain and stress 
  !**********************************************************************
  !
  subroutine allocate_post_update_n(status)
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
    allocate(pela_kk_bar(el_sub1:el_sup1), &
             psig_vec_n(0:TVEC1, 0:ngrain1, el_sub1:el_sup1), &
             stat=status)
    ! initialize variables
    pela_kk_bar = 0.0d0
    psig_vec_n = 0.0d0
    !
  end subroutine allocate_post_update_n
end module post_update_n_mod
