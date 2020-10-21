! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE MaterialMatrixEvpsModule

  USE IntrinsicTypesModule, RK=>REAL_KIND
  USE DimsModule
  USE units_mod

  use KinematicsModule
  USE READ_INPUT_MOD
  USE microstructure_mod
  USE MATRIX_OPERATIONS_MOD
  USE PolycrystalResponseEvpsModule
  USE ANISO_EVPS_MOD

  IMPLICIT NONE

  PRIVATE 
  PUBLIC :: material_matrix_evps, vel_gradient, eff_def

CONTAINS

      SUBROUTINE material_matrix_evps(&
  &   stif, tan_stif, fe, detv, dndx, dndy, dndz, gvel,&
  &   c0_angs, c_angs, sig_vec_n, sig_vec, crss_n,&
  &   crss, rstar_n, rstar, keinv,&
  &   e_bar_vec, e_elas_kk_bar, e_elas_kk, sig_kk, jiter_state,&
  &   wts, epseff, &
  &   dtime, incr,&
  &   converged_solution, auto_time, NR&
  &   )
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      LOGICAL, INTENT(INOUT) :: converged_solution
!
      INTEGER :: incr, auto_time
      INTEGER, INTENT(OUT)  :: jiter_state(0:ngrain1, el_sub1:el_sup1)
!
      REAL(RK) :: dtime
!      
      REAL(RK), INTENT(IN)    ::  keinv(0:TVEC1,1:numphases)
      REAL(RK), INTENT(OUT)   ::  stif(TVEC, TVEC, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT)   ::  tan_stif(TVEC, TVEC, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT)   ::  fe(TVEC, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT)   ::  detv(el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    ::  dndx(0:nnpe, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    ::  dndy(0:nnpe, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    ::  dndz(0:nnpe, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    ::  gvel(0:kdim1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    ::  c0_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) ::  c_angs (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    ::  rstar_n(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) ::  rstar  (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    ::  sig_vec_n(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT)   ::  sig_vec  (0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    ::  e_bar_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    ::  crss_n(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(OUT)   ::  crss  (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    ::  wts   (0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    ::  e_elas_kk_bar(el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT)   ::  e_elas_kk    (el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT)   ::  sig_kk(el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT)   ::  epseff(el_sub1:el_sup1, 0:nqpt1)
      LOGICAL,  INTENT(IN)    :: NR
      ! tm268_M13:
!      REAL(RK), INTENT(INOUT) ::  eqplas(el_sub1:el_sup1)
!
!     Locals:
!    
      INTEGER  ::  i, j, m_el
      REAL(RK) ::  d    (0:DIMS1, 0:DIMS1, el_sub1:el_sup1)
      REAL(RK) ::  w    (0:DIMS1, 0:DIMS1, el_sub1:el_sup1)      
      REAL(RK) ::  vgrad(0:DIMS1, 0:DIMS1, el_sub1:el_sup1)
      REAL(RK) ::  w_vec(0:DIMS1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK) ::  d_vec(0:TVEC1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK) ::  d_kk(el_sub1:el_sup1, 0:nqpt1)
!
!----------------------------------------------------------------------
!
      m_el = el_sup1 - el_sub1 + 1

    do i = 0, nqpt1
      
      ! Compute velocity gradient (vgrad) and its sym. (d) and skew (w) parts.
      call vel_gradient(vgrad, dndx(:,:,i), dndy(:,:,i), dndz(:,:,i), gvel)
      ! d_kk: mean/volumetric/spherical part of the symmetric part of the velocity gradient
      ! d:    deviatoric part of the symmetric part of the velocity gradient
      ! w:    skew part of the velocity gradient
      call symm_vgr(d, d_kk(:, i), vgrad, m_el)
      call skew_vgr(w, vgrad, m_el)
      ! d [3x3] --> d_vec {5}
      ! w [3x3] --> w_vec {3}
      call mat_vec_symm(d, d_vec(:, :, i), m_el)
      call mat_vec_skew(w, w_vec(:, :, i), m_el)
      ! tm268_M13:
!      ! d,dtime --> epseff,eqplas
!      call eff_def(epseff, eqplas, d, dtime, m_el)
      call eff_def(epseff(:, i), d, dtime, m_el)
    enddo
      
      ! variables @(t+dt) ? :
      ! - rstar     [3x3]
      ! - c_angs    [3x3]
      ! - crss      (1)
      ! - sig_vec   {5}
      ! - sig_kk    (1)
      ! - e_elas_kk (1)
      call polycrystal_response_evps(d_vec, w_vec,&
  &        c0_angs, c_angs, sig_vec_n, sig_vec, crss_n, crss, rstar_n,&
  &        rstar, e_bar_vec, wts, epseff, d_kk, sig_kk,&
  &        e_elas_kk_bar, e_elas_kk, jiter_state, keinv,&
  &        incr, &
  &        dtime, converged_solution, auto_time)
     
      if (.not. converged_solution .and. auto_time .eq. 1) RETURN

    do i = 0, nqpt1
      
      call aniso_evps(stif(:,:,:,i), tan_stif(:,:,:,i), fe(:,:,i), detv(:,i), c_angs(:,:,:,:,i),&
  &        sig_vec(:,:,:,i), crss(:,:,:,i), rstar_n, rstar(:,:,:,:,i), e_bar_vec(:,:,:,i), &
  &        wts, w_vec(:,:,i), e_elas_kk(:,i), keinv, dtime, NR)

    enddo
!
      END SUBROUTINE material_matrix_evps

END MODULE MaterialMatrixEvpsModule
