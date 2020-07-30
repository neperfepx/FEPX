! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE MaterialMatrixVpModule

  USE parallel_mod
  USE DimsModule
  USE READ_INPUT_MOD
  use StressSolveVpModule
  USE KinematicsModule

  IMPLICIT  NONE

  PRIVATE
  PUBLIC :: material_matrix_vp, defrate

CONTAINS

      SUBROUTINE material_matrix_vp(&
  &   type, stif, dndx, dndy, dndz, gvel, scale,&
  &   qr5x5, wts, epseff, dtime, incr)

!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER   ::  type, incr
      REAL(RK)  ::  dtime
      REAL(RK)  ::  stif(TVEC, TVEC, el_sub1:el_sup1)
      REAL(RK)  ::  dndx(0:nnpe, el_sub1:el_sup1), dndy(0:nnpe, el_sub1:el_sup1)
      REAL(RK)  ::  dndz(0:nnpe, el_sub1:el_sup1), gvel(0:kdim1, el_sub1:el_sup1)
      REAL(RK)  ::  scale(el_sub1:el_sup1)
      REAL(RK)  ::  qr5x5(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  ::  wts(0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  ::  epseff(el_sub1:el_sup1)
      ! tm268_M13:
!      REAL(RK)  ::  eqplas(el_sub1:el_sup1)      
!
!     Locals:
!    
      INTEGER i, j, m_el
      REAL(RK)&
  &   d(0:DIMS1, 0:DIMS1, el_sub1:el_sup1),&
  &   d_vec(0:TVEC1, el_sub1:el_sup1)
      REAL(RK) cmu(el_sub1:el_sup1), temp_k(el_sub1:el_sup1)
      REAL(RK) dsdt_iso(el_sub1:el_sup1), state_iso(el_sub1:el_sup1)
      REAL(RK) cconst(TVEC)
!
!     Data:
!
      DATA cconst /1.0d0, 3.0d0, 4.0d0, 4.0d0, 4.0d0/
!
!----------------------------------------------------------------------
!
      m_el = el_sup1 - el_sub1 + 1
!
      call defrate(d, dndx, dndy, dndz, gvel)

      ! tm268_M13:
!      call eff_def(epseff, eqplas, d, dtime, m_el)
      call eff_def(epseff, d, dtime, m_el)

! hr-tm
! we aren't going through this portion to update for two-phase compatibility
      if (type .eq. ANISOTROPIC_VP) then
         call par_quit('Error  :     > ANISOTROPIC_VP is no longer implemented.')
      else if (type .eq. ISOTROPIC_VP) then
!
!dbg    Note the hardwired temperature and state variable.
!
        temp_k    = 273.0 + 400.0
        state_iso = 20.0e06

        call isotropic(epseff, temp_k, state_iso, cmu, dsdt_iso)

        stif = 0.0
        stif(1, 1, :) = cconst(1) * cmu
        stif(2, 2, :) = cconst(2) * cmu
        stif(3, 3, :) = cconst(3) * cmu
        stif(4, 4, :) = cconst(4) * cmu
        stif(5, 5, :) = cconst(5) * cmu

      endif

      scale = 0.0
      do i = 1, TVEC
        scale = scale + stif(i, i, :) / cconst(i)
      enddo
      
      scale = scale / 5.0

      END SUBROUTINE material_matrix_vp
      
END MODULE MaterialMatrixVpModule
