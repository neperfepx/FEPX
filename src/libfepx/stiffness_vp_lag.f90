! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE StiffnessVpModule

  USE parallel_matrix_mod

  USE IntrinsicTypesModule, RK=>REAL_KIND
  USE quadrature_mod
  USE READ_INPUT_MOD
  USE shape_3d_mod
  USE DimsModule
  USE ConvergenceModule, ONLY: cv_options
  USE MaterialMatrixVpModule

  IMPLICIT NONE
  
  PRIVATE 

  PUBLIC :: element_stif_vp, add_to_stiffness

CONTAINS

      SUBROUTINE element_stif_vp(&
  &   itype, gstiff, gcoords, gvel, pscale, pcnst,&
  &   qr5x5, wts, eqplas, epseff,&
  &   dtime, incr)
!
!     Form elemental stiffness matrix for viscoplastic problem.
!
!----------------------------------------------------------------------
!
!     Arguments:
!
!     itype: 0=isotropic,1=anisotropic
!     incr: current increment
!
      INTEGER  ::  itype, incr
      REAL(RK) ::  dtime
!
      REAL(RK) ::  gstiff(0:kdim1, 0:kdim1, el_sub1:el_sup1)
      REAL(RK) ::  gcoords(0:kdim1, el_sub1:el_sup1)
      REAL(RK) ::  gvel(0:kdim1, el_sub1:el_sup1)
      REAL(RK) ::  pscale(el_sub1:el_sup1), pcnst(0:kdim1, el_sub1:el_sup1)
      REAL(RK) ::  qr5x5(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK) ::  wts(0:ngrain1, el_sub1:el_sup1)
      REAL(RK) ::  eqplas(el_sub1:el_sup1), epseff(el_sub1:el_sup1)
!
!     Locals:
!
      INTEGER  ::  iqpt, jqpt, kqpt
      INTEGER  ::  i, j, i1, i2, i3, j1, j2, j3, ier
      REAL(RK) ::  wt
      REAL(RK) ::  dndx(0:nnpe, el_sub1:el_sup1), dndy(0:nnpe, el_sub1:el_sup1)
      REAL(RK) ::  dndz(0:nnpe, el_sub1:el_sup1), det(el_sub1:el_sup1)
      REAL(RK) ::  s11(el_sub1:el_sup1), s12(el_sub1:el_sup1), s13(el_sub1:el_sup1)
      REAL(RK) ::  s21(el_sub1:el_sup1), s22(el_sub1:el_sup1), s23(el_sub1:el_sup1)
      REAL(RK) ::  s31(el_sub1:el_sup1), s32(el_sub1:el_sup1), s33(el_sub1:el_sup1)
      REAL(RK) ::  t11(el_sub1:el_sup1)
      REAL(RK) ::  mmtx(el_sub1:el_sup1), sclfac(el_sub1:el_sup1)
      REAL(RK) ::  c(TVEC, TVEC, el_sub1:el_sup1)
      REAL(RK) ::  xni(3, 5, el_sub1:el_sup1), xnj(5, 3, el_sub1:el_sup1)
      REAL(RK) ::  temp1(3, 5, el_sub1:el_sup1), temp2(3, 3, el_sub1:el_sup1)
      REAL(RK) ::  loc0, loc1, loc2
!
!----------------------------------------------------------------------
!
      pcnst  = 0.0
      mmtx   = 0.0
      pscale = 0.0
!
      loc0 = 0.25
      loc1 = 0.25 
      loc2 = 0.25
      
      call sfder_hpar(loc0, loc1, loc2, gcoords, dndx, dndy, dndz,&
  &   det, s11, s12, s13, s21, s22, s23, s31, s32, s33)

      ! 
      call material_matrix_vp(itype, c, dndx, dndy, dndz, gvel,&
  &   sclfac, qr5x5, wts, epseff,&
  &   dtime, incr)



      t11 = 0.0

      do iqpt = 0, nqpt1
         loc0 = qploc(0, iqpt)
         loc1 = qploc(1, iqpt)
         loc2 = qploc(2, iqpt)
             
         wt = wtqp(0, iqpt) 
             
         call sfder_hpar(loc0, loc1, loc2, gcoords, dndx, dndy,&
  &      dndz, det, s11, s12, s13, s21, s22, s23, s31, s32,&
  &      s33)

         det = det * wt
         pscale = pscale + sclfac * wt    
         mmtx = mmtx + det


         do i = 0, nnpe
           i1 = 3 * i
           i2 = i1 + 1
           i3 = i2 + 1

           pcnst(i1, :) = pcnst(i1, :) - dndx(i, :) * det
           pcnst(i2, :) = pcnst(i2, :) - dndy(i, :) * det
           pcnst(i3, :) = pcnst(i3, :) - dndz(i, :) * det

           xni(1, 1, :) = dndx(i, :)
           xni(2, 1, :) = -dndy(i, :)
           xni(3, 1, :) = 0.0
           xni(1, 2, :) = -dndx(i, :) / 3.0
           xni(2, 2, :) = -dndy(i, :) / 3.0
           xni(3, 2, :) = 2.0 * dndz(i, :) / 3.0
           xni(1, 3, :) = 0.5 * dndy(i, :)
           xni(2, 3, :) = 0.5 * dndx(i, :)
           xni(3, 3, :) = 0.0
           xni(1, 4, :) = 0.5 * dndz(i, :)
           xni(2, 4, :) = 0.0
           xni(3, 4, :) = 0.5 * dndx(i, :)
           xni(1, 5, :) = 0.0
           xni(2, 5, :) = 0.5 * dndz(i, :)
           xni(3, 5, :) = 0.5 * dndy(i, :)

           temp1 = 0.0
           call gen_matrix_mult(temp1, xni, c, 1, 2, ier)

           do j = 0, i
             j1 = 3 * j
             j2 = j1 + 1
             j3 = j2 + 1

             xnj(1, 1, :) = dndx(j, :)
             xnj(1, 2, :) = -dndy(j, :)
             xnj(1, 3, :) = 0.0
             xnj(2, 1, :) = -dndx(j, :) / 3.0
             xnj(2, 2, :) = -dndy(j, :) / 3.0
             xnj(2, 3, :) = 2.0 * dndz(j, :)  /3.0
             xnj(3, 1, :) = 0.5 * dndy(j, :)
             xnj(3, 2, :) = 0.5 * dndx(j, :)
             xnj(3, 3, :) = 0.0
             xnj(4, 1, :) = 0.5 * dndz(j, :)
             xnj(4, 2, :) = 0.0
             xnj(4, 3, :) = 0.5 * dndx(j, :)
             xnj(5, 1, :) = 0.0
             xnj(5, 2, :) = 0.5 * dndz(j, :)
             xnj(5, 3, :) = 0.5 * dndy(j, :)

             temp2 = 0.0
             call gen_matrix_mult(temp2, temp1, xnj,1, 2, ier)

             s11 = temp2(1, 1, :) * det
             s12 = temp2(1, 2, :) * det
             s13 = temp2(1, 3, :) * det
             s21 = temp2(2, 1, :) * det
             s22 = temp2(2, 2, :) * det
             s23 = temp2(2, 3, :) * det
             s31 = temp2(3, 1, :) * det
             s32 = temp2(3, 2, :) * det
             s33 = temp2(3, 3, :) * det

             call add_to_stiffness(gstiff, s11, i1, j1)
             call add_to_stiffness(gstiff, s22, i2, j2)
             call add_to_stiffness(gstiff, s33, i3, j3)
             call add_to_stiffness(gstiff, s12, i1, j2)
             call add_to_stiffness(gstiff, s13, i1, j3)
             call add_to_stiffness(gstiff, s21, i2, j1)
             call add_to_stiffness(gstiff, s23, i2, j3)
             call add_to_stiffness(gstiff, s31, i3, j1)
             call add_to_stiffness(gstiff, s32, i3, j2)                   

           enddo
        enddo
      enddo

      pscale = pscale / mmtx
      t11 = cv_options % pacc * pscale

      do i = 0, nnpe
        i1 = 3 * i
        i2 = i1 + 1
        i3 = i2 + 1

        do j = 0, i
          j1 = 3 * j
          j2 = j1 + 1
          j3 = j2 + 1

          s11 = pcnst(i1, :) * pcnst(j1, :) * t11
          s12 = pcnst(i1, :) * pcnst(j2, :) * t11
          s13 = pcnst(i1, :) * pcnst(j3, :) * t11
          s21 = pcnst(i2, :) * pcnst(j1, :) * t11
          s22 = pcnst(i2, :) * pcnst(j2, :) * t11
          s23 = pcnst(i2, :) * pcnst(j3, :) * t11
          s31 = pcnst(i3, :) * pcnst(j1, :) * t11
          s32 = pcnst(i3, :) * pcnst(j2, :) * t11
          s33 = pcnst(i3, :) * pcnst(j3, :) * t11

          call add_to_stiffness(gstiff, s11, i1, j1)
          call add_to_stiffness(gstiff, s22, i2, j2)
          call add_to_stiffness(gstiff, s33, i3, j3)
          call add_to_stiffness(gstiff, s12, i1, j2)
          call add_to_stiffness(gstiff, s13, i1, j3)
          call add_to_stiffness(gstiff, s21, i2, j1)
          call add_to_stiffness(gstiff, s23, i2, j3)
          call add_to_stiffness(gstiff, s31, i3, j1)
          call add_to_stiffness(gstiff, s32, i3, j2)                   
        enddo
      enddo


! 8    FORMAT(/10x,'===SUBROUTINE stiffness===')
! 9    FORMAT(/'c(5,5,',i2,') at the center :' / (5(2x, e12.5)))
! 10   FORMAT(/'elemental stiffness (24,24,0) :')
! 11   FORMAT (/(6(2x,e12.5)))

      END SUBROUTINE element_stif_vp

      SUBROUTINE add_to_stiffness(gstiff, sij, i, j)
!
!     Add component to stiffness matrix.
!
!----------------------------------------------------------------------
!
!     Arguments: (self-explanatory)
!
      INTEGER  i, j
      REAL(RK), INTENT(INOUT) ::  gstiff(0:kdim1, 0:kdim1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    ::  sij(el_sub1:el_sup1)
!
!----------------------------------------------------------------------
!
      if (i .ge. j) then
        if (i .eq. j) then
          gstiff(i, j, :) = gstiff(i, j, :) + sij
        else
          gstiff(i, j, :) = gstiff(i, j, :) + sij
          gstiff(j, i, :) = gstiff(j, i, :) + sij
        endif
      endif

      END SUBROUTINE add_to_stiffness
      
END MODULE StiffnessVpModule
