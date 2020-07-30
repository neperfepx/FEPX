! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
!**********************************************************************
!
      SUBROUTINE isotropic(dii, t, s_pa, visc, dsdt)
!
!     Return viscosity for isotropic problem.
!
!     Apparently, this routine was disabled and was replaced by a
!     linear viscous problem with unit viscosity.
!
!----------------------------------------------------------------------
!
      USE READ_INPUT_MOD
      USE IntrinsicTypesModule, RK=>REAL_KIND
!
      IMPLICIT NONE
!
!
!     compute the viscosity for an element based on the effective
!     strain rate.  uses the fit of dawson to hart's model for pure
!     aluminum.
!
!     input arguments -
!     
!     dii    - the effective strain rate for the element
!     t      - element temperature (deg f)
!     s_pa   - the value of the state variable (in Pascals)
!
      REAL(RK) dii(el_sub1:el_sup1), t(el_sub1:el_sup1)
      REAL(RK) s(el_sub1:el_sup1), s_pa(el_sub1:el_sup1)
!
!     output arguments -
!
!     visc   - the resulting viscosity
!     dsdt   - the derivative of the state variable
!
      REAL(RK) visc(el_sub1:el_sup1), dsdt(el_sub1:el_sup1)
!
!     modifications -
!
!     the effective viscosity computation has been changed to
!     visc = (2nd invariant of sigma) / 3*(2nd invariant of strain rate)
!     so as to be compatible with definition of strain rate invariant
!     by ajb (9/25/87).
!
!     note - units are mpa, kjoule/mole-k, etc.
!
!     constants -
!
!     REAL(RK) a0
!     PARAMETER ( a0     = 9.64d52     )
!
      REAL(RK) log_a0
      PARAMETER ( log_a0 = 122.0       )

      REAL(RK) c0, qpr, ge, m, f0, qr, smallm, lambda, mp, n

      PARAMETER ( c0     = 6.19d-9     )
      PARAMETER ( qpr    = 1.45d4      )
      PARAMETER ( ge     = 24.20d3     )
      PARAMETER ( m      = 7.80        )
      PARAMETER ( f0     = 2.12d19     )
      PARAMETER ( qr     = qpr         )
      PARAMETER ( smallm = 5.0         )
      PARAMETER ( lambda = 0.14        )
      PARAMETER ( mp     = 3.5         )
      PARAMETER ( n      = 6.0         )
!
      REAL(RK) a(el_sub1:el_sup1)
      REAL(RK) dstar(el_sub1:el_sup1)
      REAL(RK) sigma(el_sub1:el_sup1)
      REAL(RK) sigmap(el_sub1:el_sup1)
      REAL(RK) sigmav(el_sub1:el_sup1)
!
!     Locals:
!
      INTEGER i
!
!----------------------------------------------------------------------
!
! 
!     Place the state variable in Pascals.
!
      visc = 1.0
!
      RETURN

      s = s_pa / 1.0d6

!     compute flow stress in viscous (friction) element.
 
!     a = a0 * exp(-qpr/t)

      a = log_a0 + (-qpr/t)
      a = exp(a)

      sigmav = ge * (dii/a)**(1.0/m)
!
!     compute flow stress in the plastic element.
!
      dstar  = f0 * ((s/ge)**smallm) * exp(-qr/t)
      sigmap = s * exp( - (dstar/dii)**lambda )
 
!     total flow stress 
 
      sigma = sigmav + sigmap

!     compute the derivative of the state vaiable
 
      dsdt = c0*s*dii*((ge/s)**mp)*((sigmap/s)**n)

!     Convert to Pascals

!     sigma = s * dii**0.10
      sigma = sigma * 1.0e6

!     compute viscosity
 
      visc = sigma / (3.d0*dii)

      print *,'sigma', minval(sigma), maxval(sigma)
      print *,'visc', minval(visc), maxval(visc)
 
      RETURN
      END
!      
