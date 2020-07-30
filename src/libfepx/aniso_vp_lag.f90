! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
!
!     Routines for  anisotropic stiffness.
!
      SUBROUTINE invert5x5(a, n, m)
!
!     Invert array of 5x5 matrices.
!
!------------------------------------------------------------------------
!
      USE  DimsModule
!
      IMPLICIT  NONE
!
!     Arguments:
!
!     a: array of matrices on input and array of inverses on output
!     n: number of grains
!     m: number of elements
!
      INTEGER   ::  n, m
      REAL(RK)  ::  a(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER   ::  i, j
      REAL(RK)  ::  a11(0:(n - 1), 0:(m - 1)), a21(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  a22(0:(n - 1), 0:(m - 1)), a31(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  a32(0:(n - 1), 0:(m - 1)), a33(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  a41(0:(n - 1), 0:(m - 1)), a42(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  a43(0:(n - 1), 0:(m - 1)), a44(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  a51(0:(n - 1), 0:(m - 1)), a52(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  a53(0:(n - 1), 0:(m - 1)), a54(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  a55(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  ai11(0:(n - 1), 0:(m - 1)), ai21(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  ai22(0:(n - 1), 0:(m - 1)), ai31(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  ai32(0:(n - 1), 0:(m - 1)), ai33(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  ai41(0:(n - 1), 0:(m - 1)), ai42(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  ai43(0:(n - 1), 0:(m - 1)), ai44(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  ai51(0:(n - 1), 0:(m - 1)), ai52(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  ai53(0:(n - 1), 0:(m - 1)), ai54(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  ai55(0:(n - 1), 0:(m - 1))
!
!------------------------------------------------------------------------

      a11 = a(0, 0, :, :)
      a21 = a(1, 0, :, :)
      a22 = a(1, 1, :, :)
      a31 = a(2, 0, :, :)
      a32 = a(2, 1, :, :)
      a33 = a(2, 2, :, :)
      a41 = a(3, 0, :, :)
      a42 = a(3, 1, :, :)
      a43 = a(3, 2, :, :)
      a44 = a(3, 3, :, :)
      a51 = a(4, 0, :, :)
      a52 = a(4, 1, :, :)
      a53 = a(4, 2, :, :)
      a54 = a(4, 3, :, :)
      a55 = a(4, 4, :, :)

! **  A = LDL'.
!  j = 1
!
      ai55 = 1.0 / a11
      a21 = a21 * ai55
      a31 = a31 * ai55
      a41 = a41 * ai55
      a51 = a51 * ai55
!
!  j = 2
!
      ai11 = a21 * a11
      a22 = a22 - a21 * ai11
         
      a32 = a32 - a31 * ai11
      a42 = a42 - a41 * ai11
      a52 = a52 - a51 * ai11
      ai55 = 1.0 / a22
      a32 = a32 * ai55
      a42 = a42 * ai55
      a52 = a52 * ai55
!
!  j = 3
!
      ai11 = a31 * a11
      ai22 = a32 * a22
      ai55 = a31 * ai11 + a32 * ai22

      a33 = a33 - ai55
         
      ai55 = 1.0 / a33
      a43 = a43 - a41 * ai11 - a42 * ai22
      a53 = a53 - a51 * ai11 - a52 * ai22
      a43 = a43 * ai55
      a53 = a53 * ai55
!
!  j = 4
!
      ai11 = a41 * a11
      ai22 = a42 * a22
      ai33 = a43 * a33
      ai55 = a41 * ai11 + a42 * ai22 + a43 * ai33

      a44 = a44 - ai55
         
      a54 = a54 - a51 * ai11 - a52 * ai22 - a53 * ai33
      a54 = a54 / a44
!
!  j = 5
!
      ai11 = a51 * a11
      ai22 = a52 * a22
      ai33 = a53 * a33
      ai44 = a54 * a44
      ai55 = a51 * ai11 + a52 * ai22 + a53 * ai33 + a54 * ai44

      a55 = a55 - ai55
!
!  Column 1 of inverse 
!  Ly = b
!
      ai21 = - a21
      ai31 = - a31 - a32 * ai21
      ai41 = - a41 - a42 * ai21 - a43 * ai31
      ai51 = - a51 - a52 * ai21 - a53 * ai31 - a54 * ai41
!
!  Dz = y
!
      ai11 = 1.0 / a11
      ai21 = ai21 / a22
      ai31 = ai31 / a33
      ai41 = ai41 / a44
      ai51 = ai51 / a55
!
!  L'x = z
!
      ai41 = ai41 - a54 * ai51
      ai31 = ai31 - a43 * ai41 - a53 * ai51
      ai21 = ai21 - a32 * ai31 - a42 * ai41 - a52 * ai51
      ai11 = ai11 - a21 * ai21 - a31 * ai31 - a41 * ai41 - a51 * ai51
!
!  Column 2 of inverse 
!  Ly = b
!
      ai32 = - a32
      ai42 = - a42 - a43 * ai32
      ai52 = - a52 - a53 * ai32 - a54 * ai42
!
!  Dz = y
!
      ai22 = 1.0 / a22
      ai32 = ai32 / a33
      ai42 = ai42 / a44
      ai52 = ai52 / a55
!
!  L'x = z
!
      ai42 = ai42 - a54 * ai52
      ai32 = ai32 - a43 * ai42 - a53 * ai52
      ai22 = ai22 - a32 * ai32 - a42 * ai42 - a52 * ai52
!
!  Column 3 of inverse 
!  Ly = b
!
      ai43 = - a43
      ai53 = - a53 - a54 * ai43
!
!  Dz = y
!
      ai33 = 1.0 / a33
      ai43 = ai43 / a44
      ai53 = ai53 / a55
!
!  L'x = z
!
      ai43 = ai43 - a54 * ai53
      ai33 = ai33 - a43 * ai43 - a53 * ai53
!
!  Column 4 of inverse 
!  Ly = b
!
      ai54 = - a54
!
!  Dz = y
!
      ai44 = 1.0 / a44
      ai54 = ai54 / a55
!
!  L'x = z
!
      ai44 = ai44 - a54 * ai54
!
!  Column 5 of inverse 
!  Dz = y
!
      ai55 = 1.0 / a55
!
!  Recover Array
!
      a(0, 0, :, :) = ai11
      a(1, 0, :, :) = ai21
      a(2, 0, :, :) = ai31
      a(3, 0, :, :) = ai41
      a(4, 0, :, :) = ai51
      a(1, 1, :, :) = ai22
      a(2, 1, :, :) = ai32
      a(3, 1, :, :) = ai42
      a(4, 1, :, :) = ai52
      a(2, 2, :, :) = ai33
      a(3, 2, :, :) = ai43
      a(4, 2, :, :) = ai53
      a(3, 3, :, :) = ai44
      a(4, 3, :, :) = ai54
      a(4, 4, :, :) = ai55

      do i = 0, TVEC1
         do j = i + 1, TVEC1
            a(i, j, :, :) = a(j, i, :, :)
         enddo
      enddo

      RETURN
      END
!
      
