! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE KinematicsModule
  !
  !  kinematic calculations
  !
  !  ==================== Other Modules (USE statements)
  !
  USE LIBF95, RK=>REAL_KIND
  
  USE DimsModule
  USE READ_INPUT_MOD

  IMPLICIT NONE
  !
  !  ==================== Public Entities
  !
  !  variables, procedures, constants, derived types and namelist groups
  !
  PRIVATE   ! all objects are private unless declared otherwise

  PUBLIC :: vel_gradient, eff_def, defrate
  !
CONTAINS ! ============================================= MODULE PROCEDURES
  !
  SUBROUTINE vel_gradient(vgrad, dndx, dndy, dndz, gvel)
    !
    !     Compute velocity gradient.
    !
    !----------------------------------------------------------------------
    !
    !     Arguments:
    !
    REAL(RK), INTENT(OUT) :: vgrad(0:DIMS1, 0:DIMS1, el_sub1:el_sup1)
    REAL(RK), INTENT(IN)  :: dndx(0:nnpe,  el_sub1:el_sup1)
    REAL(RK), INTENT(IN)  :: dndy(0:nnpe,  el_sub1:el_sup1)
    REAL(RK), INTENT(IN)  :: dndz(0:nnpe,  el_sub1:el_sup1)
    REAL(RK), INTENT(IN)  :: gvel(0:kdim1, el_sub1:el_sup1)
    !
    !     Locals:
    !    
    INTEGER :: i, i1, i2, i3
    !
    !----------------------------------------------------------------------
    !
    vgrad = 0.0_RK

    do i = 0, nnpe
      i1 = 3 * i
      i2 = i1 + 1
      i3 = i2 + 1

      vgrad(0, 0, :) = vgrad(0, 0, :) + dndx(i, :) * gvel(i1, :)
      vgrad(0, 1, :) = vgrad(0, 1, :) + dndy(i, :) * gvel(i1, :)
      vgrad(0, 2, :) = vgrad(0, 2, :) + dndz(i, :) * gvel(i1, :)

      vgrad(1, 0, :) = vgrad(1, 0, :) + dndx(i, :) * gvel(i2, :)
      vgrad(1, 1, :) = vgrad(1, 1, :) + dndy(i, :) * gvel(i2, :)
      vgrad(1, 2, :) = vgrad(1, 2, :) + dndz(i, :) * gvel(i2, :)

      vgrad(2, 0, :) = vgrad(2, 0, :) + dndx(i, :) * gvel(i3, :)
      vgrad(2, 1, :) = vgrad(2, 1, :) + dndy(i, :) * gvel(i3, :)
      vgrad(2, 2, :) = vgrad(2, 2, :) + dndz(i, :) * gvel(i3, :)
    enddo

  END SUBROUTINE vel_gradient

  SUBROUTINE eff_def(epseff, d, dtime, m)
    !
    !     Compute the effective deformation rate and accumulated deformation
    !     (strain) from the deformation rate tensor.
    !
    !----------------------------------------------------------------------
    !
    INTEGER, INTENT(IN)     :: m
    REAL(RK), INTENT(IN)    :: dtime
    REAL(RK), INTENT(IN)    :: d(0:DIMS1, 0:DIMS1, 0:(m - 1))
    REAL(RK), INTENT(OUT)   :: epseff(0:(m - 1))
    !
    !     Locals.
    !
    REAL(RK), PARAMETER :: TWOTHIRDS=2.0_RK/3.0_RK

    INTEGER :: i
    !
    !----------------------------------------------------------------------
    !
    ! effective deformation rate
    do i=0, m-1
      epseff(i) = sqrt(TWOTHIRDS*&
           &     ( d(0, 0, i)*d(0, 0, i)&
           &     + d(1, 1, i)*d(1, 1, i)&
           &     + d(2, 2, i)*d(2, 2, i) + 2.0 * (&
           &         d(0, 1, i)*d(0, 1, i)&
           &       + d(0, 2, i)*d(0, 2, i)&
           &       + d(1, 2, i)*d(1, 2, i)&
           &       )&
           &     )&
           &     )
    enddo

  END SUBROUTINE eff_def

  SUBROUTINE defrate(d, dndx, dndy, dndz, gvel)
    !
    !----------------------------------------------------------------------
    !
    !     Arguments:
    !
    REAL(RK)    d(0:DIMS1, 0:DIMS1, el_sub1:el_sup1)
    REAL(RK)    dndx(0:nnpe, el_sub1:el_sup1), dndy(0:nnpe, el_sub1:el_sup1)
    REAL(RK)    dndz(0:nnpe, el_sub1:el_sup1), gvel(0:kdim1, el_sub1:el_sup1)
    !
    !     Locals:
    !
    INTEGER   i, i1, i2, i3
    REAL(RK)    divv(el_sub1:el_sup1)
    !
    !----------------------------------------------------------------------
    !
    d = 0.d0

    do i = 0, nnpe
      i1 = 3 * i
      i2 = i1 + 1
      i3 = i2 + 1

      d(0, 0, :) = d(0, 0, :) + dndx(i, :) * gvel(i1, :)
      d(1, 1, :) = d(1, 1, :) + dndy(i, :) * gvel(i2, :)
      d(2, 2, :) = d(2, 2, :) + dndz(i, :) * gvel(i3, :)
      !
      d(1, 0, :) = d(1, 0, :) + dndx(i, :) * gvel(i2, :)&
           &      + dndy(i, :) * gvel(i1, :)
      d(2, 0, :) = d(2, 0, :) + dndx(i, :) * gvel(i3, :)&
           &      + dndz(i, :) * gvel(i1, :)
      d(2, 1, :) = d(2, 1, :) + dndy(i, :) * gvel(i3, :)&
           &      + dndz(i, :) * gvel(i2, :)
      !
    enddo

    d(1, 0, :) = 0.5 * d(1, 0, :)
    d(2, 0, :) = 0.5 * d(2, 0, :)
    d(2, 1, :) = 0.5 * d(2, 1, :)

    d(0, 1, :) = d(1, 0, :)
    d(0, 2, :) = d(2, 0, :)
    d(1, 2, :) = d(2, 1, :)

    divv = d(0, 0, :) + d(1, 1, :) + d(2, 2, :)
    divv = divv / 3.0

    d(0, 0, :) = d(0, 0, :) - divv
    d(1, 1, :) = d(1, 1, :) - divv
    d(2, 2, :) = d(2, 2, :) - divv

  END SUBROUTINE defrate
  
END MODULE KinematicsModule
