! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module hardening_mod3

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use crys_type_mod2
  use hardening_mod4
  use matrix_operations_mod
  use parallel_mod, only: par_quit
  use gather_scatter_mod

  implicit none

  public

contains

  subroutine isotropic_hardening(mesh, crys, phase, func, dfunc, crss, crss_sat, shrate, icode, &
      & mysign, my_phase)

    ! Isotropic hardening assumption

    !---------------------------------------------------------------------------

    ! Arguments:
    ! func: Array of computed hardening rates
    ! dfunc: Array of computed hardening rate derivatives
    ! crss: Array of critical resolved shear stress (hardness)
    ! crss_sat: Array of saturation crss
    ! shrate: Shear rates
    ! icode: Hardening law
    ! mysign:
    ! my_phase:

    type(mesh_type) :: mesh
    type(crys_type) :: crys
    integer, intent(in) :: phase
    real(rk), intent(out) :: func(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(out) :: dfunc(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: crss(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: crss_sat(elt_sub:elt_sup)
    real(rk), intent(in) :: shrate(elt_sub:elt_sup)
    integer, intent(in) :: icode
    real(rk), intent(in) :: mysign(mesh%maxnumslip, elt_sub:elt_sup)
    integer, intent(in) :: my_phase(elt_sub:elt_sup)

    ! Locals:

    real(rk) :: ratio(mesh%maxnumslip, elt_sub:elt_sup)
    integer :: num_ind
    integer :: n_slip
    integer, pointer :: indices(:) => null()

    !---------------------------------------------------------------------------

    n_slip = crys%numslip

    call find_indices(my_phase, phase, indices, num_ind, elt_sub - 1)

    ratio(1, indices) = (crss_sat(indices) - &
        & crss(1, indices))/(crss_sat(indices) - crys%g_0)

    select case (icode)

    case (1) ! eval_rate
      func(1, indices) = shrate(indices) * mysign(1, indices) * crys%h_0 * &
          & ratio(1, indices)**crys%n

      func(2:n_slip, indices) = spread(func(1, indices), dim=1, &
          & ncopies=(n_slip - 2))

    case (2) ! eval_rate_der
      dfunc(1, indices) = -crys%h_0/ &
          & (crss_sat(indices) - crys%g_0)* &
          & shrate(indices)*mysign(1, indices)*crys%n

      dfunc(2:n_slip, indices) = spread(dfunc(1, indices), &
          & dim=1, ncopies=(n_slip - 2))

    case default
    end select

    deallocate (indices)

  end subroutine isotropic_hardening

  !===========================================================================

  subroutine anisotropic_hardening(mesh, crys, phase, sliprate, func, dfunc, crss, crss_sat, icode, &
      & mysign, my_phase)

    ! Anisotropic hardening assumption

    !---------------------------------------------------------------------------

    ! Arguments:
    ! func: Array of computed hardening rates
    ! dfunc: Array of computed hardening rate derivatives
    ! crss: Array of critical resolved shear stress (hardness)
    ! crss_sat: Array of saturation crss
    ! icode: Hardening law
    ! mysign:
    ! my_phase:

    type(mesh_type) :: mesh
    real(rk), intent(in) :: sliprate(mesh%maxnumslip, elt_sub:elt_sup)
    integer, intent(in) :: phase
    type(crys_type) :: crys
    real(rk), intent(out) :: func(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(out) :: dfunc(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: crss(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: crss_sat(elt_sub:elt_sup)
    integer, intent(in) :: icode
    real(rk), intent(in) :: mysign(mesh%maxnumslip, elt_sub:elt_sup)
    integer, intent(in) :: my_phase(elt_sub:elt_sup)

    ! Locals:

    real(rk) :: ratio(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: shrate(mesh%maxnumslip)
    integer :: i

    !---------------------------------------------------------------------------

    shrate = 0.0d0

    ! Proceed where the element is the right phase
    do i = elt_sub, elt_sup
      if (my_phase(i) .eq. phase) then

        ratio(1:crys%numslip, i) = (crss_sat(i) - &
            & crss(1:crys%numslip, i))/(crss_sat(i) - &
            & crys%g_0)

        call calculate_shrate(mesh, sliprate, shrate, crys%structure, i, crys)

        select case (icode)

        case (1) ! eval_rate
          ! All of the slip system func and dfunc values for the
          !   element are calculated in one go

          func(1:crys%numslip, i) = shrate(1:crys%numslip)* &
              & mysign(1:crys%numslip, i)* &
              & crys%h_0*ratio(1:crys%numslip, i) &
              & **crys%n

        case (2) ! eval_rate_der
          dfunc(1:crys%numslip, i) = -crys%h_0/ &
              & (crss_sat(i) - crys%g_0)* &
              & shrate(1:crys%numslip)*mysign(1:crys%numslip, i)* &
              & crys%n

        case default
        end select

      end if
    end do !num_phases

  end subroutine anisotropic_hardening

  !===========================================================================

  subroutine cyclic_hardening(mesh, crys, phase, acmslip, func, dfunc, crss, crss_sat, shear, epseff, &
      & icode, mysign, my_phase)

    ! Cyclic hardening assumption

    !---------------------------------------------------------------------------

    ! Arguments:
    ! func: Array of computed hardening rates
    ! dfunc: Array of computed hardening rate derivatives
    ! crss: Array of critical resolved shear stress (hardness)
    ! crss_sat: Array of satuRATIOn crss
    ! shear: Shear values
    ! epseff:
    ! icode: Hardening law
    ! mysign:
    ! my_phase:

    type(mesh_type) :: mesh
    type(crys_type) :: crys
    integer, intent(in) :: phase
    real(rk), intent(in) :: acmslip(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(out) :: func(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(out) :: dfunc(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: crss(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: crss_sat(elt_sub:elt_sup)
    real(rk), intent(in) :: shear(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: epseff(elt_sub:elt_sup)
    integer, intent(in) :: icode
    real(rk), intent(in) :: mysign(mesh%maxnumslip, elt_sub:elt_sup)
    integer, intent(in) :: my_phase(elt_sub:elt_sup)

    ! Locals:

    real(rk) :: ratio(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: ratio_sh(elt_sub:elt_sup)
    real(rk) :: active_shrate(elt_sub:elt_sup)
    real(rk) :: shear_crit(elt_sub:elt_sup)
    integer :: num_ind
    integer :: n_slip
    integer :: n_gt1
    integer :: n_lt0
    integer, pointer :: indices(:) => null()
    integer :: islip
    real(rk) :: shr_max(elt_sub:elt_sup)
    real(rk) :: shr_min(elt_sub:elt_sup)
    real(rk), parameter :: mach_eps = 2.22d-16

    !---------------------------------------------------------------------------

    active_shrate = 0.0d0
    shear_crit = 0.0d0
    n_lt0 = 0

    n_slip = crys%numslip

    call find_indices(my_phase, phase, indices, num_ind, elt_sub - 1)

    ! rc: All of the slip systems are just looped over instead of just the
    !   one. Since only the isotropic case is being examined here it should
    !   be possible to do only one computation of the slip system and then
    !   assign all the other slip systems that value. Future update should
    !   look into that implementation.

    ratio_sh(indices) = dabs(crss(1, indices)/crss_sat(indices))

    shear_crit(indices) = crys%cyclic_a* &
        & (ratio_sh(indices)**crys%cyclic_c)

    where (shear_crit(indices) .lt. mach_eps)
      shear_crit(indices) = 0.0d0
      !n_lt0 = n_lt0 + 1
    end where

    if (any(dabs(shear_crit) .ge. 1.0d0)) then
      n_gt1 = count((dabs(shear_crit) .ge. 1.0d0))
      !print *, 'Number greater than 1: ', n_gt1
      call par_quit('Error  :     > `shear_crit` is greater than 1.')
    end if

    do islip = 1, n_slip
      where (acmslip(islip, indices) .ge. shear_crit(indices))

        active_shrate(indices) = active_shrate(indices) + &
            & dabs(shear(islip, indices))
      end where
    end do

    shr_min = 1.0d-6*epseff
    shr_max = 1.0d1*epseff

    where (active_shrate .le. shr_min)
      active_shrate = shr_min
    end where

    where (active_shrate .ge. shr_max)
      active_shrate = shr_max
    end where

    do islip = 1, n_slip
      if (islip .eq. 1) then
        ratio(islip, indices) = (crss_sat(indices) - &
            & crss(islip, indices))/(crss_sat(indices) - &
            & crys%g_0)

        select case (icode)

        case (1) ! eval_rate
          func(islip, indices) = active_shrate(indices)* &
              & mysign(islip, indices)*crys%h_0* &
              & ratio(islip, indices)**crys%n

        case (2) ! eval_rate_der
          dfunc(islip, indices) = -crys%h_0/ &
              & (crss_sat(indices) - crys%g_0)* &
              & active_shrate(indices)*mysign(islip, indices)* &
              & crys%n

        case default
        end select

      else
        select case (icode)

        case (1) ! eval_rate
          func(islip, indices) = func(1, indices)

        case (2) ! eval_rate_der
          dfunc(islip, indices) = dfunc(1, indices)

        case default
        end select
      end if
    end do !n_slip

    deallocate (indices)

  end subroutine cyclic_hardening

  !===========================================================================

end module hardening_mod3
