! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module hardening_mod

! Module containing definition of hardening evolution equations

! Contains subroutines:
! hard_law: Evaluate hardening rate or derivative of hardening rate
! isotropic_hardening: Isotropic hardening assumption
! cyclic_hardening: Cyclic hardening assumption
! anisotropic_hardening: Anisotropic hardening assumption

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use crys_type_mod2
  use matrix_operations_mod
  use hardening_mod2
  use parallel_mod, only: par_quit
  use gather_scatter_mod

  implicit none

  public

contains

  subroutine update_hardening(mesh, crys, exec, results_prev, qpt, done, &
                            & results, dtime, euler_method)

    ! Add description here.

    !---------------------------------------------------------------------------

    ! Arguments:

    ! dtime: Timestep
    ! done: Mask of elements not to process

    type(mesh_type) :: mesh
    type(crys_type) :: crys (:)
    type(exec_type), intent(in) :: exec
    type(results_type), intent(in) :: results_prev
    type(results_type), intent(inout) :: results
    integer, intent(in) :: qpt
    character(len = *) :: euler_method
    real(rk), intent(in) :: dtime
    logical, intent(in)    :: done(elt_sub:elt_sup)

    ! Locals:

    integer, pointer :: indices(:) => null()
    integer  :: iter_hard, inewton, iphase, islip
    integer  :: num_ind, n_slip, i, is, num_slip
    real(rk) :: shrate(elt_sub:elt_sup)

    logical  :: newton_ok(mesh%maxnumslip, elt_sub:elt_sup)
    logical  :: newton_ok_all(elt_sub:elt_sup)
    logical  :: converged(mesh%maxnumslip, elt_sub:elt_sup)
    logical  :: converged_all(elt_sub:elt_sup)

    real(rk) :: hard_rate(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: dhard_rate(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: crss_sat(elt_sub:elt_sup)
    real(rk) :: crss_tmp(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: del_crss(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: res_n(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: res(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: ratio_res(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: fj31(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: fjac(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: xlam(mesh%maxnumslip, elt_sub:elt_sup)
    integer  :: my_phase(elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    do iphase = 1, mesh%num_phases

      call find_indices(my_phase, iphase, indices, num_ind, elt_sub - 1)

      num_slip = crys(iphase)%numslip

      shrate(indices) = 0.0d0
      do is = 1, num_slip
        shrate(indices) = shrate(indices) + abs(results%sliprate(is, indices, qpt))
      end do

      if (crys(iphase)%saturation_evolution .eqv. .true.) then
        where (my_phase .eq. iphase)
          where (shrate .gt. 0.0d0)
            crss_sat = crys(iphase)%g_s0 &
                      & *((shrate/crys(iphase)%gammadot_s0) &
                         & **crys(iphase)%m_prime)

          else where
            crss_sat = crys(iphase)%g_s0
          end where
        end where

      else if (crys(iphase)%saturation_evolution .eqv. .false.) then
        where (my_phase .eq. iphase)
          crss_sat = crys(iphase)%g_s0
        end where
      end if

      if (crys(iphase)%h_0 .eq. 0.d0) then
        results%crss(:, indices, qpt) = results_prev%crss(:, indices, qpt)
      end if

      deallocate (indices)
    end do !num_phases

    if (euler_method .eq. 'euler_forward') then
      ! Initial guess via Forward Euler

      call hard_law(mesh, crys, results%acmslip(:, :, qpt), &
          & results%sliprate(:, :, qpt), hard_rate, dhard_rate, &
          & results_prev%crss(:, :, qpt), crss_sat, results%sliprate(:, :, qpt), &
          & shrate, results%defrate_eq(:, qpt), 1)

      results%crss(:, :, qpt) = results_prev%crss(:, :, qpt) + dtime*hard_rate

    else ! if (euler_method .eq. 'euler_backward') then

      ! Start Newton iteration - Backward Euler approx. to hardening law

      ! The newton_ok and converged variables are for each individual slip system.
      ! The newton_ok_all and converged_all variables are for each individual
      ! element. It was done this way to reduce the number of changes needed to
      ! be made to the legacy code.

      newton_ok = .true.
      newton_ok_all = .true.
      converged = .true.
      converged_all = .true.

      ! Initializing all the converged indices slip systems for each individual
      ! phase to false. In a multi-phase system where the different phases have
      ! different number of slip systems. This allows the following case where
      ! line 1 is a system with 3 slip systems and line 2 only has one for the
      ! converged variable:

      ! Line 1: fff
      ! Line 2: ftt

      ! So, in the convergence variable, only the indices that are within the
      ! range of that phases number of slip systems are false.

      do iphase = 1, mesh%num_phases
        n_slip = crys(iphase)%numslip

        do islip = 1, n_slip
          where (my_phase .eq. iphase) converged(islip, :) = .false.
        end do
      end do

      ! Finding where done is true across the element setting converged
      ! true for that entire element

      do islip = 1, mesh%maxnumslip
        where (done) converged(islip, :) = .true.
      end do

      do iter_hard = 1, exec%max_iter_hard_limit
        ! Should we add "where (.not. converged) clause here?
        ! Yes, but we need to pass mask argument to hard_law

        crss_tmp = results%crss(:, :, qpt)

        call hard_law(mesh, crys, results%acmslip(:, :, qpt), &
            & results%sliprate(:, :, qpt), hard_rate, dhard_rate, &
            & crss_tmp, crss_sat, results%sliprate(:, :, qpt), &
            & shrate, results%defrate_eq(:, qpt), 1)

        res_n = -(crss_tmp - results_prev%crss(:, :, qpt) - dtime*hard_rate)

        ! Patch to avoid divide-by-zero exception

        where (abs(res_n) .eq. 0.0d0)
          converged = .true.
        end where

        call hard_law(mesh, crys, results%acmslip(:, :, qpt), &
            & results%sliprate(:, :, qpt), hard_rate, dhard_rate, &
            & crss_tmp, crss_sat, results%sliprate(:, :, qpt), &
            & shrate, results%defrate_eq(:, qpt), 2)

        fjac = 1.0 - dtime*dhard_rate
        del_crss = res_n/fjac

        ! Find where the newton_ok criterion is not satisfied

        where ((fjac .lt. vtiny) .and. (.not. converged)) newton_ok = .false.
        xlam = crss_tmp + del_crss

        call hard_law(mesh, crys, results%acmslip(:, :, qpt), &
            & results%sliprate(:, :, qpt), hard_rate, dhard_rate, &
            & xlam, crss_sat, results%sliprate(:, :, qpt), &
            & shrate, results%defrate_eq(:, qpt), 1)

        res = -(xlam - results_prev%crss(:, :, qpt) - dtime*hard_rate)

        where ((.not. converged))
          ratio_res = abs(res)/abs(res_n)

        else where
          ratio_res = 0.0d0
        end where

        ! Line search

        fj31 = 1.0d0

        do while (any(ratio_res .gt. 1.0d0 &
            & .and. newton_ok .and. .not. converged))

          where (ratio_res .gt. 1.0d0 &
              & .and. newton_ok .and. .not. converged) fj31 = fj31*0.5d0

          ! Find where the newton_ok criterion is not satisfied

          where (fj31 .lt. exec%ls_cutoff) newton_ok = .false.
          xlam = crss_tmp + fj31*del_crss

          call hard_law(mesh, crys, results%acmslip(:, :, qpt), &
              & results%sliprate(:, :, qpt), hard_rate, dhard_rate, &
              & xlam, crss_sat, results%sliprate(:, :, qpt), &
              & shrate, results%defrate_eq(:, qpt), 1)

          where (ratio_res .gt. 1.0d0 .and. newton_ok .and. .not. converged)
            res = -(xlam - results_prev%crss(:, :, qpt) - dtime*hard_rate)
            ratio_res = abs(res)/abs(res_n)
          end where
        end do

        where ((newton_ok) .and. (.not. converged)) &
            & results%crss(:, :, qpt) = crss_tmp + fj31*del_crss

        ! Find slip systems that have reached convergence tolerance

        do iphase = 1, mesh%num_phases
          n_slip = crys(iphase)%numslip

          do islip = 1, n_slip
            where (my_phase .eq. iphase)
              where ((dabs(res(islip, :)) .lt. &
                  & exec%toler_hard*crys(iphase)%g_0) .and. &
                  & (newton_ok(islip, :)) .and. (.not. done)) &
                  & converged(islip, :) = .true.
            end where
          end do
        end do

        ! Finds out which areas have not converged and which
        !   newton_ok are not okay

        do iphase = 1, mesh%num_phases
          n_slip = crys(iphase)%numslip

          do islip = 1, n_slip
            ! Converged_all set to false only when one slip system is false

            where (.not. converged(islip, :)) converged_all = .false.
            ! Newton_ok_all set to false only when one slip system is false

            where (.not. newton_ok(islip, :)) newton_ok_all = .false.
          end do

          ! If all slip systems for a given element are converged, then
          ! update the per-element arrays as needed. We limit the first index
          ! to the slip systems available for the given phase as above.

          do i = elt_sub, elt_sup
            if (all(converged(1:n_slip, i))) converged_all(i) = .true.
            if (all(newton_ok(1:n_slip, i))) newton_ok_all(i) = .true.
          end do
        end do

        inewton = count(.not. newton_ok_all)
        if ((count(converged_all) + inewton) .eq. elt_sup - elt_sub + 1) exit
      end do

      if (iter_hard > exec%max_iter_hard_limit) then
        ! All processes with this condition will write this message. This could
        ! likely be cleaned up to avoid

        write (*, '(a)') 'Warning:       . Update hardening iteration limit reached.'
      end if

      ! Perhaps a count of "not converged" is a better approach?
      where (.not. converged) results%crss(:, :, qpt) = results_prev%crss(:, :, qpt)

    end if

  end subroutine update_hardening

end module hardening_mod
