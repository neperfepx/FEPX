! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module solve_evp_rstar_mod2

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use hardening_mod
  use solve_evp_vpstress_mod
  use crys_type_mod
  use crys_type_mod2
  use matrix_operations_mod
  use orientation_conversion_mod

  implicit none

  public

contains

  !>
  subroutine compute_sliprate(mesh, crys, qpt, results)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys (:)
    type(results_type), intent(inout) :: results
    integer, intent(in) :: qpt

    integer :: num_slip
    integer :: phase, is, num_ind
    real(rk), pointer :: plocal(:, :) => null()
    integer, pointer :: indices(:) => null()
    real(rk) :: rss(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: aniso_m_temp(mesh%maxnumslip)
    integer  :: my_phase(elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    do phase = 1, mesh%num_phases
      call crys_get(crys(phase), dev=plocal)

      num_slip = crys(phase)%numslip

      ! Compute slip rate on each slip system

      call find_indices(my_phase, phase, indices, num_ind, elt_sub - 1)

      do is = 1, num_slip
        call ss_project(plocal(:, is), results%sig_vec(:, :, qpt), &
            & num_ind, indices, rss(is, :))
        rss(is, indices) = rss(is, indices)/results%crss(is, indices, qpt)

        where (abs(rss(is, indices)) .lt. crys(phase)%t_min)
          rss(is, indices) = 0.0d0
        end where

        if (crys(phase)%use_aniso_m .eqv. .false.) then
          call power_law(rss(is, :), &
                       & crys(phase)%m, crys(phase)%gammadot_0, &
                       & crys(phase)%t_min, num_ind, indices, &
                       & results%sliprate(is, :, qpt))

        else if (crys(phase)%use_aniso_m .eqv. .true.) then
          if (crys(phase)%structure .eq. "hcp") then
            aniso_m_temp(1:3) = crys(phase)%aniso_m(1)
            aniso_m_temp(4:6) = crys(phase)%aniso_m(2)
            aniso_m_temp(7:18) = crys(phase)%aniso_m(3)

            call power_law(rss(is, :), &
                         & aniso_m_temp(is), crys(phase)%gammadot_0, &
                         & crys(phase)%t_min, num_ind, indices, &
                         & results%sliprate(is, :, qpt))

          else if (crys(phase)%structure .eq. "bct") &
              & then
            aniso_m_temp(1:2) = crys(phase)%aniso_m(1)
            aniso_m_temp(3:4) = crys(phase)%aniso_m(2)
            aniso_m_temp(5:6) = crys(phase)%aniso_m(3)
            aniso_m_temp(7:10) = crys(phase)%aniso_m(4)
            aniso_m_temp(11:12) = crys(phase)%aniso_m(5)
            aniso_m_temp(13:16) = crys(phase)%aniso_m(6)
            aniso_m_temp(17:18) = crys(phase)%aniso_m(7)
            aniso_m_temp(19:20) = crys(phase)%aniso_m(8)
            aniso_m_temp(21:24) = crys(phase)%aniso_m(9)
            aniso_m_temp(25:32) = crys(phase)%aniso_m(10)

            call power_law(rss(is, :), &
                         & aniso_m_temp(is), crys(phase)%gammadot_0, &
                         & crys(phase)%t_min, num_ind, indices, &
                         & results%sliprate(is, :, qpt))
          end if
        end if
      end do !num_slip

      ! should probably be moved to write_res
      where (abs(results%sliprate(1:num_slip, indices, qpt)) .le. 1d-30)
        results%sliprate(1:num_slip, indices, qpt) = 0.0d0
      end where

      ! Compute plastic spin

      deallocate (indices)
      deallocate (plocal)
    end do !num_phases


  end subroutine compute_sliprate

  !===========================================================================

  !> Integrate the evolution equation for r* (modified rq, 02/2023)
  subroutine update_rstar(crys, mesh, results_prev, qpt, done, dtime, results)

    type(crys_type), intent(in) :: crys(:)
    type(mesh_type), intent(in) :: mesh
    type(results_type), intent(in) :: results_prev
    integer,  intent(in) :: qpt
    logical,  intent(in) :: done(elt_sub:elt_sup)
    real(rk), intent(in) :: dtime
    type(results_type), intent(inout) :: results

    integer  :: i, j, k, l, is, phase, num_ind, num_slip
    real(rk) :: rotvec(3)
    real(rk), pointer :: qlocal(:, :) => null()
    integer, pointer :: indices(:) => null()
    integer  :: my_phase(elt_sub:elt_sup)

    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    ! compute results%wp_ss from slip rates, expressed in crystal frame

    results%wp_ss(:, :, qpt) = 0.0d0

    do phase = 1, mesh%num_phases
      call crys_get(crys(phase), skw=qlocal)

      call find_indices(my_phase, phase, indices, num_ind, elt_sub - 1)

      num_slip = crys(phase)%numslip

      do is = 1, num_slip
        do i = 1, 3
          results%wp_ss(i, indices, qpt) = results%wp_ss(i, indices, qpt) + &
              & qlocal(i, is)*results%sliprate(is, indices, qpt)
        end do
      end do !num_slip

      deallocate (qlocal)
    end do

    ! do not compute results%wp_hat (already done)

    ! compute rstar from results%wp_ss and results%wp_hat

    do i = elt_sub, elt_sup
      if (.not. done(i)) then

        ! Compute the rotation vector (rotvec = axis x angle) from wp_hat and wp_ss
        ! Note: we could simply use rotvec = (wp_hat - wp_ss) * dtime, but
        !       wp_hat and wp_ss are not in the right order
        rotvec(1) = (results%wp_hat(3, i, qpt) - results%wp_ss(3, i, qpt)) * dtime
        rotvec(2) = (results%wp_hat(2, i, qpt) - results%wp_ss(2, i, qpt)) * dtime * (-1.0d0)
        rotvec(3) = (results%wp_hat(1, i, qpt) - results%wp_ss(1, i, qpt)) * dtime

        ! Compute the rotation matrix from the rotation vector
        ! The rotation matrix must be in passive convention, so we pass the opposite
        ! or the rotation vector below
        call rotvec_to_rotmat_(-rotvec, results%d_rstar(:, :, i, qpt))

        if (allocated(results%d_rstar_spin)) then
          ! repeat for wp_hat alone, and record as d_rstar_spin
          rotvec(1) = results%wp_hat(3, i, qpt) * dtime
          rotvec(2) = results%wp_hat(2, i, qpt) * dtime * (-1.0d0)
          rotvec(3) = results%wp_hat(1, i, qpt) * dtime
          call rotvec_to_rotmat_(-rotvec, results%d_rstar_spin(:, :, i, qpt))
        end if

        if (allocated(results%d_rstar_slip)) then
          ! repeat for wp_ss alone, and record as d_rstar_slip
          rotvec(1) = results%wp_ss(3, i, qpt) * dtime
          rotvec(2) = results%wp_ss(2, i, qpt) * dtime * (-1.0d0)
          rotvec(3) = results%wp_ss(1, i, qpt) * dtime
          call rotvec_to_rotmat_(-rotvec, results%d_rstar_slip(:, :, i, qpt))
        end if

        results%rstar(:, :, i, qpt) = 0.0d0
        do k = 1, 3
          do l = 1, 3
            do j = 1, 3
              results%rstar(j, k, i, qpt) = results%rstar(j, k, i, qpt) &
                                   & + results_prev%rstar(j, l, i, qpt) &
                                   & * results%d_rstar(l, k, i, qpt)
            end do
          end do
        end do
      end if
    end do

  end subroutine update_rstar

  !=============================================================================

  !> Compute new orientation: c = c_0 * rstar
  ! rq 04/2023: we now compute from the previous orientation (using d_rstar)
  ! instead of from the mesh (using rstar), which is equivalent but feels more
  ! natural (incremental, as the rest of the code)
  subroutine update_ori(results_prev, qpt, done, results)

    logical,  intent(in)  :: done(elt_sub:elt_sup)
    integer, intent(in) :: qpt
    type(results_type), intent(in) :: results_prev
    type(results_type), intent(inout) :: results

    integer  :: i, j, k, l

    ! --------------------------------------------------------------------------

    do i = elt_sub, elt_sup
      if (.not. done(i)) then

        results%ori(:, :, i, qpt) = 0.0d0
        do k = 1, 3
          do l = 1, 3
            do j = 1, 3
              results%ori(j, k, i, qpt) = results%ori(j, k, i, qpt) &
                                      & + results_prev%ori(j, l, i, qpt) &
                                      & * results%d_rstar(l, k, i, qpt)
            end do
          end do
        end do

      end if
    end do

  end subroutine update_ori

end module solve_evp_rstar_mod2
