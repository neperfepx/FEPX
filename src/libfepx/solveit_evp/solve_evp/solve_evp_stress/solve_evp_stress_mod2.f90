! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module solve_stress_evp_mod2

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use kinematics_mod
  use matrix_operations_mod
  use solve_stress_evp_mod3
  use solve_evp_vpstress_mod
  use units_mod
  use parallel_mod
  use crys_type_mod
  use crys_type_mod2

  implicit none

  public

contains

  integer function solve_stress_evp_newton(mesh, crys, exec, results, qpt, &
                                         & dt, converged, done) result(irc)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys (:)
    type(results_type), intent(inout) :: results
    integer, intent(in) :: qpt
    type(exec_type), intent(in) :: exec
    logical, intent(inout) :: converged(elt_sub:elt_sup)
    logical, intent(in)    :: done(elt_sub:elt_sup)
    real(rk), intent(in)    :: dt

    real(rk) :: wp_hat_matx(5, 5, elt_sub:elt_sup)
    integer  :: iter_newton, islip, i, inewton, n_slip
    real(rk) :: c1, xn(mesh%num_phases), xnn(mesh%num_phases)
    real(rk) :: sig_0(5, elt_sub:elt_sup)
    real(rk) :: del_s(5, elt_sub:elt_sup)
    real(rk) :: fact(elt_sub:elt_sup)
    real(rk) :: fjac(5, 5, elt_sub:elt_sup)
    real(rk) :: fj21(elt_sub:elt_sup)
    real(rk) :: fj(5, elt_sub:elt_sup)
    logical  :: newton_ok(elt_sub:elt_sup)
    real(rk) :: dgdot(mesh%maxnumslip, elt_sub:elt_sup)
    integer :: jiter(elt_sub:elt_sup)
    real(rk) :: wp_x_e(5, elt_sub:elt_sup)
    real(rk) :: res(elt_sub:elt_sup), res_n(elt_sub:elt_sup)
    real(rk) :: xlambda(5, elt_sub:elt_sup)
    real(rk) :: rss_abs(elt_sub:elt_sup)
    real(rk) :: ratio_res(elt_sub:elt_sup)
    real(rk) :: res_aux(elt_sub:elt_sup)

    real(rk), pointer :: p_hat_vec(:, :) => null()
    real(rk), pointer :: ppt(:, :, :) => null()
    real(rk), pointer :: tmp2d(:, :) => null()
    real(rk), pointer :: tmp3d(:, :, :) => null()

    real(rk), parameter :: one_dp = 1.0d0
    real(rk) :: uflow = tiny(one_dp)
    real(rk) :: toosmall(mesh%num_phases)

    integer  :: iphase, num_ind
    integer, pointer :: indices(:) => null()

    real(rk) :: aniso_m_temp(mesh%maxnumslip)
    real(rk) :: axnn(mesh%maxnumslip, mesh%num_phases)
    real(rk) :: axn(mesh%maxnumslip, mesh%num_phases)
    real(rk) :: atoosmall(mesh%maxnumslip, mesh%num_phases)

    !---------------------------------------------------------------------------

    real(rk) :: qr5x5(5, 5, elt_sub:elt_sup)
    real(rk) :: qr3x3(3, 3, elt_sub:elt_sup)
    real(rk) :: d_vec_lat(5, elt_sub:elt_sup)
    real(rk) :: e_bar_vec_r(5, elt_sub:elt_sup)
    real(rk) :: e_elas(3, 3, elt_sub:elt_sup)
    real(rk) :: e_bar(3, 3, elt_sub:elt_sup)

    ! results%ori(:, :, :, qpt) [3x3] --> qr5x5 [5x5]
    ! d_vec(sample coords) --> d_vec_lat(crystal coords)
    call rotmat_symm(results%ori(:, :, :, qpt), qr5x5, elt_sup - elt_sub + 1)
    call lattice_deform(qr5x5, results%d_vec(:, :, qpt), d_vec_lat)
    call rotmat_skew(results%ori(:, :, :, qpt), qr3x3)
    call lattice_spin(qr3x3, results%w_vec(:, :, qpt), results%w_vec_lat(:, :, qpt))

    ! d_rstar [3x3] --> qr5x5 [5x5]
    ! Apply dR* to e_bar_vec --> e_bar_vec_r
    call rotmat_symm(results%d_rstar(:, :, :, qpt), qr5x5, elt_sup - elt_sub + 1)
    call lattice_deform(qr5x5, results%e_bar_vec(:, :, qpt), e_bar_vec_r)
    call vec_mat_symm(e_bar_vec_r, e_bar)

    irc = 0
    jiter = 0
    sig_0 = 0.0d0
    newton_ok = .true.

    c1 = 1.0d0/dt

    dgdot = 0.0d0

    ! Begin iterations

    do iter_newton = 1, exec%sx_max_iters_newton
      sig_0 = results%sig_vec(:, :, qpt)

      ! Elastic strains

      do iphase = 1, mesh%num_phases
        call find_indices(mesh%elt_phase(elt_sub:elt_sup), iphase, indices, num_ind, elt_sub - 1)
        call crys_get(crys(iphase), dev=p_hat_vec, pptrans=ppt)
        n_slip = crys(iphase)%numslip

        allocate (tmp2d(5, num_ind))
        call vec_d_vec5(crys(iphase)%keinv, results%sig_vec(:, indices, qpt), tmp2d, &
            & num_ind)
        results%e_vec(:, indices, qpt) = tmp2d
        deallocate (tmp2d)

        allocate (tmp3d(3, 3, num_ind))
        call vec_mat_symm(results%e_vec(:, indices, qpt), tmp3d)
        e_elas(:, :, indices) = tmp3d
        deallocate (tmp3d)

        ! Power law viscoplastic model

        ! Added parameter check for anisotropic rate sensitivity
        if (crys(iphase)%use_aniso_m .eqv. .false.) then
          xnn(iphase) = 1.0d0/crys(iphase)%m
          xn(iphase) = xnn(iphase) - 1.0d0
          toosmall(iphase) = uflow**crys(iphase)%m

        else if (crys(iphase)%use_aniso_m .eqv. .true.) then
          if (crys(iphase)%structure .eq. "hcp") then
            aniso_m_temp(1:3) = crys(iphase)%aniso_m(1)
            aniso_m_temp(4:6) = crys(iphase)%aniso_m(2)
            aniso_m_temp(7:18) = crys(iphase)%aniso_m(3)

            axnn(:, iphase) = 1.0d0/aniso_m_temp(:)
            axn(:, iphase) = axnn(:, iphase) - 1.0d0
            atoosmall(:, iphase) = uflow**aniso_m_temp(:)

          else if (crys(iphase)%structure .eq. "bct") &
              & then
            aniso_m_temp(1:2) = crys(iphase)%aniso_m(1)
            aniso_m_temp(3:4) = crys(iphase)%aniso_m(2)
            aniso_m_temp(5:6) = crys(iphase)%aniso_m(3)
            aniso_m_temp(7:10) = crys(iphase)%aniso_m(4)
            aniso_m_temp(11:12) = crys(iphase)%aniso_m(5)
            aniso_m_temp(13:16) = crys(iphase)%aniso_m(6)
            aniso_m_temp(17:18) = crys(iphase)%aniso_m(7)
            aniso_m_temp(19:20) = crys(iphase)%aniso_m(8)
            aniso_m_temp(21:24) = crys(iphase)%aniso_m(9)
            aniso_m_temp(25:32) = crys(iphase)%aniso_m(10)

            axnn(:, iphase) = 1.0d0/aniso_m_temp(:)
            axn(:, iphase) = axnn(:, iphase) - 1.0d0
            atoosmall(:, iphase) = uflow**aniso_m_temp(:)
          end if
        end if

        do islip = 1, n_slip
          results%rss(islip,indices,qpt) = 0.0d0

          do i = 1, 5
            results%rss(islip,indices,qpt) = results%rss(islip,indices,qpt) + &
                & p_hat_vec(i, islip)*results%sig_vec(i, indices, qpt)/ &
                & results%crss(islip, indices, qpt)
          end do

          rss_abs(indices) = dabs(results%rss(islip,indices,qpt))

          if (crys(iphase)%use_aniso_m .eqv. .false.) then
            where (rss_abs(indices) .le. toosmall(iphase))
              rss_abs(indices) = 0.0d0
            end where

            fj21(indices) = crys(iphase)%gammadot_0 * rss_abs(indices)**xn(iphase)

            dgdot(islip, indices) = fj21(indices)*xnn(iphase) / results%crss(islip, indices, qpt)

            results%sliprate(islip, indices, qpt) = fj21(indices)*results%rss(islip,indices,qpt)

          else if (crys(iphase)%use_aniso_m .eqv. .true.) then
            where (rss_abs(indices) .le. atoosmall(islip, iphase)) &
                & rss_abs(indices) = 0.0d0

            fj21(indices) = crys(iphase)%gammadot_0 * rss_abs(indices)**axn(islip, iphase)

            dgdot(islip, indices) = fj21(indices) &
                                & * axnn(islip, iphase)/results%crss(islip, indices, qpt)

            results%sliprate(islip, indices, qpt) = fj21(indices)*results%rss(islip,indices,qpt)
          end if
        end do !n_slip

        do islip = 1, n_slip
          results%rss(islip,indices,qpt) = results%rss(islip,indices,qpt) &
                                       & * results%crss(islip,indices,qpt)
        end do

        ! Set up the Jacobian

        call dp_wp_hat_phase(crys(iphase), mesh, qpt, dt, &
                           & results, e_elas, e_bar, results%sliprate(:, :, qpt), &
                           & results%wp_hat(:, :, qpt), &
                           & num_ind, indices)

        call wp_hat_mat5x5(results%wp_hat(:, :, qpt), wp_hat_matx, num_ind, indices)

        allocate (tmp2d(5, num_ind))
        call mat_x_vec5(wp_hat_matx(:, :, indices), results%e_vec(:, indices, qpt), tmp2d)
        wp_x_e(:, indices) = tmp2d
        deallocate (tmp2d)

        allocate (tmp3d(5, 5, num_ind))
        call matrix_fjac(mesh, tmp3d, crys(iphase)%keinv, &
                       & wp_hat_matx(:, :, indices), c1, dgdot(:, indices), &
                       & ppt, n_slip)
        fjac(:, :, indices) = tmp3d
        deallocate (tmp3d)

        ! Set up the system function (rhs)

        call residual(res_n, del_s, d_vec_lat, results%e_vec(:, :, qpt), e_bar_vec_r, &
                    & results%dp_hat(:, :, qpt), wp_x_e, c1, num_ind, indices)

        deallocate (indices)
        deallocate (ppt)
        deallocate (p_hat_vec)
      end do !mesh%num_phases

      ! Compute new iteration

      call symmetrize_jac(fjac, del_s)
      call solvit(fjac, del_s)
      call check_diagonals_evp(fjac, newton_ok, done)

      ! Trial residual

      xlambda = sig_0 + del_s

      do iphase = 1, mesh%num_phases
        call find_indices(mesh%elt_phase(elt_sub:elt_sup), iphase, indices, num_ind, elt_sub - 1)
        call crys_get(crys(iphase), dev=p_hat_vec)
        n_slip = crys(iphase)%numslip

        allocate (tmp2d(5, num_ind))
        call vec_d_vec5(crys(iphase)%keinv, xlambda(:, indices), &
            & tmp2d, num_ind)
        results%e_vec(:, indices, qpt) = tmp2d
        deallocate (tmp2d)

        allocate (tmp3d(3, 3, num_ind))
        call vec_mat_symm(results%e_vec(:, indices, qpt), tmp3d)
        e_elas(:, :, indices) = tmp3d
        deallocate (tmp3d)

        ! Added parameter check for anisotropic rate sensitivity
        do islip = 1, n_slip
          results%rss(islip,indices,qpt) = 0.0d0

          do i = 1, 5
            results%rss(islip,indices,qpt) = results%rss(islip,indices,qpt) + &
                & p_hat_vec(i, islip)*xlambda(i, indices) &
                & /results%crss(islip, indices, qpt)
          end do

          rss_abs(indices) = dabs(results%rss(islip,indices,qpt))

          if (crys(iphase)%use_aniso_m .eqv. .false.) then
            where (rss_abs(indices) .le. toosmall(iphase)) &
                & rss_abs(indices) = 0.0d0
            results%sliprate(islip, :, qpt) = crys(iphase)%gammadot_0* &
                & results%rss(islip, :, qpt)*rss_abs**xn(iphase)

          else if (crys(iphase)%use_aniso_m .eqv. .true.) then
            where (rss_abs(indices) .le. atoosmall(islip, iphase))&
                & rss_abs(indices) = 0.0d0
            results%sliprate(islip, :, qpt) = crys(iphase)%gammadot_0* &
                & results%rss(islip, :, qpt)*rss_abs**axn(islip, iphase)
          end if
        end do !n_slip

        do islip = 1, n_slip
          results%rss(islip,indices,qpt) = results%rss(islip,indices,qpt) &
                                       & * results%crss(islip,indices,qpt)
        end do

        call dp_wp_hat_phase(crys(iphase), mesh, qpt, dt, &
                           & results, e_elas, e_bar, results%sliprate(:, :, qpt), &
                           & results%wp_hat(:, :, qpt), &
                           & num_ind, indices)

        call wp_hat_mat5x5(results%wp_hat(:, :, qpt), wp_hat_matx, num_ind, indices)

        allocate (tmp2d(5, num_ind))
        call mat_x_vec5(wp_hat_matx(:, :, indices), results%e_vec(:, indices, qpt), tmp2d)
        wp_x_e(:, indices) = tmp2d
        deallocate (tmp2d)

        call residual(res, fj, d_vec_lat, results%e_vec(:, :, qpt), e_bar_vec_r, &
                    & results%dp_hat(:, :, qpt), wp_x_e, c1, num_ind, indices)

        deallocate (indices)
        deallocate (p_hat_vec)
      end do !num_phases

      ! Line Search.

      fact = 1.0d0
      ratio_res = res/res_n

      do while (any(ratio_res .gt. 1.0 .and. newton_ok .and. .not. converged))
        where (ratio_res .gt. 1.0 .and. newton_ok .and. .not. converged) &
            & fact = fact*0.5d0

        if (any(fact .lt. 0.001)) then
          write (*, '(a,i0,i0)') 'Warning:       . Error in line &
              &search'
          !write(*, '(a,i0,i0)') 'Warning:       . Error in line &
          !    &search ', count(fact .lt. 0.1), iter_newton
          where (fact .lt. 0.001) newton_ok = .false.
        end if

        xlambda(1, :) = sig_0(1, :) + fact*del_s(1, :)
        xlambda(2, :) = sig_0(2, :) + fact*del_s(2, :)
        xlambda(3, :) = sig_0(3, :) + fact*del_s(3, :)
        xlambda(4, :) = sig_0(4, :) + fact*del_s(4, :)
        xlambda(5, :) = sig_0(5, :) + fact*del_s(5, :)

        do iphase = 1, mesh%num_phases
          call find_indices(mesh%elt_phase(elt_sub:elt_sup), iphase, indices, num_ind, elt_sub - 1)
          call crys_get(crys(iphase), dev=p_hat_vec)
          n_slip = crys(iphase)%numslip

          allocate (tmp2d(5, num_ind))
          call vec_d_vec5(crys(iphase)%keinv, xlambda(:, indices), tmp2d, num_ind)
          results%e_vec(:, indices, qpt) = tmp2d
          deallocate (tmp2d)

          allocate (tmp3d(3, 3, num_ind))
          call vec_mat_symm(results%e_vec(:, indices, qpt), tmp3d)
          e_elas(:, :, indices) = tmp3d
          deallocate (tmp3d)

          do islip = 1, n_slip
            results%rss(islip, indices, qpt) = 0.0d0

            do i = 1, 5
              results%rss(islip, indices, qpt) = results%rss(islip, indices, qpt) + &
                  & p_hat_vec(i, islip)* &
                  & xlambda(i, indices)/results%crss(islip, indices, qpt)
            end do

            rss_abs(indices) = dabs(results%rss(islip, indices, qpt))

            ! Added parameter check for anisotropic rate sensitivity
            if (crys(iphase)%use_aniso_m .eqv. .false.) then
              where (rss_abs(indices) .le. toosmall(iphase)) &
                  & rss_abs(indices) = 0.0d0
              results%sliprate(islip, :, qpt) = crys(iphase)%gammadot_0*&
                  &results%rss(islip, :, qpt)*rss_abs**xn(iphase)

            else if (crys(iphase)%use_aniso_m .eqv. .true.) then
              where (rss_abs(indices) .le. atoosmall(islip, iphase)) &
                  & rss_abs(indices) = 0.0d0
              results%sliprate(islip, :, qpt) = crys(iphase)%gammadot_0* &
                  &results%rss(islip, :, qpt)*rss_abs**axn(islip, iphase)
            end if
          end do !n_slip

          do islip = 1, n_slip
            results%rss(islip,indices,qpt) = results%rss(islip,indices,qpt) &
                                         & * results%crss(islip,indices,qpt)
          end do

          call dp_wp_hat_phase(crys(iphase), mesh, qpt, dt, &
                             & results, e_elas, e_bar, results%sliprate(:, :, qpt), &
                             & results%wp_hat(:, :, qpt), &
                             & num_ind, indices)

          call wp_hat_mat5x5(results%wp_hat(:, :, qpt), wp_hat_matx, num_ind, indices)

          allocate (tmp2d(5, num_ind))
          call mat_x_vec5(wp_hat_matx(:, :, indices), results%e_vec(:, indices, qpt), tmp2d)
          wp_x_e(:, indices) = tmp2d
          deallocate (tmp2d)

          call residual(res_aux, fj, d_vec_lat, results%e_vec(:, :, qpt), e_bar_vec_r, &
                      & results%dp_hat(:, :, qpt), wp_x_e, c1, num_ind, indices)

          deallocate (indices)
          deallocate (p_hat_vec)
        end do !num_phases

        where ((ratio_res .gt. 1.0) .and. (newton_ok) .and. &
            & (.not. converged))

          res = res_aux
          ratio_res = res/res_n
        end where
      end do !do while

      ! Update stresses and check convergence.

      where ((newton_ok) .and. (.not. converged))
        results%sig_vec(1, :, qpt) = sig_0(1, :) + fact*del_s(1, :)
        results%sig_vec(2, :, qpt) = sig_0(2, :) + fact*del_s(2, :)
        results%sig_vec(3, :, qpt) = sig_0(3, :) + fact*del_s(3, :)
        results%sig_vec(4, :, qpt) = sig_0(4, :) + fact*del_s(4, :)
        results%sig_vec(5, :, qpt) = sig_0(5, :) + fact*del_s(5, :)
      end where

      where ((res .le. exec%sx_tol) .and. (newton_ok) .and. (.not. done)&
          & .and. (.not. converged)) jiter = iter_newton
      where ((res .le. exec%sx_tol) .and. (newton_ok) .and. (.not. done))&
          & converged = .true.

      ! Return if all grains have converged.

      inewton = count(.not. newton_ok)

      if ((count(converged) + inewton) .eq. size(converged)) then
        if (inewton .gt. 0) then
          write (*, '(a)') 'Warning:       . solve_stress_evp_newton: &
              &Some elements did not converge'
          !write(*, '(a,i0,a,i0)') 'Warning:       . &
          !    &solve_stress_evp_newton: Converged = ', count(converged), &
          !    & ' Remaining = ', inewton
          !write(*, '(a,g12.5,2x,g12.5)') '                 Residual &
          !    &for converged grains = ', minval(res, mask = converged), &
          !    & maxval(res, mask = converged)
        end if

        irc = inewton

        return
      end if
    end do

    irc = -2

    return

  end function solve_stress_evp_newton

end module solve_stress_evp_mod2
