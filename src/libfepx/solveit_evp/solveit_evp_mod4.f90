! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module solveit_evp_mod4

  use general_mod
  use types_mod
  use kinematics_mod_bis
  use stiffness_mod
  use units_mod
  use crys_type_mod
  use solve_stress_evp_mod3

  implicit none

  public

contains

  !> Compute the anisotropic elasto-viscoplastic solution.
  subroutine material_matrix_evp_aniso(mesh, crys, exec, qpt, &
                                     & results, c, c_tan, f, dtime)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys (:)
    type(exec_type), intent(in) :: exec
    type(results_type), intent(inout) :: results
    integer, intent(in) :: qpt
    real(rk), intent(out) :: c(5, 5, elt_sub:elt_sup)
    real(rk), intent(out) :: c_tan(5, 5, elt_sub:elt_sup)
    real(rk), intent(out) :: f(5, elt_sub:elt_sup)
    real(rk), intent(in)  :: dtime

    integer  :: islip, i, j, m_elt, n_slip
    integer, pointer :: indices(:) => null()
    real(rk) :: dtimei, sqr2, sqr32
    real(rk) :: alpha(5)
    real(rk) :: qr5x5(5, 5, elt_sub:elt_sup)
    real(rk) :: e_bar_vec_r(5, elt_sub:elt_sup)
    real(rk) :: e_bar_sm(5, elt_sub:elt_sup)
    real(rk) :: e_vec_sm(5, elt_sub:elt_sup)
    real(rk) :: e_bar(3, 3, elt_sub:elt_sup)
    real(rk) :: e_elas(3, 3, elt_sub:elt_sup)
    real(rk) :: stif_vp(5, 5, elt_sub:elt_sup)
    real(rk) :: stif_evp(5, 5, elt_sub:elt_sup)
    real(rk) :: tan_stif_vp(5, 5, elt_sub:elt_sup)
    real(rk) :: tan_stif_evp(5, 5, elt_sub:elt_sup)
    real(rk) :: keinv_all(5, 5, elt_sub:elt_sup)
    real(rk) :: rss(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: gdot(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: comp(elt_sub:elt_sup)
    real(rk) :: temp(5, 5, elt_sub:elt_sup)
    real(rk) :: v_tensor(3, 3, elt_sub:elt_sup)
    real(rk) :: wp_hat(3, elt_sub:elt_sup)
    real(rk) :: wp_hat_mat(5, 5, elt_sub:elt_sup)
    real(rk) :: wp_x_e(5, elt_sub:elt_sup)
    real(rk) :: f_elas(5, elt_sub:elt_sup)
    real(rk), pointer :: p_hat_vec(:, :) => null()
    integer  :: iphase, num_ind
    integer  :: my_phase(elt_sub:elt_sup)
    real(rk) :: aniso_m_temp(mesh%maxnumslip)
    real(rk) :: aniso_m_min

    !---------------------------------------------------------------------------

    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    m_elt = elt_sup - elt_sub + 1

    tan_stif_evp = 0.0d0
    c_tan = 0.0d0

    ! Initialize factors to correct [c] & {f}

    sqr2 = dsqrt(2.0d0)
    sqr32 = dsqrt(1.5d0)

    alpha(1) = 1.0d0/sqr2
    alpha(2) = sqr32
    alpha(3) = sqr2
    alpha(4) = sqr2
    alpha(5) = sqr2

    ! -------------------------------------------------------------------------

    ! Rotate e_bar_vec to current configuration (lattice axes).
    call rotmat_symm(results%d_rstar(:, :, :, qpt), qr5x5, m_elt)
    call lattice_deform(qr5x5, results%e_bar_vec(:, :, qpt), e_bar_vec_r)

    ! Transform e_bar_vec_r & e_vec to sample axes: {.} = [q] {.}.

    call rotmat_symm(results%ori(:, :, :, qpt), qr5x5, m_elt)

    call mat_x_vec5(qr5x5, e_bar_vec_r, e_bar_sm)

    call mat_x_vec5(qr5x5, results%e_vec(:, :, qpt), e_vec_sm)

    ! Tensor form of e_bar_sm, e_vec_sm: {.} -> [ ].

    call vec_mat_symm(e_bar_sm, e_bar)

    call vec_mat_symm(e_vec_sm, e_elas)

    ! Compute elasto-visco-plastic crystal stiffness and
    !   visco-plastic crystal compliance.

    tan_stif_vp = 0.0d0

    dtimei = 1.0d0/dtime

    do iphase = 1, mesh%num_phases
      call crys_get(crys(iphase), dev=p_hat_vec)

      n_slip = crys(iphase)%numslip

      call find_indices(my_phase, iphase, indices, num_ind, elt_sub - 1)

      do islip = 1, n_slip
        ! Compute the rss
        call ss_project(p_hat_vec(:, islip), results%sig_vec(:, :, qpt), &
            & num_ind, indices, rss(islip, :))

        rss(islip, indices) = rss(islip, indices)/results%crss(islip, indices, qpt)

        where (abs(rss(islip, indices)) .lt. crys(iphase)%t_min)
          rss(islip, indices) = 0.0d0
        end where

        ! Introduced `anisotropic' strain rate sensitivity option here
        ! for hcp material phases - jc

        if (crys(iphase)%use_aniso_m .eqv. .false.) then
          call power_law(rss(islip, :), &
              & crys(iphase)%m, crys(iphase)%gammadot_0, &
              & crys(iphase)%t_min, num_ind, indices, &
              & gdot(islip, :))

          call compliance(rss(islip, :), gdot(islip, :), &
              & results%crss(islip, :, qpt), crys(iphase)%m, crys(iphase)%t_min, &
              & num_ind, indices, comp)

        else if (crys(iphase)%use_aniso_m .eqv. .true.) then
          if (crys(iphase)%structure .eq. "hcp") then
            aniso_m_temp(1:3) = crys(iphase)%aniso_m(1)
            aniso_m_temp(4:6) = crys(iphase)%aniso_m(2)
            aniso_m_temp(7:18) = crys(iphase)%aniso_m(3)

            call power_law(rss(islip, :), &
                & aniso_m_temp(islip), &
                & crys(iphase)%gammadot_0, crys(iphase)%t_min, &
                & num_ind, indices, &
                & gdot(islip, :))

            call compliance(rss(islip, :), gdot(islip, :), &
                & results%crss(islip, :, qpt), aniso_m_temp(islip), crys(iphase)%t_min, &
                & num_ind, indices, comp)

          else if (crys(iphase)%structure .eq. "bct") &
            then
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

            call power_law(rss(islip, :), &
                & aniso_m_temp(islip), &
                & crys(iphase)%gammadot_0, crys(iphase)%t_min, &
                & num_ind, indices, &
                & gdot(islip, :))

            call compliance(rss(islip, :), gdot(islip, :), &
                & results%crss(islip, :, qpt), aniso_m_temp(islip), crys(iphase)%t_min, &
                & num_ind, indices, comp)
          end if
        end if

        ! In order to remove one dimension from my_phase (the dimension over
        ! `ngrain'), we add a loop over that dimension here. - hritz 9/15/05

        ! Reordered loop ordering for better memory striding. - rc 6/24/16
        do j = 1, 5
          do i = 1, 5
            where (my_phase .eq. iphase)
              tan_stif_vp(i, j, :) = tan_stif_vp(i, j, :) &
                  & + comp*p_hat_vec(i, islip)* &
                  & p_hat_vec(j, islip)
            end where
          end do
        end do
      end do

      deallocate (p_hat_vec)
      deallocate (indices)
    end do

    ! Transform tan_stif_vp to sample axes: [s]_sm = [q] [s]_lat [q]'.

    call mat_x_mat5(qr5x5, tan_stif_vp, temp)

    call mat_x_matt5(temp, qr5x5, tan_stif_vp)

    ! Calculate secant moduli

    do iphase = 1, mesh%num_phases
      if (crys(iphase)%use_aniso_m .eqv. .false.) then
        do j = 1, 5
          do i = 1, 5
            where (my_phase .eq. iphase)
              stif_vp(i, j, :) = &
                  & crys(iphase)%m*tan_stif_vp(i, j, :)
            end where
          end do
        end do

      else if (crys(iphase)%use_aniso_m .eqv. .true.) then
        if (crys(iphase)%structure .eq. "hcp") then
          aniso_m_min = minval(crys(iphase)%aniso_m(1:3))

        else if (crys(iphase)%structure .eq. "bct") then
          aniso_m_min = minval(crys(iphase)%aniso_m(:))
        end if

        do j = 1, 5
          do i = 1, 5
            where (my_phase .eq. iphase)
              stif_vp(i, j, :) = &
                  & aniso_m_min*tan_stif_vp(i, j, :)
            end where
          end do
        end do
      end if
    end do

    ! Spread keinv to all crystals and transform it to sample axes.
    ! deb - 6/11/2000

    keinv_all = 0.0d0

    ! In order to remove one dimension from my_phase (the dimension over
    ! `ngrain'), we add a loop over that dimension here. - hritz 9/15/05

    do iphase = 1, mesh%num_phases
      do i = 1, 5
        where (my_phase .eq. iphase)
          keinv_all(i, i, :) = crys(iphase)%keinv(i)
        end where
      end do
    end do

    call mat_x_mat5(qr5x5, keinv_all, temp)

    call mat_x_matt5(temp, qr5x5, keinv_all)

    ! Determinant of tensor v* .

    v_tensor = e_elas

    do i = 1, 3
      v_tensor(i, i, :) = v_tensor(i, i, :) + results%e_elas_kk_bar(:, qpt)/3.0d0 + 1.0d0
    end do

    call determinant(v_tensor, results%detv(:, qpt))

    ! Compute elasto-visco-plastic crystal compliance.

    do j = 1, 5
      do i = 1, 5
        stif_evp(i, j, :) = results%detv(:, qpt)*stif_vp(i, j, :) + &
            & results%detv(:, qpt)*keinv_all(i, j, :)*dtimei
      end do
    end do

    if (exec%itmethod .eq. "NR") then
      do i = 1, 5
        do j = 1, 5
          tan_stif_evp(i, j, :) = results%detv(:, qpt)*tan_stif_vp(i, j, :) &
              & + results%detv(:, qpt)*keinv_all(i, j, :)*dtimei
        end do
      end do
    end if

    ! Compute elasto-visco-plastic crystal stiffness.

    call invert5x5(stif_evp)

    if (exec%itmethod .eq. "NR") call invert5x5(tan_stif_evp)

    ! Force vector due to elastic terms

    do iphase = 1, mesh%num_phases
      call find_indices(mesh%elt_phase(elt_sub:elt_sup), iphase, indices, num_ind, elt_sub - 1)

      call find_wp_hat_phase(crys(iphase), mesh, qpt, dtime, results, e_elas, e_bar, &
                   & gdot, wp_hat, indices)

      deallocate (indices)
    end do

    call wp_hat_mat5x5_all(wp_hat, wp_hat_mat)

    call mat_x_vec5(wp_hat_mat, e_vec_sm, wp_x_e)

    wp_x_e = wp_x_e - 1.0d0/dtime*e_bar_sm

    call mat_x_vec5(stif_evp, wp_x_e, f_elas)

    ! Averaged values of stif_evp, f_elas

    do j = 1, 5
      f(j, :) = f_elas(j, :)

      do i = 1, 5
        c(i, j, :) = stif_evp(i, j, :)
      end do
    end do

    if (exec%itmethod .eq. "NR") then
      do j = 1, 5
        do i = 1, 5
          c_tan(i, j, :) = tan_stif_evp(i, j, :)
        end do
      end do
    end if

    ! Fix [c] & {f} to be consistent with fe equations.

    do j = 1, 5
      f(j, :) = alpha(j)*f(j, :)

      do i = 1, 5
        c(i, j, :) = alpha(i)*c(i, j, :)*alpha(j)
      end do
    end do

    if (exec%itmethod .eq. "NR") then
      do j = 1, 5
        do i = 1, 5
          c_tan(i, j, :) = alpha(i)*c_tan(i, j, :)*alpha(j)
        end do
      end do
    end if

    return

  end subroutine material_matrix_evp_aniso

end module solveit_evp_mod4
