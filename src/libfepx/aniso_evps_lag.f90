! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE AnisoEvpsModule

  use units_mod
  USE DimsModule

  use microstructure_mod
  USE READ_INPUT_MOD
  USE UtilsCrystalModule

  USE StressSolveVpModule
  use MaterialMatrixVpModule
  USE StressSolveEvpsModule
  !
  IMPLICIT  NONE
  !
  PRIVATE
  PUBLIC :: aniso_evps

CONTAINS

SUBROUTINE aniso_evps(&
     &   c, c_tan, f, detv, c_angs, sig_vec, crss, rstar_n, rstar, e_bar_vec,&
     &   wts, w_vec, e_elas_kk, Keinv, dtime, NR )
  !
  !----------------------------------------------------------------------
  !
  !     Arguments:
  !
  REAL(RK), INTENT(OUT) ::  c(TVEC, TVEC, el_sub1:el_sup1)
  REAL(RK), INTENT(OUT) ::  c_tan(TVEC, TVEC, el_sub1:el_sup1)
  REAL(RK), INTENT(OUT) ::  f(TVEC, el_sub1:el_sup1)
  REAL(RK), INTENT(OUT) ::  detv(el_sub1:el_sup1)
  REAL(RK), INTENT(IN)  ::  dtime
  REAL(RK), INTENT(IN)  ::  keinv(0:TVEC1,1:numphases)
  REAL(RK), INTENT(IN)  ::  c_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK), INTENT(IN)  ::  sig_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK), INTENT(IN)  ::  crss(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK), INTENT(IN)  ::  rstar_n(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK), INTENT(IN)  ::  rstar(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK), INTENT(IN)  ::  e_bar_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK), INTENT(IN)  ::  wts(0:ngrain1, el_sub1:el_sup1)
  REAL(RK), INTENT(IN)  ::  w_vec(0:DIMS1, el_sub1:el_sup1)
  REAL(RK), INTENT(IN)  ::  e_elas_kk(el_sub1:el_sup1)
  LOGICAL, INTENT(IN)   ::  NR
  !
  !     Locals:
  !
  INTEGER  ::  islip, i, j, m_el, k, n_slip
  INTEGER, pointer ::  indices(:) => NULL()
  !
  REAL(RK)  ::  dtimei, sqr2, sqr32
  REAL(RK)  ::  alpha(TVEC)
  REAL(RK)  ::  e_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  d_rstar(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  qr5x5(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  e_bar_vec_r(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  e_bar_sm(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  e_vec_sm(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  e_bar(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  e_elas(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  w_vec_grn(0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
  !             
  REAL(RK)  ::  stif_vp(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  stif_evp(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  tan_stif_vp(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  tan_stif_evp(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  keinv_all(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  !             
  REAL(RK)  ::  rss(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  gdot(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  comp(0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  temp(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  v_tensor(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  e_elas_kk_grn(0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  determ_v(0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  wp_hat_vec(0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  wp_hat_mat(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  wp_x_e(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  REAL(RK)  ::  f_elas(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
  !
  REAL(RK), pointer  ::  e_vec_tmp(:,:,:) => NULL()
  REAL(RK), POINTER  ::  p_hat_vec(:,:) => NULL()
  !
  integer  ::  iphase, numind
  integer  ::  my_phase(0:(el_sup1-el_sub1))
  !
  !----------------------------------------------------------------------

  my_phase(:) = phase(el_sub1:el_sup1)
  !
  m_el = el_sup1 - el_sub1 + 1
  !
  tan_stif_evp = 0.0_RK
  c_tan = 0.0_RK
  !
  !     Preliminary computations
  !     ------------------------
  !
  !     Factor to correct [c] & {f}
  !
  sqr2  = dsqrt(2.0_RK)
  sqr32 = dsqrt(1.5_RK)
  !
  alpha(1) = 1.0_RK / sqr2
  alpha(2) = sqr32
  alpha(3) = sqr2
  alpha(4) = sqr2
  alpha(5) = sqr2
  !
  !     Elastic strains.
  !
  ! compute elastic strain: sig_vec --> e_vec
  do iphase=1,numphases
     call find_indices(numind, iphase, my_phase, indices)
     indices=indices+el_sub1
     if (associated(e_vec_tmp)) then
        deallocate(e_vec_tmp)
     endif
     allocate(e_vec_tmp(0:TVEC1,0:ngrain1,0:numind-1))
     call vec_d_vec5(keinv(:,iphase), sig_vec(:,:,indices), e_vec_tmp, ngrain, numind)
     e_vec(:,:,indices)=e_vec_tmp
     deallocate(indices)
     deallocate(e_vec_tmp)
  enddo !numphases
  !
  !
  !     Rotate e_bar_vec to current configuration (lattice axes).
  !
  call matt_x_mat3(rstar_n, rstar, d_rstar, ngrain, m_el)
  !
  call rot_mat_symm(d_rstar, qr5x5, ngrain, m_el)
  !
  call lattice_deform(qr5x5, e_bar_vec, e_bar_vec_r, ngrain, m_el)
  !
  !     Transform e_bar_vec_r & e_vec to sample axes: {.}=[Q] {.}.
  !
  call rot_mat_symm(c_angs, qr5x5, ngrain, m_el)
  !
  call mat_x_vec5(qr5x5, e_bar_vec_r, e_bar_sm, ngrain, m_el)
  !
  call mat_x_vec5(qr5x5, e_vec, e_vec_sm, ngrain, m_el)
  !
  !     Tensor form of e_bar_sm, e_vec_sm: {.} -> [ ].
  !
  call vec_mat_symm_grn(e_bar_sm, e_bar, ngrain, m_el)
  !
  call vec_mat_symm_grn(e_vec_sm, e_elas, ngrain, m_el)
  !
  !     {w}_sm to all grains in aggregate (Taylor Assumption).
  !
  do i = 0, DIMS1
     w_vec_grn(i, :, :) = spread(w_vec(i, :), dim = 1,  ncopies = ngrain)
  enddo
  !
  !     Elasto-visco-plastic crystal stiffness
  !     --------------------------------------
  !
  !     Visco-plastic crystal compliance.
  !
  tan_stif_vp = 0.0_RK
  !
  dtimei = 1. / dtime
  !
  do iphase=1,numphases
     call CrystalTypeGet(ctype(iphase), DEV=p_hat_vec)
     n_slip=ctype(iphase)%numslip
     call find_indices(numind, iphase, my_phase, indices)
     indices=indices+el_sub1

     do islip = 0, n_slip - 1
        ! compute rss
        call ss_project(rss(islip,:,:),p_hat_vec(:,islip+1), sig_vec, ngrain, m_el, numind, indices-el_sub1)
        rss(islip,:,indices)=rss(islip,:,indices) /crss(islip,:,indices)
        !
        ! tsh
        where (abs(rss(islip,:,indices)) .lt. t_min(iphase)) 
           rss(islip,:,indices) = 0.0_RK
        endwhere
        ! 
        call power_law(gdot(islip, :, :), rss(islip, :, :),&
             &        crystal_parm(0,iphase), crystal_parm(1,iphase),&
             &        t_min(iphase), ngrain, m_el, numind, indices-el_sub1)
        !
        call compliance(comp, rss(islip, :, :),&
             &         gdot(islip, :, :), crss(islip,:,:), crystal_parm(0,iphase),&
             &         t_min(iphase), ngrain, m_el, numind, indices-el_sub1)
        ! 
        ! hritz 9/15/05
        ! in order to remove one dimension from my_phase (the dimension over
        ! ngrain) we added a loop over that dimension here.
  ! RC 6/24/16: Reordered loop ordering for better memory striding.
        do k=0,ngrain1
           do j = 0, TVEC1
              do i = 0, TVEC1
                 where (my_phase .eq. iphase)
                    tan_stif_vp(i, j, k, :) = tan_stif_vp(i, j, k, :) +&
                         &   comp(k,:) *&
                         &   p_hat_vec(i+1, islip+1) * p_hat_vec(j+1, islip+1)
                 endwhere
              enddo
           enddo
        enddo !ngrain
     enddo !n_slip

     deallocate(p_hat_vec)
     deallocate(indices)

  enddo !numphases
  !
  !     Transform tan_stif_vp to sample axes: [S]_sm = [Q] [S]_lat [Q]'.
  !
  call mat_x_mat5(qr5x5, tan_stif_vp, temp, ngrain, m_el)
  !
  call mat_x_matt5(temp, qr5x5, tan_stif_vp, ngrain, m_el)
  !
  !     Calculate secant moduli
  !
  ! RC 6/24/16: Reordered loop ordering for better memory striding
  do iphase=1,numphases
     do k=0,ngrain1
        do j = 0, TVEC1
           do i = 0, TVEC1
              where (my_phase .eq. iphase)
                 stif_vp(i, j, k, :) = crystal_parm(0,iphase)*tan_stif_vp(i, j, k, :)
              endwhere
           enddo
        enddo
     enddo
  enddo !numphases  
  !
  !     Spread keinv to all crystals and transform it to sample axes.  
  !     (deb 6/11/2000)
  !     
  ! tsh
  keinv_all = 0.0_RK
  ! hritz 9/15/05
  ! in order to remove one dimension from my_phase (the dimension over
  ! ngrain) we added a loop over that dimension here.
  do iphase=1,numphases
     do k=0,ngrain1
        do i = 0, TVEC1
           where (my_phase .eq. iphase)
              keinv_all(i, i, k, :) = keinv(i,iphase)
           endwhere
        enddo
     enddo
  enddo !numphases
  !
  call mat_x_mat5(qr5x5, keinv_all, temp, ngrain, m_el)
  !
  call mat_x_matt5(temp, qr5x5, keinv_all, ngrain, m_el)
  !
  !     Determinant of tensor V* .
  !
  v_tensor = e_elas
  e_elas_kk_grn = spread(e_elas_kk, dim = 1, ncopies = ngrain)
  !
  do i = 0, DIMS1
     v_tensor(i, i, :, :) = v_tensor(i, i, :, :) +  e_elas_kk_grn / 3. + 1.
  enddo
  !
  call determinant_grn(v_tensor, determ_v, ngrain, m_el)
  !
  !     Elasto-visco-plastic crystal compliance.
  !
  !
  ! RC 6/24/16: Reordered loop ordering for better memory striding
  do j = 0, TVEC1
     do i = 0, TVEC1
        stif_evp(i, j, :, :) = determ_v * stif_vp(i, j, :, :) +&
             &   determ_v * keinv_all(i, j, :, :) * dtimei 
     enddo
  enddo

  if (NR) then
     do i = 0, TVEC1
        do j = 0, TVEC1
           tan_stif_evp(i, j, :, :) = determ_v * tan_stif_vp(i, j, :, :) +&
                &   determ_v * keinv_all(i, j, :, :) * dtimei 
        enddo
     enddo
  end if

  !     Elasto-visco-plastic crystal stiffness.
  !
  call invert5x5(stif_evp, ngrain, m_el)
  if (NR) call invert5x5(tan_stif_evp, ngrain, m_el)
  !
  !     Force vector due to elastic terms
  !     ---------------------------------
  !
  call find_wp_hat(wp_hat_vec, e_elas, e_bar, w_vec_grn, gdot, qr5x5, dtime, ngrain, m_el)

  call wp_hat_mat5x5_all(wp_hat_vec, wp_hat_mat, ngrain, m_el)  

  call mat_x_vec5(wp_hat_mat, e_vec_sm, wp_x_e, ngrain, m_el)
  !
  wp_x_e = wp_x_e - 1. / dtime * e_bar_sm
  !
  call mat_x_vec5(stif_evp, wp_x_e, f_elas, ngrain, m_el)
  ! 
  !     Averaged values of stif_evp, f_elas
  !     -----------------------------------
  !
  ! RC 6/24/16: Reordered loop ordering for better memory striding
  do j = 0, TVEC1
     f(j + 1, :) = sum(f_elas(j, :, :) * wts, dim = 1)
     do i= 0, TVEC1
        c(i + 1, j + 1, :) = sum(stif_evp(i, j, :, :) * wts,  dim = 1)
     enddo
  enddo

  if (NR) then
  ! RC 6/24/16: Reordered loop ordering for better memory striding
     do j = 0, TVEC1
        do i= 0, TVEC1
           c_tan(i + 1, j + 1, :) = sum(tan_stif_evp(i, j, :, :) * wts,  dim = 1)
        enddo
     enddo
  end if
  !
  detv = sum(determ_v * wts, dim = 1)
  !
  !     Fix [c] & {f} to be consistent with FE equations.
  !     ------------------------------------------------
  !
  ! RC 6/24/16: Reordered loop ordering for better memory striding
  do j = 1, TVEC
     f(j, :) = alpha(j) * f(j, :)
     do i = 1, TVEC
        c(i, j, :) = alpha(i) * c(i, j, :) * alpha(j)
     enddo
  enddo

  if (NR) then
  ! RC 6/24/16: Reordered loop ordering for better memory striding
     do j = 1, TVEC
        do i = 1, TVEC
           c_tan(i, j, :) = alpha(i) * c_tan(i, j, :) * alpha(j)
        enddo
     enddo
  end if
  !
  RETURN
END SUBROUTINE aniso_evps
!
!**********************************************************************
!
SUBROUTINE find_wp_hat(&
     &   wp_hat, e_elas, e_bar, w_vec_grn, gdot, qr5x5,&
     &   dt, n, m)
  !
  !----------------------------------------------------------------------
  !
  !     Arguments:
  !
  INTEGER, INTENT(IN)    :: n, m
  !
  REAL(RK), INTENT(OUT)  ::  wp_hat(0:DIMS1, 0:(n - 1), 0:(m - 1))  
  REAL(RK), INTENT(IN)   ::  dt
  REAL(RK), INTENT(IN)   ::  e_elas(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
  REAL(RK), INTENT(IN)   ::  e_bar(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
  REAL(RK), INTENT(IN)   ::  w_vec_grn(0:DIMS1, 0:(n - 1), 0:(m - 1))
  REAL(RK), INTENT(IN)   ::  gdot(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
  REAL(RK), INTENT(IN)   ::  qr5x5(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
  !
  !     Locals:
  !
  INTEGER  ::  i, islip, iphase, numind, n_slip
  INTEGER, POINTER  ::  indices(:) => NULL()
  !
  REAL(RK), POINTER  ::  p_hat_vec(:,:) => NULL()
  REAL(RK)  ::  ee(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
  REAL(RK)  ::  dp_hat(0:TVEC1, 0:(n - 1), 0:(m - 1))
  REAL(RK)  ::  temp(0:TVEC1, 0:(n - 1), 0:(m - 1))
  REAL(RK)  ::  dp_hat_tens(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
  REAL(RK)  ::  x(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
  integer   ::  my_phase(0:(m-1))
  !
  !----------------------------------------------------------------------
  !
  my_phase(:) = phase(el_sub1:el_sup1)
  !
  call mat_x_mat3(e_elas, e_bar, ee, n, m)
  !
  wp_hat(0, :, :) = w_vec_grn(0, :, :) + 0.5 / dt * (ee(1, 0, :, :) - ee(0, 1, :, :))
  wp_hat(1, :, :) = w_vec_grn(1, :, :) + 0.5 / dt * (ee(2, 0, :, :) - ee(0, 2, :, :))
  wp_hat(2, :, :) = w_vec_grn(2, :, :) + 0.5 / dt * (ee(2, 1, :, :) - ee(1, 2, :, :))
  !
  dp_hat = 0.0_RK
  !
  do iphase=1,numphases
     call CrystalTypeGet(ctype(iphase), DEV=p_hat_vec)
     n_slip = ctype(iphase)%numslip
     call find_indices(numind, iphase, my_phase, indices)
     do islip = 0, (n_slip - 1)
        do i = 0, TVEC1
           dp_hat(i, :, indices) = dp_hat(i, :, indices) +&
                &                  gdot(islip, :, indices) *&
                &                  p_hat_vec(i + 1, islip + 1)
        enddo
     enddo
     deallocate(p_hat_vec)
     deallocate(indices)
  enddo
  !
  call mat_x_vec5(qr5x5, dp_hat, temp, n, m)
  !
  call vec_mat_symm_grn(temp, dp_hat_tens, n, m)
  !
  call mat_x_mat3(e_elas, dp_hat_tens, x, n, m)
  !
  wp_hat(0, :, :) = wp_hat(0, :, :) - x(1, 0, :, :) + x(0, 1, :, :)
  wp_hat(1, :, :) = wp_hat(1, :, :) - x(2, 0, :, :) + x(0, 2, :, :)
  wp_hat(2, :, :) = wp_hat(2, :, :) - x(2, 1, :, :) + x(1, 2, :, :)
  !
  RETURN
  !
END SUBROUTINE find_wp_hat

END MODULE AnisoEvpsModule
