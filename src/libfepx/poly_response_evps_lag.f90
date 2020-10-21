! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE PolycrystalResponseEvpsModule

  USE IntrinsicTypesModule, RK=>REAL_KIND
  USE units_mod
  USE DimsModule

  USE READ_INPUT_MOD
  USE microstructure_mod
  USE MATRIX_OPERATIONS_MOD
  USE rstarn_solve_lag_mod
  USE StressSolveVpModule
  USE STRESS_SOLVE_EVPS_MOD, ONLY: STRESS_SOLVE_EVPS
  USE ConvergenceModule, ONLY: cv_options

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: polycrystal_response_evps, polycrystal_response_evps_qp

CONTAINS

      SUBROUTINE polycrystal_response_evps(&
  &     d_vec, w_vec, c0_angs, c_angs,&
  &     sig_vec_n, sig_vec, crss_n, crss, rstar_n, rstar, e_bar_vec,&
  &     wts, epseff, d_kk,&
  &     sig_kk, e_elas_kk_bar, e_elas_kk, jiter_state, keinv,&
  &     incr,&
  &     dtime, converged_solution, auto_time&
  &     )
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      LOGICAL, INTENT(INOUT) :: converged_solution
!
      INTEGER  :: incr, auto_time
      INTEGER, INTENT(OUT) :: jiter_state(0:ngrain1, el_sub1:el_sup1)
!
      REAL(RK)  ::  dtime

      REAL(RK), INTENT(IN) ::  keinv(0:TVEC1,1:numphases)
      REAL(RK), INTENT(IN) ::  d_vec(0:TVEC1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN) ::  w_vec(0:DIMS1, el_sub1:el_sup1, 0:nqpt1)

      REAL(RK), INTENT(IN)    :: c0_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) :: c_angs (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    :: rstar_n(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) :: rstar  (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)

      REAL(RK), INTENT(IN)  :: sig_vec_n(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT) :: sig_vec  (0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)  :: e_bar_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)

      REAL(RK), INTENT(IN)  :: crss_n(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(OUT) :: crss  (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)

      REAL(RK), INTENT(IN)  :: wts(0:ngrain1, el_sub1:el_sup1)

      !tm268_M11 (added INTENT):
      REAL(RK), INTENT(IN)  :: d_kk          (el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)  :: e_elas_kk_bar (el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT) :: sig_kk        (el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT) :: e_elas_kk     (el_sub1:el_sup1, 0:nqpt1)

      !tm268_M12 (added INTENT):
      REAL(RK), INTENT(IN)  :: epseff(el_sub1:el_sup1, 0:nqpt1)
!
!
!     Locals:
!
      INTEGER :: i, j, m_el
!
      REAL(RK) :: d_vec_grn(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK) :: w_vec_grn(0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK) :: epseff_lat(0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
!
!----------------------------------------------------------------------
!
      m_el = el_sup1 - el_sub1 + 1
!
!     (this part should be cleaned when ngrain1 stuff is removed...)
!     Spread {d}_sm & {w}_sm to all grains in aggregate (Taylor Assumption)
    do j = 0, nqpt1
      do i = 0, TVEC1
         d_vec_grn(i, :, :, j) = spread(d_vec(i, :, j), DIM = 1, NCOPIES = ngrain)
      enddo

      do i = 0, DIMS1
         w_vec_grn(i, :, :, j) = spread(w_vec(i, :, j), DIM = 1, NCOPIES = ngrain)
      enddo

      ! spread over grains: epseff --> epseff_lat
      epseff_lat(:,:,j) = spread(epseff(:,j), DIM = 1, NCOPIES = ngrain)
    enddo
!
!--------------------------------
!     Solve for State.
!--------------------------------
      ! deviatoric
      call solve_state_dev_evps(d_vec_grn, w_vec_grn,&
  &        c0_angs, c_angs, sig_vec_n, sig_vec, crss_n, crss, rstar_n,&
  &        rstar, e_bar_vec, keinv, wts, epseff_lat, jiter_state,&
  &        incr, dtime,&
  &        converged_solution, auto_time)

      if (.not. converged_solution .and. auto_time .eq. 1) RETURN

      ! volumetric
    do j = 0, nqpt1
      call solve_state_vol_evps(e_elas_kk_bar(:,j), e_elas_kk(:,j), d_kk(:,j), sig_kk(:,j), dtime, m_el)
    enddo

      END SUBROUTINE polycrystal_response_evps
!
!**********************************************************************
!
      SUBROUTINE solve_state_dev_evps(&
  &     d_vec, w_vec, c0_angs, c_angs, sig_lat_n,&
  &     sig_lat, crss_n, crss, rstar_n, rstar, e_bar_vec,&
  &     keinv, wts, epseff, jiter_state,&
  &     incr, &
  &     dtime, converged_solution,&
  &     auto_time&
  &     )
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      LOGICAL, INTENT(INOUT) :: converged_solution
!
      INTEGER  incr, auto_time

      INTEGER, INTENT(OUT)  :: jiter_state(0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: dtime

      REAL(RK), INTENT(OUT)   :: crss  (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    :: sig_lat_n(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT)   :: sig_lat  (0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
!      REAL(RK), INTENT(OUT)   :: rstar  (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) :: rstar  (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(INOUT) :: c_angs (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    :: keinv(0:TVEC1,1:numphases)
      REAL(RK), INTENT(IN)    :: d_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    :: w_vec(0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    :: c0_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: rstar_n(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: crss_n(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: e_bar_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(IN)    :: wts   (0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: epseff(0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
!
!     Locals:
!
      LOGICAL, PARAMETER :: vp_log = .false. ! write viscoplastic convergence output to log files

      LOGICAL :: done(0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      LOGICAL :: converged_newton(0:nqpt1), converged_state(0:nqpt1)
!
      INTEGER :: iter_state, m_el, islip, n_slip
!
      REAL(RK)  :: qr5x5(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: qr3x3(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: w_vec_lat  (0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK)  :: wp_hat     (0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK)  :: d_vec_lat  (0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK)  :: e_bar_vec_r(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: crss_0     (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK)  :: norm_s_0   (0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK)  :: norm_s     (0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK)  :: diff_norm_s(0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK)  :: diff_crss  (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK)  :: d_rstar(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
!
      integer :: iphase, k, i, j

      integer :: my_phase(0:(el_sup1-el_sub1))
!
!----------------------------------------------------------------------
!
      my_phase(:) = phase(el_sub1:el_sup1)
      m_el = el_sup1 - el_sub1 + 1
!
      converged_newton = .true.
      converged_state  = .true.

      jiter_state  = 0

      crss = spread(crss_n, 4, nqpt)

      done = .false.
!
!--------------------------------------------------------------------------------------------------
!     Estimate for the stresses :
!     For incr = 1 --> based on viscoplastic solution
!     For incr > 1 --> based on previous solution (unrotated)
!
    do i = 0, nqpt1

      if (incr .eq. 1) then
         ! c_angs [3x3] --> qr5x5 [5x5]
         call rot_mat_symm(c_angs(:,:,:,:,i), qr5x5, ngrain, m_el)
         ! d_vec(sample coords) --> d_vec_lat(crystal coords)
         ! {d_vec_lat} = [qr5x5]'{d_vec}
         call lattice_deform(qr5x5, d_vec(:,:,:,i), d_vec_lat(:,:,:,i), ngrain, m_el)
         ! estimate of sig_lat (visco-plastic solution)
         call stress_solve_vp(sig_lat(:,:,:,i), d_vec_lat(:,:,:,i), crss(:,:,:,i), epseff(:,:,i), vp_log)
         ! use a fraction of the viscoplastic solution as initial guess
         sig_lat(:,:,:,i) = 0.5*sig_lat(:,:,:,i)
      else
         ! use value from previous increment
         sig_lat(:,:,:,i) = sig_lat_n(:,:,:,i)
      endif
!--------------------------------------------------------------------------------------------------
!
!     First (Forward Euler) estimate for crss and rstar (c).
!
      ! c_angs [3x3] --> qr3x3 [3x3]
      call rot_mat_skew(c_angs(:,:,:,:,i), qr3x3, ngrain, m_el)
      ! w_vec(sample coords) --> w_vec_lat(crystal coords)
      call lattice_spin(qr3x3, w_vec(:,:,:,i), w_vec_lat(:,:,:,i), ngrain, m_el)
      ! estimate crss (g), rstar (R*) and d_rstar (dR*)
      ! calculate also c_angs: [c]=[c_0]*[R*]

      accumshear = gaccumshear(:,:,:,i)

      call rstarn_solve(crss_n, crss(:,:,:,i), rstar_n, rstar(:,:,:,:,i), c0_angs,&
  &        c_angs(:,:,:,:,i), sig_lat(:,:,:,i), w_vec_lat(:,:,:,i), dtime, epseff(:,:,i),&
  &        done(:,:,i), d_rstar(:,:,:,:,i), UPD_EULER_FWD)
      !
      ! Compute 2-norm for array of 5-vectors
      call norm_vec(norm_s_0(:,:,i), sig_lat(:,:,:,i), ngrain, m_el)
      !
      ! update crss_0
      crss_0(:,:,:,i) = crss(:,:,:,i)

    enddo
!
!
!     Iterate for the material state. --------------------------------------------
!
    do i = 0, nqpt1

      iter_state = 1
      accumshear = gaccumshear(:,:,:,i)

      do while ((any(.not. done(:,:,i))) .and. (iter_state .le. cv_options % sx_max_iters_state))
         !
         ! c_angs [3x3] --> qr5x5 [5x5]
         call rot_mat_symm(c_angs(:,:,:,:,i), qr5x5, ngrain, m_el)
         ! d_vec(sample coords) --> d_vec_lat(crystal coords)
         call lattice_deform(qr5x5, d_vec(:,:,:,i), d_vec_lat(:,:,:,i), ngrain, m_el)
         ! d_rstar [3x3] --> qr5x5 [5x5]
         call rot_mat_symm(d_rstar(:,:,:,:,i), qr5x5, ngrain, m_el)
         ! apply dR* to e_bar_vec --> e_bar_vec_r
         call lattice_deform(qr5x5, e_bar_vec(:,:,:,i), e_bar_vec_r, ngrain, m_el)
         ! --> sig_lat
         call stress_solve_evps(sig_lat(:,:,:,i), d_vec_lat(:,:,:,i), w_vec_lat(:,:,:,i),&
  &           e_bar_vec_r, crss(:,:,:,i), keinv,&
  &           dtime, wp_hat(:,:,:,i), iter_state,&
  &           done(:,:,i), converged_newton(i))
         !
         if (.not. converged_newton(i) .and. auto_time .eq. 1) GO TO 100
         !
         ! c_angs [3x3] --> qr3x3 [3x3]
         call rot_mat_skew(c_angs(:,:,:,:,i), qr3x3, ngrain, m_el)
         ! w_vec(sample coords) --> w_vec_lat(crystal coords)
         call lattice_spin(qr3x3, w_vec(:,:,:,i), w_vec_lat(:,:,:,i), ngrain, m_el)
         ! calculate crss (g), rstar (R*) and d_rstar (dR*)
         ! calculate also c_angs: [c]=[c_0]*[R*]
         call rstarn_solve(crss_n, crss(:,:,:,i), rstar_n, rstar(:,:,:,:,i), c0_angs,&
  &           c_angs(:,:,:,:,i), sig_lat(:,:,:,i), wp_hat(:,:,:,i), dtime, epseff(:,:,i),&
  &           done(:,:,i), d_rstar(:,:,:,:,i), UPD_EULER_BWD)
         !
         !
         call norm_vec(norm_s(:,:,i), sig_lat(:,:,:,i), ngrain, m_el)
!        deb
!        This section was originally dones with nested `where'
!        constructs, but was changed because the AIX compiler
!        rejected them, although the CM compiler had no problem.
!        deb
        do islip=0,MAXSLIP1
         where (.not. done(:,:,i))
            diff_crss(islip,:,:,i)   = dabs(crss(islip,:,:,i)   -   crss_0(islip,:,:,i))
         endwhere
        ENDDO

        where (.not. done(:,:,i))
            diff_norm_s(:,:,i) = dabs(norm_s(:,:,i) - norm_s_0(:,:,i))
        endwhere
         !
         !
         do k=0,ngrain-1
            do iphase=1,numphases
                call CrystalTypeGet(ctype(iphase))
                n_slip=ctype(iphase)%numslip
               do islip=0,n_slip-1
                    !currently have an initial crss_0 scaling factor from the simple latent hardening model but might end up getting rid of it later versions
                   where((.not. done(k,:,i))  .and.  (my_phase.eq.iphase) .and.&
      &                  (diff_norm_s(k,:,i) .lt.  (TOLER_STATE * crystal_parm(3,iphase))) .and.&
      &                  (diff_crss(islip,k,:,i)   .lt.  (TOLER_STATE * crystal_parm(3,iphase))) )
                      done(k,:,i) = .true.
                      jiter_state(k,:) = iter_state
                   endwhere
                enddo !nslips
            enddo !numphases
         enddo ! ngrain
         !
        do islip=0,MAXSLIP1
         where (.not. done(:,:,i))
            norm_s_0(:,:,i) = norm_s(:,:,i)
            crss_0(islip,:,:,i)   = crss(islip,:,:,i)
         endwhere
        enddo
         !
         iter_state = iter_state + 1
         !
         !
      enddo ! do-while loop ------------------------------------------
    enddo !nqpt do loop------------------------------------



    do i = 0,nqpt1
      if (any(.not. done(:,:,i))) then
         converged_state(i) = .false.
        WRITE(DFLT_U, *) 'Warning:       . Not all crystals converged.'
        WRITE(DFLT_U, *) 'Warning:       . Crystals = ', ngrain * m_el, ', &
            &converged = ', count(done(:,:,i)) 
!         write(*, 1000) ngrain * m_el, count(done(:,:,i))
! 1000    FORMAT (&
!  &        ' solve_state: Not all crystals converged ', /,&
!  &        ' solve_state: Crystals =', i7, ' Converged =', i7)
      endif
    enddo


 100  if ( any(.not. converged_newton) .or. any(.not. converged_state) ) converged_solution = .false.

      END SUBROUTINE solve_state_dev_evps
!
!**********************************************************************
!
      SUBROUTINE solve_state_vol_evps(&
  &       e_elas_kk_bar, e_elas_kk, d_kk, sig_kk, dtime, m)
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER :: m
!
      REAL(RK), INTENT(IN)  :: dtime
      REAL(RK), INTENT(IN)  :: e_elas_kk_bar(0:(m - 1))
      REAL(RK), INTENT(IN)  :: d_kk  (0:(m - 1))
      REAL(RK), INTENT(OUT) :: e_elas_kk(0:(m - 1))
      REAL(RK), INTENT(OUT) :: sig_kk(0:(m - 1))
!
!     Locals:
!
      integer :: my_phase(0:m-1), iphase
!
!----------------------------------------------------------------------
!
      my_phase(:) = phase(el_sub1:el_sup1)
!
!     Compute volumetric response in sample reference frame

      e_elas_kk = e_elas_kk_bar + dtime * d_kk

      do iphase=1,numphases
         where (my_phase .eq. iphase)
            sig_kk    = 3.0_RK * crystal_parm(8,iphase) * e_elas_kk
         endwhere
      enddo !numphases

      END SUBROUTINE solve_state_vol_evps

      SUBROUTINE polycrystal_response_evps_qp(&
  &     d_vec, w_vec, c0_angs, c_angs,&
  &     sig_vec_n, sig_vec, crss_n, crss, rstar_n, rstar, e_bar_vec,&
  &     wts, epseff, d_kk,&
  &     sig_kk, e_elas_kk_bar, e_elas_kk, jiter_state, keinv,&
  &     incr,&
  &     dtime, converged_solution, auto_time&
  &     )
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      LOGICAL, INTENT(INOUT) :: converged_solution
!
      INTEGER  :: incr, auto_time
      INTEGER, INTENT(OUT) :: jiter_state(0:ngrain1, el_sub1:el_sup1)
!
      REAL(RK)  ::  dtime

      REAL(RK), INTENT(IN) ::  keinv(0:TVEC1,1:numphases)
      REAL(RK), INTENT(IN) ::  d_vec(0:TVEC1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN) ::  w_vec(0:DIMS1, el_sub1:el_sup1)

      REAL(RK), INTENT(IN)    :: c0_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) :: c_angs (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: rstar_n(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) :: rstar  (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)

      REAL(RK), INTENT(IN)  :: sig_vec_n(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(OUT) :: sig_vec  (0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)  :: e_bar_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)

      REAL(RK), INTENT(IN)  :: crss_n(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(OUT) :: crss  (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)

      REAL(RK), INTENT(IN)  :: wts(0:ngrain1, el_sub1:el_sup1)

      !tm268_M11 (added INTENT):
      REAL(RK), INTENT(IN)  :: d_kk          (el_sub1:el_sup1)
      REAL(RK), INTENT(IN)  :: e_elas_kk_bar (el_sub1:el_sup1)
      REAL(RK), INTENT(OUT) :: sig_kk        (el_sub1:el_sup1)
      REAL(RK), INTENT(OUT) :: e_elas_kk     (el_sub1:el_sup1)

      !tm268_M12 (added INTENT):
      REAL(RK), INTENT(IN)  :: epseff(el_sub1:el_sup1)
!
!
!     Locals:
!
      INTEGER :: i, m_el
!
      REAL(RK) :: d_vec_grn(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK) :: w_vec_grn(0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK) :: epseff_lat(0:ngrain1, el_sub1:el_sup1)
!
!----------------------------------------------------------------------
!
      m_el = el_sup1 - el_sub1 + 1
!
!     (this part should be cleaned when ngrain1 stuff is removed...)
!     Spread {d}_sm & {w}_sm to all grains in aggregate (Taylor Assumption)
      do i = 0, TVEC1
         d_vec_grn(i, :, :) = spread(d_vec(i, :), DIM = 1, NCOPIES = ngrain)
      enddo

      do i = 0, DIMS1
         w_vec_grn(i, :, :) = spread(w_vec(i, :), DIM = 1, NCOPIES = ngrain)
      enddo

      ! spread over grains: epseff --> epseff_lat
      epseff_lat = spread(epseff, DIM = 1, NCOPIES = ngrain)
!
!--------------------------------
!     Solve for State.
!--------------------------------
      ! deviatoric
      call solve_state_dev_evps_qp(d_vec_grn, w_vec_grn,&
  &        c0_angs, c_angs, sig_vec_n, sig_vec, crss_n, crss, rstar_n,&
  &        rstar, e_bar_vec, keinv, wts, epseff_lat, jiter_state,&
  &        incr, dtime,&
  &        converged_solution, auto_time)

      if (.not. converged_solution .and. auto_time .eq. 1) RETURN

      ! volumetric
      call solve_state_vol_evps_qp(e_elas_kk_bar, e_elas_kk, d_kk, sig_kk, dtime, m_el)

      END SUBROUTINE polycrystal_response_evps_qp
!
!**********************************************************************
!
      SUBROUTINE solve_state_dev_evps_qp(&
  &     d_vec, w_vec, c0_angs, c_angs, sig_lat_n,&
  &     sig_lat, crss_n, crss, rstar_n, rstar, e_bar_vec,&
  &     keinv, wts, epseff, jiter_state,&
  &     incr, &
  &     dtime, converged_solution,&
  &     auto_time&
  &     )
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      LOGICAL, INTENT(INOUT) :: converged_solution
!
      INTEGER  incr, auto_time

      INTEGER, INTENT(OUT)  :: jiter_state(0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: dtime

      REAL(RK), INTENT(OUT)   :: crss  (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: sig_lat_n(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(OUT)   :: sig_lat  (0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
!      REAL(RK), INTENT(OUT)   :: rstar  (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) :: rstar  (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) :: c_angs (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: keinv(0:TVEC1,1:numphases)
      REAL(RK), INTENT(IN)    :: d_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: w_vec(0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: c0_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: rstar_n(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: crss_n(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: e_bar_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: wts   (0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: epseff(0:ngrain1, el_sub1:el_sup1)
!
!     Locals:
!
      LOGICAL, PARAMETER :: vp_log = .false. ! write viscoplastic convergence output to log files

      LOGICAL :: done(0:ngrain1, el_sub1:el_sup1)
      LOGICAL :: converged_newton, converged_state      
!
      INTEGER :: iter_state, m_el, islip, n_slip
!
      REAL(RK)  :: qr5x5(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: qr3x3(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: w_vec_lat  (0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: wp_hat     (0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: d_vec_lat  (0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: e_bar_vec_r(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: crss_0     (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: norm_s_0   (0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: norm_s     (0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: diff_norm_s(0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: diff_crss  (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  :: d_rstar(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
!
      integer :: iphase, k

      integer :: my_phase(0:(el_sup1-el_sub1))
!
!----------------------------------------------------------------------
!
      my_phase(:) = phase(el_sub1:el_sup1)
      m_el = el_sup1 - el_sub1 + 1
!
      converged_newton = .true.
      converged_state  = .true.

      jiter_state  = 0

      crss = crss_n
      done = .false.
!
!--------------------------------------------------------------------------------------------------
!     Estimate for the stresses :
!     For incr = 1 --> based on viscoplastic solution
!     For incr > 1 --> based on previous solution (unrotated)
!
      if (incr .eq. 1) then
         ! c_angs [3x3] --> qr5x5 [5x5]
         call rot_mat_symm(c_angs, qr5x5, ngrain, m_el)
         ! d_vec(sample coords) --> d_vec_lat(crystal coords)
         ! {d_vec_lat} = [qr5x5]'{d_vec}
         call lattice_deform(qr5x5, d_vec, d_vec_lat, ngrain, m_el)
         ! estimate of sig_lat (visco-plastic solution)
         call stress_solve_vp(sig_lat, d_vec_lat, crss, epseff, vp_log)
         ! use a fraction of the viscoplastic solution as initial guess
         sig_lat = 0.5*sig_lat
      else
         ! use value from previous increment
         sig_lat = sig_lat_n
      endif
!--------------------------------------------------------------------------------------------------
!
!     First (Forward Euler) estimate for crss and rstar (c).
!
      ! c_angs [3x3] --> qr3x3 [3x3]
      call rot_mat_skew(c_angs, qr3x3, ngrain, m_el)
      ! w_vec(sample coords) --> w_vec_lat(crystal coords)
      call lattice_spin(qr3x3, w_vec, w_vec_lat, ngrain, m_el)
      ! estimate crss (g), rstar (R*) and d_rstar (dR*)
      ! calculate also c_angs: [c]=[c_0]*[R*]
      call rstarn_solve(crss_n, crss, rstar_n, rstar, c0_angs, c_angs,&
  &        sig_lat, w_vec_lat, dtime, epseff, done, d_rstar, UPD_EULER_FWD)
      !
      ! Compute 2-norm for array of 5-vectors
      call norm_vec(norm_s_0, sig_lat, ngrain, m_el)
      !
      ! update crss_0
      crss_0 = crss
!
!
!     Iterate for the material state. --------------------------------------------
!
      iter_state = 1

      do while ((any(.not. done)) .and. (iter_state .le. cv_options % sx_max_iters_state))
         !
         ! c_angs [3x3] --> qr5x5 [5x5]
         call rot_mat_symm(c_angs, qr5x5, ngrain, m_el)
         ! d_vec(sample coords) --> d_vec_lat(crystal coords)
         call lattice_deform(qr5x5, d_vec, d_vec_lat, ngrain, m_el)
         ! d_rstar [3x3] --> qr5x5 [5x5]
         call rot_mat_symm(d_rstar, qr5x5, ngrain, m_el)
         ! apply dR* to e_bar_vec --> e_bar_vec_r
         call lattice_deform(qr5x5, e_bar_vec, e_bar_vec_r, ngrain, m_el)
         ! --> sig_lat
         call stress_solve_evps(sig_lat, d_vec_lat, w_vec_lat,&
  &           e_bar_vec_r, crss, keinv,&
  &           dtime, wp_hat, iter_state,&
  &           done, converged_newton)
         !
         if (.not. converged_newton .and. auto_time .eq. 1) GO TO 100
         !
         ! c_angs [3x3] --> qr3x3 [3x3]
         call rot_mat_skew(c_angs, qr3x3, ngrain, m_el)
         ! w_vec(sample coords) --> w_vec_lat(crystal coords)
         call lattice_spin(qr3x3, w_vec, w_vec_lat, ngrain, m_el)
         ! calculate crss (g), rstar (R*) and d_rstar (dR*)
         ! calculate also c_angs: [c]=[c_0]*[R*]
         call rstarn_solve(crss_n, crss, rstar_n, rstar, c0_angs, c_angs,&
  &           sig_lat, wp_hat, dtime, epseff, done, d_rstar, UPD_EULER_BWD)
         !
         !
         call norm_vec(norm_s, sig_lat, ngrain, m_el)
!        deb
!        This section was originally dones with nested `where'
!        constructs, but was changed because the AIX compiler
!        rejected them, although the CM compiler had no problem.
!        deb
        do islip=0,MAXSLIP1
         where (.not. done)
            diff_crss(islip,:,:)   = dabs(crss(islip,:,:)   -   crss_0(islip,:,:))
         endwhere
        ENDDO

        where (.not. done)
            diff_norm_s = dabs(norm_s - norm_s_0)
        endwhere
         !
         !
         do k=0,ngrain-1
            do iphase=1,numphases
                call CrystalTypeGet(ctype(iphase))
                n_slip=ctype(iphase)%numslip
               do islip=0,n_slip-1
                    !currently have an initial crss_0 scaling factor from the simple latent hardening model but might end up getting rid of it later versions
                   where((.not. done(k,:))  .and.  (my_phase.eq.iphase) .and.&
      &                  (diff_norm_s(k,:) .lt.  (TOLER_STATE * crystal_parm(3,iphase))) .and.&
      &                  (diff_crss(islip,k,:)   .lt.  (TOLER_STATE * crystal_parm(3,iphase))) )
                      done(k,:) = .true.
                      jiter_state(k,:) = iter_state
                   endwhere
                enddo !nslips
            enddo !numphases
         enddo ! ngrain
         !
        do islip=0,MAXSLIP1
         where (.not. done)
            norm_s_0 = norm_s
            crss_0(islip,:,:)   = crss(islip,:,:)
         endwhere
        enddo
         !
         iter_state = iter_state + 1
         !
         !
      enddo ! do-while loop ------------------------------------------




      if (any(.not. done)) then
         converged_state = .false.
        WRITE(DFLT_U, *) 'Warning:       . Not all crystals converged.'
        WRITE(DFLT_U, *) 'Warning:       . Crystals = ', ngrain * m_el, ', &
            &converged = ', count(done) 
!         write(*, 1000) ngrain * m_el, count(done)
! 1000    FORMAT (&
!  &        ' solve_state: Not all crystals converged ', /,&
!  &        ' solve_state: Crystals =', i7, ' Converged =', i7)
      endif


 100  if ( (.not. converged_newton) .or. (.not. converged_state) ) converged_solution = .false.

      END SUBROUTINE solve_state_dev_evps_qp
!
!**********************************************************************
!
      SUBROUTINE solve_state_vol_evps_qp(&
  &       e_elas_kk_bar, e_elas_kk, d_kk, sig_kk, dtime, m)
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER :: m
!
      REAL(RK), INTENT(IN)  :: dtime
      REAL(RK), INTENT(IN)  :: e_elas_kk_bar(0:(m - 1))
      REAL(RK), INTENT(IN)  :: d_kk  (0:(m - 1))
      REAL(RK), INTENT(OUT) :: e_elas_kk(0:(m - 1))
      REAL(RK), INTENT(OUT) :: sig_kk(0:(m - 1))
!
!     Locals:
!
      integer :: my_phase(0:m-1), iphase
!
!----------------------------------------------------------------------
!
      my_phase(:) = phase(el_sub1:el_sup1)
!
!     Compute volumetric response in sample reference frame

      e_elas_kk = e_elas_kk_bar + dtime * d_kk

      do iphase=1,numphases
         where (my_phase .eq. iphase)
            sig_kk    = 3.0_RK * crystal_parm(8,iphase) * e_elas_kk
         endwhere
      enddo !numphases

      END SUBROUTINE solve_state_vol_evps_qp


END MODULE PolycrystalResponseEvpsModule
