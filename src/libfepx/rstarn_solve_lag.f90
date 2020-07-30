! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE rstarn_solve_lag_mod

! This module is concerned with updating the orientations (rstar).
!
! TO DO
!
! * Handle masks (e.g. done/converged array) in helper subroutines, e.g. get_c
! * Make use of MATMUL in get_c

USE LibF95, ONLY: RK => REAL_KIND, RK_ZERO, RK_ONE, &
   &            STDERR => STANDARD_ERROR
USE parallel_mod, ONLY: par_message
USE units_mod, ONLY: ounits, LOG_U, DFLT_U

USE StressSolveVpModule
USE StressSolveEvpsModule
USE HardeningModule

IMPLICIT NONE
!
! MODULE DATA
!
PRIVATE

PUBLIC :: rstarn_solve, slip_qnt
PUBLIC :: UPD_EULER_FWD, UPD_EULER_BWD
!
! * Flags and Identifiers
!
INTEGER, PARAMETER :: UPD_EULER_FWD=1, UPD_EULER_BWD=2
!
! * Temporary storage
!
CHARACTER(LEN=128) :: message
!
! * Hardening model iteration
!
!   The equation should be quadratic for the usual hardening model,
!   and so it should only require one iteration.
!
REAL(RK), PARAMETER :: TOLER_HARD=1.0e-4_RK, LS_CUTOFF=0.001_RK
INTEGER, PARAMETER :: MAX_ITER_HARD=10

CONTAINS
!
!**********************************************************************
!
      SUBROUTINE rstarn_solve(&
  &       crss_n, crss, rstar_n, rstar, c_0, c, sig_lat,&
  &       wp_hat, dtime, epseff, done, d_rstar, iflag )
!
!----------------------------------------------------------------------
!
      USE READ_INPUT_MOD
      USE microstructure_mod
      USE DimsModule
!
      IMPLICIT  NONE
!
!     Arguments:
!

      INTEGER n_slip, iflag
      REAL(RK) :: dtime
      
      LOGICAL, INTENT(IN)    ::  done(0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)     ::  rstar_n(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT)  ::  rstar  (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)      
      REAL(RK), INTENT(INOUT)  ::  d_rstar(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT)  ::  crss  (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)     ::  crss_n(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(OUT)    ::  c  (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)     ::  c_0(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)     ::  sig_lat(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)     ::  wp_hat (0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)     ::  epseff(0:ngrain1, el_sub1:el_sup1)      
!
!     Locals:
!
      REAL(RK)  ::  shrate(0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  ::  wp_ss(0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK)  ::  shear(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
!
      INTEGER   ::  m_el,i
!
!----------------------------------------------------------------------
!
      m_el = el_sup1 - el_sub1 + 1

      ! compute gammadots (shrate) and wp_ss
      call slip_qnt(wp_ss, shear, shrate, sig_lat, crss, epseff, ngrain, m_el)
      ! update crss     
      call upd_hard(crss, crss_n, shear, shrate, epseff, dtime, iflag, done, ngrain, m_el)
      ! update rstar (R*) and d_rstar(dR*)
      call upd_rotn(rstar, rstar_n, d_rstar, dtime, wp_hat, wp_ss, done, ngrain, m_el)
      ! compute c = c_0 * R*
      call get_c(c, rstar, c_0, done, ngrain, m_el)

      RETURN
      END SUBROUTINE
!
!**********************************************************************
!
      SUBROUTINE slip_qnt(&
  &      wp_ss, shear, shrate, sig_lat, crss, epseff, n, m )
!
!----------------------------------------------------------------------
!
      USE UtilsCrystalModule
      USE READ_INPUT_MOD
      use microstructure_mod
!
      IMPLICIT  NONE
!
!
!     Arguments:
!
      INTEGER   ::  n_slip, n, m
      REAL(RK), INTENT(OUT)  ::    wp_ss(0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT)  ::    shrate(0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT)  ::  shear(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)   ::    sig_lat(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)   ::    crss(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)   ::    epseff(0:(n - 1), 0:(m - 1))
!
!
!     Locals:
!
      INTEGER  ::  islip, i, iphase, is, numind,j
! 
!
!
      REAL(RK), POINTER :: plocal(:,:) => NULL()
      REAL(RK), POINTER :: qlocal(:,:) => NULL()
      INTEGER,  POINTER :: indices(:) => NULL()
      REAL(RK)  ::  rss(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  shr_min(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  shr_max(0:(n - 1), 0:(m - 1))

  
      INTEGER  ::  my_phase(0:(m - 1))

!----------------------------------------------------------------------

! tsh, 1/27/03
      ! crystal_type of the elements (on this node)
      my_phase(:) = phase(el_sub1:el_sup1)

      ! initialize the outputs
      shrate(:,:) = 0.0_RK
      wp_ss(:,:,:) = 0.0_RK
      shear(:,:,:) = 0.0_RK

      ! for each phase
      do iphase = 1,numphases
        
         call CrystalTypeGet(ctype(iphase), DEV=plocal, SKW=qlocal)
         n_slip=ctype(iphase)%numslip
          
         ! Compute shear on each slip system.
         call find_indices(numind, iphase, my_phase, indices)
                  
         do islip = 0, n_slip - 1

            call ss_project(rss(islip,:,:), plocal(:,islip+1), sig_lat, n, m, numind, indices)            
            rss(islip,:,indices)=rss(islip,:,indices)/crss(islip,:,indices)
            where (abs(rss(islip,:,indices)) .lt. t_min(iphase)) 
               rss(islip,:,indices) = 0.0d0
            endwhere
            call power_law(shear(islip, :, :), rss(islip, :, :),&
  &              crystal_parm(0,iphase), crystal_parm(1,iphase),&
  &              t_min(iphase), n, m, numind, indices)
         enddo !n_slip

         gammadot(0:n_slip-1,:,el_sub1+indices)=shear(0:n_slip-1,:,indices)     
         where( abs(gammadot(0:n_slip-1,:,el_sub1+indices)) .le. 1d-30 )
               gammadot(0:n_slip-1,:,el_sub1+indices)=0.0d0
         endwhere
    
         ! Compute net shearing rate on all slip systems.
         do is = 0, n_slip - 1
            shrate(:, indices) = shrate(:, indices) + abs(shear(is, :, indices))
         end do
 
         ! Compute Plastic Spin.
         do is =  0, n_slip - 1
            do i = 0, DIMS1
               wp_ss(i, :, indices) = wp_ss(i, :, indices) + qlocal(i+1, is+1) * shear(is, :, indices)
            enddo
         enddo !n_slip
         !
         deallocate(indices)
         deallocate(plocal)
         deallocate(qlocal)
         !
      enddo !numphases
      
          
      shr_min = 1.0d-6 * epseff
      shr_max = 1.0d1 * epseff
      
      where (shrate .le. shr_min)   shrate = shr_min
      where (shrate .ge. shr_max)   shrate = shr_max
!
      RETURN
      END SUBROUTINE
!
!**********************************************************************
!
      SUBROUTINE upd_hard(crss, crss_0, shear, shrate, epseff, dtime, iflag, done, n, m )
!
!----------------------------------------------------------------------
! 
      USE READ_INPUT_MOD
      use microstructure_mod
      USE UtilsCrystalModule
!
      IMPLICIT  NONE
!
!     Arguments:
!
!     crss - update hardnesses [should it be in/out?]
!     crss_0 - hardnesses before timestep
!     shrate - shearing rates
!     dtime - timestep 
!     iflag - indicates whether to use forward or backward Euler
!     done - mask of elements not to process
!     n, m - number of grains and elements
!
      INTEGER, INTENT(IN)  ::  iflag, n, m
      REAL(RK), INTENT(IN) ::  dtime
      REAL(RK), INTENT(INOUT) ::  crss  (0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    ::  crss_0(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    ::  shrate(0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    ::  shear(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    ::  epseff(0:(n - 1), 0:(m - 1))
      LOGICAL, INTENT(IN)     ::  done(0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER  ::  iter_hard, nm, inewton, ihard, iphase, islip
      INTEGER  ::  numind, n_slip
      INTEGER, pointer  ::  indices(:) => NULL()
!
      LOGICAL  ::  newton_ok(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      LOGICAL  ::  newton_ok_all(0:(n - 1), 0:(m - 1))
      LOGICAL  ::  converged(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      LOGICAL  ::  converged_all(0:(n - 1), 0:(m - 1))
!
      REAL(RK)  ::  hard_rate(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  dhard_rate(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  crss_sat(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  crss_tmp(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  del_crss(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  res_n(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  res(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  ratio_res(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  fj31(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  fjac(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  xlam(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))


      integer  :: my_phase(0:(n - 1), 0:(m - 1)), i
!
!----------------------------------------------------------------------
!
!      h_0    = crystal_parm(2)
!      tausi  = crystal_parm(3)
!      taus0  = crystal_parm(4) ->g_1 or g_s0
!      xms    = crystal_parm(5) ->m'
!      gamss0 = crystal_parm(6) ->gammadot_s0

       my_phase(n-1,:) = phase(el_sub1:el_sup1)     

      ! hardwired
      ihard = 2
      do iphase=1,numphases
         !       
        WHERE (my_phase.EQ.iphase)
          WHERE(shrate .GT. RK_ZERO)
            crss_sat = crystal_parm(4,iphase) &
                 &   * ((shrate / crystal_parm(6,iphase)) &
                 &   ** crystal_parm(5,iphase))
          ELSEWHERE
            crss_sat = crystal_parm(4,iphase)
          END WHERE
        END WHERE
         !
         call find_indices(numind, iphase, my_phase(n-1,:), indices)
         if (crystal_parm(2,iphase) .eq. 0.d0)  then          
            crss(:,:,indices) = crss_0(:,:,indices)
         endif
         !
         DEALLOCATE(indices)
         !         
      enddo !numphases
      
      IF (iflag == UPD_EULER_FWD) THEN

        !Initial guess - Forward Euler

        CALL hard_law(hard_rate, dhard_rate, crss_0, crss_sat, shear, shrate, epseff, 1, ihard, n, m)
    
        crss = crss_0 + dtime * hard_rate

        RETURN 

      END IF
!
!     Start Newton iteration - Backward Euler approx. to hardening law
!
      !The newton_ok and converged variables are for each individual slip system
      !The newton_ok_all and converged_all variables are for each individual element
      !It was done this way to reduce the number of changes needed to be made to the legacy code.
      newton_ok = .TRUE.
      newton_ok_all= .TRUE.
      converged = .TRUE.
      converged_all= .TRUE.

      !Initializing all the converged indices slip systems for each individual phase to false
      !In a dual or more phase system where the different phases have different 
      !number of slip systems this allows the following case where line 1 is a
      !system with 3 slip systems and line 2 only has one for the converged variable
      !FFF
      !FTT
      !So in the convergence variable only the indices that are within the range of that phases number of slip systems are false
      do iphase=1,numphases
           call CrystalTypeGet(ctype(iphase))
           n_slip=ctype(iphase)%numslip
           do islip=0,n_slip-1
                where (my_phase .eq. iphase) converged(islip,:,:)=.FALSE.
           enddo
      enddo

      !Finding where done is true across the element setting converged true for that entire element
      do islip=0,MAXSLIP1
        WHERE(done) converged(islip,:,:) = .TRUE.
      enddo
!
!     This section is wrong for model with limiting hardness.
!
!      where ((crss_0 .eq. crss_sat) .and. (.not. done))
!        crss = crss_0
!        converged = .true.
!      endwhere

      nm = n * m

      DO iter_hard = 1, MAX_ITER_HARD
        !
        ! Should we add "where (.not. converged) clause here?
        ! .   Yes, but we need to pass mask argument to hard_law
        !
        crss_tmp = crss

        call hard_law(hard_rate, dhard_rate, crss_tmp, crss_sat, shear, shrate, epseff, 1, ihard, n, m)
          
        res_n = - (crss_tmp - crss_0 - dtime * hard_rate)

        ! next line to fix 0/0 FPE below
        WHERE (ABS(res_n) == RK_ZERO)
          converged = .TRUE.
        END WHERE
        
        call hard_law(hard_rate, dhard_rate, crss_tmp, crss_sat, shear, shrate, epseff, 2, ihard, n, m)
        fjac = 1.0 - dtime * dhard_rate
        del_crss = res_n / fjac
        !Find where the newton_ok criterion is not satisfied
        WHERE ((fjac .LT. VTINY) .AND. (.NOT. converged))  newton_ok = .FALSE.
  
        xlam = crss_tmp + del_crss        

        call hard_law(hard_rate, dhard_rate, xlam, crss_sat, shear, shrate, epseff, 1, ihard, n, m)
        
        res = - (xlam - crss_0 - dtime * hard_rate)

        WHERE ((.NOT. converged)) 
          ratio_res = ABS(res) / ABS(res_n)
        ELSEWHERE
          ratio_res = RK_ZERO
        END WHERE
!
!       Line search
!
        fj31 = RK_ONE

        DO WHILE(ANY(ratio_res .GT. RK_ONE .AND. newton_ok .AND. .NOT. converged))

          WHERE (ratio_res .GT. RK_ONE .AND. newton_ok .AND. .NOT. converged) fj31 = fj31*0.5_RK
          !Find where the newton_ok criterion is not satisfied
          WHERE(fj31 .LT. LS_CUTOFF) newton_ok = .FALSE.

          xlam = crss_tmp + fj31 * del_crss
          
          CALL hard_law(hard_rate, dhard_rate, xlam, crss_sat, shear, shrate, epseff, 1, ihard, n, m)
          
          WHERE (ratio_res .GT. RK_ONE .AND. newton_ok .AND. .NOT. converged)
            res = - (xlam - crss_0 - dtime * hard_rate)
            ratio_res = ABS(res) / ABS(res_n)
          END WHERE

        END DO

        WHERE( (newton_ok) .AND. (.NOT.converged) ) crss = crss_tmp + fj31 * del_crss
        !Find where
        do iphase=1,numphases
           call CrystalTypeGet(ctype(iphase))
           n_slip=ctype(iphase)%numslip
           do islip=0,n_slip-1
             where (my_phase .eq. iphase)
                WHERE ( (dabs(res(islip,:,:)) .LT. TOLER_HARD * crystal_parm(3,iphase)) .AND.&
                    & (newton_ok(islip,:,:)) .AND. (.NOT. done) ) converged(islip,:,:) = .TRUE.
            endwhere
           enddo
        enddo

        !Finds out which areas have not converged and which newton_ok are not okay
        !Set
        do iphase=1,numphases
            call CrystalTypeGet(ctype(iphase))
            n_slip=ctype(iphase)%numslip
            do islip=0,n_slip-1
                    !Converged_all set to false only when one slip system is false
                    where(.NOT. converged(islip,:,:)) converged_all=.FALSE.
                    !Newton_ok_all set to false only when one slip system is false
                    where(.NOT. newton_ok(islip,:,:)) newton_ok_all=.FALSE.
            enddo
        enddo

!        Write(*,*) 'Converged_all:',converged_all
!        WRITE(*,*) 'Newton_ok_all:',newton_ok_all

        inewton = COUNT( .NOT. newton_ok_all )
        IF ((COUNT(converged_all) + inewton) .EQ. nm) THEN
          !
          ! all converged or not OK: print message if any not OK, else exit loop
          !
!          IF (inewton .GT. 0) THEN
!            WRITE(ounits(LOG_U), 1000) COUNT(converged), inewton,&
!                 &  MINVAL(dabs(res), MASK=converged),&
!                 &  MAXVAL(dabs(res), MASK=converged)
!          END IF

          EXIT
        ENDIF

      enddo
            
      IF (iter_hard > MAX_ITER_HARD) THEN
        !
        ! NOTE: only processes with this condition will write this message
        !
        WRITE(DFLT_U, *) 'Warning:       . Update hardening iteration limit reached.'
        !message = 'Warning:       . Update hardening iteration limit reached.'
        !CALL par_message(STDERR, message, ALLWRITE=.TRUE.)
      END IF


      ! Should we do this? maybe log the count of not converged
      where(.not. converged) crss = crss_0
      
      RETURN

 1000 FORMAT(' upd_hard: converged =', i6, ' remaining = ', i6 /&
  &   ' min/max residual for converged grains =', g12.5, 2x, g12.5)

      END SUBROUTINE
!
      SUBROUTINE upd_rotn(rstar, rstar_n, dr, time_step, wp_hat, wp_ss, done, n, m) 
!
!------------------------------------------------------------------------
!
      USE DimsModule
      IMPLICIT NONE
!
!
!     Arguments:
!
      INTEGER n, m
      REAL(RK)  time_step
      REAL(RK), INTENT(INOUT) ::  rstar(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      !REAL(RK), INTENT(OUT) ::  rstar(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  ::  rstar_n(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(INOUT) ::  dr(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      !REAL(RK), INTENT(OUT) ::  dr(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))      
      REAL(RK), INTENT(IN)  ::  wp_hat(0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  ::  wp_ss(0:DIMS1, 0:(n - 1), 0:(m - 1))
      LOGICAL, INTENT(IN) ::  done(0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER   ::  i, j, k
      REAL(RK)  ::  th1(0:(n - 1), 0:(m - 1)), th2(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  th3(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  tau(0:(n - 1), 0:(m - 1)), taua(0:(n - 1), 0:(m - 1))
      REAL(RK)  ::  r(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
!
!------------------------------------------------------------------------
!
!     Integrate the evolution equation for r*.
!
      th1 = (wp_hat(0, :, :) - wp_ss(0, :, :)) * time_step
      th2 = (wp_hat(1, :, :) - wp_ss(1, :, :)) * time_step
      th3 = (wp_hat(2, :, :) - wp_ss(2, :, :)) * time_step
    
      tau = sqrt(th1 * th1 + th2 * th2 + th3 * th3)
      where(tau .ne. 0.0)
        taua = tan(tau / 2.0) / tau
        th1 = taua * th1
        th2 = taua * th2
        th3 = taua * th3
        tau = taua * tau
      end where
    
      tau = 2.0 / (1.0 + tau * tau)

      where (.not. done)
    
        dr(0, 0, :, :) = 1.0 - tau * (th1 * th1 + th2 * th2)
        dr(0, 1, :, :) = - tau * (th1 + th2 * th3)
        dr(0, 2, :, :) = tau * ( - th2 + th1 * th3)
        dr(1, 0, :, :) = tau * (th1 - th2 * th3)
        dr(1, 1, :, :) = 1.0 - tau * (th1 * th1 + th3 * th3)
        dr(1, 2, :, :) = - tau * (th3 + th1 * th2)
        dr(2, 0, :, :) = tau * (th2 + th1 * th3)
        dr(2, 1, :, :) = tau * (th3 - th1 * th2)
        dr(2, 2, :, :) = 1.0 - tau * (th2 * th2 + th3 * th3)
    
      endwhere
      
      r = 0.0d0
!  RC 3/24/2016: Reordered for better memory striding
      do j = 0, DIMS1
        do k = 0, DIMS1
          do i = 0, DIMS1
            r(i, j, :, :) = r(i, j, :, :) + rstar_n(i, k, :, :) * dr(k, j, :, :)
          enddo
        enddo
      enddo

      where (.not. done)
    
        rstar(0, 0, :, :) = r(0, 0, :, :) 
        rstar(0, 1, :, :) = r(0, 1, :, :) 
        rstar(0, 2, :, :) = r(0, 2, :, :) 
        rstar(1, 0, :, :) = r(1, 0, :, :) 
        rstar(1, 1, :, :) = r(1, 1, :, :) 
        rstar(1, 2, :, :) = r(1, 2, :, :) 
        rstar(2, 0, :, :) = r(2, 0, :, :) 
        rstar(2, 1, :, :) = r(2, 1, :, :) 
        rstar(2, 2, :, :) = r(2, 2, :, :) 
    
      endwhere

      RETURN
      END SUBROUTINE
!
!**********************************************************************
!
      SUBROUTINE get_c(c, r, c_0, done, n, m)
!
!------------------------------------------------------------------------
!
      USE DimsModule
      IMPLICIT NONE
!
!
!     Arguments:
! 
      INTEGER   ::  n, m
      REAL(RK), INTENT(IN)  ::  r(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m -1))
      REAL(RK), INTENT(IN)  ::  c_0(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m -1))
      REAL(RK), INTENT(OUT) ::  c(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m -1))
      LOGICAL, INTENT(IN) ::  done(0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER   ::  i, j, k
      REAL(RK)  ::  c_aux(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m -1))
!
!------------------------------------------------------------------------
!
!     Compute : c = c_0 * rstar

      c_aux = 0.0_RK
!  RC 3/24/2016: Reordered for better memory striding
      do j = 0, DIMS1
        do k = 0, DIMS1
          do i = 0, DIMS1
            c_aux(i, j, :, :) = c_aux(i, j, :, :) + c_0(i, k, :, :) * r(k, j, :, :)
          enddo
        enddo
      enddo

      where (.not. done)
    
        c(0, 0, :, :) = c_aux(0, 0, :, :) 
        c(0, 1, :, :) = c_aux(0, 1, :, :) 
        c(0, 2, :, :) = c_aux(0, 2, :, :) 
        c(1, 0, :, :) = c_aux(1, 0, :, :) 
        c(1, 1, :, :) = c_aux(1, 1, :, :) 
        c(1, 2, :, :) = c_aux(1, 2, :, :) 
        c(2, 0, :, :) = c_aux(2, 0, :, :) 
        c(2, 1, :, :) = c_aux(2, 1, :, :) 
        c(2, 2, :, :) = c_aux(2, 2, :, :) 
 
      endwhere

      RETURN
      END SUBROUTINE
!
      END MODULE rstarn_solve_lag_mod
!
      
