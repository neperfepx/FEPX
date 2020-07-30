! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE StressSolveEvpsModule

  USE parallel_mod
  USE IntrinsicTypesModule, RK=>REAL_KIND
  USE DimsModule
  USE units_mod

  USE READ_INPUT_MOD
  USE microstructure_mod
  USE UtilsCrystalModule
  USE ConvergenceModule, ONLY: cv_options
  USE StressSolveVpModule

  IMPLICIT  NONE

  PRIVATE
  PUBLIC :: stress_solve_evps, wp_hat_mat5x5_all

CONTAINS

      SUBROUTINE stress_solve_evps(&
  &      sig, d_vec_lat, w_vec_lat, e_bar_vec,&
  &      crss, keinv,&
  &      dtime, wp_hat, iter_state, done,&
  &      converged_newton&
  &      )
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER, INTENT(IN)   :: iter_state      
      LOGICAL, INTENT(INOUT):: converged_newton
      LOGICAL, INTENT(IN)   :: done(0:ngrain1, el_sub1:el_sup1)
!
      REAL(RK), INTENT(IN)    :: dtime
      REAL(RK), INTENT(IN)    :: keinv(0:TVEC1,1:numphases)
      REAL(RK), INTENT(INOUT) :: sig(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: d_vec_lat(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: w_vec_lat(0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: e_bar_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: crss(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(OUT)   :: wp_hat(0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
!
!     Locals:
!
      INTEGER :: ier , num, m_el
      INTEGER ::  jiter(0:ngrain1, el_sub1:el_sup1)
!
      REAL(RK) :: sig0_avg, sig1_avg, sig2_avg, sig3_avg, sig4_avg
!
      LOGICAL  :: converged(0:ngrain1, el_sub1:el_sup1)
!
!----------------------------------------------------------------------
!
      m_el = el_sup1 - el_sub1 + 1
!
      converged = .false.
      where(done) converged = .true.
!
!     Solve NLE for stresses.
!
!      write(6,*)'         STRESS_SOLVE_EVPS: entering SOLVE_NEWTON_EVPS'

      call solve_newton_evps(sig, d_vec_lat, e_bar_vec, w_vec_lat,&
  &        crss, cv_options % sx_max_iters_newton, ier,&
  &        keinv, dtime, wp_hat, jiter,&
  &        converged, done, ngrain, m_el)
!
      if (ier .eq. 0) GO TO 50

!      PRINT 100, count(.not. converged), iter_state
!100   FORMAT('stress_solve_evps: ', i7,&
!  &   ' grains did not converge in iter_stat =', i4 /&
!  &   ' These grains will receive these average stress values:' /&
!  &   ' (correct for constant time step, i.e. auto_time = 0)' )

      WRITE(DFLT_U, '(A,I0,A,I0)') 'Warning:       . STRESS_SOLVE_EVPS: ', &
        & count(.not. converged), ' grains did not converge in iter_state = ',&
        & iter_state
      WRITE(DFLT_U, '(A)') 'Warning:       . These grains will receive &
        &average stress values (correct for constant time step).'

      num = count(converged)
!deb  Put in to avoid infinities.
      if (num .eq. 0) num = 1
!
      sig0_avg = sum(sig(0, :, :), mask = converged) / num
      sig1_avg = sum(sig(1, :, :), mask = converged) / num
      sig2_avg = sum(sig(2, :, :), mask = converged) / num
      sig3_avg = sum(sig(3, :, :), mask = converged) / num
      sig4_avg = sum(sig(4, :, :), mask = converged) / num

!      PRINT 150, sig0_avg, sig1_avg, sig2_avg, sig3_avg, sig4_avg 
!150   FORMAT(5x, 5(e11.4, 2x))

      where (.not. converged) 
         sig(0, :, :) = sig0_avg
         sig(1, :, :) = sig1_avg
         sig(2, :, :) = sig2_avg
         sig(3, :, :) = sig3_avg
         sig(4, :, :) = sig4_avg
      endwhere

      converged_newton = .false.

50    continue

      END SUBROUTINE stress_solve_evps
!
!**********************************************************************
!
      SUBROUTINE solve_newton_evps(&
  &      sig, d_vec_lat, e_bar_vec, w_vec_lat,&
  &      crss, max_iter_newton, irc, keinv,&
  &      dt, wp_hat_vec, jiter, converged, done,&
  &      n, m&
  &      )
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER, INTENT(IN)   :: max_iter_newton, n, m
      INTEGER, INTENT(OUT)  :: jiter(0:(n - 1), 0:(m - 1))
!
      INTEGER, INTENT(OUT)  :: irc
      LOGICAL, INTENT(INOUT):: converged(0:(n - 1), 0:(m - 1))
      LOGICAL, INTENT(IN)   :: done(0:(n - 1), 0:(m - 1))
!
      REAL(RK), INTENT(IN)    :: dt
      REAL(RK), INTENT(IN)    :: keinv(0:TVEC1, 1:numphases)
      REAL(RK), INTENT(INOUT) :: sig(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    :: d_vec_lat(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    :: e_bar_vec(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    :: w_vec_lat(0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    :: crss(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT)   :: wp_hat_vec(0:DIMS1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      LOGICAL  :: newton_ok(0:(n - 1), 0:(m - 1))
!
      INTEGER  :: iter_newton, islip, i, nm, inewton, n_slip
!
      REAL(RK) :: c1, xn(numphases), xnn(numphases)
      REAL(RK) :: sig_0(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK) :: e_bar(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK) :: e_vec(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK) :: e_elas(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
               
      REAL(RK) :: gdot(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK) :: dgdot(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK) :: dp_hat_vec(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK) :: wp_x_e(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK) :: del_s(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK) :: fj(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK) :: res(0:(n - 1), 0:(m - 1)), res_n(0:(n - 1), 0:(m - 1))
      REAL(RK) :: fact(0:(n - 1), 0:(m - 1))
      REAL(RK) :: xlambda(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK) :: fjac(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
               
               
      REAL(RK) :: wp_hat_matx(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK) :: tau(0:(n - 1), 0:(m - 1)), taua(0:(n - 1), 0:(m - 1))
      REAL(RK) :: fj21(0:(n - 1), 0:(m - 1))
      REAL(RK) :: ratio_res(0:(n - 1), 0:(m - 1))
      REAL(RK) :: res_aux(0:(n - 1), 0:(m - 1))

      REAL(RK), POINTER :: p_hat_vec(:,:) => NULL()
      REAL(RK), POINTER :: ppt(:,:,:) => NULL()
      REAL(RK), pointer :: e_vec_tmp(:,:,:) => NULL()
      REAL(RK), pointer :: e_elas_tmp(:,:,:,:) => NULL()
      REAL(RK), pointer :: wp_x_e_tmp(:,:,:) => NULL()
      REAL(RK), pointer :: fjac_tmp(:,:,:,:) => NULL()
      
            
      REAL(RK), PARAMETER :: ONE_DP=1.0_RK
!      REAL(RK), PARAMETER :: ONE_DP=1.0d0      
      REAL(RK) :: UFLOW = TINY(ONE_DP)
      REAL(RK) :: toosmall(numphases)

      integer  :: my_phase(0:(m - 1)), iphase, numind, in
      integer, pointer :: indices(:) => NULL()

!----------------------------------------------------------------------
!
      !
      my_phase(:) = phase(el_sub1:el_sup1)
      ! 
      irc   = 0
      jiter = 0
      sig_0 = 0.0_RK
      newton_ok = .true.

      nm  = n * m
      c1  = 1.0_RK / dt
!      c1  = 1.d0 / dt

      ! e_bar_vec {5} --> e_bar [3x3]sym
      call vec_mat_symm_grn(e_bar_vec, e_bar, n, m)

      !slw
      ! g95 compiler gives random values of gdot and dgdot without the initialization
      gdot = 0.0_RK
      dgdot = 0.0_RK

!      
!---- start iterations -----------------------------------------------
!
      do 100 iter_newton = 1, max_iter_newton
        
         sig_0 = sig
!
!        Elastic strains.
!
         do iphase=1,numphases
           
            call find_indices(numind, iphase, my_phase, indices)
            call CrystalTypeGet(ctype(iphase), DEV=p_hat_vec,PPTRANS=ppt)
            n_slip=ctype(iphase)%numslip
         
            if(associated(e_vec_tmp)) then
               deallocate(e_vec_tmp)
            endif
            allocate(e_vec_tmp(0:TVEC1,0:(n-1),0:(numind-1)))
            call vec_d_vec5(keinv(:,iphase), sig(:,:,indices),e_vec_tmp, n, numind)
            e_vec(:,:,indices)=e_vec_tmp
         
            if(associated(e_elas_tmp)) then
               deallocate(e_elas_tmp)
            endif
            allocate(e_elas_tmp(0:DIMS1,0:DIMS1,0:(n-1),0:(numind-1)))
            call vec_mat_symm_grn(e_vec(:,:,indices),e_elas_tmp, n, numind)
            e_elas(:,:,:,indices)=e_elas_tmp
         
!           Power law viscoplastic model.
!        
            xnn(iphase) = 1.0_RK / crystal_parm(0,iphase)
            xn(iphase)  = xnn(iphase) - 1.0_RK
            toosmall(iphase) = UFLOW**crystal_parm(0,iphase)
         
            do islip = 0, n_slip - 1          
               tau(:,indices) = 0.0_RK
               do  i = 0, TVEC1
                  tau(:,indices)=tau(:,indices)+p_hat_vec(i+1,islip+1)*sig(i,:,indices)/crss(islip,:,indices)
               enddo 
               
               taua(:,indices) = dabs(tau(:,indices))!     
               WHERE (taua(:,indices) .LE. toosmall(iphase))  taua(:,indices) = 0.0_RK
                
               fj21(:,indices) = crystal_parm(1,iphase) * taua(:,indices)**xn(iphase)  
               dgdot(islip, :,indices) = fj21(:,indices) * xnn(iphase) / crss(islip,:,indices)
               gdot (islip, :, indices) = fj21(:,indices) * tau(:,indices) 
            enddo !n_slip
         
!           Set up the jacobian.
!        
            call dp_wp_hat(p_hat_vec, dp_hat_vec,&
  &              wp_hat_vec, e_elas, e_bar, w_vec_lat, gdot,&
  &              n_slip, dt, n, m, numind, indices)
         
            call wp_hat_mat5x5(wp_hat_vec,wp_hat_matx, n, m, numind, indices)
         
            if(associated(wp_x_e_tmp)) then
               deallocate(wp_x_e_tmp)
            endif
            allocate(wp_x_e_tmp(0:TVEC1, 0:(n - 1), 0:(numind - 1)))
            call mat_x_vec5(wp_hat_matx(:, :, :, indices),e_vec(:, :, indices), wp_x_e_tmp, n, numind)
            wp_x_e(:, :, indices)=wp_x_e_tmp
         
            if(associated(fjac_tmp)) then
               deallocate(fjac_tmp)
            endif
            allocate(fjac_tmp(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(numind-1)))
         
            call matrix_fjac(fjac_tmp,&
  &              keinv(:,iphase), wp_hat_matx(:, :, :, indices), c1,&
  &              dgdot(:, :, indices), ppt, n_slip,&
  &              n, numind)
            fjac(:, :, :, indices)=fjac_tmp
!        
!           Set up the system function (RHS)      
         
            call residual(res_n, del_s,&
  &              d_vec_lat, e_vec, e_bar_vec,&
  &              dp_hat_vec, wp_x_e, c1, n, m, numind, indices)
            !
            deallocate(indices)
            deallocate(fjac_tmp)
            deallocate(ppt)
            deallocate(wp_x_e_tmp)
            deallocate(e_vec_tmp)
            deallocate(e_elas_tmp)
            deallocate(p_hat_vec)
            !
         enddo !numphases

!        Compute new iterate

         call symmetrize_jac(fjac, del_s, n, m)
         call solvit(fjac, del_s, n, m)
         call check_diagonals_evps(fjac, newton_ok, done, n, m)

!        Trial residual.

         xlambda = sig_0 + del_s

!
         do iphase=1,numphases
            call CrystalTypeGet(ctype(iphase), DEV=p_hat_vec)
            call find_indices(numind, iphase, my_phase, indices)
            n_slip=ctype(iphase)%numslip
            if(associated(e_vec_tmp)) then
               deallocate(e_vec_tmp)
            endif
            allocate(e_vec_tmp(0:TVEC1,0:(n-1),0:(numind-1)))
            call vec_d_vec5(keinv(:,iphase), xlambda(:,:,indices),e_vec_tmp, n, numind)
            e_vec(:,:,indices)=e_vec_tmp

            if(associated(e_elas_tmp)) then
               deallocate(e_elas_tmp)
            endif
            allocate(e_elas_tmp(0:DIMS1,0:DIMS1,0:(n-1),0:(numind-1)))
            call vec_mat_symm_grn(e_vec(:, :, indices),e_elas_tmp, n, numind)
            e_elas(:,:,:,indices)=e_elas_tmp
            
            do islip = 0, n_slip - 1
               tau(:,indices) = 0.0_RK
               do  i = 0, TVEC1
                   tau(:,indices)=tau(:,indices)+p_hat_vec(i+1,islip+1)*xlambda(i,:,indices)/crss(islip,:,indices)
               enddo                
               taua(:,indices) = dabs(tau(:,indices))
               do in=0,(n-1)
                  WHERE (taua(in,indices) .LE. toosmall(iphase)) taua(in,indices) = 0.0_RK
                  gdot(islip, in, :) = crystal_parm(1,iphase)*tau(in,:)*taua(in,:)**xn(iphase) 
               enddo 
            enddo !n_slip
                
            call dp_wp_hat(p_hat_vec, dp_hat_vec,&
  &              wp_hat_vec, e_elas, e_bar, w_vec_lat, gdot,&
  &              n_slip, dt, n, m, numind, indices)

            call wp_hat_mat5x5(wp_hat_vec,wp_hat_matx, n, m, numind, indices)

            if(associated(wp_x_e_tmp)) then
               deallocate(wp_x_e_tmp)
            endif
            allocate(wp_x_e_tmp(0:TVEC1, 0:(n - 1), 0:(numind - 1)))
            call mat_x_vec5(wp_hat_matx(:, :, :, indices),e_vec(:, :, indices), wp_x_e_tmp, n, numind)
            wp_x_e(:, :, indices)=wp_x_e_tmp
          
            call residual(res, fj,&
  &              d_vec_lat, e_vec, e_bar_vec,&
  &              dp_hat_vec, wp_x_e, c1, n, m, numind, indices)
  
            !
            deallocate(indices)
            deallocate(wp_x_e_tmp)
            deallocate(e_vec_tmp)
            deallocate(e_elas_tmp)
            deallocate(p_hat_vec)
            !         
         enddo !numphases
!
!        Line Search.
!
         fact = 1.0_RK
!         fact = 1.0d0
         ratio_res = res / res_n

         do while(any(&
  &               ratio_res .gt. 1.0&
  &               .and. newton_ok&
  &               .and. .not. converged ))

            where(ratio_res .gt. 1.0 .and. newton_ok .and. .not. converged)  fact = fact*0.5_RK

            if( any(fact .lt. 0.001) ) then
               print *, ' error in line search',count(fact .lt. 0.1), iter_newton
               where(fact .lt. 0.001) newton_ok = .false.
            endif

            xlambda(0, :, :) = sig_0(0, :, :) + fact * del_s(0, :, :) 
            xlambda(1, :, :) = sig_0(1, :, :) + fact * del_s(1, :, :) 
            xlambda(2, :, :) = sig_0(2, :, :) + fact * del_s(2, :, :) 
            xlambda(3, :, :) = sig_0(3, :, :) + fact * del_s(3, :, :) 
            xlambda(4, :, :) = sig_0(4, :, :) + fact * del_s(4, :, :) 
         
!
            do iphase=1,numphases
               call CrystalTypeGet(ctype(iphase), DEV=p_hat_vec)
               call find_indices(numind, iphase, my_phase, indices)
               n_slip=ctype(iphase)%numslip
               if(associated(e_vec_tmp)) then
                  deallocate(e_vec_tmp)
               endif
               allocate(e_vec_tmp(0:TVEC1,0:(n-1),0:(numind-1)))
               call vec_d_vec5(keinv(:,iphase), xlambda(:,:,indices), e_vec_tmp, n, numind)
               e_vec(:,:,indices)=e_vec_tmp

               if(associated(e_elas_tmp)) then
                  deallocate(e_elas_tmp)
               endif
               allocate(e_elas_tmp(0:DIMS1,0:DIMS1,0:(n-1),0:(numind-1)))
               call vec_mat_symm_grn(e_vec(:,:,indices),e_elas_tmp, n,numind)
               e_elas(:,:,:,indices)=e_elas_tmp
               
               do islip = 0, n_slip - 1
                  tau(:,indices) = 0.0_RK
                  do  i = 0, TVEC1
                     tau(:,indices)=tau(:,indices)+p_hat_vec(i+1,islip+1)*&
  &                                 xlambda(i, :, indices) / crss(islip,:,indices)
                  enddo 
                  taua(:,indices) = dabs(tau(:,indices))
                  do in=0,n-1
                
!-tm seahag doesn't like this where statement               
!                     where (my_phase .eq. iphase)
                        WHERE (taua(in,indices) .LE. toosmall(iphase)) taua(in,indices) = 0.0
                        gdot(islip, in, :) = crystal_parm(1,iphase)*tau(in,:)*taua(in,:)**xn(iphase) 
!                     endwhere
                  enddo
               enddo !n_slip
               
               call dp_wp_hat(p_hat_vec, dp_hat_vec,&
  &                 wp_hat_vec, e_elas, e_bar, w_vec_lat, gdot,&
  &                 n_slip, dt, n, m, numind, indices)

               call wp_hat_mat5x5(wp_hat_vec,wp_hat_matx, n, m, numind, indices)

!
               if(associated(wp_x_e_tmp)) then
                  deallocate(wp_x_e_tmp)
               endif
               allocate(wp_x_e_tmp(0:TVEC1, 0:(n - 1), 0:(numind - 1)))
               call mat_x_vec5(wp_hat_matx(:, :, :, indices),e_vec(:, :, indices), wp_x_e_tmp, n, numind)
               wp_x_e(:, :, indices)=wp_x_e_tmp
!
               call residual(res_aux, fj,&
  &                 d_vec_lat, e_vec, e_bar_vec, dp_hat_vec, wp_x_e,&
  &                 c1, n, m, numind, indices)
               !
               deallocate(indices)
               deallocate(wp_x_e_tmp)
               deallocate(e_vec_tmp)
               deallocate(e_elas_tmp)
               deallocate(p_hat_vec)
               !
            enddo !numphases

            where((ratio_res .gt. 1.0) .and. (newton_ok) .and. (.not. converged))
               res = res_aux
               ratio_res = res / res_n
            endwhere

        enddo !do while
!
!       Update stresses and check convergence.
!
        where( (newton_ok) .and. (.not. converged) )
           sig(0, :, :) = sig_0(0, :, :) + fact * del_s(0, :, :)
           sig(1, :, :) = sig_0(1, :, :) + fact * del_s(1, :, :)
           sig(2, :, :) = sig_0(2, :, :) + fact * del_s(2, :, :)
           sig(3, :, :) = sig_0(3, :, :) + fact * del_s(3, :, :)
           sig(4, :, :) = sig_0(4, :, :) + fact * del_s(4, :, :)
        endwhere

!
!deb    The following section was rewritten because the AIX compiler does not
!       like the nested `where' statements, even though they
!       seem to be syntactically correct.
 
        where( (res .le. cv_options % sx_tol).and.(newton_ok).and.(.not. done).and.(.not. converged) ) jiter = iter_newton
        where( (res .le. cv_options % sx_tol).and.(newton_ok).and.(.not. done) )    converged = .true.
!
!deb    End of rewritten segment.
!
!         Return if all grains have converged.

        inewton = count(.not. newton_ok)
        if( (count(converged) + inewton) .eq. nm ) then
           if (inewton .gt. 0)&
  &           PRINT 1000, count(converged), inewton,&
  &           minval(res, mask=converged),&
  &           maxval(res, mask=converged)
              irc = inewton
           RETURN
        endif
 
  100 CONTINUE  ! -------- stop iterations --------------------------

      irc = -2

 1000 FORMAT('solve_newton_evps: converged =', i7, ' remaining =',i7/&
  &       'residual for converged grains =', g12.5, 2x, g12.5)

      END SUBROUTINE solve_newton_evps
!
!**********************************************************************
!
      SUBROUTINE matrix_fjac(&
  &      fjac, Keinv, w_matx, dti, dgdot, ppt, n_slip, n, m&
  &      )
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER, INTENT(IN)   :: n_slip, n, m
!
      REAL(RK), INTENT(OUT) :: fjac(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: ppt(0:TVEC1, 0:TVEC1, 0:MAXSLIP1)
      REAL(RK), INTENT(IN)  :: dti
      REAL(RK), INTENT(IN)  :: Keinv(0:TVEC1)
      REAL(RK), INTENT(IN)  :: w_matx(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: dgdot(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER :: islip, i, j, in, im
!
!----------------------------------------------------------------------
!
      DO im = 0, m-1
         DO in = 0, n-1
!
            DO j = 0, TVEC1
               DO i = 0, TVEC1          
                    fjac(i, j, in, im) = w_matx(i, j, in, im) * Keinv(j)
               ENDDO
               fjac(j, j, in, im) = fjac(j, j, in, im) + Keinv(j) * dti
            ENDDO
!
            DO islip = 0, n_slip-1
               fjac(:,:,in,im) = fjac(:,:,in,im) + dgdot(islip,in,im) * ppt(:,:,islip)
            ENDDO
!
         ENDDO
      ENDDO
!
      END SUBROUTINE matrix_fjac
!
!**********************************************************************
!
      SUBROUTINE jacob5x5(jac, a, n, m)
        
        ! doesn't appear to be used !
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER, INTENT(IN)   :: n, m
!
      REAL(RK), INTENT(IN)  :: a(0:DIMS1, 0:DIMS1, 0:DIMS1, 0:DIMS1, 0:(n-1), 0:(m-1))
      REAL(RK), INTENT(OUT) :: jac(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      REAL(RK) :: sqr3, sqr3b2
!
!----------------------------------------------------------------------
!
      sqr3   = dsqrt(3.d0)
      sqr3b2 = sqr3 / 2.d0 

      jac(0, 0, :, :) = (a(0, 0, 0, 0, :, :) + a(1, 1, 1, 1, :, :) -&
  &                   a(0, 0, 1, 1, :, :) - a(1, 1, 0, 0, :, :)) / 2.
      jac(0, 1, :, :) = sqr3b2 * (a(0, 0, 2, 2, :, :) -&
  &                                             a(1, 1, 2, 2, :, :))
      jac(0, 2, :, :) = a(0, 0, 0, 1, :, :) - a(1, 1, 0, 1, :, :)
      jac(0, 3, :, :) = a(0, 0, 0, 2, :, :) - a(1, 1, 0, 2, :, :)
      jac(0, 4, :, :) = a(0, 0, 1, 2, :, :) - a(1, 1, 1, 2, :, :)
      
      jac(1, 0, :, :) = sqr3b2 * (a(2, 2, 0, 0, :, :) -&
  &                                             a(2, 2, 1, 1, :, :))
      jac(1, 1, :, :) = 1.5 * a(2, 2, 2, 2, :, :) -&
  &   0.5*(a(2, 2, 0, 0, :, :) + a(2, 2, 1, 1, :, :) + a(2, 2, 2, 2, :, :))
      jac(1, 2, :, :) = sqr3 * a(2, 2, 0, 1, :, :)
      jac(1, 3, :, :) = sqr3 * a(2, 2, 0, 2, :, :)
      jac(1, 4, :, :) = sqr3 * a(2, 2, 1, 2, :, :)

      jac(2, 0, :, :) = a(0, 1, 0, 0, :, :) - a(0, 1, 1, 1, :, :)
      jac(2, 1, :, :) = sqr3 * a(0, 1, 2, 2, :, :)
      jac(2, 2, :, :) = 2. * a(0, 1, 0, 1, :, :)
      jac(2, 3, :, :) = 2. * a(0, 1, 0, 2, :, :)
      jac(2, 4, :, :) = 2. * a(0, 1, 1, 2, :, :)

      jac(3, 0, :, :) = a(0, 2, 0, 0, :, :) - a(0, 2, 1, 1, :, :)
      jac(3, 1, :, :) = sqr3 * a(0, 2, 2, 2, :, :)
      jac(3, 2, :, :) = 2. * a(0, 2, 0, 1, :, :)
      jac(3, 3, :, :) = 2. * a(0, 2, 0, 2, :, :)
      jac(3, 4, :, :) = 2. * a(0, 2, 1, 2, :, :)

      jac(4, 0, :, :) = a(1, 2, 0, 0, :, :) - a(1, 2, 1, 1, :, :)
      jac(4, 1, :, :) = sqr3 * a(1, 2, 2, 2, :, :)
      jac(4, 2, :, :) = 2. * a(1, 2, 0, 1, :, :)
      jac(4, 3, :, :) = 2. * a(1, 2, 0, 2, :, :)
      jac(4, 4, :, :) = 2. * a(1, 2, 1, 2, :, :)

      END SUBROUTINE jacob5x5
!
!**********************************************************************
!
      SUBROUTINE dp_wp_hat(&
  &       p_hat_vec, dp_hat, wp_hat, e_elas, e_bar, w_vec_lat,&
  &       gdot, n_slip, dt, n, m, numind, indices&
  &       )
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER, INTENT(IN)   :: n_slip, n, m, numind, indices(1:numind)
      REAL(RK), INTENT(OUT) :: dp_hat(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT) :: wp_hat(0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: p_hat_vec(0:TVEC1,0:MAXSLIP1)
      REAL(RK), INTENT(IN)  :: dt
      REAL(RK), INTENT(IN)  :: e_elas(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: e_bar(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: w_vec_lat(0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: gdot(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER  :: i, islip
      REAL(RK) :: dp_hat_tmp(0:TVEC1, 0:(n - 1), 0:(numind - 1))
      REAL(RK) :: wp_hat_tmp(0:DIMS1, 0:(n - 1), 0:(numind - 1))
      REAL(RK) :: e_elas_tmp(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(numind - 1))
      REAL(RK) :: e_bar_tmp(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(numind - 1))
      REAL(RK) :: w_vec_lat_tmp(0:DIMS1, 0:(n - 1), 0:(numind - 1))
      REAL(RK) :: gdot_tmp(0:MAXSLIP1, 0:(n - 1), 0:(numind - 1))
      REAL(RK) :: p_hat(0:DIMS1, 0:DIMS1, 0:MAXSLIP1)
      REAL(RK) :: x (0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(numind - 1))
      REAL(RK) :: ee(0:DIMS1, 0:DIMS1, 0:(n - 1), 0:(numind - 1))
!
!----------------------------------------------------------------------
!
               
      e_elas_tmp=e_elas(:, :, :, indices)
      e_bar_tmp=e_bar(:, :, :, indices)
      w_vec_lat_tmp=w_vec_lat(:, :, indices)
      gdot_tmp=gdot(:, :, indices)
      call vec_mat_symm(p_hat_vec, p_hat, n_slip)
      dp_hat_tmp = 0.0_RK

      call mat_x_mat3(e_elas_tmp, e_bar_tmp, ee, n, numind)

      wp_hat_tmp(0, :, :) = w_vec_lat_tmp(0, :, :) +&
  &                     0.5 / dt * (ee(1, 0, :, :) - ee(0, 1, :, :))
      wp_hat_tmp(1, :, :) = w_vec_lat_tmp(1, :, :) +&
  &                     0.5 / dt * (ee(2, 0, :, :) - ee(0, 2, :, :))
      wp_hat_tmp(2, :, :) = w_vec_lat_tmp(2, :, :) +&
  &                     0.5 / dt * (ee(2, 1, :, :) - ee(1, 2, :, :))

      do islip = 0, n_slip - 1

         call mat_x_mats3(e_elas_tmp, p_hat(0, 0, islip), x, n, numind)

         wp_hat_tmp(0, :, :) = wp_hat_tmp(0, :, :) -&
  &         gdot_tmp(islip, :, :) * (x(1, 0, :, :) - x(0, 1, :, :))
         wp_hat_tmp(1, :, :) = wp_hat_tmp(1, :, :) -&
  &         gdot_tmp(islip, :, :) * (x(2, 0, :, :) - x(0, 2, :, :))
         wp_hat_tmp(2, :, :) = wp_hat_tmp(2, :, :) -&
  &         gdot_tmp(islip, :, :) * (x(2, 1, :, :) - x(1, 2, :, :))

         do i = 0, TVEC1
            dp_hat_tmp(i, :, :) = dp_hat_tmp(i, :, :) +&
  &               gdot_tmp(islip, :, :) * p_hat_vec(i, islip)
         enddo

      enddo
      dp_hat(:,:,indices) = dp_hat_tmp
      wp_hat(:, :, indices) = wp_hat_tmp

      END SUBROUTINE dp_wp_hat
!
!**********************************************************************
!
      SUBROUTINE wp_hat_mat5x5(wp_hat, wp_hat_matx, n, m, numind, indices)
!
!----------------------------------------------------------------------
!
      USE  DimsModule
      USE IntrinsicTypesModule, RK=>REAL_KIND
      !
      IMPLICIT  NONE
!
!
!     Arguments:
!
      INTEGER, INTENT(IN)   :: n, m, numind, indices(1:numind)
      REAL(RK), INTENT(IN)  :: wp_hat(0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT) :: wp_hat_matx(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER :: i, j
      REAL(RK) :: sqr3
!
!----------------------------------------------------------------------
      
      sqr3 = dsqrt(3.d0)

      wp_hat_matx(:,:,:,indices) = 0.0

      wp_hat_matx(0, 2, :, indices) =   wp_hat(0, :, indices) * 2.0d0
      wp_hat_matx(0, 3, :, indices) =   wp_hat(1, :, indices)
      wp_hat_matx(0, 4, :, indices) = - wp_hat(2, :, indices)
                                              
      wp_hat_matx(1, 3, :, indices) = - wp_hat(1, :, indices) * sqr3
      wp_hat_matx(1, 4, :, indices) = - wp_hat(2, :, indices) * sqr3
                                              
      wp_hat_matx(2, 0, :, indices) = - wp_hat(0, :, indices) * 2.0d0
      wp_hat_matx(2, 3, :, indices) =   wp_hat(2, :, indices)
      wp_hat_matx(2, 4, :, indices) =   wp_hat(1, :, indices)
                                              
      wp_hat_matx(3, 0, :, indices) = - wp_hat(1, :, indices)
      wp_hat_matx(3, 1, :, indices) =   wp_hat(1, :, indices) * sqr3
      wp_hat_matx(3, 2, :, indices) = - wp_hat(2, :, indices)
      wp_hat_matx(3, 4, :, indices) =   wp_hat(0, :, indices)
                                              
      wp_hat_matx(4, 0, :, indices) =   wp_hat(2, :, indices)
      wp_hat_matx(4, 1, :, indices) =   wp_hat(2, :, indices) * sqr3
      wp_hat_matx(4, 2, :, indices) = - wp_hat(1, :, indices)
      wp_hat_matx(4, 3, :, indices) = - wp_hat(0, :, indices)
     
      END SUBROUTINE wp_hat_mat5x5
!
!**********************************************************************
!
      SUBROUTINE wp_hat_mat5x5_all(wp_hat, wp_hat_matx, n, m)
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER, INTENT(IN)   :: n, m
      REAL(RK), INTENT(IN)  :: wp_hat(0:DIMS1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT) :: wp_hat_matx(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER :: i, j
      REAL(RK) :: sqr3
!
!----------------------------------------------------------------------
!
      
      sqr3 = dsqrt(3.d0)

      wp_hat_matx(:,:,:,:) = 0.0

      wp_hat_matx(0, 2, :, :) =   wp_hat(0, :, :) * 2.0d0
      wp_hat_matx(0, 3, :, :) =   wp_hat(1, :, :)
      wp_hat_matx(0, 4, :, :) = - wp_hat(2, :, :)
                                              
      wp_hat_matx(1, 3, :, :) = - wp_hat(1, :, :) * sqr3
      wp_hat_matx(1, 4, :, :) = - wp_hat(2, :, :) * sqr3
                                              
      wp_hat_matx(2, 0, :, :) = - wp_hat(0, :, :) * 2.0d0
      wp_hat_matx(2, 3, :, :) =   wp_hat(2, :, :)
      wp_hat_matx(2, 4, :, :) =   wp_hat(1, :, :)
                                              
      wp_hat_matx(3, 0, :, :) = - wp_hat(1, :, :)
      wp_hat_matx(3, 1, :, :) =   wp_hat(1, :, :) * sqr3
      wp_hat_matx(3, 2, :, :) = - wp_hat(2, :, :)
      wp_hat_matx(3, 4, :, :) =   wp_hat(0, :, :)
                                              
      wp_hat_matx(4, 0, :, :) =   wp_hat(2, :, :)
      wp_hat_matx(4, 1, :, :) =   wp_hat(2, :, :) * sqr3
      wp_hat_matx(4, 2, :, :) = - wp_hat(1, :, :)
      wp_hat_matx(4, 3, :, :) = - wp_hat(0, :, :)
     
      END SUBROUTINE wp_hat_mat5x5_all
!      
!**********************************************************************
!
      SUBROUTINE residual(&
  &      res, rhs, d_vec_lat, e_vec, e_bar_vec, dp_hat_vec, wp_x_e,&
  &      c1, n, m, numind, indices)
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER, INTENT(IN)   :: n, m, numind, indices(1:numind)
!
      REAL(RK), INTENT(OUT)  :: rhs(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT)  :: res(0:(n - 1), 0:(numind - 1))
      REAL(RK), INTENT(IN)   :: c1
      REAL(RK), INTENT(IN)   :: d_vec_lat(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)   :: e_vec(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)   :: e_bar_vec(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)   :: dp_hat_vec(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)   :: wp_x_e(0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER :: i
!
!----------------------------------------------------------------------
!
      rhs = d_vec_lat - c1 * (e_vec - e_bar_vec) - dp_hat_vec - wp_x_e

      res(:,indices) = 0.0_RK
      do i = 0, TVEC1
         res(:,indices) = res(:,indices) + rhs(i, :, indices)* rhs(i, :, indices)
      enddo
      res(:,indices) = dsqrt(res(:,indices))

      END SUBROUTINE residual
!
!
!**********************************************************************
!
      SUBROUTINE check_diagonals_evps(stif, newton_ok, done, n, m)
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER, INTENT(IN)  :: n, m
!
      LOGICAL, INTENT(OUT) :: newton_ok(0:(n - 1), 0:(m - 1))
      LOGICAL, INTENT(IN)  :: done(0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN) :: stif(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER :: i
!
!----------------------------------------------------------------------
!
      do i = 0, TVEC1
         where ( (abs(stif(i, i, :, :)) .lt. VTINY) .and. (.not. done) )  newton_ok = .false.
      end do
      
      END SUBROUTINE check_diagonals_evps
!
!
!**********************************************************************
!
      SUBROUTINE symmetrize_jac(fjac, del_s, n, m)
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER, INTENT(IN) :: n, m
!
      REAL(RK), INTENT(INOUT) :: fjac(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(INOUT) :: del_s(0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER :: i, j, k
!
      REAL(RK) :: tempj(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK) :: temps(0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!----------------------------------------------------------------------
!
      tempj = 0.0_RK
      temps = 0.0_RK
!  RC 6/24/2016: Reordered loops for better memory striding      
      do j = 0, TVEC1
         do i = 0, TVEC1
            do k = 0, TVEC1
               tempj(i, j, :, :) = tempj(i, j, :, :) + fjac(k, i, :, :) * fjac(k, j, :, :)
            enddo
         enddo
      enddo

      do i = 0, TVEC1
         do j = 0, TVEC1
            temps(i, :, :) = temps(i, :, :) + fjac(j, i, :, :) * del_s(j, :, :)
         enddo
      enddo

      fjac = tempj
      del_s = temps
  
      END SUBROUTINE symmetrize_jac
      
END MODULE StressSolveEvpsModule
