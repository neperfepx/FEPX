! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE StressSolveVpModule

  !     Routines for solving visco-plastic crystal stress equations.
  !
  !     stress_solve_vp: Driver routine which takes care of scaling and 
  !     *              : initial guesses.

  !     scale_down_defr: Rescale the deformation rate to unit size.

  !     compute_work: Compute virtual plastic work for array of deformation rates
  !     *           : applied to vertex stresses.

  !     find_vertex: Select the vertex which maximizes the plastic work.

  !     vertex_stress: Set initial guess (vertex stress) for nonlinear solver.

  !     scale_stress: Rescale stress initial guess according to hardness.

  !     ss_project: Compute inner product of array of tensors with a fixed tensor.

  !     solve_newton_vp: Nonlinear solution of viscoplastic crystal stress equations.

  !     get_res: Compute residual for nonlinear VP crystal stress equation.

  !     form_crystif: Form single crystal stiffness matrix.

  !     power_law: Power law for VP single crystal.

  !     compliance: Form crystal compliance matrices.

  !     solvit: Solve an array of symmetric positive definite 5x5 systems.

  !     check_diagonals: Determine where diagonal elements are small.

  !     scale_up_sigm: Rescale the stress after solution is found.

  USE parallel_mod
  USE IntrinsicTypesModule, RK=>REAL_KIND
  USE READ_INPUT_MOD
  USE units_mod
  USE DimsModule
  USE microstructure_mod
  USE ConvergenceModule, ONLY: cv_options
  USE UtilsCrystalModule

  IMPLICIT  NONE

  PRIVATE
  PUBLIC :: stress_solve_vp, ss_project, power_law, compliance, solvit, &
       &    check_diagonals

CONTAINS
!
!**********************************************************************
!
      SUBROUTINE stress_solve_vp(sig, d_vec_lat, crss, epseff, vp_log)
!
!     Driver routine which takes care of scaling and initial guesses.
!
!----------------------------------------------------------------------
!
      REAL(RK), INTENT(OUT)   :: sig(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) :: d_vec_lat(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: crss  (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: epseff(0:ngrain1, el_sub1:el_sup1)
      LOGICAL,  INTENT(IN)    :: vp_log
!
!     Locals:
!
      LOGICAL  :: converged(0:ngrain1, el_sub1:el_sup1)
!
      INTEGER  :: i_edge, ier, m_el, ntrials
      INTEGER  :: vertex(0:ngrain1, el_sub1:el_sup1)
!
      REAL(RK) :: direction(0:ngrain1, el_sub1:el_sup1)
      REAL(RK) :: plwork(0:MAX_VERT1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK) :: sig_t(0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK) :: crss_avg  (0:ngrain1, el_sub1:el_sup1)
!
!----------------------------------------------------------------------
!
!
      m_el = el_sup1 - el_sub1 + 1
      crss_avg = 0.0_RK

      call scale_down_defr(d_vec_lat, epseff, ngrain, m_el)
      call compute_work(plwork, d_vec_lat, ngrain, m_el)
      call compute_avg_crss(crss,crss_avg,m_el)

      converged = .false.

!
!     This loop used to be from 1 to n_edge, in order to try all vertices
!     as initial guesses.  However, this can be very time-consuming and
!     doesn't usually help.
!
      ntrials = 4

      do i_edge = 1, ntrials
         
         call find_vertex(vertex, direction, plwork,  ngrain, m_el)
                  
         call vertex_stress(sig_t, vertex, direction, ngrain, m_el)
 
         call scale_stress(sig_t, crss_avg, ngrain, m_el)

         where(.not. converged)
            sig(0, :, :) = sig_t(0, :, :)
            sig(1, :, :) = sig_t(1, :, :)
            sig(2, :, :) = sig_t(2, :, :)
            sig(3, :, :) = sig_t(3, :, :)
            sig(4, :, :) = sig_t(4, :, :)
         end where

         call solve_newton_vp(sig, d_vec_lat, crss_avg,&
  &           ier, cv_options % sx_tol, converged,&
  &           ngrain, m_el, vp_log)

         if (ier .eq. 0) GO TO 50

      enddo
      
      if ((myid .eq. 0) .and. vp_log) then
            WRITE(DFLT_U, '(A)') 'Warning:       . STRESS_SOLVE_VP: '
            WRITE(DFLT_U, '(A)') 'Warning:       . ',count(.not. converged),'&
                & grains did not converge'
            WRITE(DFLT_U, '(A)') 'Warning:       . after ', ntrials,' trials.'
      endif

      where(.not. converged) 
         sig(0, :, :) = sig_t(0, :, :)
         sig(1, :, :) = sig_t(1, :, :)
         sig(2, :, :) = sig_t(2, :, :)
         sig(3, :, :) = sig_t(3, :, :)
         sig(4, :, :) = sig_t(4, :, :)
      endwhere


50    CONTINUE

      call scale_up_sigm(sig, epseff, ngrain, m_el)

      END SUBROUTINE stress_solve_vp

!**********************************************************************
!
      SUBROUTINE scale_down_defr(d_vec, epseff, n, m)
!
!     Rescale the deformation rate to unit size.
!
!----------------------------------------------------------------------
!
!     Arguments  :
!
!     d_vec: array of deformation rates as 5-vectors (input/output)
!     epseff: effective deformation of these tensors (input)
!     n: number of grains
!     m: number of elements
!  
      INTEGER, INTENT(IN)     :: n, m
      REAL(RK), INTENT(INOUT) :: d_vec(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    :: epseff(0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER :: i, j, k
!
!----------------------------------------------------------------------
!
      !-tm from donald's verion:
      do i=0, TVEC1
         where (epseff > 0.0_RK)
            d_vec(i, :, :) = d_vec(i, :, :)/ epseff
         end where
      end do
      
      END SUBROUTINE scale_down_defr
!
!**********************************************************************
!
      SUBROUTINE compute_work(plwork, d_vec, n, m)
!
!     Compute virtual plastic work for array of deformation rates
!     applied to vertex stresses.
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER, INTENT(IN)   :: n, m
      REAL(RK), INTENT(OUT) :: plwork(0:MAX_VERT1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: d_vec(0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      REAL(RK), pointer :: sig_fs(:,:) => NULL()
      integer  :: my_phase(0:(n-1),0:(m-1))
      INTEGER  :: n_edge, i, j, iphase
!
!----------------------------------------------------------------------
!
! tsh, 1/26/03
      my_phase(n-1,:) = phase(el_sub1:el_sup1)
      
      plwork = 0.0_RK
      do iphase=1,numphases
           
         call CrystalTypeGet(ctype(iphase), VERTICES=sig_fs)
                         
         n_edge=ctype(iphase)%numvertices
         sig_fs = sig_fs * crystal_parm(9,iphase)

         do i = 0, (n_edge - 1)
            do j = 0, TVEC1
               where (my_phase .eq. iphase)
                    plwork(i, :, :) = plwork(i, :, :) + sig_fs(j+1,i+1) * d_vec(j, :, :)
               endwhere
            end do
         end do
         !
         deallocate(sig_fs)
         !
      end do !numphases
      
      END SUBROUTINE compute_work
!
!**********************************************************************
!
      SUBROUTINE find_vertex(vertex, direction, plwork, n, m)
!
!     Select the vertex which maximizes the plastic work.
!
!----------------------------------------------------------------------
!
!     Arguments:
!
!     vertex: list of optimal vertex numbers for each grain (output)
!     direction: sign of the vertex (1 or -1) (output)
!     plwork: array of virtual plastic work for each grain and all vertices (in)
!     n-edge: number of vertices
!     n: number of grains
!     m: number of elements
!
      INTEGER, INTENT(IN)     :: n, m
      INTEGER, INTENT(OUT)    :: vertex(0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT)   :: direction(0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(INOUT) :: plwork(0:MAX_VERT1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER, POINTER :: indices(:) => NULL()
      REAL(RK) :: pa(0:(n - 1), 0:(m - 1))
      INTEGER  :: i, j, iphase, numind, n_edge

      integer  :: my_phase(0:(n-1),0:(m-1))
!
!----------------------------------------------------------------------
!

      my_phase(n-1,:) = phase(el_sub1:el_sup1)
!
      vertex = 0
      direction = -1.0_RK
      do iphase=1,numphases 
         call find_indices(numind, iphase, my_phase(n-1,:), indices)
         n_edge=ctype(iphase)%numvertices
         do i = 0, (n_edge - 1)
            pa(:,indices) = abs(plwork(i, :, indices))
            where (pa .gt. direction)
               direction = pa
               vertex = i
            end where
         end do
         !
         deallocate(indices)
         !
      enddo !numphases

!  RC 6/24/2016: Reordered loops for better memory striding
      do j = 0,(m - 1)
        do i = 0,(n - 1)
          direction(i, j) = plwork(vertex(i, j), i, j) / direction(i, j)
        enddo
      enddo
!
!     To keep from the vertex being reused
!
!  RC 6/24/2016: Reordered loops for better memory striding
      do j = 0,(m - 1)
        do i = 0,(n - 1)
          plwork(vertex(i, j), i, j) = 0.0_RK
        enddo
      enddo

      END SUBROUTINE find_vertex
!
!**********************************************************************
!
      SUBROUTINE vertex_stress(sig, vertex, direction, n, m)
!
!     Set initial guess (vertex stress) for nonlinear solver.
!
!----------------------------------------------------------------------
!
!     Arguments:
!
!     sig: initial guesses for stresses (output)
!     sig_fs: vertex stresses (input)
!     vertex: optimal vertex numbers
!     direction: sign to multiply vertex stress by
!     n_edge: number of vertices
!     n: number of grains
!     m: number of elements
!
      INTEGER, INTENT(IN)   :: n, m
!      
      REAL(RK), INTENT(OUT) :: sig(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: direction(0:(n - 1), 0:(m - 1))
      INTEGER,INTENT(IN)    :: vertex(0:(n - 1), 0:(m - 1))
!

      integer :: my_phase(0:(n-1),0:(m-1))
!
!     Locals:
!
!
      REAL(RK), pointer ::  sig_fs(:,:) => NULL()
!
      INTEGER :: n_edge, i, j, iphase
!
!----------------------------------------------------------------------
!
      my_phase(n-1,:) = phase(el_sub1:el_sup1)
            
      do iphase=1,numphases
 
         call CrystalTypeGet(ctype(iphase), VERTICES=sig_fs)
         
         n_edge=ctype(iphase)%numvertices
! RC 6/24/2016 reordered the where loops as well so it only runs the parts needed

         sig_fs = sig_fs * crystal_parm(9,iphase)
!  RC 6/24/2016: Reordered loops for better memory striding
        do j = 0, (n_edge - 1)
            do i = 0, TVEC1
               where (j .eq. vertex)
                where (my_phase .eq. iphase)
                    sig(i, :, :) = sig_fs(i+1,j+1) * direction(:,:)
                end where
               end where
            enddo
        enddo
         !
         deallocate(sig_fs)
         !
      enddo !numphases
 
      END SUBROUTINE vertex_stress
!
!**********************************************************************
!
      SUBROUTINE scale_stress(sig, crss, n, m)

!
!     Rescale stress initial guess according to hardness.
!
!----------------------------------------------------------------------
!
!     Arguments:
!
!     sig: initial guess for stress (input/output)
!     crss: crystal hardnesses
!     p_hat_vec (ctype%schmid_sym): array of Schmid tensors
!     n_slip: number of slip systems
!     n: number of grains
!     m: number of elements
!
      INTEGER, INTENT(IN)     :: n, m
!
      REAL(RK), INTENT(INOUT) :: sig(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    :: crss(0:(n - 1), 0:(m - 1))
!
!     Locals:
!

      INTEGER :: my_phase(0:(n - 1), 0:(m - 1))
      INTEGER :: islip, i, j, k, iphase
      INTEGER :: numind, n_slip
      INTEGER, POINTER :: indices(:)  => NULL()
!
      REAL(RK), POINTER  :: p_hat_vec(:,:) => NULL()
      REAL(RK) :: taumax(0:(n - 1), 0:(m - 1))
      REAL(RK) :: tau   (0:(n - 1), 0:(m - 1))
!
!----------------------------------------------------------------------
!

      my_phase(n-1,:) = phase(el_sub1:el_sup1)  

      !-tm  from donald's verion:      
      taumax = 1.0d-8  ! deb/ prevent div by 0
      !taumax = 0.0d0
      
      
      do iphase=1,numphases
 
         call CrystalTypeGet(ctype(iphase), DEV=p_hat_vec)
         
         n_slip=ctype(iphase)%numslip
    
         call find_indices(numind, iphase, my_phase(n-1,:), indices)
         
         do islip = 0, (n_slip - 1)   
            call ss_project(tau, p_hat_vec(:, islip+1), sig, n, m, numind, indices)
            tau(:,indices) = dabs(tau(:,indices))
            where ((my_phase .EQ. iphase) .and.(tau .gt. taumax))
                  taumax = tau
            end where
            ! maybe:
            ! where ((tau(:,indices) .gt. taumax) taumax=tau            
         enddo !n_slip
         !
         deallocate(p_hat_vec)
         deallocate(indices)
         !
      enddo !numphases
      
      tau = crss / taumax

      do k = 0,(m - 1)
        do j = 0,(n - 1)
          do i = 0,TVEC1
            sig(i, j, k) = sig(i, j, k) * tau(j, k)
          enddo
        enddo
      enddo
      
      END SUBROUTINE scale_stress
!
!**********************************************************************
!
      SUBROUTINE ss_project(proj, plocal, tensor, n, m, numind, indices)
!
!     Compute inner product of array of tensors with a fixed tensor.
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER, INTENT(IN)   :: n, m, numind, indices(1:numind)
!
      REAL(RK), INTENT(OUT) :: proj(0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: plocal(0:TVEC1)
      REAL(RK), INTENT(IN)  :: tensor(0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER  :: i
      REAL(RK) :: proj_tmp(0:(n-1),0:(numind-1))
!
!----------------------------------------------------------------------
!
      
      proj_tmp = 0.0_RK
      do i = 0, TVEC1
         proj_tmp = proj_tmp + plocal(i) * tensor(i, :, indices)
      end do
      proj(:,indices)=proj_tmp

      END SUBROUTINE ss_project
!
!**********************************************************************
!
      SUBROUTINE solve_newton_vp(&
  &       sig, d_vec, crss, irc, eps,&
  &       converged, n, m, vp_log)
!
!     Nonlinear solution of viscoplastic crystal stress equations.
!
!----------------------------------------------------------------------
!
!     Arguments:
!
!     sig: initial guess for stresses and final solution (input/output)
!     d_vec: deformation rate for which to solve
!     crss: crystal hardnesses
!     irc: return flag
!     eps: error tolerance for nonlinear solution
!     converged: array telling what crystals have already converged
!     n: number of grains
!     m: number of elements
!     vp_log: write viscoplastic convergence output to log files
!
      INTEGER, INTENT(IN)     :: n, m
      REAL(RK), INTENT(IN)    :: eps
      INTEGER, INTENT(OUT)    :: irc
      LOGICAL, INTENT(INOUT)  :: converged(0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(INOUT) :: sig  (0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    :: d_vec(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    :: crss(0:(n - 1), 0:(m - 1))
      LOGICAL, INTENT(IN)     :: vp_log
!
!     Locals:
!
!
      LOGICAL&
  &   newton_ok(0:(n - 1), 0:(m - 1))
!
      INTEGER&
  &   i, iter, nm, inewton, k, w
!
      REAL(RK)&
  &   res      (0:(n - 1), 0:(m - 1)),&
  &   res_n    (0:(n - 1), 0:(m - 1)),&
  &   fact     (0:(n - 1), 0:(m - 1)),&
  &   ratio_res(0:(n - 1), 0:(m - 1)),&
  &   res_aux  (0:(n - 1), 0:(m - 1))
      REAL(RK)&
  &   sig_0(0:TVEC1, 0:(n - 1), 0:(m - 1)),&
  &   rhs  (0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK)&
  &   rss  (0:MAXSLIP1, 0:(n - 1), 0:(m - 1)),&
  &   shear(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK)&
  &   stif(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1)),&
  &   del_s        (0:TVEC1, 0:(n - 1), 0:(m - 1)),&
  &   xlambda      (0:TVEC1, 0:(n - 1), 0:(m - 1))
! 
!----------------------------------------------------------------------
!
      nm    = n * m
      irc   = 0
      sig_0 = 0.0_RK

      newton_ok = .true.

      do iter = 1, cv_options % sx_max_iters_newton 
        
         sig_0 = sig

         call get_res(res_n, rhs, rss, shear, sig, d_vec, crss, n, m) 

         call form_crystif(stif, rss, shear, crss, n, m)

         del_s = rhs

         call solvit(stif, del_s, n, m)

         call check_diagonals(stif, newton_ok, n, m)

         xlambda = sig_0 + del_s
      
         call get_res(res, rhs, rss, shear, xlambda, d_vec, crss, n, m) 

         fact = 1.0_RK
         ratio_res = res / res_n
 
         do while(any(ratio_res .gt. 1.0_RK .and. newton_ok .and. .not. converged))
 
            where    (ratio_res .gt. 1.0_RK .and. newton_ok .and. .not. converged)   fact = fact*0.5_RK

            if( any(fact .lt. 0.001_RK) ) then
               if (vp_log .and. (myid .eq. 0)) then
                  WRITE(DFLT_U, '(A)') 'Warning:       . SOLVE_NEWTON_VP: Line &
                    &search failure for ', count(fact .lt. 0.001_RK), ' grains.'
                  !PRINT *, ' solve_newton_vp: line search failure for ',&
                  !     &           count(fact .lt. 0.001_RK), ' grains'
               endif
               where(fact .lt. 0.001_RK) newton_ok = .false.
            endif
              
            xlambda(0, :, :) = sig_0(0, :, :) + fact * del_s(0, :, :) 
            xlambda(1, :, :) = sig_0(1, :, :) + fact * del_s(1, :, :) 
            xlambda(2, :, :) = sig_0(2, :, :) + fact * del_s(2, :, :) 
            xlambda(3, :, :) = sig_0(3, :, :) + fact * del_s(3, :, :) 
            xlambda(4, :, :) = sig_0(4, :, :) + fact * del_s(4, :, :) 

            call get_res(res_aux, rhs, rss, shear, xlambda, d_vec, crss, n, m) 

            where(ratio_res .gt. 1.0_RK .and. newton_ok .and. .not. converged)
               res = res_aux
               ratio_res = res / res_n
            endwhere

         enddo !do while
 
         where(newton_ok .and. .not. converged)
            sig(0, :, :) = sig_0(0, :, :) + fact * del_s(0, :, :)
            sig(1, :, :) = sig_0(1, :, :) + fact * del_s(1, :, :)
            sig(2, :, :) = sig_0(2, :, :) + fact * del_s(2, :, :)
            sig(3, :, :) = sig_0(3, :, :) + fact * del_s(3, :, :)
            sig(4, :, :) = sig_0(4, :, :) + fact * del_s(4, :, :)
         endwhere

         where(res .le. eps .and. newton_ok) converged = .true.
!
!        Return if all grains have converged or solution is only an estimate.
!
         inewton = count(.not. newton_ok)
         if( ((count(converged) + inewton) .eq. nm) .and. (myid .eq. 0) ) then
           if ((inewton .gt. 0) .and. vp_log) then
             WRITE(DFLT_U,'(A)') 'Info   :     > SOLVE_NEWTON_VP: Converged = ', &
                & count(converged), ' remaining = ', inewton, ' minval res = ',&
                & minval(res, mask=converged), ' maxval res = ', &
                & maxval(res, mask=converged)
            end if
             irc = inewton
           RETURN
         endif
     
      enddo !newton_iterations

      irc = -2

!1000  FORMAT(' solve_newton_vp: converged =', i6, ' remaining =',i6,' minval res.=',e15.5,' maxval res.=',e15.5)

      END SUBROUTINE solve_newton_vp
!
!**********************************************************************
!
      SUBROUTINE get_res(res, rhs, rss, shear, sig, d, crss, n, m)
!
!     Compute residual for nonlinear VP crystal stress equation.
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER   n_slip, n, m
      REAL(RK), INTENT(OUT) :: res(0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT) :: rhs(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT) :: rss(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(OUT) :: shear(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: sig(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: d(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)  :: crss(0:(n - 1), 0:(m - 1))
!
      integer :: my_phase(0:(n-1),0:(m-1))
!
!     Locals:
!
      INTEGER, POINTER  :: indices(:) => NULL()
      REAL(RK), POINTER :: p(:,:) => NULL()
      INTEGER :: islip, j, k, iphase, numind, i
      REAL(RK) :: xm_fake 

!
!----------------------------------------------------------------------
!
! tsh, 1/26/03
      my_phase(n-1,:) = phase(el_sub1:el_sup1)      
      
      rhs = - d
      res = 0.0_RK

      do iphase=1,numphases
   
         call CrystalTypeGet(ctype(iphase), DEV=p)
   
         n_slip=ctype(iphase)%numslip

         call find_indices(numind, iphase, my_phase(n-1,:), indices)

         do islip = 0, (n_slip - 1)

            call ss_project(rss(islip,:,:),p(:,islip+1), sig, n, m, numind, indices)
            rss(islip,:,indices)=rss(islip,:,indices)/crss(:,indices) 

            where (abs(rss(islip,:,indices)) .lt. t_min(iphase)) 
               rss(islip,:,indices) = 0.0_RK
            endwhere

            xm_fake=0.4_RK
            !xm_fake=0.02d0
            
            call power_law(shear(islip, :, :), rss(islip, :, :),&
  &                xm_fake, crystal_parm(1,iphase), t_min(iphase),&
  &                n, m, numind, indices)
 
            do j=0,TVEC1
               rhs(j, :, indices) = rhs(j, :, indices) + p(j+1, islip+1) * shear(islip, :, indices)
            enddo
            
         enddo !n_slip
         !
         deallocate(p)
         deallocate(indices)
         !
      enddo !numphases
         
      do j = 0, DIMS1
         res = res + rhs(j, :, :)**2.0_RK
      end do
      res = sqrt(res)
      
      END SUBROUTINE get_res
!
!**********************************************************************
!
      SUBROUTINE form_crystif(stif, rss, shear, crss, n, m)
!
!     Form single crystal stiffness matrix.
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER   n_slip, n, m
      REAL(RK), pointer :: p(:,:) => NULL()
      REAL(RK), INTENT(OUT)  ::   stif(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)   ::   rss(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)   ::   shear(0:MAXSLIP1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)   ::   crss(0:(n - 1), 0:(m - 1))
!
!     Locals:

      integer   my_phase(0:(n-1),0:(m-1))
      INTEGER   i, j, k, islip, iphase, numind
      INTEGER, POINTER :: indices(:) => NULL()
      REAL(RK)    comp(0:(n - 1), 0:(m - 1))
      
      REAL(RK) xm_fake
!
!----------------------------------------------------------------------
!
      my_phase(n-1,:) = phase(el_sub1:el_sup1)
            
      stif = 0.0_RK
      do iphase=1,numphases

         call CrystalTypeGet(ctype(iphase), DEV=p)
         
         n_slip=ctype(iphase)%numslip             

         call find_indices(numind, iphase, my_phase(n-1,:), indices)

         do islip = 0, (n_slip - 1)
  
            xm_fake=0.4_RK
            !xm_fake=0.02d0
             
            call compliance(comp, rss(islip, :, :),&
  &                   shear(islip, :, :), crss, xm_fake, t_min(iphase),&
  &                   n, m, numind, indices)
  
            do j = 0, TVEC1    
               do k = 0, TVEC1
                  stif(j, k, :, indices) = stif(j, k, :, indices) -&
  &                 comp(:,indices) * p(j+1, islip+1) * p(k+1, islip+1)
               enddo
            enddo
         enddo !num_slip
         !
         deallocate(p)
         deallocate(indices)
         !
      enddo !numphases

      END SUBROUTINE form_crystif
!
!**********************************************************************
!
      SUBROUTINE power_law(power, t, xm, a_0, t_min, n, m, numind, indices)
!
!     Power law for VP single crystal.
!
!----------------------------------------------------------------------
!
!     Arguments:
!        
      INTEGER   n, m, numind, indices(1:numind)
      REAL(RK)  :: t_min, xm, a_0 
      REAL(RK), INTENT(OUT)  ::  power(0:(n - 1), 0:(m - 1)) 
      REAL(RK), INTENT(IN)   ::  t(0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      REAL(RK)    p
      REAL(RK)    power_tmp(0:(n-1),0:(numind-1))
      REAL(RK)    at(0:(n - 1), 0:(numind - 1))
      REAL(RK)    alog(0:(n - 1), 0:(numind - 1))
      REAL(RK)    blog(0:(n - 1), 0:(numind - 1))
      INTEGER   i
!
!
!----------------------------------------------------------------------

      p  = 1.0_RK / xm - 1.0_RK
      at = abs(t(:,indices))
      
      where (at .gt. t_min)
          alog = dlog(at)
          blog = p * alog
          power_tmp = a_0 * t(:,indices) * dexp(blog)
      elsewhere
          power_tmp = 0.0_RK
      endwhere  
      
      power(:,indices)=power_tmp 
     
      END SUBROUTINE power_law
!
!**********************************************************************
!
      SUBROUTINE compliance(comp, t, shear, crss, xm, t_min,&
  &     n, m, numind, indices)
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER  :: n, m, numind, indices(1:numind)
      REAL(RK) :: xm, t_min
      REAL(RK), INTENT(OUT)  ::   comp(0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)   ::   t(0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)   ::   shear(0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)   ::   crss(0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      REAL(RK) :: comp_tmp(0:(n-1),0:(numind-1))
!
!----------------------------------------------------------------------
!
      comp_tmp = 0.0_RK

      WHERE (ABS(t(:,indices)) .GT. t_min)&
           &           comp_tmp = shear(:, indices)/(t(:, indices)*crss(:, indices)*xm)
      comp(:,indices)=comp_tmp

      END SUBROUTINE compliance
!
!**********************************************************************
!
      SUBROUTINE solvit(a, x, n, m)
!
!     Solve an array of symmetric positive definite 5x5 systems.
!
!-----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER   n, m
      REAL(RK), INTENT(IN)    ::  a(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(INOUT) ::  x(0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      REAL(RK), DIMENSION(0:(n - 1), 0:(m - 1)) :: a11, a21, a22,&
  &                                       a31, a32, a33,&
  &                                       a41, a42, a43, a44,&
  &                                       a51, a52, a53, a54, a55,&
  &                                       x1, x2, x3, x4, x5,&
  &                                       v1, v2, v3, v4
!
!-----------------------------------------------------------------------
!
      a11 = a(0, 0, :, :)
      a21 = a(1, 0, :, :)
      a31 = a(2, 0, :, :)
      a41 = a(3, 0, :, :)
      a51 = a(4, 0, :, :)
      a22 = a(1, 1, :, :)
      a32 = a(2, 1, :, :)
      a42 = a(3, 1, :, :)
      a52 = a(4, 1, :, :)
      a33 = a(2, 2, :, :)
      a43 = a(3, 2, :, :)
      a53 = a(4, 2, :, :)
      a44 = a(3, 3, :, :)
      a54 = a(4, 3, :, :)
      a55 = a(4, 4, :, :)
 
      x1 = x(0, :, :)
      x2 = x(1, :, :)
      x3 = x(2, :, :)
      x4 = x(3, :, :)
      x5 = x(4, :, :)
 
! **  A = LDL'.

! **  j = 1.

      a21 = a21 / a11
      a31 = a31 / a11
      a41 = a41 / a11
      a51 = a51 / a11

! **  j = 2.

      v1 = a21 * a11

      a22 = a22 - a21 * v1

      a32 = (a32 - a31 * v1) / a22
      a42 = (a42 - a41 * v1) / a22
      a52 = (a52 - a51 * v1) / a22

! **  j = 3.

      v1 = a31 * a11
      v2 = a32 * a22

      a33 = a33 - a31 * v1 - a32 * v2
         
      a43 = (a43 - a41 * v1 - a42 * v2) / a33
      a53 = (a53 - a51 * v1 - a52 * v2) / a33

! **  j = 4.

      v1 = a41 * a11
      v2 = a42 * a22
      v3 = a43 * a33

      a44 = a44 - a41 * v1 - a42 * v2 - a43 * v3
         
      a54 = (a54 - a51 * v1 - a52 * v2 - a53 * v3) / a44

! **  j = 5.

      v1 = a51 * a11
      v2 = a52 * a22
      v3 = a53 * a33
      v4 = a54 * a44

      a55 = a55 - a51 * v1 - a52 * v2 - a53 * v3 - a54 * v4
         
! **  Ly=b.

      x2 = x2 - a21 * x1
      x3 = x3 - a31 * x1 - a32 * x2
      x4 = x4 - a41 * x1 - a42 * x2 - a43 * x3
      x5 = x5 - a51 * x1 - a52 * x2 - a53 * x3 - a54 * x4

! **  Dz=y.

      x1 = x1 / a11
      x2 = x2 / a22
      x3 = x3 / a33
      x4 = x4 / a44
      x5 = x5 / a55

! **  L'x=z.

      x4 = x4 - a54 * x5
      x3 = x3 - a43 * x4 - a53 * x5
      x2 = x2 - a32 * x3 - a42 * x4 - a52 * x5
      x1 = x1 - a21 * x2 - a31 * x3 - a41 * x4 - a51 * x5

      x(0, :, :) = x1
      x(1, :, :) = x2
      x(2, :, :) = x3
      x(3, :, :) = x4
      x(4, :, :) = x5

      END SUBROUTINE solvit
!
!**********************************************************************
!
      SUBROUTINE check_diagonals(stif, newton_ok, n, m)
!
!     Determine where diagonal elements are small.
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      INTEGER   n, m
      LOGICAL, INTENT(INOUT) :: newton_ok(0:(n - 1), 0:(m - 1))
      REAL(RK),  INTENT(IN)    :: stif(0:TVEC1, 0:TVEC1, 0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER   i
!
!----------------------------------------------------------------------
!
      do i = 0, TVEC1
         where (abs(stif(i, i, :, :)) .lt. VTINY) newton_ok = .false.
      end do
      
      END SUBROUTINE check_diagonals
!
!**********************************************************************
!
      SUBROUTINE scale_up_sigm(sig, epseff, n, m)
!
!     Rescale the stress after solution is found.
!
!----------------------------------------------------------------------
!
!     Arguments:
!
!     sig: stress (input/output)
!     xm: rate dependence
!     epseff: effective deformation rate
!     n: number of grains
!     m: number of elements
!
      INTEGER   n, m
      REAL(RK), INTENT(INOUT) :: sig(0:TVEC1, 0:(n - 1), 0:(m - 1))
      REAL(RK), INTENT(IN)    :: epseff(0:(n - 1), 0:(m - 1))
!
!     Locals:
!
      INTEGER   i, j, k, numind, iphase
      INTEGER, pointer :: indices(:)
      REAL(RK)    scale(0:(n - 1), 0:(m - 1))
      integer   my_phase(0:(m - 1))
      
      REAL(RK) xm_fake
!
!----------------------------------------------------------------------
!
      my_phase(:) = phase(el_sub1:el_sup1)
      
      do iphase=1,numphases
         !-tm        
         ! scale = epseff**crystal_parm(0,iphase)         
         xm_fake=0.4_RK
         !xm_fake=0.02d0         
         scale = epseff**xm_fake
         
         call find_indices(numind, iphase, my_phase, indices)
         
         do j = 0,(n - 1)
            do i = 0,TVEC1
               sig(i,j,indices) = sig(i,j,indices) * scale(j,indices)
            enddo
         enddo
         !
         deallocate(indices)
         !
      enddo !numphases
      
      END SUBROUTINE scale_up_sigm

!**********************************************************************
!
    SUBROUTINE compute_avg_crss(crss,crss_avg,m_el)
!
!     Computes the average strength of the crystal slip systems
!
!----------------------------------------------------------------------
!

        REAL(RK), INTENT(IN)    :: crss  (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
        REAL(RK), INTENT(INOUT)   :: crss_avg  (0:ngrain1, el_sub1:el_sup1)
        INTEGER, INTENT(IN) :: m_el

        INTEGER,  POINTER :: indices(:) => NULL()

        INTEGER :: islip, n_slip, iphase, numind
        INTEGER :: my_phase(0:m_el-1)
        my_phase(:) = phase(el_sub1:el_sup1)
        !Goes through the total number of phases in the material
        do iphase=1,numphases

            !Finds the numbers of slip systems
            call CrystalTypeGet(ctype(iphase))
            n_slip=ctype(iphase)%numslip
            !Finds the indices corresponding to the current phase the loop is on
            call find_indices(numind, iphase, my_phase, indices)
            !Sums up the crystal slip strengths corresponding to the given indices
            do islip=0, n_slip-1
                    crss_avg(:,indices+el_sub1)=crss_avg(:,indices+el_sub1)+crss(islip,:,indices+el_sub1)
            enddo
            !Calculates the average of these slip systems
            crss_avg(:,indices+el_sub1)=crss_avg(:,indices+el_sub1)/n_slip
            !WRITE(*,*) numind
            !write(*,*) size(indices)
            !write(*,*) my_phase
            ! write(*,*) crss_avg(:,indices)
            deallocate (indices)

        enddo

    END SUBROUTINE compute_avg_crss

END MODULE StressSolveVpModule
