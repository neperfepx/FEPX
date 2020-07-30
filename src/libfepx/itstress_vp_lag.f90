! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE ItMethodVpModule

  USE parallel_mod

  USE IntrinsicTypesModule, RK=>REAL_KIND
  USE gather_scatter
  USE units_mod
  USE READ_INPUT_MOD
  USE microstructure_mod
  USE surf_info_mod
  USE DimsModule
  USE ConvergenceModule, ONLY: cv_options
  USE StiffnessVpModule

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: itmethod_vp

CONTAINS

      INTEGER FUNCTION itmethod_vp(&
  &   itype, bcs, pforce, vel, elpress, evel,&
  &   dof_trace, np_trace, qr5x5,&
  &   wts, epseff, dtime, incr)
!
!     Driver for the viscoplastic iteration required for a single
!     time increment.
!
!----------------------------------------------------------------------
!
!     Arguments:
!
      TYPE(trace) dof_trace, np_trace
!
      LOGICAL :: bcs(dof_sub1:dof_sup1)
!
      INTEGER  :: itype
      INTEGER  :: n_slip, n_edge
      INTEGER  :: incr
!
      REAL(RK) :: dtime
      REAL(RK) :: pforce(dof_sub1:dof_sup1)
      REAL(RK) :: vel(dof_sub1:dof_sup1)
      REAL(RK) :: elpress(el_sub1:el_sup1), evel(0:kdim1, el_sub1:el_sup1)
      REAL(RK) :: qr5x5(0:TVEC1, 0:TVEC1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK) :: crss(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK) :: wts(0:ngrain1, el_sub1:el_sup1)
      REAL(RK) :: epseff(el_sub1:el_sup1)
      ! tm268_M13:
!      REAL(RK) :: eqplas(el_sub1:el_sup1)
!
!     Locals:
!
      INTEGER :: iter, idiv, i, j, cg_solver_ebe, cg_iter_out
      INTEGER :: cg_max_iters
      REAL(RK) :: cg_tol
!
      REAL(RK) :: &
  &   eps_o, d_norm,&
  &   part_delomax, delomax,&
  &   part_u_norm, u_norm,&
  &   part_x_norm, x_norm,&
  &   part_p_norm, p_norm
!
      REAL(RK) :: estiff(0:kdim1, 0:kdim1, el_sub1:el_sup1)
      REAL(RK) :: ecoords(0:kdim1, el_sub1:el_sup1)
      REAL(RK) :: force(dof_sub1:dof_sup1), eforce(0:kdim1, el_sub1:el_sup1)
      REAL(RK) :: vel_o(dof_sub1:dof_sup1), del(dof_sub1:dof_sup1)
      REAL(RK) :: gdiag(dof_sub1:dof_sup1)
      REAL(RK) :: pscale(el_sub1:el_sup1), pcnst(0:kdim1, el_sub1:el_sup1)
      REAL(RK) :: eqplas_tr(el_sub1:el_sup1)
!
!----------------------------------------------------------------------
!
!      write(ounits(LOG_U), *) 'INCREMENT: ', incr
!      write(ounits(LOG_U), *)&
!  &      'NL-iter CG-iters       CG-resid       L2-diff      MAX-diff        NL-tol'
!      write(ounits(LOG_U), *)&
!  &      '------- --------       --------       -------      --------        ------'
!
      itmethod_vp = 1
!
      cg_max_iters = cv_options % cg_max_iters
      cg_tol = cv_options % cg_tol

      eps_o = 1.0e+30
      vel_o = vel
!
      call part_gather(evel, vel, nodes, dof_trace)
      call part_gather(ecoords, coords, nodes, dof_trace)
!

! ----------- NL-iteration loop -----------------------------------------
      nonlinear_iteration : do iter = 1, cv_options % nl_max_iters
! -----------------------------------------------------------------------

         if (myid .eq. 0) then
            !
            write(DFLT_U,'(A,I0)') 'Info   :     > ITMETHOD_VP: Iteration ', iter
            write(DFLT_U,'(A)',ADVANCE='NO') 'Info   :       . Solving NL iteration... '
            !
         end if

         estiff = 0.d0
         eforce = 0.d0
         force  = pforce !pforce=0
         eqplas_tr = 0.d0

         call element_stif_vp(itype, estiff, ecoords, evel, pscale,&
  &      pcnst, qr5x5, wts, eqplas_tr, epseff,&
  &      dtime, incr)


         do i = 0, kdim1
           eforce(i, :) = eforce(i, :) - elpress * pcnst(i, :)
         end do

         ! eforce-->force
         call part_scatter(force, eforce, nodes, .false., dof_trace)

         ! zero the forces where velocities are specified
         where (bcs) force = 0.0

         !** preconditioned conjugate gradient solver **
         !
         ! Form the diagonal part of the stiffness matrix
         call assemble_diagonals(gdiag, estiff, kdim,dof_sub1, dof_sup1, el_sub1, el_sup1, dof_trace, nodes)

         ! compute the velocity field (vel) using the conjugate gradient method
         cg_iter_out = cg_solver_ebe(vel, d_norm, force, gdiag,&
  &                    estiff, bcs, kdim, dof_sub1, dof_sup1, el_sub1, el_sup1,&
  &                    cg_max_iters, cg_tol, dof_trace, nodes)

         ! vel-->evel
         call part_gather(evel,vel,nodes,dof_trace)

         call recover_pressure_vp(evel, pscale, elpress, pcnst)
!
         part_u_norm = sum((vel_o - vel) * (vel_o - vel))
         call par_sum(part_u_norm, u_norm)
         u_norm = sqrt(u_norm)
!
         part_x_norm = maxval(abs(vel_o - vel))
         call par_max(part_x_norm, x_norm)

!         write(ounits(LOG_U), 100) iter, cg_iter_out, d_norm, u_norm, x_norm
! 100     FORMAT(i8,1x,i8,1x,4d14.4)

         if (myid .eq. 0) then
            !
            write(DFLT_U,'(A,E10.4,A,I0,A)', ADVANCE='YES') 'R = ', u_norm, &
                & ' (', cg_iter_out, ' iters)'
!            write(DFLT_U,'(A,E10.4,A,I0)') 'Info   :       . Final residual = ', u_norm, &
!                & ', No iterations ', cg_iter_out
            !
         end if

         vel_o = vel

         if (iter .gt. 1) GO TO 20

         if (u_norm .gt. eps_o) then
            idiv = idiv + 1
         else
            idiv = 0
         endif

         if (idiv .ge. 5) then
!            write(ounits(LOG_U), *) 'Iterations are diverging.'
            itmethod_vp = -1
            RETURN
         endif

         eps_o = u_norm

      enddo nonlinear_iteration
      ! ----------- NL-iteration loop -----------------------------------------



!      write(ounits(LOG_U), 110) iter, d_norm
! 110  FORMAT(' No convergence for initial guess', /,&
!  &       ' iter =', i4, ' d_norm =', g10.3)

      itmethod_vp = -1
      RETURN

 20   CONTINUE
!deb
!     I don't know why one more iteration is required here,
!     but it seems like it can't hurt.  Note that the stiffness
!     is not recomputed.
!deb
      force = pforce
      ! tm268_M13:
!      eqplas = eqplas + eqplas_tr

      do i = 0, kdim1
        eforce(i, :) = -elpress * pcnst(i, :)
      end do

      ! eforce-->force
      call part_scatter(force, eforce, nodes, .false., dof_trace)

      ! zero the forces where velocities are specified
      where (bcs) force = 0.0

      ! compute the velocity field (vel) using the conjugate gradient method
      cg_iter_out = cg_solver_ebe(vel, d_norm, force, gdiag,&
  &                 estiff, bcs, kdim, dof_sub1, dof_sup1, el_sub1, el_sup1,&
  &                 cg_max_iters, cg_tol, dof_trace, nodes)

      ! vel-->evel
      call part_gather(evel, vel, nodes, dof_trace)

      call recover_pressure_vp(evel, pscale, elpress, pcnst)

      part_u_norm = sum((vel_o - vel) * (vel_o - vel))
      call par_sum(part_u_norm, u_norm)
      u_norm = sqrt(u_norm)

      del = abs(vel_o)
      part_delomax = maxval(del)
      call par_max(part_delomax, delomax)
      if (delomax .eq. 0.0) delomax = 1.0d0

      del = abs((vel_o - vel) / delomax)

      part_p_norm = maxval(del)
      call par_max(part_p_norm, p_norm)

      part_x_norm = maxval(abs(vel_o - vel))
      call par_max(part_x_norm, x_norm)

      iter = iter + 1
!      write(ounits(LOG_U), 100) iter, cg_iter_out, d_norm, u_norm, x_norm
!      write(ounits(LOG_U), '(a58,e16.6)')&
!  &   'relative maximum norm of difference for final iteration: ', p_norm

      vel_o = vel

      ! tm08: in the next line i think it should be iter-1
      ! original line:
      ! if (myid .eq. 0) write(DFLT_U,*) '   completed in ', iter, ' iterations.'
      ! new line:
      if (myid .eq. 0) write(DFLT_U,'(A,I0,A)') 'Info   :     > Converged in ',&
          & iter-1, ' iterations'

      RETURN
    END FUNCTION itmethod_vp
!
!**********************************************************************
!
      SUBROUTINE recover_pressure_vp(evel, pscale, elpress, pcnst)
!
!     Compute pressure from velocity field.
!
!----------------------------------------------------------------------
!
!     Arguments.
!
      REAL(RK) :: &
  &   evel (0:kdim1, el_sub1:el_sup1),&
  &   pcnst(0:kdim1, el_sub1:el_sup1)
      REAL(RK) :: &
  &   pscale (el_sub1:el_sup1),&
  &   elpress(el_sub1:el_sup1)
!
!     Locals:
!
      INTEGER i
      REAL(RK) ::   sum(el_sub1:el_sup1)
!
!----------------------------------------------------------------------
!
      sum = 0.0d0
      do i = 0, kdim1
         sum = sum + pcnst(i, :) * evel(i, :)
      end do

      elpress = elpress + cv_options % pacc * pscale * sum

      RETURN
    END SUBROUTINE recover_pressure_vp
!
  END MODULE ItMethodVpModule
