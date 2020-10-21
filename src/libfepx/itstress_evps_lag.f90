! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE ItMethodEvpsModule

  USE parallel_mod
  USE parallel_matrix_mod
  USE gather_scatter
  
  USE DimsModule
  USE units_mod

  USE READ_INPUT_MOD
  USE microstructure_mod
  USE StiffnessEvpsModule
  USE MATRIX_OPERATIONS_MOD
  USE ConvergenceModule, ONLY: cv_options
  USE WRITE_OUTPUT_MOD

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: itmethod_evps

CONTAINS

      INTEGER FUNCTION itmethod_evps(&
  &       bcs, pforce, vel, elpress, evel,&
  &       dtrace, ntrace,&
  &       c0_angs, c_angs, sig_vec_n, sig_vec, crss_n,&
  &       crss, rstar_n, rstar, keinv,&
  &       e_elas_kk_bar, e_elas_kk, sig_kk, jiter_state,&
  &       wts, epseff, &
  &       dtime, incr, e_bar_vec,&
  &       converged_solution, auto_time, iter)
!
!----------------------------------------------------------------------
!     itmethod_evps is set:
!     >0 (1)  if the solution converged
!     <0 (-1) if the solution did not converged
!             (either the maximum number of iterations was exceeded
!              or the solution was badly diverging)
!
!     Modules.
!
!
!     Arguments:
!
!     dtrace: gather/scatter trace for degrees of freedom
!     ntrace: gather/scatter trace for nodal point arrays
!
      TYPE(trace) :: dtrace, ntrace

      LOGICAL, INTENT(INOUT)  :: converged_solution
      LOGICAL, INTENT(IN)     :: bcs(dof_sub1:dof_sup1)
      REAL(RK), INTENT(INOUT) :: vel(dof_sub1:dof_sup1)

      INTEGER, INTENT(IN) :: incr   ! current load increment

      REAL(RK)  :: dtime
      REAL(RK), INTENT(IN)  :: keinv(0:TVEC1, 1:numphases)

      REAL(RK), INTENT(IN)    :: pforce(dof_sub1:dof_sup1)
      REAL(RK)                :: evel(0:kdim1, el_sub1:el_sup1)

      REAL(RK), INTENT(IN)    :: c0_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: c_angs (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) :: rstar  (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(IN)    :: rstar_n(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(INOUT) :: sig_vec_n(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)  ! enters=0 at the first increment for all the grains
                                                                                 ! exits =/0 only for grain=0
      REAL(RK), INTENT(OUT)   :: sig_vec  (0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)

      REAL(RK), INTENT(OUT)   :: e_bar_vec(0:TVEC1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)

      REAL(RK), INTENT(IN)    :: crss_n(0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1)
      REAL(RK), INTENT(OUT)   :: crss  (0:MAXSLIP1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)

      REAL(RK), INTENT(IN)    :: wts   (0:ngrain1, el_sub1:el_sup1)

      REAL(RK), INTENT(OUT)   :: epseff (el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT)   :: elpress(el_sub1:el_sup1, 0:nqpt1)

      REAL(RK), INTENT(OUT)   :: e_elas_kk_bar(el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT)   :: e_elas_kk    (el_sub1:el_sup1, 0:nqpt1)
      REAL(RK), INTENT(OUT)   :: sig_kk (el_sub1:el_sup1, 0:nqpt1)

      INTEGER, INTENT(OUT)    :: jiter_state(0:ngrain1, el_sub1:el_sup1)
      INTEGER, INTENT(OUT)    :: iter

!
!     Locals:
!
      CHARACTER(LEN=128) :: message

      INTEGER, PARAMETER :: NR_SLOPE_START = 2

      INTEGER  :: itmethod, idiv, i, j, k, cg_solver_ebe, cg_iter_out, ier
      INTEGER :: cg_max_iters, auto_time
      REAL(RK) :: cg_tol
      REAL(RK) :: nl_tol_strict, nl_tol_loose, nl_tol_min
      !
      REAL(RK) :: part_r_norm, r_norm, r_norm_o, r_norm_n
      REAL(RK) :: part_rx_norm, rx_norm
      REAL(RK) :: part_f_norm, f_norm
      REAL(RK) :: part_delu_norm, delu_norm, delu_norm_o
      REAL(RK) :: part_delux_norm, delux_norm
      REAL(RK) :: part_u_norm, u_norm
      REAL(RK) :: d_norm
      !
      REAL(RK) :: vel_o(dof_sub1:dof_sup1) ! velocity for previous iteration
      REAL(RK) :: delta_vel(dof_sub1:dof_sup1)
      REAL(RK) :: vel_SA(dof_sub1:dof_sup1)
      !
      REAL(RK) :: estiff(0:kdim1, 0:kdim1, el_sub1:el_sup1)
      REAL(RK) :: etanstiff(0:kdim1, 0:kdim1, el_sub1:el_sup1)
      REAL(RK) :: f_vec(0:kdim1, el_sub1:el_sup1)
      REAL(RK) :: tcoords(dof_sub1:dof_sup1), ecoords(0:kdim1, el_sub1:el_sup1)
      REAL(RK) :: force(dof_sub1:dof_sup1)
      REAL(RK) :: eforce(0:kdim1, el_sub1:el_sup1)
      REAL(RK) :: resid(dof_sub1:dof_sup1)
      REAL(RK) :: eresid(0:kdim1, el_sub1:el_sup1)
      REAL(RK) :: tc_angs(0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK) :: trstar  (0:DIMS1, 0:DIMS1, 0:ngrain1, el_sub1:el_sup1, 0:nqpt1)
      REAL(RK) :: gdiag(dof_sub1:dof_sup1)
      !
      INTEGER  :: m_el
      !
      REAL(RK) :: velcur(3, 10)
      integer :: io, ielem
      integer :: j1, j2, j3
      !
      ! local node number for mid-nodes correction
      ! e1, e2: end node
      ! m : mid-node
      integer :: e1, e2, m
      !
      REAL(RK) :: e_ones(0:kdim1, el_sub1:el_sup1)
      REAL(RK) :: g_ones(dof_sub1:dof_sup1)

      LOGICAL :: NR, NR_attempt, NR_slow, iter_converged
      INTEGER :: NR_iter
      REAL(RK) :: NR_tol_switch_ref, NR_tol_conv
      REAL(RK) :: NR_tol_switch, NR_conv, NR_conv_o
!
!----------------------------------------------------------------------
!
      NR_tol_switch_ref = cv_options % nr_tol_switch_ref
      NR_tol_conv = cv_options % nr_tol_conv

      cg_max_iters = cv_options % cg_max_iters
      cg_tol = cv_options % cg_tol

      nl_tol_strict = cv_options % nl_tol_strict
      nl_tol_loose = cv_options % nl_tol_loose
      nl_tol_min = cv_options % nl_tol_min

      m_el = el_sup1 - el_sub1 + 1

!      write(ounits(LOG_U), *) 'INCREMENT: ', incr
!      write(ounits(LOG_U), *) 'NL-iter CG-iters       CG-resid     delU-norm        U-norm'
!      write(ounits(LOG_U), *) '------- --------       --------     ---------        ------'


      ! tsh, initialize ones arrays
      e_ones = 1.0_RK
      g_ones = 0.0_RK

      ! scatter ones arrays and store multiplicity
      call part_scatter(g_ones, e_ones, nodes, .false., dtrace)

      ! initialization
      r_norm_o = 0.0_RK
      delu_norm_o = 0.0_RK
      vel_o = vel
      vel_SA = vel

      ! set itmethod_evps=1 at the beginning
      itmethod_evps = 1
      cg_iter_out = 0
      d_norm = 0.0_RK

      iter_converged = .false.
      NR = .false.
      NR_attempt = .false.
      NR_slow = .false.
      NR_iter = 0
      NR_tol_switch = NR_tol_switch_ref
      NR_conv = 1.0_RK
      NR_conv_o = 1.0_RK
      idiv = 0

! ----------- NL-iteration loop -----------------------------------------
      nonlinear_iteration : DO iter = 1, cv_options % nl_max_iters
! -----------------------------------------------------------------------

!tsh
! adjust mid-nodes
! re-position the mid-node to the middle of the edge-nodes

         if (myid .eq. 0) then
            !
            write(DFLT_U,'(A,I0)') 'Info   :     > ITMETHOD_EVPS: Iteration ', iter
            !
         end if

         call part_gather(ecoords, coords, nodes, dtrace)

         m  = 1
         e1 = 0
         e2 = 2
         ecoords(3*m,   el_sub1:el_sup1) =&
  &     (ecoords(3*e1,   el_sub1:el_sup1) +&
  &        ecoords(3*e2,   el_sub1:el_sup1)) / 2.0_RK
         ecoords(3*m+1, el_sub1:el_sup1) =&
  &     (ecoords(3*e1+1, el_sub1:el_sup1) +&
  &        ecoords(3*e2+1, el_sub1:el_sup1)) / 2.0_RK
         ecoords(3*m+2, el_sub1:el_sup1) =&
  &     (ecoords(3*e1+2, el_sub1:el_sup1) +&
  &        ecoords(3*e2+2, el_sub1:el_sup1)) / 2.0_RK
!
         m  = 3
         e1 = 2
         e2 = 4
         ecoords(3*m,   el_sub1:el_sup1) =&
  &     (ecoords(3*e1,   el_sub1:el_sup1) +&
  &        ecoords(3*e2,   el_sub1:el_sup1)) / 2.0_RK
         ecoords(3*m+1, el_sub1:el_sup1) =&
  &     (ecoords(3*e1+1, el_sub1:el_sup1) +&
  &        ecoords(3*e2+1, el_sub1:el_sup1)) / 2.0_RK
         ecoords(3*m+2, el_sub1:el_sup1) =&
  &     (ecoords(3*e1+2, el_sub1:el_sup1) +&
  &        ecoords(3*e2+2, el_sub1:el_sup1)) / 2.0_RK
!
         m  = 5
         e1 = 0
         e2 = 4
         ecoords(3*m,   el_sub1:el_sup1) =&
  &     (ecoords(3*e1,   el_sub1:el_sup1) +&
  &        ecoords(3*e2,   el_sub1:el_sup1)) / 2.0_RK
         ecoords(3*m+1, el_sub1:el_sup1) =&
  &     (ecoords(3*e1+1, el_sub1:el_sup1) +&
  &        ecoords(3*e2+1, el_sub1:el_sup1)) / 2.0_RK
         ecoords(3*m+2, el_sub1:el_sup1) =&
  &     (ecoords(3*e1+2, el_sub1:el_sup1) +&
  &        ecoords(3*e2+2, el_sub1:el_sup1)) / 2.0_RK
!
         m  = 6
         e1 = 0
         e2 = 9
         ecoords(3*m,   el_sub1:el_sup1) =&
  &     (ecoords(3*e1,   el_sub1:el_sup1) +&
  &        ecoords(3*e2,   el_sub1:el_sup1)) / 2.0_RK
         ecoords(3*m+1, el_sub1:el_sup1) =&
  &     (ecoords(3*e1+1, el_sub1:el_sup1) +&
  &        ecoords(3*e2+1, el_sub1:el_sup1)) / 2.0_RK
         ecoords(3*m+2, el_sub1:el_sup1) =&
  &     (ecoords(3*e1+2, el_sub1:el_sup1) +&
  &        ecoords(3*e2+2, el_sub1:el_sup1)) / 2.0_RK
!
         m  = 7
         e1 = 2
         e2 = 9
         ecoords(3*m,   el_sub1:el_sup1) =&
  &     (ecoords(3*e1,   el_sub1:el_sup1) +&
  &        ecoords(3*e2,   el_sub1:el_sup1)) / 2.0_RK
         ecoords(3*m+1, el_sub1:el_sup1) =&
  &     (ecoords(3*e1+1, el_sub1:el_sup1) +&
  &        ecoords(3*e2+1, el_sub1:el_sup1)) / 2.0_RK
         ecoords(3*m+2, el_sub1:el_sup1) =&
  &     (ecoords(3*e1+2, el_sub1:el_sup1) +&
  &        ecoords(3*e2+2, el_sub1:el_sup1)) / 2.0_RK
!
         m  = 8
         e1 = 4
         e2 = 9
         ecoords(3*m,   el_sub1:el_sup1) =&
  &     (ecoords(3*e1,   el_sub1:el_sup1) +&
  &        ecoords(3*e2,   el_sub1:el_sup1)) / 2.0_RK
         ecoords(3*m+1, el_sub1:el_sup1) =&
  &     (ecoords(3*e1+1, el_sub1:el_sup1) +&
  &        ecoords(3*e2+1, el_sub1:el_sup1)) / 2.0_RK
         ecoords(3*m+2, el_sub1:el_sup1) =&
  &     (ecoords(3*e1+2, el_sub1:el_sup1) +&
  &        ecoords(3*e2+2, el_sub1:el_sup1)) / 2.0_RK

!        reset coordinates
         coords = 0.0_RK
!        scatter e_coordinates
         call part_scatter(coords, ecoords, nodes, .false., dtrace)

!        divide the coordinates by the multiplicity
         coords = coords / g_ones

         ! estimate the new coordinates based on the velocity vel
         ! (vel is the guess for the velocity field)
         tcoords = coords + dtime * vel

         call part_gather(ecoords, tcoords, nodes, dtrace)
         call part_gather(evel, vel, nodes, dtrace)

         estiff = 0.0_RK
         etanstiff = 0.0_RK
         eforce = 0.0_RK
         force  = pforce !pforce=0
         tc_angs = spread(c_angs, 5, nqpt)
         trstar = spread(rstar, 5, nqpt)

         ! enter the subroutine that, through iterations at the constitutive level, computes:
         ! - the state of each element:
         !     e*_kk  : e_elas_kk
         !     tau_kk : sig_kk
         !     tau'   : sig_vec
         !     R*     : rstar
         !     g      : crss
         ! - matrices [S]: estiff
         ! - vectors [S]{h}: f_vec

         ! e_bar_vec: OUT
         ! sig_vec_n: INOUT (enters=0 at the first increment, only grain=0 is updated))

         call element_stif_evps(&
  &          estiff, etanstiff, f_vec, ecoords, evel,&
  &          c0_angs, tc_angs, sig_vec_n, sig_vec, crss_n,&
  &          crss, rstar_n, trstar,&
  &          e_bar_vec, e_elas_kk_bar, e_elas_kk, sig_kk, jiter_state,&
  &          keinv, wts, epseff,&
  &          dtime, incr,&
  &          converged_solution, auto_time, NR)

         ! essentially renames f_vec eforce, consider removing
         do i = 0, kdim1 ! 29
            eforce(i, :) = eforce(i, :) + f_vec(i, :)
         end do

         ! compute elemental norms
         eresid = 0.0_RK
         resid = 0.0_RK
         part_r_norm = 0.0_RK
         r_norm = 0.0_RK
         part_rx_norm = 0.0_RK
         rx_norm = 0.0_RK
         part_f_norm = 0.0_RK
         f_norm = 0.0_RK

         call gen_matrix_vector_mult(eresid, estiff, evel, 1, 2, 3, 4, ier)
         do i = 0, kdim1 ! 29
            eresid(i,:) = eforce(i,:) - eresid(i,:)
         end do

         ! scatter elemental residuals to all nodes
         call part_scatter(resid, eresid, nodes, .false., dtrace)

         ! calculate residual and force magnitudes
         part_r_norm = sum(resid * resid)
         call par_sum(part_r_norm, r_norm)
         r_norm = dsqrt(r_norm)
         if (iter .eq. 1) r_norm_n = r_norm

         part_rx_norm = maxval(abs(resid))
         call par_max(part_rx_norm, rx_norm)

         call part_scatter(force, eforce, nodes, .false., dtrace)

         part_f_norm = sum(force * force)
         call par_sum(part_f_norm, f_norm)
         f_norm = dsqrt(f_norm)

         ! solve for new velocity field

         ! zero the residual where velocities are specified
         where (bcs) resid = 0.0_RK
         delta_vel = 0.0_RK

         if (NR) then ! use Newton-Raphson

            if (myid .eq. 0) then
               write(DFLT_U,'(A)', ADVANCE='NO') 'Info   :       . Solving NR iteration... '
            endif

            itmethod = 1
            NR_iter = NR_iter+1

            ! form the diagonal part of the stiffness matrix,
            call assemble_diagonals(gdiag, etanstiff, kdim, dof_sub1, dof_sup1, el_sub1, el_sup1, dtrace, nodes)

            ! compute the velocity field (vel) using the conjugate gradient method
            cg_iter_out = cg_solver_ebe(delta_vel, d_norm, resid, gdiag,&
                 & etanstiff, bcs, kdim, dof_sub1, dof_sup1, el_sub1, el_sup1,&
                 & cg_max_iters, cg_tol, dtrace, nodes)

         else ! use Successive Approximations

            if (myid .eq. 0) then
               write(DFLT_U,'(A)', ADVANCE='NO') 'Info   :       . Solving SA iteration... '
            endif

            itmethod = 0

            call assemble_diagonals(gdiag, estiff, kdim, dof_sub1, dof_sup1, el_sub1, el_sup1, dtrace, nodes)

            ! compute the velocity field (vel) using the conjugate gradient method
            cg_iter_out = cg_solver_ebe(delta_vel, d_norm, resid, gdiag,&
                 & estiff, bcs, kdim, dof_sub1, dof_sup1, el_sub1, el_sup1,&
                 & cg_max_iters, cg_tol, dtrace, nodes)

         endif

         vel = vel_o + delta_vel

         ! calculate velocity norm
         part_delu_norm = 0.0_RK
         delu_norm = 0.0_RK
         part_delux_norm = 0.0_RK
         delux_norm = 0.0_RK
         part_u_norm = 0.0_RK
         u_norm = 0.0_RK

         part_delu_norm = sum(delta_vel * delta_vel)
         call par_sum(part_delu_norm, delu_norm)
         delu_norm = dsqrt(delu_norm)

         part_delux_norm = maxval(abs(delta_vel))
         call par_max(part_delux_norm, delux_norm)

         part_u_norm = sum(vel_o * vel_o)
         call par_sum(part_u_norm, u_norm)
         u_norm = dsqrt(u_norm)

         ! normalize  velocity norms
         delu_norm = delu_norm/u_norm
         delux_norm = delux_norm/u_norm

         ! Newton-Raphson convergence parameters
         if (NR_iter .ge. NR_SLOPE_START) then
            NR_conv = log(delu_norm_o) - log(delu_norm)
         endif

         if (NR_iter .eq. NR_SLOPE_START) NR_conv_o = NR_conv

         ! Successive Approximation convergence parameter
         if (.not. NR) then
            if (delu_norm .gt. delu_norm_o) then
               idiv = idiv + 1
            else
               idiv = 0
            end if
         end if

         ! output convergence parameters
!         write(ounits(LOG_U), 100)  iter, cg_iter_out, d_norm, delu_norm, delux_norm
!100      FORMAT(i8,1x,i8,1x,3d14.4)

         if ((myid .eq. 0) .AND. PRINT_OPTIONS%PRINT_CONV) then
            CALL WRITE_CONV_FILE_DATA(INCR, ITER, ITMETHOD, R_NORM, RX_NORM, &
                & F_NORM, DELU_NORM, DELUX_NORM, U_NORM, CG_ITER_OUT)
         endif

         if (myid .eq. 0) then
            !
            write(DFLT_U,'(A,E10.4,A,I0,A)', ADVANCE='YES') 'R = ', delu_norm, &
                & ' (', cg_iter_out, ' iters)'
            !
         end if

         ! check convergence of the velocity solution
         !
         ! cases 1-4:  convergence, strict or loose
         ! cases 5-6:  switch back to SA due to problems with NR
         ! case 7:  failure

         ! case 1
         if ((delu_norm < nl_tol_strict) .and. (iter .gt. 1)) then
            ! solution converged
            iter_converged = .TRUE.; EXIT

         ! case 2
         else if ((delu_norm*u_norm < nl_tol_min*maxdof) .and. (iter .gt. 1)) then
            ! solution converged
            if (myid .eq. 0) then
                !
                write(DFLT_U,'(A)') 'Info   :       . Change in velocity is below &
                    &threshold value.'
                !
            end if
            !
            iter_converged = .TRUE.; EXIT

         ! case 3
         else if ((NR_iter .gt. NR_SLOPE_START) .and. (NR_conv .lt. NR_tol_conv*NR_conv_o) &
              & .and. (delu_norm < nl_tol_loose)) then
            ! Newton-Raphson slow to converge, but converged satisfactorily
            if (myid .eq. 0) then
                !
                write(DFLT_U,'(A)') 'Info   :       . Newton-Raphson is slow to &
                    &converge, but converged satisfactorily.'
                !
            end if
            !
           iter_converged = .TRUE.; EXIT
           
         ! case 4
         else if ((.not. NR) .and. (NR_attempt) .and. (delu_norm < nl_tol_loose)) then
            ! SA slow to converge, but converged satisfactorily
            if (myid .eq. 0) then
                !
                write(DFLT_U,'(A)') 'Info   :       . Successive Approximations is &
                    &slow to converge, but converged satisfactorily.'
                !
            end if
            iter_converged = .TRUE.; EXIT
         
         ! case 5
         else if (NR .and. (r_norm .gt. 1.1*r_norm_n)) then
            ! Newton-Raphson solution diverging
            if (myid .eq. 0) then
               write(DFLT_U, '(A)') 'Warning:     > Newton-Raphson solution diverging.'
            end if
            ! revert to previous SA solution and use SA
            NR = .false.
            vel = vel_SA
            ! enforce stricter switch-over tolerance between methods
            NR_tol_switch = NR_tol_switch/10.0
            ! set flag that NR had been attempted
            NR_attempt = .true.
            ! reset number of NR iterations
            NR_iter = 0

         ! case 6
         else if ((NR_iter .gt. NR_SLOPE_START) .and. (NR_conv .lt. NR_tol_conv*NR_conv_o)) then
            ! Newton-Raphson converging slowly
            if (myid .eq. 0) then
               write(DFLT_U, '(A)') 'Warning:     > Newton-Raphson converging slowly.'
            end if
            ! switch to SA
            NR = .false.
            ! set flag that NR had been attempted
            NR_attempt = .true.
            NR_slow = .true.
            ! reset number of NR iterations
            NR_iter = 0

         ! case 7
         else if (idiv .ge. 5) then
            ! Successive approximations diverging
            if (myid .eq. 0) then
               write(DFLT_U, '(A)') 'Warning:     > Iterations are diverging.'
            end if
            ! failure to converge
            itmethod_evps = -1
            RETURN
         end if

         ! switch from SA to NR
         ! removed requirement for iter .gt. 1
         if ((delu_norm .lt. NR_tol_switch) &
              & .and. (.not. NR_slow) .and. (.not. NR)) then
            NR = .true.
            idiv = 0
            ! save velocity field
            vel_SA = vel
         endif

         vel_o = vel
         r_norm_o = r_norm
         delu_norm_o = delu_norm

      ENDDO nonlinear_iteration
      ! ----------- NL-iteration loop -----------------------------------------
      IF (iter_converged) THEN
            if (myid .eq. 0) then
               write(DFLT_U, '(A,I0,A)') 'Info   :     > Converged in ', iter, ' iterations'
            end if
        !WRITE(message, *) '** INC: ', incr, ' ITSTRESS converged in ', iter, ' iterations.'
        !CALL par_message(DFLT_U, message)

        ! compute pressure
        CALL recover_pressure_evps(elpress, e_elas_kk)

      ELSE
        ! Maximum number of iterations exceeded.
!        WRITE(ounits(LOG_U), 120) iter, d_norm
!120     FORMAT(' No convergence for initial guess',/,' iter =', i4, ' d_norm =', g10.3)
        ! failure to converge
        itmethod_evps = -1
      END IF

      END FUNCTION itmethod_evps
!
!
!**********************************************************************
!
      SUBROUTINE recover_pressure_evps(elpress, e_elas_kk)
!
!----------------------------------------------------------------------
!
!     Arguments.
!
      REAL(RK), INTENT(IN)  :: e_elas_kk(el_sub1:el_sup1)
      REAL(RK), INTENT(OUT) :: elpress(el_sub1:el_sup1)
!
!     Locals:
!
      INTEGER :: iphase
      INTEGER :: my_phase(el_sub1:el_sup1)
!
!----------------------------------------------------------------------
!
      my_phase(:) = phase(el_sub1:el_sup1)

      do iphase = 1,numphases
         where (my_phase .EQ. iphase)
            elpress = - crystal_parm(8,iphase) * e_elas_kk
         end where
      enddo

    END SUBROUTINE recover_pressure_evps

END MODULE ItMethodEvpsModule
