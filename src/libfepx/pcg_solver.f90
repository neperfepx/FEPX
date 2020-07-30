! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
!**********************************************************************
!
!     Routines for the preconditioned conjugate gradient solver.
!
!     assemble_diagonals : Form the diagonal part of the stiffness matrix,
!     *                  : for use in preconditioning.
!     cg_solver_ebe      : Driver for element by element conjugate gradient
!     .                  : solver.
!
!**********************************************************************
!
      SUBROUTINE assemble_diagonals(&
  &      diagonals, gstif, nnpe, nsub, nsup, esub, esup, dtrace, np)
!
!     Form the diagonal part of the stiffness matrix, for use
!     in preconditioning.
!
!----------------------------------------------------------------------
!
      USE gather_scatter
      USE IntrinsicTypesModule, RK=>REAL_KIND
!
      IMPLICIT NONE
!
!     Arguments.
!
!     diagonals  : diagonal matrix (result of this routine)
!     gstif      : global stiffness matrix
!     nnpe       : number of nodes per element
!     nsub, nsup : node range for this process
!     esub, esup : element range for this process
!     dtrace     : commmunication information for gather/scatter
!     np         : connectivity
!
      TYPE(trace)  dtrace
!
      INTEGER nnpe, nsub, nsup, esub, esup
      INTEGER np(0:(nnpe-1), esub:esup)
!
      REAL(RK) diagonals(nsub:nsup)
      REAL(RK) gstif(0:(nnpe - 1), 0:(nnpe - 1), esub:esup)
!
!     Locals.
!
      INTEGER i, j
      REAL(RK)  ediagonals(0:(nnpe - 1), esub:esup)
      REAL(RK)  diagmin, diagmax
!
!----------------------------------------------------------------------
!
!  RC 3/24/2016: Reordered for better memory striding
      do j = esub, esup
         do i = 0,(nnpe - 1)
            ediagonals(i, j) = gstif(i, i, j)
         enddo
      end do

      diagonals = 0.0
      
      call part_scatter(diagonals, ediagonals, np, .false., dtrace)

      diagonals = 1.0 / diagonals

      RETURN
      END
!
!**********************************************************************
!
      INTEGER FUNCTION cg_solver_ebe(&
  &   sol, res_norm, rhs, diagonals, gstif, bcs, nnpe, nsub, nsup,&
  &   esub, esup, max_iter, tol, dtrace, np)
!
!     Driver for element by element conjugate gradient solver.
!
!----------------------------------------------------------------------
!
!     Modules:
!
      USE gather_scatter
      USE parallel_mod
      USE parallel_matrix_mod
!
      IMPLICIT NONE
!
!     Arguments.
!
      INTEGER   max_iter, nnpe, nsub, nsup, esub, esup
      INTEGER   np(0:(nnpe-1), esub:esup)
!
      TYPE(trace) dtrace
!
      LOGICAL bcs(nsub:nsup)
!
      REAL(RK)   part_res_norm, res_norm, tol
      REAL(RK)   sol(nsub:nsup), rhs(nsub:nsup), diagonals(nsub:nsup)
      REAL(RK)   gstif(0:(nnpe - 1), 0:(nnpe - 1), esub:esup)
!
!     Locals.
!
      INTEGER iter, n_iter, ier
! tsh for interation
      INTEGER i, j
!
      REAL(RK)     elapsed
!
      REAL(RK)&
  &   zu (nsub:nsup), ru(nsub:nsup),&
  &   apu(nsub:nsup), pu(nsub:nsup)
      REAL(RK)&
  &   temp1(0:(nnpe - 1), esub:esup),&
  &   temp2(0:(nnpe - 1), esub:esup)
      REAL(RK)&
  &   part_alpha, alpha,&
  &   beta,&
  &   part_error, error,&
  &   part_xnumer, xnumer, xnumer_o,&
  &   part_mag, mag
!
!----------------------------------------------------------------------
!
!     CG initializations.
!
      call sparse_matvec_ebe(ru, sol, temp1, temp2, gstif, bcs, nnpe,&
  &        nsub, nsup, esub, esup, dtrace, np)

      ru  = rhs - ru
      rhs = rhs * diagonals
      zu  = ru * diagonals
      
      part_mag = sum(zu * zu)
      call par_sum(part_mag, mag)
      mag = dsqrt(mag)
!
!     CG iterations.
!
      xnumer_o = 1.0d0
      pu = 0.0d0
      
      error = 1.0d0
      
      n_iter = 0
      do while (error .gt. tol)

         n_iter = n_iter + 1

         if (n_iter .gt. max_iter) then
            call par_quit('Error  :       . CG_SOLVER_EBE: Convergence failure.')        
         end if
         
         part_xnumer = sum(ru * diagonals * ru)
         call par_sum(part_xnumer, xnumer)
         beta = xnumer / xnumer_o
         
         if (n_iter .eq. 1) beta = 0.0d0
         
         xnumer_o = xnumer

         pu = zu + beta * pu
         
         call sparse_matvec_ebe(apu, pu, temp1, temp2, gstif, bcs,&
  &           nnpe, nsub, nsup, esub, esup, dtrace, np)

         part_alpha = sum(pu * apu)
         call par_sum(part_alpha, alpha)
         alpha = xnumer / alpha

         sol = sol + alpha * pu
         ru = ru - alpha * apu
         zu = ru * diagonals

         part_error = sum(zu * zu)
         call par_sum(part_error, error)
         error = sqrt(error) / mag

      enddo

      part_res_norm = sum(sol * sol)
      call par_sum(part_res_norm, res_norm)
      res_norm = (res_norm)**0.5

      cg_solver_ebe = n_iter
 
      RETURN
      END
!      
