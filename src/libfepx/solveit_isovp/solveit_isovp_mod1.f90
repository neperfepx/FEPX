! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module solveit_isovp_mod

! Module handling the iteration for a single time increment for the viscoplastic
!   solution.

! Contains subroutines:
! recover_pressure_vp: Compute pressure from vel field.

! Contains functions:
! solveit_vp: Driver for the viscoplastic iteration required for a single time
!   increment.

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use solveit_isovp_mod2
  use types_mod
  use conjugate_gradient_mod
  use surface_mod
  use units_mod
  use gather_scatter_mod
  use parallel_mod

  implicit none

  public

contains

  !> Driver for the isotropic viscoplastic solution
  subroutine solveit_vp(mesh, loading, exec, results)

    type(mesh_type), intent(in) :: mesh
    type(loading_type), intent(in) :: loading
    type(exec_type), intent(in) :: exec
    type(results_type), intent(inout) :: results

    integer :: iter
    integer :: i
    integer :: cg_iter_out
    real(rk) :: part_u_norm
    real(rk) :: u_norm
    real(rk) :: estiff(kdim, kdim, elt_sub:elt_sup)
    real(rk) :: ecoos(kdim, elt_sub:elt_sup)
    real(rk) :: my_force(dof_sub:dof_sup)
    real(rk) :: eforce(kdim, elt_sub:elt_sup)
    real(rk) :: vel_o(dof_sub:dof_sup)
    real(rk) :: gdiag(dof_sub:dof_sup)
    real(rk) :: pscale(elt_sub:elt_sup)
    real(rk) :: pcnst(kdim, elt_sub:elt_sup)
    real(rk) :: elpress(elt_sub:elt_sup)
    real(rk) :: evel(kdim, elt_sub:elt_sup)
    real(rk) :: epseff(elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    elpress = 0.0d0

    vel_o = results%vel

    call part_gather(evel, results%vel, mesh%elt_dofs, exec%dof_trace)
    call part_gather(ecoos, results%coo, mesh%elt_dofs, exec%dof_trace)

    ! Non linear iteration loop
    do iter = 1, 3

      if (myid .eq. 0) then
        write (*, '(a,i0)') 'Info   :     > solveit_vp: Iteration ', &
            & iter
        write (*, '(a)', advance='no') 'Info   :       . Solving NL &
            &iteration... '
      end if

      call elt_stif_vp("isotropic_vp", estiff, ecoos, evel, pscale, pcnst, epseff)

      do i = 1, kdim
        eforce(i, :) = -elpress*pcnst(i, :)
      end do

      ! eforce --> my_force
      call part_scatter(my_force, eforce, mesh%elt_dofs, exec%dof_trace)

      ! Zero the forces where vel are specified
      where (loading%bcs_vel_defined)
        my_force = 0.0d0
      end where

      ! Preconditioned conjugate gradient solver

      ! If no MPC or PBC
      if ((loading%mpc_status .eqv. .false.) .and. (mesh%num_periodicity .eq. 0)) then

      ! Form the diagonal part of the stiffness matrix
        call assemble_diagonals("general", loading, mesh, exec, estiff, gdiag)
      ! Compute the vel field (vel) using the conjugate gradient method
        cg_iter_out = cg_solver_ebe("general", "v", gdiag, estiff, results%vel, my_force, loading, exec, mesh)
      
      ! If multi-point constraints, solve the reduced linear system only on primary dofs
      else if ((loading%mpc_status .eqv. .true.) .and. (mesh%num_periodicity .eq. 0)) then
        ! Form the diagonal part of the primary stiffness matrix (Jacobi's preconditioning)
        call assemble_diagonals("general", loading, mesh, exec, estiff, gdiag)
        ! Solve the linear system using the congujate gradient method
        cg_iter_out = cg_solver_ebe("general", "v", gdiag, estiff, results%vel, my_force, loading, exec, mesh)
    
      ! If PBC
      else if ((loading%mpc_status .eqv. .false.) .and. (mesh%num_periodicity .gt. 0)) then
          
        call assemble_diagonals("PBC", loading, mesh, exec, estiff, gdiag)
        cg_iter_out = cg_solver_ebe("PBC", "v", gdiag, estiff, results%vel, my_force, loading, exec, mesh)

      end if

      ! vel --> evel
      call part_gather(evel, results%vel, mesh%elt_dofs, exec%dof_trace)

      call recover_pressure_vp(evel, pscale, elpress, pcnst)

      part_u_norm = sum((vel_o - results%vel)*(vel_o - results%vel))
      call par_sum(part_u_norm, u_norm)
      u_norm = dsqrt(u_norm)

      if (myid .eq. 0) then
        write (*, '(a,e10.4,a,i0,a)', advance='yes') 'R = ', u_norm, &
            & ' (', cg_iter_out, ' iters)'
          if (u_norm .ne. u_norm) then
            write (*,*) "Error  :         Residual is NaN."
            call abort ()
          end if
      end if

      vel_o = results%vel

    end do

    if (myid .eq. 0) write (*, '(a,i0,a)') 'Info   :     > Converged in ',&
        & iter - 1, ' iterations'

    ! if (status .lt. 0) then
    !   call par_quit('Error  :     > solveit_vp failed to converge.')
    ! end if

  end subroutine solveit_vp

end module solveit_isovp_mod
