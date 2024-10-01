! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module solveit_evp_mod

! Module handling the iteration for a single time increment for the evp soln.

  use general_mod
  use types_mod
  use conjugate_gradient_mod
  use solveit_evp_mod2
  use units_mod
  use write_res_mod
  use gather_scatter_mod
  use parallel_mod
  use finalize_res_mod

  implicit none

  public

contains

  !> Driver for the evp iteration required for a single time increment
  !> Fails if solution did not converged (either the maximum number of
  !> iterations was exceeded or the solution was badly diverging)
  subroutine solveit_evp(mesh, crys, loading, exec, printing, &
                               & results_prev, &
                               & results, dtime, incr, &
                               & load, area, wildcard)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys (:)
    type(loading_type), intent(in) :: loading
    type(exec_type), intent(inout) :: exec
    type(printing_type), intent(in) :: printing
    type(results_type), intent(in) :: results_prev
    type(results_type), intent(inout) :: results
    real(rk), intent(in) :: dtime
    integer, intent(in) :: incr
    real(rk), intent(out), optional :: load(mesh%num_fasets, 3) ! optional for triax
    real(rk), intent(out), optional :: area(mesh%num_fasets) ! optional for triax
    character(len=*), optional :: wildcard ! for triax, remove when possible

    integer :: status
    integer :: i, j, iter
    integer :: cg_iter_out
    real(rk) :: r_norm
    real(rk) :: r_norm_prev
    real(rk) :: r_norm_ini
    real(rk) :: rx_norm
    real(rk) :: f_norm
    real(rk) :: delu_norm
    real(rk) :: delu_norm_prev
    real(rk) :: delux_norm
    real(rk) :: u_norm
    real(rk) :: vel_o(dof_sub:dof_sup) ! Velocity for previous iteration
    real(rk) :: delta_vel(dof_sub:dof_sup)
    real(rk) :: vel_sa(dof_sub:dof_sup)
    real(rk), target :: estiff(kdim, kdim, elt_sub:elt_sup)
    real(rk), target :: etanstiff(kdim, kdim, elt_sub:elt_sup)
    real(rk), pointer :: stiff(:,:,:)
    real(rk) :: tcoos(dof_sub:dof_sup), ecoos(kdim, elt_sub:elt_sup)
    real(rk) :: eforce(kdim, elt_sub:elt_sup)
    real(rk) :: resid(dof_sub:dof_sup)
    real(rk) :: eresid(kdim, elt_sub:elt_sup)
    real(rk) :: gdiag(dof_sub:dof_sup)
    real(rk) :: e_ones(kdim, elt_sub:elt_sup)
    real(rk) :: g_ones(dof_sub:dof_sup)
    real(rk) :: evel(kdim, elt_sub:elt_sup)
    integer :: midnodes(3, 6)

    ! O-indexed for efficiency
    midnodes(:, 1) = [0, 2, 1]
    midnodes(:, 2) = [2, 4, 3]
    midnodes(:, 3) = [0, 4, 5]
    midnodes(:, 4) = [0, 9, 6]
    midnodes(:, 5) = [2, 9, 7]
    midnodes(:, 6) = [4, 9, 8]

    !---------------------------------------------------------------------------

    ! Scatter ones arrays and store multiplicity
    e_ones = 1.0d0
    call part_scatter(g_ones, e_ones, mesh%elt_dofs, exec%dof_trace)

    ! Initialization
    r_norm_prev = 0.0d0
    delu_norm_prev = 0.0d0
    vel_o = results%vel
    vel_sa = results%vel

    status = 0
    cg_iter_out = 0
    exec%iter_converged = .false.
    exec%itmethod = "SA"
    exec%nr_attempted = .false.
    exec%nr_isslow = .false.
    exec%num_nr_iters = 0
    exec%nr_tol_switch = exec%nr_tol_switch_ref
    exec%nr_conv = 1.0d0
    exec%nr_conv_prev = 1.0d0
    exec%num_sa_iters = 0

    nonlinear_iteration: do iter = 1, exec%nl_max_iters

      if (myid .eq. 0) then
        write (*, '(a,i0)') 'Info   :     > solveit_evp: Iteration ', iter
      end if

      ! tsh: Re-position the mid-node to the middle of the edge-nodes
      call part_gather(ecoos, results%coo, mesh%elt_dofs, exec%dof_trace)

      do i = 1, 6
        do j = 1, 3
          ecoos(3*midnodes(3, i) + j, :) = &
            & 0.5d0 * (ecoos(3*midnodes(1, i) + j, :) + ecoos(3*midnodes(2, i) + j, :))
        end do
      end do

      ! Scatter e_coordinates
      call part_scatter(results%coo, ecoos, mesh%elt_dofs, exec%dof_trace)

      ! Divide the coordinates by the multiplicity
      results%coo = results%coo/g_ones

      ! Estimate the new coordinates based on the vel
      tcoos = results%coo + dtime*results%vel

      call part_gather(ecoos, tcoos, mesh%elt_dofs, exec%dof_trace)
      call part_gather(evel, results%vel, mesh%elt_dofs, exec%dof_trace)

      ! Dnter the subroutine that, through iterations at the constitutive
      !   level , computes:
      ! The state of each element:
      ! e*_kk: e_elas_kk_qpt
      ! tau': sig_vec_qpt
      ! r*: rstar
      ! Matrices [s]: estiff
      ! Vectors [s]{h}: eforce

      ! sig_vec_n_qpt: inout (enters=0 at the first increment, only grain=0 is
      !   updated)

      call elt_stif_evp(mesh, crys, exec, results_prev, results, ecoos, evel, &
                       & estiff, etanstiff, eforce, incr, dtime)

      ! Compute elemental norms
      call gen_matrix_vector_mult(estiff, evel, eresid)
      do i = 1, kdim
        eresid(i, :) = eforce(i, :) - eresid(i, :)
      end do

      ! Scatter elemental residuals to all nodes
      call part_scatter(resid, eresid, mesh%elt_dofs, exec%dof_trace)

      ! Calculate residual and force magnitudes
      call par_sqrtsum(sum(resid*resid), r_norm)
      if (iter .eq. 1) then
        r_norm_ini = r_norm
      end if

      call par_max(maxval(abs(resid)), rx_norm)

      call part_scatter(results%force, eforce, mesh%elt_dofs, exec%dof_trace)

      call par_sqrtsum(sum(results%force*results%force), f_norm)

      ! Solve for new vel field (zero the residual where vel are specified)
      where (loading%bcs_vel_defined)
        resid = 0.0d0
      end where

      delta_vel = 0.0d0

      if (myid .eq. 0) then
        write (*, '(a,a,a)', advance='no') &
          & 'Info   :       . Solving ', exec%itmethod, ' iteration... '
      end if

      ! Newton-Raphson (etanstiff) =============================================
      if (exec%itmethod .eq. "NR") then
        stiff => etanstiff
        exec%num_nr_iters = exec%num_nr_iters + 1
      else
        stiff => estiff
      end if

      ! Form the diagonal part of the stiffness matrix
      ! If no MPC or PBC
      if ((loading%mpc_status .eqv. .false.) .and. (mesh%num_periodicity .eq. 0)) then

        call assemble_diagonals("general", loading, mesh, exec, stiff, gdiag)
        cg_iter_out = cg_solver_ebe("general", "dv", gdiag, stiff, delta_vel, resid, loading, exec, mesh)

      ! If multi-point constraints, solve the reduced linear system only on primary dofs
      else if ((loading%mpc_status .eqv. .true.) .and. (mesh%num_periodicity .eq. 0)) then
        ! Form the diagonal part of the primary stiffness matrix (Jacobi's preconditioning)
        call assemble_diagonals("MPC", loading, mesh, exec, stiff, gdiag)
        ! Solve the linear system using the congujate gradient method
        cg_iter_out = cg_solver_ebe("MPC", "dv", gdiag, stiff, delta_vel, resid, loading, exec, mesh)

      ! If PBC
      else if ((loading%mpc_status .eqv. .false.) .and. (mesh%num_periodicity .gt. 0)) then

        call assemble_diagonals("PBC", loading, mesh, exec, stiff, gdiag)
        cg_iter_out = cg_solver_ebe("PBC", "dv", gdiag, stiff, delta_vel, resid, loading, exec, mesh)

      end if

      results%vel = vel_o + delta_vel

      ! Calculate vel norm
      call par_sqrtsum(sum(delta_vel*delta_vel), delu_norm)

      call par_max(maxval(abs(delta_vel)), delux_norm)

      call par_sqrtsum(sum(vel_o*vel_o), u_norm)

      ! Normalize vel norms
      delu_norm = delu_norm / u_norm
      delux_norm = delux_norm / u_norm

      ! Newton-Raphson convergence parameters
      if (exec%itmethod .eq. "NR") then
        if (exec%num_nr_iters .eq. exec%nr_slope_start) then
          exec%nr_conv_prev = exec%nr_conv
        else if (exec%num_nr_iters .gt. exec%nr_slope_start) then
          exec%nr_conv = log(delu_norm_prev) - log(delu_norm)
        end if
      end if

      ! Successive Approximation convergence parameter
      if (exec%itmethod .eq. "SA") then
        if (delu_norm .gt. delu_norm_prev) then
          exec%num_sa_iters = exec%num_sa_iters + 1
        else
          exec%num_sa_iters = 0
        end if
      end if

      ! Output convergence parameters
      if ((myid .eq. 0) .and. printing%print_conv) then
        call write_conv_file_data(printing, incr, iter, exec%itmethod, r_norm, rx_norm, &
            & f_norm, delu_norm, delux_norm, u_norm, cg_iter_out)
      end if

      if (myid .eq. 0) then
        write (*, '(a,e10.4,a,i0,a)', advance='yes') 'R = ', delu_norm, &
            & ' (', cg_iter_out, ' iters)'
      end if

      ! Check convergence of the vel solution
      ! Cases 1-4:  convergence, strict or loose
      ! Cases 5-6:  switch back to SA due to problems with NR
      ! Case 7:  failure

      ! Case 1: solution converged
      if ((delu_norm < exec%nl_tol_strict) .and. (iter .gt. 1)) then
        exec%iter_converged = .true.

        exit

      else if ((delu_norm*u_norm < exec%nl_tol_min*3*mesh%num_nodes) .and. &
          & (iter .gt. 1)) then
        ! Case 2
        ! Solution converged

        if (myid .eq. 0) then
          write (*, '(a)') 'Info   :       . Change in vel is &
              &below threshold value.'
        end if

        exec%iter_converged = .true.

        exit

      else if ((exec%num_nr_iters .gt. exec%nr_slope_start) .and. (exec%nr_conv .lt. &
          & exec%nr_tol_conv*exec%nr_conv_prev) .and. (delu_norm < exec%nl_tol_loose)) then
        ! Case 3
        ! Newton-Raphson slow to converge, but converged satisfactorily

        if (myid .eq. 0) then
          write (*, '(a)') 'Info   :       . Newton-Raphson is slow &
              &to converge, but converged satisfactorily.'
        end if

        exec%iter_converged = .true.

        exit

      else if ((exec%itmethod .eq. "SA") .and. (exec%nr_attempted) .and. &
          & (delu_norm < exec%nl_tol_loose)) then
        ! Case 4
        ! SA slow to converge, but converged satisfactorily

        if (myid .eq. 0) then
          write (*, '(a)') 'Info   :       . Successive Approximations&
              & is slow to converge, but converged satisfactorily.'
        end if

        exec%iter_converged = .true.

        exit

      ! Case 5: Newton-Raphson solution diverging
      else if (exec%itmethod .eq. "NR" .and. (r_norm .gt. 1.1*r_norm_ini)) then

        if (myid .eq. 0) then
          write (*, '(a)') 'Warning:     > Newton-Raphson solution &
              &diverging.'
        end if

        ! Revert to previous SA solution and use SA
        results%vel = vel_sa
        exec%itmethod = "SA"

        ! Enforce stricter switch-over tolerance between methods
        exec%nr_tol_switch = exec%nr_tol_switch/10.0

        ! Set flag that NR had been attempted
        exec%nr_attempted = .true.

        ! Reset number of NR iterations
        exec%num_nr_iters = 0

      else if ((exec%num_nr_iters .gt. exec%nr_slope_start) &
       & .and. (exec%nr_conv .lt. exec%nr_tol_conv*exec%nr_conv_prev)) then
        ! Case 6
        ! Newton-Raphson converging slowly

        if (myid .eq. 0) then
          write (*, '(a)') 'Warning:     > Newton-Raphson converging slowly.'
        end if

        ! Switch to SA
        exec%itmethod = "SA"

        ! Set flag that NR had been attempted
        exec%nr_attempted = .true.
        exec%nr_isslow = .true.

        ! Reset number of NR iterations

        exec%num_nr_iters = 0

      else if (exec%num_sa_iters .ge. 5) then
        ! Case 7
        ! Successive approximations diverging

        if (myid .eq. 0) then
          write (*, '(a)') 'Warning:     > Iterations are diverging.'
        end if

        ! Failure to converge

        status = -1

        return
      end if

      ! Switch from SA to NR
      ! Removed requirement for iter .gt. 1
      if ((delu_norm .lt. exec%nr_tol_switch) .and. (.not. exec%nr_isslow) .and. &
          & (exec%itmethod .eq. "SA")) then
        exec%itmethod = "NR"
        exec%num_sa_iters = 0

        ! Save vel field
        vel_sa = results%vel
      end if

      vel_o = results%vel

      r_norm_prev = r_norm
      delu_norm_prev = delu_norm

    end do nonlinear_iteration

    if (exec%iter_converged) then
      if (myid .eq. 0) then
        write (*, '(a,i0,a)') 'Info   :     > Converged in ', iter, ' iterations'
      end if

      ! rq 03/2023: commenting, as this as no effect
      ! call recover_pressure_evp(mesh, crys, elpress, e_elas_kk_qpt)

    ! Maximum number of iterations exceeded. Failure to converge.
    else
      status = -1

    end if

    if (status .ne. 0) then
      call par_quit('Error  :     > Failure to converge.',clock_start=exec%clock_start)
    end if

    ! Update state (at center of each element), geometry,
    ! coo @(t+dt), using dt and vel @(t+dt)
    call update_state_evp(mesh, crys, exec, results_prev, results, dtime)

    if (.not. (present(wildcard) .and. wildcard .eq. "nofinalize")) then
      call finalize_res (mesh, crys, dtime, results_prev, results, load, area)
    end if

  end subroutine solveit_evp

end module solveit_evp_mod
