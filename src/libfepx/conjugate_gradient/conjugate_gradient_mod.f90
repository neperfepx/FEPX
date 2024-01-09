! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module conjugate_gradient_mod

! Routines for the preconditioned conjugate gradient solver

! Contains subroutines:
! assemble_diagonals: Form the diagonal part of the stiffness matrix for use in
!   preconditioning
! assemble_diagonals_ps: Form the diagonal part of the stiffness matrix with mpc for use in
!   preconditioning

! Contains functions:
! cg_solver_ebe: Driver for element by element conjugate gradient solver

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod

! From libfepx:

  use matrix_operations_mod

  use gather_scatter_mod
  use parallel_mod
  use conjugate_gradient_mod2

  implicit none

  public

contains

  subroutine assemble_diagonals(type_bc, loading, mesh, exec, gstif, gdiag)

    ! Form the diagonal part of the stiffness matrix, for use in
    !   preconditioning.
    !---------------------------------------------------------------------------

    character(len=*), intent(in) :: type_bc
    type(loading_type), intent(in) :: loading
    type(mesh_type), intent(in) :: mesh
    type(exec_type), intent(in) :: exec
    real(rk), intent(in) :: gstif(kdim, kdim, elt_sub:elt_sup)
    real(rk), intent(out) :: gdiag(dof_sub:dof_sup)

    ! Locals:

    integer :: i, j, ii

    real(rk) :: ediagonals(kdim, elt_sub:elt_sup)

    real(rk) :: diagonals_ebe(kdim, elt_sub:elt_sup)

    real(rk) :: diagonals_ps_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: diagonals_primary_drive_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: diagonals_secondary_drive_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: diagonals_ps(dof_sub:dof_sup)
    real(rk) :: diagonals_primary_drive(dof_sub:dof_sup)
    real(rk) :: diagonals_secondary_drive(dof_sub:dof_sup)

    logical :: glogic(kdim, kdim, elt_sub:elt_sup)

    glogic = .false.


    select case(type_bc)

    ! System without MPC or PBC ------------------------------------------------------------------------
    case("general")

    ! rc 3/24/2016: Reordered for better memory striding

    do j = elt_sub, elt_sup
      do i = 1, kdim
        ediagonals(i, j) = gstif(i, i, j)
      end do
    end do

    call part_scatter(gdiag, ediagonals, mesh%elt_dofs, exec%dof_trace)

    gdiag = 1.0d0/gdiag


    ! System with MPC ----------------------------------------------------------------------------------
    case("MPC")

    ! Form the diagonal part of the stiffness matrix considering MPCs, for use in
    !   preconditioning.

    diagonals_ebe = 0.0d0

    do j = elt_sub, elt_sup
      do i = 1, kdim
        if (loading%conn_mpc(i,j) .eq. mesh%elt_dofs(i,j)) then
          diagonals_ebe(i,j) = gstif(i,i,j)
        else if (loading%conn_mpc(i,j) .ne. mesh%elt_dofs(i,j)) then
          diagonals_ebe(i,j) = gstif(i,i,j)

          do ii = 1, kdim
            if (mesh%elt_dofs(ii,j) .eq. loading%conn_mpc(i,j)) then
              diagonals_ebe(i,j) = diagonals_ebe(i,j) + gstif(i,ii,j) + gstif(ii,i,j)
            else if ((loading%conn_mpc(ii,j) .eq. loading%conn_mpc(i,j)) .and. (ii .ne. i)) then
              diagonals_ebe(i,j) = diagonals_ebe(i,j) + (gstif(i,ii,j) + gstif(ii,i,j)) / 2
            end if
          end do
        end if
      end do
    end do

    call part_scatter(gdiag, diagonals_ebe, int(loading%conn_mpc), loading%mpc_trace)

    where (gdiag .ne. 0.0d0)
      gdiag = 1.0d0/gdiag
    end where


    ! System with PBC ----------------------------------------------------------------------------------
    case("PBC")

    diagonals_ps_ebe = 0.0d0
    diagonals_primary_drive_ebe = 0.0d0
    diagonals_secondary_drive_ebe = 0.0d0
    diagonals_ps = 0.0d0
    diagonals_primary_drive = 0.0d0
    diagonals_secondary_drive = 0.0d0

    do j = elt_sub, elt_sup
      do i = 1, kdim

      !> Contribution of primary-secondary relations (conn_mpc)

        ! Prevent multiple summation of the same component
        if (glogic(i,i,j) .eqv. .false.) then
          diagonals_ps_ebe(i,j) = gstif(i,i,j)  ! for all dofs (primary and secondary)
          glogic(i,i,j) = .true.
        end if

        ! if secondary dof
        if (loading%conn_mpc(i,j) .ne. mesh%elt_dofs(i,j)) then
          ! parse element dofs to find others dofs which share same p-s relation
          do ii = 1, kdim
            ! if primary dof is on element j, add cross components
            if (mesh%elt_dofs(ii,j) .eq. loading%conn_mpc(i,j)) then
              diagonals_ps_ebe(i,j) = diagonals_ps_ebe(i,j) + gstif(i,ii,j) + gstif(ii,i,j)
            ! if others secondary elements, add half cross components (du to symmetric roles of secondary dofs)
            else if ((loading%conn_mpc(ii,j) .eq. loading%conn_mpc(i,j)) .and. (ii .ne. i)) then
              diagonals_ps_ebe(i,j) = diagonals_ps_ebe(i,j) + (gstif(i,ii,j) + gstif(ii,i,j)) / 2
            end if
          end do
        end if

      !> Contribution of primary drive dofs (primary_drive_ebe)

        ! if primary not drive and not imposed dof
        if ((loading%primary_drive_ebe(i,j) .ne. loading%conn_mpc(i,j)) .and. (loading%imposed_state_ebe(i,j) .eq. 0)) then

          ! Prevent multiple summation of the same component
          if (glogic(i,i,j) .eqv. .false.) then
            diagonals_ps_ebe(i,j) = gstif(i,i,j)
            glogic(i,i,j) = .true.
          end if

          ! parse element dofs to find others dofs which share a relation
          do ii = 1, 20
            ! if other primary not drive and not imposed dofs
            if ((ii .ne. i) .and. &
                &(loading%primary_drive_ebe(ii,j) .ne. loading%conn_mpc(ii,j)) .and. (loading%imposed_state_ebe(ii,j) .eq. 0)) then
              diagonals_primary_drive_ebe(i,j) = diagonals_primary_drive_ebe(i,j) + &
                & sign(1.0d0,loading%label_pv_ebe(i,j)*loading%label_pv_drive_ebe(i,j))* &
                & sign(1.0d0,loading%label_pv_ebe(ii,j)*loading%label_pv_drive_ebe(ii,j))*(gstif(i,ii,j) + gstif(ii,i,j)) / 2

            ! if primary drive dof or its secondary non drive dofs
            else if (loading%conn_mpc(ii,j) .eq. (loading%primary_drive_ebe(i,j)) .and. (loading%label_pv_ebe(ii,j) .eq. 0)) then
              diagonals_primary_drive_ebe(i,j) = diagonals_primary_drive_ebe(i,j) + &
                & sign(1.0d0,loading%label_pv_ebe(ii,j)*loading%label_pv_drive_ebe(ii,j))*(gstif(i,ii,j) + gstif(ii,i,j))

            end if

          end do
        end if

      !> Contribution of secondary drive dofs (secondary_drive_ebe)

        ! if secondary not drive and not imposed dofs
        if ((loading%secondary_drive_ebe(i,j) .ne. mesh%elt_dofs(i,j)) .and. (loading%imposed_state_ebe(i,j) .eq. 0)) then

          ! Prevent multiple summation of the same component
          if (glogic(i,i,j) .eqv. .false.) then
            diagonals_ps_ebe(i,j) = gstif(i,i,j)
            glogic(i,i,j) = .true.
          end if

          ! parse element dofs to find others dofs which share a relation
          do ii = 1, kdim
            ! if other secondary not drive and not imposed dofs
            if ((ii .ne. i) .and. &
                &(loading%secondary_drive_ebe(ii,j) .ne. mesh%elt_dofs(ii,j)) .and. (loading%imposed_state_ebe(ii,j) .eq. 0)) then
              diagonals_secondary_drive_ebe(i,j) = diagonals_secondary_drive_ebe(i,j) + &
                & sign(1.0d0,loading%label_pv_ebe(i,j)*loading%label_pv_drive_ebe(i,j))* &
                & sign(1.0d0,loading%label_pv_ebe(ii,j)*loading%label_pv_drive_ebe(ii,j))*(gstif(i,ii,j) + gstif(ii,i,j)) / 2

            ! if primary drive dof or its secondary non drive dofs
            else if (loading%conn_mpc(ii,j) .eq. (loading%secondary_drive_ebe(i,j)) .and. (loading%label_pv_ebe(ii,j) .eq. 0)) then
              diagonals_secondary_drive_ebe(i,j) = diagonals_secondary_drive_ebe(i,j) + &
                & sign(1.0d0,loading%label_pv_ebe(ii,j)*loading%label_pv_drive_ebe(ii,j))*(gstif(i,ii,j) + gstif(ii,i,j))
            end if

          end do
        end if

      end do
    end do

    call part_scatter(diagonals_ps, diagonals_ps_ebe, int(loading%conn_mpc), loading%mpc_trace)
    call part_scatter(diagonals_primary_drive, diagonals_primary_drive_ebe, &
                    & int(loading%primary_drive_ebe), loading%primary_drive_trace)
    call part_scatter(diagonals_secondary_drive, diagonals_secondary_drive_ebe, &
                    & int(loading%secondary_drive_ebe), loading%secondary_drive_trace)

   gdiag = diagonals_ps + diagonals_secondary_drive - diagonals_primary_drive

    where (gdiag .ne. 0.0d0)
      gdiag = 1.0d0/gdiag
    end where

  end select

    return

  end subroutine assemble_diagonals

  integer function cg_solver_ebe(type_bc, type_system, gdiag, gstif, delta_vel, rhs, loading, exec, mesh, res_norm)

    type(loading_type), intent(in) :: loading
    type(exec_type), intent(in) :: exec
    type(mesh_type), intent(in) :: mesh

    ! Driver for element by element conjugate gradient solver considering multi-point constraints
    ! (mpc_status = .true.) or not.

    ! Note for the solver strategy for problem with mpcs:
    ! Objective -> find the increment of velocity dv (=delta_vel) leading to df (=rhs) the residual.
    ! i.e. solving K'.dv'=df'   with K'=transpose(T).K.T the reduced stifness matrix obtained from
    ! the global stifness matrix K (=gstif) and T (=mpc_mat),
    ! the primary-secondary connectivity matrix,
    ! dv=T.dv' and df'=transpose(T).df.

    !---------------------------------------------------------------------------

    character(len=*), intent(in) :: type_bc
    character(len=*), intent(in) :: type_system
    real(rk), intent(inout) :: delta_vel(dof_sub:dof_sup)
    real(rk), intent(inout) :: rhs(dof_sub:dof_sup)
    real(rk), intent(in) :: gdiag(dof_sub:dof_sup)
    real(rk), intent(in) :: gstif(kdim, kdim, elt_sub:elt_sup)
    real(rk), intent(out), optional :: res_norm

    ! Locals:

    real(rk) :: part_res_norm
    integer :: num_iters
    real(rk) :: ru(dof_sub:dof_sup)
    real(rk) :: apu(dof_sub:dof_sup)
    real(rk) :: pu(dof_sub:dof_sup)
    real(rk) :: zu(dof_sub:dof_sup)
    real(rk) :: temp1(kdim, elt_sub:elt_sup)
    real(rk) :: temp2(kdim, elt_sub:elt_sup)
    real(rk) :: part_alpha
    real(rk) :: alpha
    real(rk) :: beta
    real(rk) :: part_error
    real(rk) :: error
    real(rk) :: part_xnumer
    real(rk) :: xnumer
    real(rk) :: xnumer_o
    real(rk) :: part_mag
    real(rk) :: mag
    ! reduced arrays at primary node dofs

    real(rk) :: sol_p(dof_sub:dof_sup)
    real(rk) :: rhs_p(dof_sub:dof_sup)
    real(rk) :: zu_p(dof_sub:dof_sup)
    real(rk) :: ru_p(dof_sub:dof_sup)
    real(rk) :: apu_p(dof_sub:dof_sup)
    real(rk) :: pu_p(dof_sub:dof_sup)
    real(rk) :: rhs_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: ru_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: apu_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: pu_p_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: buff_rhs(dof_sub:dof_sup)
    real(rk) :: buff_ru(dof_sub:dof_sup)
    real(rk) :: buff_apu(dof_sub:dof_sup)
    real(rk) :: sol_p_ebe(kdim, elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    if (present(res_norm)) res_norm = 0.0d0

      sol_p = 0.0d0
      ru_p = 0.0d0
      rhs_p = 0.0d0
      zu_p = 0.0d0
      apu_p = 0.0d0
      pu_p = 0.0d0
      rhs_ebe = 0.0d0
      ru_ebe = 0.0d0
      apu_ebe = 0.0d0
      pu_p_ebe = 0.0d0
      buff_rhs = 0.0d0
      buff_ru = 0.0d0
      buff_apu = 0.0d0
      sol_p_ebe = 0.0d0

    ! cg initializations
    ! delta_vel is normally equal to zero as input.
    call sparse_matvec_ebe(ru, delta_vel, temp1, temp2, gstif, loading%bcs_vel_defined, kdim, dof_sub, &
        & dof_sup, elt_sub, elt_sup, exec%dof_trace, mesh%elt_dofs)

    select case(type_bc)

    ! System without MPC or PBC ------------------------------------------------------------------------------
    case("general")

    ru = rhs - ru
    rhs = rhs*gdiag
    zu = ru*gdiag

    part_mag = sum(zu*zu)

    call par_sum(part_mag, mag)

    mag = dsqrt(mag)

    ! cg iterations

    xnumer_o = 1.0d0
    pu = 0.0d0
    error = 1.0d0
    num_iters = 0

    do while (error .gt. exec%cg_tol)
      num_iters = num_iters + 1

      if (num_iters .gt. exec%cg_max_iters) then
        call PAR_quit('error  : cg_solver_ebe: Max. iterations reached.')
      end if

      part_xnumer = sum(ru*gdiag*ru)

      call par_sum(part_xnumer, xnumer)

      beta = xnumer/xnumer_o

      if (num_iters .eq. 1) beta = 0.0d0

      xnumer_o = xnumer

      pu = zu + beta*pu

      call sparse_matvec_ebe(apu, pu, temp1, temp2, gstif, loading%bcs_vel_defined, kdim, dof_sub, &
          & dof_sup, elt_sub, elt_sup, exec%dof_trace, mesh%elt_dofs)

      part_alpha = sum(pu*apu)

      call par_sum(part_alpha, alpha)

      alpha = xnumer/alpha

      delta_vel = delta_vel + alpha*pu
      ru = ru - alpha*apu
      zu = ru*gdiag

      part_error = sum(zu*zu)

      call par_sum(part_error, error)

      error = dsqrt(error)/mag
    end do

    part_res_norm = sum(delta_vel*delta_vel)


    ! System with MPC ----------------------------------------------------------------------------------
    case("MPC")

      where (mesh%g_ones .ne. 0)
        buff_rhs=rhs/mesh%g_ones
        buff_ru=ru/mesh%g_ones
      end where

      call part_gather(rhs_ebe, buff_rhs, mesh%elt_dofs, exec%dof_trace)
      call part_gather(ru_ebe, buff_ru, mesh%elt_dofs, exec%dof_trace)
      call part_scatter(rhs_p,rhs_ebe, int(loading%conn_mpc), loading%mpc_trace)
      call part_scatter(ru_p,ru_ebe, int(loading%conn_mpc), loading%mpc_trace)

    ru_p = rhs_p - ru_p
    rhs_p = rhs_p*gdiag
    zu_p = ru_p*gdiag

    part_mag = sum(zu_p*zu_p)

    call par_sum(part_mag, mag)

    mag = dsqrt(mag)

    ! cg iterations

    xnumer_o = 1.0d0
    pu = 0.0d0
    error = 1.0d0
    num_iters = 0

    do while (error .gt. exec%cg_tol)
      num_iters = num_iters + 1

      if (num_iters .gt. exec%cg_max_iters) then
        call PAR_quit('error  : cg_solver_ebe: Max. iterations reached.')
      end if

      part_xnumer = sum(ru_p*gdiag*ru_p)

      call par_sum(part_xnumer, xnumer)

      beta = xnumer/xnumer_o

      if (num_iters .eq. 1) beta = 0.0d0

      xnumer_o = xnumer

      pu_p = zu_p + beta*pu_p

        call part_gather(pu_p_ebe, pu_p, int(loading%conn_mpc), loading%mpc_trace)
        call part_scatter(pu, pu_p_ebe, mesh%elt_dofs, exec%dof_trace)

        pu = loading%coeff_ps*pu/mesh%g_ones

      call sparse_matvec_ebe(apu, pu, temp1, temp2, gstif, loading%bcs_vel_defined, kdim, dof_sub, &
        & dof_sup, elt_sub, elt_sup, exec%dof_trace, mesh%elt_dofs)

        apu_p = 0.0d0

        where (mesh%g_ones .ne. 0)
          buff_apu=apu/mesh%g_ones
        end where

        call part_gather(apu_ebe, buff_apu, mesh%elt_dofs, exec%dof_trace)
        call part_scatter(apu_p,apu_ebe, int(loading%conn_mpc), loading%mpc_trace)

      part_alpha = sum(pu_p*apu_p)

      call par_sum(part_alpha, alpha)

      alpha = xnumer/alpha

      sol_p = sol_p + alpha*pu_p
      ru_p = ru_p - alpha*apu_p
      zu_p = ru_p*gdiag

      part_error = sum(zu_p*zu_p)

      call par_sum(part_error, error)

      error = dsqrt(error)/mag

    end do

      call part_gather(sol_p_ebe, sol_p, int(loading%conn_mpc), loading%mpc_trace)
      call part_scatter(delta_vel, sol_p_ebe, mesh%elt_dofs, exec%dof_trace)

      delta_vel = loading%coeff_ps*delta_vel/mesh%g_ones
      part_res_norm = sum(sol_p*sol_p)


    ! System with PBC ----------------------------------------------------------------------------------
    case("PBC")

      call all_to_prim("PBC", rhs_p, rhs, loading, mesh, exec)
      call all_to_prim("PBC", ru_p, ru, loading, mesh, exec)

      ru_p = rhs_p - ru_p
      rhs_p = rhs_p*gdiag

      where (loading%bcs_vel_defined .eqv. .true.)
        ru_p = 0.0d0
      end where

      zu_p = ru_p*gdiag

      part_mag = sum(zu_p*zu_p)

      call par_sum(part_mag, mag)

      mag = dsqrt(mag)

      ! cg iterations

      xnumer_o = 1.0d0
      pu = 0.0d0
      error = 1.0d0
      num_iters = 0

      do while (error .gt. exec%cg_tol)
        num_iters = num_iters + 1

        if (num_iters .gt. exec%cg_max_iters) then
          call PAR_quit('error  : cg_solver_ebe: Max. iterations reached.')
        end if

        part_xnumer = sum(ru_p*gdiag*ru_p)
        call par_sum(part_xnumer, xnumer)

        beta = xnumer/xnumer_o

        if (num_iters .eq. 1) beta = 0.0d0

        xnumer_o = xnumer

        pu_p = zu_p + beta*pu_p

        call prim_to_all("PBC", type_system, pu, pu_p, loading, mesh, exec)

        if (type_system .eq. "v") then
          pu = pu - loading%offset_ps
        end if

        call sparse_matvec_ebe(apu, pu, temp1, temp2, gstif, loading%bcs_vel_defined, kdim, dof_sub, &
              & dof_sup, elt_sub, elt_sup, exec%dof_trace, mesh%elt_dofs)

        call all_to_prim("PBC", apu_p, apu, loading, mesh, exec)

        part_alpha = sum(pu_p*apu_p)

        call par_sum(part_alpha, alpha)

        alpha = xnumer/alpha

        sol_p = sol_p + alpha*pu_p
        ru_p = ru_p - alpha*apu_p
        zu_p = ru_p*gdiag

        part_error = sum(zu_p*zu_p)

        call par_sum(part_error, error)

        error = dsqrt(error)/mag

      end do

      call prim_to_all("PBC", type_system, delta_vel, sol_p, loading, mesh, exec)

      part_res_norm = sum(sol_p*sol_p)

  end select

    if (present(res_norm)) then
      call par_sum(part_res_norm, res_norm)

      res_norm = (res_norm)**0.5d0
    end if

    cg_solver_ebe = num_iters

    return

  end function cg_solver_ebe

end module conjugate_gradient_mod
