! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module multi_point_constraints_mod2

! Module to calculate essential multi-point constraints.

! Contains subroutines:
! calc_mpcs: Coding dependancies between nodes dof into a specific array
! propagate_bcs_mpcs: propagate bcs through mpcs

  use general_mod
  use types_mod
  use, intrinsic :: iso_fortran_env, only: rk => real64
  use parallel_mod
  use gather_scatter_mod

  implicit none

  public

contains

  subroutine calc_mpcs(mesh, exec, loading)

    ! Defines the primary-secondary relation between node dofs:
    ! The primary and secondary dof velocity are linked by:
    ! [velocity of secondary dof] = loading%coeff_ps * [velocity of primary dof] + loading%offset_ps

    type(mesh_type), intent(inout) :: mesh
    type(exec_type), intent(in) :: exec
    type(loading_type), intent(inout) :: loading

    ! Local:
    integer :: i  ! loop variable
    real(rk) :: buff_conn_mpc(kdim,elt_sub:elt_sup) ! buffer for the convergence of conn_mpc
    real(rk) :: buff_coeff_ps_ebe(kdim,elt_sub:elt_sup)
    real(rk) :: buff_coeff_psi(dof_sub:dof_sup)
    real(rk) :: buff_coeff_psf(dof_sub:dof_sup)
    real(rk) :: buff_offset_ps_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: buff_offset_psi(dof_sub:dof_sup)
    real(rk) :: buff_offset_psf(dof_sub:dof_sup)
    real(rk) :: buff_real_dof(dof_sub:dof_sup)
    real(rk) :: e_ones(kdim, elt_sub:elt_sup)
    real(rk) :: buff_ps_dof1(dof_sub:dof_sup)
    real(rk) :: buff_ps_dof2(dof_sub:dof_sup)
    type(trace) :: buff_dof_trace

    real(rk) :: part_norm, norm

    allocate (loading%primary_dof(dof_sub:dof_sup))
    allocate (loading%conn_mpc(kdim,elt_sub:elt_sup))
    allocate (loading%count_mpc(dof_sub:dof_sup))
    allocate (loading%mask_mpc(dof_sub:dof_sup))
    allocate  (mesh%g_ones(dof_sub:dof_sup))

    e_ones = 1.0d0
    loading%count_mpc = 0.0d0
    loading%primary_dof = 0
    loading%mask_mpc = 0.0d0

    buff_coeff_ps_ebe = 1.0d0
    buff_coeff_psi = 1.0d0
    buff_coeff_psf = 1.0d0
    buff_offset_ps_ebe = 0.0d0
    buff_offset_psi = 0.0d0
    buff_offset_psf = 0.0d0

    buff_ps_dof1 = 0.0d0
    buff_ps_dof2 = 0.0d0

    buff_real_dof = 0.0d0

    ! loop to go to the top primary dof od each dof (under the form of a mpc connectivity matrix)
    buff_conn_mpc = real(mesh%elt_dofs,rk)

    part_norm = 0.0d0
    norm = 1.0d0

    call part_scatter(mesh%g_ones, e_ones, mesh%elt_dofs, exec%dof_trace)   ! multiplicity

    call part_scatter_setup(1, kdim, dof_sub, dof_sup, elt_sub, elt_sup, &
      & int(buff_conn_mpc), buff_dof_trace)

    do while (int(norm) .ne. 0)
      loading%conn_mpc = buff_conn_mpc
      buff_ps_dof1 = buff_ps_dof2

      call part_gather(buff_conn_mpc, real(loading%ps_dof,rk), int(loading%conn_mpc), buff_dof_trace)

      call part_scatter_setup(1, kdim, dof_sub, dof_sup, elt_sub, elt_sup, &
      & int(buff_conn_mpc), buff_dof_trace)

      buff_coeff_psi=loading%coeff_ps
      buff_offset_psi=loading%offset_ps
      call part_gather(buff_coeff_ps_ebe, buff_coeff_psi, int(buff_conn_mpc), buff_dof_trace)
      call part_gather(buff_offset_ps_ebe, buff_offset_psi, int(buff_conn_mpc), buff_dof_trace)

      call part_scatter(buff_coeff_psf, buff_coeff_ps_ebe, mesh%elt_dofs, exec%dof_trace)
      call part_scatter(buff_offset_psf, buff_offset_ps_ebe, mesh%elt_dofs, exec%dof_trace)
      buff_coeff_psf=buff_coeff_psf/mesh%g_ones
      buff_offset_psf=buff_offset_psf/mesh%g_ones

      loading%coeff_ps=buff_coeff_psi*buff_coeff_psf
      loading%offset_ps=buff_coeff_psi*buff_offset_psf+buff_offset_psi

      call part_scatter(buff_ps_dof2, buff_conn_mpc, mesh%elt_dofs, exec%dof_trace)

      buff_ps_dof2 = buff_ps_dof2 / mesh%g_ones

      part_norm = sum((buff_ps_dof2 - buff_ps_dof1)*(buff_ps_dof2 - buff_ps_dof1))
      call par_sum(part_norm, norm)
    end do

    loading%conn_mpc=buff_conn_mpc
    loading%mpc_trace=buff_dof_trace

    ! top primary dof array
    call part_scatter(buff_real_dof, buff_conn_mpc, mesh%elt_dofs, exec%dof_trace)
    loading%primary_dof=int(buff_real_dof/mesh%g_ones)

    ! multiplicity of the mpc connectivity matrix
    buff_real_dof = 0.0d0
    call part_scatter(buff_real_dof, e_ones, int(loading%conn_mpc), loading%mpc_trace)
    loading%count_mpc=int(buff_real_dof)

    ! mask of top primary dof (=1)
    i=1
    do i=dof_sub, dof_sup
      if (loading%count_mpc(i) .ne. 0) then
        loading%mask_mpc(i)=1
      end if
    end do

  end subroutine calc_mpcs

  subroutine propagate_bcs_mpcs(mesh, exec, loading)

    type(mesh_type), intent(in) :: mesh
    type(exec_type), intent(in) :: exec
    type(loading_type), intent(inout) :: loading

    ! Local:
    integer :: i, j
    real(rk) :: bcs_real(dof_sub:dof_sup)
    real(rk) :: buff_bcs_ebe(kdim,elt_sub:elt_sup)
    real(rk) :: buff_bcs(dof_sub:dof_sup)
    real(rk) :: buff_vel_ebe(kdim,elt_sub:elt_sup)
    real(rk) :: coeff_ps_ebe(kdim,elt_sub:elt_sup)
    real(rk) :: offset_ps_ebe(kdim,elt_sub:elt_sup)
    real(rk) :: buff_vel(dof_sub:dof_sup)
    real(rk) :: part_compare_vel
    real(rk) :: compare_vel

    ! initialization of local variables
    bcs_real = 0.0d0
    buff_bcs_ebe = 0.0d0
    buff_bcs = 0.0d0
    buff_vel_ebe = 0.0d0
    coeff_ps_ebe = 1.0d0
    offset_ps_ebe = 0.0d0
    buff_vel = 0.0d0
    part_compare_vel = 0.0d0
    compare_vel = 0.0d0

    where (loading%bcs_vel_defined)
      bcs_real = 1.0d0
    end where

    call part_gather(coeff_ps_ebe, loading%coeff_ps, mesh%elt_dofs, exec%dof_trace)
    call part_gather(offset_ps_ebe, loading%offset_ps, mesh%elt_dofs, exec%dof_trace)

    call part_gather(buff_bcs_ebe, bcs_real, mesh%elt_dofs, exec%dof_trace)
    call part_gather(buff_vel_ebe, loading%bcs_vel, mesh%elt_dofs, exec%dof_trace)

    do i=1, kdim
      do j=elt_sub, elt_sup
        if ((buff_bcs_ebe(i,j) .eq. 1.0d0) .and. (coeff_ps_ebe(i,j) .ne. 0.0d0)) then
          buff_vel_ebe(i,j) = (buff_vel_ebe(i,j) - offset_ps_ebe(i,j)) / coeff_ps_ebe(i,j)
        else if ((buff_bcs_ebe(i,j) .eq. 1.0d0) .and. (coeff_ps_ebe(i,j) .eq. 0.0d0)) then
          call par_quit('Error  :     > Conflict between MPCs and BCs')
        end if
      end do
    end do

    call part_scatter(buff_bcs, buff_bcs_ebe, int(loading%conn_mpc), loading%mpc_trace)
    call part_scatter(buff_vel, buff_vel_ebe, int(loading%conn_mpc), loading%mpc_trace)

    where (buff_bcs .ne. 0.0d0)
      buff_vel = buff_vel / buff_bcs
      buff_bcs = 1.0d0
    end where

    buff_bcs_ebe = 0.0d0
    call part_gather(buff_bcs_ebe, buff_bcs, int(loading%conn_mpc), loading%mpc_trace)
    call part_gather(buff_vel_ebe, buff_vel, int(loading%conn_mpc), loading%mpc_trace)

    where (buff_bcs_ebe .eq. 1.0d0)
      buff_vel_ebe = coeff_ps_ebe*buff_vel_ebe+offset_ps_ebe
    end where

    call part_scatter(buff_bcs, buff_bcs_ebe, mesh%elt_dofs, exec%dof_trace)
    call part_scatter(buff_vel, buff_vel_ebe, mesh%elt_dofs, exec%dof_trace)

    loading%bcs_vel_defined = .false.

    where (buff_bcs .ne. 0.0d0)
      buff_vel = buff_vel / buff_bcs
      buff_bcs = 1.0d0
      loading%bcs_vel_defined = .true.

    elsewhere
      buff_vel = buff_vel + loading%offset_ps
    end where

    if (all(buff_vel .ne. loading%bcs_vel)) then
      call par_quit('Error  :     > Conflict between MPCs and BCs')
    end if

  end subroutine propagate_bcs_mpcs

end module multi_point_constraints_mod2
