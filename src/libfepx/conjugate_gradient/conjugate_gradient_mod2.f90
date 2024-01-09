! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module conjugate_gradient_mod2

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

  implicit none

  public

contains

  subroutine prim_to_all(type_bc, type_system, allarray, primarray, loading, mesh, exec)
    
    character(len=*), intent(in) :: type_bc
    character(len=*), intent(in) :: type_system
    type(loading_type), intent(in) :: loading
    type(mesh_type), intent(in) :: mesh
    type(exec_type), intent(in) :: exec

    real(rk), intent(inout) :: allarray(dof_sub:dof_sup)
    real(rk), intent(inout) :: primarray(dof_sub:dof_sup)

    integer :: i, j
    
    real(rk) :: primarray_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: primary_drive_primarray_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: secondary_drive_primarray_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: allarray_primary_drive(dof_sub:dof_sup)
    real(rk) :: allarray_secondary_drive(dof_sub:dof_sup)
    
    primarray_ebe = 0.0d0
    primary_drive_primarray_ebe = 0.0d0
    secondary_drive_primarray_ebe = 0.0d0
    allarray_primary_drive = 0.0d0
    allarray_secondary_drive = 0.0d0

    select case(type_bc)
   
    ! System with MPC ---------------------------------------------------------------------------------- 
    case("MPC")

    ! System with PBC ---------------------------------------------------------------------------------- 
    case("PBC")
      call part_gather(primarray_ebe, primarray, int(loading%conn_mpc), loading%mpc_trace)
      call part_gather(primary_drive_primarray_ebe, primarray, int(loading%primary_drive_ebe), loading%primary_drive_trace)
      call part_gather(secondary_drive_primarray_ebe, primarray, int(loading%secondary_drive_ebe), loading%secondary_drive_trace)
      
     where ((loading%conn_mpc .eq. loading%primary_drive_ebe) .or. (loading%imposed_state_ebe .eq. 1.0d0))
      primary_drive_primarray_ebe = 0.0d0
     end where
     
     where ((loading%conn_mpc .eq. loading%secondary_drive_ebe) .or. (loading%imposed_state_ebe .eq. 1.0d0))
      secondary_drive_primarray_ebe = 0.0d0
     end where
     
      do j = elt_sub, elt_sup
        do i = 1, kdim
          primary_drive_primarray_ebe(i,j) = primary_drive_primarray_ebe(i,j) * loading%label_sg_ebe(i,j)
          secondary_drive_primarray_ebe(i,j) = secondary_drive_primarray_ebe(i,j) * loading%label_sg_ebe(i,j)
        end do
     end do
     
     allarray = 0.0d0
     allarray_primary_drive = 0.0d0
     allarray_secondary_drive = 0.0d0
     
     call part_scatter(allarray, primarray_ebe, mesh%elt_dofs, exec%dof_trace)
     call part_scatter(allarray_primary_drive, primary_drive_primarray_ebe, mesh%elt_dofs, exec%dof_trace)
     call part_scatter(allarray_secondary_drive, secondary_drive_primarray_ebe, mesh%elt_dofs, exec%dof_trace)
    
     allarray=allarray/mesh%g_ones
     allarray_primary_drive=allarray_primary_drive/mesh%g_ones
     allarray_secondary_drive=allarray_secondary_drive/mesh%g_ones

     if (type_system .eq. "dv") then
       allarray = loading%coeff_ps*(allarray + allarray_secondary_drive - allarray_primary_drive)
     else 
       allarray = loading%coeff_ps*(allarray + allarray_secondary_drive - allarray_primary_drive) &
         &+ loading%offset_ps
     end if

    end select

  end subroutine prim_to_all
  
  ! transfer quantities to primary dofs
  subroutine all_to_prim(type_bc, primarray, allarray, loading, mesh, exec)
    
    character(len=*), intent(in) :: type_bc
    type(loading_type), intent(in) :: loading
    type(mesh_type), intent(in) :: mesh
    type(exec_type), intent(in) :: exec
    
    real(rk), intent(inout) :: allarray(dof_sub:dof_sup)
    real(rk), intent(inout) :: primarray(dof_sub:dof_sup)

    integer :: i, j
    
    real(rk) :: array_norm(dof_sub:dof_sup)
    real(rk) :: array_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: array_primary_drive_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: array_secondary_drive_ebe(kdim, elt_sub:elt_sup)
    real(rk) :: primarray_primary_drive(dof_sub:dof_sup)
    real(rk) :: primarray_secondary_drive(dof_sub:dof_sup)

    array_norm = 0.0d0
    array_ebe = 0.0d0
    array_primary_drive_ebe = 0.0d0
    array_secondary_drive_ebe = 0.0d0
    primarray_primary_drive = 0.0d0
    primarray_secondary_drive = 0.0d0
    primarray = 0.0d0
  
    select case(type_bc)
   
    ! System with MPC ---------------------------------------------------------------------------------- 
    case("MPC")

    ! System with PBC ---------------------------------------------------------------------------------- 
    case("PBC")

      where (mesh%g_ones .ne. 0)
        array_norm=allarray/mesh%g_ones
      end where

      call part_gather(array_ebe, array_norm, mesh%elt_dofs, exec%dof_trace)

      array_primary_drive_ebe = array_ebe
      array_secondary_drive_ebe = array_ebe

      where ((loading%conn_mpc .eq. loading%primary_drive_ebe) .or. (loading%imposed_state_ebe .eq. 1.0d0))
        array_primary_drive_ebe = 0.0d0
      end where
      
      where ((loading%conn_mpc .eq. loading%secondary_drive_ebe) .or. (loading%imposed_state_ebe .eq. 1.0d0))
        array_secondary_drive_ebe = 0.0d0
      end where

      do j = elt_sub, elt_sup
        do i = 1, kdim
          array_primary_drive_ebe(i,j) = array_primary_drive_ebe(i,j)* loading%label_sg_ebe(i,j)
          array_secondary_drive_ebe(i,j) = array_secondary_drive_ebe(i,j) * loading%label_sg_ebe(i,j)
        end do
      end do

      call part_scatter(primarray, array_ebe, int(loading%conn_mpc), loading%mpc_trace)
      
      call part_scatter(primarray_primary_drive, array_primary_drive_ebe, int(loading%primary_drive_ebe),&
       & loading%primary_drive_trace)
      
      call part_scatter(primarray_secondary_drive, array_secondary_drive_ebe, &
              & int(loading%secondary_drive_ebe), loading%secondary_drive_trace)

      primarray = primarray + primarray_secondary_drive - primarray_primary_drive

    end select

  end subroutine all_to_prim


end module conjugate_gradient_mod2


