! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Module relative to periodic boundary conditions.
module periodicity_mod1

  use general_mod
  use types_mod
  use, intrinsic :: iso_fortran_env, only: rk => real64
  use parallel_mod
  use gather_scatter_mod
  use periodicity_mod2
  use multi_point_constraints_mod2
  use boundary_conditions_mod2

  implicit none

  public

contains

  subroutine set_periodicity(mesh, loading, loading_options, exec)

    type(mesh_type), intent(inout) :: mesh
    type(loading_options_type), intent(inout) :: loading_options
    type(loading_type), intent(inout) :: loading
    type(exec_type), intent(inout) :: exec

    ! Local:
    integer :: i  ! loop variables

    if (myid .eq. 0) then
      write (*, '(a)') "Info   :     > set periodic relations on nodes"
    end if

    ! Allocate dof arrays
    allocate (loading%ps_dof(dof_sub:dof_sup))
    allocate (loading%coeff_ps(dof_sub:dof_sup))
    allocate (loading%offset_ps(dof_sub:dof_sup))
    allocate (loading%primary_drive(dof_sub:dof_sup))
    allocate (loading%secondary_drive(dof_sub:dof_sup))
    allocate (loading%label_pv(dof_sub:dof_sup))
    allocate (loading%label_pv_drive(dof_sub:dof_sup))
    allocate (loading%label_sg(dof_sub:dof_sup))
    allocate (loading%imposed_state(dof_sub:dof_sup))
    allocate (loading%gage_length_array(dof_sub:dof_sup))
    
    allocate (loading%primary_drive_ebe(kdim, elt_sub:elt_sup))
    allocate (loading%secondary_drive_ebe(kdim, elt_sub:elt_sup))
    allocate (loading%label_pv_ebe(kdim, elt_sub:elt_sup))
    allocate (loading%label_pv_drive_ebe(kdim, elt_sub:elt_sup))
    allocate (loading%label_sg_ebe(kdim, elt_sub:elt_sup))
    allocate (loading%imposed_state_ebe(kdim, elt_sub:elt_sup))
    allocate (loading%offset_ps_ebe(kdim, elt_sub:elt_sup))

    ! Initialisation of variables
    do i = dof_sub, dof_sup
      loading%ps_dof(i) = i         ! initialisation of array to node dof index
      loading%primary_drive(i) = i
      loading%secondary_drive(i) = i
    end do

    loading%coeff_ps = 1.0d0
    loading%offset_ps = 0.0d0
    loading%offset_ps_ebe = 0.0d0
    loading%label_pv = 0
    loading%label_pv_drive = 0
    loading%label_sg = 0
    loading%label_sg_ebe = 0.0d0
    loading%imposed_state = 0
    loading%nb_secondary_dof = 0      ! at initiation, all dofs are considered as primary

    loading%pbc_status = .true.
      
    if (loading_options%imposed_strain_rate_state%nb_comp .ne. 6) then
      call par_quit("Error  : set_bc: the strain rate state is not fully imposed.")
    end if

    ! Loop on periodic relation defined in .msh
    do i = 1, mesh%num_periodicity
      ! Define the direction where mpc are applied
      select case (loading_options%periodicity)

      case ('x')
        ! check if periodic relation has to be take innto account in this case.
        if(mesh%periodicity(i)%pvect(1) .ne. 0) then
         call set_periodicity_mpc_array(mesh, loading, i) 
        end if

      case ('y')
        if(mesh%periodicity(i)%pvect(2) .ne. 0) then
         call set_periodicity_mpc_array(mesh, loading, i) 
        end if

      case ('z')
        if(mesh%periodicity(i)%pvect(3) .ne. 0) then
         call set_periodicity_mpc_array(mesh, loading, i)
        end if

      case ('xy')
        if((mesh%periodicity(i)%pvect(1) .ne. 0) .or. &
                  &(mesh%periodicity(i)%pvect(2) .ne. 0)) then
        call set_periodicity_mpc_array(mesh, loading, i) 
        end if

      case ('xz')
        if((mesh%periodicity(i)%pvect(1) .ne. 0) .or. &
                  &(mesh%periodicity(i)%pvect(3) .ne. 0)) then
         call set_periodicity_mpc_array(mesh, loading, i) 
        end if

      case ('yz')
        if((mesh%periodicity(i)%pvect(2) .ne. 0) .or. &
                  &(mesh%periodicity(i)%pvect(3) .ne. 0)) then
         call set_periodicity_mpc_array(mesh, loading, i) 
        end if

      case ('xyz', '1', 'all')
         call set_periodicity_mpc_array(mesh, loading, i) 

      end select

    end do
      
    if (myid .eq. 0) then
      write (*, '(a)') "Info   :     > set imposed strain rate state"
    end if
  
    call set_imposed_strain_rate_state(loading, loading_options)
  
    call calc_bcs_steps(loading_options, mesh, loading)
    
    call split_diagonal_periodicity_vectors(loading)
  
    call free_non_imposed_drive_nodes(loading)
    
    call calc_mpcs(mesh, exec, loading)

    call calc_driving_nodes_ebe(mesh, loading, exec)
    
    call derive_gage_length(loading, mesh, exec)

    loading%offset_ps = loading%offset_ps*loading%gage_length
    loading%step_velocity = loading%step_velocity*loading%gage_length
    
    call calc_bcs_steps(loading_options, mesh, loading)
    
    if (myid .eq. 0) then
      write (*, '(a)') "Info   :     > fix rigid body motions"
    end if
  
    call fix_rigid_body_motion(loading_options, loading, mesh)

    do i = dof_sub, dof_sup
      !loading%label_sg(i) = loading%label_pv(i)*loading%label_pv_drive(i)
      loading%label_sg(i) = loading%label_pv(i)

      if (loading%label_sg(i) .ne. 0.0d0) then
        loading%label_sg(i) = loading%label_sg(i) / abs(loading%label_sg(i))
      end if
    end do
    
    call part_gather(loading%label_sg_ebe, real(loading%label_sg,rk), mesh%elt_dofs, exec%dof_trace)

    !loading%offset_ps = loading%offset_ps*loading%label_sg
    call part_gather(loading%offset_ps_ebe, real(loading%offset_ps,rk), mesh%elt_dofs, exec%dof_trace)

    call propagate_bcs_mpcs(mesh, exec, loading)

        where (loading%bcs_vel_defined)
          loading%bcs_vel = loading%offset_ps
        end where

    if (myid .eq. 0) then
      write (*, '(a)') "Info   :     > set PBC done"
    end if
    
  end subroutine set_periodicity

 end module periodicity_mod1
