! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Module relative to periodic boundary conditions.
module periodicity_mod2

  use general_mod
  use types_mod
  use, intrinsic :: iso_fortran_env, only: rk => real64
  use parallel_mod
  use gather_scatter_mod
  use periodicity_mod3
  use multi_point_constraints_mod2

  implicit none

  public

contains

  subroutine set_periodicity_mpc_array(mesh, loading, i)
    
    type(mesh_type), intent(in) :: mesh
    type(loading_type), intent(inout) :: loading
    integer, intent(in) :: i
       
    ! Local:
    integer :: node_dof

      if (loading%pv_table%flag_label(abs(mesh%periodicity(i)%pvect_label)) .eqv. .false.) then
        loading%pv_table%flag_label(abs(mesh%periodicity(i)%pvect_label)) = .true.
        loading%pv_table%label(abs(mesh%periodicity(i)%pvect_label)) = mesh%periodicity(i)%pvect_label
        loading%pv_table%primary_drive(abs(mesh%periodicity(i)%pvect_label)) = mesh%periodicity(i)%primary
        loading%pv_table%secondary_drive(abs(mesh%periodicity(i)%pvect_label)) = mesh%periodicity(i)%secondary
      end if
     
      do node_dof = 1, 3
        call set_int_node_dof(loading%primary_drive, mesh%periodicity(i)%secondary, node_dof, &
                            & ((loading%pv_table%primary_drive(abs(mesh%periodicity(i)%pvect_label))) - 1)*3 + node_dof)
      
        call set_int_node_dof(loading%secondary_drive, mesh%periodicity(i)%secondary, node_dof, &
                            & ((loading%pv_table%secondary_drive(abs(mesh%periodicity(i)%pvect_label))) - 1)*3 + node_dof)

        
        call set_int_node_dof(loading%label_pv_drive, mesh%periodicity(i)%secondary, node_dof, &
                            & loading%pv_table%label(abs(mesh%periodicity(i)%pvect_label)))
        
        call set_int_node_dof(loading%ps_dof, mesh%periodicity(i)%secondary, node_dof, &
                            & (mesh%periodicity(i)%primary - 1)*3 + node_dof)
                          
        call set_int_node_dof(loading%label_pv, mesh%periodicity(i)%secondary, node_dof, &
                            & mesh%periodicity(i)%pvect_label)
        
        call set_real_node_dof(loading%coeff_ps, mesh%periodicity(i)%secondary, node_dof, 1.0d0)
     
        call set_real_node_dof(loading%offset_ps, mesh%periodicity(i)%secondary, node_dof, 0.0d0)
      end do

  end subroutine set_periodicity_mpc_array

  subroutine split_diagonal_periodicity_vectors(loading)

    type(loading_type), intent(inout) :: loading

    integer i, j
    integer :: label_list_xp(7)
    integer :: label_list_xo(2)
    integer :: label_list_yp(5)
    integer :: label_list_yo(2)
    integer :: label_list_yn(2)
    integer :: label_list_zp(4)
    integer :: label_list_zo(2)
    integer :: label_list_zn(3)

    label_list_xp = [4, 5, 6, 7, 10, 11, 12]
    label_list_xo = [8, 9]
    label_list_yp = [4,8,9,10,11]
    label_list_yo = [6,7]
    label_list_yn = [5,12]
    label_list_zp = [6,8,10,12]
    label_list_zo = [4,5]
    label_list_zn = [7,9,11]

    do i = dof_sub, dof_sup
      
      j = mod(i,3)

        if ((j .eq. 1) .and. (any(abs(loading%label_pv(i)) .eq. label_list_xp))) then
          loading%label_pv(i) = sign(1,loading%label_pv(i))
          loading%label_pv_drive(i) = loading%pv_table%label(1)
          loading%primary_drive(i) = (loading%pv_table%primary_drive(1) - 1)*3 + 1
          loading%secondary_drive(i) = (loading%pv_table%secondary_drive(1)- 1)*3 + 1
        
        else if ((j .eq. 1) .and. (any(abs(loading%label_pv(i)) .eq. label_list_xo))) then
          loading%label_pv(i) = sign(2,loading%label_pv(i))
          loading%label_pv_drive(i) = loading%pv_table%label(2)
          loading%primary_drive(i) = (loading%pv_table%primary_drive(2)- 1)*3 + 1
          loading%secondary_drive(i) = (loading%pv_table%secondary_drive(2)- 1)*3 + 1
        
        else if ((j .eq. 2) .and. (any(abs(loading%label_pv(i)) .eq. label_list_yp))) then
          loading%label_pv(i) = sign(2,loading%label_pv(i))
          loading%label_pv_drive(i) = loading%pv_table%label(2)
          loading%primary_drive(i) = (loading%pv_table%primary_drive(2)- 1)*3 + 2
          loading%secondary_drive(i) = (loading%pv_table%secondary_drive(2)- 1)*3 + 2

        else if ((j .eq. 2) .and. (any(abs(loading%label_pv(i)) .eq. label_list_yo))) then
          loading%label_pv(i) = sign(1,loading%label_pv(i))
          loading%label_pv_drive(i) = loading%pv_table%label(1)
          loading%primary_drive(i) = (loading%pv_table%primary_drive(1)- 1)*3 + 2
          loading%secondary_drive(i) = (loading%pv_table%secondary_drive(1)- 1)*3 + 2
        
        else if ((j .eq. 2) .and. (any(abs(loading%label_pv(i)) .eq. label_list_yn))) then
          loading%label_pv(i) = -sign(2,loading%label_pv(i))
          loading%label_pv_drive(i) = loading%pv_table%label(2)
          loading%primary_drive(i) = (loading%pv_table%primary_drive(2)- 1)*3 + 2
          loading%secondary_drive(i) = (loading%pv_table%secondary_drive(2)- 1)*3 + 2
        
        else if ((j .eq. 0) .and. (any(abs(loading%label_pv(i)) .eq. label_list_zp))) then
          loading%label_pv(i) = sign(3,loading%label_pv(i))
          loading%label_pv_drive(i) = loading%pv_table%label(3)
          loading%primary_drive(i) = (loading%pv_table%primary_drive(3)- 1)*3 + 3
          loading%secondary_drive(i) = (loading%pv_table%secondary_drive(3)- 1)*3 + 3
        
        else if ((j .eq. 0) .and. (any(abs(loading%label_pv(i)) .eq. label_list_zo))) then
          loading%label_pv(i) = sign(1,loading%label_pv(i))
          loading%label_pv_drive(i) = loading%pv_table%label(1)
          loading%primary_drive(i) = (loading%pv_table%primary_drive(1)- 1)*3 + 3
          loading%secondary_drive(i) = (loading%pv_table%secondary_drive(1)- 1)*3 + 3
        
        else if ((j .eq. 0) .and. (any(abs(loading%label_pv(i)) .eq. label_list_zn))) then
          loading%label_pv(i) = -sign(3,loading%label_pv(i))
          loading%label_pv_drive(i) = loading%pv_table%label(3)
          loading%primary_drive(i) = (loading%pv_table%primary_drive(3)- 1)*3 + 3
          loading%secondary_drive(i) = (loading%pv_table%secondary_drive(3)- 1)*3 + 3

        end if

    end do

  end subroutine split_diagonal_periodicity_vectors

  subroutine set_imposed_strain_rate_state(loading, loading_options, mesh, exec)
  
    type(loading_type), intent(inout) :: loading
    type(loading_options_type), intent(inout) :: loading_options
    type(mesh_type), intent(inout) :: mesh
    type(exec_type), intent(inout) :: exec
       
    ! Local:
    real(rk) :: max_vel
    integer :: i, curr_dir, max_dir
    integer :: label_list_x(8)
    integer :: label_list_y(8)
    integer :: label_list_z(8)

    label_list_x = [1, 4, 5, 6, 7, 10, 11, 12]
    label_list_y = [2, 4, 5, 8, 9, 10, 11, 12]
    label_list_z = [3, 6, 7, 8, 9, 10, 11, 12]
    
    max_vel = 0.0d0
    max_dir = 0
    curr_dir = 0

    do i = 1, loading_options%imposed_strain_rate_state%nb_comp
      select case (loading_options%imposed_strain_rate_state%comp(i))
        case ('xx', '11')
          call apply_loading_as_offset(loading, loading_options, label_list_x, 1, i, 1)
          curr_dir = 1

        case ('yy', '22')
          call apply_loading_as_offset(loading, loading_options, label_list_y, 2, i, 2)
          curr_dir = 2

        case ('zz', '33')
          call apply_loading_as_offset(loading, loading_options, label_list_z, 3, i, 3)
          curr_dir = 3

        case ('xy', 'yx', '12', '21')
          call apply_loading_as_offset(loading, loading_options, label_list_x, 2, i, -3)
          call apply_loading_as_offset(loading, loading_options, label_list_y, 1, i, -3)
          curr_dir = 2

        case ('xz','zx', '13', '31')
          call apply_loading_as_offset(loading, loading_options, label_list_x, 3, i, -2)
          call apply_loading_as_offset(loading, loading_options, label_list_z, 1, i, -2)
          curr_dir = 1

        case ('yz','zy', '23', '32')
          call apply_loading_as_offset(loading, loading_options, label_list_y, 3, i, -1)
          call apply_loading_as_offset(loading, loading_options, label_list_z, 2, i, -1)
          curr_dir = 3

      end select
      
      if (max_vel .lt. abs (loading_options%imposed_strain_rate_state%val(i))) then
        max_vel = abs (loading_options%imposed_strain_rate_state%val(i))
        max_dir = curr_dir
      end if

    end do

    loading%step_velocity = max_vel
    loading%loading_direction = max_dir

    do i = 1, loading%num_steps
      ! FIXME disp
      if (i .eq. 1) then
        loading%vel_factor(i) = 1.0d0
      else
        loading%vel_factor(i) = loading%step_velocity(i) / loading%step_velocity(i - 1)
        if (loading%target_sign(i) .ne. loading%target_sign(i - 1)) then
          loading%vel_factor(i) = -loading%vel_factor(i)
        end if
      end if
    end do

  end subroutine set_imposed_strain_rate_state

  subroutine free_non_imposed_drive_nodes(loading)
  
    type(loading_type), intent(inout) :: loading

    integer :: i

    do i = dof_sub, dof_sup
      if ((loading%imposed_state(i) .eq. 0) .and. (i .eq. loading%secondary_drive(i)) .and. (i .ne. loading%ps_dof(i))) then
        loading%ps_dof(i) = i
        loading%primary_drive(i) = i
      end if
    end do

  end subroutine free_non_imposed_drive_nodes

  subroutine initiate_driving_matrix_parallel(mesh, loading, exec)

    type(loading_type), intent(inout) :: loading
    type(exec_type), intent(in) :: exec
    type(mesh_type), intent(in) :: mesh

    call part_gather(loading%primary_drive_ebe, real(loading%primary_drive,rk), mesh%elt_dofs, exec%dof_trace)
    call part_gather(loading%secondary_drive_ebe, real(loading%secondary_drive,rk), mesh%elt_dofs, exec%dof_trace)
    
    call part_scatter_setup(1, kdim, dof_sub, dof_sup, elt_sub, elt_sup, &
      & int(loading%primary_drive_ebe), loading%primary_drive_trace)
    call part_scatter_setup(1, kdim, dof_sub, dof_sup, elt_sub, elt_sup, &
      & int(loading%secondary_drive_ebe), loading%secondary_drive_trace)

  end subroutine initiate_driving_matrix_parallel

  subroutine calc_driving_nodes_ebe(mesh, loading, exec)

    type(loading_type), intent(inout) :: loading
    type(exec_type), intent(in) :: exec
    type(mesh_type), intent(in) :: mesh

    real(rk) :: buff_imp_st(dof_sub:dof_sup)
    
    buff_imp_st = 0.0d0
    
    call part_gather(loading%label_pv_ebe, real(loading%label_pv,rk), mesh%elt_dofs, exec%dof_trace)
    call part_gather(loading%label_pv_drive_ebe, real(loading%label_pv_drive,rk), mesh%elt_dofs, exec%dof_trace)
    call part_gather(loading%imposed_state_ebe, real(loading%imposed_state,rk), mesh%elt_dofs, exec%dof_trace)
    
    call part_gather(loading%imposed_state_ebe, real(loading%imposed_state,rk), int(loading%secondary_drive_ebe), &
          & loading%secondary_drive_trace)
    call part_gather(loading%offset_ps_ebe, real(loading%offset_ps,rk), int(loading%secondary_drive_ebe), &
          & loading%secondary_drive_trace)
    
  end subroutine calc_driving_nodes_ebe
  
  subroutine derive_gage_length(loading, mesh, exec)

    type(loading_type), intent(inout) :: loading
    type(mesh_type), intent(in) :: mesh
    type(exec_type), intent(in) :: exec

    real(rk) :: length_ps(dof_sub:dof_sup)
    real(rk) :: coo_primary_drive(dof_sub:dof_sup)
    real(rk) :: coo_secondary_drive(dof_sub:dof_sup)
    real(rk) :: g_ones(dof_sub:dof_sup)
    real(rk) :: coo_primary_drive_ebe(kdim,elt_sub:elt_sup)
    real(rk) :: coo_secondary_drive_ebe(kdim,elt_sub:elt_sup)
    real(rk) :: e_ones(kdim,elt_sub:elt_sup)
    real(rk) :: buff
    integer :: i

    length_ps = 0.0d0
    coo_primary_drive = 0.0d0
    coo_secondary_drive = 0.0d0
    coo_primary_drive_ebe = 0.0d0
    coo_secondary_drive_ebe = 0.0d0

    loading%gage_length = 0.0d0
    buff =0.0d0

    e_ones = 1.0d0
    g_ones = 0.0d0


    call part_gather(coo_primary_drive_ebe, mesh%coo, int(loading%primary_drive_ebe), loading%primary_drive_trace)
    call part_gather(coo_secondary_drive_ebe, mesh%coo, int(loading%secondary_drive_ebe), loading%secondary_drive_trace)
      
    call part_scatter(coo_primary_drive, coo_primary_drive_ebe, mesh%elt_dofs, exec%dof_trace)
    call part_scatter(coo_secondary_drive, coo_secondary_drive_ebe, mesh%elt_dofs, exec%dof_trace)
    
    call part_scatter(g_ones, e_ones, mesh%elt_dofs, exec%dof_trace)

    length_ps = abs(coo_primary_drive - coo_secondary_drive)/g_ones

    do i = dof_sub, dof_sup
        
      if (i .eq. ((loading%pv_table%secondary_drive(1) - 1)*3 + 1)) then
        loading%gage_length_alldir(1) = length_ps((loading%pv_table%secondary_drive(1) - 1)*3 + 1)
        
      else if (i .eq. ((loading%pv_table%secondary_drive(2) - 1)*3 + 2)) then
        loading%gage_length_alldir(2) = length_ps((loading%pv_table%secondary_drive(2) - 1)*3 + 2)
        
      else if (i .eq. ((loading%pv_table%secondary_drive(3) - 1)*3 + 3)) then
        loading%gage_length_alldir(3) = length_ps((loading%pv_table%secondary_drive(3) - 1)*3 + 3)
      
      end if
      
    end do

    do i = 1, 3
      call par_max(loading%gage_length_alldir(i),buff)
      loading%gage_length_alldir(i) = buff
    end do

    do i = dof_sub, dof_sup
      loading%gage_length_array(i) = loading%gage_length_alldir(mod(i,3) + 1)
    end do

    select case (loading%loading_direction)
    
    case (1)

      loading%gage_length = loading%gage_length_alldir(1)
    
    case (2)

      loading%gage_length = loading%gage_length_alldir(2)
    
    case (3)

      loading%gage_length = loading%gage_length_alldir(3)
    
    end select


  end subroutine derive_gage_length

  subroutine fix_rigid_body_motion(loading_options, loading, mesh)
    
    type(loading_type), intent(inout) :: loading
    type(loading_options_type), intent(inout) :: loading_options
    type(mesh_type), intent(in) :: mesh

    integer :: i, j, ii

    if (loading_options%imposed_strain_rate_state%nb_comp .eq. 1) then
      select case (loading_options%imposed_strain_rate_state%comp(1))
        case ('xx')
          if (loading%pv_table%flag_label(2) .eqv. .true.) then
            call apply_bc_drive_nodes(loading, 2, 3)
          else if ((loading%pv_table%flag_label(2) .eqv. .false. ) .and. (loading%pv_table%flag_label(3) .eqv. .true.)) then
            call apply_bc_drive_nodes(loading, 3, 2)

          else if (size(mesh%nsets) .ne. 0) then
            do j = 1, size(mesh%nsets)
               if ((mesh%nsets(j)%nset_label .eq. 'z0') .and. (mesh%nsets(j)%num_nset_nodes .ne. 0)) then
                do i = 1, 3
                  call set_logical_node_dof(loading%bcs_vel_defined, mesh%nsets(j)%nset_nodes(1), i, .true.)
                  call set_real_node_dof(loading%bcs_vel, mesh%nsets(j)%nset_nodes(1), i, 0.0d0)
                end do
                
                ii = 2
                do while (mesh%coo(3*(mesh%nsets(j)%nset_nodes(ii) - 1) + 2) .eq. &
                  & mesh%coo(3*(mesh%nsets(j)%nset_nodes(1) - 1) + 2))
                  ii = ii + 1
                end do
                
                call set_logical_node_dof(loading%bcs_vel_defined, mesh%nsets(j)%nset_nodes(ii), 3, .true.)
                call set_real_node_dof(loading%bcs_vel, mesh%nsets(j)%nset_nodes(ii), 3, 0.0d0)

               end if
            end do

          else
            call par_quit('Error  :     > Automatic node pinning failed: explicit definition of BC required')

          end if

        case ('yy')
          if (loading%pv_table%flag_label(1) .eqv. .true.) then
            call apply_bc_drive_nodes(loading, 1, 3)
          else if ((loading%pv_table%flag_label(1) .eqv. .false. ) .and. (loading%pv_table%flag_label(3) .eqv. .true.)) then
            call apply_bc_drive_nodes(loading, 3, 1)

          else if (size(mesh%nsets) .ne. 0) then
            do j = 1, size(mesh%nsets)
              if ((mesh%nsets(j)%nset_label .eq. 'x0') .and. (mesh%nsets(j)%num_nset_nodes .ne. 0)) then
                do i = 1, 3
                  call set_logical_node_dof(loading%bcs_vel_defined, mesh%nsets(j)%nset_nodes(1), i, .true.)
                  call set_real_node_dof(loading%bcs_vel, mesh%nsets(j)%nset_nodes(1), i, 0.0d0)
                end do
                
                ii = 2
                do while (mesh%coo(3*(mesh%nsets(j)%nset_nodes(ii) - 1) + 3) .eq. &
                  &mesh%coo(3*(mesh%nsets(j)%nset_nodes(1) - 1) + 3))
                  ii = ii + 1
                end do
                
                call set_logical_node_dof(loading%bcs_vel_defined, mesh%nsets(j)%nset_nodes(ii), 1, .true.)
                call set_real_node_dof(loading%bcs_vel, mesh%nsets(j)%nset_nodes(ii), 1, 0.0d0)

              end if
            end do

          else
                call par_quit('Error  :     > Automatic node pinning failed: explicit definition of BC required')

          end if

        case ('zz')
          if (loading%pv_table%flag_label(1) .eqv. .true.) then
            call apply_bc_drive_nodes(loading, 1, 2)
          else if ((loading%pv_table%flag_label(1) .eqv. .false. ) .and. (loading%pv_table%flag_label(2) .eqv. .true.)) then
            call apply_bc_drive_nodes(loading, 2, 1)

          else if (size(mesh%nsets) .ne. 0) then
            do j = 1, size(mesh%nsets)
               if ((mesh%nsets(j)%nset_label .eq. 'y0') .and. (mesh%nsets(j)%num_nset_nodes .ne. 0)) then
                do i = 1, 3
                  call set_logical_node_dof(loading%bcs_vel_defined, mesh%nsets(j)%nset_nodes(1), i, .true.)
                  call set_real_node_dof(loading%bcs_vel, mesh%nsets(j)%nset_nodes(1), i, 0.0d0)
                end do

                ii = 2
                do while (mesh%coo(3*(mesh%nsets(j)%nset_nodes(ii) - 1) + 1) .eq. &
                  &mesh%coo(3*(mesh%nsets(j)%nset_nodes(1) - 1) + 1))
                  ii = ii + 1
                end do

                call set_logical_node_dof(loading%bcs_vel_defined, mesh%nsets(j)%nset_nodes(ii), 2, .true.)
                call set_real_node_dof(loading%bcs_vel, mesh%nsets(j)%nset_nodes(ii), 2, 0.0d0)

               end if
            end do

          else
            call par_quit('Error  :     > Automatic node pinning failed: explicit definition of BC required')

          end if

        end select
    
    else
      
      ! arbitrarly pinning the first node
      do i = 1, 3
        call set_logical_node_dof(loading%bcs_vel_defined, 1, i, .true.)
        call set_real_node_dof(loading%bcs_vel, 1, i, 0.0d0)
      end do

    end if

  end subroutine fix_rigid_body_motion

end module periodicity_mod2
