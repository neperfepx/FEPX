! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Module relative to periodic boundary conditions.
module periodicity_mod3

  use general_mod
  use types_mod
  use, intrinsic :: iso_fortran_env, only: rk => real64
  use parallel_mod
  use gather_scatter_mod
  use array_operations_mod

  implicit none

  public

contains

  subroutine apply_loading_as_offset(loading, loading_options, label_list, node_dof, i, type_comp)
    
    type(loading_type), intent(inout) :: loading
    type(loading_options_type), intent(inout) :: loading_options
    integer :: label_list(8)
    integer :: node_dof, i, j, ii
    integer :: type_comp
    integer :: sign_comp
    real(rk) :: length_factor 


    do j = 1, 12

      if ((loading%pv_table%flag_label(j) .eqv. .true.) .and. (any(j .eq. abs(label_list)))) then
        do ii = dof_sub, dof_sup

          sign_comp = 1
          length_factor = 1.0d0

          if ((type_comp .eq. -2) .and. (mod(ii,3) .eq. 1)) then
            length_factor = loading%gage_length_alldir(3)/loading%gage_length_alldir(1)
            select case (j)
            case (5,7,9,11)
              sign_comp = -1
            end select
          
          else if ((type_comp .eq. -2) .and. (mod(ii,3) .eq. 0)) then
            length_factor = loading%gage_length_alldir(1)/loading%gage_length_alldir(3)
            select case (j)
            case (7,11)
              sign_comp = -1
            end select

         else if ((type_comp .eq. -3) .and. (mod(ii,3) .eq. 1)) then
            length_factor = loading%gage_length_alldir(2)/loading%gage_length_alldir(1)
            select case (j)
            case (5,12)
              sign_comp = -1
            end select
         else if ((type_comp .eq. -3) .and. (mod(ii,3) .eq. 2)) then
            length_factor = loading%gage_length_alldir(1)/loading%gage_length_alldir(2)
            select case (j)
            case (5,9,12)
              sign_comp = -1
            end select

         else if ((type_comp .eq. -1) .and. (mod(ii,3) .eq. 2)) then
            length_factor = loading%gage_length_alldir(3)/loading%gage_length_alldir(2)
            select case (j)
            case (7,9,11,12)
              sign_comp = -1
            end select
         else if ((type_comp .eq. -1) .and. (mod(ii,3) .eq. 0)) then
            length_factor = loading%gage_length_alldir(2)/loading%gage_length_alldir(3)
            select case (j)
            case (5,9,11,12)
              sign_comp = -1
            end select

         else if ((type_comp .eq. 1) .and. (mod(ii,3) .eq. 2)) then
            select case (j)
            case (7,12)
              sign_comp = -1
            end select
          
          else if ((type_comp .eq. 2) .and. (mod(ii,3) .eq. 1)) then
            select case (j)
            case (12)
              sign_comp = -1
            end select
         else if ((type_comp .eq. 2) .and. (mod(ii,3) .eq. 2)) then
            select case (j)
            case (7)
              sign_comp = -1
            end select
         
          else if ((type_comp .eq. 3) .and. (mod(ii,3) .eq. 1)) then
            select case (j)
            case (12)
              sign_comp = -1
            end select
         else if ((type_comp .eq. 3) .and. (mod(ii,3) .eq. 2)) then
            select case (j)
            case (7)
              sign_comp = -1
            end select


          end if

          if (loading%secondary_drive(ii) .eq. ((loading%pv_table%secondary_drive(j) - 1)*3 + node_dof)) then
            call add_real_node_dof(loading%offset_ps, int((ii - node_dof)/3 + 1), node_dof, &
     & loading_options%imposed_strain_rate_state%val(i)*loading%gage_length_alldir(node_dof)*sign_comp*length_factor)
            call set_int_node_dof(loading%imposed_state,int((ii - node_dof)/3 + 1), node_dof, 1)

          end if
        end do
      end if
    end do

  end subroutine apply_loading_as_offset

  subroutine apply_bc_drive_nodes(loading, label_id, pinned_secondary_dof)
    type(loading_type), intent(inout) :: loading
    integer, intent(in) :: label_id
    integer, intent(in) :: pinned_secondary_dof

    integer :: i

    do i = 1, 3
      call set_logical_node_dof(loading%bcs_vel_defined, loading%pv_table%primary_drive(label_id), i, .true.)
      call set_real_node_dof(loading%bcs_vel, loading%pv_table%primary_drive(label_id), i, 0.0d0)
    end do

    call set_logical_node_dof(loading%bcs_vel_defined, loading%pv_table%secondary_drive(label_id), &
                              & pinned_secondary_dof, .true.)
    call set_real_node_dof(loading%bcs_vel, loading%pv_table%secondary_drive(label_id), pinned_secondary_dof, 0.0d0)

  end subroutine apply_bc_drive_nodes

end module periodicity_mod3
