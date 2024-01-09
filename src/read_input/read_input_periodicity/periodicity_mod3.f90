! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
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

  subroutine apply_loading_as_offset_to_drive_nodes(loading, loading_options, label_list, node_dof, i)
    
    type(loading_type), intent(inout) :: loading
    type(loading_options_type), intent(inout) :: loading_options
    integer :: label_list(8)
    integer :: node_dof, i, j, ii

    do j = 1, 12
      if ((loading%pv_table%flag_label(j) .eqv. .true.) .and. (any(j .eq. abs(label_list)))) then
        do ii = dof_sub, dof_sup
          if (loading%secondary_drive(ii) .eq. ((loading%pv_table%secondary_drive(j) - 1)*3 + node_dof)) then
            call add_real_node_dof(loading%offset_ps, int((ii - node_dof)/3 + 1), node_dof, &
             & loading_options%imposed_strain_rate_state%val(i))
            call set_int_node_dof(loading%imposed_state,int((ii - node_dof)/3 + 1), node_dof, 1)
          end if
        end do
      end if
    end do

  end subroutine apply_loading_as_offset_to_drive_nodes

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
