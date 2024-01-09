! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Module relative to array operations.
module array_operations_mod

  use general_mod
  use types_mod
  use, intrinsic :: iso_fortran_env, only: rk => real64
  use parallel_mod
  use gather_scatter_mod

  implicit none

  public

contains

  subroutine set_int_node_dof(array, node, node_dof, int_value)

    integer, intent(inout) :: array(dof_sub:dof_sup)
    integer, intent(in) :: node, node_dof, int_value

    if ((((node - 1)*3 + node_dof) .ge. dof_sub) .and. (((node - 1)*3 + node_dof) .le. dof_sup)) then
      array((node - 1)*3 + node_dof) = int_value
    end if
  
  end subroutine set_int_node_dof

  subroutine set_real_node_dof(array, node, node_dof, real_value)

    real(rk), intent(inout) :: array(dof_sub:dof_sup)
    integer, intent(in) :: node, node_dof
    real(rk), intent(in) :: real_value

    if ((((node - 1)*3 + node_dof) .ge. dof_sub) .and. (((node - 1)*3 + node_dof) .le. dof_sup)) then
      array((node - 1)*3 + node_dof) = real_value
    end if
  
  end subroutine set_real_node_dof

  subroutine set_logical_node_dof(array, node, node_dof, logical_value)

    logical, intent(inout) :: array(dof_sub:dof_sup)
    integer, intent(in) :: node, node_dof
    logical, intent(in) :: logical_value

    if ((((node - 1)*3 + node_dof) .ge. dof_sub) .and. (((node - 1)*3 + node_dof) .le. dof_sup)) then
      array((node - 1)*3 + node_dof) = logical_value
    end if
  
  end subroutine set_logical_node_dof
  
  subroutine add_int_node_dof(array, node, node_dof, int_value)

    integer, intent(inout) :: array(dof_sub:dof_sup)
    integer, intent(in) :: node, node_dof, int_value

    if ((((node - 1)*3 + node_dof) .ge. dof_sub) .and. (((node - 1)*3 + node_dof) .le. dof_sup)) then
      array((node - 1)*3 + node_dof) = array((node - 1)*3 + node_dof) + int_value
    end if
  
  end subroutine add_int_node_dof

  subroutine add_real_node_dof(array, node, node_dof, real_value)

    real(rk), intent(inout) :: array(dof_sub:dof_sup)
    integer, intent(in) :: node, node_dof
    real(rk), intent(in) :: real_value

    if ((((node - 1)*3 + node_dof) .ge. dof_sub) .and. (((node - 1)*3 + node_dof) .le. dof_sup)) then
      array((node - 1)*3 + node_dof) = array((node - 1)*3 + node_dof) + real_value
    end if
  
  end subroutine add_real_node_dof
  
end module array_operations_mod

