! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module kinematics_mod_bis

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use matrix_operations_mod
  use types_mod
  use crys_type_mod
  use crys_type_mod2

  implicit none

  public

contains

subroutine find_wp_hat(crys, mesh, qpt, dtime, results, e_elas, e_bar, sliprate, &
      & wp_hat)

    type(crys_type), intent(in) :: crys (:)
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: qpt
    real(rk), intent(in) :: dtime
    type(results_type), intent(inout) :: results
    real(rk), intent(in) :: e_elas(3, 3, elt_sub:elt_sup)
    real(rk), intent(in) :: e_bar(3, 3, elt_sub:elt_sup)
    real(rk), intent(in) :: sliprate(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(out) :: wp_hat(3, elt_sub:elt_sup)
    integer :: iphase, num_ind
    integer, pointer :: indices(:) => null()

    do iphase = 1, mesh%num_phases
      call find_indices(mesh%elt_phase(elt_sub:elt_sup), iphase, indices, num_ind, elt_sub - 1)

      call find_wp_hat_phase(crys(iphase), mesh, qpt, dtime, results, e_elas, &
                           & e_bar, sliprate, wp_hat, indices)

      deallocate (indices)
    end do

  end subroutine find_wp_hat

subroutine find_wp_hat_phase(crys, mesh, qpt, dtime, results, e_elas, e_bar, sliprate, &
                           & wp_hat, indices)

    type(crys_type), intent(in) :: crys
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: qpt
    real(rk), intent(in) :: dtime
    type(results_type), intent(inout) :: results
    integer, intent(in) :: indices(:)
    real(rk), intent(in) :: e_elas(3, 3, elt_sub:elt_sup)
    real(rk), intent(in) :: e_bar(3, 3, elt_sub:elt_sup)
    real(rk), intent(in) :: sliprate(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(out) :: wp_hat(3, elt_sub:elt_sup)

    integer :: i, islip
    real(rk) :: x(3, 3, elt_sub:elt_sup)
    real(rk) :: ee(3, 3, elt_sub:elt_sup)
    real(rk) :: dp_hat(5, elt_sub:elt_sup)
    real(rk), pointer :: p_hat_vec(:, :) => null()
    real(rk) :: qr5x5(5, 5, elt_sub:elt_sup)

    real(rk) :: dp_hat_tens(3, 3, elt_sub:elt_sup)
    real(rk) :: temp(5, elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    call crys_get(crys, dev=p_hat_vec)

    dp_hat(:, indices) = 0.0d0
    do islip = 1, crys%numslip
      do i = 1, 5
        dp_hat(i, indices) = dp_hat(i, indices) &
                           + sliprate(islip, indices) * p_hat_vec(i, islip)
      end do
    end do

    call mat_x_mat3(e_elas, e_bar, ee)

    ! using w_vec (reference frame)
    wp_hat(1, indices) = results%w_vec(1, indices, qpt) &
                     & + 0.5/dtime*(ee(2, 1, indices) - ee(1, 2, indices))
    wp_hat(2, indices) = results%w_vec(2, indices, qpt) &
                     & + 0.5/dtime*(ee(3, 1, indices) - ee(1, 3, indices))
    wp_hat(3, indices) = results%w_vec(3, indices, qpt) &
                     & + 0.5/dtime*(ee(3, 2, indices) - ee(2, 3, indices))

    call rotmat_symm(results%ori(:, :, :, qpt), qr5x5, elt_sup - elt_sub + 1)
    call mat_x_vec5(qr5x5, dp_hat, temp)
    call vec_mat_symm(temp, dp_hat_tens)
    call mat_x_mat3(e_elas, dp_hat_tens, x)

    wp_hat(1, indices) = wp_hat(1, indices) &
                     & - (x(2, 1, indices) - x(1, 2, indices))
    wp_hat(2, indices) = wp_hat(2, indices) &
                     & - (x(3, 1, indices) - x(1, 3, indices))
    wp_hat(3, indices) = wp_hat(3, indices) &
                     & - (x(3, 2, indices) - x(2, 3, indices))

    deallocate (p_hat_vec)

  end subroutine find_wp_hat_phase

end module kinematics_mod_bis
