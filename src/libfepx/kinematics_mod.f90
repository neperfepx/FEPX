! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module kinematics_mod

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

  !> Calculate deformation rate tensor
  subroutine defrate(gvel, dndx, dndy, dndz, d)

    real(rk), intent(in) :: gvel(kdim, elt_sub:elt_sup)
    real(rk), intent(in) :: dndx(ndim, elt_sub:elt_sup)
    real(rk), intent(in) :: dndy(ndim, elt_sub:elt_sup)
    real(rk), intent(in) :: dndz(ndim, elt_sub:elt_sup)
    real(rk), intent(out) :: d(3, 3, elt_sub:elt_sup)

    integer :: i
    integer :: i1
    integer :: i2
    integer :: i3
    real(rk) :: divv(elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    d = 0.0d0

    do i = 1, ndim
      i1 = 3*(i - 1) + 1
      i2 = i1 + 1
      i3 = i2 + 1

      d(1, 1, :) = d(1, 1, :) + dndx(i, :)*gvel(i2, :)
      d(2, 2, :) = d(2, 2, :) + dndy(i, :)*gvel(i3, :)
      d(3, 3, :) = d(3, 3, :) + dndz(i, :)*gvel(i3, :)
      d(2, 1, :) = d(2, 1, :) + dndx(i, :)*gvel(i3, :) + dndy(i, :)*gvel(i2, :)
      d(3, 1, :) = d(3, 1, :) + dndx(i, :)*gvel(i3, :) + dndz(i, :)*gvel(i2, :)
      d(3, 2, :) = d(3, 2, :) + dndy(i, :)*gvel(i3, :) + dndz(i, :)*gvel(i3, :)
    end do

    d(2, 1, :) = 0.5d0*d(2, 1, :)
    d(3, 1, :) = 0.5d0*d(3, 1, :)
    d(3, 2, :) = 0.5d0*d(3, 2, :)

    d(1, 2, :) = d(2, 1, :)
    d(1, 3, :) = d(3, 1, :)
    d(2, 3, :) = d(3, 2, :)

    divv = (d(1, 1, :) + d(2, 2, :) + d(3, 3, :)) / 3.0d0

    d(1, 1, :) = d(1, 1, :) - divv
    d(2, 2, :) = d(2, 2, :) - divv
    d(3, 3, :) = d(3, 3, :) - divv

  end subroutine defrate

  !===========================================================================

  subroutine dp_wp_hat_phase(crys, mesh, qpt, dtime, results, e_elas, e_bar, sliprate, &
                           & wp_hat, num_ind, indices)

    type(crys_type), intent(in) :: crys
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: qpt
    real(rk), intent(in)  :: dtime
    type(results_type), intent(inout) :: results
    integer, intent(in) :: num_ind
    integer, intent(in) :: indices(:)
    real(rk), intent(in) :: e_elas(3, 3, elt_sub:elt_sup) ! expressed in crystal frame
    real(rk), intent(in) :: e_bar(3, 3, elt_sub:elt_sup)  ! expressed in sample frame
    real(rk), intent(in) :: sliprate(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(out) :: wp_hat(3, elt_sub:elt_sup)

    integer  :: i, islip
    real(rk) :: x(3, 3, elt_sub:elt_sup)
    real(rk) :: ee(3, 3, elt_sub:elt_sup)
    real(rk) :: dp_hat(5, elt_sub:elt_sup)
    real(rk), pointer :: p_hat_vec(:, :) => null()

    real(rk) :: p_hat(3, 3, mesh%maxnumslip)

    !---------------------------------------------------------------------------

    call crys_get(crys, dev=p_hat_vec)

    dp_hat(:, indices) = 0.0d0
    do islip = 1, crys%numslip
      do i = 1, 5
        dp_hat(i, indices) = dp_hat(i, indices) &
                         & + sliprate(islip, indices)*p_hat_vec(i, islip)
      end do
    end do

    call mat_x_mat3(e_elas, e_bar, ee)

    ! using w_vec_lat (crystal frame)
    wp_hat(1, indices) = results%w_vec_lat(1, indices, qpt) &
                     & + 0.5/dtime*(ee(2, 1, indices) - ee(1, 2, indices))
    wp_hat(2, indices) = results%w_vec_lat(2, indices, qpt) &
                     & + 0.5/dtime*(ee(3, 1, indices) - ee(1, 3, indices))
    wp_hat(3, indices) = results%w_vec_lat(3, indices, qpt) &
                     & + 0.5/dtime*(ee(3, 2, indices) - ee(2, 3, indices))

    call vec_mat_symm(p_hat_vec, p_hat)

    do islip = 1, crys%numslip
      call mat_x_mats3(e_elas(:, :, indices), p_hat(1, 1, islip), x, num_ind)

      wp_hat(1, indices) = wp_hat(1, indices) &
                       & - (x(2, 1, indices) - x(1, 2, indices)) * sliprate(islip, indices)
      wp_hat(2, indices) = wp_hat(2, indices) &
                       & - (x(3, 1, indices) - x(1, 3, indices)) * sliprate(islip, indices)
      wp_hat(3, indices) = wp_hat(3, indices) &
                       & - (x(3, 2, indices) - x(2, 3, indices)) * sliprate(islip, indices)
    end do

    results%dp_hat(:, indices, qpt) = dp_hat(:, indices)

    deallocate (p_hat_vec)

  end subroutine dp_wp_hat_phase

  !===========================================================================

  !> Compute vel gradient
  subroutine nodevel_velgrad(dndx, dndy, dndz, gvel, velgrad)

    real(rk), intent(in) :: dndx(ndim, elt_sub:elt_sup)
    real(rk), intent(in) :: dndy(ndim, elt_sub:elt_sup)
    real(rk), intent(in) :: dndz(ndim, elt_sub:elt_sup)
    real(rk), intent(in) :: gvel(kdim, elt_sub:elt_sup)
    real(rk), intent(out) :: velgrad(3, 3, elt_sub:elt_sup)

    integer :: i
    integer :: i1
    integer :: i2
    integer :: i3

    !---------------------------------------------------------------------------

    velgrad = 0.0d0

    do i = 1, ndim
      i1 = 3*(i - 1) + 1
      i2 = i1 + 1
      i3 = i2 + 1

      velgrad(1, 1, :) = velgrad(1, 1, :) + dndx(i, :)*gvel(i1, :)
      velgrad(1, 2, :) = velgrad(1, 2, :) + dndy(i, :)*gvel(i1, :)
      velgrad(1, 3, :) = velgrad(1, 3, :) + dndz(i, :)*gvel(i1, :)
      velgrad(2, 1, :) = velgrad(2, 1, :) + dndx(i, :)*gvel(i2, :)
      velgrad(2, 2, :) = velgrad(2, 2, :) + dndy(i, :)*gvel(i2, :)
      velgrad(2, 3, :) = velgrad(2, 3, :) + dndz(i, :)*gvel(i2, :)
      velgrad(3, 1, :) = velgrad(3, 1, :) + dndx(i, :)*gvel(i3, :)
      velgrad(3, 2, :) = velgrad(3, 2, :) + dndy(i, :)*gvel(i3, :)
      velgrad(3, 3, :) = velgrad(3, 3, :) + dndz(i, :)*gvel(i3, :)
    end do

  end subroutine nodevel_velgrad

end module kinematics_mod
