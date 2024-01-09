! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module solve_stress_evp_mod3

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use matrix_operations_mod
  use solve_evp_vpstress_mod
  use units_mod
  use parallel_mod
  use crys_type_mod2

  implicit none

  public

contains

  subroutine matrix_fjac(mesh, fjac, keinv, w_matx, dti, dgdot, ppt, n_slip)

    type(mesh_type), intent(in) :: mesh
    integer, intent(in)   :: n_slip
    real(rk), intent(out) :: fjac(:, :, :)
    real(rk), intent(in)  :: ppt(5, 5, mesh%maxnumslip)
    real(rk), intent(in)  :: dti
    real(rk), intent(in)  :: keinv(5)
    real(rk), intent(in)  :: w_matx(:, :, :)
    real(rk), intent(in)  :: dgdot(:, :)

    integer :: islip, i, j, im, n

    !---------------------------------------------------------------------------

    n = size(fjac, 3)

    do im = 1, n
      do j = 1, 5
        do i = 1, 5
          fjac(i, j, im) = w_matx(i, j, im)*keinv(j)
        end do

        fjac(j, j, im) = fjac(j, j, im) + keinv(j)*dti
      end do

      do islip = 1, n_slip
        fjac(:, :, im) = fjac(:, :, im) + &
            & dgdot(islip, im)*ppt(:, :, islip)
      end do
    end do

  end subroutine matrix_fjac

  !===========================================================================

  subroutine wp_hat_mat5x5(wp_hat, wp_hat_matx, num_ind, indices)

    ! Add descriptions here.

    !---------------------------------------------------------------------------

    ! Arguments:

    integer, intent(in)   :: num_ind, indices(num_ind)
    real(rk), intent(in)  :: wp_hat(3, elt_sub:elt_sup)
    real(rk), intent(out) :: wp_hat_matx(5, 5, elt_sub:elt_sup)

    ! Locals:

    real(rk) :: sqr3

    !---------------------------------------------------------------------------

    sqr3 = dsqrt(3.d0)

    wp_hat_matx(:, :, indices) = 0.0d0

    wp_hat_matx(1, 3, indices) = wp_hat(1, indices)*2.0d0
    wp_hat_matx(1, 4, indices) = wp_hat(2, indices)
    wp_hat_matx(1, 5, indices) = -wp_hat(3, indices)

    wp_hat_matx(2, 4, indices) = -wp_hat(2, indices)*sqr3
    wp_hat_matx(2, 5, indices) = -wp_hat(3, indices)*sqr3

    wp_hat_matx(3, 1, indices) = -wp_hat(1, indices)*2.0d0
    wp_hat_matx(3, 4, indices) = wp_hat(3, indices)
    wp_hat_matx(3, 5, indices) = wp_hat(2, indices)

    wp_hat_matx(4, 1, indices) = -wp_hat(2, indices)
    wp_hat_matx(4, 2, indices) = wp_hat(2, indices)*sqr3
    wp_hat_matx(4, 3, indices) = -wp_hat(3, indices)
    wp_hat_matx(4, 5, indices) = wp_hat(1, indices)

    wp_hat_matx(5, 1, indices) = wp_hat(3, indices)
    wp_hat_matx(5, 2, indices) = wp_hat(3, indices)*sqr3
    wp_hat_matx(5, 3, indices) = -wp_hat(2, indices)
    wp_hat_matx(5, 4, indices) = -wp_hat(1, indices)

  end subroutine wp_hat_mat5x5

  !===========================================================================

  subroutine wp_hat_mat5x5_all(wp_hat, wp_hat_matx)

    ! Add descriptions here.

    !---------------------------------------------------------------------------

    ! Arguments:

    real(rk), intent(in)  :: wp_hat(3, elt_sub:elt_sup)
    real(rk), intent(out) :: wp_hat_matx(5, 5, elt_sub:elt_sup)

    ! Locals:

    real(rk) :: sqr3

    !---------------------------------------------------------------------------

    sqr3 = dsqrt(3.d0)

    wp_hat_matx(:, :, :) = 0.0d0

    wp_hat_matx(1, 3, :) = wp_hat(1, :)*2.0d0
    wp_hat_matx(1, 4, :) = wp_hat(2, :)
    wp_hat_matx(1, 5, :) = -wp_hat(3, :)

    wp_hat_matx(2, 4, :) = -wp_hat(2, :)*sqr3
    wp_hat_matx(2, 5, :) = -wp_hat(3, :)*sqr3

    wp_hat_matx(3, 1, :) = -wp_hat(1, :)*2.0d0
    wp_hat_matx(3, 4, :) = wp_hat(3, :)
    wp_hat_matx(3, 5, :) = wp_hat(2, :)

    wp_hat_matx(4, 1, :) = -wp_hat(2, :)
    wp_hat_matx(4, 2, :) = wp_hat(2, :)*sqr3
    wp_hat_matx(4, 3, :) = -wp_hat(3, :)
    wp_hat_matx(4, 5, :) = wp_hat(1, :)

    wp_hat_matx(5, 1, :) = wp_hat(3, :)
    wp_hat_matx(5, 2, :) = wp_hat(3, :)*sqr3
    wp_hat_matx(5, 3, :) = -wp_hat(2, :)
    wp_hat_matx(5, 4, :) = -wp_hat(1, :)

  end subroutine wp_hat_mat5x5_all

  !===========================================================================

  subroutine residual(res, rhs, d_vec_lat, e_vec, e_bar_vec, dp_hat_vec, &
      & wp_x_e, c1, num_ind, indices)

    integer, intent(in)   :: num_ind, indices(num_ind)
    real(rk), intent(out) :: rhs(5, elt_sub:elt_sup)
    real(rk), intent(out) :: res(elt_sub:elt_sup)
    real(rk), intent(in)  :: c1
    real(rk), intent(in)  :: d_vec_lat(5, elt_sub:elt_sup)
    real(rk), intent(in)  :: e_vec(5, elt_sub:elt_sup)
    real(rk), intent(in)  :: e_bar_vec(5, elt_sub:elt_sup)
    real(rk), intent(in)  :: dp_hat_vec(5, elt_sub:elt_sup)
    real(rk), intent(in)  :: wp_x_e(5, elt_sub:elt_sup)

    integer :: i

    !---------------------------------------------------------------------------

    rhs = d_vec_lat - c1*(e_vec - e_bar_vec) - dp_hat_vec - wp_x_e
    res(indices) = 0.0d0

    do i = 1, 5
      res(indices) = res(indices) + rhs(i, indices)*rhs(i, indices)
    end do

    res(indices) = dsqrt(res(indices))

  end subroutine residual

  !===========================================================================

  subroutine check_diagonals_evp(stif, newton_ok, done)

    logical, intent(out) :: newton_ok(elt_sub:elt_sup)
    logical, intent(in)  :: done(elt_sub:elt_sup)
    real(rk), intent(in) :: stif(5, 5, elt_sub:elt_sup)

    integer :: i

    !---------------------------------------------------------------------------

    do i = 1, 5
      where ((abs(stif(i, i, :)) .lt. vtiny) .and. (.not. done)) &
          & newton_ok = .false.
    end do

  end subroutine check_diagonals_evp

  !===========================================================================

  subroutine symmetrize_jac(fjac, del_s)

    real(rk), intent(inout) :: fjac(5, 5, elt_sub:elt_sup)
    real(rk), intent(inout) :: del_s(5, elt_sub:elt_sup)

    integer  :: i, j, k
    real(rk) :: tempj(5, 5, elt_sub:elt_sup)
    real(rk) :: temps(5, elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    tempj = 0.0d0
    temps = 0.0d0

    ! rc 6/24/2016: Reordered loops for better memory striding
    do j = 1, 5
      do i = 1, 5
        do k = 1, 5
          tempj(i, j, :) = tempj(i, j, :) + &
              & fjac(k, i, :)*fjac(k, j, :)
        end do
      end do
    end do

    do i = 1, 5
      do j = 1, 5
        temps(i, :) = temps(i, :) + fjac(j, i, :)*del_s(j, :)
      end do
    end do

    fjac = tempj
    del_s = temps

  end subroutine symmetrize_jac

end module solve_stress_evp_mod3
