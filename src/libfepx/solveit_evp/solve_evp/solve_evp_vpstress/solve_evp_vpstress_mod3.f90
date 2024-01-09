! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module solve_evp_vpstress_mod3

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use matrix_operations_mod
  use units_mod
  use parallel_mod
  use crys_type_mod
  use crys_type_mod2

  implicit none

! Private

  public

contains

  subroutine ss_project(plocal, tensor, num_ind, indices, proj)

    ! Compute inner product of array of tensors with a fixed tensor.

    ! Arguments:
    ! proj:
    ! plocal:
    ! tensor:
    ! num_ind:
    ! indices:

    real(rk), intent(in) :: plocal(5)
    real(rk), intent(in) :: tensor(5, elt_sub:elt_sup)
    integer, intent(in) :: num_ind
    integer, intent(in) :: indices(num_ind)
    real(rk), intent(out) :: proj(elt_sub:elt_sup)

    ! Locals:

    integer :: i
    real(rk) :: proj_tmp(num_ind)

    proj_tmp = 0.0d0

    do i = 1, 5
      proj_tmp = proj_tmp + plocal(i)*tensor(i, indices)
    end do

    proj(indices) = proj_tmp

  end subroutine ss_project

  !===========================================================================

  !> Compute residual for nonlinear vp crystal stress equation.
  subroutine get_res(mesh, crys, res, rhs, rss, shear, sig, d, crss)

    type(mesh_type) :: mesh
    type(crys_type) :: crys (:)
    real(rk), intent(out) :: res(elt_sub:elt_sup)
    real(rk), intent(out) :: rhs(5, elt_sub:elt_sup)
    real(rk), intent(out) :: rss(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(out) :: shear(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: sig(5, elt_sub:elt_sup)
    real(rk), intent(in) :: d(5, elt_sub:elt_sup)
    real(rk), intent(in) :: crss(elt_sub:elt_sup)

    integer :: my_phase(elt_sub:elt_sup)
    integer :: n_slip
    integer, pointer :: indices(:) => null()
    real(rk), pointer :: p(:, :) => null()
    integer :: islip
    integer :: j
    integer :: iphase
    integer :: num_ind
    real(rk) :: xm_fake

    ! tsh, 1/26/03
    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    rhs = -d
    res = 0.0d0

    do iphase = 1, mesh%num_phases
      call crys_get(crys(iphase), dev=p)

      n_slip = crys(iphase)%numslip

      call find_indices(my_phase, iphase, indices, num_ind, elt_sub - 1)

      do islip = 1, n_slip
        call ss_project(p(:, islip), sig, num_ind, indices, rss(islip, :))

        rss(islip, indices) = rss(islip, indices)/crss(indices)

        where (abs(rss(islip, indices)) .lt. crys(iphase)%t_min)
          rss(islip, indices) = 0.0d0
        end where

        xm_fake = 0.4d0
        !xm_fake=0.02d0

        call power_law(rss(islip, :), xm_fake, &
            & crys(iphase)%gammadot_0, crys(iphase)%t_min, num_ind, indices, shear(islip, :))

        do j = 1, 5
          rhs(j, indices) = rhs(j, indices) + p(j, islip) &
                           & *shear(islip, indices)
        end do
      end do !n_slip

      deallocate (p)
      deallocate (indices)
    end do !num_phases

    do j = 1, 3
      res = res + rhs(j, :)**2.0d0
    end do

    res = sqrt(res)

  end subroutine get_res

  !===========================================================================

  !> Form single crystal stiffness matrix.
  subroutine form_crystif(mesh, crys, stif, rss, shear, crss)

    type(mesh_type) :: mesh
    type(crys_type) :: crys (:)
    real(rk), intent(out) :: stif(5, 5, elt_sub:elt_sup)
    real(rk), intent(in) :: rss(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: shear(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk), intent(in) :: crss(elt_sub:elt_sup)

    integer :: n_slip
    real(rk), pointer :: p(:, :) => null()
    integer :: my_phase(elt_sub:elt_sup)
    integer :: islip
    integer :: j
    integer :: k
    integer :: iphase
    integer :: num_ind
    integer, pointer :: indices(:) => null()
    real(rk) :: comp(elt_sub:elt_sup)
    real(rk) :: xm_fake

    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    stif = 0.0d0

    do iphase = 1, mesh%num_phases
      call crys_get(crys(iphase), dev=p)

      n_slip = crys(iphase)%numslip

      call find_indices(my_phase, iphase, indices, num_ind, elt_sub - 1)

      do islip = 1, n_slip
        xm_fake = 0.4d0

        call compliance(rss(islip, :), shear(islip, :), crss, &
            & xm_fake, crys(iphase)%t_min, num_ind, indices, comp)

        do j = 1, 5
          do k = 1, 5
            stif(j, k, indices) = stif(j, k, indices) - &
                & comp(indices)*p(j, islip)* &
                & p(k, islip)
          end do
        end do
      end do !n_slip

      deallocate (p)
      deallocate (indices)
    end do !num_phases

  end subroutine form_crystif

  !===========================================================================

  !> Power law for vp single crystal
  ! note that rss enters normalized by crss
  ! rss_min is therefore considered as a relative value
  subroutine power_law(rss, m, gammadot_0, rss_min, num_ind, indices, gammadot)

    real(rk), intent(in) :: rss(elt_sub:elt_sup)
    real(rk), intent(in) :: m
    real(rk), intent(in) :: gammadot_0
    real(rk), intent(in) :: rss_min
    integer, intent(in) :: num_ind
    integer, intent(in) :: indices(num_ind)
    real(rk), intent(out) :: gammadot(elt_sub:elt_sup)

    real(rk) :: p
    real(rk) :: gammadot_tmp(num_ind)
    real(rk) :: abs_rss(num_ind)
    real(rk) :: log_abs_rss(num_ind)
    real(rk) :: p_log_abs_rss(num_ind)

    p = 1.0d0 / m - 1.0d0
    abs_rss = abs(rss(indices))

    where (abs_rss .gt. rss_min)
      log_abs_rss = log(abs_rss)
      p_log_abs_rss = p * log_abs_rss
      gammadot_tmp = gammadot_0*rss(indices)*exp(p_log_abs_rss)

    else where
      gammadot_tmp = 0.0d0
    end where

    gammadot(indices) = gammadot_tmp

  end subroutine power_law

  !===========================================================================

  !> Form crystal compliance matrices
  subroutine compliance(t, shear, crss, xm, t_min, num_ind, indices, comp)

    real(rk), intent(in) :: t(elt_sub:elt_sup)
    real(rk), intent(in) :: shear(elt_sub:elt_sup)
    real(rk), intent(in) :: crss(elt_sub:elt_sup)
    real(rk), intent(in) :: t_min
    integer, intent(in) :: num_ind
    integer, intent(in) :: indices(num_ind)
    real(rk), intent(out) :: comp(elt_sub:elt_sup)

    real(rk) :: xm
    real(rk) :: comp_tmp(num_ind)

    comp_tmp = 0.0d0

    where (abs(t(indices)) .gt. t_min)
      comp_tmp = shear(indices)/(t(indices)*crss(indices)*xm)
    end where

    comp(indices) = comp_tmp

  end subroutine compliance

  !===========================================================================

  !> Solve an array of symmetric positive definite 5x5 systems.
  subroutine solvit(a, x)

    real(rk), intent(in) :: a(:, :, :)
    real(rk), intent(inout) :: x(:, :)
    integer :: m

    real(rk), allocatable :: a11(:)
    real(rk), allocatable :: a21(:)
    real(rk), allocatable :: a22(:)
    real(rk), allocatable :: a31(:)
    real(rk), allocatable :: a32(:)
    real(rk), allocatable :: a33(:)
    real(rk), allocatable :: a41(:)
    real(rk), allocatable :: a42(:)
    real(rk), allocatable :: a43(:)
    real(rk), allocatable :: a44(:)
    real(rk), allocatable :: a51(:)
    real(rk), allocatable :: a52(:)
    real(rk), allocatable :: a53(:)
    real(rk), allocatable :: a54(:)
    real(rk), allocatable :: a55(:)
    real(rk), allocatable :: x1(:)
    real(rk), allocatable :: x2(:)
    real(rk), allocatable :: x3(:)
    real(rk), allocatable :: x4(:)
    real(rk), allocatable :: x5(:)
    real(rk), allocatable :: v1(:)
    real(rk), allocatable :: v2(:)
    real(rk), allocatable :: v3(:)
    real(rk), allocatable :: v4(:)

    m = size(a, 3)

    allocate (a11(m))
    allocate (a21(m))
    allocate (a22(m))
    allocate (a31(m))
    allocate (a32(m))
    allocate (a33(m))
    allocate (a41(m))
    allocate (a42(m))
    allocate (a43(m))
    allocate (a44(m))
    allocate (a51(m))
    allocate (a52(m))
    allocate (a53(m))
    allocate (a54(m))
    allocate (a55(m))
    allocate (x1(m))
    allocate (x2(m))
    allocate (x3(m))
    allocate (x4(m))
    allocate (x5(m))
    allocate (v1(m))
    allocate (v2(m))
    allocate (v3(m))
    allocate (v4(m))

    a11 = a(1, 1, :)
    a21 = a(2, 1, :)
    a31 = a(3, 1, :)
    a41 = a(4, 1, :)
    a51 = a(5, 1, :)
    a22 = a(2, 2, :)
    a32 = a(3, 2, :)
    a42 = a(4, 2, :)
    a52 = a(5, 2, :)
    a33 = a(3, 3, :)
    a43 = a(4, 3, :)
    a53 = a(5, 3, :)
    a44 = a(4, 4, :)
    a54 = a(5, 4, :)
    a55 = a(5, 5, :)
    x1 = x(1, :)
    x2 = x(2, :)
    x3 = x(3, :)
    x4 = x(4, :)
    x5 = x(5, :)

    ! **  a = ldl'.
    ! **  j = 1.

    a21 = a21/a11
    a31 = a31/a11
    a41 = a41/a11
    a51 = a51/a11

    ! **  j = 2.

    v1 = a21*a11
    a22 = a22 - a21*v1
    a32 = (a32 - a31*v1)/a22
    a42 = (a42 - a41*v1)/a22
    a52 = (a52 - a51*v1)/a22

    ! **  j = 3.

    v1 = a31*a11
    v2 = a32*a22
    a33 = a33 - a31*v1 - a32*v2
    a43 = (a43 - a41*v1 - a42*v2)/a33
    a53 = (a53 - a51*v1 - a52*v2)/a33

    ! **  j = 4.

    v1 = a41*a11
    v2 = a42*a22
    v3 = a43*a33
    a44 = a44 - a41*v1 - a42*v2 - a43*v3
    a54 = (a54 - a51*v1 - a52*v2 - a53*v3)/a44

    ! **  j = 5.

    v1 = a51*a11
    v2 = a52*a22
    v3 = a53*a33
    v4 = a54*a44
    a55 = a55 - a51*v1 - a52*v2 - a53*v3 - a54*v4

    ! **  Ly=b.

    x2 = x2 - a21*x1
    x3 = x3 - a31*x1 - a32*x2
    x4 = x4 - a41*x1 - a42*x2 - a43*x3
    x5 = x5 - a51*x1 - a52*x2 - a53*x3 - a54*x4

    ! **  Dz=y.

    x1 = x1/a11
    x2 = x2/a22
    x3 = x3/a33
    x4 = x4/a44
    x5 = x5/a55

    ! **  l'x=z.

    x4 = x4 - a54*x5
    x3 = x3 - a43*x4 - a53*x5
    x2 = x2 - a32*x3 - a42*x4 - a52*x5
    x1 = x1 - a21*x2 - a31*x3 - a41*x4 - a51*x5
    x(1, :) = x1
    x(2, :) = x2
    x(3, :) = x3
    x(4, :) = x4
    x(5, :) = x5

  end subroutine solvit

  !===========================================================================

  !> Determine where diagonal elements are small.
  subroutine check_diagonals(stif, newton_ok)

    real(rk), intent(in) :: stif(5, 5, elt_sub:elt_sup)
    logical, intent(inout) :: newton_ok(elt_sub:elt_sup)

    integer :: i

    do i = 1, 5
      where (abs(stif(i, i, :)) .lt. vtiny)
        newton_ok = .false.
      end where
    end do

  end subroutine check_diagonals

  !===========================================================================

  !> Rescale the stress after solution is found
  subroutine scale_up_sigm(mesh, sig, epseff)

    type(mesh_type) :: mesh
    real(rk), intent(inout) :: sig(5, elt_sub:elt_sup)
    real(rk), intent(in) :: epseff(elt_sub:elt_sup)

    integer :: i
    integer :: num_ind
    integer :: iphase
    integer, pointer :: indices(:)
    real(rk) :: scale(elt_sub:elt_sup)
    integer :: my_phase(elt_sub:elt_sup)
    real(rk) :: xm_fake

    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    do iphase = 1, mesh%num_phases
      !-tm
      xm_fake = 0.4d0
      !xm_fake=0.02d0

      scale = epseff**xm_fake
      !scale = epseff**crystal_parm(1, iphase)

      call find_indices(my_phase, iphase, indices, num_ind, elt_sub - 1)

      do i = 1, 5
        sig(i, indices) = sig(i, indices)*scale(indices)
      end do

      deallocate (indices)
    end do !num_phases

  end subroutine scale_up_sigm

end module solve_evp_vpstress_mod3
