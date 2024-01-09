! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module solveit_isovp_mod2

! Contains subroutines:
! recover_pressure_vp: Compute pressure from vel field.

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use solveit_isovp_mod3
  use conjugate_gradient_mod
  use quadrature_mod
  use shape_3d_mod
  use surface_mod
  use units_mod
  use gather_scatter_mod
  use parallel_mod
  use stiffness_mod

  implicit none

  public

contains

  !> Compute pressure from velocity field
  subroutine recover_pressure_vp(evel, pscale, elpress, pcnst)

    real(rk) :: evel(kdim, elt_sub:elt_sup)
    real(rk) :: pscale(elt_sub:elt_sup)
    real(rk) :: elpress(elt_sub:elt_sup)
    real(rk) :: pcnst(kdim, elt_sub:elt_sup)

    integer :: i
    real(rk) :: sum(elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    sum = 0.0d0

    do i = 1, kdim
      sum = sum + pcnst(i, :)*evel(i, :)
    end do

    elpress = elpress + pscale*sum

    return

  end subroutine recover_pressure_vp

  !> Form elemental stiffness matrix for viscoplastic problem.
  subroutine elt_stif_vp(itype, gstiff, gcoos, gvel, pscale, pcnst, epseff)

    character(len = *) :: itype
    real(rk), intent(out) :: gstiff(kdim, kdim, elt_sub:elt_sup)
    real(rk) :: gcoos(kdim, elt_sub:elt_sup)
    real(rk) :: gvel(kdim, elt_sub:elt_sup)
    real(rk) :: pscale(elt_sub:elt_sup)
    real(rk) :: pcnst(kdim, elt_sub:elt_sup)
    real(rk) :: epseff(elt_sub:elt_sup)

    integer :: iqpt
    integer :: i
    integer :: j
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: j1
    integer :: j2
    integer :: j3
    real(rk) :: wt
    real(rk) :: dndx(ndim, elt_sub:elt_sup)
    real(rk) :: dndy(ndim, elt_sub:elt_sup)
    real(rk) :: dndz(ndim, elt_sub:elt_sup)
    real(rk) :: det(elt_sub:elt_sup)
    real(rk) :: s11(elt_sub:elt_sup)
    real(rk) :: s12(elt_sub:elt_sup)
    real(rk) :: s13(elt_sub:elt_sup)
    real(rk) :: s21(elt_sub:elt_sup)
    real(rk) :: s22(elt_sub:elt_sup)
    real(rk) :: s23(elt_sub:elt_sup)
    real(rk) :: s31(elt_sub:elt_sup)
    real(rk) :: s32(elt_sub:elt_sup)
    real(rk) :: s33(elt_sub:elt_sup)
    real(rk) :: t11(elt_sub:elt_sup)
    real(rk) :: mmtx(elt_sub:elt_sup)
    real(rk) :: sclfac(elt_sub:elt_sup)
    real(rk) :: c(5, 5, elt_sub:elt_sup)
    real(rk) :: xni(3, 5, elt_sub:elt_sup)
    real(rk) :: xnj(5, 3, elt_sub:elt_sup)
    real(rk) :: temp1(3, 5, elt_sub:elt_sup)
    real(rk) :: temp2(3, 3, elt_sub:elt_sup)
    real(rk) :: loc0, loc1, loc2

    !---------------------------------------------------------------------------

    gstiff = 0.0d0

    pcnst = 0.0d0
    mmtx = 0.0d0
    pscale = 0.0d0

    loc0 = 0.25d0
    loc1 = 0.25d0
    loc2 = 0.25d0

    call sfder_hpar(loc0, loc1, loc2, gcoos, dndx, dndy, dndz, det, s11, &
        & s12, s13, s21, s22, s23, s31, s32, s33)

    call material_matrix_vp(itype, c, dndx, dndy, dndz, gvel, sclfac, epseff)

    t11 = 0.0d0

    do iqpt = 1, nqpt
      loc0 = qploc(1, iqpt)
      loc1 = qploc(2, iqpt)
      loc2 = qploc(3, iqpt)

      wt = wtqp(1, iqpt)

      call sfder_hpar(loc0, loc1, loc2, gcoos, dndx, dndy, dndz, det, s11, &
          & s12, s13, s21, s22, s23, s31, s32, s33)

      det = det*wt
      pscale = pscale + sclfac*wt
      mmtx = mmtx + det

      do i = 1, ndim
        i1 = 3*(i - 1) + 1
        i2 = i1 + 1
        i3 = i2 + 1

        pcnst(i1, :) = pcnst(i1, :) - dndx(i, :)*det
        pcnst(i2, :) = pcnst(i2, :) - dndy(i, :)*det
        pcnst(i3, :) = pcnst(i3, :) - dndz(i, :)*det

        xni(1, 1, :) = dndx(i, :)
        xni(2, 1, :) = -dndy(i, :)
        xni(3, 1, :) = 0.0d0
        xni(1, 2, :) = -dndx(i, :)/3.0d0
        xni(2, 2, :) = -dndy(i, :)/3.0d0
        xni(3, 2, :) = 2.0d0*dndz(i, :)/3.0d0
        xni(1, 3, :) = 0.5d0*dndy(i, :)
        xni(2, 3, :) = 0.5d0*dndx(i, :)
        xni(3, 3, :) = 0.0d0
        xni(1, 4, :) = 0.5d0*dndz(i, :)
        xni(2, 4, :) = 0.0d0
        xni(3, 4, :) = 0.5d0*dndx(i, :)
        xni(1, 5, :) = 0.0d0
        xni(2, 5, :) = 0.5d0*dndz(i, :)
        xni(3, 5, :) = 0.5d0*dndy(i, :)

        temp1 = 0.0d0

        call gen_matrix_mult(xni, c, temp1)

        do j = 1, i
          j1 = 3*(j - 1) + 1
          j2 = j1 + 1
          j3 = j2 + 1

          xnj(1, 1, :) = dndx(j, :)
          xnj(1, 2, :) = -dndy(j, :)
          xnj(1, 3, :) = 0.0d0
          xnj(2, 1, :) = -dndx(j, :)/3.0d0
          xnj(2, 2, :) = -dndy(j, :)/3.0d0
          xnj(2, 3, :) = 2.0d0*dndz(j, :)/3.0d0
          xnj(3, 1, :) = 0.5d0*dndy(j, :)
          xnj(3, 2, :) = 0.5d0*dndx(j, :)
          xnj(3, 3, :) = 0.0d0
          xnj(4, 1, :) = 0.5d0*dndz(j, :)
          xnj(4, 2, :) = 0.0d0
          xnj(4, 3, :) = 0.5d0*dndx(j, :)
          xnj(5, 1, :) = 0.0d0
          xnj(5, 2, :) = 0.5d0*dndz(j, :)
          xnj(5, 3, :) = 0.5d0*dndy(j, :)

          temp2 = 0.0d0

          call gen_matrix_mult(temp1, xnj, temp2)

          s11 = temp2(1, 1, :)*det
          s12 = temp2(1, 2, :)*det
          s13 = temp2(1, 3, :)*det
          s21 = temp2(2, 1, :)*det
          s22 = temp2(2, 2, :)*det
          s23 = temp2(2, 3, :)*det
          s31 = temp2(3, 1, :)*det
          s32 = temp2(3, 2, :)*det
          s33 = temp2(3, 3, :)*det

          call add_to_stiffness(gstiff, s11, i1, j1)
          call add_to_stiffness(gstiff, s22, i2, j2)
          call add_to_stiffness(gstiff, s33, i3, j3)
          call add_to_stiffness(gstiff, s12, i1, j2)
          call add_to_stiffness(gstiff, s13, i1, j3)
          call add_to_stiffness(gstiff, s21, i2, j1)
          call add_to_stiffness(gstiff, s23, i2, j3)
          call add_to_stiffness(gstiff, s31, i3, j1)
          call add_to_stiffness(gstiff, s32, i3, j2)
        end do
      end do
    end do !nqpt

    pscale = pscale/mmtx
    t11 = pscale

    do i = 1, ndim
      i1 = 3*(i - 1) + 1
      i2 = i1 + 1
      i3 = i2 + 1

      do j = 1, i
        j1 = 3*(j - 1) + 1
        j2 = j1 + 1
        j3 = j2 + 1

        s11 = pcnst(i1, :)*pcnst(j1, :)*t11
        s12 = pcnst(i1, :)*pcnst(j2, :)*t11
        s13 = pcnst(i1, :)*pcnst(j3, :)*t11
        s21 = pcnst(i2, :)*pcnst(j1, :)*t11
        s22 = pcnst(i2, :)*pcnst(j2, :)*t11
        s23 = pcnst(i2, :)*pcnst(j3, :)*t11
        s31 = pcnst(i3, :)*pcnst(j1, :)*t11
        s32 = pcnst(i3, :)*pcnst(j2, :)*t11
        s33 = pcnst(i3, :)*pcnst(j3, :)*t11

        call add_to_stiffness(gstiff, s11, i1, j1)
        call add_to_stiffness(gstiff, s22, i2, j2)
        call add_to_stiffness(gstiff, s33, i3, j3)
        call add_to_stiffness(gstiff, s12, i1, j2)
        call add_to_stiffness(gstiff, s13, i1, j3)
        call add_to_stiffness(gstiff, s21, i2, j1)
        call add_to_stiffness(gstiff, s23, i2, j3)
        call add_to_stiffness(gstiff, s31, i3, j1)
        call add_to_stiffness(gstiff, s32, i3, j2)
      end do
    end do

  end subroutine elt_stif_vp

end module solveit_isovp_mod2
