! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module solveit_evp_mod2

! Module for the elemental stiffness matrix for the elastic viscoplastic problem

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use solveit_evp_mod3
  use quadrature_mod
  use shape_3d_mod
  use stiffness_mod
  use units_mod
  use parallel_mod

  implicit none

  public

contains

  !> Construct elemental stiffness matrix for evp problem
  subroutine elt_stif_evp(mesh, crys, exec, results_prev, results, ecoos, &
                        & evel, estiff, etanstiff, eforce, incr, dtime)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys (:)
    type(exec_type), intent(in) :: exec
    type(results_type), intent(in) :: results_prev
    type(results_type), intent(inout) :: results

    real(rk), intent(in) :: ecoos(kdim, elt_sub:elt_sup)
    real(rk), intent(in) :: evel(kdim, elt_sub:elt_sup)
    real(rk), intent(out) :: estiff(kdim, kdim, elt_sub:elt_sup)
    real(rk), intent(out) :: etanstiff(kdim, kdim, elt_sub:elt_sup)
    real(rk), intent(out) :: eforce(kdim, elt_sub:elt_sup)
    integer, intent(in) :: incr
    real(rk), intent(in) :: dtime

    integer :: iqpt
    integer :: phase
    integer :: i, j, k, l
    integer :: i1
    integer :: i2
    integer :: i3
    integer :: j1
    integer :: j2
    integer :: j3
    real(rk) :: dndx(ndim, elt_sub:elt_sup, nqpt)
    real(rk) :: dndy(ndim, elt_sub:elt_sup, nqpt)
    real(rk) :: dndz(ndim, elt_sub:elt_sup, nqpt)
    real(rk) :: det(elt_sub:elt_sup, nqpt)
    real(rk) :: s11(elt_sub:elt_sup)
    real(rk) :: s12(elt_sub:elt_sup)
    real(rk) :: s13(elt_sub:elt_sup)
    real(rk) :: s21(elt_sub:elt_sup)
    real(rk) :: s22(elt_sub:elt_sup)
    real(rk) :: s23(elt_sub:elt_sup)
    real(rk) :: s31(elt_sub:elt_sup)
    real(rk) :: s32(elt_sub:elt_sup)
    real(rk) :: s33(elt_sub:elt_sup)
    real(rk) :: c(5, 5, elt_sub:elt_sup, nqpt)
    real(rk) :: c_tan(5, 5, elt_sub:elt_sup, nqpt)
    real(rk) :: f(5, elt_sub:elt_sup, nqpt)
    real(rk) :: xni(3, 5, elt_sub:elt_sup)
    real(rk) :: xnj(5, 3, elt_sub:elt_sup)
    real(rk) :: temp1(3, 5, elt_sub:elt_sup)
    real(rk) :: temp2(3, 3, elt_sub:elt_sup)
    real(rk) :: temp3(3, 5, elt_sub:elt_sup)
    real(rk) :: temp4(3, 3, elt_sub:elt_sup)
    real(rk) :: bulk_fac1(elt_sub:elt_sup)
    real(rk) :: bulk_fac2(elt_sub:elt_sup)
    real(rk) :: ftemp(3, elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    estiff = 0.0d0
    etanstiff = 0.0d0
    eforce = 0.0d0

    call material_matrix_evp(mesh, crys, exec, results_prev, results, c, &
                            & c_tan, f, dndx, dndy, dndz, evel, &
                            & incr, dtime, det, ecoos)

    do iqpt = 1, nqpt
      do phase = 1, mesh%num_phases
        where (mesh%elt_phase(elt_sub:elt_sup) .eq. phase)
          bulk_fac1 = crys(phase)%bulk_mod/results%detv(:, iqpt)* &
              & dtime
          bulk_fac2 = crys(phase)%bulk_mod/results%detv(:, iqpt)* &
              & results_prev%e_elas_kk_bar(:, iqpt)
        end where
      end do

      if (minval(det) .lt. 0.0d0) then
        do i = elt_sub, elt_sup
          if (det(i, iqpt) .lt. 0.0d0) then
            write (*, *) 'Error  :       . Element: ', i, ', &
                &determinant: ', det(i, iqpt)
          end if
        end do

        call par_quit('Error  :       . elt_stif_evp: Negative &
            &Jacobian(s)', abort=.true.)
      end if

      det(:, iqpt) = det(:, iqpt)*wtqp(1, iqpt)

      do i = 1, ndim
        i1 = 3*(i - 1) + 1
        i2 = i1 + 1
        i3 = i2 + 1

        xni(1, 1, :) = dndx(i, :, iqpt)
        xni(2, 1, :) = -dndy(i, :, iqpt)
        xni(3, 1, :) = 0.0d0
        xni(1, 2, :) = -dndx(i, :, iqpt)/3.0d0
        xni(2, 2, :) = -dndy(i, :, iqpt)/3.0d0
        xni(3, 2, :) = 2.0d0*dndz(i, :, iqpt)/3.0d0
        xni(1, 3, :) = 0.5d0*dndy(i, :, iqpt)
        xni(2, 3, :) = 0.5d0*dndx(i, :, iqpt)
        xni(3, 3, :) = 0.0d0
        xni(1, 4, :) = 0.5d0*dndz(i, :, iqpt)
        xni(2, 4, :) = 0.0d0
        xni(3, 4, :) = 0.5d0*dndx(i, :, iqpt)
        xni(1, 5, :) = 0.0d0
        xni(2, 5, :) = 0.5d0*dndz(i, :, iqpt)
        xni(3, 5, :) = 0.5d0*dndy(i, :, iqpt)

        ftemp = 0.0d0

        !  rc 6/24/2016: Reordered for better memory striding

        do l = 1, 5
          do k = 1, 3
            ftemp(k, :) = ftemp(k, :) + xni(k, l, :)*f(l, :, iqpt)
          end do
        end do

        ftemp(1, :) = ftemp(1, :) - dndx(i, :, iqpt)*bulk_fac2
        ftemp(2, :) = ftemp(2, :) - dndy(i, :, iqpt)*bulk_fac2
        ftemp(3, :) = ftemp(3, :) - dndz(i, :, iqpt)*bulk_fac2

        eforce(i1, :) = eforce(i1, :) + ftemp(1, :)*det(:, iqpt)
        eforce(i2, :) = eforce(i2, :) + ftemp(2, :)*det(:, iqpt)
        eforce(i3, :) = eforce(i3, :) + ftemp(3, :)*det(:, iqpt)

        call gen_matrix_mult(xni, c(:, :, :, iqpt), temp1)

        if (exec%itmethod .eq. "NR") call gen_matrix_mult(xni, c_tan(:, :, :, iqpt), temp3)

        do j = 1, i
          j1 = 3*(j - 1) + 1
          j2 = j1 + 1
          j3 = j2 + 1

          xnj(1, 1, :) = dndx(j, :, iqpt)
          xnj(1, 2, :) = -dndy(j, :, iqpt)
          xnj(1, 3, :) = 0.0d0
          xnj(2, 1, :) = -dndx(j, :, iqpt)/3.0d0
          xnj(2, 2, :) = -dndy(j, :, iqpt)/3.0d0
          xnj(2, 3, :) = 2.0d0*dndz(j, :, iqpt)/3.0d0
          xnj(3, 1, :) = 0.5d0*dndy(j, :, iqpt)
          xnj(3, 2, :) = 0.5d0*dndx(j, :, iqpt)
          xnj(3, 3, :) = 0.0d0
          xnj(4, 1, :) = 0.5d0*dndz(j, :, iqpt)
          xnj(4, 2, :) = 0.0d0
          xnj(4, 3, :) = 0.5d0*dndx(j, :, iqpt)
          xnj(5, 1, :) = 0.0d0
          xnj(5, 2, :) = 0.5d0*dndz(j, :, iqpt)
          xnj(5, 3, :) = 0.5d0*dndy(j, :, iqpt)

          call gen_matrix_mult(temp1, xnj, temp2)

          if (exec%itmethod .eq. "NR") call gen_matrix_mult(temp3, xnj, temp4)

          ! Assemble secant + volumentric stiffness

          s11 = temp2(1, 1, :)*det(:, iqpt) + dndx(i, :, iqpt)* &
              & dndx(j, :, iqpt)*bulk_fac1*det(:, iqpt)
          s12 = temp2(1, 2, :)*det(:, iqpt) + dndx(i, :, iqpt)* &
              & dndy(j, :, iqpt)*bulk_fac1*det(:, iqpt)
          s13 = temp2(1, 3, :)*det(:, iqpt) + dndx(i, :, iqpt)* &
              & dndz(j, :, iqpt)*bulk_fac1*det(:, iqpt)
          s21 = temp2(2, 1, :)*det(:, iqpt) + dndy(i, :, iqpt)* &
              & dndx(j, :, iqpt)*bulk_fac1*det(:, iqpt)
          s22 = temp2(2, 2, :)*det(:, iqpt) + dndy(i, :, iqpt)* &
              & dndy(j, :, iqpt)*bulk_fac1*det(:, iqpt)
          s23 = temp2(2, 3, :)*det(:, iqpt) + dndy(i, :, iqpt)* &
              & dndz(j, :, iqpt)*bulk_fac1*det(:, iqpt)
          s31 = temp2(3, 1, :)*det(:, iqpt) + dndz(i, :, iqpt)* &
              & dndx(j, :, iqpt)*bulk_fac1*det(:, iqpt)
          s32 = temp2(3, 2, :)*det(:, iqpt) + dndz(i, :, iqpt)* &
              & dndy(j, :, iqpt)*bulk_fac1*det(:, iqpt)
          s33 = temp2(3, 3, :)*det(:, iqpt) + dndz(i, :, iqpt)* &
              & dndz(j, :, iqpt)*bulk_fac1*det(:, iqpt)

          call add_to_stiffness(estiff, s11, i1, j1)
          call add_to_stiffness(estiff, s22, i2, j2)
          call add_to_stiffness(estiff, s33, i3, j3)
          call add_to_stiffness(estiff, s12, i1, j2)
          call add_to_stiffness(estiff, s13, i1, j3)
          call add_to_stiffness(estiff, s21, i2, j1)
          call add_to_stiffness(estiff, s23, i2, j3)
          call add_to_stiffness(estiff, s31, i3, j1)
          call add_to_stiffness(estiff, s32, i3, j2)

          ! Assemble tangent + volumentric stiffness

          if (exec%itmethod .eq. "NR") then
            s11 = temp4(1, 1, :)*det(:, iqpt) + dndx(i, :, iqpt)* &
                & dndx(j, :, iqpt)*bulk_fac1*det(:, iqpt)
            s12 = temp4(1, 2, :)*det(:, iqpt) + dndx(i, :, iqpt)* &
                & dndy(j, :, iqpt)*bulk_fac1*det(:, iqpt)
            s13 = temp4(1, 3, :)*det(:, iqpt) + dndx(i, :, iqpt)* &
                & dndz(j, :, iqpt)*bulk_fac1*det(:, iqpt)
            s21 = temp4(2, 1, :)*det(:, iqpt) + dndy(i, :, iqpt)* &
                & dndx(j, :, iqpt)*bulk_fac1*det(:, iqpt)
            s22 = temp4(2, 2, :)*det(:, iqpt) + dndy(i, :, iqpt)* &
                & dndy(j, :, iqpt)*bulk_fac1*det(:, iqpt)
            s23 = temp4(2, 3, :)*det(:, iqpt) + dndy(i, :, iqpt)* &
                & dndz(j, :, iqpt)*bulk_fac1*det(:, iqpt)
            s31 = temp4(3, 1, :)*det(:, iqpt) + dndz(i, :, iqpt)* &
                & dndx(j, :, iqpt)*bulk_fac1*det(:, iqpt)
            s32 = temp4(3, 2, :)*det(:, iqpt) + dndz(i, :, iqpt)* &
                & dndy(j, :, iqpt)*bulk_fac1*det(:, iqpt)
            s33 = temp4(3, 3, :)*det(:, iqpt) + dndz(i, :, iqpt)* &
                & dndz(j, :, iqpt)*bulk_fac1*det(:, iqpt)

            call add_to_stiffness(etanstiff, s11, i1, j1)
            call add_to_stiffness(etanstiff, s22, i2, j2)
            call add_to_stiffness(etanstiff, s33, i3, j3)
            call add_to_stiffness(etanstiff, s12, i1, j2)
            call add_to_stiffness(etanstiff, s13, i1, j3)
            call add_to_stiffness(etanstiff, s21, i2, j1)
            call add_to_stiffness(etanstiff, s23, i2, j3)
            call add_to_stiffness(etanstiff, s31, i3, j1)
            call add_to_stiffness(etanstiff, s32, i3, j2)
          end if ! exec%nr_use
        end do
      end do
    end do

  end subroutine elt_stif_evp

  !> Compute pressure from vel field
  subroutine recover_pressure_evp(mesh, crys, elpress, e_elas_kk_qpt)

    type(mesh_type) :: mesh
    type(crys_type) :: crys (:)
    real(rk), intent(out) :: elpress(elt_sub:elt_sup)
    real(rk), intent(in) :: e_elas_kk_qpt(elt_sub:elt_sup)

    integer :: phase

    !---------------------------------------------------------------------------

    do phase = 1, mesh%num_phases
      where (mesh%elt_phase(elt_sub:elt_sup) .eq. phase)
        elpress = -crys(phase)%bulk_mod*e_elas_kk_qpt
      end where
    end do

  end subroutine recover_pressure_evp

  !> update_state_evp: Update crystal states for entire mesh
  subroutine update_state_evp(mesh, crys, exec, results_prev, results, dtime)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys(:)
    type(exec_type), intent(in) :: exec
    type(results_type), intent(in) :: results_prev
    type(results_type), intent(inout) :: results
    real(rk), intent(in) :: dtime

    integer  :: i
    real(rk) :: dndx(ndim, elt_sub:elt_sup, nqpt)
    real(rk) :: dndy(ndim, elt_sub:elt_sup, nqpt)
    real(rk) :: dndz(ndim, elt_sub:elt_sup, nqpt)
    real(rk) :: det(elt_sub:elt_sup, nqpt)
    real(rk) :: ecoos(kdim, elt_sub:elt_sup)
    real(rk) :: evel(kdim, elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    ! Compute coordinates @(t+dt) using the vel @(t+dt)
    results%coo = results%coo + results%vel*dtime

    ! coo --> ecoos [30 x m]
    ! vel --> evel  [30 x m]
    call part_gather(ecoos, results%coo, mesh%elt_dofs, exec%dof_trace)
    call part_gather(evel, results%vel, mesh%elt_dofs, exec%dof_trace)

    ! Looping over gauss quadrature points
    do i = 1, nqpt

      call solve_evp(mesh, crys, exec, results_prev, i, results, 9999, &
                     & dtime, ecoos, evel, dndx, dndy, dndz, det)

      call compute_sliprate(mesh, crys, i, results)

      ! FIXME: doing this from vel_iteration (which runs this subroutine
      ! iteratively) will lead to wrong acmslip values
      results%acmslip(:, :, i) = results%acmslip(:, :, i) + dabs(results%sliprate(:, :, i))*dtime
    end do

  end subroutine update_state_evp

end module solveit_evp_mod2
