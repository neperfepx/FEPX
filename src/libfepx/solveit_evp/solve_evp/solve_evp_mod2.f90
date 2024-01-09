! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module solve_evp_mod2

! Module hangling elastic-viscoplastic response for polycrystals.

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use crys_type_mod2
  use solve_evp_rstar_mod
  use solve_stress_evp_mod
  use solve_evp_vpstress_mod
  use matrix_operations_mod
  use units_mod
  use kinematics_mod

  implicit none

  public

contains

  !> solve_evp pre-processing
  subroutine solve_evp_pre(mesh, crys, results_prev, qpt, results, ecoos, &
                         & evel, dndx, dndy, dndz, det)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys (:)
    type(results_type), intent(in) :: results_prev
    type(results_type), intent(inout) :: results
    integer, intent(in) :: qpt
    real(rk), intent(in) :: ecoos(kdim, elt_sub:elt_sup)
    real(rk), intent(in) :: evel(kdim, elt_sub:elt_sup)
    real(rk), intent(inout) :: dndx(ndim, elt_sub:elt_sup, nqpt)
    real(rk), intent(inout) :: dndy(ndim, elt_sub:elt_sup, nqpt)
    real(rk), intent(inout) :: dndz(ndim, elt_sub:elt_sup, nqpt)
    real(rk), intent(out) :: det(elt_sub:elt_sup, nqpt)

    integer :: phase
    integer :: num_ind
    integer, pointer :: indices(:) => null()
    real(rk), pointer :: tmp2d(:, :) => null()
    real(rk) :: s11(elt_sub:elt_sup)
    real(rk) :: s12(elt_sub:elt_sup)
    real(rk) :: s13(elt_sub:elt_sup)
    real(rk) :: s21(elt_sub:elt_sup)
    real(rk) :: s22(elt_sub:elt_sup)
    real(rk) :: s23(elt_sub:elt_sup)
    real(rk) :: s31(elt_sub:elt_sup)
    real(rk) :: s32(elt_sub:elt_sup)
    real(rk) :: s33(elt_sub:elt_sup)

      do phase = 1, mesh%num_phases
        call find_indices(mesh%elt_phase(elt_sub:elt_sup), phase, indices, num_ind, elt_sub - 1)
        allocate (tmp2d(5, num_ind))

        call vec_d_vec5(crys(phase)%keinv, results_prev%sig_vec(:, indices, qpt), &
                      & tmp2d, num_ind)

        results%e_bar_vec(:, indices, qpt) = tmp2d

        deallocate (tmp2d)
        deallocate (indices)
      end do

      ! Compute quadrature quantities given a set of local coordinates
      call sfder_hpar(qploc(1, qpt), qploc(2, qpt), qploc(3, qpt), ecoos, &
          & dndx(:, :, qpt), dndy(:, :, qpt), dndz(:, :, qpt), det(:, qpt), s11, &
          & s12, s13, s21, s22, s23, s31, s32, s33)

      ! Compute d_vec, deff and w_vec via the vel gradient and its
      ! sym (d/d_kk) and skew (w) parts.
      call nodevel_velgrad(dndx(:, :, qpt), dndy(:, :, qpt), dndz(:, :, qpt), evel, &
                         & results%velgrad(:, :, :, qpt))
      call velgrad_sympart(results%velgrad(:, :, :, qpt), results%d(:, :, :, qpt), results%d_kk(:, qpt))
      call mat_vec_symm(results%d(:, :, :, qpt), results%d_vec(:, :, qpt))
      call defrate_eq(results%d(:, :, :, qpt), results%defrate_eq(:, qpt))
      call velgrad_skewpart(results%velgrad(:, :, :, qpt), results%w(:, :, :, qpt))
      call mat_vec_skew(results%w(:, :, :, qpt), results%w_vec(:, :, qpt))

  end subroutine solve_evp_pre

  subroutine solve_evp_dev (mesh, crys, exec, results_prev, qpt, results, &
                          & incr, dtime)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys (:)
    type(exec_type), intent(in) :: exec
    type(results_type), intent(in) :: results_prev
    type(results_type), intent(inout) :: results
    integer, intent(in) :: qpt
    integer, intent(in) :: incr
    real(rk), intent(in) :: dtime

    integer :: jiter_state(elt_sub:elt_sup)
    logical, parameter :: vp_log = .false.
    logical :: done(elt_sub:elt_sup)
    logical :: converged_newton
    logical :: converged_state
    integer :: iter_state
    integer :: m_elt
    integer :: islip
    integer :: n_slip
    real(rk) :: crss_0(mesh%maxnumslip, elt_sub:elt_sup)
    real(rk) :: norm_s_0(elt_sub:elt_sup)
    real(rk) :: norm_s(elt_sub:elt_sup)
    real(rk) :: diff_norm_s(elt_sub:elt_sup)
    real(rk) :: diff_crss(mesh%maxnumslip, elt_sub:elt_sup)
    integer :: iphase
    integer :: my_phase(elt_sub:elt_sup)
    logical :: converged_solution
    real(rk) :: qr3x3(3, 3, elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    my_phase = mesh%elt_phase(elt_sub:elt_sup)
    m_elt = elt_sup - elt_sub + 1

    converged_newton = .true.
    converged_state = .true.

    jiter_state = 0

    results%crss(:, :, qpt) = results_prev%crss(:, :, qpt)
    done = .false.

    ! Estimate for the stresses :
    ! For incr = 1 --> based on viscoplastic solution
    ! For incr > 1 --> based on previous solution (unrotated)

    if (incr .eq. 1) then
      ! Estimate of sig_vec (visco-plastic solution) and use a fraction
      ! as initial guess
      call solve_evp_vpstress(mesh, crys, exec, qpt, results, vp_log)

      results%sig_vec(:, :, qpt) = 0.5*results%sig_vec(:, :, qpt)
    else
      ! Use value from previous increment
      results%sig_vec(:, :, qpt) = results_prev%sig_vec(:, :, qpt)
    end if

    ! First (Forward Euler) estimate for crss and rstar (C)

    ! Estimate crss (g), rstar (r*) and d_rstar (dR*)
    ! Calculate also results%ori(:, :, :, qpt): [c]=[c_0]*[r*]
    call rotmat_skew(results%ori(:, :, :, qpt), qr3x3)
    call lattice_spin(qr3x3, results%w_vec(:, :, qpt), results%wp_hat(:, :, qpt))
    call solve_evp_rstar(mesh, crys, exec, results_prev, qpt, done, results, &
                       & dtime, 'euler_forward')

    ! Compute 2-norm for array of 5-vectors
    call vec5_norm(results%sig_vec(:, :, qpt), norm_s_0)

    ! Update crss_0
    crss_0 = results%crss(:, :, qpt)

    ! Iterate for the material state

    iter_state = 1
    converged_solution = .true.

    do while ((any(.not. done)) .and. (iter_state .le. &
        & exec%sx_max_iters_state))

      ! --> sig_vec
      call solve_evp_stress(mesh, crys, exec, qpt, done, results, &
          & dtime, converged_newton)

      if (.not. converged_newton .and. exec%auto_time .eq. 1) then
        if ((.not. converged_newton) .or. (.not. converged_state)) then
          converged_solution = .false.
        end if
      end if

      ! Calculate crss (g), rstar (r*) and d_rstar (dR*)
      ! Calculate also results%ori(:, :, :, qpt): [c]=[c_0]*[r*]
      ! no apparent reason to pass results%wp_hat(:, :, qpt)
      call solve_evp_rstar(mesh, crys, exec, results_prev, qpt, done, results, &
                         & dtime, 'euler_backward')

      call vec5_norm(results%sig_vec(:, :, qpt), norm_s)

      ! deb This section was originally dones with nested `where' constructs,
      !   but was changed because the aix compiler rejected them, although the
      !   cm compiler had no problem.

      do islip = 1, mesh%maxnumslip
        where (.not. done)
          diff_crss(islip, :) = dabs(results%crss(islip, :, qpt) - crss_0(islip, :))
        end where
      end do

      where (.not. done)
        diff_norm_s = dabs(norm_s - norm_s_0)
      end where

      do iphase = 1, mesh%num_phases
        n_slip = crys(iphase)%numslip

        do islip = 1, n_slip
          ! Currently have an initial crss_0 scaling factor from the
          !   simple anisotropic hardening model but might end up getting
          !   rid of it later versions
          where ((.not. done) .and. (my_phase .eq. iphase) &
          & .and. (diff_norm_s .lt. (exec%toler_state*crys(iphase)%g_0)) &
          & .and. (diff_crss(islip, :) .lt. (exec%toler_state*crys(iphase)%g_0)))

            done = .true.
            jiter_state = iter_state
          end where
        end do !n_slip
      end do

      do islip = 1, mesh%maxnumslip
        where (.not. done)
          norm_s_0 = norm_s
          crss_0(islip, :) = results%crss(islip, :, qpt)
        end where
      end do

      iter_state = iter_state + 1
    end do ! do while

    if (any(.not. done)) then
      converged_state = .false.
      write (*, '(a)') 'Warning:       . Not all crystals converged.'
    end if

    if (.not. converged_solution .and. exec%auto_time .eq. 1) then
      call par_quit('Error  :     > No converged solution found.')
    end if

  end subroutine solve_evp_dev

  subroutine solve_evp_vol (mesh, crys, results_prev, qpt, results, dtime)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys(:)
    type(results_type), intent(in) :: results_prev
    type(results_type), intent(inout) :: results
    integer, intent(in) :: qpt
    real(rk), intent(in) :: dtime

    integer :: m
    integer, allocatable :: my_phase(:)
    integer :: iphase

    ! somehow we have to do this instead of proceeding assumed-shape, otherwise,
    ! the code fails in debug mode
    m = elt_sup - elt_sub + 1
    allocate (my_phase(m))

    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    results%e_elas_kk_bar(:, qpt) = results_prev%e_elas_kk_bar(:, qpt) + dtime*results%d_kk(:, qpt)

    do iphase = 1, mesh%num_phases
      where (my_phase .eq. iphase)
        results%sig_kk(:, qpt) = 3.0d0*crys(iphase)%bulk_mod*results%e_elas_kk_bar(:, qpt)
      end where
    end do

  end subroutine solve_evp_vol

end module solve_evp_mod2
