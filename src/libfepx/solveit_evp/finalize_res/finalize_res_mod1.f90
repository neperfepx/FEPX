! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module finalize_res_mod

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use matrix_operations_mod
  use solve_evp_mod
  use quadrature_mod
  use shape_3d_mod
  use surface_mod
  use units_mod
  use gather_scatter_mod
  use parallel_mod

  implicit none

  public

contains

  ! oris may be results%ori or tmp oris
  subroutine finalize_res (mesh, crys, dtime, results_prev, results, load, area)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys(:)
    real(rk), intent(in) :: dtime
    type(results_type), intent(in) :: results_prev
    type(results_type), intent(inout) :: results
    real(rk), intent(out), optional :: load(mesh%num_fasets, 3)
    real(rk), intent(out), optional :: area(mesh%num_fasets)

    integer :: i, j
    real(rk) :: tmp

    call finalize_res_stressstrain(mesh, crys, results, load, area)

    call finalize_res_dphat_wphat(mesh, crys, results)

    if (allocated(results%strain_pl)) then
      results%strain_pl = results%strain_pl + results%dp_hat(:, :, cqpt)*dtime
    end if

    if (allocated(results%slip)) then
      results%slip = results%slip + results%sliprate(:, :, cqpt)*dtime
    end if

    if (allocated(results%defrate)) then
      results%defrate = results%d(:, :, :, cqpt)

      do i = elt_sub, elt_sup
        tmp = (results%velgrad(1, 1, i, cqpt) + results%velgrad(2, 2, i, cqpt) &
             + results%velgrad(3, 3, i, cqpt)) / 3.0d0
        do j = 1, 3
          results%defrate(j, j, i) = results%defrate(j, j, i) + tmp
        end do
      end do
    end if

    if (allocated(results%strain)) then
      results%strain = results%strain + results%defrate*dtime
    end if

    ! Calculate work, plastic work

    if (allocated(results%work)) then
      call finalize_res_totalwork(results_prev, dtime, results)
    end if

    if (allocated(results%work_pl)) then
      call finalize_res_plasticwork(results_prev, dtime, results)
    end if

  end subroutine finalize_res

  !> Calculate elemental stress/strain, mesh surface areas, macroscopic loads.
  subroutine finalize_res_stressstrain(mesh, crys, results, load, area)


    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys (:)
    type(results_type), intent(inout) :: results
    real(rk), intent(out) :: load(mesh%num_fasets, 3)
    real(rk), intent(out) :: area(mesh%num_fasets)

    real(rk), pointer :: e_elas_vec_dev_tmp(:, :) => null()
    real(rk) :: e_elas_vec_dev(5, elt_sub:elt_sup)
    real(rk) :: elas_t3x3(3, 3, elt_sub:elt_sup)
    real(rk) :: v_tensor(3, 3, elt_sub:elt_sup)
    real(rk) :: determ_v(elt_sub:elt_sup)
    real(rk) :: qr5x5(5, 5, elt_sub:elt_sup)
    real(rk) :: sig_sm(5, elt_sub:elt_sup)
    real(rk) :: s_avg_kk(elt_sub:elt_sup)
    real(rk), allocatable :: sig_avg_all(:)
    integer :: my_phase(elt_sub:elt_sup)
    integer :: elt_dof_min, elt_dof_max, idof, i, k
    integer :: iphase, num_ind
    integer, pointer :: indices(:) => null()
    real(rk) :: sig_vec(5, elt_sub:elt_sup)

    !---------------------------------------------------------------------------

    sig_vec = results%sig_vec(:, :, cqpt)

    my_phase = mesh%elt_phase(elt_sub:elt_sup)

    ! Deviatoric elastic lattice strain @(t+dt):
    ! sig_vec --> e_elas_vec_dev
    do iphase = 1, mesh%num_phases
      call find_indices(my_phase, iphase, indices, num_ind, elt_sub - 1)

      if (associated(e_elas_vec_dev_tmp)) then
        deallocate (e_elas_vec_dev_tmp)
      end if

      allocate (e_elas_vec_dev_tmp(5, num_ind))

      call vec_d_vec5(crys(iphase)%keinv, sig_vec(:, indices), &
          & e_elas_vec_dev_tmp, num_ind)

      e_elas_vec_dev(:, indices) = e_elas_vec_dev_tmp

      deallocate (e_elas_vec_dev_tmp)
      deallocate (indices)
    end do

    ! Deviatoric elastic strain:
    ! 5-vector --> symmetric matrix (3x3)
    ! e_elas_vec_dev --> elas_t3x3

    call vec_mat_symm(e_elas_vec_dev, elas_t3x3)

    ! Elastic strain tensor: deviatoric + volumetric

    elas_t3x3(1, 1, :) = elas_t3x3(1, 1, :) + results%e_elas_kk_bar(:, cqpt)/3.0d0
    elas_t3x3(2, 2, :) = elas_t3x3(2, 2, :) + results%e_elas_kk_bar(:, cqpt)/3.0d0
    elas_t3x3(3, 3, :) = elas_t3x3(3, 3, :) + results%e_elas_kk_bar(:, cqpt)/3.0d0

    ! Elastic strain tensor in 6-vec format

    results%elas_tot6(1, :, cqpt) = elas_t3x3(1, 1, :)
    results%elas_tot6(2, :, cqpt) = elas_t3x3(1, 2, :)
    results%elas_tot6(3, :, cqpt) = elas_t3x3(1, 3, :)
    results%elas_tot6(4, :, cqpt) = elas_t3x3(2, 2, :)
    results%elas_tot6(5, :, cqpt) = elas_t3x3(2, 3, :)
    results%elas_tot6(6, :, cqpt) = elas_t3x3(3, 3, :)

    ! Determinant of tensor v* in lattice axes.
    ! v*=i+e*

    v_tensor = elas_t3x3

    do i = 1, 3
      v_tensor(i, i, :) = elas_t3x3(i, i, :) + 1.0d0
    end do

    ! det(v*)

    call determinant(v_tensor, determ_v)

    ! Cauchy Stress at current configuration (crystal coo).
    ! Calculate deviatoric part: divide by det(v*)
    ! rc 6/24/16: Reordered loop ordering for better memory striding

    do k = elt_sub, elt_sup
      do i = 1, 5
        sig_vec(i, k) = sig_vec(i, k)/determ_v(k)
      end do
    end do

    ! Cauchy stress in sample coo
    ! oris [3x3] --> qr5x5 [5x5]

    call rotmat_symm(results%ori(:, :, :, cqpt), qr5x5, elt_sup - elt_sub + 1)

    ! deviatoric cauchy stress: sig_vec (crystal) --> sig_sm (sample)

    call mat_x_vec5(qr5x5, sig_vec, sig_sm)

    ! Average Cauchy stress (sample coo)

    ! Consider volumetric part: divide by det(v*)
    s_avg_kk = results%sig_kk(:, cqpt)/determ_v

    ! Deviatoric stress:
    ! 5-vector --> symmetric matrix (3x3)
    ! s_avg    --> stress

    call vec_mat_symm(sig_sm, results%stress)

    ! Total stress tensor: deviatoric + volumetric

    results%stress(1, 1, :) = results%stress(1, 1, :) + s_avg_kk/3.0d0
    results%stress(2, 2, :) = results%stress(2, 2, :) + s_avg_kk/3.0d0
    results%stress(3, 3, :) = results%stress(3, 3, :) + s_avg_kk/3.0d0

    ! Update surface information

    elt_dof_min = 6*elt_sub
    elt_dof_max = 6*elt_sup + 5

    allocate (sig_avg_all(elt_dof_min:elt_dof_max))

    do i = elt_sub, elt_sup
      idof = 6*i
      sig_avg_all(idof) = results%stress(1, 1, i)
      sig_avg_all(idof + 1) = results%stress(2, 1, i)
      sig_avg_all(idof + 2) = results%stress(3, 1, i)
      sig_avg_all(idof + 3) = results%stress(2, 2, i)
      sig_avg_all(idof + 4) = results%stress(3, 2, i)
      sig_avg_all(idof + 5) = results%stress(3, 3, i)
    end do

    ! Update areas and loads

    call surf_update(mesh, results%coo, sig_avg_all, load, area)

  end subroutine finalize_res_stressstrain

  subroutine finalize_res_totalwork(results_prev, dtime, results)

    ! Calculates elemental total work for a given time step.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! dtime: Time step for current increment

    type(results_type), intent(in) :: results_prev
    type(results_type), intent(inout) :: results
    real(rk), intent(in) :: dtime

    ! Locals:
    ! work_rate: Elemental total work rate
    ! work_step: Elemental total work for the current step
    ! i: Generic looping index

    real(rk) :: work_rate(elt_sub:elt_sup)
    real(rk) :: work_step(elt_sub:elt_sup)
    integer  :: i

    !---------------------------------------------------------------------------

    ! Initialize

    work_rate = 0.0d0
    work_step = 0.0d0
    results%work = 0.0d0

    do i = elt_sub, elt_sup
      ! Calculate total work (tensor inner product of cauchy stress and
      !   deformation rate tensor, here both 3x3 matrices)

      work_rate(i) = (results%defrate(1, 1, i)*results%stress(1, 1, i)) + &
          & (results%defrate(1, 2, i)*results%stress(1, 2, i)) + &
          & (results%defrate(1, 3, i)*results%stress(1, 3, i)) + &
          & (results%defrate(2, 1, i)*results%stress(2, 1, i)) + &
          & (results%defrate(2, 2, i)*results%stress(2, 2, i)) + &
          & (results%defrate(2, 3, i)*results%stress(2, 3, i)) + &
          & (results%defrate(3, 1, i)*results%stress(3, 1, i)) + &
          & (results%defrate(3, 2, i)*results%stress(3, 2, i)) + &
          & (results%defrate(3, 3, i)*results%stress(3, 3, i))

      ! Calculate work over step (trapezoidal time integration)
      !   Use previous work rate workrate

      work_step(i) = dtime*0.5d0* &
          & (work_rate(i) + (results%workrate(i)))

      ! Calculate cumulative work at current step

      results%work(i) = work_step(i) + results_prev%work(i)
    end do

    ! Update previous variables

    results%workrate = work_rate

  end subroutine finalize_res_totalwork

  subroutine finalize_res_plasticwork(results_prev, dtime, results)

    ! Calculates elemental plastic work for a given time step.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! dtime: Time step for current increment

    type(results_type), intent(in) :: results_prev
    real(rk), intent(in) :: dtime
    type(results_type), intent(inout) :: results

    ! Locals:
    ! workp_rate: Elemental plastic work rate
    ! workp_step: Elemental plastic work for the current step
    ! i: Looping indices
    ! qr5x5: 5x5 rotation matrix
    ! sigm: Mean stress
    ! sigdev: Deviatoric stress
    ! sigdev5: Deviatoric stress (5 vector)
    ! dp_hat_sam: Pl. def. rate tensor in sample basis

    real(rk) :: workp_rate(elt_sub:elt_sup)
    real(rk) :: workp_step(elt_sub:elt_sup)
    integer  :: i
    real(rk) :: qr5x5(5, 5)
    real(rk) :: sigm
    real(rk) :: sigdev(3, 3)
    real(rk) :: sigdev5(5)
    real(rk) :: dp_hat_sam(5)

    !---------------------------------------------------------------------------

    ! Initialize

    workp_rate = 0.0d0
    workp_step = 0.0d0
    results%work_pl = 0.0d0

    do i = elt_sub, elt_sup
      ! Initialize looped variables

      sigm = 0.0d0
      sigdev = 0.0d0
      sigdev5 = 0.0d0
      dp_hat_sam = 0.0d0

      ! Use transpose of crys-to-sample transformation:
      ! Lattice deform (below) transposes input (usually intented to go
      ! sample-to-crys). Transpose will let lattice_deform go crys-to-sample.

      call rotmat_symm(transpose(results%ori(:, :, i, cqpt)), qr5x5, 1)

      ! Calculate plastic work rate (tensor inner product of deviatoric cauchy
      !   stress and plastic deformation rate tensor). Since dp tensor is in 5
      !   vector form, easiest to convert 3x3 stress to deviatoric 6 vector,
      !   then to deviatoric 5 vector

      ! First, construct deviatoric stress tensor (3x3)

      sigm = (results%stress(1, 1, i) + results%stress(2, 2, i) + results%stress(3, 3, i)) &
          & /3.0d0
      sigdev = results%stress(:, :, i)
      sigdev(1, 1) = sigdev(1, 1) - sigm
      sigdev(2, 2) = sigdev(2, 2) - sigm
      sigdev(3, 3) = sigdev(3, 3) - sigm

      ! Next, construct deviatoric 5 vector for stress (sigdev5)
      !   Ordering (11-22), 33, 12, 13, 23 (with proper scalings)

      call mat_vec_symm_(sigdev, sigdev5)

      ! Rotate dp_hat to sample reference frame, find plastic work rate

      call lattice_deform_(qr5x5, results%dp_hat(:, i, cqpt), dp_hat_sam)
      workp_rate(i) = (dp_hat_sam(1)*sigdev5(1)) + &
          & (dp_hat_sam(2)*sigdev5(2)) + (dp_hat_sam(3)*sigdev5(3)) + &
          & (dp_hat_sam(4)*sigdev5(4)) + (dp_hat_sam(5)*sigdev5(5))

      ! Calculate work over step (trapezoidal time integration)
      !   Use previous work rate results%workrate_pl

      workp_step(i) = dtime*0.5d0* &
          & (workp_rate(i) + (results%workrate_pl(i)))

      ! Calculate cumulative work at current step

      results%work_pl(i) = workp_step(i) + results_prev%work_pl(i)
    end do

    results%workrate_pl = workp_rate

  end subroutine finalize_res_plasticwork

  subroutine finalize_res_dphat_wphat(mesh, crys, results)

    type(mesh_type), intent(in) :: mesh
    type(crys_type), intent(in) :: crys (:)
    type(results_type), intent(inout) :: results

    integer :: islip, i, iphase, num_id
    integer, pointer  :: indices(:) => null()
    real(rk), pointer :: p_hat_vec(:, :) => null()
    real(rk), pointer :: q_hat_vec(:, :) => null()

    results%dp_hat = 0.0d0
    results%wp_hat = 0.0d0

    do iphase = 1, mesh%num_phases
      call crys_get(crys(iphase), dev=p_hat_vec, skw=q_hat_vec)

      call find_indices(mesh%elt_phase(elt_sub:elt_sup), iphase, &
           & indices, num_id, elt_sub - 1)

      do islip = 1, crys(iphase)%numslip
        do i = 1, 5
          results%dp_hat(i, indices, cqpt) = results%dp_hat(i, indices, cqpt) &
            + results%sliprate(islip, indices, cqpt)*p_hat_vec(i, islip)
        end do

        do i = 1, 3
          results%wp_hat(i, indices, cqpt) = results%wp_hat(i, indices, cqpt) + &
                & results%sliprate(islip, indices, cqpt)*q_hat_vec(i, islip)
        end do
      end do

      deallocate (p_hat_vec, q_hat_vec)
      deallocate (indices)
    end do

  end subroutine finalize_res_dphat_wphat

end module finalize_res_mod
