! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module write_res_mod

! Module for printing variables to file at the end of a step.

! Contains subroutines:

! General printing and handling of field variable data:
! write_res: Print variables on a given load step

! Handling the .sim, forces/X, and post.conv files
! write_force_file_headers: Writes the file headers for table formatting
! write_force_file_data: Writes the surface forces [x y z] for all surfaces
! write_conv_file_headers: Writes the file headers for the table formatting
! write_conv_file_data: Writes the various convergence statistics

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use matrix_operations_mod
  use orientation_conversion_mod
  use units_mod
  use parallel_mod
  use write_res_mod2
  use kinematics_mod
  use utils_mod

  implicit none

! Public

  public

contains

  subroutine write_res(step, mesh, crys, results, dtime, printing)

    ! Print variables on a given load step

    !---------------------------------------------------------------------------

    ! Arguments:
    ! step: Step number

    integer, intent(in) :: step
    type(mesh_type), intent(in) :: mesh
    type(crys_type) :: crys(:)
    type(results_type), intent(in) :: results
    real(rk), intent(in) :: dtime
    type(printing_type), intent(in) :: printing

    include 'mpif-handles.h'

    ! Locals:
    ! io: File control integer
    ! i, j, igrain: Looping indices
    ! p: Passive (1) or active (-1) orientation convention
    ! aa: Angle-axis arrays
    ! eb: Euler-Bunge angles
    ! ek: Euler-Kocks angles
    ! rod: Rodrigues' vectors
    ! quat: Quaternions
    ! elas_sam: Elastic strain in the sample basis
    ! stress_eq: Equivalent (von Mises) stress
    integer :: io
    integer :: i, j, phase, numdim
    real(rk) :: tmp_oris(3, 3, elt_sub:elt_sup)
    real(rk) :: aa(4, mesh%num_elts)
    real(rk) :: eb(3, mesh%num_elts)
    real(rk) :: ek(3, mesh%num_elts)
    real(rk) :: rod(3, mesh%num_elts)
    real(rk) :: quat(4, mesh%num_elts)
    real(rk) :: elas_sam(6, elt_sub:elt_sup)
    real(rk) :: plstrain6(6, elt_sub:elt_sup)
    real(rk) :: plstrain6_sam(6, elt_sub:elt_sup)
    real(rk) :: stress_eq(elt_sub:elt_sup)
    real(rk) :: dp_hat_mat(6, elt_sub:elt_sup)
    real(rk) :: dp_hat_sam(6)
    real(rk) :: wp_hat_sam(3, elt_sub:elt_sup)
    real(rk) :: qr5x5(5, 5)
    real(rk) :: qr3x3(3, 3)
    real(rk) :: elt_val(elt_sub:elt_sup)
    real(rk) :: tmp_mat3x3(3, 3, elt_sub:elt_sup)
    real(rk) :: tmp_vec3(3, elt_sub:elt_sup)
    real(rk) :: tmp_vec4(4, elt_sub:elt_sup)

    real(rk), allocatable  :: result_buffer(:), result_buffer_arr(:, :)
    real(rk), allocatable  :: result_buffer_tensor(:, :, :)
    character(len=10000):: filename, dir_name, cwd
    character(len=8):: step_num
    character(len=8) :: charid ! assumes less than 10,000 processes
    real(rk), parameter :: pi_over_180 = 4.0d0*datan(1.0d0)/180.0d0
    character(len=50) :: fmt

    real(rk) :: eps = 1e-15

    !---------------------------------------------------------------------------

    ! Allocate results string arrays with "headers"
    call getcwd(cwd)
    dir_name = trim(cwd)//"/simulation.sim"
    write (step_num, '(a,i0)') "step", step
    write (charid, '(i0)') myid + 1
    ! Begin writing output

    ! Nodal values

    ! Nodal coordinates
    if (printing%print_coo) then
      allocate (result_buffer(1:mesh%num_nodes*3))
      call par_gatherv(results%coo(dof_sub:dof_sup), result_buffer, 3*mesh%global_info(2, :), mesh%global_info(3, :))
      filename = trim(dir_name)//'/results/nodes/coo/coo.'//trim(step_num)
      if (myid == 0) then
        open (printing%coo_u, file=filename)
        call write_res_vector(result_buffer(:), printing%coo_u, 3, step, 1, mesh%num_nodes*3)
      end if
      deallocate (result_buffer)
    end if

    ! Nodal displacements
    if (printing%print_disp) then
      allocate (result_buffer(1:mesh%num_nodes*3))
      call par_gatherv((results%coo - mesh%coo), result_buffer, 3*mesh%global_info(2, :), mesh%global_info(3, :))
      filename = trim(dir_name)//'/results/nodes/disp/disp.'//trim(step_num)
      if (myid == 0) then
        open (printing%disp_u, file=filename)
        call write_res_vector(result_buffer(:), printing%disp_u, 3, step, 1, mesh%num_nodes*3)
      end if
      deallocate (result_buffer)
    end if

    ! Nodal velocities
    if (printing%print_vel) then
      allocate (result_buffer(1:mesh%num_nodes*3))
      call par_gatherv(results%vel, result_buffer, 3*mesh%global_info(2, :), mesh%global_info(3, :))
      filename = trim(dir_name)//'/results/nodes/vel/vel.'//trim(step_num)
      if (myid == 0) then
        open (printing%vel_u, file=filename)
        call write_res_vector(result_buffer(:), printing%vel_u, 3, step, 1, mesh%num_nodes*3)
      end if
      deallocate (result_buffer)
    end if

    ! Elemental values

    ! Orientations
    if (printing%print_ori) then
      !
      io = printing%ori_u
      filename = trim(dir_name)//'/results/elts/ori/ori.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      !
      ! If passive required in output, transposing
      if (mesh%orientation_convention .eq. 'active') then
        tmp_oris = results%ori(:, :, :, cqpt)
      else if (mesh%orientation_convention .eq. 'passive') then
        do j = elt_sub, elt_sup
          tmp_oris(:, :, j) = transpose(results%ori(:, :, j, cqpt))
        end do
      end if
      !
      numdim = 9
      allocate (result_buffer_tensor(3, 3, mesh%num_elts))
      call par_gatherv_tensor(tmp_oris, result_buffer_tensor, numdim*mesh%global_info(1, :),&
        &numdim*mesh%global_info(4, :))
      !
      if (myid .eq. 0) then
        ! Determine parameterization from orientation options
        if (mesh%orientation_parameterization .eq. 'axis-angle') then
          call rotmat_to_axis_angle(result_buffer_tensor, aa)
          call write_res_array(aa, io, 4, step, 1, mesh%num_elts)

        else if (mesh%orientation_parameterization .eq. 'euler-bunge') then
          call rotmat_to_euler_bunge(result_buffer_tensor, eb)
          call write_res_array(eb, io, 3, step, 1, mesh%num_elts)

        else if (mesh%orientation_parameterization .eq. 'euler-kocks') then
          call rotmat_to_euler_kocks(result_buffer_tensor, ek)
          call write_res_array(ek, io, 3, step, 1, mesh%num_elts)

        else if (mesh%orientation_parameterization .eq. 'rodrigues') then
          call rotmat_to_rodrigues(result_buffer_tensor, rod)
          call write_res_array(rod, io, 3, step, 1, mesh%num_elts)

        else if (mesh%orientation_parameterization .eq. 'quaternion') then
          call rotmat_to_quat(result_buffer_tensor, quat)
          call write_res_array(quat, io, 4, step, 1, mesh%num_elts)
        end if
      end if
      deallocate (result_buffer_tensor)

    end if

    !!! Need to figure out how to handle phase
    ! crss Values
    if (printing%print_crss) then
      io = printing%crss_u
      filename = trim(dir_name)//'/results/elts/crss/crss.'//trim(step_num)
      if (myid == 0) open (io, file=filename)

      numdim = mesh%maxnumslip
      allocate (result_buffer_arr(numdim, mesh%num_elts))
      call par_gatherv_array(results%crss(:, :, cqpt), result_buffer_arr, numdim*mesh%global_info(1, :), &
                             numdim*mesh%global_info(4, :))

      if (myid .eq. 0) then
        do j = 1, mesh%num_elts
          phase = mesh%elt_phase(j)
          ! Check if the element is fcc
          if (crys(phase)%structure .eq. "fcc") then
            ! Next, check if the hardening is isotropic or anisotropic.
            write (io, '(12(e13.7,1x))') (result_buffer_arr(:, j))
            !
            ! Else, check if the element is bcc
          else if (crys(phase)%structure .eq. "bcc") then
            ! Next, check if 112 slip is considered
            if (crys(phase)%g_0_112 .eqv. .false.) then ! Not considered
              ! Next, check if the hardening is isotropic or anisotropic.
              write (io, '(12(e13.7,1x))') (result_buffer_arr(:, j))
            else if (crys(phase)%g_0_112 .eqv. .true.) then ! Considered
              ! Next, check if the hardening is isotropic or anisotropic.
              write (io, '(24(e13.7,1x))') (result_buffer_arr(:, j))
            end if
            ! Else, check if the element is hcp.
          else if (crys(phase)%structure .eq. "hcp") then
            ! Next, check if the hardening is isotropic or anisotropic.
            write (io, '(18(e13.7,1x))') (result_buffer_arr(:, j))
            ! Else, check if the element is bct.
          else if (crys(phase)%structure .eq. "bct") then
            ! Next, check if the hardening is isotropic or anisotropic.
            write (io, '(32(e13.7,1x))') (result_buffer_arr(:, j))
          end if
        end do
      end if
      deallocate (result_buffer_arr)
      close (io)
    end if

    ! Elemental elastic strain tensors (6-vector)
    ! Written as 11, 22, 33, 23, 13, 12 in sample basis
    if (printing%print_strain_el) then
      !
      io = printing%strain_el_u
      filename = trim(dir_name)//'/results/elts/strain_el/strain_el.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      !
      numdim = 6
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, numdim, 1, mesh%num_elts)
      else
        allocate (result_buffer_arr(numdim, 1:mesh%num_elts))
        call vec6_crystosam(results%elas_tot6(:, :, cqpt), results%ori(:, :, :, cqpt), elas_sam)
        !
        call par_gatherv_array(elas_sam, result_buffer_arr, numdim*mesh%global_info(1, :), &
             & numdim*mesh%global_info(4, :))
        if (myid == 0) call write_res_array(result_buffer_arr, io, numdim, step, 1, mesh%num_elts)
        deallocate (result_buffer_arr)
      end if
      close (io)
    end if

    ! Elemental plastic strain tensors
    ! Written as 11, 22, 33, 23, 13, 12 in sample basis
    if (printing%print_strain_pl) then
      !
      io = printing%strain_pl_u
      filename = trim(dir_name)//'/results/elts/strain_pl/strain_pl.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      !
      numdim = 6
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, numdim, 1, mesh%num_elts)
      else
        allocate (result_buffer_arr(numdim, 1:mesh%num_elts))
        call vec_mat_symm(results%strain_pl, tmp_mat3x3)
        call mat_vec6(tmp_mat3x3, plstrain6)
        call vec6_crystosam(plstrain6, results%ori(:, :, :, cqpt), plstrain6_sam)
        !
        call par_gatherv_array(plstrain6_sam, result_buffer_arr, numdim*mesh%global_info(1, :), &
             & numdim*mesh%global_info(4, :))
        if (myid == 0) call write_res_array(result_buffer_arr, io, numdim, step, 1, mesh%num_elts)
        !
        deallocate (result_buffer_arr)
      end if
      close (io)
    end if

    ! Elemental total strain tensors
    ! Written as 11, 22, 33, 23, 13, 12 in sample basis
    if (printing%print_strain) then
      !
      io = printing%strain_u
      filename = trim(dir_name)//'/results/elts/strain/strain.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      !
      numdim = 9
      if (step .eq. 0) then
        numdim = 6
        if (myid == 0) call write_res_zeros(io, numdim, 1, mesh%num_elts)
      else
        allocate (result_buffer_tensor(3, 3, mesh%num_elts))
        call par_gatherv_tensor(results%strain, result_buffer_tensor, numdim*mesh%global_info(1, :)&
          &, numdim*mesh%global_info(4, :))
        if (myid == 0) call write_res_tensor(result_buffer_tensor, io, 6, step, 1, mesh%num_elts)
        deallocate (result_buffer_tensor)
      end if
      close (io)
    end if

    ! Elemental stress tensors (3x3 matrix)
    ! Written as 11, 22, 33, 23, 13, 12 in sample basis
    if (printing%print_stress) then
      !
      io = printing%stress_u
      filename = trim(dir_name)//'/results/elts/stress/stress.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      !
      numdim = 9
      if (step .eq. 0) then
        numdim = 6
        if (myid == 0) call write_res_zeros(io, numdim, 1, mesh%num_elts)
      else
        allocate (result_buffer_tensor(3, 3, mesh%num_elts))
        call par_gatherv_tensor(results%stress, result_buffer_tensor, numdim*mesh%global_info(1, :)&
          &, numdim*mesh%global_info(4, :))
        if (myid == 0) call write_res_tensor(result_buffer_tensor, io, 6, step, 1, mesh%num_elts)
        deallocate (result_buffer_tensor)
      end if
      close (io)
    end if

    !!! Needs thought on how to handle phase
    ! Elemental shear rates (sliprate)
    if (printing%print_sliprate) then
      io = printing%sliprate_u
      filename = trim(dir_name)//'/results/elts/sliprate/sliprate.'//trim(step_num)
      if (myid == 0) open (io, file=filename)

      numdim = mesh%maxnumslip
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, numdim, 1, mesh%num_elts)
      else
        allocate (result_buffer_arr(numdim, mesh%num_elts))
        call par_gatherv_array(results%sliprate(:, :, cqpt), result_buffer_arr, numdim*mesh%global_info(1, :)&
          &, numdim*mesh%global_info(4, :))
        if (myid .eq. 0) then
          do i = 1, mesh%num_elts
            phase = mesh%elt_phase(i)
            write (fmt, *) '(', crys(phase)%numslip, '(e13.7,1x))'

            write (io, fmt) (anint(result_buffer_arr(j, i)/eps)*eps, j=1, crys(phase)%numslip)
          end do
        end if
        deallocate (result_buffer_arr)
      end if
      close (io)
    end if

    !!! Needs thought on how to handle phase
    ! Elemental slip system shear stresses (rss)
    if (printing%print_rss) then
      io = printing%rss_u
      filename = trim(dir_name)//'/results/elts/rss/rss.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      numdim = mesh%maxnumslip
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, numdim, 1, mesh%num_elts)
      else
        allocate (result_buffer_arr(numdim, mesh%num_elts))
        call par_gatherv_array(results%rss(:, :, cqpt), result_buffer_arr, numdim*mesh%global_info(1, :), &
             & numdim*mesh%global_info(4, :))
        if (myid .eq. 0) then
          do i = 1, mesh%num_elts
            phase = mesh%elt_phase(i)
            write (fmt, *) '(', crys(phase)%numslip, '(e13.7,1x))'

            write (io, fmt) (anint(result_buffer_arr(j, i)/eps)*eps, j=1, crys(phase)%numslip)
          end do
        end if
        deallocate (result_buffer_arr)
      end if
      close (io)
    end if

    ! Elemental equivalent stress (scalar)
    if (printing%print_stress_eq) then
      io = printing%stress_eq_u
      filename = trim(dir_name)//'/results/elts/stress_eq/stress_eq.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, 1, 1, mesh%num_elts)
      else
        allocate (result_buffer(1:mesh%num_elts))
        call stress_equiv_3x3(results%stress, stress_eq)
        call par_gatherv(stress_eq, result_buffer, mesh%global_info(1, :), mesh%global_info(4, :))
        if (myid == 0) call write_res_vector(result_buffer, io, 1, step, 1, mesh%num_elts)
        deallocate (result_buffer)
      end if
      close (io)
    end if

    ! Elemental effective deformation rate (scalar)
    if (printing%print_defrate_eq) then
      io = printing%defrate_eq_u
      filename = trim(dir_name)//'/results/elts/defrate_eq/defrate_eq.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, 1, 1, mesh%num_elts)
      else
        allocate (result_buffer(1:mesh%num_elts))
        call par_gatherv(results%defrate_eq(:, cqpt), result_buffer, mesh%global_info(1, :), mesh%global_info(4, :))
        if (myid == 0) call write_res_vector(result_buffer, printing%defrate_eq_u, 1, step, 1, mesh%num_elts)
        deallocate (result_buffer)
      end if
      close (io)
    end if

    ! Elemental equivalent strain (scalar)
    if (printing%print_strain_eq) then
      io = printing%strain_eq_u
      filename = trim(dir_name)//'/results/elts/strain_eq/strain_eq.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, 1, 1, mesh%num_elts)
      else
        allocate (result_buffer(1:mesh%num_elts))
        call strain_equiv_3x3(results%strain, elt_val)
        call par_gatherv(elt_val, result_buffer, mesh%global_info(1, :), mesh%global_info(4, :))
        if (myid == 0) call write_res_vector(result_buffer, io, 1, step, 1, mesh%num_elts)
        deallocate (result_buffer)
      end if
      close (io)
    end if

    ! Elemental effective plastic deformation rate (scalar)
    if (printing%print_defrate_pl_eq) then
      io = printing%defrate_pl_eq_u
      filename = trim(dir_name)//'/results/elts/defrate_pl_eq/defrate_pl_eq.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, 1, 1, mesh%num_elts)
      else
        allocate (result_buffer(1:mesh%num_elts))
        call strain_equiv(results%dp_hat(:, :, cqpt), elt_val)
        call par_gatherv(elt_val, result_buffer, mesh%global_info(1, :), mesh%global_info(4, :))
        if (myid == 0) call write_res_vector(result_buffer, io, 1, step, 1, mesh%num_elts)
        deallocate (result_buffer)
      end if
      close (io)
    end if

    ! Elemental equivalent plastic strain (scalar)
    if (printing%print_strain_pl_eq) then
      io = printing%strain_pl_eq_u
      filename = trim(dir_name)//'/results/elts/strain_pl_eq/strain_pl_eq.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, 1, 1, mesh%num_elts)
      else
        allocate (result_buffer(1:mesh%num_elts))
        call vec_mat_symm(results%strain_pl, tmp_mat3x3)
        call strain_equiv_3x3(tmp_mat3x3, elt_val)
        call par_gatherv(elt_val, result_buffer, mesh%global_info(1, :), mesh%global_info(4, :))
        if (myid == 0) call write_res_vector(result_buffer, io, 1, step, 1, mesh%num_elts)
        deallocate (result_buffer)
      end if
      close (io)
    end if

    ! Elemental equivalent elastic strain (scalar)
    if (printing%print_strain_el_eq) then
      io = printing%strain_el_eq_u
      filename = trim(dir_name)//'/results/elts/strain_el_eq/strain_el_eq.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, 1, 1, mesh%num_elts)
      else
        allocate (result_buffer(1:mesh%num_elts))
        call vec6_mat_symm(results%elas_tot6(:, :, cqpt), tmp_mat3x3, elt_sup - elt_sub + 1)
        call strain_equiv_3x3(tmp_mat3x3, elt_val)
        call par_gatherv(elt_val, result_buffer, mesh%global_info(1, :), mesh%global_info(4, :))
        if (myid == 0) call write_res_vector(result_buffer, io, 1, step, 1, mesh%num_elts)
        deallocate (result_buffer)
      end if
      close (io)
    end if

    ! Elemental vel gradient (3x3 matrix)
    ! Written as 11, 12, 13, 21, 22, 23, 31, 32, 33 in sample basis
    if (printing%print_velgrad) then
      io = printing%velgrad_u
      filename = trim(dir_name)//'/results/elts/velgrad/velgrad.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      numdim = 9
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, numdim, 1, mesh%num_elts)
      else
        allocate (result_buffer_tensor(3, 3, mesh%num_elts))
        call par_gatherv_tensor(results%velgrad(:, :, :, cqpt), result_buffer_tensor,&
          &numdim*mesh%global_info(1, :), numdim*mesh%global_info(4, :))
        if (myid == 0) call write_res_tensor(result_buffer_tensor, io, numdim, step, 1, mesh%num_elts)
        deallocate (result_buffer_tensor)
      end if
      close (io)
    end if

    !!! Needs work
    ! Elemental plastic deformation rate (6-vector)
    if (printing%print_defrate_pl) then
      io = printing%defrate_pl_u
      filename = trim(dir_name)//'/results/elts/defrate_pl/defrate_pl.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      numdim = 6
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, numdim, 1, mesh%num_elts)
      else
        do i = elt_sub, elt_sup
          ! Convert the transpose of the elemental rotation matrix into
          ! an operator on 5-vectors. The transpose outputs qr5x5 in
          ! a crystal-to-sample form.
          call rotmat_symm(transpose(results%ori(:, :, i, cqpt)), qr5x5, 1)
          ! lattice_deform transforms dp_hat from crystal to sample frame.
          ! This is because oris was transposed (above).
          ! 5-vector order is maintained here as (11-22), 33, 12, 13, 23
          ! with the proper scalings.
          call lattice_deform_(qr5x5, results%dp_hat(:, i, cqpt), dp_hat_sam)
          ! Convert 5-vector with scaling into proper 6-vector form.
          call vec5_vec6(dp_hat_sam, dp_hat_mat(:, i))
        end do
        allocate (result_buffer_arr(numdim, 1:mesh%num_elts))
        ! Written as 12, 13, 23 in sample basis.
        call par_gatherv_array(dp_hat_mat, result_buffer_arr, numdim*mesh%global_info(1, :), &
             & numdim*mesh%global_info(4, :))
        if (myid == 0) call write_res_array(result_buffer_arr, io, numdim, step, 1, mesh%num_elts)
        deallocate (result_buffer_arr)
      end if
      close (io)
    end if

    !!! Needs work
    ! Elemental plastic spin rate (3-vector)
    if (printing%print_spinrate) then
      io = printing%spinrate_u
      filename = trim(dir_name)//'/results/elts/spinrate/spinrate.'//trim(step_num)
      numdim = 3
      if (myid == 0) open (io, file=filename)
      allocate (result_buffer_arr(numdim, 1:mesh%num_elts))
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, numdim, 1, mesh%num_elts)
      else
        do i = elt_sub, elt_sup
          ! Convert the transpose of the elemental rotation matrix into
          ! an operator on skew 3-vectors. The transpose outputs qr3x3 in
          ! a crystal-to-sample form.
          call rotmat_skew_(transpose(results%ori(:, :, i, cqpt)), qr3x3)
          ! lattice_spin transforms wp_hat from crystal to sample frame.
          ! This is because oris was transposed (above).
          ! Skew 3-vector order is maintained here as 21 31 32.
          call lattice_spin_(qr3x3, results%wp_hat(:, i, cqpt), wp_hat_sam(:, i))
          ! Negative values are output here in order to obtain desired
          ! order while maintaining skew-symmetry constraints.
        end do
        ! Written as 12, 13, 23 in sample basis.
        call par_gatherv_array(-wp_hat_sam, result_buffer_arr, numdim*mesh%global_info(1, :), numdim*mesh%global_info(4, :))
        if (myid == 0) call write_res_array(result_buffer_arr, io, numdim, step, 1, mesh%num_elts)
        deallocate (result_buffer_arr)
      end if
      close (io)
    end if

    !!! Needs work
    ! Elemental plastic spin rate (3-vector)
    if (printing%print_rotrate) then
      io = printing%rotrate_u
      filename = trim(dir_name)//'/results/elts/rotrate/rotrate.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      numdim = 3
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, numdim, 1, mesh%num_elts)
      else
        call rotmat_to_axis_angle(results%d_rstar(:, :, :, cqpt), tmp_vec4)
        tmp_vec3(1, :) = tmp_vec4(1, :)*tmp_vec4(4, :)*pi_over_180
        tmp_vec3(2, :) = tmp_vec4(2, :)*tmp_vec4(4, :)*pi_over_180
        tmp_vec3(3, :) = tmp_vec4(3, :)*tmp_vec4(4, :)*pi_over_180
        tmp_vec3 = -tmp_vec3
        tmp_vec3 = tmp_vec3/dtime

        do i = elt_sub, elt_sup
          ! Convert the transpose of the elemental rotation matrix into
          ! an operator on rotation vectors. The transpose outputs qr3x3 in
          ! a crystal-to-sample form.
          call rotmat_skew2(transpose(results%ori(:, :, i, cqpt)), qr3x3, 1)
          call lattice_spin_(qr3x3, tmp_vec3(:, i), tmp_vec3(:, i))
          ! lattice_rot transforms wp_hat from crystal to sample frame.
          ! This is because oris was transposed (above).
          ! Skew 3-vector order is maintained here as 21 31 32.
          ! Negative values are output here in order to obtain desired
          ! order while maintaining skew-symmetry constraints.
        end do
        allocate (result_buffer_arr(numdim, 1:mesh%num_elts))
        ! Written as 12, 13, 23 in sample basis.
        ! write (io, '(3(e13.7,1x))') tmp_vec3(1, i), tmp_vec3(2, i), tmp_vec3(3, i)
        call par_gatherv_array(tmp_vec3, result_buffer_arr, numdim*mesh%global_info(1, :), numdim*mesh%global_info(4, :))
        if (myid == 0) call write_res_array(result_buffer_arr, io, numdim, step, 1, mesh%num_elts)
        deallocate (result_buffer_arr)
      end if
    end if

    ! Elemental plastic spin rate (3-vector)
    if (printing%print_rotrate_spin) then

      io = printing%rotrate_spin_u
      filename = trim(dir_name)//'/results/elts/rotrate_spin/rotrate_spin.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      numdim = 3
      !
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, numdim, 1, mesh%num_elts)
      else
        call rotmat_to_axis_angle(results%d_rstar_spin(:, :, :, cqpt), tmp_vec4)

        tmp_vec3(1, :) = tmp_vec4(1, :)*tmp_vec4(4, :)*pi_over_180
        tmp_vec3(2, :) = tmp_vec4(2, :)*tmp_vec4(4, :)*pi_over_180
        tmp_vec3(3, :) = tmp_vec4(3, :)*tmp_vec4(4, :)*pi_over_180
        tmp_vec3 = -tmp_vec3
        tmp_vec3 = tmp_vec3/dtime

        do i = elt_sub, elt_sup
          ! Convert the transpose of the elemental rotation matrix into
          ! an operator on rotation vectors. The transpose outputs qr3x3 in
          ! a crystal-to-sample form.
          call rotmat_skew2(transpose(results%ori(:, :, i, cqpt)), qr3x3, 1)
          call lattice_spin_(qr3x3, tmp_vec3(:, i), tmp_vec3(:, i))
          ! lattice_rot transforms wp_hat from crystal to sample frame.
          ! This is because oris was transposed (above).
          ! Skew 3-vector order is maintained here as 21 31 32.
          ! Negative values are output here in order to obtain desired
          ! order while maintaining skew-symmetry constraints.
        end do
        allocate (result_buffer_arr(numdim, 1:mesh%num_elts))
        ! Written as 12, 13, 23 in sample basis.
        call par_gatherv_array(tmp_vec3, result_buffer_arr, numdim*mesh%global_info(1, :),&
          & numdim*mesh%global_info(4, :))
        if (myid == 0) call write_res_array(result_buffer_arr, io, numdim, step, 1, mesh%num_elts)
        deallocate (result_buffer_arr)
      end if
      close (io)
    end if

    ! Elemental plastic slip rate (3-vector)
    if (printing%print_rotrate_slip) then
      io = printing%rotrate_slip_u
      filename = trim(dir_name)//'/results/elts/rotrate_slip/rotrate_slip.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      numdim = 3
      !
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, numdim, 1, mesh%num_elts)
      else
        call rotmat_to_axis_angle(results%d_rstar_slip(:, :, :, cqpt), tmp_vec4)

        tmp_vec3(1, :) = tmp_vec4(1, :)*tmp_vec4(4, :)*pi_over_180
        tmp_vec3(2, :) = tmp_vec4(2, :)*tmp_vec4(4, :)*pi_over_180
        tmp_vec3(3, :) = tmp_vec4(3, :)*tmp_vec4(4, :)*pi_over_180
        tmp_vec3 = -tmp_vec3
        tmp_vec3 = tmp_vec3/dtime

        do i = elt_sub, elt_sup
          ! Convert the transpose of the elemental rotation matrix into
          ! an operator on rotation vectors. The transpose outputs qr3x3 in
          ! a crystal-to-sample form.
          call rotmat_skew2(transpose(results%ori(:, :, i, cqpt)), qr3x3, 1)
          call lattice_spin_(qr3x3, tmp_vec3(:, i), tmp_vec3(:, i))
          ! lattice_rot transforms wp_hat from crystal to sample frame.
          ! This is because oris was transposed (above).
          ! Skew 3-vector order is maintained here as 21 31 32.
          ! Negative values are output here in order to obtain desired
          ! order while maintaining skew-symmetry constraints.
          ! Written as 12, 13, 23 in sample basis.
        end do
        allocate (result_buffer_arr(numdim, 1:mesh%num_elts))
        ! Written as 12, 13, 23 in sample basis.
        call par_gatherv_array(tmp_vec3, result_buffer_arr, numdim*mesh%global_info(1, :), numdim*mesh%global_info(4, :))
        if (myid == 0) call write_res_array(result_buffer_arr, io, numdim, step, 1, mesh%num_elts)
        deallocate (result_buffer_arr)
      end if
      close (io)
    end if

    !!! Needs thought on how to handle phase
    ! Elemental accumulated shears
    if (printing%print_slip) then
      io = printing%slip_u
      filename = trim(dir_name)//'/results/elts/slip/slip.'//trim(step_num)
      if (myid == 0) open (io, file=filename)

      numdim = mesh%maxnumslip
      !
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, numdim, 1, mesh%num_elts)
      else
        allocate (result_buffer_arr(numdim, mesh%num_elts))
        call par_gatherv_array(results%slip(:, :), result_buffer_arr, numdim*mesh%global_info(1, :),&
         & numdim*mesh%global_info(4, :))
        !
        if (myid .eq. 0) then
          do i = 1, mesh%num_elts
            phase = mesh%elt_phase(i)
            write (fmt, *) '(', crys(phase)%numslip, '(e13.7,1x))'

            write (io, fmt) (anint(result_buffer_arr(j, i)/eps)*eps, j=1, crys(phase)%numslip)
          end do
        end if
        deallocate (result_buffer_arr)
      end if
      close (io)
    end if

    ! Elemental total work (scalar)
    if (printing%print_work) then
      io = printing%work_u
      filename = trim(dir_name)//'/results/elts/work/work.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      !
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, 1, 1, mesh%num_elts)
      else
        allocate (result_buffer(1:mesh%num_elts))
        call par_gatherv(results%work, result_buffer, mesh%global_info(1, :), mesh%global_info(4, :))
        if (myid == 0) call write_res_vector(result_buffer, io, 1, step, 1, mesh%num_elts)
        deallocate (result_buffer)
      end if
      close (io)
    end if

    ! Elemental plastic work (scalar)
    if (printing%print_work_pl) then
      io = printing%work_pl_u
      filename = trim(dir_name)//'/results/elts/work_pl/work_pl.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      !
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, 1, 1, mesh%num_elts)
      else
        allocate (result_buffer(1:mesh%num_elts))
        call par_gatherv(results%work_pl, result_buffer, mesh%global_info(1, :), mesh%global_info(4, :))
        if (myid == 0) call write_res_vector(result_buffer, io, 1, step, 1, mesh%num_elts)
        deallocate (result_buffer)
      end if
      close (io)
    end if

    ! Elemental total deformation rate tensor (3x3 matrix)
    ! Written as 11, 22, 33, 23, 13, 12 in sample basis
    if (printing%print_defrate) then
      io = printing%defrate_u
      filename = trim(dir_name)//'/results/elts/defrate/defrate.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      numdim = 9
      !
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, 6, 1, mesh%num_elts)
      else
        allocate (result_buffer_tensor(3, 3, mesh%num_elts))
        call par_gatherv_tensor(results%defrate, result_buffer_tensor, numdim*mesh%global_info(1, :), &
             & numdim*mesh%global_info(4, :))
        if (myid == 0) call write_res_tensor(result_buffer_tensor, io, 6, step, 1, mesh%num_elts)
        deallocate (result_buffer_tensor)
      end if
      close (io)
    end if

    ! Elemental total work rate (scalar)
    if (printing%print_workrate) then
      io = printing%workrate_u
      filename = trim(dir_name)//'/results/elts/workrate/workrate.'//trim(step_num)
      if (myid == 0) open (io, file=filename)
      !
      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, 1, 1, mesh%num_elts)
      else
        allocate (result_buffer(1:mesh%num_elts))
        call par_gatherv(results%workrate, result_buffer, mesh%global_info(1, :), mesh%global_info(4, :))
        if (myid == 0) call write_res_vector(result_buffer, io, 1, step, 1, mesh%num_elts)
        deallocate (result_buffer)
      end if
      close (io)
    end if

    ! Elemental plastic work rate (scalar)
    if (printing%print_workrate_pl) then
      io = printing%workrate_pl_u
      filename = trim(dir_name)//'/results/elts/workrate_pl/workrate_pl.'//trim(step_num)
      if (myid == 0) open (io, file=filename)

      if (step .eq. 0) then
        if (myid == 0) call write_res_zeros(io, 1, 1, mesh%num_elts)
      else
        allocate (result_buffer(1:mesh%num_elts))
        call par_gatherv(results%workrate_pl, result_buffer, mesh%global_info(1, :), mesh%global_info(4, :))
        if (myid == 0) call write_res_vector(result_buffer, io, 1, step, 1, mesh%num_elts)
        deallocate (result_buffer)
      end if
      close (io)
    end if

    call write_dot_sim_file(printing, mesh, step)

  end subroutine write_res

  !===========================================================================

  subroutine write_force_file_headers(mesh, printing, driver_type)

    ! This writes the file headers for surface force files. Consistent format
    ! across all drivers is maintained here. Assumes that the call is wrapped
    ! in `if (printing%print_forces) then` logic from where it is called.
    !---------------------------------------------------------------------------

    ! Arguments:
    ! driver_type: Flag that denotes if the calling driver is uni- or triaxial
    !   (1) = loading, (2) = triaxial

    type(mesh_type), intent(in) :: mesh
    type(printing_type), intent(in) :: printing
    integer :: driver_type

    integer :: i

    ! Locals:
    ! force_headeru1/2: Headers for the uniaxial drivers
    ! force_headert1/2: Headers for the triaxial drivers

    character(len=*), parameter :: force_headeru1 = &
        &   '% step     incr     Fx             Fy            &
        & Fz             area a         time'
    character(len=*), parameter :: force_headeru2 = &
        &   '% -------- -------- -------------- --------------&
        & -------------- -------------- --------------'

    character(len=*), parameter :: force_headert1 = &
        &   '% step     incr     Fx             Fy            &
        & Fz             area a         time           length'
    character(len=*), parameter :: force_headert2 = &
        &   '% -------- -------- -------------- --------------&
        & -------------- -------------- -------------- -----------------'

    !---------------------------------------------------------------------------

    ! Check the input driver_type and print the correct headers
    do i = 1, mesh%num_fasets

      if (driver_type .eq. 1) then ! uniaxial
        write (printing%force_u + i - 1, '(a)') force_headeru1
        write (printing%force_u + i - 1, '(a)') force_headeru2

      else if (driver_type .eq. 2) then ! triaxial
        write (printing%force_u + i - 1, '(a)') force_headert1
        write (printing%force_u + i - 1, '(a)') force_headert2

      end if

    end do

  end subroutine write_force_file_headers

  !===========================================================================

  subroutine write_force_file_data(mesh, printing, istep, incr, load, &
      & area, time, length)

    ! This writes the actual data for surface force files. Consistent format
    ! across all drivers is maintained here. Assumes that the call is wrapped
    ! in `if (printing%print_forces) then` logic from where it is called.
    !---------------------------------------------------------------------------

    ! Arguments:
    !   (1) = loading, (2) = triaxial
    ! istep: Current timestep being printed
    ! incr: Current total increment being printed
    ! load: Contains the [x y z] loads on all surfaces
    ! area: Current surface face areas bring printed
    ! time: Current time value being printed
    ! length: Current length of the mesh edges (triaxial only)

    type(mesh_type), intent(in) :: mesh
    type(printing_type), intent(in) :: printing
    integer :: istep
    integer :: incr
    real(rk) :: load(:, :)
    real(rk) :: area(:)
    real(rk) :: time
    real(rk), optional :: length(:)

    integer :: i

    !---------------------------------------------------------------------------

    if (.not. present(length)) then
      ! Print the data to the files
      do i = 1, mesh%num_fasets
        write (printing%force_u + i - 1, '(2(i8),5(e15.5))') istep, incr, load(i, :), area(i), time
      end do

    else
      ! Print the data to the files
      do i = 1, mesh%num_fasets
        write (printing%force_u + i - 1, '(2(i8),5(e15.5),e18.8)') istep, incr, &
               & load(i, :), area(i), time, length(i)
      end do
    end if

  end subroutine write_force_file_data

  !===========================================================================

  subroutine write_conv_file_headers(printing)

    ! This writes the file headers for convergence report. Consistent format
    ! across all drivers is maintained here. Assumes that the call is wrapped
    ! in `if (printing%print_conv) then` logic from where it is called.
    !---------------------------------------------------------------------------
    type(printing_type), intent(in) :: printing

    ! Print the headers
    write (printing%conv_u, '(a)') '%   incr     iter       NR     r_norm&
        &        rx_norm       f_norm        delu_norm     delux_norm&
        &    u_norm      cg_iter'

  end subroutine write_conv_file_headers

  !===========================================================================

  subroutine write_conv_file_data(printing, incr, iter, itmethod, r_norm, rx_norm, &
      & f_norm, delu_norm, delux_norm, u_norm, cg_iter_out)

    ! This writes the file data for convergence report. This is only called
    ! from the solveit_evp subroutine. Assumes that the call is wrapped in
    ! `if (printing%print_conv) then` logic from where it is called.
    ! Arguments:
    ! incr: Current total increment value
    ! iter: Current iteration within step
    ! itmethod: "SA" for successive approx., NR for newton-raphson
    ! r_norm: Residual norm -> sqrt(sum(resid * resid))
    ! rx_norm: Maximal absolute residual -> maxval(abs(resid))
    ! f_norm: Force norm -> sqrt(sum(force * force))
    ! delu_norm: Change in vel norm -> sqrt(sum(delta_vel * delta_vel))
    ! delux_norm: Maximal absolute change in vel -> maxval(abs(delta_vel))
    ! u_norm: Velocity norm -> sqrt(sum(vel_o * vel_o))
    ! cg_iter_out: Number of conjugate-gradient (cg) iterations for this incr

    type(printing_type), intent(in) :: printing
    integer, intent(in)  :: incr, iter, cg_iter_out
    character(len=2), intent(in) :: itmethod
    real(rk), intent(in) :: r_norm, rx_norm, f_norm, delu_norm, delux_norm, u_norm
    !---------------------------------------------------------------------------

    ! Print the data to the files
    write (printing%conv_u, '(i8,1x,i8,1x,a,6d14.4,1x,i8)') incr, iter, &
        & itmethod, r_norm, rx_norm, f_norm, delu_norm, delux_norm, u_norm, &
        & cg_iter_out

  end subroutine write_conv_file_data

  subroutine write_dot_sim_file (printing, mesh, step)

    type(printing_type), intent(in) :: printing
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: step

    character(len=256):: filename
    character(len=256):: dir_name

    dir_name = "simulation.sim"
    filename = trim(dir_name)//'/.sim'

    if (myid .eq. 0) then
      open (printing%dot_sim_u, file=filename)

      call write_dot_sim_file_header(printing, mesh)
      call write_dot_sim_file_output_files(printing)
      call write_dot_sim_file_complete_steps(printing, step)

      close (printing%dot_sim_u)
    end if

  end subroutine write_dot_sim_file

end module write_res_mod
