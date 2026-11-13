! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of orientation_mod: descriptor names, Neper's documented 30°
!> example, Neper-written fixtures in every descriptor and convention, and round trips.
!>
!> Source : FEPX (ROADMAP step 1.3; conventions.md §2).
!>
!> Usage: test_orientation <fixture directory> (tests/math/fixtures/neper_ori).
!>
!> Checks:
!>   - names <-> codes, parameter counts, conventions;
!>   - Neper's documentation example (doc/exprskeys.rst, Neper 5.0): a rotation of 30° about X_s,
!>     passive values of all 7 descriptors -> rot_c2s whose columns are the crystal axes of the
!>     documentation figure: X_c = X_s, Y_c = (0, cos 30°, sin 30°), Z_c = (0, -sin 30°,
!>     cos 30°); the opposite rotation in the active convention gives the same rot_c2s;
!>     tolerance 1e-9 (the documented values have 9 digits);
!>   - Neper fixtures: 20 orientations x 7 descriptors x 2 conventions all give the same
!>     rot_c2s (= transpose of the passive rotmat) to 1e-10 (12-decimal files), and
!>     descriptor_from_rot_c2s writes Neper's values back (Euler angles modulo 360°) to 1e-9;
!>   - random rotations (500) and edge cases (identity, 180° about general and coordinate axes,
!>     Phi = 0 and 180°, axes in a coordinate plane, 1e-6° rotation): rot_c2s -> parameters ->
!>     rot_c2s for every descriptor and convention to 1e-12, except Euler angles near
!>     Phi = 0/180° (not exactly 0/180°), where the individual angles are determined only to
!>     about eps/sin(Phi): tolerance max(1e-12, 1e-15/sin(Phi)) (measured 3e-9 for the 1e-6°
!>     rotation, sin(Phi) = 1.7e-8); quaternion q0 >= 0; axis-angle angle in [0, 180]; small
!>     angles with full relative precision;
!>   - normalization: a quaternion or axis given with length 2 gives the same rotation.
!>
!> Locals:
!>   fixture_dir  character(512)  fixture directory (command-line argument)
!>   des, conv    integer         descriptor and convention codes
!>   i, j, k      integer         loop indices
!>   num          integer         number of descriptor parameters
!>   params       real(rk)(9)     descriptor parameters
!>   back         real(rk)(9)     parameters written back
!>   neper        real(rk)(9,20)  parameters of one fixture file
!>   rot_ref      real(rk)(3,3,20) reference rot_c2s of the fixture orientations
!>   rot, rot2    real(rk)(3,3)   rotations
!>   c30, s30     real(rk)        cos 30°, sin 30°
!>   rot_doc      real(rk)(3,3)   rot_c2s of Neper's documentation example
!>   edge         real(rk)(3,3,12) edge-case rotations
!>   label        character(96)   check label
!>   sin_phi      real(rk)        sin(Phi) of a rotation (Phi = Bunge/Roe Phi = Kocks Theta)
!>   tol_rt       real(rk)        round-trip tolerance
program test_orientation

  use kinds_mod, only: rk
  use constants_mod, only: identity_3, deg_to_rad
  use orientation_mod, only: ori_des_rodrigues, ori_des_euler_bunge, ori_des_euler_kocks, &
    & ori_des_euler_roe, ori_des_rotmat, ori_des_axis_angle, ori_des_quaternion, &
    & ori_conv_passive, ori_conv_active, ori_descriptor_from_name, ori_descriptor_name, &
    & ori_num_params, ori_convention_from_name, ori_convention_name, &
    & rot_c2s_from_descriptor, descriptor_from_rot_c2s, g_from_axis_angle, g_from_euler_bunge
  use tensor_test_utils_mod, only: seed_random, random_rotation, determinant3
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary

  implicit none

  character(len=512) :: fixture_dir
  integer :: des, conv, i, j, k, num
  real(rk) :: params(9), back(9)
  real(rk) :: neper(9, 20)
  real(rk) :: rot_ref(3, 3, 20)
  real(rk) :: rot(3, 3), rot2(3, 3)
  real(rk) :: c30, s30
  real(rk) :: rot_doc(3, 3)
  real(rk) :: edge(3, 3, 12)
  character(len=96) :: label
  real(rk) :: sin_phi, tol_rt

  call get_command_argument(1, fixture_dir)
  call seed_random(20261008)

  ! ------------------------------------------------------------------------------------------
  ! Names and codes.
  ! ------------------------------------------------------------------------------------------

  do des = 1, 7
    call assert_equal(ori_descriptor_from_name(ori_descriptor_name(des)), des, &
      & 'name_round_trip_' // ori_descriptor_name(des))
  end do
  call assert_equal(ori_descriptor_from_name('euler-bunge'), ori_des_euler_bunge, 'name_bunge')
  call assert_equal(ori_descriptor_from_name('rotmat'), ori_des_rotmat, 'name_rotmat')
  call assert_equal(ori_descriptor_from_name('euler_bunge'), 0, 'name_unknown')
  call assert_true(all([(ori_num_params(des), des=1, 7)] == [3, 3, 3, 3, 9, 4, 4]), 'num_params')
  call assert_equal(ori_num_params(8), 0, 'num_params_invalid')
  call assert_equal(ori_convention_from_name('passive'), ori_conv_passive, 'conv_passive')
  call assert_equal(ori_convention_from_name('active'), ori_conv_active, 'conv_active')
  call assert_equal(ori_convention_from_name('Passive'), 0, 'conv_unknown')
  call assert_equal(ori_convention_name(ori_conv_active), 'active', 'conv_name')

  ! ------------------------------------------------------------------------------------------
  ! Neper's documentation example: 30° about X_s (passive), and the opposite rotation (active).
  ! ------------------------------------------------------------------------------------------

  c30 = cos(30.0_rk*deg_to_rad)
  s30 = sin(30.0_rk*deg_to_rad)
  ! Columns: crystal axes in the sample frame, as in Neper's figure.
  rot_doc(:, 1) = [1.0_rk, 0.0_rk, 0.0_rk]
  rot_doc(:, 2) = [0.0_rk, c30, s30]
  rot_doc(:, 3) = [0.0_rk, -s30, c30]

  call check_doc(ori_des_rodrigues, ori_conv_passive, [0.267949192_rk])
  call check_doc(ori_des_euler_bunge, ori_conv_passive, [0.0_rk, 30.0_rk, 0.0_rk])
  call check_doc(ori_des_euler_kocks, ori_conv_passive, [270.0_rk, 30.0_rk, 90.0_rk])
  call check_doc(ori_des_euler_roe, ori_conv_passive, [270.0_rk, 30.0_rk, 90.0_rk])
  call check_doc(ori_des_rotmat, ori_conv_passive, [1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
    & 0.866025404_rk, 0.5_rk, 0.0_rk, -0.5_rk, 0.866025404_rk])
  call check_doc(ori_des_axis_angle, ori_conv_passive, [1.0_rk, 0.0_rk, 0.0_rk, 30.0_rk])
  call check_doc(ori_des_quaternion, ori_conv_passive, [0.965925826_rk, 0.258819045_rk, &
    & 0.0_rk, 0.0_rk])

  ! Active convention: the opposite rotation (inverse) gives the same orientation. Euler
  ! angles of the inverse of Bunge (0, 30, 0) are (180, 30, 180); Kocks Psi = phi1 - 90,
  ! phi = 90 - phi2; Roe psi = phi1 - 90, phi = phi2 + 90.
  call check_doc(ori_des_rodrigues, ori_conv_active, [-0.267949192_rk])
  call check_doc(ori_des_euler_bunge, ori_conv_active, [180.0_rk, 30.0_rk, 180.0_rk])
  call check_doc(ori_des_euler_kocks, ori_conv_active, [90.0_rk, 30.0_rk, -90.0_rk])
  call check_doc(ori_des_euler_roe, ori_conv_active, [90.0_rk, 30.0_rk, 270.0_rk])
  call check_doc(ori_des_rotmat, ori_conv_active, [1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
    & 0.866025404_rk, -0.5_rk, 0.0_rk, 0.5_rk, 0.866025404_rk])
  call check_doc(ori_des_axis_angle, ori_conv_active, [1.0_rk, 0.0_rk, 0.0_rk, -30.0_rk])
  call check_doc(ori_des_quaternion, ori_conv_active, [0.965925826_rk, -0.258819045_rk, &
    & 0.0_rk, 0.0_rk])

  ! ------------------------------------------------------------------------------------------
  ! Neper fixtures: every descriptor and convention gives the same rot_c2s, and is written back.
  ! ------------------------------------------------------------------------------------------

  ! Reference: transpose of the passive rotmat (rot_c2s = g^T).
  call read_fixture('rotmat', 'passive', 9, neper)
  do i = 1, 20
    rot_ref(:, :, i) = reshape(neper(1:9, i), [3, 3])
  end do

  do des = 1, 7
    num = ori_num_params(des)
    do conv = ori_conv_passive, ori_conv_active
      call read_fixture(ori_descriptor_name(des), ori_convention_name(conv), num, neper)
      do i = 1, 20
        write (label, '(a,":",a,"[",i0,"]")') ori_descriptor_name(des), &
          & ori_convention_name(conv), i
        ! Descriptor -> rot_c2s.
        rot = rot_c2s_from_descriptor(des, conv, neper(1:num, i))
        call assert_close(rot, rot_ref(:, :, i), 'fixture_read_' // trim(label), abs_tol=1.0e-10_rk)
        ! rot_c2s -> descriptor: Neper's values back (Euler angles modulo 360°).
        back = descriptor_from_rot_c2s(des, conv, rot_ref(:, :, i))
        if (des == ori_des_euler_bunge .or. des == ori_des_euler_kocks &
          & .or. des == ori_des_euler_roe) then
          do k = 1, 3
            call assert_close(modulo(back(k) - neper(k, i) + 180.0_rk, 360.0_rk) - 180.0_rk, &
              & 0.0_rk, 'fixture_write_' // trim(label), abs_tol=1.0e-9_rk)
          end do
        else
          call assert_close(back(1:num), neper(1:num, i), 'fixture_write_' // trim(label), &
            & abs_tol=1.0e-9_rk)
        end if
      end do
    end do
  end do

  ! ------------------------------------------------------------------------------------------
  ! Round trips: random rotations and edge cases.
  ! ------------------------------------------------------------------------------------------

  ! Edge cases (as rot_c2s).
  edge(:, :, 1) = identity_3
  edge(:, :, 2) = transpose(g_from_axis_angle([1.0_rk, 2.0_rk, 2.0_rk]/3.0_rk, 180.0_rk))
  edge(:, :, 3) = transpose(g_from_axis_angle([1.0_rk, 0.0_rk, 0.0_rk], 180.0_rk))
  edge(:, :, 4) = transpose(g_from_axis_angle([0.0_rk, 0.0_rk, 1.0_rk], 180.0_rk))
  edge(:, :, 5) = transpose(g_from_euler_bunge([30.0_rk, 0.0_rk, 40.0_rk]))
  edge(:, :, 6) = transpose(g_from_euler_bunge([30.0_rk, 180.0_rk, 40.0_rk]))
  edge(:, :, 7) = transpose(g_from_axis_angle([0.6_rk, 0.8_rk, 0.0_rk], 77.0_rk))
  edge(:, :, 8) = transpose(g_from_axis_angle([0.0_rk, 0.6_rk, -0.8_rk], 123.0_rk))
  edge(:, :, 9) = transpose(g_from_axis_angle([1.0_rk, 2.0_rk, 2.0_rk]/3.0_rk, 1.0e-6_rk))
  edge(:, :, 10) = transpose(g_from_euler_bunge([30.0_rk, 1.0e-6_rk, 40.0_rk]))
  edge(:, :, 11) = transpose(g_from_euler_bunge([0.0_rk, 90.0_rk, 0.0_rk]))
  edge(:, :, 12) = transpose(g_from_axis_angle([0.0_rk, 1.0_rk, 0.0_rk], 90.0_rk))

  do j = 1, 500 + 12
    if (j <= 500) then
      rot = random_rotation()
    else
      rot = edge(:, :, j - 500)
    end if
    do des = 1, 7
      ! A 180° rotation has no finite Rodrigues vector.
      if (des == ori_des_rodrigues .and. (j == 502 .or. j == 503 .or. j == 504)) cycle
      num = ori_num_params(des)
      do conv = ori_conv_passive, ori_conv_active
        write (label, '(a,":",a,"[",i0,"]")') ori_descriptor_name(des), &
          & ori_convention_name(conv), j
        params = descriptor_from_rot_c2s(des, conv, rot)
        rot2 = rot_c2s_from_descriptor(des, conv, params(1:num))
        ! Euler angles near (not at) gimbal lock: tolerance scaled by 1/sin(Phi).
        tol_rt = 1.0e-12_rk
        sin_phi = sqrt(rot(3, 1)**2 + rot(3, 2)**2)
        if ((des == ori_des_euler_bunge .or. des == ori_des_euler_kocks &
          & .or. des == ori_des_euler_roe) .and. sin_phi > 0.0_rk) then
          tol_rt = max(tol_rt, 1.0e-15_rk/sin_phi)
        end if
        call assert_close(rot2, rot, 'round_trip_' // trim(label), abs_tol=tol_rt)
        if (des == ori_des_quaternion) then
          call assert_true(params(1) >= 0.0_rk, 'quaternion_q0_nonnegative_' // trim(label))
        end if
        if (des == ori_des_axis_angle) then
          call assert_true(params(4) >= 0.0_rk .and. params(4) <= 180.0_rk, &
            & 'axis_angle_range_' // trim(label))
        end if
      end do
    end do
  end do

  ! Small angles keep full relative precision (edge cases 9 and 10).
  params = descriptor_from_rot_c2s(ori_des_axis_angle, ori_conv_passive, edge(:, :, 9))
  call assert_close(params(4), 1.0e-6_rk, 'small_rotation_angle', rel_tol=1.0e-9_rk)
  params = descriptor_from_rot_c2s(ori_des_euler_bunge, ori_conv_passive, edge(:, :, 10))
  call assert_close(params(2), 1.0e-6_rk, 'small_bunge_phi', rel_tol=1.0e-9_rk)

  ! 180° about a general axis: the quaternion has q0 = 0 and the first nonzero component > 0;
  ! the axis is +-(1,2,2)/3.
  params = descriptor_from_rot_c2s(ori_des_quaternion, ori_conv_passive, edge(:, :, 2))
  call assert_close(params(1:4), [0.0_rk, 1.0_rk/3.0_rk, 2.0_rk/3.0_rk, 2.0_rk/3.0_rk], &
    & 'quaternion_180', abs_tol=1.0e-15_rk)

  ! ------------------------------------------------------------------------------------------
  ! Normalization of axis-angle and quaternion inputs.
  ! ------------------------------------------------------------------------------------------

  rot = rot_c2s_from_descriptor(ori_des_quaternion, ori_conv_passive, &
    & 2.0_rk*[0.965925826_rk, 0.258819045_rk, 0.0_rk, 0.0_rk])
  call assert_close(rot, rot_doc, 'quaternion_normalized', abs_tol=1.0e-9_rk)
  rot = rot_c2s_from_descriptor(ori_des_axis_angle, ori_conv_passive, &
    & [2.0_rk, 0.0_rk, 0.0_rk, 30.0_rk])
  call assert_close(rot, rot_doc, 'axis_normalized', abs_tol=1.0e-12_rk)
  call assert_close(determinant3(rot), 1.0_rk, 'axis_normalized_det', abs_tol=1.0e-14_rk)

  ! Report; error stop on any failure.
  call test_summary()

contains

  !> Check that the documented descriptor values give rot_doc (and that the first parameter of
  !> a Rodrigues vector is passed with zeros for the other two).
  !>
  !> Arguments:
  !>   descriptor  (in)  integer      descriptor code
  !>   convention  (in)  integer      convention code
  !>   values      (in)  real(rk)(:)  documented values (a single value for Rodrigues: x only)
  !> Locals:
  !>   vals        real(rk)(9)  parameters padded with zeros
  subroutine check_doc(descriptor, convention, values)

    integer, intent(in) :: descriptor
    integer, intent(in) :: convention
    real(rk), intent(in) :: values(:)

    real(rk) :: vals(9)

    vals = 0.0_rk
    vals(1:size(values)) = values
    call assert_close(rot_c2s_from_descriptor(descriptor, convention, &
      & vals(1:ori_num_params(descriptor))), rot_doc, 'neper_doc_30deg_' // &
      & ori_descriptor_name(descriptor) // ':' // ori_convention_name(convention), &
      & abs_tol=1.0e-9_rk)

  end subroutine check_doc

  !> Read the 20 orientations of one fixture file r_<descriptor>_<convention>.ori.
  !>
  !> Arguments:
  !>   des_name   (in)  character        descriptor name
  !>   conv_name  (in)  character        convention name
  !>   num        (in)  integer          number of values per line
  !>   values     (out) real(rk)(9,20)   values (rows 1:num)
  !> Locals:
  !>   unit       integer         file unit
  !>   iostat     integer         I/O status
  !>   path       character(640)  file path
  !>   line       integer         line index
  subroutine read_fixture(des_name, conv_name, num, values)

    character(len=*), intent(in) :: des_name
    character(len=*), intent(in) :: conv_name
    integer, intent(in) :: num
    real(rk), intent(out) :: values(9, 20)

    integer :: unit, iostat, line
    character(len=640) :: path

    values = 0.0_rk
    path = trim(fixture_dir) // '/r_' // des_name // '_' // conv_name // '.ori'
    open (newunit=unit, file=trim(path), status='old', action='read', iostat=iostat)
    call assert_equal(iostat, 0, 'open ' // trim(path))
    if (iostat /= 0) return
    do line = 1, 20
      read (unit, *, iostat=iostat) values(1:num, line)
      call assert_equal(iostat, 0, 'read ' // trim(path))
    end do
    close (unit)

  end subroutine read_fixture

end program test_orientation
