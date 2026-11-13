! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Orientation descriptors: conversion between Neper's orientation descriptors (with their
!> active/passive convention) and the internal crystal-to-sample rotation rot_c2s.
!>
!> Theory : dev/conventions.md §2. Neper's documentation (doc/exprskeys.rst, "Rotations and
!>          Orientations"; Neper 5.0): the passive convention "is based on the rotation of the
!>          sample coordinate system into the crystal coordinate system"; under it, a rotation of
!>          30° about X_s has the rotation matrix g = [[1,0,0],[0,c,s],[0,-s,c]] (c = cos 30°,
!>          s = sin 30°), and the documentation figure shows the crystal axes as the sample axes
!>          rotated by +30° about X_s. The columns of rot_c2s are the crystal axes in the sample
!>          frame, so rot_c2s = g^T for the passive convention. "The values of the orientation
!>          descriptors under the active convention are obtained by taking the opposite
!>          rotation", so rot_c2s = g for the active convention. In the user's terms (ROADMAP
!>          D11): passive = the transformation from sample to crystal (v_c = g v_s), active =
!>          the transformation from crystal to sample (v_s = g v_c); internal storage is the
!>          active form rot_c2s, so passive input is transposed once when read and once more
!>          when written.
!> Source : descriptor -> g: port of v2 orientation_conversion_mod.f90 (axis_angle_to_rotmat,
!>          euler_bunge_to_rotmat, euler_kocks_to_rotmat, rodrigues_to_rotmat, quat_to_rotmat),
!>          including v2's zeroing of |sin|, |cos| < 1e-15 (`vtiny`). Euler-Roe (not in v2): via
!>          Bunge with phi1 = psi + 90°, Phi = theta, phi2 = phi - 90° (Neper orilib
!>          ol_er_e_rad; Kocks, Tomé & Wenk (1998), Texture and Anisotropy, Ch. 2). g -> descriptor:
!>          port of v2-fixed after CHANGES.md C18 (quaternion by Shepperd's method, Shepperd,
!>          S.W. (1978), J. Guidance and Control 1(3):223-224; angles with atan2), and v2's
!>          rotmat_to_euler_bunge / rotmat_to_euler_kocks with Phi = atan2(sin, cos) (C18) and
!>          the Kocks Theta = 180° case (C19). Near Phi = 0 or 180° the individual Euler angles
!>          are determined only to about eps/sin(Phi) (inherent to Euler angles; as v2).
!>          The conventions are those of Neper's msh file version 2.4, the only version FEPX
!>          reads (ROADMAP D33; older files had the active/passive labels swapped).
!>
!> Descriptor parameters (params(1:num), real(rk)), as in Neper's files:
!>   rodrigues    3  Rodrigues vector, axis tan(angle/2)
!>   euler-bunge  3  phi1, Phi, phi2 [deg]
!>   euler-kocks  3  Psi, Theta, phi [deg]
!>   euler-roe    3  psi, theta, phi [deg]
!>   rotmat       9  g11 g12 g13 g21 g22 g23 g31 g32 g33 (the matrix g, row by row)
!>   axis-angle   4  axis (3 components), angle [deg]
!>   quaternion   4  q0 (scalar part), q1, q2, q3
!> The descriptor matrix g is the matrix of v2's conversions and Neper's `rotmat` descriptor
!> (e.g. g(1,2) = 2(q1 q2 + q0 q3) for a quaternion).
!>
!> Output ranges of descriptor_from_rot_c2s: Euler angles from atan2, in (-180°, 180°] for
!> phi1, phi2, Psi, phi (v2), Phi and Theta in [0°, 180°], Roe psi = phi1 - 90° and
!> phi = phi2 + 90°; quaternion with q0 >= 0 (for q0 = 0, the first nonzero of q1..q3 > 0);
!> axis-angle with angle in [0°, 180°] (axis (0,0,1) for the identity); Rodrigues vector
!> large but finite for a 180° rotation (as v2).
!>
!> Public parameters:
!>   ori_des_rodrigues, ori_des_euler_bunge, ori_des_euler_kocks, ori_des_euler_roe,
!>   ori_des_rotmat, ori_des_axis_angle, ori_des_quaternion   integer  descriptor codes 1..7
!>   ori_conv_passive, ori_conv_active                        integer  convention codes 1, 2
!>   ori_max_params                                           integer  9
!>
!> Public routines:
!>   ori_descriptor_from_name   descriptor code of a Neper descriptor name (0 if unknown)
!>   ori_descriptor_name        Neper name of a descriptor code
!>   ori_num_params             number of parameters of a descriptor
!>   ori_convention_from_name   convention code of "passive" / "active" (0 if unknown)
!>   ori_convention_name        name of a convention code
!>   rot_c2s_from_descriptor    descriptor + convention + parameters -> rot_c2s (validated)
!>   descriptor_from_rot_c2s    rot_c2s -> parameters of a descriptor in a convention
!>   g_from_rodrigues, g_from_euler_bunge, g_from_euler_kocks, g_from_euler_roe,
!>   g_from_axis_angle, g_from_quaternion        descriptor -> g (pure, no validation)
!>   quaternion_from_g          g -> unit quaternion (Shepperd's method, pure)
module orientation_mod

  use kinds_mod, only: rk
  use constants_mod, only: pi, deg_to_rad, rad_to_deg, identity_3
  use error_mod, only: fatal_error, internal_error

  implicit none

  private

  public :: ori_des_rodrigues, ori_des_euler_bunge, ori_des_euler_kocks, ori_des_euler_roe
  public :: ori_des_rotmat, ori_des_axis_angle, ori_des_quaternion
  public :: ori_conv_passive, ori_conv_active, ori_max_params
  public :: ori_descriptor_from_name, ori_descriptor_name, ori_num_params
  public :: ori_convention_from_name, ori_convention_name
  public :: rot_c2s_from_descriptor, descriptor_from_rot_c2s
  public :: g_from_rodrigues, g_from_euler_bunge, g_from_euler_kocks, g_from_euler_roe
  public :: g_from_axis_angle, g_from_quaternion, quaternion_from_g

  ! Descriptor codes.
  integer, parameter :: ori_des_rodrigues = 1
  integer, parameter :: ori_des_euler_bunge = 2
  integer, parameter :: ori_des_euler_kocks = 3
  integer, parameter :: ori_des_euler_roe = 4
  integer, parameter :: ori_des_rotmat = 5
  integer, parameter :: ori_des_axis_angle = 6
  integer, parameter :: ori_des_quaternion = 7

  ! Convention codes.
  integer, parameter :: ori_conv_passive = 1
  integer, parameter :: ori_conv_active = 2

  ! Largest number of parameters of a descriptor (rotmat).
  integer, parameter :: ori_max_params = 9

  ! Neper names and parameter counts, indexed by descriptor code.
  character(len=11), parameter :: descriptor_names(7) = [character(len=11) :: 'rodrigues', &
    & 'euler-bunge', 'euler-kocks', 'euler-roe', 'rotmat', 'axis-angle', 'quaternion']
  integer, parameter :: descriptor_num_params(7) = [3, 3, 3, 3, 9, 4, 4]

  ! v2's `vtiny`: computed sines and cosines below this magnitude are set to zero, so that
  ! angles such as 90° give exact zeros (v2 orientation_conversion_mod.f90).
  real(rk), parameter :: trig_zero_tol = 1.0e-15_rk

  ! Largest |g^T g - I| entry and smallest det(g) accepted for a `rotmat` input (input files
  ! carry about 12 significant digits).
  real(rk), parameter :: rotmat_tol = 1.0e-6_rk

contains

  ! ------------------------------------------------------------------------------------------
  ! Names and codes.
  ! ------------------------------------------------------------------------------------------

  !> Descriptor code of a Neper descriptor name (exact, lowercase), or 0 if unknown.
  !>
  !> Arguments:
  !>   name  (in)  character  descriptor name, e.g. "euler-bunge"
  !>   code  (out, result) integer  descriptor code, 0 if the name is not a descriptor
  !> Locals:
  !>   i     integer  loop index over descriptors
  pure function ori_descriptor_from_name(name) result(code)

    character(len=*), intent(in) :: name
    integer :: code

    integer :: i

    ! Look the name up in the table.
    code = 0
    do i = 1, size(descriptor_names)
      if (trim(name) == trim(descriptor_names(i))) code = i
    end do

  end function ori_descriptor_from_name

  !> Neper name of a descriptor code ("unknown" for an invalid code).
  !>
  !> Arguments:
  !>   code  (in)  integer  descriptor code
  !>   name  (out, result) character(:), allocatable  descriptor name
  pure function ori_descriptor_name(code) result(name)

    integer, intent(in) :: code
    character(len=:), allocatable :: name

    ! Table lookup with a guard.
    if (code >= 1 .and. code <= size(descriptor_names)) then
      name = trim(descriptor_names(code))
    else
      name = 'unknown'
    end if

  end function ori_descriptor_name

  !> Number of parameters of a descriptor (0 for an invalid code).
  !>
  !> Arguments:
  !>   code  (in)  integer  descriptor code
  !>   num   (out, result) integer  number of parameters
  pure function ori_num_params(code) result(num)

    integer, intent(in) :: code
    integer :: num

    ! Table lookup with a guard.
    num = 0
    if (code >= 1 .and. code <= size(descriptor_num_params)) num = descriptor_num_params(code)

  end function ori_num_params

  !> Convention code of a convention name ("passive" or "active"), or 0 if unknown.
  !>
  !> Arguments:
  !>   name  (in)  character  convention name
  !>   code  (out, result) integer  convention code, 0 if unknown
  pure function ori_convention_from_name(name) result(code)

    character(len=*), intent(in) :: name
    integer :: code

    ! Compare with the two names.
    select case (trim(name))
    case ('passive')
      code = ori_conv_passive
    case ('active')
      code = ori_conv_active
    case default
      code = 0
    end select

  end function ori_convention_from_name

  !> Name of a convention code ("unknown" for an invalid code).
  !>
  !> Arguments:
  !>   code  (in)  integer  convention code
  !>   name  (out, result) character(:), allocatable  "passive", "active" or "unknown"
  pure function ori_convention_name(code) result(name)

    integer, intent(in) :: code
    character(len=:), allocatable :: name

    ! Two names and a guard.
    select case (code)
    case (ori_conv_passive)
      name = 'passive'
    case (ori_conv_active)
      name = 'active'
    case default
      name = 'unknown'
    end select

  end function ori_convention_name

  ! ------------------------------------------------------------------------------------------
  ! Descriptor -> rot_c2s.
  ! ------------------------------------------------------------------------------------------

  !> Crystal-to-sample rotation of an orientation given by a descriptor, a convention and its
  !> parameters. The axis of an axis-angle pair and a quaternion are normalized; a `rotmat`
  !> input must be orthogonal with det +1 within rotmat_tol. Invalid input is a fatal error.
  !>
  !> Theory : module header (passive: rot_c2s = g^T; active: rot_c2s = g).
  !>
  !> Arguments:
  !>   descriptor  (in)  integer             descriptor code
  !>   convention  (in)  integer             convention code
  !>   params      (in)  real(rk)(:)         parameters (at least ori_num_params(descriptor))
  !>   rot_c2s     (out, result) real(rk)(3,3)  crystal-to-sample rotation
  !> Locals:
  !>   g           real(rk)(3,3)  descriptor matrix
  !>   norm        real(rk)       norm of an axis or quaternion
  !>   det         real(rk)       determinant of a rotmat input
  !>   detail      character(160) error message
  function rot_c2s_from_descriptor(descriptor, convention, params) result(rot_c2s)

    integer, intent(in) :: descriptor
    integer, intent(in) :: convention
    real(rk), intent(in) :: params(:)
    real(rk) :: rot_c2s(3, 3)

    character(len=*), parameter :: routine = 'rot_c2s_from_descriptor'
    real(rk) :: g(3, 3)
    real(rk) :: norm
    real(rk) :: det
    character(len=160) :: detail

    ! Check the descriptor, the convention and the number of parameters.
    if (ori_num_params(descriptor) == 0) then
      write (detail, '(a,i0)') 'invalid orientation descriptor code ', descriptor
      call internal_error(routine, trim(detail))
    end if
    if (convention /= ori_conv_passive .and. convention /= ori_conv_active) then
      write (detail, '(a,i0)') 'invalid orientation convention code ', convention
      call internal_error(routine, trim(detail))
    end if
    if (size(params) < ori_num_params(descriptor)) then
      write (detail, '(a,i0,a,a,a,i0)') 'expected ', ori_num_params(descriptor), &
        & ' parameters for ', ori_descriptor_name(descriptor), ', got ', size(params)
      call internal_error(routine, trim(detail))
    end if

    ! Descriptor matrix g.
    select case (descriptor)
    case (ori_des_rodrigues)
      g = g_from_rodrigues(params(1:3))
    case (ori_des_euler_bunge)
      g = g_from_euler_bunge(params(1:3))
    case (ori_des_euler_kocks)
      g = g_from_euler_kocks(params(1:3))
    case (ori_des_euler_roe)
      g = g_from_euler_roe(params(1:3))
    case (ori_des_rotmat)
      ! Row-by-row input; check that it is a rotation.
      g = transpose(reshape(params(1:9), [3, 3]))
      det = g(1, 1)*(g(2, 2)*g(3, 3) - g(2, 3)*g(3, 2)) &
        & - g(1, 2)*(g(2, 1)*g(3, 3) - g(2, 3)*g(3, 1)) &
        & + g(1, 3)*(g(2, 1)*g(3, 2) - g(2, 2)*g(3, 1))
      if (maxval(abs(matmul(transpose(g), g) - identity_3)) > rotmat_tol &
        & .or. det < 1.0_rk - rotmat_tol) then
        call fatal_error(routine, 'rotmat orientation is not a rotation matrix (orthogonal, '// &
          & 'det +1) to 1e-6')
      end if
    case (ori_des_axis_angle)
      ! Normalize the axis.
      norm = norm2(params(1:3))
      if (.not. (norm > 0.0_rk)) call fatal_error(routine, 'axis-angle orientation with a zero axis')
      g = g_from_axis_angle(params(1:3)/norm, params(4))
    case (ori_des_quaternion)
      ! Normalize the quaternion.
      norm = norm2(params(1:4))
      if (.not. (norm > 0.0_rk)) call fatal_error(routine, 'zero quaternion orientation')
      g = g_from_quaternion(params(1:4)/norm)
    end select

    ! Convention: passive descriptors give rot_c2s = g^T, active ones rot_c2s = g.
    if (convention == ori_conv_passive) then
      rot_c2s = transpose(g)
    else
      rot_c2s = g
    end if

  end function rot_c2s_from_descriptor

  !> Descriptor matrix of a Rodrigues vector r = n tan(angle/2).
  !>
  !> Source : port of v2 orientation_conversion_mod.f90::rodrigues_to_rotmat (via axis-angle).
  !>
  !> Arguments:
  !>   rod  (in)  real(rk)(3)    Rodrigues vector
  !>   g    (out, result) real(rk)(3,3)  descriptor matrix
  !> Locals:
  !>   angle  real(rk)     rotation angle, 2 atan |r| [rad]
  !>   axis   real(rk)(3)  rotation axis (e1 for a zero vector)
  pure function g_from_rodrigues(rod) result(g)

    real(rk), intent(in) :: rod(3)
    real(rk) :: g(3, 3)

    real(rk) :: angle
    real(rk) :: axis(3)

    ! Axis-angle form of the Rodrigues vector.
    angle = 2.0_rk*atan(norm2(rod))
    if (angle > trig_zero_tol) then
      axis = rod/norm2(rod)
    else
      axis = [1.0_rk, 0.0_rk, 0.0_rk]
    end if

    ! Matrix of the axis-angle pair (angle converted to degrees for the shared routine).
    g = g_from_axis_angle(axis, angle*rad_to_deg)

  end function g_from_rodrigues

  !> Descriptor matrix of an axis-angle pair (unit axis, angle in degrees).
  !>
  !> Source : port of v2 orientation_conversion_mod.f90::axis_angle_to_rotmat.
  !>
  !> Arguments:
  !>   axis   (in)  real(rk)(3)  unit rotation axis
  !>   angle  (in)  real(rk)     rotation angle [deg]
  !>   g      (out, result) real(rk)(3,3)  descriptor matrix
  !> Locals:
  !>   com, som  real(rk)  cosine and sine of the angle (zeroed below trig_zero_tol)
  pure function g_from_axis_angle(axis, angle) result(g)

    real(rk), intent(in) :: axis(3)
    real(rk), intent(in) :: angle
    real(rk) :: g(3, 3)

    real(rk) :: com, som

    ! Cosine and sine, with v2's zeroing of tiny values.
    com = cos(angle*deg_to_rad)
    som = sin(angle*deg_to_rad)
    if (abs(com) < trig_zero_tol) com = 0.0_rk
    if (abs(som) < trig_zero_tol) som = 0.0_rk

    ! v2's matrix, row by row.
    g(1, :) = [(1.0_rk - com)*axis(1)**2 + com, (1.0_rk - com)*axis(1)*axis(2) + som*axis(3), &
      & (1.0_rk - com)*axis(1)*axis(3) - som*axis(2)]
    g(2, :) = [(1.0_rk - com)*axis(1)*axis(2) - som*axis(3), (1.0_rk - com)*axis(2)**2 + com, &
      & (1.0_rk - com)*axis(2)*axis(3) + som*axis(1)]
    g(3, :) = [(1.0_rk - com)*axis(1)*axis(3) + som*axis(2), &
      & (1.0_rk - com)*axis(2)*axis(3) - som*axis(1), (1.0_rk - com)*axis(3)**2 + com]

  end function g_from_axis_angle

  !> Descriptor matrix of Euler angles in Bunge's convention (phi1, Phi, phi2) [deg].
  !>
  !> Source : port of v2 orientation_conversion_mod.f90::euler_bunge_to_rotmat.
  !>
  !> Arguments:
  !>   angles  (in)  real(rk)(3)  phi1, Phi, phi2 [deg]
  !>   g       (out, result) real(rk)(3,3)  descriptor matrix
  !> Locals:
  !>   sps1, cps1, sph, cph, sps2, cps2  real(rk)  sines and cosines (zeroed below trig_zero_tol)
  pure function g_from_euler_bunge(angles) result(g)

    real(rk), intent(in) :: angles(3)
    real(rk) :: g(3, 3)

    real(rk) :: sps1, cps1, sph, cph, sps2, cps2

    ! Sines and cosines, with v2's zeroing of tiny values.
    sps1 = trig_zeroed(sin(angles(1)*deg_to_rad))
    cps1 = trig_zeroed(cos(angles(1)*deg_to_rad))
    sph = trig_zeroed(sin(angles(2)*deg_to_rad))
    cph = trig_zeroed(cos(angles(2)*deg_to_rad))
    sps2 = trig_zeroed(sin(angles(3)*deg_to_rad))
    cps2 = trig_zeroed(cos(angles(3)*deg_to_rad))

    ! v2's matrix, row by row.
    g(1, :) = [cps1*cps2 - sps1*cph*sps2, sps1*cps2 + cps1*cph*sps2, sph*sps2]
    g(2, :) = [-cps1*sps2 - sps1*cph*cps2, -sps1*sps2 + cps1*cph*cps2, sph*cps2]
    g(3, :) = [sps1*sph, -cps1*sph, cph]

  end function g_from_euler_bunge

  !> Descriptor matrix of Euler angles in Kocks' convention (Psi, Theta, phi) [deg].
  !>
  !> Source : port of v2 orientation_conversion_mod.f90::euler_kocks_to_rotmat.
  !>
  !> Arguments:
  !>   angles  (in)  real(rk)(3)  Psi, Theta, phi [deg]
  !>   g       (out, result) real(rk)(3,3)  descriptor matrix
  !> Locals:
  !>   sps, cps, sth, cth, sph, cph  real(rk)  sines and cosines (zeroed below trig_zero_tol)
  pure function g_from_euler_kocks(angles) result(g)

    real(rk), intent(in) :: angles(3)
    real(rk) :: g(3, 3)

    real(rk) :: sps, cps, sth, cth, sph, cph

    ! Sines and cosines, with v2's zeroing of tiny values.
    sps = trig_zeroed(sin(angles(1)*deg_to_rad))
    cps = trig_zeroed(cos(angles(1)*deg_to_rad))
    sth = trig_zeroed(sin(angles(2)*deg_to_rad))
    cth = trig_zeroed(cos(angles(2)*deg_to_rad))
    sph = trig_zeroed(sin(angles(3)*deg_to_rad))
    cph = trig_zeroed(cos(angles(3)*deg_to_rad))

    ! v2's matrix, row by row.
    g(1, :) = [-sps*sph - cps*cph*cth, cps*sph - sps*cph*cth, cph*sth]
    g(2, :) = [cph*sps - sph*cps*cth, -cps*cph - sps*sph*cth, sph*sth]
    g(3, :) = [cps*sth, sps*sth, cth]

  end function g_from_euler_kocks

  !> Descriptor matrix of Euler angles in Roe's convention (psi, theta, phi) [deg], through
  !> Bunge's: phi1 = psi + 90°, Phi = theta, phi2 = phi - 90°.
  !>
  !> Source : Neper orilib ol_er_e_rad; Kocks, Tomé & Wenk (1998), Texture and Anisotropy, Ch. 2.
  !>
  !> Arguments:
  !>   angles  (in)  real(rk)(3)  psi, theta, phi [deg]
  !>   g       (out, result) real(rk)(3,3)  descriptor matrix
  pure function g_from_euler_roe(angles) result(g)

    real(rk), intent(in) :: angles(3)
    real(rk) :: g(3, 3)

    ! Convert to Bunge angles and build the matrix.
    g = g_from_euler_bunge([angles(1) + 90.0_rk, angles(2), angles(3) - 90.0_rk])

  end function g_from_euler_roe

  !> Descriptor matrix of a unit quaternion (q0 scalar part).
  !>
  !> Source : port of v2 orientation_conversion_mod.f90::quat_to_rotmat.
  !>
  !> Arguments:
  !>   quat  (in)  real(rk)(4)    unit quaternion q0, q1, q2, q3
  !>   g     (out, result) real(rk)(3,3)  descriptor matrix
  !> Locals:
  !>   qbar  real(rk)  q0^2 - (q1^2 + q2^2 + q3^2)
  pure function g_from_quaternion(quat) result(g)

    real(rk), intent(in) :: quat(4)
    real(rk) :: g(3, 3)

    real(rk) :: qbar

    ! v2's matrix, row by row.
    qbar = quat(1)**2 - (quat(2)**2 + quat(3)**2 + quat(4)**2)
    g(1, :) = [qbar + 2.0_rk*quat(2)**2, 2.0_rk*(quat(2)*quat(3) + quat(1)*quat(4)), &
      & 2.0_rk*(quat(2)*quat(4) - quat(1)*quat(3))]
    g(2, :) = [2.0_rk*(quat(2)*quat(3) - quat(1)*quat(4)), qbar + 2.0_rk*quat(3)**2, &
      & 2.0_rk*(quat(3)*quat(4) + quat(1)*quat(2))]
    g(3, :) = [2.0_rk*(quat(2)*quat(4) + quat(1)*quat(3)), &
      & 2.0_rk*(quat(3)*quat(4) - quat(1)*quat(2)), qbar + 2.0_rk*quat(4)**2]

  end function g_from_quaternion

  ! ------------------------------------------------------------------------------------------
  ! rot_c2s -> descriptor.
  ! ------------------------------------------------------------------------------------------

  !> Parameters of a descriptor, in a convention, of the orientation rot_c2s (for output in the
  !> descriptor and convention that were read; conventions §2). Invalid codes are fatal errors.
  !>
  !> Theory : module header; output ranges listed there.
  !>
  !> Arguments:
  !>   descriptor  (in)  integer        descriptor code
  !>   convention  (in)  integer        convention code
  !>   rot_c2s     (in)  real(rk)(3,3)  crystal-to-sample rotation
  !>   params      (out, result) real(rk)(ori_max_params)  parameters (unused entries 0)
  !> Locals:
  !>   g           real(rk)(3,3)  descriptor matrix
  !>   quat        real(rk)(4)    unit quaternion of g
  !>   vec_norm    real(rk)       norm of the vector part of quat
  !>   bunge       real(rk)(3)    Bunge angles [deg]
  !>   detail      character(80)  error message
  function descriptor_from_rot_c2s(descriptor, convention, rot_c2s) result(params)

    integer, intent(in) :: descriptor
    integer, intent(in) :: convention
    real(rk), intent(in) :: rot_c2s(3, 3)
    real(rk) :: params(ori_max_params)

    character(len=*), parameter :: routine = 'descriptor_from_rot_c2s'
    real(rk) :: g(3, 3)
    real(rk) :: quat(4)
    real(rk) :: vec_norm
    real(rk) :: bunge(3)
    character(len=80) :: detail

    params = 0.0_rk

    ! Descriptor matrix of the convention (inverse of rot_c2s_from_descriptor).
    select case (convention)
    case (ori_conv_passive)
      g = transpose(rot_c2s)
    case (ori_conv_active)
      g = rot_c2s
    case default
      write (detail, '(a,i0)') 'invalid orientation convention code ', convention
      call internal_error(routine, trim(detail))
    end select

    ! Parameters of the descriptor.
    select case (descriptor)
    case (ori_des_rodrigues)
      ! n tan(angle/2) = q_vec/q0, written as v2 does (large but finite for 180°).
      quat = quaternion_from_g(g)
      vec_norm = norm2(quat(2:4))
      if (vec_norm > trig_zero_tol) then
        params(1:3) = (quat(2:4)/vec_norm)*tan(atan2(vec_norm, quat(1)))
      end if
    case (ori_des_euler_bunge)
      params(1:3) = euler_bunge_from_g(g)
    case (ori_des_euler_kocks)
      params(1:3) = euler_kocks_from_g(g)
    case (ori_des_euler_roe)
      ! Roe from Bunge: psi = phi1 - 90°, theta = Phi, phi = phi2 + 90°.
      bunge = euler_bunge_from_g(g)
      params(1:3) = [bunge(1) - 90.0_rk, bunge(2), bunge(3) + 90.0_rk]
    case (ori_des_rotmat)
      ! g row by row.
      params(1:9) = reshape(transpose(g), [9])
    case (ori_des_axis_angle)
      ! Axis q_vec/|q_vec| and angle 2 atan2(|q_vec|, q0) in [0°, 180°].
      quat = quaternion_from_g(g)
      vec_norm = norm2(quat(2:4))
      if (vec_norm > 0.0_rk) then
        params(1:3) = quat(2:4)/vec_norm
        params(4) = 2.0_rk*atan2(vec_norm, quat(1))*rad_to_deg
      else
        params(1:4) = [0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk]
      end if
    case (ori_des_quaternion)
      params(1:4) = quaternion_from_g(g)
    case default
      write (detail, '(a,i0)') 'invalid orientation descriptor code ', descriptor
      call internal_error(routine, trim(detail))
    end select

  end function descriptor_from_rot_c2s

  !> Unit quaternion (q0, q1, q2, q3) of a descriptor matrix g, by Shepperd's method: the
  !> largest of 4 q0^2 = 1 + tr g and 4 qi^2 = 1 + 2 g(i,i) - tr g gives one component
  !> (magnitude >= 1/2) by a square root; the other three are off-diagonal sums or differences
  !> divided by 4 times it. Sign rule: q0 >= 0; for q0 = 0, the first nonzero of q1..q3 > 0.
  !>
  !> Source : port of v2-fixed orientation_conversion_mod.f90::rotmat_to_quat_shepperd
  !>          (CHANGES.md C18); Shepperd, S.W. (1978), "Quaternion from rotation matrix",
  !>          J. Guidance and Control 1(3):223-224. g as built by g_from_quaternion.
  !>
  !> Arguments:
  !>   g     (in)  real(rk)(3,3)  descriptor matrix (a rotation)
  !>   quat  (out, result) real(rk)(4)  unit quaternion
  !> Locals:
  !>   t     real(rk)(4)  4 q0^2, 4 q1^2, 4 q2^2, 4 q3^2 from the diagonal
  !>   k     integer      index of the largest entry of t
  !>   i     integer      loop index over the vector part
  pure function quaternion_from_g(g) result(quat)

    real(rk), intent(in) :: g(3, 3)
    real(rk) :: quat(4)

    real(rk) :: t(4)
    integer :: k, i

    ! Squared components (times 4) from the diagonal; pick the largest.
    t(1) = 1.0_rk + g(1, 1) + g(2, 2) + g(3, 3)
    t(2) = 1.0_rk + g(1, 1) - g(2, 2) - g(3, 3)
    t(3) = 1.0_rk - g(1, 1) + g(2, 2) - g(3, 3)
    t(4) = 1.0_rk - g(1, 1) - g(2, 2) + g(3, 3)
    k = maxloc(t, 1)

    ! Largest component by a square root, the others from the off-diagonal entries.
    select case (k)
    case (1)
      quat(1) = 0.5_rk*sqrt(t(1))
      quat(2) = (g(2, 3) - g(3, 2))/(4.0_rk*quat(1))
      quat(3) = (g(3, 1) - g(1, 3))/(4.0_rk*quat(1))
      quat(4) = (g(1, 2) - g(2, 1))/(4.0_rk*quat(1))
    case (2)
      quat(2) = 0.5_rk*sqrt(t(2))
      quat(1) = (g(2, 3) - g(3, 2))/(4.0_rk*quat(2))
      quat(3) = (g(1, 2) + g(2, 1))/(4.0_rk*quat(2))
      quat(4) = (g(1, 3) + g(3, 1))/(4.0_rk*quat(2))
    case (3)
      quat(3) = 0.5_rk*sqrt(t(3))
      quat(1) = (g(3, 1) - g(1, 3))/(4.0_rk*quat(3))
      quat(2) = (g(1, 2) + g(2, 1))/(4.0_rk*quat(3))
      quat(4) = (g(2, 3) + g(3, 2))/(4.0_rk*quat(3))
    case default
      quat(4) = 0.5_rk*sqrt(t(4))
      quat(1) = (g(1, 2) - g(2, 1))/(4.0_rk*quat(4))
      quat(2) = (g(1, 3) + g(3, 1))/(4.0_rk*quat(4))
      quat(3) = (g(2, 3) + g(3, 2))/(4.0_rk*quat(4))
    end select

    ! Sign rule (q0 = 0 is tested as "neither < 0 nor > 0", a nonzero component as |qi| > 0).
    if (quat(1) < 0.0_rk) then
      quat = -quat
    else if (.not. (quat(1) > 0.0_rk)) then
      do i = 2, 4
        if (abs(quat(i)) > 0.0_rk) then
          if (quat(i) < 0.0_rk) quat = -quat
          exit
        end if
      end do
    end if

    ! Normalize.
    quat = quat/norm2(quat)

  end function quaternion_from_g

  !> Bunge Euler angles (phi1, Phi, phi2) [deg] of a descriptor matrix g.
  !>
  !> Source : port of v2-fixed orientation_conversion_mod.f90::rotmat_to_euler_bunge (C18:
  !>          Phi = atan2(sqrt(g13^2 + g23^2), g33) instead of acos(g33); the positive scale
  !>          factor 1/sin(Phi) of v2's atan2 arguments cancels and is not computed).
  !>
  !> Arguments:
  !>   g       (in)  real(rk)(3,3)  descriptor matrix
  !>   angles  (out, result) real(rk)(3)  phi1, Phi, phi2 [deg]
  pure function euler_bunge_from_g(g) result(angles)

    real(rk), intent(in) :: g(3, 3)
    real(rk) :: angles(3)

    ! General case, and the degenerate case Phi = 0 or 180° (only phi1 + phi2 or phi1 - phi2
    ! is defined; phi2 = 0 as v2). v2 tests |g33| /= 1; |g33| < 1 is the same test for a
    ! rotation (|g33| <= 1) and also sends a rounded |g33| > 1 to the degenerate case.
    if (abs(g(3, 3)) < 1.0_rk) then
      angles(1) = atan2(g(3, 1), -g(3, 2))
      angles(2) = atan2(sqrt(g(1, 3)**2 + g(2, 3)**2), g(3, 3))
      angles(3) = atan2(g(1, 3), g(2, 3))
    else
      angles(1) = atan2(g(1, 2), g(1, 1))
      angles(2) = 0.5_rk*pi*(1.0_rk - g(3, 3))
      angles(3) = 0.0_rk
    end if
    angles = angles*rad_to_deg

  end function euler_bunge_from_g

  !> Kocks Euler angles (Psi, Theta, phi) [deg] of a descriptor matrix g.
  !>
  !> Source : port of v2-fixed orientation_conversion_mod.f90::rotmat_to_euler_kocks (C18:
  !>          Theta = atan2(sqrt(g13^2 + g23^2), g33) instead of acos(g33); the positive scale
  !>          factor 1/sin(Theta) of v2's atan2 arguments cancels and is not computed; C19:
  !>          sign of the degenerate case Theta = 180°).
  !>
  !> Arguments:
  !>   g       (in)  real(rk)(3,3)  descriptor matrix
  !>   angles  (out, result) real(rk)(3)  Psi, Theta, phi [deg]
  pure function euler_kocks_from_g(g) result(angles)

    real(rk), intent(in) :: g(3, 3)
    real(rk) :: angles(3)

    ! General case, and the degenerate case Theta = 0 or 180° (Psi = 0 as v2; |g33| < 1 as in
    ! euler_bunge_from_g).
    angles(2) = atan2(sqrt(g(1, 3)**2 + g(2, 3)**2), g(3, 3))
    if (abs(g(3, 3)) < 1.0_rk) then
      angles(1) = atan2(g(3, 2), g(3, 1))
      angles(3) = atan2(g(2, 3), g(1, 3))
    else
      ! Theta = 0: g11 = -cos(Psi - phi), g21 = sin(Psi - phi); Theta = 180°: g11 =
      ! cos(Psi + phi), g21 = sin(Psi + phi); with Psi = 0, phi = atan2(-g21 g33, -g11 g33)
      ! (v2-fixed, CHANGES.md C19).
      angles(1) = 0.0_rk
      angles(3) = atan2(-g(2, 1)*g(3, 3), -g(1, 1)*g(3, 3))
    end if
    angles = angles*rad_to_deg

  end function euler_kocks_from_g

  !> A sine or cosine, set to zero below trig_zero_tol (v2's `vtiny` treatment).
  !>
  !> Arguments:
  !>   val     (in)  real(rk)  sine or cosine value
  !>   zeroed  (out, result) real(rk)  val, or 0 if |val| < trig_zero_tol
  pure function trig_zeroed(val) result(zeroed)

    real(rk), intent(in) :: val
    real(rk) :: zeroed

    zeroed = val
    if (abs(val) < trig_zero_tol) zeroed = 0.0_rk

  end function trig_zeroed

end module orientation_mod
