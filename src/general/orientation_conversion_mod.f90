! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module orientation_conversion_mod

! Module containing orientation conversions for input and printing.

! Note: The conversions here assume no convention ("active" or "passive"). The
! default assumption in FEPX is a "passive" convention. If an "active"
! convention is chosen by the user, the internal rotation matrix is transposed
! when converted from it's original orientation convention, and again before
! output to files. FEPX assumes "active" to mean "sample-to-crystal"
! transformation, and "passive" to mean "crystal-to-sample" transformation.

! All conversions are summarized in a paper by Rowenhorst et. al.:
!     doi:10.1088/0965-0393/23/8/083501, 2015.

! Note: The following subroutines rely on other conversions:
! rodrigues_to_rotmat: First converts to axis-angle, then to rotation matrix.
! rotmat_to_axis_angle: First converts to quaternions, then to axis-angle.
! rotmat_to_rodrigues: First converts to quaternions, then to Rodrigues.

! Contains subroutines:
! axis_angle_to_rotmat: Axis-angle (degrees) to rotation matrices.
! euler_bunge_to_rotmat: Euler-Bunge angles (degrees) to rotation matrices.
! euler_kocks_to_rotmat: Euler-Kocks angles (degrees) to rotation matrices.
! rodrigues_to_rotmat: Rodrigues vectors to rotation matrices.
! rotmat_to_axis_angle: Rotation matrices to axis-angle (degrees)
! rotmat_to_euler_bunge: Rotation matrices to Euler-Bunge angles (degrees).
! rotmat_to_euler_kocks: Rotation matrices to Euler-Kocks angles (degrees).
! rotmat_to_rodrigues: Rotation matrices to rodrigues vectors.
! rotmat_to_quat: Rotation matrices to quaternions.
! quat_to_rotmat: Quaternions to rotation matrices.

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use matrix_operations_mod

  implicit none

  real(rk), parameter, private :: pi = 4.0d0*datan(1.0d0)
  real(rk), parameter, private :: pi_over_180 = pi/180.0d0

  public

contains

  !===========================================================================

  !> Convert axis-angle pairs (degrees) to rotation matrices.
  subroutine axis_angle_to_rotmat(axis, angle, r)

    real(rk), intent(in) :: axis(:, :)
    real(rk), intent(in) :: angle(:)
    real(rk), intent(inout) :: r(:, :, :)

    ! com: Cosine of angle
    ! som: Sine of angle
    integer  :: m, i
    real(rk) :: anglerad
    real(rk) :: com
    real(rk) :: som

    !-----------------------------------------------------------------------

    m = size(axis, 2)

    r = 0.0d0

    do i = 1, m
      anglerad = 0.0d0
      com = 0.0d0
      som = 0.0d0

      anglerad = angle(i)*pi_over_180
      com = cos(anglerad)
      som = sin(anglerad)

      ! Check if returned cos/sin values are near machine epsilon
      if (abs(com) .lt. vtiny) com = 0.0d0
      if (abs(som) .lt. vtiny) som = 0.0d0

      r(1, 1, i) = (1.0d0 - com)*(axis(1, i)**2) + com
      r(1, 2, i) = (1.0d0 - com)*axis(1, i)*axis(2, i) + som*axis(3, i)
      r(1, 3, i) = (1.0d0 - com)*axis(1, i)*axis(3, i) - som*axis(2, i)
      r(2, 1, i) = (1.0d0 - com)*axis(1, i)*axis(2, i) - som*axis(3, i)
      r(2, 2, i) = (1.0d0 - com)*(axis(2, i)**2) + com
      r(2, 3, i) = (1.0d0 - com)*axis(2, i)*axis(3, i) + som*axis(1, i)
      r(3, 1, i) = (1.0d0 - com)*axis(1, i)*axis(3, i) + som*axis(2, i)
      r(3, 2, i) = (1.0d0 - com)*axis(2, i)*axis(3, i) - som*axis(1, i)
      r(3, 3, i) = (1.0d0 - com)*(axis(3, i)**2) + com
    end do

  end subroutine axis_angle_to_rotmat

  !===========================================================================

  !> Convert rotation vectors to rotation matrices.
  subroutine rotvec_to_rotmat(rotvec, r)

    ! rotvec: Array of rotation vectors (axis * angle, radians)
    ! r: Array of rotation matrices
    real(rk), intent(in) :: rotvec(:, :)
    real(rk), intent(inout) :: r(:, :, :)

    real(rk), allocatable :: axis(:, :)
    real(rk), allocatable :: angle(:)
    integer :: i, m

    !-----------------------------------------------------------------------

    m = size(rotvec, 2)

    allocate(axis(3, m))
    allocate(angle(m))

    r = 0.0d0
    axis(1, :) = 1.0d0

    call vec3_norm(rotvec, angle)

    do i = 1, m
      if (angle(i) .gt. vtiny) then
        axis(:, i) = rotvec(:, i)/angle(i)
      end if
    end do

    angle = angle/pi_over_180

    call axis_angle_to_rotmat(axis, angle, r)

  end subroutine rotvec_to_rotmat

  !===========================================================================

  !> Same as rotvec_to_mat, but operates on one value only
  subroutine rotvec_to_rotmat_(rotvec, r)

    real(rk), intent(in) :: rotvec(3)
    real(rk), intent(out) :: r(3, 3)

    real(rk) :: rotvec2(3, 1)
    real(rk) :: r2(3, 3, 1)

    rotvec2(:, 1) = rotvec

    call rotvec_to_rotmat(rotvec2, r2)

    r = r2(:, :, 1)

  end subroutine rotvec_to_rotmat_

  !===========================================================================

  !> Convert Euler-Bunge angles (degrees) to rotation matrices
  subroutine euler_bunge_to_rotmat(psi1, phi, psi2, r)

    ! psi1, phi, psi2: Euler-Bunge angles
    ! r: Array of rotation matrices
    real(rk), intent(in) :: psi1(:)
    real(rk), intent(in) :: phi(:)
    real(rk), intent(in) :: psi2(:)
    real(rk), intent(inout) :: r(:, :, :)

    ! Locals:
    ! psi1rad, phirad, psi2rad: Angles in radians
    ! sps1: Sine of psi1
    ! cps1: Cosine of psi1
    ! sph: Sine of phi
    ! cph: Cosine of phi
    ! sps2: Sine of psi2
    ! cps2: Cosine of psi2
    ! j: Looping indices

    real(rk), allocatable :: psi1rad(:)
    real(rk), allocatable :: phirad(:)
    real(rk), allocatable :: psi2rad(:)
    real(rk), allocatable :: sps1(:)
    real(rk), allocatable :: cps1(:)
    real(rk), allocatable :: sph(:)
    real(rk), allocatable :: cph(:)
    real(rk), allocatable :: sps2(:)
    real(rk), allocatable :: cps2(:)
    integer :: j, m

    !-----------------------------------------------------------------------

    m = size(psi1)

    allocate(psi1rad(m))
    allocate(phirad(m))
    allocate(psi2rad(m))
    allocate(sps1(m))
    allocate(cps1(m))
    allocate(sph(m))
    allocate(cph(m))
    allocate(sps2(m))
    allocate(cps2(m))

    r = 0.0d0

    psi1rad = 0.0d0
    phirad = 0.0d0
    psi2rad = 0.0d0
    sps1 = 0.0d0
    cps1 = 0.0d0
    sph = 0.0d0
    cph = 0.0d0
    sps2 = 0.0d0
    cps2 = 0.0d0

    ! Convert from degrees to radians

    psi1rad = psi1*pi_over_180
    phirad = phi*pi_over_180
    psi2rad = psi2*pi_over_180

    sps1 = sin(psi1rad)
    cps1 = cos(psi1rad)
    sph = sin(phirad)
    cph = cos(phirad)
    sps2 = sin(psi2rad)
    cps2 = cos(psi2rad)

    ! Check if returned cos/sin values are near machine epsilon
    do j = 1, m
      if (abs(sps1(j)) .lt. vtiny) sps1(j) = 0.0d0
      if (abs(cps1(j)) .lt. vtiny) cps1(j) = 0.0d0
      if (abs(sph(j)) .lt. vtiny) sph(j) = 0.0d0
      if (abs(cph(j)) .lt. vtiny) cph(j) = 0.0d0
      if (abs(sps2(j)) .lt. vtiny) sps2(j) = 0.0d0
      if (abs(cps2(j)) .lt. vtiny) cps2(j) = 0.0d0
    end do

    r(1, 1, :) = cps1*cps2 - sps1*cph*sps2
    r(1, 2, :) = sps1*cps2 + cps1*cph*sps2
    r(1, 3, :) = sph*sps2
    r(2, 1, :) = -cps1*sps2 - sps1*cph*cps2
    r(2, 2, :) = -sps1*sps2 + cps1*cph*cps2
    r(2, 3, :) = sph*cps2
    r(3, 1, :) = sps1*sph
    r(3, 2, :) = -cps1*sph
    r(3, 3, :) = cph

    return

  end subroutine euler_bunge_to_rotmat

  !===========================================================================

  !> Convert Euler-Kocks angles (degrees) to rotation matrices.
  subroutine euler_kocks_to_rotmat(psi, the, phi, r)

    ! Arguments:
    ! psi, the, phi: Euler-Kocks angles (degrees)
    ! r: Array of rotation matrices

    real(rk), intent(in) :: psi(:)
    real(rk), intent(in) :: the(:)
    real(rk), intent(in) :: phi(:)
    real(rk), intent(inout) :: r(:, :, :)

    ! Locals:
    ! psirad, therad, phirad: Angles in radians
    ! sps: Sine of psi
    ! cps: Cosine of psi
    ! sth: Sine of theta
    ! cth: Cosine of theta
    ! sph: Sine of phi
    ! cph: Cosine of phi
    ! i/j: Looping indices

    integer :: m

    real(rk), allocatable :: psirad(:)
    real(rk), allocatable :: therad(:)
    real(rk), allocatable :: phirad(:)
    real(rk), allocatable :: sps(:)
    real(rk), allocatable :: cps(:)
    real(rk), allocatable :: sth(:)
    real(rk), allocatable :: cth(:)
    real(rk), allocatable :: sph(:)
    real(rk), allocatable :: cph(:)
    integer  :: j

    m = size(psi)

    allocate(psirad(m))
    allocate(therad(m))
    allocate(phirad(m))
    allocate(sps(m))
    allocate(cps(m))
    allocate(sth(m))
    allocate(cth(m))
    allocate(sph(m))
    allocate(cph(m))

    !-----------------------------------------------------------------------

    r = 0.0d0

    psirad = 0.0d0
    therad = 0.0d0
    phirad = 0.0d0
    sps = 0.0d0
    cps = 0.0d0
    sth = 0.0d0
    cth = 0.0d0
    sph = 0.0d0
    cph = 0.0d0

    ! Convert from degrees to radians

    psirad = psi*pi_over_180
    therad = the*pi_over_180
    phirad = phi*pi_over_180

    sps = sin(psirad)
    cps = cos(psirad)
    sth = sin(therad)
    cth = cos(therad)
    sph = sin(phirad)
    cph = cos(phirad)

    ! Check if returned cos/sin values are near machine epsilon
    do j = 1, m
      if (abs(sps(j)) .lt. vtiny) sps(j) = 0.0d0
      if (abs(cps(j)) .lt. vtiny) cps(j) = 0.0d0
      if (abs(sth(j)) .lt. vtiny) sth(j) = 0.0d0
      if (abs(cth(j)) .lt. vtiny) cth(j) = 0.0d0
      if (abs(sph(j)) .lt. vtiny) sph(j) = 0.0d0
      if (abs(cph(j)) .lt. vtiny) cph(j) = 0.0d0
    end do

    r(1, 1, :) = -sps*sph - cps*cph*cth
    r(1, 2, :) = cps*sph - sps*cph*cth
    r(1, 3, :) = cph*sth
    r(2, 1, :) = cph*sps - sph*cps*cth
    r(2, 2, :) = -cps*cph - sps*sph*cth
    r(2, 3, :) = sph*sth
    r(3, 1, :) = cps*sth
    r(3, 2, :) = sps*sth
    r(3, 3, :) = cth

  end subroutine euler_kocks_to_rotmat

  !===========================================================================

  subroutine rodrigues_to_rotmat(rods, r)

    ! Convert Rodrigues vectors to rotation matrices.
    ! Note: first converts Rodrigues vector to axis-angle parameterization,
    ! then to rotation matrix.

    !-----------------------------------------------------------------------

    ! Arguments:
    ! m: Number of elements
    ! rods: Array of Rodrigues vectors
    ! r: Array of rotation matrices

    real(rk), intent(in) :: rods(:, :)
    real(rk), intent(inout) :: r(:, :, :)

    ! Locals:
    ! j: Looping index
    ! angle: Angle of rotation about axis
    ! axis: Axis of rotation
    ! com: Cosine of angle
    ! som: Sine of angle

    integer :: j, m
    real(rk) :: angle
    real(rk) :: axis(3)
    real(rk) :: com
    real(rk) :: som

    m = size(rods, 2)

    !-----------------------------------------------------------------------

    r = 0.0d0

    do j = 1, m
      angle = 0.0d0
      axis = 0.0d0
      com = 0.0d0
      som = 0.0d0

      ! Convert to axis-angle representation

      angle = 2.0d0*atan(norm2(rods(:, j)))

      if (angle .gt. vtiny) then
        axis = rods(:, j)/norm2(rods(:, j))

      else
        axis(1) = 1.0d0
      end if

      com = cos(angle)
      som = sin(angle)

      ! Check if returned cos/sin values are near machine epsilon
      if (abs(com) .lt. vtiny) com = 0.0d0
      if (abs(som) .lt. vtiny) som = 0.0d0

      ! Then from axis-angle to rotation matrix

      r(1, 1, j) = com + (1.0d0 - com)*(axis(1)**2)
      r(1, 2, j) = (1.0d0 - com)*axis(1)*axis(2) + som*axis(3)
      r(1, 3, j) = (1.0d0 - com)*axis(1)*axis(3) - som*axis(2)
      r(2, 1, j) = (1.0d0 - com)*axis(1)*axis(2) - som*axis(3)
      r(2, 2, j) = com + (1.0d0 - com)*(axis(2)**2)
      r(2, 3, j) = (1.0d0 - com)*axis(2)*axis(3) + som*axis(1)
      r(3, 1, j) = (1.0d0 - com)*axis(1)*axis(3) + som*axis(2)
      r(3, 2, j) = (1.0d0 - com)*axis(2)*axis(3) - som*axis(1)
      r(3, 3, j) = com + (1.0d0 - com)*(axis(3)**2)
    end do

  end subroutine rodrigues_to_rotmat

  !===========================================================================

  subroutine rotmat_to_axis_angle(r, aa)

    ! Convert rotation matrices to axis-angle (degrees)

    !-----------------------------------------------------------------------

    ! Arugments:
    ! m: Number of elements
    ! r: Array of rotation matrices
    ! aa: Array of axes (indices 1:3) and angles/rotations (degrees, index 4)
    real(rk), intent(in) :: r(:, :, :)
    real(rk), intent(inout) :: aa(:, :)

    integer :: m
    ! Locals:
    ! q0, q1, q2, q3: Individual components of quaternion
    ! s: Normalization factor for axis
    ! j: Looping index
    real(rk) :: q0
    real(rk) :: q1
    real(rk) :: q2
    real(rk) :: q3
    real(rk) :: s
    integer :: j

    m = size(aa, 2)

    !-----------------------------------------------------------------------

    aa = 0.0d0

    do j = 1, m
      ! Convert first to quaternion

      q0 = 0.0d0
      q1 = 0.0d0
      q2 = 0.0d0
      q3 = 0.0d0

      q0 = 0.5d0*dsqrt(1 + r(1, 1, j) + r(2, 2, j) + &
          & r(3, 3, j))
      q1 = (-0.5d0)*dsqrt(1 + r(1, 1, j) - r(2, 2, j) - &
          & r(3, 3, j))
      q2 = (-0.5d0)*dsqrt(1 - r(1, 1, j) + r(2, 2, j) - &
          & r(3, 3, j))
      q3 = (-0.5d0)*dsqrt(1 - r(1, 1, j) - r(2, 2, j) + &
          & r(3, 3, j))

      if (r(3, 2, j) .lt. r(2, 3, j)) then
        q1 = -q1
      end if

      if (r(1, 3, j) .lt. r(3, 1, j)) then
        q2 = -q2
      end if

      if (r(2, 1, j) .lt. r(1, 2, j)) then
        q3 = -q3
      end if

      ! Then from quaternion to axis-angle

      aa(4, j) = 2.0d0*acos(q0)

      if (aa(4, j) .ne. 0.0d0) then
        if (q0 .ne. 0.0d0) then
          s = dsqrt((q1**2) + (q2**2) + (q3**2))
          aa(1, j) = q1/s
          aa(2, j) = q2/s
          aa(3, j) = q3/s

          if (q0 .lt. 0.0d0) then
            aa(1:3, j) = -aa(1:3, j)
          end if

        else if (q0 .eq. 0.0d0) then
          aa(4, j) = pi
          aa(1, j) = q1
          aa(2, j) = q2
          aa(3, j) = q3
        end if

      else if (aa(4, j) .eq. 0.0d0) then
        aa(1, j) = 0.0d0
        aa(2, j) = 0.0d0
        aa(3, j) = 1.0d0
      end if

      ! Convert from radians to degrees

      aa(4, j) = aa(4, j)/pi_over_180
    end do

  end subroutine rotmat_to_axis_angle

  !===========================================================================

  subroutine rotmat_to_euler_bunge(r, eb)

    ! Convert rotation matrices to Euler-Bunge angles (degrees).

    !-----------------------------------------------------------------------

    ! Arugments:
    ! m: Number of elements
    ! r: Array of rotation matrices
    ! eb: Array of Euler-Bunge angles (in order psi1, phi, psi2)
    real(rk), intent(in) :: r(:, :, :)
    real(rk), intent(inout) :: eb(:, :)
    ! Locals:
    ! xi: Internal scaling factor
    real(rk), allocatable :: xi(:)

    integer :: m

    m = size(eb, 2)

    allocate(xi(m))

    !-----------------------------------------------------------------------

    xi = 0.0d0
    eb = 0.0d0

    xi = 1.0d0/dsqrt(1 - (r(3, 3, :)**2))

    where (abs(r(3, 3, :)) .ne. 1.0d0)
      eb(1, :) = atan2(r(3, 1, :)*xi, -r(3, 2, :)*xi)
      eb(2, :) = acos(r(3, 3, :))
      eb(3, :) = atan2(r(1, 3, :)*xi, r(2, 3, :)*xi)

    else where
      eb(1, :) = atan2(r(1, 2, :), r(1, 1, :))
      eb(2, :) = (pi/2.0d0)*(1.0d0 - r(3, 3, :))
      eb(3, :) = 0.0d0
    end where

    eb(1, :) = eb(1, :)/pi_over_180
    eb(2, :) = eb(2, :)/pi_over_180
    eb(3, :) = eb(3, :)/pi_over_180

  end subroutine rotmat_to_euler_bunge

  !===========================================================================

  subroutine rotmat_to_euler_kocks(r, ek)

    ! Convert rotation matrices to Euler-Kocks angles (degrees).

    !-----------------------------------------------------------------------

    ! Arugments:
    ! m: Number of elements
    ! r: Array of rotation matrices
    ! ek: Array of Euler-Kocks angles (in order psi, the, phi)
    real(rk), intent(in) :: r(:, :, :)
    real(rk), intent(inout) :: ek(:, :)
    ! Locals:
    ! sth: Sine of theta
    real(rk), allocatable :: sth(:)

    integer :: m

    m = size(ek, 2)

    allocate(sth(m))

    !-----------------------------------------------------------------------

    sth = 0.0d0
    ek = 0.0d0

    ek(2, :) = acos(r(3, 3, :))

    where (abs(r(3, 3, :)) .ne. 1.)
      sth = sin(ek(2, :))
      ek(1, :) = atan2(r(3, 2, :)/sth, r(3, 1, :)/sth)
      ek(3, :) = atan2(r(2, 3, :)/sth, r(1, 3, :)/sth)

    else where
      ek(1, :) = 0.0d0
      ek(3, :) = atan2(-r(2, 1, :), -r(1, 1, :))
    end where

    ek(1, :) = ek(1, :)/pi_over_180
    ek(2, :) = ek(2, :)/pi_over_180
    ek(3, :) = ek(3, :)/pi_over_180

  end subroutine rotmat_to_euler_kocks

  !===========================================================================

  !> Convert rotation matrices to Rodrigues vectors.
  !> Note: first converts rotation matrices to quaternions and then to
  !> Rodrigues vectors.
  subroutine rotmat_to_rodrigues(r, rods)

    real(rk), intent(in) :: r(:, :, :)
    real(rk), intent(inout) :: rods(:, :)

    ! q0, q1, q2, q3: Individual components of quaternion
    ! s: Normalization factor for axis
    ! i, j: Looping indices

    real(rk) :: q0
    real(rk) :: q1
    real(rk) :: q2
    real(rk) :: q3
    real(rk) :: s
    real(rk) :: t
    integer :: j, m

    m = size(rods, 2)

    !-----------------------------------------------------------------------

    rods = 0.0d0

    do j = 1, m
      ! Convert first to quaternions

      q0 = 0.0d0
      q1 = 0.0d0
      q2 = 0.0d0
      q3 = 0.0d0

      ! Take abs() of quantity within dsqrt() to avoid NaNs
      q0 = 0.5d0*dsqrt(abs(1 + r(1, 1, j) + r(2, 2, j) + &
          & r(3, 3, j)))
      q1 = (-0.5d0)*dsqrt(abs(1 + r(1, 1, j) - r(2, 2, j) - &
          & r(3, 3, j)))
      q2 = (-0.5d0)*dsqrt(abs(1 - r(1, 1, j) + r(2, 2, j) - &
          & r(3, 3, j)))
      q3 = (-0.5d0)*dsqrt(abs(1 - r(1, 1, j) - r(2, 2, j) + &
          & r(3, 3, j)))

      if (r(3, 2, j) .lt. r(2, 3, j)) then
        q1 = -q1
      end if

      if (r(1, 3, j) .lt. r(3, 1, j)) then
        q2 = -q2
      end if

      if (r(2, 1, j) .lt. r(1, 2, j)) then
        q3 = -q3
      end if

      ! Then from quaternions to Rodrigues

      s = dsqrt((q1**2) + (q2**2) + (q3**2))
      t = tan(acos(q0))

      if (s .le. vtiny) then
        rods(1, j) = 0.0d0
        rods(2, j) = 0.0d0
        rods(3, j) = 0.0d0

      else
        rods(1, j) = (q1/s)*t
        rods(2, j) = (q2/s)*t
        rods(3, j) = (q3/s)*t
      end if
    end do

  end subroutine rotmat_to_rodrigues

  !===========================================================================

  !> Convert rotation matrices to quaternions
  subroutine rotmat_to_quat(r, quat)

    ! r: Array of rotation matrices
    ! quat: Array of quaternions

    real(rk), intent(in) :: r(:, :, :)
    real(rk), intent(inout) :: quat(:, :)

    ! Locals:
    ! q0, q1, q2, q3: Individual components of quaternion
    real(rk) :: q0
    real(rk) :: q1
    real(rk) :: q2
    real(rk) :: q3
    integer :: i, m

    m = size(quat, 2)

    !-----------------------------------------------------------------------

    quat = 0.0d0

    do i = 1, m
      q0 = 0.0d0
      q1 = 0.0d0
      q2 = 0.0d0
      q3 = 0.0d0

      q0 = 0.5d0*dsqrt(1 + r(1, 1, i) + r(2, 2, i) + &
          & r(3, 3, i))
      q1 = (-0.5d0)*dsqrt(1 + r(1, 1, i) - r(2, 2, i) - &
          & r(3, 3, i))
      q2 = (-0.5d0)*dsqrt(1 - r(1, 1, i) + r(2, 2, i) - &
          & r(3, 3, i))
      q3 = (-0.5d0)*dsqrt(1 - r(1, 1, i) - r(2, 2, i) + &
          & r(3, 3, i))

      if (r(3, 2, i) .lt. r(2, 3, i)) then
        q1 = -q1
      end if

      if (r(1, 3, i) .lt. r(3, 1, i)) then
        q2 = -q2
      end if

      if (r(2, 1, i) .lt. r(1, 2, i)) then
        q3 = -q3
      end if

      quat(1, i) = q0
      quat(2, i) = q1
      quat(3, i) = q2
      quat(4, i) = q3
      quat(:, i) = quat(:, i)/norm2(quat(:, i))
    end do

  end subroutine rotmat_to_quat

  !===========================================================================

  !> Convert quaternions to rotation matrices
  subroutine quat_to_rotmat(quat, r)

    ! quat: Array of quaternions
    ! r: Array of rotation matrices
    real(rk), intent(in) :: quat(:, :)
    real(rk), intent(inout) :: r(:, :, :)

    ! qbar: Quaternion magnitude
    ! j: Looping indices
    real(rk) :: qbar
    integer :: j, m

    !-----------------------------------------------------------------------

    m = size(quat, 2)

    r = 0.0d0

    do j = 1, m
      qbar = 0.0d0
      qbar = (quat(1, j)**2) - ((quat(2, j)**2) + &
          & (quat(3, j)**2) + (quat(4, j)**2))

      r(1, 1, j) = qbar + (2*(quat(2, j)**2))
      r(1, 2, j) = 2*((quat(2, j)*quat(3, j)) + (quat(1, j)*quat(4, j)))
      r(1, 3, j) = 2*((quat(2, j)*quat(4, j)) - (quat(1, j)*quat(3, j)))
      r(2, 1, j) = 2*((quat(2, j)*quat(3, j)) - (quat(1, j)*quat(4, j)))
      r(2, 2, j) = qbar + (2*(quat(3, j)**2))
      r(2, 3, j) = 2*((quat(3, j)*quat(4, j)) + (quat(1, j)*quat(2, j)))
      r(3, 1, j) = 2*((quat(2, j)*quat(4, j)) + (quat(1, j)*quat(3, j)))
      r(3, 2, j) = 2*((quat(3, j)*quat(4, j)) - (quat(1, j)*quat(2, j)))
      r(3, 3, j) = qbar + (2*(quat(4, j)**2))
    end do

  end subroutine quat_to_rotmat

end module orientation_conversion_mod
