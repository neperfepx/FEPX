! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Frame changes between the crystal frame and the sample frame.
!>
!> Theory : dev/conventions.md §2. Each material point stores one rotation matrix rot_c2s
!>          (crystal-to-sample): v_s = rot_c2s v_c; its columns are the crystal axes expressed
!>          in the sample frame. A second-order tensor transforms as A_s = rot_c2s A_c rot_c2s^T,
!>          a Mandel 6-vector with [Q6(rot_c2s)] (recoupled_formulation.md Eq. (R4)), a 6x6
!>          Mandel matrix as [Q6][M_c][Q6]^T. The routines below are the only places where a
!>          transpose of rot_c2s may appear (conventions §2: no ad hoc transposes elsewhere).
!>          (The deviatoric 5-vector transforms, [Q5(rot_c2s)] of KD25 Eq. 31, were removed in
!>          step 9.4: nothing used them.)
!> Source : FEPX (ROADMAP step 1.3). Uses tensor_mod for the Mandel operators.
!>
!> All routines are pure functions of one material point. The 3x3 forms (vec3, mat3, mandel,
!> skew) are computed with 3x3 products, which is cheaper than building [Q6] for one vector;
!> callers that rotate many quantities with one rotation may build [Q6] once with
!> tensor_mod::mandel_rotation.
!>
!> Public routines (rot_c2s real(rk)(3,3) is the first argument of each):
!>   transform_crystal_to_sample_vec3    / transform_sample_to_crystal_vec3     3-vector
!>   transform_crystal_to_sample_mat3    / transform_sample_to_crystal_mat3     3x3 tensor
!>   transform_crystal_to_sample_mandel  / transform_sample_to_crystal_mandel   Mandel 6-vector
!>   transform_crystal_to_sample_mat6    / transform_sample_to_crystal_mat6     6x6 Mandel matrix
!>   transform_crystal_to_sample_skew    / transform_sample_to_crystal_skew     skew 3-vector
!> Rotation updates (step 5.2):
!>   rotation_exp_skew   exact exponential of a skew matrix (Rodrigues' formula)
!>   rot_c2s_spin_update the lattice rotation over an increment from the lattice spin (v2)
!>   rot_c2s_spin_update_crystal  the same, from a spin in the crystal frame (step 6.2)
!>   rotation_log_skew   skew 3-vector of the logarithm of a rotation (step 6.5)
!>   rot_c2s_rotation_rate  sample-frame rotation rate of an increment, log(R R_prev^T)/dt
module rotation_mod

  use kinds_mod, only: rk
  use orientation_mod, only: quaternion_from_g
  use tensor_mod, only: mandel_from_matrix, matrix_from_mandel, skew_from_matrix, &
    & matrix_from_skew, mandel_rotation

  implicit none

  private

  public :: transform_crystal_to_sample_vec3, transform_sample_to_crystal_vec3
  public :: transform_crystal_to_sample_mat3, transform_sample_to_crystal_mat3
  public :: transform_crystal_to_sample_mandel, transform_sample_to_crystal_mandel
  public :: transform_crystal_to_sample_mat6, transform_sample_to_crystal_mat6
  public :: transform_crystal_to_sample_skew, transform_sample_to_crystal_skew
  public :: rotation_exp_skew, rotation_log_skew, rot_c2s_spin_update, &
    & rot_c2s_spin_update_crystal, rot_c2s_rotation_rate

  ! Rotation angle below which Rodrigues' coefficients are taken from their Taylor series (the
  ! truncation error, about angle^4/120, is below the round-off of the closed forms there).
  real(rk), parameter :: small_angle = 1.0e-4_rk

contains

  !> 3-vector from the crystal frame to the sample frame, v_s = rot_c2s v_c.
  !>
  !> Arguments:
  !>   rot_c2s  (in)  real(rk)(3,3)  crystal-to-sample rotation
  !>   vec_c    (in)  real(rk)(3)    vector, crystal frame
  !>   vec_s    (out, result) real(rk)(3)  vector, sample frame
  pure function transform_crystal_to_sample_vec3(rot_c2s, vec_c) result(vec_s)

    real(rk), intent(in) :: rot_c2s(3, 3)
    real(rk), intent(in) :: vec_c(3)
    real(rk) :: vec_s(3)

    ! v_s = R v_c.
    vec_s = matmul(rot_c2s, vec_c)

  end function transform_crystal_to_sample_vec3

  !> 3-vector from the sample frame to the crystal frame, v_c = rot_c2s^T v_s.
  !>
  !> Arguments:
  !>   rot_c2s  (in)  real(rk)(3,3)  crystal-to-sample rotation
  !>   vec_s    (in)  real(rk)(3)    vector, sample frame
  !>   vec_c    (out, result) real(rk)(3)  vector, crystal frame
  pure function transform_sample_to_crystal_vec3(rot_c2s, vec_s) result(vec_c)

    real(rk), intent(in) :: rot_c2s(3, 3)
    real(rk), intent(in) :: vec_s(3)
    real(rk) :: vec_c(3)

    ! v_c = R^T v_s (R^T v = v^T R).
    vec_c = matmul(vec_s, rot_c2s)

  end function transform_sample_to_crystal_vec3

  !> Second-order tensor (3x3, any symmetry) from the crystal to the sample frame,
  !> A_s = rot_c2s A_c rot_c2s^T.
  !>
  !> Arguments:
  !>   rot_c2s  (in)  real(rk)(3,3)  crystal-to-sample rotation
  !>   mat_c    (in)  real(rk)(3,3)  tensor, crystal frame
  !>   mat_s    (out, result) real(rk)(3,3)  tensor, sample frame
  pure function transform_crystal_to_sample_mat3(rot_c2s, mat_c) result(mat_s)

    real(rk), intent(in) :: rot_c2s(3, 3)
    real(rk), intent(in) :: mat_c(3, 3)
    real(rk) :: mat_s(3, 3)

    ! A_s = R A_c R^T.
    mat_s = matmul(rot_c2s, matmul(mat_c, transpose(rot_c2s)))

  end function transform_crystal_to_sample_mat3

  !> Second-order tensor (3x3, any symmetry) from the sample to the crystal frame,
  !> A_c = rot_c2s^T A_s rot_c2s.
  !>
  !> Arguments:
  !>   rot_c2s  (in)  real(rk)(3,3)  crystal-to-sample rotation
  !>   mat_s    (in)  real(rk)(3,3)  tensor, sample frame
  !>   mat_c    (out, result) real(rk)(3,3)  tensor, crystal frame
  pure function transform_sample_to_crystal_mat3(rot_c2s, mat_s) result(mat_c)

    real(rk), intent(in) :: rot_c2s(3, 3)
    real(rk), intent(in) :: mat_s(3, 3)
    real(rk) :: mat_c(3, 3)

    ! A_c = R^T A_s R.
    mat_c = matmul(transpose(rot_c2s), matmul(mat_s, rot_c2s))

  end function transform_sample_to_crystal_mat3

  !> Mandel 6-vector from the crystal to the sample frame, {a}_s = [Q6(rot_c2s)] {a}_c.
  !>
  !> Arguments:
  !>   rot_c2s  (in)  real(rk)(3,3)  crystal-to-sample rotation
  !>   vec_c    (in)  real(rk)(6)    Mandel vector, crystal frame
  !>   vec_s    (out, result) real(rk)(6)  Mandel vector, sample frame
  pure function transform_crystal_to_sample_mandel(rot_c2s, vec_c) result(vec_s)

    real(rk), intent(in) :: rot_c2s(3, 3)
    real(rk), intent(in) :: vec_c(6)
    real(rk) :: vec_s(6)

    ! Through the 3x3 form: {R A_c R^T}.
    vec_s = mandel_from_matrix(transform_crystal_to_sample_mat3(rot_c2s, &
      & matrix_from_mandel(vec_c)))

  end function transform_crystal_to_sample_mandel

  !> Mandel 6-vector from the sample to the crystal frame, {a}_c = [Q6(rot_c2s)]^T {a}_s.
  !>
  !> Arguments:
  !>   rot_c2s  (in)  real(rk)(3,3)  crystal-to-sample rotation
  !>   vec_s    (in)  real(rk)(6)    Mandel vector, sample frame
  !>   vec_c    (out, result) real(rk)(6)  Mandel vector, crystal frame
  pure function transform_sample_to_crystal_mandel(rot_c2s, vec_s) result(vec_c)

    real(rk), intent(in) :: rot_c2s(3, 3)
    real(rk), intent(in) :: vec_s(6)
    real(rk) :: vec_c(6)

    ! Through the 3x3 form: {R^T A_s R}.
    vec_c = mandel_from_matrix(transform_sample_to_crystal_mat3(rot_c2s, &
      & matrix_from_mandel(vec_s)))

  end function transform_sample_to_crystal_mandel

  !> 6x6 Mandel matrix (e.g. stiffness) from the crystal to the sample frame,
  !> [M_s] = [Q6] [M_c] [Q6]^T, [Q6] = [Q6(rot_c2s)] (Eq. (R4)).
  !>
  !> Arguments:
  !>   rot_c2s  (in)  real(rk)(3,3)  crystal-to-sample rotation
  !>   mat_c    (in)  real(rk)(6,6)  Mandel matrix, crystal frame
  !>   mat_s    (out, result) real(rk)(6,6)  Mandel matrix, sample frame
  !> Locals:
  !>   q6       real(rk)(6,6)  [Q6(rot_c2s)]
  pure function transform_crystal_to_sample_mat6(rot_c2s, mat_c) result(mat_s)

    real(rk), intent(in) :: rot_c2s(3, 3)
    real(rk), intent(in) :: mat_c(6, 6)
    real(rk) :: mat_s(6, 6)

    real(rk) :: q6(6, 6)

    ! [Q6] [M_c] [Q6]^T.
    q6 = mandel_rotation(rot_c2s)
    mat_s = matmul(q6, matmul(mat_c, transpose(q6)))

  end function transform_crystal_to_sample_mat6

  !> 6x6 Mandel matrix from the sample to the crystal frame, [M_c] = [Q6]^T [M_s] [Q6].
  !>
  !> Arguments:
  !>   rot_c2s  (in)  real(rk)(3,3)  crystal-to-sample rotation
  !>   mat_s    (in)  real(rk)(6,6)  Mandel matrix, sample frame
  !>   mat_c    (out, result) real(rk)(6,6)  Mandel matrix, crystal frame
  !> Locals:
  !>   q6       real(rk)(6,6)  [Q6(rot_c2s)]
  pure function transform_sample_to_crystal_mat6(rot_c2s, mat_s) result(mat_c)

    real(rk), intent(in) :: rot_c2s(3, 3)
    real(rk), intent(in) :: mat_s(6, 6)
    real(rk) :: mat_c(6, 6)

    real(rk) :: q6(6, 6)

    ! [Q6]^T [M_s] [Q6].
    q6 = mandel_rotation(rot_c2s)
    mat_c = matmul(transpose(q6), matmul(mat_s, q6))

  end function transform_sample_to_crystal_mat6

  !> Skew 3-vector {w23, w13, w12} (conventions §1.4) from the crystal to the sample frame:
  !> W_s = rot_c2s W_c rot_c2s^T.
  !>
  !> Arguments:
  !>   rot_c2s  (in)  real(rk)(3,3)  crystal-to-sample rotation
  !>   skw_c    (in)  real(rk)(3)    skew 3-vector, crystal frame
  !>   skw_s    (out, result) real(rk)(3)  skew 3-vector, sample frame
  pure function transform_crystal_to_sample_skew(rot_c2s, skw_c) result(skw_s)

    real(rk), intent(in) :: rot_c2s(3, 3)
    real(rk), intent(in) :: skw_c(3)
    real(rk) :: skw_s(3)

    ! Through the 3x3 skew matrix.
    skw_s = skew_from_matrix(transform_crystal_to_sample_mat3(rot_c2s, matrix_from_skew(skw_c)))

  end function transform_crystal_to_sample_skew

  !> Skew 3-vector from the sample to the crystal frame: W_c = rot_c2s^T W_s rot_c2s.
  !>
  !> Arguments:
  !>   rot_c2s  (in)  real(rk)(3,3)  crystal-to-sample rotation
  !>   skw_s    (in)  real(rk)(3)    skew 3-vector, sample frame
  !>   skw_c    (out, result) real(rk)(3)  skew 3-vector, crystal frame
  pure function transform_sample_to_crystal_skew(rot_c2s, skw_s) result(skw_c)

    real(rk), intent(in) :: rot_c2s(3, 3)
    real(rk), intent(in) :: skw_s(3)
    real(rk) :: skw_c(3)

    ! Through the 3x3 skew matrix.
    skw_c = skew_from_matrix(transform_sample_to_crystal_mat3(rot_c2s, matrix_from_skew(skw_s)))

  end function transform_sample_to_crystal_skew

  !> Exact exponential of the skew matrix W of a skew 3-vector {w23, w13, w12}
  !> (conventions §1.4): the rotation by the angle theta = |omega| about omega/theta, omega the
  !> axial vector of W (W v = omega x v; omega = (-w23, w13, -w12)), in active form.
  !>
  !> Theory : Rodrigues' formula, exp(W) = I + (sin theta/theta) W + ((1 - cos theta)/theta^2) W^2,
  !>          the closed form of the axis-angle rotation matrix: the same rotation as v2's
  !>          orientation_conversion_mod.f90::rotvec_to_rotmat (axis_angle_to_rotmat of the
  !>          rotation vector, cos/sin form). Below small_angle the two coefficients are
  !>          1 - theta^2/6 and 1/2 - theta^2/24 (their Taylor series).
  !> Source : port of v2 rotvec_to_rotmat as used by solve_evp_rstar_mod2.f90::update_rstar
  !>          (ROADMAP step 5.2), written with the skew matrix (no axis normalization).
  !>
  !> Arguments:
  !>   skw  (in)  real(rk)(3)    skew 3-vector of W (e.g. a spin times a time increment) [-]
  !>   rot  (out, result) real(rk)(3,3)  exp(W), a proper rotation [-]
  !> Locals:
  !>   w_mat   real(rk)(3,3)  W
  !>   theta   real(rk)       rotation angle |omega| = |skw| [rad]
  !>   c1, c2  real(rk)       sin(theta)/theta and (1 - cos(theta))/theta^2
  !>   i       integer        diagonal index
  pure function rotation_exp_skew(skw) result(rot)

    real(rk), intent(in) :: skw(3)
    real(rk) :: rot(3, 3)

    real(rk) :: w_mat(3, 3), theta, c1, c2
    integer :: i

    ! W and its rotation angle.
    w_mat = matrix_from_skew(skw)
    theta = norm2(skw)

    ! Rodrigues' coefficients (series for small angles).
    if (theta < small_angle) then
      c1 = 1.0_rk - theta**2/6.0_rk
      c2 = 0.5_rk - theta**2/24.0_rk
    else
      c1 = sin(theta)/theta
      c2 = (1.0_rk - cos(theta))/theta**2
    end if

    ! I + c1 W + c2 W^2.
    rot = c1*w_mat + c2*matmul(w_mat, w_mat)
    do i = 1, 3
      rot(i, i) = rot(i, i) + 1.0_rk
    end do

  end function rotation_exp_skew

  !> Crystal-to-sample rotation at the end of an increment from the lattice spin: the spin is
  !> expressed in the crystal frame of the current estimate rot_c2s_cur, the increment
  !> rotation is d_rot = exp(W_c dt), and rot_c2s = rot_c2s_prev d_rot (equal to
  !> exp(W_s dt) rot_c2s_prev when rot_c2s_cur = rot_c2s_prev).
  !>
  !> Theory : recoupled_formulation.md §7 (v2 conventions: ori_new = ori_prev d_rstar); with no
  !>          slip the lattice spin is w^p of (R2).
  !> Source : port of v2 solve_evp_rstar_mod2.f90::update_rstar (rotation vector
  !>          (w^p - w_slip) dt in the lattice frame of the current orientation, lattice_spin
  !>          with results%ori) and ::update_ori (ori = ori_prev d_rstar).
  !>
  !> Arguments:
  !>   rot_c2s_prev  (in)  real(rk)(3,3)  rotation at the start of the increment [-]
  !>   rot_c2s_cur   (in)  real(rk)(3,3)  current estimate of the end rotation (gives the
  !>                                      crystal frame of the spin) [-]
  !>   spin_s        (in)  real(rk)(3)    lattice spin {w23, w13, w12}, sample frame [1/s]
  !>   dt            (in)  real(rk)       time increment [s]
  !>   rot_c2s       (out, result) real(rk)(3,3)  rotation at the end of the increment [-]
  !> Locals:
  !>   d_rot         real(rk)(3,3)  increment rotation exp(W_c dt), crystal frame [-]
  pure function rot_c2s_spin_update(rot_c2s_prev, rot_c2s_cur, spin_s, dt) result(rot_c2s)

    real(rk), intent(in) :: rot_c2s_prev(3, 3)
    real(rk), intent(in) :: rot_c2s_cur(3, 3)
    real(rk), intent(in) :: spin_s(3)
    real(rk), intent(in) :: dt
    real(rk) :: rot_c2s(3, 3)

    real(rk) :: d_rot(3, 3)

    ! Increment rotation in the crystal frame (a local first: gfortran -O3 warns spuriously on
    ! matmul of a function result), applied to the start rotation (v2 order).
    d_rot = rotation_exp_skew(dt*transform_sample_to_crystal_skew(rot_c2s_cur, spin_s))
    rot_c2s = matmul(rot_c2s_prev, d_rot)

  end function rot_c2s_spin_update

  !> Crystal-to-sample rotation at the end of an increment from a lattice spin given in the
  !> crystal frame of the current estimate: rot_c2s = rot_c2s_prev exp(W_c dt) (the crystal-frame
  !> form of rot_c2s_spin_update, for spins computed in that frame, e.g. w^p minus the slip
  !> spin of crystal plasticity).
  !>
  !> Theory : as rot_c2s_spin_update.
  !> Source : port of v2 solve_evp_rstar_mod2.f90::update_rstar and ::update_ori (rotation
  !>          vector (w^p - w_slip) dt in the lattice frame; ori = ori_prev d_rstar).
  !>
  !> Arguments:
  !>   rot_c2s_prev  (in)  real(rk)(3,3)  rotation at the start of the increment [-]
  !>   spin_c        (in)  real(rk)(3)    lattice spin {w23, w13, w12}, crystal frame of the
  !>                                      current estimate [1/s]
  !>   dt            (in)  real(rk)       time increment [s]
  !>   rot_c2s       (out, result) real(rk)(3,3)  rotation at the end of the increment [-]
  !> Locals:
  !>   d_rot         real(rk)(3,3)  increment rotation exp(W_c dt), crystal frame [-]
  pure function rot_c2s_spin_update_crystal(rot_c2s_prev, spin_c, dt) result(rot_c2s)

    real(rk), intent(in) :: rot_c2s_prev(3, 3)
    real(rk), intent(in) :: spin_c(3)
    real(rk), intent(in) :: dt
    real(rk) :: rot_c2s(3, 3)

    real(rk) :: d_rot(3, 3)

    ! Increment rotation (a local first, as above), applied to the start rotation.
    d_rot = rotation_exp_skew(dt*spin_c)
    rot_c2s = matmul(rot_c2s_prev, d_rot)

  end function rot_c2s_spin_update_crystal

  !> Skew 3-vector {w23, w13, w12} of the logarithm of a rotation: the inverse of
  !> rotation_exp_skew for rotation angles in [0, pi) (at pi either sign of the axis).
  !> rotation_exp_skew(rotation_log_skew(rot)) = rot for every rotation.
  !>
  !> Theory : the axis-angle form of a rotation (inverse of Rodrigues' formula, see
  !>          rotation_exp_skew): the quaternion (q0, q_vec) of rot^T (quaternion_from_g, which
  !>          returns the quaternion of the passive matrix it is given; q0 >= 0) has
  !>          q_vec = sin(theta/2) n for the active rotation rot by theta about n (checked by
  !>          the round trips of test_rotation), so the axial vector is
  !>          omega = (theta/|q_vec|) q_vec with theta = 2 atan2(|q_vec|, q0), and
  !>          the skew 3-vector is (-omega1, omega2, -omega3) (conventions §1.4).
  !> Source : port of v2-fixed orientation_conversion_mod.f90::rotmat_to_axis_angle (C18:
  !>          Shepperd's quaternion, angle by atan2; accurate for small angles and near pi), as
  !>          used for the rotation-rate outputs (write_res_mod1.f90 rotrate: rotation vector of
  !>          d_rstar). Below small_angle, theta/|q_vec| = 2 (1 + theta^2/24) (Taylor series of
  !>          theta/sin(theta/2)).
  !>
  !> Arguments:
  !>   rot  (in)  real(rk)(3,3)  rotation (proper orthogonal) [-]
  !>   skw  (out, result) real(rk)(3)  skew 3-vector of log(rot) [-]
  !> Locals:
  !>   rot_t   real(rk)(3,3)  rot^T
  !>   quat    real(rk)(4)  quaternion (q0, q_vec) of rot^T, q0 >= 0
  !>   s       real(rk)     |q_vec| = sin(theta/2)
  !>   theta   real(rk)     rotation angle [rad]
  !>   fac     real(rk)     theta/|q_vec|
  !>   omega   real(rk)(3)  axial vector of log(rot) [-]
  pure function rotation_log_skew(rot) result(skw)

    real(rk), intent(in) :: rot(3, 3)
    real(rk) :: skw(3)

    real(rk) :: rot_t(3, 3), quat(4), s, theta, fac, omega(3)

    ! Quaternion of rot^T (quaternion_from_g takes a passive matrix; the transpose is the
    ! passive form of the active rotation rot; a local avoids an array temporary).
    rot_t = transpose(rot)
    quat = quaternion_from_g(rot_t)

    ! Angle and the factor theta/|q_vec| (series for small angles).
    s = norm2(quat(2:4))
    theta = 2.0_rk*atan2(s, quat(1))
    if (theta < small_angle) then
      fac = 2.0_rk*(1.0_rk + theta**2/24.0_rk)
    else
      fac = theta/s
    end if

    ! Axial vector, then the skew 3-vector {w23, w13, w12} = (-omega1, omega2, -omega3).
    omega = fac*quat(2:4)
    skw = [-omega(1), omega(2), -omega(3)]

  end function rotation_log_skew

  !> Lattice rotation rate of an increment in the sample frame: the skew 3-vector of
  !> log(rot_c2s rot_c2s_prev^T) divided by dt. For the updates of this module
  !> (rot_c2s = rot_c2s_prev exp(W_c dt)) it is the spin W_c in the sample frame (R W_c R^T, with
  !> either rotation).
  !>
  !> Theory : rot_c2s rot_c2s_prev^T is the rotation that takes the start lattice to the end
  !>          lattice, in the sample frame; its rotation vector is that of rot_c2s_prev^T rot_c2s
  !>          (v2's d_rstar, crystal frame) transformed to the sample frame.
  !> Source : port of v2-fixed write_res_mod1.f90 rotrate (rotation vector of d_rstar over
  !>          dtime, transformed to the sample frame; CHANGES.md C41: the increment's dtime).
  !>
  !> Arguments:
  !>   rot_c2s_prev  (in)  real(rk)(3,3)  rotation at the start of the increment [-]
  !>   rot_c2s       (in)  real(rk)(3,3)  rotation at the end of the increment [-]
  !>   dt            (in)  real(rk)       time increment [s]
  !>   rate_s        (out, result) real(rk)(3)  rotation rate {w23, w13, w12}, sample [1/s]
  !> Locals:
  !>   rot_incr      real(rk)(3,3)  rot_c2s rot_c2s_prev^T [-]
  pure function rot_c2s_rotation_rate(rot_c2s_prev, rot_c2s, dt) result(rate_s)

    real(rk), intent(in) :: rot_c2s_prev(3, 3)
    real(rk), intent(in) :: rot_c2s(3, 3)
    real(rk), intent(in) :: dt
    real(rk) :: rate_s(3)

    real(rk) :: rot_incr(3, 3)

    ! Rotation of the increment, then its logarithm over dt.
    rot_incr = matmul(rot_c2s, transpose(rot_c2s_prev))
    rate_s = rotation_log_skew(rot_incr)/dt

  end function rot_c2s_rotation_rate

end module rotation_mod
