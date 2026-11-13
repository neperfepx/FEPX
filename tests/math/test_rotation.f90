! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of rotation_mod (crystal <-> sample transforms) and a physical check of the
!> orientation convention with an HCP crystal.
!>
!> Source : FEPX (ROADMAP step 1.3; conventions.md §2).
!>
!> Checks, for 500 random rotations rot_c2s and random inputs (entries in [-1, 1)):
!>   - every crystal_to_sample transform against explicit algebra: R v, R A R^T, {R A R^T},
!>     [Q6][M][Q6]^T, R W R^T, and every sample_to_crystal transform as
!>     its inverse (round trip);
!>   - invariants: A_s:B_s = A_c:B_c, [M_s]{a_s} = {(M_c a_c)_s} (a material law commutes with
!>     the frame change), v_s.u_s = v_c.u_c;
!>   - the columns of rot_c2s are the crystal axes in the sample frame:
!>     transform_crystal_to_sample_vec3(rot_c2s, e_i) = rot_c2s(:, i).
!> Physical HCP check (material point; the finite-element version follows in Stage 4): an HCP
!> crystal (Ti-like constants, MPa) with the Voigt compliances of a transversely isotropic
!> crystal (closed form, checked against [C]^-1). Its uniaxial modulus along a sample direction
!> n is 1/({n n}^T [S_s] {n n}), [S_s] the sample-frame Mandel compliance, and must equal the
!> analytic directional modulus 1/E = S11 (1 - c^2)^2 + S33 c^4 + (2 S13 + S44) c^2 (1 - c^2),
!> c = cos(angle between n and the c-axis). Orientations:
!>   - axis-angle (1,1,1)/sqrt3, 120° (maps x -> y -> z -> x), passive: c-axis along X_s, so
!>     E_x = 1/S33 and E_y = E_z = 1/S11; active: c-axis along Y_s (the inverse rotation), so
!>     E_y = 1/S33. This distinguishes the conventions (the inverse rotation moves the c-axis
!>     to a different sample axis);
!>   - axis-angle (0,1,0), 30°, passive: c-axis (sin 30°, 0, cos 30°): E along X_s, Y_s, Z_s
!>     against the formula with c = sin 30°, 0, cos 30°.
!> Rotation updates (step 5.2): rotation_exp_skew against the Taylor series of exp(W) (40
!> terms), orthogonality, det = 1, the axis omega fixed (W v = omega x v), agreement with v2's
!> axis-angle matrix (orientation_mod::g_from_axis_angle, passive, so exp(W) = g^T), and
!> rotation angles on both sides of the small-angle switch (1e-4) and down to 1e-12;
!> rot_c2s_spin_update with rot_c2s_cur = rot_c2s_prev equals exp(W_s dt) rot_c2s_prev;
!> rot_c2s_spin_update_crystal equals rot_c2s_prev exp(W_c dt) and the sample-frame form with the
!> spin transformed by the current estimate (step 6.2).
!> rotation_log_skew (step 6.5): log(exp(W)) = W for the random skew vectors above (angles up
!> to 3 rad, both sides of the small-angle switch), exp(log(R)) = R for random rotations and for
!> rotations within 1e-6 and 1e-10 rad of pi, and log(I) = 0; rot_c2s_rotation_rate of
!> rot_c2s_spin_update_crystal(R_prev, w_c, dt) is R_prev w_c in the sample frame (round-off
!> 1e-15 absolute from forming R R_prev^T).
!> Tolerances: 1e-13 (transforms, inputs of order 1), 1e-12 relative (moduli).
!>
!> Locals:
!>   num_samples   integer         number of random samples
!>   tol           real(rk)        tolerance of the transform checks
!>   i_sample, i   integer         sample index, loop index
!>   rot           real(rk)(3,3)   random rot_c2s
!>   vec_c, vec_u  real(rk)(3)     random 3-vectors, crystal frame
!>   mat_c, mat_b  real(rk)(3,3)   random symmetric tensors, crystal frame
!>   mat_g         real(rk)(3,3)   random general tensor, crystal frame (later a rotation)
!>   mat_r         real(rk)(3,3)   random rotation (current estimate of the spin-update tests)
!>   mat_w         real(rk)(3,3)   random skew tensor, crystal frame
!>   stiff_c       real(rk)(6,6)   random symmetric Mandel matrix, crystal frame
!>   vec_a, vec_b  real(rk)(6)     Mandel vectors of mat_c, mat_b
!>   q6            real(rk)(6,6)   [Q6(rot)]
!>   eye3          real(rk)(3,3)   identity
!>   comp_c        real(rk)(6,6)   HCP Mandel compliance, crystal frame
!>   stiff_hcp     real(rk)(6,6)   HCP Mandel stiffness, crystal frame
!>   s11, s12, s13, s33, s44   real(rk)  HCP Voigt compliances [1/MPa]
!>   skw, omega    real(rk)(3)     random skew 3-vector, its axial vector
!>   theta         real(rk)        rotation angle [rad]
!>   exp_series, term  real(rk)(3,3)  Taylor sum of exp(W) and its current term
!>   k             integer         series term
!>   scales        real(rk)(8)     magnitudes of the random skew vectors
program test_rotation

  use kinds_mod, only: rk
  use constants_mod, only: identity_3, deg_to_rad
  use tensor_mod, only: mandel_from_matrix, matrix_from_mandel, skew_from_matrix, &
    & matrix_from_skew, mandel_rotation
  use rotation_mod, only: transform_crystal_to_sample_vec3, transform_sample_to_crystal_vec3, &
    & transform_crystal_to_sample_mat3, transform_sample_to_crystal_mat3, &
    & transform_crystal_to_sample_mandel, transform_sample_to_crystal_mandel, &
    & transform_crystal_to_sample_mat6, transform_sample_to_crystal_mat6, &
    & transform_crystal_to_sample_skew, transform_sample_to_crystal_skew, &
    & rotation_exp_skew, &
    & rotation_log_skew, rot_c2s_spin_update, rot_c2s_spin_update_crystal, rot_c2s_rotation_rate
  use orientation_mod, only: rot_c2s_from_descriptor, ori_des_axis_angle, ori_conv_passive, &
    & ori_conv_active, g_from_axis_angle
  use tensor_test_utils_mod, only: seed_random, random_matrix, random_symmetric, random_skew, &
    & random_rotation
  use assert_mod, only: assert_close, test_summary

  implicit none

  integer, parameter :: num_samples = 500
  real(rk), parameter :: tol = 1.0e-13_rk

  integer :: i_sample, i
  real(rk) :: rot(3, 3), vec_c(3), vec_u(3)
  real(rk) :: mat_c(3, 3), mat_b(3, 3), mat_g(3, 3), mat_w(3, 3), mat_r(3, 3)
  real(rk) :: stiff_c(6, 6), vec_a(6), vec_b(6), q6(6, 6), eye3(3, 3)
  real(rk) :: comp_c(6, 6), stiff_hcp(6, 6)
  real(rk) :: s11, s12, s13, s33, s44
  real(rk) :: skw(3), omega(3), theta, exp_series(3, 3), term(3, 3)
  integer :: k
  real(rk), parameter :: scales(8) = [1.0e-12_rk, 1.0e-6_rk, 0.99e-4_rk, 1.01e-4_rk, 1.0e-2_rk, &
    & 0.3_rk, 1.0_rk, 3.0_rk]

  call seed_random(20261008)
  eye3 = identity_3

  ! ------------------------------------------------------------------------------------------
  ! Random property tests.
  ! ------------------------------------------------------------------------------------------

  do i_sample = 1, num_samples

    rot = random_rotation()
    q6 = mandel_rotation(rot)
    call random_number(vec_c)
    call random_number(vec_u)
    mat_c = random_symmetric()
    mat_b = random_symmetric()
    mat_g = random_matrix()
    mat_w = random_skew()
    call random_number(stiff_c)
    stiff_c = stiff_c + transpose(stiff_c) - 1.0_rk
    vec_a = mandel_from_matrix(mat_c)
    vec_b = mandel_from_matrix(mat_b)

    ! Columns of rot_c2s are the crystal axes in the sample frame.
    do i = 1, 3
      call assert_close(transform_crystal_to_sample_vec3(rot, eye3(:, i)), rot(:, i), &
        & 'columns_are_crystal_axes', abs_tol=tol)
    end do

    ! 3-vectors.
    call assert_close(transform_crystal_to_sample_vec3(rot, vec_c), matmul(rot, vec_c), &
      & 'vec3_c2s', abs_tol=tol)
    call assert_close(transform_sample_to_crystal_vec3(rot, &
      & transform_crystal_to_sample_vec3(rot, vec_c)), vec_c, 'vec3_round_trip', abs_tol=tol)
    call assert_close(dot_product(transform_crystal_to_sample_vec3(rot, vec_c), &
      & transform_crystal_to_sample_vec3(rot, vec_u)), dot_product(vec_c, vec_u), &
      & 'vec3_dot_invariant', abs_tol=tol)

    ! General 3x3 tensors.
    call assert_close(transform_crystal_to_sample_mat3(rot, mat_g), &
      & matmul(matmul(rot, mat_g), transpose(rot)), 'mat3_c2s', abs_tol=tol)
    call assert_close(transform_sample_to_crystal_mat3(rot, &
      & transform_crystal_to_sample_mat3(rot, mat_g)), mat_g, 'mat3_round_trip', abs_tol=tol)

    ! Mandel 6-vectors: 3x3 algebra and [Q6]; round trip; inner product invariant.
    call assert_close(transform_crystal_to_sample_mandel(rot, vec_a), &
      & mandel_from_matrix(matmul(matmul(rot, mat_c), transpose(rot))), 'mandel_c2s', &
      & abs_tol=tol)
    call assert_close(transform_crystal_to_sample_mandel(rot, vec_a), matmul(q6, vec_a), &
      & 'mandel_c2s_q6', abs_tol=tol)
    call assert_close(transform_sample_to_crystal_mandel(rot, &
      & transform_crystal_to_sample_mandel(rot, vec_a)), vec_a, 'mandel_round_trip', abs_tol=tol)
    call assert_close(dot_product(transform_crystal_to_sample_mandel(rot, vec_a), &
      & transform_crystal_to_sample_mandel(rot, vec_b)), dot_product(vec_a, vec_b), &
      & 'mandel_inner_invariant', abs_tol=tol)

    ! 6x6 matrices: [Q6][M][Q6]^T; round trip; a material law commutes with the frame change.
    call assert_close(transform_crystal_to_sample_mat6(rot, stiff_c), &
      & matmul(q6, matmul(stiff_c, transpose(q6))), 'mat6_c2s', abs_tol=tol)
    call assert_close(transform_sample_to_crystal_mat6(rot, &
      & transform_crystal_to_sample_mat6(rot, stiff_c)), stiff_c, 'mat6_round_trip', abs_tol=tol)
    call assert_close(matmul(transform_crystal_to_sample_mat6(rot, stiff_c), &
      & transform_crystal_to_sample_mandel(rot, vec_a)), &
      & transform_crystal_to_sample_mandel(rot, matmul(stiff_c, vec_a)), &
      & 'mat6_law_commutes', abs_tol=1.0e-12_rk)

    ! Skew 3-vectors.
    call assert_close(transform_crystal_to_sample_skew(rot, skew_from_matrix(mat_w)), &
      & skew_from_matrix(matmul(matmul(rot, mat_w), transpose(rot))), 'skew_c2s', abs_tol=tol)
    call assert_close(transform_sample_to_crystal_skew(rot, &
      & transform_crystal_to_sample_skew(rot, skew_from_matrix(mat_w))), &
      & skew_from_matrix(mat_w), 'skew_round_trip', abs_tol=tol)

  end do

  ! ------------------------------------------------------------------------------------------
  ! Physical HCP check of the orientation convention.
  ! ------------------------------------------------------------------------------------------

  ! HCP stiffness (Ti-like, MPa; C33 independent), Mandel form: shear diagonal 2 C44, 2 C44,
  ! 2 C66 with C66 = (C11 - C12)/2.
  stiff_hcp = 0.0_rk
  stiff_hcp(1:3, 1:3) = reshape([162400.0_rk, 92000.0_rk, 69000.0_rk, &
    & 92000.0_rk, 162400.0_rk, 69000.0_rk, 69000.0_rk, 69000.0_rk, 180700.0_rk], [3, 3])
  stiff_hcp(4, 4) = 2.0_rk*46700.0_rk
  stiff_hcp(5, 5) = 2.0_rk*46700.0_rk
  stiff_hcp(6, 6) = 162400.0_rk - 92000.0_rk

  ! Voigt compliances of a transversely isotropic crystal (closed form), and the Mandel
  ! compliance (shear diagonal S44/2, S44/2, S66/2 with S66 = 2 (S11 - S12)).
  call hcp_compliance(162400.0_rk, 92000.0_rk, 69000.0_rk, 180700.0_rk, 46700.0_rk, &
    & s11, s12, s13, s33, s44)
  comp_c = 0.0_rk
  comp_c(1:3, 1:3) = reshape([s11, s12, s13, s12, s11, s13, s13, s13, s33], [3, 3])
  comp_c(4, 4) = 0.5_rk*s44
  comp_c(5, 5) = 0.5_rk*s44
  comp_c(6, 6) = s11 - s12
  call assert_close(matmul(comp_c, stiff_hcp), reshape([(merge(1.0_rk, 0.0_rk, &
    & mod(i, 7) == 1), i=1, 36)], [6, 6]), 'hcp_compliance_inverse', abs_tol=1.0e-12_rk)

  ! Passive (1,1,1)/sqrt3, 120°: c-axis along X_s.
  rot = rot_c2s_from_descriptor(ori_des_axis_angle, ori_conv_passive, &
    & [1.0_rk, 1.0_rk, 1.0_rk, 120.0_rk])
  call assert_close(rot(:, 3), [1.0_rk, 0.0_rk, 0.0_rk], 'hcp_passive_c_axis_x', abs_tol=1.0e-15_rk)
  call check_modulus(rot, [1.0_rk, 0.0_rk, 0.0_rk], 1.0_rk, 'hcp_passive_Ex')
  call check_modulus(rot, [0.0_rk, 1.0_rk, 0.0_rk], 0.0_rk, 'hcp_passive_Ey')
  call check_modulus(rot, [0.0_rk, 0.0_rk, 1.0_rk], 0.0_rk, 'hcp_passive_Ez')

  ! Active, same values: the inverse rotation puts the c-axis along Y_s.
  rot = rot_c2s_from_descriptor(ori_des_axis_angle, ori_conv_active, &
    & [1.0_rk, 1.0_rk, 1.0_rk, 120.0_rk])
  call assert_close(rot(:, 3), [0.0_rk, 1.0_rk, 0.0_rk], 'hcp_active_c_axis_y', abs_tol=1.0e-15_rk)
  call check_modulus(rot, [1.0_rk, 0.0_rk, 0.0_rk], 0.0_rk, 'hcp_active_Ex')
  call check_modulus(rot, [0.0_rk, 1.0_rk, 0.0_rk], 1.0_rk, 'hcp_active_Ey')

  ! Passive (0,1,0), 30°: c-axis (sin 30°, 0, cos 30°).
  rot = rot_c2s_from_descriptor(ori_des_axis_angle, ori_conv_passive, &
    & [0.0_rk, 1.0_rk, 0.0_rk, 30.0_rk])
  call check_modulus(rot, [1.0_rk, 0.0_rk, 0.0_rk], sin(30.0_rk*deg_to_rad), 'hcp_tilt30_Ex')
  call check_modulus(rot, [0.0_rk, 1.0_rk, 0.0_rk], 0.0_rk, 'hcp_tilt30_Ey')
  call check_modulus(rot, [0.0_rk, 0.0_rk, 1.0_rk], cos(30.0_rk*deg_to_rad), 'hcp_tilt30_Ez')

  ! Report; error stop on any failure.
  ! ------------------------------------------------------------------------------------------
  ! Rotation updates (step 5.2).
  ! ------------------------------------------------------------------------------------------

  do i_sample = 1, num_samples
    ! Random skew vector of a given angle (the scales cross the small-angle switch).
    call random_number(skw)
    skw = 2.0_rk*skw - 1.0_rk
    skw = scales(mod(i_sample - 1, size(scales)) + 1)*skw/norm2(skw)
    theta = norm2(skw)
    mat_w = matrix_from_skew(skw)
    omega = [-skw(1), skw(2), -skw(3)]

    ! Against the Taylor series of exp(W).
    exp_series = eye3
    term = eye3
    do k = 1, 40
      term = matmul(term, mat_w)/real(k, rk)
      exp_series = exp_series + term
    end do
    rot = rotation_exp_skew(skw)
    call assert_close(rot, exp_series, 'exp_skew_series', abs_tol=tol)

    ! A proper rotation about omega.
    call assert_close(matmul(transpose(rot), rot), eye3, 'exp_skew_orthogonal', abs_tol=tol)
    call assert_close(det3(rot), 1.0_rk, 'exp_skew_det', abs_tol=tol)
    call assert_close(matmul(rot, omega), omega, 'exp_skew_axis', abs_tol=tol)
    call assert_close(matmul(mat_w, vec_c), cross(omega, vec_c), 'skew_axial_vector', &
      & abs_tol=tol)

    ! v2's axis-angle matrix (passive): exp(W) is its transpose (angles above 1e-9 rad, where
    ! the axis is well defined).
    if (theta > 1.0e-9_rk) call assert_close(rot, transpose(g_from_axis_angle(omega/theta, &
      & theta/deg_to_rad)), 'exp_skew_v2_axis_angle', abs_tol=tol)

    ! The lattice update with the start rotation as estimate: exp(W_s dt) rot_c2s_prev.
    mat_g = random_rotation()
    call assert_close(rot_c2s_spin_update(mat_g, mat_g, skw/0.5_rk, 0.5_rk), &
      & matmul(rot, mat_g), 'spin_update_sample_frame', abs_tol=tol)

    ! Crystal-frame form (step 6.2): equal to the sample-frame form with the spin transformed
    ! by the current estimate, and to rot_c2s_prev exp(W_c dt).
    mat_r = random_rotation()
    call assert_close(rot_c2s_spin_update_crystal(mat_g, skw/0.5_rk, 0.5_rk), &
      & rot_c2s_spin_update(mat_g, mat_r, transform_crystal_to_sample_skew(mat_r, skw/0.5_rk), &
      & 0.5_rk), 'spin_update_crystal_frame', abs_tol=tol)
    call assert_close(rot_c2s_spin_update_crystal(mat_g, skw/0.5_rk, 0.5_rk), &
      & matmul(mat_g, rotation_exp_skew(skw)), 'spin_update_crystal_exp', abs_tol=tol)

    ! Logarithm (step 6.5): the inverse of the exponential (angles below pi), relative to the
    ! angle for the smallest ones.
    call assert_close(rotation_log_skew(rot), skw, 'log_skew_of_exp', &
      & abs_tol=tol*max(theta, 1.0e-3_rk))

    ! Rotation rate of a crystal-frame update: the spin in the sample frame (relative to the
    ! angle, plus the absolute round-off 1e-15 of forming R R_prev^T, over dt).
    call assert_close(rot_c2s_rotation_rate(mat_g, rot_c2s_spin_update_crystal(mat_g, &
      & skw/0.5_rk, 0.5_rk), 0.5_rk), transform_crystal_to_sample_skew(mat_g, skw/0.5_rk), &
      & 'rotation_rate_of_update', abs_tol=(tol*theta + 1.0e-15_rk)/0.5_rk)

    ! exp(log(R)) = R for a random rotation.
    call assert_close(rotation_exp_skew(rotation_log_skew(mat_g)), mat_g, 'exp_of_log_skew', &
      & abs_tol=tol)
  end do

  ! Logarithm near pi (1e-6 and 1e-10 rad below, random axes) and of the identity.
  do i_sample = 1, 2*num_samples
    call random_number(skw)
    skw = 2.0_rk*skw - 1.0_rk
    theta = acos(-1.0_rk) - merge(1.0e-6_rk, 1.0e-10_rk, mod(i_sample, 2) == 0)
    rot = rotation_exp_skew(theta*skw/norm2(skw))
    call assert_close(rotation_exp_skew(rotation_log_skew(rot)), rot, 'exp_of_log_near_pi', &
      & abs_tol=1.0e-12_rk)
  end do
  call assert_close(rotation_log_skew(eye3), [0.0_rk, 0.0_rk, 0.0_rk], 'log_skew_identity', &
    & abs_tol=0.0_rk)

  call test_summary()

contains

  !> Determinant of a 3x3 matrix.
  !>
  !> Arguments:
  !>   a  (in)  real(rk)(3,3)  matrix
  pure function det3(a) result(det)

    real(rk), intent(in) :: a(3, 3)
    real(rk) :: det

    det = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) - a(1, 2)*(a(2, 1)*a(3, 3) - &
      & a(2, 3)*a(3, 1)) + a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))

  end function det3

  !> Cross product of two 3-vectors.
  !>
  !> Arguments:
  !>   a, b  (in)  real(rk)(3)  vectors
  pure function cross(a, b) result(c)

    real(rk), intent(in) :: a(3), b(3)
    real(rk) :: c(3)

    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]

  end function cross

  !> Compare the uniaxial modulus along a sample direction with the analytic directional
  !> modulus of the transversely isotropic crystal.
  !>
  !> Arguments:
  !>   rot_c2s   (in)  real(rk)(3,3)  crystal-to-sample rotation
  !>   dir_s     (in)  real(rk)(3)    unit loading direction, sample frame
  !>   cos_c     (in)  real(rk)       cosine of the angle between dir_s and the c-axis
  !>   label     (in)  character      check label
  !> Locals:
  !>   comp_s    real(rk)(6,6)  Mandel compliance, sample frame
  !>   nn        real(rk)(6)    Mandel vector of dir_s (x) dir_s
  !>   inv_e     real(rk)       1/E from the rotated compliance [1/MPa]
  !>   inv_e_ref real(rk)       1/E from the analytic formula [1/MPa]
  subroutine check_modulus(rot_c2s, dir_s, cos_c, label)

    real(rk), intent(in) :: rot_c2s(3, 3)
    real(rk), intent(in) :: dir_s(3)
    real(rk), intent(in) :: cos_c
    character(len=*), intent(in) :: label

    real(rk) :: comp_s(6, 6), nn(6), inv_e, inv_e_ref

    comp_s = transform_crystal_to_sample_mat6(rot_c2s, comp_c)
    nn = mandel_from_matrix(spread(dir_s, 2, 3)*spread(dir_s, 1, 3))
    inv_e = dot_product(nn, matmul(comp_s, nn))
    inv_e_ref = s11*(1.0_rk - cos_c**2)**2 + s33*cos_c**4 &
      & + (2.0_rk*s13 + s44)*cos_c**2*(1.0_rk - cos_c**2)
    call assert_close(1.0_rk/inv_e, 1.0_rk/inv_e_ref, label, rel_tol=1.0e-12_rk)

  end subroutine check_modulus

  !> Voigt compliances of a transversely isotropic (hexagonal) crystal, closed form (e.g. Nye,
  !> J.F. (1985), Physical Properties of Crystals, Oxford, Ch. VIII):
  !>   S11 = (C11 C33 - C13^2)/((C11 - C12) D),  S12 = -(C12 C33 - C13^2)/((C11 - C12) D),
  !>   S13 = -C13/D,  S33 = (C11 + C12)/D,  S44 = 1/C44,  D = (C11 + C12) C33 - 2 C13^2.
  !> Checked against the matrix inverse by the caller.
  !>
  !> Arguments:
  !>   c11, c12, c13, c33, c44  (in)   real(rk)  Voigt stiffnesses [MPa]
  !>   s11, s12, s13, s33, s44  (out)  real(rk)  Voigt compliances [1/MPa]
  !> Locals:
  !>   den  real(rk)  (C11 + C12) C33 - 2 C13^2
  subroutine hcp_compliance(c11, c12, c13, c33, c44, s11, s12, s13, s33, s44)

    real(rk), intent(in) :: c11, c12, c13, c33, c44
    real(rk), intent(out) :: s11, s12, s13, s33, s44

    real(rk) :: den

    den = (c11 + c12)*c33 - 2.0_rk*c13**2
    s11 = (c11*c33 - c13**2)/((c11 - c12)*den)
    s12 = -(c12*c33 - c13**2)/((c11 - c12)*den)
    s13 = -c13/den
    s33 = (c11 + c12)/den
    s44 = 1.0_rk/c44

  end subroutine hcp_compliance

end program test_rotation
