! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Property tests of tensor_mod against explicit 3x3 (and fourth-order) tensor algebra.
!>
!> Source : FEPX (ROADMAP step 1.2; conventions.md §1.5).
!>
!> Fixed examples with hand-computed values pin the component order (11 22 33 23 13 12), the
!> sqrt(2) factors and the skew-vector order (23 13 12). Then, for num_samples random inputs
!> (entries uniform in [-1, 1), fixed seed), every routine is compared with the same operation
!> written in 3x3 or fourth-order index form (tensor_test_utils_mod):
!>   - conversions: round trips, inner product A:B = {a}.{b}, symmetric part of a general
!>     matrix, I/O components, decomposition A = sym(A) + skew(A);
!>   - trace, deviatoric projection (traceless, idempotent), Frobenius norm;
!>   - equivalent stress and strain vs v2's 3x3 formulas, plus uniaxial, shear and hydrostatic
!>     values;
!>   - [Q6(R)]: {R A R^T} = [Q6]{A}; column q = {R E_q R^T} (Eq. R4); orthogonality;
!>     [Q6(R^T)] = [Q6(R)]^T; [Q6(R1 R2)] = [Q6(R1)][Q6(R2)]; [Q6(I)] = I; fourth-order
!>     rotation [C'] = [Q6][C][Q6]^T vs R_ip R_jq R_kr R_ls C_pqrs; Hooke's law {C:e} = [C]{e};
!>   - [W^]: [W^]{A} = {A W - W A} (Eq. R8); antisymmetry; traceless output; consistency with
!>     [Q6]: d/dt [Q6(exp(t W))] at t = 0 is -[W^] (central difference);
!>   - deviatoric 5-vectors (dev/theory/addendum_5vectors.md; KD25 = Kasemer & Dawson 2025,
!>     arXiv:2502.19531): (F1) componentwise and as projections on the basis B1..B5; inner
!>     product (F2); [P][P]^T = I5 and [P]^T[P] = deviatoric projector (F5); round trips;
!>     [P] vs KD25 Eqs. 28-29; the 5x5 form [P][Q6][P]^T of the rotation (mat5_from_mat6),
!>     orthogonal, transpose, composition, action; [W^'] action and antisymmetry; [M*] round
!>     trip through the 5x5 form ([P]^T[M'][P] written out);
!>     mapping to v2's 5-vector order (12 13 23).
!>
!> Tolerances (absolute, for inputs of order 1): tol = 1e-14 for single operations of a few
!> terms; tol_sum = 1e-13 for products of 6x6 or 3x3 matrices; tol_t4 = 1e-12 for the 81-term,
!> 5-factor sums of the fourth-order rotation; tol_fd = 1e-8 for the central difference with
!> step h = 1e-5 (truncation error O(h^2) ~ 1e-10, rounding O(eps/h) ~ 1e-11).
!>
!> Locals:
!>   num_samples    integer        number of random samples per property
!>   tol, tol_sum, tol_t4, tol_fd  real(rk)  tolerances (see above)
!>   h              real(rk)       step of the central difference
!>   i_sample       integer        sample index
!>   p, q, i        integer        Mandel component, Mandel component, loop index
!>   mat_a, mat_b   real(rk)(3,3)  random symmetric matrices A, B
!>   mat_g          real(rk)(3,3)  random general matrix
!>   mat_w          real(rk)(3,3)  random skew matrix W
!>   rot, rot2      real(rk)(3,3)  random rotation R1; R1^T, then a second random rotation R2
!>   basis          real(rk)(3,3)  Mandel basis tensor E_q
!>   vec_a, vec_b   real(rk)(6)    Mandel vectors of A, B
!>   vec_dev        real(rk)(6)    deviatoric part of A
!>   q6, q6_t, q6_2 real(rk)(6,6)  [Q6(R1)], [Q6(R1^T)], [Q6(R2)]
!>   stiff, stiff_rot  real(rk)(6,6)  random symmetric Mandel matrix [C], its rotation
!>   spin           real(rk)(6,6)  [W^]
!>   tensor4        real(rk)(3,3,3,3)  fourth-order form of [C]
!>   eye6           real(rk)(6,6)  6x6 identity
!>   eye5           real(rk)(5,5)  5x5 identity
!>   vec5_a, vec5_b real(rk)(5)    deviatoric 5-vectors of A, B
!>   q5             real(rk)(5,5)  [P][Q6(R1)][P]^T
!>   mat_v          real(rk)(5,6)  KD25 Eq. 28 [V] (tensor components -> 5-vector)
!>   mat_vinv       real(rk)(6,5)  KD25 Eq. 29 [V]^-1 (5-vector -> tensor components)
!>   mat_d          real(rk)(6,6)  diag(1,1,1,sqrt2,sqrt2,sqrt2): Mandel = [D] tensor components
!>   mstar          real(rk)(6,6)  sum of {p}{p}^T over random traceless {p} (an [M*])
!>   vec_p          real(rk)(6)    random traceless Mandel vector
!>   dev_basis      real(rk)(3,3,5) deviatoric basis tensors B1..B5
!>   mat_dev        real(rk)(3,3)  deviatoric part of A, 3x3 form
!>   expected       real(rk)       expected scalar value
program test_tensor

  use kinds_mod, only: rk
  use constants_mod, only: sqrt2, sqrt3, sqrt6, sqrt_3_2, identity_3, identity_mandel
  use tensor_mod, only: mandel_index_i, mandel_index_j, mandel_from_matrix, matrix_from_mandel, &
    & mandel_from_components, components_from_mandel, mandel_trace, mandel_dev, mandel_norm, &
    & equiv_stress, equiv_strain, skew_from_matrix, matrix_from_skew, mandel_rotation, &
    & spin_commutator, dev5_projection, dev5_from_mandel, mandel_from_dev5, mat5_from_mat6
  use tensor_test_utils_mod, only: seed_random, random_matrix, random_symmetric, random_skew, &
    & random_rotation, exp_skew, tensor4_from_mandel, mandel_from_tensor4, rotate_tensor4, &
    & determinant3
  use assert_mod, only: assert_true, assert_close, test_summary

  implicit none

  integer, parameter :: num_samples = 500
  real(rk), parameter :: tol = 1.0e-14_rk
  real(rk), parameter :: tol_sum = 1.0e-13_rk
  real(rk), parameter :: tol_t4 = 1.0e-12_rk
  real(rk), parameter :: tol_fd = 1.0e-8_rk
  real(rk), parameter :: h = 1.0e-5_rk

  integer :: i_sample
  integer :: p, q, i
  real(rk) :: mat_a(3, 3), mat_b(3, 3), mat_g(3, 3), mat_w(3, 3)
  real(rk) :: rot(3, 3), rot2(3, 3), basis(3, 3), mat_dev(3, 3)
  real(rk) :: vec_a(6), vec_b(6), vec_dev(6)
  real(rk) :: q6(6, 6), q6_t(6, 6), q6_2(6, 6)
  real(rk) :: stiff(6, 6), stiff_rot(6, 6), spin(6, 6), eye6(6, 6)
  real(rk) :: eye5(5, 5), vec5_a(5), vec5_b(5), q5(5, 5)
  real(rk) :: mat_v(5, 6), mat_vinv(6, 5), mat_d(6, 6), mstar(6, 6), vec_p(6)
  real(rk) :: dev_basis(3, 3, 5)
  real(rk) :: tensor4(3, 3, 3, 3)
  real(rk) :: expected

  call seed_random(20261008)

  eye6 = 0.0_rk
  do i = 1, 6
    eye6(i, i) = 1.0_rk
  end do
  eye5 = eye6(1:5, 1:5)

  ! KD25 Eq. 28 [V] and Eq. 29 [V]^-1, typed from the paper (rows written out, stored by
  ! columns with reshape and transposed), and the Mandel scaling [D].
  mat_v = transpose(reshape([ &
    & 1.0_rk/sqrt2, -1.0_rk/sqrt2, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
    & 0.0_rk, 0.0_rk, sqrt_3_2, 0.0_rk, 0.0_rk, 0.0_rk, &
    & 0.0_rk, 0.0_rk, 0.0_rk, sqrt2, 0.0_rk, 0.0_rk, &
    & 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, sqrt2, 0.0_rk, &
    & 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, sqrt2], [6, 5]))
  mat_vinv = transpose(reshape([ &
    & 1.0_rk/sqrt2, -1.0_rk/sqrt6, 0.0_rk, 0.0_rk, 0.0_rk, &
    & -1.0_rk/sqrt2, -1.0_rk/sqrt6, 0.0_rk, 0.0_rk, 0.0_rk, &
    & 0.0_rk, sqrt(2.0_rk/3.0_rk), 0.0_rk, 0.0_rk, 0.0_rk, &
    & 0.0_rk, 0.0_rk, 1.0_rk/sqrt2, 0.0_rk, 0.0_rk, &
    & 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk/sqrt2, 0.0_rk, &
    & 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk/sqrt2], [5, 6]))
  mat_d = eye6
  do i = 4, 6
    mat_d(i, i) = sqrt2
  end do

  ! Deviatoric basis B1..B5 (addendum_5vectors.md §1).
  dev_basis = 0.0_rk
  dev_basis(1, 1, 1) = 1.0_rk/sqrt2
  dev_basis(2, 2, 1) = -1.0_rk/sqrt2
  dev_basis(1, 1, 2) = -1.0_rk/sqrt6
  dev_basis(2, 2, 2) = -1.0_rk/sqrt6
  dev_basis(3, 3, 2) = 2.0_rk/sqrt6
  dev_basis(2, 3, 3) = 1.0_rk/sqrt2
  dev_basis(3, 2, 3) = 1.0_rk/sqrt2
  dev_basis(1, 3, 4) = 1.0_rk/sqrt2
  dev_basis(3, 1, 4) = 1.0_rk/sqrt2
  dev_basis(1, 2, 5) = 1.0_rk/sqrt2
  dev_basis(2, 1, 5) = 1.0_rk/sqrt2

  ! ------------------------------------------------------------------------------------------
  ! Fixed examples (hand-computed).
  ! ------------------------------------------------------------------------------------------

  ! Index tables.
  call assert_true(all(mandel_index_i == [1, 2, 3, 2, 1, 1]), 'index_i')
  call assert_true(all(mandel_index_j == [1, 2, 3, 3, 3, 2]), 'index_j')

  ! A = [[1, 6, 5], [6, 2, 4], [5, 4, 3]] -> {1, 2, 3, 4 sqrt2, 5 sqrt2, 6 sqrt2}.
  mat_a = reshape([1.0_rk, 6.0_rk, 5.0_rk, 6.0_rk, 2.0_rk, 4.0_rk, 5.0_rk, 4.0_rk, 3.0_rk], &
    & [3, 3])
  vec_a = [1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk*sqrt2, 5.0_rk*sqrt2, 6.0_rk*sqrt2]
  call assert_close(mandel_from_matrix(mat_a), vec_a, 'fixed_mandel_from_matrix', abs_tol=tol)
  call assert_close(matrix_from_mandel(vec_a), mat_a, 'fixed_matrix_from_mandel', abs_tol=tol)
  call assert_close(components_from_mandel(vec_a), &
    & [1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk, 5.0_rk, 6.0_rk], 'fixed_components', abs_tol=tol)
  call assert_close(mandel_trace(vec_a), 6.0_rk, 'fixed_trace')
  call assert_close(mandel_dev(vec_a), &
    & [-1.0_rk, 0.0_rk, 1.0_rk, 4.0_rk*sqrt2, 5.0_rk*sqrt2, 6.0_rk*sqrt2], 'fixed_dev', &
    & abs_tol=tol)
  ! |A|^2 = 1 + 4 + 9 + 2 (16 + 25 + 36) = 168.
  call assert_close(mandel_norm(vec_a), sqrt(168.0_rk), 'fixed_norm', abs_tol=tol)

  ! G = [[0, 1, 2], [3, 0, 4], [5, 6, 0]]: skew part (G - G^T)/2 has
  ! W23 = (4 - 6)/2 = -1, W13 = (2 - 5)/2 = -1.5, W12 = (1 - 3)/2 = -1.
  mat_g = reshape([0.0_rk, 3.0_rk, 5.0_rk, 1.0_rk, 0.0_rk, 6.0_rk, 2.0_rk, 4.0_rk, 0.0_rk], &
    & [3, 3])
  call assert_close(skew_from_matrix(mat_g), [-1.0_rk, -1.5_rk, -1.0_rk], 'fixed_skew')
  ! Skew vector {w23, w13, w12} = {1, 2, 3} -> W = [[0, 3, 2], [-3, 0, 1], [-2, -1, 0]].
  call assert_close(matrix_from_skew([1.0_rk, 2.0_rk, 3.0_rk]), &
    & reshape([0.0_rk, -3.0_rk, -2.0_rk, 3.0_rk, 0.0_rk, -1.0_rk, 2.0_rk, 1.0_rk, 0.0_rk], &
    & [3, 3]), 'fixed_matrix_from_skew')

  ! Equivalent measures: uniaxial stress sigma -> sigma; pure shear tau -> sqrt(3) tau;
  ! hydrostatic -> 0; incompressible uniaxial strain diag(e, -e/2, -e/2) -> e.
  call assert_close(equiv_stress([250.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]), &
    & 250.0_rk, 'fixed_equiv_stress_uniaxial', abs_tol=1.0e-12_rk)
  call assert_close(equiv_stress(mandel_from_components( &
    & [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, -40.0_rk])), sqrt3*40.0_rk, &
    & 'fixed_equiv_stress_shear', abs_tol=1.0e-12_rk)
  call assert_close(equiv_stress([-7.0_rk, -7.0_rk, -7.0_rk, 0.0_rk, 0.0_rk, 0.0_rk]), &
    & 0.0_rk, 'fixed_equiv_stress_hydrostatic', abs_tol=1.0e-12_rk)
  call assert_close(equiv_strain([0.002_rk, -0.001_rk, -0.001_rk, 0.0_rk, 0.0_rk, 0.0_rk]), &
    & 0.002_rk, 'fixed_equiv_strain_uniaxial', abs_tol=tol)

  ! [Q6] of a 90-degree rotation about e3 (R e1 = e2, R e2 = -e1): R A R^T for
  ! A = [[1,6,5],[6,2,4],[5,4,3]] is [[2,-6,-4],[-6,1,5],[-4,5,3]] (by hand, checked with NumPy).
  rot = reshape([0.0_rk, 1.0_rk, 0.0_rk, -1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk], &
    & [3, 3])
  call assert_close(matmul(mandel_rotation(rot), mandel_from_matrix(mat_a)), &
    & [2.0_rk, 1.0_rk, 3.0_rk, 5.0_rk*sqrt2, -4.0_rk*sqrt2, -6.0_rk*sqrt2], 'fixed_q6_rot90', &
    & abs_tol=tol)

  ! [W^] of W = {1, 2, 3}, applied to A (Mandel vector of A as above): with
  ! W = [[0, 3, 2], [-3, 0, 1], [-2, -1, 0]], A W - W A = [[-56, -16, -10], [-16, 28, 26],
  ! [-10, 26, 28]] (by hand, checked with NumPy) -> Mandel {-56, 28, 28, 26 sqrt2, -10 sqrt2,
  ! -16 sqrt2}.
  vec_a = [1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk*sqrt2, 5.0_rk*sqrt2, 6.0_rk*sqrt2]
  call assert_close(matmul(spin_commutator([1.0_rk, 2.0_rk, 3.0_rk]), vec_a), &
    & [-56.0_rk, 28.0_rk, 28.0_rk, 26.0_rk*sqrt2, -10.0_rk*sqrt2, -16.0_rk*sqrt2], &
    & 'fixed_spin_commutator', abs_tol=1.0e-13_rk)

  ! Deviatoric 5-vectors of A = [[1,6,5],[6,2,4],[5,4,3]] (a'33 = 3 - 2 = 1):
  ! {(1 - 2)/sqrt2, sqrt(3/2), 4 sqrt2, 5 sqrt2, 6 sqrt2}; and of the 90-degree rotated
  ! R A R^T = [[2,-6,-4],[-6,1,5],[-4,5,3]] (a'33 = 1): {1/sqrt2, sqrt(3/2), 5 sqrt2, -4 sqrt2,
  ! -6 sqrt2}, which [P][Q6(R)][P]^T must produce from the former.
  vec_a = [1.0_rk, 2.0_rk, 3.0_rk, 4.0_rk*sqrt2, 5.0_rk*sqrt2, 6.0_rk*sqrt2]
  call assert_close(dev5_from_mandel(vec_a), [-1.0_rk/sqrt2, sqrt_3_2, 4.0_rk*sqrt2, &
    & 5.0_rk*sqrt2, 6.0_rk*sqrt2], 'fixed_dev5', abs_tol=tol)
  call assert_close(matmul(mat5_from_mat6(mandel_rotation(rot)), dev5_from_mandel(vec_a)), [1.0_rk/sqrt2, &
    & sqrt_3_2, 5.0_rk*sqrt2, -4.0_rk*sqrt2, -6.0_rk*sqrt2], 'fixed_q5_rot90', abs_tol=tol)
  call assert_close(dev5_from_mandel(identity_mandel), [0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, &
    & 0.0_rk], 'fixed_dev5_of_identity', abs_tol=tol)
  call assert_close(mat5_from_mat6(mandel_rotation(identity_3)), eye5, 'q5_identity', &
    & abs_tol=tol)

  ! [Q6(I)] = I and [W^(0)] = 0.
  call assert_close(mandel_rotation(identity_3), eye6, 'q6_identity', abs_tol=tol)
  call assert_close(spin_commutator([0.0_rk, 0.0_rk, 0.0_rk]), 0.0_rk*eye6, 'spin_zero')

  ! ------------------------------------------------------------------------------------------
  ! Random property tests.
  ! ------------------------------------------------------------------------------------------

  do i_sample = 1, num_samples

    mat_a = random_symmetric()
    mat_b = random_symmetric()
    mat_g = random_matrix()
    mat_w = random_skew()
    rot = random_rotation()
    vec_a = mandel_from_matrix(mat_a)
    vec_b = mandel_from_matrix(mat_b)

    ! The random rotations are proper rotations (checks the test helper itself).
    call assert_close(matmul(transpose(rot), rot), identity_3, &
      & labelled('helper_rotation_orthogonal'), abs_tol=tol_sum)
    call assert_close(determinant3(rot), 1.0_rk, labelled('helper_rotation_det'), &
      & abs_tol=tol_sum)

    ! --- Conversions ---------------------------------------------------------------------

    ! Round trips.
    call assert_close(matrix_from_mandel(vec_a), mat_a, labelled('round_trip_matrix'), &
      & abs_tol=tol)
    call assert_close(mandel_from_matrix(matrix_from_mandel(vec_a)), vec_a, &
      & labelled('round_trip_mandel'), abs_tol=tol)
    call assert_close(mandel_from_components(components_from_mandel(vec_a)), vec_a, &
      & labelled('round_trip_components'), abs_tol=tol)

    ! I/O components are the plain tensor components in order 11 22 33 23 13 12.
    call assert_close(components_from_mandel(vec_a), [(mat_a(mandel_index_i(p), &
      & mandel_index_j(p)), p=1, 6)], labelled('components_are_tensor_components'), &
      & abs_tol=tol)

    ! Inner product preserved: A:B = {a}.{b}.
    call assert_close(dot_product(vec_a, vec_b), sum(mat_a*mat_b), labelled('inner_product'), &
      & abs_tol=tol_sum)

    ! A general matrix converts as its symmetric part, and G = sym(G) + skew(G).
    call assert_close(mandel_from_matrix(mat_g), &
      & mandel_from_matrix(0.5_rk*(mat_g + transpose(mat_g))), labelled('symmetric_part'), &
      & abs_tol=tol)
    call assert_close(matrix_from_mandel(mandel_from_matrix(mat_g)) &
      & + matrix_from_skew(skew_from_matrix(mat_g)), mat_g, labelled('sym_skew_decomposition'), &
      & abs_tol=tol)

    ! Skew vectors: round trip, skew part of a symmetric matrix is zero.
    call assert_close(matrix_from_skew(skew_from_matrix(mat_w)), mat_w, &
      & labelled('round_trip_skew'), abs_tol=tol)
    call assert_close(skew_from_matrix(mat_a), [0.0_rk, 0.0_rk, 0.0_rk], &
      & labelled('skew_of_symmetric'), abs_tol=tol)

    ! --- Trace, deviatoric part, norm -----------------------------------------------------

    call assert_close(mandel_trace(vec_a), mat_a(1, 1) + mat_a(2, 2) + mat_a(3, 3), &
      & labelled('trace'), abs_tol=tol)
    mat_dev = mat_a - (sum([(mat_a(i, i), i=1, 3)])/3.0_rk)*identity_3
    vec_dev = mandel_dev(vec_a)
    call assert_close(matrix_from_mandel(vec_dev), mat_dev, labelled('dev_vs_3x3'), abs_tol=tol)
    call assert_close(mandel_trace(vec_dev), 0.0_rk, labelled('dev_traceless'), abs_tol=tol)
    call assert_close(mandel_dev(vec_dev), vec_dev, labelled('dev_idempotent'), abs_tol=tol)
    call assert_close(mandel_norm(vec_a), sqrt(sum(mat_a*mat_a)), labelled('norm'), &
      & abs_tol=tol)

    ! --- Equivalent measures (v2 3x3 formulas) -------------------------------------------

    expected = sqrt(1.5_rk*sum(mat_dev*mat_dev))
    call assert_close(equiv_stress(vec_a), expected, labelled('equiv_stress'), abs_tol=tol)
    expected = sqrt((2.0_rk/3.0_rk)*sum(mat_dev*mat_dev))
    call assert_close(equiv_strain(vec_a), expected, labelled('equiv_strain'), abs_tol=tol)

    ! --- [Q6(R)] ---------------------------------------------------------------------------

    q6 = mandel_rotation(rot)

    ! Action: {R A R^T} = [Q6]{A}.
    call assert_close(matmul(q6, vec_a), &
      & mandel_from_matrix(matmul(rot, matmul(mat_a, transpose(rot)))), &
      & labelled('q6_action'), abs_tol=tol_sum)

    ! Column definition of Eq. (R4): column q is {R E_q R^T}.
    do q = 1, 6
      ! E_q = e_i (x) e_i (q <= 3) or (e_i (x) e_j + e_j (x) e_i)/sqrt(2) (q >= 4).
      basis = 0.0_rk
      if (q <= 3) then
        basis(mandel_index_i(q), mandel_index_i(q)) = 1.0_rk
      else
        basis(mandel_index_i(q), mandel_index_j(q)) = 1.0_rk/sqrt2
        basis(mandel_index_j(q), mandel_index_i(q)) = 1.0_rk/sqrt2
      end if
      call assert_close(q6(:, q), mandel_from_matrix(matmul(rot, matmul(basis, &
        & transpose(rot)))), labelled('q6_column_definition'), abs_tol=tol_sum)
    end do

    ! Orthogonality, transpose, composition.
    call assert_close(matmul(transpose(q6), q6), eye6, labelled('q6_orthogonal'), &
      & abs_tol=tol_sum)
    rot2 = transpose(rot)
    q6_t = mandel_rotation(rot2)
    call assert_close(q6_t, transpose(q6), labelled('q6_transpose'), abs_tol=tol)
    rot2 = random_rotation()
    q6_2 = mandel_rotation(rot2)
    call assert_close(mandel_rotation(matmul(rot, rot2)), matmul(q6, q6_2), &
      & labelled('q6_composition'), abs_tol=tol_sum)

    ! Inverse rotation back: [Q6(R)]^T [Q6(R)] {A} = {A}, i.e. sample -> crystal -> sample.
    call assert_close(matmul(transpose(q6), matmul(q6, vec_a)), vec_a, &
      & labelled('q6_round_trip'), abs_tol=tol_sum)

    ! Fourth-order: random symmetric Mandel matrix (major and minor symmetries);
    ! [Q6][C][Q6]^T equals the index-form rotation of C.
    call random_number(stiff)
    stiff = stiff + transpose(stiff) - 1.0_rk
    tensor4 = tensor4_from_mandel(stiff)
    call assert_close(mandel_from_tensor4(tensor4), stiff, labelled('helper_tensor4_round_trip'), &
      & abs_tol=tol)
    stiff_rot = matmul(q6, matmul(stiff, transpose(q6)))
    call assert_close(stiff_rot, mandel_from_tensor4(rotate_tensor4(rot, tensor4)), &
      & labelled('q6_fourth_order'), abs_tol=tol_t4)

    ! Hooke's law in Mandel form: {C:A} = [C]{A}, with (C:A)_ij = C_ijkl A_kl.
    call assert_close(matmul(stiff, vec_a), &
      & mandel_from_matrix(reshape(matmul(reshape(tensor4, [9, 9]), reshape(mat_a, [9])), &
      & [3, 3])), labelled('hooke_mandel'), abs_tol=tol_sum)

    ! --- [W^] -------------------------------------------------------------------------------

    spin = spin_commutator(skew_from_matrix(mat_w))

    ! Action, Eq. (R8): [W^]{A} = {A W - W A}.
    call assert_close(matmul(spin, vec_a), &
      & mandel_from_matrix(matmul(mat_a, mat_w) - matmul(mat_w, mat_a)), &
      & labelled('spin_action'), abs_tol=tol_sum)

    ! Antisymmetric; output traceless.
    call assert_close(spin, -transpose(spin), labelled('spin_antisymmetric'), abs_tol=tol)
    call assert_close(mandel_trace(matmul(spin, vec_a)), 0.0_rk, labelled('spin_traceless'), &
      & abs_tol=tol_sum)

    ! Consistency with [Q6]: d/dt [Q6(exp(t W))] at t = 0 is -[W^], since
    ! d/dt (R A R^T) = W A - A W at R = I.
    call assert_close((mandel_rotation(exp_skew(h*mat_w)) &
      & - mandel_rotation(exp_skew(-h*mat_w)))/(2.0_rk*h), -spin, &
      & labelled('spin_q6_derivative'), abs_tol=tol_fd)

    ! --- Deviatoric 5-vectors ---------------------------------------------------------------

    vec5_a = dev5_from_mandel(vec_a)
    vec5_b = dev5_from_mandel(vec_b)

    ! (F1) componentwise from the 3x3 deviatoric part (mat_dev of A, computed above).
    call assert_close(vec5_a, [(mat_a(1, 1) - mat_a(2, 2))/sqrt2, sqrt_3_2*mat_dev(3, 3), &
      & sqrt2*mat_a(2, 3), sqrt2*mat_a(1, 3), sqrt2*mat_a(1, 2)], labelled('dev5_components'), &
      & abs_tol=tol)

    ! Components are the projections A:B_k on the deviatoric basis.
    call assert_close(vec5_a, [(sum(mat_a*dev_basis(:, :, i)), i=1, 5)], &
      & labelled('dev5_basis_projection'), abs_tol=tol)

    ! (F2): A':B' = {a'}.{b'}; norm and equivalent stress.
    call assert_close(dot_product(vec5_a, vec5_b), sum(mat_dev*(matrix_from_mandel( &
      & mandel_dev(vec_b)))), labelled('dev5_inner_product'), abs_tol=tol_sum)
    call assert_close(norm2(vec5_a), mandel_norm(mandel_dev(vec_a)), labelled('dev5_norm'), &
      & abs_tol=tol)
    call assert_close(sqrt_3_2*norm2(vec5_a), equiv_stress(vec_a), labelled('dev5_equiv'), &
      & abs_tol=tol)

    ! Round trips: 5 -> 6 -> 5 is the identity; 6 -> 5 -> 6 is the deviatoric part.
    call assert_close(dev5_from_mandel(mandel_from_dev5(vec5_a)), vec5_a, &
      & labelled('dev5_round_trip_5'), abs_tol=tol)
    call assert_close(mandel_from_dev5(vec5_a), mandel_dev(vec_a), &
      & labelled('dev5_round_trip_6'), abs_tol=tol)

    ! Functions agree with the matrix [P]; (F5).
    call assert_close(matmul(dev5_projection, vec_a), vec5_a, labelled('dev5_matrix_form'), &
      & abs_tol=tol)
    call assert_close(matmul(transpose(dev5_projection), vec5_a), mandel_from_dev5(vec5_a), &
      & labelled('dev5_matrix_form_back'), abs_tol=tol)
    if (i_sample == 1) then
      call assert_close(matmul(dev5_projection, transpose(dev5_projection)), eye5, &
        & 'dev5_p_pt_identity', abs_tol=tol)
      call assert_close(matmul(transpose(dev5_projection), dev5_projection), &
        & eye6 - spread(identity_mandel, 2, 6)*spread(identity_mandel, 1, 6)/3.0_rk, &
        & 'dev5_pt_p_projector', abs_tol=tol)
      ! KD25 Eq. 29: [P]^T = [D][V]^-1; KD25 Eq. 28: [V][V]^-1 = I5.
      call assert_close(transpose(dev5_projection), matmul(mat_d, mat_vinv), &
        & 'dev5_kd25_eq29', abs_tol=tol)
      call assert_close(matmul(mat_v, mat_vinv), eye5, 'dev5_kd25_eq28_right_inverse', &
        & abs_tol=tol)
    end if

    ! KD25 Eq. 28 on the tensor components of A' gives the same 5-vector (deviatoric input).
    call assert_close(matmul(mat_v, components_from_mandel(mandel_dev(vec_a))), vec5_a, &
      & labelled('dev5_kd25_eq28'), abs_tol=tol)

    ! The 5x5 rotation [P][Q6][P]^T: action; orthogonal; transpose; composition.
    q5 = mat5_from_mat6(q6)
    call assert_close(matmul(q5, vec5_a), &
      & dev5_from_mandel(mandel_from_matrix(matmul(rot, matmul(mat_a, transpose(rot))))), &
      & labelled('q5_action'), abs_tol=tol_sum)
    call assert_close(matmul(transpose(q5), q5), eye5, labelled('q5_orthogonal'), &
      & abs_tol=tol_sum)
    rot2 = transpose(rot)
    call assert_close(mat5_from_mat6(mandel_rotation(rot2)), transpose(q5), labelled('q5_transpose'), &
      & abs_tol=tol)
    rot2 = random_rotation()
    call assert_close(mat5_from_mat6(mandel_rotation(matmul(rot, rot2))), &
      & matmul(q5, mat5_from_mat6(mandel_rotation(rot2))), &
      & labelled('q5_composition'), abs_tol=tol_sum)

    ! [W^'] = [P][W^][P]^T: action on A', antisymmetry.
    call assert_close(matmul(mat5_from_mat6(spin), vec5_a), dev5_from_mandel(mandel_from_matrix( &
      & matmul(mat_dev, mat_w) - matmul(mat_w, mat_dev))), labelled('spin5_action'), &
      & abs_tol=tol_sum)
    call assert_close(mat5_from_mat6(spin), -transpose(mat5_from_mat6(spin)), &
      & labelled('spin5_antisymmetric'), abs_tol=tol)

    ! [M*] = sum of {p}{p}^T over traceless {p} survives the 5x5 round trip exactly.
    mstar = 0.0_rk
    do i = 1, 4
      vec_p = mandel_dev(mandel_from_matrix(random_symmetric()))
      mstar = mstar + spread(vec_p, 2, 6)*spread(vec_p, 1, 6)
    end do
    call assert_close(matmul(transpose(dev5_projection), matmul(mat5_from_mat6(mstar), &
      & dev5_projection)), mstar, &
      & labelled('mstar_round_trip_5x5'), abs_tol=tol_sum)

    ! v2's 5-vector (order 12 13 23, mat_vec_symm) on deviatoric input = ours with components
    ! 3 and 5 swapped.
    call assert_close([(mat_dev(1, 1) - mat_dev(2, 2))/sqrt2, sqrt_3_2*mat_dev(3, 3), &
      & sqrt2*mat_dev(1, 2), sqrt2*mat_dev(1, 3), sqrt2*mat_dev(2, 3)], &
      & [vec5_a(1), vec5_a(2), vec5_a(5), vec5_a(4), vec5_a(3)], labelled('dev5_v2_order'), &
      & abs_tol=tol)

  end do

  ! Report; error stop on any failure.
  call test_summary()

contains

  !> Check label with the current sample number appended, e.g. "q6_action[17]".
  !>
  !> Arguments:
  !>   name   (in)  character  property name
  !>   label  (out, result) character(64)  name[i_sample]
  function labelled(name) result(label)

    character(len=*), intent(in) :: name
    character(len=64) :: label

    write (label, '(a,"[",i0,"]")') name, i_sample

  end function labelled

end program test_tensor
