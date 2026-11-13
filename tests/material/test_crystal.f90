! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of crystal_mod: slip systems, Schmid tensors, families, stiffness matrices.
!>
!> Source : FEPX (ROADMAP step 1.4).
!>
!> Usage: test_crystal <fixture file> (tests/material/fixtures/schmid_reference.txt).
!>
!> Families (step 6.1, ROADMAP D48):
!>   - tokens: round trip token -> family -> token for every structure; tokens of another
!>     structure, near misses (spaces, other case of "all") and "all" itself are not tokens;
!>   - family sizes add up to the full tables (12, 48, 18, 32); defaults: all families, bcc
!>     {110}<111> only; the default build equals the explicit all-families build;
!>   - order: families given out of canonical order (bcc {112} first, bcc {123} then {110}, hcp
!>     pyramidal then basal,
!>     bct reversed) give the same systems as the canonical build, family after family in the
!>     given order (directions, normals, Schmid tensors bit-identical), with slip_family the
!>     position in the given list and slip_table the full-table index;
!>   - Miller labels of known systems (fcc 1, bcc 13 and 25, hcp 7, bct 25).
!> Slip systems, for fcc, bcc (12, 24 and 48), hcp (c/a = 1.587 and 1.633), bct (c/a = 0.5456):
!>   - counts and families (v2 index ranges);
!>   - unit directions and normals, d.n = 0 (tolerance 1e-15);
!>   - p = sym(d (x) n) and q = skw(d (x) n) against explicit 3x3 algebra; p traceless,
!>     |p|^2 = 1/2; matrix(p) + matrix(q) = d (x) n;
!>   - no two systems with equal or opposite d (x) n (p and q together); for bct, exactly the six
!>     conjugate pairs of v2's table (plane and direction swapped: 1-17, 2-18, 3-19, 4-20, 5-6,
!>     11-12) share p up to sign while their skew parts q differ (same resolved shear stress,
!>     different lattice spin); fcc, bcc and hcp have no such pairs;
!>   - span: the p of fcc, bcc and hcp span the 5-dimensional deviatoric space (det of the
!>     5x5 form of sum p p^T > 0);
!>   - crystallographic families: fcc normals {111}, directions <110>; bcc {110}<111>,
!>     {112}<111> and {123}<111> (step 7.3; the 24 distinct {123} planes); hcp basal normals along c, prismatic normals and basal/prismatic directions
!>     in the basal plane, system 1 direction a1 = X_c (Neper >= 4.10 convention);
!>   - reference: d (x) n equals FEPX's own (fixture schmid_reference.txt, written by
!>     dump_schmid.f90 in step 9.13; 1e-15; until step 9.13 v2-fixed's Schmid tensors, in v2's
!>     order, which is FEPX's, also 1e-15).
!> Stiffness:
!>   - symmetric and positive definite (Cholesky), for Cu, alpha-Fe, Ti and beta-Sn constants;
!>   - Mandel Hooke's law equals the Voigt definition sigma_i = C_ij gamma_j (engineering shear)
!>     on random strains;
!>   - symmetry: cubic with C44 = (C11 - C12)/2 is isotropic (invariant under random rotations);
!>     hexagonal is invariant under any rotation about c; tetragonal under 90° about c and,
!>     with C66 /= (C11 - C12)/2, not under 45°;
!>   - v2 kelas cross-check: with v2's constraint C33 = C11 + C12 - C13, the 5x5 deviatoric
!>     block [P][C][P]^T, in v2's 5-vector order (components 3 and 5 swapped), is the diagonal
!>     kelas of v2-fixed crys_type_mod1.f90 (cubic, hcp, bct), and [C]{I} = 3K{I} with v2's
!>     bulk modulus (no deviatoric-volumetric coupling); with an independent C33 the coupling is
!>     present ([C]{I} not parallel to {I}).
!> Hill shear modulus (step 7.1, D59): shear_modulus_hill equals v2-fixed's per-structure
!>   Voigt-Reuss-Hill formulas (crys_type_mod1.f90, written out here): cubic (v2 test 02
!>   constants), hcp with v2's C33 = C11 + C12 - C13, bct (C13 fixed formulas, v2 test 43
!>   constants), relative 1e-13; isotropic C44 = (C11 - C12)/2 gives the shear modulus C44.
!>
!> Locals:
!>   fixture        character(512)  fixture path (command-line argument)
!>   crys           type(crystal_t) crystal under test
!>   i, j, k        integer         loop indices
!>   stiff          real(rk)(6,6)   stiffness
!>   rot            real(rk)(3,3)   rotation
!>   kelas          real(rk)(5)     v2's kelas
!>   mat5           real(rk)(5,5)   5x5 deviatoric block in v2's order
!>   vec            real(rk)(6)     work vector
!>   c11 ... c66    real(rk)        elastic constants [MPa]
!>   g_v, g_r, c_sqr, bulk, s11, s12, s13, s33  real(rk)  v2's Hill formula terms
program test_crystal

  use kinds_mod, only: rk
  use constants_mod, only: identity_mandel
  use tensor_mod, only: matrix_from_mandel, mandel_from_matrix, matrix_from_skew, &
    & skew_from_matrix, mandel_trace, components_from_mandel, mat5_from_mat6
  use rotation_mod, only: transform_crystal_to_sample_mat6
  use orientation_mod, only: g_from_axis_angle
  use crystal_mod, only: crystal_t, crystal_init, crystal_fcc, crystal_bcc, crystal_hcp, &
    & crystal_bct, crystal_structure_from_name, crystal_structure_name, stiffness_cubic, &
    & stiffness_hexagonal, stiffness_tetragonal, crystal_num_families, crystal_family_token, &
    & crystal_family_from_token, crystal_family_size, crystal_default_families, &
    & crystal_system_label, shear_modulus_hill
  use tensor_test_utils_mod, only: seed_random, random_rotation, random_symmetric
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary

  implicit none

  character(len=512) :: fixture
  type(crystal_t) :: crys
  integer :: i, j, k
  real(rk) :: stiff(6, 6), rot(3, 3), kelas(5), mat5(5, 5), vec(6)
  real(rk) :: c11, c12, c13, c33, c44, c66
  real(rk) :: g_v, g_r, c_sqr, bulk, s11, s12, s13, s33

  call get_command_argument(1, fixture)
  call seed_random(20261008)

  ! ------------------------------------------------------------------------------------------
  ! Names.
  ! ------------------------------------------------------------------------------------------

  do i = 1, 4
    call assert_equal(crystal_structure_from_name(crystal_structure_name(i)), i, &
      & 'name_round_trip_' // crystal_structure_name(i))
  end do
  call assert_equal(crystal_structure_from_name('FCC'), 0, 'name_unknown')

  ! ------------------------------------------------------------------------------------------
  ! Slip systems.
  ! ------------------------------------------------------------------------------------------

  ! fcc: 12 systems, 1 family; {111} normals, <110> directions.
  call crystal_init(crys, crystal_fcc)
  call check_geometry('fcc', 12, 1, .true.)
  do i = 1, 12
    call assert_close(abs(crys%slip_normal_c(:, i)), spread(1.0_rk/sqrt(3.0_rk), 1, 3), &
      & 'fcc_normal_111', abs_tol=1.0e-15_rk)
    call assert_equal(count(abs(crys%slip_dir_c(:, i)) < 1.0e-15_rk), 1, 'fcc_dir_110_zero')
  end do
  call check_reference('fcc')

  ! bcc: 12 systems; with {112}: 24, 2 families.
  call crystal_init(crys, crystal_bcc)
  call check_geometry('bcc12', 12, 1, .true.)
  call crystal_init(crys, crystal_bcc, families=[1, 2])
  call check_geometry('bcc24', 24, 2, .true.)
  call assert_true(all(crys%slip_family(1:12) == 1) .and. all(crys%slip_family(13:24) == 2), &
    & 'bcc24_families')
  do i = 1, 24
    call assert_close(abs(crys%slip_dir_c(:, i)), spread(1.0_rk/sqrt(3.0_rk), 1, 3), &
      & 'bcc_dir_111', abs_tol=1.0e-15_rk)
  end do
  do i = 13, 24
    call assert_close(sum(abs(crys%slip_normal_c(:, i))), 4.0_rk/sqrt(6.0_rk), &
      & 'bcc_normal_112', abs_tol=1.0e-15_rk)
  end do
  ! bcc with {123} (step 7.3, not in v2): 48 systems, 3 families; the {123} normals are the 24
  ! distinct planes (components 1, 2, 3 in some order and sign, up to the sign of the plane).
  call crystal_init(crys, crystal_bcc, families=[1, 2, 3])
  call check_geometry('bcc48', 48, 3, .true.)
  call assert_true(all(crys%slip_family(25:48) == 3), 'bcc48_families')
  do i = 25, 48
    call assert_close(abs(crys%slip_dir_c(:, i)), spread(1.0_rk/sqrt(3.0_rk), 1, 3), &
      & 'bcc_123_dir_111', abs_tol=1.0e-15_rk)
    call assert_close(sort3(abs(crys%slip_normal_c(:, i)))*sqrt(14.0_rk), &
      & [1.0_rk, 2.0_rk, 3.0_rk], 'bcc_normal_123', abs_tol=1.0e-14_rk)
    do j = 25, i - 1
      call assert_true(norm2(crys%slip_normal_c(:, i) - crys%slip_normal_c(:, j)) > 0.1_rk &
        & .and. norm2(crys%slip_normal_c(:, i) + crys%slip_normal_c(:, j)) > 0.1_rk, &
        & 'bcc_123_distinct_planes')
    end do
  end do
  ! v2's table has the first two families only.
  call crystal_init(crys, crystal_bcc, families=[1, 2])
  call check_reference('bcc')

  ! hcp: 18 systems, 3 families; c/a = 1.633 (ideal) and 1.587 (Ti, fixture).
  call crystal_init(crys, crystal_hcp, c_over_a=1.633_rk)
  call check_geometry('hcp1.633', 18, 3, .true.)
  call crystal_init(crys, crystal_hcp, c_over_a=1.587_rk)
  call check_geometry('hcp1.587', 18, 3, .true.)
  call assert_true(all(crys%slip_family(1:3) == 1) .and. all(crys%slip_family(4:6) == 2) &
    & .and. all(crys%slip_family(7:18) == 3), 'hcp_families')
  call assert_close(crys%slip_dir_c(:, 1), [1.0_rk, 0.0_rk, 0.0_rk], 'hcp_a1_along_xc', &
    & abs_tol=1.0e-15_rk)
  do i = 1, 3
    call assert_close(crys%slip_normal_c(:, i), [0.0_rk, 0.0_rk, 1.0_rk], 'hcp_basal_normal', &
      & abs_tol=1.0e-15_rk)
  end do
  do i = 1, 6
    call assert_close(crys%slip_dir_c(3, i), 0.0_rk, 'hcp_a_direction_in_basal_plane', &
      & abs_tol=1.0e-15_rk)
  end do
  do i = 4, 6
    call assert_close(crys%slip_normal_c(3, i), 0.0_rk, 'hcp_prismatic_normal_in_basal_plane', &
      & abs_tol=1.0e-15_rk)
  end do
  do i = 7, 18
    call assert_true(abs(crys%slip_dir_c(3, i)) > 0.5_rk, 'hcp_pyramidal_c_plus_a')
  end do
  call check_reference('hcp')

  ! bct: 32 systems, 10 families (v2 ranges); c/a = 0.5456 (beta-Sn, fixture).
  call crystal_init(crys, crystal_bct, c_over_a=0.5456_rk)
  call check_geometry('bct', 32, 10, .false.)
  call assert_true(all(crys%slip_family == [1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, &
    & 7, 7, 8, 8, 9, 9, 9, 9, (10, i=1, 8)]), 'bct_families')
  call check_reference('bct')

  ! ------------------------------------------------------------------------------------------
  ! Families: tokens, sizes, defaults, order (D48), labels.
  ! ------------------------------------------------------------------------------------------

  ! Tokens round trip; sizes add up to the full tables.
  do i = 1, 4
    do j = 1, crystal_num_families(i)
      call assert_equal(crystal_family_from_token(i, crystal_family_token(i, j)), j, &
        & 'token_round_trip_' // crystal_structure_name(i))
    end do
    call assert_equal(sum([(crystal_family_size(i, j), j = 1, crystal_num_families(i))]), &
      & full_size(i), 'family_sizes_' // crystal_structure_name(i))
  end do
  call assert_equal(crystal_num_families(crystal_bct), 10, 'bct_num_families')
  call assert_equal(crystal_family_token(crystal_hcp, 3), '{10-11}<11-23>', 'hcp_token_3')
  call assert_equal(crystal_family_from_token(crystal_bcc, '{112}<111>'), 2, 'bcc_112_token')
  call assert_equal(crystal_family_from_token(crystal_bcc, '{111}<110>'), 0, 'fcc_token_in_bcc')
  call assert_equal(crystal_family_from_token(crystal_fcc, '{111} <110>'), 0, 'token_with_space')
  call assert_equal(crystal_family_from_token(crystal_fcc, 'all'), 0, 'all_not_a_token')
  call assert_equal(crystal_family_from_token(crystal_hcp, '{10-1 1}<11-23>'), 0, 'near_miss')
  call assert_equal(crystal_family_token(crystal_fcc, 2), '', 'token_out_of_range')
  call assert_equal(crystal_family_size(crystal_hcp, 4), 0, 'size_out_of_range')

  ! Defaults: every family, except bcc {110}<111> only; the default build is the explicit one.
  call assert_true(all(crystal_default_families(crystal_fcc) == [1]), 'default_fcc')
  call assert_true(all(crystal_default_families(crystal_bcc) == [1]), 'default_bcc')
  call assert_true(all(crystal_default_families(crystal_hcp) == [1, 2, 3]), 'default_hcp')
  call assert_true(all(crystal_default_families(crystal_bct) == [(i, i = 1, 10)]), 'default_bct')
  call check_order(crystal_hcp, 1.587_rk, [1, 2, 3], 'hcp_explicit_all')
  call crystal_init(crys, crystal_hcp, c_over_a=1.587_rk)
  call assert_true(all(crys%families == [1, 2, 3]) .and. &
    & all(crys%slip_table == [(i, i = 1, 18)]), 'hcp_default_is_canonical')

  ! Order as given (D48).
  call check_order(crystal_bcc, 0.0_rk, [2, 1], 'bcc_112_first')
  call check_order(crystal_bcc, 0.0_rk, [2], 'bcc_112_only')
  call check_order(crystal_bcc, 0.0_rk, [3, 1], 'bcc_123_first')
  call check_order(crystal_hcp, 1.587_rk, [3, 1], 'hcp_pyramidal_basal')
  call check_order(crystal_bct, 0.5456_rk, [(11 - i, i = 1, 10)], 'bct_reversed')
  call check_order(crystal_bct, 0.5456_rk, [4, 10], 'bct_two')
  call crystal_init(crys, crystal_bcc, families=[2, 1])
  call assert_true(all(crys%slip_family(1:12) == 1) .and. all(crys%slip_family(13:24) == 2), &
    & 'bcc_112_first_slip_family')
  call assert_true(all(crys%slip_table == [(i, i = 13, 24), (i, i = 1, 12)]), &
    & 'bcc_112_first_slip_table')
  do i = 1, 12
    call assert_close(sum(abs(crys%slip_normal_c(:, i))), 4.0_rk/sqrt(6.0_rk), &
      & 'bcc_112_first_normals_112', abs_tol=1.0e-15_rk)
  end do

  ! Miller labels (v2's tables).
  call assert_equal(crystal_system_label(crystal_fcc, 1), '(111)[01-1]', 'label_fcc_1')
  call assert_equal(crystal_system_label(crystal_bcc, 13), '(-1-12)[111]', 'label_bcc_13')
  call assert_equal(crystal_system_label(crystal_bcc, 25), '(12-3)[111]', 'label_bcc_25')
  call assert_equal(crystal_family_from_token(crystal_bcc, '{123}<111>'), 3, 'bcc_123_token')
  call assert_equal(crystal_system_label(crystal_hcp, 7), '(10-11)[-2113]', 'label_hcp_7')
  call assert_equal(crystal_system_label(crystal_bct, 25), '(121)[-101]', 'label_bct_25')
  call assert_equal(crystal_system_label(crystal_fcc, 13), '', 'label_out_of_range')

  ! ------------------------------------------------------------------------------------------
  ! Stiffness.
  ! ------------------------------------------------------------------------------------------

  ! Cu (cubic), alpha-Fe (cubic), Ti (hexagonal), beta-Sn (tetragonal), MPa.
  call check_stiffness('cubic_Cu', stiffness_cubic(168400.0_rk, 121400.0_rk, 75400.0_rk), &
    & [168400.0_rk, 121400.0_rk, 121400.0_rk, 168400.0_rk, 75400.0_rk, 75400.0_rk])
  call check_stiffness('cubic_Fe', stiffness_cubic(231400.0_rk, 134700.0_rk, 116400.0_rk), &
    & [231400.0_rk, 134700.0_rk, 134700.0_rk, 231400.0_rk, 116400.0_rk, 116400.0_rk])
  call check_stiffness('hexagonal_Ti', stiffness_hexagonal(162400.0_rk, 92000.0_rk, &
    & 69000.0_rk, 180700.0_rk, 46700.0_rk), [162400.0_rk, 92000.0_rk, 69000.0_rk, &
    & 180700.0_rk, 46700.0_rk, 0.5_rk*(162400.0_rk - 92000.0_rk)])
  call check_stiffness('tetragonal_Sn', stiffness_tetragonal(72300.0_rk, 59400.0_rk, &
    & 35800.0_rk, 88400.0_rk, 22000.0_rk, 24000.0_rk), [72300.0_rk, 59400.0_rk, 35800.0_rk, &
    & 88400.0_rk, 22000.0_rk, 24000.0_rk])

  ! Isotropic cubic (C44 = (C11 - C12)/2): invariant under random rotations.
  stiff = stiffness_cubic(250000.0_rk, 100000.0_rk, 75000.0_rk)
  do i = 1, 50
    call assert_close(transform_crystal_to_sample_mat6(random_rotation(), stiff), stiff, &
      & 'isotropic_cubic_invariant', abs_tol=1.0e-9_rk)
  end do

  ! Hexagonal: invariant under any rotation about c (crystal 3).
  stiff = stiffness_hexagonal(162400.0_rk, 92000.0_rk, 69000.0_rk, 180700.0_rk, 46700.0_rk)
  do i = 1, 12
    rot = transpose(g_from_axis_angle([0.0_rk, 0.0_rk, 1.0_rk], 7.5_rk + 29.0_rk*i))
    call assert_close(transform_crystal_to_sample_mat6(rot, stiff), stiff, &
      & 'hexagonal_transversely_isotropic', abs_tol=1.0e-9_rk)
  end do

  ! Tetragonal: invariant under 90° about c; not under 45° (C66 /= (C11 - C12)/2).
  stiff = stiffness_tetragonal(72300.0_rk, 59400.0_rk, 35800.0_rk, 88400.0_rk, 22000.0_rk, &
    & 24000.0_rk)
  rot = transpose(g_from_axis_angle([0.0_rk, 0.0_rk, 1.0_rk], 90.0_rk))
  call assert_close(transform_crystal_to_sample_mat6(rot, stiff), stiff, &
    & 'tetragonal_invariant_90', abs_tol=1.0e-9_rk)
  rot = transpose(g_from_axis_angle([0.0_rk, 0.0_rk, 1.0_rk], 45.0_rk))
  call assert_true(maxval(abs(transform_crystal_to_sample_mat6(rot, stiff) - stiff)) &
    & > 1000.0_rk, 'tetragonal_not_invariant_45')

  ! v2 kelas cross-check (v2 order: (11-22), 33', 12, 13, 23), with C33 = C11 + C12 - C13.
  ! Cubic.
  c11 = 168400.0_rk
  c12 = 121400.0_rk
  c44 = 75400.0_rk
  kelas = [c11 - c12, c11 - c12, 2.0_rk*c44, 2.0_rk*c44, 2.0_rk*c44]
  call check_kelas('cubic', stiffness_cubic(c11, c12, c44), kelas, (c11 + 2.0_rk*c12)/3.0_rk)
  ! hcp.
  c11 = 162400.0_rk
  c12 = 92000.0_rk
  c13 = 69000.0_rk
  c44 = 46700.0_rk
  c33 = c11 + c12 - c13
  kelas = [c11 - c12, c11 + c12 - 2.0_rk*c13, c11 - c12, 2.0_rk*c44, 2.0_rk*c44]
  call check_kelas('hcp', stiffness_hexagonal(c11, c12, c13, c33, c44), kelas, &
    & (c11 + c12 + c13)/3.0_rk)
  ! bct.
  c11 = 72300.0_rk
  c12 = 59400.0_rk
  c13 = 35800.0_rk
  c44 = 22000.0_rk
  c66 = 24000.0_rk
  c33 = c11 + c12 - c13
  kelas = [c11 - c12, c11 + c12 - 2.0_rk*c13, 2.0_rk*c66, 2.0_rk*c44, 2.0_rk*c44]
  call check_kelas('bct', stiffness_tetragonal(c11, c12, c13, c33, c44, c66), kelas, &
    & (c11 + c12 + c13)/3.0_rk)

  ! Independent C33 (Ti, 180700 /= 185400): deviatoric-volumetric coupling present.
  stiff = stiffness_hexagonal(162400.0_rk, 92000.0_rk, 69000.0_rk, 180700.0_rk, 46700.0_rk)
  vec = matmul(stiff, identity_mandel)
  call assert_true(maxval(abs(vec - (sum(vec(1:3))/3.0_rk)*identity_mandel)) > 1000.0_rk, &
    & 'hcp_independent_c33_couples')

  ! Hill shear modulus against v2-fixed's per-structure formulas (crys_type_mod1.f90).
  ! Cubic (v2 test 02).
  c11 = 245.0e3_rk
  c12 = 155.0e3_rk
  c44 = 62.5e3_rk
  g_v = (c11 - c12 + 3.0_rk*c44)/5.0_rk
  g_r = 5.0_rk*(c11 - c12)*c44/(3.0_rk*c11 - 3.0_rk*c12 + 4.0_rk*c44)
  call assert_close(shear_modulus_hill(stiffness_cubic(c11, c12, c44)), 0.5_rk*(g_v + g_r), &
    & 'hill_cubic_v2', rel_tol=1.0e-13_rk)
  ! Isotropic cubic: G = C44.
  call assert_close(shear_modulus_hill(stiffness_cubic(c11, c12, 0.5_rk*(c11 - c12))), &
    & 0.5_rk*(c11 - c12), 'hill_isotropic', rel_tol=1.0e-13_rk)
  ! hcp with v2's C33 (Ti constants).
  c11 = 162400.0_rk
  c12 = 92000.0_rk
  c13 = 69000.0_rk
  c44 = 46700.0_rk
  c33 = c11 + c12 - c13
  c66 = 0.5_rk*(c11 - c12)
  bulk = (c11 + c12 + c13)/3.0_rk
  g_v = (7.0_rk*c11 - 5.0_rk*c12 + 12.0_rk*c44 + 2.0_rk*c33 - 4.0_rk*c13)/30.0_rk
  c_sqr = (c11 + c12)*c33 - 2.0_rk*c13**2
  g_r = 5.0_rk*c_sqr*c44*c66/(6.0_rk*bulk*c44*c66 + 2.0_rk*c_sqr*(c44 + c66))
  call assert_close(shear_modulus_hill(stiffness_hexagonal(c11, c12, c13, c33, c44)), &
    & 0.5_rk*(g_v + g_r), 'hill_hcp_v2', rel_tol=1.0e-13_rk)
  ! bct (v2 test 43; v2-fixed C13), with v2's C33.
  c11 = 72.3e3_rk
  c12 = 59.4e3_rk
  c13 = 35.8e3_rk
  c44 = 22.0e3_rk
  c66 = 24.0e3_rk
  c33 = c11 + c12 - c13
  g_v = ((2.0_rk*c11 + c33) - (c12 + 2.0_rk*c13) + 3.0_rk*(2.0_rk*c44 + c66))/15.0_rk
  c_sqr = (c11 + c12)*c33 - 2.0_rk*c13**2
  s11 = 0.5_rk*(c33/c_sqr + 1.0_rk/(c11 - c12))
  s12 = 0.5_rk*(c33/c_sqr - 1.0_rk/(c11 - c12))
  s13 = -c13/c_sqr
  s33 = (c11 + c12)/c_sqr
  g_r = 15.0_rk/(4.0_rk*(2.0_rk*s11 + s33) - 4.0_rk*(s12 + 2.0_rk*s13) + &
    & 3.0_rk*(2.0_rk/c44 + 1.0_rk/c66))
  call assert_close(shear_modulus_hill(stiffness_tetragonal(c11, c12, c13, c33, c44, c66)), &
    & 0.5_rk*(g_v + g_r), 'hill_bct_v2', rel_tol=1.0e-13_rk)

  ! Report; error stop on any failure.
  call test_summary()

contains

  !> The three components of a vector in ascending order (bcc {123} normal check).
  !>
  !> Arguments:
  !>   v      (in)  real(rk)(3)  vector
  !>   w      (out, result) real(rk)(3)  its components, ascending
  pure function sort3(v) result(w)

    real(rk), intent(in) :: v(3)
    real(rk) :: w(3)

    ! Smallest, middle and largest component.
    w(1) = minval(v)
    w(3) = maxval(v)
    w(2) = sum(v) - w(1) - w(3)

  end function sort3

  !> Number of systems of the full table of a structure (12, 24, 18, 32).
  !>
  !> Arguments:
  !>   structure  (in)  integer  structure code
  pure function full_size(structure) result(num)

    integer, intent(in) :: structure
    integer :: num

    integer, parameter :: sizes(4) = [12, 48, 18, 32]

    num = sizes(structure)

  end function full_size

  !> Families given in some order give the full-table systems of those families, family after
  !> family in that order, bit-identical to the all-families canonical build; slip_family is the
  !> position in the given list, slip_table the full-table index.
  !>
  !> Arguments:
  !>   structure  (in)  integer     structure code
  !>   ratio      (in)  real(rk)    c/a (hcp, bct) [-]
  !>   families   (in)  integer(:)  canonical families, in the order to test
  !>   name       (in)  character   label
  !> Locals:
  !>   full       type(crystal_t)   canonical build with every family
  !>   built      type(crystal_t)   build with the given families
  !>   i, n, k    integer           given family, full-table system, built system
  subroutine check_order(structure, ratio, families, name)

    integer, intent(in) :: structure
    real(rk), intent(in) :: ratio
    integer, intent(in) :: families(:)
    character(len=*), intent(in) :: name

    type(crystal_t) :: full, built
    integer :: i, k, n

    ! Canonical build and the build under test.
    if (structure == crystal_hcp .or. structure == crystal_bct) then
      call crystal_init(full, structure, c_over_a=ratio, &
        & families=[(i, i = 1, crystal_num_families(structure))])
      call crystal_init(built, structure, c_over_a=ratio, families=families)
    else
      call crystal_init(full, structure, families=[(i, i = 1, crystal_num_families(structure))])
      call crystal_init(built, structure, families=families)
    end if
    call assert_true(all(full%slip_table == [(i, i = 1, full_size(structure))]), &
      & name // '_full_canonical')
    call assert_equal(built%num_families, size(families), name // '_num_families')
    call assert_true(all(built%families == families), name // '_families')
    call assert_equal(built%num_slip, sum([(crystal_family_size(structure, families(i)), &
      & i = 1, size(families))]), name // '_num_slip')

    ! Family after family, in the full table's order within a family.
    k = 0
    do i = 1, size(families)
      do n = 1, full%num_slip
        if (full%families(full%slip_family(n)) /= families(i)) cycle
        k = k + 1
        call assert_equal(built%slip_family(k), i, name // '_slip_family')
        call assert_equal(built%slip_table(k), n, name // '_slip_table')
        ! Bit-identical (written as an ordered comparison, -Wcompare-reals).
        call assert_true(.not. any(abs([built%schmid_sym_c(:, k) - full%schmid_sym_c(:, n), &
          & built%schmid_skew_c(:, k) - full%schmid_skew_c(:, n), &
          & built%slip_dir_c(:, k) - full%slip_dir_c(:, n), &
          & built%slip_normal_c(:, k) - full%slip_normal_c(:, n)]) > 0.0_rk), &
          & name // '_same_system')
      end do
    end do
    call assert_equal(k, built%num_slip, name // '_all_systems')

  end subroutine check_order

  !> Geometric checks of the current crystal (crys).
  !>
  !> Arguments:
  !>   name          (in)  character  case name (labels)
  !>   num_slip      (in)  integer    expected number of systems
  !>   num_families  (in)  integer    expected number of families
  !>   spans         (in)  logical    the p must span the deviatoric space (not checked for bct)
  !> Locals:
  !>   dn            real(rk)(3,3)  d (x) n
  !>   sum_ppt       real(rk)(6,6)  sum of p p^T
  !>   num_pairs     integer        number of system pairs with p equal or opposite
  !>   conjugate_pairs  integer(6)  bct conjugate pairs, encoded 100*j + k
  subroutine check_geometry(name, num_slip, num_families, spans)

    integer, parameter :: conjugate_pairs(6) = [117, 218, 319, 420, 506, 1112]
    character(len=*), intent(in) :: name
    integer, intent(in) :: num_slip
    integer, intent(in) :: num_families
    logical, intent(in) :: spans

    real(rk) :: dn(3, 3), sum_ppt(6, 6)
    integer :: num_pairs

    num_pairs = 0
    call assert_equal(crys%num_slip, num_slip, name // '_num_slip')
    call assert_equal(crys%num_families, num_families, name // '_num_families')
    call assert_equal(size(crys%schmid_sym_c, 2), num_slip, name // '_array_size')
    call assert_true(minval(crys%slip_family) == 1 .and. maxval(crys%slip_family) == &
      & num_families, name // '_family_range')

    sum_ppt = 0.0_rk
    do j = 1, crys%num_slip
      ! Unit vectors, orthogonal.
      call assert_close(norm2(crys%slip_dir_c(:, j)), 1.0_rk, name // '_dir_unit', abs_tol=1.0e-15_rk)
      call assert_close(norm2(crys%slip_normal_c(:, j)), 1.0_rk, name // '_normal_unit', &
        & abs_tol=1.0e-15_rk)
      call assert_close(dot_product(crys%slip_dir_c(:, j), crys%slip_normal_c(:, j)), 0.0_rk, &
        & name // '_dir_normal_orthogonal', abs_tol=1.0e-15_rk)

      ! Schmid tensors against 3x3 algebra.
      dn = spread(crys%slip_dir_c(:, j), 2, 3)*spread(crys%slip_normal_c(:, j), 1, 3)
      call assert_close(crys%schmid_sym_c(:, j), mandel_from_matrix(0.5_rk*(dn + transpose(dn))), &
        & name // '_schmid_sym', abs_tol=1.0e-15_rk)
      call assert_close(crys%schmid_skew_c(:, j), skew_from_matrix(dn), name // '_schmid_skew', &
        & abs_tol=1.0e-15_rk)
      call assert_close(matrix_from_mandel(crys%schmid_sym_c(:, j)) &
        & + matrix_from_skew(crys%schmid_skew_c(:, j)), dn, name // '_sym_plus_skew', &
        & abs_tol=1.0e-15_rk)
      call assert_close(mandel_trace(crys%schmid_sym_c(:, j)), 0.0_rk, name // '_traceless', &
        & abs_tol=1.0e-15_rk)
      call assert_close(dot_product(crys%schmid_sym_c(:, j), crys%schmid_sym_c(:, j)), 0.5_rk, &
        & name // '_p_norm', abs_tol=1.0e-15_rk)

      ! Distinct systems: no other d (x) n equal or opposite; pairs with p equal or opposite
      ! must be the expected conjugate pairs (bct only), with different q.
      do k = j + 1, crys%num_slip
        call assert_true(maxval(abs([crys%schmid_sym_c(:, j) - crys%schmid_sym_c(:, k), &
          & crys%schmid_skew_c(:, j) - crys%schmid_skew_c(:, k)])) > 1.0e-6_rk &
          & .and. maxval(abs([crys%schmid_sym_c(:, j) + crys%schmid_sym_c(:, k), &
          & crys%schmid_skew_c(:, j) + crys%schmid_skew_c(:, k)])) > 1.0e-6_rk, &
          & name // '_distinct')
        if (maxval(abs(crys%schmid_sym_c(:, j) - crys%schmid_sym_c(:, k))) < 1.0e-6_rk &
          & .or. maxval(abs(crys%schmid_sym_c(:, j) + crys%schmid_sym_c(:, k))) < 1.0e-6_rk) then
          num_pairs = num_pairs + 1
          call assert_true(crys%structure == crystal_bct .and. any(j*100 + k == conjugate_pairs), &
            & name // '_same_p_only_bct_conjugate_pairs')
        end if
      end do
      sum_ppt = sum_ppt + spread(crys%schmid_sym_c(:, j), 2, 6)*spread(crys%schmid_sym_c(:, j), 1, 6)
    end do

    ! Number of pairs sharing p (up to sign): 6 for bct, 0 otherwise.
    call assert_equal(num_pairs, merge(6, 0, crys%structure == crystal_bct), name // '_num_same_p_pairs')

    ! Span of the deviatoric space.
    if (spans) then
      call assert_true(determinant(mat5_from_mat6(sum_ppt)) > 1.0e-6_rk, name // '_spans_5d')
    end if

  end subroutine check_geometry

  !> Compare d (x) n of the current crystal with the reference fixture.
  !>
  !> Arguments:
  !>   name   (in)  character  structure name of the fixture case
  !> Locals:
  !>   unit, iostat   integer         file unit and status
  !>   head           character(3)    case name in a header line
  !>   num            integer         number of systems in the case
  !>   ratio          real(rk)        c/a of the case
  !>   reference      real(rk)(9)     one reference tensor, row by row
  !>   dn             real(rk)(3,3)   FEPX d (x) n
  !>   found          logical         the case was found
  subroutine check_reference(name)

    character(len=*), intent(in) :: name

    integer :: unit, iostat, num
    character(len=3) :: head
    real(rk) :: ratio, reference(9), dn(3, 3)
    logical :: found

    found = .false.
    open (newunit=unit, file=trim(fixture), status='old', action='read', iostat=iostat)
    call assert_equal(iostat, 0, 'open ' // trim(fixture))
    if (iostat /= 0) return
    do
      read (unit, *, iostat=iostat) head, num, ratio
      if (iostat /= 0) exit
      if (head == name) then
        found = .true.
        call assert_equal(num, crys%num_slip, 'reference_' // name // '_num_slip')
        call assert_close(ratio, crys%c_over_a, 'reference_' // name // '_c_over_a')
        do j = 1, num
          read (unit, *) reference
          dn = spread(crys%slip_dir_c(:, j), 2, 3)*spread(crys%slip_normal_c(:, j), 1, 3)
          call assert_close(reshape(transpose(dn), [9]), reference, 'reference_' // name // '_schmid', &
            & abs_tol=1.0e-15_rk)
        end do
        exit
      else
        do j = 1, num
          read (unit, *)
        end do
      end if
    end do
    close (unit)
    call assert_true(found, 'reference_' // name // '_case_found')

  end subroutine check_reference

  !> Stiffness checks: symmetric, positive definite, Mandel Hooke's law = Voigt definition.
  !>
  !> Arguments:
  !>   name    (in)  character      case name
  !>   stiff6  (in)  real(rk)(6,6)  Mandel stiffness [MPa]
  !>   consts  (in)  real(rk)(6)    C11, C12, C13, C33, C44, C66 (Voigt) [MPa]
  !> Locals:
  !>   voigt   real(rk)(6,6)  Voigt stiffness (engineering shear strains)
  !>   strain  real(rk)(3,3)  random strain
  !>   gamma   real(rk)(6)    Voigt strain (e11 e22 e33 2e23 2e13 2e12)
  subroutine check_stiffness(name, stiff6, consts)

    character(len=*), intent(in) :: name
    real(rk), intent(in) :: stiff6(6, 6)
    real(rk), intent(in) :: consts(6)

    real(rk) :: voigt(6, 6), strain(3, 3), gamma(6)

    call assert_close(stiff6, transpose(stiff6), name // '_symmetric')
    call assert_true(positive_definite(stiff6), name // '_positive_definite')

    ! Voigt matrix from the constants (C22 = C11, C23 = C13, C55 = C44).
    voigt = 0.0_rk
    voigt(1, 1:3) = [consts(1), consts(2), consts(3)]
    voigt(2, 1:3) = [consts(2), consts(1), consts(3)]
    voigt(3, 1:3) = [consts(3), consts(3), consts(4)]
    voigt(4, 4) = consts(5)
    voigt(5, 5) = consts(5)
    voigt(6, 6) = consts(6)
    do j = 1, 20
      strain = random_symmetric()
      gamma = [strain(1, 1), strain(2, 2), strain(3, 3), 2.0_rk*strain(2, 3), &
        & 2.0_rk*strain(1, 3), 2.0_rk*strain(1, 2)]
      call assert_close(components_from_mandel(matmul(stiff6, mandel_from_matrix(strain))), &
        & matmul(voigt, gamma), name // '_hooke_voigt', abs_tol=1.0e-9_rk)
    end do

  end subroutine check_stiffness

  !> v2 kelas cross-check (constrained C33): [P][C][P]^T in v2's order is diag(kelas), and
  !> [C]{I} = 3K{I}.
  !>
  !> Arguments:
  !>   name      (in)  character      case name
  !>   stiff6    (in)  real(rk)(6,6)  Mandel stiffness [MPa]
  !>   kelas_v2  (in)  real(rk)(5)    v2's kelas [MPa]
  !>   bulk      (in)  real(rk)       v2's bulk modulus [MPa]
  !> Locals:
  !>   perm      integer(5)     FEPX component of each v2 component (3 and 5 swapped)
  !>   expected  real(rk)(5,5)  diag(kelas)
  subroutine check_kelas(name, stiff6, kelas_v2, bulk)

    character(len=*), intent(in) :: name
    real(rk), intent(in) :: stiff6(6, 6)
    real(rk), intent(in) :: kelas_v2(5)
    real(rk), intent(in) :: bulk

    integer, parameter :: perm(5) = [1, 2, 5, 4, 3]
    real(rk) :: expected(5, 5)

    mat5 = mat5_from_mat6(stiff6)
    mat5 = mat5(perm, perm)
    expected = 0.0_rk
    do j = 1, 5
      expected(j, j) = kelas_v2(j)
    end do
    call assert_close(mat5, expected, 'kelas_' // name, abs_tol=1.0e-9_rk)
    call assert_close(matmul(stiff6, identity_mandel), 3.0_rk*bulk*identity_mandel, &
      & 'bulk_decoupled_' // name, abs_tol=1.0e-9_rk)

  end subroutine check_kelas

  !> Positive definiteness of a symmetric matrix by Cholesky factorization (all pivots > 0).
  !>
  !> Arguments:
  !>   mat  (in)  real(rk)(:,:)  symmetric matrix
  !>   pd   (out, result) logical  true if positive definite
  !> Locals:
  !>   low  real(rk)(:,:), allocatable  Cholesky factor
  !>   n, p, q  integer  size and indices
  !>   pivot    real(rk) diagonal pivot
  function positive_definite(mat) result(pd)

    real(rk), intent(in) :: mat(:, :)
    logical :: pd

    real(rk), allocatable :: low(:, :)
    integer :: n, p, q
    real(rk) :: pivot

    n = size(mat, 1)
    allocate (low(n, n))
    low = 0.0_rk
    pd = .true.
    do q = 1, n
      pivot = mat(q, q) - sum(low(q, 1:q - 1)**2)
      if (.not. (pivot > 0.0_rk)) then
        pd = .false.
        return
      end if
      low(q, q) = sqrt(pivot)
      do p = q + 1, n
        low(p, q) = (mat(p, q) - sum(low(p, 1:q - 1)*low(q, 1:q - 1)))/low(q, q)
      end do
    end do

  end function positive_definite

  !> Determinant of a square matrix by Gaussian elimination with partial pivoting.
  !>
  !> Arguments:
  !>   mat  (in)  real(rk)(:,:)  matrix
  !>   det  (out, result) real(rk)  determinant
  !> Locals:
  !>   work     real(rk)(:,:), allocatable  eliminated copy
  !>   n, p, q, r  integer  size, pivot column, pivot row, row
  !>   tmp      real(rk)(:), allocatable  row swap buffer
  function determinant(mat) result(det)

    real(rk), intent(in) :: mat(:, :)
    real(rk) :: det

    real(rk), allocatable :: work(:, :), tmp(:)
    integer :: n, p, q, r

    n = size(mat, 1)
    work = mat
    det = 1.0_rk
    do p = 1, n
      q = p - 1 + maxloc(abs(work(p:n, p)), 1)
      if (q /= p) then
        tmp = work(p, :)
        work(p, :) = work(q, :)
        work(q, :) = tmp
        det = -det
      end if
      det = det*work(p, p)
      if (.not. (abs(work(p, p)) > 0.0_rk)) return
      do r = p + 1, n
        work(r, p:n) = work(r, p:n) - (work(r, p)/work(p, p))*work(p, p:n)
      end do
    end do

  end function determinant

end program test_crystal
