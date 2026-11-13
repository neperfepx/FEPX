! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of the material-solver contract (material_state_mod, material_model_mod) and of
!> the elastic-only model (elastic_model_mod).
!>
!> Source : FEPX (ROADMAP steps 4.1 and 5.2).
!>
!> Point kernel elastic_point_solve (cubic Cu and hexagonal Ti constants with an independent
!> C33, random orientations; 200 random points each with |d| ~ 1e-3/s, |w| ~ 1/s, dt = 0.1 s,
!> |e^e_0| ~ 1e-3):
!>   - converged (status 0); [s] = (dt/beta)[C_s] with beta = det(I + e^e) and [C_s] of the
!>     returned orientation, symmetric;
!>   - orientation (step 5.2): a rotation (1e-14), the fixed point of the backward-Euler update
!>     rot = rot_prev exp(W^p_c dt) (W^p_c in the frame of rot) to the loop tolerance (1e-8:
!>     the loop stops on the stress norm, 1e-5 relative, and the last update used the previous
!>     estimate's frame; |w dt| ~ 0.1 rad here);
!>   - (R2) holds exactly at the result: w^p = w + 1/2 skw(e dot e - e e dot) (1e-14);
!>   - (R1) holds to the loop tolerance: {d} - ({e^e} - {e^e_0})/dt - [W^]{e^e} relative to |d|
!>     (1e-8; the last strain solve used the previous iterate of w^p);
!>   - (R10) sigma = [s]({d} - {h}) and sigma = [C_s]{e^e}/beta, relative 1e-6 (the stress is
!>     that of the last solve, [C_s] that of the final orientation, as v2; they differ by the last
!>     orientation update, about 5e-8 relative at |w dt| ~ 0.1 rad);
!>   - no spin and e^e_0 = 0: e^e = dt d exactly (w^p = 0 since d commutes with itself), the
!>     orientation unchanged;
!>   - pure spin (d = 0, e^e_0 = 0): zero strain and stress, rot = exp(W dt) rot_prev (1e-14);
!>   - frame invariance: rotating the inputs by Q (rot_c2s -> Q rot_c2s, {d}, {e^e_0} -> [Q6]
!>     {.}, w -> skw(Q W Q^T)) rotates every output the same way (1e-12), rot -> Q rot;
!>   - iteration limit: max_iters_state = 1 with a nonzero spin gives status 1.
!> Contract (solve_points through class(material_model_t)):
!>   - material_state_alloc: sizes, zero fields, identity orientations, crss with 0 rows;
!>   - only the points of the given element list change; tangent = secant; num_failed counts
!>     the points that hit the iteration limit.
!>
!> Locals:
!>   stiff_cu, stiff_ti  real(rk)(6,6)  crystal-frame stiffnesses [MPa]
!>   zero6               real(rk)(6)    zero 6-vector
program test_elastic_model

  use kinds_mod, only: rk, ik
  use tensor_mod, only: mandel_from_matrix, matrix_from_mandel, matrix_from_skew, &
    & skew_from_matrix, spin_commutator, mandel_rotation
  use rotation_mod, only: transform_crystal_to_sample_mat6, rot_c2s_spin_update, &
    & rotation_exp_skew
  use crystal_mod, only: stiffness_cubic, stiffness_hexagonal
  use material_state_mod, only: material_state_t, material_state_alloc
  use material_model_mod, only: material_model_t
  use elastic_model_mod, only: elastic_model_t, elastic_model_init, elastic_point_solve
  use tensor_test_utils_mod, only: seed_random, random_rotation, random_symmetric, random_skew, &
    & determinant3
  use assert_mod, only: assert_true, assert_equal, assert_close, test_summary

  implicit none

  real(rk) :: stiff_cu(6, 6), stiff_ti(6, 6)
  real(rk), parameter :: zero6(6) = 0.0_rk

  call seed_random(20261009)
  stiff_cu = stiffness_cubic(168.4e3_rk, 121.4e3_rk, 75.4e3_rk)
  stiff_ti = stiffness_hexagonal(162.4e3_rk, 92.0e3_rk, 69.0e3_rk, 180.7e3_rk, 46.7e3_rk)

  call check_point_solve(stiff_cu, 'Cu')
  call check_point_solve(stiff_ti, 'Ti')
  call check_no_spin(stiff_cu)
  call check_invariance(stiff_ti)
  call check_iteration_limit(stiff_cu)
  call check_contract(stiff_ti)

  call test_summary()

contains

  !> Random inputs of one point.
  !>
  !> Arguments:
  !>   rot, defrate, spin, strain_prev  (out)  orientation, {d}, w, {e^e_0}
  subroutine random_point(rot, defrate, spin, strain_prev)

    real(rk), intent(out) :: rot(3, 3), defrate(6), spin(3), strain_prev(6)

    rot = random_rotation()
    defrate = 1e-3_rk*mandel_from_matrix(random_symmetric())
    spin = skew_from_matrix(random_skew())
    strain_prev = 1e-3_rk*mandel_from_matrix(random_symmetric())

  end subroutine random_point

  !> Equations (R1), (R2), (R10), (R11) at the converged point.
  !>
  !> Arguments:
  !>   stiff  (in)  real(rk)(6,6)  crystal-frame stiffness [MPa]
  !>   label  (in)  character      check label prefix
  !> Locals: inputs and outputs of elastic_point_solve; stiff_s [C_s]; beta; residuals.
  subroutine check_point_solve(stiff, label)

    real(rk), intent(in) :: stiff(6, 6)
    character(len=*), intent(in) :: label

    real(rk), parameter :: dt = 0.1_rk
    real(rk) :: rot(3, 3), defrate(6), spin(3), strain_prev(6)
    real(rk) :: stress(6), strain(6), spin_lat(3), secant(6, 6), vec_h(6), stiff_s(6, 6)
    real(rk) :: strain_mat(3, 3), rate_mat(3, 3), beta, rot_end(3, 3)
    integer :: i_pt, status

    do i_pt = 1, 200
      call random_point(rot, defrate, spin, strain_prev)
      call elastic_point_solve(stiff, rot, strain_prev, defrate, spin, dt, 1e-5_rk, 100, stress, &
        & strain, spin_lat, rot_end, secant, vec_h, status)
      call assert_equal(status, 0, label // ' converged')

      ! Orientation: a rotation, the fixed point of the backward-Euler update.
      call assert_close(reshape(matmul(transpose(rot_end), rot_end), [9]), [1.0_rk, 0.0_rk, &
        & 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk], label // ' rotation', &
        & abs_tol=1e-14_rk)
      call assert_close(reshape(rot_end, [9]), reshape(rot_c2s_spin_update(rot, rot_end, &
        & spin_lat, dt), [9]), label // ' backward-Euler orientation', abs_tol=1e-8_rk)

      ! (R11) and beta.
      stiff_s = transform_crystal_to_sample_mat6(rot_end, stiff)
      strain_mat = matrix_from_mandel(strain)
      beta = determinant3(strain_mat + reshape([1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk, &
        & 0.0_rk, 0.0_rk, 1.0_rk], [3, 3]))
      call assert_close(reshape(secant, [36]), reshape(dt/beta*stiff_s, [36]), label // &
        & ' [s] = (dt/beta)[C_s]', rel_tol=1e-14_rk)
      call assert_close(reshape(secant, [36]), reshape(transpose(secant), [36]), label // &
        & ' [s] symmetric', abs_tol=1e-12_rk*maxval(abs(secant)))

      ! (R2) exactly at the result.
      rate_mat = matrix_from_mandel((strain - strain_prev)/dt)
      call assert_close(spin_lat, spin + skew_from_matrix(0.5_rk*(matmul(rate_mat, strain_mat) - &
        & matmul(strain_mat, rate_mat))), label // ' (R2) lattice spin', abs_tol=1e-14_rk)

      ! (R1) to the loop tolerance.
      call assert_close(defrate - (strain - strain_prev)/dt - matmul(spin_commutator(spin_lat), &
        & strain), zero6, label // ' (R1) residual', &
        & abs_tol=1e-8_rk*maxval(abs(defrate)))

      ! Stress: Hooke's law and (R10).
      call assert_close(stress, matmul(stiff_s, strain)/beta, label // ' sigma = [C_s]{e}/beta', &
        & abs_tol=1e-6_rk*maxval(abs(stress)))
      call assert_close(stress, matmul(secant, defrate - vec_h), label // &
        & ' (R10) sigma = [s]({d} - {h})', abs_tol=1e-6_rk*maxval(abs(stress)))
    end do

  end subroutine check_point_solve

  !> No spin and zero start strain: e^e = dt d exactly, w^p = 0, sigma = [C_s]{e^e}/beta.
  !>
  !> Arguments:
  !>   stiff  (in)  real(rk)(6,6)  crystal-frame stiffness [MPa]
  subroutine check_no_spin(stiff)

    real(rk), intent(in) :: stiff(6, 6)

    real(rk), parameter :: dt = 0.1_rk
    real(rk) :: rot(3, 3), defrate(6), spin(3), strain_prev(6)
    real(rk) :: stress(6), strain(6), spin_lat(3), secant(6, 6), vec_h(6), rot_end(3, 3)
    integer :: i_pt, status

    do i_pt = 1, 50
      call random_point(rot, defrate, spin, strain_prev)
      call elastic_point_solve(stiff, rot, zero6, defrate, [0.0_rk, 0.0_rk, 0.0_rk], dt, &
        & 1e-5_rk, 100, stress, strain, spin_lat, rot_end, secant, vec_h, status)
      call assert_equal(status, 0, 'no spin converged')
      call assert_close(strain, dt*defrate, 'no spin: e = dt d', abs_tol=1e-18_rk)
      call assert_close(spin_lat, [0.0_rk, 0.0_rk, 0.0_rk], 'no spin: w^p = 0', abs_tol=1e-18_rk)
      call assert_close(reshape(rot_end, [9]), reshape(rot, [9]), &
        & 'no spin: orientation unchanged', abs_tol=1e-15_rk)

      ! Pure spin: no strain or stress; the lattice turns with the spin, exp(W dt) rot.
      call elastic_point_solve(stiff, rot, zero6, zero6, spin, dt, 1e-5_rk, 100, stress, &
        & strain, spin_lat, rot_end, secant, vec_h, status)
      call assert_equal(status, 0, 'pure spin converged')
      call assert_close(stress, zero6, 'pure spin: no stress', abs_tol=1e-18_rk)
      call assert_close(spin_lat, spin, 'pure spin: w^p = w', abs_tol=1e-18_rk)
      call assert_close(reshape(rot_end, [9]), reshape(matmul(rotation_exp_skew(dt*spin), rot), &
        & [9]), 'pure spin: rot = exp(W dt) rot_prev', abs_tol=1e-14_rk)
    end do

  end subroutine check_no_spin

  !> Frame invariance of the point solve.
  !>
  !> Arguments:
  !>   stiff  (in)  real(rk)(6,6)  crystal-frame stiffness [MPa]
  subroutine check_invariance(stiff)

    real(rk), intent(in) :: stiff(6, 6)

    real(rk), parameter :: dt = 0.1_rk
    real(rk) :: rot(3, 3), defrate(6), spin(3), strain_prev(6), qrot(3, 3), q6(6, 6)
    real(rk) :: stress(6), strain(6), spin_lat(3), secant(6, 6), vec_h(6), rot_end(3, 3)
    real(rk) :: stress2(6), strain2(6), spin_lat2(3), secant2(6, 6), vec_h2(6), rot_end2(3, 3)
    integer :: i_pt, status

    do i_pt = 1, 100
      call random_point(rot, defrate, spin, strain_prev)
      call elastic_point_solve(stiff, rot, strain_prev, defrate, spin, dt, 1e-12_rk, 100, stress, &
        & strain, spin_lat, rot_end, secant, vec_h, status)
      qrot = random_rotation()
      q6 = mandel_rotation(qrot)
      call elastic_point_solve(stiff, matmul(qrot, rot), matmul(q6, strain_prev), &
        & matmul(q6, defrate), skew_from_matrix(matmul(qrot, matmul(matrix_from_skew(spin), &
        & transpose(qrot)))), dt, 1e-12_rk, 100, stress2, strain2, spin_lat2, rot_end2, secant2, &
        & vec_h2, status)
      call assert_close(reshape(rot_end2, [9]), reshape(matmul(qrot, rot_end), [9]), &
        & 'invariance: orientation', abs_tol=1e-13_rk)
      call assert_close(stress2, matmul(q6, stress), 'invariance: stress', &
        & abs_tol=1e-12_rk*maxval(abs(stress)))
      call assert_close(strain2, matmul(q6, strain), 'invariance: elastic strain', &
        & abs_tol=1e-15_rk)
      call assert_close(spin_lat2, skew_from_matrix(matmul(qrot, matmul(matrix_from_skew( &
        & spin_lat), transpose(qrot)))), 'invariance: lattice spin', abs_tol=1e-13_rk)
      call assert_close(reshape(secant2, [36]), reshape(matmul(q6, matmul(secant, &
        & transpose(q6))), [36]), 'invariance: [s]', abs_tol=1e-12_rk*maxval(abs(secant)))
      call assert_close(vec_h2, matmul(q6, vec_h), 'invariance: {h}', &
        & abs_tol=1e-12_rk*maxval(abs(vec_h)))
    end do

  end subroutine check_invariance

  !> The iteration limit is reported.
  !>
  !> Arguments:
  !>   stiff  (in)  real(rk)(6,6)  crystal-frame stiffness [MPa]
  subroutine check_iteration_limit(stiff)

    real(rk), intent(in) :: stiff(6, 6)

    real(rk) :: rot(3, 3), defrate(6), spin(3), strain_prev(6)
    real(rk) :: stress(6), strain(6), spin_lat(3), secant(6, 6), vec_h(6), rot_end(3, 3)
    integer :: status

    call random_point(rot, defrate, spin, strain_prev)
    call elastic_point_solve(stiff, rot, strain_prev, defrate, spin, 0.1_rk, 1e-5_rk, 1, stress, &
      & strain, spin_lat, rot_end, secant, vec_h, status)
    call assert_equal(status, 1, 'iteration limit: status 1')

  end subroutine check_iteration_limit

  !> The contract: state allocation, dispatch through class(material_model_t), element lists,
  !> tangent, failure count.
  !>
  !> Arguments:
  !>   stiff  (in)  real(rk)(6,6)  crystal-frame stiffness [MPa]
  !> Locals:
  !>   model            class(material_model_t), allocatable  the model, used polymorphically
  !>   state_prev, state  type(material_state_t)  states
  !>   defrate, spin    real(rk)(:,:,:)   per-point inputs
  !>   secant, tangent  real(rk)(:,:,:,:) per-point matrices
  !>   vec_h            real(rk)(:,:,:)   per-point vectors
  !>   num_failed       integer(ik)       failed points
  !>   e, q             integer           element, point
  subroutine check_contract(stiff)

    real(rk), intent(in) :: stiff(6, 6)

    class(material_model_t), allocatable :: model
    type(elastic_model_t) :: elastic
    type(material_state_t) :: state_prev, state
    real(rk), allocatable :: defrate(:, :, :), spin(:, :, :), secant(:, :, :, :), &
      & tangent(:, :, :, :), vec_h(:, :, :)
    integer(ik) :: num_failed
    integer :: e, q
    real(rk) :: rot(3, 3), d1(6), w1(3), e1(6), stress1(6), strain1(6), spin1(3), sec1(6, 6), &
      & h1(6), rot1(3, 3)
    integer :: status

    ! State allocation.
    call material_state_alloc(state_prev, 15_ik, 4_ik, 0_ik)
    call assert_equal(size(state_prev%stress_s, 2), 15, 'state qpts')
    call assert_equal(size(state_prev%stress_s, 3), 4, 'state elts')
    call assert_equal(size(state_prev%crss, 1), 0, 'state crss has 0 rows (elastic)')
    call assert_true(all(abs(state_prev%strain_el_s) < tiny(1.0_rk)), 'state zero strain')
    call assert_close(reshape(state_prev%rot_c2s(:, :, 7, 3), [9]), [1.0_rk, 0.0_rk, 0.0_rk, &
      & 0.0_rk, 1.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 1.0_rk], 'state identity orientation')

    ! Random start state and inputs.
    allocate (defrate(6, 15, 4), spin(3, 15, 4), secant(6, 6, 15, 4), tangent(6, 6, 15, 4), &
      & vec_h(6, 15, 4))
    do e = 1, 4
      do q = 1, 15
        call random_point(rot, defrate(:, q, e), spin(:, q, e), state_prev%strain_el_s(:, q, e))
        state_prev%rot_c2s(:, :, q, e) = rot
      end do
    end do
    state = state_prev
    secant = -1.0_rk
    tangent = -1.0_rk
    vec_h = -1.0_rk

    ! Polymorphic model; solve the points of elements 2 and 4 only.
    call elastic_model_init(elastic, stiff)
    allocate (model, source=elastic)
    call model%init_state([2_ik, 4_ik], state)
    call model%solve_points([2_ik, 4_ik], 0.1_rk, defrate, spin, state_prev, state, secant, &
      & vec_h, num_failed, tangent)
    call assert_equal(num_failed, 0_ik, 'contract: no failed points')

    ! Elements 1 and 3 untouched; elements 2 and 4 equal the point kernel; tangent = secant.
    call assert_true(all(abs(secant(:, :, :, [1, 3]) + 1.0_rk) < tiny(1.0_rk)), &
      & 'contract: other elements untouched')
    call assert_close(reshape(state%strain_el_s(:, :, [1, 3]), [180]), &
      & reshape(state_prev%strain_el_s(:, :, [1, 3]), [180]), 'contract: other states untouched')
    call assert_close(reshape(tangent(:, :, :, [2, 4]), [1080]), &
      & reshape(secant(:, :, :, [2, 4]), [1080]), 'contract: tangent = secant')
    do e = 2, 4, 2
      do q = 1, 15
        rot = state_prev%rot_c2s(:, :, q, e)
        d1 = defrate(:, q, e)
        w1 = spin(:, q, e)
        e1 = state_prev%strain_el_s(:, q, e)
        call elastic_point_solve(stiff, rot, e1, d1, w1, 0.1_rk, 1e-5_rk, 100, stress1, strain1, &
          & spin1, rot1, sec1, h1, status)
        call assert_close(state%stress_s(:, q, e), stress1, 'contract: stress = kernel')
        call assert_close(reshape(state%rot_c2s(:, :, q, e), [9]), reshape(rot1, [9]), &
          & 'contract: orientation = kernel')
        call assert_close(vec_h(:, q, e), h1, 'contract: {h} = kernel')
      end do
    end do

    ! Failure count: a one-iteration limit fails every point with a nonzero spin.
    call elastic_model_init(elastic, stiff, max_iters_state=1)
    deallocate (model)
    allocate (model, source=elastic)
    call model%solve_points([1_ik, 3_ik], 0.1_rk, defrate, spin, state_prev, state, secant, &
      & vec_h, num_failed)
    call assert_equal(num_failed, 30_ik, 'contract: num_failed counts the points')

  end subroutine check_contract

end program test_elastic_model
