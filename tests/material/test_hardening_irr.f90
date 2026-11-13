! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of the irradiation damage term of hardening_mod (step 7.3, D61,
!> dev/theory/addendum_irradiation.md §10 items 1-3):
!>   - the average loop tensor (1/N) sum (I - n (x) n) against explicit 3x3 algebra for random
!>     normals, and (2/3) I for the four <111> and the three <100> normals;
!>   - N^a = n^a (x) n^a in Mandel form: N^a . H = n^T H n for random symmetric H;
!>   - H_0 (I4) with Barton's loops (1.63e22 m^-3, d = 100 b): 8.08e8 I [1/mm^2], and the defect
!>     term b mu_0 sqrt(h_d N:H_0) = 498.2 MPa with Barton's mu_0 = 86 GPa (addendum §5);
!>   - strengths (I1): equal to Kocks-Mecking for H = 0; the defect term 0 where N^a : H < 0
!>     (<x>_+, errata T10);
!>   - the rate (I3) against explicit 3x3 algebra, and its Jacobian (J H = H dot);
!>   - forward Euler and its guard (dt eta Gammadot >= 1 gives H_prev);
!>   - backward Euler: residual at round-off after one Newton iteration, slip increments up to
!>     eta Delta Gamma = 10; H_prev = 0 gives 0 with no iteration;
!>   - one active system from isotropic H_0: N^a : H to exp(-eta gamma) and a latent plane to
!>     1 - (1 - exp(-eta gamma)) (n^a . n^b)^2 at first order in dt (errors halve); the active
!>     contact density decreases monotonically;
!>   - constant slip-rate ratios (three active systems): H = exp(-eta t M) H_0 with
!>     M = sum |gammadot^a| N^a N^a^T, from the eigen-decomposition of M (cyclic Jacobi, checked:
!>     M V = V Lambda, V^T V = I), at first order in dt (errors halve) (addendum §10 item 2).
!>
!> Locals: see the declarations (parameters, normals, tensors, errors).
program test_hardening_irr
  use kinds_mod, only: rk
  use tensor_mod, only: mandel_from_matrix, matrix_from_mandel
  use hardening_mod, only: hardening_params_t, hardening_irr_t, hardening_kocks_mecking, hardening_km_strengths, &
    & hardening_irr_loop_average, hardening_irr_planes, hardening_irr_dmg_0, &
    & hardening_irr_strengths, hardening_irr_rate, hardening_irr_forward, hardening_irr_backward, &
    & hardening_work_t, hardening_work_alloc
  use crystal_mod, only: crystal_loop_normals, crystal_fcc, crystal_bcc, crystal_hcp
  use tensor_test_utils_mod, only: seed_random, symmetric_eigen
  use assert_mod, only: assert_true, assert_close, test_summary

  implicit none

  integer, parameter :: num = 3
  real(rk), parameter :: dt_list(4) = [1.0e-4_rk, 1.0e-2_rk, 2.0e-2_rk, 1.0e-1_rk]
  type(hardening_params_t) :: hard
  type(hardening_work_t) :: work
  type(hardening_irr_t) :: irr
  real(rk) :: normals(3, 5), avg(6), avg_ref(3, 3), eye(3, 3), mat(3, 3), sym(3, 3), unit(3)
  real(rk) :: slip_normal(3, num), plane(6, num), dmg(6), dmg_prev(6), dmg_est(6), rate(6)
  real(rk) :: rate_ref(6), jac(6, 6), tau_0(num), crss(num), sliprate(num), errors_a(5)
  real(rk), allocatable :: family(:, :)
  real(rk) :: errors_b(5), dt, gamma_end, cos2, contact, contact_prev, value, residual(6)
  logical :: conv, monotonic
  real(rk) :: rate_mat(6, 6), eig_values(6), eig_vectors(6, 6), dmg_exact(6), eye6(6, 6), t_end
  integer :: i, k, a, iters, num_incr, level

  call seed_random(20261012)
  eye = 0.0_rk
  do k = 1, 3
    eye(k, k) = 1.0_rk
  end do

  ! Average loop tensor against explicit 3x3 algebra, random normals of any length.
  call random_number(normals)
  normals = normals - 0.5_rk
  avg_ref = 0.0_rk
  do i = 1, 5
    unit = normals(:, i)/norm2(normals(:, i))
    do k = 1, 3
      avg_ref(:, k) = avg_ref(:, k) + eye(:, k) - unit*unit(k)
    end do
  end do
  avg_ref = avg_ref/5.0_rk
  avg = hardening_irr_loop_average(normals)
  call assert_close(matrix_from_mandel(avg), avg_ref, 'loop average, random normals', &
    & abs_tol=1.0e-14_rk)
  call assert_close(matrix_from_mandel(hardening_irr_loop_average(reshape([1, 1, 1, 1, 1, -1, &
    & 1, -1, 1, -1, 1, 1]*1.0_rk, [3, 4]))), 2.0_rk/3.0_rk*eye, 'loop average <111> = 2/3 I', &
    & abs_tol=1.0e-15_rk)
  call assert_close(matrix_from_mandel(hardening_irr_loop_average(reshape([1, 0, 0, 0, 1, 0, &
    & 0, 0, 1]*1.0_rk, [3, 3]))), 2.0_rk/3.0_rk*eye, 'loop average <100> = 2/3 I', &
    & abs_tol=1.0e-15_rk)

  ! Loop-normal families (crystal_mod): 4, 3, 6 unit normals, distinct up to sign, each family
  ! averaging to (2/3) I; unknown tokens and hcp give none.
  call crystal_loop_normals(crystal_bcc, '<111>', family)
  call assert_true(size(family, 2) == 4, 'loop normals <111>: 4')
  call check_family(family, '<111>')
  call crystal_loop_normals(crystal_fcc, '<100>', family)
  call assert_true(size(family, 2) == 3, 'loop normals <100>: 3')
  call check_family(family, '<100>')
  call crystal_loop_normals(crystal_bcc, '<110>', family)
  call assert_true(size(family, 2) == 6, 'loop normals <110>: 6')
  call check_family(family, '<110>')
  call crystal_loop_normals(crystal_bcc, '<112>', family)
  call assert_true(size(family, 2) == 0, 'loop normals: unknown token')
  call crystal_loop_normals(crystal_hcp, '<111>', family)
  call assert_true(size(family, 2) == 0, 'loop normals: hcp not defined')

  ! Slip planes (110), (101), (112) (bcc), and N^a . H = n^T H n for a random symmetric H.
  slip_normal(:, 1) = [1.0_rk, 1.0_rk, 0.0_rk]/sqrt(2.0_rk)
  slip_normal(:, 2) = [1.0_rk, 0.0_rk, 1.0_rk]/sqrt(2.0_rk)
  slip_normal(:, 3) = [1.0_rk, 1.0_rk, 2.0_rk]/sqrt(6.0_rk)
  plane = hardening_irr_planes(slip_normal)
  call random_number(sym)
  sym = sym + transpose(sym)
  do a = 1, num
    call assert_close(dot_product(plane(:, a), mandel_from_matrix(sym)), &
      & dot_product(slip_normal(:, a), matmul(sym, slip_normal(:, a))), &
      & 'N^a : H = n^T H n', rel_tol=1.0e-14_rk)
  end do

  ! Barton's parameters in FEPX units (Table 1, §3.2), mu_0 = 86 GPa for the check of §5.
  hard%law = hardening_kocks_mecking
  hard%km_b = 2.48e-7_rk
  hard%km_mu_0 = 86000.0_rk
  hard%km_h_n = 0.125_rk
  irr%h_d = 0.675_rk
  irr%eta = 100.0_rk
  irr%loop_avg_c = hardening_irr_loop_average(reshape([1, 1, 1, 1, 1, -1, 1, -1, 1, -1, &
    & 1, 1]*1.0_rk, [3, 4]))
  irr%plane_c = plane
  hard%tol_hardening = 1.0e-4_rk
  hard%max_iters = 10
  call hardening_work_alloc(work, hard, num, irradiation=.true.)
  tau_0 = [90.0_rk, 95.0_rk, 100.0_rk]

  ! H_0 (I4) and the initial defect term (addendum §5).
  dmg = hardening_irr_dmg_0(irr, 1.63e13_rk, 100.0_rk*hard%km_b)
  call assert_close(matrix_from_mandel(dmg), 1.63e13_rk*3.0_rk*2.48e-5_rk*2.0_rk/3.0_rk*eye, &
    & 'H_0 Barton', rel_tol=1.0e-14_rk)
  crss = hardening_irr_strengths(hard, irr, tau_0, 2.0e7_rk, dmg) - &
    & hardening_km_strengths(hard, tau_0, 2.0e7_rk)
  call assert_close(crss, spread(hard%km_b*hard%km_mu_0*sqrt(0.675_rk*8.0848e8_rk), 1, num), &
    & 'defect term at H_0', rel_tol=1.0e-12_rk)
  call assert_close(crss(1), 498.2_rk, 'defect term about 498 MPa (Barton Fig. 3: ~480)', &
    & abs_tol=0.1_rk)

  ! H = 0: Kocks-Mecking exactly; N^a : H < 0 on a system: its defect term is 0.
  call assert_close(hardening_irr_strengths(hard, irr, tau_0, 2.0e7_rk, spread(0.0_rk, 1, 6)), &
    & hardening_km_strengths(hard, tau_0, 2.0e7_rk), 'H = 0: Kocks-Mecking', rel_tol=0.0_rk)
  sym = 0.0_rk
  sym(3, 3) = -1.0e8_rk
  sym(1, 1) = 1.0e8_rk
  crss = hardening_irr_strengths(hard, irr, tau_0, 2.0e7_rk, mandel_from_matrix(sym)) - &
    & hardening_km_strengths(hard, tau_0, 2.0e7_rk)
  call assert_true(dot_product(plane(:, 3), mandel_from_matrix(sym)) < 0.0_rk .and. &
    & .not. abs(crss(3)) > 0.0_rk, 'negative contact density: no defect term')
  call assert_close(crss(1), hard%km_b*hard%km_mu_0*sqrt(0.675_rk*0.5e8_rk), &
    & 'positive contact density kept', rel_tol=1.0e-13_rk)

  ! Rate (I3) against explicit 3x3 algebra; J H = H dot.
  call random_number(sym)
  sym = 1.0e8_rk*(sym + transpose(sym))
  sliprate = [0.3_rk, -0.7_rk, 0.0_rk]
  call hardening_irr_rate(irr, sliprate, mandel_from_matrix(sym), rate, jac)
  mat = 0.0_rk
  do a = 1, num
    mat = mat - irr%eta*abs(sliprate(a))*dot_product(slip_normal(:, a), &
      & matmul(sym, slip_normal(:, a)))*spread(slip_normal(:, a), 2, 3)* &
      & spread(slip_normal(:, a), 1, 3)
  end do
  rate_ref = mandel_from_matrix(mat)
  call assert_close(rate, rate_ref, 'H dot (I3)', abs_tol=1.0e-13_rk*maxval(abs(rate_ref)))
  call assert_close(matmul(jac, mandel_from_matrix(sym)), rate, 'J H = H dot', &
    & abs_tol=1.0e-13_rk*maxval(abs(rate_ref)))
  call assert_close(jac, transpose(jac), 'J symmetric', abs_tol=0.0_rk)

  ! Forward Euler, and its guard.
  dmg_prev = hardening_irr_dmg_0(irr, 1.63e13_rk, 2.48e-5_rk)
  dt = 1.0e-3_rk
  call hardening_irr_forward(irr, dmg_prev, sliprate, dt, dmg_est)
  call hardening_irr_rate(irr, sliprate, dmg_prev, rate)
  call assert_close(dmg_est, dmg_prev + dt*rate, 'forward Euler', rel_tol=1.0e-15_rk)
  call hardening_irr_forward(irr, dmg_prev, sliprate, 1.0_rk/irr%eta, dmg_est)
  call assert_close(dmg_est, dmg_prev, 'forward Euler guard: H_prev', rel_tol=0.0_rk)

  ! Backward Euler: one Newton iteration, residual at round-off, slip increments to 0.1.
  do level = 1, 4
    dt = dt_list(level)
    call hardening_irr_forward(irr, dmg_prev, sliprate, dt, dmg)
    call hardening_irr_backward(hard, irr, dmg_prev, sliprate, dt, dmg, conv, work, iters)
    call hardening_irr_rate(irr, sliprate, dmg, rate)
    residual = dmg - dmg_prev - dt*rate
    call assert_true(conv .and. iters == 1, 'backward Euler: one iteration')
    call assert_true(maxval(abs(residual)) < 1.0e-12_rk*maxval(abs(dmg_prev)), &
      & 'backward Euler: residual at round-off')
  end do
  dmg = 1.0_rk
  call hardening_irr_backward(hard, irr, spread(0.0_rk, 1, 6), sliprate, dt, dmg, conv, work, &
    & iters)
  call assert_true(conv .and. iters == 0 .and. .not. any(abs(dmg) > 0.0_rk), &
    & 'no loops: H stays 0')

  ! One active system (110) at slip rate 1/s from isotropic H_0 to gamma = 0.02: N^1 : H against
  ! exp(-eta gamma), the latent (112) plane against 1 - (1 - exp(-eta gamma)) cos^2.
  sliprate = [1.0_rk, 0.0_rk, 0.0_rk]
  gamma_end = 0.02_rk
  cos2 = dot_product(slip_normal(:, 1), slip_normal(:, 3))**2
  value = dot_product(plane(:, 1), dmg_prev)
  monotonic = .true.
  do level = 1, 5
    num_incr = 10*2**(level - 1)
    dt = gamma_end/num_incr
    dmg = dmg_prev
    contact_prev = value
    do i = 1, num_incr
      dmg_est = dmg
      call hardening_irr_forward(irr, dmg, sliprate, dt, dmg_est)
      call hardening_irr_backward(hard, irr, dmg, sliprate, dt, dmg_est, conv, work)
      dmg = dmg_est
      contact = dot_product(plane(:, 1), dmg)
      if (.not. contact < contact_prev) monotonic = .false.
      contact_prev = contact
    end do
    errors_a(level) = abs(dot_product(plane(:, 1), dmg)/value - exp(-irr%eta*gamma_end))
    errors_b(level) = abs(dot_product(plane(:, 3), dmg)/value - (1.0_rk - (1.0_rk - &
      & exp(-irr%eta*gamma_end))*cos2))
  end do
  call assert_true(monotonic, 'one active system: contact density decreases monotonically')
  do level = 2, 5
    call assert_close(errors_a(level - 1)/errors_a(level), 2.0_rk, &
      & 'active plane: first order', abs_tol=0.15_rk)
    call assert_close(errors_b(level - 1)/errors_b(level), 2.0_rk, &
      & 'latent plane: first order', abs_tol=0.15_rk)
  end do
  call assert_true(errors_a(5) < 1.0e-2_rk .and. errors_b(5) < 1.0e-2_rk, &
    & 'finest level close to the analytic solution')

  ! Constant slip-rate ratios, three active systems: H(t) = exp(-eta t M) H_0 with M the
  ! (positive semidefinite) sum of |gammadot^a| N^a N^a^T, from the eigen-decomposition of M.
  sliprate = [1.0_rk, -0.5_rk, 0.25_rk]
  t_end = 0.02_rk
  rate_mat = 0.0_rk
  do a = 1, num
    rate_mat = rate_mat + abs(sliprate(a))*spread(plane(:, a), 2, 6)*spread(plane(:, a), 1, 6)
  end do
  call symmetric_eigen(rate_mat, eig_values, eig_vectors)
  eye6 = 0.0_rk
  do k = 1, 6
    eye6(k, k) = 1.0_rk
  end do
  call assert_close(matmul(rate_mat, eig_vectors), matmul(eig_vectors, &
    & eye6*spread(eig_values, 1, 6)), 'Jacobi: M V = V Lambda', abs_tol=1.0e-14_rk)
  call assert_close(matmul(transpose(eig_vectors), eig_vectors), eye6, 'Jacobi: V^T V = I', &
    & abs_tol=1.0e-14_rk)
  dmg_exact = matmul(eig_vectors, exp(-irr%eta*t_end*eig_values)* &
    & matmul(transpose(eig_vectors), dmg_prev))
  do level = 1, 5
    num_incr = 10*2**(level - 1)
    dt = t_end/num_incr
    dmg = dmg_prev
    do i = 1, num_incr
      dmg_est = dmg
      call hardening_irr_backward(hard, irr, dmg, sliprate, dt, dmg_est, conv, work)
      dmg = dmg_est
    end do
    errors_a(level) = maxval(abs(dmg - dmg_exact))/maxval(abs(dmg_prev))
  end do
  do level = 2, 5
    call assert_close(errors_a(level - 1)/errors_a(level), 2.0_rk, &
      & 'constant ratios: first order', abs_tol=0.15_rk)
  end do
  call assert_true(errors_a(5) < 1.0e-2_rk, &
    & 'constant ratios: finest level close to exp(-eta t M) H_0')

  call test_summary()

contains

  !> A loop-normal family: unit normals, no two parallel, average loop tensor (2/3) I.
  !>
  !> Arguments:
  !>   normals_c  (in)  real(rk)(3, :)  the family's normals
  !>   label      (in)  character       family token
  !> Locals:
  !>   i, j  integer  normals
  subroutine check_family(normals_c, label)

    real(rk), intent(in) :: normals_c(:, :)
    character(len=*), intent(in) :: label

    integer :: i, j

    do i = 1, size(normals_c, 2)
      call assert_close(norm2(normals_c(:, i)), 1.0_rk, 'unit normal ' // label, &
        & rel_tol=1.0e-15_rk)
      do j = 1, i - 1
        call assert_true(abs(dot_product(normals_c(:, i), normals_c(:, j))) < 0.99_rk, &
          & 'distinct normals ' // label)
      end do
    end do
    call assert_close(matrix_from_mandel(hardening_irr_loop_average(normals_c)), &
      & 2.0_rk/3.0_rk*eye, 'loop average (2/3) I ' // label, abs_tol=1.0e-15_rk)

  end subroutine check_family

end program test_hardening_irr
