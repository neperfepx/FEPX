! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Unit tests of the output quantities of the material state (steps 6.5, 6.5b):
!> material_state_mod::material_state_update_derived (point kinematics, total and sample-frame
!> plastic strain, work rate, work, plastic work rate, plastic work, lattice rotation rate) and
!> the slip quantities set by crystal_plasticity_model_mod's solve_points.
!>
!> Source : FEPX (ROADMAP steps 6.5, 6.5b).
!>
!> Checks (fixed random seed; 2 elements x 3 points; random states):
!>   - derived quantities: the work rate equals sigma : d with explicit 3x3 tensors (the full
!>     contraction), the work equals work_prev + dt/2 (workrate_prev + workrate), and the
!>     rotation rate of rot = exp(W dt) rot_prev (W a random sample-frame spin, angles up to
!>     0.3 rad) equals W (relative 1e-12); step 6.5b: d and w are stored as given, the total
!>     strain equals strain_prev + dt d and the sample-frame plastic strain strain_pl_prev +
!>     dt d^p (explicit 3x3 sums), the plastic work rate equals sigma' : d^p with the explicit
!>     deviatoric stress (v2 finalize_res_plasticwork; d^p a random deviatoric tensor), and the
!>     plastic work follows the trapezoidal rule (relative 1e-12);
!>   - crystal plasticity (fcc and hcp, Cu and Ti constants, strengths at the stress level so
!>     that systems slip, random start slip and plastic strain): after solve_points, at every
!>     point, the slip rates equal those of cp_point_solve on the same inputs; the resolved
!>     shear stresses equal s^a . tau m^a with tau = beta R^T sigma R (explicit 3x3, s and m the
!>     slip direction and normal of crystal_t), and the power law of slip_kinetics_mod gives the
!>     slip rates from them; slip = slip_prev + gammadot dt; the crystal-frame plastic strain
!>     equals its start value plus dt sum gammadot^a sym(s^a (x) m^a); the plastic deformation
!>     rate and spin equal R (sum gammadot^a sym / skw(s^a (x) m^a)) R^T (explicit 3x3); the
!>     rotation rate of material_state_update_derived equals the w^p minus the slip spin of the
!>     orientation update (spin_hat_update_s - spin_slip_update_s, step 6.5b); at least
!>     one system slips (the check is not vacuous); relative 1e-10 of the largest entry.
!>
!> Locals:
!>   work   type(cp_work_t)  point-solve work arrays (cp_work_alloc before each call; D100)
!>   crys     type(crystal_t)  slip geometry
!>   stiff    real(rk)(6, 6)   stiffness [MPa]
program test_material_outputs

  use kinds_mod, only: rk, ik
  use tensor_mod, only: matrix_from_mandel, mandel_from_matrix, skew_from_matrix, &
    & matrix_from_skew
  use rotation_mod, only: rotation_exp_skew
  use material_state_mod, only: material_state_t, material_state_alloc, &
    & material_state_update_derived
  use crystal_mod, only: crystal_t, crystal_init, crystal_fcc, crystal_hcp, stiffness_cubic, &
    & stiffness_hexagonal
  use crystal_plasticity_model_mod, only: cp_model_t, cp_model_init, cp_point_solve
  use crystal_plasticity_model_mod, only: cp_work_t, cp_work_alloc
  use slip_kinetics_mod, only: slip_power_law
  use tensor_test_utils_mod, only: seed_random, random_rotation, random_symmetric
  use assert_mod, only: assert_true, assert_close, test_summary

  implicit none

  type(cp_work_t) :: work
  type(crystal_t) :: crys
  real(rk) :: stiff(6, 6)

  call seed_random(20261011)

  call check_derived()

  call crystal_init(crys, crystal_fcc)
  stiff = stiffness_cubic(168.4e3_rk, 121.4e3_rk, 75.4e3_rk)
  call check_cp_outputs(stiff, 100.0_rk, 'fcc')

  call crystal_init(crys, crystal_hcp, c_over_a=1.587_rk)
  stiff = stiffness_hexagonal(162.4e3_rk, 92.0e3_rk, 69.0e3_rk, 180.7e3_rk, 46.7e3_rk)
  call check_cp_outputs(stiff, 150.0_rk, 'hcp')

  call test_summary()

contains

  !> Derived quantities of material_state_update_derived (header).
  !>
  !> Locals:
  !>   state_prev, state  type(material_state_t)  start and current states (2 elements,
  !>                                              3 points)
  !>   defrate_s     real(rk)(6, 3, 2)   {d} [1/s]
  !>   spin          real(rk)(3, 3, 2)   random sample-frame spins of the rotation [1/s]
  !>   spin_s        real(rk)(3, 3, 2)   random continuum spins w (stored as given) [1/s]
  !>   dt            real(rk)            time increment [s]
  !>   i_elt, i_qpt  integer(ik)         element, point
  !>   sk            real(rk)(3)         random numbers
  !>   expected      real(rk)            sigma : d by explicit 3x3 algebra [MPa/s]
  !>   expected_pl   real(rk)            sigma' : d^p by explicit 3x3 algebra [MPa/s]
  !>   sig_dev, dp   real(rk)(3, 3)      deviatoric stress [MPa], d^p [1/s]
  !>   mat_expected  real(rk)(3, 3)      expected strain (3x3) [-]
  subroutine check_derived()

    type(material_state_t) :: state_prev, state
    real(rk) :: defrate_s(6, 3, 2), spin(3, 3, 2), spin_s(3, 3, 2), dt, sk(3), expected
    real(rk) :: expected_pl, sig_dev(3, 3), dp(3, 3), mat_expected(3, 3)
    integer(ik) :: i_elt, i_qpt, i

    call material_state_alloc(state_prev, 3_ik, 2_ik, 0_ik)
    dt = 0.7_rk
    do i_elt = 1_ik, 2_ik
      do i_qpt = 1_ik, 3_ik
        ! Random start: rotation, work and work rate; current stress and {d}; a spin.
        state_prev%rot_c2s(:, :, i_qpt, i_elt) = random_rotation()
        call random_number(sk)
        state_prev%work(i_qpt, i_elt) = 10.0_rk*sk(1)
        state_prev%workrate(i_qpt, i_elt) = sk(2) - 0.5_rk
        state_prev%work_pl(i_qpt, i_elt) = 5.0_rk*sk(3)
        state_prev%workrate_pl(i_qpt, i_elt) = sk(1)*sk(2)
        state_prev%strain_s(:, i_qpt, i_elt) = 1.0e-2_rk*mandel_from_matrix(random_symmetric())
        state_prev%strain_pl_s(:, i_qpt, i_elt) = 1.0e-2_rk* &
          & mandel_from_matrix(random_symmetric())
        call random_number(spin(:, i_qpt, i_elt))
        spin(:, i_qpt, i_elt) = 0.3_rk/dt*(2.0_rk*spin(:, i_qpt, i_elt) - 1.0_rk)/sqrt(3.0_rk)
      end do
    end do
    state = state_prev
    do i_elt = 1_ik, 2_ik
      do i_qpt = 1_ik, 3_ik
        state%stress_s(:, i_qpt, i_elt) = 100.0_rk*mandel_from_matrix(random_symmetric())
        defrate_s(:, i_qpt, i_elt) = 1.0e-3_rk*mandel_from_matrix(random_symmetric())
        state%rot_c2s(:, :, i_qpt, i_elt) = matmul(rotation_exp_skew(dt* &
          & spin(:, i_qpt, i_elt)), state_prev%rot_c2s(:, :, i_qpt, i_elt))
        ! A deviatoric d^p (set by the models) and a continuum spin.
        dp = 1.0e-3_rk*random_symmetric()
        sk(1) = (dp(1, 1) + dp(2, 2) + dp(3, 3))/3.0_rk
        do i = 1, 3
          dp(i, i) = dp(i, i) - sk(1)
        end do
        state%defrate_pl_s(:, i_qpt, i_elt) = mandel_from_matrix(dp)
        call random_number(spin_s(:, i_qpt, i_elt))
        spin_s(:, i_qpt, i_elt) = 1.0e-3_rk*(2.0_rk*spin_s(:, i_qpt, i_elt) - 1.0_rk)
      end do
    end do

    call material_state_update_derived(state_prev, defrate_s, spin_s, dt, state)

    do i_elt = 1_ik, 2_ik
      do i_qpt = 1_ik, 3_ik
        ! sigma : d, full contraction of the 3x3 tensors.
        expected = sum(matrix_from_mandel(state%stress_s(:, i_qpt, i_elt))* &
          & matrix_from_mandel(defrate_s(:, i_qpt, i_elt)))
        call assert_close(state%workrate(i_qpt, i_elt), expected, 'derived_workrate', &
          & rel_tol=1.0e-12_rk)
        call assert_close(state%work(i_qpt, i_elt), state_prev%work(i_qpt, i_elt) + 0.5_rk* &
          & dt*(state_prev%workrate(i_qpt, i_elt) + expected), 'derived_work', rel_tol=1.0e-12_rk)

        ! Rotation rate = the spin of the update.
        call assert_close(state%rotrate_s(:, i_qpt, i_elt), spin(:, i_qpt, i_elt), &
          & 'derived_rotrate', abs_tol=1.0e-12_rk*maxval(abs(spin(:, i_qpt, i_elt))))

        ! Kinematics stored as given.
        call assert_close(state%defrate_s(:, i_qpt, i_elt), defrate_s(:, i_qpt, i_elt), &
          & 'derived_defrate', abs_tol=0.0_rk)
        call assert_close(state%spin_s(:, i_qpt, i_elt), spin_s(:, i_qpt, i_elt), &
          & 'derived_spin', abs_tol=0.0_rk)

        ! Total strain and sample-frame plastic strain: 3x3 sums of the rates times dt.
        mat_expected = matrix_from_mandel(state_prev%strain_s(:, i_qpt, i_elt)) + &
          & dt*matrix_from_mandel(defrate_s(:, i_qpt, i_elt))
        call assert_close(matrix_from_mandel(state%strain_s(:, i_qpt, i_elt)), mat_expected, &
          & 'derived_strain', abs_tol=1.0e-12_rk*maxval(abs(mat_expected)))
        dp = matrix_from_mandel(state%defrate_pl_s(:, i_qpt, i_elt))
        mat_expected = matrix_from_mandel(state_prev%strain_pl_s(:, i_qpt, i_elt)) + dt*dp
        call assert_close(matrix_from_mandel(state%strain_pl_s(:, i_qpt, i_elt)), &
          & mat_expected, 'derived_strain_pl_sample', &
          & abs_tol=1.0e-12_rk*maxval(abs(mat_expected)))

        ! Plastic work rate sigma' : d^p (explicit deviatoric stress) and plastic work.
        sig_dev = matrix_from_mandel(state%stress_s(:, i_qpt, i_elt))
        expected = (sig_dev(1, 1) + sig_dev(2, 2) + sig_dev(3, 3))/3.0_rk
        do i = 1, 3
          sig_dev(i, i) = sig_dev(i, i) - expected
        end do
        expected_pl = sum(sig_dev*dp)
        call assert_close(state%workrate_pl(i_qpt, i_elt), expected_pl, 'derived_workrate_pl', &
          & rel_tol=1.0e-12_rk)
        call assert_close(state%work_pl(i_qpt, i_elt), state_prev%work_pl(i_qpt, i_elt) + &
          & 0.5_rk*dt*(state_prev%workrate_pl(i_qpt, i_elt) + expected_pl), 'derived_work_pl', &
          & rel_tol=1.0e-12_rk)
      end do
    end do

  end subroutine check_derived

  !> Slip quantities set by the crystal plasticity solve_points (header).
  !>
  !> Arguments:
  !>   stiff_c  (in)  real(rk)(6, 6)  stiffness [MPa]
  !>   crss     (in)  real(rk)        strength of every system [MPa]
  !>   name     (in)  character       label
  !> Locals:
  !>   model        type(cp_model_t)        model (m 0.02, gammadot_0 1)
  !>   state_prev, state  type(material_state_t)  start and current states
  !>   n            integer                 slip systems
  !>   elts         integer(ik)(2)          elements of the phase
  !>   defrate_s, spin_s  real(rk)          {d} (6, 3, 2) and w (3, 3, 2) [1/s]
  !>   sec, h, tan  real(rk)                [s], {h}, [s_tan] of every point
  !>   dt           real(rk)                time increment [s]
  !>   num_failed   integer(ik)             failed point solves
  !>   i_elt, i_qpt integer(ik)             element, point
  !>   a, status    integer                 system, point-solve status
  !>   tau_c        real(rk)(6)             start Kirchhoff stress, crystal frame [MPa]
  !>   e_c          real(rk)(6)             its elastic strain [-]
  !>   rot          real(rk)(3, 3)          end orientation
  !>   tau_mat      real(rk)(3, 3)          end Kirchhoff stress beta R^T sigma R, crystal [MPa]
  !>   schmid       real(rk)(3, 3)          s (x) m of a system
  !>   dp_mat, wp_mat  real(rk)(3, 3)       sum gammadot sym / skw(s (x) m), crystal [1/s]
  !>   sig_n, e_n, wp  real(rk)             point-solve outputs (sliprate reference)
  !>   rot_n        real(rk)(3, 3)          point-solve orientation
  !>   rate_ref     real(rk)(:)             slip rates of cp_point_solve [1/s]
  !>   crss_n       real(rk)(:)             strengths of cp_point_solve [MPa]
  !>   rate_law, coef_sec, coef_tan  real(rk)(:)  power law of the resolved shear stresses
  !>   rss_ref      real(rk)(:)             s.tau m [MPa]
  !>   sec1, tan1   real(rk)(6, 6)          point matrices of cp_point_solve
  !>   h1           real(rk)(6)             point vector of cp_point_solve
  !>   beta         real(rk)                det(I + e^e) [-]
  !>   slipping     logical                 a system slips somewhere
  subroutine check_cp_outputs(stiff_c, crss, name)

    real(rk), intent(in) :: stiff_c(6, 6)
    real(rk), intent(in) :: crss
    character(len=*), intent(in) :: name

    type(cp_model_t) :: model
    type(material_state_t) :: state_prev, state
    integer :: n, a, status
    integer(ik) :: elts(2), num_failed, i_elt, i_qpt
    real(rk) :: defrate_s(6, 3, 2), spin_s(3, 3, 2), sec(6, 6, 3, 2), h(6, 3, 2), tan(6, 6, 3, 2)
    real(rk) :: dt, tau_c(6), e_c(6), rot(3, 3), tau_mat(3, 3), schmid(3, 3), dp_mat(3, 3)
    real(rk) :: wp_mat(3, 3), sig_n(6), e_n(6), wp(3), rot_n(3, 3), sec1(6, 6), tan1(6, 6)
    real(rk) :: h1(6), beta
    real(rk), allocatable :: rate_ref(:), crss_n(:), rate_law(:), coef_sec(:), coef_tan(:)
    real(rk), allocatable :: rss_ref(:)
    logical :: slipping

    ! Model: every system at strength crss, m 0.02, gammadot_0 1, tight tolerances.
    n = crys%num_slip
    call cp_model_init(model, stiff_c, crys%schmid_sym_c, crys%schmid_skew_c, &
      & spread(0.02_rk, 1, n), spread(1.0_rk, 1, n), spread(crss, 1, n), &
      & tol_newton=1.0e-12_rk, &
      & max_iters_newton=100, max_iters_state=200, tol_state=1.0e-14_rk)
    allocate (rate_ref(n), crss_n(n), rate_law(n), coef_sec(n), coef_tan(n), rss_ref(n))

    ! Start states: random orientation, a stress whose largest resolved shear stress is crss, its
    ! elastic strain, random accumulated slip and plastic strain; random {d} and w.
    call material_state_alloc(state_prev, 3_ik, 2_ik, int(n, ik))
    elts = [1_ik, 2_ik]
    call model%init_state(elts, state_prev)
    dt = 0.2_rk
    do i_elt = 1_ik, 2_ik
      do i_qpt = 1_ik, 3_ik
        state_prev%rot_c2s(:, :, i_qpt, i_elt) = random_rotation()
        tau_c = mandel_from_matrix(random_symmetric())
        tau_c = crss*tau_c/maxval(abs(matmul(tau_c, model%schmid_sym_c)))
        e_c = matmul(model%compl_el_c, tau_c)
        beta = det3(matrix_from_mandel(e_c))
        rot = state_prev%rot_c2s(:, :, i_qpt, i_elt)
        state_prev%stress_s(:, i_qpt, i_elt) = mandel_from_matrix(matmul(rot, &
          & matmul(matrix_from_mandel(tau_c), transpose(rot))))/beta
        state_prev%strain_el_s(:, i_qpt, i_elt) = mandel_from_matrix(matmul(rot, &
          & matmul(matrix_from_mandel(e_c), transpose(rot))))
        call random_number(state_prev%slip(:, i_qpt, i_elt))
        state_prev%slip(:, i_qpt, i_elt) = 1.0e-2_rk*(state_prev%slip(:, i_qpt, i_elt) - 0.5_rk)
        state_prev%strain_pl_c(:, i_qpt, i_elt) = 1.0e-2_rk*mandel_from_matrix( &
          & random_symmetric())
        defrate_s(:, i_qpt, i_elt) = 2.0e-3_rk*mandel_from_matrix(random_symmetric())
        call random_number(spin_s(:, i_qpt, i_elt))
        spin_s(:, i_qpt, i_elt) = 2.0e-3_rk*(spin_s(:, i_qpt, i_elt) - 0.5_rk)
      end do
    end do
    state = state_prev

    ! The phase's point solves.
    call model%solve_points(elts, dt, defrate_s, spin_s, state_prev, state, sec, h, num_failed, &
      & tan)
    call assert_true(num_failed == 0_ik, name // '_cp_outputs_converged')
    ! The derived quantities (rotation rate of the update).
    call material_state_update_derived(state_prev, defrate_s, spin_s, dt, state)

    slipping = .false.
    do i_elt = 1_ik, 2_ik
      do i_qpt = 1_ik, 3_ik
        ! Slip rates of cp_point_solve on the same inputs.
        call cp_work_alloc(work, model)
        call cp_point_solve(model, state_prev%crss(:, i_qpt, i_elt), &
          & state_prev%rot_c2s(:, :, i_qpt, i_elt), state_prev%stress_s(:, i_qpt, i_elt), &
          & state_prev%strain_el_s(:, i_qpt, i_elt), defrate_s(:, i_qpt, i_elt), &
          & spin_s(:, i_qpt, i_elt), dt, sig_n, e_n, wp, rot_n, crss_n, sec1, h1, tan1, &
          & rate_ref, status, work=work)
        call assert_close(state%sliprate(:, i_qpt, i_elt), rate_ref, name // '_sliprate', &
          & abs_tol=1.0e-10_rk*maxval(abs(rate_ref)))
        slipping = slipping .or. maxval(abs(rate_ref)) > 1.0e-6_rk

        ! Resolved shear stresses s . tau m with tau = beta R^T sigma R (explicit 3x3), and the
        ! power law of them gives the slip rates.
        rot = state%rot_c2s(:, :, i_qpt, i_elt)
        beta = det3(matrix_from_mandel(state%strain_el_s(:, i_qpt, i_elt)))
        tau_mat = beta*matmul(transpose(rot), matmul(matrix_from_mandel( &
          & state%stress_s(:, i_qpt, i_elt)), rot))
        do a = 1, n
          rss_ref(a) = dot_product(crys%slip_dir_c(:, a), matmul(tau_mat, &
            & crys%slip_normal_c(:, a)))
        end do
        call assert_close(state%rss(:, i_qpt, i_elt), rss_ref, name // '_rss', &
          & abs_tol=1.0e-10_rk*maxval(abs(rss_ref)))
        call slip_power_law(state%rss(:, i_qpt, i_elt), state%crss(:, i_qpt, i_elt), model%m, &
          & model%gammadot_0, model%rss_cutoff, rate_law, coef_sec, coef_tan)
        where (.not. (abs(rate_law) > 1.0e-30_rk)) rate_law = 0.0_rk
        call assert_close(state%sliprate(:, i_qpt, i_elt), rate_law, name // '_sliprate_law', &
          & abs_tol=1.0e-10_rk*maxval(abs(rate_law)))

        ! Accumulated slip.
        call assert_close(state%slip(:, i_qpt, i_elt), state_prev%slip(:, i_qpt, i_elt) + &
          & dt*state%sliprate(:, i_qpt, i_elt), name // '_slip', rel_tol=1.0e-14_rk)

        ! Plastic deformation rate and spin from the slip directions and normals.
        dp_mat = 0.0_rk
        wp_mat = 0.0_rk
        do a = 1, n
          schmid = spread(crys%slip_dir_c(:, a), 2, 3)*spread(crys%slip_normal_c(:, a), 1, 3)
          dp_mat = dp_mat + state%sliprate(a, i_qpt, i_elt)*0.5_rk*(schmid + transpose(schmid))
          wp_mat = wp_mat + state%sliprate(a, i_qpt, i_elt)*0.5_rk*(schmid - transpose(schmid))
        end do
        call assert_close(state%strain_pl_c(:, i_qpt, i_elt), &
          & state_prev%strain_pl_c(:, i_qpt, i_elt) + dt*mandel_from_matrix(dp_mat), &
          & name // '_strain_pl_c', abs_tol=1.0e-10_rk*maxval(abs(state_prev%strain_pl_c)))
        call assert_close(state%defrate_pl_s(:, i_qpt, i_elt), mandel_from_matrix(matmul(rot, &
          & matmul(dp_mat, transpose(rot)))), name // '_defrate_pl_s', &
          & abs_tol=1.0e-10_rk*maxval(abs(dp_mat)) + tiny(1.0_rk))
        call assert_close(state%spin_pl_s(:, i_qpt, i_elt), skew_from_matrix(matmul(rot, &
          & matmul(wp_mat, transpose(rot)))), name // '_spin_pl_s', &
          & abs_tol=1.0e-10_rk*maxval(abs(wp_mat)) + tiny(1.0_rk))

        ! The two parts of the orientation update give the rotation rate (step 6.5b).
        call assert_close(state%rotrate_s(:, i_qpt, i_elt), &
          & state%spin_hat_update_s(:, i_qpt, i_elt) - state%spin_slip_update_s(:, i_qpt, i_elt), &
          & name // '_rotrate_parts', &
          & abs_tol=1.0e-10_rk*maxval(abs(state%rotrate_s(:, i_qpt, i_elt))) + tiny(1.0_rk))
      end do
    end do
    call assert_true(slipping, name // '_cp_outputs_slipping')

  end subroutine check_cp_outputs

  !> det(I + A) of a 3 x 3 matrix.
  !>
  !> Arguments:
  !>   a  (in)  real(rk)(3,3)  matrix
  !> Locals:
  !>   m  real(rk)(3,3)        I + A
  pure function det3(a) result(det)

    real(rk), intent(in) :: a(3, 3)
    real(rk) :: det

    real(rk) :: m(3, 3)

    m = a
    m(1, 1) = m(1, 1) + 1.0_rk
    m(2, 2) = m(2, 2) + 1.0_rk
    m(3, 3) = m(3, 3) + 1.0_rk
    det = m(1, 1)*(m(2, 2)*m(3, 3) - m(2, 3)*m(3, 2)) - m(1, 2)*(m(2, 1)*m(3, 3) - &
      & m(2, 3)*m(3, 1)) + m(1, 3)*(m(2, 1)*m(3, 2) - m(2, 2)*m(3, 1))

  end function det3

end program test_material_outputs
