! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Driver for addendum_hardening §5 test 3 (time-integration accuracy of the strengths and the
!> stress; step 6.4), extended by addendum_point_verification.md (step 7.4a: lattice spin, the
!> modifiers, bcc, hcp, the Kocks-Mecking limits). A single crystal point under a constant
!> deformation rate {d} and spin w, from a zero state, is integrated by cp_point_solve
!> (staggered forward/backward-Euler strengths in the state loop, orientation by the exponential
!> map) to t = 40 s with dt = 0.8, 0.4, 0.2, 0.1, 0.05 s, with tight point and hardening
!> tolerances, so that the time discretization dominates the error;
!> tests/verification/check_hardening_integration.py integrates the rate equations
!> independently (RK4, the orientation as a Rodrigues vector by arXiv Eq. 31) and compares
!> (first order in dt, error bound).
!>
!> Cases 1-6 (steps 6.4-7.3): fcc (copper), cube orientation, uniaxial
!> {d} = 1e-3 [-1/2, -1/2, 1, 0, 0, 0] 1/s, w = 0 (no lattice spin, by symmetry):
!>   1: isotropic; 2: anisotropic, H off-diagonal 1.4; 3: anisotropic, off-diagonal 5 (h_0 200
!>   MPa, tau_s 330 MPa, n 1.5, tau_0 100 MPa, m 0.02, gammadot_0 1 1/s); 4, 5: Kocks-Mecking
!>   (addendum_kocks_mecking.md §9 item 4; Barton's parameters in FEPX units: b 2.48e-7 mm,
!>   h_n 0.125, k_1 4.5e5 /mm, k_20 14, gammadot_k0 1e10 /s, n_km 100; mu_0 Hill's of the
!>   stiffness) from rho_0 = 2e7 and 1e4 /mm^2; 6: Kocks-Mecking from 2e7 /mm^2 with
!>   irradiation softening (addendum_irradiation.md §10 item 4: h_d 0.675, eta 100, <111> loops
!>   of 3.61e12 /mm^3 and diameter 100 b; H in the state, damped state loop, D62).
!> Cases 7-18 (addendum_point_verification.md §5): start orientation Rodrigues vector
!> (0.1, 0.2, 0.3) (rotation of 2 atan|r| about r/|r|, crystal to sample), the case 1 Voce
!> parameters unless stated:
!>   7: fcc, isotropic, uniaxial; 8: fcc, isotropic, simple shear l = gammadot e1 (x) e2 with
!>   gammadot = sqrt(3) 1e-3 1/s (d_12 = w_12 = gammadot/2); 9: fcc, anisotropic (1.4), simple
!>   shear; 10: fcc, isotropic with saturation evolution (gammadot_s0 1e-2 1/s, m' 0.1);
!>   11: fcc, isotropic, cyclic (a 5e-3, c 2), uniaxial reversed at t = 20 s (accumulated slip
!>   reset there); 12: fcc, isotropic, precipitate cutting (a_p 100 MPa, f_p 0.01, r_p 1e-6 mm,
!>   b_p 2.5e-7 mm: tau_p = 20 MPa added to tau_0); 13: bcc {110}<111> + {112}<111> (iron),
!>   isotropic; 14: hcp basal + prismatic + pyramidal <c+a> (titanium, c/a 1.587), isotropic,
!>   per-family tau_0 150, 100, 250 MPa and m 0.02, 0.02, 0.05; 15: fcc, Kocks-Mecking with
!>   h_n = 0 (crss = tau_0 for all time); 16-18: fcc, Kocks-Mecking with k_1 = k_20 = 0 from
!>   rho_0 = 1e4, 1e2, 1 /mm^2 (crss = tau_0 + b mu_0 sqrt(h_n rho_0) for all time). Cases 7,
!>   10-18 are uniaxial.
!>
!> Output (file given as the first command argument), one line per case, dt and output time
!> (t = 0, 4, 40 s; case 11 also 20 s):
!>   case  dt  t  its  n  e(1:6)  rot(1:9)  crss(1:24)
!> its: the largest number of state iterations so far; n: number of slip systems; e: elastic
!> strain (Mandel 11 22 33 23 13 12, sample frame); rot: rot_c2s, column by column; crss: the
!> strengths [MPa] of the n systems, zero-padded to 24.
!>
!> Locals:
!>   crys          type(crystal_t)            slip systems of the case
!>   model         type(cp_model_t)           point model
!>   hard          type(hardening_params_t)   the law
!>   work          type(cp_work_t)            point-solve work arrays, with the per-point
!>                                            hardening values (cyclic: accumulated slip, d_eq)
!>   irr           type(hardening_irr_t)      irradiation parameters (case 6)
!>   dmg           real(rk)(6)                damage tensor H, Mandel, crystal frame [1/mm^2]
!>   stiff         real(rk)(6, 6)             crystal-frame stiffness [MPa]
!>   dts           real(rk)(5)                time increments [s]
!>   t_out         real(rk)(4)                output times [s]
!>   t_rev         real(rk)                   reversal time (0: none) [s]
!>   rot0          real(rk)(3, 3)             start orientation rot_c2s [-]
!>   d0, w0        real(rk)(6), real(rk)(3)   deformation rate (Mandel) and spin (w23 w13 w12)
!>                                            of the first loading segment, sample frame [1/s]
!>   d, w          the current segment's (sign -1 after the reversal) [1/s]
!>   m, gdot0, tau_0  real(rk)(:)             rate sensitivity [-], reference slip rate [1/s],
!>                                            initial strengths [MPa] per system
!>   accum         real(rk)(:)                accumulated slip since the last reversal [-]
!>   rot, sig, e   orientation, Cauchy stress [MPa], elastic strain carried between increments
!>   rot_n, sig_n, e_n, wp, sec, h, tan  point-solve output
!>   crss_prev, crss, rate  strengths [MPa], slip rates [1/s]
!>   rho, rho_0    Kocks-Mecking density (current, start) [1/mm^2]
!>   cyclic, is_km, irradiation  logical      the case's modifier, law and irradiation flags
!>   file_name     character(256)             output file
!>   icase, i_dt, k, num_incr, n, it, max_its, status, unit  indices, counts, status
program drive_hardening_point

  use kinds_mod, only: rk
  use crystal_mod, only: crystal_t, crystal_init, crystal_fcc, crystal_bcc, crystal_hcp, &
    & stiffness_cubic, stiffness_hexagonal, shear_modulus_hill
  use crystal_plasticity_model_mod, only: cp_model_t, cp_model_init, cp_point_solve, &
    & cp_converged
  use crystal_plasticity_model_mod, only: cp_work_t, cp_work_alloc
  use hardening_mod, only: hardening_params_t, hardening_point_set, &
    & hardening_isotropic, hardening_anisotropic, hardening_kocks_mecking, &
    & hardening_km_strengths, hardening_irr_t, hardening_irr_loop_average, &
    & hardening_irr_planes, hardening_irr_dmg_0, hardening_irr_strengths
  use orientation_mod, only: rot_c2s_from_descriptor, ori_des_rodrigues, ori_conv_passive
  use tensor_mod, only: equiv_strain

  implicit none

  type(cp_work_t) :: work
  integer, parameter :: num_cases = 18
  integer, parameter :: max_out_slip = 24
  type(crystal_t) :: crys
  type(cp_model_t) :: model
  type(hardening_params_t) :: hard
  type(hardening_irr_t) :: irr
  real(rk) :: dmg(6)
  real(rk), parameter :: dts(5) = [0.8_rk, 0.4_rk, 0.2_rk, 0.1_rk, 0.05_rk]
  real(rk), parameter :: t_out(4) = [0.0_rk, 4.0_rk, 20.0_rk, 40.0_rk]
  real(rk) :: t_rev, rot0(3, 3), d0(6), w0(3), d(6), w(3)
  real(rk) :: stiff(6, 6), rot(3, 3), sig(6), e(6), rot_n(3, 3), sig_n(6), e_n(6), wp(3)
  real(rk) :: sec(6, 6), h(6), tan(6, 6)
  real(rk) :: rho, rho_0
  real(rk), allocatable :: m(:), gdot0(:), tau_0(:), accum(:)
  real(rk), allocatable :: crss_prev(:), crss(:), rate(:)
  logical :: cyclic, is_km, irradiation
  character(len=256) :: file_name
  integer :: icase, i_dt, k, num_incr, n, it, max_its, status, unit

  call get_command_argument(1, file_name)
  open (newunit=unit, file=trim(file_name), status='replace', action='write')

  do icase = 1, num_cases
    ! The case: crystal, stiffness, law, start orientation and loading.
    call case_setup(icase)
    n = crys%num_slip
    if (allocated(crss_prev)) deallocate (crss_prev, crss, rate, accum)
    allocate (crss_prev(n), crss(n), rate(n), accum(n))
    call cp_model_init(model, stiff, crys%schmid_sym_c, crys%schmid_skew_c, m, gdot0, tau_0, &
      & tol_newton=1.0e-14_rk, max_iters_newton=100, max_iters_state=200, &
      & tol_state=1.0e-14_rk, hard=hard)
    ! Irradiation softening (case 6): the four <111> loop normals, H_0 by (I4).
    if (irradiation) then
      irr%h_d = 0.675_rk
      irr%eta = 100.0_rk
      irr%loop_density = 3.61e12_rk
      irr%loop_diameter = 100.0_rk*hard%km_b
      irr%loop_avg_c = hardening_irr_loop_average(reshape([1, 1, 1, -1, 1, 1, 1, -1, 1, &
        & 1, 1, -1]*1.0_rk, [3, 4]))
      irr%plane_c = hardening_irr_planes(crys%slip_normal_c)
      model%irradiation = .true.
      model%irr = irr
    end if

    call cp_work_alloc(work, model)
    do i_dt = 1, size(dts)
      ! Zero start state at the case's orientation; strengths from the law's start state.
      rot = rot0
      sig = 0.0_rk
      e = 0.0_rk
      accum = 0.0_rk
      crss_prev = model%tau_0
      rho = rho_0
      if (is_km) crss_prev = hardening_km_strengths(hard, model%tau_0, rho)
      dmg = 0.0_rk
      if (irradiation) then
        dmg = hardening_irr_dmg_0(irr, irr%loop_density, irr%loop_diameter)
        crss_prev = hardening_irr_strengths(hard, irr, model%tau_0, rho, dmg)
      end if
      max_its = 0
      num_incr = nint(t_out(4)/dts(i_dt))
      call write_line(0.0_rk)

      ! Increments; output at the output times.
      do k = 1, num_incr
        ! Loading segment: reversed (and the accumulated slip reset) from t_rev on.
        d = d0
        w = w0
        if (t_rev > 0.0_rk .and. (k - 1)*dts(i_dt) > t_rev - 1.0e-9_rk) then
          d = -d0
          w = -w0
          if (abs((k - 1)*dts(i_dt) - t_rev) < 1.0e-9_rk) accum = 0.0_rk
        end if
        ! Per-point values: the start-of-increment accumulated slip and d_eq (cyclic).
        call hardening_point_set(work%hardening%point, hard, slip_accum=accum, &
          & defrate_eq=equiv_strain(d))
        call cp_point_solve(model, crss_prev, rot, sig, e, d, w, dts(i_dt), sig_n, e_n, wp, &
          & rot_n, crss, sec, h, tan, rate, status, num_state_iters=it, km_rho=rho, dmg_c=dmg, &
          & work=work)
        if (status /= cp_converged) error stop 'drive_hardening_point: point solve failed'
        max_its = max(max_its, it)
        ! Accumulated slip with the end-of-increment rates (as cp_solve_points).
        if (cyclic) accum = accum + abs(rate)*dts(i_dt)
        crss_prev = crss
        rot = rot_n
        sig = sig_n
        e = e_n
        if (any(abs(k*dts(i_dt) - t_out(2:)) < 1.0e-9_rk)) then
          ! The 20 s output only for the case with a reversal.
          if (abs(k*dts(i_dt) - t_out(3)) > 1.0e-9_rk .or. t_rev > 0.0_rk) &
            & call write_line(k*dts(i_dt))
        end if
      end do
    end do
  end do
  close (unit)

contains

  !> Sets the case's crystal, stiffness, slip parameters, law, start orientation and loading
  !> (program header) in the program's variables crys, stiff, m, gdot0, tau_0, hard, rho_0,
  !> cyclic, is_km, irradiation, rot0, d0, w0, t_rev.
  !>
  !> Arguments:
  !>   icase  (in)  integer  case number 1..num_cases
  !> Locals:
  !>   ns       integer           number of slip systems
  !>   b        integer           system
  !>   latent   real(rk)          off-diagonal of the interaction matrix [-]
  !>   gdot     real(rk)          simple-shear rate [1/s]
  !>   tau_p    real(rk)          precipitate strength (cutting) [MPa]
  !>   rod      real(rk)(3)       start Rodrigues vector of cases 7-18 [-]
  !>   uniaxial real(rk)(6)       the uniaxial deformation rate, Mandel [1/s]
  !>   hcp_tau_0, hcp_m  real(rk)(3)  hcp initial strength [MPa] and m [-] per family (case 14)
  !>   km_rho_limit   real(rk)(3)  rho_0 of cases 16-18 [1/mm^2]
  subroutine case_setup(icase)

    integer, intent(in) :: icase

    integer :: ns, b
    real(rk) :: latent, gdot, tau_p
    real(rk), parameter :: rod(3) = [0.1_rk, 0.2_rk, 0.3_rk]
    real(rk), parameter :: uniaxial(6) = 1.0e-3_rk*[-0.5_rk, -0.5_rk, 1.0_rk, 0.0_rk, &
      & 0.0_rk, 0.0_rk]
    real(rk), parameter :: hcp_tau_0(3) = [150.0_rk, 100.0_rk, 250.0_rk]
    real(rk), parameter :: hcp_m(3) = [0.02_rk, 0.02_rk, 0.05_rk]
    real(rk), parameter :: km_rho_limit(3) = [1.0e4_rk, 1.0e2_rk, 1.0_rk]

    ! Crystal and stiffness: fcc copper unless stated.
    select case (icase)
    case (13)
      call crystal_init(crys, crystal_bcc, families=[1, 2])
      stiff = stiffness_cubic(231400.0_rk, 134700.0_rk, 116400.0_rk)
    case (14)
      call crystal_init(crys, crystal_hcp, c_over_a=1.587_rk)
      stiff = stiffness_hexagonal(162400.0_rk, 92000.0_rk, 69000.0_rk, 180700.0_rk, &
        & 46700.0_rk)
    case default
      call crystal_init(crys, crystal_fcc)
      stiff = stiffness_cubic(168400.0_rk, 121400.0_rk, 75400.0_rk)
    end select
    ns = crys%num_slip

    ! Slip parameters per system: m 0.02, gammadot_0 1 1/s, tau_0 100 MPa unless stated.
    if (allocated(m)) deallocate (m, gdot0, tau_0)
    allocate (m(ns), gdot0(ns), tau_0(ns))
    m = 0.02_rk
    gdot0 = 1.0_rk
    tau_0 = 100.0_rk
    if (icase == 14) then
      ! hcp, per family (basal, prismatic, pyramidal <c+a>).
      tau_0 = hcp_tau_0(crys%slip_family)
      m = hcp_m(crys%slip_family)
    end if
    if (icase == 12) then
      ! Precipitate cutting (H-P1), added to tau_0 (as fepx.f90 does for the phase).
      tau_p = 100.0_rk*sqrt(0.01_rk*1.0e-6_rk/2.5e-7_rk)
      tau_0 = tau_0 + tau_p
    end if

    ! The law: isotropic unless stated; Voce parameters of case 1.
    latent = 0.0_rk
    if (icase == 2 .or. icase == 9) latent = 1.4_rk
    if (icase == 3) latent = 5.0_rk
    hard = hardening_params_t()
    hard%h_0 = spread(200.0_rk, 1, ns)
    hard%tau_s = spread(330.0_rk, 1, ns)
    hard%n = spread(1.5_rk, 1, ns)
    hard%tol_hardening = 1.0e-13_rk
    hard%max_iters = 50
    hard%interaction = spread(spread(latent, 1, ns), 1, ns)
    do b = 1, ns
      hard%interaction(b, b) = 1.0_rk
    end do
    hard%law = hardening_isotropic
    if (latent > 0.0_rk) hard%law = hardening_anisotropic
    ! Modifiers.
    if (icase == 10) then
      hard%saturation_evolution = .true.
      hard%gammadot_s0 = 1.0e-2_rk
      hard%m_prime = 0.1_rk
    end if
    cyclic = icase == 11
    if (cyclic) then
      hard%cyclic = .true.
      hard%cyclic_a = 5.0e-3_rk
      hard%cyclic_c = 2.0_rk
    end if
    ! Kocks-Mecking (cases 4-6, 15-18).
    is_km = (icase >= 4 .and. icase <= 6) .or. icase >= 15
    irradiation = icase == 6
    rho_0 = 0.0_rk
    if (is_km) then
      hard%law = hardening_kocks_mecking
      hard%km_b = 2.48e-7_rk
      hard%km_mu_0 = shear_modulus_hill(stiff)
      hard%km_h_n = 0.125_rk
      hard%km_k1 = 4.5e5_rk
      hard%km_k20 = 14.0_rk
      hard%km_gammadot_k0 = 1.0e10_rk
      hard%km_n = 100.0_rk
      rho_0 = 2.0e7_rk
      if (icase == 5) rho_0 = 1.0e4_rk
      if (icase == 15) hard%km_h_n = 0.0_rk
      if (icase >= 16) then
        hard%km_k1 = 0.0_rk
        hard%km_k20 = 0.0_rk
        rho_0 = km_rho_limit(icase - 15)
      end if
      hard%km_rho_0 = rho_0
    end if

    ! Start orientation: cube (cases 1-6) or the Rodrigues vector rod (the passive descriptor
    ! gives rot_c2s = the rotation by 2 atan|rod| about rod/|rod|, conventions.md §2; the
    ! check script verifies it against its own (P5) at t = 0).
    rot0 = 0.0_rk
    do b = 1, 3
      rot0(b, b) = 1.0_rk
    end do
    if (icase >= 7) rot0 = rot_c2s_from_descriptor(ori_des_rodrigues, ori_conv_passive, rod)

    ! Loading: uniaxial, or simple shear (cases 8, 9); reversal at 20 s (case 11).
    d0 = uniaxial
    w0 = 0.0_rk
    if (icase == 8 .or. icase == 9) then
      gdot = sqrt(3.0_rk)*1.0e-3_rk
      d0 = 0.0_rk
      d0(6) = sqrt(2.0_rk)*0.5_rk*gdot
      w0 = [0.0_rk, 0.0_rk, 0.5_rk*gdot]
    end if
    t_rev = 0.0_rk
    if (icase == 11) t_rev = 20.0_rk

  end subroutine case_setup

  !> Writes one output line of the current case and increment (program header).
  !>
  !> Arguments:
  !>   t  (in)  real(rk)  time [s]
  !> Locals:
  !>   crss_out  real(rk)(max_out_slip)  strengths, zero-padded [MPa]
  subroutine write_line(t)

    real(rk), intent(in) :: t

    real(rk) :: crss_out(max_out_slip)

    ! Strengths of the start state at t = 0, else of the last increment.
    crss_out = 0.0_rk
    if (t > 0.0_rk) then
      crss_out(1:n) = crss
    else
      crss_out(1:n) = crss_prev
    end if
    write (unit, '(i3, 2es12.4, 2i5, 39es25.16)') icase, dts(i_dt), t, max_its, n, e, &
      & reshape(rot, [9]), crss_out

  end subroutine write_line

end program drive_hardening_point
