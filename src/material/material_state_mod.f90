! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Material state at the quadrature points of the rank's elements: one non-polymorphic type,
!> shared by every material model, holding named component-first arrays (ROADMAP D39).
!>
!> Fields common to all models are always allocated; the slip-system fields (crss, sliprate,
!> rss, slip, slip_accum) have num_slip rows (zero for the elastic-only model, so they are allocated with
!> zero size and no code tests whether they exist). The output quantities of step 6.5 (slip
!> rates, resolved shear stresses, accumulated slip and plastic strain, plastic deformation rate
!> and spin, lattice rotation rate, work) are state fields too: the models set them in their
!> point solves from the start-of-increment state, like every other field, and
!> material_state_update_derived sets the model-independent ones (the point kinematics, total
!> strain, sample-frame plastic strain, work and plastic work, rotation rate; step 6.5b).
!> Variables that only some models evolve in addition (e.g. dislocation density, irradiation
!> damage, twin fractions) are added here as named fields with those models (Stages 7 and 11):
!> the network dislocation density km_rho and the recovery coefficient km_k2 and the density rate km_rhorate of the Kocks-Mecking law
!> (step 7.2, D60; one value per point, zero for the other models). An increment is accepted by
!> assignment: state_prev = state.
!>
!> Source : FEPX (ROADMAP step 4.1, D38-D39). v2 keeps the corresponding values in
!>          results_type (sig_vec, e_elas, ori, wp_hat, crss), with 5-vector deviatoric parts and
!>          separate volumetric scalars; FEPX stores full Mandel 6-vectors (D2).
!>
!> Public types:
!>   material_state_t  state of every quadrature point (components documented at the type)
!>
!> Public routines:
!>   material_state_alloc           allocate (and zero) a state for given sizes
!>   material_state_update_derived  kinematics, strains, work and rotation rate of every point
module material_state_mod

  use kinds_mod, only: rk, ik
  use rotation_mod, only: rot_c2s_rotation_rate

  implicit none

  private

  public :: material_state_t, material_state_alloc, material_state_update_derived

  !> State of every quadrature point of the rank's elements.
  !>
  !> Components (sizes; frames; units):
  !>   num_qpts     integer(ik)    quadrature points per element
  !>   num_elts     integer(ik)    elements
  !>   num_slip     integer(ik)    slip systems per point (0: elastic-only)
  !>   stress_s     real(rk)(6, num_qpts, num_elts)     Cauchy stress, Mandel, sample [MPa]
  !>   strain_el_s  real(rk)(6, num_qpts, num_elts)     elastic strain e^e, Mandel, sample [-]
  !>   rot_c2s      real(rk)(3, 3, num_qpts, num_elts)  crystal-to-sample rotation [-]
  !>   spin_lat_s   real(rk)(3, num_qpts, num_elts)     lattice spin w^p (w23, w13, w12),
  !>                                                     sample [1/s]
  !>   crss         real(rk)(num_slip, num_qpts, num_elts)  slip-system strengths [MPa]
  !>   sliprate     real(rk)(num_slip, num_qpts, num_elts)  slip rates gammadot [1/s]
  !>   rss          real(rk)(num_slip, num_qpts, num_elts)  resolved shear stresses
  !>                                                     tau^a = p^a : tau (Kirchhoff) [MPa]
  !>   slip         real(rk)(num_slip, num_qpts, num_elts)  accumulated slip, sum of
  !>                                                     gammadot dt (signed) [-]
  !>   slip_accum   real(rk)(num_slip, num_qpts, num_elts)  accumulated slip magnitude since
  !>                                                     the last reversal of the control
  !>                                                     strain rate, sum of |gammadot| dt
  !>                                                     (cyclic hardening; v2 acmslip; reset
  !>                                                     to zero by the driver) [-]
  !>   strain_pl_c  real(rk)(6, num_qpts, num_elts)     accumulated plastic strain, sum of
  !>                                                     d^p dt in the crystal frame, Mandel [-]
  !>   defrate_pl_s real(rk)(6, num_qpts, num_elts)     plastic deformation rate d^p, Mandel,
  !>                                                     sample [1/s]
  !>   spin_pl_s    real(rk)(3, num_qpts, num_elts)     plastic (slip) spin sum gammadot q
  !>                                                     (w23, w13, w12), sample [1/s]
  !>   rotrate_s    real(rk)(3, num_qpts, num_elts)     lattice rotation rate of the increment,
  !>                                                     log(R R_prev^T)/dt (w23, w13, w12),
  !>                                                     sample [1/s]
  !>   workrate     real(rk)(num_qpts, num_elts)        work rate sigma : d [MPa/s]
  !>   work         real(rk)(num_qpts, num_elts)        work, trapezoidal integral of the work
  !>                                                     rate [MPa]
  !>   defrate_s    real(rk)(6, num_qpts, num_elts)     deformation rate d = sym(grad v) of the
  !>                                                     current velocities, Mandel, sample [1/s]
  !>   spin_s       real(rk)(3, num_qpts, num_elts)     continuum spin w = skw(grad v)
  !>                                                     (w23, w13, w12), sample [1/s]
  !>   strain_s     real(rk)(6, num_qpts, num_elts)     total strain, sum of d dt in the sample
  !>                                                     frame, Mandel [-]
  !>   strain_pl_s  real(rk)(6, num_qpts, num_elts)     plastic strain accumulated in the sample
  !>                                                     frame, sum of d^p dt, Mandel [-]
  !>   spin_hat_update_s   real(rk)(3, num_qpts, num_elts)  the spin w^p of the increment's
  !>                                                     orientation update (w23, w13, w12),
  !>                                                     sample [1/s] (output rotrate_spin)
  !>   spin_slip_update_s  real(rk)(3, num_qpts, num_elts)  the slip spin sum gammadot q of
  !>                                                     that update, sample [1/s] (output
  !>                                                     rotrate_slip; zero for elastic-only)
  !>   start_from_vertex  logical                    the point solves of this increment start
  !>                                                     from v2's viscoplastic vertex estimate
  !>                                                     (initialization
  !>                                                     isotropic_viscoplastic,
  !>                                                     increment 1; set by the driver) [-]
  !>   workrate_pl  real(rk)(num_qpts, num_elts)        plastic work rate sigma : d^p [MPa/s]
  !>   work_pl      real(rk)(num_qpts, num_elts)        plastic work, trapezoidal integral of
  !>                                                     the plastic work rate [MPa]
  !>   km_rho       real(rk)(num_qpts, num_elts)        network dislocation density of the
  !>                                                     Kocks-Mecking law (the law's stored
  !>                                                     state; zero for other models) [1/mm^2]
  !>   km_k2        real(rk)(num_qpts, num_elts)        Kocks-Mecking recovery coefficient at
  !>                                                     the slip rates of the final stress
  !>                                                     (zero where Gammadot = 0 and for other
  !>                                                     models) [-]
  !>   dmg_c        real(rk)(6, num_qpts, num_elts)     irradiation damage tensor H (Mandel,
  !>                                                     crystal frame; step 7.3, D61; zero
  !>                                                     without irradiation) [1/mm^2]
  !>   dmgrate_c    real(rk)(6, num_qpts, num_elts)     its rate H dot at the slip rates of
  !>                                                     the final stress [1/(mm^2 s)]
  !>   km_rhorate   real(rk)(num_qpts, num_elts)        rate of km_rho, (K2) at the slip rates
  !>                                                     of the final stress and the updated
  !>                                                     density (zero where Gammadot = 0 and
  !>                                                     for other models) [1/(mm^2 s)]
  type :: material_state_t
    integer(ik) :: num_qpts = 0_ik
    integer(ik) :: num_elts = 0_ik
    integer(ik) :: num_slip = 0_ik
    logical :: start_from_vertex = .false.
    real(rk), allocatable :: stress_s(:, :, :)
    real(rk), allocatable :: strain_el_s(:, :, :)
    real(rk), allocatable :: rot_c2s(:, :, :, :)
    real(rk), allocatable :: spin_lat_s(:, :, :)
    real(rk), allocatable :: crss(:, :, :)
    real(rk), allocatable :: sliprate(:, :, :)
    real(rk), allocatable :: rss(:, :, :)
    real(rk), allocatable :: slip(:, :, :)
    real(rk), allocatable :: slip_accum(:, :, :)
    real(rk), allocatable :: strain_pl_c(:, :, :)
    real(rk), allocatable :: defrate_pl_s(:, :, :)
    real(rk), allocatable :: spin_pl_s(:, :, :)
    real(rk), allocatable :: rotrate_s(:, :, :)
    real(rk), allocatable :: workrate(:, :)
    real(rk), allocatable :: work(:, :)
    real(rk), allocatable :: defrate_s(:, :, :)
    real(rk), allocatable :: spin_s(:, :, :)
    real(rk), allocatable :: strain_s(:, :, :)
    real(rk), allocatable :: strain_pl_s(:, :, :)
    real(rk), allocatable :: spin_hat_update_s(:, :, :)
    real(rk), allocatable :: spin_slip_update_s(:, :, :)
    real(rk), allocatable :: workrate_pl(:, :)
    real(rk), allocatable :: work_pl(:, :)
    real(rk), allocatable :: km_rho(:, :)
    real(rk), allocatable :: km_k2(:, :)
    real(rk), allocatable :: km_rhorate(:, :)
    real(rk), allocatable :: dmg_c(:, :, :)
    real(rk), allocatable :: dmgrate_c(:, :, :)
  end type material_state_t

contains

  !> Allocate a state and set every field to zero (rot_c2s to the identity; the orientations
  !> are set from the mesh by the caller).
  !>
  !> Arguments:
  !>   state     (out) type(material_state_t)  state
  !>   num_qpts  (in)  integer(ik)             quadrature points per element
  !>   num_elts  (in)  integer(ik)             elements
  !>   num_slip  (in)  integer(ik)             slip systems per point (0: elastic-only)
  !> Locals:
  !>   i         integer                       diagonal index
  subroutine material_state_alloc(state, num_qpts, num_elts, num_slip)

    type(material_state_t), intent(out) :: state
    integer(ik), intent(in) :: num_qpts
    integer(ik), intent(in) :: num_elts
    integer(ik), intent(in) :: num_slip

    integer :: i

    ! Sizes.
    state%num_qpts = num_qpts
    state%num_elts = num_elts
    state%num_slip = num_slip

    ! Common fields, zero; orientations the identity.
    allocate (state%stress_s(6, num_qpts, num_elts), source=0.0_rk)
    allocate (state%strain_el_s(6, num_qpts, num_elts), source=0.0_rk)
    allocate (state%spin_lat_s(3, num_qpts, num_elts), source=0.0_rk)
    allocate (state%rot_c2s(3, 3, num_qpts, num_elts), source=0.0_rk)
    do i = 1, 3
      state%rot_c2s(i, i, :, :) = 1.0_rk
    end do

    ! Slip-system fields (zero rows for the elastic-only model).
    allocate (state%crss(num_slip, num_qpts, num_elts), source=0.0_rk)
    allocate (state%sliprate(num_slip, num_qpts, num_elts), source=0.0_rk)
    allocate (state%rss(num_slip, num_qpts, num_elts), source=0.0_rk)
    allocate (state%slip(num_slip, num_qpts, num_elts), source=0.0_rk)
    allocate (state%slip_accum(num_slip, num_qpts, num_elts), source=0.0_rk)

    ! Plastic and kinematic output quantities (zero at the initial state).
    allocate (state%strain_pl_c(6, num_qpts, num_elts), source=0.0_rk)
    allocate (state%defrate_pl_s(6, num_qpts, num_elts), source=0.0_rk)
    allocate (state%spin_pl_s(3, num_qpts, num_elts), source=0.0_rk)
    allocate (state%rotrate_s(3, num_qpts, num_elts), source=0.0_rk)
    allocate (state%workrate(num_qpts, num_elts), source=0.0_rk)
    allocate (state%work(num_qpts, num_elts), source=0.0_rk)

    ! Point kinematics, total and sample-frame plastic strain, plastic work (step 6.5b).
    allocate (state%defrate_s(6, num_qpts, num_elts), source=0.0_rk)
    allocate (state%spin_s(3, num_qpts, num_elts), source=0.0_rk)
    allocate (state%strain_s(6, num_qpts, num_elts), source=0.0_rk)
    allocate (state%strain_pl_s(6, num_qpts, num_elts), source=0.0_rk)
    allocate (state%spin_hat_update_s(3, num_qpts, num_elts), source=0.0_rk)
    allocate (state%spin_slip_update_s(3, num_qpts, num_elts), source=0.0_rk)
    allocate (state%workrate_pl(num_qpts, num_elts), source=0.0_rk)
    allocate (state%work_pl(num_qpts, num_elts), source=0.0_rk)

    ! Kocks-Mecking fields (step 7.2; set by that law's model, zero otherwise).
    allocate (state%km_rho(num_qpts, num_elts), source=0.0_rk)
    allocate (state%km_k2(num_qpts, num_elts), source=0.0_rk)
    allocate (state%km_rhorate(num_qpts, num_elts), source=0.0_rk)

    ! Irradiation damage (step 7.3; set by that model, zero otherwise).
    allocate (state%dmg_c(6, num_qpts, num_elts), source=0.0_rk)
    allocate (state%dmgrate_c(6, num_qpts, num_elts), source=0.0_rk)

  end subroutine material_state_alloc

  !> Model-independent derived quantities of every point at the current iterate: the point
  !> kinematics (d, w), the total strain and the sample-frame plastic strain, the work rate
  !> sigma : d and the plastic work rate sigma : d^p with their time integrals (trapezoidal rule
  !> over the increment), and the lattice rotation rate of the increment. Called after the point
  !> solves (increment_mod::element_system), so the values at the converged velocities are those
  !> of the last call. The models set d^p (defrate_pl_s) in their point solves; it is zero for
  !> the elastic-only model.
  !>
  !> Theory : work rate = sigma : d (Cauchy stress, total deformation rate; the Mandel dot product
  !>          is the full contraction); plastic work rate = sigma : d^p (= sigma' : d^p, d^p is
  !>          deviatoric); work = work_prev + dt/2 (workrate_prev + workrate), the same for the
  !>          plastic work; strain = strain_prev + d dt and strain_pl_s = strain_pl_s_prev +
  !>          d^p dt, both summed in the sample frame (backward Euler with the end-of-increment
  !>          rates); rotation rate = log(R R_prev^T)/dt, the sample-frame rotation of the
  !>          increment divided by dt (the lattice spin when the orientation update is
  !>          exp(W dt)).
  !> Source : port of v2-fixed finalize_res_mod1.f90: finalize_res (defrate = the deviatoric d
  !>          plus the volumetric part of the velocity gradient; strain = strain + defrate dtime
  !>          in the sample frame), finalize_res_totalwork (sigma : d, trapezoidal rule with the
  !>          previous increment's rate), finalize_res_plasticwork (sigma' : d^p in the sample
  !>          frame, the same rule), and write_res_mod1.f90 rotrate (rotation vector of
  !>          d_rstar = R_prev^T R, divided by dtime, in the sample frame; with CHANGES.md C41:
  !>          the increment's own dtime). The rotation vector of R_prev^T R in the sample frame
  !>          is that of R R_prev^T. strain_pl_s is not a v2 quantity: the sample-frame variant
  !>          of v2's strain_pl (crystal-frame sum), ROADMAP D53 (2).
  !>
  !> Arguments:
  !>   state_prev  (in)    type(material_state_t)  state at the start of the increment
  !>   defrate_s   (in)    real(rk)(6, :, :)       {d}, Mandel, sample frame, every point [1/s]
  !>   spin_s      (in)    real(rk)(3, :, :)       w (w23, w13, w12), sample frame, every point
  !>                                               [1/s]
  !>   dt          (in)    real(rk)                time increment [s]
  !>   state       (inout) type(material_state_t)  current-iterate state (stress_s, rot_c2s and
  !>                                               defrate_pl_s in; defrate_s, spin_s,
  !>                                               strain_s, strain_pl_s, workrate, work,
  !>                                               workrate_pl, work_pl, rotrate_s out)
  !> Locals:
  !>   i_elt, i_qpt  integer(ik)    element, point
  subroutine material_state_update_derived(state_prev, defrate_s, spin_s, dt, state)

    type(material_state_t), intent(in) :: state_prev
    real(rk), intent(in) :: defrate_s(:, :, :)
    real(rk), intent(in) :: spin_s(:, :, :)
    real(rk), intent(in) :: dt
    type(material_state_t), intent(inout) :: state

    integer(ik) :: i_elt, i_qpt

    do i_elt = 1_ik, state%num_elts
      do i_qpt = 1_ik, state%num_qpts
        ! Point kinematics of the current velocities.
        state%defrate_s(:, i_qpt, i_elt) = defrate_s(:, i_qpt, i_elt)
        state%spin_s(:, i_qpt, i_elt) = spin_s(:, i_qpt, i_elt)

        ! Total strain and sample-frame plastic strain (sums of the rates times dt).
        state%strain_s(:, i_qpt, i_elt) = state_prev%strain_s(:, i_qpt, i_elt) + &
          & dt*defrate_s(:, i_qpt, i_elt)
        state%strain_pl_s(:, i_qpt, i_elt) = state_prev%strain_pl_s(:, i_qpt, i_elt) + &
          & dt*state%defrate_pl_s(:, i_qpt, i_elt)

        ! Work rate and work (trapezoidal rule over the increment).
        state%workrate(i_qpt, i_elt) = dot_product(state%stress_s(:, i_qpt, i_elt), &
          & defrate_s(:, i_qpt, i_elt))
        state%work(i_qpt, i_elt) = state_prev%work(i_qpt, i_elt) + 0.5_rk*dt* &
          & (state_prev%workrate(i_qpt, i_elt) + state%workrate(i_qpt, i_elt))

        ! Plastic work rate and plastic work (the same rule).
        state%workrate_pl(i_qpt, i_elt) = dot_product(state%stress_s(:, i_qpt, i_elt), &
          & state%defrate_pl_s(:, i_qpt, i_elt))
        state%work_pl(i_qpt, i_elt) = state_prev%work_pl(i_qpt, i_elt) + 0.5_rk*dt* &
          & (state_prev%workrate_pl(i_qpt, i_elt) + state%workrate_pl(i_qpt, i_elt))

        ! Lattice rotation rate: the rotation of the increment in the sample frame over dt.
        state%rotrate_s(:, i_qpt, i_elt) = rot_c2s_rotation_rate( &
          & state_prev%rot_c2s(:, :, i_qpt, i_elt), state%rot_c2s(:, :, i_qpt, i_elt), dt)
      end do
    end do

  end subroutine material_state_update_derived

end module material_state_mod
