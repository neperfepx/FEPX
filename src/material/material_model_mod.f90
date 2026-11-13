! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> The material-solver contract (ROADMAP D38): the abstract material model type that every
!> constitutive model extends, one instance per phase.
!>
!> Division of work:
!>   - the FE solver does all kinematics: at every quadrature point it computes the deformation
!>     rate {d}_s (Mandel 6-vector) and the spin w_s (skew 3-vector, w23 w13 w12) from the nodal
!>     velocities (element_mod), and passes them with the time increment; the material never
!>     sees geometry or velocities;
!>   - the material model, per point of its phase (pure per-point kernels; polymorphic dispatch
!>     happens once per phase, never per point, conventions §6), solves for the current-iterate
!>     state from the start-of-increment state and returns the matrices of
!>     recoupled_formulation.md (R10)-(R12): the secant [s]_s, the vector {h}_s, and, when the
!>     solver asks for Newton-Raphson, the tangent [s_tan]_s; the solver then forms the element
!>     matrices (R14) and the stress is {sigma} = [s]({d} - {h});
!>   - the state lives in material_state_t (D39); accepting an increment (state_prev = state)
!>     is done by the solver.
!>
!> Source : FEPX (ROADMAP step 4.1, D38). v2 computes {d} and w inside its material routines
!>          (solve_evp_mod2.f90::solve_evp_pre) and returns its 5 x 5 matrices from
!>          solveit_evp_mod4.f90::material_matrix_evp_aniso.
!>
!> Public types:
!>   material_model_t   abstract model (components and deferred procedures at the type)
module material_model_mod

  use kinds_mod, only: rk, ik
  use material_state_mod, only: material_state_t

  implicit none

  private

  public :: material_model_t

  !> Abstract material model (one instance per phase).
  !>
  !> Components:
  !>   num_slip  integer(ik)  slip systems per point of this model (0: no slip)
  !>
  !> Deferred procedures:
  !>   init_state    set the model's initial state at the points of the given elements
  !>   solve_points  point solve at the points of the given elements: current-iterate state,
  !>                 [s]_s, {h}_s, and optionally [s_tan]_s
  type, abstract :: material_model_t
    integer(ik) :: num_slip = 0_ik
  contains
    procedure(init_state_iface), deferred :: init_state
    procedure(solve_points_iface), deferred :: solve_points
  end type material_model_t

  abstract interface

    !> Set the model's initial state (e.g. slip-system strengths) at every quadrature point of
    !> the given elements. Orientations are set by the caller from the mesh.
    !>
    !> Arguments:
    !>   self   (in)    class(material_model_t)  model
    !>   elts   (in)    integer(ik)(:)           local elements of this model's phase
    !>   state  (inout) type(material_state_t)   state (allocated by the caller)
    subroutine init_state_iface(self, elts, state)
      import :: material_model_t, material_state_t, ik
      class(material_model_t), intent(in) :: self
      integer(ik), intent(in) :: elts(:)
      type(material_state_t), intent(inout) :: state
    end subroutine init_state_iface

    !> Point solve at every quadrature point of the given elements.
    !>
    !> Arguments:
    !>   self           (in)    class(material_model_t)  model
    !>   elts           (in)    integer(ik)(:)           local elements of this model's phase
    !>   dt             (in)    real(rk)                 time increment [s]
    !>   defrate_s      (in)    real(rk)(6, :, :)        {d}, Mandel, sample frame [1/s]
    !>   spin_s         (in)    real(rk)(3, :, :)        w (w23 w13 w12), sample frame [1/s]
    !>   state_prev     (in)    type(material_state_t)   state at the start of the increment
    !>   state          (inout) type(material_state_t)   current-iterate state (out at the points
    !>                                                   of elts)
    !>   mat_secant_s   (inout) real(rk)(6, 6, :, :)     [s], sample frame [MPa s] (out at elts)
    !>   vec_h_s        (inout) real(rk)(6, :, :)        {h}, sample frame [1/s] (out at elts)
    !>   num_failed     (out)   integer(ik)              points whose solve did not converge
    !>   mat_tangent_s  (inout, optional) real(rk)(6, 6, :, :)  [s_tan], sample frame [MPa s]
    !>                                                   (out at elts; Newton-Raphson only)
    subroutine solve_points_iface(self, elts, dt, defrate_s, spin_s, state_prev, state, &
      & mat_secant_s, vec_h_s, num_failed, mat_tangent_s)
      import :: material_model_t, material_state_t, rk, ik
      class(material_model_t), intent(in) :: self
      integer(ik), intent(in) :: elts(:)
      real(rk), intent(in) :: dt
      real(rk), intent(in) :: defrate_s(:, :, :)
      real(rk), intent(in) :: spin_s(:, :, :)
      type(material_state_t), intent(in) :: state_prev
      type(material_state_t), intent(inout) :: state
      real(rk), intent(inout) :: mat_secant_s(:, :, :, :)
      real(rk), intent(inout) :: vec_h_s(:, :, :)
      integer(ik), intent(out) :: num_failed
      real(rk), intent(inout), optional :: mat_tangent_s(:, :, :, :)
    end subroutine solve_points_iface

  end interface

end module material_model_mod
