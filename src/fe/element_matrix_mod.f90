! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Element-level finite element operations of the recoupled formulation: the kinematics at the
!> quadrature points (the solver's side of the material contract, D38) and the element matrix
!> and force vector (R14).
!>
!> Theory : recoupled_formulation.md §6: {d} = [B]{V} (R13), w = skw(grad v);
!>            [k_e] = sum_q w_q det J_q [B_q]^T [s_q] [B_q],
!>            {f_h,e} = sum_q w_q det J_q [B_q]^T [s_q] {h_q}                        (R14)
!>          (the full stiffness, volumetric part included, at every quadrature point).
!> Source : port of v2 solveit_evp_mod2.f90::elt_stif_evp (element matrices, sum over the
!>          quadrature points of B^T c B det J w, full symmetric storage as v2's
!>          stiffness_mod.f90::add_to_stiffness) and solve_evp_mod2.f90::solve_evp_pre (velocity
!>          gradient, its symmetric and skew parts), in Mandel form (v2 assembles a deviatoric
!>          5 x 3 form plus a separate volumetric term).
!>
!> Work arrays (ROADMAP D100): the routines take an element_work_t with their arrays at the
!> element's size (shape derivatives, B, [s][B]), allocated once by
!> element_work_alloc outside the element loops, so that nothing is allocated per element or
!> point.
!>
!> Public types:
!>   element_work_t          work arrays of the routines below (components at the type)
!>
!> Public routines:
!>   element_work_alloc      allocate the work arrays for an element type
!>   element_kinematics_qpt  {d}_s and w_s at one quadrature point from the element velocities
!>   element_stiffness       element matrix [k_e] (R14)
!>   element_force_h         element vector {f_h,e} (R14)
module element_matrix_mod

  use kinds_mod, only: rk, ik
  use tensor_mod, only: mandel_from_matrix, skew_from_matrix
  use element_mod, only: element_t, element_geometry_qpt, mandel_b_matrix

  implicit none

  private

  public :: element_work_t, element_work_alloc, element_kinematics_qpt, element_stiffness, &
    & element_force_h

  !> Work arrays of the element routines, at the size of one element type (element_work_alloc;
  !> ROADMAP D100); their values are not kept from one call to the next.
  !>
  !> Components:
  !>   dshape_s  real(rk)(3, num_nodes)        d N_a / d x_i at a point [1/mm]
  !>   b         real(rk)(6, num_dofs)         Mandel B matrix at a point [1/mm]
  !>   sb        real(rk)(6, num_dofs)         [s][B] [MPa s/mm]
  type :: element_work_t
    real(rk), allocatable :: dshape_s(:, :)
    real(rk), allocatable :: b(:, :)
    real(rk), allocatable :: sb(:, :)
  end type element_work_t

contains

  !> Allocate the work arrays of the element routines for an element type, at its size (once,
  !> outside the element loops; ROADMAP D100); arrays already of that size are kept.
  !>
  !> Arguments:
  !>   work  (inout) type(element_work_t)  work arrays
  !>   elt   (in)    type(element_t)       element descriptor
  !> Locals:
  !>   nn, nd  integer                     nodes and dofs of the element
  pure subroutine element_work_alloc(work, elt)

    type(element_work_t), intent(inout) :: work
    type(element_t), intent(in) :: elt

    integer :: nn, nd

    ! Nothing to do when the arrays have the element's size.
    nn = int(elt%num_nodes)
    nd = int(elt%num_dofs)
    if (allocated(work%dshape_s)) then
      if (size(work%dshape_s, 2) == nn .and. size(work%b, 2) == nd) return
      deallocate (work%dshape_s, work%b, work%sb)
    end if

    ! The arrays at the element's size.
    allocate (work%dshape_s(3, nn), work%b(6, nd), work%sb(6, nd))

  end subroutine element_work_alloc

  !> Deformation rate and spin at one quadrature point of one element.
  !>
  !> Theory : grad v(i, j) = sum_a V_a(i) d N_a / d x_j; {d} = {sym grad v} = [B]{V} (R13);
  !>          w = skw(grad v) as (w23, w13, w12).
  !> Source : port of v2 solve_evp_mod2.f90::solve_evp_pre (nodevel_velgrad, velgrad_sympart,
  !>          velgrad_skewpart).
  !>
  !> Arguments:
  !>   elt        (in)  type(element_t)  element descriptor
  !>   coo        (in)  real(rk)(:, :)   (3, num_nodes) node coordinates, sample frame [mm]
  !>   vel        (in)  real(rk)(:, :)   (3, num_nodes) node velocities, sample frame [mm/s]
  !>   i_qpt      (in)  integer(ik)      quadrature point
  !>   defrate_s  (out) real(rk)(6)      {d}, Mandel, sample frame [1/s]
  !>   spin_s     (out) real(rk)(3)      w (w23 w13 w12), sample frame [1/s]
  !>   work       (inout) type(element_work_t)  work arrays (element_work_alloc): dshape_s,
  !>                                     d N_a / d x_i [1/mm]
  !> Locals:
  !>   det_jac    real(rk)               det J [mm^3]
  !>   velgrad    real(rk)(3, 3)         grad v [1/s]
  pure subroutine element_kinematics_qpt(elt, coo, vel, i_qpt, defrate_s, spin_s, work)

    type(element_t), intent(in) :: elt
    real(rk), intent(in) :: coo(:, :)
    real(rk), intent(in) :: vel(:, :)
    integer(ik), intent(in) :: i_qpt
    real(rk), intent(out) :: defrate_s(6)
    real(rk), intent(out) :: spin_s(3)
    type(element_work_t), intent(inout) :: work

    real(rk) :: det_jac
    real(rk) :: velgrad(3, 3)

    ! Velocity gradient, then its symmetric (Mandel) and skew parts.
    call element_geometry_qpt(elt, coo, i_qpt, work%dshape_s, det_jac)
    velgrad = matmul(vel, transpose(work%dshape_s))
    defrate_s = mandel_from_matrix(velgrad)
    spin_s = skew_from_matrix(velgrad)

  end subroutine element_kinematics_qpt

  !> Element matrix [k_e] = sum_q w_q det J_q [B_q]^T [s_q] [B_q] (R14).
  !>
  !> Arguments:
  !>   elt     (in)  type(element_t)  element descriptor
  !>   coo     (in)  real(rk)(:, :)   (3, num_nodes) node coordinates, sample frame [mm]
  !>   mat_s   (in)  real(rk)(:, :, :)  (6, 6, num_qpts) [s] (secant or tangent), sample frame
  !>                                    [MPa s]
  !>   k_e     (out) real(rk)(:, :)   (num_dofs, num_dofs) element matrix [N s/mm]
  !>   work    (inout) type(element_work_t)  work arrays (element_work_alloc): dshape_s, b, sb
  !>                                   at each point
  !> Locals:
  !>   i_qpt     integer(ik)            quadrature point
  !>   det_jac   real(rk)               det J [mm^3]
  !>   wdet      real(rk)               w_q det J_q [mm^3]
  !>   acc       real(rk)               a sum over the six Mandel components
  !>   i, j, k, l  integer              element dofs (i, j), Mandel components (k, l)
  pure subroutine element_stiffness(elt, coo, mat_s, k_e, work)

    type(element_t), intent(in) :: elt
    real(rk), intent(in) :: coo(:, :)
    real(rk), intent(in) :: mat_s(:, :, :)
    real(rk), intent(out) :: k_e(:, :)
    type(element_work_t), intent(inout) :: work

    integer(ik) :: i_qpt
    real(rk) :: det_jac, wdet, acc
    integer :: i, j, k, l

    ! Sum over the quadrature points of w det J B^T [s] B, by explicit loops: [s][B], then each
    ! entry of B^T([s][B]) added to k_e as it is formed (the products and the order of the sums
    ! of the former matmul expression, without its intermediate matrix; twice as fast).
    k_e = 0.0_rk
    do i_qpt = 1_ik, elt%num_qpts
      call element_geometry_qpt(elt, coo, i_qpt, work%dshape_s, det_jac)
      call mandel_b_matrix(work%dshape_s, work%b)
      ! [s][B].
      do j = 1, size(work%b, 2)
        do k = 1, 6
          acc = 0.0_rk
          do l = 1, 6
            acc = acc + mat_s(k, l, i_qpt)*work%b(l, j)
          end do
          work%sb(k, j) = acc
        end do
      end do
      ! k_e += w det J B^T ([s][B]).
      wdet = elt%qpt_weight(i_qpt)*det_jac
      do j = 1, size(work%b, 2)
        do i = 1, size(work%b, 2)
          acc = 0.0_rk
          do k = 1, 6
            acc = acc + work%b(k, i)*work%sb(k, j)
          end do
          k_e(i, j) = k_e(i, j) + wdet*acc
        end do
      end do
    end do

  end subroutine element_stiffness

  !> Element vector {f_h,e} = sum_q w_q det J_q [B_q]^T [s_q] {h_q} (R14).
  !>
  !> Arguments:
  !>   elt     (in)  type(element_t)    element descriptor
  !>   coo     (in)  real(rk)(:, :)     (3, num_nodes) node coordinates, sample frame [mm]
  !>   mat_s   (in)  real(rk)(:, :, :)  (6, 6, num_qpts) [s] (secant), sample frame [MPa s]
  !>   vec_h   (in)  real(rk)(:, :)     (6, num_qpts) {h}, sample frame [1/s]
  !>   f_e     (out) real(rk)(:)        (num_dofs) element vector [N]
  !>   work    (inout) type(element_work_t)  work arrays (element_work_alloc): dshape_s, b at
  !>                                     each point
  !> Locals:
  !>   i_qpt     integer(ik)            quadrature point
  !>   det_jac   real(rk)               det J [mm^3]
  !>   wdet      real(rk)               w_q det J_q [mm^3]
  !>   sh        real(rk)(6)            [s]{h} [MPa]
  !>   j         integer                element dof
  pure subroutine element_force_h(elt, coo, mat_s, vec_h, f_e, work)

    type(element_t), intent(in) :: elt
    real(rk), intent(in) :: coo(:, :)
    real(rk), intent(in) :: mat_s(:, :, :)
    real(rk), intent(in) :: vec_h(:, :)
    real(rk), intent(out) :: f_e(:)
    type(element_work_t), intent(inout) :: work

    integer(ik) :: i_qpt
    real(rk) :: det_jac, wdet
    real(rk) :: sh(6)
    integer :: j

    ! Sum over the quadrature points of w det J B^T [s] {h}.
    f_e = 0.0_rk
    do i_qpt = 1_ik, elt%num_qpts
      call element_geometry_qpt(elt, coo, i_qpt, work%dshape_s, det_jac)
      call mandel_b_matrix(work%dshape_s, work%b)
      sh = matmul(mat_s(:, :, i_qpt), vec_h(:, i_qpt))
      ! f_e += w det J B^T ([s]{h}), one dot product per dof (the summation order of the former
      ! expression: gfortran's matmul(transpose(b), v) sums in another order when b is an array).
      wdet = elt%qpt_weight(i_qpt)*det_jac
      do j = 1, size(work%b, 2)
        f_e(j) = f_e(j) + wdet*dot_product(work%b(:, j), sh)
      end do
    end do

  end subroutine element_force_h

end module element_matrix_mod
