! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Type-independent element descriptor and the element geometry kernels: Jacobian, its
!> determinant and inverse, shape-function derivatives in the sample frame, the Mandel B
!> matrix, element volume, and the area-weighted normal of a face element.
!>
!> The descriptor (element_t) is built once from the tables of the chosen element module
!> (tet10_mod or hex20_mod) and holds everything the solver loops need, including the shape
!> functions and their reference derivatives precomputed at the quadrature points (as v2
!> precomputes them), so that no loop over elements or quadrature points branches on the
!> element type (ROADMAP D36; conventions §6: no per-point dispatch).
!>
!> Theory : standard isoparametric mapping x(xi) = sum_a N_a(xi) x_a:
!>            J(i, j) = d x_i / d xi_j = sum_a x_a(i) d N_a / d xi_j,
!>            d N_a / d x_i = sum_j d N_a / d xi_j (J^-1)(j, i),
!>            integral over the element of f = sum_q w_q f(xi_q) det J(xi_q);
!>          Mandel B matrix: recoupled_formulation.md §6, Eq. (R13), {d} = [B]{V};
!>          face element: dx/dxi1 x dx/dxi2 is the normal of the face's node order, its length
!>          the surface Jacobian (area element dA = |n| dxi1 dxi2).
!> Source : Jacobian, determinant, inverse and global derivatives: port of v2
!>          shape_3d_mod.f90::sfder_hpar (the same cofactor formulas; v2 stores the transpose
!>          of J, jac_rc = d x_c / d xi_r, which has the same determinant). B matrix: new
!>          (R13; v2 builds a deviatoric 5 x 3 form, solveit_evp_mod2.f90::elt_stif_evp).
!>
!> Public types:
!>   element_t   descriptor of one element type (components documented at the type)
!>
!> Public routines:
!>   element_init          build the descriptor of a Gmsh element type (11: tet10, 17: hex20)
!>   element_jacobian      J = dx/dxi from node coordinates and reference derivatives
!>   jacobian_det_inv      determinant and inverse of J (cofactors, as v2)
!>   shape_deriv_sample    d N / d x (sample frame) from reference derivatives and J^-1
!>   element_geometry_qpt  d N / d x and det J at one quadrature point of one element
!>   mandel_b_matrix       Mandel B matrix [B] (6 x 3 num_nodes), Eq. (R13)
!>   element_volume        element volume (quadrature of det J)
!>   face_area_normal      dx/dxi1 x dx/dxi2 of a face element at one point
module element_mod

  use kinds_mod, only: rk, ik
  use constants_mod, only: inv_sqrt2
  use error_mod, only: fatal_error
  use tet10_mod, only: tet10_gmsh_type, tet10_num_nodes, tet10_num_qpts, tet10_qpt_coo, &
    & tet10_qpt_weight, tet10_shape, tet10_shape_deriv, tet10_face_num_nodes, &
    & tet10_face_num_qpts, tet10_face_qpt_coo, tet10_face_qpt_weight, tet10_face_shape, &
    & tet10_face_shape_deriv, tet10_qpt_centroid
  use hex20_mod, only: hex20_gmsh_type, hex20_num_nodes, hex20_num_qpts, hex20_qpt_coo, &
    & hex20_qpt_weight, hex20_shape, hex20_shape_deriv, hex20_face_num_nodes, &
    & hex20_face_num_qpts, hex20_face_qpt_coo, hex20_face_qpt_weight, hex20_face_shape, &
    & hex20_face_shape_deriv, hex20_qpt_centroid

  implicit none

  private

  public :: element_t
  public :: element_init, element_jacobian, jacobian_det_inv, shape_deriv_sample
  public :: element_geometry_qpt, mandel_b_matrix, element_volume, face_area_normal

  !> Descriptor of one element type (reference element; no mesh data).
  !>
  !> Components (sizes; units):
  !>   gmsh_type        integer(ik)    Gmsh element type code (11 or 17)
  !>   name             character(:)   "tet10" or "hex20"
  !>   num_nodes        integer(ik)    nodes per element (Gmsh order, D35)
  !>   num_dofs         integer(ik)    3 num_nodes (dof 3 (a - 1) + i: node a, component i)
  !>   num_qpts         integer(ik)    volume quadrature points
  !>   qpt_centroid     integer(ik)    the point at the element centroid (elemental output)
  !>   qpt_coo          real(rk)(3, num_qpts)            reference coordinates [-]
  !>   qpt_weight       real(rk)(num_qpts)               weights (sum: reference volume) [-]
  !>   shape_qpt        real(rk)(num_nodes, num_qpts)    N_a at each quadrature point [-]
  !>   dshape_qpt       real(rk)(3, num_nodes, num_qpts) d N_a / d xi_i at each point [-]
  !>   face_num_nodes   integer(ik)    nodes per face element (6 or 8)
  !>   face_num_qpts    integer(ik)    face quadrature points (7 or 9)
  !>   face_qpt_coo     real(rk)(2, face_num_qpts)       face reference coordinates [-]
  !>   face_qpt_weight  real(rk)(face_num_qpts)          face weights (sum: reference area) [-]
  !>   face_shape_qpt   real(rk)(face_num_nodes, face_num_qpts)     face N_a [-]
  !>   face_dshape_qpt  real(rk)(2, face_num_nodes, face_num_qpts)  face d N_a / d xi_i [-]
  type :: element_t
    integer(ik) :: gmsh_type = 0_ik
    character(len=:), allocatable :: name
    integer(ik) :: num_nodes = 0_ik
    integer(ik) :: num_dofs = 0_ik
    integer(ik) :: num_qpts = 0_ik
    integer(ik) :: qpt_centroid = 0_ik
    real(rk), allocatable :: qpt_coo(:, :)
    real(rk), allocatable :: qpt_weight(:)
    real(rk), allocatable :: shape_qpt(:, :)
    real(rk), allocatable :: dshape_qpt(:, :, :)
    integer(ik) :: face_num_nodes = 0_ik
    integer(ik) :: face_num_qpts = 0_ik
    real(rk), allocatable :: face_qpt_coo(:, :)
    real(rk), allocatable :: face_qpt_weight(:)
    real(rk), allocatable :: face_shape_qpt(:, :)
    real(rk), allocatable :: face_dshape_qpt(:, :, :)
  end type element_t

contains

  !> Build the descriptor of an element type from its element module's tables, with the shape
  !> functions and reference derivatives evaluated at every (volume and face) quadrature point.
  !>
  !> Arguments:
  !>   elt        (out) type(element_t)  descriptor
  !>   gmsh_type  (in)  integer(ik)      Gmsh element type code (11: tet10, 17: hex20)
  !> Locals:
  !>   i_qpt      integer(ik)            quadrature point
  !>   code_text  character(24)          gmsh_type as text, for the error message
  subroutine element_init(elt, gmsh_type)

    type(element_t), intent(out) :: elt
    integer(ik), intent(in) :: gmsh_type

    integer(ik) :: i_qpt
    character(len=24) :: code_text

    ! Sizes and quadrature rules of the element type.
    elt%gmsh_type = gmsh_type
    select case (gmsh_type)
    case (tet10_gmsh_type)
      elt%name = 'tet10'
      elt%num_nodes = tet10_num_nodes
      elt%num_qpts = tet10_num_qpts
      elt%qpt_centroid = tet10_qpt_centroid
      elt%qpt_coo = tet10_qpt_coo
      elt%qpt_weight = tet10_qpt_weight
      elt%face_num_nodes = tet10_face_num_nodes
      elt%face_num_qpts = tet10_face_num_qpts
      elt%face_qpt_coo = tet10_face_qpt_coo
      elt%face_qpt_weight = tet10_face_qpt_weight
    case (hex20_gmsh_type)
      elt%name = 'hex20'
      elt%num_nodes = hex20_num_nodes
      elt%num_qpts = hex20_num_qpts
      elt%qpt_centroid = hex20_qpt_centroid
      elt%qpt_coo = hex20_qpt_coo
      elt%qpt_weight = hex20_qpt_weight
      elt%face_num_nodes = hex20_face_num_nodes
      elt%face_num_qpts = hex20_face_num_qpts
      elt%face_qpt_coo = hex20_face_qpt_coo
      elt%face_qpt_weight = hex20_face_qpt_weight
    case default
      write (code_text, '(i0)') gmsh_type
      call fatal_error('element_init', 'unsupported element type (Gmsh type code ' // &
        & trim(code_text) // ', supported: 11 and 17)')
    end select
    elt%num_dofs = 3_ik*elt%num_nodes

    ! Shape functions and reference derivatives at the volume quadrature points.
    allocate (elt%shape_qpt(elt%num_nodes, elt%num_qpts))
    allocate (elt%dshape_qpt(3, elt%num_nodes, elt%num_qpts))
    do i_qpt = 1_ik, elt%num_qpts
      if (gmsh_type == tet10_gmsh_type) then
        elt%shape_qpt(:, i_qpt) = tet10_shape(elt%qpt_coo(:, i_qpt))
        elt%dshape_qpt(:, :, i_qpt) = tet10_shape_deriv(elt%qpt_coo(:, i_qpt))
      else
        elt%shape_qpt(:, i_qpt) = hex20_shape(elt%qpt_coo(:, i_qpt))
        elt%dshape_qpt(:, :, i_qpt) = hex20_shape_deriv(elt%qpt_coo(:, i_qpt))
      end if
    end do

    ! Face shape functions and reference derivatives at the face quadrature points.
    allocate (elt%face_shape_qpt(elt%face_num_nodes, elt%face_num_qpts))
    allocate (elt%face_dshape_qpt(2, elt%face_num_nodes, elt%face_num_qpts))
    do i_qpt = 1_ik, elt%face_num_qpts
      if (gmsh_type == tet10_gmsh_type) then
        elt%face_shape_qpt(:, i_qpt) = tet10_face_shape(elt%face_qpt_coo(:, i_qpt))
        elt%face_dshape_qpt(:, :, i_qpt) = tet10_face_shape_deriv(elt%face_qpt_coo(:, i_qpt))
      else
        elt%face_shape_qpt(:, i_qpt) = hex20_face_shape(elt%face_qpt_coo(:, i_qpt))
        elt%face_dshape_qpt(:, :, i_qpt) = hex20_face_shape_deriv(elt%face_qpt_coo(:, i_qpt))
      end if
    end do

  end subroutine element_init

  !> Jacobian of the isoparametric map at one point.
  !>
  !> Theory : module header, J(i, j) = sum_a x_a(i) d N_a / d xi_j.
  !> Source : port of v2 shape_3d_mod.f90::sfder_hpar (Jacobian sums; v2 stores J^T).
  !>
  !> Arguments:
  !>   coo     (in)  real(rk)(:, :)  (3, num_nodes) node coordinates, sample frame [mm]
  !>   dshape  (in)  real(rk)(:, :)  (3, num_nodes) d N_a / d xi_i at the point [-]
  !>   jac     (out) real(rk)(3, 3)  J(i, j) = d x_i / d xi_j [mm]
  pure function element_jacobian(coo, dshape) result(jac)

    real(rk), intent(in) :: coo(:, :)
    real(rk), intent(in) :: dshape(:, :)
    real(rk) :: jac(3, 3)

    ! J = coo dshape^T: sum over the nodes.
    jac = matmul(coo, transpose(dshape))

  end function element_jacobian

  !> Determinant and inverse of a Jacobian by cofactors. The inverse is computed only for a
  !> positive determinant (a valid element); otherwise it is set to zero and the caller, which
  !> checks det_jac > 0, reports the invalid element.
  !>
  !> Source : port of v2 shape_3d_mod.f90::sfder_hpar (determinant by the rule of Sarrus and
  !>          inverse = cofactor matrix transposed / determinant, the same products).
  !>
  !> Arguments:
  !>   jac      (in)  real(rk)(3, 3)  Jacobian J [mm]
  !>   det_jac  (out) real(rk)        det J [mm^3]
  !>   jac_inv  (out) real(rk)(3, 3)  J^-1 [1/mm] (zero if det J <= 0)
  pure subroutine jacobian_det_inv(jac, det_jac, jac_inv)

    real(rk), intent(in) :: jac(3, 3)
    real(rk), intent(out) :: det_jac
    real(rk), intent(out) :: jac_inv(3, 3)

    ! Determinant (rule of Sarrus).
    det_jac = jac(1, 1)*jac(2, 2)*jac(3, 3) + jac(1, 2)*jac(2, 3)*jac(3, 1) &
      & + jac(1, 3)*jac(2, 1)*jac(3, 2) - jac(1, 1)*jac(2, 3)*jac(3, 2) &
      & - jac(1, 2)*jac(2, 1)*jac(3, 3) - jac(1, 3)*jac(2, 2)*jac(3, 1)

    ! Inverse: transposed cofactor matrix divided by the determinant (valid elements only).
    if (det_jac > 0.0_rk) then
      jac_inv(1, 1) = jac(2, 2)*jac(3, 3) - jac(2, 3)*jac(3, 2)
      jac_inv(1, 2) = jac(1, 3)*jac(3, 2) - jac(1, 2)*jac(3, 3)
      jac_inv(1, 3) = jac(1, 2)*jac(2, 3) - jac(1, 3)*jac(2, 2)
      jac_inv(2, 1) = jac(2, 3)*jac(3, 1) - jac(2, 1)*jac(3, 3)
      jac_inv(2, 2) = jac(1, 1)*jac(3, 3) - jac(1, 3)*jac(3, 1)
      jac_inv(2, 3) = jac(1, 3)*jac(2, 1) - jac(1, 1)*jac(2, 3)
      jac_inv(3, 1) = jac(2, 1)*jac(3, 2) - jac(2, 2)*jac(3, 1)
      jac_inv(3, 2) = jac(1, 2)*jac(3, 1) - jac(1, 1)*jac(3, 2)
      jac_inv(3, 3) = jac(1, 1)*jac(2, 2) - jac(1, 2)*jac(2, 1)
      jac_inv = jac_inv/det_jac
    else
      jac_inv = 0.0_rk
    end if

  end subroutine jacobian_det_inv

  !> Shape-function derivatives with respect to the sample coordinates at one point.
  !>
  !> Theory : chain rule, d N_a / d x_i = sum_j d N_a / d xi_j (J^-1)(j, i).
  !> Source : port of v2 shape_3d_mod.f90::sfder_hpar (dndx, dndy, dndz).
  !>
  !> Arguments:
  !>   dshape    (in)  real(rk)(:, :)          (3, num_nodes) d N_a / d xi_i [-]
  !>   jac_inv   (in)  real(rk)(3, 3)          J^-1 [1/mm]
  !>   dshape_s  (out) real(rk)(3, num_nodes)  d N_a / d x_i, sample frame [1/mm]
  pure function shape_deriv_sample(dshape, jac_inv) result(dshape_s)

    real(rk), intent(in) :: dshape(:, :)
    real(rk), intent(in) :: jac_inv(3, 3)
    real(rk) :: dshape_s(3, size(dshape, 2))

    ! dshape_s = J^-T dshape.
    dshape_s = matmul(transpose(jac_inv), dshape)

  end function shape_deriv_sample

  !> Sample-frame shape-function derivatives and Jacobian determinant at one quadrature point
  !> of one element (the per-point geometry of the element loops).
  !>
  !> Arguments:
  !>   elt       (in)  type(element_t)         element descriptor
  !>   coo       (in)  real(rk)(:, :)          (3, num_nodes) node coordinates, sample frame [mm]
  !>   i_qpt     (in)  integer(ik)             quadrature point
  !>   dshape_s  (out) real(rk)(:, :)          (3, num_nodes) d N_a / d x_i, sample frame [1/mm]
  !>   det_jac   (out) real(rk)                det J at the point [mm^3]
  !> Locals:
  !>   jac       real(rk)(3, 3)                Jacobian [mm]
  !>   jac_inv   real(rk)(3, 3)                its inverse [1/mm]
  pure subroutine element_geometry_qpt(elt, coo, i_qpt, dshape_s, det_jac)

    type(element_t), intent(in) :: elt
    real(rk), intent(in) :: coo(:, :)
    integer(ik), intent(in) :: i_qpt
    real(rk), intent(out) :: dshape_s(:, :)
    real(rk), intent(out) :: det_jac

    real(rk) :: jac(3, 3)
    real(rk) :: jac_inv(3, 3)

    ! Jacobian, its determinant and inverse, then the chain rule.
    jac = element_jacobian(coo, elt%dshape_qpt(:, :, i_qpt))
    call jacobian_det_inv(jac, det_jac, jac_inv)
    dshape_s = shape_deriv_sample(elt%dshape_qpt(:, :, i_qpt), jac_inv)

  end subroutine element_geometry_qpt

  !> Mandel B matrix of an element at one point: {d} = [B]{V} with {d} the Mandel 6-vector
  !> (11 22 33 23 13 12) of the deformation rate and {V} the element velocities, node by node
  !> (column 3 (a - 1) + i: node a, component i). For node a (N,x = d N_a / d x, ...):
  !>   | N,x      0        0      |
  !>   | 0        N,y      0      |
  !>   | 0        0        N,z    |
  !>   | 0        N,z/rt2  N,y/rt2|
  !>   | N,z/rt2  0        N,x/rt2|
  !>   | N,y/rt2  N,x/rt2  0      |
  !>
  !> Theory : recoupled_formulation.md §6, Eq. (R13).
  !> Source : new (R13).
  !>
  !> Arguments:
  !>   dshape_s  (in)  real(rk)(:, :)              (3, num_nodes) d N_a / d x_i [1/mm]
  !>   b         (out) real(rk)(6, 3 num_nodes)    Mandel B matrix [1/mm] (the caller's array: a
  !>                                               subroutine, so that no result temporary is
  !>                                               allocated per point, ROADMAP D100)
  !> Locals:
  !>   a         integer                           node
  !>   col       integer                           first column of node a
  pure subroutine mandel_b_matrix(dshape_s, b)

    real(rk), intent(in) :: dshape_s(:, :)
    real(rk), intent(out) :: b(:, :)

    integer :: a
    integer :: col

    b = 0.0_rk
    do a = 1, size(dshape_s, 2)
      col = 3*(a - 1)
      ! Normal components 11, 22, 33.
      b(1, col + 1) = dshape_s(1, a)
      b(2, col + 2) = dshape_s(2, a)
      b(3, col + 3) = dshape_s(3, a)
      ! Shear components 23, 13, 12 (Mandel: sqrt2 times the tensor component).
      b(4, col + 2) = inv_sqrt2*dshape_s(3, a)
      b(4, col + 3) = inv_sqrt2*dshape_s(2, a)
      b(5, col + 1) = inv_sqrt2*dshape_s(3, a)
      b(5, col + 3) = inv_sqrt2*dshape_s(1, a)
      b(6, col + 1) = inv_sqrt2*dshape_s(2, a)
      b(6, col + 2) = inv_sqrt2*dshape_s(1, a)
    end do

  end subroutine mandel_b_matrix

  !> Volume of one element: sum over the quadrature points of w_q det J (exact for affine
  !> elements; the quadrature value otherwise).
  !>
  !> Arguments:
  !>   elt     (in)  type(element_t)  element descriptor
  !>   coo     (in)  real(rk)(:, :)   (3, num_nodes) node coordinates, sample frame [mm]
  !>   vol     (out) real(rk)         volume [mm^3]
  !> Locals:
  !>   i_qpt    integer(ik)           quadrature point
  !>   det_jac  real(rk)              det J [mm^3]
  !>   jac_inv  real(rk)(3, 3)        unused inverse [1/mm]
  pure function element_volume(elt, coo) result(vol)

    type(element_t), intent(in) :: elt
    real(rk), intent(in) :: coo(:, :)
    real(rk) :: vol

    integer(ik) :: i_qpt
    real(rk) :: det_jac
    real(rk) :: jac_inv(3, 3)

    ! Quadrature of det J.
    vol = 0.0_rk
    do i_qpt = 1_ik, elt%num_qpts
      call jacobian_det_inv(element_jacobian(coo, elt%dshape_qpt(:, :, i_qpt)), det_jac, jac_inv)
      vol = vol + elt%qpt_weight(i_qpt)*det_jac
    end do

  end function element_volume

  !> Area-weighted normal of a face element at one point: n = dx/dxi1 x dx/dxi2, whose
  !> direction is the normal of the face's node order (out of the body for faces read from
  !> $Fasets, D37) and whose length is the surface Jacobian (dA = |n| dxi1 dxi2).
  !>
  !> Theory : module header (face element).
  !> Source : as v2 surf_mod.f90::surf_update (tangents and their cross product).
  !>
  !> Arguments:
  !>   coo     (in)  real(rk)(:, :)  (3, face_num_nodes) face node coordinates, sample frame [mm]
  !>   dshape  (in)  real(rk)(:, :)  (2, face_num_nodes) face d N_a / d xi_i at the point [-]
  !>   normal  (out) real(rk)(3)     dx/dxi1 x dx/dxi2 [mm^2]
  !> Locals:
  !>   tang    real(rk)(3, 2)        tangents dx/dxi1, dx/dxi2 [mm]
  pure function face_area_normal(coo, dshape) result(normal)

    real(rk), intent(in) :: coo(:, :)
    real(rk), intent(in) :: dshape(:, :)
    real(rk) :: normal(3)

    real(rk) :: tang(3, 2)

    ! Tangents, then their cross product.
    tang = matmul(coo, transpose(dshape))
    normal = [tang(2, 1)*tang(3, 2) - tang(3, 1)*tang(2, 2), &
      & tang(3, 1)*tang(1, 2) - tang(1, 1)*tang(3, 2), &
      & tang(1, 1)*tang(2, 2) - tang(2, 1)*tang(1, 2)]

  end function face_area_normal

end module element_mod
