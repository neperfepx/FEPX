! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Forces and areas of the face sets: the integral of the traction sigma . n over each face set,
!> with each face's stress taken as the stress of its element (at the element's centroid
!> quadrature point, the stress that is written to the output).
!>
!> Theory : F = sum over the faces f of the set of int_f sigma_e(f) . n dA, A = sum of int_f dA,
!>          with n dA = (dx/dxi1 x dx/dxi2) dxi1 dxi2 (element_mod::face_area_normal; faces
!>          oriented out of the body, D37) and the face quadrature of the element type (D36).
!> Source : port of v2 surf_mod.f90::surf_update (constant element stress per face, face
!>          quadrature, sums over the ranks with par_sum). Differences: the stress is a Mandel
!>          6-vector (v2: the 6 components 11 12 13 22 23 33); the face orientation comes from the
!>          mesh reader (D37; v2 orients each face from its element's centroid, CHANGES.md C24).
!>
!> Public routines:
!>   face_forces  forces and areas of every face set of the mesh (collective)
module face_force_mod

  use kinds_mod, only: rk, ik
  use tensor_mod, only: matrix_from_mandel
  use element_mod, only: element_t, face_area_normal
  use mesh_mod, only: mesh_t
  use parallel_mod, only: par_sum

  implicit none

  private

  public :: face_forces

contains

  !> Forces and areas of every face set, summed over the ranks (each face belongs to one local
  !> element of one rank). Collective.
  !>
  !> Arguments:
  !>   mesh      (in)  type(mesh_t)       rank-local mesh (face sets with local node numbers)
  !>   elt       (in)  type(element_t)    element descriptor (face quadrature)
  !>   coo_s     (in)  real(rk)(:, :)     (3, num_nodes_loc) node coordinates, sample frame [mm]
  !>   stress_s  (in)  real(rk)(:, :)     (6, num_elts) element stress, Mandel, sample frame [MPa]
  !>   loads     (out) real(rk)(:, :)     (3, num_fasets) force on each face set, sample frame [N]
  !>   areas     (out) real(rk)(:)        (num_fasets) area of each face set [mm^2]
  !> Locals:
  !>   i_set, i_face, i_qpt  integer(ik)  face set, face and face point
  !>   coo_face  real(rk)(:, :)           (3, face_num_nodes) face node coordinates [mm]
  !>   sigma     real(rk)(3, 3)           element stress as a matrix [MPa]
  !>   normal    real(rk)(3)              area-weighted normal at a point [mm^2]
  !>   load      real(rk)(3)              force of this rank's faces of the set [N]
  !>   area      real(rk)                 area of this rank's faces of the set [mm^2]
  subroutine face_forces(mesh, elt, coo_s, stress_s, loads, areas)

    type(mesh_t), intent(in) :: mesh
    type(element_t), intent(in) :: elt
    real(rk), intent(in) :: coo_s(:, :)
    real(rk), intent(in) :: stress_s(:, :)
    real(rk), intent(out) :: loads(:, :)
    real(rk), intent(out) :: areas(:)

    integer(ik) :: i_set, i_face, i_qpt
    real(rk) :: coo_face(3, elt%face_num_nodes)
    real(rk) :: sigma(3, 3)
    real(rk) :: normal(3)
    real(rk) :: load(3)
    real(rk) :: area

    do i_set = 1_ik, size(mesh%fasets, kind=ik)
      ! This rank's faces of the set: sum of sigma . n dA and of dA over the face points.
      load = 0.0_rk
      area = 0.0_rk
      do i_face = 1_ik, size(mesh%fasets(i_set)%elts, kind=ik)
        coo_face = coo_s(:, mesh%fasets(i_set)%nodes(:, i_face))
        sigma = matrix_from_mandel(stress_s(:, mesh%fasets(i_set)%elts(i_face)))
        do i_qpt = 1_ik, elt%face_num_qpts
          normal = face_area_normal(coo_face, elt%face_dshape_qpt(:, :, i_qpt))
          load = load + elt%face_qpt_weight(i_qpt)*matmul(sigma, normal)
          area = area + elt%face_qpt_weight(i_qpt)*norm2(normal)
        end do
      end do

      ! Sums over the ranks.
      call par_sum(load, loads(:, i_set))
      call par_sum(area, areas(i_set))
    end do

  end subroutine face_forces

end module face_force_mod
