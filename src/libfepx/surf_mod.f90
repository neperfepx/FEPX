! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module surface_mod

! Module for using surface information.

! Contains function:
! surf_alloc: Allocate space to enough elements
! surf_init: Initialize surface arrays
! surf_update: Update surface information

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use gather_scatter_mod
  use quadrature_mod
  use shape_2d_mod
  use types_mod
  use parallel_mod

  implicit none

  private

  public:: surf_alloc, surf_init, surf_update

! Public

! The type `surface_section' is for use in parallel codes, and contains a
!   section of the surface mesh.

contains

  function surf_alloc(type, semin, semax, surf) result(status)

    ! Allocate space to enough elements.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! type: Surface element type, now only q6 is available
    ! semin, semax: Range of surface element numbers
    ! surf: The surface to be allocated

    integer :: type
    integer :: semin
    integer :: semax
    type(surface_section) :: surf
    integer :: status

    !---------------------------------------------------------------------------

    status = 0

    !if (type .ne. 4) then
    ! tsh for 3-node triangle
    !if (type .ne. 3) then
    ! tsh for 6-node triangle
    if (type .ne. 6) then
      status = 1

      return
    end if
    surf%type = type ! type=6

    if (semax .ge. semin) then
      surf%semin = semin
      surf%semax = semax

    else
      status = 2

      return
    end if

    allocate (surf%conn(type, semin:semax), stat=status)
    allocate (surf%conn3d(3*type, semin:semax), stat=status)
    allocate (surf%econn(6, semin:semax), stat=status)
    allocate (surf%coo(3*type, semin:semax), stat=status)
    allocate (surf%elt(semin:semax), stat=status)

    return

  end function surf_alloc

  !===========================================================================

  subroutine surf_init()

    ! Initialize surface arrays

    !---------------------------------------------------------------------------

    !   Locals:

    integer :: istat

    !---------------------------------------------------------------------------

    call sf2d(sftype, 7, qp2d, sfqp2d, sftype, istat)
    call sf2dg(sftype, 7, qp2d, sfgqp2d, sftype, istat)

    return

  end subroutine surf_init

  !===========================================================================

  subroutine surf_update(mesh, coo, sig_all, load_out, area_out)

    ! Update surface information.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! coo: Coordinates
    ! sig_all: Stresses
    ! load_out: Output loads (enters as 0)
    ! area_out: Output area (enters as 0)

    type(mesh_type), intent(in) :: mesh
    real(rk), intent(in) :: coo(:)
    real(rk), intent(in) :: sig_all(:)
    real(rk), intent(inout) :: load_out(mesh%num_fasets, 3)
    real(rk), intent(inout) :: area_out(mesh%num_fasets)

    ! Locals:

    integer :: i
    integer :: iel
    integer :: jq
    integer :: jn
    integer :: jdof
    integer :: iload
    real(rk) :: tangent(3, 2, 7)
    real(rk) :: normal(3, 7)
    real(rk) :: nmag
    real(rk) :: sjac(7)
    real(rk), allocatable:: sig_se(:, :)
    real(rk) :: load(3)
    real(rk) :: p_total_load(3)
    real(rk) :: total_load(3)
    real(rk) :: p_area
    real(rk) :: area

    !----------------------------------------------------------------------

    load_out = 0.0d0
    area_out = 0.0d0

    do i = 1, mesh%num_fasets
      call part_gather(mesh%fasets(i)%coo, coo, mesh%fasets(i)%conn3d, &
          & mesh%fasets(i)%tr)

      allocate (sig_se(6, mesh%fasets(i)%semin:mesh%fasets(i)%semax))

      call part_gather(sig_se, sig_all, mesh%fasets(i)%econn, mesh%fasets(i)%etr)

      ! Now compute normal vectors at quadrature points.

      p_total_load = 0.0d0
      p_area = 0.0d0

      do iel = mesh%fasets(i)%semin, mesh%fasets(i)%semax
        tangent = 0.0d0

        do jq = 1, 7
          do jn = 1, sftype !number of nodes in the surface element
            jdof = 3*(jn - 1) + 1

            ! First tangent vector.

            tangent(1, 1, jq) = tangent(1, 1, jq) + &
                & mesh%fasets(i)%coo(jdof, iel)*sfgqp2d(1, jn, jq)
            tangent(2, 1, jq) = tangent(2, 1, jq) + &
                & mesh%fasets(i)%coo(jdof + 1, iel)*sfgqp2d(1, jn, jq)
            tangent(3, 1, jq) = tangent(3, 1, jq) + &
                & mesh%fasets(i)%coo(jdof + 2, iel)*sfgqp2d(1, jn, jq)

            ! Second tangent vector.

            tangent(1, 2, jq) = tangent(1, 2, jq) + &
                 & mesh%fasets(i)%coo(jdof, iel)*sfgqp2d(2, jn, jq)
            tangent(2, 2, jq) = tangent(2, 2, jq) + &
                 & mesh%fasets(i)%coo(jdof + 1, iel)*sfgqp2d(2, jn, jq)
            tangent(3, 2, jq) = tangent(3, 2, jq) + &
                 & mesh%fasets(i)%coo(jdof + 2, iel)*sfgqp2d(2, jn, jq)
          end do

          ! Now compute normals

          normal(1, jq) = tangent(2, 1, jq)*tangent(3, 2, jq) - &
              & tangent(3, 1, jq)*tangent(2, 2, jq)
          normal(2, jq) = tangent(3, 1, jq)*tangent(1, 2, jq) - &
              & tangent(1, 1, jq)*tangent(3, 2, jq)
          normal(3, jq) = tangent(1, 1, jq)*tangent(2, 2, jq) - &
              & tangent(2, 1, jq)*tangent(1, 2, jq)
          nmag = dsqrt(normal(1, jq)*normal(1, jq) + &
              & normal(2, jq)*normal(2, jq) + &
              & normal(3, jq)*normal(3, jq))

          if (nmag > 0) then
            normal(:, jq) = normal(:, jq)/nmag
            sjac(jq) = nmag

          else
            call par_quit('Error  :     > Surface normal zero &
                &magnitude.')
          end if
        end do

        load = 0.0d0

        do jq = 1, 7
          load(1) = load(1) + wt2d(jq)*sjac(jq)*(sig_se(1, iel)* &
              & normal(1, jq) + sig_se(2, iel)*normal(2, jq) + &
              & sig_se(3, iel)*normal(3, jq))
          load(2) = load(2) + wt2d(jq)*sjac(jq)*(sig_se(2, iel)* &
              & normal(1, jq) + sig_se(4, iel)*normal(2, jq) + &
              & sig_se(5, iel)*normal(3, jq))
          load(3) = load(3) + wt2d(jq)*sjac(jq)*(sig_se(3, iel)* &
              & normal(1, jq) + sig_se(5, iel)*normal(2, jq) + &
              & sig_se(6, iel)*normal(3, jq))
          p_area = p_area + wt2d(jq)*sjac(jq)
        end do

        p_total_load = p_total_load + load
      end do !iel

      do iload = 1, 3
        call par_sum(p_total_load(iload), total_load(iload))
      end do

      call par_sum(p_area, area)

      load_out(i, :) = total_load
      area_out(i) = area

      deallocate (sig_se)
    end do !mesh%num_fasets

    return

  end subroutine surf_update

end module surface_mod
