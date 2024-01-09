! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module mesh_type_mod

! Module to define types

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use parallel_mod
  use gather_scatter_mod
  use types_mod

  implicit none

  public

  contains

  subroutine mesh_surfaceareas(mesh, coo, areas)

    ! Compute the areas of the surfaces of a mesh

    !---------------------------------------------------------------------------

    type(mesh_type), intent(in) :: mesh
    real(rk), intent(in), optional :: coo(:)
    real(rk), intent(out) :: areas(mesh%num_fasets)

    ! Locals:

    integer :: i
    integer :: iel
    integer :: jq
    integer :: jn
    integer :: jdof
    real(rk) :: tangent(3, 2, 7)
    real(rk) :: normal(3, 7)
    real(rk) :: nmag
    real(rk) :: sjac(7)
    real(rk) :: p_area, area

    !---------------------------------------------------------------------------

    areas = 0.0d0

    do i = 1, mesh%num_fasets
      if (present(coo)) then
        call part_gather(mesh%fasets(i)%coo, coo, mesh%fasets(i)%conn3d, &
          & mesh%fasets(i)%tr)

      else
        call part_gather(mesh%fasets(i)%coo, mesh%coo, mesh%fasets(i)%conn3d, &
          & mesh%fasets(i)%tr)

      end if

      ! Now compute normal vectors at quadrature points

      p_area = 0.0d0

      do iel = mesh%fasets(i)%semin, mesh%fasets(i)%semax
        tangent = 0.0d0

        do jq = 1, 7
          do jn = 1, sftype
            jdof = 3*(jn - 1) + 1

            ! First tangent vector

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

          ! Now compute normals.

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

        do jq = 1, 7
          p_area = p_area + wt2d(jq)*sjac(jq)
        end do
      end do !iel

      call par_sum(p_area, area)

      areas(i) = area
    end do !mesh%num_fasets

    return

  end subroutine mesh_surfaceareas

  subroutine mesh_bbox(mesh, bbox, results)

    type(mesh_type), intent(in) :: mesh
    real(rk), intent(out) :: bbox(2, 3)
    type(results_type), intent(in), optional :: results

    integer :: i, j, dir
    real(rk), allocatable :: temp(:)

    allocate (temp(size(mesh%coo)/3 + 1))

    do dir = 1, 3
      temp = 0.0d0
      j = 1

      if (.not. present (results)) then
        do i = dof_sub, dof_sup, 3
          temp(j) = mesh%coo(i + dir - 1)
          j = j + 1
        end do

      else
        do i = dof_sub, dof_sup, 3
          temp(j) = results%coo(i + dir - 1)
          j = j + 1
        end do

      end if

      call par_min(minval(temp), bbox(1,dir))
      call par_max(maxval(temp), bbox(2,dir))
    end do

    deallocate (temp)

  end subroutine mesh_bbox

  subroutine mesh_lengths(mesh, length, results)

    type(mesh_type), intent(in) :: mesh
    real(rk), intent(out) :: length(3)
    type(results_type), intent(in), optional :: results

    integer :: i
    real(rk) :: bbox (2, 3)

    call mesh_bbox(mesh, bbox, results)

    do i = 1, 3
      length(i) = bbox(2, i) - bbox(1, i)
    end do

  end subroutine mesh_lengths

  subroutine mesh_length(mesh, dir, length, results)

    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: dir
    real(rk), intent(out) :: length
    type(results_type), intent(in), optional :: results

    real(rk) :: bbox (2, 3)

    if (dir .lt. 1 .or. dir .gt. 3) then
      write (*,*) "dir undefined"
      call par_quit ("Error")
    end if

    call mesh_bbox(mesh, bbox, results)

    length = bbox(2, dir) - bbox(1, dir)

  end subroutine mesh_length

  subroutine calc_mesh_dim(results, length, indx, indy, indz)

    ! Calculate mesh dimensions for triaxial loading.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! length:
    ! indx:
    ! indy:
    ! indz:

    type(results_type), intent(in) :: results
    real(rk), intent(out) :: length(3)
    integer, intent(in) :: indx((dof_sup - dof_sub + 1)/3)
    integer, intent(in) :: indy((dof_sup - dof_sub + 1)/3)
    integer, intent(in) :: indz((dof_sup - dof_sub + 1)/3)

    ! Locals:
    ! part_xmax,ymax,zmax:
    ! length_x,y,z:

    real(rk) :: part_xmax, part_ymax, part_zmax
    real(rk) :: length_x, length_y, length_z

    !---------------------------------------------------------------------------

    ! Locate maximum coordinate values of the domain.

    part_xmax = maxval(results%coo(indx))
    part_ymax = maxval(results%coo(indy))
    part_zmax = maxval(results%coo(indz))

    call par_max(part_xmax, length_x)
    call par_max(part_ymax, length_y)
    call par_max(part_zmax, length_z)

    ! Store returned maximums into output array.

    length(1) = length_x
    length(2) = length_y
    length(3) = length_z

    return

  end subroutine calc_mesh_dim

  end module mesh_type_mod
