! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module boundary_conditions_mod3

  use general_mod
  use mesh_type_mod
  use types_mod
  use, intrinsic :: iso_fortran_env, only: rk => real64

  use parallel_mod, only: par_min, par_max, par_quit
  use types_mod
  use utils_mod

  implicit none

contains

  ! ============================================================================

  !> Calculate general boundary conditions / velocity
  subroutine calc_bcs_general_vel(mesh, loading_options, id, loading)

    type(mesh_type), intent(in) :: mesh
    type(loading_options_type), intent(inout) :: loading_options
    integer, intent(in) :: id
    type(loading_type), intent(inout) :: loading

    integer :: j, k, num_nodes, status, uniqnode, dof, offset
    integer, allocatable :: nodes(:)

    read (loading_options%bc_nset(id), *, iostat=status) uniqnode

    if (status == 0) then
      num_nodes = 1
      allocate (nodes(num_nodes))
      nodes(1) = uniqnode

    else
      do j = 1, size(mesh%nsets)
        if (loading_options%bc_nset(id) .eq. mesh%nsets(j)%nset_label) then
          num_nodes = mesh%nsets(j)%num_nset_nodes
          allocate (nodes(num_nodes))
          nodes = mesh%nsets(j)%nset_nodes
          exit
        end if
      end do
    end if

    if (.not. allocated (nodes)) then
      write (*,'(a,a,a)') "Error  : ", trim(loading_options%bc_nset(id)), ":"
      call par_quit("         Failed to find node or node set.")
    end if

    ! Define the direction where bc are applied
    select case (loading_options%bc_dir(id))

    case ('x', 'vx')
      offset = 1
    case ('y', 'vy')
      offset = 2
    case ('z', 'vz')
      offset = 3

    end select

    ! Affect the bc defined value to each nodes of the current nset
    do k = 1, num_nodes

      dof = (nodes(k) - 1) * 3 + offset
      if (dof .ge. dof_sub .and. dof .le. dof_sup) then
        if ((loading%bcs_vel_defined(dof) .eqv. .true.) .and. &
          & (loading%bcs_vel(dof) .ne. loading_options%bc_vel(id))) then
          write (*,*) "Error  :     >  Node ", nodes(k), loading%bcs_vel(dof), loading_options%bc_vel(id)
          call par_quit("Error  :     >  At least one node has conflicting boundary conditions")

        else
          loading%bcs_vel_defined(dof) = .true.
          loading%bcs_vel(dof) = loading_options%bc_vel(id)
        end if
      end if
    end do

    deallocate(nodes)

  end subroutine calc_bcs_general_vel

  !> Calculate general boundary conditions / strain rate
  subroutine calc_bcs_general_velgrad(mesh, loading_options, L, L_defined, L_defined_bcindex, loading)

    type(mesh_type), intent(in) :: mesh
    type(loading_options_type), intent(inout) :: loading_options
    real(rk) :: L(3, 3)
    logical :: L_defined(3,3)
    integer :: L_defined_bcindex(3,3)
    type(loading_type), intent(inout) :: loading
    character(len=2) :: nset1, nset2
    character(len=1) :: dir

    integer :: i, j, k, num_components
    real(rk) :: lengths(3), strainrate

    num_components = 0
    do i = 1, 3
      do j = 1, 3
        if (L_defined(i, j)) then
          loading%loading_direction = i
          loading%loading_face = j

          strainrate = L(i,j)

          nset1 = achar (loading%loading_face + ichar ('x') - 1) // '0'
          nset2 = achar (loading%loading_face + ichar ('x') - 1) // '1'

          if (loading_options%bc_type(L_defined_bcindex(i, j)) .eq. "grip") then

            call mesh_lengths (mesh, lengths)
            loading%gage_length = lengths(loading%loading_face)

            loading_options%bc_var = "vel"
            loading_options%num_bcs = 6
            loading_options%bc_vel(1:6)  = 0.0d0

            loading_options%bc_nset(1:3) = nset1
            loading_options%bc_dir(1)  = "x"
            loading_options%bc_dir(2)  = "y"
            loading_options%bc_dir(3)  = "z"

            loading_options%bc_nset(4:6) = nset2
            loading_options%bc_dir(4)  = "x"
            loading_options%bc_dir(5)  = "y"
            loading_options%bc_dir(6)  = "z"
            loading_options%bc_vel(3 + loading%loading_direction) =  strainrate*loading%gage_length

          else if (loading_options%bc_type(L_defined_bcindex(i, j)) .eq. "minimal") then

            call mesh_lengths (mesh, lengths)
            loading%gage_length = lengths(loading%loading_face)
            dir = achar (loading%loading_direction + ichar ('x') - 1)

            loading_options%bc_var = "vel"
            loading_options%num_bcs = 2

            loading_options%bc_nset(1) = nset1
            loading_options%bc_dir(1)  = dir
            loading_options%bc_vel(1)  = 0.0d0

            loading_options%bc_nset(2) = nset2
            loading_options%bc_dir(2)  = dir
            loading_options%bc_vel(2) =  strainrate*loading%gage_length

          else
            call par_quit ("Error  :   - Failed to apply bcs.")

          end if

          do k = 1, loading_options%num_bcs
            call calc_bcs_general_vel(mesh, loading_options, k, loading)
          end do

        end if
      end do
    end do

  end subroutine calc_bcs_general_velgrad

  ! ============================================================================

end module boundary_conditions_mod3
