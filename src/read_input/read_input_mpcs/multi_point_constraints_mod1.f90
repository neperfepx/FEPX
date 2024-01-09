! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

!> Module to calculate essential multi-point constraints.
module multi_point_constraints_mod

  use general_mod
  use types_mod
  use, intrinsic :: iso_fortran_env, only: rk => real64
  use parallel_mod
  use gather_scatter_mod
  use multi_point_constraints_mod2
  use periodicity_mod1

  implicit none

  public

contains

  subroutine read_mpcs(mesh, loading_options, loading)

    ! Defines the primary-secondary relation between node dofs:
    ! The primary and secondary dof velocity are linked by:
    ! [velocity of secondary dof] = loading%coeff_ps * [velocity of primary dof] + loading%offset_ps

    type(mesh_type), intent(in) :: mesh
    type(loading_options_type), intent(inout) :: loading_options
    type(loading_type), intent(inout) :: loading

    ! Local:
    integer :: i, j, ii, jj   ! loop variables
    integer :: nb_tot_dof, node_dof

    ! Allocate dof arrays to mpc
    allocate (loading%ps_dof(dof_sub:dof_sup))
    allocate (loading%coeff_ps(dof_sub:dof_sup))
    allocate (loading%offset_ps(dof_sub:dof_sup))

    ! Determine the number of total dofs (probably redundant with an other variable in the code)
    nb_tot_dof = mesh%num_nodes * 3

    ! Initialisation of variables
    do i = dof_sub, dof_sup
      loading%ps_dof(i) = i         ! initialisation of array to node dof index
    end do

    loading%coeff_ps = 1.0d0
    loading%offset_ps = 0.0d0
    loading%nb_secondary_dof = 0      ! at initiation, all dofs are considered as primary

    ! Loop on mpc defined in .cfg
    do i = 1, loading_options%num_mpcs
      ! Loop on nsets defined in .msh
      do j = 1, size(mesh%nsets)
        ! If a nset has a defined mpc
        if (loading_options%general_mpc(i)%nset .eq. mesh%nsets(j)%nset_label) then
          ! Parsing all the components of the mpc array
          do ii = 1, loading_options%general_mpc(i)%nbdof
            ! Define the direction where mpc are applied
            select case (loading_options%general_mpc(i)%dir(ii))

            case ('x', 'vx')
              node_dof = 1
            case ('y', 'vy')
              node_dof = 2
            case ('z', 'vz')
              node_dof = 3

            end select

            ! Affect the mpc relations between dof based on nodes of nset:
            !     - first one (the primary) does not change the value in loading%ps_dof
            !     - the others (secondaries) take the opposite value of the primary dof node id
            ! Set of the corresponding linear relation between primary and secondary dof
            ! velocity.
            do jj = 2, mesh%nsets(j)%num_nset_nodes
              ! if in dof_sub, dof_sup range
              if ((((mesh%nsets(j)%nset_nodes(jj)-1)*3 + node_dof) .ge. dof_sub) .and. &
                  &(((mesh%nsets(j)%nset_nodes(jj)-1)*3 + node_dof) .le. dof_sup)) then

                loading%ps_dof((mesh%nsets(j)%nset_nodes(jj) - 1)*3 + node_dof) &
                  &= ((mesh%nsets(j)%nset_nodes(1) - 1)*3 + node_dof)
                loading%coeff_ps((mesh%nsets(j)%nset_nodes(jj) - 1)*3 + node_dof) &
                  &= loading_options%general_mpc(i)%coeff(ii)
                loading%offset_ps((mesh%nsets(j)%nset_nodes(jj) - 1)*3 + node_dof) &
                  &= loading_options%general_mpc(i)%offset(ii)

                ! increment the nb of secondary dofs
                loading%nb_secondary_dof = loading%nb_secondary_dof + 1

                if (loading_options%general_mpc(i)%coeff(ii) .eq. 0) then
                  call par_quit('Error  :     > Multi-point constraints &
                    &must have nonzero coefficient.')
                end if
              end if
            end do
          end do
        end if
      end do
    end do

    ! loop to go to the top primary dof od each dof (under the form of a mpc connectivity matrix)
        ! Determine the number or primary dofs
    loading%nb_primary_dof = nb_tot_dof - loading%nb_secondary_dof


  end subroutine read_mpcs

  !> launch_mpcs: set mpcs
  subroutine launch_mpcs(mesh, exec, loading)

    type(mesh_type), intent(inout) :: mesh
    type(exec_type), intent(inout) :: exec
    type(loading_type), intent(inout) :: loading

    ! Calculate multi-point constraints.
  if (loading%mpc_status .eqv. .true.) then

    if (myid .eq. 0) then
      write (*, '(a)') "Info   :   Set multi-point constraints"
    end if

    call calc_mpcs(mesh, exec, loading)
    call propagate_bcs_mpcs(mesh, exec, loading)

  end if

  end subroutine launch_mpcs

 end module multi_point_constraints_mod
