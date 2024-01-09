! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
module read_input_mod

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use read_input_cfg_mod
  use read_input_msh_mod
  use units_mod
  use crys_type_mod
  use loading_options_type_mod
  use boundary_conditions_mod
  use multi_point_constraints_mod
  use multi_point_constraints_mod2
  use periodicity_mod1
  use parallel_mod
  use gather_scatter_mod

  implicit none

  public

contains

  subroutine read_input (mesh, crys, loading, printing, exec)

    type(mesh_type), intent(inout) :: mesh
    type(crys_type), allocatable, intent(inout) :: crys (:)
    type(loading_type), intent(inout) :: loading
    type(printing_type), intent(inout) :: printing
    type(exec_type), intent(inout) :: exec

    integer :: i
    character(len=20) :: mesh_file = "", ori_file = "", phase_file = ""
    type(loading_options_type) :: loading_options

  call loading_options_set_default(loading_options)

  call read_config(mesh_file, ori_file, phase_file, printing, loading_options, &
                  & crys, exec)

  call crys_initparams(crys)

  ! read mesh and initialize partition sizes -----------------------------------
  call read_meshsize(mesh_file, mesh)

  call par_partition(mesh%num_elts, num_procs, myid, elt_sub, elt_sup)
  call par_partition(mesh%num_nodes, num_procs, myid, node_sub, node_sup)
  dof_sub = 3*(node_sub - 1) + 1
  dof_sup = 3*node_sup

  call read_mesh(mesh_file, ori_file, phase_file, mesh)
  !  ---------------------------------------------------------------------------

  if (mesh%num_phases .ne. size(crys)) then
    call par_quit('Error  :     > Number of phases in mesh and config file&
        & do not match.')
      end if

  mesh%maxnumslip = 0
  do i = 1, size(crys)
    mesh%maxnumslip = max (mesh%maxnumslip, crys(i)%numslip)
  end do

  if (myid .eq. 0) then
    write (*, '(a)') 'Info   : Initializing simulation...'
  end if

  ! compute boundary conditions ------------------------------------------------
  call calc_bcs(mesh, loading_options, loading)
  !  ---------------------------------------------------------------------------

  ! compute mpcs, if necessary -------------------------------------------------
  ! If mpcs defined in .cfg, read them

  if (loading%mpc_status .eqv. .true.) then
    call read_mpcs(mesh, loading_options, loading)
  end if

  ! Setting parallel (Mika's MPI gather/scatter routines)
  call part_scatter_setup(1, kdim, dof_sub, dof_sup, elt_sub, elt_sup, &
      & mesh%elt_dofs, exec%dof_trace)
  call part_scatter_setup(1, ndim, node_sub, node_sup, elt_sub, elt_sup, &
      & mesh%elt_nodes, exec%node_trace)

  ! Setting multi-point constraints
  if (loading%mpc_status .eqv. .true.) then
    call launch_mpcs(mesh, exec, loading)
  end if

  ! Loading of periodic relations.
  ! Note: periodic relation erase other MPCs if specified on same nodes
  if (mesh%num_periodicity .ne. 0) then
     call set_periodicity(mesh, loading, loading_options, exec)
  end if

end subroutine read_input

end module read_input_mod
