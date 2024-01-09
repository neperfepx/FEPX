! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.

module fepx_mod

! Full-field 3d polycrystal fem analysis with anisotropic elasticity and
! viscoplasticity. Primary arrays are automatically allocated here to be
! used in the subsequent driver routines which perform the simulation.

!-------------------------------------------------------------------------------

  use parallel_mod
  use gather_scatter_mod
  use write_res_mod
  use fepx_config_mod
  use surface_mod

  implicit none

  public

  contains

subroutine fepx_init (exec)

  type(exec_type), intent(inout) :: exec

  call cpu_time (exec%clock_start)

  call par_init
  call quadrature_init
  call surf_init

end subroutine fepx_init

subroutine fepx_header(num_procs, myid)

  integer, intent(in) :: num_procs, myid
  integer :: timevalues(8)

  if (myid .eq. 0) then
    call date_and_time(values=timevalues)
    write (*, '(a)') '==========================    &
        & F   E   P   X   =========================='
    write (*, '(a)') 'Info   : A finite element software package for &
        & polycrystal plasticity.'
    write (*, '(a,a)') 'Info   : Version ', version
    write (*, '(a,i0,a)') 'Info   : Running on ', num_procs, ' cores.'
    write (*, '(a)') 'Info   : <https://fepx.info>'
    write (*, '(a)') 'Info   : Copyright (C) 1996-2023, DPLab, ACME Lab.'
    write (*, '(a)') 'Info   : &
        &---------------------------------------------------------------'
    write (*, '(a,i0,a,i0,a,i0,a,i0,a,i2.2)') 'Info   : Start time: ',&
        & timevalues(1), '-', timevalues(2), '-', timevalues(3), ' at ',&
        & timevalues(5), ':', timevalues(6)
  end if

end subroutine fepx_header

! Gather partition information for the .sim file.
subroutine fepx_print_partinfo (printing, mesh, num_procs, myid)

  type(printing_type), intent(in) :: printing
  type(mesh_type), intent(inout) :: mesh
  integer, intent(in) :: num_procs, myid

  integer, allocatable :: part_info(:, :)
  integer, allocatable :: global_info(:, :)
  integer :: num_elt_part, num_node_part

  allocate (part_info(4, num_procs))
  allocate (global_info(4, num_procs))
  allocate (mesh%global_info(4, num_procs))
  part_info = 0
  global_info = 0

  num_elt_part = (elt_sup - elt_sub) + 1
  num_node_part = (node_sup - node_sub) + 1
  part_info(1, myid + 1) = num_elt_part
  part_info(2, myid + 1) = num_node_part
  part_info(3, myid + 1) = (node_sub-1)*3
  part_info(4, myid + 1) = (elt_sub-1)

  call par_gather(part_info(:, myid + 1), global_info, 4)

  mesh%global_info(:,:) = global_info(:,:)
  
  if (myid .eq. 0) then
    call write_dot_sim_file_header(printing, mesh)
  end if

  deallocate (part_info)
  deallocate (global_info)

end subroutine fepx_print_partinfo

end module fepx_mod
