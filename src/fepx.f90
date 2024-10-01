! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

program FEPX

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use types_mod
  use printing_type_mod
  use exec_type_mod
  use read_input_mod
  use multi_point_constraints_mod
  use driver_triaxcsr_mod
  use driver_triaxclr_mod
  use driver_uniaxial_control_mod
  use units_mod
  use parallel_mod
  use gather_scatter_mod
  use res_init_mod
  use fepx_mod

  implicit none

  type(mesh_type) :: mesh
  type(crys_type), allocatable :: crys (:)
  type(loading_type) :: loading
  type(exec_type) :: exec
  type(results_type) :: results
  type(printing_type) :: printing

  ! Managing command-line arguments ----------------------------------------------

  if (command_argument_count() .ge. 2) then
    call fepx_runonarguments
  end if

  ! Initializing -----------------------------------------------------------------

  call fepx_init(exec)

  call fepx_header(num_procs, myid)

  ! call mesh_set_default(mesh)
  ! call loading_set_default(mesh)
  call exec_set_default (exec)
  call printing_set_default(printing)

  ! Reading and processing input -----------------------------------------------

  if (myid .eq. 0) then
    write (*, '(a)') 'Info   : Loading simulation...'
  end if

  call read_input (mesh, crys, loading, printing, exec)

  ! Allocating result variables --------------------------------------------------

  call res_init (mesh, crys, loading, printing, results)

  ! Beginning the deformation simulation proper ----------------------------------
  call open_output_files(mesh, printing, myid)

  call fepx_print_partinfo (printing, mesh, num_procs, myid)

  if (loading%def_control_by(1:8) .eq. 'uniaxial') then
    call driver_uniaxial_control(mesh, crys, loading, exec, results, printing)

  else if (loading%def_control_by .eq. 'triaxial_constant_strain_rate') then
    call driver_triaxcsr(mesh, crys, loading, exec, results, printing)

  else if (loading%def_control_by .eq. 'triaxial_constant_load_rate') then
    call driver_triaxclr(mesh, crys, loading, exec, results, printing)

  else
    write (*,*) loading%def_control_by
    call par_quit('Error  :     > Invalid deformation control option.')

  end if

  ! Exiting ----------------------------------------------------------------------

  call par_quit('Info   : Completed simulation successfully',clock_start=exec%clock_start)

end program FEPX
