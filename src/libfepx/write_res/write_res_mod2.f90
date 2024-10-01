! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module write_res_mod2

! Module for printing variables to file at the end of a step.

! Contains subroutines:

! General printing and handling of field variable data:
! write_res: Print variables on a given load step

! Handling the .sim, post.forceX, and post.conv files
! write_force_file_headers: Writes the file headers for table formatting
! write_force_file_data: Writes the surface forces [x y z] for all surfaces
! write_conv_file_headers: Writes the file headers for the table formatting
! write_conv_file_data: Writes the various convergence statistics

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use matrix_operations_mod
  use orientation_conversion_mod
  use units_mod
  use parallel_mod
  use fepx_config_mod

  implicit none

! Public

  public

contains

  subroutine write_res_vector(var, ounit, numdim, step, start, end)

    real(rk), intent(in) :: var(start:end)
    integer, intent(in) :: ounit
    integer, intent(in) :: numdim
    integer, intent(in) :: step
    integer, intent(in) :: start
    integer, intent(in) :: end
    ! Locals:
    integer :: i
    character(len=50) :: fmt

    !---------------------------------------------------------------------------

    ! Assemble format string
    write (fmt, *) '(', numdim, '(e13.7,1x))'

    ! Write to file
    write (ounit, fmt) (anint(var(i)/output_precision)*output_precision, i=start, end)

  end subroutine write_res_vector

  !=============================================================================

  subroutine write_res_array(var, ounit, numdim, step, start, end)

    real(rk), intent(in) :: var(numdim, start:end)
    integer, intent(in) :: ounit
    integer, intent(in) :: numdim
    integer, intent(in) :: step
    integer, intent(in) :: start
    integer, intent(in) :: end
    ! Locals:
    integer :: i
    character(len=50) :: fmt

    !---------------------------------------------------------------------------

    ! Assemble format string
    write (fmt, *) '(', numdim, '(e13.7,1x))'

    ! Write to file
    write (ounit, fmt) (anint(var(:, i)/output_precision)*output_precision, i=start, end)

  end subroutine write_res_array

  !=============================================================================

  subroutine write_res_tensor(var, ounit, numdim, step, start, end)

    real(rk), intent(in) :: var(3, 3, start:end)
    integer, intent(in) :: ounit
    integer, intent(in) :: numdim
    integer, intent(in) :: step
    integer, intent(in) :: start
    integer, intent(in) :: end
    ! Locals:
    integer :: i
    character(len=50) :: fmt

    !---------------------------------------------------------------------------

    ! Assemble format string
    write (fmt, *) '(', numdim, '(e13.7,1x))'

    ! Write to file
    if (numdim .eq. 6) then
      write (ounit, fmt) ((/anint(var(1, 1, i)/output_precision)*output_precision, &
                            anint(var(2, 2, i)/output_precision)*output_precision, &
                            anint(var(3, 3, i)/output_precision)*output_precision, &
                            anint(var(2, 3, i)/output_precision)*output_precision, &
                            anint(var(1, 3, i)/output_precision)*output_precision, &
                            anint(var(1, 2, i)/output_precision)*output_precision/), i=start, end)
    else if (numdim .eq. 9) then
      write (ounit, fmt) ((/anint(var(1, 1, i)/output_precision)*output_precision, &
                            anint(var(1, 2, i)/output_precision)*output_precision, &
                            anint(var(1, 3, i)/output_precision)*output_precision, &
                            anint(var(2, 1, i)/output_precision)*output_precision, &
                            anint(var(2, 2, i)/output_precision)*output_precision, &
                            anint(var(2, 3, i)/output_precision)*output_precision, &
                            anint(var(3, 1, i)/output_precision)*output_precision, &
                            anint(var(3, 2, i)/output_precision)*output_precision, &
                            anint(var(3, 3, i)/output_precision)*output_precision/), i=start, end)
    end if

  end subroutine write_res_tensor

  !=============================================================================

  subroutine write_res_zeros(ounit, numdim, start, end)

    integer, intent(in) :: ounit
    integer, intent(in) :: numdim
    integer, intent(in) :: start
    integer, intent(in) :: end
    ! Locals:
    integer :: i, var(numdim, start:end)
    character(len=50) :: fmt

    !---------------------------------------------------------------------------
    var = 0
    ! Assemble format string
    write (fmt, *) '(', numdim, '(i0,1x))'

    ! Write to file
    write (ounit, fmt) (var(:, i), i=start, end)

  end subroutine write_res_zeros

  !===========================================================================

  subroutine write_dot_sim_file_header(printing, mesh)

    ! Write a hidden .sim file containing required information to facilitate
    ! automated post-processing via Neper.

    ! This prints the number of nodes, elements, partitions, elements by part,
    ! nodes by part, number of slip systems, and orientation definition.

    ! Arugments:
    ! part_array: Gathered array from fepx.f90 that contains partition info.

    type(printing_type), intent(in) :: printing
    type(mesh_type), intent(in) :: mesh

    write (printing%dot_sim_u, '(a)') '***sim'
    write (printing%dot_sim_u, '(a)') ' **format'
    write (printing%dot_sim_u, '(a)') '   1.1'
    write (printing%dot_sim_u, '(a)') ' **input'
    if (ut_file_exists("simulation.tess")) then
      write (printing%dot_sim_u, '(a)') '  *tess'
      write (printing%dot_sim_u, '(a)') '   simulation.tess'
    end if
    if (ut_file_exists("simulation.cfg")) then
      write (printing%dot_sim_u, '(a)') '  *cfg'
      write (printing%dot_sim_u, '(a)') '   simulation.cfg'
    end if
    if (ut_file_exists("simulation.msh")) then
      write (printing%dot_sim_u, '(a)') '  *msh'
      write (printing%dot_sim_u, '(a)') '   simulation.msh'
    end if
    if (ut_file_exists("simulation.ori")) then
      write (printing%dot_sim_u, '(a)') '  *ori'
      write (printing%dot_sim_u, '(a)') '   simulation.ori'
    end if
    if (ut_file_exists("simulation.phase")) then
      write (printing%dot_sim_u, '(a)') '  *phase'
      write (printing%dot_sim_u, '(a)') '   simulation.phase'
    end if
    if (ut_file_exists("simulation.opt")) then
      write (printing%dot_sim_u, '(a)') '  *opt'
      write (printing%dot_sim_u, '(a)') '   simulation.opt'
    end if
    if (ut_file_exists("simulation.tesr")) then
      write (printing%dot_sim_u, '(a)') '  *tesr'
      write (printing%dot_sim_u, '(a)') '   simulation.tesr'
    end if
    write (printing%dot_sim_u, '(a)') ' **general'
    write (printing%dot_sim_u, '(a,5(i0,1x))') '   ', mesh%num_cell, mesh%num_nodes,&
        & mesh%num_elts, mesh%num_elsets, num_procs
    write (printing%dot_sim_u, '(a)') '  *orides'
    write (printing%dot_sim_u, '(a,a,a,a)') '   ', trim(mesh%orientation_parameterization), ':', trim(mesh%orientation_convention)

  end subroutine write_dot_sim_file_header

  ! ============================================================================

  !> write_dot_sim_file_output_files: Writes the files requested to be printed
  subroutine write_dot_sim_file_output_files(printing)

    type(printing_type), intent(in) :: printing

    integer :: i

    write (printing%dot_sim_u, '(a)') '**entity node'
    write (printing%dot_sim_u, '(a)') '  *result'
    write (printing%dot_sim_u, '(a,i0)') '   ', size(printing%node_results)
    write (printing%dot_sim_u, '(25(a,1x))') '  ', (trim(printing%node_results(i)), i=1, size(printing%node_results))

    write (printing%dot_sim_u, '(a)') '**entity elt'
    write (printing%dot_sim_u, '(a)') '  *result'
    write (printing%dot_sim_u, '(a,i0)') '   ', size(printing%elt_results)
    write (printing%dot_sim_u, '(25(a,1x))') '  ', (trim(printing%elt_results(i)), i=1, size(printing%elt_results))

  end subroutine write_dot_sim_file_output_files

  !===========================================================================

  subroutine write_dot_sim_file_complete_steps(printing, step)
    !
    ! This writes the completed number of steps from the driver as the current
    ! simulation is either finishing successfully or terminated early.

    type(printing_type), intent(in) :: printing
    integer, intent(in) :: step

    ! Print the completed number of steps

    write (printing%dot_sim_u, '(a)') ' **step'
    write (printing%dot_sim_u, '(a,i0)') '   ', step

    write (printing%dot_sim_u, '(a)') '***end'
    ! Close the .sim file before ending the process

    close (printing%dot_sim_u)

  end subroutine write_dot_sim_file_complete_steps

  !===========================================================================

end module write_res_mod2
