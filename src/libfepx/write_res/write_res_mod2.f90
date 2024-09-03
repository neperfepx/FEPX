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
    write (ounit, fmt) (anint (var(i) / output_precision) * output_precision, i=start, end)

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
    write (ounit, fmt) (anint (var(:, i) / output_precision) * output_precision, i=start, end)

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
      write (ounit, fmt) ((/anint (var(1, 1, i) / output_precision) * output_precision, &
                            anint (var(2, 2, i) / output_precision) * output_precision, &
                            anint (var(3, 3, i) / output_precision) * output_precision, &
                            anint (var(2, 3, i) / output_precision) * output_precision, &
                            anint (var(1, 3, i) / output_precision) * output_precision, &
                            anint (var(1, 2, i) / output_precision) * output_precision/), i=start, end)
    else if (numdim .eq. 9) then
      write (ounit, fmt) ((/anint (var(1, 1, i) / output_precision) * output_precision, &
                            anint (var(1, 2, i) / output_precision) * output_precision, &
                            anint (var(1, 3, i) / output_precision) * output_precision, &
                            anint (var(2, 1, i) / output_precision) * output_precision, &
                            anint (var(2, 2, i) / output_precision) * output_precision, &
                            anint (var(2, 3, i) / output_precision) * output_precision, &
                            anint (var(3, 1, i) / output_precision) * output_precision, &
                            anint (var(3, 2, i) / output_precision) * output_precision, &
                            anint (var(3, 3, i) / output_precision) * output_precision/), i=start, end)
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

  ! ============================================================================

  !> write_dot_sim_file_output_files: Writes the files requested to be printed
  subroutine write_dot_sim_file_output_files(printing, curr_step, &
    & input_string, pflag)

    ! This writes the file names of the user-defined output files from the
    ! simulation.cfg.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! input_string: Character array denoting file name to print.
    ! curr_step: Current step at which we are printing information for.
    ! pflag: Determine if printing should happen or not
    type(printing_type), intent(in) :: printing
    character(len=*), intent(in) :: input_string(:)
    integer, intent(in) :: curr_step
    logical, intent(in) :: pflag
    ! Locals:
    ! i: Looping variable for printing of input_string
    integer :: i

    !---------------------------------------------------------------------------
    ! Only print to the .sim file from the root processor
    if (myid .eq. 0) then
      ! If the first step is being printed then print to the .sim file.
      ! If the "first" step after restarting a simulation, also print.
      if ((pflag .eqv. .true.) .or. ((printing%restart .eqv. .true.) &
        & .and. (curr_step .eq. printing%restart_initial_step))) then
        write (printing%dot_sim_u, '(a)') trim(input_string(1))
        write (printing%dot_sim_u, '(a)') '  *result'
        write (printing%dot_sim_u, '(a,i0)') '   ', size(input_string)-1
        write (printing%dot_sim_u, '(25(a,1x))') '  ',&
          & (trim(input_string(i)), i=2, size(input_string))
      end if
    end if

  end subroutine write_dot_sim_file_output_files

  !===========================================================================

  subroutine add_to_output_files_list(printing, curr_step, &
    & input_string, string_array, pflag)

    ! This writes the file names of the user-defined output files from the
    ! simulation.cfg.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! input_string: Character array denoting file name to print.
    ! curr_step: Current step at which we are printing information for.

    type(printing_type), intent(in) :: printing
    character(len=*), intent(in) :: input_string
    character(len=16), intent(inout), allocatable :: string_array(:)
    integer, intent(in) :: curr_step
    logical, intent(in) :: pflag

    character(len=16), allocatable :: temp_array(:)
    integer :: arraysize, i

    !---------------------------------------------------------------------------

    ! Only print to the .sim file from the root processor
    if (myid .eq. 0) then
      ! If the first step is being printed then print to the .sim file.
      ! If the "first" step after restarting a simulation, also print.
      if ((pflag .eqv. .true.) .or. ((printing%restart .eqv. .true.) &
          & .and. (curr_step .eq. printing%restart_initial_step))) then
        ! Append the input_string to the present array and reallocate.
        if (allocated(string_array)) then
          arraysize = size(string_array)
          allocate (temp_array(arraysize + 1))

          do i = 1, arraysize
            temp_array(i) = string_array(i)
          end do

          temp_array(arraysize + 1) = input_string
          deallocate (string_array)
          call move_alloc(temp_array, string_array)
        end if
      end if
    end if

  end subroutine add_to_output_files_list

end module write_res_mod2
