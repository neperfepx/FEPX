! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
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

  subroutine fepx_init(exec)

    type(exec_type), intent(inout) :: exec

    call cpu_time(exec%clock_start)

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
      write (*, '(a)') 'Info   : Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.'
      write (*, '(a)') 'Info   : &
          &---------------------------------------------------------------'
      write (*, '(a,i0,a,i0,a,i0,a,i0,a,i2.2)') 'Info   : Start time: ',&
          & timevalues(1), '-', timevalues(2), '-', timevalues(3), ' at ',&
          & timevalues(5), ':', timevalues(6)
    end if

  end subroutine fepx_header

! Gather partition information for the .sim file.
  subroutine fepx_print_partinfo(printing, mesh, num_procs, myid)

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
    part_info(3, myid + 1) = (node_sub - 1)*3
    part_info(4, myid + 1) = (elt_sub - 1)

    call par_gather(part_info(:, myid + 1), global_info, 4)

    mesh%global_info(:, :) = global_info(:, :)

    if (myid .eq. 0) then
      call write_dot_sim_file_header(printing, mesh)
    end if

    deallocate (part_info)
    deallocate (global_info)

  end subroutine fepx_print_partinfo

  subroutine fepx_runonarguments

    ! Variables
    integer :: i, diff
    real(rk) :: eps, reps
    character(len=1000) :: option, file1, file2, eps_string

    ! Default epsilon value
    eps = 1.0e-6
    reps = 1.0e-4

    call get_command_argument(1, option) ! Get the option

    if (option .eq. "--diff") then
      if (command_argument_count() .lt. 3) then
        stop - 1
      else
        call get_command_argument(2, file1)
        call get_command_argument(3, file2)

        if (command_argument_count() .ge. 4) then
          call get_command_argument(4, eps_string)
          read (eps_string, *) eps
        end if
        if (command_argument_count() .ge. 5) then
          call get_command_argument(5, eps_string)
          read (eps_string, *) reps
        end if
      end if

      call fepx_diff(trim(file1), trim(file2), eps, reps, diff)

      stop diff
    end if

  end subroutine fepx_runonarguments

  subroutine fepx_diff(file1, file2, eps, reps, diff)
    implicit none

    character(len=*), intent(in) :: file1, file2
    real(8), intent(in) :: eps, reps
    integer, intent(out) :: diff

    character(len=1000) :: line1, line2
    integer :: ios1, ios2, status
    real(rk) :: val1, val2

    ! Initialize the diff flag to 0 (no difference)
    diff = 0

    ! Open the files
    open (unit=10, file=trim(file1), status='old', action='read', iostat=ios1)
    if (ios1 .ne. 0) then
      print *, 'Error opening file: ', trim(file1)
      stop - 1
    end if

    open (unit=20, file=trim(file2), status='old', action='read', iostat=ios2)
    if (ios2 .ne. 0) then
      print *, 'Error opening file: ', trim(file2)
      stop - 1
    end if

    ! Compare files line by line
    do
      read (10, '(A)', iostat=ios1) line1
      read (20, '(A)', iostat=ios2) line2

      ! Check for end-of-file
      if (ios1 .ne. 0 .or. ios2 .ne. 0) exit

      ! Try to convert lines to floating-point numbers
      read (line1, *, iostat=status) val1
      if (status .eq. 0) then
        read (line2, *, iostat=status) val2
        if (status .eq. 0) then
          if (abs(val1 - val2) > eps .and. abs(val1 - val2) > reps * max(abs(val1), abs(val2))) then
            diff = 1
            exit
          end if
        else
          diff = 1
          exit
        end if
      else
        ! If not numbers, compare strings
        if (trim(line1) .ne. trim(line2)) then
          diff = 1
          exit
        end if
      end if
    end do

    ! Close the files
    close (10)
    close (20)

  end subroutine fepx_diff

end module fepx_mod
