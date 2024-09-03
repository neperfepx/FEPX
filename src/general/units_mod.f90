! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module units_mod

! Module to handle fortran unit numbers for printing.

! Contains subroutines:
! open_log_files: Opens log files.
! open_output_files: Opens output files.

! From libfepx:

  use general_mod
  use types_mod
  use parallel_mod
  use utils_mod

  implicit none

!  Output units.


contains

  subroutine open_output_files(mesh, printing, myid)

    ! Open files for output.  Append process identifier string to file name.

    !---------------------------------------------------------------------------

    ! Arugments
    ! printing: Type containing print options
    ! myid: Processor number
    type(mesh_type), intent(in) :: mesh
    type(printing_type), intent(in) :: printing
    integer, intent(in) :: myid
    ! Locals:
    integer :: iostatus
    character(len=10000):: filename,dir_name,cwd
    character(len=8) :: charid ! assumes less than 10,000 processes
    character(len=256) :: message
    logical :: file_exists
    integer :: rst_num, i
    character(len=8) :: rst_num_str

    !---------------------------------------------------------------------------

    ! Print 1-indexed process numbers
    write (charid, '(i0)') myid + 1

    ! If restarting simulation, get new restart cycle number
    ! Find max value n of rstN.control
    if (printing%restart) then
      rst_num = 1000
      file_exists = .false.
      do while (.not. file_exists)
        write (rst_num_str, '(i0)') rst_num
        filename = 'rst'//trim(rst_num_str)//'.control'
        inquire (file=filename, exist=file_exists)
        rst_num = rst_num - 1
        if (rst_num .eq. -2) then
          call par_quit('Error  :     > Restart control file not found.')
        end if
      end do
      ! Increase n to write new files for next simulation cycle
      rst_num = rst_num + 2
      write (rst_num_str, '(i0)') rst_num
    end if

    ! Only for main process
    if (myid .eq. 0) then
      !  .sim
      call getcwd(cwd)
      dir_name = trim(cwd)//"/simulation.sim"
      call system('mkdir '//trim(dir_name))
      call system('cd '//trim(dir_name))
      call system('mkdir '//trim(dir_name)//'/inputs')
      call system('cp simulation.tess simulation.cfg simulation.msh simulation.ori&
                   & simulation.phase simulation.opt simulation.tesr '//trim(dir_name)//'/inputs 2> /dev/null')
      call system('cp *.sh '//trim(dir_name)//'/inputs'//' 2> /dev/null')
      call system('mkdir '//trim(dir_name)//'/results')
      call system('mkdir '//trim(dir_name)//'/results/nodes')
      call system('mkdir '//trim(dir_name)//'/results/elts')
      if (printing%print_coo) then
        call system('mkdir '//trim(dir_name)//'/results/nodes/coo')
      end if
      !
      if (printing%print_disp) then
        call system('mkdir '//trim(dir_name)//'/results/nodes/disp')
      end if
      !
      if (printing%print_crss) then
        call system('mkdir '//trim(dir_name)//'/results/elts/crss')
      end if
      !
      if (printing%print_defrate_eq) then
        call system('mkdir '//trim(dir_name)//'/results/elts/defrate_eq')
      end if
      !
      if (printing%print_defrate_pl_eq) then
        call system('mkdir '//trim(dir_name)//'/results/elts/defrate_pl_eq')
      end if
      !
      if (printing%print_defrate_pl) then
        call system('mkdir '//trim(dir_name)//'/results/elts/defrate_pl')
      end if
      !
      if (printing%print_strain_pl_eq) then
        call system('mkdir '//trim(dir_name)//'/results/elts/strain_pl_eq')
      end if
      !
      if (printing%print_strain_el_eq) then
        call system('mkdir '//trim(dir_name)//'/results/elts/strain_el_eq')
      end if
      !
      if (printing%print_strain_eq) then
        call system('mkdir '//trim(dir_name)//'/results/elts/strain_eq')
      end if
      !
      if (printing%print_stress_eq) then
        call system('mkdir '//trim(dir_name)//'/results/elts/stress_eq')
      end if
      !
      if (printing%print_slip) then
        call system('mkdir '//trim(dir_name)//'/results/elts/slip')
      end if
      !
      if (printing%print_sliprate) then
        call system('mkdir '//trim(dir_name)//'/results/elts/sliprate')
      end if
      !
      if (printing%print_rss) then
        call system('mkdir '//trim(dir_name)//'/results/elts/rss')
      end if
      !
      if (printing%print_ori) then
        call system('mkdir '//trim(dir_name)//'/results/elts/ori')
      end if
      !
      if (printing%print_strain) then
        call system('mkdir '//trim(dir_name)//'/results/elts/strain')
      end if
      !
      if (printing%print_strain_el) then
        call system('mkdir '//trim(dir_name)//'/results/elts/strain_el')
      end if
      !
      if (printing%print_strain_pl) then
        call system('mkdir '//trim(dir_name)//'/results/elts/strain_pl')
      end if
      !
      if (printing%print_stress) then
        call system('mkdir '//trim(dir_name)//'/results/elts/stress')
      end if
      !
      if (printing%print_vel) then
        call system('mkdir '//trim(dir_name)//'/results/nodes/vel')
      end if
      !
      if (printing%print_velgrad) then
        call system('mkdir '//trim(dir_name)//'/results/elts/velgrad')
      end if
      !
      if (printing%print_spinrate) then
        call system('mkdir '//trim(dir_name)//'/results/elts/spinrate')
      end if
      !
      if (printing%print_rotrate) then
        call system('mkdir '//trim(dir_name)//'/results/elts/rotrate')
      end if
      !
      if (printing%print_rotrate_spin) then
        call system('mkdir '//trim(dir_name)//'/results/elts/rotrate_spin')
      end if
      !
      if (printing%print_rotrate_slip) then
        call system('mkdir '//trim(dir_name)//'/results/elts/rotrate_slip')
      end if
      !
      if (printing%print_work) then
        call system('mkdir '//trim(dir_name)//'/results/elts/work')
      end if
      !
      if (printing%print_work_pl) then
        call system('mkdir '//trim(dir_name)//'/results/elts/work_pl')
      end if
      !
      if (printing%print_defrate) then
        call system('mkdir '//trim(dir_name)//'/results/elts/defrate')
      end if
      !
      if (printing%print_workrate) then
        call system('mkdir '//trim(dir_name)//'/results/elts/workrate')
      end if
      !
      if (printing%print_workrate_pl) then
        call system('mkdir '//trim(dir_name)//'/results/elts/workrate_pl')
      end if
      !
      ! post.force.#
      if (printing%print_forces) then          
        call system('mkdir '//trim(dir_name)//'/results/forces')
        select case (printing%restart_file_handling)
        case ("normal") ! Normal simulation, no restart
          do i = 1, mesh%num_fasets
            ! Loop over all force units
            filename = trim(dir_name)//'/results/forces/'//mesh%faset_labels(i)
            open (printing%force_u + i - 1, &
              & file=filename, iostat=iostatus,status='new')
            if (iostatus .ne. 0) then
              write (message, '(a,a,a)') 'Error  :     > io &
                &Failure to open ', trim(filename), ' file'
              call par_quit(trim(adjustl(message)))
            end if
          end do
        case ("restart") ! Restart simulation, new file
          do i = 1, mesh%num_fasets
            filename = 'post.force.rst'//trim(rst_num_str)//'.'//mesh%faset_labels(i)
            open (printing%force_u + i - 1, &
              & file=filename, iostat=iostatus,status='new')
            if (iostatus .ne. 0) then
              write (message, '(a,a,a)') 'Error  :     > io &
                &Failure to open ', trim(filename), ' file.'
              call par_quit(trim(adjustl(message)))
            end if
          end do
        case default ! unsupported value
          call par_quit('Error  :     > Invalid restart file handling&
            & option provided.')
        end select
      end if

      ! post.conv
      if (printing%print_conv) then
        select case (printing%restart_file_handling)
        case ("normal") ! Normal simulation, no restart
          filename = 'post.conv'
          open (printing%conv_u, file=filename, iostat=iostatus)
          if (iostatus .ne. 0) then
            write (message, '(a,a,a)') 'Error  :     > io &
              &Failure to open ', trim(filename), ' file.'
            call par_quit(trim(adjustl(message)))
          end if
        case ("restart") ! Restart simulation, new file
          filename = 'post.conv.rst'//trim(rst_num_str)
          open (printing%conv_u, file=filename, iostat=iostatus)
          if (iostatus .ne. 0) then
            write (message, '(a,a,a)') 'Error  :     > io &
              &Failure to open ', trim(filename), ' file.'
            call par_quit(trim(adjustl(message)))
          end if
        case default ! unsupported value
          call par_quit('Error  :     > Invalid restart file handling&
            & option provided.')
        end select
      end if
    end if

  end subroutine open_output_files

end module units_mod
