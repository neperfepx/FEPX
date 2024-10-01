! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module read_input_msh_mod

! Module to handle all mesh parsing.

! Contains subroutines:
! allocate_msh: Allocates mesh arrays based on the problem size.
! get_msh_size: Scrapes the msh file and retrieve the number of nodes and elems.
! read_spatial_msh: Parent subroutine to handle msh file parsing process.

! Helper subroutines:
! For parsing the mesh (and related external) file:
! get_node_info: Gets the number of nodes in the mesh.
! read_mesh_format: Reads in gmsh file format - 2.2 0 8 only.
! read_mesh_version: Reads in gmsh file version.
! read_nodes: Read in node id and [x y z] spatial coordinates.
! read_elts: Read in all elements and only store 3d elements.
! read_nsets: Read in node sets - currently not stored for read-in BCs.
! read_fasets: Read in surface face node sets (2d elements).
! read_nodepartitions: Read in per-node partition distribution. -optional
! read_physicalnames: Read in physical names and do not store (gmsh data only).
! read_elsetorientations: Read in elset ori.
! read_eltorientations: Read in element ori.
! read_groups: Read in elset phases (assumes single-phase by default). -optional

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use gather_scatter_mod
  use parallel_mod
  use read_input_msh_mod2

  implicit none

  public

contains

  !===========================================================================

  !> @brief Extract the number of elements and nodes from the mesh
  !! @param[in]   mesh_file    file name
  subroutine read_meshsize (mesh_file, mesh)

    ! Scrapes $Nodes and $Elements fields for the number of nodes and 3d elems
    ! to determine problem size. Also, attempt to retrieve embedded mesh
    ! partitioning info from $Elements and $NodePartitions (if available).

    character(len=*), intent(in) :: mesh_file
    type(mesh_type), intent(in) :: mesh
    integer :: iostatus
    integer :: eofstat
    integer :: file_id
    character(len=256) :: line

    open (newunit=file_id, file=mesh_file, status='old', action='read', &
        & iostat=iostatus)

    if (iostatus .ne. 0) then
      call par_quit("Error  :     > Failure to open `simulation.msh'.")
    end if

    eofstat = 0

    do while (eofstat .eq. 0)
      ! Read in line and determine which section is to be parsed.
      ! The below line is a temporary fix to avoid program hangs - jc
      line = ""

      read (file_id, '(a)', iostat=eofstat) line

      select case (line)

      case ('$Nodes')
        backspace (file_id)
        call read_meshsize_nodes(file_id, mesh)

      case ('$Elements')
        backspace (file_id)
        call read_meshsize_elts(file_id, mesh)

      case ('$NodePartitions')
        backspace (file_id)
        call read_meshsize_nodeparts(file_id)

      end select
    end do

    close (file_id)

    return

  end subroutine read_meshsize

  !===========================================================================

  subroutine read_mesh(mesh_file, opt_file, ori_file, phase_file, mesh)

    ! Parse the .msh file and have each processor store only what it needs.

    ! Locals:
    ! eofstat: Value that confirms if fortran is at the end of a file.
    ! i: Generic looping index.
    ! line: Input line on current record to be parsed.

    character(len=*), intent(in) :: mesh_file
    character(len=*), intent(in) :: opt_file
    character(len=*), intent(in) :: ori_file
    character(len=*), intent(in) :: phase_file
    type(mesh_type), intent(inout) :: mesh
    

    integer :: iostatus
    integer :: eofstat
    character(len=10) :: version
    character(len=256) :: line
    character(len=32)  :: iarray(16)
    integer :: status, nlines, i
    integer, allocatable :: elset_ids_inv(:)
    integer :: nb_elsets
    integer, allocatable :: elt_elset(:)
    real(rk), allocatable :: elset_oris(:, :)
    real(rk), allocatable :: elset_crss(:, :)
    real(rk), allocatable :: elset_sat_str(:, :)
    integer :: file_id

    logical :: eltoris_defined = .false.
    logical :: eltcrss_defined = .false.
    logical :: eltsat_str_defined = .false.

    open (newunit=file_id, file=mesh_file, status='old', action='read', &
        & iostat=iostatus)

    if (myid .eq. 0) then
      write (*, '(a,a,a)') "Info   :   [i] Parsing file `", trim(mesh_file), "'..."
      write (*, '(a)') 'Info   :   - Mesh parameters:'
      write (*, '(a,i0)') 'Info   :     > Node number: ', mesh%num_nodes
      write (*, '(a,i0)') 'Info   :     > Elt  number: ', mesh%num_elts
    end if

    ! Notes:
    ! The .msh file follows the Gmsh format 2.2 and is parsed in chunks
    ! and the field headers ($FieldName) are read to determine which
    ! helper subroutines to call.

    !---------------------------------------------------------------------------
    status = allocate_msh (mesh, elt_elset)

    if (status .ne. 0) then
      call par_quit('Error  :     > Failed to allocate mesh.')
    end if

    ! Initialize eofstat and loop over the entire file until eof is reached.
    eofstat = 0

    do while (eofstat .eq. 0)
      ! Read in line and determine which section is to be parsed.
      ! The below line is a temporary fix to avoid program hangs - jc
      line = ""

      read (file_id, '(a)', iostat=eofstat) line

      select case (line)

      case ('$MeshFormat')
        backspace (file_id)
        call read_mesh_format(file_id)

      case ('$MeshVersion')
        backspace (file_id)
        call read_mesh_version(file_id, version)

      case ('$Nodes')
        backspace (file_id)
        call read_nodes(file_id, mesh)

      case ('$Elements')
        backspace (file_id)
        call read_elts(file_id, mesh, elt_elset, elset_ids_inv, nb_elsets)

      case ('$NSets') ! Nodal values for BCs assignments (optional).
        backspace (file_id)
        call read_nsets(file_id, mesh)
      
      case ('$Periodicity') ! Periodicity relations for periodic BCs (optional).
        backspace (file_id)
        call read_periodicity(file_id, mesh)

      case ('$Fasets')
        backspace (file_id)
        call read_fasets(file_id, mesh)

      case ('$PhysicalNames')
        backspace (file_id)
        call read_physicalnames(file_id)

      case ('$ElsetOrientations') ! Per-elset ori.
        backspace (file_id)
        call read_elsetorientations(file_id, mesh, eltoris_defined, elset_ids_inv, elset_oris)

      case ('$ElementOrientations') ! Per-element ori (optional).
        backspace (file_id)
        call read_eltorientations(file_id, mesh, eltoris_defined, elset_oris)

      case ('$Groups') ! Grain/phase assignment for multiphase (optional).
        backspace (file_id)
        call read_groups(file_id, mesh, elt_elset, elset_ids_inv)

      end select
    end do

    if (myid .eq. 0) then
      write (*, '(a,a,a)') "Info   :   [i] Parsed file `", trim(mesh_file), "'."
    end if

    close (file_id)

    if ((myid .eq. 0) .and. num_procs .gt. 1 .and. (mesh%num_nodes/num_procs .le. 50)) then
      write (*, '(a)') 'Warning:     > The average number of nodes per &
          &processor is small. mpi '
      write (*, '(a)') '               issues may cause this simulation &
          &to hang.'
    end if

    
    if (opt_file .ne. "") then
      open (newunit=file_id, file=opt_file, status='old', action='read', &
          & iostat=iostatus)

      ! Initialize eofstat and loop over the entire file until eof is reached.
      eofstat = 0

      if (iostatus .ne. 0) then
        call par_quit("Error  :     > Failure to open `simulation.opt' file.")
      end if

      if (myid .eq. 0) then
        write (*, '(a,a,a)') "Info   :   [i] Parsing file `", trim(opt_file), "'..."
      end if

      do while (eofstat .eq. 0)
        ! Read in line and determine which section is to be parsed.
        ! The below line is a temporary fix to avoid program hangs - jc
        line = ""
  
        read (file_id, '(a)', iostat=eofstat) line
        
        select case (line)

        case ('$ElsetCrss') ! Per-elset crss (optional).
          backspace (file_id)
          if ((mesh%g_0_from_file .eqv. .true.) .and. (eltcrss_defined .eqv. .false.)) then
            call read_elsetflag(file_id, eltcrss_defined, elset_crss,'ElsetCrss', &
            & mesh)
            mesh%crss_defined = .true.
          else ! if not reading skip lines
            read (file_id, '(a)', iostat=eofstat) line
            read (line, *) iarray(1:2)
            read (iarray(1), *) nlines
            do i = 1, nlines
              read (file_id, *) line
            end do
          end if
        case ('$ElementCrss') ! Per-element crss (optional).
          backspace (file_id)
          if (mesh%g_0_from_file .eqv. .true.) then
            call read_eltflag(file_id, eltcrss_defined, elset_crss,'ElementCrss', &
            & mesh)
            mesh%crss_defined = .true.
          else ! if not reading skip lines 
            read (file_id, '(a)', iostat=eofstat) line
            read (line, *) iarray(1:2)
            read (iarray(1), *) nlines
            do i = 1, nlines
              read (file_id, *) line
            end do
          end if

        case ('$ElsetCrssSat') ! Per-elset saturation strength (optional).
          backspace (file_id)
          if ((mesh%g_s_from_file .eqv. .true.) .and. (eltsat_str_defined .eqv. .false.)) then
            call read_elsetflag(file_id, eltsat_str_defined, elset_sat_str,'ElsetCrssSat', &
                & mesh)
            mesh%sat_str_defined = .true.
          else ! if not reading skip lines
            read (file_id, '(a)', iostat=eofstat) line
            read (line, *) iarray(1:2)
            read (iarray(1), *) nlines
            do i = 1, nlines
              read (file_id, *) line
            end do
          end if

        case ('$ElementCrssSat') ! Per-element saturation strength (optional).
          backspace (file_id)
          if (mesh%g_s_from_file .eqv. .true.) then
            call read_eltflag(file_id, eltsat_str_defined, elset_sat_str,'ElementCrssSat', &
                & mesh)
            mesh%sat_str_defined = .true.
          else ! if not reading skip lines
            read (file_id, '(a)', iostat=eofstat) line
            read (line, *) iarray(1:2)
            read (iarray(1), *) nlines
            do i = 1, nlines
              read (file_id, *) line
            end do
          end if
        end select

      end do  

      if (myid .eq. 0) then
        write (*, '(a,a,a)') "Info   :   [i] Parsed file `", trim(opt_file), "'."
      end if
      close (file_id) 

    end if

    if (ori_file .ne. "") then

      open (newunit=file_id, file=ori_file, status='old', action='read', &
          & iostat=iostatus)

      read (file_id, '(a)', iostat=eofstat) line

      select case (line)

      case ('$ElsetOrientations') ! Per-elset ori.
        backspace (file_id)
        call read_elsetorientations(file_id, mesh, eltoris_defined, elset_ids_inv, elset_oris)

      case ('$ElementOrientations') ! Per-element ori (optional).
        backspace (file_id)
        call read_eltorientations(file_id, mesh, eltoris_defined, elset_oris)

      end select

      close (file_id)
    end if

    if (phase_file .ne. "") then

      open (newunit=file_id, file=phase_file, status='old', &
          & action='read', iostat=iostatus)

      if (iostatus .ne. 0) then
        call par_quit("Error  :     > Failure to open `simulation.phase' file.")
      end if

      if (myid .eq. 0) then
        write (*, '(a)') "Info   :   [i] Parsing file `simulation.phase'..."
      end if

      call read_groups(file_id, mesh, elt_elset, elset_ids_inv)

      if (myid .eq. 0) then
        write (*, '(a)') "Info   :   [i] Parsed file `simulation.phase'."
      end if

      close (file_id)
    end if
    !
    if (mesh%crss_defined .eqv. .true.) then

      allocate (mesh%crss(size(elset_crss,1), elt_sub:elt_sup,nqpt))

      call read_mesh_init_eltcrss(mesh, eltcrss_defined, elset_crss, elt_elset, nb_elsets)

    end if
    !
    if (mesh%sat_str_defined .eqv. .true.) then

      allocate (mesh%sat_str(size(elset_sat_str,1), elt_sub:elt_sup))

      call read_mesh_init_eltsat_str(mesh, eltsat_str_defined, elset_sat_str, elt_elset, nb_elsets)

    end if
    !
    allocate (mesh%ori(3, 3, elt_sub:elt_sup))

    call read_mesh_init_eltoris(mesh, version, eltoris_defined, elt_elset, elset_oris, nb_elsets)

    return

  end subroutine read_mesh

end module read_input_msh_mod
