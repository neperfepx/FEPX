! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module read_input_msh_mod2

! allocate_msh: Allocates mesh arrays based on the problem size.

  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod
  use types_mod
  use crys_type_mod2
  use orientation_conversion_mod
  use gather_scatter_mod
  use surface_mod
  use parallel_mod

  implicit none

  public

contains

  function allocate_msh(mesh, elt_elset) result(status)

    ! Allocate mesh arrays according to problem size.

    !---------------------------------------------------------------------------

    implicit none

    ! Arguments:
    ! status: Returns status of allocation to main program.

    type(mesh_type), intent(inout) :: mesh
    integer, allocatable :: elt_elset(:)

    integer :: status

    !---------------------------------------------------------------------------

    status = 0

    allocate (mesh%elt_nodes(ndim, elt_sub:elt_sup), &
        & mesh%elt_dofs(kdim, elt_sub:elt_sup), &
        & elt_elset(mesh%num_elts), &
        & mesh%elt_phase(mesh%num_elts), &
        & stat=status)

    ! Initialize these arrays here before the spatial mesh is parsed.

    ! Note: elt_elset stores the elset for an element
    ! and we set phase to 1 to assume a single-phase simulation by default.

    elt_elset = 0
    mesh%elt_phase = 1

    return

  end function allocate_msh

  subroutine read_mesh_format(io)

    ! Read the mesh file format and confirm it is of correct type.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Input unit for .msh file.

    integer, intent(in) :: io

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! nspace/nvals: Storage values used for dividing a line into substrings.
    ! file_type: Gmsh defined file type of mesh - must be '0'.
    ! data_size: Gmsh defined data size of mesh - must be '8'.
    ! mesh_format: Gmsh defined mesh file format - must be '2.2'.
    ! iarray: Substring array for internal read parsing.
    ! line: Input line on current record to be parsed.

    integer :: ierr, nspace, nvals
    integer :: file_type, data_size
    character(len=3)   :: mesh_format
    character(len=12)  :: iarray(16)
    character(len=256) :: line

    !---------------------------------------------------------------------------

    ! Read the first line and confirm it is the correct record.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$MeshFormat')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &mesh format.')
    end if

    ! Read in the mesh format string.
    read (io, '(a)', iostat=ierr) line

    ! Trim the full string into nvals number of substrings to parse.
    nspace = count(transfer(line, 'a', len_trim(line)) .eq. " ")
    nvals = nspace + 1
    iarray = ""

    ! Internal read of line to store substrings.
    read (line, *) iarray(1:nvals)
    read (iarray(1), *) mesh_format
    read (iarray(2), *) file_type
    read (iarray(3), *) data_size

    if (mesh_format .ne. '2.2') &
        &call par_quit('Error  :     > Incorrect mesh format version provided.')
    if (file_type .ne. 0) &
        &call par_quit('Error  :     > Incorrect mesh file type provided.')
    if (data_size .ne. 8) &
        &call par_quit('Error  :     > Incorrect mesh data size provided.')

    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndMeshFormat')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &mesh format.')
    end if

    return

  end subroutine read_mesh_format

  !===========================================================================

  subroutine read_mesh_version(io)

    ! Read the mesh file version

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Input unit for .msh file.

    integer, intent(in) :: io

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! line: Input line on current record to be parsed.

    integer :: ierr
    character(len=256) :: line
    character(len=256) :: version

    !---------------------------------------------------------------------------

    ! Read the first line and confirm it is the correct record.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$MeshVersion')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &mesh version.')
    end if

    ! Read in the mesh version string.
    read (io, '(a)', iostat=ierr) version

    if (version(1:3) .ne. '2.2') &
        &call par_quit('Error  :     > Incorrect mesh version provided.')

    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndMeshVersion')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &mesh version.')
    end if

    return

  end subroutine read_mesh_version

  !===========================================================================

  subroutine read_nodes(io, mesh)

    ! Read in node id and coordinates and store per processor.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Input unit for .msh file.

    integer, intent(in) :: io
    type(mesh_type) :: mesh

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! inode: Node id value.
    ! nlines: Read-in value for the number of lines that will be parsed.
    ! i: Generic loop index.
    ! k1/k2/k3: dof mapping values for the node coordinates.
    ! x/y/z: Cartesian coordinate value of the node.
    ! line: Input line on current record to be parsed.

    integer  :: ierr, inode, nlines
    integer  :: i, k1, k2, k3
    real(rk) :: x, y, z
    character(len=256) :: line

    !---------------------------------------------------------------------------

    allocate (mesh%coo(dof_sub:dof_sup))

    ! Read the first line and confirm it is the correct record.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$Nodes')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &nodes.')
    end if

    ! Read in the number of nodes in the section (and set nlines).
    read (io, *) nlines

    do i = 1, nlines
      read (io, *) inode, x, y, z

      ! Store the nodal coords only if within local processor range.
      if ((inode .ge. node_sub) .and. (inode .le. node_sup)) then
        k1 = 3*(inode - 1) + 1
        k2 = k1 + 1
        k3 = k2 + 1

        mesh%coo(k1) = x
        mesh%coo(k2) = y
        mesh%coo(k3) = z
      end if
    end do

    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndNodes')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &nodes.')
    end if

    return

  end subroutine read_nodes

  !===========================================================================

  subroutine read_elts(io, mesh, elt_elset, elset_ids_inv, nb_elsets)

    ! Read in all elements [0d -> 3d] and only store 3d elements within
    ! local processor range.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Input unit for .msh file.

    integer, intent(in) :: io
    type(mesh_type) :: mesh
    integer :: elt_elset(:)
    integer, allocatable, intent(out) :: elset_ids_inv(:)
    integer, intent(out) :: nb_elsets

    integer, allocatable :: elset_ids(:)

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! i/j: Generic loop index.
    ! j1/j2/k3/k1/k2/k3: dof mapping values for the element nodes.
    ! nlines: Read-in value for the number of lines that will be parsed.
    ! nspace/nvals: Storage values used for dividing a line into substrings.
    ! elttype: Gmsh defined value that determines the element type.
    ! nelt: Locally determined number of 3d elements in mesh.
    ! ntags: Number of tags in the 3d element line - must be 3.
    ! elset_id: Grain (or elset) id
    ! nodes_fe: Local element nodal connectivity array.
    ! iarray: Substring array for internal read parsing.
    ! line: Input line on current record to be parsed.

    integer :: ierr, i, j, j1, j2, j3, k1, k2, k3, elset_id_max
    integer :: nlines, nspace, nvals, elttype, nelt, ntags, elset_id, elset_pos
    integer :: nodes_fe(ndim)
    character(len=12) :: iarray(16)
    character(len=256) :: line

    !---------------------------------------------------------------------------

    ! Read the first line and confirm it is the correct record.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$Elements')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &elements.')
    end if

    ! Read in the number of elements in the section (and set nlines).
    read (io, *) nlines
    allocate (elset_ids(nlines))

    nelt = 0
    nb_elsets = 0

    do i = 1, nlines
      read (io, '(a)') line

      ! Trim the full string into nvals number of substrings to parse.
      nspace = count(transfer(line, 'a', len_trim(line)) .eq. " ")
      nvals = nspace + 1
      iarray = ""

      ! Internal read of line to store substrings.
      read (line, *) iarray(1:nvals)
      read (iarray(2), *) elttype

      if (elttype .eq. 11) then
        ! Increment nelt by 1
        nelt = nelt + 1

        ! If 10-node tetrahedral element, extract information from the line.
        ! The line has format of:
        ! ielt elttype ntags elset_id tag2 tag3 node0 ... node9
        read (iarray(3), *) ntags
        read (iarray(4), *) elset_id ! elset

        ! Do we know this elset?
        elset_pos = -1
        do j = 1, nb_elsets
          if (elset_ids(j) .eq. elset_id) then
            elset_pos = j
            ! escape here
          end if
        end do

        if (elset_pos .eq. -1) then
          nb_elsets = nb_elsets + 1
          elset_pos = nb_elsets
          elset_ids(elset_pos) = elset_id
        end if

        ! Set the elt_elset value for this element (assigns elset/elset).
        elt_elset(nelt) = elset_pos

        if (ntags .eq. 3) then
          ! The .msh format modifies the local node order of an element.
          ! The remapping follows from Gmsh to FEPX:
          ! gmsh order: 1 2 3 4 5 6 7 8 9 10
          ! FEPX order: 1 3 5 10 2 4 6 7 9 8
          read (iarray(7), *) nodes_fe(1)
          read (iarray(8), *) nodes_fe(3)
          read (iarray(9), *) nodes_fe(5)
          read (iarray(10), *) nodes_fe(10)
          read (iarray(11), *) nodes_fe(2)
          read (iarray(12), *) nodes_fe(4)
          read (iarray(13), *) nodes_fe(6)
          read (iarray(14), *) nodes_fe(7)
          read (iarray(15), *) nodes_fe(9)
          read (iarray(16), *) nodes_fe(8)

          ! Store elements only on local processor range.
          if ((nelt .ge. elt_sub) .and. (nelt .le. elt_sup)) then
            mesh%elt_nodes(:, nelt) = nodes_fe

            ! Map the local dof to global dof.
            do j = 1, ndim
              j1 = 3*(j - 1) + 1
              j2 = j1 + 1
              j3 = j2 + 1

              k1 = 3*(nodes_fe(j) - 1) + 1
              k2 = k1 + 1
              k3 = k2 + 1

              mesh%elt_dofs(j1, nelt) = k1
              mesh%elt_dofs(j2, nelt) = k2
              mesh%elt_dofs(j3, nelt) = k3
            end do
          end if
        end if
      end if
    end do

    elset_id_max = 0

    do i = 1, nb_elsets
      elset_id_max = max(elset_ids(i), elset_id_max)
    end do

    allocate (elset_ids_inv(elset_id_max))
    elset_ids_inv = 0

    do i = 1, nb_elsets
      elset_ids_inv(elset_ids(i)) = i
    end do
    
    mesh%num_elsets = nb_elsets
    ! Confirm that nelt (local) and mesh%num_elts (global) values match
    if (nelt .ne. mesh%num_elts) call par_quit('Error  :     > &
        &Number of elements in $Elements does not match problem size.')

    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndElements')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &elements.')
    end if

    return

  end subroutine read_elts

  !===========================================================================

  subroutine read_nsets(io, mesh)

    ! Read in node sets and do not store - currently unused, but will be
    ! integrated into custom boundary conditions in the future.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Input unit for .msh file.

    integer, intent(in) :: io
    type(mesh_type) :: mesh

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! i/j: Generic loop index.
    ! num_nsets: Number of NSets within the field.
    ! line: Input line on current record to be parsed.

    integer :: ierr, i, j, num_nsets
    character(len=256) :: line

    !---------------------------------------------------------------------------

    ! Read the first line and confirm it is the correct record.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$NSets')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &nsets.')
    end if

    ! Read in the number of nsets in the section (and set nsets).
    read (io, *) num_nsets

    ! Allocate the overall type, nsets, to accomodate num_nsets (total num.)
    allocate (mesh%nsets(num_nsets))

    ! Loop over the number of nsets and extract information per set.
    do i = 1, num_nsets
      ! Read in the NSet label.
      read (io, *) mesh%nsets(i)%nset_label

      ! Read in the number of nodes in this set, allocate nset_nodes
      read (io, *) mesh%nsets(i)%num_nset_nodes
      allocate (mesh%nsets(i)%nset_nodes(mesh%nsets(i)%num_nset_nodes))

      ! Loop over the NSet and store all associated nodes
      do j = 1, mesh%nsets(i)%num_nset_nodes
        ! Read to vector
        read (io, *) mesh%nsets(i)%nset_nodes(j)
      end do
    end do

    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndNSets')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &nsets.')
    end if

    return

  end subroutine read_nsets

  !===========================================================================

  subroutine read_periodicity(io, mesh)

    ! Read periodicity relations

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Input unit for .msh file.

    integer, intent(in) :: io
    type(mesh_type) :: mesh

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! i/j: Generic loop index.
    ! line: Input line on current record to be parsed.

    integer :: ierr, i, prim, sec, pvectx, pvecty, pvectz, sum_pvect
    character(len=256) :: line

    !---------------------------------------------------------------------------

    ! Read the first line and confirm it is the correct record.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$Periodicity')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &periodicity.')
    end if

    ! Turn on flag for periodic mesh
    mesh%periodic_mesh = .true.

    ! Read in the number of periodicity relations in the section (and set periodicity relation).
    read (io, *) mesh%num_periodicity

    ! Allocate the overall type, periodicity, to accomodate mun_periodicity (total num.)
    allocate (mesh%periodicity(mesh%num_periodicity))

    ! Loop over the number of periodicity relations and extract information per set.
    do i = 1, mesh%num_periodicity
      ! Read in the NSet label.
      read (io, *) sec, prim, pvectx, pvecty, pvectz

      mesh%periodicity(i)%primary=prim
      mesh%periodicity(i)%secondary=sec
      mesh%periodicity(i)%pvect(1)=pvectx
      mesh%periodicity(i)%pvect(2)=pvecty
      mesh%periodicity(i)%pvect(3)=pvectz

      sum_pvect=100*pvectx+10*pvecty+pvectz

      ! label periodic vector directions
      select case (sum_pvect)
      
      case (100)
        mesh%periodicity(i)%pvect_label = 1

      case (10)
        mesh%periodicity(i)%pvect_label = 2

      case (1)
        mesh%periodicity(i)%pvect_label = 3

      case (110)
        mesh%periodicity(i)%pvect_label = 4

      case (90)
        mesh%periodicity(i)%pvect_label = 5

      case (101)
        mesh%periodicity(i)%pvect_label = 6

      case (99)
        mesh%periodicity(i)%pvect_label = 7

      case (11)
        mesh%periodicity(i)%pvect_label = 8

      case (9)
        mesh%periodicity(i)%pvect_label = 9

      case (111)
        mesh%periodicity(i)%pvect_label = 10

      case (109)
        mesh%periodicity(i)%pvect_label = 11

      case (91)
        mesh%periodicity(i)%pvect_label = 12
      
      case (-100)
        mesh%periodicity(i)%pvect_label = -1

      case (-10)
        mesh%periodicity(i)%pvect_label = -2

      case (-1)
        mesh%periodicity(i)%pvect_label = -3

      case (-110)
        mesh%periodicity(i)%pvect_label = -4

      case (-90)
        mesh%periodicity(i)%pvect_label = -5

      case (-101)
        mesh%periodicity(i)%pvect_label = -6

      case (-99)
        mesh%periodicity(i)%pvect_label = -7

      case (-11)
        mesh%periodicity(i)%pvect_label = -8

      case (-9)
        mesh%periodicity(i)%pvect_label = -9

      case (-111)
        mesh%periodicity(i)%pvect_label = -10

      case (-109)
        mesh%periodicity(i)%pvect_label = -11

      case (-91)
        mesh%periodicity(i)%pvect_label = -12

      end select

    end do

    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndPeriodicity')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &nsets.')
    end if

    return

  end subroutine read_periodicity

  !===========================================================================

  subroutine read_fasets(io, mesh)

    ! Read in facesets and store the surface 2d elements. This subroutine also
    ! prepares the 3d connectivity for the surface nodes which is later used
    ! to integrate the loads on each surface.

    ! Notes:
    ! This information should be readily extracted from the read_elts
    ! subroutine by searching for elt_type (9) and extracting, however, this
    ! would not provide faset labels or IDs.

    ! Legacy note: type = 6 for triangle - tsh

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Input unit for .msh file.

    integer, intent(in) :: io
    type(mesh_type), intent(inout):: mesh

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! i/j/k: Generic loop index.
    ! type: Parameters to define the internal element type/bounds.
    ! elt_id: Surface element id - 1-indexed.
    ! elt_nodes/elt_nodes_tmp: Surface element connectivity.
    ! num_elts: Number of elements on a given surface.
    ! semin/semax: Partitioned array bounds to distribute surface sections.
    ! status: Confirms if allocation of surface section arrays was successful.
    ! j3/n3: 3d connectivity dof mapping values.
    ! ierr: Value that confirms if a read() fails.
    ! elt_dof_min/elt_dof_max: dof mapping to scatter surface elements.
    ! fasets: Number of fasets in the field to be parsed.
    ! line: Input line on current record to be parsed.

    integer, parameter :: type = 6
    integer :: elt_id, elt_nodes(type), elt_nodes_tmp(type)
    integer :: ii, i, j, k, num_elts, semin, semax, status, j3, n3, ierr
    integer :: elt_dof_min, elt_dof_max
    character(len=256) :: line

    !---------------------------------------------------------------------------

    ! Read the first line and confirm it is the correct record.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$Fasets')) then
      call par_quit('Error  :     > Parse error attempting to read in mesh%num_fasets.')
    end if

    ! Read in the number of fasets in the section (and set fasets).
    read (io, *) mesh%num_fasets

    allocate (mesh%faset_labels(mesh%num_fasets))
    allocate (mesh%fasets(mesh%num_fasets))

    ii = 1

    do i = 1, mesh%num_fasets
      ! Read the faset label
      read (io, *) mesh%faset_labels(ii)
      ! number of elements on the surface.
      read (io, *) num_elts

      ! skip empty fasets
      if (num_elts .ne. 0) then
        
        ! Should semin/semax use the the surface elements from the 3d elements
        ! that are stored on the same processor range instead of something else?
        call par_partition(num_elts, num_procs, myid, semin, semax)
        status = surf_alloc(type, semin, semax, mesh%fasets(ii))

        if (status .ne. 0) then
          call par_quit('Error  :     > Surface allocation error.')
        end if

        do j = 1, num_elts
          ! For each element, read in the id and 3d connectivity.
          read (io, *) elt_id, elt_nodes_tmp

          ! Shift the id for msh file version > 2.2
          ! if (mesh_version .ne. '2.2') then  -- commented rq 03/13/23
          elt_id = elt_id - mesh%elt_startid + 1
          ! end if

          if (elt_id .lt. 1 .or. elt_id .gt. mesh%num_elts) then
            call par_quit('Error  :     > Element index out of bounds')
          end if

          ! The node values need to be reordered.
          ! gmsh order: 1 2 3 4 5 6
          ! FEPX order: 6 4 2 5 3 1
          elt_nodes(1) = elt_nodes_tmp(6)
          elt_nodes(2) = elt_nodes_tmp(3)
          elt_nodes(3) = elt_nodes_tmp(5)
          elt_nodes(4) = elt_nodes_tmp(2)
          elt_nodes(5) = elt_nodes_tmp(4)
          elt_nodes(6) = elt_nodes_tmp(1)

          if ((j .le. semax) .and. (j .ge. semin)) then
            ! Create the 2d connectivity first.
            mesh%fasets(ii)%econn(1, j) = 6*(elt_id - 1) + 1
            mesh%fasets(ii)%econn(2, j) = 6*(elt_id - 1) + 2
            mesh%fasets(ii)%econn(3, j) = 6*(elt_id - 1) + 3
            mesh%fasets(ii)%econn(4, j) = 6*(elt_id - 1) + 4
            mesh%fasets(ii)%econn(5, j) = 6*(elt_id - 1) + 5
            mesh%fasets(ii)%econn(6, j) = 6*(elt_id - 1) + 6

          ! Store into nodes and id into global connectivity?
          mesh%fasets(ii)%conn(:, j) = elt_nodes(:)
          mesh%fasets(ii)%elt(j) = elt_id
          end if
        end do

        ! Create a 3d connectivity for the coordinate gather.

        do k = semin, semax
          do j = 1, type
            j3 = 3*(j - 1) + 1
            n3 = 3*(mesh%fasets(ii)%conn(j, k) - 1) + 1
            mesh%fasets(ii)%conn3d(j3, k) = n3
            mesh%fasets(ii)%conn3d(j3 + 1, k) = n3 + 1
            mesh%fasets(ii)%conn3d(j3 + 2, k) = n3 + 2
          end do
        end do

        call part_scatter_setup(1, 3 * type, dof_sub, dof_sup,&
            & semin, semax, mesh%fasets(ii)%conn3d, mesh%fasets(ii)%tr)
        elt_dof_min = 6*(elt_sub - 1) + 1
        elt_dof_max = 6*elt_sup
        call part_scatter_setup(1, 6, elt_dof_min, elt_dof_max,&
            & semin, semax, mesh%fasets(ii)%econn, mesh%fasets(ii)%etr)

        ii = ii + 1

      end if
      end do

      mesh%num_fasets = ii - 1
    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndFasets')) then
      call par_quit('Error  :     > Parse error attempting to read in mesh%num_fasets.')
    end if

    return

  end subroutine read_fasets

  !===========================================================================

  subroutine read_meshsize_nodes(io, mesh)

    ! Read the $Nodes field and retrieves the total number of nodes in mesh.

    integer, intent(in) :: io
    type(mesh_type) :: mesh

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! i: Generic looping index.
    ! line: Input line on current record to be parsed.

    integer :: ierr, i
    character(len=256) :: line

    !---------------------------------------------------------------------------

    ! Read the first line and confirm it is the correct record.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$Nodes')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &nodes.')
    end if

    ! Read in the number of nodes in the section (and set mesh%num_nodes).
    read (io, *) mesh%num_nodes

    do i = 1, mesh%num_nodes
      read (io, *)
    end do

    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndNodes')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &nodes.')
    end if

    return

  end subroutine read_meshsize_nodes

  !===========================================================================

  ! read_meshsize_elts: Gets the number of elements in the mesh and partitions.

  subroutine read_meshsize_elts(io, mesh)

    ! Read the $Elements field and retrieves the total number of nodes in mesh.
    ! Also, retrieves elemental partition information if present in the mesh.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Input unit for .msh file.

    integer, intent(in) :: io
    type(mesh_type) :: mesh

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! i: Generic looping index.
    ! nlines: Read-in value for the number of lines that will be parsed.
    ! nspace/nvals: Storage values used for dividing a line into substrings.
    ! elttype: Gmsh defined value that determines the element type.
    ! nelt: Locally determined number of 3d elements in mesh.
    ! minbnd/maxbnd: Array bounds to assist in trimming element partition array.
    ! nparts: Number of partitions read in from the mesh.
    ! ipart: Partition loop index to check which partition is being checked.
    ! line: Input line on current record to be parsed.
    ! iarray: Substring array for internal read parsing.
    ! temp: Temporary array storage for element paritions including non-3d.
    ! elt_parts: Array storing the trimmed element paritions.
    ! tempval: Array of elt_sub/elt_sup bounds for each partition read in.
    ! chg_index: Index in array immediately after a change in value occurs.
    ! mask: Masking array to trim temp into elt_parts.

    integer :: ierr, i
    integer :: nlines, nspace, nvals, elttype, nelt
    integer :: minbnd, maxbnd, nparts, ipart
    character(len=256) :: line
    character(len=12)  :: iarray(16)
    integer, allocatable :: temp(:), elt_parts(:), tempval(:, :), chg_index(:)
    logical, allocatable :: mask(:)
    logical :: read_elt_part = .false.

    !---------------------------------------------------------------------------

    ! Read the first line and confirm it is the correct record.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$Elements')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &elements.')
    end if

    ! Read in the number of elements in the section (and set nlines).
    read (io, *) nlines

    ! Initialize the 3d elem counter and partition array.
    nelt = 0

    allocate (temp(nlines))
    allocate (mask(nlines))
    temp = -1

    do i = 1, nlines
      read (io, '(a)') line

      ! Trim the full string into nvals number of substrings to parse.
      nspace = count(transfer(line, 'a', len_trim(line)) .eq. " ")
      nvals = nspace + 1
      iarray = ""

      ! Internal read of line to store substrings.
      read (line, *) iarray(1:nvals)
      read (iarray(2), *) elttype

      if (elttype .eq. 11) then
        if (mesh%elt_startid .eq. -1) then
          read (iarray(1), *) mesh%elt_startid
        end if

        ! Assign data from iarray(6) the tag3 value which contain partition.
        read (iarray(6), *) temp(i)

        ! Increment nelt by 1 to ensure mesh%num_elts is assigned correctly.
        nelt = nelt + 1
      end if
    end do

    ! Check to see if any elttype=11 have been parsed (if not, quit)
    if (nelt .eq. 0) then
      call par_quit('Error  :     > Mesh requires second-order elements.')
    end if

    ! Set the global num_elts value.
    mesh%num_elts = nelt

    ! Reallocate temp array to properly sized array (1:num_elts).
    mask = temp .eq. -1
    minbnd = count(mask) + 1
    maxbnd = size(temp)

    allocate (elt_parts(maxbnd-minbnd+1))
    elt_parts = temp(minbnd:maxbnd)
    deallocate (mask)
    deallocate (temp)

    ! Set a flag if we read in paritions and strip the elt_parts array down.
    if (maxval(elt_parts) .gt. 0) then
      read_elt_part = .true.
      nparts = maxval(elt_parts)
      allocate (tempval(nparts, 2))
      allocate (chg_index(nparts - 1))
      ipart = 1

      ! We can split elt_parts
      do i = 1, size(elt_parts)
        ! a change has been detected so do something.
        if (elt_parts(i) .ne. ipart) then
          chg_index(ipart) = i
          if (ipart .ne. nparts) ipart = ipart + 1
        end if
      end do

      ! Build bound array for nparts -> (1,:) is elt_sub, (2,:) is elt_sup.
      tempval = 0
      do i = 1, nparts
        ! Need to handle the edge cases differently.
        if (i .eq. 1) then
          tempval(1, 1) = 1
          tempval(1, 2) = chg_index(i) - 1

        else if (i .eq. nparts) then
          tempval(i, 1) = chg_index(i - 1)
          tempval(i, 2) = size(elt_parts)

        else
          tempval(i, 1) = chg_index(i - 1)
          tempval(i, 2) = chg_index(i) - 1
        end if
      end do
    end if

    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndElements')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &elements.')
    end if

    ! Do not need to maintain element partitions for now so deallocate.
    if (allocated(elt_parts)) deallocate (elt_parts)
    if (allocated(chg_index)) deallocate (chg_index)
    if (allocated(tempval)) deallocate (tempval)

    return

  end subroutine read_meshsize_elts

  subroutine read_meshsize_nodeparts(io)

    ! Read in node paritions and do not store - currently unused, but will be
    ! integrated in the future for utilizing imprted scotch partitions.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Input unit for .msh file.

    integer, intent(in) :: io

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! i: Generic loop index.
    ! nlines: Read-in value for the number of lines that will be parsed.
    ! nparts: Number of partitions read in from the mesh.
    ! ipart: Partition loop index to check which partition is being checked.
    ! temp: Temporary storage for current node id - not stored.
    ! node_parts: Array storing the nodal paritions.
    ! tempval: Array of node_sub/node_sup bounds for each partition read in.
    ! chg_index: Index in array immediately after a change in value occurs.
    ! line: Input line on current record to be parsed.

    integer :: ierr, i, nlines, nparts, ipart, temp
    integer, allocatable :: node_parts(:), tempval(:, :), chg_index(:)
    character(len=256) :: line
    logical :: read_node_part = .false.

    ! Read the first line and confirm it is the correct record.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$NodePartitions')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &node partitions.')
    end if

    ! Set the flag if this section is being read in
    read_node_part = .true.

    ! Read in the number of node partitions in the section (and set nlines).
    read (io, *) nlines

    ! Allocate the nodal partition array -> n(1,:) is 1-index id, n(2,:) is part.
    allocate (node_parts(nlines))

    ! Loop over the number of nlines and do not store information.
    do i = 1, nlines
      read (io, *) temp, node_parts(i)
    end do

    nparts = maxval(node_parts)
    allocate (tempval(nparts, 2))
    allocate (chg_index(nparts - 1))
    ipart = 1

    ! We can split elt_parts
    do i = 1, size(node_parts)
      ! a change has been detected so do something
      if (node_parts(i) .ne. ipart) then
        chg_index(ipart) = i
        if (ipart .ne. nparts) ipart = ipart + 1
      end if
    end do

    ! Build the bound array for nparts -> (1,:) is node_sub, (2,:) is node_sup.
    tempval = 0
    do i = 1, nparts
      ! Handle the edge cases differently
      if (i .eq. 1) then
        tempval(1, 1) = 1
        tempval(1, 2) = chg_index(i) - 1

      else if (i .eq. nparts) then
        tempval(i, 1) = chg_index(i - 1)
        tempval(i, 2) = size(node_parts)

      else
        tempval(i, 1) = chg_index(i - 1)
        tempval(i, 2) = chg_index(i) - 1
      end if
    end do

    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndNodePartitions')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &node partitions.')
    end if

    ! Do not need to maintain node partitions for now so deallocate.
    deallocate (node_parts)
    deallocate (chg_index)
    deallocate (tempval)

    return

  end subroutine read_meshsize_nodeparts

  !===========================================================================

  subroutine read_physicalnames(io)

    ! Read in physical names and do not store - will never be used as these
    ! are for gmsh internal use only and must be stored to successfully open
    ! the msh file in any version of gmsh.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Input unit for .msh file.

    integer, intent(in) :: io

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! i: Generic loop index.
    ! nlines: Read-in value for the number of lines that will be parsed.
    ! line: Input line on current record to be parsed.

    integer :: ierr, i, nlines
    character(len=256) :: line

    !---------------------------------------------------------------------------

    ! Read the first line and confirm it is the correct record.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$PhysicalNames')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &physical names.')
    end if

    ! Read in the number of physical names in the section (and set nlines).
    read (io, *) nlines

    ! Loop over the number of nlines and do not store information.
    do i = 1, nlines
      read (io, *)
    end do

    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndPhysicalNames')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &physical names.')
    end if

    return

  end subroutine read_physicalnames

  !===========================================================================

  subroutine read_elsetorientations(io, mesh, eltoris_defined, elset_ids_inv, elset_oris)

    ! Read elset orientation information from either the msh file or an external
    ! simulation.ori file. If we are reading in ori from an external
    ! file we must deallocate the current storage array and start over.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Input unit for .msh file.

    integer, intent(in) :: io
    type(mesh_type) :: mesh
    logical :: eltoris_defined
    integer :: elset_ids_inv(:)
    real(rk), allocatable, intent(inout) :: elset_oris(:, :)

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! i/j: Generic loop index.
    ! s: Status value that read in string is a valid option.
    ! nlines: Read-in value for the number of lines that will be parsed.
    ! nspace/nvals: Storage values used for dividing a line into substrings.
    ! elset_id: Grain id - 1-indexed.
    ! delim_pos: Line position where the ':' delimiter is found.

    ! parm_string: String describing input orientation parameterization.
    ! conv_string: String describing input orientation convention.
    ! ori_string: Temp string used to split the 'descriptor:convention' line.
    ! iarray: Substring array for internal read parsing.
    ! line: Input line on current record to be parsed.

    integer :: ierr, i, j, s
    integer :: nlines, nspace, nvals, elset_id, delim_pos
    real(rk) :: ori(4)
    character(len=50)  :: parm_string, conv_string, ori_string
    character(len=32)  :: iarray(16)
    character(len=256) :: line

    !---------------------------------------------------------------------------

    ! If per-element ori have been read-in previously then skip this.
    if (eltoris_defined .eqv. .true.) then
      ! Read the first line and confirm it is the correct record.
      read (io, '(a)', iostat=ierr) line

      if ((ierr .lt. 0) .or. (line .ne. '$ElsetOrientations')) then
        call par_quit('Error  :     > Parse error attempting to read in &
            &elset ori.')
      end if

      ! Extract the number of lines to skip
      read (io, '(a)', iostat=ierr) line
      nspace = count(transfer(line, 'a', len_trim(line)) .eq. " ")
      nvals = nspace + 1
      iarray = ""

      ! Internal read of the line into the primary substring arrays.
      read (line, *) iarray(1:nvals)
      read (iarray(1), *) nlines

      ! Loop over nlines and do not store read in values.
      do i = 1, nlines
        read (io, *) line
      end do

      ! Read the end of section footer.
      read (io, '(a)', iostat=ierr) line

      if ((ierr .lt. 0) .or. (line .ne. '$EndElsetOrientations')) then
        call par_quit('Error  :     > Parse error attempting to read in &
            &elset ori.')
      end if

      return
    end if

    ! Check if the array is already allocated from reading the mesh file.
    if (allocated(elset_oris) .eqv. .true.) then
      deallocate (elset_oris)
    end if

    ! Read the first line and confirm it is the correct record.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$ElsetOrientations')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &elset ori.')
    end if

    ! Read in the section information string and prepare to parse.

    read (io, '(a)', iostat=ierr) line

    ! Trim the full string into two substring arrays split by the space delim.
    nspace = count(transfer(line, 'a', len_trim(line)) .eq. " ")
    nvals = nspace + 1
    iarray = ""

    ! Internal read of the line into the primary substring arrays.
    read (line, *) iarray(1:nvals)

    ! Internal read of the primary substring arrays to get the number of
    ! lines to set to parse and prepare to parse the orientation description.
    read (iarray(1), *) nlines
    read (iarray(2), '(a)') ori_string

    ! Find the position of the ":" character in ori_string as fortran
    ! refuses to allow the read() function to handle this automatically.
    delim_pos = index(ori_string, ":")

    ! Internal read to store the orientation parameterization and convention.
    read (ori_string(1:delim_pos - 1), *) parm_string
    read (ori_string(delim_pos + 1:len_trim(ori_string)), *) conv_string

    ! Set a per-element logical to false.
    eltoris_defined = .false.

    ! Parse the parameterization and confirm it is a valid option.

    s = 0

    read (parm_string, '(a)') &
        & mesh%orientation_parameterization

    ! Confirm that the rest of the string is a valid input.
    if (mesh%orientation_parameterization &
        & .eq. 'axis-angle') then
      s = 0

    else if (mesh%orientation_parameterization &
        & .eq. 'euler-bunge') then
      s = 0

    else if (mesh%orientation_parameterization &
        & .eq. 'euler-kocks') then
      s = 0

    else if (mesh%orientation_parameterization &
        & .eq. 'rodrigues') then
      s = 0

    else if (mesh%orientation_parameterization &
        & .eq. 'quaternion') then
      s = 0

    else ! Orientation parameterization defined is incorrect.
      s = 1
    end if

    if (s .eq. 1) then
      call par_quit("Error  :     > orientation_parameterization &
          & contains an error or unexpected input type.")
    end if

    ! Parse the convention and confirm it is a valid option.

    s = 0

    read (conv_string, '(a)') &
        & mesh%orientation_convention

    ! Confirm that the rest of the string is a valid input.
    if (mesh%orientation_convention &
        & .eq. 'active') then
      s = 0

    else if (mesh%orientation_convention &
        & .eq. 'passive') then
      s = 0

    else ! Orientation convention defined is incorrect.
      s = 1
    end if

    if (s .eq. 1) then
      call par_quit("Error  :     > orientation_convention &
          & contains an error or unexpected input type.")
    end if

    ! Check orientation parameterization to set column number.
    if ((mesh%orientation_parameterization &
        & .eq. 'axis-angle') .or. &
        & (mesh%orientation_parameterization &
        & .eq. 'quaternion')) then
      allocate (elset_oris(4, nlines))

    else ! All other parameterization have 3 values per orientation.
      allocate (elset_oris(3, nlines))
    end if

    ! Only quaternion or axis-angle.
    if (size(elset_oris, 1) .eq. 3 .or. size(elset_oris, 1) .eq. 4) then
      do i = 1, nlines
        read (io, *) elset_id, ori(1:size(elset_oris, 1))

        do j = 1, size(elset_oris, 1)
          elset_oris(j, elset_ids_inv(elset_id)) = ori(j)
        end do
      end do

    else
      call par_quit("Error  :     > elset_oris array not initialized&
          & correctly.")
    end if

    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndElsetOrientations')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &elset ori.')
    end if

    return

  end subroutine read_elsetorientations

  !===========================================================================

  subroutine read_eltorientations(io, mesh, eltoris_defined, elset_oris)

    ! Read elem orientation information from either the msh file or an external
    ! simulation.ori file. If we are reading in ori from an external
    ! file we must deallocate the current storage array and start over.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Input unit for .msh file.

    integer :: io
    type(mesh_type) :: mesh
    logical :: eltoris_defined
    real(rk), allocatable, intent(inout) :: elset_oris(:, :)

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! i/j: Generic loop index.
    ! s: Status value that read in string is a valid option.
    ! nlines: Read-in value for the number of lines that will be parsed.
    ! nspace/nvals: Storage values used for dividing a line into substrings.
    ! elt: Element id - 1-indexed.
    ! delim_pos: Line position where the ':' delimiter is found.
    ! parm_string: String describing input orientation parameterization.
    ! conv_string: String describing input orientation convention.
    ! ori_string: Temp string used to split the 'descriptor:convention' line.
    ! iarray: Substring array for internal read parsing.
    ! line: Input line on current record to be parsed.

    integer :: ierr, i, j, s
    integer :: nlines, nspace, nvals, elt, delim_pos
    character(len=50)  :: parm_string, conv_string, ori_string
    character(len=32)  :: iarray(16)
    character(len=256) :: line

    !---------------------------------------------------------------------------

    ! Check if the array is already allocated from reading the mesh file.
    if (allocated(elset_oris) .eqv. .true.) then
      deallocate (elset_oris)
    end if

    ! Read the first line and confirm it is the correct record.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$ElementOrientations')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &element ori.')
    end if

    ! Read in the section information string and prepare to parse.

    read (io, '(a)', iostat=ierr) line

    ! Trim the full string into two substring arrays split by the space delim.
    nspace = count(transfer(line, 'a', len_trim(line)) .eq. " ")
    nvals = nspace + 1
    iarray = ""

    ! Internal read of the line into the primary substring arrays.
    read (line, *) iarray(1:nvals)

    ! Internal read of the primary substring arrays to get the number of
    ! lines to set to parse and prepare to parse the orientation description.
    read (iarray(1), *) nlines
    read (iarray(2), '(a)') ori_string

    ! Find the position of the ":" character in ori_string as fortran
    ! refuses to allow the read() function to handle this automatically.
    delim_pos = index(ori_string, ":")

    ! Internal read to store the orientation parameterization and convention.
    read (ori_string(1:delim_pos - 1), *) parm_string
    read (ori_string(delim_pos + 1:len_trim(ori_string)), *) conv_string

    ! Set a per-element logical to true.
    eltoris_defined = .true.

    ! Parse the parameterization and confirm it is a valid option.

    s = 0

    read (parm_string, '(a)') &
        & mesh%orientation_parameterization

    ! Confirm that the rest of the string is a valid input.
    if (mesh%orientation_parameterization &
        & .eq. 'axis-angle') then
      s = 0

    else if (mesh%orientation_parameterization &
        & .eq. 'euler-bunge') then
      s = 0

    else if (mesh%orientation_parameterization &
        & .eq. 'euler-kocks') then
      s = 0

    else if (mesh%orientation_parameterization &
        & .eq. 'rodrigues') then
      s = 0

    else if (mesh%orientation_parameterization &
        & .eq. 'quaternion') then
      s = 0

    else ! Orientation parameterization defined is incorrect.
      s = 1
    end if

    if (s .eq. 1) then
      call par_quit("Error  :     > orientation_parameterization &
          & contains an error or unexpected input type.")
    end if

    ! Parse the convention and confirm it is a valid option.

    s = 0

    read (conv_string, '(a)') &
        & mesh%orientation_convention

    ! Confirm that the rest of the string is a valid input.
    if (mesh%orientation_convention &
        & .eq. 'active') then
      s = 0

    else if (mesh%orientation_convention &
        & .eq. 'passive') then
      s = 0

    else ! Orientation convention defined is incorrect.
      s = 1
    end if

    if (s .eq. 1) then
      call par_quit("Error  :     > orientation_convention &
          & contains an error or unexpected input type.")
    end if

    ! Check orientation parameterization to set column number.
    if ((mesh%orientation_parameterization &
        & .eq. 'axis-angle') .or. &
        & (mesh%orientation_parameterization &
        & .eq. 'quaternion')) then
      allocate (elset_oris(4, nlines))

    else ! All other parameterization have 3 values per orientation.
      allocate (elset_oris(3, nlines))
    end if

    ! Only quaternion or axis-angle.
    if (size(elset_oris, 1) .eq. 3 .or. size(elset_oris, 1) .eq. 4) then
      do i = 1, nlines
        read (io, *) elt, (elset_oris(j, i), j=1, size(elset_oris, 1))
      end do

    else
      call par_quit("Error  :     > elset_oris array not initialized&
          & correctly.")
    end if

    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndElementOrientations')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &ori.')
    end if

    return

  end subroutine read_eltorientations

  !===========================================================================

  subroutine read_groups(io, mesh, elt_elset, elset_ids_inv)

    ! Read elset and phase information from either the .msh file or an external
    ! simulation.phase file. If we are reading in phases from an external
    ! file we must deallocate the current storage array and start over.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Input unit for .msh file.

    integer, intent(in) :: io
    type(mesh_type), intent(inout):: mesh
    integer, intent(in) :: elt_elset(:)
    integer, allocatable, intent(in) :: elset_ids_inv(:)

    ! Locals:
    ! ierr: Value that confirms if a read() fails.
    ! i: Generic loop index.
    ! nlines: Read-in value for the number of lines that will be parsed.
    ! elset_id: Grain id - 1-indexed.
    ! elset_phase: Storage array to hold which phase is assigned to which elset.
    ! line: Input line on current record to be parsed.

    integer :: i, ierr, nlines, elset_id, phase_tmp
    integer, allocatable :: elset_phase(:)
    character(len=256)   :: line

    ! Notes:
    ! If the simulation is single-phase then we don't input any information
    ! on the elset phases anymore. Therefore we need to attempt to read in
    ! the next record on the mesh unit and see if it fails by eof. This only
    ! works as $Groups will always be the last section in the mesh file iff
    ! it is present and generated via Neper.

    ! Todo:
    ! It seems to be possible to define a number of elset/phase assignments
    ! that is out of bounds of the previously allocated array since we don't
    ! explicity enforce this.
    !----------------------------------------------------------------------

    ! Attempt to read the next line (may be eof, but unknown)
    read (io, '(a)', iostat=ierr) line

    if (ierr .lt. 0) then
      if (myid .eq. 0) then
        write (*, '(a)') 'Info   :   [i] No phase assignments &
            &found in `simulation.msh`.'
      end if

      return

    else if (ierr .eq. 0) then
      ! Check that the line read in previously was $Groups.
      if (line .ne. '$Groups') then
        call par_quit('Error  :     > Phase assignment section header &
            &should be `$Groups`.')
      end if

      ! Phase assignments are per elset so we need to do some mapping with
      ! temporary arrays.
      read (io, *) ! This skips the `elset` line
      read (io, *) nlines

      ! Check if the array is already allocated from the mesh read in
      if (allocated(elset_phase) .eqv. .true.) then
        deallocate (elset_phase)
      end if

      ! Allocate and fill an array with elset-based `elset/phase` pairs
      ! (1,:) is phase assign to elset id (row index)
      allocate (elset_phase(nlines))

      do i = 1, nlines
        read (io, *) elset_id, phase_tmp

        elset_phase(elset_ids_inv(elset_id)) = phase_tmp
      end do

      ! Testing elset_phase
      if (minval(elset_phase) .le. 0) then
        call par_quit('Error  :     > elset_phase is out of bounds.')
      end if

      ! Now, loop over elt_elset and reassign phase accordingly
      do i = 1, mesh%num_elts
        ! Assign phase per element based on elset value of elt_elset per element
        ! elt_elset(i) returns an elset value that is 0-indexed that will grab
        ! the corresponding elset row in elset_phase and, in turn, the phase
        mesh%elt_phase(i) = elset_phase(elt_elset(i))
      end do

      ! Testing phase
      if (minval(mesh%elt_phase) .le. 0) then
        call par_quit('Error  :     > Phase index is out of bounds.')
      end if

      ! Recording number of phases as the max phase id (note that this assumes that the
      ! phases are numbered contiguously from 1, which should always be the case.
      ! A test could be added here.
      mesh%num_phases = maxval(elset_phase)

    else
      ! If something goes wrong terminate the simulation.
      call par_quit('Error  :     > Phase assignment parsing failed.')
    end if

    ! Read the end of section footer.
    read (io, '(a)', iostat=ierr) line

    if ((ierr .lt. 0) .or. (line .ne. '$EndGroups')) then
      call par_quit('Error  :     > Parse error attempting to read in &
          &groups.')
    end if

    return

  end subroutine read_groups

  subroutine read_mesh_init_eltoris(mesh, eltoris_defined, elt_elset, elset_oris, nb_elsets)

    ! Assign ori to the elements and compute rotation matrices

    !---------------------------------------------------------------------------

    ! Arguments:
    ! io: Fortran unit number
    !   orientation, initialized to identity in this routine

    !integer :: io
    type(mesh_type), intent(inout):: mesh
    logical :: eltoris_defined
    integer, intent(in):: elt_elset(:)
    real(rk), intent(in) :: elset_oris (:, :)
    integer, intent(in) :: nb_elsets

    ! Locals:

    integer  :: my_phase(elt_sub:elt_sup)
    integer  :: my_elt_elset(elt_sub:elt_sup)

    real(rk) :: angle(elt_sub:elt_sup)
    real(rk) :: axis(3, elt_sub:elt_sup)
    real(rk) :: phi1(elt_sub:elt_sup)
    real(rk) :: phi(elt_sub:elt_sup)
    real(rk) :: phi2(elt_sub:elt_sup)
    real(rk) :: psi(elt_sub:elt_sup)
    real(rk) :: the(elt_sub:elt_sup)
    real(rk) :: rods(3, elt_sub:elt_sup)
    real(rk) :: quat(4, elt_sub:elt_sup)

    integer  :: i, j, m

    !---------------------------------------------------------------------------

    my_phase = mesh%elt_phase(elt_sub:elt_sup)
    my_elt_elset = elt_elset(elt_sub:elt_sup)

    m = elt_sup - elt_sub + 1

    ! Check whether or not element ori are to be used. Assign
    ! accordingly either per-elset or per-element.
    if (eltoris_defined .eqv. .false.) then
      ! Assign the ori to each element from elsets

      do i = 1, nb_elsets
        if (mesh%orientation_parameterization &
            & .eq. 'axis-angle') then
          where (my_elt_elset .eq. i)
            axis(1, :) = elset_oris(1, i)
            axis(2, :) = elset_oris(2, i)
            axis(3, :) = elset_oris(3, i)
            angle(:) = elset_oris(4, i)
          end where

        else if (mesh%orientation_parameterization &
            & .eq. 'euler-bunge') then
          where (my_elt_elset .eq. i)
            phi1(:) = elset_oris(1, i)
            phi(:) = elset_oris(2, i)
            phi2(:) = elset_oris(3, i)
          end where

        else if (mesh%orientation_parameterization &
            & .eq. 'euler-kocks') then
          where (my_elt_elset .eq. i)
            psi(:) = elset_oris(1, i)
            the(:) = elset_oris(2, i)
            phi(:) = elset_oris(3, i)
          end where

        else if (mesh%orientation_parameterization &
            & .eq. 'rodrigues') then
          where (my_elt_elset .eq. i)
            rods(1, :) = elset_oris(1, i)
            rods(2, :) = elset_oris(2, i)
            rods(3, :) = elset_oris(3, i)
          end where

        else if (mesh%orientation_parameterization &
            & .eq. 'quaternion') then
          where (my_elt_elset .eq. i)
            quat(1, :) = elset_oris(1, i)
            quat(2, :) = elset_oris(2, i)
            quat(3, :) = elset_oris(3, i)
            quat(4, :) = elset_oris(4, i)
          end where
        end if
      end do

    else if (eltoris_defined .eqv. .true.) then
      ! Assign the ori to each element directly

      do i = 1, mesh%num_elts
        if (mesh%orientation_parameterization &
            & .eq. 'axis-angle') then
          if ((i .ge. elt_sub) .and. (i .le. elt_sup)) then
            axis(1, i) = elset_oris(1, i)
            axis(2, i) = elset_oris(2, i)
            axis(3, i) = elset_oris(3, i)
            angle(i) = elset_oris(4, i)
          end if

        else if (mesh%orientation_parameterization &
            & .eq. 'euler-bunge') then
          if ((i .ge. elt_sub) .and. (i .le. elt_sup)) then
            phi1(i) = elset_oris(1, i)
            phi(i) = elset_oris(2, i)
            phi2(i) = elset_oris(3, i)
          end if

        else if (mesh%orientation_parameterization &
            & .eq. 'euler-kocks') then
          if ((i .ge. elt_sub) .and. (i .le. elt_sup)) then
            psi(i) = elset_oris(1, i)
            the(i) = elset_oris(2, i)
            phi(i) = elset_oris(3, i)
          end if

        else if (mesh%orientation_parameterization &
            & .eq. 'rodrigues') then
          if ((i .ge. elt_sub) .and. (i .le. elt_sup)) then
            rods(1, i) = elset_oris(1, i)
            rods(2, i) = elset_oris(2, i)
            rods(3, i) = elset_oris(3, i)
          end if

        else if (mesh%orientation_parameterization &
            & .eq. 'quaternion') then
          if ((i .ge. elt_sub) .and. (i .le. elt_sup)) then
            quat(1, i) = elset_oris(1, i)
            quat(2, i) = elset_oris(2, i)
            quat(3, i) = elset_oris(3, i)
            quat(4, i) = elset_oris(4, i)
          end if
        end if
      end do
    end if

    ! Find initial orientation matrix, c0_angs

    ! Determine parameterization from orientation options

    if (mesh%orientation_parameterization .eq. &
        & 'axis-angle') then
      call axis_angle_to_rotmat(axis, angle, mesh%ori)

    else if (mesh%orientation_parameterization .eq. &
        & 'euler-bunge') then
      call euler_bunge_to_rotmat(phi1, phi, phi2, &
          & mesh%ori)

    else if (mesh%orientation_parameterization .eq. &
        & 'euler-kocks') then
      call euler_kocks_to_rotmat(psi, the, phi, &
          & mesh%ori)

    else if (mesh%orientation_parameterization .eq. &
        & 'rodrigues') then
      call rodrigues_to_rotmat(rods, mesh%ori)

    else if (mesh%orientation_parameterization .eq. &
        & 'quaternion') then
      call quat_to_rotmat(quat, mesh%ori)
    end if

    ! Determine passive (c2s) or active (s2c) from orientation options

    if (mesh%orientation_convention .eq. &
        & 'passive') then
      ! Don't do anything - FEPX assumes passive convention!!

    else if (mesh%orientation_convention .eq. &
        & 'active') then
      do j = elt_sub, elt_sup
        mesh%ori(:, :, j) = transpose(mesh%ori(:, :, j))
      end do
    end if

    ! Initialize hardnesses

    return

  end subroutine read_mesh_init_eltoris

end module read_input_msh_mod2
