! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Distribution of the mesh to the ranks: root_rank reads the mesh files (mesh_file_mod), decides
!> the rank of every element and the owner of every node, and sends each rank its part; the
!> result on each rank is a mesh_t (mesh_mod).
!>
!> Source : FEPX (ROADMAP step 2.2; D5, D32, D108). Contains no numerical method.
!>
!> Policy:
!>   - Partitioning (D5, D108): the mesh's partition table gives the ranks (element partition p ->
!>     rank p - 1, node partition p -> owner rank p - 1); the table need not be sequential. The
!>     mesh must have as many partitions as there are ranks, a fatal error otherwise, except on
!>     1 rank, which takes any mesh whole (the warning for a partitioned mesh is written by
!>     fepx's setup summary, from mesh_t%num_parts_file).
!>   - Each rank receives its elements (with nodes, elset, phase, rot_c2s and per-element values),
!>     its owned nodes and its ghost nodes (nodes of its elements owned by other ranks, with their
!>     owners) with coordinates (mesh_mod layout). Small whole-mesh data (elsets, phases,
!>     orientation descriptor, nsets, fasets, periodicity, per-elset values) are broadcast;
!>     nsets and fasets are then restricted to the rank's local nodes and elements (the nsets are
!>     also kept whole, nsets_glob, step 8.3a).
!>   - Only root_rank opens the files (D32); a fatal error there ends every rank (MPI abort
!>     handler, step 2.1).
!>   - The distributed mesh is then localized (mesh_local_mod::mesh_localize: local numbering
!>     and communication graph), so mesh_load returns a mesh ready for use.
!>
!> Public entities:
!>   mesh_load        read the mesh files on root_rank and distribute (collective)
!>   mesh_distribute  distribute a parsed mesh held by root_rank (collective)
module mesh_distribute_mod

  use kinds_mod, only: rk, ik
  use error_mod, only: fatal_error, internal_error
  use parallel_mod, only: num_ranks, my_rank, root_rank, par_bcast, par_scatterv
  use mesh_file_mod, only: mesh_file_t, mesh_file_load, mesh_nset_t, mesh_faset_t, &
    & mesh_values_t, num_value_sections, value_section_is_elset
  use mesh_mod, only: mesh_t, find_sorted
  use mesh_local_mod, only: mesh_localize

  implicit none

  private

  public :: mesh_load, mesh_distribute

contains

  ! ==========================================================================================
  ! Public routines
  ! ==========================================================================================

  !> Read the mesh file (and auxiliary files) on root_rank and distribute it. Collective.
  !>
  !> Arguments:
  !>   mesh          (out)           type(mesh_t)     this rank's mesh
  !>   mesh_name     (in)            character        mesh file name (used on root_rank)
  !>   aux_names     (in, optional)  character(*)(:)  auxiliary file names (root_rank)
  !> Locals:
  !>   file_mesh  type(mesh_file_t)  whole mesh (root_rank only)
  subroutine mesh_load(mesh, mesh_name, aux_names)

    type(mesh_t), intent(out) :: mesh
    character(len=*), intent(in) :: mesh_name
    character(len=*), intent(in), optional :: aux_names(:)

    type(mesh_file_t) :: file_mesh

    ! Parse on the root only (D32), then distribute.
    if (my_rank == root_rank) call mesh_file_load(file_mesh, mesh_name, aux_names)
    call mesh_distribute(file_mesh, mesh)

  end subroutine mesh_load

  !> Distribute a parsed mesh held by root_rank (module policy). Collective.
  !>
  !> Arguments:
  !>   file_mesh     (in)   type(mesh_file_t)  whole mesh, finalized (significant on root_rank)
  !>   mesh          (out)  type(mesh_t)       this rank's mesh
  !> Locals:
  !>   elt_rank     integer(ik)(:)  rank of each element (root; empty elsewhere)
  !>   node_rank    integer(ik)(:)  owner rank of each node (root; empty elsewhere)
  !>   is_root      logical         whether this is root_rank
  subroutine mesh_distribute(file_mesh, mesh)

    type(mesh_file_t), intent(in) :: file_mesh
    type(mesh_t), intent(out) :: mesh

    integer(ik), allocatable :: elt_rank(:)
    integer(ik), allocatable :: node_rank(:)
    logical :: is_root

    is_root = (my_rank == root_rank)

    ! Root: rank of every element and owner of every node (empty on the other ranks, so that
    ! the arrays are defined everywhere: GCC's -O3 otherwise warns that their bounds may be
    ! used uninitialized, step 9.7).
    if (is_root) then
      call choose_ranks(file_mesh, elt_rank, node_rank)
    else
      allocate (elt_rank(0), node_rank(0))
    end if
    ! The mesh file's partition count (the setup summary, step 9.7).
    if (is_root) mesh%num_parts_file = file_mesh%num_parts
    call par_bcast(mesh%num_parts_file)

    ! Whole-mesh data, then each rank's part.
    call broadcast_global(file_mesh, mesh)
    call scatter_parts(file_mesh, elt_rank, node_rank, mesh)
    call broadcast_sets(file_mesh, mesh)

    ! Local numbering and communication graph (step 2.3).
    call mesh_localize(mesh)

  end subroutine mesh_distribute

  ! ==========================================================================================
  ! Partitioning (root_rank)
  ! ==========================================================================================

  !> The rank of every element and the owner of every node (module policy).
  !>
  !> Arguments:
  !>   file_mesh  (in)   type(mesh_file_t)  whole mesh
  !>   elt_rank   (out)  integer(ik)(:)     rank of each element
  !>   node_rank  (out)  integer(ik)(:)     owner rank of each node
  !> Locals:
  !>   parts  character(:)  the mesh's partitions, for the error message
  subroutine choose_ranks(file_mesh, elt_rank, node_rank)

    type(mesh_file_t), intent(in) :: file_mesh
    integer(ik), allocatable, intent(out) :: elt_rank(:)
    integer(ik), allocatable, intent(out) :: node_rank(:)

    character(len=:), allocatable :: parts

    allocate (elt_rank(file_mesh%num_elts), node_rank(file_mesh%num_nodes))

    if (num_ranks == 1_ik) then
      ! One rank takes any mesh whole.
      elt_rank = 0_ik
      node_rank = 0_ik
    else if (file_mesh%num_parts == num_ranks) then
      ! The partition table.
      elt_rank = file_mesh%elt_part - 1_ik
      node_rank = file_mesh%node_part - 1_ik
    else
      ! Another number of partitions (D108).
      if (file_mesh%num_parts == 0_ik) then
        parts = 'The mesh is not partitioned'
      else
        parts = 'The mesh has ' // int_text(file_mesh%num_parts) // ' partition(s)'
      end if
      call fatal_error('mesh_distribute', parts // ', but ' // int_text(int(num_ranks, ik)) // &
        & ' ranks were started: partition it for ' // int_text(int(num_ranks, ik)) // &
        & ' ranks (neper -M -loadmesh simulation.msh -part ' // &
        & int_text(int(num_ranks, ik)) // ' -o simulation) or start as many ranks as ' // &
        & 'it has partitions')
    end if

  end subroutine choose_ranks

  ! ==========================================================================================
  ! Communication
  ! ==========================================================================================

  !> Broadcast the whole-mesh scalars and small arrays (sizes, elsets, phases, orientation
  !> descriptor, domain, per-elset value sections).
  !>
  !> Arguments:
  !>   file_mesh  (in)     type(mesh_file_t)  whole mesh (significant on root_rank)
  !>   mesh       (inout)  type(mesh_t)       this rank's mesh
  !> Locals:
  !>   is_root  logical  whether this is root_rank
  !>   i_sec    integer  value-section index
  subroutine broadcast_global(file_mesh, mesh)

    type(mesh_file_t), intent(in) :: file_mesh
    type(mesh_t), intent(inout) :: mesh

    logical :: is_root
    integer :: i_sec

    is_root = (my_rank == root_rank)

    ! Sizes and element type.
    if (is_root) then
      mesh%num_nodes_glob = file_mesh%num_nodes
      mesh%num_elts_glob = file_mesh%num_elts
      mesh%elt_type = file_mesh%elt_type
      mesh%num_elt_nodes = file_mesh%num_elt_nodes
      mesh%num_elsets = file_mesh%num_elsets
      mesh%num_phases = file_mesh%num_phases
      mesh%has_ori = file_mesh%has_ori
      mesh%ori_from_elts = file_mesh%ori_from_elts
      mesh%ori_descriptor = file_mesh%ori_descriptor
      mesh%ori_convention = file_mesh%ori_convention
      mesh%domain = ''
      if (allocated(file_mesh%domain)) mesh%domain = file_mesh%domain
      mesh%elset_id = file_mesh%elset_id
      mesh%elset_phase = file_mesh%elset_phase
    end if
    call par_bcast(mesh%num_nodes_glob)
    call par_bcast(mesh%num_elts_glob)
    call par_bcast(mesh%elt_type)
    call par_bcast(mesh%num_elt_nodes)
    call par_bcast(mesh%num_elsets)
    call par_bcast(mesh%num_phases)
    call par_bcast(mesh%has_ori)
    call par_bcast(mesh%ori_from_elts)
    call par_bcast(mesh%ori_descriptor)
    call par_bcast(mesh%ori_convention)
    call bcast_string(mesh%domain)
    call bcast_int_array(mesh%elset_id)
    call bcast_int_array(mesh%elset_phase)

    ! Per-elset value sections (replicated); per-element ones: header data only here.
    do i_sec = 1, num_value_sections
      if (is_root) then
        mesh%values(i_sec)%defined = file_mesh%values(i_sec)%defined
        if (mesh%values(i_sec)%defined) then
          mesh%values(i_sec)%source = file_mesh%values(i_sec)%source
          mesh%values(i_sec)%num_values_hdr = file_mesh%values(i_sec)%num_values_hdr
          if (value_section_is_elset(i_sec)) then
            mesh%values(i_sec)%num_values = file_mesh%values(i_sec)%num_values
            mesh%values(i_sec)%values = file_mesh%values(i_sec)%values
          end if
        end if
      end if
      call par_bcast(mesh%values(i_sec)%defined)
      if (.not. mesh%values(i_sec)%defined) cycle
      call bcast_string(mesh%values(i_sec)%source)
      call par_bcast(mesh%values(i_sec)%num_values_hdr)
      if (value_section_is_elset(i_sec)) then
        call bcast_int_array(mesh%values(i_sec)%num_values)
        call bcast_real_matrix(mesh%values(i_sec)%values)
      end if
    end do

  end subroutine broadcast_global

  !> Send each rank its elements and nodes (module policy and mesh_mod layout).
  !>
  !> Arguments:
  !>   file_mesh  (in)     type(mesh_file_t)  whole mesh (significant on root_rank)
  !>   elt_rank   (in)     integer(ik)(:)     rank of each element (root)
  !>   node_rank  (in)     integer(ik)(:)     owner rank of each node (root)
  !>   mesh       (inout)  type(mesh_t)       this rank's mesh
  !> Locals (root unless stated):
  !>   rank_elts    integer(ik)(:)  element numbers grouped by rank (increasing within a rank)
  !>   count_elts   integer(ik)(:)  elements of each rank
  !>   rank_nodes   integer(ik)(:)  local node lists of all ranks (owned then ghosts)
  !>   count_nodes  integer(ik)(:)  local nodes of each rank
  !>   ghost_owner  integer(ik)(:)  owners of the ghost nodes of all ranks
  !>   count_ghost  integer(ik)(:)  ghost nodes of each rank
  !>   counts_2     integer(ik)(:)  2 per rank (for the owned/ghost counts)
  !>   own_ghost    integer(ik)(:)  owned and ghost counts of each rank, interleaved
  !>   recv_int     integer(ik)(:)  received integers (all ranks)
  !>   recv_real    real(rk)(:)     received reals (all ranks)
  !>   send_int     integer(ik)(:)  integer send buffer (root)
  !>   send_real    real(rk)(:)     real send buffer (root)
  !>   i_sec        integer         value-section index
  !>   width        integer(ik)     values per element of a value section
  !>   is_root      logical         whether this is root_rank
  subroutine scatter_parts(file_mesh, elt_rank, node_rank, mesh)

    type(mesh_file_t), intent(in) :: file_mesh
    integer(ik), allocatable, intent(in) :: elt_rank(:)
    integer(ik), allocatable, intent(in) :: node_rank(:)
    type(mesh_t), intent(inout) :: mesh

    integer(ik), allocatable :: rank_elts(:)
    integer(ik), allocatable :: count_elts(:)
    integer(ik), allocatable :: rank_nodes(:)
    integer(ik), allocatable :: count_nodes(:)
    integer(ik), allocatable :: ghost_owner(:)
    integer(ik), allocatable :: count_ghost(:)
    integer(ik), allocatable :: counts_2(:)
    integer(ik), allocatable :: own_ghost(:)
    integer(ik), allocatable :: recv_int(:)
    real(rk), allocatable :: recv_real(:)
    integer(ik), allocatable :: send_int(:)
    real(rk), allocatable :: send_real(:)
    integer :: i_sec
    integer(ik) :: width
    logical :: is_root

    is_root = (my_rank == root_rank)

    ! Root: element and node lists of every rank.
    if (is_root) then
      call group_elements(elt_rank, rank_elts, count_elts)
      call build_node_lists(file_mesh, node_rank, rank_elts, count_elts, rank_nodes, &
        & count_nodes, ghost_owner, count_ghost, own_ghost)
      allocate (counts_2(num_ranks))
      counts_2 = 2_ik
    else
      allocate (rank_elts(0), count_elts(0), rank_nodes(0), count_nodes(0), ghost_owner(0), &
        & count_ghost(0), counts_2(0), own_ghost(0))
    end if

    ! Elements: global numbers, nodes, elsets, phases.
    call par_scatterv(rank_elts, count_elts, mesh%elt_id_glob)
    mesh%num_elts = int(size(mesh%elt_id_glob), ik)
    if (is_root) then
      send_int = reshape(file_mesh%elt_nodes(:, rank_elts), [size(rank_elts)* &
        & file_mesh%num_elt_nodes])
    else
      allocate (send_int(0))
    end if
    call par_scatterv(send_int, count_elts*mesh_nodes_per_elt(file_mesh, is_root), recv_int)
    mesh%elt_nodes_glob = reshape(recv_int, [mesh%num_elt_nodes, mesh%num_elts])
    if (is_root) send_int = file_mesh%elt_elset(rank_elts)
    call par_scatterv(send_int, count_elts, mesh%elt_elset)
    if (is_root) send_int = file_mesh%elt_phase(rank_elts)
    call par_scatterv(send_int, count_elts, mesh%elt_phase)

    ! Orientations (zero when none were given).
    if (mesh%has_ori) then
      if (is_root) then
        send_real = reshape(file_mesh%elt_rot_c2s(:, :, rank_elts), [9*size(rank_elts)])
      else
        allocate (send_real(0))
      end if
      call par_scatterv(send_real, 9_ik*count_elts, recv_real)
      mesh%rot_c2s = reshape(recv_real, [3, 3, int(mesh%num_elts)])
    else
      allocate (mesh%rot_c2s(3, 3, mesh%num_elts))
      mesh%rot_c2s = 0.0_rk
    end if

    ! Nodes: owned and ghost counts, global ids, ghost owners, coordinates.
    call par_scatterv(own_ghost, counts_2, recv_int)
    mesh%num_nodes_owned = recv_int(1)
    mesh%num_nodes_ghost = recv_int(2)
    mesh%num_nodes_loc = mesh%num_nodes_owned + mesh%num_nodes_ghost
    call par_scatterv(rank_nodes, count_nodes, mesh%node_id_glob)
    call par_scatterv(ghost_owner, count_ghost, mesh%ghost_owner)
    if (is_root) then
      send_real = reshape(file_mesh%coo(:, rank_nodes), [3*size(rank_nodes)])
    else
      if (allocated(send_real)) deallocate (send_real)
      allocate (send_real(0))
    end if
    call par_scatterv(send_real, 3_ik*count_nodes, recv_real)
    mesh%coo_s = reshape(recv_real, [3, int(size(mesh%node_id_glob))])

    ! Per-element value sections: values of the local elements.
    do i_sec = 1, num_value_sections
      if (value_section_is_elset(i_sec) .or. .not. mesh%values(i_sec)%defined) cycle
      width = 0_ik
      if (is_root) width = int(size(file_mesh%values(i_sec)%values, 1), ik)
      call par_bcast(width)
      if (is_root) send_int = file_mesh%values(i_sec)%num_values(rank_elts)
      call par_scatterv(send_int, count_elts, mesh%values(i_sec)%num_values)
      if (is_root) then
        send_real = reshape(file_mesh%values(i_sec)%values(:, rank_elts), &
          & [int(width)*size(rank_elts)])
      end if
      call par_scatterv(send_real, width*count_elts, recv_real)
      mesh%values(i_sec)%values = reshape(recv_real, [int(width), int(mesh%num_elts)])
    end do

  end subroutine scatter_parts

  !> Broadcast the node sets, face sets and periodicity, and restrict the sets to this rank's
  !> local nodes and elements; the node sets are also kept whole (mesh_t%nsets_glob, step 8.3a).
  !>
  !> Arguments:
  !>   file_mesh  (in)     type(mesh_file_t)  whole mesh (significant on root_rank)
  !>   mesh       (inout)  type(mesh_t)       this rank's mesh (elements and nodes received)
  !> Locals:
  !>   is_root    logical          whether this is root_rank
  !>   num_sets   integer(ik)      number of nsets or fasets
  !>   i_set      integer(ik)      set index
  !>   label      character(:)     set label
  !>   nodes      integer(ik)(:)   nset nodes (global)
  !>   keep       logical(:)       which entries are local
  !>   elts       integer(ik)(:)   faset element numbers (global)
  !>   faces      integer(ik)(:,:) faset face nodes
  !>   loc_elt    integer(ik)(:)   local index of each faset element (0: not local)
  !>   num_face_nodes  integer(ik)  nodes per face
  !>   i          integer(ik)      entry index
  subroutine broadcast_sets(file_mesh, mesh)

    type(mesh_file_t), intent(in) :: file_mesh
    type(mesh_t), intent(inout) :: mesh

    logical :: is_root
    integer(ik) :: num_sets
    integer(ik) :: i_set
    character(len=:), allocatable :: label
    integer(ik), allocatable :: nodes(:)
    logical, allocatable :: keep(:)
    integer(ik), allocatable :: elts(:)
    integer(ik), allocatable :: faces(:, :)
    integer(ik), allocatable :: loc_elt(:)
    integer(ik) :: num_face_nodes
    integer(ik) :: i

    is_root = (my_rank == root_rank)

    ! Node sets: broadcast each, keep its local nodes (file order), and the whole set (global
    ! ids, file order; for the multi-point constraints, step 8.3a).
    num_sets = 0_ik
    if (is_root) num_sets = int(size(file_mesh%nsets), ik)
    call par_bcast(num_sets)
    allocate (mesh%nsets(num_sets), mesh%nsets_glob(num_sets))
    do i_set = 1_ik, num_sets
      if (is_root) then
        label = file_mesh%nsets(i_set)%label
        nodes = file_mesh%nsets(i_set)%nodes
      end if
      call bcast_string(label)
      call bcast_int_array(nodes)
      keep = [(is_local_node(mesh, nodes(i)), i = 1_ik, int(size(nodes), ik))]
      mesh%nsets(i_set)%label = label
      mesh%nsets(i_set)%nodes = pack(nodes, keep)
      mesh%nsets_glob(i_set)%label = label
      mesh%nsets_glob(i_set)%nodes = nodes
    end do

    ! Face sets: broadcast each, keep the faces of local elements (local element index).
    num_sets = 0_ik
    if (is_root) num_sets = int(size(file_mesh%fasets), ik)
    call par_bcast(num_sets)
    allocate (mesh%fasets(num_sets))
    do i_set = 1_ik, num_sets
      if (is_root) then
        label = file_mesh%fasets(i_set)%label
        num_face_nodes = file_mesh%fasets(i_set)%num_face_nodes
        elts = file_mesh%fasets(i_set)%elts
        faces = file_mesh%fasets(i_set)%nodes
      end if
      call bcast_string(label)
      call par_bcast(num_face_nodes)
      call bcast_int_array(elts)
      call bcast_int_matrix(faces)
      loc_elt = [(find_sorted(mesh%elt_id_glob, elts(i)), i = 1_ik, int(size(elts), ik))]
      keep = loc_elt > 0_ik
      mesh%fasets(i_set)%label = label
      mesh%fasets(i_set)%num_face_nodes = num_face_nodes
      mesh%fasets(i_set)%elts = pack(loc_elt, keep)
      mesh%fasets(i_set)%nodes = faces(:, pack([(i, i = 1_ik, int(size(elts), ik))], keep))
    end do

    ! Periodicity: replicated.
    if (is_root) then
      mesh%num_periodicity = file_mesh%num_periodicity
      mesh%per_secondary = file_mesh%per_secondary
      mesh%per_primary = file_mesh%per_primary
      mesh%per_vect = file_mesh%per_vect
    end if
    call par_bcast(mesh%num_periodicity)
    call bcast_int_array(mesh%per_secondary)
    call bcast_int_array(mesh%per_primary)
    call bcast_int_matrix(mesh%per_vect)

  end subroutine broadcast_sets

  ! ==========================================================================================
  ! Root-side list building
  ! ==========================================================================================

  !> Group the element numbers by rank (counting sort; increasing within a rank).
  !>
  !> Arguments:
  !>   elt_rank    (in)   integer(ik)(:)  rank of each element
  !>   rank_elts   (out)  integer(ik)(:), allocatable  element numbers grouped by rank
  !>   count_elts  (out)  integer(ik)(:), allocatable  elements of each rank (rank r at r + 1)
  !> Locals:
  !>   next   integer(ik)(:)  next free position of each rank's group
  !>   i_elt  integer(ik)     element number (also rank index in the prefix sum)
  subroutine group_elements(elt_rank, rank_elts, count_elts)

    integer(ik), intent(in) :: elt_rank(:)
    integer(ik), allocatable, intent(out) :: rank_elts(:)
    integer(ik), allocatable, intent(out) :: count_elts(:)

    integer(ik), allocatable :: next(:)
    integer(ik) :: i_elt

    ! Count the elements of each rank (count_elts(r + 1) for rank r).
    allocate (count_elts(num_ranks), rank_elts(size(elt_rank)))
    count_elts = 0_ik
    do i_elt = 1_ik, int(size(elt_rank), ik)
      count_elts(elt_rank(i_elt) + 1_ik) = count_elts(elt_rank(i_elt) + 1_ik) + 1_ik
    end do

    ! First position of each rank's group, then place the elements in increasing order.
    allocate (next(num_ranks))
    next(1) = 1_ik
    do i_elt = 2_ik, num_ranks
      next(i_elt) = next(i_elt - 1) + count_elts(i_elt - 1)
    end do
    do i_elt = 1_ik, int(size(elt_rank), ik)
      rank_elts(next(elt_rank(i_elt) + 1_ik)) = i_elt
      next(elt_rank(i_elt) + 1_ik) = next(elt_rank(i_elt) + 1_ik) + 1_ik
    end do

  end subroutine group_elements

  !> Local node lists of every rank: owned nodes (increasing), then ghost nodes (nodes of the
  !> rank's elements owned by other ranks, increasing), with the ghosts' owners.
  !>
  !> Arguments:
  !>   file_mesh    (in)   type(mesh_file_t)  whole mesh
  !>   node_rank    (in)   integer(ik)(:)     owner rank of each node
  !>   rank_elts    (in)   integer(ik)(:)     element numbers grouped by rank
  !>   count_elts   (in)   integer(ik)(:)     elements of each rank
  !>   rank_nodes   (out)  integer(ik)(:), allocatable  local nodes of all ranks
  !>   count_nodes  (out)  integer(ik)(:), allocatable  local nodes of each rank
  !>   ghost_owner  (out)  integer(ik)(:), allocatable  owners of the ghosts of all ranks
  !>   count_ghost  (out)  integer(ik)(:), allocatable  ghosts of each rank
  !>   own_ghost    (out)  integer(ik)(:), allocatable  owned and ghost counts, interleaved
  !> Locals:
  !>   mark      integer(ik)(:)  last rank whose elements touched each node
  !>   i_rank    integer(ik)     rank
  !>   first     integer(ik)     first position of the rank's elements in rank_elts
  !>   i_elt     integer(ik)     position in rank_elts
  !>   i_node    integer(ik)     node id
  !>   owned     integer(ik)(:)  owned nodes of the rank
  !>   ghosts    integer(ik)(:)  ghost nodes of the rank
  subroutine build_node_lists(file_mesh, node_rank, rank_elts, count_elts, rank_nodes, &
    & count_nodes, ghost_owner, count_ghost, own_ghost)

    type(mesh_file_t), intent(in) :: file_mesh
    integer(ik), intent(in) :: node_rank(:)
    integer(ik), intent(in) :: rank_elts(:)
    integer(ik), intent(in) :: count_elts(:)
    integer(ik), allocatable, intent(out) :: rank_nodes(:)
    integer(ik), allocatable, intent(out) :: count_nodes(:)
    integer(ik), allocatable, intent(out) :: ghost_owner(:)
    integer(ik), allocatable, intent(out) :: count_ghost(:)
    integer(ik), allocatable, intent(out) :: own_ghost(:)

    integer(ik), allocatable :: mark(:)
    integer(ik) :: i_rank
    integer(ik) :: first
    integer(ik) :: i_elt
    integer(ik) :: i_node
    integer(ik), allocatable :: owned(:)
    integer(ik), allocatable :: ghosts(:)

    allocate (mark(file_mesh%num_nodes), rank_nodes(0), ghost_owner(0))
    allocate (count_nodes(num_ranks), count_ghost(num_ranks), own_ghost(2*num_ranks))
    mark = -1_ik
    first = 1_ik
    do i_rank = 0_ik, num_ranks - 1_ik
      ! Mark the nodes of the rank's elements.
      do i_elt = first, first + count_elts(i_rank + 1) - 1_ik
        mark(file_mesh%elt_nodes(:, rank_elts(i_elt))) = i_rank
      end do
      first = first + count_elts(i_rank + 1)

      ! Owned nodes, then marked nodes owned elsewhere (both increasing).
      owned = pack([(i_node, i_node = 1_ik, file_mesh%num_nodes)], node_rank == i_rank)
      ghosts = pack([(i_node, i_node = 1_ik, file_mesh%num_nodes)], &
        & mark == i_rank .and. node_rank /= i_rank)
      rank_nodes = [rank_nodes, owned, ghosts]
      ghost_owner = [ghost_owner, node_rank(ghosts)]
      count_nodes(i_rank + 1) = int(size(owned) + size(ghosts), ik)
      count_ghost(i_rank + 1) = int(size(ghosts), ik)
      own_ghost(2*i_rank + 1:2*i_rank + 2) = [int(size(owned), ik), int(size(ghosts), ik)]
    end do

  end subroutine build_node_lists

  ! ==========================================================================================
  ! Helpers
  ! ==========================================================================================

  !> Nodes per element on the root (for scatter counts; 0 elsewhere, where counts are ignored).
  !>
  !> Arguments:
  !>   file_mesh  (in)  type(mesh_file_t)  whole mesh
  !>   is_root    (in)  logical            whether this is root_rank
  !>   num        (out, result) integer(ik)
  pure function mesh_nodes_per_elt(file_mesh, is_root) result(num)

    type(mesh_file_t), intent(in) :: file_mesh
    logical, intent(in) :: is_root
    integer(ik) :: num

    num = 0_ik
    if (is_root) num = file_mesh%num_elt_nodes

  end function mesh_nodes_per_elt

  !> Whether a global node is local (owned or ghost) on this rank (binary search in the two
  !> increasing groups of node_id_glob).
  !>
  !> Arguments:
  !>   mesh     (in)  type(mesh_t)  this rank's mesh
  !>   node_id  (in)  integer(ik)   global node id
  !>   found    (out, result) logical
  pure function is_local_node(mesh, node_id) result(found)

    type(mesh_t), intent(in) :: mesh
    integer(ik), intent(in) :: node_id
    logical :: found

    found = find_sorted(mesh%node_id_glob(1:mesh%num_nodes_owned), node_id) > 0_ik
    if (.not. found) found = find_sorted(mesh%node_id_glob(mesh%num_nodes_owned + 1:), &
      & node_id) > 0_ik

  end function is_local_node

  !> Broadcast a string from root_rank (length, then characters).
  !>
  !> Arguments:
  !>   text  (inout)  character(:), allocatable  significant on root_rank; set on every rank
  !> Locals:
  !>   length  integer(ik)  string length
  subroutine bcast_string(text)

    character(len=:), allocatable, intent(inout) :: text

    integer(ik) :: length

    length = 0_ik
    if (my_rank == root_rank) length = int(len(text), ik)
    call par_bcast(length)
    if (my_rank /= root_rank) then
      if (allocated(text)) deallocate (text)
      allocate (character(len=length) :: text)
    end if
    if (length > 0_ik) call par_bcast(text)

  end subroutine bcast_string

  !> Broadcast an integer array from root_rank (size, then values).
  !>
  !> Arguments:
  !>   values  (inout)  integer(ik)(:), allocatable  significant on root_rank; set on every rank
  !> Locals:
  !>   n  integer(ik)  size
  subroutine bcast_int_array(values)

    integer(ik), allocatable, intent(inout) :: values(:)

    integer(ik) :: n

    n = 0_ik
    if (my_rank == root_rank) n = int(size(values), ik)
    call par_bcast(n)
    if (my_rank /= root_rank) then
      if (allocated(values)) deallocate (values)
      allocate (values(n))
    end if
    if (n > 0_ik) call par_bcast(values)

  end subroutine bcast_int_array

  !> Broadcast an integer matrix from root_rank (shape, then values).
  !>
  !> Arguments:
  !>   values  (inout)  integer(ik)(:,:), allocatable  significant on root_rank
  !> Locals:
  !>   shape_  integer(ik)(2)   shape
  !>   flat    integer(ik)(:)   values in array element order
  subroutine bcast_int_matrix(values)

    integer(ik), allocatable, intent(inout) :: values(:, :)

    integer(ik) :: shape_(2)
    integer(ik), allocatable :: flat(:)

    shape_ = 0_ik
    if (my_rank == root_rank) then
      shape_ = int(shape(values), ik)
      flat = reshape(values, [size(values)])
    else
      allocate (flat(0))
    end if
    call par_bcast(shape_)
    if (my_rank /= root_rank) then
      deallocate (flat)
      allocate (flat(shape_(1)*shape_(2)))
    end if
    if (size(flat) > 0) call par_bcast(flat)
    values = reshape(flat, int(shape_))

  end subroutine bcast_int_matrix

  !> Broadcast a real matrix from root_rank (shape, then values).
  !>
  !> Arguments:
  !>   values  (inout)  real(rk)(:,:), allocatable  significant on root_rank
  !> Locals:
  !>   shape_  integer(ik)(2)  shape
  !>   flat    real(rk)(:)     values in array element order
  subroutine bcast_real_matrix(values)

    real(rk), allocatable, intent(inout) :: values(:, :)

    integer(ik) :: shape_(2)
    real(rk), allocatable :: flat(:)

    shape_ = 0_ik
    if (my_rank == root_rank) then
      shape_ = int(shape(values), ik)
      flat = reshape(values, [size(values)])
    else
      allocate (flat(0))
    end if
    call par_bcast(shape_)
    if (my_rank /= root_rank) then
      deallocate (flat)
      allocate (flat(shape_(1)*shape_(2)))
    end if
    if (size(flat) > 0) call par_bcast(flat)
    values = reshape(flat, int(shape_))

  end subroutine bcast_real_matrix

  !> An integer as text, without blanks.
  !>
  !> Arguments:
  !>   value  (in)  integer(ik)  the integer
  !>   text   (out, result) character(:), allocatable
  !> Locals:
  !>   buffer  character(24)  formatted value
  pure function int_text(value) result(text)

    integer(ik), intent(in) :: value
    character(len=:), allocatable :: text

    character(len=24) :: buffer

    write (buffer, '(i0)') value
    text = trim(buffer)

  end function int_text

end module mesh_distribute_mod
