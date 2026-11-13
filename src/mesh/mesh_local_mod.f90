! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Local numbering and communication graph of a distributed mesh (ROADMAP step 2.3).
!>
!> Source : FEPX (D5: owned + ghost nodes, local 1-based numbering, halo exchange). The graph
!>          replaces v2 gather_scatter_mod.f90::part_scatter_setup (the "trace" of v2's element
!>          gather/scatter) for the owned + ghost layout. Contains no numerical method.
!>
!> Policy (layout: mesh_mod):
!>   - The local node numbers are the positions in node_id_glob (owned nodes, then ghosts);
!>     elt_nodes_loc and the nodes of the nsets and fasets are converted to them.
!>   - Communication graph: the ghosts of this rank are grouped by owner (receive lists); each
!>     rank tells each owner which of its nodes it holds as ghosts (par_alltoall of the counts,
!>     par_alltoallv of the global ids, both in increasing global id), and the owner converts
!>     them to its local numbers (send lists). Neighbors are the ranks this rank sends to or
!>     receives from. The all-to-all exchanges cost O(num_ranks) per rank, once per mesh.
!>   - Inconsistencies (an element node that is not local, a requested node that is not owned)
!>     are fatal errors; they would be errors of the distribution.
!>
!> Public entities:
!>   mesh_localize  build the local numbering and the communication graph (collective)
module mesh_local_mod

  use kinds_mod, only: ik
  use error_mod, only: internal_error
  use parallel_mod, only: num_ranks, my_rank, par_alltoall, par_alltoallv
  use mesh_mod, only: mesh_t, node_loc_from_glob

  implicit none

  private

  public :: mesh_localize

contains

  !> Build the local numbering and the communication graph of a distributed mesh (module
  !> policy). Collective: every rank must call it.
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_t)  this rank's mesh (from mesh_distribute_mod)
  subroutine mesh_localize(mesh)

    type(mesh_t), intent(inout) :: mesh

    ! Local numbers in the connectivity and the sets, then the graph.
    call localize_connectivity(mesh)
    call localize_sets(mesh)
    call build_comm_graph(mesh)

  end subroutine mesh_localize

  !> Local node numbers of the element nodes, and the element dof table (elt_dofs_loc).
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_t)  this rank's mesh
  !> Locals:
  !>   i_elt   integer(ik)  local element
  !>   i_node  integer(ik)  node of the element
  subroutine localize_connectivity(mesh)

    type(mesh_t), intent(inout) :: mesh

    integer(ik) :: i_elt
    integer(ik) :: i_node

    allocate (mesh%elt_nodes_loc(mesh%num_elt_nodes, mesh%num_elts))
    do i_elt = 1_ik, mesh%num_elts
      do i_node = 1_ik, mesh%num_elt_nodes
        mesh%elt_nodes_loc(i_node, i_elt) = node_loc_from_glob(mesh, &
          & mesh%elt_nodes_glob(i_node, i_elt))
      end do
    end do

    ! Every node of a local element is local (owned or ghost).
    if (any(mesh%elt_nodes_loc == 0_ik)) call internal_error('mesh_localize', &
      & 'a node of a local element is neither owned nor a ghost on rank ' // int_text(my_rank))

    ! Local dofs of every element (node a, component k: 3 (a_loc - 1) + k).
    allocate (mesh%elt_dofs_loc(3*mesh%num_elt_nodes, mesh%num_elts))
    do i_elt = 1_ik, mesh%num_elts
      do i_node = 1_ik, mesh%num_elt_nodes
        mesh%elt_dofs_loc(3*i_node - 2:3*i_node, i_elt) = &
          & 3_ik*(mesh%elt_nodes_loc(i_node, i_elt) - 1_ik) + [1_ik, 2_ik, 3_ik]
      end do
    end do

  end subroutine localize_connectivity

  !> Local node numbers of the nsets and fasets (already restricted to the local nodes and
  !> elements by mesh_distribute_mod).
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_t)  this rank's mesh
  !> Locals:
  !>   i_set   integer      set index
  !>   i       integer(ik)  entry index
  !>   j       integer(ik)  face node index
  subroutine localize_sets(mesh)

    type(mesh_t), intent(inout) :: mesh

    integer :: i_set
    integer(ik) :: i
    integer(ik) :: j

    ! Node sets.
    do i_set = 1, size(mesh%nsets)
      do i = 1_ik, int(size(mesh%nsets(i_set)%nodes), ik)
        mesh%nsets(i_set)%nodes(i) = node_loc_from_glob(mesh, mesh%nsets(i_set)%nodes(i))
      end do
      if (any(mesh%nsets(i_set)%nodes == 0_ik)) call internal_error('mesh_localize', &
        & 'nset ' // mesh%nsets(i_set)%label // ' has a node that is not local')
    end do

    ! Face sets (the faces of local elements: all their nodes are local).
    do i_set = 1, size(mesh%fasets)
      do i = 1_ik, int(size(mesh%fasets(i_set)%nodes, 2), ik)
        do j = 1_ik, mesh%fasets(i_set)%num_face_nodes
          mesh%fasets(i_set)%nodes(j, i) = node_loc_from_glob(mesh, &
            & mesh%fasets(i_set)%nodes(j, i))
        end do
      end do
      if (any(mesh%fasets(i_set)%nodes == 0_ik)) call internal_error('mesh_localize', &
        & 'faset ' // mesh%fasets(i_set)%label // ' has a node that is not local')
    end do

  end subroutine localize_sets

  !> Build the communication graph (module policy).
  !>
  !> Arguments:
  !>   mesh  (inout)  type(mesh_t)  this rank's mesh (local numbering done)
  !> Locals:
  !>   req_counts    integer(ik)(num_ranks)  ghosts of this rank owned by each rank
  !>   need_counts   integer(ik)(num_ranks)  owned nodes of this rank held as ghosts by each rank
  !>   req_ids       integer(ik)(:)  ghost global ids grouped by owner (increasing in a group)
  !>   req_loc       integer(ik)(:)  the same ghosts as local numbers
  !>   need_ids      integer(ik)(:)  global ids requested from this rank, grouped by requester
  !>   next          integer(ik)(:)  next free position of each owner's group
  !>   i_ghost       integer(ik)     ghost index (1..num_nodes_ghost)
  !>   owner         integer(ik)     owner rank of a ghost
  !>   i_rank        integer(ik)     rank (0-based)
  !>   k             integer(ik)     neighbor index
  !>   i             integer(ik)     entry index
  !>   recv_pos      integer(ik)     position in req_loc of the current owner's group
  !>   send_pos      integer(ik)     position in need_ids of the current requester's group
  subroutine build_comm_graph(mesh)

    type(mesh_t), intent(inout) :: mesh

    integer(ik) :: req_counts(num_ranks)
    integer(ik) :: need_counts(num_ranks)
    integer(ik), allocatable :: req_ids(:)
    integer(ik), allocatable :: req_loc(:)
    integer(ik), allocatable :: need_ids(:)
    integer(ik), allocatable :: next(:)
    integer(ik) :: i_ghost
    integer(ik) :: owner
    integer(ik) :: i_rank
    integer(ik) :: k
    integer(ik) :: i
    integer(ik) :: recv_pos
    integer(ik) :: send_pos

    ! Ghosts per owner (an owner is another rank).
    req_counts = 0_ik
    do i_ghost = 1_ik, mesh%num_nodes_ghost
      owner = mesh%ghost_owner(i_ghost)
      if (owner < 0_ik .or. owner >= num_ranks .or. owner == my_rank) then
        call internal_error('mesh_localize', 'invalid owner of a ghost node on rank ' // &
          & int_text(my_rank))
      end if
      req_counts(owner + 1_ik) = req_counts(owner + 1_ik) + 1_ik
    end do

    ! Ghosts grouped by owner, keeping the increasing global order within each group.
    allocate (next(num_ranks), req_ids(mesh%num_nodes_ghost), req_loc(mesh%num_nodes_ghost))
    next(1) = 1_ik
    do i_rank = 2_ik, num_ranks
      next(i_rank) = next(i_rank - 1_ik) + req_counts(i_rank - 1_ik)
    end do
    do i_ghost = 1_ik, mesh%num_nodes_ghost
      owner = mesh%ghost_owner(i_ghost)
      req_ids(next(owner + 1_ik)) = mesh%node_id_glob(mesh%num_nodes_owned + i_ghost)
      req_loc(next(owner + 1_ik)) = mesh%num_nodes_owned + i_ghost
      next(owner + 1_ik) = next(owner + 1_ik) + 1_ik
    end do

    ! Tell each owner which of its nodes this rank needs.
    call par_alltoall(req_counts, need_counts)
    call par_alltoallv(req_ids, req_counts, need_ids, need_counts)

    ! Neighbors: ranks this rank receives from or sends to.
    mesh%comm%num_neighbors = int(count(req_counts > 0_ik .or. need_counts > 0_ik), ik)
    allocate (mesh%comm%neighbor(mesh%comm%num_neighbors))
    allocate (mesh%comm%recv_offset(mesh%comm%num_neighbors + 1_ik))
    allocate (mesh%comm%send_offset(mesh%comm%num_neighbors + 1_ik))
    allocate (mesh%comm%recv_nodes(mesh%num_nodes_ghost))
    allocate (mesh%comm%send_nodes(size(need_ids)))

    ! Receive and send lists of each neighbor, in rank order.
    k = 0_ik
    recv_pos = 1_ik
    send_pos = 1_ik
    mesh%comm%recv_offset(1) = 1_ik
    mesh%comm%send_offset(1) = 1_ik
    do i_rank = 0_ik, num_ranks - 1_ik
      if (req_counts(i_rank + 1_ik) == 0_ik .and. need_counts(i_rank + 1_ik) == 0_ik) cycle
      k = k + 1_ik
      mesh%comm%neighbor(k) = i_rank

      ! Ghosts received from this neighbor (local numbers).
      mesh%comm%recv_nodes(recv_pos:recv_pos + req_counts(i_rank + 1_ik) - 1_ik) = &
        & req_loc(recv_pos:recv_pos + req_counts(i_rank + 1_ik) - 1_ik)
      recv_pos = recv_pos + req_counts(i_rank + 1_ik)
      mesh%comm%recv_offset(k + 1_ik) = recv_pos

      ! Owned nodes sent to this neighbor: its requests, converted to local numbers.
      do i = send_pos, send_pos + need_counts(i_rank + 1_ik) - 1_ik
        mesh%comm%send_nodes(i) = node_loc_from_glob(mesh, need_ids(i))
        if (mesh%comm%send_nodes(i) < 1_ik .or. &
          & mesh%comm%send_nodes(i) > mesh%num_nodes_owned) then
          call internal_error('mesh_localize', 'rank ' // int_text(i_rank) // ' needs node ' // &
            & int_text(need_ids(i)) // ', which rank ' // int_text(my_rank) // ' does not own')
        end if
      end do
      send_pos = send_pos + need_counts(i_rank + 1_ik)
      mesh%comm%send_offset(k + 1_ik) = send_pos
    end do

  end subroutine build_comm_graph

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

end module mesh_local_mod
