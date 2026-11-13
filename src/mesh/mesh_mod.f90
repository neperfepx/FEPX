! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Rank-local mesh: the elements of this rank, its owned and ghost nodes, the mesh data it
!> needs, the local numbering and the communication graph (filled by mesh_distribute_mod and
!> mesh_local_mod).
!>
!> Source : FEPX (ROADMAP steps 2.2-2.3; D5: owned + ghost nodes, local 1-based numbering).
!>          The communication graph replaces v2 gather_scatter_mod.f90::part_scatter_setup (the
!>          "trace" of v2's element gather/scatter) for the owned + ghost layout of D5.
!>          Contains no numerical method.
!>
!> Layout:
!>   - Elements: local element i is global element elt_id_glob(i) (increasing).
!>   - Nodes: local node i (1..num_nodes_loc) is global node node_id_glob(i); the owned nodes
!>     come first (1..num_nodes_owned, increasing global id), then the ghost nodes (nodes of
!>     local elements owned by another rank, increasing global id). An owned node need not belong
!>     to a local element (Neper's partitioner can give a node to a partition none of whose
!>     elements contains it).
!>   - Degrees of freedom: local node i has the local dofs 3(i - 1) + 1..3 (x, y, z), so the
!>     owned dofs are 1..3 num_nodes_owned, followed by the ghost dofs.
!>   - After mesh_local_mod::mesh_localize (step 2.3): elt_nodes_loc holds local node numbers,
!>     and the nodes of the nsets and fasets are local node numbers (before: global ids).
!>     Periodicity relations keep global node ids (they may refer to nodes of other ranks).
!>   - Communication graph (comm_graph_t): for each neighbor rank, the ghost nodes received from
!>     it (its owned nodes) and the owned nodes sent to it (its ghosts), both in increasing global
!>     id, so that a send list of rank p for rank q and the receive list of q from p match entry
!>     by entry.
!>
!> Public entities:
!>   mesh_t                      rank-local mesh
!>   comm_graph_t                communication graph of the owned/ghost nodes
!>   node_loc_from_glob          local node number of a global node id (0 if not local)
!>   elt_loc_from_glob           local element number of a global element number (0 if not local)
!>   node_dofs                   the 3 local dofs of a local node
!>   find_sorted                 position of a value in an increasing array (0 if absent)
module mesh_mod

  use kinds_mod, only: rk, ik
  use mesh_file_mod, only: mesh_nset_t, mesh_faset_t, mesh_values_t, num_value_sections

  implicit none

  private

  public :: mesh_t, comm_graph_t
  public :: node_loc_from_glob, elt_loc_from_glob, node_dofs, find_sorted

  !> Communication graph of the owned and ghost nodes (see the module layout).
  !>
  !> Components:
  !>   num_neighbors  integer(ik)     number of neighbor ranks (ranks this rank sends to or
  !>                                  receives from)
  !>   neighbor       integer(ik)(num_neighbors)      neighbor ranks (0-based), increasing
  !>   recv_offset    integer(ik)(num_neighbors + 1)  ghosts from neighbor k: recv_nodes(
  !>                                  recv_offset(k) : recv_offset(k + 1) - 1)
  !>   recv_nodes     integer(ik)(num_nodes_ghost)    local node numbers of the ghosts, grouped by
  !>                                  neighbor, increasing global id within a group
  !>   send_offset    integer(ik)(num_neighbors + 1)  owned nodes for neighbor k: send_nodes(
  !>                                  send_offset(k) : send_offset(k + 1) - 1)
  !>   send_nodes     integer(ik)(:)  local node numbers of the owned nodes each neighbor holds as
  !>                                  ghosts, grouped by neighbor, increasing global id in a group
  type :: comm_graph_t
    integer(ik) :: num_neighbors = 0_ik
    integer(ik), allocatable :: neighbor(:)
    integer(ik), allocatable :: recv_offset(:)
    integer(ik), allocatable :: recv_nodes(:)
    integer(ik), allocatable :: send_offset(:)
    integer(ik), allocatable :: send_nodes(:)
  end type comm_graph_t

  !> Rank-local mesh.
  !>
  !> Components (sizes; units):
  !>   num_parts_file   integer(ik)    partitions of the mesh file (0: not partitioned; the
  !>                                   setup summary, step 9.7)
  !>   num_nodes_glob   integer(ik)    number of nodes of the whole mesh
  !>   num_elts_glob    integer(ik)    number of 3D elements of the whole mesh
  !>   elt_type         integer(ik)    Gmsh type of the 3D elements (11 or 17)
  !>   num_elt_nodes    integer(ik)    nodes per element (10 or 20)
  !>   domain           character(:)   $Domain of the mesh file (e.g. "cube"), or empty
  !>   num_elts         integer(ik)    elements of this rank
  !>   elt_id_glob      integer(ik)(num_elts)  global element number, increasing
  !>   elt_nodes_glob   integer(ik)(num_elt_nodes, num_elts)  global node ids, Gmsh order
  !>                                   (the internal element node order, ROADMAP D35)
  !>   elt_nodes_loc    integer(ik)(num_elt_nodes, num_elts)  local node numbers, same order
  !>   elt_dofs_loc     integer(ik)(3 num_elt_nodes, num_elts)  local dofs 3 (a_loc - 1) + k of
  !>                                   node a (local number a_loc) and component k, element node
  !>                                   order (step 4.4: built once for the element-by-element
  !>                                   operations)
  !>                                   (step 2.3)
  !>   elt_elset        integer(ik)(num_elts)  elset index (1..num_elsets)
  !>   elt_phase        integer(ik)(num_elts)  phase
  !>   rot_c2s          real(rk)(3, 3, num_elts)  crystal-to-sample rotation (zero if has_ori
  !>                                   is false)
  !>   num_nodes_owned  integer(ik)    nodes owned by this rank
  !>   num_nodes_ghost  integer(ik)    ghost nodes
  !>   num_nodes_loc    integer(ik)    local nodes, num_nodes_owned + num_nodes_ghost
  !>   node_id_glob     integer(ik)(owned + ghost)  global node ids (owned, then ghosts; each
  !>                                   group increasing)
  !>   ghost_owner      integer(ik)(num_nodes_ghost)  rank (0-based) owning each ghost node
  !>   coo_s            real(rk)(3, owned + ghost)  node coordinates, sample frame [mm]
  !>   num_elsets       integer(ik)    number of elsets (whole mesh)
  !>   elset_id         integer(ik)(num_elsets)  elset ids of the mesh file
  !>   num_phases       integer(ik)    number of phases
  !>   elset_phase      integer(ik)(num_elsets)  phase of each elset
  !>   has_ori          logical        whether orientations were given
  !>   ori_from_elts    logical        element (.true.) or elset orientations
  !>   ori_descriptor, ori_convention  integer  descriptor and convention codes (for output)
  !>   nsets            mesh_nset_t(:)   node sets, restricted to the local (owned and ghost)
  !>                                   nodes, file order; local node numbers (global ids before
  !>                                   mesh_localize)
  !>   nsets_glob       mesh_nset_t(:)   node sets of the whole mesh, replicated: global ids, file
  !>                                   order (never localized; for the replicated constraint table
  !>                                   of multi-point constraints, step 8.3a, ROADMAP D74)
  !>   fasets           mesh_faset_t(:)  face sets, restricted to the faces of local elements
  !>                                   (elts: local element number); local node numbers (global
  !>                                   ids before mesh_localize)
  !>   num_periodicity  integer(ik)    periodicity relations (whole mesh, replicated)
  !>   per_secondary, per_primary  integer(ik)(num_periodicity)  global node ids
  !>   per_vect         integer(ik)(3, num_periodicity)  scaled vector primary -> secondary
  !>   values           mesh_values_t(num_value_sections)  value sections: per elset
  !>                                   (replicated, by elset index) or per element (local
  !>                                   elements, by local element index)
  !>   comm             comm_graph_t     communication graph (step 2.3)
  type :: mesh_t
    integer(ik) :: num_parts_file = 0_ik
    integer(ik) :: num_nodes_glob = 0_ik
    integer(ik) :: num_elts_glob = 0_ik
    integer(ik) :: elt_type = 0_ik
    integer(ik) :: num_elt_nodes = 0_ik
    character(len=:), allocatable :: domain
    integer(ik) :: num_elts = 0_ik
    integer(ik), allocatable :: elt_id_glob(:)
    integer(ik), allocatable :: elt_nodes_glob(:, :)
    integer(ik), allocatable :: elt_nodes_loc(:, :)
    integer(ik), allocatable :: elt_dofs_loc(:, :)
    integer(ik), allocatable :: elt_elset(:)
    integer(ik), allocatable :: elt_phase(:)
    real(rk), allocatable :: rot_c2s(:, :, :)
    integer(ik) :: num_nodes_owned = 0_ik
    integer(ik) :: num_nodes_ghost = 0_ik
    integer(ik) :: num_nodes_loc = 0_ik
    integer(ik), allocatable :: node_id_glob(:)
    integer(ik), allocatable :: ghost_owner(:)
    real(rk), allocatable :: coo_s(:, :)
    integer(ik) :: num_elsets = 0_ik
    integer(ik), allocatable :: elset_id(:)
    integer(ik) :: num_phases = 0_ik
    integer(ik), allocatable :: elset_phase(:)
    logical :: has_ori = .false.
    logical :: ori_from_elts = .false.
    integer :: ori_descriptor = 0
    integer :: ori_convention = 0
    type(mesh_nset_t), allocatable :: nsets(:)
    type(mesh_nset_t), allocatable :: nsets_glob(:)
    type(mesh_faset_t), allocatable :: fasets(:)
    integer(ik) :: num_periodicity = 0_ik
    integer(ik), allocatable :: per_secondary(:)
    integer(ik), allocatable :: per_primary(:)
    integer(ik), allocatable :: per_vect(:, :)
    type(mesh_values_t) :: values(num_value_sections)
    type(comm_graph_t) :: comm
  end type mesh_t

contains

  !> Local node number of a global node id: binary search in the owned and in the ghost part of
  !> node_id_glob (both increasing); 0 if the node is not local.
  !>
  !> Arguments:
  !>   mesh     (in)  type(mesh_t)  rank-local mesh
  !>   id_glob  (in)  integer(ik)   global node id
  !>   i_loc    (out, result) integer(ik)  local node number, or 0
  pure function node_loc_from_glob(mesh, id_glob) result(i_loc)

    type(mesh_t), intent(in) :: mesh
    integer(ik), intent(in) :: id_glob
    integer(ik) :: i_loc

    ! Owned nodes first, then ghosts (offset by the owned count).
    i_loc = find_sorted(mesh%node_id_glob(1:mesh%num_nodes_owned), id_glob)
    if (i_loc > 0_ik) return
    i_loc = find_sorted(mesh%node_id_glob(mesh%num_nodes_owned + 1:mesh%num_nodes_loc), id_glob)
    if (i_loc > 0_ik) i_loc = i_loc + mesh%num_nodes_owned

  end function node_loc_from_glob

  !> Local element number of a global element number (0 if not local).
  !>
  !> Arguments:
  !>   mesh     (in)  type(mesh_t)  rank-local mesh
  !>   id_glob  (in)  integer(ik)   global element number
  !>   i_loc    (out, result) integer(ik)  local element number, or 0
  pure function elt_loc_from_glob(mesh, id_glob) result(i_loc)

    type(mesh_t), intent(in) :: mesh
    integer(ik), intent(in) :: id_glob
    integer(ik) :: i_loc

    ! The local elements are in increasing global order.
    i_loc = find_sorted(mesh%elt_id_glob, id_glob)

  end function elt_loc_from_glob

  !> The three local degrees of freedom (x, y, z) of a local node: 3(i - 1) + 1..3.
  !>
  !> Arguments:
  !>   i_node  (in)  integer(ik)  local node number
  !>   dofs    (out, result) integer(ik)(3)  its local dofs
  pure function node_dofs(i_node) result(dofs)

    integer(ik), intent(in) :: i_node
    integer(ik) :: dofs(3)

    dofs = 3_ik*(i_node - 1_ik) + [1_ik, 2_ik, 3_ik]

  end function node_dofs

  !> Position of a value in an increasing array (0 if absent), by bisection.
  !>
  !> Arguments:
  !>   sorted  (in)  integer(ik)(:)  increasing values
  !>   value   (in)  integer(ik)     value searched
  !>   pos     (out, result) integer(ik)
  !> Locals:
  !>   low, high, mid  integer(ik)  bisection bounds
  pure function find_sorted(sorted, value) result(pos)

    integer(ik), intent(in) :: sorted(:)
    integer(ik), intent(in) :: value
    integer(ik) :: pos

    integer(ik) :: low
    integer(ik) :: high
    integer(ik) :: mid

    pos = 0_ik
    low = 1_ik
    high = int(size(sorted), ik)
    do while (low <= high)
      mid = (low + high)/2_ik
      if (sorted(mid) == value) then
        pos = mid
        return
      else if (sorted(mid) < value) then
        low = mid + 1_ik
      else
        high = mid - 1_ik
      end if
    end do

  end function find_sorted

end module mesh_mod
