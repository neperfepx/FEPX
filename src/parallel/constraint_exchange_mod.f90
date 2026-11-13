! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Transfer of the values of independent dofs between ranks for the constraint relations
!> v = T v_hat + g (multi-point constraints, periodicity; ROADMAP step 8.3a, D74;
!> dev/theory/addendum_periodic_mpc.md §3.1).
!>
!> Source : FEPX (D74, user decision Q-P14: the least intrusive way for FEPX's parallel
!>          implementation, without changing its existing parallel routines: halo_mod,
!>          parallel_mod and the communication graph of mesh_local_mod are used as they are and
!>          not modified). The message pattern follows halo_mod (non-blocking point-to-point
!>          messages, MPI_Irecv / MPI_Isend / MPI_Waitall, one exchange at a time); the set-up
!>          follows mesh_local_mod::build_comm_graph (the requests sent to their owners with
!>          par_alltoall / par_alltoallv). Contains no numerical method. v2 applies T and T^T
!>          through its element gather/scatter with a division by the dof multiplicity
!>          (conjugate_gradient_mod2.f90::prim_to_all / all_to_prim; R-P4 §3.2.5); this module
!>          replaces that transfer for FEPX's owned + ghost layout (D5).
!>
!> Policy:
!>   - A rank's constrained dofs (owned or ghost) depend on independent dofs (columns of T) that
!>     may be owned by other ranks, and that need not be local at all. Each such dof that the
!>     rank does not own is a "slot": the rank holds one value per slot in a buffer, numbered
!>     1..num_slots, grouped by owner rank (increasing) and in increasing global dof within a
!>     group.
!>   - Set-up (collective, once): the owner of every candidate column is found from a replicated
!>     list of all columns (every rank marks those it owns; a sum over the ranks); each rank sends
!>     each owner the global dofs of its slots (par_alltoall, par_alltoallv); the owner converts
!>     them to its local (owned) dofs, its send lists. Neighbors: ranks this rank sends to or
!>     receives from.
!>   - constraint_forward: each owner sends the values of its send lists; each rank stores them
!>     in its slot buffer.
!>   - constraint_reverse_add: each rank sends its slot buffer to the owners; each owner adds the
!>     received values to its owned entries, neighbor by neighbor in increasing rank order and in
!>     list order (fixed order: the result is reproducible for a given partitioning).
!>   - One message tag (constraint_tag, distinct from halo_mod's). Neighbors with an empty list
!>     in the exchange's direction get no message.
!>
!> Public entities:
!>   constraint_exchange_t      the slots and lists of one rank (components at the type)
!>   constraint_exchange_setup  build the slots and lists (collective)
!>   constraint_forward         owners' values -> slot buffers (collective over the neighbors)
!>   constraint_reverse_add     slot buffers added to the owners' values
module constraint_exchange_mod

  use mpi_f08, only: MPI_Irecv, MPI_Isend, MPI_Waitall, MPI_Request, MPI_COMM_WORLD, &
    & MPI_STATUSES_IGNORE
  use kinds_mod, only: rk, ik
  use error_mod, only: internal_error
  use parallel_mod, only: num_ranks, my_rank, par_sum, par_alltoall, par_alltoallv, &
    & par_check_mpi, mpi_type_rk
  use mesh_mod, only: mesh_t, node_loc_from_glob

  implicit none

  private

  public :: constraint_exchange_t, constraint_exchange_setup, constraint_forward, &
    & constraint_reverse_add

  !> Message tag of the constraint exchanges (halo_mod uses 2401).
  integer, parameter :: constraint_tag = 2402

  !> The slots and lists of one rank.
  !>
  !> Components:
  !>   num_slots      integer(ik)     independent dofs this rank needs and does not own
  !>   slot_dof_glob  integer(ik)(num_slots)  global dof 3 (node id - 1) + k of each slot,
  !>                                  grouped by owner (increasing rank), increasing in a group
  !>   num_neighbors  integer(ik)     ranks this rank sends to or receives from
  !>   neighbor       integer(ik)(num_neighbors)      neighbor ranks (0-based), increasing
  !>   recv_offset    integer(ik)(num_neighbors + 1)  slots served by neighbor k: recv_offset(k)
  !>                                  .. recv_offset(k + 1) - 1
  !>   send_offset    integer(ik)(num_neighbors + 1)  entries of send_dofs for neighbor k
  !>   send_dofs      integer(ik)(:)  local owned dofs whose values neighbor k needs (its slots,
  !>                                  in its slot order)
  type :: constraint_exchange_t
    integer(ik) :: num_slots = 0_ik
    integer(ik), allocatable :: slot_dof_glob(:)
    integer(ik) :: num_neighbors = 0_ik
    integer(ik), allocatable :: neighbor(:)
    integer(ik), allocatable :: recv_offset(:)
    integer(ik), allocatable :: send_offset(:)
    integer(ik), allocatable :: send_dofs(:)
  end type constraint_exchange_t

contains

  !> Build the slots and lists (module policy). Collective: every rank must call it, with the
  !> same all_cols.
  !>
  !> Arguments:
  !>   cx        (out)  type(constraint_exchange_t)  this rank's slots and lists
  !>   mesh      (in)   type(mesh_t)    this rank's mesh (owned nodes)
  !>   all_cols  (in)   integer(ik)(:)  every column (independent global dof) of the replicated
  !>                                    constraint table, increasing, without repetition (the
  !>                                    same on every rank)
  !>   needed    (in)   logical(:)      (size(all_cols)) columns this rank needs (its local rows
  !>                                    refer to them); those it does not own become slots
  !> Locals:
  !>   owner_mark   integer(ik)(:)  my_rank + 1 at the columns this rank owns, 0 elsewhere
  !>   owner_p1     integer(ik)(:)  owner rank + 1 of each column (sum over the ranks)
  !>   req_counts   integer(ik)(num_ranks)  slots of this rank owned by each rank
  !>   need_counts  integer(ik)(num_ranks)  dofs of this rank requested by each rank
  !>   need_ids     integer(ik)(:)  global dofs requested from this rank, grouped by requester
  !>   i_col        integer(ik)     column index
  !>   i_rank       integer(ik)     rank (0-based)
  !>   k            integer(ik)     neighbor index
  !>   i            integer(ik)     entry index
  !>   i_loc        integer(ik)     local node of a requested dof
  !>   n_slot       integer(ik)     slots placed so far
  !>   send_pos     integer(ik)     position in need_ids of the current requester's group
  subroutine constraint_exchange_setup(cx, mesh, all_cols, needed)

    type(constraint_exchange_t), intent(out) :: cx
    type(mesh_t), intent(in) :: mesh
    integer(ik), intent(in) :: all_cols(:)
    logical, intent(in) :: needed(:)

    integer(ik), allocatable :: owner_mark(:)
    integer(ik), allocatable :: owner_p1(:)
    integer(ik) :: req_counts(num_ranks)
    integer(ik) :: need_counts(num_ranks)
    integer(ik), allocatable :: need_ids(:)
    integer(ik) :: i_col
    integer(ik) :: i_rank
    integer(ik) :: k
    integer(ik) :: i
    integer(ik) :: i_loc
    integer(ik) :: n_slot
    integer(ik) :: send_pos

    if (size(needed) /= size(all_cols)) call internal_error('constraint_exchange_setup', &
      & 'needed and all_cols differ in size')

    ! Owner of every column: each rank marks the columns of its owned nodes; the sum over the
    ! ranks is owner + 1 (every node has exactly one owner).
    allocate (owner_mark(size(all_cols)), owner_p1(size(all_cols)))
    owner_mark = 0_ik
    do i_col = 1_ik, int(size(all_cols), ik)
      i_loc = node_loc_from_glob(mesh, (all_cols(i_col) - 1_ik)/3_ik + 1_ik)
      if (i_loc >= 1_ik .and. i_loc <= mesh%num_nodes_owned) owner_mark(i_col) = my_rank + 1_ik
    end do
    call par_sum(owner_mark, owner_p1)
    if (any(owner_p1 < 1_ik .or. owner_p1 > num_ranks)) call internal_error( &
      & 'constraint_exchange_setup', 'a constraint column has no single owner')

    ! Slots: the needed columns owned elsewhere, grouped by owner (all_cols is increasing, so
    ! each group is increasing too).
    req_counts = 0_ik
    do i_col = 1_ik, int(size(all_cols), ik)
      if (.not. needed(i_col) .or. owner_p1(i_col) - 1_ik == my_rank) cycle
      req_counts(owner_p1(i_col)) = req_counts(owner_p1(i_col)) + 1_ik
    end do
    cx%num_slots = sum(req_counts)
    allocate (cx%slot_dof_glob(cx%num_slots))
    n_slot = 0_ik
    do i_rank = 0_ik, num_ranks - 1_ik
      if (i_rank == my_rank) cycle
      do i_col = 1_ik, int(size(all_cols), ik)
        if (.not. needed(i_col) .or. owner_p1(i_col) - 1_ik /= i_rank) cycle
        n_slot = n_slot + 1_ik
        cx%slot_dof_glob(n_slot) = all_cols(i_col)
      end do
    end do

    ! Tell each owner which of its dofs this rank needs.
    call par_alltoall(req_counts, need_counts)
    call par_alltoallv(cx%slot_dof_glob, req_counts, need_ids, need_counts)

    ! Neighbors, receive ranges (slots) and send lists (local owned dofs), in rank order.
    cx%num_neighbors = int(count(req_counts > 0_ik .or. need_counts > 0_ik), ik)
    allocate (cx%neighbor(cx%num_neighbors), cx%recv_offset(cx%num_neighbors + 1_ik))
    allocate (cx%send_offset(cx%num_neighbors + 1_ik), cx%send_dofs(size(need_ids)))
    cx%recv_offset(1) = 1_ik
    cx%send_offset(1) = 1_ik
    k = 0_ik
    send_pos = 1_ik
    do i_rank = 0_ik, num_ranks - 1_ik
      if (req_counts(i_rank + 1_ik) == 0_ik .and. need_counts(i_rank + 1_ik) == 0_ik) cycle
      k = k + 1_ik
      cx%neighbor(k) = i_rank
      cx%recv_offset(k + 1_ik) = cx%recv_offset(k) + req_counts(i_rank + 1_ik)
      do i = send_pos, send_pos + need_counts(i_rank + 1_ik) - 1_ik
        ! The requested global dof as a local owned dof of this rank.
        i_loc = node_loc_from_glob(mesh, (need_ids(i) - 1_ik)/3_ik + 1_ik)
        if (i_loc < 1_ik .or. i_loc > mesh%num_nodes_owned) call internal_error( &
          & 'constraint_exchange_setup', 'a requested constraint dof is not owned here')
        cx%send_dofs(i) = 3_ik*(i_loc - 1_ik) + mod(need_ids(i) - 1_ik, 3_ik) + 1_ik
      end do
      send_pos = send_pos + need_counts(i_rank + 1_ik)
      cx%send_offset(k + 1_ik) = send_pos
    end do

  end subroutine constraint_exchange_setup

  !> Owners' values -> slot buffers (module policy). Collective over the neighbors.
  !>
  !> Arguments:
  !>   cx   (in)   type(constraint_exchange_t)  slots and lists
  !>   x    (in)   real(rk)(:)  dof vector (3 num_nodes_loc); its owned entries are sent
  !>   buf  (out)  real(rk)(:)  (num_slots) the values of the slots
  !> Locals:
  !>   send_buf  real(rk)(:)  packed outgoing values (ASYNCHRONOUS: buffer of non-blocking calls)
  !>   recv_buf  real(rk)(:)  incoming values (ASYNCHRONOUS)
  !>   requests  type(MPI_Request)(:)  receive requests, then send requests
  !>   num_req   integer      requests posted
  !>   k         integer(ik)  neighbor index
  !>   first, last  integer(ik)  list range of a neighbor
  !>   ierror    integer      MPI error code
  subroutine constraint_forward(cx, x, buf)

    type(constraint_exchange_t), intent(in) :: cx
    real(rk), intent(in) :: x(:)
    real(rk), intent(out) :: buf(:)

    real(rk), allocatable, asynchronous :: send_buf(:)
    real(rk), allocatable, asynchronous :: recv_buf(:)
    type(MPI_Request), allocatable :: requests(:)
    integer :: num_req
    integer(ik) :: k
    integer(ik) :: first
    integer(ik) :: last
    integer :: ierror

    if (size(buf) /= cx%num_slots) call internal_error('constraint_forward', &
      & 'buffer size differs from the number of slots')
    if (cx%num_neighbors == 0_ik) return
    allocate (recv_buf(cx%num_slots), send_buf(size(cx%send_dofs)))
    allocate (requests(2*cx%num_neighbors))
    num_req = 0

    ! Post the receives (the slots served by each neighbor).
    do k = 1_ik, cx%num_neighbors
      first = cx%recv_offset(k)
      last = cx%recv_offset(k + 1_ik) - 1_ik
      if (last < first) cycle
      num_req = num_req + 1
      call MPI_Irecv(recv_buf(first:last), int(last - first + 1_ik), mpi_type_rk, &
        & int(cx%neighbor(k)), constraint_tag, MPI_COMM_WORLD, requests(num_req), ierror)
      call par_check_mpi(ierror, 'constraint exchange', 'MPI_Irecv')
    end do

    ! Pack and send the requested owned values.
    send_buf = x(cx%send_dofs)
    do k = 1_ik, cx%num_neighbors
      first = cx%send_offset(k)
      last = cx%send_offset(k + 1_ik) - 1_ik
      if (last < first) cycle
      num_req = num_req + 1
      call MPI_Isend(send_buf(first:last), int(last - first + 1_ik), mpi_type_rk, &
        & int(cx%neighbor(k)), constraint_tag, MPI_COMM_WORLD, requests(num_req), ierror)
      call par_check_mpi(ierror, 'constraint exchange', 'MPI_Isend')
    end do

    ! Wait for all messages; the slot values.
    call MPI_Waitall(num_req, requests(1:num_req), MPI_STATUSES_IGNORE, ierror)
    call par_check_mpi(ierror, 'constraint exchange', 'MPI_Waitall')
    buf = recv_buf

  end subroutine constraint_forward

  !> Slot buffers added to the owners' values (module policy). Collective over the neighbors.
  !>
  !> Arguments:
  !>   cx   (in)     type(constraint_exchange_t)  slots and lists
  !>   buf  (in)     real(rk)(:)  (num_slots) contributions of this rank to its slots
  !>   x    (inout)  real(rk)(:)  dof vector (3 num_nodes_loc); the contributions of the other
  !>                              ranks are added to its owned entries
  !> Locals:
  !>   send_buf  real(rk)(:)  outgoing contributions (ASYNCHRONOUS)
  !>   recv_buf  real(rk)(:)  incoming contributions, in send-list order (ASYNCHRONOUS)
  !>   requests  type(MPI_Request)(:)  receive requests, then send requests
  !>   num_req   integer      requests posted
  !>   k         integer(ik)  neighbor index
  !>   i         integer(ik)  entry index
  !>   first, last  integer(ik)  list range of a neighbor
  !>   ierror    integer      MPI error code
  subroutine constraint_reverse_add(cx, buf, x)

    type(constraint_exchange_t), intent(in) :: cx
    real(rk), intent(in) :: buf(:)
    real(rk), intent(inout) :: x(:)

    real(rk), allocatable, asynchronous :: send_buf(:)
    real(rk), allocatable, asynchronous :: recv_buf(:)
    type(MPI_Request), allocatable :: requests(:)
    integer :: num_req
    integer(ik) :: k
    integer(ik) :: i
    integer(ik) :: first
    integer(ik) :: last
    integer :: ierror

    if (size(buf) /= cx%num_slots) call internal_error('constraint_reverse_add', &
      & 'buffer size differs from the number of slots')
    if (cx%num_neighbors == 0_ik) return
    allocate (recv_buf(size(cx%send_dofs)), send_buf(cx%num_slots))
    allocate (requests(2*cx%num_neighbors))
    num_req = 0

    ! Post the receives (the contributions to this rank's send lists).
    do k = 1_ik, cx%num_neighbors
      first = cx%send_offset(k)
      last = cx%send_offset(k + 1_ik) - 1_ik
      if (last < first) cycle
      num_req = num_req + 1
      call MPI_Irecv(recv_buf(first:last), int(last - first + 1_ik), mpi_type_rk, &
        & int(cx%neighbor(k)), constraint_tag, MPI_COMM_WORLD, requests(num_req), ierror)
      call par_check_mpi(ierror, 'constraint exchange', 'MPI_Irecv')
    end do

    ! Send the contributions of each owner's slots.
    send_buf = buf
    do k = 1_ik, cx%num_neighbors
      first = cx%recv_offset(k)
      last = cx%recv_offset(k + 1_ik) - 1_ik
      if (last < first) cycle
      num_req = num_req + 1
      call MPI_Isend(send_buf(first:last), int(last - first + 1_ik), mpi_type_rk, &
        & int(cx%neighbor(k)), constraint_tag, MPI_COMM_WORLD, requests(num_req), ierror)
      call par_check_mpi(ierror, 'constraint exchange', 'MPI_Isend')
    end do

    ! Wait for all messages, then add in neighbor order and list order.
    call MPI_Waitall(num_req, requests(1:num_req), MPI_STATUSES_IGNORE, ierror)
    call par_check_mpi(ierror, 'constraint exchange', 'MPI_Waitall')
    do i = 1_ik, int(size(cx%send_dofs), ik)
      x(cx%send_dofs(i)) = x(cx%send_dofs(i)) + recv_buf(i)
    end do

  end subroutine constraint_reverse_add

end module constraint_exchange_mod
