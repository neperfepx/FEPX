! This file is part of the FEPX software package.
! Copyright (C) 1996-2026, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

!> Halo exchange of nodal vectors (owned <-> ghost nodes) and distributed vector operations
!> (ROADMAP step 2.4).
!>
!> Source : FEPX (D5: owned + ghost nodes, halo exchange). Replaces v2
!>          gather_scatter_mod.f90::part_gather / part_scatter (element gather/scatter with
!>          non-blocking point-to-point messages) for the owned + ghost layout; the dot product
!>          and norm follow v2 parallel_mod.f90::par_sum / par_sqrtsum (local sum, then an
!>          allreduce). Contains no numerical method.
!>
!> Policy:
!>   - A nodal vector holds num_comp values per local node, node after node:
!>     x((i - 1) num_comp + c) is component c of local node i (mesh_mod layout; num_comp = 3 for
!>     a dof vector, whose owned entries are 1..3 num_nodes_owned). Its size must be
!>     num_comp num_nodes_loc.
!>   - halo_forward: each owner sends the values of the nodes in its send lists; each rank
!>     overwrites its ghost values with them. Owned values are unchanged.
!>   - halo_reverse_add: each rank sends the values of its ghost nodes; each owner adds them to
!>     its owned values, neighbor by neighbor in increasing rank order (fixed order: the result
!>     is reproducible from run to run for a given partitioning). Ghost values are unchanged
!>     (a following halo_forward makes them equal to the owners' sums).
!>   - Messages: MPI_Irecv for every neighbor, then MPI_Isend, then MPI_Waitall (as v2's
!>     part_gather_start / _stop); one tag; one exchange at a time. Neighbors with an empty
!>     list in the exchange's direction get no message.
!>   - dist_dot / dist_norm: sum over the owned entries only (each node counted once), then
!>     par_sum; the same value on every rank.
!>
!> Public entities:
!>   halo_forward      owners' values -> ghosts (collective over the neighbors)
!>   halo_reverse_add  ghosts' values added to the owners' values
!>   dist_dot          dot product of two distributed vectors (owned entries)
!>   dist_norm         Euclidean norm of a distributed vector (owned entries)
module halo_mod

  use mpi_f08, only: MPI_Irecv, MPI_Isend, MPI_Waitall, MPI_Request, MPI_COMM_WORLD, &
    & MPI_STATUSES_IGNORE
  use kinds_mod, only: rk, ik
  use error_mod, only: internal_error
  use parallel_mod, only: par_sum, par_check_mpi, mpi_type_rk
  use mesh_mod, only: mesh_t

  implicit none

  private

  public :: halo_forward, halo_reverse_add, dist_dot, dist_norm

  !> Message tag of the halo exchanges.
  integer, parameter :: halo_tag = 2401

contains

  !> Overwrite the ghost values of a nodal vector with the owners' values (module policy).
  !>
  !> Arguments:
  !>   mesh      (in)     type(mesh_t)  this rank's mesh (communication graph)
  !>   x         (inout)  real(rk)(:)   nodal vector, num_comp num_nodes_loc entries
  !>   num_comp  (in)     integer(ik)   values per node
  subroutine halo_forward(mesh, x, num_comp)

    type(mesh_t), intent(in) :: mesh
    real(rk), intent(inout) :: x(:)
    integer(ik), intent(in) :: num_comp

    ! Owners send their send-list values; ghosts receive into their receive-list positions.
    call check_size(mesh, x, num_comp, 'halo_forward')
    call exchange(mesh, x, num_comp, .false.)

  end subroutine halo_forward

  !> Add the ghost values of a nodal vector to the owners' values (module policy).
  !>
  !> Arguments:
  !>   mesh      (in)     type(mesh_t)  this rank's mesh (communication graph)
  !>   x         (inout)  real(rk)(:)   nodal vector, num_comp num_nodes_loc entries
  !>   num_comp  (in)     integer(ik)   values per node
  subroutine halo_reverse_add(mesh, x, num_comp)

    type(mesh_t), intent(in) :: mesh
    real(rk), intent(inout) :: x(:)
    integer(ik), intent(in) :: num_comp

    ! Ghosts send their values; owners add them, neighbor by neighbor in rank order.
    call check_size(mesh, x, num_comp, 'halo_reverse_add')
    call exchange(mesh, x, num_comp, .true.)

  end subroutine halo_reverse_add

  !> Dot product of two distributed vectors over their owned entries (v2 par_sum pattern).
  !>
  !> Arguments:
  !>   x          (in)  real(rk)(:)  first vector
  !>   y          (in)  real(rk)(:)  second vector
  !>   num_owned  (in)  integer(ik)  number of owned entries (the first num_owned of x and y;
  !>                                 3 num_nodes_owned for a dof vector)
  !>   dot        (out, result) real(rk)  sum over all ranks of x . y on the owned entries
  !> Locals:
  !>   local_dot  real(rk)  this rank's part
  function dist_dot(x, y, num_owned) result(dot)

    real(rk), intent(in) :: x(:)
    real(rk), intent(in) :: y(:)
    integer(ik), intent(in) :: num_owned
    real(rk) :: dot

    real(rk) :: local_dot

    ! Owned entries only (each node once), then the sum over the ranks.
    if (num_owned < 0_ik .or. num_owned > size(x) .or. num_owned > size(y)) then
      call internal_error('dist_dot', 'num_owned exceeds the vector size')
    end if
    local_dot = dot_product(x(1:num_owned), y(1:num_owned))
    call par_sum(local_dot, dot)

  end function dist_dot

  !> Euclidean norm of a distributed vector over its owned entries (v2 par_sqrtsum pattern).
  !>
  !> Arguments:
  !>   x          (in)  real(rk)(:)  the vector
  !>   num_owned  (in)  integer(ik)  number of owned entries
  !>   norm       (out, result) real(rk)  sqrt of the sum over all ranks of x . x (owned)
  function dist_norm(x, num_owned) result(norm)

    real(rk), intent(in) :: x(:)
    integer(ik), intent(in) :: num_owned
    real(rk) :: norm

    norm = sqrt(dist_dot(x, x, num_owned))

  end function dist_norm

  !> One halo exchange in either direction: forward sends the send lists and receives into the
  !> receive lists; reverse the opposite (the lists are passed, not copied).
  !>
  !> Arguments:
  !>   mesh      (in)     type(mesh_t)  this rank's mesh
  !>   x         (inout)  real(rk)(:)   nodal vector
  !>   num_comp  (in)     integer(ik)   values per node
  !>   reverse   (in)     logical       .false.: owners -> ghosts (overwrite);
  !>                                    .true.: ghosts -> owners (add)
  subroutine exchange(mesh, x, num_comp, reverse)

    type(mesh_t), intent(in) :: mesh
    real(rk), intent(inout) :: x(:)
    integer(ik), intent(in) :: num_comp
    logical, intent(in) :: reverse

    associate (comm => mesh%comm)
      if (comm%num_neighbors == 0_ik) return
      if (reverse) then
        call exchange_lists(x, num_comp, comm%neighbor, comm%recv_offset, comm%recv_nodes, &
          & comm%send_offset, comm%send_nodes, .true.)
      else
        call exchange_lists(x, num_comp, comm%neighbor, comm%send_offset, comm%send_nodes, &
          & comm%recv_offset, comm%recv_nodes, .false.)
      end if
    end associate

  end subroutine exchange

  !> Send the values of the outgoing node lists and receive into the incoming node lists, one
  !> list pair per neighbor; incoming values overwrite or are added (module policy).
  !>
  !> Arguments:
  !>   x           (inout)  real(rk)(:)     nodal vector
  !>   num_comp    (in)     integer(ik)     values per node
  !>   neighbor    (in)     integer(ik)(:)  neighbor ranks
  !>   out_offset  (in)     integer(ik)(:)  offsets of the outgoing lists (num_neighbors + 1)
  !>   out_nodes   (in)     integer(ik)(:)  outgoing node lists (local numbers)
  !>   in_offset   (in)     integer(ik)(:)  offsets of the incoming lists
  !>   in_nodes    (in)     integer(ik)(:)  incoming node lists (local numbers)
  !>   add         (in)     logical         add (.true.) or overwrite (.false.) incoming values
  !> Locals:
  !>   send_buf, recv_buf  real(rk)(:)     packed values (ASYNCHRONOUS: buffers of
  !>                                       non-blocking calls; contiguous sections are passed)
  !>   requests            type(MPI_Request)(:)  receive requests, then send requests
  !>   num_req             integer         requests posted
  !>   k                   integer(ik)     neighbor index
  !>   first, last         integer(ik)     list range of a neighbor
  !>   ierror              integer         MPI error code
  subroutine exchange_lists(x, num_comp, neighbor, out_offset, out_nodes, in_offset, in_nodes, &
    & add)

    real(rk), intent(inout) :: x(:)
    integer(ik), intent(in) :: num_comp
    integer(ik), intent(in) :: neighbor(:)
    integer(ik), intent(in) :: out_offset(:)
    integer(ik), intent(in) :: out_nodes(:)
    integer(ik), intent(in) :: in_offset(:)
    integer(ik), intent(in) :: in_nodes(:)
    logical, intent(in) :: add

    real(rk), allocatable, asynchronous :: send_buf(:)
    real(rk), allocatable, asynchronous :: recv_buf(:)
    type(MPI_Request), allocatable :: requests(:)
    integer :: num_req
    integer(ik) :: k
    integer(ik) :: first
    integer(ik) :: last
    integer :: ierror

    allocate (recv_buf(num_comp*size(in_nodes)), send_buf(num_comp*size(out_nodes)))
    allocate (requests(2*size(neighbor)))
    num_req = 0

    ! Post the receives (one per neighbor with a non-empty incoming list).
    do k = 1_ik, int(size(neighbor), ik)
      first = in_offset(k)
      last = in_offset(k + 1_ik) - 1_ik
      if (last < first) cycle
      num_req = num_req + 1
      call MPI_Irecv(recv_buf(num_comp*(first - 1_ik) + 1_ik:num_comp*last), &
        & int(num_comp*(last - first + 1_ik)), mpi_type_rk, int(neighbor(k)), halo_tag, &
        & MPI_COMM_WORLD, requests(num_req), ierror)
      call par_check_mpi(ierror, 'halo exchange', 'MPI_Irecv')
    end do

    ! Pack and send the outgoing values.
    call pack_values(x, num_comp, out_nodes, send_buf)
    do k = 1_ik, int(size(neighbor), ik)
      first = out_offset(k)
      last = out_offset(k + 1_ik) - 1_ik
      if (last < first) cycle
      num_req = num_req + 1
      call MPI_Isend(send_buf(num_comp*(first - 1_ik) + 1_ik:num_comp*last), &
        & int(num_comp*(last - first + 1_ik)), mpi_type_rk, int(neighbor(k)), halo_tag, &
        & MPI_COMM_WORLD, requests(num_req), ierror)
      call par_check_mpi(ierror, 'halo exchange', 'MPI_Isend')
    end do

    ! Wait for all messages, then unpack (overwrite or add, in neighbor order).
    call MPI_Waitall(num_req, requests(1:num_req), MPI_STATUSES_IGNORE, ierror)
    call par_check_mpi(ierror, 'halo exchange', 'MPI_Waitall')
    call unpack_values(recv_buf, num_comp, in_nodes, add, x)

  end subroutine exchange_lists

  !> Pack the values of a node list: buf((j - 1) num_comp + c) = component c of node nodes(j).
  !>
  !> Arguments:
  !>   x         (in)   real(rk)(:)     nodal vector
  !>   num_comp  (in)   integer(ik)     values per node
  !>   nodes     (in)   integer(ik)(:)  local node numbers
  !>   buf       (out)  real(rk)(:)     packed values, num_comp size(nodes)
  !> Locals:
  !>   j  integer(ik)  list position
  !>   c  integer(ik)  component
  pure subroutine pack_values(x, num_comp, nodes, buf)

    real(rk), intent(in) :: x(:)
    integer(ik), intent(in) :: num_comp
    integer(ik), intent(in) :: nodes(:)
    real(rk), intent(out) :: buf(:)

    integer(ik) :: j
    integer(ik) :: c

    do j = 1_ik, int(size(nodes), ik)
      do c = 1_ik, num_comp
        buf(num_comp*(j - 1_ik) + c) = x(num_comp*(nodes(j) - 1_ik) + c)
      end do
    end do

  end subroutine pack_values

  !> Unpack received values into a node list, overwriting or adding (in list order, i.e.
  !> neighbor by neighbor in increasing rank order).
  !>
  !> Arguments:
  !>   buf       (in)     real(rk)(:)     packed values
  !>   num_comp  (in)     integer(ik)     values per node
  !>   nodes     (in)     integer(ik)(:)  local node numbers
  !>   add       (in)     logical         add (.true.) or overwrite (.false.)
  !>   x         (inout)  real(rk)(:)     nodal vector
  !> Locals:
  !>   j  integer(ik)  list position
  !>   c  integer(ik)  component
  !>   p  integer(ik)  position in x
  pure subroutine unpack_values(buf, num_comp, nodes, add, x)

    real(rk), intent(in) :: buf(:)
    integer(ik), intent(in) :: num_comp
    integer(ik), intent(in) :: nodes(:)
    logical, intent(in) :: add
    real(rk), intent(inout) :: x(:)

    integer(ik) :: j
    integer(ik) :: c
    integer(ik) :: p

    do j = 1_ik, int(size(nodes), ik)
      do c = 1_ik, num_comp
        p = num_comp*(nodes(j) - 1_ik) + c
        if (add) then
          x(p) = x(p) + buf(num_comp*(j - 1_ik) + c)
        else
          x(p) = buf(num_comp*(j - 1_ik) + c)
        end if
      end do
    end do

  end subroutine unpack_values

  !> Fatal error if a nodal vector does not have num_comp num_nodes_loc entries.
  !>
  !> Arguments:
  !>   mesh      (in)  type(mesh_t)  this rank's mesh
  !>   x         (in)  real(rk)(:)   nodal vector
  !>   num_comp  (in)  integer(ik)   values per node
  !>   routine   (in)  character     calling routine
  !> Locals:
  !>   text  character(80)  sizes as text
  subroutine check_size(mesh, x, num_comp, routine)

    type(mesh_t), intent(in) :: mesh
    real(rk), intent(in) :: x(:)
    integer(ik), intent(in) :: num_comp
    character(len=*), intent(in) :: routine

    character(len=80) :: text

    if (num_comp >= 1_ik .and. size(x) == num_comp*mesh%num_nodes_loc) return
    write (text, '(a,i0,a,i0,a,i0)') 'vector size ', size(x), ', expected ', num_comp, &
      & ' x ', mesh%num_nodes_loc
    call internal_error(routine, trim(text) // ' (values per node x local nodes)')

  end subroutine check_size

end module halo_mod
