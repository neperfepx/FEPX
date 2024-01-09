! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module gather_scatter_mod

! Parallel version of gather/scatter routines. Originally developed by Dave
!   Mika (~1998). Modified by Don Boyce 9/2000---see end of file.

! deb changes (9/2000):
! Put routines in module structure.
! Fixed bug in PART_SCATTER_setup. Nodes to be passed which were part of only
!   one element on the local process were ignored. See section of code setting
!   up copy_address_1 and _2.

! Contains subroutines:
! part_gather: Gather takes global arrays of nodal values to arrays on each
!   element.
! part_gather_start: Starts part_gather
! part_gather_stop: Stops part_gather
! part_scatter: Accumulate elemental values into a global nodal array.
! part_scatter_setup: Prepare data structures for later scatter.
! int_sort:

! Notes:
! To use these routines, you need first call part_scatter_setup to set up the
!   trace structure associated with a connectivity array. then you can freely
!   use the part_scatter and part_gather.
! part_scatter has the possibly unexpected side effect of also modifying the
!   elemental array. This is because the elemental array is used as the send
!   buffer for communications and values are accumulated here before sending.
!   Also, the nodal array is not zeroed before use.

  use, intrinsic :: iso_fortran_env, only: rk => real64, real_kind_d => real64
  use general_mod
  use types_mod, only: trace
  use parallel_mod, only: par_quit

  implicit none

! On vel, there is a bug in mpi_waitany which returns 0-based indices
!   instead of 1-based indices for fortran CALLs.

! Public

  public :: part_gather
  public :: part_scatter
  public :: part_scatter_setup

! Private

  private :: part_gather_start
  private :: part_gather_stop
  private :: int_sort

contains

  subroutine part_gather(elts, y, connectivity, tr)

    !  Gather takes global arrays of nodal values to arrays on each
    !  element.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! elts: Elemental array
    ! y: Nodal array
    ! connectivity: Connectivity
    ! tr: Data structure interface with mpi

    type(trace) :: tr
    real(rk) :: elts(tr%dim1_sub:tr%dim1_sup, tr%elt_sub:tr%elt_sup)
    real(rk) :: y(tr%dim2_sub:tr%dim2_sup)
    integer :: connectivity(tr%dim1_sub:tr%dim1_sup, tr%elt_sub:tr%elt_sup)

    !---------------------------------------------------------------------------

    call part_gather_start(elts, y, connectivity, tr)
    call part_gather_stop(elts, y, connectivity, tr)

    return

  end subroutine part_gather

  !===========================================================================

  subroutine part_gather_start(elts, y, connectivity, tr)

    ! Starts part_gather

    !---------------------------------------------------------------------------

    ! Arguments:
    ! elts:
    ! y:
    ! connectivity:
    ! tr:

    type(trace) :: tr
    real(rk) :: elts(tr%dim1_sub:tr%dim1_sup, tr%elt_sub:tr%elt_sup)
    real(rk) :: y(tr%dim2_sub:tr%dim2_sup)
    integer :: connectivity(tr%dim1_sub:tr%dim1_sup, tr%elt_sub:tr%elt_sup)

    ! Locals:

    integer :: i
    integer :: j
    integer :: k
    integer :: ierr

    !---------------------------------------------------------------------------

    tr%nreceived = 0
    tr%ntotal = tr%nelt_traces ! this is number we are going to receive.

    do i = 1, tr%num_procs
      if (tr%n_to_move(i) .ne. 0) then
        tr%nreceived = tr%nreceived + 1

        call mpi_irecv(elts(tr%dim1_sub, tr%elt_sub), 1, tr%elt_trace(i), &
            & i - 1, i - 1, tr%gs_comm, tr%req(tr%nreceived), ierr)
      end if

      if (tr%nodes_to_move(i) .ne. 0) then
        tr%ntotal = tr%ntotal + 1
        call mpi_isend(y(tr%dim2_sub), 1, tr%node_trace(i), i - 1, &
            & tr%myid, tr%gs_comm, tr%req(tr%ntotal), ierr)
      end if
    end do

    ! Gather local quantities

    do j = 1, tr%nlocals
      k = tr%locals(j)

      do i = tr%dim1_sub, tr%dim1_sup
        elts(i, k) = y(connectivity(i, k))
      end do
    end do

    return

  end subroutine part_gather_start

  !===========================================================================

  subroutine part_gather_stop(elts, y, connectivity, tr)

    ! Stops part_gather

    !---------------------------------------------------------------------------

    include 'mpif-config.h'
    include 'mpif-constants.h'

    ! Arguments:
    ! elts:
    ! y:
    ! connectivity:
    ! tr

    type(trace) :: tr
    real(rk) :: elts(tr%dim1_sub:tr%dim1_sup, tr%elt_sub:tr%elt_sup)
    real(rk) :: y(tr%dim2_sub:tr%dim2_sup)
    integer :: connectivity(tr%dim1_sub:tr%dim1_sup, tr%elt_sub:tr%elt_sup)

    ! Locals:

    integer :: i
    integer :: j
    integer :: jj
    integer :: k
    integer :: l
    integer :: m
    integer :: n
    integer :: ierr
    integer :: mpistatus(mpi_status_size)

    !---------------------------------------------------------------------------

    ! Pick up the local nodes of the boundary elements
    ! These elements have some off-procesor nodes

    do j = tr%nlocals + 1, tr%elt_sup - tr%elt_sub + 1
      k = tr%locals(j)

      do i = tr%dim1_sub, tr%dim1_sup
        l = connectivity(i, k)

        if ((l .ge. tr%dim2_sub) .and. (l .le. tr%dim2_sup)) then
          elts(i, k) = y(l)
        end if
      end do
    end do

    ! Disseminate nodal quantities after checking for completion of sends

    do i = 1, tr%ntotal
      call mpi_waitany(tr%ntotal, tr%req, jj, mpistatus, ierr)

      if (jj .le. tr%nreceived) then
        k = mpistatus(mpi_source) + 1 ! Message from proc k is confirmed

        m = 1

        do j = 1, tr%min_n_to_move(k)
          n = m + tr%copy_ptr(k)

          do l = 1, tr%copies(j + tr%min_ptr(k)) - 1
            elts(tr%copy_address_1(l + n) + 1, tr%copy_address_2(l + n)) &
                & = elts(tr%copy_address_1(n) + 1, tr%copy_address_2(n))
          end do

          m = m + tr%copies(j + tr%min_ptr(k))
        end do
      end if
    end do

    return

  end subroutine part_gather_stop

  !===========================================================================

  subroutine part_scatter(y, elts, connectivity, tr)

    ! Accumulate elemental values into a global nodal array.

    ! note: The elemental array `elts' is altered in this routine, and this may
    !   not be expected behavior.
    !        edit F. Villette (06/2023) : the spurious behavior was fixed using a buffer elts_buff

    !---------------------------------------------------------------------------

    include 'mpif-handles.h'
    include 'mpif-config.h'
    include 'mpif-constants.h'

    ! Arguments:
    ! y: Nodal array
    ! elts: Elemental array to be scattered
    ! connectivity: Connectivity array
    ! tr: Trace structure

    type(trace) :: tr
    real(rk) :: y(tr%dim2_sub:tr%dim2_sup)
    real(rk) :: elts(tr%dim1_sub:tr%dim1_sup, tr%elt_sub:tr%elt_sup)
    integer :: connectivity(tr%dim1_sub:tr%dim1_sup, tr%elt_sub:tr%elt_sup)

    ! Locals:

    integer :: i
    integer :: j
    integer :: jj
    integer :: k
    integer :: l
    integer :: m
    integer :: ierr
    integer :: mpistatus(mpi_status_size)

    
    real(rk) :: elts_buff(tr%dim1_sub:tr%dim1_sup, tr%elt_sub:tr%elt_sup)

    !---------------------------------------------------------------------------

    y = 0.0d0

    tr%nreceived = 0

    elts_buff = elts

    ! Post accounts receivable.

    do i = 1, tr%num_procs
      if (tr%nodes_to_move(i) .ne. 0) then
        tr%nreceived = tr%nreceived + 1

        call mpi_irecv(tr%buffer(1 + tr%node_ptr(i)), &
            & tr%nodes_to_move(i), mpi_double_precision, i - 1, i - 1, &
            & tr%gs_comm, tr%req(tr%nreceived), ierr)
      end if
    end do

    ! Assemble before sending.

    tr%ntotal = tr%nreceived
    do i = 1, tr%num_procs
      if (tr%n_to_move(i) .ne. 0) then
        m = 1

        do jj = 1, tr%min_n_to_move(i)
          l = m + tr%copy_ptr(i)

          do k = 1, tr%copies(jj + tr%min_ptr(i)) - 1
            elts(tr%copy_address_1(l) + 1, tr%copy_address_2(l)) = &
                & elts(tr%copy_address_1(l) + 1, tr%copy_address_2(l)) + &
                & elts(tr%copy_address_1(k + l) + 1, &
                & tr%copy_address_2(k + l))
          end do

          m = m + tr%copies(jj + tr%min_ptr(i))
        end do

        tr%ntotal = tr%ntotal + 1
        call mpi_isend(elts(tr%dim1_sub, tr%elt_sub), 1, tr%elt_trace(i), &
            & i - 1, tr%myid, tr%gs_comm, tr%req(tr%ntotal), ierr)
      end if
    end do

    ! While waiting for sends to complete, scatter local quantities

    ! The guaranteed local elements:

    do k = 1, tr%nlocals
      j = tr%locals(k)

      do i = tr%dim1_sub, tr%dim1_sup
        y(connectivity(i, j)) = y(connectivity(i, j)) + elts(i, j)
      end do
    end do

    ! ------- cut here for part_scatter start (above) and stop (below) ---------

    ! You are guaranteed all the locals of y are updated

    ! Now scatter the local part of elements are on the boundary:

    do l = tr%nlocals + 1, tr%elt_sup - tr%elt_sub + 1
      j = tr%locals(l)

      do i = tr%dim1_sub, tr%dim1_sup
        k = connectivity(i, j)

        if ((k .ge. tr%dim2_sub) .and. (k .le. tr%dim2_sup)) then
          y(k) = y(k) + elts(i, j)
        end if
      end do
    end do

    ! Check for completion of sends

    do i = 1, tr%ntotal
      call mpi_waitany(tr%ntotal, tr%req, j, mpistatus, ierr)

      if (j .le. tr%nreceived) then
        k = mpistatus(mpi_source) + 1 ! received message from processor k-1

        do l = 1, tr%nodes_to_move(k)
          y(tr%node_address(l + tr%node_ptr(k))) = &
              & y(tr%node_address(l + tr%node_ptr(k))) + &
              & tr%buffer(l + tr%node_ptr(k))
        end do
      end if
    end do

    elts = elts_buff
    return

  end subroutine part_scatter

  !===========================================================================

  subroutine part_scatter_setup(dim1_sub, dim1_sup, dim2_sub, dim2_sup, &
      & elt_sub, elt_sup, connectivity, tr)

    ! Prepare data structures for later scatter.

    !---------------------------------------------------------------------------

    include 'mpif-handles.h'

    ! Arguments:
    ! dim1_sub:
    ! dim1_sup:
    ! dim2_sub:
    ! dim2_sup:
    ! elt_sub:
    ! elt_sup:
    ! connectivity:
    ! tr:

    integer, intent(in) :: dim1_sub
    integer, intent(in) :: dim1_sup
    integer, intent(in) :: dim2_sub
    integer, intent(in) :: dim2_sup
    integer, intent(in) :: elt_sub
    integer, intent(in) :: elt_sup
    integer, intent(in) :: connectivity(dim1_sub:dim1_sup, elt_sub:elt_sup)
    type(trace), intent(out) :: tr

    ! Locals:

    integer, allocatable :: buffer(:)
    integer, allocatable :: recvbuff(:)
    integer, allocatable :: addresses(:)
    integer, allocatable :: eltsizes(:)
    integer, allocatable :: dim2_to_processors(:)
    integer, allocatable :: order(:)
    integer, allocatable :: aux(:)
    integer, allocatable :: min_node_offsets(:)
    integer, allocatable :: min_elt_offsets(:)
    integer, allocatable :: displace(:)
    integer, allocatable :: dim2_subs(:)
    integer, allocatable :: dim2_sups(:)
    integer, allocatable :: recvbufsizes(:)
    integer, allocatable :: temp(:)
    integer, allocatable :: bufsizes(:)
    integer, allocatable :: sdispls(:)
    integer, allocatable :: node_offsets(:, :)
    integer, allocatable :: elt_offsets_1(:, :)
    integer, allocatable :: elt_offsets_2(:, :)
    integer :: ipoint
    integer :: position
    integer :: bufsize
    integer :: max_send
    integer :: local_nodes
    integer :: nonlocal_nodes
    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: max_nodes
    integer :: max_elts
    integer :: naux
    integer :: ierr
    integer :: np2
    integer :: temp_send(2)

    !----------------------------------------------------------------------

    call mpi_comm_size(mpi_comm_world, tr%num_procs, ierr)

    ! Get a new communicator for the gather/scatters.

    call mpi_cart_create(mpi_comm_world, 1, (/tr%num_procs/), (/.false./), &
        & .true., tr%gs_comm, ierr)

    ! Get my position in this communicator.

    call mpi_comm_rank(tr%gs_comm, tr%myid, ierr)

    ! Allocate storage for structure components and temporary local variables.

    allocate (displace(tr%num_procs))
    allocate (dim2_subs(tr%num_procs))
    allocate (dim2_sups(tr%num_procs))
    allocate (recvbufsizes(tr%num_procs))
    allocate (temp(2*tr%num_procs))
    allocate (bufsizes(tr%num_procs))
    allocate (sdispls(tr%num_procs))

    allocate (tr%n_to_move(tr%num_procs), stat=ierr)
    allocate (tr%nodes_to_move(tr%num_procs), stat=ierr)
    allocate (tr%min_n_to_move(tr%num_procs), stat=ierr)
    allocate (tr%locals(elt_sup - elt_sub + 1), stat=ierr)
    allocate (tr%node_trace(tr%num_procs), stat=ierr)
    allocate (tr%elt_trace(tr%num_procs), stat=ierr)

    if (ierr .ne. 0) then
      call par_quit('allocation failure')
    end if

    allocate (tr%min_ptr(tr%num_procs), stat=ierr)
    allocate (tr%copy_ptr(tr%num_procs), stat=ierr)
    allocate (tr%node_ptr(tr%num_procs), stat=ierr)

    np2 = 2*(tr%num_procs - 1)

    allocate (tr%req(np2), stat=ierr)

    tr%dim1_sub = dim1_sub
    tr%dim1_sup = dim1_sup
    tr%dim2_sub = dim2_sub
    tr%dim2_sup = dim2_sup
    tr%elt_sub = elt_sub
    tr%elt_sup = elt_sup

    ! This section appears very limiting. Each node allocates storage for all of
    !   the arrays on all of the nodes.

    i = (elt_sup - elt_sub + 1)*(tr%dim1_sup - tr%dim1_sub + 1)
    allocate (node_offsets(i, tr%num_procs))
    allocate (elt_offsets_1(i, tr%num_procs))
    allocate (elt_offsets_2(i, tr%num_procs))

    ! Get the decomposition boundaries from all processes.

    temp_send(1) = tr%dim2_sub
    temp_send(2) = tr%dim2_sup
    call mpi_allgather(temp_send, 2, mpi_integer, temp, 2, mpi_integer,&
        & tr%gs_comm, ierr)

    j = 1
    do i = 1, tr%num_procs
      dim2_subs(i) = temp(j)
      dim2_sups(i) = temp(j + 1)
      j = j + 2
    end do

    ! Find the traces for gather/scatter.

    ! Again, all nodes in the entire problem are allocated here, and then set to
    !   a constant value.

    allocate (dim2_to_processors(dim2_subs(1):dim2_sups(tr%num_procs)))

    do i = 1, tr%num_procs
      dim2_to_processors(dim2_subs(i):dim2_sups(i)) = i
    end do

    tr%n_to_move = 0
    tr%nlocals = 0
    nonlocal_nodes = 0

    do i = elt_sub, elt_sup
      local_nodes = 0

      do j = tr%dim1_sub, tr%dim1_sup
        ipoint = connectivity(j, i)

        if (ipoint .lt. 0) call par_quit('Error  :     > &
          &ipoint out of bounds')

        if (ipoint > dim2_sups(tr%num_procs)) then
          call par_quit('segv', .true.)
        end if

        if (dim2_to_processors(ipoint) - 1 .ne. tr%myid) then ! moving
          k = dim2_to_processors(ipoint)
          tr%n_to_move(k) = tr%n_to_move(k) + 1
          elt_offsets_1(tr%n_to_move(k), k) = j
          elt_offsets_2(tr%n_to_move(k), k) = i
          node_offsets(tr%n_to_move(k), k) = ipoint

        else ! node is local
          local_nodes = local_nodes + 1
        end if
      end do

      ! Check if all nodes were local.  In that case, the element is local

      if (local_nodes .eq. (tr%dim1_sup - tr%dim1_sub + 1)) then
        ! Element is local

        tr%nlocals = tr%nlocals + 1
        tr%locals(tr%nlocals) = i

      else
        ! Part of element is on another processor.
        ! Save these too, but starting from the top and counting down

        tr%locals(elt_sup - elt_sub + 1 - nonlocal_nodes) = i
        nonlocal_nodes = nonlocal_nodes + 1
      end if
    end do

    deallocate (dim2_to_processors)

    ! Each processor has tagged the nodes that it does not have

    ! To save space, use a new storage pattern from here out  (i.e. n_to_move
    ! could be large for one processor and zero for the next)

    i = sum(tr%n_to_move)

    allocate (min_elt_offsets(i))
    allocate (min_node_offsets(i))
    allocate (tr%copies(i))
    allocate (tr%copy_address_1(i))
    allocate (tr%copy_address_2(i))

    naux = maxval(tr%n_to_move)
    allocate (order(naux))
    allocate (aux(naux))

    ! Find the minimal set (we do not want to send duplicates)
    ! First sort the nodes to make finding the minimal set more efficient

    tr%min_n_to_move = 0
    tr%copies = 0
    k = 0
    l = 0

    do i = 1, tr%num_procs
      tr%min_ptr(i) = k
      tr%copy_ptr(i) = l
      if (tr%n_to_move(i) .ne. 0) then
        do j = 1, tr%n_to_move(i)
          order(j) = j
        end do

        ! This is not a "stable" sort

        call int_sort(node_offsets(:, i), tr%n_to_move(i), order)

        tr%min_n_to_move(i) = tr%min_n_to_move(i) + 1
        min_node_offsets(1 + tr%min_ptr(i)) = node_offsets(1, i)

        ! Switching min_elt_offsets to an absolute offset from the
        !   beginning of the buffer instead of a 1 and 2 index into the
        !   array. The routine mpi_type_indexed below will need it this way.

        min_elt_offsets(1 + tr%min_ptr(i)) = (tr%dim1_sup - tr%dim1_sub &
            & + 1)*(elt_offsets_2(order(1), i) - elt_sub) + &
            & elt_offsets_1(order(1), i) - 1

        j = 1
        l = l + 1
        tr%copy_address_1(l) = elt_offsets_1(order(j), i) - 1
        tr%copy_address_2(l) = elt_offsets_2(order(j), i)

        tr%copies(tr%min_n_to_move(i) + tr%min_ptr(i)) = &
          tr%copies(tr%min_n_to_move(i) + tr%min_ptr(i)) + 1

        do j = 2, tr%n_to_move(i)
          if (node_offsets(j, i) .eq. node_offsets(j - 1, i)) then
            ! The connectivities are the same.  In a scatter, the
            ! unassembled elems with the same global nodes will need to
            ! be added before sending, and in a gather, these will need
            ! to be assigned the same value after that value is received
            ! from another processor.

          else
            ! This is the first copy. The value in this address will be
            !   sent or assigned in the communication.

            tr%min_n_to_move(i) = tr%min_n_to_move(i) + 1
            min_node_offsets(tr%min_n_to_move(i) + tr%min_ptr(i)) = &
                & node_offsets(j, i)
            min_elt_offsets(tr%min_n_to_move(i) + tr%min_ptr(i)) = &
                & (tr%dim1_sup - tr%dim1_sub + 1)* &
                & (elt_offsets_2(order(j), i) - elt_sub) + &
                & elt_offsets_1(order(j), i) - 1
          end if

          tr%copies(tr%min_n_to_move(i) + tr%min_ptr(i)) = &
              & tr%copies(tr%min_n_to_move(i) + tr%min_ptr(i)) + 1

          l = l + 1
          tr%copy_address_1(l) = elt_offsets_1(order(j), i) - 1
          tr%copy_address_2(l) = elt_offsets_2(order(j), i)
        end do
      end if

      k = k + tr%min_n_to_move(i)
    end do

    deallocate (order, aux, elt_offsets_1, elt_offsets_2)

    ! Distribute the traces we just calculated but first determine size of
    !   needed buffer

    max_send = sum(tr%min_n_to_move) + tr%num_procs

    call mpi_pack_size(max_send, mpi_integer, tr%gs_comm, bufsize, ierr)

    allocate (buffer(bufsize/4)) ! bufsize is in bytes

    ! Pack data into buffer

    position = 0  ! Start at beginning of buffer

    do i = 1, tr%num_procs
      sdispls(i) = position
      j = tr%min_n_to_move(i)

      call mpi_pack(j, 1, mpi_integer, buffer, bufsize, position, &
          & tr%gs_comm, ierr)

      if (j .ne. 0) then
        call mpi_pack(min_node_offsets(1 + tr%min_ptr(i)), j, mpi_integer, &
            & buffer, bufsize, position, tr%gs_comm, ierr)
      end if

      bufsizes(i) = position - sdispls(i)
    end do

    deallocate (min_node_offsets)

    ! Distribute the buffersizes in preperation for data transfer

    call mpi_alltoall(bufsizes, 1, mpi_integer, recvbufsizes, 1, mpi_integer, &
        & tr%gs_comm, ierr)

    ! Allocate space for incomming data

    bufsize = sum(recvbufsizes)   ! bufsize is in bytes
    allocate (recvbuff(bufsize/4))

    j = 0

    do i = 1, tr%num_procs
      displace(i) = j
      j = j + recvbufsizes(i)
    end do

    ! distribute! (finally)

    call mpi_alltoallv(buffer, bufsizes, sdispls, mpi_packed, recvbuff, &
        & recvbufsizes, displace, mpi_packed, tr%gs_comm, ierr)

    deallocate (buffer)

    ! unpack data...
    ! get sizes first

    position = 0

    do i = 1, tr%num_procs
      call mpi_unpack(recvbuff, bufsize, position, tr%nodes_to_move(i), &
          & 1, mpi_integer, tr%gs_comm, ierr)

      position = position + tr%nodes_to_move(i)*4 ! position in bytes
    end do

    max_nodes = maxval(tr%nodes_to_move)
    max_elts = maxval(tr%min_n_to_move)
    deallocate (node_offsets)
    allocate (node_offsets(max_nodes, tr%num_procs))
    allocate (eltsizes(Max(max_nodes, max_elts)))
    allocate (addresses(max_nodes))

    position = 0

    do i = 1, tr%num_procs
      position = position + 4  ! position is in bytes
      j = tr%nodes_to_move(i)

      if (j .ne. 0) then
        call mpi_unpack(recvbuff, bufsize, position, node_offsets(1, i), &
            & j, mpi_integer, tr%gs_comm, ierr)
      end if
    end do

    deallocate (recvbuff)

    ! Commit the traces

    eltsizes = 1

    tr%nelt_traces = 0
    tr%nnode_traces = 0

    do i = 1, tr%num_procs
      if (tr%nodes_to_move(i) .ne. 0) then
        tr%nnode_traces = tr%nnode_traces + 1

        do j = 1, tr%nodes_to_move(i)
          addresses(j) = node_offsets(j, i) - tr%dim2_sub
        end do

        call mpi_type_indexed(tr%nodes_to_move(i), eltsizes, &
            & addresses, mpi_double_precision, tr%node_trace(i), ierr)

        call mpi_type_commit(tr%node_trace(i), ierr)
      end if

      if (tr%min_n_to_move(i) .ne. 0) then
        tr%nelt_traces = tr%nelt_traces + 1

        call mpi_type_indexed(tr%min_n_to_move(i), eltsizes, &
            & min_elt_offsets(1 + tr%min_ptr(i)), &
            & mpi_double_precision, tr%elt_trace(i), ierr)

        call mpi_type_commit(tr%elt_trace(i), ierr)
      end if
    end do

    ! Since 1 does not provide a recv with an add, it is inefficient to use
    !   mpi_type_indexed with the part_scatter.  Instead, use a buffer to hold
    !   the incomming message and add message to appropriate locations.
    !   Hopefully the needed function will be added to the 1 standard - then the
    !   above node_trace can be used.

    bufsize = sum(tr%nodes_to_move)
    allocate (tr%node_address(bufsize))
    allocate (tr%buffer(bufsize))

    j = 0

    do i = 1, tr%num_procs
      tr%node_ptr(i) = j
      j = j + tr%nodes_to_move(i)
    end do

    do i = 1, tr%num_procs
      do j = 1, tr%nodes_to_move(i)
        tr%node_address(j + tr%node_ptr(i)) = node_offsets(j, i)
      end do
    end do

    deallocate (eltsizes, addresses, node_offsets, min_elt_offsets)

    deallocate (displace, dim2_subs, dim2_sups, recvbufsizes, temp, bufsizes, &
        & sdispls)

    tr%setup = 0

    return

  end subroutine part_scatter_setup

  !===========================================================================

  subroutine int_sort(a, n, t)

    ! Version 1.00, 14 July 1995
    !Author: Alan Miller
    !        csiro Division of Mathematics & Statistics
    !        Private Bag 10, Rosebank mdc
    !        Clayton 3169, Victoria, Australia
    !        Phone: (+61) 3 9545-8036      Fax: (+61) 3 9545-8080
    !        e-mail: Alan.Miller @ mel.dms.csiro.au

    ! non-recursive stack version of quicksort from n.wirth's pascal
    ! book, 'algorithms + data structures = programs'.

    ! single precision, also changes the order of the associated array t.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! a:
    ! n:
    ! t:

    integer :: a(n)
    integer, intent(in) :: n
    integer, intent(inout) :: t(n)

    ! Locals:

    integer :: i
    integer :: j
    integer :: k
    integer :: l
    integer :: r
    integer :: s
    integer :: stackl(15)
    integer :: stackr(15)
    integer :: ww
    integer :: w
    integer :: x

    !---------------------------------------------------------------------------

    s = 1
    stackl(1) = 1
    stackr(1) = n

    ! keep taking the top request from the stack until s = 0.

10  continue

    l = stackl(s)
    R = stackr(s)
    s = s - 1

    ! keep splitting a(l), ... , a(r) until l >= r.

20  continue

    i = l
    j = r
    k = (l + r)/2
    x = a(k)

    !     repeat until i > j.

    do
      do
        if (a(i) .lt. x) then ! Search from lower end
          i = i + 1

          cycle

        else
          exit
        end if
      end do

      do
        if (x .lt. a(j)) then ! Search from upper end
          j = j - 1

          cycle

        else
          exit
        end if
      end do

      if (i .le. j) then ! Swap positions i & j
        w = a(i)
        ww = t(i)
        a(i) = a(j)
        t(i) = t(j)
        a(j) = w
        t(j) = ww
        i = i + 1
        j = j - 1

        if (i .gt. j) exit

      else
        exit
      end if
    end do

    if (j - l .ge. r - i) then
      if (l .lt. j) then
        s = s + 1
        stackl(s) = l
        stackr(s) = j
      end if

      l = i

    else
      if (i .lt. r) then
        s = s + 1
        stackl(s) = i
        stackr(s) = r
      end if

      R = j
    end if

    if (l .lt. r) go to 20
    if (s .ne. 0) go to 10

    return

  end subroutine int_sort

end module gather_scatter_mod
