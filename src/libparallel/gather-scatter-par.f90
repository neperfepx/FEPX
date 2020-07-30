! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE gather_scatter
  !
  !  Parallel version of gather/scatter routines.
  !  Originally developed by Dave Mika.
  !
  !  Modified by Don Boyce 9/2000---see end of file.
  !
  !
  USE IntrinsicTypesModule, RK=>REAL_KIND
  USE parallel_mod,  ONLY: par_quit, myid, myidstr
  !
  IMPLICIT NONE
  !
  TYPE trace
     !
     !  buffer : data to be received
     !
     REAL(RK),  ALLOCATABLE, DIMENSION(:) :: buffer
     !
     INTEGER, ALLOCATABLE, DIMENSION(:) :: copies, copy_address_1, copy_address_2
     INTEGER, ALLOCATABLE, DIMENSION(:) :: node_address
     INTEGER, ALLOCATABLE, DIMENSION(:) :: locals
     INTEGER, ALLOCATABLE, DIMENSION(:) :: n_to_move, nodes_to_move, min_n_to_move
     INTEGER, ALLOCATABLE, DIMENSION(:) :: node_trace, elm_trace
     INTEGER, ALLOCATABLE, DIMENSION(:) :: min_ptr, copy_ptr, node_ptr
     INTEGER, ALLOCATABLE, DIMENSION(:) :: req
     !
     INTEGER  dim1_sub1, dim1_sup1
     INTEGER  dim2_sub1, dim2_sup1
     INTEGER  el_sub1,   el_sup1
     INTEGER  nlocals, nreceived
     INTEGER  nelem_traces, nnode_traces
     INTEGER  ntotal, setup
     !
     !  numprocs : number of processes
     !  myid     : rank of this process
     !  gs_comm  : communicator for gather/scatter
     !
     INTEGER  numprocs, myid, gs_comm
     !
  END TYPE trace
  !
  !  On velocity, there is a bug in MPI_Waitany which returns 0-based indices
  !  instead of 1-based indices for FORTRAN calls.
  !
  INTEGER, PRIVATE, PARAMETER :: WAITANY_START = 1   ! should be one 
  !
  !---------------------ROUTINES-----------------------------------------
  !
  !  Accessible Routines.
  !
  !  To use these routines, you need first call part_scatter_setup to set up
  !  the trace structure associated with a connectivity array.  Then you 
  !  you can freely use the part_scatter and part_gather.
  !
  !  NOTES:
  !  
  !  (*)  part_scatter has the possibly unexpected side effect of also modifying
  !  .    the elemental array; this is because the elemental array is used as the
  !  .    send buffer for communications and values are accumulated here before
  !  .    sending;  also, the nodal array is not zeroed before use;
  !
  PUBLIC part_gather, part_scatter, part_scatter_setup
  !
  !  part_scatter_setup -- sets up trace data structure for use in gather/scatter
  !  part_gather        -- assigns nodal values to an elemental array
  !  part_scatter       -- accumulates elemental values into a nodal array
  !
  !  Inaccessible Routines.
  !
  PRIVATE part_gather_start, part_gather_stop, int_sort
  !
  !  part_gather_start  -- 
  !  part_gather_stop   --
  !  int_sort           -- NOTE: Dave must have gotten this from the web; there
  !  .                  --       is no information about right to use/license, etc.
  !
CONTAINS
  !
  ! *********************************************************************
  !
  SUBROUTINE part_gather(elms, y, connectivity, tr)
    !
    !  Gather takes global arrays of nodal values to arrays on each
    !  element.
    !
    IMPLICIT NONE
    !
    !  Arguments:
    !
    !  elms         : elemental array
    !  y            : nodal array
    !  connectivity : connectivity
    !  trace        : data structure interface with MPI
    !
    TYPE (trace) tr
    !
    REAL(RK)  elms(tr%dim1_sub1:tr%dim1_sup1,tr%el_sub1:tr%el_sup1)
    REAL(RK)  y(tr%dim2_sub1:tr%dim2_sup1)
    !
    INTEGER connectivity(tr%dim1_sub1:tr%dim1_sup1,tr%el_sub1:tr%el_sup1)

    call  part_gather_start(elms, y, connectivity, .False., tr)
    call  part_gather_stop (elms, y, connectivity, .False., tr)

    RETURN
    !  
  END SUBROUTINE part_gather
  !
  ! *********************************************************************
  !
  SUBROUTINE part_gather_start(elms, y, connectivity, persist, tr)
    !
    !     Author: David Mika (current with GE)
    !
    IMPLICIT NONE
    !
    !  Arguments.
    !
    TYPE (trace) tr
    !
    REAL(RK)  elms(tr%dim1_sub1:tr%dim1_sup1,tr%el_sub1:tr%el_sup1)
    REAL(RK)  y(tr%dim2_sub1:tr%dim2_sup1)
    !
    INTEGER connectivity(tr%dim1_sub1:tr%dim1_sup1,tr%el_sub1:tr%el_sup1)
    !
    LOGICAL persist
    !
    !  Locals.
    !
    INTEGER i, j, k, ierr
    !
    !----------------------------------------------------------------------
    !
    tr%nreceived = 0
    tr%ntotal = tr%nelem_traces ! this is number we are going to receive.
    Do i = 1, tr%numprocs
       If(tr%n_to_move(i) .Ne. 0) Then
          tr%nreceived = tr%nreceived + 1
          Call MPI_IRecv(elms(tr%dim1_sub1,tr%el_sub1), 1, tr%elm_trace(i), &
               i-1, i-1, tr%gs_comm, tr%req(tr%nreceived), ierr)
       Endif

       If(tr%nodes_to_move(i) .Ne. 0) Then
          tr%ntotal = tr%ntotal + 1
          Call MPI_ISend(y(tr%dim2_sub1), 1, tr%node_trace(i), i-1, &
               tr%myid, tr%gs_comm, tr%req(tr%ntotal), ierr)
       Endif
    Enddo

    ! gather local quantities 

    Do j = 1, tr%nlocals
       k = tr%locals(j)
       Do i = tr%dim1_sub1, tr%dim1_sup1
          elms(i,k) = y(connectivity(i,k))
       Enddo
    Enddo

    RETURN
  END SUBROUTINE part_gather_start
  !
  ! *********************************************************************
  !
  SUBROUTINE part_gather_stop(elms, y, connectivity, persist, tr)
    !
    !     Author: David Mika (currently with GE)
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INCLUDE 'mpif-config.h'
    INCLUDE 'mpif-constants.h'
    !
    !  Arguments.
    !
    TYPE (trace) tr
    !
    REAL(RK)  elms(tr%dim1_sub1:tr%dim1_sup1,tr%el_sub1:tr%el_sup1)
    REAL(RK)  y(tr%dim2_sub1:tr%dim2_sup1)
    !
    INTEGER connectivity(tr%dim1_sub1:tr%dim1_sup1,tr%el_sub1:tr%el_sup1)
    !
    LOGICAL persist
    !
    !  Locals.
    !
    INTEGER  i, j, jj, k, l, m, n, ierr
    INTEGER  mpistatus(MPI_STATUS_SIZE)
    !
    !----------------------------------------------------------------------
    !
    ! pick up the local nodes of the boundary elements
    ! These elements have some off-procesor nodes

    Do j = tr%nlocals+1, tr%el_sup1-tr%el_sub1+1
       k = tr%locals(j)
       Do i = tr%dim1_sub1, tr%dim1_sup1
          l = connectivity(i,k)
          If((l.Ge.tr%dim2_sub1).And.(l.Le.tr%dim2_sup1)) Then
             elms(i,k) = y(l)
          Endif
       Enddo
    Enddo

    ! disseminate nodal quantities after checking for completion of sends
    Do i = 1, tr%ntotal

       Call MPI_Waitany(tr%ntotal, tr%req, jj, mpistatus, ierr)
       jj = jj + 1 - WAITANY_START            !  fixes bug on velocity

       If( jj .Le. tr%nreceived ) Then

          k = mpistatus(MPI_SOURCE) + 1 ! message from proc k is confirmed

          m = 1
          Do j = 1, tr%min_n_to_move(k)
             n = m + tr%copy_ptr(k)
             Do l = 1, tr%copies(j + tr%min_ptr(k))-1
                elms(tr%copy_address_1( l + n ), tr%copy_address_2( l + n ))  &
                     = elms(tr%copy_address_1( n ), tr%copy_address_2( n ))
             End Do
             m = m + tr%copies(j + tr%min_ptr(k))
          End Do

       Endif
    Enddo

    RETURN
  END SUBROUTINE part_gather_stop
  !
  ! *********************************************************************
  !
  SUBROUTINE part_scatter(y, elms, connectivity, persist, tr)
    !
    ! Accumulate elemental values into a global nodal array.
    !
    ! NOTE:  the elemental array `elms' is altered in this routine, and
    ! .   :  this may not be expected behavior.
    !
    ! Author: David Mika (currently with GE)
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INCLUDE 'mpif-handles.h'
    INCLUDE 'mpif-config.h'
    INCLUDE 'mpif-constants.h'
    !
    !  Arguments.
    !
    !  y            : nodal array                                  (output)
    !  elms         : elemental array to be scattered        (input/output)
    !  connectivity : connectivity array                            (input)
    !  persist      : flag for persistent communication             (input)
    !  tr           : trace structure                        (input/output)
    !        
    TYPE(trace) tr
    !
    REAL(RK)  elms(tr%dim1_sub1:tr%dim1_sup1,tr%el_sub1:tr%el_sup1)
    REAL(RK)  y(tr%dim2_sub1:tr%dim2_sup1)
    !
    INTEGER connectivity(tr%dim1_sub1:tr%dim1_sup1,tr%el_sub1:tr%el_sup1)
    !
    LOGICAL persist
    !
    !  Locals.
    !
    INTEGER i, j, jj, k, l, m, ierr
    !
    INTEGER  mpistatus(MPI_STATUS_SIZE)
    INTEGER  status_array(MPI_STATUS_SIZE,tr%numprocs*2)
    !
    !----------------------------------------------------------------------
    !
    tr%nreceived = 0
    !
    !  Post accounts receivable.
    !  
    Do i = 1, tr%numprocs
       If(tr%nodes_to_move(i) .Ne. 0) Then
          tr%nreceived = tr%nreceived + 1
          Call MPI_IRecv(tr%buffer(1+tr%node_ptr(i)), tr%nodes_to_move(i), &
               MPI_DOUBLE_PRECISION, i-1, i-1, tr%gs_comm, tr%req(tr%nreceived), ierr)
       Endif
    Enddo

    !
    !  Assemble before sending.
    !  
    tr%ntotal = tr%nreceived
    Do i = 1, tr%numprocs
       If(tr%n_to_move(i) .Ne. 0) Then

          m = 1
          Do jj = 1, tr%min_n_to_move(i)
             l = m + tr%copy_ptr(i)
             Do k = 1, tr%copies(jj + tr%min_ptr(i))-1
                elms(tr%copy_address_1( l ),  tr%copy_address_2( l ) ) = &
                     elms(tr%copy_address_1( l ), tr%copy_address_2( l ) ) + & 
                     elms(tr%copy_address_1( k + l ), tr%copy_address_2( k + l ) )
             End Do
             m = m + tr%copies(jj + tr%min_ptr(i))
          End Do

          tr%ntotal = tr%ntotal + 1
          Call MPI_ISend(elms(tr%dim1_sub1,tr%el_sub1),1,tr%elm_trace(i), &
               i-1,tr%myid,tr%gs_comm,tr%req(tr%ntotal),ierr) 
       Endif
    Enddo

    ! While waiting for sends to complete, scatter local quantities  

    ! the guaranteed local elements:

    Do k = 1, tr%nlocals
       j = tr%locals(k) 
       Do i = tr%dim1_sub1, tr%dim1_sup1
          y(connectivity(i,j)) = y(connectivity(i,j)) + elms(i,j)
       Enddo
    Enddo

    ! ------- cut here for part_scatter start (above) and stop (below) ----------
    ! you are guaranteed all the locals of y are updated

    ! now scatter the local part of elements are on the boundary: 
    
    Do l = tr%nlocals+1, tr%el_sup1-tr%el_sub1+1
       j = tr%locals(l)
       Do i = tr%dim1_sub1, tr%dim1_sup1
          k = connectivity(i,j)
          If ((k .Ge. tr%dim2_sub1) .And. (k .Le. tr%dim2_sup1)) Then
             y(k) = y(k) + elms(i,j)
          Endif
       Enddo
    Enddo

    ! check for completion of sends
    Do i = 1, tr%ntotal

       Call MPI_Waitany(tr%ntotal,tr%req,j,mpistatus,ierr)
       j = j + 1 - WAITANY_START            !  fixes bug on velocity

       If( j .Le. tr%nreceived ) Then
          k = mpistatus(MPI_SOURCE) + 1 ! received message from processor k-1

          Do l = 1, tr%nodes_to_move(k)
             y(tr%node_address(l+tr%node_ptr(k))) = &
                  y(tr%node_address(l+tr%node_ptr(k))) + &
                  tr%buffer(l+tr%node_ptr(k))
          End Do

       Endif
    Enddo

    RETURN
  END SUBROUTINE part_scatter
  !
  ! *********************************************************************
  !
  SUBROUTINE part_scatter_setup(dim1_sub1, dim1_sup1, dim2_sub1, dim2_sup1, &
       el_sub1, el_sup1, connectivity, tr)
    !
    ! Prepare data structures for later scatter.
    !
    ! Author: David Mika (now with GE)
    !
    !----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INCLUDE 'mpif-handles.h'
    !
    !  Arguments.
    !
    TYPE (trace) tr
    !
    INTEGER  dim1_sub1, dim1_sup1, dim2_sub1, dim2_sup1, el_sub1, el_sup1
    INTEGER  connectivity(dim1_sub1:dim1_sup1,el_sub1:el_sup1)
    !
    INTENT(IN) :: dim1_sub1, dim1_sup1, dim2_sub1, dim2_sup1, el_sub1, el_sup1
    INTENT(IN) :: connectivity
    INTENT(OUT) :: tr
    !
    !  Locals.
    !
    INTEGER, ALLOCATABLE, DIMENSION(:) :: buffer, recvbuff, addresses, &
         elmsizes, dim2_to_processors, order, aux,  min_node_offsets, &
         min_elem_offsets, displace, dim2_sub1s, dim2_sup1s, &
         recvbufsizes, temp, bufsizes, sdispls

    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: node_offsets, elem_offsets_1, &
         elem_offsets_2

    INTEGER  ipoint, position, bufsize, max_send, &
         local_nodes, nonlocal_nodes, &
         i, j, k, l, max_nodes, max_elems, naux, &
         el_sub, el_sup, ierr, np2, temp_send(2)
    !
    !----------------------------------------------------------------------
    !
    call MPI_COMM_SIZE( MPI_COMM_WORLD, tr%numprocs, ierr )
    !
    !  Get a new communicator for the gather/scatters.
    !  
    CALL MPI_CART_CREATE( MPI_COMM_WORLD, 1, (/ tr%numprocs /), (/ .FALSE. /), &
         .true., tr % gs_comm, ierr )       
    !
    !  Get my position in this communicator.
    !  
    call MPI_COMM_RANK( tr%gs_comm, tr%myid, ierr )
    !
    !  Allocate storage for structure components and temporary local variables.
    !
    el_sub = el_sub1 + 1
    el_sup = el_sup1 + 1
    !
    Allocate (&
         &  displace  (tr%numprocs),        dim2_sub1s  (tr%numprocs), &
         &  dim2_sup1s(tr%numprocs),        recvbufsizes(tr%numprocs), &
         &  temp    (2*tr%numprocs),        bufsizes    (tr%numprocs), &
         &  sdispls   (tr%numprocs) &
         &  )

    Allocate (&
         &  tr%n_to_move    (tr%numprocs),      tr%nodes_to_move(tr%numprocs),&
         &  tr%min_n_to_move(tr%numprocs),      tr%locals    (el_sup-el_sub1),&
         &  tr%node_trace   (tr%numprocs),      tr%elm_trace    (tr%numprocs),&
         &  STAT=ierr )
    IF (ierr /= 0) THEN
      call par_quit('allocation failure')
    END IF

    Allocate (&
         &  tr%min_ptr      (tr%numprocs),      tr%copy_ptr     (tr%numprocs),&
         &  STAT=ierr )

    Allocate (&
         &  tr%node_ptr     (tr%numprocs),&
         &  STAT=ierr )



    np2 = 2*(tr%numprocs-1)

    ALLOCATE (tr%req(np2), STAT=ierr)

    tr%dim1_sub1 = dim1_sub1
    tr%dim1_sup1 = dim1_sup1
    tr%dim2_sub1 = dim2_sub1
    tr%dim2_sup1 = dim2_sup1
    tr%el_sub1   = el_sub1
    tr%el_sup1   = el_sup1
    !
    !  This section appears very limiting.  Each node allocates storage for
    !  all of the arrays on all of the nodes.
    !
    i = (el_sup-el_sub1)*(tr%dim1_sup1-tr%dim1_sub1+1)
    Allocate (node_offsets  (i,tr%numprocs))
    Allocate (elem_offsets_1(i,tr%numprocs))
    Allocate (elem_offsets_2(i,tr%numprocs))
    !
    !  Get the decomposition boundaries from all processes.
    !
    temp_send(1) = tr%dim2_sub1
    temp_send(2) = tr%dim2_sup1 
    Call MPI_Allgather(&
         &  temp_send, 2, MPI_INTEGER,&      ! send
         &  temp, 2, MPI_INTEGER,&      ! receive
         &  tr%gs_comm, ierr)
    !
    j = 1
    Do i = 1, tr%numprocs
       dim2_sub1s(i) = temp(j)
       dim2_sup1s(i) = temp(j+1)
       j = j + 2
    Enddo
    !
    !  Find the traces for gather/scatter.
    !
    !  Again, all nodes in the entire problem are allocated here,
    !  and then set to a constant value.
    !  
    Allocate (dim2_to_processors(dim2_sub1s(1):dim2_sup1s(tr%numprocs)))
    Do i = 1, tr%numprocs
       dim2_to_processors(dim2_sub1s(i):dim2_sup1s(i)) = i
    Enddo
    tr%n_to_move   = 0
    tr%nlocals     = 0
    nonlocal_nodes = 0
    !
    Do i = el_sub1, el_sup1
       local_nodes = 0
       Do j = tr%dim1_sub1, tr%dim1_sup1
          ipoint = connectivity(j,i)
          IF (ipoint > dim2_sup1s(tr%numprocs) ) THEN
            CALL par_quit('segv', .TRUE.)
          END IF
          If (dim2_to_processors(ipoint)-1 .Ne. tr%myid) Then       ! moving
             k = dim2_to_processors(ipoint)
             tr%n_to_move(k) = tr%n_to_move(k) + 1
             elem_offsets_1(tr%n_to_move(k),k) = j
             elem_offsets_2(tr%n_to_move(k),k) = i
             node_offsets(tr%n_to_move(k),k) = ipoint
          Else                                                      ! node is local
             local_nodes = local_nodes + 1
          Endif
       Enddo
       ! check if all nodes were local.  In that case, the element is local
       If(local_nodes .Eq. (tr%dim1_sup1 - tr%dim1_sub1 + 1)) Then  
                                        ! element is local
          tr%nlocals = tr%nlocals + 1
          tr%locals(tr%nlocals) = i
       Else                             ! part of element is on another processor.  
                         ! Save these too, but starting from the top and counting down
          tr%locals(el_sup-el_sub1-nonlocal_nodes) = i
          nonlocal_nodes = nonlocal_nodes + 1
       Endif
    Enddo
    Deallocate(dim2_to_processors)

    ! Each processor has tagged the nodes that it does not have

    ! To save space, use a new storage pattern from here out  (i.e. n_to_move 
    ! could be large for one processor and zero for the next)

    i = Sum(tr%n_to_move)

    Allocate (min_elem_offsets(i))
    Allocate (min_node_offsets(i))
    Allocate (tr%copies(i))
    Allocate (tr%copy_address_1(i))
    Allocate (tr%copy_address_2(i))

    naux = Maxval(tr%n_to_move)
    Allocate (order(naux))
    Allocate (aux(naux))

    ! Find the minimal set (we do not want to send duplicates)
    ! First sort the nodes to make finding the minimal set more efficient

    tr%min_n_to_move = 0
    tr%copies = 0
    k = 0
    l = 0
    Do i = 1, tr%numprocs

       tr%min_ptr(i) = k
       tr%copy_ptr(i) = l
       If(tr%n_to_move(i) .Ne. 0) Then

          Do j = 1, tr%n_to_move(i)
             order(j) = j
          Enddo
          ! This is not a "stable" sort
          Call int_sort(node_offsets(:, i), tr%n_to_move(i), order) 

          tr%min_n_to_move(i) = tr%min_n_to_move(i) + 1
          min_node_offsets(1+tr%min_ptr(i)) = node_offsets(1,i)

          ! switching min_elem_offsets to an absolute offset from the beginning
          ! of the buffer instead of a 1 and 2 index into the array.  The routine
          ! MPI_Type_indexed below will need it this way.

          min_elem_offsets(1+tr%min_ptr(i)) = (tr%dim1_sup1 - tr%dim1_sub1 + 1) &
               * (elem_offsets_2(order(1),i) - el_sub1) + elem_offsets_1(order(1),i)
          
          j = 1
          l = l + 1
          tr%copy_address_1(l) = elem_offsets_1(order(j),i)
          tr%copy_address_2(l) = elem_offsets_2(order(j),i)

          tr%copies(tr%min_n_to_move(i)+tr%min_ptr(i)) = &
               tr%copies(tr%min_n_to_move(i)+tr%min_ptr(i)) + 1

          Do j = 2, tr%n_to_move(i)

             If(node_offsets(j,i) .Eq. node_offsets(j-1,i)) Then

                ! The connectivities are the same.  In a scatter, the unassembled
                ! elems with the same global nodes will need to be added before
                ! sending, and in a gather, these will need to be assigned the
                ! same value after that value is received from another processor.

             Else 

                ! This is the first copy.  
                ! The value in this address will be sent or assigned in
                ! the communication.

                tr%min_n_to_move(i) = tr%min_n_to_move(i) + 1
                min_node_offsets(tr%min_n_to_move(i)+tr%min_ptr(i)) = node_offsets(j,i)
                min_elem_offsets(tr%min_n_to_move(i)+tr%min_ptr(i)) = &
                     (tr%dim1_sup1 - tr%dim1_sub1 + 1) * &
                     (elem_offsets_2(order(j),i) - el_sub1) + &
                     elem_offsets_1(order(j),i)

             Endif

             tr%copies(tr%min_n_to_move(i)+tr%min_ptr(i)) = &
                  tr%copies(tr%min_n_to_move(i)+tr%min_ptr(i)) + 1

             l = l + 1
             tr%copy_address_1(l) = elem_offsets_1(order(j),i)
             tr%copy_address_2(l) = elem_offsets_2(order(j),i)

          Enddo
       Endif
       k = k + tr%min_n_to_move(i)
    Enddo


    Deallocate (order, aux, elem_offsets_1, elem_offsets_2)

    ! distribute the traces we just calculated
    ! but first determine size of needed buffer

    max_send = Sum(tr%min_n_to_move)+tr%numprocs
    Call MPI_Pack_Size(max_send, MPI_INTEGER, tr%gs_comm, bufsize, ierr)
    Allocate(buffer(bufsize/4)) ! bufsize is in bytes

    ! Pack data into buffer
    position = 0  ! Start at beginning of buffer

    Do i = 1, tr%numprocs
       sdispls(i) = position
       j = tr%min_n_to_move(i)
       Call MPI_Pack(j, 1, MPI_INTEGER, buffer, bufsize, &
            position, tr%gs_comm, ierr)
       If( j .Ne. 0 ) Then
          Call MPI_Pack(min_node_offsets(1+tr%min_ptr(i)), j, MPI_INTEGER, &
               buffer, bufsize, position, tr%gs_comm, ierr)
       Endif
       bufsizes(i) = position - sdispls(i)
    Enddo

    Deallocate (min_node_offsets)

    ! distribute the buffersizes in preperation for data transfer
    Call MPI_Alltoall(bufsizes, 1, MPI_INTEGER, recvbufsizes, &
         1, MPI_INTEGER, tr%gs_comm, ierr)

    ! allocate space for incomming data
    bufsize = Sum(recvbufsizes)   ! bufsize is in bytes
    Allocate(recvbuff(bufsize/4))

    j = 0
    Do i = 1, tr%numprocs
       displace(i) = j
       j = j + recvbufsizes(i)
    Enddo

    ! distribute! (finally)

    Call MPI_Alltoallv(buffer, bufsizes, sdispls, MPI_PACKED, recvbuff, &
         recvbufsizes, displace, MPI_PACKED, tr%gs_comm, ierr) 
    Deallocate(buffer)

    ! unpack data...
    ! get sizes first

    position = 0
    Do i = 1, tr%numprocs
       Call MPI_UNPack(recvbuff, bufsize, position, tr%nodes_to_move(i), &
            1, MPI_INTEGER, tr%gs_comm, ierr)
       position = position + tr%nodes_to_move(i)*4 ! position in bytes
    Enddo

    max_nodes = Maxval(tr%nodes_to_move)
    max_elems = Maxval(tr%min_n_to_move)
    Deallocate(node_offsets)
    Allocate(node_offsets(max_nodes, tr%numprocs))
    Allocate(elmsizes(Max(max_nodes, max_elems)))
    Allocate(addresses(max_nodes))

    position = 0
    Do i = 1, tr%numprocs
       position = position + 4  ! position is in bytes
       j = tr%nodes_to_move(i)
       If( j .Ne. 0 ) Then
          Call MPI_UNPack(recvbuff, bufsize, position, node_offsets(1,i), &
               j, MPI_INTEGER, tr%gs_comm, ierr)
       Endif
    Enddo
    Deallocate(recvbuff)

    ! commit the traces

    elmsizes = 1

    tr%nelem_traces = 0; tr%nnode_traces = 0
    Do i = 1, tr%numprocs

       If(tr%nodes_to_move(i) .Ne. 0) Then 
          tr%nnode_traces = tr%nnode_traces + 1
          Do j = 1, tr%nodes_to_move(i)
             addresses(j) = node_offsets(j,i) - tr%dim2_sub1
          Enddo
          Call MPI_Type_indexed(tr%nodes_to_move(i), elmsizes, addresses, &
               MPI_DOUBLE_PRECISION, tr%node_trace(i), ierr)
          Call MPI_Type_commit( tr%node_trace(i), ierr )
       Endif

       If(tr%min_n_to_move(i) .Ne. 0) Then
          tr%nelem_traces = tr%nelem_traces + 1
          Call MPI_Type_indexed(tr%min_n_to_move(i), elmsizes, &
               min_elem_offsets(1+tr%min_ptr(i)), MPI_DOUBLE_PRECISION, &
               tr%elm_trace(i), ierr )

          Call MPI_Type_commit( tr%elm_trace(i), ierr )
       Endif

    Enddo

    ! Since 1 does not provide a Recv with an ADD, it is inefficient to use 
    ! MPI_Type_indexed with the part_scatter.  Instead, use a buffer to hold
    ! the incomming message and add message to appropriate locations.  Hopefully
    ! the needed function will be added to the 1 standard - then the above 
    ! node_trace can be used. 

    bufsize = Sum(tr%nodes_to_move)
    Allocate(tr%node_address(bufsize))
    Allocate(tr%buffer(bufsize))

    j = 0
    Do i = 1, tr%numprocs
       tr%node_ptr(i) = j
       j = j + tr%nodes_to_move(i)
    End Do

    Do i = 1, tr%numprocs
       Do j = 1, tr%nodes_to_move(i)
          tr%node_address(j+tr%node_ptr(i)) = node_offsets(j,i)
       End Do
    End Do

    Deallocate(elmsizes, addresses, node_offsets, min_elem_offsets)

    DEALLOCATE(displace, dim2_sub1s, dim2_sup1s, recvbufsizes, temp, bufsizes, sdispls)

    tr%setup = 0

    RETURN
  END SUBROUTINE part_scatter_setup
  !
  ! *********************************************************************
  !
  SUBROUTINE int_sort(a, n, t)

    !     Version 1.00, 14 July 1995
    !     Author: Alan Miller
    !             CSIRO Division of Mathematics & Statistics
    !             Private Bag 10, Rosebank MDC
    !             Clayton 3169, Victoria, Australia
    !     Phone: (+61) 3 9545-8036      Fax: (+61) 3 9545-8080
    !     e-mail: Alan.Miller @ mel.dms.csiro.au

    !     NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
    !     BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.

    !     SINGLE PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.

    IMPLICIT NONE

    INTEGER, INTENT(IN)    :: n
    INTEGER                   a(n)
    INTEGER, INTENT(INOUT) :: t(n)

    !     Local Variables

    INTEGER                :: i, j, k, l, r, s, stackl(15), stackr(15), ww
    INTEGER                :: w, x

    s = 1
    stackl(1) = 1
    stackr(1) = n

    !     KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

10  CONTINUE
    l = stackl(s)
    r = stackr(s)
    s = s - 1

    !     KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

20  CONTINUE
    i = l
    j = r
    k = (l+r) / 2
    x = a(k)

    !     REPEAT UNTIL I > J.

    DO
       DO
          IF (a(i).LT.x) THEN                ! Search from lower end
             i = i + 1
             CYCLE
          ELSE
             EXIT
          END IF
       END DO

       DO
          IF (x.LT.a(j)) THEN                ! Search from upper end
             j = j - 1
             CYCLE
          ELSE
             EXIT
          END IF
       END DO

       IF (i.LE.j) THEN                     ! Swap positions i & j
          w = a(i)
          ww = t(i)
          a(i) = a(j)
          t(i) = t(j)
          a(j) = w
          t(j) = ww
          i = i + 1
          j = j - 1
          IF (i.GT.j) EXIT
       ELSE
          EXIT
       END IF
    END DO

    IF (j-l.GE.r-i) THEN
       IF (l.LT.j) THEN
          s = s + 1
          stackl(s) = l
          stackr(s) = j
       END IF
       l = i
    ELSE
       IF (i.LT.r) THEN
          s = s + 1
          stackl(s) = i
          stackr(s) = r
       END IF
       r = j
    END IF

    IF (l.LT.r) GO TO 20
    IF (s.NE.0) GO TO 10

    RETURN
  END SUBROUTINE int_sort
  !
  !  *** CHANGES ***
  !
  !  9/2000, don boyce
  !
  !  (*)  put routines in module structure;
  !
  !  (*)  fixed bug in part_scatter_setup; nodes to be passed which were
  !  .    part of only one element on the local process were ignored; see
  !  .    section of code setting up copy_address_1 and _2;
  !
  !  (*)  put in WAITANY_START to fix MPI bug on velocity;
  !
END MODULE gather_scatter
!
