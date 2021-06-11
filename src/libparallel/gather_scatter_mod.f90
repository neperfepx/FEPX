! This file is part of the FEPX software package.
! Copyright (C) 1996-2021, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE GATHER_SCATTER_MOD
!
! Parallel version of gather/scatter routines. Originally developed by Dave
!   Mika (~1998). Modified by Don Boyce 9/2000---see end of file.
!
! DEB changes (9/2000):
! Put routines in module structure.
! Fixed bug in PART_SCATTER_setup. Nodes to be passed which were part of only
!   one element on the local process were ignored. See section of code setting
!   up COPY_ADDRESS_1 and _2.
! Put in WAITANY_START to fix MPI bug on velocity.
!
! Contains subroutines:
! PART_GATHER: Gather takes global arrays of nodal values to arrays on each
!   element.
! PART_GATHER_START: Starts PART_GATHER
! PART_GATHER_STOP: Stops PART_GATHER
! PART_SCATTER: Accumulate elemental values into a global nodal array.
! PART_SCATTER_SETUP: Prepare data structures for later scatter.
! INT_SORT:
!
! Notes:
! To use these routines, you need first call PART_SCATTER_SETUP to set up the
!   trace structure associated with a connectivity array. THEN you can freely
!   use the PART_SCATTER and PART_GATHER.
! PART_SCATTER has the possibly unexpected side effect of also modifying the
!   elemental array. This is because the elemental array is used as the send
!   buffer for communications and values are accumulated here before sending.
!   Also, the nodal array is not zeroed before use.
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, RK=>REAL_KIND
!
! From libparallel:
!
USE PARALLEL_MOD, ONLY: PAR_QUIT, MYID, MYIDSTR
!
IMPLICIT NONE
!
TYPE TRACE
    !
    ! Buffer : data to be received
    !
    REAL(RK),  ALLOCATABLE, DIMENSION(:) :: BUFFER
    INTEGER, ALLOCATABLE, DIMENSION(:) :: COPIES
    INTEGER, ALLOCATABLE, DIMENSION(:) :: COPY_ADDRESS_1
    INTEGER, ALLOCATABLE, DIMENSION(:) :: COPY_ADDRESS_2
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NODE_ADDRESS
    INTEGER, ALLOCATABLE, DIMENSION(:) :: LOCALS
    INTEGER, ALLOCATABLE, DIMENSION(:) :: N_TO_MOVE
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NODES_TO_MOVE
    INTEGER, ALLOCATABLE, DIMENSION(:) :: MIN_N_TO_MOVE
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NODE_TRACE
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ELM_TRACE
    INTEGER, ALLOCATABLE, DIMENSION(:) :: MIN_PTR
    INTEGER, ALLOCATABLE, DIMENSION(:) :: COPY_PTR
    INTEGER, ALLOCATABLE, DIMENSION(:) :: NODE_PTR
    INTEGER, ALLOCATABLE, DIMENSION(:) :: REQ
    !
    INTEGER :: DIM1_SUB1
    INTEGER :: DIM1_SUP1
    INTEGER :: DIM2_SUB1
    INTEGER :: DIM2_SUP1
    INTEGER :: EL_SUB1
    INTEGER :: EL_SUP1
    INTEGER :: NLOCALS
    INTEGER :: NRECEIVED
    INTEGER :: NELEM_TRACES
    INTEGER :: NNODE_TRACES
    INTEGER :: NTOTAL
    INTEGER :: SETUP
    !
    ! NUMPROCS: Number of processes
    ! MYID: Rank of this process
    ! GS_COMM: Communicator for gather/scatter
    !
    INTEGER :: NUMPROCS
    INTEGER :: MYID
    INTEGER :: GS_COMM
    !
END TYPE TRACE
!
! On velocity, there is a bug in MPI_WAITANY which returns 0-based indices
!   instead of 1-based indices for FORTRAN CALLs.
!
INTEGER, PRIVATE, PARAMETER :: WAITANY_START = 1   ! should be one
!
! Public
!
PUBLIC :: PART_GATHER
PUBLIC :: PART_SCATTER
PUBLIC :: PART_SCATTER_SETUP
!
! Private
!
PRIVATE :: PART_GATHER_START
PRIVATE :: PART_GATHER_STOP
PRIVATE :: INT_SORT
!
CONTAINS
    !
    SUBROUTINE PART_GATHER(ELMS, Y, CONNECTIVITY, TR)
    !
    !  Gather takes global arrays of nodal values to arrays on each
    !  element.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! ELMS: Elemental array
    ! Y: Nodal array
    ! CONNECTIVITY: Connectivity
    ! TR: Data structure interface with MPI
    !
    TYPE(TRACE) :: TR
    REAL(RK) :: ELMS(TR%DIM1_SUB1:TR%DIM1_SUP1, TR%EL_SUB1:TR%EL_SUP1)
    REAL(RK) :: Y(TR%DIM2_SUB1:TR%DIM2_SUP1)
    INTEGER :: CONNECTIVITY(TR%DIM1_SUB1:TR%DIM1_SUP1, TR%EL_SUB1:TR%EL_SUP1)
    !
    !---------------------------------------------------------------------------
    !
    CALL  PART_GATHER_START(ELMS, Y, CONNECTIVITY, TR)
    CALL  PART_GATHER_STOP(ELMS, Y, CONNECTIVITY, TR)
    !
    RETURN
    !  
    END SUBROUTINE PART_GATHER
    !
    !===========================================================================
    !
    SUBROUTINE PART_GATHER_START(ELMS, Y, CONNECTIVITY, TR)
    !
    ! Starts PART_GATHER
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! ELMS:
    ! Y:
    ! CONNECTIVITY:
    ! TR:
    !
    TYPE (TRACE) :: TR
    REAL(RK) :: ELMS(TR%DIM1_SUB1:TR%DIM1_SUP1, TR%EL_SUB1:TR%EL_SUP1)
    REAL(RK) :: Y(TR%DIM2_SUB1:TR%DIM2_SUP1)
    INTEGER :: CONNECTIVITY(TR%DIM1_SUB1:TR%DIM1_SUP1, TR%EL_SUB1:TR%EL_SUP1)
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    INTEGER :: IERR
    !
    !---------------------------------------------------------------------------
    !
    TR%NRECEIVED = 0
    TR%NTOTAL = TR%NELEM_TRACES ! this is number we are going to receive.
    !
    DO I = 1, TR%NUMPROCS
        !
        IF (TR%N_TO_MOVE(I) .NE. 0) THEN
            !
            TR%NRECEIVED = TR%NRECEIVED + 1
            !
            CALL MPI_IRECV(ELMS(TR%DIM1_SUB1,TR%EL_SUB1), 1, TR%ELM_TRACE(I), &
                & I - 1, I - 1, TR%GS_COMM, TR%REQ(TR%NRECEIVED), IERR)
            !
        END IF
        !
        IF (TR%NODES_TO_MOVE(I) .NE. 0) THEN
            !
            TR%NTOTAL = TR%NTOTAL + 1
            CALL MPI_ISEND(Y(TR%DIM2_SUB1), 1, TR%NODE_TRACE(I), I - 1, &
                & TR%MYID, TR%GS_COMM, TR%REQ(TR%NTOTAL), IERR)
            !
        END IF
        !
    END DO
    !
    ! Gather local quantities
    !
    DO J = 1, TR%NLOCALS
        !
        K = TR%LOCALS(J)
        !
        DO I = TR%DIM1_SUB1, TR%DIM1_SUP1
            !
            ELMS(I, K) = Y(CONNECTIVITY(I, K))
            !
        END DO
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE PART_GATHER_START
    !
    !===========================================================================
    !
    SUBROUTINE PART_GATHER_STOP(ELMS, Y, CONNECTIVITY, TR)
    !
    ! Stops PART_GATHER
    !
    !---------------------------------------------------------------------------
    !
    INCLUDE 'mpif-config.h'
    INCLUDE 'mpif-constants.h'
    !
    ! Arguments:
    ! ELMS:
    ! Y:
    ! CONNECTIVITY:
    ! TR
    !
    TYPE (TRACE) :: TR
    REAL(RK) :: ELMS(TR%DIM1_SUB1:TR%DIM1_SUP1, TR%EL_SUB1:TR%EL_SUP1)
    REAL(RK) :: Y(TR%DIM2_SUB1:TR%DIM2_SUP1)
    INTEGER :: CONNECTIVITY(TR%DIM1_SUB1:TR%DIM1_SUP1, TR%EL_SUB1:TR%EL_SUP1)
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: JJ
    INTEGER :: K
    INTEGER :: L
    INTEGER :: M
    INTEGER :: N
    INTEGER :: IERR
    INTEGER :: MPISTATUS(MPI_STATUS_SIZE)
    !
    !---------------------------------------------------------------------------
    !
    ! Pick up the local nodes of the boundary elements
    ! These elements have some off-procesor nodes
    !
    DO J = TR%NLOCALS + 1, TR%EL_SUP1 - TR%EL_SUB1 + 1
        !
        K = TR%LOCALS(J)
        !
        DO I = TR%DIM1_SUB1, TR%DIM1_SUP1
            !
            L = CONNECTIVITY(I, K)
            !
            IF ((L .GE. TR%DIM2_SUB1) .AND. (L .LE. TR%DIM2_SUP1)) THEN
                !
                ELMS(I, K) = Y(L)
                !
            END IF
            !
        END DO
        !
    END DO
    !
    ! Disseminate nodal quantities after checking for completion of sends
    !
    DO I = 1, TR%NTOTAL
        !
        CALL MPI_WAITANY(TR%NTOTAL, TR%REQ, JJ, MPISTATUS, IERR)
        !
        JJ = JJ + 1 - WAITANY_START ! Fixes bug on velocity
        !
        IF (JJ .LE. TR%NRECEIVED ) THEN
            !
            K = MPISTATUS(MPI_SOURCE) + 1 ! Message from proc k is confirmed
            !
            M = 1
            !
            DO J = 1, TR%MIN_N_TO_MOVE(K)
                !
                N = M + TR%COPY_PTR(K)
                !
                DO L = 1, TR%COPIES(J + TR%MIN_PTR(K)) - 1
                    !
                    ELMS(TR%COPY_ADDRESS_1(L + N), TR%COPY_ADDRESS_2(L + N)) &
                        & = ELMS(TR%COPY_ADDRESS_1(N), TR%COPY_ADDRESS_2(N))
                    !
                END DO
                !
                M = M + TR%COPIES(J + TR%MIN_PTR(K))
                !
            END DO
            !
        END IF
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE PART_GATHER_STOP
    !
    !===========================================================================
    !
    SUBROUTINE PART_SCATTER(Y, ELMS, CONNECTIVITY, TR)
    !
    ! Accumulate elemental values into a global nodal array.
    !
    ! NOTE: The elemental array `ELMS' is altered in this routine, and this may
    !   not be expected behavior.
    !
    !---------------------------------------------------------------------------
    !
    INCLUDE 'mpif-handles.h'
    INCLUDE 'mpif-config.h'
    INCLUDE 'mpif-constants.h'
    !
    ! Arguments:
    ! Y: Nodal array
    ! ELMS: Elemental array to be scattered
    ! CONNECTIVITY: Connectivity array
    ! TR: Trace structure
    !
    TYPE(TRACE) :: TR
    REAL(RK) :: Y(TR%DIM2_SUB1:TR%DIM2_SUP1)
    REAL(RK) :: ELMS(TR%DIM1_SUB1:TR%DIM1_SUP1, TR%EL_SUB1:TR%EL_SUP1)
    INTEGER :: CONNECTIVITY(TR%DIM1_SUB1:TR%DIM1_SUP1, TR%EL_SUB1:TR%EL_SUP1)
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: JJ
    INTEGER :: K
    INTEGER :: L
    INTEGER :: M
    INTEGER :: IERR
    INTEGER :: MPISTATUS(MPI_STATUS_SIZE)
    !
    !---------------------------------------------------------------------------
    !
    TR%NRECEIVED = 0
    !
    ! Post accounts receivable.
    !  
    DO I = 1, TR%NUMPROCS
        !
        IF (TR%NODES_TO_MOVE(I) .NE. 0) THEN
            !
            TR%NRECEIVED = TR%NRECEIVED + 1
            !
            IF (RK .EQ. REAL_KIND_D) THEN
                !
                CALL MPI_IRECV(TR%BUFFER(1 + TR%NODE_PTR(I)), &
                    & TR%NODES_TO_MOVE(I), MPI_DOUBLE_PRECISION, I - 1, I - 1, &
                    & TR%GS_COMM, TR%REQ(TR%NRECEIVED), IERR)
                !
            ELSE
                !
                CALL PAR_QUIT("Error  :     > Unsupported MPI data type set &
                    &in `INTRINSIC_TYPES_MOD'.")
                !
            END IF
            !
        END IF
        !
    END DO
    !
    ! Assemble before sending.
    !  
    TR%NTOTAL = TR%NRECEIVED
    DO I = 1, TR%NUMPROCS
        !
        IF (TR%N_TO_MOVE(I) .NE. 0) THEN
            !
            M = 1
            !
            DO JJ = 1, TR%MIN_N_TO_MOVE(I)
                !
                L = M + TR%COPY_PTR(I)
                !
                DO K = 1, TR%COPIES(JJ + TR%MIN_PTR(I)) - 1
                    !
                    ELMS(TR%COPY_ADDRESS_1(L),  TR%COPY_ADDRESS_2(L) ) = &
                        & ELMS(TR%COPY_ADDRESS_1(L), TR%COPY_ADDRESS_2(L) ) + &
                        & ELMS(TR%COPY_ADDRESS_1(K + L), &
                        & TR%COPY_ADDRESS_2(K + L) )
                    !
                END DO
                !
                M = M + TR%COPIES(JJ + TR%MIN_PTR(I))
                !
            END DO
            !
            TR%NTOTAL = TR%NTOTAL + 1
            CALL MPI_ISEND(ELMS(TR%DIM1_SUB1,TR%EL_SUB1), 1, TR%ELM_TRACE(I), &
                & I - 1, TR%MYID, TR%GS_COMM, TR%REQ(TR%NTOTAL), IERR)
            !
        END IF
        !
    END DO
    !
    ! While waiting for sends to complete, scatter local quantities  
    !
    ! The guaranteed local elements:
    !
    DO K = 1, TR%NLOCALS
        !
        J = TR%LOCALS(K)
        !
        DO I = TR%DIM1_SUB1, TR%DIM1_SUP1
            !
            Y(CONNECTIVITY(I, J)) = Y(CONNECTIVITY(I, J)) + ELMS(I, J)
            !
        END DO
        !
    END DO
    !
    ! ------- cut here for PART_SCATTER start (above) and stop (below) ---------
    !
    ! You are guaranteed all the locals of Y are updated
    !
    ! Now scatter the local part of elements are on the boundary:
    !
    DO L = TR%NLOCALS + 1, TR%EL_SUP1 - TR%EL_SUB1 + 1
        !
        J = TR%LOCALS(L)
        !
        DO I = TR%DIM1_SUB1, TR%DIM1_SUP1
            !
            K = CONNECTIVITY(I, J)
            !
            IF ((k .GE. TR%DIM2_SUB1) .AND. (K .LE. TR%DIM2_SUP1)) THEN
                !
                Y(K) = Y(K) + ELMS(I, J)
                !
            END IF
            !
        END DO
        !
    END DO
    !
    ! Check for completion of sends
    !
    DO I = 1, TR%NTOTAL
        !
        CALL MPI_WAITANY(TR%NTOTAL, TR%REQ, J, MPISTATUS, IERR)
        !
        J = J + 1 - WAITANY_START ! Fixes bug on velocity
        !
        IF ( J .LE. TR%NRECEIVED ) THEN
            !
            K = MPISTATUS(MPI_SOURCE) + 1 ! received message from processor k-1
            !
            DO L = 1, TR%NODES_TO_MOVE(K)
                !
                Y(TR%NODE_ADDRESS(L + TR%NODE_PTR(K))) = &
                    & Y(TR%NODE_ADDRESS(L + TR%NODE_PTR(K))) + &
                    & TR%BUFFER(L + TR%NODE_PTR(K))
                !
            END DO
            !
        END IF
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE PART_SCATTER
    !
    !===========================================================================
    !
    SUBROUTINE PART_SCATTER_SETUP(DIM1_SUB1, DIM1_SUP1, DIM2_SUB1, DIM2_SUP1, &
        & EL_SUB1, EL_SUP1, CONNECTIVITY, TR)
    !
    ! Prepare data structures for later scatter.
    !
    !---------------------------------------------------------------------------
    !
    INCLUDE 'mpif-handles.h'
    !
    ! Arguments:
    ! DIM1_SUB1:
    ! DIM1_SUP1:
    ! DIM2_SUB1:
    ! DIM2_SUP1:
    ! EL_SUB1:
    ! EL_SUP1:
    ! CONNECTIVITY:
    ! TR:
    !
    INTEGER, INTENT(IN) :: DIM1_SUB1
    INTEGER, INTENT(IN) :: DIM1_SUP1
    INTEGER, INTENT(IN) :: DIM2_SUB1
    INTEGER, INTENT(IN) :: DIM2_SUP1
    INTEGER, INTENT(IN) :: EL_SUB1
    INTEGER, INTENT(IN) :: EL_SUP1
    INTEGER, INTENT(IN) :: CONNECTIVITY(DIM1_SUB1:DIM1_SUP1, EL_SUB1:EL_SUP1)
    TYPE (TRACE), INTENT(OUT) :: TR
    !
    ! Locals:
    !
    INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFFER
    INTEGER, ALLOCATABLE, DIMENSION(:) :: RECVBUFF
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ADDRESSES
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ELMSIZES
    INTEGER, ALLOCATABLE, DIMENSION(:) :: DIM2_TO_PROCESSORS
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ORDER
    INTEGER, ALLOCATABLE, DIMENSION(:) :: AUX
    INTEGER, ALLOCATABLE, DIMENSION(:) :: MIN_NODE_OFFSETS
    INTEGER, ALLOCATABLE, DIMENSION(:) :: MIN_ELEM_OFFSETS
    INTEGER, ALLOCATABLE, DIMENSION(:) :: DISPLACE
    INTEGER, ALLOCATABLE, DIMENSION(:) :: DIM2_SUB1S
    INTEGER, ALLOCATABLE, DIMENSION(:) :: DIM2_SUP1S
    INTEGER, ALLOCATABLE, DIMENSION(:) :: RECVBUFSIZES
    INTEGER, ALLOCATABLE, DIMENSION(:) :: TEMP
    INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFSIZES
    INTEGER, ALLOCATABLE, DIMENSION(:) :: SDISPLS
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NODE_OFFSETS
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ELEM_OFFSETS_1
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ELEM_OFFSETS_2
    INTEGER :: IPOINT
    INTEGER :: POSITION
    INTEGER :: BUFSIZE
    INTEGER :: MAX_SEND
    INTEGER :: LOCAL_NODES
    INTEGER :: NONLOCAL_NODES
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    INTEGER :: L
    INTEGER :: MAX_NODES
    INTEGER :: MAX_ELEMS
    INTEGER :: NAUX
    INTEGER :: EL_SUB
    INTEGER :: EL_SUP
    INTEGER :: IERR
    INTEGER :: NP2
    INTEGER :: TEMP_SEND(2)
    !
    !----------------------------------------------------------------------
    !
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, TR%NUMPROCS, IERR)
    !
    ! Get a new communicator for the gather/scatters.
    !  
    CALL MPI_CART_CREATE(MPI_COMM_WORLD, 1, (/ TR%NUMPROCS /), (/ .FALSE. /), &
        & .TRUE., TR%GS_COMM, IERR )
    !
    ! Get my position in this communicator.
    !  
    CALL MPI_COMM_RANK(TR%GS_COMM, TR%MYID, IERR)
    !
    ! Allocate storage for structure components and temporary local variables.
    !
    EL_SUB = EL_SUB1 + 1
    EL_SUP = EL_SUP1 + 1
    !
    ALLOCATE(DISPLACE(TR%NUMPROCS))
    ALLOCATE(DIM2_SUB1S(TR%NUMPROCS))
    ALLOCATE(DIM2_SUP1S(TR%NUMPROCS))
    ALLOCATE(RECVBUFSIZES(TR%NUMPROCS))
    ALLOCATE(TEMP(2 * TR%NUMPROCS))
    ALLOCATE(BUFSIZES(TR%NUMPROCS))
    ALLOCATE(SDISPLS(TR%NUMPROCS))
    !
    ALLOCATE(TR%N_TO_MOVE(TR%NUMPROCS), STAT = IERR)
    ALLOCATE(TR%NODES_TO_MOVE(TR%NUMPROCS), STAT = IERR)
    ALLOCATE(TR%MIN_N_TO_MOVE(TR%NUMPROCS), STAT = IERR)
    ALLOCATE(TR%LOCALS(EL_SUP-EL_SUB1), STAT = IERR)
    ALLOCATE(TR%NODE_TRACE(TR%NUMPROCS), STAT = IERR)
    ALLOCATE(TR%ELM_TRACE(TR%NUMPROCS), STAT = IERR)
    !
    IF (IERR /= 0) THEN
        !
        CALL PAR_QUIT('allocation failure')
        !
    END IF
    !
    ALLOCATE(TR%MIN_PTR(TR%NUMPROCS), STAT = IERR)
    ALLOCATE(TR%COPY_PTR(TR%NUMPROCS), STAT = IERR)
    ALLOCATE(TR%NODE_PTR(TR%NUMPROCS), STAT = IERR)
    !
    NP2 = 2 * (TR%NUMPROCS - 1)
    !
    ALLOCATE(TR%REQ(NP2), STAT=IERR)
    !
    TR%DIM1_SUB1 = DIM1_SUB1
    TR%DIM1_SUP1 = DIM1_SUP1
    TR%DIM2_SUB1 = DIM2_SUB1
    TR%DIM2_SUP1 = DIM2_SUP1
    TR%EL_SUB1 = EL_SUB1
    TR%EL_SUP1 = EL_SUP1
    !
    ! This section appears very limiting. Each node allocates storage for all of
    !   the arrays on all of the nodes.
    !
    I = (EL_SUP - EL_SUB1) * (TR%DIM1_SUP1 - TR%DIM1_SUB1 + 1)
    ALLOCATE(NODE_OFFSETS (I, TR%NUMPROCS))
    ALLOCATE(ELEM_OFFSETS_1(I, TR%NUMPROCS))
    ALLOCATE(ELEM_OFFSETS_2(I, TR%NUMPROCS))
    !
    ! Get the decomposition boundaries from all processes.
    !
    TEMP_SEND(1) = TR%DIM2_SUB1
    TEMP_SEND(2) = TR%DIM2_SUP1
    CALL MPI_ALLGATHER( TEMP_SEND, 2, MPI_INTEGER, TEMP, 2, MPI_INTEGER,&
        & TR%GS_COMM, IERR)
    !
    J = 1
    DO I = 1, TR%NUMPROCS
        !
        DIM2_SUB1S(I) = TEMP(J)
        DIM2_SUP1S(I) = TEMP(J + 1)
        J = J + 2
        !
    END DO
    !
    ! Find the traces for gather/scatter.
    !
    ! Again, all nodes in the entire problem are allocated here, and then set to
    !   a constant value.
    !  
    ALLOCATE (DIM2_TO_PROCESSORS(DIM2_SUB1S(1):DIM2_SUP1S(TR%NUMPROCS)))
    !
    DO I = 1, TR%NUMPROCS
        !
        DIM2_TO_PROCESSORS(DIM2_SUB1S(I):DIM2_SUP1S(I)) = I
        !
    END DO
    !
    TR%N_TO_MOVE = 0
    TR%NLOCALS = 0
    NONLOCAL_NODES = 0
    !
    DO I = EL_SUB1, EL_SUP1
        !
        LOCAL_NODES = 0
        !
        DO J = TR%DIM1_SUB1, TR%DIM1_SUP1
            !
            IPOINT = CONNECTIVITY(J, I)
            !
            IF (IPOINT .LT. 0) CALL PAR_QUIT('Error  :     > &
              &IPOINT out of bounds')
            !
            IF (IPOINT > DIM2_SUP1S(TR%NUMPROCS) ) THEN
                !
                CALL PAR_QUIT('segv', .TRUE.)
                !
            END IF
            !
            IF (DIM2_TO_PROCESSORS(IPOINT) - 1 .NE. TR%MYID) THEN ! moving
                !
                K = DIM2_TO_PROCESSORS(IPOINT)
                TR%N_TO_MOVE(K) = TR%N_TO_MOVE(K) + 1
                ELEM_OFFSETS_1(TR%N_TO_MOVE(K), K) = J
                ELEM_OFFSETS_2(TR%N_TO_MOVE(K), K) = I
                NODE_OFFSETS(TR%N_TO_MOVE(K), K) = IPOINT
                !
            ELSE ! node is local
                !
                LOCAL_NODES = LOCAL_NODES + 1
                !
            END IF
            !
        END DO
        !
        ! Check if all nodes were local.  In that case, the element is local
        !
        IF (LOCAL_NODES .EQ. (TR%DIM1_SUP1 - TR%DIM1_SUB1 + 1)) THEN
            !
            ! Element is local
            !
            TR%NLOCALS = TR%NLOCALS + 1
            TR%LOCALS(TR%NLOCALS) = I
            !
        ELSE
            !
            ! Part of element is on another processor.
            ! Save these too, but starting from the top and counting down
            !
            TR%LOCALS(EL_SUP - EL_SUB1 - NONLOCAL_NODES) = I
            NONLOCAL_NODES = NONLOCAL_NODES + 1
            !
        END IF
        !
    END DO
    !
    DEALLOCATE(DIM2_TO_PROCESSORS)
    !
    ! Each processor has tagged the nodes that it does not have
    !
    ! To save space, use a new storage pattern from here out  (i.e. N_TO_MOVE
    ! could be large for one processor and zero for the next)
    !
    I = SUM(TR%N_TO_MOVE)
    !
    ALLOCATE(MIN_ELEM_OFFSETS(I))
    ALLOCATE(MIN_NODE_OFFSETS(I))
    ALLOCATE(TR%COPIES(I))
    ALLOCATE(TR%COPY_ADDRESS_1(I))
    ALLOCATE(TR%COPY_ADDRESS_2(I))
    !
    NAUX = MAXVAL(TR%N_TO_MOVE)
    ALLOCATE (ORDER(NAUX))
    ALLOCATE (AUX(NAUX))
    !
    ! Find the minimal set (we do not want to send duplicates)
    ! First sort the nodes to make finding the minimal set more efficient
    !
    TR%MIN_N_TO_MOVE = 0
    TR%COPIES = 0
    K = 0
    L = 0
    !
    DO I = 1, TR%NUMPROCS
        !
        TR%MIN_PTR(I) = K
        TR%COPY_PTR(I) = L
        IF (TR%N_TO_MOVE(I) .NE. 0) THEN
            !
            DO J = 1, TR%N_TO_MOVE(I)
                !
                ORDER(J) = J
                !
            END DO
            !
            ! This is not a "stable" sort
            !
            CALL INT_SORT(NODE_OFFSETS(:, I), TR%N_TO_MOVE(I), ORDER)
            !
            TR%MIN_N_TO_MOVE(I) = TR%MIN_N_TO_MOVE(I) + 1
            MIN_NODE_OFFSETS(1 + TR%MIN_PTR(I)) = NODE_OFFSETS(1, I)
            !
            ! Switching MIN_ELEM_OFFSETS to an absolute offset from the
            !   beginning of the buffer instead of a 1 and 2 index into the
            !   array. The routine MPI_TYPE_INDEXED below will need it this way.
            !
            MIN_ELEM_OFFSETS(1 + TR%MIN_PTR(I)) = (TR%DIM1_SUP1 - TR%DIM1_SUB1 &
                & + 1) * (ELEM_OFFSETS_2(ORDER(1), I) - EL_SUB1) + &
                & ELEM_OFFSETS_1(ORDER(1), I)
            !
            J = 1
            L = L + 1
            TR%COPY_ADDRESS_1(L) = ELEM_OFFSETS_1(ORDER(J), I)
            TR%COPY_ADDRESS_2(L) = ELEM_OFFSETS_2(ORDER(J), I)
            !
            TR%COPIES(TR%MIN_N_TO_MOVE(I) + TR%MIN_PTR(I)) = &
            TR%COPIES(TR%MIN_N_TO_MOVE(I) + TR%MIN_PTR(I)) + 1
            !
            DO J = 2, TR%N_TO_MOVE(I)
                !
                IF (NODE_OFFSETS(J, I) .EQ. NODE_OFFSETS(J - 1, I)) THEN
                    !
                    ! The connectivities are the same.  In a scatter, the
                    ! unassembled elems with the same global nodes will need to
                    ! be added before sending, and in a gather, these will need
                    ! to be assigned the same value after that value is received
                    ! from another processor.
                    !
                ELSE
                    !
                    ! This is the first copy. The value in this address will be
                    !   sent or assigned in the communication.
                    !
                    TR%MIN_N_TO_MOVE(I) = TR%MIN_N_TO_MOVE(I) + 1
                    MIN_NODE_OFFSETS(TR%MIN_N_TO_MOVE(I) + TR%MIN_PTR(I)) = &
                        & NODE_OFFSETS(J, I)
                    MIN_ELEM_OFFSETS(TR%MIN_N_TO_MOVE(I) + TR%MIN_PTR(I)) = &
                        & (TR%DIM1_SUP1 - TR%DIM1_SUB1 + 1) * &
                        & (ELEM_OFFSETS_2(ORDER(J), I) - EL_SUB1) + &
                        & ELEM_OFFSETS_1(ORDER(J), I)
                    !
                END IF
                !
                TR%COPIES(TR%MIN_N_TO_MOVE(I)  +TR%MIN_PTR(I)) = &
                    & TR%COPIES(TR%MIN_N_TO_MOVE(I) + TR%MIN_PTR(I)) + 1
                !
                L = L + 1
                TR%COPY_ADDRESS_1(l) = ELEM_OFFSETS_1(ORDER(J), I)
                TR%COPY_ADDRESS_2(l) = ELEM_OFFSETS_2(ORDER(J), I)
                !
            END DO
            !
        END IF
        !
        K = K + TR%MIN_N_TO_MOVE(I)
        !
    END DO
    !
    DEALLOCATE (ORDER, AUX, ELEM_OFFSETS_1, ELEM_OFFSETS_2)
    !
    ! Distribute the traces we just calculated but first determine size of
    !   needed buffer
    !
    MAX_SEND = SUM(TR%MIN_N_TO_MOVE) + TR%NUMPROCS
    !
    CALL MPI_PACK_SIZE(MAX_SEND, MPI_INTEGER, TR%GS_COMM, BUFSIZE, IERR)
    !
    ALLOCATE(BUFFER(BUFSIZE/4)) ! BUFSIZE is in bytes
    !
    ! Pack data into buffer
    !
    POSITION = 0  ! Start at beginning of buffer
    !
    DO I = 1, TR%NUMPROCS
        !
        SDISPLS(I) = POSITION
        J = TR%MIN_N_TO_MOVE(I)
        !
        CALL MPI_PACK(J, 1, MPI_INTEGER, BUFFER, BUFSIZE, POSITION, &
            & TR%GS_COMM, IERR)
        !
        IF (J .NE. 0) THEN
            !
            CALL MPI_PACK(MIN_NODE_OFFSETS(1 + TR%MIN_PTR(I)), J, MPI_INTEGER, &
                & BUFFER, BUFSIZE, POSITION, TR%GS_COMM, IERR)
            !
        END IF
        !
        BUFSIZES(I) = POSITION - SDISPLS(I)
        !
    END DO
    !
    DEALLOCATE (MIN_NODE_OFFSETS)
    !
    ! Distribute the buffersizes in preperation for data transfer
    !
    CALL MPI_ALLTOALL(BUFSIZES, 1, MPI_INTEGER, RECVBUFSIZES, 1, MPI_INTEGER, &
        & TR%GS_COMM, IERR)
    !
    ! Allocate space for incomming data
    !
    BUFSIZE = SUM(RECVBUFSIZES)   ! BUFSIZE is in bytes
    ALLOCATE(RECVBUFF(BUFSIZE/4))
    !
    J = 0
    !
    DO I = 1, TR%NUMPROCS
        !
        DISPLACE(I) = J
        J = J + RECVBUFSIZES(I)
        !
    END DO
    !
    ! distribute! (finally)
    !
    CALL MPI_ALLTOALLV(BUFFER, BUFSIZES, SDISPLS, MPI_PACKED, RECVBUFF, &
        & RECVBUFSIZES, DISPLACE, MPI_PACKED, TR%GS_COMM, IERR)
    !
    DEALLOCATE(BUFFER)
    !
    ! unpack data...
    ! get sizes first
    !
    POSITION = 0
    !
    DO I = 1, TR%NUMPROCS
        !
        CALL MPI_UNPACK(RECVBUFF, BUFSIZE, POSITION, TR%NODES_TO_MOVE(I), &
            & 1, MPI_INTEGER, TR%GS_COMM, IERR)
        !
        POSITION = POSITION + TR%NODES_TO_MOVE(I) * 4 ! POSITION in bytes
        !
    END DO
    !
    MAX_NODES = MAXVAL(TR%NODES_TO_MOVE)
    MAX_ELEMS = MAXVAL(TR%MIN_N_TO_MOVE)
    DEALLOCATE(NODE_OFFSETS)
    ALLOCATE(NODE_OFFSETS(MAX_NODES, TR%NUMPROCS))
    ALLOCATE(ELMSIZES(Max(MAX_NODES, MAX_ELEMS)))
    ALLOCATE(ADDRESSES(MAX_NODES))
    !
    POSITION = 0
    !
    DO I = 1, TR%NUMPROCS
        !
        POSITION = POSITION + 4  ! POSITION is in bytes
        J = TR%NODES_TO_MOVE(I)
        !
        IF (J .NE. 0) THEN
            !
            CALL MPI_UNPACK(RECVBUFF, BUFSIZE, POSITION, NODE_OFFSETS(1, I), &
                & J, MPI_INTEGER, TR%GS_COMM, IERR)
            !
        END IF
        !
    END DO
    !
    DEALLOCATE(RECVBUFF)
    !
    ! Commit the traces
    !
    ELMSIZES = 1
    !
    TR%NELEM_TRACES = 0
    TR%NNODE_TRACES = 0
    !
    DO I = 1, TR%NUMPROCS
        !
        IF (TR%NODES_TO_MOVE(I) .NE. 0) THEN
            !
            TR%NNODE_TRACES = TR%NNODE_TRACES + 1
            !
            DO J = 1, TR%NODES_TO_MOVE(I)
                !
                ADDRESSES(J) = NODE_OFFSETS(J, I) - TR%DIM2_SUB1
                !
            END DO
            !
            IF (RK .EQ. REAL_KIND_D) THEN
                !
                CALL MPI_TYPE_INDEXED(TR%NODES_TO_MOVE(I), ELMSIZES, &
                    & ADDRESSES, MPI_DOUBLE_PRECISION, TR%NODE_TRACE(I), IERR)
                !
            ELSE
                !
                CALL PAR_QUIT("Error  :     > Unsupported MPI data type set &
                    &in `INTRINSIC_TYPES_MOD'.")
                !
            END IF
            !
            CALL MPI_TYPE_COMMIT(TR%NODE_TRACE(I), IERR)
            !
        END IF
        !
        IF (TR%MIN_N_TO_MOVE(I) .NE. 0) THEN
            !
            TR%NELEM_TRACES = TR%NELEM_TRACES + 1
            !
            IF (RK .EQ. REAL_KIND_D) THEN
                !
                CALL MPI_TYPE_INDEXED(TR%MIN_N_TO_MOVE(I), ELMSIZES, &
                    & MIN_ELEM_OFFSETS(1 + TR%MIN_PTR(I)), &
                    & MPI_DOUBLE_PRECISION, TR%ELM_TRACE(I), IERR )
                !
            ELSE
                !
                CALL PAR_QUIT("Error  :     > Unsupported MPI data type set &
                    &in `INTRINSIC_TYPES_MOD'.")
                !
            END IF
            !
            CALL MPI_TYPE_COMMIT(TR%ELM_TRACE(I), IERR)
            !
        END IF
        !
    END DO
    !
    ! Since 1 does not provide a RECV with an ADD, it is inefficient to use
    !   MPI_TYPE_INDEXED with the PART_SCATTER.  Instead, use a buffer to hold
    !   the incomming message and add message to appropriate locations.
    !   Hopefully the needed function will be added to the 1 standard - then the
    !   above NODE_TRACE can be used.
    !
    BUFSIZE = SUM(TR%NODES_TO_MOVE)
    ALLOCATE(TR%NODE_ADDRESS(BUFSIZE))
    ALLOCATE(TR%BUFFER(BUFSIZE))
    !
    J = 0
    !
    DO I = 1, TR%NUMPROCS
        !
        TR%NODE_PTR(I) = J
        J = J + TR%NODES_TO_MOVE(I)
        !
    END DO
    !
    DO I = 1, TR%NUMPROCS
        !
        DO J = 1, TR%NODES_TO_MOVE(I)
            !
            TR%NODE_ADDRESS(J + TR%NODE_PTR(I)) = NODE_OFFSETS(J, I)
            !
        END DO
        !
    END DO
    !
    DEALLOCATE(ELMSIZES, ADDRESSES, NODE_OFFSETS, MIN_ELEM_OFFSETS)
    !
    DEALLOCATE(DISPLACE, DIM2_SUB1S, DIM2_SUP1S, RECVBUFSIZES, TEMP, BUFSIZES, &
        & SDISPLS)
    !
    TR%SETUP = 0
    !
    RETURN
    !
    END SUBROUTINE PART_SCATTER_SETUP
    !
    !===========================================================================
    !
    SUBROUTINE INT_SORT(A, N, T)
    !
    ! Version 1.00, 14 July 1995
    !Author: Alan Miller
    !        CSIRO Division of Mathematics & Statistics
    !        Private Bag 10, Rosebank MDC
    !        Clayton 3169, Victoria, Australia
    !        Phone: (+61) 3 9545-8036      Fax: (+61) 3 9545-8080
    !        e-mail: Alan.Miller @ mel.dms.csiro.au
    !
    ! NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
    ! BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
    !
    ! SINGLE PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! A:
    ! N:
    ! T:
    !
    INTEGER :: A(N)
    INTEGER, INTENT(IN) :: N
    INTEGER, INTENT(INOUT) :: T(N)
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: J
    INTEGER :: K
    INTEGER :: L
    INTEGER :: R
    INTEGER :: S
    INTEGER :: STACKL(15)
    INTEGER :: STACKR(15)
    INTEGER :: WW
    INTEGER :: W
    INTEGER :: X
    !
    !---------------------------------------------------------------------------
    !
    S = 1
    STACKL(1) = 1
    STACKR(1) = N
    !
    ! KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.
    !
    10  CONTINUE
    !
    L = STACKL(S)
    R = STACKR(S)
    S = S - 1
    !
    ! KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.
    !
    20  CONTINUE
    !
    I = L
    J = R
    K = (L + R) / 2
    X = A(K)
    !
    !     REPEAT UNTIL I > J.
    !
    DO
        !
        DO
            !
            IF (A(I) .LT. X) THEN ! Search from lower end
                !
                I = I + 1
                !
                CYCLE
                !
            ELSE
                !
                EXIT
                !
            END IF
            !
        END DO
        !
        DO
            !
            IF (X .LT. A(J)) THEN ! Search from upper end
                !
                J = J - 1
                !
                CYCLE
                !
            ELSE
                !
                EXIT
                !
            END IF
            !
        END DO
        !
        IF (I .LE. J) THEN ! Swap positions i & j
            !
            W = A(I)
            WW = T(I)
            A(I) = A(J)
            T(I) = T(J)
            A(J) = W
            T(J) = WW
            I = I + 1
            J = J - 1
            !
            IF (I .GT. J) EXIT
            !
        ELSE
            !
            EXIT
            !
        END IF
        !
    END DO
    !
    IF (J - L .GE. R - I) THEN
        !
        IF (L .LT. J) THEN
            !
            S = S + 1
            STACKL(S) = L
            STACKR(S) = J
            !
        END IF
        !
        L = I
        !
    ELSE
        !
        IF (I .LT. R) THEN
            !
            S = S + 1
            STACKL(S) = I
            STACKR(S) = R
            !
        END IF
        !
        R = J
        !
    END IF
    !
    IF (L .LT. R) GO TO 20
    IF (S .NE. 0) GO TO 10
    !
    RETURN
    !
    END SUBROUTINE INT_SORT
    !
END MODULE GATHER_SCATTER_MOD
