! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE parallel_mod
  !
  !  Module for parallel variables and operations. (MPI version)
  !
  !  This module provides information and elementary
  !  operations for the primary communicator in an MPI-based
  !  parallel program.
  !
  IMPLICIT NONE
  INCLUDE 'mpif.h'
  !
  ! ==================== Public Entities
  !
  PRIVATE
  !
  PUBLIC :: numprocs, myid, myidstr
  PUBLIC :: par_init, par_partition, par_sum, par_max, par_min, par_message, par_quit, par_gather
  !
  ! ==================== Communicator Info
  ! 
  !  numprocs : number of processes
  !  myid     : number of the current process
  !  myidstr  : id string for current process
  !
  INTEGER :: numprocs, myid
  CHARACTER(LEN=8) :: myidstr
  !
CONTAINS
  !
  !  PUBLIC PROCEDURES:
  !  -----------------
  !
  !  par_init      : initialize for parallel execution.
  !  par_partition : partition arrays across processors.
  !  par_sum       : sum scalars across partitions.
  !  par_message   : prints a message from one or all processes
  !  par_quit      : clean up and exit.
  !
  !  MRE_DECOMP1D  : decompose an array across processors.
  !
  SUBROUTINE par_init()
    !
    !  Performs MPI initializations.  Sets global variables
    !  `numprocs' and `myid'.  This routine has no arguments.
    !
    IMPLICIT NONE
    !
    !  Arguments.
    !
    !  *** NONE ***
    !
    !  Locals.
    !
    INTEGER ierr
    !
    !----------------------------------------------------------------------
    !
    call MPI_INIT( ierr )

    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
    !
    !  Get my position in this communicator.
    !
    call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
    !
    write(myidstr, '(i5)') myid
    myidstr = ADJUSTL(myidstr)
    !
    RETURN
  END SUBROUTINE par_init
  !
  !**********************************************************************
  !
  SUBROUTINE par_partition(na, np, me, dl, du)
    !
    !  Partition arrays for multiple processes.  Ranges go
    !  from 0 to (na - 1).
    !
    IMPLICIT NONE
    !
    !  Arguments
    !
    INTEGER, INTENT(IN)  :: na, np, me
    INTEGER, INTENT(OUT) :: dl, du
    !
    !  na : size of array
    !  np : number of processes
    !  me : my process number
    !  dl : lower dimension of my array segment
    !  du : upper dimension of my array segment
    !
    !----------------------------------------------------------------------
    !
    call MRE_DECOMP1D(na, np, me, dl, du)
    dl = dl - 1
    du = du - 1
    !
    RETURN
  END SUBROUTINE par_partition
  !
  !**********************************************************************
  !
  SUBROUTINE par_sum(part, whole)
    !
    !  Sum scalars across partitions.
    !
    IMPLICIT NONE
    !
    !  Arguments.
    !
    !  part  : this processor's part
    !  whole : the sum of all parts
    !
    DOUBLE PRECISION part, whole
    !
    !  Locals.
    !
    INTEGER ierr
    !
    !----------------------------------------------------------------------
    !
    call MPI_AllReduce(part, whole, 1, MPI_DOUBLE_PRECISION, &
         &    MPI_SUM, MPI_COMM_WORLD, ierr)
    !
    RETURN
  END SUBROUTINE par_sum
  !
  !**********************************************************************
  !
  SUBROUTINE par_max(part, maxpart)
    !
    !  Find maximum of parts.
    !
    IMPLICIT NONE
    !
    !  Arguments.
    !
    !  part    : this processor's part
    !  maxpart : the maximum of all parts
    !
    DOUBLE PRECISION part, maxpart
    !
    !  Locals.
    !
    INTEGER ierr
    !
    !----------------------------------------------------------------------
    !
    call MPI_AllReduce(part, maxpart, 1, MPI_DOUBLE_PRECISION, &
         &    MPI_MAX, MPI_COMM_WORLD, ierr)
    !
    RETURN
  END SUBROUTINE par_max
  !
  !**********************************************************************
  !
  SUBROUTINE par_min(part, minpart)
    !
    !  Find minimum of parts.
    !
    IMPLICIT NONE
    !
    !  Arguments.
    !
    !  part    : this processor's part
    !  minpart : the minimum of all parts
    !
    DOUBLE PRECISION part, minpart
    !
    !  Locals.
    !
    INTEGER ierr
    !
    !----------------------------------------------------------------------
    !
    call MPI_AllReduce(part, minpart, 1, MPI_DOUBLE_PRECISION, &
         &    MPI_MIN, MPI_COMM_WORLD, ierr)
    !
    RETURN
  END SUBROUTINE par_min
  !
  !**********************************************************************
  !
  SUBROUTINE par_gather(part, full, n)
    !
    !  Gather array data across partitions.
    !
    IMPLICIT NONE
    !
    !  Arguments.
    !
    !  part  : this processor's part
    !  full  : the array containing all parts
    !  n     : the amount of data being sent and recieved (length of array)
    !
    INTEGER :: n
    INTEGER, DIMENSION(:) :: part
    INTEGER, DIMENSION(:,:) :: full
    !
    !  Locals.
    !
    INTEGER ierr
    !
    !----------------------------------------------------------------------
    !
    call MPI_Gather(part, n, MPI_INTEGER, full, n, MPI_INTEGER, 0, &
        & MPI_COMM_WORLD, ierr)
    !
    RETURN
  END SUBROUTINE par_gather
  !
  !**********************************************************************
  !
  SUBROUTINE par_quit(message, ABORT)
    !
    !  Clean up and exit.
    !
    IMPLICIT NONE
    !
    !  Arguments.
    !
    CHARACTER*(*), INTENT(IN) :: message
    LOGICAL, INTENT(IN), OPTIONAL :: ABORT
    !
    !  Locals.
    !
    INTEGER :: ierr, enverr = 1
    LOGICAL :: use_abort
    CHARACTER(LEN=72), PARAMETER :: footer = "============================&
        &============================================"
    !
    !----------------------------------------------------------------------
    !
    use_abort = .FALSE.
    if (PRESENT(ABORT)) then
       if (ABORT) then
          use_abort = .TRUE.
       end if
    end if
    !
    call par_message(6, message, ALLWRITE=.FALSE.)
    call par_message(6, footer, ALLWRITE=.FALSE.)
    !call par_message(6, 'Exiting.', ALLWRITE=.FALSE.)
    !
    if (use_abort) then
       call MPI_Abort(MPI_COMM_WORLD, enverr, ierr)
    else
       call MPI_Finalize(ierr)
    end if
    !
    STOP
    !
  END SUBROUTINE par_quit
  !
  !**********************************************************************
  !
  SUBROUTINE par_message(unit, message, ALLWRITE)
    !
    !  Write message, either from master process or from all.
    !
    !  Arguments:
    !
    !  unit -- unit number to write to
    !  message -- message to write
    !  ALLWRITE -- flags whether all processes write or just master 
    !  .        -- default=FALSE [only master node writes message]
    !
    INTEGER, INTENT(IN) :: unit
    !
    CHARACTER(LEN=*), INTENT(IN) :: message
    !
    LOGICAL, INTENT(IN), OPTIONAL :: ALLWRITE
    !
    !  Locals:
    !
    LOGICAL :: aw
    !
    !--------------*-------------------------------------------------------
    !
    if (PRESENT(ALLWRITE)) then
       aw = ALLWRITE
    else
       aw = .FALSE.
    end if
    !
    if (aw) then
       write(unit, '(a)') 'Info   : [Rank '//TRIM(myidstr(1:5))//'] - '//message
    else
       if (myid == 0) then
          !write(unit, '(a)') 'Info   : [Rank '//TRIM(myidstr(1:5))//'] - '//message
          write(unit, '(a)') message
       end if
    end if
    !
  END SUBROUTINE par_message
  !
  !**********************************************************************
  !
  SUBROUTINE MRE_DECOMP1D( n, numprocs, myid, s, e )
    !
    !  This file contains a routine for producing a decomposition of a 1-d array
    !  when given a number of processors.  It may be used in "direct" product
    !  decomposition.  The values returned assume a "global" domain in [1:n]
    !
    !  Arguments.
    !
    !  n        : size of array
    !  numprocs : number of processes
    !  myid     : my process number (0 to numprocs-1 )
    !  s, e     : start and end locations for this process
    !
    INTEGER n, numprocs, myid, s, e
    !
    !  Locals.
    !
    INTEGER nlocal, deficit
    !
    !----------------------------------------------------------------------
    !
    nlocal  = n / numprocs
    s       = myid * nlocal + 1
    deficit = mod(n,numprocs)
    s       = s + min(myid,deficit)
    if (myid .lt. deficit) then
       nlocal = nlocal + 1
    endif
    e = s + nlocal - 1
    if (e .gt. n .or. myid .eq. numprocs-1) e = n

    RETURN
  END SUBROUTINE MRE_DECOMP1D

END MODULE parallel_mod
