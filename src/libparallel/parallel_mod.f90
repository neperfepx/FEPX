! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE PARALLEL_MOD
!
! Module for parallel variables and operations (MPI version). This module
!   provides information and elementary operations for the primary communicator
!   in an MPI-based parallel program.
!
! Contains subroutines:
! PAR_INIT: Performs MPI initializations. Sets variables `NUMPROCS' and `MYID'
! PAR_PARTITION: Partition arrays for multiple processes
! PAR_SUM: Sum scalars across partitions
! PAR_MAX: Find maximum of parts
! PAR_MIN: Find minimum of parts
! PAR_GATHER: Gather array data across partitions
! PAR_QUIT: Clean up and exit
! PAR_MESSAGE: Write message, either from master process or from all
! MRE_DECOMP1D:  This file contains a routine for producing a decomposition of a
!   1-d array when given a number of processors.  It may be used in "direct"
!   product decomposition. The values returned assume a "global" domain in [1:n]
!
IMPLICIT NONE
!
INCLUDE 'mpif.h'
!
! Private
!
PRIVATE
!
! Public
!
PUBLIC :: NUMPROCS
PUBLIC :: MYID
PUBLIC :: MYIDSTR
PUBLIC :: PAR_INIT
PUBLIC :: PAR_PARTITION
PUBLIC :: PAR_SUM
PUBLIC :: PAR_MAX
PUBLIC :: PAR_MIN
PUBLIC :: PAR_MESSAGE
PUBLIC :: PAR_QUIT
PUBLIC :: PAR_GATHER
!
! NUMPROCS: Number of processes
! MYID: Number of the current process
! MYIDSTR: ID string for current process
!
INTEGER :: NUMPROCS
INTEGER :: MYID
CHARACTER(LEN=8) :: MYIDSTR
  !
CONTAINS
    !
    SUBROUTINE PAR_INIT()
    !
    !  Performs MPI initializations.  Sets global variables `NUMPROCS' and
    !   `MYID'.  This routine has no arguments.
    !
    !---------------------------------------------------------------------------
    !
    ! Locals:
    !
    INTEGER :: IERR
    !
    !---------------------------------------------------------------------------
    !
    CALL MPI_INIT(IERR)

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NUMPROCS, IERR)
    !
    ! Get my position in this communicator.
    !
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)
    !
    WRITE(MYIDSTR, '(i5)') MYID
    MYIDSTR = ADJUSTL(MYIDSTR)
    !
    RETURN
    !
    END SUBROUTINE PAR_INIT
    !
    !===========================================================================
    !
    SUBROUTINE PAR_PARTITION(NA, NP, ME, DL, DU)
    !
    ! Partition arrays for multiple processes. Ranges go from 0 to (NA - 1).
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! NA: Size of array
    ! NP: Number of processes
    ! ME: My process number
    ! DL: Lower dimension of my array segment
    ! DU: Upper dimension of my array segment
    !
    INTEGER, INTENT(IN) :: NA
    INTEGER, INTENT(IN) :: NP
    INTEGER, INTENT(IN) :: ME
    INTEGER, INTENT(OUT) :: DL
    INTEGER, INTENT(OUT) :: DU
    !
    !---------------------------------------------------------------------------
    !
    CALL MRE_DECOMP1D(NA, NP, ME, DL, DU)
    !
    DL = DL - 1
    DU = DU - 1
    !
    RETURN
    !
    END SUBROUTINE PAR_PARTITION
    !
    !===========================================================================
    !
    SUBROUTINE PAR_SUM(PART, WHOLE)
    !
    ! Sum scalars across partitions.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! PART: This processor's part
    ! WHOLE: The sum of all parts
    !
    DOUBLE PRECISION :: PART
    DOUBLE PRECISION :: WHOLE
    !
    ! Locals:
    !
    INTEGER :: IERR
    !
    !---------------------------------------------------------------------------
    !
    CALL MPI_ALLREDUCE(PART, WHOLE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
        & MPI_COMM_WORLD, IERR)
    !
    RETURN
    !
    END SUBROUTINE PAR_SUM
    !
    !===========================================================================
    !
    SUBROUTINE PAR_MAX(PART, MAXPART)
    !
    ! Find maximum of parts.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! PART: This processor's part
    ! MAXPART: The maximum of all parts
    !
    DOUBLE PRECISION :: PART
    DOUBLE PRECISION :: MAXPART
    !
    ! Locals:
    !
    INTEGER :: IERR
    !
    !---------------------------------------------------------------------------
    !
    CALL MPI_ALLREDUCE(PART, MAXPART, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
        & MPI_COMM_WORLD, IERR)
    !
    RETURN
    !
    END SUBROUTINE PAR_MAX
    !
    !===========================================================================
    !
    SUBROUTINE PAR_MIN(PART, MINPART)
    !
    ! Find minimum of parts.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! PART: This processor's part
    ! MINPART: The minimum of all parts
    !
    DOUBLE PRECISION :: PART
    DOUBLE PRECISION :: MINPART
    !
    ! Locals:
    !
    INTEGER :: IERR
    !
    !---------------------------------------------------------------------------
    !
    CALL MPI_ALLREDUCE(PART, MINPART, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
        & MPI_COMM_WORLD, IERR)
    !
    RETURN
    !
    END SUBROUTINE PAR_MIN
    !
    !===========================================================================
    !
    SUBROUTINE PAR_GATHER(PART, FULL, N)
    !
    !  Gather array data across partitions.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! PART: This processor's part
    ! FULL: The array containing all parts
    ! N: The amount of data being sent and recieved (length of array)
    !
    INTEGER, DIMENSION(:) :: PART
    INTEGER, DIMENSION(:, :) :: FULL
    INTEGER :: N
    !
    ! Locals:
    !
    INTEGER :: IERR
    !
    !---------------------------------------------------------------------------
    !
    CALL MPI_GATHER(PART, N, MPI_INTEGER, FULL, N, MPI_INTEGER, 0, &
        & MPI_COMM_WORLD, IERR)
    !
    RETURN
    !
    END SUBROUTINE PAR_GATHER
    !
    !===========================================================================
    !
    SUBROUTINE PAR_QUIT(MESSAGE, ABORT)
    !
    !  Clean up and exit.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! MESSAGE: Message to write
    ! ABORT: Flag that forces MPI to abort over finalize
    !
    CHARACTER*(*), INTENT(IN) :: MESSAGE
    LOGICAL, INTENT(IN), OPTIONAL :: ABORT
    !
    ! Locals:
    !
    INTEGER :: IERR, ENVERR = 1
    LOGICAL :: USE_ABORT
    CHARACTER(LEN=72), PARAMETER :: FOOTER = "============================&
        &============================================"
    !
    !---------------------------------------------------------------------------
    !
    USE_ABORT = .FALSE.
    !
    IF (PRESENT(ABORT)) THEN
        !
        IF (ABORT) THEN
            !
            USE_ABORT = .TRUE.
            !
        END IF
        !
    END IF
    !
    CALL PAR_MESSAGE(6, MESSAGE, ALLWRITE = .FALSE.)
    CALL PAR_MESSAGE(6, FOOTER, ALLWRITE = .FALSE.)
    !
    IF (USE_ABORT) THEN
        !
        CALL MPI_ABORT(MPI_COMM_WORLD, ENVERR, IERR)
        !
    ELSE
        !
        CALL MPI_FINALIZE(IERR)
        !
    END IF
    !
    STOP
    !
    END SUBROUTINE PAR_QUIT
    !
    !===========================================================================
    !
    SUBROUTINE PAR_MESSAGE(UNIT, MESSAGE, ALLWRITE)
    !
    ! Write message, either from master process or from all.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! UNIT: Unit number to write to
    ! MESSAGE: Message to write
    ! ALLWRITE: Flags whether all processes write or just master
    !         : default = FALSE [only master node writes message]
    !
    INTEGER, INTENT(IN) :: UNIT
    CHARACTER(LEN=*), INTENT(IN) :: MESSAGE
    LOGICAL, INTENT(IN), OPTIONAL :: ALLWRITE
    !
    ! Locals:
    !
    LOGICAL :: AW
    !
    !---------------------------------------------------------------------------
    !
    IF (PRESENT(ALLWRITE)) THEN
        !
        AW = ALLWRITE
        !
    ELSE
        !
        AW = .FALSE.
        !
    END IF
    !
    IF (AW) THEN
        !
        WRITE(UNIT, '(a)') 'Info   : [Rank '//TRIM(MYIDSTR(1:5))//'] - '//&
            &MESSAGE
        !
    ELSE
        !
        IF (MYID == 0) THEN
            !
            WRITE(UNIT, '(a)') MESSAGE
            !
        END IF
        !
    END IF
    !
    END SUBROUTINE PAR_MESSAGE
    !
    !===========================================================================
    !
    SUBROUTINE MRE_DECOMP1D(N, NUMPROCS, MYID, S, E)
    !
    ! This file contains a routine for producing a decomposition of a 1-d array
    !   when given a number of processors.  It may be used in "direct" product
    !   decomposition. The values returned assume a "global" domain in [1:n]
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! N: Size of array
    ! NUMPROCS: Number of processes
    ! MYID: My process number (0 to NUMPROCS-1 )
    ! S, E: Start and end locations for this process
    !
    INTEGER :: N
    INTEGER :: NUMPROCS
    INTEGER :: MYID
    INTEGER :: S
    INTEGER :: E
    !
    ! Locals:
    !
    INTEGER :: NLOCAL
    INTEGER :: DEFICIT
    !
    !---------------------------------------------------------------------------
    !
    NLOCAL = N / NUMPROCS
    S = MYID * NLOCAL + 1
    DEFICIT = MOD(N, NUMPROCS)
    S = S + MIN(MYID, DEFICIT)
    !
    IF (MYID .LT. DEFICIT) THEN
        !
        NLOCAL = NLOCAL + 1
        !
    END IF
    !
    E = S + NLOCAL - 1
    !
    IF (E .GT. N .OR. MYID .EQ. NUMPROCS - 1) THEN
        !
        E = N
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE MRE_DECOMP1D
    !
END MODULE PARALLEL_MOD
