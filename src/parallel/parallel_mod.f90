! This file is part of the FEPX software package.
! Copyright (C) 1996-2023, DPLab, ACME Lab, CNRS.
! See the COPYING file in the top-level directory.

module parallel_mod

! Module for parallel variables and operations (mpi version). This module
!   provides information and elementary operations for the primary communicator
!   in an mpi-based parallel program.

! Contains subroutines:
! par_init: Performs mpi initializations. Sets variables `num_procs' and `myid'
! par_partition: Partition arrays for multiple processes
! par_sum: Sum scalars across partitions
! par_max: Find maximum of parts
! par_min: Find minimum of parts
! par_gather: Gather array data across partitions
! par_quit: Clean up and exit
! par_message: Write message, either from master process or from all
! mre_decomp1d:  This file contains a routine for producing a decomposition of a
!   1-d array when given a number of processors.  It may be used in "direct"
!   product decomposition. The values returned assume a "global" domain in [1:n]

  ! Do not change rk; otherwise also change mpi_double_precision
  use, intrinsic :: iso_fortran_env, only: rk => real64
  use general_mod

  implicit none
  include 'mpif-handles.h'

  private

! Public

  public :: par_init
  public :: par_partition
  public :: par_sum
  public :: par_sqrtsum
  public :: par_max
  public :: par_min
  public :: par_message
  public :: par_quit
  public :: par_gather
  public :: par_gatherv
  public :: par_gatherv_array
  public :: par_gatherv_tensor

contains

  subroutine par_init()

    !  Performs mpi initializations.  Sets global variables `num_procs' and
    !   `myid'.  This routine has no arguments.

    !---------------------------------------------------------------------------

    ! Locals:

    integer :: ierr

    !---------------------------------------------------------------------------

    call mpi_init(ierr)

    call mpi_comm_size(mpi_comm_world, num_procs, ierr)

    ! Get my position in this communicator.

    call mpi_comm_rank(mpi_comm_world, myid, ierr)

    return

  end subroutine par_init

  !===========================================================================

  subroutine par_partition(na, np, me, dl, du)

    ! Partition arrays for multiple processes. Ranges go from 0 to (na - 1).

    ! Note moved here from fepx.f90 (rq 03/13)
    ! Note: the dof array was partitioned to keep the degrees of freedom on the same
    !   processor as the corresponding nodes. This shouldn't be necessary, but there
    !   was some kind of bug when i did it the other way. (deb)

    !---------------------------------------------------------------------------

    ! Arguments:
    ! na: Size of array
    ! np: Number of processes
    ! me: My process number
    ! dl: Lower dimension of my array segment
    ! du: Upper dimension of my array segment

    integer, intent(in) :: na
    integer, intent(in) :: np
    integer, intent(in) :: me
    integer, intent(out) :: dl
    integer, intent(out) :: du

    !---------------------------------------------------------------------------

    call mre_decomp1d(na, np, me, dl, du)

    return

  end subroutine par_partition

  !===========================================================================

  subroutine par_sqrtsum(part, whole)

    real(rk), intent(in) :: part
    real(rk), intent(out) :: whole

    call par_sum(part, whole)

    whole = dsqrt(whole)

  end subroutine par_sqrtsum

  subroutine par_sum(part, whole)

    ! Sum scalars across partitions.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! part: This processor's part
    ! whole: The sum of all parts

    real(rk), intent(in) :: part
    real(rk), intent(out) :: whole

    ! Locals:

    integer :: ierr

    !---------------------------------------------------------------------------

    call mpi_allreduce(part, whole, 1, mpi_double_precision, mpi_sum, &
        & mpi_comm_world, ierr)

    return

  end subroutine par_sum

  !===========================================================================

  subroutine par_max(part, maxpart)

    ! Find maximum of parts.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! part: This processor's part
    ! maxpart: The maximum of all parts

    real(rk) :: part
    real(rk) :: maxpart

    ! Locals:

    integer :: ierr

    !---------------------------------------------------------------------------

    call mpi_allreduce(part, maxpart, 1, mpi_double_precision, mpi_max, &
        & mpi_comm_world, ierr)

    return

  end subroutine par_max

  !===========================================================================

  subroutine par_min(part, minpart)

    ! Find minimum of parts.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! part: This processor's part
    ! minpart: The minimum of all parts

    real(rk) :: part
    real(rk) :: minpart

    ! Locals:

    integer :: ierr

    !---------------------------------------------------------------------------

    call mpi_allreduce(part, minpart, 1, mpi_double_precision, mpi_min, &
        & mpi_comm_world, ierr)

    return

  end subroutine par_min

  !===========================================================================

  subroutine par_gather(part, full, n)

    !  Gather array data across partitions.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! part: This processor's part
    ! full: The array containing all parts
    ! n: The amount of data being sent and recieved (length of array)

    integer :: part(:)
    integer :: full(:, :)
    integer :: n

    ! Locals:

    integer :: ierr

    !---------------------------------------------------------------------------

    call mpi_gather(part, n, mpi_integer, full, n, mpi_integer, 0, &
        & mpi_comm_world, ierr)

    return

  end subroutine par_gather

  !===========================================================================

  subroutine par_gatherv(part, full,n,m)

    !  Gather array data across partitions.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! part: This processor's part
    ! full: The array containing all parts

    real(rk) :: part(:)
    real(rk) :: full(:)
    integer :: n(:),m(:)
    ! Locals:
    ! n: The amount of data being sent and recieved (length of array)
    ! m: The starting index of data being sent and recieved 

    integer :: len,ierr
    len = size(part(:))
    
    !---------------------------------------------------------------------------
    call mpi_gatherv(part, len,mpi_double_precision, full,n,m &
      & , mpi_double_precision, 0,mpi_comm_world, ierr)

    return

  end subroutine par_gatherv

  !===========================================================================

  subroutine par_gatherv_array(part, full,n,m)

    !  Gather array data across partitions.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! part: This processor's part
    ! full: The array containing all parts

    real(rk) :: part(:,:)
    real(rk) :: full(:,:)
    integer :: n(:),m(:)
    real(rk) :: flat_part(size(part))
    real(rk) :: flat_full(size(full))
    ! Locals:
    ! n: The amount of data being sent and recieved (length of array)

    integer :: len,ierr
    len = size(part)
    flat_part = reshape(part,shape(flat_part))
    !---------------------------------------------------------------------------
    call mpi_gatherv(flat_part, len,mpi_double_precision, flat_full,n,m &
      & , mpi_double_precision, 0,mpi_comm_world, ierr)
    
    full = reshape(flat_full, shape(full))
    return

  end subroutine par_gatherv_array

  !===========================================================================

  subroutine par_gatherv_tensor(part, full,n,m)

    !  Gather array data across partitions.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! part: This processor's part
    ! full: The array containing all parts

    real(rk) :: part(:,:,:)
    real(rk) :: full(:,:,:)
    real(rk) :: flat_part(size(part))
    real(rk) :: flat_full(size(full))

    integer :: n(:),m(:)
    ! Locals:
    ! n: The amount of data being sent and recieved (length of array)
    integer :: len,ierr
    len = size(part)
    flat_part = reshape(part,shape(flat_part))
    !---------------------------------------------------------------------------
    call mpi_gatherv(flat_part, len,mpi_double_precision, flat_full,n,m &
      & , mpi_double_precision, 0,mpi_comm_world, ierr)
    
    full = reshape(flat_full, shape(full))
    return

  end subroutine par_gatherv_tensor

  !===========================================================================

  subroutine par_quit(message, abort)

    !  Clean up and exit.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! message: Message to write
    ! abort: Flag that forces mpi to abort over finalize

    character*(*), intent(in) :: message
    logical, intent(in), optional :: abort

    ! Locals:

    integer :: ierr, enverr = 1
    logical :: use_abort
    character(len=72), parameter :: footer = "============================&
        &============================================"

    !---------------------------------------------------------------------------

    use_abort = .false.

    if (present(abort)) then
      if (abort) then
        use_abort = .true.
      end if
    end if

    call par_message(6, message, allwrite=.false.)
    call par_message(6, footer, allwrite=.false.)

    if (use_abort) then
      call mpi_abort(mpi_comm_world, enverr, ierr)

    else
      call mpi_finalize(ierr)
    end if

    stop

  end subroutine par_quit

  !===========================================================================

  subroutine par_message(unit, message, allwrite)

    ! Write message, either from master process or from all.

    !---------------------------------------------------------------------------

    ! Arguments:
    ! unit: Unit number to write to
    ! message: Message to write
    ! allwrite: Flags whether all processes write or just master
    !         : default = false [only master node writes message]

    integer, intent(in) :: unit
    character(len=*), intent(in) :: message
    logical, intent(in), optional :: allwrite

    ! Locals:

    logical :: aw
    character(len=8) :: myidstr

    write (myidstr, '(i5)') myid
    myidstr = adjustl(myidstr)

    !---------------------------------------------------------------------------

    if (present(allwrite)) then
      aw = allwrite

    else
      aw = .false.
    end if

    if (aw) then
      write (unit, '(a)') 'Info   : [Rank '//trim(myidstr(1:5))//'] - '//&
          &message

    else
      if (myid .eq. 0) then
        write (unit, '(a)') message
      end if
    end if

  end subroutine par_message

  !===========================================================================

  subroutine mre_decomp1d(n, num_procs, myid, s, e)

    ! This file contains a routine for producing a decomposition of a 1-d array
    !   when given a number of processors.  It may be used in "direct" product
    !   decomposition. The values returned assume a "global" domain in [1:n]

    !---------------------------------------------------------------------------

    ! Arguments:
    ! n: Size of array
    ! num_procs: Number of processes
    ! myid: My process number (0 to num_procs-1 )
    ! s, e: Start and end locations for this process

    integer :: n
    integer :: num_procs
    integer :: myid
    integer :: s
    integer :: e

    ! Locals:

    integer :: nlocal
    integer :: deficit

    !---------------------------------------------------------------------------

    nlocal = n/num_procs
    s = myid*nlocal + 1
    deficit = mod(n, num_procs)
    s = s + min(myid, deficit)

    if (myid .lt. deficit) then
      nlocal = nlocal + 1
    end if

    e = s + nlocal - 1

    if (e .gt. n .or. myid .eq. num_procs - 1) then
      e = n
    end if

    return

  end subroutine mre_decomp1d

end module parallel_mod
