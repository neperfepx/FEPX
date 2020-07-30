! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE surface_mod
  !
  ! Module for using surface information.
  !
  USE gather_scatter
  USE IntrinsicTypesModule, RK=>REAL_KIND
  !
  IMPLICIT NONE
  !
  INTEGER, PARAMETER :: NDIM3=3
  !
  ! The type `surface_section' is for use in parallel codes,
  ! and contains a section of the surface mesh.
  !
  TYPE surface_section
     !
     ! type: currently a number indicating nodes per element
     ! nsel: number of surface elements in this section
     ! semin to semax: range of surface element numbers
     ! 
     INTEGER :: type, semin, semax
     !
     ! conn: connectivity
     ! conn3d: dof connectivity
     ! econn: elemental connectivity to be used for gather/scatter stress
     ! 
     INTEGER, POINTER, DIMENSION(:,:)  :: conn, conn3d, econn
     !
     ! crds: coordinate array
     ! 
     REAL(RK), POINTER, DIMENSION(:,:) :: crds
     ! tm08:
     INTEGER, POINTER, DIMENSION(:)    :: elem
     !
     ! tr: the trace data structure for gather/scatter operations
     ! etr: elemental trace for stress
     ! 
     TYPE(trace) :: tr, etr
     !
  END TYPE surface_section
  !
CONTAINS
  !
  !
  FUNCTION allocate_surface_section(type, semin, semax, surf) RESULT(status) 
    !
    ! Allocate space to enough elements.
    !
    IMPLICIT NONE
    !
    !   Arguments:
    !
    !   type : surface element type, now only q6 is available
    !   semin, semax: range of surface element numbers
    !   surf : the surface to be allocated
    !
    INTEGER :: type, semin, semax
    TYPE(surface_section) :: surf
    INTEGER :: status 
    !
    !   Locals:
    !

    !
    !----------------------------------------------------------------------
    !
    status = 0
    !
    !if (type /= 4) then
    ! tsh for 3-node triangle
    !if (type /= 3) then
    ! tsh for 6-node triangle
    if (type /= 6) then
       status = 1
       RETURN
    endif
    surf%type = type  !type=6
    !
    if (semax >= semin) then
       surf%semin = semin
       surf%semax = semax
    else
       status = 2 
       RETURN
    endif
    !
    ALLOCATE(surf%conn(0:(type-1), semin:semax),STAT=status)
    ALLOCATE(surf%conn3d(0:(3*type-1), semin:semax),STAT=status)
    ALLOCATE(surf%econn(0:5, semin:semax),STAT=status)
    ALLOCATE(surf%crds(0:(NDIM3*type-1), semin:semax),STAT=status)
    !
    ! tm08:
    ALLOCATE(surf%elem(semin:semax),STAT=status)
    !
    RETURN
    !
  END FUNCTION allocate_surface_section

  !
END MODULE surface_mod
!
