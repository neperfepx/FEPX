! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE SURFACE_MOD
!
! Module for using surface information.
!
! Contains function:
! ALLOCATE_SURFACE_SECTION: Allocate space to enough elements
! COMPUTE_AREA: Compute surface area
! INIT_SURF: Initialize surface arrays
! UPD_SURF: Update surface information
!
! From libf95:
!
USE INTRINSIC_TYPES_MOD, ONLY:RK=>REAL_KIND
!
! From libparallel
!
USE GATHER_SCATTER_MOD
USE QUADRATURE_MOD
USE SHAPE_2D_MOD
!
! From libparallel:
!
USE PARALLEL_MOD
!
IMPLICIT NONE
!
! Public
!
! The type `surface_section' is for use in parallel codes, and contains a
!   section of the surface mesh.
!
TYPE SURFACE_SECTION
    !
    ! TYPE: Currently a number indicating nodes per element
    ! NSEL: Number of surface elements in this section
    ! SEMIN to SEMAX: Range of surface element numbers
    ! CONN: Connectivity
    ! CONN3D: DOF connectivity
    ! ECONN: Elemental connectivity to be used for gather/scatter stress
    ! CRDS: Coordinate array
    ! ELEM:
    ! TR: The trace data structure for gather/scatter operations
    ! ETR: Elemental trace for stress
    !
    INTEGER :: TYPE, SEMIN, SEMAX
    INTEGER, POINTER, DIMENSION(:,:) :: CONN
    INTEGER, POINTER, DIMENSION(:,:) :: CONN3D
    INTEGER, POINTER, DIMENSION(:,:) :: ECONN
    REAL(RK), POINTER, DIMENSION(:,:) :: CRDS
    INTEGER, POINTER, DIMENSION(:) :: ELEM
    TYPE(TRACE) :: TR
    TYPE(TRACE) :: ETR
    !
END TYPE SURFACE_SECTION
!
INTEGER, PARAMETER:: MAXSURFACES = 6
TYPE(SURFACE_SECTION) :: SURFACES(MAXSURFACES)
INTEGER :: NSURFACES
CHARACTER(LEN=12) :: FASET(6)
!
! 3 node triangular element
!INTEGER, PARAMETER :: SFTYPE = 3
! 6 node triangular element
INTEGER, PARAMETER :: SFTYPE = 6
! 4 node quadralateral element
!INTEGER, PARAMETER :: SFTYPE = 4
REAL(RK) :: SFQP2D(SFTYPE, MAXQP2D)
REAL(RK) :: SFGQP2D(NDIM2, SFTYPE, MAXQP2D)
!
CONTAINS
    !
    FUNCTION ALLOCATE_SURFACE_SECTION(TYPE, SEMIN, SEMAX, SURF) RESULT(STATUS)
    !
    ! Allocate space to enough elements.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! TYPE: Surface element type, now only q6 is available
    ! SEMIN, SEMAX: Range of surface element numbers
    ! SURF: The surface to be allocated
    !
    INTEGER :: TYPE
    INTEGER :: SEMIN
    INTEGER :: SEMAX
    TYPE(SURFACE_SECTION) :: SURF
    INTEGER :: STATUS
    !
    !---------------------------------------------------------------------------
    !
    STATUS = 0
    !
    !if (type /= 4) then
    ! tsh for 3-node triangle
    !if (type /= 3) then
    ! tsh for 6-node triangle
    IF (TYPE /= 6) THEN
        !
        STATUS = 1
        !
        RETURN
    END IF
    SURF%TYPE = TYPE ! TYPE=6
    !
    IF (SEMAX >= SEMIN) THEN
        !
        SURF%SEMIN = SEMIN
        SURF%SEMAX = SEMAX
        !
    ELSE
        !
        STATUS = 2
        !
        RETURN
        !
    END IF
    !
    ALLOCATE(SURF%CONN(0:(TYPE-1), SEMIN:SEMAX), STAT = STATUS)
    ALLOCATE(SURF%CONN3D(0:(3*TYPE-1), SEMIN:SEMAX), STAT = STATUS)
    ALLOCATE(SURF%ECONN(0:5, SEMIN:SEMAX), STAT = STATUS)
    ALLOCATE(SURF%CRDS(0:(NDIM3*TYPE-1), SEMIN:SEMAX), STAT = STATUS)
    ALLOCATE(SURF%ELEM(SEMIN:SEMAX), STAT = STATUS)
    !
    RETURN
    !
    END FUNCTION ALLOCATE_SURFACE_SECTION
    !
    !===========================================================================
    !
    SUBROUTINE COMPUTE_AREA(COORDS, AREA0)
    !
    ! Compute surface area
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! COORDS: Coordinates
    ! AREA0: Initial area (enters as 0)
    !
    REAL(RK), INTENT(IN) :: COORDS(:)
    REAL(RK), INTENT(INOUT) :: AREA0(NSURFACES)
    !
    ! Locals:
    !
    INTEGER :: I
    INTEGER :: IEL
    INTEGER :: JQ
    INTEGER :: JN
    INTEGER :: JDOF
    REAL(RK) :: TANGENT(3, 2, MAXQP2D)
    REAL(RK) :: NORMAL(3, MAXQP2D)
    REAL(RK) :: NMAG
    REAL(RK) :: SJAC(MAXQP2D)
    REAL(RK) :: P_AREA, AREA
    !
    !---------------------------------------------------------------------------
    !
    AREA0 = 0.0D0
    !
    DO I = 1, NSURFACES
        !
        CALL PART_GATHER(SURFACES(I)%CRDS, COORDS,SURFACES(I)%CONN3D, &
            & SURFACES(I)%TR)
        !
        ! Now compute normal vectors at quadrature points
        !
        P_AREA = 0.0D0
        !
        DO IEL = SURFACES(I)%SEMIN, SURFACES(I)%SEMAX
            !
            TANGENT = 0.0D0
            !
            DO JQ = 1, NQP2D
                !
                DO JN = 1, SFTYPE
                    !
                    JDOF = 3 * (JN - 1)
                    !
                    ! First tangent vector
                    !
                    TANGENT(1, 1, JQ) = TANGENT(1, 1, JQ) + &
                        & SURFACES(I)%CRDS(JDOF, IEL) * SFGQP2D(1, JN, JQ)
                    TANGENT(2, 1, JQ) = TANGENT(2, 1, JQ) + &
                        & SURFACES(I)%CRDS(JDOF + 1, IEL) * SFGQP2D(1, JN, JQ)
                    TANGENT(3, 1, JQ) = TANGENT(3, 1, JQ) + &
                        & SURFACES(I)%CRDS(JDOF + 2, IEL) * SFGQP2D(1, JN, JQ)
                    !
                    ! Second tangent vector.
                    !
                    TANGENT(1, 2, JQ) = TANGENT(1, 2, JQ) + &
                        & SURFACES(I)%CRDS(JDOF, IEL) * SFGQP2D(2, JN, JQ)
                    TANGENT(2, 2, JQ) = TANGENT(2, 2, JQ) + &
                        & SURFACES(I)%CRDS(JDOF + 1, IEL) * SFGQP2D(2, JN, JQ)
                    TANGENT(3, 2, JQ) = TANGENT(3, 2, JQ) + &
                        & SURFACES(I)%CRDS(JDOF + 2, IEL) * SFGQP2D(2, JN, JQ)
                    !
                END DO
                !
                ! Now compute normals.
                !
                NORMAL(1, JQ) = TANGENT(2, 1, JQ) * TANGENT(3, 2, JQ) - &
                    & TANGENT(3, 1, JQ) * TANGENT(2, 2, JQ)
                NORMAL(2, JQ) = TANGENT(3, 1, JQ) * TANGENT(1, 2, JQ) - &
                    & TANGENT(1, 1, JQ) * TANGENT(3, 2, JQ)
                NORMAL(3, JQ) = TANGENT(1, 1, JQ) * TANGENT(2, 2, JQ) - &
                    & TANGENT(2, 1, JQ) * TANGENT(1, 2, JQ)
                NMAG = DSQRT(NORMAL(1, JQ) * NORMAL(1, JQ) + &
                    & NORMAL(2, JQ) * NORMAL(2, JQ) + &
                    & NORMAL(3, JQ) * NORMAL(3, JQ) )
                !
                IF (NMAG > 0) THEN
                    !
                    NORMAL(:, JQ) = NORMAL(:, JQ) / NMAG
                    SJAC(JQ) = NMAG
                    !
                ELSE
                    !
                    CALL PAR_QUIT('Error  :     > Surface normal zero magnitude.')
                    !
                END IF
                !
            END DO
            !
            DO JQ = 1, NQP2D
                !
                P_AREA = P_AREA + WT2D(JQ) * SJAC(JQ)
                !
            END DO
            !
        END DO !IEL
        !
        CALL PAR_SUM(P_AREA, AREA)
        !
        AREA0(I) = AREA
        !
    END DO !NSURFACES
    !
    RETURN
    !
    END SUBROUTINE COMPUTE_AREA
    !
    !===========================================================================
    !
    SUBROUTINE INIT_SURF()
    !
    ! Initialize surface arrays
    !
    !---------------------------------------------------------------------------
    !
    !   Locals:
    !
    INTEGER :: ISTAT
    !
    !---------------------------------------------------------------------------
    !
    CALL SF2D(SFTYPE, NQP2D, QP2D, SFQP2D, SFTYPE, ISTAT)
    CALL SF2DG(SFTYPE, NQP2D, QP2D, SFGQP2D, SFTYPE, ISTAT)
    !
    RETURN
    !
    END SUBROUTINE INIT_SURF
    !
    !===========================================================================
    !
    SUBROUTINE UPD_SURF(UNIT, COORDS, SIG_ALL, LOAD_OUT, AREA_OUT)
    !
    ! Update surface information.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! UNIT:
    ! COORDS: Coordinates
    ! SIG_ALL: Stresses
    ! LOAD_OUT: Output loads (enters as 0)
    ! AREA_OUT: Output area (enters as 0)
    !
    INTEGER :: UNIT !--> DEBUG_U
    REAL(RK), INTENT(IN) :: COORDS(:)
    REAL(RK), INTENT(IN) :: SIG_ALL (:)
    REAL(RK), INTENT(INOUT) :: LOAD_OUT(NSURFACES,3)
    REAL(RK), INTENT(INOUT) :: AREA_OUT(NSURFACES)
    !
    ! Locals:
    !
    INTEGER :: i
    INTEGER :: IEL
    INTEGER :: JQ
    INTEGER :: JN
    INTEGER :: JDOF
    INTEGER :: ILOAD
    INTEGER :: j
    REAL(RK) :: TANGENT(3, 2, MAXQP2D)
    REAL(RK) :: NORMAL(3, MAXQP2D)
    REAL(RK) :: NMAG
    REAL(RK) :: SJAC(MAXQP2D)
    REAL(RK), ALLOCATABLE:: SIG_SE(:,:)
    REAL(RK) :: LOAD(3)
    REAL(RK) :: P_TOTAL_LOAD(3)
    REAL(RK) :: TOTAL_LOAD(3)
    REAL(RK) :: P_AREA
    REAL(RK) :: AREA
    !
    !----------------------------------------------------------------------
    !
    ! write(UNIT,*) 'Updating surface info.'
    !
    LOAD_OUT = 0.0D0
    AREA_OUT = 0.0D0
    !
    !
    DO I = 1, NSURFACES
        !
        CALL PART_GATHER(SURFACES(I)%CRDS, COORDS, SURFACES(I)%CONN3D, &
            & SURFACES(I)%TR)
        !
        ALLOCATE(SIG_SE(0:5, SURFACES(I)%SEMIN:SURFACES(I)%SEMAX))
        !
        CALL PART_GATHER(SIG_SE, SIG_ALL, SURFACES(I)%ECONN, SURFACES(I)%ETR)
        !
        ! Now compute normal vectors at quadrature points.
        !
        P_TOTAL_LOAD = 0.0D0
        P_AREA = 0.0D0
        !
        DO IEL = SURFACES(I)%SEMIN, SURFACES(I)%SEMAX
            !
            TANGENT = 0.0D0
            !
            DO JQ = 1, NQP2D
                !
                DO JN = 1, SFTYPE !number of nodes in the surface element
                    !
                    JDOF = 3 * (JN - 1)
                    !
                    ! First tangent vector.
                    !
                    TANGENT(1, 1, JQ) = TANGENT(1, 1, JQ) + &
                        & SURFACES(I)%CRDS(JDOF, IEL) * SFGQP2D(1, JN, JQ)
                    TANGENT(2, 1, JQ) = TANGENT(2, 1, JQ) + &
                        & SURFACES(I)%CRDS(JDOF + 1, IEL) * SFGQP2D(1, JN, JQ)
                    TANGENT(3, 1, JQ) = TANGENT(3, 1, JQ) + &
                        & SURFACES(I)%CRDS(JDOF + 2, IEL) * SFGQP2D(1, JN, JQ)
                    !
                    ! Second tangent vector.
                    !
                    TANGENT(1, 2, JQ) = TANGENT(1, 2, JQ) + &
                         & SURFACES(I)%CRDS(JDOF, IEL) * SFGQP2D(2, JN, JQ)
                    TANGENT(2, 2, JQ) = TANGENT(2, 2, JQ) + &
                         & SURFACES(I)%CRDS(JDOF + 1, IEL) * SFGQP2D(2, JN, JQ)
                    TANGENT(3, 2, JQ) = TANGENT(3, 2, JQ) + &
                         & SURFACES(I)%CRDS(JDOF + 2, IEL)*SFGQP2D(2, JN, JQ)
                    !
                END DO
                !
                ! Now compute normals
                !
                NORMAL(1, JQ) = TANGENT(2, 1, JQ) * TANGENT(3, 2, JQ) - &
                    & TANGENT(3, 1, JQ) * TANGENT(2, 2, JQ)
                NORMAL(2, JQ) = TANGENT(3, 1, JQ) * TANGENT(1, 2, JQ) - &
                    & TANGENT(1, 1, JQ) * TANGENT(3, 2, JQ)
                NORMAL(3, JQ) = TANGENT(1, 1, JQ) * TANGENT(2, 2, JQ) - &
                    & TANGENT(2, 1, JQ) * TANGENT(1, 2, JQ)
                NMAG = DSQRT(NORMAL(1, JQ) * NORMAL(1, JQ) + &
                    & NORMAL(2, JQ) * NORMAL(2, JQ) + &
                    & NORMAL(3, JQ) * NORMAL(3, JQ))
                !
                IF (NMAG > 0) THEN
                    !
                    NORMAL(:, JQ) = NORMAL(:, JQ) / NMAG
                    SJAC(JQ) = NMAG
                    !
                ELSE
                    !
                    CALL PAR_QUIT('Error  :     > Surface normal zero magnitude.')
                    !
                END IF
                !
            END DO
            !
            !write(UNIT, *) 'element: ', IEL
            !write(UNIT, '(3e14.4)') NORMAL
            !
            LOAD = 0.0d0
            !
            DO JQ = 1, NQP2D
                !
                LOAD(1) = LOAD(1) + WT2D(JQ) * SJAC(JQ) * (SIG_SE(0, IEL) * &
                    & NORMAL(1, JQ) + SIG_SE(1, IEL) * NORMAL(2, JQ) + &
                    & SIG_SE(2, IEL) * NORMAL(3, JQ))
                LOAD(2) = LOAD(2) + WT2D(JQ) * SJAC(JQ) * (SIG_SE(1, IEL) * &
                    & NORMAL(1, JQ) + SIG_SE(3, IEL) * NORMAL(2, JQ) + &
                    & SIG_SE(4, IEL) * NORMAL(3, JQ))
                LOAD(3) = LOAD(3) + WT2D(JQ) * SJAC(JQ) * (SIG_SE(2, IEL) * &
                    & NORMAL(1, JQ) + SIG_SE(4, IEL) * NORMAL(2, JQ) + &
                    & SIG_SE(5, IEL) * NORMAL(3, JQ))
                P_AREA = P_AREA + WT2D(JQ) * SJAC(JQ)
                !
            END DO
            !
            P_TOTAL_LOAD = P_TOTAL_LOAD + LOAD
            !
        END DO !IEL
        !
        DO ILOAD = 1,3
            !
            CALL PAR_SUM(P_TOTAL_LOAD(ILOAD), TOTAL_LOAD(ILOAD))
            !
        END DO
        !
        CALL PAR_SUM(P_AREA, AREA)
        !
        !IF (myid .eq. 0) THEN
        !  write(UNIT, '(a9,i3,a1)', ADVANCE='NO') 'LOAD surf', i, ':'
        !  write(UNIT, '(5x,3e14.4)') TOTAL_LOAD
        !  write(UNIT, *) 'AREA = ', AREA
        !endif
        !
        LOAD_OUT(I, :) = TOTAL_LOAD
        AREA_OUT(I) = AREA
        !
        DEALLOCATE(SIG_SE)
        !
    END DO !NSURFACES
    !
    RETURN
    !
    END SUBROUTINE UPD_SURF
    !
END MODULE SURFACE_MOD
