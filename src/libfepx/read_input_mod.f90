! This file is part of the FEPX software package.
! Copyright (C) 1996-2020, DPLab, ACME Lab.
! See the COPYING file in the top-level directory.
!
MODULE READ_INPUT_MOD
!
! Module to handle all generic file read-in functionality. Specific file
! handling (restart, BCs, and powder-diffraction) are maintained in their
! own respective modules for ease of modification.
!
! Contains subroutines:
! ALLOCATE_MSH: Allocates mesh arrays based on the problem size.
! GET_MSH_SIZE: Scrapes the msh file and retrieve the number of nodes and elems.
! READ_SPATIAL_MSH: Parent subroutine to handle msh file parsing process.
!
! Contains helper subroutines (supporting READ_SPATIAL_MSH):
! GET_NODE_INFO: Gets the number of nodes in the mesh.
! GET_ELEM_INFO: Gets the number of elements in the mesh and partitions.
! READ_MESH_FORMAT: Reads in gmsh file format - 2.2 0 8 only.
! READ_NODES: Read in node ID and [X Y Z] spatial coordinates.
! READ_ELEMENTS: Read in ALL elements and only store 3D elements.
! READ_NSETS: Read in node sets - currently not stored for read-in BCs.
! READ_FASETS: Read in surface face node sets (2D elements).
! READ_NODEPARTITIONS: Read in per-node partition distribution. -OPTIONAL
! READ_PHYSICALNAMES: Read in physical names and do not store (gmsh data only).
! READ_ORIENTATIONS: Read in grain (or element) orientations.
! READ_GROUPS: Read in grain phases (assumes single-phase by default). -OPTIONAL
!
USE INTRINSICTYPESMODULE, ONLY: RK=>REAL_KIND
USE PARALLEL_MOD
! 
USE SURFACE_MOD
USE SURF_INFO_MOD
USE SIMULATION_CONFIGURATION_MOD
USE UNITS_MOD
!
IMPLICIT NONE
!
! Public
!
PUBLIC
!
! Element/DOF array bounds
!
INTEGER, PARAMETER :: NDIM  = 10 ! 10-noded tetrahedra
INTEGER, PARAMETER :: NNPE  = NDIM - 1
INTEGER, PARAMETER :: KDIM  = 3*NDIM
INTEGER, PARAMETER :: KDIM1 = KDIM - 1
!
! Array sizes and ranges.
!  
INTEGER :: NUMELM, NUMNP, MAXEL, MAXEL1, MAXNP, MAXNP1, MAXDOF, MAXDOF1
INTEGER :: EL_SUB1, EL_SUP1, NP_SUB1, NP_SUP1, DOF_SUB1, DOF_SUP1
!
! Connectivities and coordinates.
!
INTEGER,  ALLOCATABLE :: NP(:,:), NODES(:,:)
INTEGER,  ALLOCATABLE :: UESUB(:), PHASE(:)
REAL(RK), ALLOCATABLE :: COORDS(:), ELEMENT_CRDS(:,:), &
    & COORDS_ORIG(:)
!
! Orientation options
!
LOGICAL :: ELEMENT_ORIS = .FALSE.
TYPE ORIENTATION_TYPE
    !
    CHARACTER(LEN=50) :: ORIENTATION_PARAMETERIZATION
    CHARACTER(LEN=50) :: ORIENTATION_CONVENTION
    !
END TYPE ORIENTATION_TYPE
TYPE(ORIENTATION_TYPE) :: ORIENTATION_OPTIONS
!
! Microstructure options
!
INTEGER :: NUM_GRAINS
REAL(RK), ALLOCATABLE :: GRAIN_ORIENTATION(:,:)
!
! Partition flags
!
LOGICAL :: READ_ELEM_PART = .FALSE.
LOGICAL :: READ_NODE_PART = .FALSE.
!
CONTAINS
    !
    SUBROUTINE ALLOCATE_MSH(STATUS)
    !
    ! Allocate mesh arrays according to problem size.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! STATUS: Returns status of allocation to main program.
    !
    INTEGER, INTENT(INOUT) :: STATUS
    !
    !---------------------------------------------------------------------------
    !
    STATUS = 0
    !
    ALLOCATE (NP(0:NNPE, EL_SUB1:EL_SUP1), &
        & NODES(0:KDIM1, EL_SUB1:EL_SUP1), &
        & COORDS(DOF_SUB1:DOF_SUP1), &
        & COORDS_ORIG(DOF_SUB1:DOF_SUP1), &
        & ELEMENT_CRDS(0:KDIM1, EL_SUB1:EL_SUP1), &
        & UESUB(0:NUMELM-1), &
        & PHASE(0:NUMELM-1), &
        & STAT=STATUS)
    !
    ! Initialize these arrays here before the spatial mesh is parsed.
    !
    ! Note: UESUB stores the 1-indexed elset for a element of ID:(0:NUMELM-1)
    ! and we set phase to 1 to assume a single-phase simulation by default.
    !
    UESUB = 0
    PHASE = 1
    !
    RETURN
    !
    END SUBROUTINE ALLOCATE_MSH
    !
    !===========================================================================
    !
    SUBROUTINE GET_MSH_SIZE(IO)
    !
    ! Scrapes $Nodes and $Elements fields for the number of nodes and 3D elems
    ! to determine problem size. Also, attempt to retrieve embedded mesh 
    ! partitioning info from $Elements and $NodePartitions (if available).
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! EOFSTAT: Value that confirms if FORTRAN is at the end of a file.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: EOFSTAT
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Initialize EOFSTAT and loop over the entire file until EOF is reached.
    EOFSTAT = 0
    !
    DO WHILE (EOFSTAT .EQ. 0)
        !
        ! Read in line and determine which section is to be parsed.
        ! The below line is a temporary fix to avoid program hangs - JC
        LINE = ""
        !
        READ(IO, '(A)', IOSTAT = EOFSTAT) LINE
        !
        SELECT CASE (LINE)
            !
            CASE('$Nodes')
                !
                BACKSPACE(IO)
                CALL GET_NODE_INFO(IO)
                !
            CASE('$Elements')
                !
                BACKSPACE(IO)
                CALL GET_ELEM_INFO(IO)
                !
            CASE('$NodePartitions')
                !
                BACKSPACE(IO)
                CALL READ_NODEPARTITIONS(IO)
                !
            ! END CASE            
            !
        END SELECT
        !
    END DO
    !
    REWIND (IO)
    !
    RETURN
    !
    END SUBROUTINE GET_MSH_SIZE
    !
    !===========================================================================
    !
    SUBROUTINE GET_NODE_INFO(IO)
    !
    ! Read the $Nodes field and retrieves the total number of nodes in mesh.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I: Generic looping index.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, I
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$Nodes')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &nodes.')
        !
    END IF
    !
    ! Read in the number of nodes in the section (and set global NUMNP).
    READ(IO, *) NUMNP
    !
    DO I = 1, NUMNP
        !
        READ(IO, *)
        !  
    END DO
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndNodes')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &nodes.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE GET_NODE_INFO
    !
    !===========================================================================
    !
    SUBROUTINE GET_ELEM_INFO(IO)
    !
    ! Read the $Elements field and retrieves the total number of nodes in mesh.
    ! Also, retrieves elemental partition information if present in the mesh.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I: Generic looping index.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! NSPACE/NVALS: Storage values used for dividing a line into substrings.
    ! ELMTYPE: Gmsh defined value that determines the element type.
    ! NELM: Locally determined number of 3D elements in mesh.
    ! MINBND/MAXBND: Array bounds to assist in trimming element partition array.
    ! NPARTS: Number of partitions read in from the mesh.
    ! IPART: Partition loop index to check which partition is being checked.
    ! LINE: Input line on current record to be parsed.
    ! IARRAY: Substring array for internal read parsing.
    ! TEMP: Temporary array storage for element paritions including non-3D.
    ! ELEM_PARTS: Array storing the trimmed element paritions.
    ! TEMPVAL: Array of EL_SUB1/EL_SUP1 bounds for each partition read in.
    ! CHG_INDEX: Index in array immediately after a change in value occurs.
    ! MASK: Masking array to trim TEMP into ELEM_PARTS.
    !
    INTEGER :: IERR, I
    INTEGER :: NLINES, NSPACE, NVALS, ELMTYPE, NELM
    INTEGER :: MINBND, MAXBND, NPARTS, IPART
    CHARACTER(LEN=256) :: LINE
    CHARACTER(LEN=12)  :: IARRAY(16)
    INTEGER, ALLOCATABLE :: TEMP(:), ELEM_PARTS(:), TEMPVAL(:,:), CHG_INDEX(:)
    LOGICAL, ALLOCATABLE :: MASK(:)
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$Elements')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &elements.')
        !
    END IF
    !
    ! Read in the number of elements in the section (and set NLINES).
    READ(IO, *) NLINES
    !
    ! Initialize the 3D elem counter and partition array.
    NELM = 0
    !
    ALLOCATE(TEMP(1:NLINES))
    ALLOCATE(MASK(1:NLINES))
    TEMP = -1
    !
    DO I = 1, NLINES
        !
        READ(IO, '(A)') LINE
        !
        ! Trim the full string into NVALS number of substrings to parse.
        NSPACE = COUNT(TRANSFER(LINE, 'A', LEN_TRIM(LINE)) == " ")
        NVALS = NSPACE + 1
        IARRAY = ""
        !
        ! Internal read of line to store substrings.
        READ(LINE, *) IARRAY(1:NVALS)
        READ(IARRAY(2), *) ELMTYPE
        !
        IF (ELMTYPE .EQ. 11) THEN
            !
            ! Assign data from IARRAY(6) the TAG3 value which contain partition.
            READ(IARRAY(6), *) TEMP(I)
            !
            ! Increment NELM by 1 to ensure NUMELM is assigned correctly.
            NELM = NELM + 1
            !
        END IF
        !
    END DO
    !
    ! Set the global NUMELM value.
    NUMELM = NELM
    !
    ! Reallocate TEMP array to properly sized array (1:NUMELM).
    MASK = TEMP .EQ. -1
    MINBND = COUNT(MASK) + 1
    MAXBND = SIZE(TEMP)
    !
    ELEM_PARTS = TEMP(MINBND:MAXBND)
    DEALLOCATE(MASK)
    DEALLOCATE(TEMP)
    !
    ! Set a flag if we read in paritions and strip the ELEM_PARTS array down.
    IF (MAXVAL(ELEM_PARTS) .GT. 0) THEN
        !
        READ_ELEM_PART = .TRUE.
        NPARTS = MAXVAL(ELEM_PARTS)
        ALLOCATE(TEMPVAL(1:NPARTS, 2))
        ALLOCATE(CHG_INDEX(NPARTS-1))
        IPART = 1
        !
        ! We can split ELEM_PARTS in 0-indexed EL_SUB1 and EL_SUP1 here.
        DO I = 1, SIZE(ELEM_PARTS)
            !
            ! A change has been detected so do something.
            IF (ELEM_PARTS(I) .NE. IPART) THEN
                !
                CHG_INDEX(IPART) = I
                IF (IPART .NE. NPARTS) IPART = IPART + 1
                !
            END IF
            !
        END DO
        !
        ! DEBUG WRITE
!        WRITE(*,*) 'Number of elem partitions: ', NPARTS
!        WRITE(*,*) 'Element step indices: ', CHG_INDEX
        !
        ! Build bound array for NPARTS -> (1,:) is EL_SUB1, (2,:) is EL_SUP1.
        TEMPVAL = 0
        DO I = 1, NPARTS
            !
            ! Need to handle the edge cases differently.
            IF (I .EQ. 1) THEN
                !
                TEMPVAL(1,1) = 1
                TEMPVAL(1,2) = CHG_INDEX(I) - 1
                !
            ELSE IF (I .EQ. NPARTS) THEN
                !
                TEMPVAL(I,1) = CHG_INDEX(I-1)
                TEMPVAL(I,2) = SIZE(ELEM_PARTS)
                !
            ELSE
                !
                TEMPVAL(I,1) = CHG_INDEX(I-1)
                TEMPVAL(I,2) = CHG_INDEX(I) - 1
                !
            END IF
            !
        END DO
        !
        ! DEBUG WRITE
!        DO I = 1, NPARTS
!            WRITE(*,*) 'Partition ', I, ': ', TEMPVAL(I,:)
!        END DO
        !
    END IF
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndElements')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &elements.')
        !
    END IF
    !
    ! Do not need to maintain element partitions for now so deallocate.
    IF (ALLOCATED(ELEM_PARTS)) DEALLOCATE(ELEM_PARTS)
    IF (ALLOCATED(CHG_INDEX))  DEALLOCATE(CHG_INDEX)
    IF (ALLOCATED(TEMPVAL))    DEALLOCATE(TEMPVAL)
    ! 
    RETURN
    !
    END SUBROUTINE GET_ELEM_INFO
    !
    !===========================================================================
    !
    SUBROUTINE READ_SPATIAL_MSH(IO)    
    !
    ! Parse the .msh file and have each processor store only what it needs.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! EOFSTAT: Value that confirms if FORTRAN is at the end of a file.
    ! I: Generic looping index.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: EOFSTAT, I
    CHARACTER(LEN=256) :: LINE
    ! 
    ! Notes:
    ! The .msh file follows the Gmsh format 2.2 and is parsed in chunks
    ! and the field headers ($FieldName) are read to determine which
    ! helper subroutines to call.
    !
    !---------------------------------------------------------------------------
    !
    ! Initialize EOFSTAT and loop over the entire file until EOF is reached.
    EOFSTAT = 0
    !
    DO WHILE (EOFSTAT .EQ. 0)
        !
        ! Read in line and determine which section is to be parsed.
        ! The below line is a temporary fix to avoid program hangs - JC
        LINE = ""
        !
        READ(IO, '(A)', IOSTAT = EOFSTAT) LINE
        !
        SELECT CASE (LINE)
            !
            CASE('$MeshFormat')
                !
                BACKSPACE(IO)
                CALL READ_MESH_FORMAT(IO)
                !
            CASE('$Nodes')
                !
                BACKSPACE(IO)
                CALL READ_NODES(IO)
                !
            CASE('$Elements')
                !
                BACKSPACE(IO)
                CALL READ_ELEMENTS(IO)
                !
            CASE('$NSets') ! Nodal values for BCs assignments (OPTIONAL).
                !
                BACKSPACE(IO)
                CALL READ_NSETS(IO)
                !
            CASE('$Fasets')
                !
                BACKSPACE(IO)
                CALL READ_FASETS(IO)
                !
                DO I = 1, NSURFACES
                    !
                    CALL PART_GATHER(SURFACES(I)%CRDS, COORDS, &
                        &SURFACES(I)%CONN3D, SURFACES(I)%TR)
                    !
                ENDDO
                !
            CASE('$PhysicalNames')
                !
                BACKSPACE(IO)
                CALL READ_PHYSICALNAMES(IO)
                !
            CASE('$ElsetOrientations') ! Per-grain orientations.
                !
                BACKSPACE(IO)
                CALL READ_ELSETORIENTATIONS(IO)
                !
            CASE('$ElementOrientations') ! Per-element orientations (OPTIONAL).
                !
                BACKSPACE(IO)
                CALL READ_ELEMENTORIENTATIONS(IO)
                !
            CASE('$Groups') ! Grain/phase assignment for multiphase (OPTIONAL).
                !
                BACKSPACE(IO)
                CALL READ_GROUPS(IO)
                !
            ! END CASE            
            !
        END SELECT
        !
    END DO
    !
    RETURN
    !
    END SUBROUTINE READ_SPATIAL_MSH
    !
    !===========================================================================
    !
    SUBROUTINE READ_MESH_FORMAT(IO)
    !
    ! Read the mesh file format and confirm it is of correct type.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! NSPACE/NVALS: Storage values used for dividing a line into substrings.
    ! FILE_TYPE: Gmsh defined file type of mesh - must be '0'.
    ! DATA_SIZE: Gmsh defined data size of mesh - must be '8'.
    ! MESH_FORMAT: Gmsh defined mesh file format - must be '2.2'.
    ! IARRAY: Substring array for internal read parsing.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, NSPACE, NVALS
    INTEGER :: FILE_TYPE, DATA_SIZE
    CHARACTER(LEN=3)   :: MESH_FORMAT
    CHARACTER(LEN=12)  :: IARRAY(16)
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$MeshFormat')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &mesh format.')
        !
    END IF
    !
    ! Read in the mesh format string.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    ! Trim the full string into NVALS number of substrings to parse.
    NSPACE = COUNT(TRANSFER(LINE, 'A', LEN_TRIM(LINE)) == " ")
    NVALS = NSPACE + 1
    IARRAY = ""
    !
    ! Internal read of line to store substrings.
    READ(LINE, *) IARRAY(1:NVALS)
    READ(IARRAY(1), *) MESH_FORMAT
    READ(IARRAY(2), *) FILE_TYPE
    READ(IARRAY(3), *) DATA_SIZE
    !
    IF (MESH_FORMAT .NE. '2.2') &
        &CALL PAR_QUIT('Error  :     > Incorrect mesh format version provided.')
    IF (FILE_TYPE .NE. 0) &
        &CALL PAR_QUIT('Error  :     > Incorrect mesh file type provided.')
    IF (DATA_SIZE .NE. 8) &
        &CALL PAR_QUIT('Error  :     > Incorrect mesh data size provided.')
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndMeshFormat')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &mesh format.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_MESH_FORMAT
    !
    !===========================================================================
    !
    SUBROUTINE READ_NODES(IO)
    !
    ! Read in node ID and coordinates and store per processor.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! INODE: Node 1-indexed ID value.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! I: Generic loop index.
    ! K1/K2/K3: DOF mapping values for the node coordinates.
    ! X/Y/Z: Cartesian coordinate value of the node.
    ! IARRAY: Substring array for internal read parsing.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER  :: IERR, INODE, NLINES
    INTEGER  :: I, K1, K2, K3
    REAL(RK) :: X, Y, Z
    CHARACTER(LEN=12)  :: IARRAY(16)
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$Nodes')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &nodes.')
        !
    END IF
    !
    ! Read in the number of nodes in the section (and set NLINES).
    READ(IO, *) NLINES
    !
    DO I = 1, NLINES
        !
        READ(IO, *) INODE, X, Y, Z
        !
        ! Need to shift by 1 to store zero-indexed node ID for internal use.
        INODE = INODE - 1
        !
        ! Store the nodal coords only if within local processor range.
        IF ((INODE .GE. NP_SUB1) .AND. (INODE .LE. NP_SUP1)) THEN
            !
            K1 = 3 * INODE
            K2 = K1 + 1
            K3 = K2 + 1
            
            COORDS(K1) = X
            COORDS(K2) = Y
            COORDS(K3) = Z
            !
        END IF
        !  
    END DO
    !
    ! Make copy of original coordinates.
    !
    COORDS_ORIG = COORDS
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndNodes')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &nodes.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_NODES
    !
    !===========================================================================
    !
    SUBROUTINE READ_ELEMENTS(IO)
    !
    ! Read in all elements [0D -> 3D] and only store 3D elements within
    ! local processor range.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I/J: Generic loop index.
    ! J1/J2/K3/K1/K2/K3: DOF mapping values for the element nodes.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! NSPACE/NVALS: Storage values used for dividing a line into substrings.
    ! ELMTYPE: Gmsh defined value that determines the element type.
    ! NELM: Locally determined number of 3D elements in mesh.
    ! NTAGS: Number of tags in the 3D element line - must be 3.
    ! TAG1: Grain (or elset) ID - 1-indexed.
    ! NODES_FE: Local element nodal connectivity array.
    ! IARRAY: Substring array for internal read parsing.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER  :: IERR, I, J, J1, J2, J3, K1, K2, K3
    INTEGER  :: NLINES, NSPACE, NVALS, ELMTYPE, NELM, NTAGS, TAG1
    INTEGER  :: NODES_FE(0:NNPE)
    CHARACTER(LEN=12)  :: IARRAY(16)
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$Elements')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &elements.')
        !
    END IF
    !
    ! Read in the number of elements in the section (and set NLINES).
    READ(IO, *) NLINES
    !
    NELM = 0
    !
    DO I = 1, NLINES
        !
        READ(IO, '(A)') LINE
        !
        ! Trim the full string into NVALS number of substrings to parse.
        NSPACE = COUNT(TRANSFER(LINE, 'A', LEN_TRIM(LINE)) == " ")
        NVALS = NSPACE + 1
        IARRAY = ""
        !
        ! Internal read of line to store substrings.
        READ(LINE, *) IARRAY(1:NVALS)
        READ(IARRAY(2), *) ELMTYPE
        !
        IF (ELMTYPE .EQ. 11) THEN
            !
            ! If 10-node tetrahedral element, extract information from the line.
            ! The line has format of:
            ! IELEM ELMTYPE NTAGS TAG1 TAG2 TAG3 NODE0 ... NODE9
            READ(IARRAY(3), *) NTAGS
            READ(IARRAY(4), *) TAG1 ! elset
            !
            ! Set the UESUB value for this element (assigns elset or grain ID).
            UESUB(NELM) = TAG1
            !
            IF (NTAGS .EQ. 3) THEN
                !
                ! The .msh format modifies the local node order of an element.
                ! The remapping follows from Gmsh to FEPX:
                ! GMSH ORDER: 1 2 3 4 5 6 7 8 9 10
                ! FEPX ORDER: 1 3 5 10 2 4 6 7 9 8
                READ(IARRAY(7),*) NODES_FE(0)
                READ(IARRAY(8),*) NODES_FE(2)
                READ(IARRAY(9),*) NODES_FE(4)
                READ(IARRAY(10),*) NODES_FE(9)
                READ(IARRAY(11),*) NODES_FE(1)
                READ(IARRAY(12),*) NODES_FE(3)
                READ(IARRAY(13),*) NODES_FE(5)
                READ(IARRAY(14),*) NODES_FE(6)
                READ(IARRAY(15),*) NODES_FE(8)
                READ(IARRAY(16),*) NODES_FE(7)
                !
                ! Shift the array in order to maintain internal zero-indexing.
                NODES_FE(0:9) = NODES_FE(0:9) - 1
                !
                ! Store elements only on local processor range.
                IF ((NELM .GE. EL_SUB1) .AND. (NELM .LE. EL_SUP1)) THEN
                    !
                    NP(:, NELM) = NODES_FE
                    !
                    ! Map the local DOF to global DOF.
                    DO J = 0, NNPE
                        !
                        J1 = 3 * J
                        J2 = J1 + 1
                        J3 = J2 + 1
                        !
                        K1 = 3 * NODES_FE(J)
                        K2 = K1 + 1
                        K3 = K2 + 1
                        !
                        NODES(J1, NELM) = K1
                        NODES(J2, NELM) = K2
                        NODES(J3, NELM) = K3
                        !
                    END DO
                    !
                END IF
                !  
            END IF
            !
            ! Increment NELM by 1 to ensure NP is assigned correctly
            NELM = NELM + 1
            NUM_GRAINS = TAG1
            !
        END IF
        !
    END DO
    !
    ! Confirm that NELM (local) and NUMELM (global) values match
    IF (NELM .NE. NUMELM) CALL PAR_QUIT('Error  :     > &
        &Number of elements in $Elements does not match problem size.')
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndElements')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &elements.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_ELEMENTS
    !
    !===========================================================================
    !
    SUBROUTINE READ_NSETS(IO)
    !
    ! Read in node sets and do not store - currently unused, but will be 
    ! integrated into custom boundary conditions in the future.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I/J: Generic loop index.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! NSETS: Number of NSets within the field.
    ! NSETLABELS: Array of NSet labels for global use.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, I, J, NLINES, NSETS
    CHARACTER(LEN=12)  :: NSETLABELS(6)
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$NSets')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &nsets.')
        !
    END IF
    !
    ! Read in the number of nsets in the section (and set NSETS).
    READ(IO, *) NSETS
    !
    ! Loop over the number of NSETS and extract information per set.
    DO I = 1, NSETS
        !
        ! Read in the NSet label.
        READ(IO, *) NSETLABELS(I)
        !
        ! Read in the number of nodes in this set.
        READ(IO, *) NLINES
        !
        ! Loop over single set and do NOT store.
        DO J = 1, NLINES
            !
            READ(IO, *)
            !
        END DO
        !
    END DO
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndNSets')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &nsets.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_NSETS
    !
    !===========================================================================
    !
    SUBROUTINE READ_FASETS(IO)
    !
    ! Read in facesets and store the surface 2D elements. This subroutine also
    ! prepares the 3D connectivity for the surface nodes which is later used
    ! to integrate the loads on each surface. 
    !
    ! Notes:
    ! This information should be readily extracted from the READ_ELEMENTS
    ! subroutine by searching for ELM_TYPE (9) and extracting, however, this
    ! would not provide faset labels or IDs.
    !
    ! Legacy note: type = 6 for triangle - tsh
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I/J/K: Generic loop index.
    ! TYPE/TYPE1/TYPE3: Parameters to define the internal element type/bounds.
    ! ELEM_ID: Surface element ID - 1-indexed.
    ! ELEM_NODES/ELEM_NODES_TMP: Surface element connectivity.
    ! NSEL: Number of elements on a given surface.
    ! SEMIN/SEMAX: Partitioned array bounds to distribute surface sections.
    ! STATUS: Confirms if allocation of surface section arrays was successful.
    ! J3/N3: 3D connectivity DOF mapping values.
    ! IERR: Value that confirms if a READ() fails.
    ! EL_DOF_MIN/EL_DOF_MAX: DOF mapping to scatter surface elements.
    ! FASETS: Number of fasets in the field to be parsed.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER, PARAMETER :: TYPE = 6, TYPE1 = TYPE - 1, TYPE3 = 3*TYPE -1
    INTEGER :: ELEM_ID, ELEM_NODES(0:TYPE1), ELEM_NODES_TMP(0:TYPE1)
    INTEGER :: I, J, K, NSEL, SEMIN, SEMAX, STATUS, J3, N3, IERR
    INTEGER :: EL_DOF_MIN, EL_DOF_MAX, FASETS
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$Fasets')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &fasets.')
        !
    END IF
    !
    ! Read in the number of fasets in the section (and set FASETS).
    READ(IO, *) FASETS
    !
    ! This needs to be set here, but is in the scope of SURF_INFO_MOD currently.
    NSURFACES = FASETS 
    !
    ! Check if 6 surfaces are to be read in. This is hardwired for cubic geom.
    IF (FASETS .NE. 6) CALL PAR_QUIT('Error  :     > &
        &6 surfaces necessary in $Fasets.')
    !
    DO I = 1, FASETS
        !
        ! Read the faset label - unused currently.
        READ(IO,*) FASET(I)
        ! number of elements on the surface.
        READ(IO,*) NSEL
        !
        ! Should SEMIN/SEMAX use the the surface elements from the 3D elements
        ! that are stored on the same processor range instead of something else?
        CALL PAR_PARTITION(NSEL, NUMPROCS, MYID, SEMIN, SEMAX)
        STATUS = ALLOCATE_SURFACE_SECTION(TYPE, SEMIN, SEMAX, SURFACES(I))
        !
        IF (STATUS .NE. 0) THEN
            !        
            CALL PAR_QUIT('Error  :     > Surface allocation error.')
            !       
        ENDIF
        !
        DO J = 0, NSEL-1
            !
            ! For each element, read in the ID and 3D connectivity.
            READ(IO, *) ELEM_ID, ELEM_NODES_TMP
            !
            ! Shift the ID to maintain internal zero-indexing.
            ELEM_ID = ELEM_ID - 1
            ! 
            ! The node values read also need to be shifted and reordered.
            ! GMSH ORDER: 1 2 3 4 5 6
            ! FEPX ORDER: 6 4 2 5 3 1
            ELEM_NODES(0) = ELEM_NODES_TMP(5) - 1
            ELEM_NODES(1) = ELEM_NODES_TMP(2) - 1
            ELEM_NODES(2) = ELEM_NODES_TMP(4) - 1
            ELEM_NODES(3) = ELEM_NODES_TMP(1) - 1
            ELEM_NODES(4) = ELEM_NODES_TMP(3) - 1
            ELEM_NODES(5) = ELEM_NODES_TMP(0) - 1
            !            
            IF ((J <= SEMAX) .AND. (J>=SEMIN)) THEN
                !
                ! Create the 2D connectivity first.
                SURFACES(I)%ECONN(0,J) = 6 * ELEM_ID
                SURFACES(I)%ECONN(1,J) = 6 * ELEM_ID + 1
                SURFACES(I)%ECONN(2,J) = 6 * ELEM_ID + 2
                SURFACES(I)%ECONN(3,J) = 6 * ELEM_ID + 3
                SURFACES(I)%ECONN(4,J) = 6 * ELEM_ID + 4
                SURFACES(I)%ECONN(5,J) = 6 * ELEM_ID + 5
                !
                ! Store into nodes and ID into global connectivity?
                SURFACES(I)%CONN(:,J)=ELEM_NODES(:)
                SURFACES(I)%ELEM(J)=ELEM_ID
                !
            END IF
            !
        END DO
        !
        ! Create a 3d connectivity for the coordinate gather.
        ! 
        DO K = SEMIN, SEMAX
            !
            DO J = 0, TYPE1
                !
                J3 = 3 * J
                N3 = 3 * SURFACES(I)%CONN(J,K)
                SURFACES(I)%CONN3D(J3,K)   = N3
                SURFACES(I)%CONN3D(J3 + 1,K) = N3 + 1
                SURFACES(I)%CONN3D(J3 + 2,K) = N3 + 2
                !
            ENDDO
            !
        ENDDO
        !
        CALL PART_SCATTER_SETUP(0,  TYPE3, DOF_SUB1, DOF_SUP1,&
            & SEMIN, SEMAX, SURFACES(I)%CONN3D, SURFACES(I)%TR)
        EL_DOF_MIN = 6 * EL_SUB1
        EL_DOF_MAX = 6 * EL_SUP1 + 5
        CALL PART_SCATTER_SETUP(0,  5, EL_DOF_MIN, EL_DOF_MAX,&
            & SEMIN, SEMAX, SURFACES(I)%ECONN, SURFACES(I)%ETR)
        !
    ENDDO
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndFasets')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &fasets.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_FASETS
    !
    !===========================================================================
    !
    SUBROUTINE READ_NODEPARTITIONS(IO)
    !
    ! Read in node paritions and do not store - currently unused, but will be 
    ! integrated in the future for utilizing imprted SCOTCH partitions.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I: Generic loop index.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! NPARTS: Number of partitions read in from the mesh.
    ! IPART: Partition loop index to check which partition is being checked.
    ! TEMP: Temporary storage for current node ID - not stored.
    ! NODE_PARTS: Array storing the nodal paritions.
    ! TEMPVAL: Array of NP_SUB1/NP_SUP1 bounds for each partition read in.
    ! CHG_INDEX: Index in array immediately after a change in value occurs.
    ! IARRAY: Substring array for internal read parsing.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, I, NLINES, NPARTS, IPART, TEMP
    INTEGER, ALLOCATABLE :: NODE_PARTS(:), TEMPVAL(:,:), CHG_INDEX(:)
    CHARACTER(LEN=12)  :: IARRAY(16)
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$NodePartitions')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &node partitions.')
        !
    END IF
    !
    ! Set the flag if this section is being read in
    READ_NODE_PART = .TRUE.
    !
    ! Read in the number of node partitions in the section (and set NLINES).
    READ(IO, *) NLINES
    !
    ! Allocate the nodal partition array -> N(1,:) is 1-index ID, N(2,:) is part. 
    ALLOCATE(NODE_PARTS(1:NLINES))
    !
    ! Loop over the number of NLINES and do NOT store information.
    DO I = 1, NLINES
        !
        READ(IO, *) TEMP, NODE_PARTS(I)
        !
    END DO
    !
    NPARTS = MAXVAL(NODE_PARTS)
    ALLOCATE(TEMPVAL(1:NPARTS, 2))
    ALLOCATE(CHG_INDEX(NPARTS-1))
    IPART = 1
    !
    ! We can split ELEM_PARTS in 0-indexed EL_SUB1 and EL_SUP1 here.
    DO I = 1, SIZE(NODE_PARTS)
        !
        ! A change has been detected so do something
        IF (NODE_PARTS(I) .NE. IPART) THEN
            !
            CHG_INDEX(IPART) = I
            IF (IPART .NE. NPARTS) IPART = IPART + 1
            !
        END IF
        !
    END DO
    !
    ! DEBUG WRITE
!    WRITE(*,*) 'Number of node partitions: ', NPARTS
!    WRITE(*,*) 'Nodal step indices: ', CHG_INDEX
    !
    ! Build the bound array for NPARTS -> (1,:) is NP_SUB1, (2,:) is NP_SUP1.
    TEMPVAL = 0
    DO I = 1, NPARTS
        !
        ! Handle the edge cases differently
        IF (I .EQ. 1) THEN
            !
            TEMPVAL(1,1) = 1
            TEMPVAL(1,2) = CHG_INDEX(I) - 1
            !
        ELSE IF (I .EQ. NPARTS) THEN
            !
            TEMPVAL(I,1) = CHG_INDEX(I-1)
            TEMPVAL(I,2) = SIZE(NODE_PARTS)
            !
        ELSE
            !
            TEMPVAL(I,1) = CHG_INDEX(I-1)
            TEMPVAL(I,2) = CHG_INDEX(I) - 1
            !
        END IF
        !
    END DO
    !
    ! DEBUG WRITE
!    DO I = 1, NPARTS
!        WRITE(*,*) 'Partition ', I, ': ', TEMPVAL(I,:)
!    END DO
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndNodePartitions')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &node partitions.')
        !
    END IF
    !
    ! Do not need to maintain node partitions for now so deallocate.
    DEALLOCATE(NODE_PARTS)
    DEALLOCATE(CHG_INDEX)
    DEALLOCATE(TEMPVAL)
    !
    RETURN
    !
    END SUBROUTINE READ_NODEPARTITIONS
    !
    !===========================================================================
    !
    SUBROUTINE READ_PHYSICALNAMES(IO)
    !
    ! Read in physical names and do not store - will never be used as these
    ! are for gmsh internal use only and must be stored to successfully open
    ! the msh file in any version of gmsh.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I: Generic loop index.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, I, NLINES
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$PhysicalNames')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &physical names.')
        !
    END IF
    !
    ! Read in the number of physical names in the section (and set NLINES).
    READ(IO, *) NLINES
    !
    ! Loop over the number of NLINES and do NOT store information.
    DO I = 1, NLINES
        !
        READ(IO, *)
        !
    END DO
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndPhysicalNames')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &physical names.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_PHYSICALNAMES
    !
    !===========================================================================
    !
    SUBROUTINE READ_ELSETORIENTATIONS(IO)
    !
    ! Read elset orientation information from either the msh file or an external
    ! simulation.ori file. If we are reading in orientations from an external
    ! file we must deallocate the current storage array and start over.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I/J: Generic loop index.
    ! S: Status value that read in string is a valid option.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! NSPACE/NVALS: Storage values used for dividing a line into substrings.
    ! ELSET: Grain ID - 1-indexed.
    ! DELIM_POS: Line position where the ':' delimiter is found.
    ! IDEAL_DIM/IDEAL_DIM1: Array widths for storing orientations (3 or 4).
    ! PARM_STRING: String describing input orientation parameterization.
    ! CONV_STRING: String describing input orientation convention.
    ! ORI_STRING: Temp string used to split the 'descriptor:convention' line.
    ! IARRAY: Substring array for internal read parsing.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, I, J, S
    INTEGER :: NLINES, NSPACE, NVALS, ELSET, DELIM_POS
    INTEGER, PARAMETER :: IDEAL_DIM  = 4
    INTEGER, PARAMETER :: IDEAL_DIM1 = IDEAL_DIM - 1
    CHARACTER(LEN=50)  :: PARM_STRING, CONV_STRING, ORI_STRING
    CHARACTER(LEN=32)  :: IARRAY(16)
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! If per-element orientations have been read-in previously then skip this.
    IF (ELEMENT_ORIS .EQV. .TRUE.) THEN
        !
        ! Read the first line and confirm it is the correct record.
        READ(IO, '(A)', IOSTAT = IERR) LINE
        !
        IF ((IERR .LT. 0) .OR. (LINE .NE. '$ElsetOrientations')) THEN
            !
            CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
                &elset orientations.')
            !
        END IF
        !
        ! Extract the number of lines to skip
        READ(IO, '(A)', IOSTAT = IERR) LINE
        NSPACE = COUNT(TRANSFER(LINE, 'A', LEN_TRIM(LINE)) == " ")
        NVALS = NSPACE + 1
        IARRAY = ""
        !
        ! Internal read of the line into the primary substring arrays.
        READ(LINE, *) IARRAY(1:NVALS)
        READ(IARRAY(1), *) NLINES
        !
        ! Loop over NLINES and do not store read in values.
        DO I = 0, (NLINES - 1)
            !
            READ(IO, *) LINE
            !
        END DO
        !
        ! Read the end of section footer.
        READ(IO, '(A)', IOSTAT = IERR) LINE
        !
        IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndElsetOrientations')) THEN
            !
            CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
                &elset orientations.')
            !
        END IF
        !
        RETURN
        !
    END IF
    !
    ! Check if the array is already allocated from reading the mesh file.
    IF (ALLOCATED(GRAIN_ORIENTATION) .EQV. .TRUE.) THEN
        !
        DEALLOCATE(GRAIN_ORIENTATION)
        !
    END IF
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$ElsetOrientations')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &elset orientations.')
        !
    END IF
    !
    ! Read in the section information string and prepare to parse.
    !
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    ! Trim the full string into two substring arrays split by the space delim.
    NSPACE = COUNT(TRANSFER(LINE, 'A', LEN_TRIM(LINE)) == " ")
    NVALS = NSPACE + 1
    IARRAY = ""
    !
    ! Internal read of the line into the primary substring arrays.
    READ(LINE, *) IARRAY(1:NVALS)
    !
    ! Internal read of the primary substring arrays to get the number of 
    ! lines to set to parse and prepare to parse the orientation description.
    READ(IARRAY(1), *) NLINES
    READ(IARRAY(2), '(A)') ORI_STRING
    !
    ! Find the position of the ":" character in ORI_STRING as FORTRAN
    ! refuses to allow the READ() function to handle this automatically.
    DELIM_POS = INDEX(ORI_STRING, ":")
    !
    ! Internal read to store the orientation parameterization and convention.
    READ(ORI_STRING(1:DELIM_POS-1), *) PARM_STRING
    READ(ORI_STRING(DELIM_POS+1:LEN_TRIM(ORI_STRING)), *) CONV_STRING
    !
    ! Set a per-element logical to false.
    ELEMENT_ORIS = .FALSE.
    !
    ! Parse the parameterization and confirm it is a valid option.
    !
    S = 0
    !
    READ(PARM_STRING, '(A)') &
        & ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION
    !
    ! Confirm that the rest of the string is a valid input.
    IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'axis-angle') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'euler-bunge') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'euler-kocks') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'rodrigues') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'quaternion') THEN
        !
        S = 0
        !
    ELSE ! Orientation parameterization defined is incorrect.
        !
        S = 1
        !
    END IF
    !
    IF (S .EQ. 1) THEN
        !
        CALL PAR_QUIT("Error  :     > ORIENTATION_PARAMETERIZATION &
            & contains an error or unexpected input type.")
        !
    END IF
    !
    ! Parse the convention and confirm it is a valid option.
    !
    S = 0
    !
    READ(CONV_STRING, '(A)') &
        & ORIENTATION_OPTIONS%ORIENTATION_CONVENTION
    !
    ! Confirm that the rest of the string is a valid input.
    IF (ORIENTATION_OPTIONS%ORIENTATION_CONVENTION &
        & .EQ. 'active') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_CONVENTION &
        & .EQ. 'passive') THEN
        !
        S = 0
        !
    ELSE ! Orientation convention defined is incorrect.
        !
        S = 1
        !
    END IF
    !
    IF (S .EQ. 1) THEN
        !
        CALL PAR_QUIT("Error  :     > ORIENTATION_CONVENTION &
            & contains an error or unexpected input type.")
        !
    END IF
    !
    ! Check orientation parameterization to set column number.
    IF ( (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'axis-angle') .OR. &
        & (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'quaternion') ) THEN
        !
        ALLOCATE(GRAIN_ORIENTATION(0:IDEAL_DIM1, 0:NLINES - 1))
        !
    ELSE ! All other parameterization have 3 values per orientation.
        !
        ALLOCATE(GRAIN_ORIENTATION(0:IDEAL_DIM1-1, 0:NLINES - 1))
        !
    END IF
    !
    IF (SIZE(GRAIN_ORIENTATION, 1) .EQ. 4) THEN ! Only quaternion or axis-angle.
        !
        DO I = 0, (NLINES - 1)
            !
            READ(IO, *) ELSET, (GRAIN_ORIENTATION(J,I), J = 0, IDEAL_DIM1)
            !
        END DO
        !
    ELSE IF (SIZE(GRAIN_ORIENTATION,1) .EQ. 3) THEN
        !
        DO I = 0, (NLINES - 1)
            !
            READ(IO, *) ELSET, (GRAIN_ORIENTATION(J,I), J = 0, IDEAL_DIM1-1)
            !
        END DO
        !
    ELSE
        !
        CALL PAR_QUIT("Error  :     > GRAIN_ORIENTATION array not initialized&
            & correctly.")
        !
    END IF
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndElsetOrientations')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &elset orientations.')
        !
    END IF
    !
    ! DEBUG PRINT
    ! WRITE(*, '(3(E17.7,1X))') GRAIN_ORIENTATION
    ! WRITE(*,*) SIZE(GRAIN_ORIENTATION,1), SIZE(GRAIN_ORIENTATION,2)
    !
    RETURN
    !
    END SUBROUTINE READ_ELSETORIENTATIONS
    !
    !===========================================================================
    !
    SUBROUTINE READ_ELEMENTORIENTATIONS(IO)
    !
    ! Read elem orientation information from either the msh file or an external
    ! simulation.ori file. If we are reading in orientations from an external
    ! file we must deallocate the current storage array and start over.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I/J: Generic loop index.
    ! S: Status value that read in string is a valid option.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! NSPACE/NVALS: Storage values used for dividing a line into substrings.
    ! ELSET: Grain ID - 1-indexed.
    ! DELIM_POS: Line position where the ':' delimiter is found.
    ! IDEAL_DIM/IDEAL_DIM1: Array widths for storing orientations (3 or 4).
    ! PARM_STRING: String describing input orientation parameterization.
    ! CONV_STRING: String describing input orientation convention.
    ! ORI_STRING: Temp string used to split the 'descriptor:convention' line.
    ! IARRAY: Substring array for internal read parsing.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: IERR, I, J, S
    INTEGER :: NLINES, NSPACE, NVALS, ELSET, DELIM_POS
    INTEGER, PARAMETER :: IDEAL_DIM  = 4
    INTEGER, PARAMETER :: IDEAL_DIM1 = IDEAL_DIM - 1
    CHARACTER(LEN=50)  :: PARM_STRING, CONV_STRING, ORI_STRING
    CHARACTER(LEN=32)  :: IARRAY(16)
    CHARACTER(LEN=256) :: LINE
    !
    !---------------------------------------------------------------------------
    !
    ! Check if the array is already allocated from reading the mesh file.
    IF (ALLOCATED(GRAIN_ORIENTATION) .EQV. .TRUE.) THEN
        !
        DEALLOCATE(GRAIN_ORIENTATION)
        !
    END IF
    !
    ! Read the first line and confirm it is the correct record.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$ElementOrientations')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &element orientations.')
        !
    END IF  
    !
    ! Read in the section information string and prepare to parse.
    !
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    ! Trim the full string into two substring arrays split by the space delim.
    NSPACE = COUNT(TRANSFER(LINE, 'A', LEN_TRIM(LINE)) == " ")
    NVALS = NSPACE + 1
    IARRAY = ""
    !
    ! Internal read of the line into the primary substring arrays.
    READ(LINE, *) IARRAY(1:NVALS)
    !
    ! Internal read of the primary substring arrays to get the number of 
    ! lines to set to parse and prepare to parse the orientation description.
    READ(IARRAY(1), *) NLINES
    READ(IARRAY(2), '(A)') ORI_STRING
    !
    ! Find the position of the ":" character in ORI_STRING as FORTRAN
    ! refuses to allow the READ() function to handle this automatically.
    DELIM_POS = INDEX(ORI_STRING, ":")
    !
    ! Internal read to store the orientation parameterization and convention.
    READ(ORI_STRING(1:DELIM_POS-1), *) PARM_STRING
    READ(ORI_STRING(DELIM_POS+1:LEN_TRIM(ORI_STRING)), *) CONV_STRING
    !
    ! Set a per-element logical to true.
    ELEMENT_ORIS = .TRUE.
    !
    ! Parse the parameterization and confirm it is a valid option.
    !
    S = 0
    !
    READ(PARM_STRING, '(A)') &
        & ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION
    !
    ! Confirm that the rest of the string is a valid input.
    IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'axis-angle') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'euler-bunge') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'euler-kocks') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'rodrigues') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'quaternion') THEN
        !
        S = 0
        !
    ELSE ! Orientation parameterization defined is incorrect.
        !
        S = 1
        !
    END IF
    !
    IF (S .EQ. 1) THEN
        !
        CALL PAR_QUIT("Error  :     > ORIENTATION_PARAMETERIZATION &
            & contains an error or unexpected input type.")
        !
    END IF
    !
    ! Parse the convention and confirm it is a valid option.
    !
    S = 0
    !
    READ(CONV_STRING, '(A)') &
        & ORIENTATION_OPTIONS%ORIENTATION_CONVENTION
    !
    ! Confirm that the rest of the string is a valid input.
    IF (ORIENTATION_OPTIONS%ORIENTATION_CONVENTION &
        & .EQ. 'active') THEN
        !
        S = 0
        !
    ELSE IF (ORIENTATION_OPTIONS%ORIENTATION_CONVENTION &
        & .EQ. 'passive') THEN
        !
        S = 0
        !
    ELSE ! Orientation convention defined is incorrect.
        !
        S = 1
        !
    END IF
    !
    IF (S .EQ. 1) THEN
        !
        CALL PAR_QUIT("Error  :     > ORIENTATION_CONVENTION &
            & contains an error or unexpected input type.")
        !
    END IF
    !
    ! Check orientation parameterization to set column number.
    IF ( (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'axis-angle') .OR. &
        & (ORIENTATION_OPTIONS%ORIENTATION_PARAMETERIZATION &
        & .EQ. 'quaternion') ) THEN
        !
        ALLOCATE(GRAIN_ORIENTATION(0:IDEAL_DIM1, 0:NLINES - 1))
        !
    ELSE ! All other parameterization have 3 values per orientation.
        !
        ALLOCATE(GRAIN_ORIENTATION(0:IDEAL_DIM1-1, 0:NLINES - 1))
        !
    END IF
    !
    IF (SIZE(GRAIN_ORIENTATION, 1) .EQ. 4) THEN ! Only quaternion or axis-angle.
        !
        DO I = 0, (NLINES - 1)
            !
            READ(IO, *) ELSET, (GRAIN_ORIENTATION(J,I), J = 0, IDEAL_DIM1)
            !
        END DO
        !
    ELSE IF (SIZE(GRAIN_ORIENTATION,1) .EQ. 3) THEN
        !
        DO I = 0, (NLINES - 1)
            !
            READ(IO, *) ELSET, (GRAIN_ORIENTATION(J,I), J = 0, IDEAL_DIM1-1)
            !
        END DO
        !
    ELSE
        !
        CALL PAR_QUIT("Error  :     > GRAIN_ORIENTATION array not initialized&
            & correctly.")
        !
    END IF
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndElementOrientations')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &orientations.')
        !
    END IF
    !
    ! DEBUG PRINT
    ! WRITE(*, '(3(E17.7,1X))') GRAIN_ORIENTATION
    ! WRITE(*,*) SIZE(GRAIN_ORIENTATION,1), SIZE(GRAIN_ORIENTATION,2)
    !
    RETURN
    !
    END SUBROUTINE READ_ELEMENTORIENTATIONS
    !
    !===========================================================================
    !
    SUBROUTINE READ_GROUPS(IO)
    !
    ! Read grain and phase information from either the .msh file or an external
    ! simulation.phase file. If we are reading in phases from an external
    ! file we must deallocate the current storage array and start over.
    !
    !---------------------------------------------------------------------------
    !
    ! Arguments:
    ! IO: Input unit for .msh file.
    !
    INTEGER, INTENT(IN) :: IO
    !
    ! Locals:
    ! IERR: Value that confirms if a READ() fails.
    ! I: Generic loop index.
    ! NLINES: Read-in value for the number of lines that will be parsed.
    ! ELSET: Grain ID - 1-indexed.
    ! ELSET_PHASE: Storage array to hold which phase is assigned to which grain.
    ! LINE: Input line on current record to be parsed.
    !
    INTEGER :: I, IERR, NLINES, ELSET
    INTEGER, ALLOCATABLE :: ELSET_PHASE(:)
    CHARACTER(LEN=256)   :: LINE
    !
    ! Notes:
    ! If the simulation is single-phase then we don't input any information
    ! on the grain phases anymore. Therefore we need to attempt to read in
    ! the next record on the mesh unit and see if it fails by EOF. This only 
    ! works as $Groups will always be the last section in the mesh file iff 
    ! it is present and generated via Neper. 
    ! 
    ! Todo:
    ! It seems to be possible to define a number of grain/phase assignments 
    ! that is out of bounds of the previously allocated array since we don't 
    ! explicity enforce this.
    !----------------------------------------------------------------------
    !
    ! Attempt to read the next line (may be EOF, but unknown)
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF (IERR .LT. 0) THEN
        !
        IF (MYID .EQ. 0) THEN
            !
            WRITE(DFLT_U, '(A)') 'Info   :   [i] No phase assignments found &
                &in `simulation.msh`.'
            !
        END IF
        !
        RETURN
        !
    ELSE IF (IERR .EQ. 0) THEN
        !
        ! Check that the line read in previously was $Groups.
        IF (LINE .NE. '$Groups') THEN
            !
            CALL PAR_QUIT('Error  :     > Phase assignment section header &
                &should be `$Groups`.')
            !
        END IF
        !
        ! Phase assignments are per elset so we need to do some mapping with
        ! temporary arrays.
        READ(IO, *) ! This skips the `elset` line
        READ(IO, *) NLINES
        !
        ! Check if the array is already allocated from the mesh read in
        IF (ALLOCATED(ELSET_PHASE) .EQV. .TRUE.) THEN
            !
            DEALLOCATE(ELSET_PHASE)
            !
        END IF
        !
        ! Allocate and fill an array with elset-based `grain/phase` pairs
        ! (1,:) is phase assign to grain ID (row index)
        ALLOCATE(ELSET_PHASE(1:NLINES))
        !
        DO I = 1, NLINES
            !
            READ(IO, *) ELSET, ELSET_PHASE(I)
            !
        END DO
        !
        ! Now, loop over UESUB and reassign phase accordingly
        DO I = 0, NUMELM-1
            !
            ! Assign phase per element based on elset value of UESUB per element
            ! UESUB(I) returns an elset value that is 1-indexed that will grab
            ! the corresponding elset row in ELSET_PHASE and, in turn, the phase
            PHASE(I) = ELSET_PHASE(UESUB(I))
            !
        END DO
        !
        ! Now that PHASE is built there is a chance that the user supplied
        ! an incorrect `number_of_phases` parameter for what the .msh has
        IF (MAXVAL(ELSET_PHASE) .NE. CRYS_OPTIONS%NUMBER_OF_PHASES) THEN
            !
            CALL PAR_QUIT('Error  :     > Number of phases in mesh does not&
                & match the config file.')
            !
        END IF
        !
    ELSE
        ! If something goes wrong terminate the simulation. 
        CALL PAR_QUIT('Error  :     > Phase assignment parsing failed.')
        !
    END IF
    !
    !
    ! Read the end of section footer.
    READ(IO, '(A)', IOSTAT = IERR) LINE
    !
    IF ((IERR .LT. 0) .OR. (LINE .NE. '$EndGroups')) THEN
        !
        CALL PAR_QUIT('Error  :     > Parse error attempting to read in &
            &groups.')
        !
    END IF
    !
    RETURN
    !
    END SUBROUTINE READ_GROUPS
    !
END MODULE READ_INPUT_MOD
